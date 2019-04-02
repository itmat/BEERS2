import os
import collections
import numpy

from beers_utils.molecule_packet import MoleculePacket
from beers_utils.molecule import Molecule
from beers_utils.sample import Sample
from beers_utils.general_utils import GeneralUtils

class MoleculeMaker:
    """
    MoleculeMaker generates molecules based off of gene, intron, and allelic quantification files
    as well as customized genomic sequence and annotation
    """
    def __init__(self, sample, sample_directory):
        """
        :param sample: sample object to work on
        :param sample_directory: path to the directory that contains the sample's file (i.e. the CAMPAREE data directory for the sample)
        """
        self.sample = sample

        intron_quant_file_path = os.path.join(sample_directory, "intron_quantifications.txt")
        self.transcript_intron_quants, self.intron_quants = self.load_intron_quants(intron_quant_file_path)

        gene_quant_file_path = os.path.join(sample_directory, "gene_quantifications.txt")
        self.genes, self.gene_quants = self.load_gene_quants(gene_quant_file_path)
        self.gene_probabilities = self.gene_quants / numpy.sum(self.gene_quants)

        self.transcriptomes = [self.load_transcriptome(os.path.join(sample_directory, f"transcriptome_{i}.fa"))
                                    for i in [1,2]]

        self.annotations = [self.load_annotation(os.path.join(sample_directory, f"updated_annotation_{i}.txt"))
                                    for i in [1,2]]

        self.genomes = [self.load_genome(os.path.join(sample_directory, f"custom_genome_{i}.fa"))
                                    for i in [1,2]]

        indel_data = [self.load_indels(os.path.join(sample_directory, f"custom_genome_indels_{i}.txt"))
                                    for i in [1,2]]
        self.indels = [indels for indels, offset_data in indel_data]
        self.offset_data = [offset_data for indels, offset_data in indel_data]

        isoform_quant_file_path = os.path.join(sample_directory, "isoform_psi_value_quantifications.txt")
        self.isoform_quants = self.load_isoform_quants(isoform_quant_file_path)

        allelic_quant_file_path = os.path.join(sample_directory, "allelic_imbalance_quantifications.txt")
        self.allelic_quant = self.load_allelic_quants(allelic_quant_file_path)

    def load_annotation(self, file_path):
        transcripts = dict()
        with open(file_path) as annotation_file:
            for line in annotation_file:
                if line.startswith("#"):
                    continue # Comment/header line

                chrom, strand, tx_start, tx_end, exon_count, exon_starts, exon_ends, transcript_id, gene_id, gene_sybmol, *other \
                        = line.split("\t")

                transcripts[transcript_id] = (chrom, strand, int(tx_start), int(tx_end),
                                                [int(start) for start in exon_starts.split(",")],
                                                [int(end) for end in exon_ends.split(",")])
        return transcripts


    def load_transcriptome(self, file_path):
        """
        Read in a fasta file and load is a dictionary id -> sequence
        assumed one-line for the whole contig
        """
        transcripts = dict()
        with open(file_path) as transcriptome_file:
            while True:
                line = transcriptome_file.readline()
                if not line:
                    break

                assert line[0] == ">"
                transcript_id, chrom, region = line[1:].strip().split(":")
                sequence = transcriptome_file.readline().strip()
                transcripts[transcript_id] = sequence
        return transcripts

    def load_genome(self, file_path):
        """
        Read in a fasta file and load is a dictionary id -> sequence
        """
        genome = dict()
        with open(file_path) as transcriptome_file:
            while True:
                line = transcriptome_file.readline()
                if not line:
                    break

                assert line[0] == ">"
                contig = line[1:].strip()
                sequence = transcriptome_file.readline().strip()
                genome[contig] = sequence
        return genome

    def load_indels(self, file_path):
        """
        Read in the file of indel locations for a given custom genome

        Store it as a dictionary chrom -> (indel_starts, indel_data)
        where indel_starts is a numpy array of start locations of the indels
        (i.e. 1 based coordinates of the base in the custom genome where the insertion/deletion occurs immediately after)
        and indel_data is a list of tuples (indel_start, indel_type, indel_length)
        where indel_type is 'I' or 'D'

        The indel file is tab-separated with format "chrom:start type length"  and looks like the following:
        1:4897762       I       2
        1:7172141       I       2
        1:7172378       D       1

        Assumption is that the file is sorted by start and no indels overlap

        Moreover, return a dictionary chrom -> (offset_starts, offset_values)
        where offset_starts is as indel_starts, a sorted numpy array with indicating the positions
        where the offset (i.e. reference_genome_position - custom_genome_position values) change
        and offset_values is a list in the same order indicating these values.
        Note that the offset value is to be used for all bases AFTER the offset position, not on that base
        """
        indels = collections.defaultdict(lambda : (collections.deque(), collections.deque()))
        offset_values = collections.defaultdict(lambda : collections.deque())
        current_offsets = collections.defaultdict(lambda : 0)
        with open(file_path) as indel_file:
            for line in indel_file:
                loc, indel_type, length = line.split('\t')
                chrom, start = loc.split(':')

                # indel_start needs to be relative to the custom genome but the indel file has starts
                # that are relative to the reference genome
                indel_start = int(start) + current_offsets[chrom]
                indel_length = int(length)

                starts, data = indels[chrom]
                starts.append(indel_start)
                data.append((indel_start, indel_type, indel_length))

                if indel_type == "I":
                    current_offsets[chrom] += indel_length
                elif indel_type == "D":
                    current_offsets[chrom] -= indel_length

                offset_values[chrom].append(current_offsets[chrom])

        # We make these defaultdicts in case there are chromosomes without any indels
        # in which case we don't see them in this file, but we don't want to crash on them
        result = collections.defaultdict(lambda: (numpy.array([]), list()))
        offset_data = collections.defaultdict(lambda: (numpy.array([]), list()))
        for key, (starts, data) in indels.items():
            result[key] = (numpy.array(starts), list(data))
            offset_data[key] = (numpy.array(starts), offset_values[key])

        return result, offset_data

    def load_intron_quants(self, file_path):
        """
        Load an intron quantification file as two dictionaries,
        (transcript ID -> sum FPK of all introns in transcript) and
        (transcript ID -> list of FPKs of each intron in transcript)
        """
        transcript_intron_quants = dict() # Dictionary transcript -> FPK for all introns in the transcript, combined
        intron_quants = dict() # Dictioanry transcript -> array of FPKs for each intron in the transcript

        with open(file_path) as intron_quants_file:
            for line in intron_quants_file:

                if line.startswith("#"):
                    continue # Comment/header line

                gene, transcript, chrom, strand, transcript_intron_reads_FPK, intron_reads_FPK = line.strip().split("\t")

                transcript_intron_quants[transcript] = float(transcript_intron_reads_FPK)
                intron_quants[transcript] = [float(quant) for quant in intron_reads_FPK.split(",")]

        return transcript_intron_quants, intron_quants

    def load_gene_quants(self, file_path):
        """
        Read in a gene quantification file as two lists of gene IDs and of their read quantifications
        """
        genes = []
        gene_quants = []

        with open(file_path) as gene_quant_file:
            for line in gene_quant_file:
                if line.startswith("#"):
                    continue # Comment/header line

                gene, quant = line.strip().split("\t")

                genes.append(gene)
                gene_quants.append(float(quant))

        return genes, numpy.array(gene_quants)

    def load_isoform_quants(self, file_path):
        """
        Reads an isoform quant file into a dictionary gene -> (list of transcript IDs, list of psi values)
        """
        isoform_quants = dict()

        with open(file_path) as isoform_quant_file:
            for line in isoform_quant_file:

                if line.startswith("#"):
                    continue # Comment/header line

                gene, entries = line.strip().split("\t")
                isoforms = [entry.split(":") for entry in entries.split(",")]
                isoforms = [(isoform, float(psi)) for isoform, psi in isoforms]
                isoform_list = [isoform for isoform, psi in isoforms]
                psi_list = [psi for isoform, psi in isoforms]

                isoform_quants[gene] = (isoform_list, psi_list)

        return isoform_quants

    def load_allelic_quants(self, file_path):
        """
        Reads allelic quantification file into a dictionary: gene_id -> (allele 1 probability, allele 2 probability)
        """

        allelic_quant = dict()

        with open(file_path) as allele_quant_file:
            for line in allele_quant_file:
                if line.startswith("#"):
                    continue # Comment/header line

                gene, allele1, allele2 = line.split("\t")
                allele1 = float(allele1)
                allele2 = float(allele2)

                allelic_quant[gene] = (allele1, allele2)
        return allelic_quant

    def make_molecule(self):
        # Pick random gene
        gene_index = numpy.random.choice(len(self.genes), p=self.gene_probabilities)
        gene = self.genes[gene_index]
        gene_quant = self.gene_quants[gene_index]

        # Pick random transcript in gene
        transcripts, psis = self.isoform_quants[gene]
        transcript = numpy.random.choice(transcripts, p=psis)

        # Pick random allele based on the gene's allelic distribution
        allele_number = numpy.random.choice([1,2], p=self.allelic_quant[gene])

        # Read in annotation for the chosen transcript
        chrom,strand,tx_start,tx_end,starts,ends= self.annotations[allele_number - 1][transcript]

        # Determine if pre_mRNA or mature mRNA
        intron_quant = self.transcript_intron_quants[transcript]
        # TODO: check that this gives the appropriate fraction as pre_mRNA
        #       previously was using intron_quant / (intron_quant + gene_quant)
        #       but if assuming everything is either full pre_mRNA or mature mRNA then this should be
        #       the right fraction, which could happen to be greater than one (!)
        try:
            fraction_pre_mRNA = min(intron_quant / (gene_quant), 1)
        except ZeroDivisionError:
            # Should not ever get here since a gene with 0 gene_quant should have 0 chance of being chosen
            # however, if we do, we will just always give pre_mRNA
            fraction_pre_mRNA = 1.0

        pre_mRNA = numpy.random.uniform() < fraction_pre_mRNA
        if pre_mRNA:
            # If chosen to be pre_mRNA, overwrite the usual exon starts/ends with a single, big "exon"
            starts = [tx_start]
            ends = [tx_end]

        # Find cigar string relative to the reference genome (i.e. custom_genome_1 or custom_genome_2)
        gaps = [next_start - last_end - 1 for next_start,last_end in zip(starts[1:],ends[:-1])]
        cigar = ''.join( f"{end - start + 1}M{gap}N" for start,end,gap in zip(starts[:-1],ends[:-1],gaps)) \
                    + f"{ends[-1] - starts[-1] + 1}M"

        ref_cigar =  ''.join( f"{self.get_reference_cigar(start, end, chrom, allele_number)}{gap}N"
                            for start,end,gap in zip(starts[:-1],ends[:-1],gaps)) \
                        + self.get_reference_cigar(starts[-1], ends[-1], chrom, allele_number)
        ref_start = self.convert_genome_position_to_reference(starts[0], chrom, allele_number)

        transcript_id = f"{transcript}_{allele_number}{'_pre_mRNA' if pre_mRNA else ''}"


        # Build the actual sequence
        chrom_sequence = self.genomes[allele_number - 1][chrom]
        sequence = ''.join( chrom_sequence[start-1:end] for start,end in zip(starts, ends) )

        if strand == '-':
            # We always give the sequence from 5' to 3' end of the RNA molecule
            # so reverse complement this
            sequence = GeneralUtils.create_complement_strand(sequence)
            # NOTE: cigar string stays the same since that is relative to the + strand

        # TODO: for now, everything gets polyA but maybe shouldn't
        polyA_tail = True
        if polyA_tail:
            # TODO: polyA tails should vary in length
            # Add polyA tail to 3' end
            sequence = sequence + "A"*200
            # Soft-clip the polyA tail at the end since it shouldn't align
            if strand == "+":
                cigar = cigar + "200S"
                ref_cigar = ref_cigar + "200S"
            else:
                cigar = "200S" + cigar # Relative to + strand, the A's are going on the 5' end
                ref_cigar = "200S" + ref_cigar


        return sequence, starts[0], cigar, ref_start, ref_cigar, strand, chrom, transcript_id

    def get_reference_cigar(self, start, end, chrom, allele):
        """
        Returns the cigar string for the part of the custom chromosome on the segment from  start to end (inclusive, one based)
        relative to the reference genome
        """
        indel_starts, indel_data = self.indels[allele-1][chrom]

        cigar_components = []
        start_of_remaining = start

        # Skip to the first relevant indel, then continue until they go past this region
        index_idx = numpy.searchsorted(indel_starts, start)
        for indel in indel_data[index_idx:]:
            indel_start, indel_type, indel_length = indel
            if indel_start >= end:
                # The deletion or insertion begins immediately AFTER indel_start
                # so we're done, we've gone through all indels in our exon
                break

            if indel_start > start_of_remaining:
                # Match a region without indels
                cigar_components.append(f"{indel_start - start_of_remaining + 1}M")

            if indel_type == "I":
                # Custom genome has an insertion right after indel_start
                # Our region might end before the insertion does, so cap the length of the insertion
                length = min(end, indel_start + indel_length) - indel_start
                cigar_components.append(f"{length}I")
                # Increment to after the indel, which includes the inserted bases
                start_of_remaining = indel_start + indel_length + 1
            elif indel_type == "D":
                # Custom genome has a deletion right after indel_start
                cigar_components.append(f"{indel_length}D")
                # Advance to just after the deletion, but don't increment by the length of the deletion
                # since the deleted bases don't appear in our custom genome segment
                start_of_remaining = indel_start + 1

        # The last chunk of matches after the last indel
        if end >= start_of_remaining:
            cigar_components.append(f"{end - start_of_remaining + 1}M")

        return ''.join(cigar_components)

    def convert_genome_position_to_reference(self, position, chrom, allele):
        """
        Convert (1-indexed) position into the current (custom) genome into a position
        relative to the reference genome
        """
        # offset_positions is a sorted list of all offsets in this allele+chromosome
        # so find where our position fits in, so offset_index -1 is the last offset strictly before position
        offset_starts, offset_values = self.offset_data[allele-1][chrom]
        offset_index = numpy.searchsorted(offset_starts, position)

        if offset_index == 0:
            # If ours is before all offsets, then we need no offset
            return position
        else:
            offset = offset_values[offset_index-1]
            return position - offset

    def make_packet(self, id="packet0", N=10_000):
        molecules = []
        for i in range(N):
            sequence, start, cigar, strand, ref_start, ref_cigar, chrom, transcript_id = self.make_molecule()
            # TODO: use ref_start, ref_cigar when making Molecules
            #       this requires that the Molecule class include these cigar strings
            mol = Molecule(Molecule.new_id(), sequence, start=1, cigar=f"{len(sequence)}M",
                                source_start=start, source_cigar=cigar, source_strand=strand,
                                transcript_id=transcript_id)
            molecules.append(mol)
        return MoleculePacket(id, self.sample, molecules)

    def make_molecule_file(self, filepath, N=10_000):
        """
        Write out molecules to a tab-separated file

        Note: we write out a molecules start and cigar relative to the appropriate
        custom genome, either _1 or _2 as per the transcript id
        """
        with open(filepath, "w") as molecule_file:
            header = "#transcript_id\tchrom\tstart\tcigar\tref_start\tref_cigar\tstrand\tsequence\n"
            molecule_file.write(header)
            for i in range(N):
                sequence, start, cigar, ref_start, ref_cigar, strand, chrom, transcript_id = self.make_molecule()
                # NOTE: Not outputing the molecules start or cigar string since those are relative to parent
                #       which in this case is always trivial (start=1, cigar=###M) since the molecule is new
                line = "\t".join([transcript_id,
                                  chrom,
                                  str(start),
                                  cigar,
                                  str(ref_start),
                                  ref_cigar,
                                  strand,
                                  sequence]
                                  ) + "\n"

                molecule_file.write(line)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Molecule Maker")
    parser.add_argument("sample", help="sample object (as output by repr() in single quotes)")
    parser.add_argument("sample_directory", help="directory in which sample results go")

    args = parser.parse_args()

    # Remove single quotes and eval the sample
    sample_object = eval(args.sample.strip("'"))

    print("Starting to make molecules...")
    mm = MoleculeMaker(sample_object, args.sample_directory)
    packet = mm.make_packet()

    import pickle
    print("Writing packet out to molecule_packet.pickle")
    with open(os.path.join(args.sample_directory, "molecule_packet.pickle")) as out_file:
        pickle.dump(packet, out_file)

''' Example:
python molecule_maker.py 'Sample(1, "Sample 1", "/", ["ACGT", "CGCA"], gender="male")' data/_run1/sample1/
'''
