import os
import numpy

from beers.molecule_packet import MoleculePacket
from beers.molecule import Molecule
from beers.sample import Sample
from beers.utilities.general_utils import GeneralUtils

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
            else:
                cigar = "200S" + cigar # Relative to + strand, the A's are going on the 5' end


        return Molecule(Molecule.new_id(), sequence, start=1, cigar=f"{len(sequence)}M",
                            source_start=starts[0], source_cigar=cigar, source_strand=strand,
                            transcript_id=transcript_id)

    def make_packet(self, id="packet0", N=10_000):
        molecules = [self.make_molecule() for i in range(N)]
        return MoleculePacket(id, self.sample, molecules)

    def make_molecule_file(self, filepath, N=10_000):
        """
        Write out molecules to a tab-separated file

        Note: we write out a molecules start and cigar relative to the appropriate
        custom genome, either _1 or _2 as per the transcript id
        """
        with open(filepath, "w") as molecule_file:
            header = "#transcript_id\tstart\tcigar\tstrand\tsequence\n"
            molecule_file.write(header)
            for i in range(N):
                molecule = self.make_molecule()
                # NOTE: Not outputing the molecules start or cigar string since those are relative to parent
                #       which in this case is always trivial (start=1, cigar=###M) since the molecule is new
                line = "\t".join([molecule.transcript_id,
                                  str(molecule.source_start),
                                  molecule.source_cigar,
                                  molecule.source_strand,
                                  molecule.sequence]
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
