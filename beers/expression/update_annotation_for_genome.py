import argparse
import sys
import os
import collections
import bisect
from timeit import default_timer as timer
from beers.utils import Utils

class UpdateAnnotationForGenome:
    """Updates a gene annotation's coordinates to account for insertions &
    deletions (indels) introduced by GenomeFilesPreparation when it creates
    variant genomes. Note, this is designed to update annotation for a single
    variant genome at a time.

    Parameters
    ----------
    genome_indel_filename : string
        Path to file containing list of indel locations generated by the
        GenomeFilesPreparation. This file has no header and contains three,
        tab-delimited colums (example):
            1:134937	D	2
            1:138813	I	1
        First column = Chromosome and coordinate of indel in the original
                       reference genome. Note, coordinate is zero-based.
        Second column = "D" if variant is a deletion, "I" if variant is an
                        insertion.
        Third column = Length in bases of indel.
    input_annot_filename : string
        Path to file containing gene/transcript annotations using coordinates
        to the original reference genome. This file should have 11, tab-delimited
        columns and includes a header (example):
            chrom  strand  txStart  txEnd exonCount  exonStarts         exonEnds           transcriptID     geneID           geneSymbol  biotype
            1      +       11869    14409 3          11869,12613,13221  12227,12721,14409  ENST00000456328  ENSG00000223972  DDX11L1     pseudogene
        An annotation file with this format can be generated from a GTF file
        using the convert_gtf_to_annot_file_format() function in the Utils
        package. A template for this annotation format is available in the
        class variable Utils.annot_output_format.
    updated_annot_filename : string
        Path to the output file containing gene/transcript annotations with
        coordinates updated to match the variant genome.
    log_filename : string
        Path to the log file.

    """

    def __init__(self, genome_indel_filename, input_annot_filename,
                 updated_annot_filename, log_filename):
        """Short summary.

        Parameters
        ----------
        genome_indel_filename : string
            Full path to indel file generated by GenomeFilesPreparation.
        input_annot_filename : string
            Full path to annotation file with coordinates for reference genome.
        updated_annot_filename : string
            Full path to output annotation file with coordinates updated for
            variant genome.
        log_filename : string
            Full path to log file.

        """
        self.genome_indel_filename = genome_indel_filename
        self.input_annot_filename = input_annot_filename
        self.updated_annot_filename = updated_annot_filename
        self.log_filename = log_filename

        #TODO add checks to make sure given files exist. Or shift the test below
        #     to another method called after another files are tested (maybe the
        #     main?).

        #Clear the output file first. Maybe we should update this so it checks
        #for an existing file first, stops and warns the user, and then the user
        #has the option of using something like a "--force" parameter to just
        #overwrite the existing file anyway?
        try:
            os.remove(self.updated_annot_filename)
        except OSError:
            pass

    def update_annotation(self):
        """Main work-horse function that generates the updated annotation.
        """

        #Load indel offsets from the indel file
        indel_offsets = UpdateAnnotationForGenome._get_offsets_from_variant_file(self.genome_indel_filename)

        with open(self.input_annot_filename, 'r') as input_annot_file, \
                open(self.updated_annot_filename, 'w') as updated_annot_file, \
                open(self.log_filename, 'w') as log_file:

            #Copy header from annotation file
            annot_feature = input_annot_file.readline()
            updated_annot_file.write(annot_feature)

            current_chrom = ""

            for annot_feature in input_annot_file:

                annot_feature = annot_feature.rstrip('\n')
                line_data = annot_feature.split('\t')

                if current_chrom != line_data[0]:
                    current_chrom = line_data[0]
                    log_file.write(f"Processing indels and annotated features from chromosome {current_chrom}.\n")

                    if current_chrom in indel_offsets:
                        """
                        Since code below will be performing many lookups and index-
                        based references to the values and keys in current_chrom_variants,
                        it will likely be more efficient to create a list of values
                        and a list of keys from current_chrom_variants once, rather
                        than re-creating them each time the code needs to access a
                        key or value by ordered index.
                        """
                        current_chrom_variant_coords = list(indel_offsets[current_chrom].keys())
                        current_chrom_variant_offsets = list(indel_offsets[current_chrom].values())
                    else:
                        #New chromosome contains no variants
                        log_file.write(f"----No indels from chromosome {current_chrom}.\n")
                        current_chrom_variant_coords = ()
                        current_chrom_variant_coords = ()

                #Current chromosome contains variants
                if current_chrom_variant_coords:

                    tx_start = int(line_data[2])
                    tx_end = int(line_data[3])
                    #exon_count = int(line_data[4])
                    exon_starts = [int(coord) for coord in line_data[5].split(',')]
                    exon_ends = [int(coord) for coord in line_data[6].split(',')]

                    #bisect_right() finds the index at which to insert the given
                    #coordinate in sorted order. Since I'm looking for the
                    #closest coordinate <= the given coordinate, subtract 1 from
                    #the result of bisect_right() to get the correct index.
                    tx_start_offset_index = bisect.bisect_right(current_chrom_variant_coords, tx_start) - 1
                    tx_end_offset_index = bisect.bisect_right(current_chrom_variant_coords, tx_end) - 1

                    #No indels before start of current feature.
                    if tx_start_offset_index == -1:

                        updated_tx_start = tx_start

                        #No indels before end of current feature
                        if tx_end_offset_index == -1:
                            updated_tx_end = tx_end
                            updated_exon_starts = exon_starts
                            updated_exon_ends = exon_ends
                        #First indels occur before end of current feature
                        else:
                            updated_tx_end = tx_end + current_chrom_variant_offsets[tx_end_offset_index]

                            updated_exon_starts = []
                            updated_exon_ends = []
                            for coord in exon_starts:
                                ex_coord_offset_index = bisect.bisect_right(current_chrom_variant_coords, coord) - 1
                                updated_exon_coord = coord
                                if ex_coord_offset_index >= 0:
                                    updated_exon_coord += current_chrom_variant_offsets[ex_coord_offset_index]
                                updated_exon_starts.append(updated_exon_coord)
                            for coord in exon_ends:
                                ex_coord_offset_index = bisect.bisect_right(current_chrom_variant_coords, coord) - 1
                                updated_exon_coord = coord
                                if ex_coord_offset_index >= 0:
                                    updated_exon_coord += current_chrom_variant_offsets[ex_coord_offset_index]
                                updated_exon_ends.append(updated_exon_coord)
                    #No new variants between the start and stop coordinates, so
                    #apply the same offset to all coordinates in the current
                    #feature.
                    elif tx_start_offset_index == tx_end_offset_index:
                        offset = current_chrom_variant_offsets[tx_start_offset_index]
                        updated_tx_start = tx_start + offset
                        updated_tx_end = tx_end + offset
                        updated_exon_starts = [coord+offset for coord in exon_starts]
                        updated_exon_ends = [coord+offset for coord in exon_ends]
                    else:
                        updated_tx_start = tx_start + current_chrom_variant_offsets[tx_start_offset_index]
                        updated_tx_end = tx_end + current_chrom_variant_offsets[tx_end_offset_index]

                        #Update lists of exon starts/ends with correct offsets
                        updated_exon_starts = []
                        updated_exon_ends = []
                        for coord in exon_starts:
                            ex_coord_offset_index = bisect.bisect_right(current_chrom_variant_coords, coord) - 1
                            updated_exon_coord = coord + current_chrom_variant_offsets[ex_coord_offset_index]
                            updated_exon_starts.append(updated_exon_coord)
                        for coord in exon_ends:
                            ex_coord_offset_index = bisect.bisect_right(current_chrom_variant_coords, coord) - 1
                            updated_exon_coord = coord + current_chrom_variant_offsets[ex_coord_offset_index]
                            updated_exon_ends.append(updated_exon_coord)

                    #Format updated annotation data and output
                    updated_annot_file.write(
                        Utils.annot_output_format.format(
                            chrom=line_data[0],
                            strand=line_data[1],
                            txStart=updated_tx_start,
                            txEnd=updated_tx_end,
                            exonCount=line_data[4],
                            exonStarts=','.join([str(x) for x in updated_exon_starts]),
                            exonEnds=','.join([str(x) for x in updated_exon_ends]),
                            transcriptID=line_data[7],
                            geneID=line_data[8],
                            geneSymbol=line_data[9],
                            biotype=line_data[10]
                        )
                    )

                #No variants in the current chromosome, so no need to update
                #feature coordinates.
                else:
                    updated_annot_file.write(f"{annot_feature}\n")

    @staticmethod
    def _get_offsets_from_variant_file(genome_indel_filename):
        """Read indel file, calculate rolling offset at each variant position
        and return results as a dictionary of ordered dictionaries, indexed by
        chromoeomse name.

        This method requires the genome_indel_filename attribute is set and
        contains a valid filename.

        Parameters
        ----------
        genome_indel_filename : string
            Full path to indel file generated by GenomeFilesPreparation.

        Returns
        -------
        OrderedDict nested in defaultdict
            Ordered collection of rolling offsets for each chromosome.
            For outer defaultdict:
                Key = chromosome/contig name from indel file
                Value = OrderedDict (see below)
            For inner OrderedDict:
                Key = chromosomal coordinate of variant position
                Value = rolling offest at variant position
            So variant_offsets["chr1"][12345] stores the rolling offset at
            position 12345 on chromosome 1.

        """

        with open(genome_indel_filename, 'r') as genome_indel_file:

            """
            By using the defaultdict object, I can specify the default value
            used every time I create a new key. This way, I don't need to
            include code to check and instantiate keys with empty OrderedDict
            objects. It's all handed by the defaultdict
            """
            variant_offsets = collections.defaultdict(collections.OrderedDict)
            rolling_offset = 0

            for line in genome_indel_file:
                line_data = line.split('\t')
                indel_chrom, indel_position = line_data[0].split(':')
                indel_position = int(indel_position)
                indel_type = line_data[1]
                indel_offset = int(line_data[2])

                #If variant is a deletion, make offset negative so it subtracts
                #from the rolling offset.
                if indel_type == 'D':
                    indel_offset *= -1

                rolling_offset += indel_offset
                variant_offsets[indel_chrom][indel_position] = rolling_offset

            return variant_offsets

    @staticmethod
    def main():
        """Entry point into script.

        Parses arguments, gathers input and output filenames, and calls scripts
        that perform the actual operation.

        Returns
        -------
        type
            Description of returned object.

        """
        parser = argparse.ArgumentParser(description='Update annotation file with'
                                                     ' coordinates for variant genome')
        parser.add_argument('-g', '--genome_indel_filename')
        parser.add_argument('-i', '--input_annot_filename')
        parser.add_argument('-o', '--updated_annot_filename')
        parser.add_argument('-l', '--log_filename')

        args = parser.parse_args()
        updated_annotation = UpdateAnnotationForGenome(args.genome_indel_filename,
                                                       args.input_annot_filename,
                                                       args.updated_annot_filename,
                                                       args.log_filename)
        start = timer()
        updated_annotation.update_annotation()
        end = timer()
        print(f"Updated annotation: {end - start}")

if __name__ == '__main__':
    sys.exit(UpdateAnnotationForGenome.main())