import contextlib
import os
import sys
import argparse
import json

import numpy

from camparee.camparee_utils import CampareeUtils
from camparee.abstract_camparee_step import AbstractCampareeStep

class VariantsCompilationStep(AbstractCampareeStep):

    #Name of the file where script output is stored.
    VARIANTS_OUTPUT_FILENAME = "all_variants.vcf"
    #Name of the file where script logging is stored.
    VARIANTS_LOG_FILENAME = "VariantsCompilationStep.log"

    def __init__(self, log_directory_path, data_directory_path, parameters=None):
        self.data_directory_path = data_directory_path
        self.log_directory_path = log_directory_path
        self.chr_ploidy_data = None

    def validate(self):
        return True

    def execute(self, sample_id_list, chr_ploidy_data, reference_genome, seed=None):
        """
        Entry point into variants_compilation.

        Parameters
        ----------
        sample_id_list : list
            List of sample IDs
        chr_ploidy_data : dict
            Dictionary of chromosomes as keys and a dictionary of male/female
            ploidy as values.
        reference_genome : dict
            Dictionary representation of the reference genome
        seed : int
            Seed for random number generator. Used so repeated runs will produce
            the same results.

        """
        self.chr_ploidy_data = chr_ploidy_data
        contig_order = list(reference_genome.keys())
        log_file_path = os.path.join(self.log_directory_path, VariantsCompilationStep.VARIANTS_LOG_FILENAME)

        # The common_variant() method defined below uses a random number when
        # choosing which of two equally prevalent variants to keep.
        if seed is not None:
            numpy.random.seed(seed)

        with open(log_file_path, "w") as log_file:
            print("Converting variants into vcf file")
            log_file.write("Converting variants into vcf file\n")
            contigs_so_far = []
            last_chromosome = None
            # Open then process all the files
            all_variants_file_path = os.path.join(self.data_directory_path, self.VARIANTS_OUTPUT_FILENAME)
            with contextlib.ExitStack() as stack, open(all_variants_file_path, "w") as out:
                # Output the header
                out.write("##fileformat=VCFv4.0\n")
                sample_ids = '\t'.join(['sample' + str(sample_id) for sample_id in sample_id_list])
                out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_ids}\n")

                # Open variant files
                variant_file_paths = [os.path.join(self.data_directory_path,
                                                   'sample' + str(sample_id), "variants.txt")
                                      for sample_id in sample_id_list]
                variant_files = [stack.enter_context(open(file_path)) for file_path in variant_file_paths]

                # read the first line of each variant file
                next_lines = [vf.readline() for vf in variant_files]

                # proceed until we have exhausted (and removed from the list) all the variant files
                i = 0
                while any(line != '' for line in next_lines):
                    i += 1
                    parsed_lines = [CampareeUtils.parse_variant_line(line) for line in next_lines]
                    chromosome = min(chr for chr, pos, vars in parsed_lines if chr != "DONE")
                    position = min(pos for chr, pos, vars in parsed_lines if chr == chromosome)
                    if chromosome != last_chromosome:
                        if last_chromosome is not None:
                            print(f"Wrote chromosome {last_chromosome} to vcf file")
                            log_file.write(f"Wrote chromosome {last_chromosome} to vcf file\n")
                            contigs_so_far.append(last_chromosome)
                            if contig_order.index(last_chromosome) > contig_order.index(chromosome):
                                print(f"WARNING: chromosomes {last_chromosome}, {chromosome} are not in the same order as in the reference genome")
                                log_file.write(f"WARNING: chromosomes {last_chromosome}, {chromosome} are not in the same order as in the reference genome\n")
                            if sorted(contigs_so_far) != contigs_so_far:
                                print("ERROR: Need chromosomes to be in alphabetical order (i.e. 1, 10, 11, 12, .., 19, 2, 20, 21, 22, 3, .., MT, X,Y")
                                log_file.write("ERROR: Need chromosomes to be in alphabetical order (i.e. 1, 10, 11, 12, .., 19, 2, 20, 21, 22, 3, .., MT, X,Y\n")
                                # TODO: probably will need to use an alternative means of sorting them since alphabetical is terrible but what we currently use in the test files
                        last_chromosome = chromosome

                    ref_base = reference_genome[chromosome][position - 1]

                    # List of variant dictionaries (variant -> count) with empty dicts for files that have no entry for this loc
                    variants = [vars if (chr == chromosome and pos == position) else {}
                                for chr, pos, vars in parsed_lines]

                    # which samples we actually used a variant of at this location
                    used_samples = [True if (chr == chromosome and pos == position) else False
                                    for chr, pos, vars in parsed_lines]

                    # How many variants we had in each sample
                    num_variants = [len(vars) for vars in variants]

                    # Filter variants by those that are NOT just the reference
                    # and do not include "N"s
                    variants = [{var: count for var, count in vars.items() if (var != ref_base and "N" not in var)}
                                for vars in variants]

                    #TODO: note that current versions of variant_finder.py do largely the same job as this already
                    # So one could remove most of these checks at some point, but no harm is done in running these
                    # (except for slowing down the process slightly)

                    # Find each sample's most likely variants that aren't just the reference
                    def common_variant(variants):
                        if variants:
                            maximum_count = max(count for var, count in variants.items())
                            common_variants = [var for var, count in variants.items() if count == maximum_count]
                            # Take random choice of the tied maximally seen variants
                            return numpy.random.choice(common_variants)
                        else:
                            return None

                    most_common_variants = [common_variant(vars) for vars in variants]

                    # Output nothing if we've now discarded all the variants
                    if all(var == None for var in most_common_variants):
                        next_lines = [line if not used else sample_file.readline()
                                      for line, used, sample_file in zip(next_lines, used_samples, variant_files)]
                        continue

                    # Find the reference length of the longest variant (i.e. length of the longest indel, or 1 if only SNPs)
                    def variant_length(variant):
                        if variant.startswith("D"):
                            return int(variant[1:]) + 1  # +1 includes the reference base before our deletion
                        else:
                            return 1

                    length = max(variant_length(variant) for variant in most_common_variants if variant is not None)

                    # Do we have any deletions to handle?
                    # any_dels = any(var.startswith("D") for var in most_common_variants)
                    any_dels = (length > 1)

                    # accommodate the base before the deletion
                    if any_dels:
                        position -= 1

                    # Give reference of the appropriate size to accommodate the longest variant here
                    ref = reference_genome[chromosome][position - 1:position - 1 + length]

                    sample_descriptions = []
                    alts = []
                    for variant, num_vars in zip(most_common_variants, num_variants):
                        if variant is None:
                            sample_descriptions.append("0/0")  # ref-ref if no variants for this sample
                            continue

                        # With deletions, include one extra base
                        if any_dels:
                            if variant in "ACGT":  # SNP
                                alt = ref[0] + variant + ref[2:]
                            elif variant.startswith("I"):
                                alt = ref[0:2] + variant[1:] + ref[2:]
                            elif variant.startswith("D"):
                                l = int(variant[1:])
                                alt = ref[0] + ref[1 + l:]
                            else:
                                raise ValueError(f"unknown variant {variant}")
                        else:
                            if variant in "ACGT":  # SNP
                                alt = variant
                            elif variant.startswith("I"):
                                alt = ref[0] + variant[1:]
                            elif variant.startswith("D"):
                                raise ValueError(f"shouldn't have a deletion here in variant {variant}")
                            else:
                                raise ValueError(f"unknown variant {variant}")

                        if alt not in alts:
                            alts.append(alt)
                        idx = alts.index(alt) + 1  # 0 is the ref, 1 is the first alt, etc.

                        # We assume that a lone variant is going to be alt-alt
                        # but if there are multiple variants, we always use ref as one
                        # and the most common alt as the other
                        # TODO: should we ever use alt-alt with two different alts? we don't currently
                        if num_vars == 1:
                            sample_descriptions.append(f"{idx}/{idx}")  # alt-alt
                        else:
                            sample_descriptions.append(f"0/{idx}")  # ref-alt

                    line = "\t".join([chromosome,
                                      str(position),
                                      ".",
                                      ref,
                                      ",".join(alts),
                                      ".",
                                      ".",
                                      ".",
                                      "GT",
                                      "\t".join(sample_descriptions)])
                    out.write(line + "\n")

                    next_lines = [line if not used else sample_file.readline()
                                  for line, used, sample_file in zip(next_lines, used_samples, variant_files)]

            print(f"Wrote chromosome {last_chromosome} to vcf file")
            print(f"Finished creating VCF file for Beagle with {i} variant entries")
            log_file.write(f"Wrote chromosome {last_chromosome} to vcf file\n")
            log_file.write(f"Finished creating VCF file for Beagle with {i} variant entries\n")
            log_file.write("ALL DONE!\n")

    def get_commandline_call(self, samples, chr_ploidy_file_path, reference_genome_file_path, seed=None):
        """
        Prepare command to execute the VariantsCompilationStep from the command
        line, given all of the arugments used to run the execute() function.

        Parameters
        ----------
        samples : list
            List of Sample() objects for which variants have been called and need
            to be merged.
        chr_ploidy_file_path : string
            File that maps chromosome names to their male/female ploidy.
        reference_genome_file_path : string
            File that maps chromosome names in reference to nucleotide sequence.
        seed : integer
            Seed for random number generator. Used so repeated runs will produce
            the same results.

        Returns
        -------
        string
            Command to execute on the command line. It will perform the same
            operations as a call to execute() with the same parameters.
        """

        #Retrieve path to the variants_compilation.py script.
        variant_compilation_path = os.path.realpath(__file__)
        #If the above command returns a string with a "pyc" extension, instead
        #of "py", strip off "c" so it points to this script.
        variant_compilation_path = variant_compilation_path.rstrip('c')

        command = (f" python {variant_compilation_path}"
                   f" --log_directory_path {self.log_directory_path}"
                   f" --data_directory_path {self.data_directory_path}"
                   f" --sample_ids '{json.dumps(samples)}'"
                   f" --chr_ploidy_file_path {chr_ploidy_file_path}"
                   f" --reference_genome_file_path {reference_genome_file_path}")

        if seed is not None:
            command += f" --seed {seed}"

        return command

    def get_validation_attributes(self, samples, chr_ploidy_file_path, reference_genome_file_path, seed=None):
        """
        Prepare attributes required by is_output_valid() function to validate
        output generated the VariantsCompilationStep job.

        Parameters
        ----------
        samples : list
            List of Sample() objects for which variants have been called and need
            to be merged. [Note: this parameter is captured just so
            get_validation_attributes() accepts the same arguments as
            get_commandline_call(). It is not used here.]
        chr_ploidy_file_path : string
            File that maps chromosome names to their male/female ploidy. [Note:
            this parameter is captured just so get_validation_attributes() accepts
            the same arguments as get_commandline_call(). It is not used here.]
        reference_genome_file_path : string
            File that maps chromosome names in reference to nucleotide sequence.
            [Note: this parameter is captured just so get_validation_attributes()
            accepts the same arguments as get_commandline_call(). It is not used
            here.]
        seed : integer
            Seed for random number generator. Used so repeated runs will produce
            the same results. [Note: this parameter is captured just so
            get_validation_attributes() accepts the same arguments as
            get_commandline_call(). It is not used here.]

        Returns
        -------
        dict
            A VariantsCompilationStep run's data_directory and log_directory.
        """
        validation_attributes = {}
        validation_attributes['data_directory'] = self.data_directory_path
        validation_attributes['log_directory'] = self.log_directory_path
        return validation_attributes


    @staticmethod
    def main():
        """
        Entry point into script. Allows script to be executed/submitted via the
        command line.
        """

        parser = argparse.ArgumentParser(description='Command line wrapper around'
                                                     ' the variant compilation step')
        parser.add_argument('--log_directory_path')
        parser.add_argument('--data_directory_path')
        parser.add_argument('--sample_ids')
        parser.add_argument('--chr_ploidy_file_path')
        parser.add_argument('--reference_genome_file_path')
        parser.add_argument('--seed', type=int, default=None)
        args = parser.parse_args()

        variants_compiler = VariantsCompilationStep(args.log_directory_path,
                                                    args.data_directory_path)
        sample_id_list = json.loads(args.sample_ids)
        reference_genome = CampareeUtils.create_genome(args.reference_genome_file_path)
        chr_ploidy_data = CampareeUtils.create_chr_ploidy_data(args.chr_ploidy_file_path)
        variants_compiler.execute(sample_id_list,
                                  chr_ploidy_data,
                                  reference_genome,
                                  args.seed)

    @staticmethod
    def is_output_valid(validation_attributes):
        """
        Check if output of VariantsCompilationStep for a specific job/execution
        is correctly formed and valid, given the run's data and log directories.
        Prepare these attributes using the get_validation_attributes() method.

        Parameters
        ----------
        validation_attributes : dict
            A CAMPAREE run's data_directory and log_directory.

        Returns
        -------
        boolean
            True  - VariantsCompilationStep output files were created and are
                    well formed.
            False - VariantsCompilationStep output files do not exist or are
                    missing data.
        """
        data_directory = validation_attributes['data_directory']
        log_directory = validation_attributes['log_directory']

        valid_output = False

        output_file_path = os.path.join(data_directory, VariantsCompilationStep.VARIANTS_OUTPUT_FILENAME)
        log_file_path = os.path.join(log_directory, VariantsCompilationStep.VARIANTS_LOG_FILENAME)
        if os.path.isfile(output_file_path) and os.path.isfile(log_file_path):
            #Read last line in variants_finder log file
            line = ""
            with open(log_file_path, "r") as log_file:
                for line in log_file:
                    line = line.rstrip()
            if line == "ALL DONE!":
                valid_output = True

        return valid_output

if __name__ == "__main__":
    sys.exit(VariantsCompilationStep.main())
