"""
Module containing a rough draft of transcriptome preparation including:
    1) creation of transcriptome fasta files from genome fasta files and annotations
    2) creation of STAR indexes for transcriptomes
    3) alignment of fastq files to STAR indexes
    4) Quantification at transcript and allele levels
    5) Expression of molecules based on quantifications

Runs these tasks in parallel through bsub
"""

import subprocess
import os
import sys
import argparse
import time
import pickle

import numpy

from beers_utils.constants import MAX_SEED
from camparee.gene_files_preparation import GeneFilesPreparation
from camparee.transcript_gene_quant import TranscriptGeneQuantificationStep
from camparee.allelic_imbalance_quant import AllelicImbalanceQuantificationStep
from camparee.molecule_maker import MoleculeMaker
from beers_utils.sample import Sample
from camparee.abstract_camparee_step import AbstractCampareeStep

# TODO: Split this class up into each of the separate steps. This way those that
#       are independent from each other can run in parallel.
class TranscriptQuantificatAndMoleculeGenerationStep():

    MOLECULES_PER_PACKET = 10_000

    TRANSCRIPTOMES_LOG_FILENAME = "TranscriptQuantificatAndMoleculeGenerationStep.log"

    def __init__(self, log_directory_path, data_directory_path, parameters=dict()):
        self.data_directory_path = data_directory_path
        self.log_directory_path = log_directory_path
        self.num_bowtie_threads = parameters.get('num_bowtie_threads')

        if not self.num_bowtie_threads:
            self.num_bowtie_threads = 1

    def validate(self):
        return True

    def done_file_name(self, sample_id):
        return os.path.join(self.data_directory_path, f"sample{sample_id}", "transcriptome_prep_done.txt")


    def execute(self, sample, kallisto_file_path, bowtie2_dir_path, output_type, output_molecule_count, seed=None):
        if seed is not None:
            numpy.random.seed(seed)
        sample_dir = os.path.join(self.data_directory_path, f"sample{sample.sample_id}")
        transcriptome_logfile_path = os.path.join(self.log_directory_path, f"sample{sample.sample_id}",
                                                  TranscriptQuantificatAndMoleculeGenerationStep.TRANSCRIPTOMES_LOG_FILENAME)
        fastq_file_1, fastq_file_2 = sample.fastq_file_paths

        # Prepare the transcriptome fasta files
        for i in [1,2]:
            print(f"Preparing fasta files for transcriptome {i} of sample{sample.sample_id}")
            genome = os.path.join(sample_dir, f"custom_genome_{i}.fa")
            geneinfo = os.path.join(sample_dir, f"updated_annotation_{i}.txt")
            transcriptome = os.path.join(sample_dir, f"transcriptome_{i}.fa")
            file_prep = GeneFilesPreparation(genome, geneinfo, transcriptome)
            file_prep.prepare_gene_files()

        # Build Kallisto indexes
        for i in [1, 2]:
            print(f"Building Kallisto indexes for transcriptome {i} of sample{sample.sample_id}")
            transcriptome = os.path.join(sample_dir, f"transcriptome_{i}.fa")
            transcriptome_index_dir = os.path.join(sample_dir, f"transcriptome_{i}_index")
            os.mkdir(transcriptome_index_dir)
            transcriptome_index = os.path.join(transcriptome_index_dir, f"transcriptome_{i}.kallisto.index")

            command = f"{kallisto_file_path} index -i {transcriptome_index} {transcriptome}"
            subprocess.run(command, shell=True, check=True)


        # Kallisto Quantification Step
        for i in [1, 2]:
            print(f"Running Kallisto quantification for transcriptome {i} of sample{sample.sample_id}")
            transcriptome = os.path.join(sample_dir, f"transcriptome_{i}.fa")
            transcriptome_index_dir = os.path.join(sample_dir, f"transcriptome_{i}_index")
            transcriptome_index = os.path.join(transcriptome_index_dir, f"transcriptome_{i}.kallisto.index")

            kallisto_output_dir = os.path.join(sample_dir, f"transcriptome_{i}_kallisto_quant")
            os.mkdir(kallisto_output_dir)
            command = f"{kallisto_file_path} quant -i {transcriptome_index} -o {kallisto_output_dir} {fastq_file_1} {fastq_file_2}"
            subprocess.run(command, shell=True, check=True)

        # Run transcript/gene/allelic_imbalance quantification scripts
        geneinfo_filename = os.path.join(sample_dir, "updated_annotation_1.txt")

        print(f"Quantifying transcriptome reads for sample{sample.sample_id}")
        transcript_gene_quant = TranscriptGeneQuantificationStep(geneinfo_filename, sample_dir)
        transcript_gene_quant.make_transcript_gene_dist_file()

        # Build Bowtie2 indexes
        for i in [1,2]:
            print(f"Building Bowtie2 indexes for transcriptome {i} of sample{sample.sample_id}")
            transcriptome = os.path.join(sample_dir, f"transcriptome_{i}.fa")
            bowtie2_build_file_path = os.path.join(bowtie2_dir_path, f"bowtie2-build")
            transcriptome_index_dir = os.path.join(sample_dir, f"transcriptome_{i}_index")
            transcriptome_bowtie2_index = os.path.join(transcriptome_index_dir, f"transcriptome_{i}")

            command = f"{bowtie2_build_file_path} --threads {self.num_bowtie_threads} {transcriptome} {transcriptome_bowtie2_index}"
            subprocess.run(command, shell=True, check=True)

        # Align fastq files to Bowtie2 indexes
        for i in [1, 2]:
            print(f"Aligning by Bowtie2 indexes for transcriptome {i} of sample{sample.sample_id}")
            transcriptome = os.path.join(sample_dir, f"transcriptome_{i}.fa")
            bowtie2_file_path = os.path.join(bowtie2_dir_path, f"bowtie2")
            transcriptome_index_dir = os.path.join(sample_dir, f"transcriptome_{i}_index")
            transcriptome_bowtie2_index = os.path.join(transcriptome_index_dir, f"transcriptome_{i}")

            aligned_file_path = os.path.join(sample_dir, f"{i}_Aligned.out.sam")
            command = f"{bowtie2_file_path} --threads {self.num_bowtie_threads} --very-sensitive -x {transcriptome_bowtie2_index} -1 {fastq_file_1} -2 {fastq_file_2} -S {aligned_file_path}"
            subprocess.run(command, shell=True, check=True)

        print(f"Quantifying allelic imabalance for sample{sample.sample_id}")
        allelic_imbalance_quant = AllelicImbalanceQuantificationStep(sample_dir)
        allelic_imbalance_quant.quantify_allelic_imbalance()
        allelic_imbalance_quant.make_allele_imbalance_dist_file()

        # Run Molecule maker to generate a packet
        # TODO: what sample object can we use?
        print(f"Generating molecules for sample{sample.sample_id}")
        molecule_maker = MoleculeMaker(sample, sample_dir)
        if output_type == "packet":
            # TODO: potentially rounds down the number of molecules to make
            # TODO: MOLECULES_PER_PACKET maybe should not be a constant
            num_packets = output_molecule_count // TranscriptQuantificatAndMoleculeGenerationStep.MOLECULES_PER_PACKET
            for i in range(1,num_packets+1):
                print(f"Generating packet {i} of {num_packets} for sample{sample.sample_id}")
                packet = molecule_maker.make_packet(id=f"sample{sample.sample_id}.{i}")

                with open(os.path.join(sample_dir, f"molecule_packet{i}.pickle"), "wb") as out_file:
                    pickle.dump(packet, out_file)
        elif output_type == "molecule_file":
            molecule_maker.make_molecule_file(filepath=os.path.join(sample_dir, f"molecule_file"),
                                              N = output_molecule_count)

        print(f"Done with transcriptome.py for sample{sample.sample_id}")

        # Output to DONE file
        with open(transcriptome_logfile_path, "w") as done_file:
            done_file.write("\nALL DONE!\n")

    def get_commandline_call(self, sample, kallisto_file_path, bowtie2_dir_path, output_type, output_molecule_count, seed=None):
        #Retrieve path to the transcriptomes.py script.
        transcriptomes_path = os.path.realpath(__file__)
        #If the above command returns a string with a "pyc" extension, instead
        #of "py", strip off "c" so it points to this script.
        transcriptomes_path = transcriptomes_path.rstrip('c')

        command = (f" python {transcriptomes_path}"
                   f" --log_directory_path {self.log_directory_path}"
                   f" --data_directory_path {self.data_directory_path}"
                   f" --sample '{repr(sample)}'"
                   f" --kallisto_file_path {kallisto_file_path}"
                   f" --bowtie2_dir_path {bowtie2_dir_path}"
                   f" --output_type {output_type}"
                   f" --output_molecule_count {output_molecule_count}"
                   f" --num_bowtie_threads {self.num_bowtie_threads}")

        if seed is not None:
            command += f" --seed {seed}"

        return command

    def get_validation_attributes(self, sample, kallisto_file_path, bowtie2_dir_path, output_type, output_molecule_count, seed=None):
        validation_attributes = {}
        validation_attributes['data_directory'] = self.data_directory_path
        validation_attributes['log_directory'] = self.log_directory_path
        validation_attributes['sample_id'] = sample.sample_id
        validation_attributes['output_molecule_count'] = output_molecule_count
        return validation_attributes

    @staticmethod
    def main():
        parser = argparse.ArgumentParser()
        parser.add_argument("--log_directory_path", help="Directory for log files")
        parser.add_argument("--data_directory_path", help="Directory for the data (input + output)")
        parser.add_argument("--sample", help="Serialized Sample object containing FASTQ file info.")
        parser.add_argument("--kallisto_file_path", help="path to Kallisto executable")
        parser.add_argument("--bowtie2_dir_path", help="path to Bowtie2 directory")
        parser.add_argument("--output_type", help="type of output file to write", choices=["packet", "molecule_file"], default="molecule_file")
        parser.add_argument("--output_molecule_count", help="number of molecules to output", type=int)
        parser.add_argument("--num_bowtie_threads", type=int, default=None,
                            help="Number of threads to use for Bowtie alignment")
        parser.add_argument("--seed", help="seed for RNG", default=None, type=int)
        args = parser.parse_args()

        config_parameters = {}
        config_parameters['num_bowtie_threads'] = args.num_bowtie_threads
        transcriptomes = TranscriptQuantificatAndMoleculeGenerationStep(args.log_directory_path,
                                                                        args.data_directory_path,
                                                                        config_parameters)
        sample = eval(args.sample)
        transcriptomes.execute(sample,
                               args.kallisto_file_path,
                               args.bowtie2_dir_path,
                               args.output_type,
                               args.output_molecule_count,
                               args.seed)

    @staticmethod
    def is_output_valid(validation_attributes):
        data_directory = validation_attributes['data_directory']
        log_directory = validation_attributes['log_directory']
        sample_id = validation_attributes['sample_id']
        output_molecule_count = validation_attributes['output_molecule_count']

        valid_output = False

        transcriptome_logfile_path = os.path.join(log_directory, f"sample{sample_id}",
                                                  TranscriptQuantificatAndMoleculeGenerationStep.TRANSCRIPTOMES_LOG_FILENAME)
        molecule_file = os.path.join(data_directory, f"sample{sample_id}", "molecule_file")
        if os.path.isfile(transcriptome_logfile_path) and \
           os.path.isfile(molecule_file):
            #Read last line in transcriptomes log file
            line = ""
            with open(transcriptome_logfile_path, "r") as transcriptome_log_file:
                for line in transcriptome_log_file:
                    line = line.rstrip()
            if line == "ALL DONE!":
                valid_output = True

        return valid_output

if __name__ == '__main__':
    sys.exit(TranscriptQuantificatAndMoleculeGenerationStep.main())
