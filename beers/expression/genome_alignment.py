import os
import sys
import io
import shutil
import subprocess

import pysam

class GenomeAlignmentStep():

    name = "Genome Alignment Step"

    def __init__(self, log_directory_path, data_directory_path, parameters):
        self.log_directory_path = log_directory_path
        self.data_directory_path = data_directory_path

    def validate(self):
        # TODO: validate and use parameters for STAR
        return True

    def execute(self, sample, reference_genome_file_path, star_file_path, mode="serial"):
        # Right now, options are mode="serial", "parallel" and "lsf"
        # and mode="parallel" is the same as "serial"

        #TODO: module load STAR-v2.6.0c
        assert 1 <= len(sample.input_file_paths) <= 2

        # Check if they are already bam files and so we don't need to re-run the alignment
        if sample.input_file_paths[0].endswith("bam"):
            assert len(sample.input_file_paths) == 1
            print(f"sample{sample.sample_id} {sample.sample_name} is already aligned.")

            # Give the path to the already-existing bam file
            return sample.input_file_paths[0]

        read_files = ' '.join(sample.input_file_paths)

        out_file_prefix = os.path.join(self.data_directory_path, f"sample{sample.sample_id}", "genome_alignment.")

        stdout_log = os.path.join(self.log_directory_path, f"sample{sample.sample_id}", "star_bsub.out")
        stderr_log = os.path.join(self.log_directory_path, f"sample{sample.sample_id}", "star_bsub.err")
        star_index_path = os.path.join(reference_genome_file_path, "genome")

        # STAR's output bam files are always of this path
        bam_output_file = f"{out_file_prefix}Aligned.sortedByCoord.out.bam"

        #TODO: always use zcat?
        star_command = (f"{star_file_path}"
                          f" --outFileNamePrefix {out_file_prefix}"
                          f" --genomeDir {star_index_path}"
                          f" --runMode alignReads"
                          f" --runThreadN 4"
                          f" --outSAMtype BAM SortedByCoordinate"
                          f" --outFilterMismatchNmax 33"
                          f" --seedSearchStartLmax 33"
                          f" --alignSJoverhangMin 8"
                          f" --readFilesCommand zcat"
                          f" --readFilesIn {read_files}")

        bsub_command = (f"bsub"
                         f" -n 4"
                         f" -M 40000"
                         f" -R \"span[hosts=1]\""
                         f" -J STAR_{sample.sample_id}_{sample.sample_name}"
                         f" -oo {stdout_log}"
                         f" -eo {stderr_log}"
                         f" {star_command}")

        if mode == "serial" or mode == "parallel":
            print(f"Starting STAR on sample {sample.sample_name}.")
            result =  subprocess.run(star_command, shell=True, check=True)
            print(f"Finished running STAR on sample{sample.sample_id} {sample.sample_name}.")
        elif mode == "lsf":
            #TODO: implement lsf support - or move this functionality to other files
            raise NotImplementedError()

            # bsub it and capture stdout to read bsub's results (job id)
            result =  subprocess.run(bsub_command, shell=True, check=True, stdout = subprocess.PIPE, encoding="ascii")
            print(result.stdout)

        # Return the path of the output so that later steps can use it
        # Or move it to a standard location?
        return bam_output_file

    def index(self, sample, bam_file, mode = "serial"):
        if mode == "serial" or mode == "parallel":
            print(f"Creating index for sample{sample.sample_id} {sample.sample_name}")
            pysam.index(bam_file)
            print(f"Done creating index for sample{sample.sample_id} {sample.sample_name}")
        elif mode == "lsf":
            raise NotImplementedError()

            bsub_command = (f"bsub samtools index {bam_file}")
            result = subprocess.run(bsub_command, shell=True, check=True, stdout = subprocess.PIPE, encoding="ascii")
            print(result.stdout)
