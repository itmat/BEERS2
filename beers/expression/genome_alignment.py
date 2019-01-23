import shutil
import os
import subprocess

class GenomeAlignmentStep():

    name = "Genome Alignment Step"

    def __init__(self, log_directory_path, data_directory_path, parameters):
        self.log_directory_path = log_directory_path
        self.data_directory_path = data_directory_path

    def validate(self):
        # TODO: validate and use parameters for STAR
        return True

    def execute(self, sample, reference_genome, mode="serial"):
        # Right now, options are mode="serial", "parallel" and "lsf"
        # and mode="parallel" is the same as "serial"

        #TODO: module load STAR-v2.6.0c
        assert 1 <= len(self.sample.input_file_paths) <= 2
        read_files = ' '.join(self.sample.input_file_paths)

        out_file_prefix = os.path.join(self.data_directory_path, sample.sample_id, "genome_alignment")

        stdout_log = os.path.join(self.log_directory_path, sample.sample_id, "star_bsub.out")
        stderr_log = os.path.join(self.log_directory_path, sample.sample_id, "star_bsub.err")

        #TODO: always use zcat?
        star_command = "STAR \
                            --outFileNamePrefix {out_file_prefix} \
                            --genomeDir {reference_genome} \
                            --runMode alignReads \
                            --runThreadN 4 \
                            --outSAMtype BAM Unsorted \
                            --outFilterMismatchNmax 33 \
                            --seedSearchStartLmax 33 \
                            --alignSJoverhangMin 8 \
                            --readFilesCommand zcat \
                            --readFilesIn {read_files}"

        bsub_command = f"bsub \
                            -n 4 \
                            -M 40000 \
                            -R \"span[hosts=1]\" \
                            -J STAR_{sample.sample_id}_{sample.sample_name} \
                            -oo {stdout_log} \
                            -eo {stderr_log} \
                            {star_command}"

        if mode == "serial" or mode == "parallel":
            print(f"Starting STAR on sample {sample.sample_name}.")
            result =  subprocess.run(star_command, shell=True, check=True)
            print(f"Finished running STAR on sample{sample.sample_id} {sample.sample_name}.")
        elif mode == "lsf":
            #TODO: implement lsf support - or move this functionality to other files
            raise NotImplementedError()
