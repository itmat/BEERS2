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

    def execute(self, sample, reference_genome_file_path, star_file_path, mode="serial"):
        # Right now, options are mode="serial", "parallel" and "lsf"
        # and mode="parallel" is the same as "serial"

        #TODO: module load STAR-v2.6.0c
        assert 1 <= len(sample.input_file_paths) <= 2
        read_files = ' '.join(sample.input_file_paths)

        out_file_prefix = os.path.join(self.data_directory_path, f"sample{sample.sample_id}", "genome_alignment.")

        stdout_log = os.path.join(self.log_directory_path, f"sample{sample.sample_id}", "star_bsub.out")
        stderr_log = os.path.join(self.log_directory_path, f"sample{sample.sample_id}", "star_bsub.err")
        star_index_path = os.path.join(reference_genome_file_path, "genome")

        #TODO: always use zcat?
        star_command = ' '.join([f"{star_file_path}",
                           f"--outFileNamePrefix {out_file_prefix}",
                           f"--genomeDir {star_index_path}",
                           f"--runMode alignReads",
                           f"--runThreadN 4",
                           f"--outSAMtype BAM Unsorted",
                           f"--outFilterMismatchNmax 33",
                           f"--seedSearchStartLmax 33",
                           f"--alignSJoverhangMin 8",
                           f"--readFilesCommand zcat",
                           f"--readFilesIn {read_files}" ])

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
