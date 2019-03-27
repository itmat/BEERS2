import os
import sys
import io
import shutil
import subprocess
import re #Probably won't need this if we switch to an LSF API

import pysam

class GenomeAlignmentStep():

    name = "Genome Alignment Step"

    #Suffix STAR uses for naming all output BAM files.
    star_bamfile_suffix = "Aligned.sortedByCoord.out.bam"

    #Regex for recognizing and extracting lsf job IDs following submission.
    lsf_bsub_output_pattern = re.compile(r'Job <(?P<job_id>\d+?)> is submitted .*')

    def __init__(self, log_directory_path, data_directory_path, parameters = dict()):
        self.log_directory_path = log_directory_path
        self.data_directory_path = data_directory_path
        self.options = parameters

    def validate(self):
        invalid_STAR_parameters = ["--outFileNamePrefix", "--genomeDir", "--runMode", "--outSAMtype", "--readFilesIn"]
        for key, value in self.options.items():
            if not key.startswith("--"):
                print(f"Genome Alignment parameter {key} with value {value} needs to be a STAR option starting with double dashes --", sys.stderr)
                return False
            if key in invalid_STAR_parameters:
                print(f"Genome Alignment parameter {key} with value {value} cannot be used as a STAR option since the value is already determined by the genome alignment script.")
                return False

        return True

    def execute(self, sample, reference_genome_file_path, star_file_path, mode="serial"):
        # Right now, options are mode="serial", "parallel" and "lsf"
        # and mode="parallel" is the same as "serial"

        job_id = ""

        # Check if they are already bam files and so we don't need to re-run the alignment
        if sample.bam_file_path is not '':
            print(f"sample{sample.sample_id} {sample.sample_name} is already aligned.")
            job_id = "ALREADY_ALIGNED"

            # Give the path to the already-existing bam file
            return sample.bam_file_path, job_id

        read_files = ' '.join(sample.fastq_file_paths)

        out_file_prefix = os.path.join(self.data_directory_path, f"sample{sample.sample_id}", "genome_alignment.")

        stdout_log = os.path.join(self.log_directory_path, f"sample{sample.sample_id}", "Genome_Alignment.STAR.bsub.%J.out")
        stderr_log = os.path.join(self.log_directory_path, f"sample{sample.sample_id}", "Genome_Alignment.STAR.bsub.%J.err")
        star_index_path = os.path.join(reference_genome_file_path, "genome")

        # STAR's output bam files are always of this path
        bam_output_file = f"{out_file_prefix}Aligned.sortedByCoord.out.bam"

        options = ' '.join( f"{key} {value}" for key,value in self.options.items() )
        star_command = (f"{star_file_path}"
                          f" --outFileNamePrefix {out_file_prefix}"
                          f" --genomeDir {star_index_path}"
                          f" --runMode alignReads"
                          f" --outSAMtype BAM SortedByCoordinate"
                          f" {options}"
                          f" --readFilesIn {read_files}")

        bsub_command = (f"bsub"
                        f" -n 4"
                        f" -M 40000"
                        f" -R \"span[hosts=1]\""
                        f" -J Genome_Alignment.STAR.sample{sample.sample_id}_{sample.sample_name}"
                        f" -oo {stdout_log}"
                        f" -eo {stderr_log}"
                        f" {star_command}")

        if mode == "serial" or mode == "parallel":
            print(f"Starting STAR on sample {sample.sample_name}.")
            result =  subprocess.run(star_command, shell=True, check=True)
            print(f"Finished running STAR on {sample.sample_name}.")
        elif mode == "lsf":
            print(f"Submitting STAR command to {mode} for sample {sample.sample_name}.")
            result = subprocess.run(bsub_command, shell=True, check=True, stdout=subprocess.PIPE, encoding="ascii")

            #Extract job ID from LSF stdout
            job_id = GenomeAlignmentStep.lsf_bsub_output_pattern.match(result.stdout).group("job_id")
            print(f"\t{result.stdout.rstrip()}")

            print(f"Finished submitting STAR command to {mode} for sample {sample.sample_name}.")

        # Return the path of the output so that later steps can use it
        # Or move it to a standard location?
        return bam_output_file, job_id

    #TODO: We might want to remove the bam_file argument from this method and
    #      have it generate the bai filename using the same rules the BAM file
    #      was named in the code above.
    def index(self, sample, bam_file, mode = "serial"):
        ''' Build index of a given bam file '''
        # indexing is necessary for the variant finding step

        job_id = ""
        #If there are existing BAM files, BAI files will be made in the same
        #directory. The BAI files will subsequently be copied to the output path.
        out_file_path = os.path.join(self.data_directory_path, f"sample{sample.sample_id}", 'genome_alignment.Aligned.sortedByCoord.out.bam.bai')

        if mode == "serial" or mode == "parallel":
            print(f"Creating BAM index for sample {sample.sample_name}")
            pysam.index(bam_file)
            if not os.path.isfile(out_file_path):
                shutil.copyfile(f"{bam_file}.bai", out_file_path)
            print(f"Finished creating BAM index for sample {sample.sample_name}")
        elif mode == "lsf":
            print(f"Submitting BAM index command to {mode} for sample {sample.sample_name}.")
            #TODO: This code currently uses whichever samtools binary is in the
            #      system path and will fail if there isn't one present. To avoid
            #      requiring a separate samtools install, it might be better to
            #      write a small wrapper script around the pysam.index() method,
            #      so we can submit that script directly to bsub. This way we
            #      are effectively doing the same thing as in serial/parallel mode,
            #      but through an lsf submission.

            stdout_log = os.path.join(self.log_directory_path, f"sample{sample.sample_id}", "Index_Genome_BAM.bsub.%J.out")
            stderr_log = os.path.join(self.log_directory_path, f"sample{sample.sample_id}", "Index_Genome_BAM.bsub.%J.err")

            # TODO: no need to copy this file at all, except that we are currently using it's {out_file_path}
            #       to determine if the job was successful. In the future, we should check it anywhere and not need a copy
            if f"{bam_file}.bai" != out_file_path:
                # Need to copy the output file to the expected location
                # for the case the reads were already aligned (a bam file given in the config)
                copy_command = f"cp -p {bam_file}.bai {out_file_path}"
            else:
                copy_command = '' # no need to copy, it's already the right file
            bsub_command = (f"bsub"
                            f" -J Index_Genome_BAM.sample{sample.sample_id}_{sample.sample_name}"
                            f" -oo {stdout_log}"
                            f" -eo {stderr_log}"
                            f" 'python -m beers.expression.run_pysam_index {bam_file}; {copy_command}'")
            result = subprocess.run(bsub_command, shell=True, check=True, stdout = subprocess.PIPE, encoding="ascii")

            #Extract job ID from LSF stdout
            job_id = GenomeAlignmentStep.lsf_bsub_output_pattern.match(result.stdout).group("job_id")
            print(f"\t{result.stdout.rstrip()}")

            print(f"Finished submitting BAM index command to {mode} for sample {sample.sample_name}.")
        else:
            raise NotImplementedError("dispatcher mode must be serial, parallel, or lsf")

        return job_id

    @staticmethod
    def is_output_valid(job_attributes):
        """
        Check if output of GenomeAlignment for a specific job/execution is
        correctly formed and valid, given a job's data directory, log directory,
        and sample id.

        Parameters
        ----------
        job_attributes : dict
            A job's data_directory, log_directory, and sample_id.

        Returns
        -------
        boolean
            True  - GenomeAlignment output files were created and are well formed.
            False - GenomeAlignment output files do not exist or are missing data.
        """
        data_directory = job_attributes['data_directory']
        log_directory = job_attributes['log_directory']
        sample_id = job_attributes['sample_id']

        valid_output = False

        alignment_outfile_path = os.path.join(data_directory, f"sample{sample_id}", "genome_alignment.Aligned.sortedByCoord.out.bam")
        alignment_logfile_path = os.path.join(log_directory, f"sample{sample_id}", "genome_alignment.Log.progress.out")
        if os.path.isfile(alignment_outfile_path) and \
           os.path.isfile(alignment_logfile_path):
            #Read last line in variants_finder log file
            line = ""
            with open(alignment_logfile_path, "r") as alignment_log_file:
                for line in alignment_log_file:
                    line = line.rstrip()
            if line == "ALL DONE!":
                valid_output = True

        return valid_output
