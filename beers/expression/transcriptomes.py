"""
Module containing a rough draft of transcriptome preparation including:
    1) creation of transcriptome fasta files from genome fasta files and annotations
    2) creation of STAR indexes for transcriptomes
    3) alignment of fastq files to STAR indexes

Runs these tasks in parallel through bsub
"""

import subprocess
import os
import argparse
import time

from beers.personal.gene_files_preparation import GeneFilesPreparation

def done_file_name(data_directory, sample_id):
    return os.path.join(data_directory, f"sample{sample_id}", "transcriptome_prep_done.txt")

def prep_transcriptomes(samples, data_directory,  log_directory, star_file_path, dispatcher_mode="serial"):
    for sample in samples:
        command = f"python -m beers.expression.transcriptomes {sample.sample_id} {data_directory} {log_directory} {star_file_path} --fastq_files {' '.join(sample.input_file_paths)}"

        if dispatcher_mode == "serial":
            subprocess.run(command, shell=True)
        elif dispatcher_mode == "lsf":
            stdout_log = os.path.join(log_directory, f"sample{sample.sample_id}", "Transcriptome.bsub.%J.out")
            stderr_log = os.path.join(log_directory, f"sample{sample.sample_id}", "Transcriptome.bsub.%J.err")
            bsub_command = (f"bsub -M 40000"
                            f" -J Variant_Finder.sample{sample.sample_id}_{sample.sample_name}"
                            f" -oo {stdout_log}"
                            f" -eo {stderr_log}"
                            f" {command}")
            result = subprocess.run(bsub_command, shell=True, check=True, stdout = subprocess.PIPE, encoding="ascii")
            print(result)
        else:
            raise NotImplementedError("only lsf/serial dispatching")

    sample_done = {sample.sample_id: False for sample in samples}
    while True:
        print("Checking on Transcriptome preparation status")
        some_not_done = False
        for sample in samples:
            is_done = sample_done[sample.sample_id]
            if is_done:
                continue

            done_file = done_file_name(data_directory, sample.sample_id)
            is_done = os.path.isfile(done_file)
            sample_done[sample.sample_id] = is_done
            some_not_done = True
        if not some_not_done:
            return #All done!

        time.sleep(10)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_id", help="Id of the sample to process")
    parser.add_argument("data_directory", help="Directory for the data (input + output)")
    parser.add_argument("log_directory", help="Directory for log files")
    parser.add_argument("star_file_path", help="path to STAR executable")
    parser.add_argument("--fastq_files", help="fastq files to use", nargs="+")

    args = parser.parse_args()
    data_directory = args.data_directory
    sample_id = args.sample_id
    log_directory = args.log_directory
    star_file_path = args.star_file_path
    fastq_files = ' '.join(args.fastq_files)

    sample_dir = os.path.join(data_directory, f"sample{sample_id}")

    # Prepare the transcriptome fasta files
    for i in [1,2]:
        genome = os.path.join(sample_dir, f"custom_genome_{i}.fa")
        geneinfo = os.path.join(sample_dir, f"updated_annotation_{i}.txt")
        transcriptome = os.path.join(sample_dir, f"transcriptome_{i}.fa")
        file_prep = GeneFilesPreparation(genome, geneinfo, transcriptome)
        file_prep.prepare_gene_files()

    # Create STAR indexes
    for i in [1,2]:
        transcriptome = os.path.join(sample_dir, f"transcriptome_{i}.fa")
        transcriptome_dir = os.path.join(sample_dir, f"transcriptome_{i}")
        os.mkdir(transcriptome_dir)
        command = f"{star_file_path} --runThreadN 4 --runMode genomeGenerate --genomeDir {transcriptome_dir} --genomeFastaFiles {transcriptome}"
        subprocess.run(command, shell=True, check=True)

    # Align fastq files to new indexes
    for i in [1,2]:
        transcriptome_dir = os.path.join(sample_dir, f"transcriptome_{i}")
        align_file_prefix = os.path.join(sample_dir, f"{i}_")
        command = f"{star_file_path} --runThreadN 4 --genomeDir {transcriptome_dir}  --readFilesIn {fastq_files} --readFilesCommand zcat --outFileNamePrefix {align_file_prefix} --outSAMunmapped Within"
        subprocess.run(command, shell=True, check=True)

    # Output to DONE file
    with open(done_file_name(data_directory, sample_id), "w") as done_file:
        done_file.write("DONE\n")
