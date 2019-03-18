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
import argparse
import time
import pickle

from beers.personal.gene_files_preparation import GeneFilesPreparation
from beers.expression.transcript_gene_quant import TranscriptGeneQuantificationStep
from beers.expression.allelic_imbalance_quant import AllelicImbalanceQuantificationStep
from beers.expression.molecule_maker import MoleculeMaker
from beers.sample import Sample

def done_file_name(data_directory, sample_id):
    return os.path.join(data_directory, f"sample{sample_id}", "transcriptome_prep_done.txt")


def prep_transcriptomes(samples, data_directory, log_directory, kallisto_file_path, bowtie2_dir_path, dispatcher_mode="serial"):
    for sample in samples:
        command = f"python -m beers.expression.transcriptomes {sample.sample_id} {data_directory} {log_directory} {kallisto_file_path} {bowtie2_dir_path} --fastq_files {' '.join(sample.input_file_paths)}"

        if dispatcher_mode == "serial":
            subprocess.run(command, shell=True)
        elif dispatcher_mode == "lsf":
            stdout_log = os.path.join(log_directory, f"sample{sample.sample_id}", "Transcriptome.bsub.%J.out")
            stderr_log = os.path.join(log_directory, f"sample{sample.sample_id}", "Transcriptome.bsub.%J.err")
            bsub_command = (f"bsub -M 40000"
                            f" -J Prep_Transcriptomes.sample{sample.sample_id}_{sample.sample_name}"
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
    parser.add_argument("kallisto_file_path", help="path to Kallisto executable")
    parser.add_argument("bowtie2_dir_path", help="path to Bowtie2 directory")
    parser.add_argument("--fastq_files", help="fastq files to use", nargs="+")

    args = parser.parse_args()
    data_directory = args.data_directory
    sample_id = args.sample_id
    log_directory = args.log_directory
    kallisto_file_path = args.kallisto_file_path
    bowtie2_dir_path = args.bowtie2_dir_path
    fastq_file_1, fastq_file_2 = args.fastq_files

    sample_dir = os.path.join(data_directory, f"sample{sample_id}")

    # Prepare the transcriptome fasta files
    for i in [1,2]:
        print(f"Preparing fasta files for transcriptome {i} of sample{sample_id}")
        genome = os.path.join(sample_dir, f"custom_genome_{i}.fa")
        geneinfo = os.path.join(sample_dir, f"updated_annotation_{i}.txt")
        transcriptome = os.path.join(sample_dir, f"transcriptome_{i}.fa")
        file_prep = GeneFilesPreparation(genome, geneinfo, transcriptome)
        file_prep.prepare_gene_files()


    # Build Kallisto indexes
    for i in [1, 2]:
        print(f"Building Kallisto indexes for transcriptome {i} of sample{sample_id}")
        transcriptome = os.path.join(sample_dir, f"transcriptome_{i}.fa")
        transcriptome_index_dir = os.path.join(sample_dir, f"transcriptome_{i}_index")
        os.mkdir(transcriptome_index_dir)
        transcriptome_index = os.path.join(transcriptome_index_dir, f"transcriptome_{i}.kallisto.index")

        command = f"{kallisto_file_path} index -i {transcriptome_index} {transcriptome}"
        subprocess.run(command, shell=True, check=True)


    # Kallisto Quantification Step
    for i in [1, 2]:
        print(f"Running Kallisto quantification for transcriptome {i} of sample{sample_id}")
        transcriptome = os.path.join(sample_dir, f"transcriptome_{i}.fa")
        transcriptome_index_dir = os.path.join(sample_dir, f"transcriptome_{i}_index")
        transcriptome_index = os.path.join(transcriptome_index_dir, f"transcriptome_{i}.kallisto.index")

        kallisto_output_dir = os.path.join(sample_dir, f"transcriptome_{i}_kallisto_quant")
        os.mkdir(kallisto_output_dir)
        command = f"{kallisto_file_path} quant -i {transcriptome_index} -o {kallisto_output_dir} {fastq_file_1} {fastq_file_2}"
        subprocess.run(command, shell=True, check=True)


    # Build Bowtie2 indexes
    for i in [1,2]:
        print(f"Building Bowtie2 indexes for transcriptome {i} of sample{sample_id}")
        transcriptome = os.path.join(sample_dir, f"transcriptome_{i}.fa")
        bowtie2_build_file_path = os.path.join(bowtie2_dir_path, f"bowtie2-build")
        transcriptome_index_dir = os.path.join(sample_dir, f"transcriptome_{i}_index")
        transcriptome_bowtie2_index = os.path.join(transcriptome_index_dir, f"transcriptome_{i}")

        command = f"{bowtie2_build_file_path} {transcriptome} {transcriptome_bowtie2_index}"
        subprocess.run(command, shell=True, check=True)

    # Align fastq files to Bowtie2 indexes
    for i in [1, 2]:
        print(f"Aligning by Bowtie2 indexes for transcriptome {i} of sample{sample_id}")
        transcriptome = os.path.join(sample_dir, f"transcriptome_{i}.fa")
        bowtie2_file_path = os.path.join(bowtie2_dir_path, f"bowtie2")
        transcriptome_index_dir = os.path.join(sample_dir, f"transcriptome_{i}_index")
        transcriptome_bowtie2_index = os.path.join(transcriptome_index_dir, f"transcriptome_{i}")

        aligned_file_path = os.path.join(sample_dir, f"{i}_Aligned.out.sam")
        command = f"{bowtie2_file_path} -x {transcriptome_bowtie2_index} -1 {fastq_file_1} -2 {fastq_file_2} -S {aligned_file_path}"
        subprocess.run(command, shell=True, check=True)
        
    # Run transcript/gene/allelic_imbalance quantification scripts
    geneinfo_filename = os.path.join(sample_dir, "updated_annotation_1.txt")

    print(f"Quantifying transcriptome reads for sample{sample_id}")
    transcript_gene_quant = TranscriptGeneQuantificationStep(geneinfo_filename, sample_dir)
    #transcript_gene_quant.quantify_transcript()
    transcript_gene_quant.make_transcript_gene_dist_file()

    print(f"Quantifying allelic imabalance for sample{sample_id}")
    allelic_imbalance_quant = AllelicImbalanceQuantificationStep(sample_dir)
    allelic_imbalance_quant.quantify_allelic_imbalance()
    allelic_imbalance_quant.make_allele_imbalance_dist_file()

    # Run Molecule maker to generate a packet
    # TODO: what sample object can we use?
    print(f"Generating molecules for sample{sample_id}")
    sample = Sample(sample_id, "sample{sample_id}", [sample_dir], ["unkown", "unknown"])
    molecule_maker = MoleculeMaker(sample, sample_dir)
    packet = molecule_maker.make_packet()

    with open(os.path.join(sample_dir, "molecule_packet.pickle"), "wb") as out_file:
        pickle.dump(packet, out_file)

    print(f"Done with transcriptome.py for sample{sample_id}")

    # Output to DONE file
    with open(done_file_name(data_directory, sample_id), "w") as done_file:
        done_file.write("DONE\n")
