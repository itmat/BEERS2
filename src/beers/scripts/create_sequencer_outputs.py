'''
Create the sequencer output files (fastq/sam/bam)
by gathering all the sequenced cluster packets together
'''

import beers.sam
import beers.fast_q
import beers.flowcell
from beers_utils.read_fasta import read_fasta

filetype = snakemake.wildcards.sequencer_output_filetype
config = snakemake.config
sample = snakemake.wildcards.sample

flowcell = beers.flowcell.Flowcell(
        config['sequence_pipeline']['flowcell'],
        rng = None, # We don't use any of the rng features here, that's only coordinate assignment
)

barcode = (config['global']['samples'][sample]['barcodes']['i5'] + "+" + config['global']['samples'][sample]['barcodes']['i7'])

if filetype in ['sam', 'bam']:
    SAM = beers.sam.SAM(
        flowcell = flowcell,
        sample_id = sample,
        sample_barcode = barcode,
    )
    reference_genome_file = snakemake.input.reference_genome
    SAM.generate_report(
        cluster_packet_paths = snakemake.input.cluster_packets,
        output_file_paths = snakemake.output.sam_file_paths,
        bad_barcode_file_paths = snakemake.output.bad_file_paths,
        reference_seqs = read_fasta(reference_genome_file),
        BAM = (filetype == 'bam'),
        sort_by_coordinates = False, # TODO: allow configuring sorted sam/bam files
    )
elif filetype == 'fastq':
    FastQ = beers.fast_q.FastQ(
        flowcell = flowcell,
        sample_id = sample,
        sample_barcode = barcode,
    )
    FastQ.generate_report(
        cluster_packet_paths = snakemake.input.cluster_packets,
        output_file_paths = [snakemake.output.output_file_paths_R1, snakemake.output.output_file_paths_R2],
        bad_barcode_file_paths = [snakemake.output.bad_file_paths_R1, snakemake.output.bad_file_paths_R2],
        sort_by_coordinates = False, # TODO: allow sorting fastq?
    )
else:
    raise ValueError(f"Unknown sequencer output filetype of {filetype}. Should be one of bam/sam/fastq.")
