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
        config['flowcell'],
)

barcodes = {sample: (config['samples'][sample]['barcodes']['i5'] + "+" + config['samples'][sample]['barcodes']['i7'])}

print(dir(snakemake))
if filetype in ['sam', 'bam']:
    SAM = beers.sam.SAM(
        flowcell = flowcell,
        cluster_packet_directory = snakemake.params.cluster_packet_dir,
        sam_output_directory = snakemake.params.outdir,
        sample_id = sample,
        sample_barcodes = barcodes,
    )
    reference_genome_file = snakemake.input.reference_genome
    SAM.generate_report(
        reference_seqs = read_fasta(reference_genome_file),
        BAM = (filetype == 'bam'),
        sort_by_coordinates = False, # TODO: allow configuring sorted sam/bam files
    )
elif filetype == 'fastq':
    FastQ = beers.fast_q.FastQ(
        flowcell = flowcell,
        cluster_packet_directory = snakemake.params.cluster_packet_dir,
        fastq_output_directory = snakemake.params.outdir,
        sample_id = sample,
        sample_barcodes = barcodes,
    )
    FastQ.generate_report(
        sort_by_coordinates = False, # TODO: allow sorting fastq?
    )
else:
    raise ValueError(f"Unknown sequencer output filetype of {filetype}. Should be one of bam/sam/fastq.")
