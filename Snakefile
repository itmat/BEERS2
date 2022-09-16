import pathlib
import functools
import json

if not config:
    raise Exception("No config provided - must pass `--configfile path/to/config.yaml` to snakemake to run pipeline")

samples = config['samples']
lanes =  config['flowcell']['lanes_to_use']
seed = config['seed']

@functools.cache
def input_molecule_files(sample):
    path = config['library_prep_pipeline']['input'].get('directory_path', None)
    if not path:
        return []
    input_molecule_file_dir = pathlib.Path(path)/ f"sample{sample}"
    return sorted(list(input_molecule_file_dir.glob("*")))
def try_absolute_and_relative_path(path):
    # Some resources should be allowed to be relative to the pipeline
    # since we need to specify the provided resources in the example config files
    # This function enables that by searching workflow-relative paths too
    if pathlib.Path(path).exists():
        return path
    rel_path = workflow.source_path(path) # Get workflow-relative path
    if pathlib.Path(rel_path).exists():
        return rel_path
    raise ValueError(f"No file found at location {path} either as an absolute path or relative to the BEERS2 pipeline")


# Count the numbers of packets for each sample
# by source (distribution or molecule packet) and total
num_packets_from_distribution_for_sample = {
    sample: config['library_prep_pipeline']['input']['from_distribution_data'][sample]['num_packets']
        for sample in samples.keys()
}
num_packets_from_molecule_file_for_sample = {
    sample: len(input_molecule_files(sample))
        for sample in samples.keys()
}
num_total_packets_for_sample = {
    sample: num_packets_from_distribution_for_sample[sample] + num_packets_from_molecule_file_for_sample[sample]
        for sample in samples.keys()
}

wildcard_constraints:
    packet_num = "[0-9]+",
    sample = "[0-9]+",
    sequencer_output_filetype = "sam|bam|fastq",

rule all:
    # These are the default files that BEERS2 will create
    # All dependencies of them will be generated as well.
    input:
        expand(
            "results/S{sample}_L{lane}.sam",
            sample = samples.keys(),
            lane = lanes,
        ),
        expand(
            "results/S{sample}_L{lane}_R1.fastq",
            sample = samples.keys(),
            lane = lanes,
        ),

# Prep directories:
lib_prep_dir = pathlib.Path("library_prep_pipeline")
lib_prep_dir.mkdir(exist_ok=True)
seq_dir = pathlib.Path("sequence_pipeline")
seq_dir.mkdir(exist_ok=True)
for sample in samples.keys():
    sample_dir = lib_prep_dir / f"sample{sample}"
    sample_dir.mkdir(exist_ok=True)
    (sample_dir / "from_molecule_files").mkdir(exist_ok=True)
    (sample_dir / "from_distribution").mkdir(exist_ok=True)
    (sample_dir / "logs").mkdir(exist_ok=True)

    sample_dir = seq_dir / f"sample{sample}"
    sample_dir.mkdir(exist_ok=True)
    (sample_dir / "logs").mkdir(exist_ok=True)

rule run_library_prep_packet_from_molecule_file:
    input:
        molecule_file = lambda wildcards: input_molecule_files(wildcards.sample)[int(wildcards.packet_num)]
    output:
        packet_file = "library_prep_pipeline/sample{sample}/from_molecule_files/library_prep_pipeline_result_molecule_pkt{packet_num}.txt",
        quant_file = "library_prep_pipeline/sample{sample}/from_molecule_files/library_prep_pipeline_result_molecule_pkt{packet_num}.quant_file",
    params:
        outdir = "library_prep_pipeline/sample{sample}/from_molecule_files/",
        logdir = "library_prep_pipeline/sample{sample}/logs/",
        config = json.dumps(config['library_prep_pipeline']),
        global_config = json.dumps(config),
        seed = seed,
    resources:
        mem_mb = 12_000
    script:
        "scripts/run_library_prep_pipeline.py"

rule run_library_prep_packet_from_distribution:
    input:
        sample_data_dir = lambda wildcards: config['library_prep_pipeline']['input']['from_distribution_data'][wildcards.sample]['sample_data_directory'],
    output:
        packet_file ="library_prep_pipeline/sample{sample}/from_distribution/library_prep_pipeline_result_molecule_pkt{packet_num}.txt",
        quant_file = "library_prep_pipeline/sample{sample}/from_distribution/library_prep_pipeline_result_molecule_pkt{packet_num}.quant_file",
    params:
        outdir = "library_prep_pipeline/sample{sample}/from_distribution/",
        logdir = "library_prep_pipeline/sample{sample}/logs/",
        num_molecules_per_packet = lambda wildcards: config['library_prep_pipeline']['input']['from_distribution_data'][wildcards.sample]['num_molecules_per_packet'],
        config = json.dumps(config['library_prep_pipeline']),
        global_config = json.dumps(config),
        seed = seed,
    resources:
        mem_mb = 12_000
    script:
        "scripts/run_library_prep_pipeline.py"

def get_lib_prep_result_packet_path(sample, packet_num):
    from_molecule_files = num_packets_from_molecule_file_for_sample[sample]
    if int(packet_num) < from_molecule_files:
        return f"library_prep_pipeline/sample{sample}/from_molecule_files/library_prep_pipeline_result_molecule_pkt{packet_num}.txt"
    else:
        return f"library_prep_pipeline/sample{sample}/from_distribution/library_prep_pipeline_result_molecule_pkt{packet_num}.txt"

rule create_cluster_packet:
    input:
        packet = lambda wildcards: get_lib_prep_result_packet_path(wildcards.sample, wildcards.packet_num)
    output:
        "sequence_pipeline/sample{sample}/input_cluster_packets/cluster_packet_start_pkt{packet_num}.gzip",
    params:
        outdir = "sequence_pipeline/",
        configuration = json.dumps(config),
    script:
        "scripts/create_cluster_packets.py"

# TODO: this has been commented out for now but this is to to perform an aggregation
# of all the packets to distribute all the flowcell coordinates so that they are unique
# However this considerably slows down execution as its a big bottleneck that is not
# parallelizable. We should later allow both this and the rule above to operate, choosing
# them based off the config settings
#rule create_cluster_packets_distinct_coords:
#    input:
#        packet_files_from_molecule_files = [f"library_prep_pipeline/sample{sample}/from_molecule_files/library_prep_pipeline_result_molecule_pkt{packet_num}.txt"
#            for sample in samples.keys()
#            for packet_num in range(num_packets_from_molecule_file_for_sample[sample])
#        ],
#        packet_files_from_distribution = [f"library_prep_pipeline/sample{sample}/from_distribution/library_prep_pipeline_result_molecule_pkt{packet_num}.txt"
#            for sample in samples.keys()
#            # numbered starting where the molecule file ones left off - no repeated sample numbers
#            for packet_num in range(num_packets_from_molecule_file_for_sample[sample], num_total_packets_for_sample[sample])
#        ],
#    output:
#        cluster_packets = [f"sequence_pipeline/sample{sample}/input_cluster_packets/cluster_packet_start_pkt{packet_num}.gzip"
#                                for sample in samples.keys()
#                                for packet_num in range(num_total_packets_for_sample[sample])]
#    resources:
#        mem_mb = 12_000
#    params:
#        outdir = "sequence_pipeline/",
#        parameters = json.dumps(config['flowcell']),
#    script:
#        "scripts/create_cluster_packets.py"

rule sequence_cluster_packet:
    input:
        cluster_packet = "sequence_pipeline/sample{sample}/input_cluster_packets/cluster_packet_start_pkt{packet_num}.gzip"
    output:
        cluster_packet = "sequence_pipeline/sample{sample}/output_cluster_packets/sequence_cluster_packet{packet_num}.gzip"
    params:
        seed = seed,
        config = json.dumps(config),
        logdir = "sequence_pipeline/sample{sample}/logs/",
    resources:
        mem_mb = 12_000
    script:
        'scripts/run_sequence_pipeline.py'

rule create_sequencer_outputs_sam_or_bam:
    input:
        cluster_packets = lambda wildcards: [f"sequence_pipeline/sample{{sample}}/output_cluster_packets/sequence_cluster_packet{packet_num}.gzip"
                                                for packet_num in range(num_total_packets_for_sample[wildcards.sample])],
        reference_genome = try_absolute_and_relative_path(config['resources']['reference_genome_fasta']),
    output:
        expand("results/S{{sample}}_L{lane}.{{sequencer_output_filetype}}", lane=lanes),
    params:
        cluster_packet_dir = "sequence_pipeline/sample{sample}/output_cluster_packets",
        outdir = "results/",
    resources:
        mem_mb = 6_000
    script:
        'scripts/create_sequencer_outputs.py'

rule create_sequencer_outputs_fastq:
    input:
        cluster_packets = lambda wildcards: [f"sequence_pipeline/sample{{sample}}/output_cluster_packets/sequence_cluster_packet{packet_num}.gzip"
                                                for packet_num in range(num_total_packets_for_sample[wildcards.sample])],
        reference_genome = try_absolute_and_relative_path(config['resources']['reference_genome_fasta']),
    output:
        # NOTE: fastq outputs have both R1 and R2 files to create
        expand("results/S{{sample}}_L{lane}_R1.{{sequencer_output_filetype}}", lane=lanes),
        expand("results/S{{sample}}_L{lane}_R2.{{sequencer_output_filetype}}", lane=lanes),
    params:
        cluster_packet_dir = "sequence_pipeline/sample{sample}/output_cluster_packets",
        outdir = "results/",
    resources:
        mem_mb = 6_000
    script:
        'scripts/create_sequencer_outputs.py'
