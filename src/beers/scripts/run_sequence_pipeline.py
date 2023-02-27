import json
from beers.sequence.sequence_pipeline import SequencePipeline

config = json.loads(snakemake.params.config)
global_config = json.loads(snakemake.params.global_config)

SequencePipeline.main(
        snakemake.params.seed,
        config,
        global_config,
        snakemake.input.cluster_packet,
        snakemake.output.cluster_packet,
        snakemake.output.log_paths,
)
