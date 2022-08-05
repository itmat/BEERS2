import json
from beers.sequence.sequence_pipeline import SequencePipeline

config = json.loads(snakemake.params.config)

SequencePipeline.main(
        snakemake.params.seed,
        config['sequence_pipeline'],
        config,
        snakemake.input.cluster_packet,
        snakemake.output.cluster_packet,
        snakemake.params.logdir,
)
