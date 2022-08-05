from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline

if 'molecule_file' in snakemake.input.keys():
    LibraryPrepPipeline.main(
            snakemake.params.seed,
            snakemake.params.config,
            snakemake.params.global_config,
            snakemake.params.outdir,
            snakemake.params.logdir,
            snakemake.input.molecule_file,
            snakemake.wildcards.packet_num,
            sample_id = snakemake.wildcards.sample,
    )
else:
    # From distribution
    LibraryPrepPipeline.main(
            snakemake.params.seed,
            snakemake.params.config,
            snakemake.params.global_config,
            snakemake.params.outdir,
            snakemake.params.logdir,
            None,
            snakemake.wildcards.packet_num,
            sample_id = snakemake.wildcards.sample,
            distribution_directory = snakemake.input.sample_data_dir,
            molecules_per_packet_from_distribution = snakemake.params.num_molecules_per_packet,
    )
