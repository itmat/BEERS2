from beers.library_prep.library_prep_pipeline import LibraryPrepPipeline

if 'molecule_file' in snakemake.input.keys():
    LibraryPrepPipeline.main(
            snakemake.params.seed,
            snakemake.params.config,
            snakemake.params.global_config,
            snakemake.output.packet_file,
            snakemake.output.output_quant_file,
            snakemake.output.input_quant_file,
            snakemake.output.log_paths,
            snakemake.input.molecule_file,
            snakemake.wildcards.packet_num,
            sample_id = snakemake.wildcards.sample,
            full_logs = snakemake.params.full_logs,
    )
else:
    # From distribution
    LibraryPrepPipeline.main(
            snakemake.params.seed,
            snakemake.params.config,
            snakemake.params.global_config,
            snakemake.output.packet_file,
            snakemake.output.output_quant_file,
            snakemake.output.input_quant_file,
            snakemake.output.log_paths,
            None,
            snakemake.wildcards.packet_num,
            sample_id = snakemake.wildcards.sample,
            distribution_directory = snakemake.params.sample_data_dir,
            molecules_per_packet_from_distribution = snakemake.params.num_molecules_per_packet,
            full_logs = snakemake.params.full_logs,
    )
