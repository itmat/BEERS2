import pathlib
import shutil
import uuid
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

    # Construct a sample directory out of the input files
    # In case files are coming from a temporary download when
    # run in a no-shared-fs mode, they won't necessarily be in
    # the right place.
    tmp_dir = pathlib.Path("/tmp/") / str(uuid.uuid4())
    tmp_dir.mkdir()
    for file_path in snakemake.input.sample_data_files:
        file_path = pathlib.Path(file_path)
        new_path = tmp_dir / file_path.name
        new_path.symlink_to(file_path.resolve())

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
            distribution_directory = tmp_dir,
            molecules_per_packet_from_distribution = snakemake.params.num_molecules_per_packet,
            full_logs = snakemake.params.full_logs,
    )

    shutil.rmtree(tmp_dir)
