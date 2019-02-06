import sys
import os
import importlib
import time
from beers.constants import CONSTANTS,SUPPORTED_DISPATCHER_MODES
from beers.job_monitor import Monitor
from beers.utilities.expression_utils import ExpressionUtils


class ExpressionPipeline:
    def __init__(self, configuration, dispatcher_mode, resources, output_directory_path, input_samples):
        self.dispatcher_mode = dispatcher_mode
        self.samples = input_samples
        self.output_directory_path = output_directory_path
        log_directory_path = os.path.join(output_directory_path, CONSTANTS.LOG_DIRECTORY_NAME)
        data_directory_path = os.path.join(output_directory_path, CONSTANTS.DATA_DIRECTORY_NAME)
        self.create_intermediate_data_subdirectories(data_directory_path, log_directory_path)
        self.log_file_path = os.path.join(log_directory_path, "expression_pipeline.log")
        self.steps = {}
        for step in configuration['steps']:
            if "step_name" not in step:
                raise ExpressionPipelineValidationException("Every step in the configuration must have an associated"
                                                            " step name written as 'module name.class_name'.")
            module_name, step_name = step["step_name"].rsplit(".")
            parameters = step.get("parameters", dict())
            module = importlib.import_module(f'.{module_name}', package="beers.expression")
            step_class = getattr(module, step_name)
            self.steps[step_name] = step_class(log_directory_path, data_directory_path, parameters)
        valid, reference_genome_file_path, chr_ploidy_file_path, beagle_file_path, annotation_file_path, star_file_path, resources_index_files_directory_path =\
            self.validate_and_gather_resources(resources)
        if not valid:
            raise ExpressionPipelineValidationException("The resources data is not completely valid."
                                                        "  Consult the standard error file for details.")
        self.reference_genome = ExpressionUtils.create_genome(reference_genome_file_path)
        self.chr_ploidy_data = ExpressionUtils.create_chr_ploidy_data(chr_ploidy_file_path)
        self.beagle_file_path = beagle_file_path
        self.annotation_file_path = annotation_file_path
        self.star_file_path = star_file_path
        self.reference_genome_file_path = reference_genome_file_path
        self.resources_index_files_directory_path = resources_index_files_directory_path

    def create_intermediate_data_subdirectories(self, data_directory_path, log_directory_path):
        for sample in self.samples:
            os.makedirs(os.path.join(data_directory_path, f'sample{sample.sample_id}'), mode=0o0755, exist_ok=True)
            os.makedirs(os.path.join(log_directory_path, f'sample{sample.sample_id}'), mode=0o0755, exist_ok=True)

    @staticmethod
    def validate_and_gather_resources(resources):
        """
        Since the resources are input file intensive, and since information about resource paths is found in the
        configuration file, this method validates that all needed resource information is complete, consistent and all
        input data is found.
        :param resources: dictionary containing resources from the configuration file
        :return: a tuple = valid (True/False), reference genome file path, chr ploidy file path and beagle file path
        """
        # TODO a some point STAR and samtools will be in thrid party software and may require validation
        reference_genome_file_path, chr_ploidy_file_path, beagle_file_path, annotation_file_path = \
            None, None, None, None
        valid = True
        if 'species_model' not in resources:
            print("The species_model must be listed in the resources section of the configuration file.",
                  file=sys.stderr)
            valid = False
        else:
            resources_index_files_directory_path = \
                os.path.join(resources['resources_folder'], "index_files", resources['species_model'])
            if not(os.path.exists(resources_index_files_directory_path) and
                   os.path.isdir(resources_index_files_directory_path)):
                print(f"The index files directory, {resources_index_files_directory_path}, must exist as a directory",
                      file=sys.stderr)
                valid = False
            else:
                if 'reference_genome_filename' not in resources:
                    print(f"The reference genome filename must be specified in the resources section of the "
                          f"configuration file.", file=sys.stderr)
                    valid = False
                else:
                    reference_genome_file_path = \
                        os.path.join(resources_index_files_directory_path, resources['reference_genome_filename'])
                    if not(os.path.exists(reference_genome_file_path) and
                           os.path.isfile(reference_genome_file_path)):
                        print(f"The reference genome file path, {reference_genome_file_path}, must exist as"
                              f" a file.", file=sys.stderr)
                        valid = False
                if 'chr_ploidy_filename' not in resources:
                    print(f"The chr ploidy filename must be specified in the resources section of the configuration"
                          f" file.", file=sys.stderr)
                    valid = False
                else:
                    chr_ploidy_file_path = \
                        os.path.join(resources_index_files_directory_path, resources['chr_ploidy_filename'])
                    if not(os.path.exists(chr_ploidy_file_path) and
                           os.path.isfile(chr_ploidy_file_path)):
                        print(f"The chr ploidy file path, {chr_ploidy_file_path} must exist as a file",
                              file=sys.stderr)
                        valid = False
                if 'annotation_filename' not in resources:
                    print(f"The annotation filename must be specified in the resources section of the configuration"
                          f" file.", file=sys.stderr)
                    valid = False
                else:
                    annotation_file_path = \
                        os.path.join(resources_index_files_directory_path, resources['annotation_filename'])
                    if not(os.path.exists(annotation_file_path) and
                           os.path.isfile(annotation_file_path)):
                        print(f"The annotation file path, {annotation_file_path} must exist as a file",
                              file=sys.stderr)
                        valid = False
        third_party_software_directory_path = os.path.join(resources['resources_folder'], "third_party_software")
        if not (os.path.exists(third_party_software_directory_path) and
                os.path.isdir(third_party_software_directory_path)):
                    print(f"The third party software directory path , {third_party_software_directory_path}, must exist"
                          f" as a directory.", file=sys.stderr)
                    valid = False
        else:
            beagle_filenames = [filename for filename in os.listdir(third_party_software_directory_path)
                                if "beagle" in filename]

            if not beagle_filenames:
                print(f"No file is the third party software directory can be identified as the Beagle program",
                      file=sys.stderr)
                valid = False
            else:
                beagle_file_path = os.path.join(third_party_software_directory_path, beagle_filenames[0])
                if not (os.path.exists(beagle_file_path) and os.path.isfile(beagle_file_path) and
                        beagle_filenames[0].endswith('jar')):
                    print(f"The beagle file path, {beagle_file_path}, must exist as an executable jar file.",
                          file=sys.stderr)
                    valid = False

            star_filenames = [filename for filename in os.listdir(third_party_software_directory_path)
                                if "STAR" in filename]

            if not star_filenames:
                print(f"No file is the third party software directory can be identified as the STAR program",
                      file=sys.stderr)
                valid = False
            else:
                star_file_path = os.path.join(third_party_software_directory_path, star_filenames[0])
                if not (os.path.exists(star_file_path) and os.path.isfile(star_file_path)):
                    print(f"The star file path, {star_file_path}, must exist as an executable",
                          file=sys.stderr)
                    valid = False

        return valid, reference_genome_file_path, chr_ploidy_file_path, beagle_file_path, annotation_file_path, star_file_path, resources_index_files_directory_path

    def validate(self):
        valid = True
        reference_genome_chromosomes = self.reference_genome.keys()
        ploidy_chromosomes = set(self.chr_ploidy_data.keys())
        if not ploidy_chromosomes.issubset(reference_genome_chromosomes):
            missing_chromosomes = ' '.join(chrom for chrom in ploidy_chromosomes.difference(reference_genome_chromosomes))
            reference_chroms = ' '.join(chrom for chrom in reference_genome_chromosomes.keys())
            print(f"The chromosome ploidy has chromosomes `{missing_chromosomes}` not found in the reference genome file", file=sys.stderr)
            print(f"The reference genome has chromosomes {reference_chroms}", file=sys.stderr)
            valid = False
        if not all([step.validate() for step in self.steps.values()]):
            valid = False
        return valid

    def execute(self):
        print("Execution of the Expression Pipeline Started...")

        expression_pipeline_monitor = Monitor(self.output_directory_path, self.dispatcher_mode)

        genome_alignment = self.steps['GenomeAlignmentStep']

        bam_files = []
        for sample in self.samples:
            (bam_file, system_id) = genome_alignment.execute(sample, self.resources_index_files_directory_path,
                                                             self.star_file_path, self.dispatcher_mode)
            bam_files.append(bam_file)
            expression_pipeline_monitor.submit_new_job(f'GenomeAlignment.{sample.sample_id}',
                                                       sample, "GenomeAlignmentStep",
                                                       self.output_directory_path, system_id)

        for sample in self.samples:
            #system_id = genome_alignment.index(sample, self.dispatcher_mode)
            expression_pipeline_monitor.submit_new_job(f'GenomeBamIndex.{sample.sample_id}',
                                                       sample, "GenomeBamIndexStep",
                                                       self.output_directory_path, system_id=None,
                                                       dependency_list=[f'GenomeAlignment.{sample.sample_id}'])

        #TODO: Create a generalized job submitter framework for the steps in
        #      this pipeline. Given a sample object , and the step name, we should
        #      be able to call a single function which will perform the step-
        #      specific code to launch each job and return a system id (or
        #      perhaps even handles queuing to the job monitor itself). This will
        #      cut down on the code redundancy between the rest of the script and
        #      the loop below.
        #Wait here until all of the preceding steps have finished. Submit any
        #jobs that were waiting on dependencies to complete and resubmit any
        #failed jobs.
        print(f'Waiting until all samples finish processing before running Beagle.')
        while not expression_pipeline_monitor.is_processing_complete():
            #Check for jobs requiring resubmission
            resubmission_jobs = expression_pipeline_monitor.resubmission_list.copy()
            if resubmission_jobs:
                print(f'Resubmitting {len(resubmission_jobs)} jobs that failed/stalled.')
                for resub_job_id, resub_job in resubmission_jobs.items():
                    resub_sample = expression_pipeline_monitor.get_sample(resub_job.sample_id)

                    #Generate bam file name from sample info using same rules as the
                    #execute code from genome_alignment.
                    bam_filename = os.path.join(genome_alignment.data_directory_path,
                                                f'sample{resub_sample.sample_id}',
                                                f'genome_alignment.{genome_alignment.star_bamfile_suffix}')

                    if resub_job.step_name == "GenomeAlignmentStep":
                        (bam_file, system_id) = genome_alignment.execute(resub_sample,
                                                                         self.resources_index_files_directory_path,
                                                                         self.star_file_path, self.dispatcher_mode)
                    elif resub_job.step_name == "GenomeBamIndexStep":
                        system_id = genome_alignment.index(resub_sample, bam_filename, self.dispatcher_mode)
                    """
                    elif resub_job.step_name == "VariantsFinderStep":
                        print(f'Processing variants in sample {resub_sample.sample_id} ({resub_sample.sample_name})...')
                        # Use chr_ploidy as the gold std for alignment, variants, VCF, genome_maker
                        variants_finder = self.steps['VariantsFinderStep']
                        system_id = variants_finder.execute(resub_sample, bam_file,
                                                            self.chr_ploidy_data,
                                                            self.reference_genome)
                    """
                    expression_pipeline_monitor.resubmit_job(resub_job_id, system_id)
            #Check if pending jobs have satisfied their dependencies
            pending_jobs = expression_pipeline_monitor.pending_list.copy()
            if pending_jobs:
                print(f'Check {len(pending_jobs)} pending jobs for satisfied dependencies:')
                for pend_job_id, pend_job in pending_jobs.items():
                    if expression_pipeline_monitor.are_dependencies_satisfied(pend_job_id):
                        pend_sample = expression_pipeline_monitor.get_sample(pend_job.sample_id)
                        #Generate bam file name from sample info using same rules as the
                        #execute code from genome_alignment.
                        bam_filename = os.path.join(genome_alignment.data_directory_path,
                                                    f'sample{pend_sample.sample_id}',
                                                    f'genome_alignment.{genome_alignment.star_bamfile_suffix}')

                        if pend_job.step_name == "GenomeAlignmentStep":
                            (bam_file, system_id) = genome_alignment.execute(pend_sample,
                                                                             self.resources_index_files_directory_path,
                                                                             self.star_file_path, self.dispatcher_mode)
                        elif pend_job.step_name == "GenomeBamIndexStep":
                            system_id = genome_alignment.index(pend_sample, bam_filename, self.dispatcher_mode)
                        expression_pipeline_monitor.submit_pending_job(pend_job_id, system_id)
            time.sleep(60)

        for bam_file, sample in zip(bam_files, self.samples):
            print(f"Processing variants in sample {sample.sample_id} ({sample.sample_name})...")
            # Use chr_ploidy as the gold std for alignment, variants, VCF, genome_maker
            variants_finder = self.steps['VariantsFinderStep']
            variants_finder.execute(sample, bam_file, self.chr_ploidy_data, self.reference_genome)
            """
            expression_pipeline_monitor.submit_new_job(f'VariantsFinder.{sample.sample_id}',
                                                       sample, "VariantsFinderStep",
                                                       self.output_directory_path, system_id,
                                                       dependency_list=[f'GenomeBamIndex.{sample.sample_id}'])
            """

        variants_compilation = self.steps['VariantsCompilationStep']
        variants_compilation.execute(self.samples, self.chr_ploidy_data, self.reference_genome)

        print(f"Processing combined samples...")
        beagle = self.steps['BeagleStep']
        outcome = beagle.execute(self.beagle_file_path)
        if outcome != 0:
            raise ExpressionPipelineException("Beagle process failed.")

        for sample in self.samples:
            print(f"Processing sample{sample.sample_id} ({sample.sample_name}...")

            genome_builder = self.steps['GenomeBuilderStep']
            genome_builder.execute(sample, self.chr_ploidy_data, self.reference_genome)

            for suffix in [1,2]:

                annotation_updater = self.steps['UpdateAnnotationForGenomeStep']
                annotation_updater.execute(sample, suffix, self.annotation_file_path)

            #for _ in range(2):
            #    quantifier = Quantify(annotation_updates, self.alignment_filename)
            #    transcript_distributions.append(quantifier.quantify())
            #for i, item in enumerate(genomes):
            #    genome = item[0]
            #    transcript_maker = TranscriptMaker(genome, annotation_updates[i], transcript_distributions[i], self.parameters)
            #    molecules.append(transcript_maker.prepare_transcripts())
            #r_rna = RibosomalRNA(self.parameters)
            #molecules.append(r_rna.generate_rRNA_sample())
        print("Execution of the Expression Pipeline Ended")

    @staticmethod
    def main(configuration, dispatcher_mode, resources, output_directory_path, input_samples):
        pipeline = ExpressionPipeline(configuration, dispatcher_mode, resources, output_directory_path, input_samples)
        if not pipeline.validate():
            raise ExpressionPipelineValidationException("Expression Pipeline Validation Failed.  "
                                              "Consult the standard error file for details.")
        pipeline.execute()


class ExpressionPipelineValidationException(Exception):
    pass

class ExpressionPipelineException(Exception):
    pass
