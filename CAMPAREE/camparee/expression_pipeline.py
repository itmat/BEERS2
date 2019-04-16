import sys
import os
import importlib
import time
from beers_utils.constants import CONSTANTS,SUPPORTED_SCHEDULER_MODES, MAX_SEED
from beers_utils.job_monitor import JobMonitor
from camparee.camparee_utils import CampareeUtils, CampareeException
#To enable export of config parameter dictionaries to command line
import json
import subprocess
import inspect
import numpy

import camparee.transcriptomes as transcriptomes


class ExpressionPipeline:
    """
    This class represents a pipeline of steps that take user supplied fastq files through alignment, variants
    finding, parental genome construction, annotation, quantification and generation of transcripts and finally the
    generation of packets of molecules that may be used to simulate RNA sequencing.
    """

    CAMPAREE_ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    THIRD_PARTY_SOFTWARE_DIR_PATH = os.path.join(CAMPAREE_ROOT_DIR, "third_party_software")
    REQUIRED_RESOURCE_MAPPINGS = ['species_model', 'star_genome_index_directory_name',
                                  'reference_genome_filename', 'annotation_filename', 'chr_ploidy_filename']
    REQUIRED_OUTPUT_MAPPINGS = ['directory_path', 'type', 'override_sample_molecule_count', 'default_molecule_count']

    def __init__(self, configuration, scheduler_mode, output_directory_path, input_samples):
        self.scheduler_mode = scheduler_mode
        self.samples = input_samples
        self.output_directory_path = output_directory_path
        log_directory_path = os.path.join(output_directory_path, CONSTANTS.LOG_DIRECTORY_NAME)
        data_directory_path = os.path.join(output_directory_path, CONSTANTS.DATA_DIRECTORY_NAME)
        self.data_directory_path = data_directory_path
        self.log_directory_path = log_directory_path
        self.create_intermediate_data_subdirectories(data_directory_path, log_directory_path)
        self.log_file_path = os.path.join(log_directory_path, "expression_pipeline.log")
        self.steps = {}
        #Track pathes of scripts for each step. This is needed when running the
        #steps from the command line, as we do when submitting to lsf.
        self.__step_paths = {}
        for step, props in configuration['steps'].items():
            module_name, step_name = step.rsplit(".")
            parameters = props["parameters"] if props and "parameters" in props else None
            module = importlib.import_module(f'.{module_name}', package="camparee")
            step_class = getattr(module, step_name)
            self.steps[step_name] = step_class(log_directory_path, data_directory_path, parameters)
            self.__step_paths[step_name] = inspect.getfile(module)
            JobMonitor.PIPELINE_STEPS[step_name] = step_class

        # Validate the resources and set file and directory paths as needed.
        if not self.validate_and_set_resources(configuration['resources']):
            raise CampareeValidationException("The resources data is not completely valid.  "
                                              "Consult the standard error file for details.")

        # Collect the data from the ref genome and chromosome ploidy files
        self.reference_genome = CampareeUtils.create_genome(self.reference_genome_file_path)
        self.chr_ploidy_data = CampareeUtils.create_chr_ploidy_data(self.chr_ploidy_file_path)

        # Set 3rd party software paths
        self.set_third_party_software()

        # Validate output data and set
        if not self.validate_and_set_output_data(configuration['output']):
            raise CampareeValidationException("The output data is not completely valid.  "
                                              "Consult the standard error file for details.")

        self.expression_pipeline_monitor = None
        if self.scheduler_mode != "serial":
            self.expression_pipeline_monitor = JobMonitor(self.output_directory_path,
                                                          self.scheduler_mode)
            print(f"Running CAMPAREE using the {self.scheduler_mode} job scheduler.",
                  file=sys.stderr)
        else:
            print("Running CAMPAREE in serial mode.", file=sys.stderr)

    def create_intermediate_data_subdirectories(self, data_directory_path, log_directory_path):
        for sample in self.samples:
            os.makedirs(os.path.join(data_directory_path, f'sample{sample.sample_id}'), mode=0o0755, exist_ok=True)
            os.makedirs(os.path.join(log_directory_path, f'sample{sample.sample_id}'), mode=0o0755, exist_ok=True)

    def validate_and_set_output_data(self, output):
        """
        Helper method to validate and set output data.
        :param output: The output dictionary extracted from the configuration file.
        :return: True for valid output data and False otherwise
        """

        valid = True

        # Insure that all required resources keys are in place.  No point in continuing until this
        # problem is resolved.
        missing_output_keys = [item
                                 for item in ExpressionPipeline.REQUIRED_OUTPUT_MAPPINGS
                                 if item not in output]
        if missing_output_keys:
            print(f"The following required mappings were not found under 'outputs': "
                  f"{(',').join(missing_output_keys)}", file=sys.stderr)
            return False

        # Insure type mapping exists
        if "type" not in output:
            print(f"The required mapping 'type' was not found under 'output.", file=sys.stderr)
            valid = False
        else:
            self.output_type = output["type"]

        # Insure default_molecule_count exists and is an int
        # TODO: is this redundant given the check performed with missing_output_keys above?
        if "default_molecule_count" not in output:
            print(f"The required mapping 'default_molecule_count' was not found under 'output.", file=sys.stderr)
            valid = False
        else:
            self.default_molecule_count = output["default_molecule_count"]
            if not isinstance(self.default_molecule_count, int) and \
               self.default_molecule_count < 0:
                print(f"The 'default_molecule_count' must be a positive integer - not "
                      f"{output['default_molecule_count']}",
                      file=sys.stderr)
                valid = False

        #Check to make sure override_sample_molecule_count is a boolean
        self.override_sample_molecule_count = output["override_sample_molecule_count"]
        if not isinstance(self.override_sample_molecule_count, bool):
            print(f"The 'override_sample_molecule_count' must be a boolean (True/False) "
                  f"- not {output['override_sample_molecule_count']}",
                  file=sys.stderr)
            valid = False

        return valid

    def set_third_party_software(self):
        """
        Helper method to gather the names of all the 3rd party application files or directories and use them to set
        all the paths needed in the pipeline.  Since the third party software is shipped with this application, validation
        should not be necessary.  Software is identified generally by name and not specifically by filename since
        filenames may contain versioning and other artefacts.
        :return: the filenames for beagle, star, and kallisto, and the directory name for bowtie2.
        """

        beagle_filename = [filename for filename in os.listdir(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH)
                           if "beagle" in filename.lower()][0]
        self.beagle_file_path = os.path.join(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH, beagle_filename)

        star_filename = [filename for filename in os.listdir(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH)
                         if "STAR" in filename][0]
        self.star_file_path = os.path.join(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH, star_filename)

        kallisto_filename = [filename for filename in os.listdir(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH)
                              if "kallisto" in filename][0]
        self.kallisto_file_path = os.path.join(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH, kallisto_filename)

        bowtie2_dir_name = [filename for filename in os.listdir(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH)
                             if "bowtie2" in filename][0]
        self.bowtie2_dir_path = os.path.join(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH, bowtie2_dir_name)

    def validate_and_set_resources(self, resources):
        """
        Since the resources are input file intensive, and since information about resource paths is found in the
        configuration file, this method validates that all needed resource information is complete, consistent and all
        input data is found.
        :param resources: dictionary containing resources from the configuration file
        :return: True if valid and False otherwise
        """
        # TODO a some point STAR will be in third party software and may require validation

        valid = True

        # Insure that all required resources keys are in place.  No point in continuing until this
        # problem is resolved.
        missing_resource_keys = [item
                                 for item in ExpressionPipeline.REQUIRED_RESOURCE_MAPPINGS
                                 if item not in resources]
        if missing_resource_keys:
            print(f"The following required mappings were not found under 'resources': "
                  f"{(',').join(missing_resource_keys)}", file=sys.stderr)
            return False

        # If user did not provide a path to the resources directory, use the
        # directory contained in the CAMPAREE install path.
        resources_directory_path = resources.get('directory_path', None)
        if not resources_directory_path:
            resources_directory_path = os.path.join(ExpressionPipeline.CAMPAREE_ROOT_DIR, "resources")
        elif not(os.path.exists(resources_directory_path) and os.path.isdir(resources_directory_path)):
            print(f"The given resources directory, {resources_directory_path}, must exist as a directory.",
                  file=sys.stderr)
            return False

        # Insure that the species model directory exists.  No point in continuing util this problem is
        # resolved.
        species_model_directory_path = os.path.join(resources_directory_path, resources['species_model'])
        if not(os.path.exists(species_model_directory_path) and os.path.isdir(species_model_directory_path)):
                print(f"The species model directory, {species_model_directory_path}, must exist as a directory",
                      file=sys.stderr)
                return False

        # Insure that the reference genome file path exists and points to a file.
        self.reference_genome_file_path = os.path.join(species_model_directory_path, resources['reference_genome_filename'])
        if not(os.path.exists(self.reference_genome_file_path) and os.path.isfile(self.reference_genome_file_path)):
            print(f"The reference genome file path, {self.reference_genome_file_path}, must exist as"
                  f" a file.", file=sys.stderr)
            valid = False

        # Insure that the chromosome ploidy file path exists and points to a file.
        self.chr_ploidy_file_path = os.path.join(species_model_directory_path, resources['chr_ploidy_filename'])
        if not(os.path.exists(self.chr_ploidy_file_path) and os.path.isfile(self.chr_ploidy_file_path)):
            print(f"The chr ploidy file path, {self.chr_ploidy_file_path} must exist as a file", file=sys.stderr)
            valid = False

        # Insure that the annotations file path exists and points to a file.
        self.annotation_file_path = os.path.join(species_model_directory_path, resources['annotation_filename'])
        if not(os.path.exists(self.annotation_file_path) and os.path.isfile(self.annotation_file_path)):
            print(f"The annotation file path, {self.annotation_file_path} must exist as a file", file=sys.stderr)
            valid = False

        # Insure that the star index directory exists as a directory.
        self.star_index_directory_path = \
            os.path.join(species_model_directory_path, resources['star_genome_index_directory_name'])
        if not (os.path.exists(self.star_index_directory_path) and os.path.isdir(self.star_index_directory_path)):
            print(f"The star index directory, {self.star_index_directory_path}, must exist as a directory",
                  file=sys.stderr)
            valid = False

        return valid

    def validate(self):
        valid = True
        reference_genome_chromosomes = self.reference_genome.keys()
        ploidy_chromosomes = set(self.chr_ploidy_data.keys())
        if not ploidy_chromosomes.issubset(reference_genome_chromosomes):
            missing_chromosomes = ' '.join(chrom for chrom in ploidy_chromosomes.difference(reference_genome_chromosomes))
            reference_chroms = ' '.join(chrom for chrom in reference_genome_chromosomes)
            print(f"The chromosome ploidy has chromosomes `{missing_chromosomes}` not found in the reference genome file", file=sys.stderr)
            print(f"The reference genome has chromosomes {reference_chroms}", file=sys.stderr)
            valid = False
        if not all([step.validate() for step in self.steps.values()]):
            valid = False
        return valid

    def execute(self):
        print("Execution of the Expression Pipeline Started...")

        seeds = self.generate_job_seeds()

        bam_files = {}
        for sample in self.samples:

            #Retrieve name of BAM file associated with this sample. This is either
            #the path to a user provided BAM file, or the default path the
            #GenomeAlignmentStep will use to store the alignment results for this
            #sample.
            genome_alignment = self.steps['GenomeAlignmentStep']
            bam_file = genome_alignment.get_genome_bam_path(sample)
            bam_files[sample.sample_id] = bam_file
            self.run_step(step_name='GenomeAlignmentStep',
                          sample=sample,
                          execute_args=[sample, self.star_index_directory_path,
                                        self.star_file_path],
                          cmd_line_args=[sample, self.star_index_directory_path,
                                         self.star_file_path],
                          scheduler_memory_in_mb=40000,
                          scheduler_num_processors=4)

        for sample in self.samples:

            bam_filename = bam_files[sample.sample_id]
            self.run_step(step_name='GenomeBamIndexStep',
                          sample=sample,
                          execute_args=[sample, bam_filename],
                          cmd_line_args=[sample, bam_filename],
                          dependency_list=[f"GenomeAlignmentStep.{sample.sample_id}"])

        for sample_id, bam_file in bam_files.items():
            # Use chr_ploidy as the gold std for alignment, variants, VCF, genome_maker
            sample = self.expression_pipeline_monitor.get_sample(sample_id)
            seed = seeds[f"VariantsFinderStep.{sample_id}"]
            self.run_step(step_name='VariantsFinderStep',
                          sample=sample,
                          execute_args=[sample, bam_file, self.chr_ploidy_data,
                                        self.reference_genome, seed],
                          cmd_line_args=[sample, bam_file, self.chr_ploidy_file_path,
                                         self.reference_genome_file_path, seed],
                          dependency_list=[f"GenomeBamIndexStep.{sample_id}"])

        for sample_id, bam_file in bam_files.items():
            output_directory = os.path.join(self.data_directory_path, f"sample{sample.sample_id}")
            sample = self.expression_pipeline_monitor.get_sample(sample_id)
            self.run_step(step_name='IntronQuantificationStep',
                          sample=sample,
                          execute_args=[bam_file, output_directory, self.annotation_file_path],
                          cmd_line_args=[bam_file, output_directory, self.annotation_file_path],
                          dependency_list=[f"GenomeBamIndexStep.{sample_id}"])
            #TODO: do we need to depend upon the index being done? or just the alignment?
            #      I'm hypothesizing that some failures are being caused by indexing and quantification happening
            #      on the same BAM file at the same time, though I don't know why this would be a problem.

        seed = seeds["VariantsCompilationStep"]
        self.run_step(step_name='VariantsCompilationStep',
                      sample=None,
                      execute_args=[[sample.sample_id for sample in self.samples],
                                    self.chr_ploidy_data, self.reference_genome, seed],
                      cmd_line_args=[[sample.sample_id for sample in self.samples],
                                     self.chr_ploidy_file_path,
                                     self.reference_genome_file_path, seed],
                      dependency_list=[f"VariantsFinderStep.{sample.sample_id}" for sample in self.samples])

        seed = seeds["BeagleStep"]
        self.run_step(step_name='BeagleStep',
                      sample=None,
                      execute_args=[self.beagle_file_path, seed],
                      cmd_line_args=[self.beagle_file_path, seed],
                      dependency_list=[f"VariantsCompilationStep"])

        #TODO: We could load all of the steps in the entire pipeline into the queue
        #      and then just have the queue keep running until everything finishes.
        #      The only advantage of explicitly waiting here is that the user gets
        #      stdout indicating which stage is running for the pipeline.
        self.expression_pipeline_monitor.monitor_until_all_jobs_completed(queue_update_interval=10)

        for sample in self.samples:
            print(f"Processing sample{sample.sample_id} ({sample.sample_name}...")
            self.run_step(step_name='GenomeBuilderStep',
                          sample=sample,
                          execute_args=[sample, self.chr_ploidy_data, self.reference_genome],
                          cmd_line_args=[sample, self.chr_ploidy_file_path,
                                         self.reference_genome_file_path],
                          dependency_list=[f"BeagleStep"])

            for suffix in [1, 2]:

                self.run_step(step_name='UpdateAnnotationForGenomeStep',
                              sample=sample,
                              execute_args=[sample, suffix, self.annotation_file_path,
                                            self.chr_ploidy_file_path],
                              cmd_line_args=[sample, suffix, self.annotation_file_path,
                                             self.chr_ploidy_file_path],
                              dependency_list=[f"GenomeBuilderStep.{sample.sample_id}"],
                              jobname_suffix=suffix)

            seed = seeds[f"TranscriptQuantificatAndMoleculeGenerationStep.{sample.sample_id}"]
            num_molecules_to_generate = sample.molecule_count
            # If no molecule count specified for this sample, use the default count.
            if not num_molecules_to_generate or self.override_sample_molecule_count:
                num_molecules_to_generate = self.default_molecule_count
            self.run_step(step_name='TranscriptQuantificatAndMoleculeGenerationStep',
                          sample=sample,
                          execute_args=[sample, self.kallisto_file_path, self.bowtie2_dir_path,
                                        self.output_type, num_molecules_to_generate, seed],
                          cmd_line_args=[sample, self.kallisto_file_path, self.bowtie2_dir_path,
                                         self.output_type, num_molecules_to_generate, seed],
                          dependency_list=[f"UpdateAnnotationForGenomeStep.{sample.sample_id}.1",
                                           f"UpdateAnnotationForGenomeStep.{sample.sample_id}.2"],
                          scheduler_memory_in_mb=40000,
                          scheduler_num_processors=7)

        self.expression_pipeline_monitor.monitor_until_all_jobs_completed(queue_update_interval=10)

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

    def generate_job_seeds(self):
        """
        Generate one seed per job that needs a seed, returns a dictionary mapping
        job names to seeds

        We generate seeds for each job since they run on separate nodes of the cluster, potentially
        and so do not simply share Numpy seeds. We generate them all ahead of time so that if jobs need
        to be restart, they can reuse the same seed.
        """
        seeds = {}
        # Seeds for jobs that don't run per sample
        for job in ["VariantsCompilationStep", "BeagleStep"]:
            seeds[job] = numpy.random.randint(MAX_SEED)
        # Seeds for jobs that are run per sample
        for job in ["VariantsFinderStep", "TranscriptQuantificatAndMoleculeGenerationStep"]:
            for sample in self.samples:
                seeds[f"{job}.{sample.sample_id}"] = numpy.random.randint(MAX_SEED)
        return seeds

    def run_step(self, step_name, sample, execute_args, cmd_line_args, dependency_list=None,
                 jobname_suffix=None, scheduler_memory_in_mb=6000, scheduler_num_processors=1):
        """
        Helper function that runs the given step, with the given parameters. If
        CAMPAREE is configured to use a scheduler/job monitor, this helper function
        wraps submission of the step to the job monitor.

        Parameters
        ----------
        step_name : string
            Name of the CAMPAREE step to run. It should be in the list of steps
            stored in the steps dictionary.
        sample : Sample
            Sample to run through the step. For steps that aren't associated with
            specific samples, set this to None.
        execute_args : list
            List of positional paramteres to pass to the execute() method for the
            given step.
        cmd_line_args : list
            List of positional parameters to pass to the get_commandline_call()
            method for the given step.
        dependency_list : list
            List of job names (if any) the current step depends on. Default: None.
        jobname_suffix : string
            Suffix to add to job submission ID. Default: None.
        scheduler_memory_in_mb : int
            Amount of RAM (in MB) to request if submitting this step to a job scheduler.
        scheduler_num_processors : int
            Number of processors to request if submitting this step to a job scheduler.

        """
        if step_name not in list(self.steps.keys()):
            raise CampareeException(f"{step_name} not in the list of loaded steps (see config file).")

        step_class = self.steps[step_name]
        if self.scheduler_mode == "serial":
            status_msg = f"Performing {step_name}"
            status_msg += f".{jobname_suffix}" if jobname_suffix else ""
            status_msg += f" on sample{sample.sample_id}" if sample else ""
            print(status_msg)
            step_class.execute(*execute_args)
            status_msg = f"Finished {step_name}"
            status_msg += f".{jobname_suffix}" if jobname_suffix else ""
            status_msg += f" on sample{sample.sample_id}" if sample else ""
            print(status_msg + "\n")
        else:
            stdout_log = os.path.join(step_class.log_directory_path,
                                      f"sample{sample.sample_id}" if sample else "",
                                      f"{step_name}{f'.{jobname_suffix}' if jobname_suffix else ''}.bsub.%J.out")
            stderr_log = os.path.join(step_class.log_directory_path,
                                      f"sample{sample.sample_id}" if sample else "",
                                      f"{step_name}{f'.{jobname_suffix}' if jobname_suffix else ''}.bsub.%J.err")
            scheduler_job_name = (f"{step_name}{f'.sample{sample.sample_id}_{sample.sample_name}' if sample else ''}"
                                  f"{f'.{jobname_suffix}' if jobname_suffix else ''}")
            scheduler_args = {'job_name' : scheduler_job_name,
                              'stdout_logfile' : stdout_log,
                              'stderr_logfile' : stderr_log,
                              'memory_in_mb' : scheduler_memory_in_mb,
                              'num_processors' : scheduler_num_processors}
            command = step_class.get_commandline_call(*cmd_line_args)
            validation_attributes = step_class.get_validation_attributes(*cmd_line_args)
            output_directory = os.path.join(step_class.data_directory_path,
                                            f"sample{sample.sample_id}" if sample else "")
            self.expression_pipeline_monitor.submit_new_job(job_id=f"{step_name}{f'.{sample.sample_id}' if sample else ''}"
                                                                   f"{f'.{jobname_suffix}' if jobname_suffix else ''}",
                                                            job_command=command,
                                                            sample=sample,
                                                            step_name=step_name,
                                                            scheduler_arguments=scheduler_args,
                                                            validation_attributes=validation_attributes,
                                                            output_directory_path=output_directory,
                                                            system_id=None,
                                                            dependency_list=dependency_list)

    @staticmethod
    def main(configuration, scheduler_mode, output_directory_path, input_samples):
        pipeline = ExpressionPipeline(configuration, scheduler_mode, output_directory_path, input_samples)
        if not pipeline.validate():
            raise CampareeValidationException("Expression Pipeline Validation Failed.  "
                                              "Consult the standard error file for details.")
        pipeline.execute()


class CampareeValidationException(CampareeException):
    pass
