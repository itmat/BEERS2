import sys
import os
import importlib
import time
from beers_utils.constants import CONSTANTS,SUPPORTED_DISPATCHER_MODES, MAX_SEED
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
    def __init__(self, configuration, dispatcher_mode, resources, output_directory_path, input_samples):
        self.dispatcher_mode = dispatcher_mode
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
        for step in configuration['steps']:
            if "step_name" not in step:
                raise CampareeValidationException("Every step in the configuration must have an associated"
                                                  " step name written as 'module name.class_name'.")
            module_name, step_name = step["step_name"].rsplit(".")
            parameters = step.get("parameters", dict())
            module = importlib.import_module(f'.{module_name}', package="camparee")
            step_class = getattr(module, step_name)
            self.steps[step_name] = step_class(log_directory_path, data_directory_path, parameters)
            self.__step_paths[step_name] = inspect.getfile(module)
            JobMonitor.PIPELINE_STEPS[step_name] = step_class
        valid, reference_genome_file_path, chr_ploidy_file_path, beagle_file_path, annotation_file_path, star_file_path, \
            kallisto_file_path, bowtie2_dir_path, resources_index_files_directory_path = self.validate_and_gather_resources(resources)
        if not valid:
            raise CampareeValidationException("The resources data is not completely valid."
                                              "  Consult the standard error file for details.")
        self.reference_genome = CampareeUtils.create_genome(reference_genome_file_path)
        self.chr_ploidy_data = CampareeUtils.create_chr_ploidy_data(chr_ploidy_file_path)
        self.beagle_file_path = beagle_file_path
        self.annotation_file_path = annotation_file_path
        self.star_file_path = star_file_path
        self.kallisto_file_path = kallisto_file_path
        self.bowtie2_dir_path = bowtie2_dir_path
        self.chr_ploidy_file_path = chr_ploidy_file_path
        self.reference_genome_file_path = reference_genome_file_path
        self.resources_index_files_directory_path = resources_index_files_directory_path

        self.output_type = configuration["output"]["type"]
        self.output_molecule_count = configuration["output"]["molecule_count"]

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
        reference_genome_file_path, chr_ploidy_file_path, annotation_file_path, beagle_file_path, \
            star_file_path, kallisto_file_path, bowtie2_dir_path = None, None, None, None, None, None, None
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


            kallisto_filenames = [filename for filename in os.listdir(third_party_software_directory_path)
                                if "kallisto" in filename]
            if not kallisto_filenames:
                print(f"No file is the third party software directory can be identified as the Kallisto program",
                      file=sys.stderr)
                valid = False
            else:
                kallisto_file_path = os.path.join(third_party_software_directory_path, kallisto_filenames[0])
                if not (os.path.exists(kallisto_file_path) and os.path.isfile(kallisto_file_path)):
                    print(f"The Kallisto file path, {kallisto_file_path}, must exist as an executable",
                          file=sys.stderr)
                    valid = False


            bowtie2_dir_names = [filename for filename in os.listdir(third_party_software_directory_path)
                                if "bowtie2" in filename]
            if not bowtie2_dir_names:
                print(f"No file is the third party software directory can be identified as the Bowtie2 program",
                      file=sys.stderr)
                valid = False
            else:
                bowtie2_dir_path = os.path.join(third_party_software_directory_path, bowtie2_dir_names[0])
                # Make sure build and run files are there in the directory
                if not (os.path.exists(bowtie2_dir_path)):
                    print(f"The Bowtie2 directory path, {bowtie2_dir_path}, must exist as an executable",
                          file=sys.stderr)
                    valid = False

        return valid, reference_genome_file_path, chr_ploidy_file_path, beagle_file_path, annotation_file_path, star_file_path, \
            kallisto_file_path, bowtie2_dir_path, resources_index_files_directory_path

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

        expression_pipeline_monitor = JobMonitor(self.output_directory_path, self.dispatcher_mode)

        genome_alignment = self.steps['GenomeAlignmentStep']

        bam_files = {}
        for sample in self.samples:
            (bam_file, system_id) = genome_alignment.execute(sample, self.resources_index_files_directory_path,
                                                             self.star_file_path, self.dispatcher_mode)
            bam_files[sample.sample_id] = bam_file
            expression_pipeline_monitor.submit_new_job(job_id=f"GenomeAlignment.{sample.sample_id}",
                                                       job_command="", sample=sample,
                                                       step_name='GenomeAlignmentStep',
                                                       scheduler_arguments=None,
                                                       validation_attributes=None,
                                                       output_directory_path=self.output_directory_path,
                                                       system_id=system_id,
                                                       dependency_list=None)
            #TODO: This is a hack to handle the cases where the bam files were
            #      already present, or the script was run in serial/multicore mode.
            #      Need to change this to something more legit.
            if system_id == "ALREADY_ALIGNED" or self.dispatcher_mode != "lsf":
                expression_pipeline_monitor.mark_job_completed(f"GenomeAlignment.{sample.sample_id}")

        for sample in self.samples:
            if self.dispatcher_mode == "lsf":
                expression_pipeline_monitor.submit_new_job(job_id=f"GenomeBamIndex.{sample.sample_id}",
                                                           job_command="", sample=sample,
                                                           step_name='GenomeBamIndexStep',
                                                           scheduler_arguments=None,
                                                           validation_attributes=None,
                                                           output_directory_path=self.output_directory_path,
                                                           system_id=None,
                                                           dependency_list=[f"GenomeAlignment.{sample.sample_id}"])
            else:
                print(f"Running Bam Index creating in sample {sample.sample_id}")
                bam_filename = bam_files[sample.sample_id]
                genome_alignment.index(sample, bam_filename, self.dispatcher_mode)

        for sample_id, bam_file in bam_files.items():
            # Use chr_ploidy as the gold std for alignment, variants, VCF, genome_maker
            variants_finder = self.steps['VariantsFinderStep']
            sample = expression_pipeline_monitor.get_sample(sample_id)
            seed = seeds[f"variant_finder.{sample_id}"]

            if self.dispatcher_mode == "serial":
                print(f"Processing variants in sample {sample_id}...")
                variants_finder.execute(sample, bam_file, self.chr_ploidy_data, self.reference_genome, seed)
            else:
                stdout_log = os.path.join(variants_finder.log_directory_path, f"sample{sample_id}", "Variants_Finder.bsub.%J.out")
                stderr_log = os.path.join(variants_finder.log_directory_path, f"sample{sample_id}", "Variants_Finder.bsub.%J.err")
                command = variants_finder.get_commandline_call(sample, bam_file, self.chr_ploidy_file_path,
                                                               self.reference_genome_file_path, seed)
                scheduler_args = {'job_name' : f"Variant_Finder.sample{sample_id}_{sample.sample_name}",
                                  'stdout_logfile' : stdout_log,
                                  'stderr_logfile' : stderr_log,
                                  'memory_in_mb' : 6000,
                                  'num_processors' : 1}
                validation_attributes = variants_finder.get_validation_attributes(sample)
                expression_pipeline_monitor.submit_new_job(job_id=f"VariantsFinder.{sample_id}",
                                                           job_command=command, sample=sample,
                                                           step_name='VariantsFinderStep',
                                                           scheduler_arguments=scheduler_args,
                                                           validation_attributes=validation_attributes,
                                                           output_directory_path=self.output_directory_path,
                                                           system_id=None,
                                                           dependency_list=[f"GenomeBamIndex.{sample_id}"])

        for sample_id, bam_file in bam_files.items():
            output_directory = os.path.join(self.output_directory_path, f"sample{sample_id}")
            sample = expression_pipeline_monitor.get_sample(sample_id)
            if self.dispatcher_mode == "lsf":
                expression_pipeline_monitor.submit_new_job(job_id=f"IntronQuantificationStep.{sample_id}",
                                                           job_command="", sample=sample,
                                                           step_name='IntronQuantificationStep',
                                                           scheduler_arguments=None,
                                                           validation_attributes=None,
                                                           output_directory_path=output_directory,
                                                           system_id=None,
                                                           dependency_list=[f"GenomeBamIndex.{sample_id}"])
                #TODO: do we need to depend upon the index being done? or just the alignment?
                #      I'm hypothesizing that some failures are being caused by indexing and quantification happening
                #      on the same BAM file at the same time, though I don't know why this would be a problem.
            else:
                print(f"Computing Intron Quantifications in sample {sample_id}")
                intron_quant = self.steps["IntronQuantificationStep"]
                output_directory = os.path.join(intron_quant.data_directory_path, f"sample{sample.sample_id}")
                intron_quant.execute(bam_file, output_directory, self.annotation_file_path)


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
                print(f"Resubmitting {len(resubmission_jobs)} jobs that failed/stalled.")
                for resub_job_id, resub_job in resubmission_jobs.items():
                    resub_sample = expression_pipeline_monitor.get_sample(resub_job.sample_id)

                    bam_filename = bam_files[resub_job.sample_id]

                    if resub_job.step_name == "GenomeAlignmentStep":
                        (bam_file, system_id) = genome_alignment.execute(resub_sample,
                                                                         self.resources_index_files_directory_path,
                                                                         self.star_file_path, self.dispatcher_mode)
                    elif resub_job.step_name == "GenomeBamIndexStep":
                        system_id = genome_alignment.index(resub_sample, bam_filename, self.dispatcher_mode)
                    elif resub_job.step_name == "VariantsFinderStep":
                        print(f"Submitting variant finder command to {self.dispatcher_mode} for sample {resub_sample.sample_name}.")
                        # Use chr_ploidy as the gold std for alignment, variants, VCF, genome_maker
                        variants_finder = self.steps['VariantsFinderStep']
                        variant_finder_path = self.__step_paths['VariantsFinderStep']
                        stdout_log = os.path.join(variants_finder.log_directory_path, f"sample{resub_sample.sample_id}", "Variants_Finder.bsub.%J.out")
                        stderr_log = os.path.join(variants_finder.log_directory_path, f"sample{resub_sample.sample_id}", "Variants_Finder.bsub.%J.err")

                        #TODO: Modify code to maintain access to step parameters.
                        #      Need these parameters to instantiate these objects
                        #      inside the main methods (like the variant_finder.py
                        #      script below). If I store this somewhere, it means
                        #      I won't have to manually recreate it below (which
                        #      we'll need to update here in the code every time we
                        #      update the config file).

                        #Recreate parameter dictionary for VariantsFinderStep
                        variant_finder_params = {}
                        variant_finder_params['sort_by_entropy'] = variants_finder.entropy_sort
                        variant_finder_params['min_threshold'] = variants_finder.min_abundance_threshold

                        seed = seeds[f"variant_finder.{resub_sample.sample_id}"]

                        bsub_command = (f"bsub"
                                        f" -J Variant_Finder.sample{resub_sample.sample_id}_{resub_sample.sample_name}"
                                        f" -oo {stdout_log}"
                                        f" -eo {stderr_log}"
                                        f" python {variant_finder_path}"
                                        f" --log_directory_path {variants_finder.log_directory_path}"
                                        f" --data_directory_path {variants_finder.data_directory_path}"
                                        f" --config_parameters '{json.dumps(variant_finder_params)}'"
                                        f" --sample '{repr(resub_sample)}'"
                                        f" --bam_filename {bam_filename}"
                                        f" --chr_ploidy_file_path {self.chr_ploidy_file_path}"
                                        f" --reference_genome_file_path {self.reference_genome_file_path}"
                                        f" --seed {seed}")

                        result = subprocess.run(bsub_command, shell=True, check=True, stdout = subprocess.PIPE, encoding="ascii")
                        print(f"\t{result.stdout.rstrip()}")
                        #Extract job ID from LSF stdout
                        system_id = expression_pipeline_monitor.job_scheduler.LSF_BSUB_OUTPUT_PATTERN.match(result.stdout).group('job_id')

                        print(f"Finished submitting variant finder command to {self.dispatcher_mode} for sample {resub_sample.sample_name}.")

                    elif resub_job.step_name == "IntronQuantificationStep":
                        intron_quant = self.steps["IntronQuantificationStep"]

                        stdout_log = os.path.join(intron_quant.log_directory_path, f"sample{resub_sample.sample_id}", "IntronQuantification.bsub.%J.out")
                        stderr_log = os.path.join(intron_quant.log_directory_path, f"sample{resub_sample.sample_id}", "IntronQuantification.bsub.%J.err")

                        intron_quant_path = self.__step_paths["IntronQuantificationStep"]
                        output_directory = os.path.join(intron_quant.data_directory_path, f"sample{resub_sample.sample_id}")
                        intron_quant_params = {"forward_read_is_sense": intron_quant.forward_read_is_sense, "flank_size": intron_quant.flank_size}
                        params = json.dumps(intron_quant_params)

                        bsub_command = (f"bsub"
                                        f" -J IntronQuantification.sample{resub_sample.sample_id}_{resub_sample.sample_name}"
                                        f" -oo {stdout_log}"
                                        f" -eo {stderr_log}"
                                        f" python {intron_quant_path}"
                                        f" --log_directory_path {intron_quant.log_directory_path}"
                                        f" --data_directory_path {intron_quant.data_directory_path}"
                                        f" --output_directory {output_directory}"
                                        f" --parameters '{params}'"
                                        f" --bam_file {bam_file}"
                                        f" --info_file {self.annotation_file_path}")

                        result = subprocess.run(bsub_command, shell=True, check=True, stdout = subprocess.PIPE, encoding="ascii")
                        print(f"\t{result.stdout.rstrip()}")
                        #Extract job ID from LSF stdout
                        system_id = expression_pipeline_monitor.job_scheduler.LSF_BSUB_OUTPUT_PATTERN.match(result.stdout).group('job_id')
                        print(f"Finished submitting intron quantification command to {self.dispatcher_mode} for sample {resub_sample.sample_name}")

                    # Finish resubmission
                    expression_pipeline_monitor.resubmit_job(resub_job_id, system_id)

            #Check if pending jobs have satisfied their dependencies
            pending_jobs = expression_pipeline_monitor.pending_list.copy()
            if pending_jobs:
                print(f"Check {len(pending_jobs)} pending jobs for satisfied dependencies:")
                for pend_job_id, pend_job in pending_jobs.items():
                    if expression_pipeline_monitor.are_dependencies_satisfied(pend_job_id):
                        pend_sample = expression_pipeline_monitor.get_sample(pend_job.sample_id)

                        bam_filename = bam_files[pend_job.sample_id]

                        if pend_job.step_name == "GenomeAlignmentStep":
                            (bam_file, system_id) = genome_alignment.execute(pend_sample,
                                                                             self.resources_index_files_directory_path,
                                                                             self.star_file_path, self.dispatcher_mode)
                        elif pend_job.step_name == "GenomeBamIndexStep":
                            system_id = genome_alignment.index(pend_sample, bam_filename, self.dispatcher_mode)
                        elif pend_job.step_name == "VariantsFinderStep":
                            print(f"Submitting variant finder command to {self.dispatcher_mode} for sample {pend_sample.sample_name}.")
                            # Use chr_ploidy as the gold std for alignment, variants, VCF, genome_maker
                            variants_finder = self.steps['VariantsFinderStep']
                            variant_finder_path = self.__step_paths['VariantsFinderStep']
                            stdout_log = os.path.join(variants_finder.log_directory_path, f"sample{pend_sample.sample_id}", "Variants_Finder.bsub.%J.out")
                            stderr_log = os.path.join(variants_finder.log_directory_path, f"sample{pend_sample.sample_id}", "Variants_Finder.bsub.%J.err")

                            #Recreate parameter dictionary for VariantsFinderStep
                            variant_finder_params = {}
                            variant_finder_params['sort_by_entropy'] = variants_finder.entropy_sort
                            variant_finder_params['min_threshold'] = variants_finder.min_abundance_threshold

                            seed = seeds[f"variant_finder.{pend_sample.sample_id}"]

                            bsub_command = (f"bsub"
                                            f" -J Variant_Finder.sample{pend_sample.sample_id}_{pend_sample.sample_name}"
                                            f" -oo {stdout_log}"
                                            f" -eo {stderr_log}"
                                            f" python {variant_finder_path}"
                                            f" --log_directory_path {variants_finder.log_directory_path}"
                                            f" --data_directory_path {variants_finder.data_directory_path}"
                                            f" --config_parameters '{json.dumps(variant_finder_params)}'"
                                            f" --sample '{repr(pend_sample)}'"
                                            f" --bam_filename {bam_filename}"
                                            f" --chr_ploidy_file_path {self.chr_ploidy_file_path}"
                                            f" --reference_genome_file_path {self.reference_genome_file_path}"
                                            f" --seed {seed}")

                            result = subprocess.run(bsub_command, shell=True, check=True, stdout = subprocess.PIPE, encoding="ascii")
                            print(f"\t{result.stdout.rstrip()}")
                            #Extract job ID from LSF stdout
                            system_id = expression_pipeline_monitor.job_scheduler.LSF_BSUB_OUTPUT_PATTERN.match(result.stdout).group('job_id')

                            print(f"Finished submitting variant finder command to {self.dispatcher_mode} for sample {pend_sample.sample_name}.")

                        elif pend_job.step_name == "IntronQuantificationStep":
                            intron_quant = self.steps["IntronQuantificationStep"]

                            stdout_log = os.path.join(intron_quant.log_directory_path, f"sample{pend_sample.sample_id}", "IntronQuantification.bsub.%J.out")
                            stderr_log = os.path.join(intron_quant.log_directory_path, f"sample{pend_sample.sample_id}", "IntronQuantification.bsub.%J.err")

                            print(stderr_log)

                            intron_quant_path = self.__step_paths["IntronQuantificationStep"]
                            output_directory = os.path.join(intron_quant.data_directory_path, f"sample{pend_sample.sample_id}")
                            intron_quant_params = {"forward_read_is_sense": intron_quant.forward_read_is_sense, "flank_size": intron_quant.flank_size}
                            params = json.dumps(intron_quant_params)

                            bsub_command = (f"bsub"
                                            f" -J IntronQuantification.sample{pend_sample.sample_id}_{pend_sample.sample_name}"
                                            f" -oo {stdout_log}"
                                            f" -eo {stderr_log}"
                                            f" python {intron_quant_path}"
                                            f" --log_directory_path {intron_quant.log_directory_path}"
                                            f" --data_directory_path {intron_quant.data_directory_path}"
                                            f" --output_directory {output_directory}"
                                            f" --parameters '{params}'"
                                            f" --bam_file {bam_file}"
                                            f" --info_file {self.annotation_file_path}")

                            result = subprocess.run(bsub_command, shell=True, check=True, stdout = subprocess.PIPE, encoding="ascii")
                            print(f"\t{result.stdout.rstrip()}")
                            #Extract job ID from LSF stdout
                            system_id = expression_pipeline_monitor.job_scheduler.LSF_BSUB_OUTPUT_PATTERN.match(result.stdout).group('job_id')
                            print(f"Finished submitting intron quantification command to {self.dispatcher_mode} for sample {pend_sample.sample_name}")

                        # Finish submission
                        expression_pipeline_monitor.submit_pending_job(pend_job_id, system_id)
            time.sleep(10)

        variants_compilation = self.steps['VariantsCompilationStep']
        variants_compilation.execute(self.samples, self.chr_ploidy_data, self.reference_genome)

        print(f"Processing combined samples...")
        beagle = self.steps['BeagleStep']
        outcome = beagle.execute(self.beagle_file_path, seeds["beagle"])
        if outcome != 0:
            raise CampareeException("Beagle process failed.")

        for sample in self.samples:
            print(f"Processing sample{sample.sample_id} ({sample.sample_name}...")

            genome_builder = self.steps['GenomeBuilderStep']
            genome_builder.execute(sample, self.chr_ploidy_data, self.reference_genome)

            for suffix in [1,2]:

                annotation_updater = self.steps['UpdateAnnotationForGenomeStep']
                annotation_updater.execute(sample, suffix, self.annotation_file_path, self.chr_ploidy_file_path)

        transcriptomes.prep_transcriptomes(self.samples,
                                            self.data_directory_path,
                                            self.log_directory_path,
                                            self.kallisto_file_path,
                                            self.bowtie2_dir_path,
                                            self.output_type,
                                            self.output_molecule_count,
                                            self.dispatcher_mode,
                                            seeds["transcriptomes"])

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
        for job in ["beagle", "transcriptomes"]:
            seeds[job] = numpy.random.randint(MAX_SEED)
        # Seeds for jobs that are run per sample
        for job in ["variant_finder"]:
            for sample in self.samples:
                seeds[f"{job}.{sample.sample_id}"] = numpy.random.randint(MAX_SEED)
        return seeds

    @staticmethod
    def main(configuration, dispatcher_mode, resources, output_directory_path, input_samples):
        pipeline = ExpressionPipeline(configuration, dispatcher_mode, resources, output_directory_path, input_samples)
        if not pipeline.validate():
            raise CampareeValidationException("Expression Pipeline Validation Failed.  "
                                              "Consult the standard error file for details.")
        pipeline.execute()


class CampareeValidationException(CampareeException):
    pass
