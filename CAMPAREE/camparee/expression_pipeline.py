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
    def __init__(self, configuration, dispatcher_mode, output_directory_path, input_samples):
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
            kallisto_file_path, bowtie2_dir_path, resources_index_files_directory_path =\
            self.validate_and_gather_resources(configuration['resources'])
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
                os.path.join(resources['directory_path'], "index_files", resources['species_model'])
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
        third_party_software_directory_path = os.path.join(resources['directory_path'], "third_party_software")
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

        bam_files = {}
        for sample in self.samples:

            step_name = 'GenomeAlignmentStep'
            genome_alignment = self.steps[step_name]
            bam_file = genome_alignment.get_genome_bam_path(sample)
            bam_files[sample.sample_id] = bam_file

            if self.dispatcher_mode == "serial":
                genome_alignment.execute(sample, self.resources_index_files_directory_path,
                                         self.star_file_path)
            else:
                stdout_log = os.path.join(genome_alignment.log_directory_path, f"sample{sample.sample_id}", f"{step_name}.bsub.%J.out")
                stderr_log = os.path.join(genome_alignment.log_directory_path, f"sample{sample.sample_id}", f"{step_name}.bsub.%J.err")
                command = genome_alignment.get_commandline_call(sample, self.resources_index_files_directory_path,
                                                                self.star_file_path)
                scheduler_args = {'job_name' : f"{step_name}.sample{sample.sample_id}_{sample.sample_name}",
                                  'stdout_logfile' : stdout_log,
                                  'stderr_logfile' : stderr_log,
                                  'memory_in_mb' : 40000,
                                  'num_processors' : 4}
                validation_attributes = genome_alignment.get_validation_attributes(sample)
                expression_pipeline_monitor.submit_new_job(job_id=f"{step_name}.{sample.sample_id}",
                                                           job_command=command, sample=sample,
                                                           step_name=step_name,
                                                           scheduler_arguments=scheduler_args,
                                                           validation_attributes=validation_attributes,
                                                           output_directory_path=self.output_directory_path,
                                                           system_id=None,
                                                           dependency_list=None)

        for sample in self.samples:

            step_name = 'GenomeBamIndexStep'
            genome_index = self.steps[step_name]
            bam_filename = bam_files[sample.sample_id]

            if self.dispatcher_mode == "serial":
                genome_index.execute(sample, bam_filename)
            else:
                stdout_log = os.path.join(genome_index.log_directory_path, f"sample{sample.sample_id}", f"{step_name}.bsub.%J.out")
                stderr_log = os.path.join(genome_index.log_directory_path, f"sample{sample.sample_id}", f"{step_name}.bsub.%J.err")
                command = genome_index.get_commandline_call(sample, bam_filename)

                scheduler_args = {'job_name' : f"{step_name}.sample{sample.sample_id}_{sample.sample_name}",
                                  'stdout_logfile' : stdout_log,
                                  'stderr_logfile' : stderr_log,
                                  'memory_in_mb' : 6000,
                                  'num_processors' : 1}
                validation_attributes = genome_index.get_validation_attributes(bam_filename)
                expression_pipeline_monitor.submit_new_job(job_id=f"{step_name}.{sample.sample_id}",
                                                           job_command=command, sample=sample,
                                                           step_name=step_name,
                                                           scheduler_arguments=scheduler_args,
                                                           validation_attributes=validation_attributes,
                                                           output_directory_path=self.output_directory_path,
                                                           system_id=None,
                                                           dependency_list=[f"GenomeAlignmentStep.{sample.sample_id}"])


        for sample_id, bam_file in bam_files.items():
            # Use chr_ploidy as the gold std for alignment, variants, VCF, genome_maker
            step_name = 'VariantsFinderStep'
            variants_finder = self.steps[step_name]
            sample = expression_pipeline_monitor.get_sample(sample_id)
            seed = seeds[f"variant_finder.{sample_id}"]

            if self.dispatcher_mode == "serial":
                print(f"Processing variants in sample {sample_id}...")
                variants_finder.execute(sample, bam_file, self.chr_ploidy_data, self.reference_genome, seed)
            else:
                stdout_log = os.path.join(variants_finder.log_directory_path, f"sample{sample_id}", f"{step_name}.bsub.%J.out")
                stderr_log = os.path.join(variants_finder.log_directory_path, f"sample{sample_id}", f"{step_name}.bsub.%J.err")
                command = variants_finder.get_commandline_call(sample, bam_file, self.chr_ploidy_file_path,
                                                               self.reference_genome_file_path, seed)
                scheduler_args = {'job_name' : f"{step_name}.sample{sample_id}_{sample.sample_name}",
                                  'stdout_logfile' : stdout_log,
                                  'stderr_logfile' : stderr_log,
                                  'memory_in_mb' : 6000,
                                  'num_processors' : 1}
                validation_attributes = variants_finder.get_validation_attributes(sample)
                expression_pipeline_monitor.submit_new_job(job_id=f"{step_name}.{sample_id}",
                                                           job_command=command, sample=sample,
                                                           step_name=step_name,
                                                           scheduler_arguments=scheduler_args,
                                                           validation_attributes=validation_attributes,
                                                           output_directory_path=self.output_directory_path,
                                                           system_id=None,
                                                           dependency_list=[f"GenomeBamIndexStep.{sample_id}"])

        for sample_id, bam_file in bam_files.items():
            step_name = 'IntronQuantificationStep'
            intron_quant = self.steps[step_name]
            output_directory = os.path.join(intron_quant.data_directory_path, f"sample{sample.sample_id}")
            sample = expression_pipeline_monitor.get_sample(sample_id)

            if self.dispatcher_mode == "serial":
                print(f"Computing Intron Quantifications in sample {sample_id}")
                intron_quant.execute(bam_file, output_directory, self.annotation_file_path)
            else:
                stdout_log = os.path.join(intron_quant.log_directory_path, f"sample{sample.sample_id}", f"{step_name}.bsub.%J.out")
                stderr_log = os.path.join(intron_quant.log_directory_path, f"sample{sample.sample_id}", f"{step_name}.bsub.%J.err")
                command = intron_quant.get_commandline_call(bam_file, output_directory, self.annotation_file_path)

                scheduler_args = {'job_name' : f"{step_name}.sample{sample_id}_{sample.sample_name}",
                                  'stdout_logfile' : stdout_log,
                                  'stderr_logfile' : stderr_log,
                                  'memory_in_mb' : 6000,
                                  'num_processors' : 1}
                validation_attributes = intron_quant.get_validation_attributes(output_directory)
                expression_pipeline_monitor.submit_new_job(job_id=f"{step_name}.{sample_id}",
                                                           job_command=command, sample=sample,
                                                           step_name=step_name,
                                                           scheduler_arguments=scheduler_args,
                                                           validation_attributes=validation_attributes,
                                                           output_directory_path=output_directory,
                                                           system_id=None,
                                                           dependency_list=[f"GenomeBamIndexStep.{sample_id}"])
                #TODO: do we need to depend upon the index being done? or just the alignment?
                #      I'm hypothesizing that some failures are being caused by indexing and quantification happening
                #      on the same BAM file at the same time, though I don't know why this would be a problem.

        #Wait here until all of the preceding steps have finished.
        print(f'Waiting until all samples finish processing before running Beagle.')
        expression_pipeline_monitor.monitor_until_all_jobs_completed(queue_update_interval=10)

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
    def main(configuration, dispatcher_mode, output_directory_path, input_samples):
        pipeline = ExpressionPipeline(configuration, dispatcher_mode, output_directory_path, input_samples)
        if not pipeline.validate():
            raise CampareeValidationException("Expression Pipeline Validation Failed.  "
                                              "Consult the standard error file for details.")
        pipeline.execute()


class CampareeValidationException(CampareeException):
    pass
