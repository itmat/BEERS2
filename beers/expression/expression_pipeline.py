import sys
import os
import importlib
from beers.constants import CONSTANTS
from beers.expression.variants_finder import VariantsFinderStep
from beers.expression.beagle import BeagleStep
from beers.utilities.expression_utils import ExpressionUtils
from beers.sample import Sample
import pandas as pd

class ExpressionPipeline:
    def __init__(self, configuration, resources, output_directory_path, input_samples):
        self.samples = input_samples
        self.resources = resources
        self.output_directory_path = output_directory_path
        log_directory_path = os.path.join(output_directory_path, CONSTANTS.LOG_DIRECTORY_NAME)
        data_directory_path = os.path.join(output_directory_path, CONSTANTS.DATA_DIRECTORY_NAME)
        self.create_intermediate_data_subdirectories(data_directory_path, log_directory_path)
        self.log_file_path = os.path.join(log_directory_path, "expression_pipeline.log")
        self.steps = {}
        for step in configuration['steps']:
            module_name, step_name = step["step_name"].rsplit(".")
            parameters = step["parameters"]
            module = importlib.import_module(f'.{module_name}', package="beers.expression")
            step_class = getattr(module, step_name)
            self.steps[step_name] = step_class(log_directory_path, data_directory_path, parameters)
        self.reference_genome = dict()
        self.reference_genome_file_path = \
            os.path.join(resources['resources_folder'], "index_files", f"{resources['species_model']}",
                         f"{resources['reference_genome_filename']}")
        self.chr_ploidy_file_path = \
            os.path.join(resources['resources_folder'], "index_files", f"{resources['species_model']}",
                         f"{resources['chr_ploidy_filename']}")
        self.reference_genome = ExpressionUtils.create_reference_genome(self.reference_genome_file_path)
        self.create_chr_ploidy_data()

    def create_intermediate_data_subdirectories(self, data_directory_path, log_directory_path):
        for sample in self.samples:
            os.makedirs(os.path.join(data_directory_path, f'sample{sample.sample_id}'), mode=0o0755, exist_ok=True)
            os.makedirs(os.path.join(log_directory_path, f'sample{sample.sample_id}'), mode=0o0755, exist_ok=True)

    def create_chr_ploidy_data(self):
        df = pd.read_csv(self.chr_ploidy_file_path, sep='\t', index_col=False)
        self.chr_ploidy_data = df.to_dict(orient='records')


    def validate(self):
        reference_genome_chromosomes = self.reference_genome.keys()
        ploidy_chromosomes = set([item['chr'] for item in self.chr_ploidy_data])
        if not ploidy_chromosomes.issubset(reference_genome_chromosomes):
            raise BeersExpressionValidationException(
                f"The chromosome ploidy has chromosomes not found in the reference genome file")

        if not all([step.validate() for step in self.steps.values()]):
            raise BeersExpressionValidationException(f"Validation error in an expression step. See stderr for details.")

    def execute(self):
        print("Execution of the Expression Pipeline Started...")

        for sample in self.samples:

            genome_alignment = self.steps['GenomeAlignmentStep']
            genome_alignment.execute(sample, self.reference_genome)

            variants_finder = self.steps['VariantsFinderStep']
            variants_finder.execute(sample, self.reference_genome, ['19'])

        variants_compilation = self.steps['VariantsCompilationStep']
        variants_compilation.execute(self.samples, self.reference_genome)

        beagle = self.steps['BeagleStep']
        outcome = beagle.execute(self.resources)
        if outcome != 0:
            sys.stderr.write("Beagle process failed.\n")
            sys.exit(1)

        for sample in self.samples:

            genome_builder = self.steps['GenomeBuilderStep']
            genome_builder.execute(sample, self.reference_genome, ['19'])

            #annotation_updates = []
            #transcript_distributions = []
            #molecules = []
            #for item in genomes:
            #    indel_data = item[1]
            #    annotation_updater = UpdateAnnotationForGenome(indel_data, self.annotation_filename)
            #    annotation_updates.append(annotation_updater.update_annotation())
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
        #return molecules

    @staticmethod
    def main(configuration, resources, output_directory_path, input_samples):
        pipeline = ExpressionPipeline(configuration, resources, output_directory_path, input_samples)
        pipeline.validate()
        pipeline.execute()

class BeersExpressionValidationException(Exception):
    pass

class ExpressionPipelineException(Exception):
    pass
