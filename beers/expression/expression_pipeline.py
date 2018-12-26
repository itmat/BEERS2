import sys
import os
import importlib
from beers.expression.variants_finder import VariantsFinderStep
from beers.expression.beagle import BeagleStep
from beers.utilities.expression_utils import ExpressionUtils
from beers.sample import Sample

class ExpressionPipeline:
    def __init__(self, configuration, output_directory_path, input_samples):
        self.samples = input_samples
        self.output_directory_path = output_directory_path
        log_directory_path = os.path.join(output_directory_path, "logs")
        data_directory_path = os.path.join(output_directory_path, 'data')
        self.log_file_path = os.path.join(log_directory_path, "expression_pipeline.log")
        self.steps = {}
        for step in configuration['steps']:
            module_name, step_name = step["step_name"].rsplit(".")
            step_log_filename = f"{step_name}.log"
            step_log_file_path = os.path.join(log_directory_path, step_log_filename)
            parameters = step["parameters"]
            module = importlib.import_module(f'.{module_name}', package="beers.expression")
            step_class = getattr(module, step_name)
            self.steps[step_name] = step_class(step_log_file_path, data_directory_path, parameters)
        input_directory_path = configuration["input"]["directory_path"]
        self.reference_genome = dict()
        self.reference_genome_file_path = \
            os.path.join(input_directory_path, configuration["input"]["model_files"]["reference_genome"])

    def validate(self):
        if not all([step.validate() for step in self.steps.values()]):
            raise BeersExpressionValidationException(f"Validation error in an expression step. See stderr for details.")

    def execute(self):
        print("Execution of the Expression Pipeline Started...")

        self.reference_genome = ExpressionUtils.create_reference_genome(self.reference_genome_file_path)
        #for chromosome,sequence in self.reference_genome.items():
        #    print(chromosome, sequence[:25])

        for sample in self.samples:

            variants_finder = self.steps['VariantsFinderStep']
            inferred_gender = variants_finder.execute(sample, self.reference_genome, 'X')
            sample.gender = sample.gender or inferred_gender

            # Variants to VCF conversion goes here.

        for sample in self.samples:

            beagle = self.steps['BeagleStep']
            outcome = beagle.execute()
            if outcome != 0:
                sys.stderr.write("Beagle process failed.\n")
                sys.exit(1)

            #genome_builder = self.steps['GenomeBuilderStep']
            #genome_builder.execute(sample, self.reference_genome)

            #annotation_updates = []
            #transcript_distributions = []
            molecules = []
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
    def main(configuration, output_directory_path, input_samples):
        pipeline = ExpressionPipeline(configuration, output_directory_path, input_samples)
        pipeline.validate()
        pipeline.execute()

class BeersExpressionValidationException(Exception):
    pass

class ExpressionPipelineException(Exception):
    pass
