import sys
import json
import pysam
import os
import re
from io import StringIO
from beers.expression.variants_finder import VariantsFinder

class ExpressionPipeline:
    def __init__(self, configuration):
        self.reference_genome = dict()
        input_directory_path = configuration["input"]["directory_path"]
        self.reference_genome_file_path = os.path.join(input_directory_path, configuration["input"]["files"]["reference_genome"])
        self.alignment_file_paths = [os.path.join(input_directory_path, entry)
                                     for entry in configuration["input"]["files"]["alignment"]]
        try:
            self.output_directory_path = configuration["output"]["directory_path"]
        except FileExistsError:
            pass
        self.parameters = {}
        for item in configuration["processes"]:
            self.parameters[item["class_name"]] = item.get("parameters", dict())

    def create_reference_genome(self):
        fasta_chromosome_pattern = re.compile(">([^\s]*)")
        chromosome, sequence = '', None
        building_sequence = False
        with open(self.reference_genome_file_path, 'r') as reference_genome_file:
            for line in reference_genome_file:
                if line.startswith(">"):
                    if building_sequence:
                        self.reference_genome[chromosome] = sequence.getvalue()
                        sequence.close()
                    chromosome_match = re.match(fasta_chromosome_pattern, line)
                    chromosome = chromosome_match.group(1)
                    building_sequence = True
                    sequence = StringIO()
                    continue
                elif building_sequence:
                    sequence.write(line.rstrip('\n').upper())

    def validate(self):
        pass

    def execute(self):
        print("Execution of the Expression Pipeline Started...")

        self.create_reference_genome()
        #for chromosome,sequence in self.reference_genome.items():
        #    print(chromosome, sequence[:25])

        for alignment_file_path in self.alignment_file_paths:

            variants_finder = \
                VariantsFinder(None, alignment_file_path, self.reference_genome, self.parameters["VariantsFinder"], self.output_directory_path)
            self.gender = variants_finder.find_variants()

            #genome_maker = GenomeMaker(chromosome, variants, reference_sequence, self.parameters["GenomeMaker"])
            #genomes = genome_maker.make_genomes()
            #variants = None  # Should we dereference objects no longer needed to save memory?
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
    def main():
        with open("../../config/config.json", "r+") as configuration_file:
            configuration = json.load(configuration_file)
        pipeline = ExpressionPipeline(configuration["expression_pipeline"])
        pipeline.validate()
        pipeline.execute()

if __name__ == "__main__":
    sys.exit(ExpressionPipeline.main())