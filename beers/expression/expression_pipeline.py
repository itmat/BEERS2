import sys
import json
import pysam
import os
import re
from io import StringIO
from expression.variants_finder import VariantsFinder

class ExpressionPipeline:
    def __init__(self, configuration):
        self.alignment_file_path = self.get_input_file_path(configuration, "alignment")
        self.chromosomes = []
        self.alignment_file = pysam.AlignmentFile(self.alignment_file_path, "rb")
        for item in self.alignment_file.header['SQ']:
            self.chromosomes.append(item['SN'])

        self.reference_genome_file_path = self.get_input_file_path(configuration, "reference_genome")
        self.reference_genome_input = self.file_input_generator(self.reference_genome_file_path)

        self.output_directory = configuration["output"]

        self.parameters = {}
        for item in configuration["expression_pipeline"]["processes"]:
            self.parameters[item["class_name"]] = item.get("parameters", dict())

    def get_input_file_path(self, configuration, file_type):
        """
        Convenience method to obtain full path to input file of the given file type
        :param configuration: dictionary idenitifying input directory and file name of given file type
        :param file_type: alignment, reference_genome, annotation
        :return: path to the input file of the given file type
        """
        input_directory_path = configuration["input"]["directory_path"]
        filename = configuration["input"]["files"][file_type]
        return os.path.join(input_directory_path, filename)

    def file_input_generator(self, file_path):
        """
        Convenience method to obtain a generator for a given file_path
        :return: generator - .next() returns next line or None at end of file.
        """
        with open(file_path, 'r') as input_file:
            for line in input_file:
                yield line
        return None

    def retrieve_reference_chromosome_sequence(self, generator, chromosome):
        building_sequence = False
        sequence = StringIO()
        for line in generator:
            if line.startswith(">"):
                if building_sequence:
                    reference_chromosome_sequence = sequence.getvalue()
                    sequence.close()
                    return reference_chromosome_sequence
                    identifier = re.sub(r'[ \t].*\n', '', line)[1:]
                if identifier == chromosome:
                    building_sequence = True
                continue
            elif building_sequence:
                sequence.write(line.rstrip('\n').upper())
        return None

    def validate(self):
        pass

    def execute(self):
        print("Execution of the Expression Pipeline Started...")
        for chromosome in self.chromosomes:

            reference_sequence = \
                self.retrieve_reference_chromosome_sequence(self.reference_genome_input, chromosome)
            variants_finder = \
                VariantsFinder(chromosome, self.alignment_file, reference_sequence, self.parameters["VariantsFinder"])
            variants = variants_finder.collect_reads()

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
        pipeline = ExpressionPipeline(configuration)
        pipeline.validate()
        pipeline.execute()

if __name__ == "__main__":
    sys.exit(ExpressionPipeline.main())
