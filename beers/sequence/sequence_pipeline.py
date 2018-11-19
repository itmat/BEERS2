class SequencePipeline:

    def __init__(self, configuration):
        self.parameters = configuration["parameters"]

    def validate(self, **kwargs):
        pass

    def execute(self):
        print("Execution of the Sequence Pipeline Started...")

    @staticmethod
    def main(configuration):
        sequence_pipeline = SequencePipeline(configuration)
        sequence_pipeline.validate()
        sequence_pipeline.execute()


class BeersSequenceValidationException(Exception):
    pass

