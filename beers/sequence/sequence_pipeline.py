class SequencePipeline:

    def __init__(self, cluster_packet, configuration):
        self.cluster_packet = cluster_packet
        self.parameters = configuration["parameters"]

    def validate(self, **kwargs):
        pass

    def execute(self):
        print("Execution of the Sequence Pipeline Started...")

    @staticmethod
    def main(cluster_packet, configuration):
        sequence_pipeline = SequencePipeline(cluster_packet, configuration)
        sequence_pipeline.validate()
        sequence_pipeline.execute()


class BeersSequenceValidationException(Exception):
    pass

