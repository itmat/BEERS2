from beers.sequence.bridge_amplification_step import BridgeAmplificationStep
from beers.sequence.sequence_by_synthesis_step import SequenceBySynthesisStep

class SequencePipeline:

    def __init__(self, configuration):
        self.parameters = configuration["parameters"]

    def validate(self, **kwargs):
        pass

    def execute(self, cluster_packet):
        print("Execution of the Sequence Pipeline Started...")
        bridge_amplification_step = BridgeAmplificationStep(self.parameters)
        bridge_amplification_step.execute(cluster_packet)
        sequence_by_synthesis_step = SequenceBySynthesisStep(self.parameters)
        sequence_by_synthesis_step.execute(cluster_packet)

    @staticmethod
    def main(cluster_packet, configuration):
        sequence_pipeline = SequencePipeline(configuration)
        sequence_pipeline.validate()
        sequence_pipeline.execute(cluster_packet)


class BeersSequenceValidationException(Exception):
    pass

