from beers.sequence.flowcell_loading_step import FlowcellLoadingStep
from beers.cluster import Cluster
from beers.cluster_packet import ClusterPacket
import os
import pickle

class SequencePipeline:

    def __init__(self, configuration):
        input_directory_path = configuration["input"]["directory_path"]
        self.cluster_packets = []
        for molecule_packet_filename in configuration["input"]["packets"]:
            molecule_packet_file_path = os.path.join(input_directory_path, molecule_packet_filename)
            with open(molecule_packet_file_path, 'rb') as molecule_packet_file:
                molecule_packet = pickle.load(molecule_packet_file)
                self.cluster_packets.append(self.convert_molecule_pkt_to_cluster_pkt(molecule_packet))
        self.parameters = {}
        for item in configuration["steps"]:
            self.parameters[item["class_name"]] = item.get("parameters", dict())

    def validate(self, **kwargs):
        pass

    def execute(self):
        print("Execution of the Sequence Pipeline Started...")
        flowcell_loading_step = FlowcellLoadingStep(self.cluster_packets, self.parameters["FlowcellLoadingStep"])
        flowcell_loading_step.execute()

    @staticmethod
    def convert_molecule_pkt_to_cluster_pkt(molecule_packet):
        clusters = []
        for molecule in molecule_packet.molecules:
            cluster_id = Cluster.next_cluster_id
            clusters.append(Cluster(cluster_id, molecule))
            Cluster.next_cluster_id += 1
        return ClusterPacket(molecule_packet.sample, clusters)

    @staticmethod
    def main(configuration):
        sequence_pipeline = SequencePipeline(configuration)
        sequence_pipeline.validate()
        sequence_pipeline.execute()


class BeersSequenceValidationException(Exception):
    pass

