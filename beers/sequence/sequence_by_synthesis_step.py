import sys

class SequenceBySynthesisStep:

    name = "Sequence By Synthesis Step"

    def __init__(self, step_log_file_path, parameters):
        self.log_filename = step_log_file_path
        self.read_length = parameters['read_length']
        self.foward_is_5_prime = parameters.get('forward_is_5_prime', True)
        self.paired_ends = parameters.get('paired_ends', False)
        self.barcode_data = (parameters['barcode_5_prime_start'],
                             parameters['barcode_5_prime_length'],
                             parameters['barcode_3_prime_start'],
                             parameters['barcode_3_prime_length'])
        print(f"{SequenceBySynthesisStep.name} instantiated")

    def execute(self, cluster_packet):
        for cluster in cluster_packet.clusters:
            cluster.read(self.read_length, self.foward_is_5_prime, self.paired_ends, self.barcode_data)
        cluster_packet.clusters = sorted(cluster_packet.clusters, key=lambda cluster: cluster.coordinates)
        return cluster_packet

    def validate(self):
        print(f"{SequenceBySynthesisStep.name} validating parameters")
        if not self.read_length or self.read_length < 0 or not isinstance(self.read_length, int):
            print(f"The read length, {self.read_length}, is required and must be a non-negative integer,", file=sys.stderr)
            return False
        return True
