import sys

class SequenceBySynthesisStep:

    name = "Sequence By Synthesis Step"

    def __init__(self, step_log_file_path, parameters, global_config):
        self.log_filename = step_log_file_path
        self.read_length = parameters['read_length']
        self.forward_is_5_prime = parameters.get('forward_is_5_prime', True)
        self.paired_ends = parameters.get('paired_ends', False)

        # Determine where the barcodes lie
        self.i5_start = len(global_config['resources']['pre_i5_adapter']) + 1 # one past the first part of the adapter
        self.i7_start = len(global_config['resources']['post_i7_adapter']) + 1 # one past the first part of the adapter, reading from 3' end
        i5_adapter_lengths = [len(sample['barcodes']['i5']) for sample in  global_config['samples'].values()]
        assert len(set(i5_adapter_lengths)) == 1
        self.i5_length = i5_adapter_lengths[0]
        i7_adapter_lengths = [len(sample['barcodes']['i7']) for sample in  global_config['samples'].values()]
        assert len(set(i7_adapter_lengths)) == 1
        self.i7_length = i7_adapter_lengths[0]

        # Determine where the sequencing starts
        self.forward_read_start = self.i5_start + self.i5_length + len(global_config['resources']['post_i5_adapter'])
        # NOTE: extra 'A' is ligated to the 3' end of the reads prior to adapter ligation (done in the ligation step in BEERS)
        #       hence we have an additional +1 on this read start
        self.reverse_read_start = self.i7_start + self.i7_length + len(global_config['resources']['pre_i7_adapter']) + 1

        print(f"{SequenceBySynthesisStep.name} instantiated")

    def execute(self, cluster_packet):
        for cluster in cluster_packet.clusters:
            cluster.set_forward_direction(self.forward_is_5_prime)
            cluster.read(
                    self.read_length, self.paired_ends,
                    self.i5_start, self.i5_length, self.i7_start, self.i7_length,
                    self.forward_read_start, self.reverse_read_start
            )
        cluster_packet.clusters = sorted(cluster_packet.clusters, key=lambda cluster: cluster.coordinates)
        return cluster_packet

    def validate(self):
        print(f"{SequenceBySynthesisStep.name} validating parameters")
        if not self.read_length or self.read_length < 0 or not isinstance(self.read_length, int):
            print(f"The read length, {self.read_length}, is required and must be a non-negative integer,", file=sys.stderr)
            return False
        return True
