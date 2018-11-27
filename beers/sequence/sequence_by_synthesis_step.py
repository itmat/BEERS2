class SequenceBySynthesisStep:

    def __init__(self, step_log_file_path, parameters):
        self.log_filename = step_log_file_path
        self.read_length = parameters['read_length']
        self.paired_ends = parameters['paired_ends']
        self.barcode_data = (parameters['barcode_5_prime_start'],
                             parameters['barcode_5_prime_length'],
                             parameters['barcode_3_prime_start'],
                             parameters['barcode_3_prime_length'])
        print("Sequence by synthesis step instantiated")

    def execute(self, cluster_packet):
        for cluster in cluster_packet.clusters:
            cluster.compute_quality_scores(self.read_length, self.paired_ends, self.barcode_data)
        cluster_packet.clusters = sorted(cluster_packet.clusters, key=lambda cluster: cluster.coordinates)
        return cluster_packet

    def validate(self):
        return True
