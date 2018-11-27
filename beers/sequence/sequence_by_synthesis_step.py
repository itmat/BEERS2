class SequenceBySynthesisStep:

    def __init__(self, step_log_file_path, parameters):
        self.log_filename = step_log_file_path
        print("Sequence by synthesis step instantiated")

    def execute(self, cluster_packet):
        for cluster in cluster_packet.clusters:
            cluster.compute_quality_score()
        return cluster_packet

    def validate(self):
        return True



