class SequenceBySynthesisStep:

    def __init__(self, parameters):
        self.parameters = parameters

    def execute(self, cluster_packet):
        for cluster in cluster_packet.clusters:
            cluster.compute_quality_score()
        return cluster_packet



