class ClusterPacket:

    next_cluster_packet_id = 1  # Static variable for creating increasing cluster packet id's

    def __init__(self, cluster_packet_id, sample, clusters):
        self.cluster_packet_id = cluster_packet_id
        self.sample = sample
        self.clusters = clusters