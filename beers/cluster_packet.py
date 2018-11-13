class ClusterPacket:

    def __init__(self, packet_id, sample, clusters):
        self.packet_id = packet_id
        self.sample = sample
        self.clusters = clusters