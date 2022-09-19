from beers_utils.sample import Sample
from beers.cluster import Cluster
import os
import gzip

class ClusterPacket:
    """
    The cluster packet object is the medium of exchange between steps in the sequence pipeline.  Each step execution
    in the sequence pipeline accepts a cluster packet as input and returns a chuster packet, usually modified somehow,
    as output.  The cluster packet is primarily composed of clusters that are derived from molecules.  The cluster
    packet itself is derived from the molecule packet that contained those precursor molecules to begin with.  At the
    end of the pipeline, when FASTQ reports are being created the data contained in the cluster packet is used to
    populate those FASTQ reports.
    """

    next_cluster_packet_id = 0  # Static variable for creating increasing cluster packet id's

    def __init__(self, cluster_packet_id, sample, clusters):
        """
        The cluster packet is mostly a wrapper for the clusters but additionally, it does contain data from the
        sample from which the original molecules were drawn.  Note then that a cluster packet respresents a portion
        of just one sample.

        Parameters
        ---------
        cluster_packet_id: int
            Unique identifier for a cluster packet.  This is normally included in filenames to
            assure uniqueness when saving cluster data to disk.
        sample: beers_utils.sample.Sample
            A sample object representing the sample ancestor of the contained clusters
        clusters: list of Cluster
            The contained clusters
        """
        self.cluster_packet_id = cluster_packet_id
        self.sample = sample
        self.clusters = clusters

    def __str__(self):
        """
        String representation of a cluster packet that may be displayed when a cluster packet object is printed.
        This may not be a complete representation.

        Returns
        -------
        A string representing the cluster packet.
        """
        return f"cluster_packet_id: {self.cluster_packet_id}, sample_name: {self.sample.sample_name}, " \
               f"# of clusters: {len(self.clusters)}"

    def serialize(self, file_path):
        """
        Cluster packets are serialized and saved to the file system and likewise de-serialized from the file system.
        The serialized data is written, in compressed form, to a gzip file.  The first line contains the cluster packet
        id and the serialized sample data, prepended with a '#'.  Following that, each cluster is serialized and
        added to the file.

        Parameters
        ----------
        file_path: str
            location of the file into which the serialized, compressed data is to go.
        """
        with gzip.open(file_path, 'wb') as obj_file:
            obj_file.write(f"#{self.cluster_packet_id}\n#{self.sample.serialize()}\n".encode())
            for cluster in self.clusters:
                obj_file.write(cluster.serialize().encode())
                # Clusters take up a variable number of lines, so we need a separator
                obj_file.write("-\n".encode())

    @staticmethod
    def deserialize(file_path, skip_base_counts=False):
        """
        This method re-rendered that serialized, compressed data found in the gzipped file located via the given
        file path, into a fully restored object of the ClusterPacket class.

        Parameters
        ----------
        file_path: str
            The locatation of the gzipped file containing the serialized object.
        skip_base_counts: bool
            if True, don't load base counts (for memory efficiency)

        Returns
        -------
        ClusterPacket
        """
        cluster_lines = []
        clusters = []
        cluster_packet_id = 0
        sample = None
        with gzip.open(file_path, 'rb') as obj_file:
            for line_number, line in enumerate(obj_file):
                line = line.rstrip(b'\n')
                if line_number == 0:
                    cluster_packet_id = int(line[1:].decode())
                elif line_number == 1:
                    sample = Sample.deserialize(line.decode())
                else:
                    if line.decode() == '-':
                        if cluster_lines:
                            clusters.append(Cluster.deserialize("\n".join(cluster_lines), skip_base_counts))
                        cluster_lines = []
                    else:
                        cluster_lines.append(line.decode())
        return ClusterPacket(cluster_packet_id, sample, clusters)
