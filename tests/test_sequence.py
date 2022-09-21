import numpy
from beers_utils.molecule import Molecule
from beers_utils.sample import Sample
from beers.cluster import Cluster
from beers.cluster_packet import ClusterPacket
from beers.sequence.bridge_amplification_step  import BridgeAmplificationStep
from beers.sequence.sequence_by_synthesis_step  import SequenceBySynthesisStep

def make_cluster_packet(count, length, rng):
    def random_sequence(length):
        return ''.join(rng.choice(['A', 'C', 'G', 'T'], size=length))

    clusters = []
    for i in range(count):
        clusters.append(Cluster(
            cluster_id = 5,
            molecule = Molecule(
                molecule_id = "5",
                sequence = random_sequence(length),
                start = 5,
                cigar = f"{length}M",
                strand = '+',
                transcript_id = "ENSMUSG0000000000",
                source_start = rng.integers(1, 5_000_000_000),
                source_cigar =  f"{length}M",
                source_strand = rng.choice(["-", "+"]),
                source_chrom = "X",
            ),
            lane = 1,
            coordinates = (1,3,4),
        ))

    return ClusterPacket(
        clusters = clusters,
        cluster_packet_id = 5,
        sample = Sample("1", "sample1", None, None, None)
    )

def test_BridgeAmplificationStep(tmp_path):
    rng = numpy.random.default_rng(0)

    cluster_packet = make_cluster_packet(count = 10, length= 300, rng = rng)
    original_clusters = cluster_packet.clusters

    step = BridgeAmplificationStep(
        step_log_file_path = tmp_path / "log.txt",
        parameters = {
            "cycles": 10,
            "substitution_rate": 0.01,
        },
        global_config = {
        }
    )

    # Run the step
    output = step.execute(cluster_packet, rng)

    for cluster in output.clusters:
        assert cluster.molecule_count == 2**10
        assert cluster.base_counts.shape[0] == 4
        assert (cluster.base_counts.sum(axis=0) == 2**10).all()

def test_SequenceBySynthesis(tmp_path):
    rng = numpy.random.default_rng(0)

    cluster_packet = make_cluster_packet(count = 10, length= 300, rng = rng)
    for cluster in cluster_packet.clusters:
        cluster.initialize_base_counts()
        cluster.base_counts *= 1024
        cluster.molecule_count = 1024
    original_clusters = cluster_packet.clusters

    step = SequenceBySynthesisStep(
        step_log_file_path = tmp_path / "log.txt",
        parameters = {
            "forward_is_5_prime": True,
            "paired_ends": True,
            "read_length": 100,
            "skip_rate": 0.002,
            "drop_rate": 0.002,
        },
        global_config = {
            "samples": {
                '1': dict( barcodes = dict( i5 = "AGCGCTAG", i7 = "AACCGCGG") ),
                '2': dict( barcodes = dict( i5 = "GATATCGA", i7 = "TTATAACC") ),
            },
            "resources": {
                "pre_i5_adapter": "AATGATACGGCGACCACCGAGATCTACAC",
                "post_i5_adapter": "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
                "pre_i7_adapter": "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
                "post_i7_adapter": "ATCTCGTATGCCGTCTTCTGCTTG",
            }
        }
    )

    # Run the step
    output = step.execute(cluster_packet, rng)

    for cluster in output.clusters:
        assert isinstance(cluster.called_barcode,str)
        assert len(cluster.quality_scores) == 2
        assert len(cluster.read_starts) == 2
        assert len(cluster.read_cigars) == 2
        assert len(cluster.read_strands) == 2
