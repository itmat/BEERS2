from collections import namedtuple
import numpy as np
from io import StringIO
import math
from contextlib import closing
from beers.utilities.general_utils import GeneralUtils


BaseCounts = namedtuple('BaseCounts', ["G", "A", "T", "C"])

class Cluster:

    next_cluster_id = 1  # Static variable for creating increasing cluster id's

    MIN_ASCII = 33
    MAX_QUALITY = 41

    def __init__(self, run_id, cluster_id, molecule):
        self.coordinates = None
        self.cluster_id = cluster_id
        self.run_id = run_id
        self.molecule = molecule
        self.diameter = 0
        self.molecule_count = 1
        self.quality_scores = []
        self.called_sequences = []
        encoded = np.frombuffer(molecule.sequence.encode("ascii"), dtype='uint8')
        counts = {nt: (encoded == ord(nt)).astype('int32') for nt in "ACGT"} 
        self.base_counts = BaseCounts(counts['G'],counts['A'], counts['T'], counts['C'])

    def assign_coordinates(self, coordinates):
        self.coordinates = coordinates

    def generate_fasta_header(self, sample_name):
        self.header = f"{self.encode_sequence_identifier()}\t{sample_name}"

    def encode_sequence_identifier(self):
        # TODO the 1 is a placeholder for flowcell.  What should we do with this?
        return f"@BEERS:{self.run_id}:1:{self.coordinates.lane}:{self.coordinates.tile}:{self.coordinates.x}:{self.coordinates.y}"

    def read(self, read_length, forward_is_5_prime, paired_ends, barcode_data):
        self.forward_is_5_prime = forward_is_5_prime
        if forward_is_5_prime:
            self.read_in_5_prime_direction(read_length, barcode_data)
            if paired_ends:
                self.read_in_3_prime_direction(read_length, barcode_data)
        else:
            self.read_in_3_prime_direction(read_length, barcode_data)
            if paired_ends:
                self.read_in_5_prime_direction(read_length, barcode_data)

    def read_over_range(self, range_start, range_end):
        with closing(StringIO()) as called_bases:
            with closing(StringIO()) as quality_scores:
                for position in range(range_start, range_end):
                    base_counts = list(zip('GATC', self.get_base_counts_by_position(position)))
                    max_base_count = max(base_counts, key=lambda base_count: base_count[1])
                    number_max_values = len([base_count[0] for base_count in base_counts
                                        if base_count[1] == max_base_count[1]])
                    if number_max_values > 1:
                        called_bases.write('N')
                        quality_scores.write(str(chr(Cluster.MIN_ASCII)))
                    else:
                        other_bases_total_count = sum(count for base, count in base_counts if base != max_base_count[0])
                        prob = other_bases_total_count / self.molecule_count
                        quality_value = min(Cluster.MAX_QUALITY, math.floor(-10 * math.log10(prob))) \
                                        if prob != 0 else Cluster.MAX_QUALITY
                        called_bases.write(max_base_count[0])
                        quality_scores.write(str(chr(quality_value + Cluster.MIN_ASCII)))
                return called_bases.getvalue(), quality_scores.getvalue()

    def read_in_5_prime_direction(self, read_length, barcode_data):
        range_start = barcode_data[0]
        range_end = range_start + read_length
        called_bases, quality_scores = self.read_over_range(range_start, range_end)
        self.quality_scores.append(quality_scores)
        self.called_sequences.append(called_bases)

    def read_in_3_prime_direction(self, read_length, barcode_data):
        range_end = len(self.molecule.sequence) - barcode_data[1]
        range_start = range_end - read_length
        called_bases, quality_scores = self.read_over_range(range_start, range_end)
        quality_scores.reverse()
        called_bases = GeneralUtils.create_complement_strand(called_bases)
        self.quality_scores.append(quality_scores.getvalue())
        self.called_sequences.append(called_bases.getvalue())

    def get_base_counts_by_position(self, index):
        return  (self.base_counts.G[index],
                 self.base_counts.A[index],
                 self.base_counts.T[index],
                 self.base_counts.C[index])

    def __str__(self):
        with closing(StringIO()) as output:
            output.write("pos\tG\tA\tT\tC\torig\n")
            for index in range(len(self.molecule.sequence)):
                output.write(f"{index}\t")
                [output.write(f"{base_count}\t") for base_count in self.get_base_counts_by_position(index)]
                output.write(f"{self.molecule.sequence[index]}\n")
            return output.getvalue()
