import math
import numpy as np
import sys

class BridgeAmplificationStep:

    BASES = ['G','A','T','C']
    MAX_CYCLES = 16
    MAX_SNP_RATE = 1

    def __init__(self, step_log_file_path, parameters):
        self.log_filename = step_log_file_path
        self.cycles = parameters["bridge_amplification_cycles"]
        self.snp_rate = parameters["bridge_amplification_snp_percentage"]/100
        print("Bridge amplification step instantiated")


    def execute(self, cluster_packet):
        print(f"Cluster Pkt size: {len(cluster_packet.clusters)}")
        for cycle in range(1,self.cycles + 1):
            for cluster in cluster_packet.clusters:
                cluster.molecule_count *= 2
                snps = []
                if self.snp_rate > 0 and cycle < 4:
                    snps = self.determine_snps(cluster, cycle)
                for index, original_base in enumerate(cluster.molecule.sequence):
                    if not snps or index not in snps:
                        for base in 'GATC':
                            getattr(cluster.base_counts, base)[index] *= 2
                    else:
                        count = snps.count(index)
                        original_base = cluster.molecule.sequence[index]
                        alternative_bases = BridgeAmplificationStep.BASES[:]
                        alternative_bases.remove(original_base)
                        selected_bases = np.random.choice(alternative_bases, size=count)
                        getattr(cluster.base_counts, original_base)[index] *= 2
                        getattr(cluster.base_counts, original_base)[index] -= count
                        unselected_bases = [x for x in BridgeAmplificationStep.BASES if x not in selected_bases and x != original_base]
                        for base in unselected_bases:
                            getattr(cluster.base_counts, base)[index] *= 2
                        for base in selected_bases:
                            getattr(cluster.base_counts, base)[index] += 1
        return cluster_packet

    def determine_snps(self, cluster, cycle):
        snp_count = math.floor(len(cluster.molecule.sequence) * self.snp_rate * 2 ** cycle)
        snps = []
        if snp_count > 0:
            snps = np.random.choice(range(0,len(cluster.molecule.sequence)), size=snp_count)
        return sorted(snps)

    def validate(self):
        print(f"Bridge amplification step validating parameters")
        if self.cycles < 0 or self.cycles > BridgeAmplificationStep.MAX_CYCLES:
            print(f"The number of bridge amplification cycles must be a non-negative integer less than or"
                  f" equal to {BridgeAmplificationStep.MAX_CYCLES}", file=sys.stderr)
            return False
        if self.snp_rate < 0 or self.snp_rate > BridgeAmplificationStep.MAX_SNP_RATE:
            print(f"The snp rate must be a non-negative value less than or"
                  f" equal to {BridgeAmplificationStep.MAX_SNP_RATE} ", file=sys.stderr)
            return False
        return True






