import math
import numpy as np
import sys
import resource

class BridgeAmplificationStep:
    """
    This step serves to simulate the amplification of the molecule in each cluster on the flowcell.  In the idealized
    case, no snps are introduced.  The non-ideal case of introducing snps is allows via a snp_rate parameter.
    """

    BASES = ['G','A','T','C']
    MAX_CYCLES = 16
    MAX_SNP_RATE = 1
    CYCLE_AT_WHICH_TO_IGNORE_SNPS = 4

    name = "Bridge Amplification Step"

    def __init__(self, step_log_file_path, parameters, global_config):
        """
        Initializes the step with a file path to the step log and a dictionary of parameters.  Missing parameters that
        control non-idealized behavior are replaced with default values that produce idealized step behavior.
        :param step_log_file_path: location of step logfile
        :param parameters: dictionary of parameters.  Any required parameters not provided are identified by the
        validate method.
        :param global_config: dictionary of general parameters
        """
        self.log_filename = step_log_file_path
        self.cycles = parameters.get("bridge_amplification_cycles")
        self.snp_rate = parameters.get("bridge_amplification_snp_percentage", 0)/100
        self.global_config = global_config
        print(f"{BridgeAmplificationStep.name} instantiated")


    def execute(self, cluster_packet):
        """
        Amplifies the molecule in the cluster for each cluster in the cluster packet by the given number of cycles
        and keeps track of which base positions have incurred snps
        :param cluster_packet:
        :return:
        """
        for cycle in range(1,self.cycles + 1):
            for cluster in cluster_packet.clusters:
                cluster.molecule_count *= 2
                snps = []
                # For non-zero snp rate, early in the bridge amplification, identify positions receiving snps.  Later
                # occurring snps won't affect the outcome of a read as they won't amplify sufficiently.
                if self.snp_rate > 0 and cycle < BridgeAmplificationStep.CYCLE_AT_WHICH_TO_IGNORE_SNPS:
                    snps = self.determine_snps(cluster, cycle)
                for index, original_base in enumerate(cluster.molecule.sequence):
                    # If no snp occurs at this position, just double the counts for each of the four possible bases.
                    if not snps or index not in snps:
                        for base in 'GATC':
                            getattr(cluster.base_counts, base)[index] *= 2
                    else:
                        count = snps.count(index)
                        original_base = cluster.molecule.sequence[index]
                        # Remove original base from alternatives.  Assuming that a snp does not replace like with like.
                        alternative_bases = BridgeAmplificationStep.BASES[:]
                        alternative_bases.remove(original_base)
                        selected_bases = np.random.choice(alternative_bases, size=count)
                        # Adjust base counts based on snp count.
                        getattr(cluster.base_counts, original_base)[index] *= 2
                        getattr(cluster.base_counts, original_base)[index] -= count
                        unselected_bases = [x for x in BridgeAmplificationStep.BASES if x not in selected_bases and x != original_base]
                        for base in unselected_bases:
                            getattr(cluster.base_counts, base)[index] *= 2
                        for base in selected_bases:
                            getattr(cluster.base_counts, base)[index] += 1
        print(f"Molecules per cluster after amplification is {cluster_packet.clusters[0].molecule_count}.")
        print(f"Total molecules over all clusters is {len(cluster_packet.clusters)*cluster_packet.clusters[0].molecule_count}")
        return cluster_packet

    def determine_snps(self, cluster, cycle):
        """
        Provides a sorted list of snp counts for each position along the cluster molecule for the given cycle.
        :param cluster: cluster containing the molecule being amplified
        :param cycle: the amplification round
        :return: A sorted list of the number of snps at each position along the molecule.
        """
        snp_count = math.floor(len(cluster.molecule.sequence) * self.snp_rate * 2 ** cycle)
        snps = []
        if snp_count > 0:
            snps = np.random.choice(range(0,len(cluster.molecule.sequence)), size=snp_count)
        return sorted(snps)

    def validate(self):
        """
        Insures that the parameters provided are valid.  Error messages are sent to stderr.
        :return: True if the step's parameters are all valid and false otherwise.
        """
        print(f"{BridgeAmplificationStep.name} validating parameters")
        if not self.cycles or self.cycles < 0 or self.cycles > BridgeAmplificationStep.MAX_CYCLES:
            print(f"The number of bridge amplification cycles must be a non-negative integer less than or"
                  f" equal to {BridgeAmplificationStep.MAX_CYCLES}", file=sys.stderr)
            return False
        if self.snp_rate < 0 or self.snp_rate > BridgeAmplificationStep.MAX_SNP_RATE:
            print(f"The snp rate must be a non-negative value less than or"
                  f" equal to {BridgeAmplificationStep.MAX_SNP_RATE} ", file=sys.stderr)
            return False
        return True
