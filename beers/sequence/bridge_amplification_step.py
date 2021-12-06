import sys
import numpy as np

class BridgeAmplificationStep:
    """
    This step serves to simulate the amplification of the molecule in each cluster on the flowcell.  In the idealized
    case, no substitutions are introduced.  The non-ideal case of introducing substitutions is allows via a substitution_rate parameter.
    """

    BASES = ['G','A','T','C']
    MAX_CYCLES = 16
    MAX_SNP_RATE = 1
    CYCLE_AT_WHICH_TO_IGNORE_SUBS = 4

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
        self.cycles = parameters.get("cycles")
        self.substitution_rate = parameters.get("substitution_rate", 0)
        self.global_config = global_config
        print(f"{BridgeAmplificationStep.name} instantiated")


    def execute(self, cluster_packet):
        """
        Amplifies the molecule in the cluster for each cluster in the cluster packet by the given number of cycles
        and keeps track of which base positions have incurred substitutions
        :param cluster_packet:
        :return:
        """
        for cycle in range(1,self.cycles + 1):
            for cluster in cluster_packet.clusters:
                # Start with a perfect copy
                copies = cluster.base_counts.copy()

                if cycle < self.CYCLE_AT_WHICH_TO_IGNORE_SUBS:
                    # For low cycles, we produce subsitutions
                    # in later ones, they don't matter enough to include

                    # Number of substitutions in each position and of each base type
                    substitutions = np.random.binomial(copies, p = self.substitution_rate)

                    # First remove substituted bases
                    copies -= substitutions
                    # then replace them with random draws from A, C, G, T
                    copies += multinomial(substitutions.sum(axis=0), [0.25, 0.25, 0.25, 0.25])

                # Add to existing molecules
                cluster.base_counts += copies
                cluster.molecule_count *= 2
        print(f"Molecules per cluster after amplification is {cluster_packet.clusters[0].molecule_count}.")
        print(f"Total molecules over all clusters is {len(cluster_packet.clusters)*cluster_packet.clusters[0].molecule_count}")
        return cluster_packet


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
        if self.substitution_rate < 0 or self.substitution_rate > BridgeAmplificationStep.MAX_SNP_RATE:
            print(f"The substitution rate must be a non-negative value less than or"
                  f" equal to {BridgeAmplificationStep.MAX_SNP_RATE} ", file=sys.stderr)
            return False
        return True


def multinomial(n, p):
    '''
    Partially vectorized version of np.random.multinomial, see
    https://stackoverflow.com/questions/55818845/fast-vectorized-multinomial-in-python

    :param n: 1-d array of non-negative integers denoting the 'number of experiments'
    :param p: 1-d array of probabilities of each possibility being drawn

    Returns an 2-d array where the first dimension corresponds to the possible draws (length of p)
    and the second to the length of n
    '''
    n = np.array(n)
    p = np.array(p)
    count = n.copy()
    out = np.empty((len(p), len(n)), dtype=int)
    ps = p.cumsum(axis=-1)
    # Conditional probabilities
    with np.errstate(divide='ignore', invalid='ignore'):
        condp = p / ps
    condp[np.isnan(condp)] = 0.0
    for i in range(p.shape[-1]-1, 0, -1):
        binsample = np.random.binomial(count, condp[i])
        out[i] = binsample
        count -= binsample
    out[0] = count
    return out
