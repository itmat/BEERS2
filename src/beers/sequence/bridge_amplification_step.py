import numpy as np
from beers.cluster_packet import ClusterPacket

class BridgeAmplificationStep:
    """
    This step serves to simulate the amplification of the molecule in each cluster on the flowcell.
    Introducing substitutions is allowed via a substitution_rate parameter.

    Config Example::

        # Number of cycles to perform. Generates 2^N molecules in the cluster
        # via a PCR-like bridge amplification process.
        cycles: 10
        # Substitution rate per base in the cluster formation
        # No indels are generated here, for computation simplicity.
        substitution_rate: 0.01
    """

    CYCLE_AT_WHICH_TO_IGNORE_SUBS = 4

    name = "Bridge Amplification Step"

    def __init__(self, step_log_file_path: str, parameters: dict, global_config: dict):
        """
        Initializes the step with a file path to the step log and a dictionary of parameters.  Missing parameters that
        control non-idealized behavior are replaced with default values that produce idealized step behavior.

        Parameters
        ----------
        step_log_file_path:
            location of step logfile
        parameters:
            json-like object of parameters.
        global_config:
            json-like object of configuration for the overall BEERS2 run.
        """
        self.log_filepath = step_log_file_path
        self.cycles: int = parameters["cycles"]
        self.substitution_rate: float = parameters.get("substitution_rate", 0)
        self.global_config = global_config
        print(f"{BridgeAmplificationStep.name} instantiated")


    def execute(self, cluster_packet: ClusterPacket, rng: np.random.Generator) -> ClusterPacket:
        """
        Amplifies the molecule in the cluster for each cluster in the cluster packet by the given number of cycles
        and keeps track of which base positions have incurred substitutions

        Parameters
        ----------
        cluster_packet:
            The input cluster packet
        rng:
            Random number generator instance

        Returns
        ------
            The updated cluster packet
        """
        with open(self.log_filepath, "wt") as log:
            for cluster in cluster_packet.clusters:
                # Generate the initial base_count values
                cluster.initialize_base_counts()

            for cycle in range(1,self.cycles + 1):
                for cluster in cluster_packet.clusters:
                    assert cluster.base_counts is not None
                    # Start with a perfect copy
                    copies = cluster.base_counts.copy()

                    if cycle < self.CYCLE_AT_WHICH_TO_IGNORE_SUBS:
                        # For low cycles, we produce subsitutions
                        # in later ones, they don't matter enough to include

                        # Number of substitutions in each position and of each base type
                        substitutions = rng.binomial(copies, p = self.substitution_rate)

                        # First remove substituted bases
                        copies -= substitutions
                        # then replace them with random draws from A, C, G, T
                        copies += multinomial(substitutions.sum(axis=0), [0.25, 0.25, 0.25, 0.25], rng)

                    # Add to existing molecules
                    cluster.base_counts += copies
                    cluster.molecule_count *= 2
        print(f"Molecules per cluster after amplification is {cluster_packet.clusters[0].molecule_count}.")
        print(f"Total molecules over all clusters is {len(cluster_packet.clusters)*cluster_packet.clusters[0].molecule_count}")
        return cluster_packet


    @staticmethod
    def validate(parameters, global_config):
        errors = []
        if 'cycles' not in parameters:
            errors.append("Must specify 'cycles' for number of bridge amplification cycles.")
        elif not isinstance(parameters['cycles'], int) or parameters['cycles'] <= 0:
            errors.append("'cycles' must be a positive integer")

        substitution_rate = parameters.get("substitution_rate", 0)
        if substitution_rate < 0 or substitution_rate > 1:
            errors.append(f"The substitution rate must be between 0 and 1")

        return errors


def multinomial(n: np.ndarray, p: list, rng: np.random.Generator):
    '''
    Partially vectorized version of np.random.multinomial, see
    https://stackoverflow.com/questions/55818845/fast-vectorized-multinomial-in-python

    Parameters
    ----------

    n:
        1-d array of non-negative integers denoting the 'number of experiments'
    p:
        1-d array of probabilities of each possibility being drawn

    Returns
    ------
    ndarray
        A 2-d array where the first dimension corresponds to the possible draws (length of p)
        and the second to the length of n.
    '''
    n = np.array(n)
    p_array = np.array(p)
    count = n.copy()
    out = np.empty((len(p_array), len(n)), dtype=int)
    ps = p_array.cumsum(axis=-1)
    # Conditional probabilities
    with np.errstate(divide='ignore', invalid='ignore'):
        condp = p / ps
    condp[np.isnan(condp)] = 0.0
    for i in range(p_array.shape[-1]-1, 0, -1):
        binsample = rng.binomial(count, condp[i])
        out[i] = binsample
        count -= binsample
    out[0] = count
    return out
