import math
import collections
import sys

import numpy
import scipy.stats

from beers_utils.molecule import Molecule
from beers_utils.molecule_packet import MoleculePacket


METHODS = {"uniform", "beta"}

# TODO: should this be a parameter?
# For fragment2
# just makes lambda into a reasonable value for a molecule of the "typical" size
# i.e. lambda gives the rate of breakage of a molecule which has this length
TYPICAL_MOLECULE_SIZE = 1000


class FragmentStep:
    """
    Fragmentation step breaks each molecule into pieces

    Two fragmentation methods are available:
     1. 'uniform' fragments at each base equally and is the default.
        Note that fragment length distributions will not be uniform.
        Instead, the fragment locations are chosen uniformly.
     2. 'beta' biases fragment location by position within the fragment
        and can bias not fragmenting already short fragments
        NOTE: using 'beta' can generate unusual coverage plots since there are
        significant edge effects around the transcript. However, it can be used
        to give a more realistic fragment distribution.

    Configuration
    -------------

    The 'lambda' parameter is the rate parameter for an exponential distribution
    which determines the time it takes until a molecule to fragments.
    Molecules which would take longer then 'runtime' to fragment, do not fragment.
    Molecules may also fragment multiple times if lambda is high enough.

    Since fragmentation generates many very small fragments that will
    be quickly lost in the following steps, we have the option to
    drop those fragments immediately. Setting ``min_frag_size`` can significantly
    decrease runtime and memory requirements.

    If method == 'beta', then the fragmentation sites can be biased within
    the molecule, according to a beta distribution.
    If using 'beta', you must also specify the following:
   
    The parameters for the beta distribution, beta_A and beta_B.
    Setting these so that A = B > 0 to bias towards fragmentation in
    the middle of the fragment, with larger values biasing further towards the middle.
    If A > B, then would bias towards the 5' end instead.
   
    And the beta_N factor allows a non-linear fragmentation rate depending
    upon the length of the molecule. Values >1 indicate that longer molecules
    are more likely to fragment than smaller ones, biasing towards larger
    fragment sizes.

    Config example::

        parameters:
            # Fragmentation methods are available.
            # 'uniform' fragments at each base equally and is the default
            method: uniform
            # The 'lambda' parameter is the rate parameter for an exponential distribution
            # which determines the time it takes until a molecule to fragments.
            # Molecules which would take longer then 'runtime' to fragment, do not fragment.
            # Molecules may also fragment multiple times if lambda is high enough.
            lambda: 0.005
            runtime: 1
            # Since fragmentation generates many very small fragments that will
            # be quickly lost in the following steps, we have the option to
            # drop those fragments immediately. Setting this value can significantly
            # decrease runtime and memory requirements.
            min_frag_size: 20

            # If method == 'beta', then the fragmentation sites can be biased within
            # the molecule, according to a beta distribution.
            # NOTE: using 'beta' can generate unusual coverage plots since there are
            # significant edge effects around the transcript. However, it can be used
            # to give a more realistic fragment distribution.
            #
            # If using 'beta', you must also specify the following:
            #
            # The parameters for the beta distribution. Set these so that
            # A = B > 0 to bias towards fragmentation in the middle of the fragment,
            # with larger values biasing further towards the middle.
            # If A > B, then would bias towards the 5' end instead.
            # beta_A: 3.0
            # beta_B: 3.0
            #
            # And the N factor allows a non-linear fragmentation rate depending
            # upon the length of the molecule. Values >1 indicate that longer molecules
            # are more likely to fragment than smaller ones, biasing towards larger
            # fragment sizes.
            # beta_N: 2.0

    """

    name = "Fragment Step"

    def __init__(self, parameters, global_config):
        self.global_config = global_config
        self.method = parameters["method"]
        self.lambda_ = parameters["lambda"]
        self.runtime = parameters["runtime"]
        # Reject fragments that are below this size. This limits computation and memory significantly
        # and tiny fragments fail priming and/or size selection later regardless. Most generated fragments
        # are very small, so even small minimum size requirements reduce the number of fragments dramatically
        self.min_frag_size = parameters["min_frag_size"]

        # Parameters used ONLY for beta_fragmentation method
        if self.method == "beta":
            self.beta_A = parameters["beta_A"]
            self.beta_B = parameters["beta_B"]
            self.beta_N = parameters["beta_N"]

    def execute(self, molecule_packet, rng, log):
        print("Fragment Step acting on sample")
        sample = molecule_packet.molecules
        if self.method == "uniform":
            fragment_locations =  _compute_fragment_locations_uniform(sample, self.lambda_, self.runtime, rng)
        elif self.method == "beta":
            fragment_locations = _compute_fragment_locations_beta(sample, self.lambda_, self.beta_N, self.beta_A, self.beta_B, self.runtime, rng)
        else:
            raise NotImplementedError(f"Unknown fragmentation method {self.method}")

        result = [sample[k].make_fragment(start+1,end) for (start, end, k) in fragment_locations
                    if end - start >= self.min_frag_size]

        for molecule in result:
            log.write(molecule)

        molecule_packet.molecules = result
        return molecule_packet

    @staticmethod
    def validate(parameters, global_config):
        errors = []
        if 'method' not in parameters:
            errors.append(f"Must specify 'method' as one of {METHODS}")
        else:
            method = parameters['method']
            if method not in METHODS:
                errors.append(f"Fragmentation method must be one of the following: {METHODS}, instead received {method}")

            # Parameters used ONLY for beta_fragmentation method
            if method == "beta":
                if any((beta_param not in parameters) for beta_param in ['beta_A', 'beta_B', 'beta_N']):
                    errors.append("Must specify all of 'beta_A', 'beta_B', and 'beta_N' when using method == 'beta'")
                else:
                    beta_A = parameters["beta_A"]
                    beta_B = parameters["beta_B"]
                    beta_N = parameters["beta_N"]
                    if any((not isinstance(beta_param, (float, int))
                                or (beta_param < 0))
                                for beta_param in [beta_A, beta_B, beta_N]):
                        errors.append("All of 'beta_A', 'beta_B', and 'beta_N' must be positive numbers")
        if "lambda" not in parameters:
            errors.append("Must specify 'lambda' value")
        elif (not isinstance(parameters['lambda'], (int, float)) or parameters['lambda'] <= 0):
            errors.append("Rate of fragmentation must be positive number")

        if "runtime" not in parameters:
            errors.append("Must specify 'runtime' value")
        elif (not isinstance(parameters['runtime'], (int, float)) or parameters['runtime'] <= 0):
            errors.append("Fragmentation 'runtime' must be positive number")

        if "min_frag_size" not in parameters:
            errors.append("Must specify 'min_frag_size' value")
        elif (not isinstance(parameters['min_frag_size'], int) or parameters['min_frag_size'] <= 0):
            errors.append("Fragmentation 'min_frag_size' must be positive integer")

        return errors

### FRAGMENTATION METHODS
def _sample_without_replacement(n, k, rng):
    """uniformly sample k numbers from 0,...,n-1 without replacement

    Intended to be faster than numpy.random.choice(n, size=k, replace=False)
    for the case when n is large and k is small.
    (About 1000x times faster for n=1,000,000 and k=4.)
    Sample is returned in sorted order """

    sample = set(rng.choice(n, size=k, replace=True))
    while len(sample) < k:
        sample = sample.union( rng.choice(n, size=(k-len(sample)), replace=True) )

    # Sort it since list(set) will give a sort-of arbitrary but not random order
    return numpy.array(sorted(sample), dtype=int)

def _compute_fragment_locations_uniform(molecules, lambda_, runtime, rng):
    """
    uniform fragmentation with a rate lambda_ parameter
    All bonds between adjacent bases are equally likely to break (they're iid)

    Use the fact that the methods of fragment is equivalent (when lambda constant)
    to first sampling the number of break points from a binomial distribution and then
    picking those points uniformly on the molecule

    returns a list of tuples (start, end, k) of fragments
    each fragment being k'th molecule from positions start to end (zero-based, non-inclusive of end)
    """
    assert runtime > 0
    assert lambda_ > 0

    # Breaking occurs with probability `lambda_` per unit time
    # so the time to break a bond is exponentially distributed.
    # We convert that to the chance that this time-to-break
    # is at most equal to the runtime.
    # With runtime = 1, this is very nearly just `lambda_`
    probability_of_base_breaking = scipy.stats.expon(scale=1/lambda_).cdf(runtime)

    output = collections.deque()
    for k, molecule in enumerate(molecules):
        num_bonds = len(molecule) - 1

        num_breakpoints = rng.binomial(n=num_bonds, p=probability_of_base_breaking)

        breakpoints = _sample_without_replacement(num_bonds, num_breakpoints, rng)
        bps = numpy.concatenate([[0], breakpoints + 1, [len(molecule)]])
        output.extend( (bps[i], bps[i+1], k) for i in range(len(bps)-1) )

    return list(output)

### Parametric method with the following assumptions:
# the rate at which a molecule breaks is proportional to length**N
# for some  power N
# and the point at which it breaks follows a Beta distribution with parameters A,B
# (this gives the a value from 0 to 1 denoting the point on the fragment to break)
# For symmetric breakage (most useful), you'll want A=B
# If A=B > 1 then breakage tends towards the middle. If 0 < A=B < 1 then tends towards the edges
# A = B = 5 gives reasonable values
# NOTE: there is no theoretical justification for why this should be an appropriate model
#       but it gives a reasonable looking length distribution while the uniform methods do not
def _compute_fragment_locations_beta(molecules, lambda_, N, A,B, runtime, rng):
    """
    fragment molecules with varying lambas

    returns a list of tuples (start, end, k) of fragments
    each fragment being k'th molecule from positions start to end (zero-based, non-inclusive of end)

    Parameters
    ----------
    molecules:
        a list of molecules to fragment
    lambda_:
        the rate of breakage of the bond between adjacent bases
    runtime:
        length of this process, longer times means more breakage

    """
    assert runtime > 0
    assert A > 0
    assert B > 0
    assert lambda_ > 0

    todo = collections.deque()
    todo.extend(((0, len(molecule)), runtime, k) for k, molecule in enumerate(molecules))
    done = collections.deque()

    while todo:
        # start, end are 0-based python-style half-open indexes
        (start, end), time_left, k = todo.pop()
        if end - start == 1:
            done.append(((start, end), k))
            continue


        num_bonds = end - start - 1
        break_rate = lambda_ * num_bonds**N / TYPICAL_MOLECULE_SIZE**(N-1)
        time_until_break = rng.exponential(scale = 1/break_rate)

        if time_until_break < time_left:
            # Break!

            # pick breakpoint by a beta distribution times the length of the molecule and round
            beta_variable = rng.beta(A,B)
            break_point = math.floor(beta_variable*(num_bonds)) + start
            if beta_variable == 1.0:
                # Surprisingly, rng.beta CAN return a 1.0
                # Which won't round down in the floor above, so we handle it here
                break_point = start+num_bonds - 1

            todo.append(((start, break_point + 1), time_left - time_until_break, k))
            todo.append(((break_point + 1, end), time_left - time_until_break, k))
        else:
            # No break!
            done.append( ((start, end), k) )

    return [(start, end, k) for (start, end), k in done]
