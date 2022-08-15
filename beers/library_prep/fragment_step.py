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

    name = "Fragment Step"

    def __init__(self, logfile, parameters, global_config):
        """
        Make fragmentation step that with rate of fragmentation lambda_ running for run_time

        lambda_ -- either a value > 0 in which case all positions in each molecule are equally likely to break
                or a function (k, (start,end), molecule) -> rate
                that gives a rate of breaking (in breaks/unit time) of the bond between nucleotides
                at positions k and k+1 of the fragment that goes from (start,end) in molecule
        runtime -- length of running the experiment (Has the inverse units of lambda_)
        """
        self.global_config = global_config
        self.history_filename = logfile
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

    def execute(self, molecule_packet, rng):
        print("Fragment Step acting on sample")
        sample = molecule_packet.molecules
        if self.method == "uniform":
            fragment_locations =  compute_fragment_locations_uniform(sample, self.lambda_, self.runtime, rng)
        elif self.method == "beta":
            fragment_locations = compute_fragment_locations_beta(sample, self.lambda_, self.beta_N, self.beta_A, self.beta_B, self.runtime, rng)
        else:
            raise NotImplementedError(f"Unknown fragmentation method {self.method}")

        result = [sample[k].make_fragment(start+1,end) for (start, end, k) in fragment_locations
                    if end - start >= self.min_frag_size]

        # Output all the molecules to log
        with open(self.history_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in result:
                log_file.write(molecule.log_entry())

        molecule_packet.molecules = result
        return molecule_packet

    def validate(self):
        if self.lambda_ <= 0:
            print("Rate of fragmentation must be positive", file=sys.stderr)
            return False
        if self.runtime <= 0:
            print("Fragmentation runtime must be a positive number", file=sys.stderr)
            return False
        if self.method not in METHODS:
            print("Fragmentation method must be one of the following" + str(METHODS), file=sys.stderr)
            return False
        if self.min_frag_size <= 0:
            print("min_frag_size must be a positive number", file=sys.stderr)
            return False
        return True


def estimate_uniform_lambda(starting_median_length, desired_median_length):
    """Computes a (lambda_, rate) pair that give the desired amount of fragmentation"""
    #TODO: verify correctness of this
    runtime = 1 # Chosing this just rescales the lambda_ (choosing different units for time)
    lambda_ = starting_median_length/desired_median_length

    return lambda_, runtime

### FRAGMENTATION METHODS

# MOST GENERAL SOLUTION
# Much slower than uniform fragmentation
# Unlike other methods, allows an arbitrary lambda_ depending upon location
def compute_fragment_locations(molecules, lambda_, runtime, rng):
    """fragment molecules with varying lambas

    molecules -- a list of molecules
    lambda_ -- a function mapping (j, (start,end), molecule) to the rate of breakage of
            the bond at position j of the molecule if it is in a fragment (start,end) of molecule
    runtime -- length of this process, longer times means more breakage
    rng -- random number generator

    returns a list of tuples (start, end, k) of fragments
    each fragment being k'th molecule from positions start to end (zero-based, non-inclusive of end)
    """
    assert runtime > 0

    todo = collections.deque()
    todo.extend(((0, len(molecule)), runtime, k) for k, molecule in enumerate(molecules))
    done = collections.deque()

    while todo:
        # start, end are zero based python slice-style half-open
        (start, end), time_left, k = todo.pop()
        if end - start == 1:
            done.append(((start, end), k))
            continue


        num_bonds = end - start - 1
        lambdas = numpy.array([lambda_(j,start, end, molecules[k]) for j in range(start, end-1)])
        total_lambda = sum(lambdas)
        time_until_break = rng.exponential(scale = 1/total_lambda)

        if time_until_break < time_left:
            # Break!
            break_point = rng.choice(num_bonds, p=lambdas/total_lambda) + start
            todo.append(((start, break_point + 1), time_left - time_until_break, k))
            todo.append(((break_point + 1, end), time_left - time_until_break, k))
        else:
            # No break!
            done.append( ((start, end), k) )

    return [(start, end, k) for (start, end), k in done]

# Used by direct_fragment
def sample_without_replacement(n, k, rng):
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

def compute_fragment_locations_uniform(molecules, lambda_, runtime, rng):
    """uniform fragmentation with a rate lambda_ parameter
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

        breakpoints = sample_without_replacement(num_bonds, num_breakpoints, rng)
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
def compute_fragment_locations_beta(molecules, lambda_, N, A,B, runtime, rng):
    """fragment molecules with varying lambas

    molecules -- a list of molecules
    lambda_ -- the rate of breakage of the bond between adjacent bases
    runtime -- length of this process, longer times means more breakage

    returns a list of tuples (start, end, k) of fragments
    each fragment being k'th molecule from positions start to end (zero-based, non-inclusive of end)
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
