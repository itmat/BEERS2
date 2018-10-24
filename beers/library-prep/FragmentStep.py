import numpy
import math
import collections
from molecule import Molecule
import sys


METHODS = {"uniform", "beta"}

# TODO: should this be a parameter?
# For fragment2
# just makes lambda into a reasonable value for a molecule of the "typical" size
# i.e. lambda gives the rate of breakage of a molecule which has this length
TYPICAL_MOLECULE_SIZE = 1000

class FragmentStep:
    def __init__(self, logfile, parameters):
        """
        Make fragmentation step that with rate of fragmentation lambda_ running for run_time

        lambda_ -- either a value > 0 in which case all positions in each molecule are equally likely to break
                or a function (k, (start,end), molecule) -> rate
                that gives a rate of breaking (in breaks/unit time) of the bond between nucleotides
                at positions k and k+1 of the fragment that goes from (start,end) in molecule
        runtime -- length of running the experiment (Has the inverse units of lambda_)
        """
        # TODO: these parameters should be read in from a config file
        self.history_filename = logfile
        self.method = parameters["method"]
        self.lambda_ = parameters["lambda"]
        self.runtime = parameters["runtime"]

        # Parameters used ONLY for beta_fragmentation method
        if self.method == "beta":
            self.beta_A = parameters["beta_A"]
            self.beta_B = parameters["beta_B"]
            self.beta_N = parameters["beta_N"]

    def execute(self, sample):
        print("Fragment Step acting on sample")
        if self.method == "uniform":
            fragment_locations = fast_compute_fragment_locations(sample, self.lambda_, self.runtime)
        elif self.method == "beta":
            fragment_locations = compute_fragment_locations_beta(sample, self.lambda_, self.beta_N, self.beta_A, self.beta_B, self.runtime)

        result = [sample[k].make_fragment(start,end) for (start, end, k) in fragment_locations]

        # Output all the molecules to log
        with open(self.history_filename, "w+") as log_file:
            log_file.write(Molecule.header)
            for molecule in result:
                log_file.write(molecule.log_entry())

        return result

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
        return True


def estimate_uniform_lambda(starting_median_length, desired_median_length):
    """Computes a (lambda_, rate) pair that give the desired amount of fragmentation"""
    #TODO: verify correctness of this
    runtime = 1 # Chosing this just rescales the lambda_ (choosing different units for time)
    lambda_ = starting_median_length/desired_median_length

    return lambda_, runtime

### FRAGMENTATION METHODS

# BEST GENERAL SOLUTION
# Much slower than uniform fragmentation like direct_fragment
def compute_fragment_locations(molecules, lambda_, runtime):
    """fragment molecules with varying lambas

    molecules -- a list of molecules
    lambda_ -- a function mapping (j, (start,end), molecule) to the rate of breakage of
            the bond at position j of the molecule if it is in a fragment (start,end) of molecule
    runtime -- length of this process, longer times means more breakage

    returns a list of tuples (start, end, k) of fragments
    each fragment being k'th molecule from positions start to end (non-inclusive of end)
    """
    assert runtime > 0

    todo = collections.deque()
    todo.extend(((0, len(molecule)), runtime, k) for k, molecule in enumerate(molecules))
    done = collections.deque()

    while todo:
        (start, end), time_left, k = todo.pop()
        if end - start == 1:
            done.append(((start, end), k))
            continue


        num_bonds = end - start - 1
        lambdas = numpy.array([lambda_(j,start, end, molecules[k]) for j in range(start, end-1)])
        total_lambda = sum(lambdas)
        time_until_break = numpy.random.exponential(scale = 1/total_lambda)

        if time_until_break < time_left:
            # Break!
            break_point = numpy.random.choice(num_bonds, p=lambdas/total_lambda) + start
            todo.append(((start, break_point + 1), time_left - time_until_break, k))
            todo.append(((break_point + 1, end), time_left - time_until_break, k))
        else:
            # No break!
            done.append( ((start, end), k) )

    return [(start, end, k) for (start, end), k in done]

# The general idea, but slow method. Operates on one fragment at a time.
# Prefer fragment() above. Not of the same format as other functions in this module
def simple_fragment(molecule, lambda_function, runtime):
    """fragment any iterable `molecule` into some number of pieces.

    Prefer fragment() function instead
    Each inter-base-pair bond has a rate lambda of breaking determined by
    lambda = lambda_function(k, molecule)
    where k is the position of the bond in a molecule of given length"""

    assert runtime > 0

    # Can't fragment a single base
    if len(molecule) == 1:
        yield molecule
        return

    num_bonds = len(molecule) - 1
    lambda_array = numpy.array([lambda_function(k, molecule) for k in range(num_bonds)])
    total_lambda = numpy.sum(lambda_array)
    time_until_break = numpy.random.exponential(scale=1/total_lambda)

    if time_until_break < runtime:
        # Break!
        break_point = numpy.random.choice(num_bonds, p = lambda_array/total_lambda)
        yield from simple_fragment(molecule[:break_point+1], lambda_function, runtime - time_until_break)
        yield from simple_fragment(molecule[break_point+1:], lambda_function, runtime - time_until_break)
    else:
        # No break!
        yield molecule

# used in fast_fragment
def random_range(mins, maxes):
    """Return random numpy array of random integers with varying mins,maxes"""
    #Surprisingly, numpy doesn't have a built-in function for this
    # numpy.random.randomint gives each sample the same min and max
    return numpy.floor(numpy.random.random(mins.shape)*(maxes - mins) + mins)

# Fastest UNIFORM fragmentation
# But is harder to read than direct_fragment() and only very slightly faster (10%)
def fast_compute_fragment_locations(molecules, lambda_, runtime):
    """ Each molecule is fragment randomly according

    molecules -- the list of all molecules
    lambda_ -- the rate of breaking of a single bond between base adjacent base pairs
    runtime -- how long the fragmentation is run for

    The algorithm used to fragment them is the following:
    on each molecule, the time until the next base fragments is distributed as
    an exponential with rate (lambda) the sum of the rates of all of its bonds.
    So we cycle through all molecules, randomly taking the next base to fragment from each molecule
    so long as the time this takes to fragment is less than the runtime of the experiment.
    Then repeat on each half of each molecule that fragmented.
    """

    assert runtime > 0
    assert lambda_ > 0

    batch = numpy.array([(0, len(molecule), runtime, k) for k, molecule in enumerate(molecules)])
    done = numpy.array([])
    done.shape = (0,4)

    while batch.size > 0:
        # Potentially break each molecule of the batch
        # then add those that broke to the second batch
        # and those that don't break to the 'done' list

        starts = batch[:,0]
        ends = batch[:,1]
        time_lefts = batch[:,2]
        ks = batch[:,3]

        num_bonds = ends - starts - 1
        valid = (num_bonds > 0)
        total_lambdas = lambda_ * num_bonds

        times_until_break = numpy.full(shape = starts.shape, fill_value=numpy.inf)
        times_until_break[valid] = numpy.random.exponential(scale = 1/total_lambdas[valid])
        #Note: everything that is not valid (i.e. is a single base and hence can't break)
        # will have infinite breaking time and hence won't break

        remaining_time = time_lefts - times_until_break
        broke = (remaining_time > 0)


        # Any that didn't break are put on the done list
        done = numpy.append(done, batch[numpy.logical_not(broke),:], axis=0)

        # Now break the ones that did break
        break_points = random_range(starts, ends - 1)

        # The left halves of broken
        left_starts = starts[broke]
        left_ends = break_points[broke]+1
        left_time_lefts = remaining_time[broke]
        left_ks = ks[broke]
        lefts = numpy.vstack( (left_starts, left_ends, left_time_lefts, left_ks) ).T

        # The right halves
        right_starts = break_points[broke]+1
        right_ends = ends[broke]
        right_time_lefts = remaining_time[broke]
        right_ks = ks[broke]
        rights = numpy.vstack( (right_starts, right_ends, right_time_lefts, right_ks) ).T

        # New batch to potentially break
        batch = numpy.append(lefts, rights, axis=0)

    return [(int(start), int(end), int(k)) for start, end, time_left, k in done]


# Used by direct_fragment
def sample_without_replacement(n, k):
    """uniformly sample k numbers from 0,...,n-1 without replacement

    Intended to be faster than numpy.random.choice(n, size=k, replace=False)
    for the case when n is large and k is small.
    (About 1000x times faster for n=1,000,000 and k=4.)
    Sample is returned in sorted order """

    sample = set(numpy.random.choice(n, size=k, replace=True))
    while len(sample) < k:
        sample = sample.union( numpy.random.choice(n, size=(k-len(sample)), replace=True) )

    # Sort it since list(set) will give a sort-of arbitrary but not random order
    return sorted(sample)

### Simplest and very fast method for UNIFORM fragmentation
# around 10% slower than fast_fragment
# i.e. all bonds have the same chance of breaking
def uniform_direct_compute_fragment_locations(molecules, lambda_, runtime):
    """uniform lambda_ parameter fragmentation

    Use the fact that the methods of fragment is equivalent (when lambda constant)
    to first sampling the number of break points from a binomial distribution and then
    picking those points uniformly on the molecule"""
    assert runtime > 0
    assert lambda_ > 0

    output = collections.deque()
    for k, molecule in enumerate(molecules):
        num_bonds = len(molecule) - 1
        probability_of_base_breaking = lambda_ * runtime

        num_breakpoints = numpy.random.binomial(n=num_bonds, p=probability_of_base_breaking)

        breakpoints = sample_without_replacement(num_bonds, num_breakpoints)
        bps = [0]+ breakpoints + [len(molecule)]
        output.extend( (bps[i], bps[i+1], k) for i in range(len(bps)-1) )

    return list(output)

### Parametric method with the following assumptions:
# the rate at which a molecule breaks is proportional to length**N
# for some  power k
# and the point at which it breaks follows a Beta distribution with parameters A,B
# (this gives the a value from 0 to 1 denoting the point on the fragment to break)
# For symmetric breakage (most useful), you'll want A=B
# If A=B > 1 then breakage tends towards the middle. If 0 < A=B < 1 then tends towards the edges
# A = B = 5 gives reasonable values
# NOTE: there is no theoretical justification for why this should be an appropriate model
#       but it gives a reasonable looking length distribution while the uniform methods do not
def compute_fragment_locations_beta(molecules, lambda_, N, A,B, runtime):
    """fragment molecules with varying lambas

    molecules -- a list of molecules
    lambda_ -- a function mapping (j, (start,end), molecule) to the rate of breakage of
            the bond at position j of the molecule if it is in a fragment (start,end) of molecule
    runtime -- length of this process, longer times means more breakage

    returns a list of tuples (start, end, k) of fragments
    each fragment being k'th molecule from positions start to end (non-inclusive of end)
    """
    assert runtime > 0
    assert A > 0
    assert B > 0
    assert lambda_ > 0

    todo = collections.deque()
    todo.extend(((0, len(molecule)), runtime, k) for k, molecule in enumerate(molecules))
    done = collections.deque()

    while todo:
        (start, end), time_left, k = todo.pop()
        if end - start == 1:
            done.append(((start, end), k))
            continue


        num_bonds = end - start - 1
        break_rate = lambda_ * num_bonds**N / TYPICAL_MOLECULE_SIZE**(N-1)
        time_until_break = numpy.random.exponential(scale = 1/break_rate)

        if time_until_break < time_left:
            # Break!

            # pick breakpoint by a beta distribution times the length of the molecule and round
            beta_variable = numpy.random.beta(A,B)
            break_point = math.floor(beta_variable*(num_bonds)) + start
            if beta_variable == 1.0:
                # Surprisingly, numpy.random.beta CAN return a 1.0
                # Which won't round down in the floor above, so we handle it here
                break_point = start+num_bonds - 1

            todo.append(((start, break_point + 1), time_left - time_until_break, k))
            todo.append(((break_point + 1, end), time_left - time_until_break, k))
        else:
            # No break!
            done.append( ((start, end), k) )

    return [(start, end, k) for (start, end), k in done]