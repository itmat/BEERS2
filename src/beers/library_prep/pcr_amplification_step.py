import sys
from timeit import default_timer as timer
from copy import copy
import pickle

import numpy as np

from beers_utils.molecule import Molecule
from beers.utilities.gc_content import gc_content


def hypergeometric(ngood: int, nbad: int, nsamp: np.ndarray, rng: np.random.Generator):
    '''
    Sample from the hypergeometric distribution

    See np.random.hypergeometric but allows nsamp to have zeros
    '''

    out = np.zeros(shape = len(nsamp), dtype=int)
    nonzero = (nsamp != 0)
    out[nonzero] = rng.hypergeometric(ngood, nbad, nsamp[nonzero])
    return out

class PCRAmplificationStep:
    """
    PCR amplification step creates many copies of each molecule
    Since the number of molecules generated is quite large and
    since in real sequencing, only a small fraction of the molecules
    prepped end up forming clusters on the flowcell and being
    sequenced, we give the option to remove most of the PCR
    generated molecules prior to the end of the step.
    This provides a very large increase in speed and decrease
    in memory usage.

    Config Example::

        parameters:
            # Number of cycles to use
            # Generates 2^n fragments for each input fragment
            number_cycles: 10

            # Retain only this percent of the data
            # NOTE: scale is 0-100, not 0-1, so 0.08 means 0.0008 of the generated molecules are kep
            # This value can be approximated from real data with UMI tags by the following:
            # Compute the PCR duplication rate using the UMI tags.
            # Let N be the number of PCR steps. Assuming that each molecule generates 2^N
            # PCR descendants, if all were sequenced, we would expect all molecules to be duped
            # 2^N times. Instead, choose retention_percentage with the following strategy:
            #
            # '''
            # import scipy.stats
            # at_least_two = scipy.stats.binom(p=retention_percentage / 100, n=2**N).sf(1)
            # at_least_one = scipy.stats.binom(p=retention_percentage / 100, n=2**N).sf(0)
            # dupe_rate = at_least_two / at_least_one
            # # Choose retention_percentage so that dupe_rate is approximately the observed PCR dupe rate
            # '''
            #
            # For example, retention_percentage = 0.08 and number_cycles=10 gives dupe rate
            # of about 35%, which is pretty typical. Note that we are defining dupe rate as
            # the fraction of all reads after deduping that had at least one duplicate.
            # If number_cycles changes, this should be modified too to maintain dupe rate.
            retention_percentage: 0.08

            # Induce a GC bias by discarding some PCR duplicates
            # according to their GC bias. All fragments overall
            # GC content is computed and then the following three
            # parameters are used to compute a probability of retention
            # gc_bias_constant + gc_bias_linear*(gc - 0.5) + gc_bias_quadratic*(gc - 0.5)^2
            # clipped to always be within 0 and 1
            # For no bias, set to 1, 0, and 0 for const, linear, and quadratic, respectively.
            # To bias towards GC=0.5 content, set gc_bias_quadratic to be negative,
            # such as -100 for a large GC bias.
            gc_bias_constant: 1.0
            gc_bias_linear: 0.0
            gc_bias_quadratic: -100

            # During PCR, we have some chance of mis-copying the molecules
            # These are specified here as probability per-base, so 0.001 is one error per 1kb
            deletion_rate: 0.0001
            insertion_rate: 0.0001
            substitution_rate: 0.001
    """

    name = "PCR Amplification Step"

    # Putting an upper bound on the number of allowed PCR amplification cycles to avoid accidentally launching
    # a virtually runaway population of molecule objects.
    MAX_CYCLE_NUMBER = 16

    def __init__(self, parameters, global_config):
        self.number_cycles = parameters.get("number_cycles")
        self.retention_percentage = parameters.get('retention_percentage')
        self.substitution_rate = parameters.get('substitution_rate')
        self.insertion_rate = parameters.get('insertion_rate')
        self.deletion_rate = parameters.get('deletion_rate')
        self.parameters = parameters
        self.global_config = global_config

        # The sample id counter maintains a counter for each of the molecules input into this step so that
        # unique ids may be assigned to those molecule's descendants.
        self.sample_id_ctr = dict()

    def execute(self, molecule_packet, rng, log):
        """
        For each cycle, duplicates each molecule from the current list of molecules (including those created in
        prior cycles) and assigns a unique id.
        :param molecule_packet: packet containing sample info and molecules from the prior step
        :return: original molecule packet but with molecules from this step
        """
        print("PCR Amplication step starting")

        # Initialize the sample id counter for all starting molecules, using the molecule id as the key.
        self.sample_id_ctr = dict.fromkeys([molecule.molecule_id for molecule in molecule_packet.molecules], 0)

        # The sample represents the input.  The molecules variable starts with the sample input and will grow
        # with each completed cycle.
        molecules = molecule_packet.molecules

        # For each cycle, we generate copies. However, since PCR grows
        # the number of molecule exponentially and the next step after PCR
        # is always a downsampling to much lower numbers (hence why PCR must
        # be performed), we compute ahead of time how many descendants of
        # each molecule will be retained at the end. Then we only retain and
        # clone molecules as necessary, dropping off any that will end up
        # not being retained. Thereby, much computation time is saved, but
        # the tree structure of copies-and-originals with descendants that are
        # copies-of-copies is retained.

        # Compute the number of retained descendants each molecule will have
        # These are computed under the assumption of no GC-content bias
        # GC-bias will further reduce the actual numbers later on
        retention_rate = self.retention_percentage / 100
        descendants = rng.binomial(n=2**self.number_cycles, p = retention_rate, size=len(molecules))

        # Rate of success (per molecule) of each PCR cycle
        gc = np.array([gc_content(molecule) for molecule in molecules])
        gc_const = self.parameters['gc_bias_constant']
        gc_linear = self.parameters['gc_bias_linear']
        gc_quadratic = self.parameters['gc_bias_quadratic']
        success_rate = np.clip(gc_const + gc_linear*(gc - 0.5) + gc_quadratic*(gc - 0.5)**2, 0, 1)

        # Iterate over the number of cycles requested.
        for cycle in range(self.number_cycles):

            print(f"Starting cycle {cycle + 1}")


            # Each of the assigned descendants comes from either the original molecule
            # or from the copy of it we make this round.
            # We divvy those descendants up between original or copy 50/50
            possible_descendants = 2**(self.number_cycles - cycle - 1)
            original_descendants = hypergeometric(
                    possible_descendants,
                    possible_descendants,
                    descendants,
                    rng)
            copy_descendants = descendants - original_descendants

            # Account for failure in PCR
            pcr_succeeded = rng.random(success_rate.shape) < success_rate

            # Retain molecules if the original has descendants (not from the copy made this round)
            new_molecules = [molecule
                                for desc, molecule
                                in zip(original_descendants, molecules)
                                if desc > 0]

            # Create copy molecules if they will have at least one descendant
            for desc, molecule, pcr_success in zip(copy_descendants, molecules, pcr_succeeded):
                if desc == 0:
                    continue

                if not pcr_success:
                    # This molecule failed to PCR duplicate this round
                    continue

                # Copy the original molecule
                new_molecule = copy(molecule)

                # Add any errors in copying
                new_molecule.generate_errors(
                        self.substitution_rate,
                        self.insertion_rate,
                        self.deletion_rate,
                        rng,
                )

                # TODO: copies should be reverse-complemented, right?
                self.assign_id(new_molecule, molecule.molecule_id)
                new_molecules.append(new_molecule)

            # Descendents of the final molecules from this step were already computed
            # so we assemble them here: original molecules go first and copies second
            descendants = np.concatenate((
                original_descendants[original_descendants > 0],
                copy_descendants[(copy_descendants > 0) & pcr_succeeded],
            ))

            # success_rates of the copied molecules are assumed to be
            # the same as the original. This is a simplification since
            # errors in PCR copies could influence GC bias, but should
            # be a very close approximation.
            success_rate = np.concatenate((
                success_rate[original_descendants > 0],
                success_rate[(copy_descendants > 0) & pcr_succeeded],
            ))

            molecules = new_molecules

        # Output the updated list of molecules to the log file
        for molecule in molecules:
            log.write(molecule)

        # Return the new list of molecules
        molecule_packet.molecules = molecules

        return molecule_packet

    def assign_id(self, new_molecule, ancestor_id):
        """
        Molecule ids read as A.B.C.D and so on where the dot indicates a new step in which new molecules are
        generated.  Since we do not know how many former steps have created molecules and since this step will
        only add one dot, a molecule created in this step will have an id of its parent plus a dot plus a number.
        Original molecules (parent molecule) will only have its original id.  So we need to back up at most one
        ".nnn" to find the original id and we can use that to do a lookup and increment in the sample id counter
        held in this step.
        :param new_molecule: freshly created molecule
        :param ancestor_id: parent id of the new molecule (which may or may not be from the input sample)
        """

        # If the current ancestor's id is in the sample id counter, we are at the id of the input molecule.
        if ancestor_id not in self.sample_id_ctr:
            # Otherwise, we need to back up the id string to the last occurence of the dot and check again.
            index = ancestor_id.rfind(".")
            assert index > 0, f"cannot locate the ancestor of {new_molecule.molecule_id}"
            ancestor_id = ancestor_id[:index]
        # The sample id counter has the number we can assign to the new id of this molecule, after
        # which, we bump the sample id counter and return
        new_molecule.molecule_id = f"{ancestor_id}.{self.sample_id_ctr[ancestor_id]}"
        self.sample_id_ctr[ancestor_id] += 1

    @staticmethod
    def validate(parameters, global_config):
        """
        Insure that the parameters supplied for this step are appropriate.  In this case, the number of cycles
        parameter should be a counting integer no greater than the number of maximum cycles allowed.
        :return: True if validation passed, False otherwise
        """
        errors = []


        if 'number_cycles' not in parameters:
            errors.append("Must specify 'number_cycles'")
        elif not isinstance(parameters['number_cycles'],int):
            errors.append("'number_cylces' must be an integer")
        elif parameters['number_cycles'] < 0 or parameters['number_cycles'] > PCRAmplificationStep.MAX_CYCLE_NUMBER:
            errors.append(
                  f"The cycle number parameter must be greater than 0 and less"
                  f" than {PCRAmplificationStep.MAX_CYCLE_NUMBER}")

        if 'retention_percentage' not in parameters:
            errors.append("Must specify 'retention_percentage'")
        elif parameters['retention_percentage'] <= 0 or parameters['retention_percentage'] > 100:
            errors.append(
                    f"The retention_percentage, {parameters['retention_percentage']},"
                   f" must be between 0 and 100.")

        for var in ['substitution_rate', 'insertion_rate', 'deletion_rate']:
            if var not in parameters:
                errors.append(f"Must specify {var}")
            elif not (0 <= parameters[var] <= 1):
                errors.append(f"{var} must be betwen 0 and 1")

        for var in ['gc_bias_constant', 'gc_bias_linear', 'gc_bias_quadratic']:
            if var not in parameters:
                errors.append(f"Must specify {var}")
            elif not isinstance(parameters[var], (float, int)):
                errors.append(f"{var} must be a number")

        return errors
