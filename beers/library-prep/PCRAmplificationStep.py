import numpy as np
from molecule import Molecule
import sys
from utils import Utils
from timeit import default_timer as timer
from copy import copy


class PCRAmplificationStep:

    # Putting an upper bound on the number of allowed PCR amplification cycles to avoid accidentally launching
    # a virtually runaway population of molecule objects.
    MAX_CYCLE_NUMBER = 12

    def __init__(self, log_filename, parameters):
        """
        Instantiates the step with a log file name and a list of parameters.  For PCR Amplification, the only
        parameter at this point is the number of cycles.
        :param log_filename: path to log file
        :param parameters: parameter json object
        """
        self.log_filename = log_filename
        self.number_cycles = parameters.get("number_cycles")

        # The sample id counter maintains a counter for each of the molecules input into this step so that
        # unique ids may be assigned to those molecule's descendents.
        self.sample_id_ctr = dict()

    def execute(self, sample):
        """
        For each cycle, duplicates each molecule from the current list of molecules (including those created in
        prior cycles) and assigns a unique id.
        :param sample: molecules from the prior step
        :return: molecules from this step
        """
        print("PCR Amplication step starting")

        # Initialize the sample id counter for all starting molecules, using the molecule id as the key.
        self.sample_id_ctr = dict.fromkeys([molecule.molecule_id for molecule in sample], 0)

        # The sample represents the input.  The molecules variable starts with the sample input and will grow
        # with each completed cycle.
        molecules = sample

        # Open the log file for writing only (overwrite any existing log file of the same name/location)"
        with open(self.log_filename, "w") as log_file:

            # Add a header to the log file.
            log_file.write(Molecule.header)

            # Iterate over the number of cycles requested.
            for cycle in range(self.number_cycles):

                # Start a new empty list of new molecules
                new_molecules = []

                # Iterate over the current list of molecules
                for molecule in molecules:

                    # Clone the original molecule (no biases as yet being introduced)
                    new_molecule = copy(molecule)

                    # Get a unique id for the new molecule that is derived from its parent
                    self.assign_id(new_molecule, molecule.molecule_id)

                    # Add the new molecule to the list of new molecules
                    new_molecules.append(new_molecule)

                # Dump the newly created molecules into the list of current molecules so that they will
                # also be cloned during subsequent cycles.
                molecules.extend(new_molecules)

            # Output the updated list of molecules to the log file
            for molecule in molecules:
                log_file.write(molecule.log_entry())

        # Return the new list of molecules
        return molecules

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

        # At most, we will need to iterate twice (once if the new molecule has an original molecule as its
        # parent and twice if the new molecule is the offspring of a new molecule created in a prior cycle.
        for _ in range(2):

            # If the current ancestor's id is in the sample id counter, we are at the id of the input molecule.
            if ancestor_id in self.sample_id_ctr:

                # The sample id counter has the number we can assign to the new id of this molecule, after
                # which, we bump the sample id counter and return
                new_molecule.id = str(ancestor_id) + "." + str(self.sample_id_ctr[ancestor_id])
                self.sample_id_ctr[ancestor_id] += 1
                break

            # Otherwise, we need to back up the id string to the last occurence of the dot and check again.
            index = ancestor_id[:ancestor_id.rfind(".")]
            assert index < 0, f"cannot locate the ancestor of {new_molecule.molecule_id}"
            ancestor_id = ancestor_id[:index]
        assert True, f"cannot locate the ancestor of {new_molecule.molecule_id}"

    def validate(self):
        """
        Insure that the parameters supplied for this step are appropriate.  In this case, the number of cycles
        parameter should be a counting integer no greater than the number of maximum cycles allowed.
        :return: True if validation passed, False otherwise
        """
        print(f"PCR Amplification step validating parameters")
        if self.number_cycles < 0 or self.number_cycles > PCRAmplificationStep.MAX_CYCLE_NUMBER:
            print(f"The cycle number parameter must be greater than 0 and less"
                  f" than {PCRAmplificationStep.MAX_CYCLE_NUMBER}", file=sys.stderr)
            return False
        return True


if __name__ == "__main__":
    # Single step test - using a seed to allow reproducible runs.
    np.random.seed(100)

    # Taking advantage of an existing log file to grab molecules.
    input_sample = Utils.convert_log_data_into_molecules("../../data/sizing_step.log")

    # Copying these molecules into a separate log file
    input_data_log_file = "../../data/pcr_amplification_step_input_data.log"
    with open(input_data_log_file, "w+") as input_data_log:
        input_data_log.write(Molecule.header)
        for rna_molecule in input_sample:
            input_data_log.write(rna_molecule.log_entry())

    # Selecting step log file and parameter info and using both to instantiate a step
    # object (not bothering with validation)
    output_data_log_file = "../../data/pcr_amplification_step_output_data.log"
    input_parameters = {
        "number_cycles": 8
    }
    step = PCRAmplificationStep(output_data_log_file, input_parameters)

    # Executing the step and noting the time taken.
    start = timer()
    result = step.execute(input_sample)
    end = timer()
    print(f"PCRAmplification Step: {end - start} sec for {len(result)} molecules.")