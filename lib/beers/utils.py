import pandas as pd
from molecule import Molecule

class Utils:
    base_complements = {"A":"T","T":"A","G":"C","C":"G"}

    @staticmethod
    def create_complement_strand(strand):
        """
        Simple utility to provide the complement of the strand and return it in the
        5' to 3' direction.  Note that T is used rather than U even for RNA
        :param strand: RNA/DNA strand to complement.
        :return: complement strand in 5' to 3' direction
        """
        complement_strand = ''.join(Utils.base_complements[base] for base in strand)
        return complement_strand[::-1]

    @staticmethod
    def convert_log_data_into_molecules(log_filename):
        molecules = []
        log_df = pd.read_csv(log_filename)
        log_retained_df = log_df[log_df["note"] != "removed"]
        for index, row in log_retained_df.iterrows():
            molecule = Molecule(row["id"], row["sequence"], row["start"], row["cigar"])
            molecules.append(molecule)
        return molecules


if __name__ == "__main__":
    #print(Utils.create_complement_strand("AAGTGACCTAAG"))

    Utils.convert_log_data_into_molecules("../../data/sizing_step.log")
