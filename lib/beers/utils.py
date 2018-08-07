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
        complement_strand = ''
        for base in strand:
            complement_strand += Utils.base_complements[base]
        return complement_strand[::-1]


if __name__ == "__main__":
    print(Utils.create_complement_strand("AAGTGACCTAAG"))
