import time as time
import os
import glob
import numpy as np

class GeneralUtils:
    """
    Utility scripts expected to have wide applicability across the various modules in this app.
    """

    FILES_PER_DIRECTORY_LIMIT = 100

    @staticmethod
    def generate_seed():
        """
        Provides an integer of 32 bits or less using a seconds based timestamp.  If the timestamp exceeds 32 bits, the
        integer obtained from the lowest 32 bits is returned
        :return: a 32 bit integer seed for the numpy random number generator.
        """
        candidate_seed = int(time.time())

        # If the timestamp bit length exceeds 32 bits, mask out all but the lowest 32 bits.
        if candidate_seed.bit_length() > 32:
            candidate_seed = candidate_seed & 0xffffffff
        return candidate_seed

    base_complements = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    # TODO: should we allow 'N's to be complemented? - Yes when reverse complementing reads

    @staticmethod
    def create_complement_strand(strand):
        """
        Simple utility to provide the complement of the strand and return it in the
        5' to 3' direction.  Note that T is used rather than U even for RNA
        :param strand: RNA/DNA strand to complement.
        :return: complement strand in 5' to 3' direction
        """
        complement_strand = ''.join(GeneralUtils.base_complements[base] for base in strand)
        return complement_strand[::-1]

    @staticmethod
    def create_subdirectories(file_count, directory_path):
        """
        This method creates a series of subdirectories extending from the given directory path organized such as to
        avoid exceeding the FILES_PER_DIRECTORY_LIMIT in any given subdirectory.  The number of files to be housed
        under the given directory path is used to determine the number of and nesting level of the subdirectories
        needed.  A string representation of the resulting directory structure is returned.
        :param file_count: The number of files to be housed under the given directory path
        :param directory_path: the path of the parent directory to house the expected number of files
        :return: A string representation of the number of level of nesting and the number of directories on each level.
        """
        directory_structure = {}
        GeneralUtils.create_nested_subdirectories(file_count, directory_structure, directory_path, 1)
        directory_structure = [directory_structure[key] for key in sorted(directory_structure.keys())]
        return ",".join(directory_structure)

    @staticmethod
    def create_nested_subdirectories(file_count, directory_structure, directory_path, depth):
        full_subdirectories, remaining_files = divmod(file_count, GeneralUtils.FILES_PER_DIRECTORY_LIMIT)
        partial_subdirectory = 1 if remaining_files else 0
        number_subdirectories_to_create = min(GeneralUtils.FILES_PER_DIRECTORY_LIMIT,
                                              full_subdirectories + partial_subdirectory)
        if depth not in directory_structure.keys():
            directory_structure[depth] = str(number_subdirectories_to_create)
        for subdirectory_id in range(number_subdirectories_to_create):
            remaining_file_count = file_count - number_subdirectories_to_create * GeneralUtils.FILES_PER_DIRECTORY_LIMIT
            subdirectory_path = os.path.join(directory_path, str(subdirectory_id))
            os.makedirs(subdirectory_path, mode=0o0755, exist_ok=True)
            if remaining_file_count > 0:
                GeneralUtils.create_nested_subdirectories(remaining_file_count,
                                                          directory_structure,
                                                          subdirectory_path,
                                                          depth + 1)

    @staticmethod
    def get_output_subdirectories(pkt_id, directory_structure):
        """
        This returns a list of the subdirectory names intervening between the user's output directory (data, logs, etc)
        and the actual packet file data, log. etc and is used to locate where the actual file should be written.
        :param pkt_id:  The pkt id of the pkt file is used to locate the file in the nested hierarchy.  It is
        treated as the (pkt_id + 1)th file to write
        :param directory_structure: string representation of the directory structure.  The number of levels are
        represented by the number of comma delimited values and the number of directories on each level are given
        by the values themselves.
        :return: a list of the sbudirectory names to be contatenated with the user's specific output directory to
        create a full path to the packet file data, log, etc
        """
        directory_counts = [int(directory_count) for directory_count in directory_structure.split(",")]
        subdirectory_list = []
        remaining_file_count = pkt_id
        # Handles edge case where everything fits into one directory.
        level_file_count = max(np.product(directory_counts),GeneralUtils.FILES_PER_DIRECTORY_LIMIT)
        nesting_depth = len(directory_counts)
        for _ in range(nesting_depth):
            assert level_file_count > 0, f"Too many nesting levels {nesting_depth}"
            level_path, remaining_file_count = divmod(remaining_file_count, level_file_count)
            subdirectory_list.append(str(level_path))
            level_file_count //= directory_counts.pop(-1)
        return subdirectory_list

class BeersUtilsException(Exception):
    """Base class for other Utils exceptions."""
    pass

if __name__ == "__main__":
    pass
