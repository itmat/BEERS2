import time as time
import os
import glob

class GeneralUtils:

    NUMBER_OF_FILE_PER_DIRECTORY = 200

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
    def create_subdirectories(packet_count, parent_directory_path):
        max_packet_group = packet_count // GeneralUtils.NUMBER_OF_FILE_PER_DIRECTORY
        log_directory_path = os.path.join(parent_directory_path, "logs")
        data_directory_path = os.path.join(parent_directory_path, 'data')
        for packet_group in range(max_packet_group + 1):
            log_subdirectory_path = os.path.join(log_directory_path, f'pkt_grp{packet_group}')
            data_subdirectory_path = os.path.join(data_directory_path, f'pkt_grp{packet_group}')
            os.makedirs(log_subdirectory_path, mode=0o0755, exist_ok=True)
            os.makedirs(data_subdirectory_path, mode=0o0755, exist_ok=True)
        return log_subdirectory_path, data_subdirectory_path

    @staticmethod
    def get_output_subdirectories(packet_id, output_directory_path):
        log_directory_path = os.path.join(output_directory_path, "logs")
        data_directory_path = os.path.join(output_directory_path, 'data')
        packet_group = packet_id // GeneralUtils.NUMBER_OF_FILE_PER_DIRECTORY
        log_subdirectory_path = os.path.join(log_directory_path, f'pkt_grp{packet_group}')
        data_subdirectory_path = os.path.join(data_directory_path, f'pkt_grp{packet_group}')
        return log_subdirectory_path, data_subdirectory_path


if __name__ == "__main__":
    pass

