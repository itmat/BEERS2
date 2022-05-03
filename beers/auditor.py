import os
import fcntl
import time
from beers_utils.constants import CONSTANTS
from beers.beers_exception import BeersException


class Auditor:

    def __init__(self, total, output_directory):
        self.total = total
        self.audit_file_path = os.path.join(output_directory, CONSTANTS.AUDIT_FILENAME)
        with open(self.audit_file_path, 'w'):
            pass
        self.packets_submitted = []

    def submit_packet_id(self, packet_id):
        self.packets_submitted.append(packet_id)

    def is_processing_complete(self):
        for ctr in range(1000):
            try:
                with open(self.audit_file_path, 'r') as audit_file:
                    fcntl.flock(audit_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    count = len(audit_file.readlines())
                    fcntl.flock(audit_file, fcntl.LOCK_UN)
                    return True if count == self.total else False
            except IOError:
                time.sleep(0.05)
        else:
            raise BeersException(
                f"Unable to access audit file, {self.audit_file_path}, for reading within a reasonable period.")

    @staticmethod
    def note_packet_completed(packet_id, audit_directory):
        audit_file_path = os.path.join(audit_directory, CONSTANTS.AUDIT_FILENAME)
        for ctr in range(1000):
            try:
                with open(audit_file_path, 'a') as audit_file:
                    fcntl.flock(audit_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    audit_file.write(f"{packet_id}\n")
                    fcntl.flock(audit_file, fcntl.LOCK_UN)
                    break
            except IOError:
                time.sleep(0.05)
        else:
            raise BeersException(
                f"Unable to access audit file, {audit_file_path}, for writing within a reasonable period.")
