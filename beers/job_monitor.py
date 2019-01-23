import os
import sys
import time
from beers.constants import SUPPORTED_DISPATCHER_MODES
from beers.beers_exception import BeersException

class Monitor:
    """
    The class monitors the status of various subprocesses running throughout the
    pipeline. It checks for jobs that are pending, running, stalled, or halted
    (either due to success or error/failure).
    """

    def __init__(self, output_directory_path, dispatcher_mode):

        self.process_list = []
        self.output_directory = output_directory_path

        if dispatcher_mode not in SUPPORTED_DISPATCHER_MODES:
            raise BeersException(f'{dispatcher_mode} is not a supported mode.\n'
                                 'Please select one of {",".join(SUPPORTED_DISPATCHER_MODES)}.\n')
        else:
            self.dispatcher_mode = dispatcher_mode

    def is_processing_complete(self, process_list):
        pass
