from collections import namedtuple
import os

Constants = namedtuple('Constants', ['ROOT_DIR',
                                     'MALE_GENDER',
                                     'FEMALE_GENDER',
                                     'AUDIT_FILENAME',
                                     'FILES_PER_DIRECTORY_LIMIT',
                                     'DATA_DIRECTORY_NAME',
                                     'LOG_DIRECTORY_NAME',
                                     'STDOUT_SUBDIRECTORY_NAME',
                                     'STDERR_SUBDIRECTORY_NAME',
                                     'DIRECTION_CONVENTION'])
Constants.__doc__ = """
Provides a list of constants for the BEERS ecosystem.  Constants for steps are not included here since those are
subject to change by the user.
"""

_ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONSTANTS = Constants(ROOT_DIR=_ROOT_DIR,
                      MALE_GENDER='male',
                      FEMALE_GENDER='female',
                      AUDIT_FILENAME='audit.txt',
                      FILES_PER_DIRECTORY_LIMIT=100,
                      DATA_DIRECTORY_NAME='data',
                      LOG_DIRECTORY_NAME='logs',
                      STDOUT_SUBDIRECTORY_NAME='stdout',
                      STDERR_SUBDIRECTORY_NAME='stderr',
                      DIRECTION_CONVENTION=[1, 2])

SUPPORTED_SCHEDULER_MODES = ['lsf', 'serial']

# Maximum size for a job seed, generated from the controller's seed
# which is specified in the config.json file
MAX_SEED = 2_000_000_000
