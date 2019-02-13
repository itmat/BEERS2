from collections import namedtuple
import os

Constants = namedtuple('Constants', ['ROOT_DIR',
                                     'VARIANTS_FILE_NAME',
                                     'BEAGLE_DATA_FILE_NAME',
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

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONSTANTS = Constants(ROOT_DIR, "variants.txt", "beagle.vcf.vcf.gz", "male", "female",
                      "audit.txt", 100, 'data', 'logs', 'stdout', 'stderr', [1, 2])

SUPPORTED_DISPATCHER_MODES = ['multicore', 'lsf', 'serial']
