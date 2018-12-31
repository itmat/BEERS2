from collections import namedtuple

Constants = namedtuple('Constants', ["AUDIT_FILENAME",
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

CONSTANTS = Constants("audit.txt", 100, 'data', 'logs', 'stdout', 'stderr', [1, 2])
