from collections import namedtuple

Constants = namedtuple('Constants', ["AUDIT_FILENAME",
                                     'FILES_PER_DIRECTORY_LIMIT',
                                     'DATA_DIRECTORY_NAME',
                                     'LOG_DIRECTORY_NAME'])

CONSTANTS = Constants("audit.txt", 100, 'data', 'logs')
