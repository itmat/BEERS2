import gzip
from beers_utils.molecule import Molecule

class Logger:
    def __init__(self, log_file, compression, full_logs):
        '''
        Logger for BEERS steps.

        log_file:
            path to write to
        compression:
            what to use any compression, default None
        full_logs:
            if True, then write out all molecules
            otherwise, we skip that step
        '''
        self.log_filename = log_file
        self.compression = compression
        self.full_logs = full_logs

    def __enter__(self):
        if self.compression is None:
            self.log_file = open(self.log_filename, "wt+")
        elif self.compression == 'gzip':
            self.log_file = gzip.open(self.log_filename, "wt+", compresslevel=6)
        else:
            raise NotImplementedError(f"Unknown compression {repr(compression)} for log file {self.log_file}")
        self.log_file.write(Molecule.header)
        return self

    def __exit__(self, type, value, traceback):
        self.log_file.close()

    def write(self, molecule, note = ''):
        '''
        Write out an entry to the log 
        
        Skipped if full_logs not enabled
        '''
        if self.full_logs:
            self.log_file.write(molecule.log_entry(note))

