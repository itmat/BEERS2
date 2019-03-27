import math
import random

mature_rRNA_classes = ["5.8s", "18s", "28s", "5s"]
premature_rRNA_classes = ["45s pre-ribosomal RNA"]
rRNA_classes = mature_rRNA_classes + premature_rRNA_classes

class RibosomalRNA:
    '''
    Load rRNA data and allow creation of rRNA samples.
    '''
    def __init__(self, logfile, parameters):
        self.fasta = self.read_fasta(parameters["rRNA_fasta_file"])

        for rRNA in rRNA_classes:
            assert rRNA in self.fasta, f"rRNA FASTA file {parameters['rRNA_fasta_file']} must contain a '{rRNA}' entry"


        # The percentage of the rRNA that is mature (versus just the pre-ribosomal RNA
        self.percent_mature = parameters["percent_mature"]

    def generate_rRNA_sample(self, num_molecules):
        ''' returns a sample of size `num_molcules` of sequences of
        un-degraded rRNA molecules '''

        num_mature = math.floor(self.percent_mature * num_molecules)
        num_premature = num_molecules - num_mature

        sample = []
        for i in range(num_mature):
            rRNA_class = mature_rRNA_classes[random.randrange(len(mature_rRNA_classes))]
            sample.append( self.fasta[rRNA_class] )

        for j in range(num_premature):
            rRNA_class = premature_rRNA_classes[random.randrange(len(premature_rRNA_classes))]
            sample.append( self.fasta[rRNA_class] )

        return sample

    def read_fasta(self, filename):
        ''' load a fasta file as a dictionary of description:sequence pairs '''
        fasta_file = open(filename, "r")

        entries = {}
        current_entry_description = None
        current_entry = None
        for line in fasta_file:
            if line.startswith(">"):
                if current_entry is not None:
                    entries[current_entry_description] = ''.join(current_entry)

                current_entry_description = line[1:-1] #remove '>' and newline
                current_entry = []
            else:
                if current_entry is None:
                    raise Exception("FASTA file has no description line")
                
                current_entry.append(line[:-1]) #remove newline

        if current_entry is not None:
            entries[current_entry_description] = ''.join(current_entry)

        return entries

if __name__ == '__main__':
    rRNA = RibosomalRNA(logfile = None,
                        parameters = dict(percent_mature = 0.9,
                                      rRNA_fasta_file = "../../data/mm9/rRNA.fa"))

    rRNA_sample = rRNA.generate_rRNA_sample(100)
    print(f"Generated {len(rRNA_sample)} rRNA molecules")
    print(f"The first ten are: {rRNA_sample[:10]}")
