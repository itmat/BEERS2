import sys
import argparse
import os
import json
import numpy as np
import pickle
import bisect
from collections import namedtuple
from beers.library_prep.molecule import Molecule

Files = namedtuple('Files', ['designator', 'quantities_filename', 'transcriptome_filename'])


class TranscriptMaker:

    def __init__(self,
                 maternal_quantities_filename,
                 maternal_transcriptome_filename,
                 paternal_quantities_filename,
                 paternal_transcriptome_filename,
                 log_filename,
                 seed,
                 transcriptome_size):

        self.parents = []
        self.parents.append(Files("m", maternal_quantities_filename, maternal_transcriptome_filename))
        self.parents.append(Files("p", paternal_quantities_filename, paternal_transcriptome_filename))

        self.log_filename = log_filename

        # TODO - create and log a seed if none provided.
        self.seed = int(seed)
        np.random.seed(self.seed)

        self.transcriptome_size = int(transcriptome_size)
        self.polya_tail = 200 * 'A'

        # Holds molecules created from transcript
        self.molecules = set()

        # Holds the cumulative transcript distribution
        self.transcript_distribution = TranscriptCumulativeDistribution()

        # Hold the probabalistically determined transcript counts
        self.transcript_counts_map = dict()

        # Since the log file will be open for appending, we need to insure that the file doesn't
        # currently exist.
        try:
            os.remove(self.log_filename)
        except OSError:
            pass

    def create_transcript_distribution(self):

        # Iterate over both parents
        for parent in self.parents:
            with open(parent.quantities_filename, 'r') as quantities_file:
                for line in quantities_file:
                    if line.startswith(' '):
                        continue

                    fields = line.rstrip('\n').split('\t')

                    # Transcript id and counts - Tack on parent designator as prefix
                    transcript_id = parent.designator + '_' + fields[0]
                    # Convert string to float first to avoid "invalid literal for int()" error,
                    # since some entries in the dist file may have decimal values.
                    transcript_counts = int(float(fields[1]))

                    # Skip over transcripts that have no counts.
                    if transcript_counts == 0:
                        continue

                    # Accumulate
                    self.transcript_distribution.append(transcript_id, transcript_counts)

    def create_random_transcript_id(self):
        for _ in range(self.transcriptome_size):
            return self.transcript_distribution.locate_id(np.random.uniform())

    def prepare_transcripts(self):

        # Create a map of the transcript distributions
        self.create_transcript_distribution()
        self.transcript_distribution.normalize()

        # Map to stash transcript sequences keyed by transcript id.
        transcript_sequences = dict()
        for parent in self.parents:
            with open(parent.transcriptome_filename, 'r') as transcriptome_file:

                transcriptome = json.load(transcriptome_file)

                for transcript_id in transcriptome:
                    # Make sure the contributions of each parent are distinguishable
                    transcript_sequences[parent.designator + "_" + transcript_id] = transcriptome[transcript_id]

        # Dereference this variable since it can consume at lot of memory
        transcriptome = None

        # Create the molecules
        self.make_molecules(transcript_sequences)

    def make_molecules(self, transcript_sequences):

        # Generate the requested number of molecules randomly.
        for _ in range(self.transcriptome_size):

            # Protect against the possibility that we have no transcript sequence for a transcript id selected.
            for _ in range(100):
                transcript_id = self.create_random_transcript_id()
                if transcript_id in transcript_sequences:
                    break
            else:
                print("Cannot find a transcript id for which a sequence is available in 100 attempts.")
                sys.exit(1)

            # Pick an immature transcript sequence 3% of the time.
            transcript_sequence = np.random.choice(transcript_sequences[transcript_id], p = [0.97, 0.03])

            # Create an initial cigar sequence (all matches)
            cigar = str(len(transcript_sequence)) + 'M'

            # Affix a poly A tail
            transcript_sequence += self.polya_tail

            # Generate as many molecules of the transcript as dictated by the transcript's quantity.
            self.molecules.add(Molecule(Molecule.next_molecule_id, transcript_sequence, 1,
                               cigar, transcript_id))
            Molecule.next_molecule_id += 1

    def serialize_molecules(self):
        #Output pickled file
        with open("molecules.pickle", 'wb') as molecules_file:
            pickle.dump(self.molecules, molecules_file)
        #Output human-readable file
        with open("molecules.txt", 'wt') as molecules_file:
            for molecule in self.molecules:
                molecules_file.write(f'{molecule}\n')

    @staticmethod
    def main():
        parser = argparse.ArgumentParser(description='Transcript_Maker')
        # TODO accept output molecules filename or an output directory as an argument
        parser.add_argument('-q', '--maternal_quant_filename')
        parser.add_argument('-t', '--maternal_transcriptome_filename')
        parser.add_argument('-p', '--paternal_quant_filename')
        parser.add_argument('-f', '--paternal_transcriptome_filename')
        parser.add_argument('-l', '--log_filename')
        parser.add_argument('-d', '--seed')
        parser.add_argument('-s', '--transcriptome_size', type=int)
        args = parser.parse_args()
        print(args)
        transcript_maker = TranscriptMaker(args.maternal_quant_filename,
                                           args.maternal_transcriptome_filename,
                                           args.paternal_quant_filename,
                                           args.paternal_transcriptome_filename,
                                           args.log_filename,
                                           args.seed,
                                           args.transcriptome_size)
        transcript_maker.prepare_transcripts()
        transcript_maker.serialize_molecules()
        print(f'Number of molecules: {len(transcript_maker.molecules)}')


class TranscriptCumulativeDistribution:

    def __init__(self):
        self.transcript_ids = []
        self.transcript_cumulative_counts = []
        self.count_accumulator = 0

    def append(self, transcript_id, count):
        self.transcript_ids.append(transcript_id)
        self.count_accumulator += count
        self.transcript_cumulative_counts.append(self.count_accumulator)

    def normalize(self):
        self.transcript_cumulative_counts = \
            [count/self.count_accumulator for count in self.transcript_cumulative_counts]

    def locate_id(self, value):
        return self.transcript_ids[bisect.bisect(self.transcript_cumulative_counts, value)]


if __name__ == "__main__":
    sys.exit(TranscriptMaker.main())

'''
Example command
python transcript_maker.py \
    -q ../../data/expression/GRCh38/maternal_transcript_dist.txt \
    -t ../../data/expression/GRCh38/maternal.GRCh38.ensembl92.transcriptseq.json \
    -p ../../data/expression/GRCh38/maternal_transcript_dist.txt \
    -f ../../data/expression/GRCh38/maternal.GRCh38.ensembl92.transcriptseq.json \
    -l ../../data/preBEERS/mm10/transcript_maker.log \
    -d 100 \
    -s 10000
    '''
