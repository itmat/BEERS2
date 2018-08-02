import numpy as np
import argparse

def generate(filename, percent_mrna, count, mean_length, std_length, a_tail_mean_length, a_tail_std_dev_length):
    with open(filename, "w+") as data:
        sequence_lengths = np.random.normal(mean_length, std_length, count)
        poly_a_lengths = np.random.normal(a_tail_mean_length, a_tail_std_dev_length, count)
        for index in range(count):
            nucleotides = ['A','T','G','C']
            sequence = np.random.choice(nucleotides, int(sequence_lengths[index]))
            tail = ''
            mrna = np.random.choice([1, 0], 1, p=[float(percent_mrna), 1 - float(percent_mrna)])[0]
            if mrna:
                tail = int(poly_a_lengths[index])*'A'
            data.write("".join(sequence) + tail + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Starter RNA Generator')
    parser.add_argument('-f', '--filename', default='../../data/molecules.txt',
            help = 'Provide the name of the file where the sequences are stored.')
    parser.add_argument('-s', '--seed', default=None, type=int,
            help = 'Provide a random seed - default is none.')
    parser.add_argument('-c', '--count', default=100, type=int,
                        help='Number of molecules to generate')
    parser.add_argument('-m', '--mean_length', default=1000, type=int,
                        help='Mean length (gaussian distribution)')
    parser.add_argument('-d', '--std_dev_length', default=100, type=int,
                        help='Standard deviation length (gaussian distribution)')
    parser.add_argument('-a', '--a_tail_mean_length', default=20, type=int,
                        help='Mean a tail length (gaussian distribution)')
    parser.add_argument('-t', '--a_tail_std_dev_length', default=5, type=int,
                        help='Std dev a tail length (gaussian distribution)')
    parser.add_argument('-p', '--percent_mRNA', default=0.25, type=float,
                        help='Percent of sequences that represent mRNA (gaussian distribution)')
    args = parser.parse_args()
    generate(args.filename, args.percent_mRNA, args.count, args.mean_length, args.std_dev_length, args.a_tail_mean_length, args.a_tail_std_dev_length)
