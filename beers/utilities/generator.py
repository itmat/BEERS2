import numpy as np
import argparse
import functools
import sys
import pandas as pd
from plotly import graph_objs as go
import plotly.offline as pyo

def normal_model_generator(args):
    """
    Creates a file containing normally distributed lenghts of random sequences with a mean length and a
    standard deviation length as provided in the argument list.  Additionally, a poly A tail having a
    mean length and a standard deviation length as provided in the argument list, is tacked onto the 3' end
    of the sequence.
    :param args:
    """
    with open(args.filename, "w+") as data:
        sequence_lengths = np.random.normal(args.mean_length, args.std_length, args.count)
        poly_a_lengths = 0
        if args.a_tail_mean_length > 0 and args.a_tail_std_dev_length > 0:
            poly_a_lengths = np.random.normal(args.a_tail_mean_length, args.a_tail_std_dev_length, args.count)
        for index in range(args.count):
            nucleotides = ['A','T','G','C']
            sequence = np.random.choice(nucleotides, int(sequence_lengths[index]))
            tail = ''
            mrna = np.random.choice([1, 0], 1, p=[float(args.percent_mrna), 1 - float(args.percent_mrna)])[0]
            if mrna:
                tail = int(poly_a_lengths[index])*'A'
            data.write("".join(sequence) + tail + "\n")
    show_histogram(args.filename, "Normal Distribution", args.count)

def uniform_model_generator(args):
    """
    Creates a file containing uniformly distributed lengths of random sequences varying from the
    minimum to the maximum lengths provided in the argument list.  A histogram is produced showing the
    results of the sampling.
    :param args:
    """
    with open(args.filename, "w+") as data:
        sequence_lengths = np.random.uniform(args.max_length, args.min_length, args.count)
        for index in range(args.count):
            nucleotides = ['A', 'T', 'G', 'C']
            sequence = np.random.choice(nucleotides, int(sequence_lengths[index]))
            data.write("".join(sequence) + "\n")
    show_histogram(args.filename, "Uniform Distribution", args.count)

def show_histogram(filename):
    """
    Displays a histogram of the script results, to serve as a sanity check.
    :param filename:
    """
    df = pd.read_csv(filename, names = ["sequence"])
    df['seq_length'] = df.apply(lambda row: len(row['sequence']), axis=1)
    data = [
        go.Histogram(
            x=df["seq_length"],
            opacity=0.75,
            name="Molecules"
        )
    ]
    layout = go.Layout(title="Generated Sequence Lengths Histogram")
    fig = go.Figure(data=data, layout=layout)
    pyo.plot(fig, filename='../../data/generated_sequence.html')

def main():
    """
    Utility script for creating either normal or uniform RNA distributions.  Parameters differ depending upon
    the type of sample requested.
    """
    parser = argparse.ArgumentParser(description='Starter RNA Generator')
    subparsers = parser.add_subparsers(help='sequence length probability model')
    parser_normal_model = subparsers.add_parser('normal_model', help='sequence lengths normally distributed')
    parser_uniform_model = subparsers.add_parser('uniform_model', help='sequence lengths uniformly distributed')
    parser.add_argument('-f', '--filename', default='../../data/molecules.txt',
                        help = 'Provide the name of the file where the sequences are stored.')
    parser.add_argument('-s', '--seed', default=None, type=int,
                        help = 'Provide a random seed - default is none.')
    parser.add_argument('-c', '--count', default=10000, type=int,
                        help='Number of molecules to generate - default is 10000.')
    parser_normal_model.add_argument('-m', '--mean_length', default=1000, type=int,
                                     help='Mean length (gaussian distribution)')
    parser_normal_model.add_argument('-d', '--std_dev_length', default=100, type=int,
                                     help='Standard deviation length (gaussian distribution)')
    parser_normal_model.add_argument('-a', '--a_tail_mean_length', default=20, type=int,
                                     help='Mean a tail length (gaussian distribution)')
    parser_normal_model.add_argument('-t', '--a_tail_std_dev_length', default=5, type=int,
                                     help='Std dev a tail length (gaussian distribution)')
    parser_normal_model.add_argument('-p', '--percent_mRNA', default=0.25, type=float,
                                     help='Percent of sequences that represent mRNA (gaussian distribution)')
    parser_normal_model.set_defaults(func=normal_model_generator)
    parser_uniform_model.add_argument('-x', '--max_length', type=int,
                                      help='Maximum length of a sequence')
    parser_uniform_model.add_argument('-n', '--min_length', type=int,
                                      help='Minimum length of a sequece')
    parser_uniform_model.set_defaults(func=uniform_model_generator)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    sys.exit(main())


# To create 50000 rna strands having ~ 5000 bases and with 10% having ~ 200 polyA tails
# python generator.py -f "../../data/molecules5k50k.txt" -c 50000 normal_model -m 5000 -d 500 -a 200 -t 25 -p 0.10

# To create 500000 rna stands having ~5000 bases and with 10% having ~ 200 polyA tails
# python generator.py -f "../../data/molecules5k500k.txt" -c 500000 normal_model  -m 5000 -d 500 -a 200 -t 25 -p 0.10

# To create 50000 rna strands having uniform distribution from 25 - 3000 bases
# python generator.py  -f ../../data/uni_mol50k5k.txt -c 50000 uniform_model -x 3000 -n 25

# To create 100000 rna strands having uniform distribution from 25 - 1000 bases
# python generator.py -f ../../data/uni_mol100k1k.txt -c 100000 uniform_model -x 1000 -n 25
