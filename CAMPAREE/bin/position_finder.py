#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
import re
from prettytable import PrettyTable

def main():
    parser = argparse.ArgumentParser(description='Genome Position Finder')
    parser.add_argument('-p', '--path', required=True, help='Full path to mapping file.')
    args = parser.parse_args()

    position_pattern = re.compile(r"(.+):(\d+)-(\d+)")

    # Separating on tabs and dashes (so start and end are in separate cols)
    df = pd.read_csv(args.path, sep="\t|-", header=0, names=['chr','rstart', 'rend', 'gstart', 'gend'], engine='python')

    # Change stars to -1 as they represent the indel 'loser'
    df=df.replace({'*':-1})

    # Make start/end value ints
    df.rstart = df.rstart.astype(int)
    df.rend = df.rend.astype(int)
    df.gstart = df.gstart.astype(int)
    df.gend = df.gend.astype(int)

    # Display data as it looks now
    print("The data looks like this:")
    print(df.head())

    # Continuous loop - user breaks out with a CR
    while(True):

        # Note!  No validation
        position = input("Indicate positions as chr:start-end or hit CR to escape: ")
        if position == '':
            break
        position_match = re.match(position_pattern, position.strip())
        if position_match:
            chr = position_match.group(1)
            start = int(position_match.group(2))
            end = int(position_match.group(3))

            # Isolate the dataset to the chr of interest.
            chr_df = df[df['chr']==chr]
            if chr_df.empty:
                print("That chr is not found")
                continue

            # Print output in tabular format (we do the mapping both ways)
            output_table = PrettyTable()
            output_table.field_names = ['Starting from', 'Given span', 'Mapped span',
                                        "# Bases Inserted/Deleted"]
            output_table.align['Starting from'] = 'l'
            output_table.align['Given span'] = 'l'
            output_table.align['Mapped span'] = 'l'
            output_table.align['# Bases Inserted/Deleted'] = 'r'

            # Map reference to custom genome
            os, oe = mapping_span(['rstart','rend'], ['gstart', 'gend'], start, end, chr_df)
            output_table.add_row(["Reference", f"{chr}:{start}-{end}", f"{chr}:{os}-{oe}", (oe-os) - (end-start)])

            # Map custon genome to reference
            os, oe = mapping_span(['gstart', 'gend'], ['rstart', 'rend'], start, end, chr_df)
            output_table.add_row(["Genome", f"{chr}:{start}-{end}", f"{chr}:{os}-{oe}", (oe - os) - (end - start)])

            # Print results
            print(output_table)

def mapping_span(i, o, start, end, chr_df):
    """
    Do the work of locating the mapped span start and end values
    :param i:  array of start and end names for input columns in dataset
    :param o: array of start and end names for mapped columns in dataset
    :param start: input span start value
    :param end: input span end value
    :param chr_df: dataset narrowed to requested chromosome
    :return: tuple of the mapped span start and end values
    """
    os = None
    oe = None

    # Narrow the dataset to the contiguous rows containing the start and end values for the input span
    i_df = chr_df[chr_df[i[1]] >= start]
    i_df = i_df[i_df[i[0]] <= end]

    # Get the offsets relative to the input start value for the appropriate row for the input start and end values
    # The 1st row should contain the start value and the last row should contain the end value
    start_offset = start - i_df[i[0]].iloc[0]
    end_offset = end - i_df[i[0]].iloc[-1]

    # Iterate over all the rows in this more limited dataset
    for n in range(len(i_df.index)):

        # Collect the row
        row = i_df.iloc[n]

        # The row is interesting only if it contains mapped start/end values.  Otherwise get the
        # next row
        if row[o[0]] > 0:

            # If no mapped span start value is set and we are in the row containing the input span start value, the
            # mapped span start is just the start offset from the row mapped start value.  Otherwise the mapped span
            # start is the just the row mapped start value.
            if not os:
                os = row[o[0]] + start_offset if n == 0 else row[o[0]]

            # If the row contains no input start/end values or it is not the last row in this dataset, set the mapped
            # span end to the row mapped end value.  This will possibly change with subsequent iterations.
            if row[i[1]] < 0 or n != len(i_df.index) - 1:
                oe = o[1]

            # If we are in the row containing the input span end value, the mappend span end is just the end offset
            # from the row mapped start value.
            elif n == len(i_df.index) - 1:
                oe = row[o[0]] + end_offset

    # Return the mapped start and end values.  Note that both could be None if the input span is completely
    # deleted in the mapped data in the dataset.
    return os, oe



if __name__ == "__main__":
    main()
