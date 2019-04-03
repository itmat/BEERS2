#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
import re

def main():
    parser = argparse.ArgumentParser(description='Genome Position Finder')
    parser.add_argument('-p', '--path', required=True, help='Full path to mapping file.')
    args = parser.parse_args()

    position_pattern = re.compile(r"(.+):(\d+)-(\d+)")

    df = pd.read_csv(args.path, sep="\t|-", header=0, names=['chr','rstart', 'rend', 'gstart', 'gend'], engine='python')
    df=df.replace({'*':-1})
    df.rstart = df.rstart.astype(int)
    df.rend = df.rend.astype(int)
    df.gstart = df.gstart.astype(int)
    df.gend = df.gend.astype(int)
    print("The data looks like this:")
    print(df.head)g

    while(True):
        position = input("Indicate positions as chr:start-end or type 'end': ")
        if position == 'end':
            break
        position_match = re.match(position_pattern, position.strip())
        if position_match:
            chr = position_match.group(1)
            start = int(position_match.group(2))
            end = int(position_match.group(3))
            chr_df = df[df['chr']==chr]
            if chr_df.empty:
                print("That chr is not found")
                continue
            os, oe = mapping_span(['rstart','rend'], ['gstart', 'gend'], start, end, chr_df)
            print(f"Reference {chr}:{start}-{end} -> {chr}:{os}-{oe} #bases changed {(oe-os) - (end-start)}")
            os, oe = mapping_span(['gstart', 'gend'], ['rstart', 'rend'], start, end, chr_df)
            print(f"Genome {chr}:{start}-{end} -> {chr}:{os}-{oe} #bases changed {(oe-os) - (end-start)}")


def mapping_span(i,o, start, end, chr_df):
        os = None
        oe = None
        i_df = chr_df[chr_df[i[1]] >= start]
        i_df = i_df[i_df[i[0]] <= end]
        soffset = start - i_df[i[0]].iloc[0]
        eoffset = end - i_df[i[0]].iloc[-1]
        for n in range(len(i_df.index)):
            row = i_df.iloc[n]
            if row[o[0]] > 0:
                if not os:
                    os = row[o[0]] + soffset
                if row[i[1]] < end:
                    oe = o[1]
                else:
                    oe = row[o[0]] + eoffset
        return os,oe



if __name__ == "__main__":
    main()
