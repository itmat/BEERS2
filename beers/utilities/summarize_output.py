#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="Summarizes BEERS2 output for the purposes of assessing biases, alignment, and apparent expression.")

parser.add_argument("-c", "--config", help="BEERS config file to summarize", required=True)
parser.add_argument("-r", "--run_id", help="ID of the run to summarize", required=True)
parser.add_argument("-o", "--outdir", help="Directory to output the summary results to", required=True)

args = parser.parse_args()



import json
import pathlib
import matplotlib
import pylab
import numpy
import pysam
from beers_utils.read_fasta import read_fasta
from beers_utils.molecule_packet import MoleculePacket
from beers_utils.molecule import Molecule
from beers.utilities.gc_content import packet_gc_content
from beers.cluster_packet import ClusterPacket
from beers_utils import cigar

def make_pileups(reference):
    pileups = {ref: {strand: numpy.zeros(len(seq)) for strand in ["+", "-"]}
                        for ref, seq in reference.items()}
    return pileups
def pileup(molecules, pileups):
    for molecule in molecules:
        chrom = molecule.source_chrom
        start = molecule.source_start
        cigar_str = molecule.source_cigar
        strand = molecule.source_strand

        i = start - 1
        for op, num in cigar.split_cigar(cigar_str):
            if cigar.consumes[op]['ref']:
                if cigar.consumes[op]['query']:
                    pileups[chrom][strand][i:i+num] += 1
                i += num
    return pileups

def plot_pileups(pileups_dict, reference):
    pileup_dir = outdir/"pileups"
    pileup_dir.mkdir(exist_ok=True)

    for chrom in reference:
        length = len(reference[chrom])

        fig, axes = pylab.subplots(
                nrows=len(pileups_dict),
                figsize=(40, 2.5*len(pileups_dict)),
        )
        for (pileup_name, pileup), ax in zip(pileups_dict.items(), axes):
            ax.set_ylabel(pileup_name)
            for strand in ['+', '-']:
                strand_pileup = pileup[chrom][strand]
                x = numpy.arange(length)
                strand_pileup *= (1 if strand == '+' else -1)
                change = numpy.concatenate((numpy.diff(strand_pileup), [0])) != 0
                ax.fill_between(
                    x[change],
                    strand_pileup[change],
                    color = 'r' if strand == '+' else 'b',
                    step = "pre",
                    linewidth = 0,
                )
            ax.set_xlim(0, length)
        fig.savefig(pileup_dir / f"{chrom}.png", dpi=300)

def plot_gcs(gcs_dict):
    fig, ax = pylab.subplots(
            figsize = (5, 5)
    )

    for name, gcs in gcs_dict.items():
        bins = [bin for dist, bin in gcs]
        dists = [dist for dist, bin in gcs]
        dist = numpy.concatenate((numpy.nanmean(dists, axis=0), [0]))
        ax.plot(bins[0], dist, label=name)
    ax.set_ylabel("Percent Reads")
    ax.set_xlabel("Percent GC")
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    ax.legend()
    return fig

outdir = pathlib.Path(args.outdir)
outdir.mkdir(exist_ok=True)

config = json.load(open(args.config))

reference_genome = read_fasta(config['resources']['reference_genome_fasta'])
reference_genome = {ref.split()[0]: seq for ref, seq in reference_genome.items()} # Drop to just the chrom names

input_dir = pathlib.Path(config['library_prep_pipeline']['input']['directory_path'])
input_samples = input_dir.rglob("**/molecule_file.txt")

input_pileups = make_pileups(reference_genome)
input_gcs = []
input_aligned_gcs = []
for sample in input_samples:
    packet = MoleculePacket.from_CAMPAREE_molecule_file(sample, packet_id=0)
    pileup(packet.molecules, input_pileups)
    input_gcs.append(packet_gc_content(packet))
    input_aligned_gcs.append(packet_gc_content(packet, aligned_only=True))

output_dir = pathlib.Path(config['controller']['output_directory_path'] + f"_run{args.run_id}")
lib_prep_out = output_dir/"library_prep_pipeline/data/"

lib_prep_molecule_packets = lib_prep_out.rglob("**/library_prep_pipeline_result_molecule_pkt*.txt")

lib_prep_pileups = make_pileups(reference_genome)
lib_prep_gcs = []
lib_prep_aligned_gcs = []
for packet in lib_prep_molecule_packets:
    packet = MoleculePacket.deserialize(packet)
    pileup(packet.molecules, lib_prep_pileups)
    lib_prep_gcs.append(packet_gc_content(packet))
    lib_prep_aligned_gcs.append(packet_gc_content(packet, aligned_only=True))

seq_out = output_dir/"controller/data/"
sam_files = seq_out.rglob("*.sam")
seq_pileups = make_pileups(reference_genome)
seq_gcs = []
seq_aligned_gcs = []
for sam_file in sam_files:
    molecules = []
    #TODO: handle overlapping reads from the same fragment properly
    with pysam.AlignmentFile(str(sam_file), "r") as sam:
        for read in sam:
            strand = '-' if ((read.is_reverse and read.is_read1) or (not read.is_reverse and read.is_read2)) else '+'
            molecules.append(Molecule(
                start = 1,
                molecule_id = 0,
                sequence = read.query_sequence,
                source_chrom = read.reference_name,
                source_start = read.reference_start + 1,
                source_cigar = read.cigarstring,
                source_strand = strand,
            ))
    pileup(molecules, seq_pileups)
    packet = MoleculePacket(0, None, molecules)
    seq_gcs.append(packet_gc_content(packet))
    seq_aligned_gcs.append(packet_gc_content(packet, aligned_only=True))

plot_pileups(
    {
        "CAMPAREE": input_pileups,
        "Lib Prep": lib_prep_pileups,
        "Seq": seq_pileups
    },
    reference=reference_genome)

gc_dir = outdir/"gcs"
gc_dir.mkdir(exist_ok=True)
fig = plot_gcs(
    {
        "CAMPAREE": input_gcs,
        "Lib Prep": lib_prep_gcs,
        "Seq": seq_gcs
    }
)
fig.savefig(gc_dir / "gc.png", dpi=300)

fig = plot_gcs(
    {
        "CAMPAREE": input_aligned_gcs,
        "Lib Prep": lib_prep_aligned_gcs,
        "Seq": seq_aligned_gcs
    }
)
fig.savefig(gc_dir / "gc.aligned_only.png", dpi=300)
