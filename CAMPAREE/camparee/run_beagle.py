#!/usr/bin/env python
# Taking a list of variant files from BEERS2, creates a VCF file
# and runs Beagle on it to do phase calling


import gzip
import argparse
import contextlib
import collections
import random
import subprocess

BEAGLE_JAR = "beagle.28Sep18.793.jar"

#TODO: move this method to the camparee utilities (it's also defined and used in
#the variants_compilation.py script).
def parse_line(line):
    ''' reads a line of a variant file from BEERS2'''
    # sample line is: (note tabs and spaces both used)
    #1:28494 | C:1 | T:1    TOT=2   0.5,0.5 E=1.0
    if line == '':
        return "DONE", 0,{}
    entries = line.split('\t')
    loc_and_vars, total, fractions, entropy = entries
    loc, *variants = loc_and_vars.split(" | ")
    chromosome,position = loc.split(":")
    position = int(position)
    variants = {base:int(count) for base,count in [variant.split(":") for variant in variants]}
    return chromosome, position, variants

def read_maybe_gzipped(filename):
    if ".gz" in filename:
        return gzip.open(filename, mode="rt", encoding="ascii")
    else:
        return open(filename, mode="r")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert variant files from BEERS 2 to a VCF and run Beagle on it to do phasing, resulting in a phased VCF file")
    parser.add_argument("variant_files", help="list of variant files from BEERS2", nargs="+")
    parser.add_argument("--vcf", help="vcf file to create and then pass to Beagle", required=True)
    parser.add_argument("--out", help="file name to output Beagle results", required=True)
    parser.add_argument("--ref", help="reference genome file to use", required=True)

    args = parser.parse_args()

    # Read the reference file in
    with read_maybe_gzipped(args.ref) as reference_file:
        reference = {}
        contig_order = []
        working_name = None
        working_sequence_list = []
        for line in reference_file:
            if line.startswith(">"):
                if working_name is not None:
                    # Have to add to the reference the sequence we've constructed to start a new one
                    assert len(working_sequence_list) > 0 # Need some sequence
                    assert working_name not in reference
                    sequence = ''.join(working_sequence_list)
                    reference[working_name] = sequence
                    contig_order.append(working_name)

                working_name = line.split()[0][1:]
                working_sequence_list = []
            else:
                working_sequence_list.append(line.strip()) # Don't include newline character in sequences
    print(f"Read reference file with {len(reference)} chromosome")

    print("Converting variants into vcf file")
    contigs_so_far = []
    last_chromosome = None
    # Open then process all the files
    with contextlib.ExitStack() as stack, open(args.vcf, "w") as out:
        # Output the header
        out.write("##fileformat=VCFv4.0\n")
        sample_names = '\t'.join(args.variant_files) # TODO find better sample names!
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_names}\n")

        # Open variant files
        variant_files = [stack.enter_context(open(file_name)) for file_name in args.variant_files]

        # read the first line of each variant file
        next_lines = [vf.readline() for  vf in variant_files]

        # proceed until we have exhausted (and removed from the list) all the variant files
        i = 0
        while any(line != '' for line in next_lines):
            i += 1
            parsed_lines = [parse_line(line) for line in next_lines]
            chromosome = min(chr for chr, pos, vars in parsed_lines if chr != "DONE")
            position = min(pos for chr, pos, vars in parsed_lines if chr == chromosome)
            if chromosome != last_chromosome:
                print(f"Wrote chromosome {last_chromosome} to vcf file")
                if last_chromosome is not None:
                    contigs_so_far.append(last_chromosome)
                    if contig_order.index(last_chromosome) > contig_order.index(chromosome):
                        print(f"WARNING: chromosomes {last_chromosome}, {chromosome} are not in the same order as in the reference genome")
                    if sorted(contigs_so_far) != contigs_so_far:
                        print("ERROR: Need chromosomes to be in alphabetical order (i.e. 1, 10, 11, 12, .., 19, 2, 20, 21, 22, 3, .., MT, X,Y")
                        #TODO: probably will need to use an alternative means of sorting them since alphabetical is terrible but what we currently use in the test files
                last_chromosome = chromosome


            ref_base = reference[chromosome][position-1]

            # List of variant dictionaries (variant -> count) with empty dicts for files that have no entry for this loc
            variants = [vars if (chr == chromosome and pos == position) else {}
                                for chr, pos, vars in parsed_lines]

            # which samples we actually used a variant of at this location
            used_samples = [True if (chr == chromosome and pos == position) else False
                                for chr, pos, vars in parsed_lines]

            # How many variants we had in each sample
            num_variants = [len(vars) for vars in variants]

            # Filter variants by those that are NOT just the reference
            # and do not include "N"s
            variants = [{var:count for var, count in vars.items() if (var != ref_base and "N" not in var)}
                            for vars in variants]

            # Find each sample's most likely variants that aren't just the reference
            def common_variant(variants):
                if variants:
                    maximum_count = max(count for var,count in variants.items())
                    common_variants = [var for var, count in variants.items() if count == maximum_count]
                    # Take random choice of the tied maximally seen variants
                    return random.choice(common_variants)
                else:
                    return None
            most_common_variants = [common_variant(vars) for vars in variants]

            # Output nothing if we've now discarded all the variants
            if all(var == None for var in most_common_variants):
                next_lines = [line if not used else sample_file.readline()
                                for line, used, sample_file in zip(next_lines, used_samples, variant_files)]
                continue

            # Find the reference length of the longest variant (i.e. length of the longest indel, or 1 if only SNPs)
            def variant_length(variant):
                if variant.startswith("D"):
                    return int(variant[1:]) + 1 # +1 includes the reference base before our deletion
                else:
                    return 1
            length = max(variant_length(variant) for variant in most_common_variants if variant is not None)

            # Do we have any deletions to handle?
            #any_dels = any(var.startswith("D") for var in most_common_variants)
            any_dels = (length > 1)

            # accommodate the base before the deletion
            if any_dels:
                position -= 1

            # Give reference of the appropriate size to accommodate the longest variant here
            ref = reference[chromosome][position-1:position-1+length]

            sample_descriptions = []
            alts = []
            for variant, num_vars in zip(most_common_variants, num_variants):
                if variant is None:
                    sample_descriptions.append("0/0") # ref-ref if no variants for this sample
                    continue

                # With deletions, include one extra base
                if any_dels:
                    if variant in "ACGT": # SNP
                        alt = ref[0] + variant + ref[2:]
                    elif variant.startswith("I"):
                        alt = ref[0:2] + variant[1:] + ref[2:]
                    elif variant.startswith("D"):
                        l = int(variant[1:])
                        alt = ref[0] + ref[1+l:]
                    else:
                        raise ValueError(f"unknown variant {variant}")
                else:
                    if variant in "ACGT": # SNP
                        alt = variant
                    elif variant.startswith("I"):
                        alt = ref[0] + variant[1:]
                    elif variant.startswith("D"):
                        raise ValueError(f"shouldn't have a deletion here in variant {variant}")
                    else:
                        raise ValueError(f"unknown variant {variant}")

                if alt not in alts:
                    alts.append(alt)
                idx = alts.index(alt) + 1 # 0 is the ref, 1 is the first alt, etc.

                # We assume that a lone variant is going to be alt-alt
                # but if there are multiple variants, we always use ref as one
                # and the most common alt as the other
                # TODO: should we ever use alt-alt with two different alts? we don't currently
                if num_vars == 1:
                    sample_descriptions.append(f"{idx}/{idx}") # alt-alt
                else:
                    sample_descriptions.append(f"0/{idx}") # ref-alt


            line = "\t".join([chromosome,
                            str(position),
                            ".",
                            ref,
                            ",".join(alts),
                            ".",
                            ".",
                            ".",
                            "GT",
                            "\t".join(sample_descriptions)])
            out.write(line+"\n")

            next_lines = [line if not used else sample_file.readline()
                            for line, used, sample_file in zip(next_lines, used_samples, variant_files)]

    print(f"Finished creating VCF file for Beagle with {i} variant entries")

    # Next step: Run Beagle
    command = f"java -jar {BEAGLE_JAR} gt={args.vcf} out={args.out}"
    print(f"Calling beagle with command:")
    print(command)

    result =  subprocess.run(command, shell=True, check=True)
    print(f"Finished running Beagle.\nAll Done.")
