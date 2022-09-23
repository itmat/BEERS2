# BEERS2

## Release status

This is publicly funded research, therefore for the sake of full transparency we are developing BEERS2 as a public repository. Please stay tuned until our first beta release, sometime later this summer.

In the meantime please check out CAMPAREE, a member of the BEERS2 suite of tools that simulates the collection of RNA molecules in a biological sample, prior to RNA-seq library preparation.
https://github.com/itmat/CAMPAREE

## Installation

Somewhere under your home directory, clone the develop branch of the BEERS2 repository:

```bash
git clone --recurse-submodules git@github.com:itmat/BEERS2.git
```

Now set up a virtual environment using python 3.10.

```bash
cd BEERS2
python -m venv ./venv
source ./venv/bin/activate
```

You'll know the virtual environment is activated because the virtual environment path with precede
your terminal prompt.  Now you need to add the python packages/modules upon which BEERS2 depends.  You
do that by installing the packages/modules listed in the requirements.txt file like so:

```bash
pip install -r requirements.txt
```

Then install BEERS2 and its submodules, BEERS_UTILS and CAMPAREE.

```bash
pip install .
cd BEERS_UTILS
pip install .
cd ../CAMPAREE
pip install .
cd ..
```

### Snakemake Configuration

BEERS2 uses Snakemake to run its pipeline, which allows distribution across multiple cores or across a cluster environment.
Many computing environments (including LSF, SGE or SLURM) already have configuration available at the [Snakemake-profiles](https://github.com/Snakemake-Profiles/doc) project.
Choose your environment and follow the instructions there, then use the installed Snakemake profile with the `--profile` option when running BEERS2 (see below).

## Running BEERS

### Example Dataset
We will run a simplified 'baby' example, with a reduced mouse genome and provided output to test the installation and Snakemake configuration.
With the activated venv created above and in the BEERS2 directory, run:

```bash
snakemake --configfile config/baby.config.yaml --directory output/ --cores 1
```

Here BEERS2 will send all output to the `output/` directory, which can be set to where-ever one wishes.

The `--cores 1` option sets to run this on a single-core locally.
Increasing this number will allow Snakemake to run multiple processes simultaneously on the machine your execute the `snakemake` command from.
Alternatively, if you installed a Snakemake profile for your computing environment above, you can include `--profile myprofilename` (where `myprofilename` is whatever name you chose in the profile creation process).
This will allow Snakemake to execute BEERS commands on other machines in your cluster and can be configured using `--cores` and `--jobs` options to control the number of cores and jobs to run simultaneously.

If you wish to run BEERS2 in anther location (without `BEERS2/` as the current working directory), you can also supply `--snakefile path/to/BEERS2/Snakefile`.

### Configuring BEERS

The command above used the example config file `config/baby.config.yaml`.
To configure BEERS, setting the input files and parameters for generating the data, please copy the template config from `BEERS2/config/template.config.yaml` and manually edit it.
This config is commented to guide you through selecting or changing values.
Values required to be set are marked with a TODO comment.

Configuration files are recommended to use absolute file paths to input files.
Otherwise, paths are interpretted relative to the set `--directory` output directory or relative to the BEERS2 directory.

## Useful Snakemake commands

In addition to the already mentioned `--profile`, `--cores`, `--jobs` and `--directory` values, the following are useful to run BEERS.
To check that everything is in place, use `--dry-run`, which will report all steps that will be run.
To re-run a portion of the pipeline, use `--force-rerun {rule name or output file}` where `rule name` could be any of the names of rules in `Snakefile`, for example `run_library_prep_packet_from_molecule_file`,
or you can pass any file produced by BEERS2 and the last rule used to create it will be rerun (as well as any downstream rules).

Depending upon settings like the number of molecules in an input packet, BEERS2 may run out of memory on some steps.
This can be mitigated either by splitting into more packets or by increasing memory limits by passing `--set-resources {rule_name}:mem_mb=32000` or similar with the name of the rule that ran out of memory.

## Input

The recommended source for input files to BEERS2 is [CAMPAREE](https://github.com/itmat/CAMPAREE).
Alternatively, an easy way to provide input RNA to sequence is as a molecule file.
These are tab-separated files with the following structure:

```
transcript_id chrom parental_start parental_cigar ref_start ref_cigar strand sequence
```
The `parent_start, parental_cigar` are 1-indexed alignment start position and [cigar string](https://genome.sph.umich.edu/wiki/SAM#What_is_a_CIGAR.3F) relative to the molecule's parent molecule, while the `ref_start` and `ref_cigar` are alignments relative to the reference genome.
CAMPAREE considers the simulated genome to be the parental molecule, which differs from the reference genome by any variants that the simulated genome may have that are not present in the reference.
In other cases, it is safe to set `parental_start` and `parent_cigar` to the same as `ref_start` and `ref_cigar`.
Strand should be either "+" or "-" and the sequence is a string of ACGT values (N is not allowed and should be replaced with an A, C, G, or T prior to running BEERS).


## Output

BEERS2 produces output through the entire process, so that individual aspects of library preparation and sequencing can be investigated.
However, for most uses, the final outputs will be all that is needed.
These are placed in the `results/` subdirectory of whatever output directory was set with `--directory`.
Output types can be set in the config but are either:

 - FASTQ file: these contain the reads much like a real sequencing run output.
 - SAM/BAM: these contain the reads and quality scores as well as a ideal alignment, giving the true simulated origin of the molecule, even if such an origin can no longer be recreated from its sequence.
         This causes some unusual alignments, such as "100S" alignments that skip all bases but still have a position,
         indicating that the read is entirely adapter sequence or PolyA or other things that do not originate in the genome, but that the original source fragment originated at the specified location.

For all types of output, the read ids contain additional useful information.
An example below is annotated:

```
@BEERS:1:1:1103:32736:10931:ENSMUST00000055032_1.74.941.cdna.cdna.7	1:N:0:AGCGCTAG+AACCGCGG
@BEERS:{flowcell id}:{lane id}:{x}:{y}:{input_molecule ID}.{BEERS additional fields, separated by .}  {read number}:N:0:{i5 barcode}+{i7 barcode}
```

Read numbers are either 1 or 2 to denote forward or reverse reads if paired end.
Here the additional fields are added by most steps of the pipeline when producing a new molecule, so that they can be distinguished from other molecules.
For example, the 'cdna' markers are added by the FirstStrandSynthesisStep and SecondStrandSynthesisStep steps.
The final number (7 in the example), is added by PCR amplitifcation, meaning that one can identify all PCR duplicates by looking for ids that differ only in that final number.
If using input from CAMPAREE, the transcript name will have a `_1` or `_2`  at the end, indicating which allele it derived from.

## Funding

Work on the BEERS2 project is supported by R21-LM012763-01A1: “The Next Generation of RNA-Seq Simulators for Benchmarking Analyses” (PI: Gregory R. Grant).
