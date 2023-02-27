# BEERS2

## Release status

This is publicly funded research, therefore for the sake of full transparency we are developing BEERS2 as a public repository. Please stay tuned until our first beta release, sometime later this summer.

In the meantime please check out CAMPAREE, a member of the BEERS2 suite of tools that simulates the collection of RNA molecules in a biological sample, prior to RNA-seq library preparation.
https://github.com/itmat/CAMPAREE

## Installation

First set up a virtual environment using python 3.10.

```bash
python -m venv ./venv_beers
source ./venv_beers/bin/activate
```

You'll know the virtual environment is activated because the virtual environment path with precede
your terminal prompt.
Now we need to install the BEERS2 package.
Note that we do not recommend doing this outside of the newly created virtual environment.

```bash
pip install git+https://github.com/itmat/BEERS2
```

## Running BEERS

### Example Dataset
We will run a simplified 'baby' example, with a reduced mouse genome and provided output to test the installation and Snakemake configuration.
We first download the example data and config file.

```bash
wget -c https://s3.amazonaws.com/itmat.data/BEERS2/examples/baby_mouse_example.tar.gz -O - | tar -xz
cd baby_mouse_example/
```

Now we are ready to run BEERS2 using this dataset.

```bash
run_beers --configfile baby.config.yaml --jobs 1
```

Verify that the run has been successful by examing `results/` which should contain output FASTQ and SAM files.
To confirm that the expected results were produced, compare md5 hashes to this reference:

```bash
$ md5sum results/*sam
c02704029f6e16a77ea6890ebf92375a  results/S1_L1.sam
35dca3ee057e7a8bac8b8a5d790b054f  results/S1_L2.sam
7c963445bab5153c63c87d507ba589aa  results/S1_unidentifiedL1.sam
00dcaa96464e2e207aa2989650d0b8f3  results/S1_unidentifiedL2.sam
e1b0b04a7e3ca9274b5fa59b03479f97  results/S2_L1.sam
4a7e526ea6b333b19ff40689f659e680  results/S2_L2.sam
19937368990e8bff3ed6b86a8a6718fc  results/S2_unidentifiedL1.sam
bef563a2384e92961bd9701e5f5e2720  results/S2_unidentifiedL2.sam
```

The `--jobs 1` option sets to run this on a single-core locally.
Increasing this number will allow Snakemake to run multiple processes simultaneously on the machine your execute the `run_beers` command from.

### Configuring BEERS

The command above used the example config file `baby.config.yaml`.
To configure BEERS, setting the input files and parameters for generating the data, please copy the template config from [config/template.config.yaml](config/template.config.yaml) and manually edit it.
This config is commented to guide you through selecting or changing values.
Values required to be set are marked with a TODO comment.

Configuration files are recommended to use absolute file paths to input files.
Otherwise, paths are interpreted relative to the set `--directory` output directory or relative to the BEERS2 directory.

### Cluster Execution

BEERS2 uses Snakemake to run its pipeline, which allows distribution across multiple cores or across a cluster environment.
Many computing environments (including LSF, SGE or SLURM) already have configuration available at the [Snakemake-profiles](https://github.com/Snakemake-Profiles/doc) project.
Choose your environment and follow the instructions there, then use the installed Snakemake profile with the `--profile` option when running BEERS2 (see below).
After configuring the Snakemake profile, our example run command could change to:

```bash
run_beers --configfile baby.config.yaml --jobs 4 --profile MYPROFILENAME
```
Where `MYPROFILENAME` should be replaced with the name of the profile you chose during profile configuration. This will tell Snakemake to execute BEERS2 commands on other machines in your cluster.
The `--jobs` options should be changed according to the number of jobs to use simultaneously; at most one per molecule packet provided as input will be used (4 for the example we ran).
BEERS2 will use up to one job per packet at a time, so increasing `jobs` above that number will not cause further parallelism.

### Useful Snakemake commands

The `run_beers` command is a thin wrapper around `snakemake` command, pointing it to the correct `Snakefile` to run, so any commands passed to `run_beers` will be given to `snakemake`.
In addition to the already mentioned `--profile`, `--jobs` and `--directory` values, the following are useful to run BEERS.
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
The final number (7 in the example), is added by PCR amplification, meaning that one can identify all PCR duplicates by looking for ids that differ only in that final number.
If using input from CAMPAREE, the transcript name will have a `_1` or `_2`  at the end, indicating which allele it derived from.

### Quantifications

In addition to providing the alignments, BEERS2 generates quantification counts of both the input molecules and the output molecules as `.quant_file` files.
The format for each of these is a tab-separated file with two columns, labelled `transcript_id` and `counts`.
The `transcript_id` is whatever ID was provided to the BEERS in the input molecule packet or distribution (see above).

The output sequenced fragment counts are available in the `library_prep_pipeline/sample*/from_distribution/` and `library_prep_pipeline/sample*/from_molecule_file/` directories and the input quant files are available in the `input_quants` subdirectories of those directories.
One quant file is produced for each molecule packet.
To obtain quantifications for a sample, sum all quant files from all packets (from both `from_distribution` and `from_molecule_file`).

Input quant files quantify the transcript counts that were in the original sequence and are therefore useful for computing TPM values.
The output quant files count the number of fragments sequenced and therefore are useful for computing FPKM values.
Since one input transcript may generate zero, one, or many sequenced fragments, the input and output quantifications can differ substantially depending upon the transcript sequence and the configuration of BEERS.
The output of the library prep pipeline exactly corresponds one-to-one with the sequenced fragments, so there are no additional quantification files provided for the sequencing pipeline (just use the library prep output quants).

### Logs

Logs store partially-completed information about the BEERS2 run and may be useful for diagnosing where in the pipeline various sequencing artefacts arise.
For example, the `library_prep_pipeline/sample1/logs/FragmentStep` directory contains one log per molecule packet that was processed by the Fragment Step (which simulates fragmentation of transcripts).
The log files contain the output molecules from the step, sometimes with additional 'notes' on the molecule about what happened in the step, in a tab-separated file.

## Development

To contribute to or modify BEERS2, we recommend the following setup. First, clone the git repository including the submodules:

```bash
git clone --recurse-submodules https://github.com/itmat/BEERS2.git
```

Then create and activate a virtual environment for BEERS2 development using Python 3.10:

```bash
cd BEERS2
python -m venv ./venv_beers
source ./venv_beers/bin/activate
```

Install BEERS for local editing with development requirements:

```bash
pip install -e .[dev]
cd BEERS_UTILS
pip install -e .
cd ../CAMPAREE
pip install -e .
cd ..
```

Now any edits made to BEERS2 will be reflected when running it.
Also, we can now run some useful development commands: `pytest tests` to run the unit tests, `mypy src` to check typing errors, and `make html` in the `doc` directory to build the documentation.

## Funding

Work on the BEERS2 project is supported by R21-LM012763-01A1: “The Next Generation of RNA-Seq Simulators for Benchmarking Analyses” (PI: Gregory R. Grant).
