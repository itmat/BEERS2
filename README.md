# BEERS2

**BEERS2** is an RNA-seq simulator that combines a flexible and [highly configurable](#development) design with detailed simulation of the entire library preparation and sequencing pipeline.

## Installation

BEERS2 can be installed locally, on an HPC cluster, or in the AWS cloud.

### Local and cluster installation

Python version 3.11 is required to install and run the latest version of BEERS2 locally or on a cluster. It is strongly recommended to install BEERS2 inside of a virtual environment:

```bash
python -m venv venv_beers
source venv_beers/bin/activate
```

Make sure the virtual environment is activated – the virtual environment name has to precede your terminal prompt. Then install BEERS2:

```bash
pip install git+https://github.com/itmat/BEERS2
```

If you wish to install a specific version of BEERS2, add `@` followed by the version name to the previous line.

### Cloud installation

BEERS2 software and the cloud infrastructure on which it runs are packaged together in a [stack](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/stacks.html). Deploying this stack in your AWS account is an easy, wizard-driven process. The stack can also be easily [removed](#removal) from your account when it is not needed anymore. To deploy the stack with default configuration please click the following button:

[![Launch Stack](https://cdn.rawgit.com/buildkite/cloudformation-launch-stack-button-svg/master/launch-stack.svg)](https://console.aws.amazon.com/cloudformation#/stacks/quickcreate?stackName=Beers&templateURL=https://s3.amazonaws.com/bioinf.itmat.upenn.edu/BEERS2/stack/templates/default.yaml)

After providing your AWS credentials, you will be redirected to the stack deployment page. Here you can specify the **stack name** or change the deployment region (using the dropdown button in the upper right corner). Scroll to the bottom of the page, acknowledge the creation of IAM resources, and click the create stack button.

The deployment process takes about 5 minutes. The event table will be refreshed periodically every minute so you can monitor the progress. The first event will have the stack name as the logical ID and <span style="color:#44b9d6">CREATE_IN_PROGRESS</span> as status. The deployment will be complete when the event whose logical ID is the stack name, having status <span style="color:#6aaf35">CREATE_COMPLETE</span>, appears.

It is assumed that you have sufficient privileges to deploy the stack.

#### Configuring stack

During the deployment of the stack a new S3 bucket will be created. If you prefer to use an existing bucket for storing the input and output of BEERS2 pipeline, please deploy the stack by clicking [here](https://s3.amazonaws.com/bioinf.itmat.upenn.edu/BEERS2/stack/templates/bucket.yaml) instead of the button above, and provide the name of your bucket. In this case make sure you deploy the stack in the same region in which your bucket is located, otherwise additional data transfer charges may apply.

The stack can be further configured and tuned using [CDK](https://aws.amazon.com/cdk/) – you can find the source code defining the stack [here](https://github.com/itmat/cask).

#### Pricing

The stack incurs practically no charges when not in use. Charges are incurred for storage of data in the bucket and for the computation resources used during the pipeline execution.

## Running the simulation pipeline

### Local execution

We will run a simplified “baby” example, with a reduced mouse genome to test the installation.
We first download the example input data and the configuration file:

```bash
curl https://s3.amazonaws.com/itmat.data/BEERS2/examples/baby_mouse_example.tar.gz | tar -xz
```

Now we are ready to run BEERS2 using this dataset:

```bash
cd baby_mouse_example
run_beers --configfile baby.config.yaml --jobs 1
```

Verify that the run has been successful by examining `results` directory which should contain output FASTQ and SAM files.
To confirm that the expected results were produced, compare md5 hashes to this reference:

```bash
md5sum results/*sam

c02704029f6e16a77ea6890ebf92375a  results/S1_L1.sam
35dca3ee057e7a8bac8b8a5d790b054f  results/S1_L2.sam
7c963445bab5153c63c87d507ba589aa  results/S1_unidentifiedL1.sam
00dcaa96464e2e207aa2989650d0b8f3  results/S1_unidentifiedL2.sam
e1b0b04a7e3ca9274b5fa59b03479f97  results/S2_L1.sam
4a7e526ea6b333b19ff40689f659e680  results/S2_L2.sam
19937368990e8bff3ed6b86a8a6718fc  results/S2_unidentifiedL1.sam
bef563a2384e92961bd9701e5f5e2720  results/S2_unidentifiedL2.sam
```

### Cluster execution

BEERS2 uses Snakemake to run its pipeline, which facilitates distribution of the workload across multiple CPU cores or cluster nodes. To run BEERS2 in your cluster environment, please go [here](https://github.com/Snakemake-Profiles/doc) to locate the Snakemake profile for your environment (Slurm, LSF, SGE..) and then install the corresponding profile by following the provided installation instructions. Then you can use the installed profile by passing it via `profile` argument:

```bash
run_beers --configfile baby.config.yaml --jobs 4 --profile PROFILE
```

`PROFILE` should be replaced with the name of the profile chosen during profile configuration.

### Cloud execution

An [S3 bucket](https://docs.aws.amazon.com/AmazonS3/latest/userguide) created during the stack deployment will be used for storing the input and output of simulations. Please locate this bucket [here](https://console.aws.amazon.com/s3/buckets) – you can sort buckets by creation time and the newest one should be the one you are looking for, its name should start with the name of the stack.

Open the bucket by clicking on its name. Inside you will find an HTML file which holds detailed instructions on how to execute and monitor the simulation pipeline. This instruction HTML file is almost identical to [this one](http://bioinf.itmat.upenn.edu/BEERS2/stack/instructions/template.html) but the links and commands will be specific to your stack.

### Execution arguments

Execution arguments can be passed to `run_beers` command if executing locally or on a cluster, or given in the input JSON if executing in the cloud. Ultimately, these arguments will be passed down to `snakemake` command which orchestrates the pipeline execution. The full list of arguments can be found [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html), below we list some of the more useful ones.

The `jobs` argument specifies the maximum number of jobs which can be run simultaneously. BEERS2 workload distribution unit is molecule package, so increasing the number of `jobs` beyond the total number of molecule packets (4 in the example above) will not yield further parallelization.

The `directory` argument defines how relative paths in the configuration file are interpreted – they are interpreted relative to `directory`. Moreover, the pipeline results will be placed inside the `directory`. Default directory is the current working directory in local and cluster executions, and the root of the data bucket in cloud executions.

To check that everything is in place, you can use `dry-run`, which will report all steps that will be run.

To re-run a portion of the pipeline, use `force` argument, providing it with a rule name or output file. Rule name could be any of the names of rules in [Snakefile](src/beers/Snakefile), for example `run_library_prep_packet_from_molecule_file`, or you can pass any file produced by BEERS2 and the last rule used to create it will be rerun (as well as any downstream rules).

Depending upon settings like the number of molecules in an input packet, BEERS2 may run out of memory on some steps.
This can be mitigated either by splitting the workload into more packets or by increasing memory limits for specific rules. To accomplish this use the `set-resources` argument which accepts (a sequence of) `RULE:RESOURCE=VALUE` declarations as its value. These declarations specify resource requirements for desired rules. For memory, use `mem_mb` as the `RESOURCE`. Here is an example for local and cluster executions:

```bash
run_beers --configfile ... --set-resources sequence_cluster_packet:mem_mb=32000
```

Analogous input JSON example for cloud execution:

```json
{
  "configfile": "...",
  "set-resources": "sequence_cluster_packet:mem_mb=32000"
}
```

## Removal

To remove local or cluster installation of BEERS2, simply delete the directory which holds the virtual environment in which BEERS2 was installed.

To remove the stack deployed in your AWS account, please follow the [instructions](https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/cfn-console-delete-stack.html). To prevent accidental data loss, the S3 bucket used for storing the input and output will not be deleted in the stack removal process. To delete the bucket after removing the stack, please follow these [instructions](https://docs.aws.amazon.com/AmazonS3/latest/userguide/delete-bucket.html). You may also want to delete logs.

## Simulation configuration

To configure BEERS2 simulation parameters and to specify locations of input files, please copy the configuration template located [here](config/template.config.yaml) and manually edit it.
This template is commented to guide you through selecting or changing values.
Values required to be set are marked with a `TODO` comment.

For locations of input files it is recommended to use absolute paths. Otherwise, paths are interpreted relative to the value of the `directory` argument.

### Input data

The recommended source for input files to BEERS2 is [CAMPAREE](https://github.com/itmat/CAMPAREE).
Alternatively, an easy way to provide input RNA for sequencing is as a molecule file.
These are tab-separated files with the following structure:

```
transcript_id chrom parental_start parental_cigar ref_start ref_cigar strand sequence
```

The `parent_start, parental_cigar` are 1-indexed alignment start position and [CIGAR string](https://genome.sph.umich.edu/wiki/SAM#What_is_a_CIGAR.3F) relative to the molecule's parent molecule, while the `ref_start` and `ref_cigar` are alignments relative to the reference genome.

CAMPAREE considers the simulated genome to be the parental molecule, which differs from the reference genome by any variants that the simulated genome may have that are not present in the reference.
In other cases, it is safe to set `parental_start` and `parent_cigar` to the same as `ref_start` and `ref_cigar`.
Strand should be either `+` or `-` and the sequence is a string of A, C, G, T values. Note that N is not allowed and should be replaced with an A, C, G, or T prior to running BEERS2.

## Output

BEERS2 produces output through the entire process, so that individual aspects of library preparation and sequencing can be investigated.
However, for most uses, the final outputs will be all that is needed.
The final outputs are placed in the `results` subdirectory of `directory` (value of the `directory` argument). Output types can be set in the configuration YAML file and can be:

- FASTQ files – these contain the reads much like a real sequencing run output.
- SAM/BAM files – these contain the reads and quality scores as well as a ideal alignment, giving the true simulated origin of the molecule, even if such an origin can no longer be recreated from its sequence.
  This causes some unusual alignments, such as `100S` alignments that skip all bases but still have a position,
  indicating that the read is entirely adapter sequence or PolyA or other things that do not originate in the genome, but that the original source fragment originated at the specified location.

For all types of output, the read IDs contain additional useful information.
An example below is annotated:

```
@BEERS:1:1:1103:32736:10931:ENSMUST00000055032_1.74.941.cdna.cdna.7	1:N:0:AGCGCTAG+AACCGCGG
@BEERS:{flowcell ID}:{lane ID}:{x}:{y}:{input molecule ID}.{BEERS additional fields, separated by .}  {read number}:N:0:{i5 barcode}+{i7 barcode}
```

Read numbers are either 1 or 2 to denote forward or reverse reads if paired end.
Here the additional fields are added by most steps of the pipeline when producing a new molecule, so that they can be distinguished from other molecules.
For example, the `cdna` markers are added by the first and second strand synthesis steps.
The final number (7 in the example), is added by PCR amplification, meaning that one can identify all PCR duplicates by looking for IDs that differ only in that final number.
If using input from CAMPAREE, the transcript name will have a `_1` or `_2` at the end, indicating which allele it derived from.

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
The log files contain the output molecules from the step, sometimes with additional “notes” on the molecule about what happened in the step, in a tab-separated file.

## Development

BEERS2 has a modular and highly configurable design which allows addition of new steps to the simulation pipeline. If you wish to add new steps, or otherwise modify the code, you need to setup BEERS2 in a slightly different way than described previously. First, clone the desired branch, including the submodules:

```bash
git clone --branch ... --recurse-submodule https://github.com/itmat/BEERS2.git
```

If you have [Visual Studio Code](https://code.visualstudio.com/) with [remote development extension pack](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack), use it to open the cloned BEERS2 directory. A prompt will ask you if you wish to reopen the directory to develop in a container. Click to affirm. This will install all necessary development tools and requirements.

Alternatively, you can execute the first three lines of the [setup script](.devcontainer/setup.sh) inside the BEERS2 directory, after creating and activating virtual environment as described [previously](#local-and-cluster-installation).

With this setup, any edits made to BEERS2 will be reflected when running it. If you wish to deploy your changes in the cloud, please follow the instructions [here](https://github.com/itmat/cask).

## Funding

Work on the BEERS2 project is supported by R21-LM012763-01A1: “The Next Generation of RNA-Seq Simulators for Benchmarking Analyses” (PI: Gregory R. Grant).
