Configuring
==================

To configure BEERS, setting the input files and parameters for generating the data, please copy the template config from `config/template.config.yaml <config/template.config.yaml>`_ and manually edit it.
This config is commented to guide you through selecting or changing values.
Values required to be set are marked with a TODO comment.

Configuration files are recommended to use absolute file paths to input files.
Otherwise, paths are interpreted relative to the set :code:`--directory` output directory or relative to the BEERS2 directory.

Cluster Execution
-----------------

BEERS2 uses Snakemake to run its pipeline, which allows distribution across multiple cores or across a cluster environment.
Many computing environments (including LSF, SGE or SLURM) already have configuration available at the `Snakemake-profiles <https://github.com/Snakemake-Profiles/doc>`_ project.
Choose your environment and follow the instructions there, then use the installed Snakemake profile with the :code:`--profile` option when running BEERS2 (see below).
After configuring the Snakemake profile, our example run command could change to::

    run_beers --configfile baby.config.yaml --jobs 4 --profile MYPROFILENAME

Where :code:`MYPROFILENAME` should be replaced with the name of the profile you chose during profile configuration. This will tell Snakemake to execute BEERS2 commands on other machines in your cluster.
The :code:`--jobs` options should be changed according to the number of jobs to use simultaneously; at most one per molecule packet provided as input will be used (4 for the example we ran).
BEERS2 will use up to one job per packet at a time, so increasing :code:`jobs` above that number will not cause further parallelism.

Useful Snakemake commands
^^^^^^^^^^^^^^^^^^^^^^^^^

The :code:`run_beers` command is a thin wrapper around :code:`snakemake` command, pointing it to the correct :code:`Snakefile` to run, so any commands passed to :code:`run_beers` will be given to :code:`snakemake`.
In addition to the already mentioned :code:`--profile`, :code:`--jobs` and :code:`--directory` values, the following are useful to run BEERS.
To check that everything is in place, use :code:`--dry-run`, which will report all steps that will be run.
To re-run a portion of the pipeline, use :code:`--force-rerun {rule name or output file}` where :code:`rule name` could be any of the names of rules in :code:`Snakefile`, for example :code:`run_library_prep_packet_from_molecule_file`,
or you can pass any file produced by BEERS2 and the last rule used to create it will be rerun (as well as any downstream rules).

Depending upon settings like the number of molecules in an input packet, BEERS2 may run out of memory on some steps.
This can be mitigated either by splitting into more packets or by increasing memory limits by passing :code:`--set-resources {rule_name}:mem_mb=32000` or similar with the name of the rule that ran out of memory.



Input
===============

The recommended source for input files to BEERS2 is `CAMPAREE <https://github.com/itmat/CAMPAREE>`_.
Alternatively, an easy way to provide input RNA to sequence is as a molecule file.
These are tab-separated files with the following structure::

    transcript_id chrom parental_start parental_cigar ref_start ref_cigar strand sequence

The :code:`parent_start, parental_cigar` are 1-indexed alignment start position and `cigar string <https://genome.sph.umich.edu/wiki/SAM#What_is_a_CIGAR.3F>`_ relative to the molecule's parent molecule, while the :code:`ref_start` and :code:`ref_cigar` are alignments relative to the reference genome.
CAMPAREE considers the simulated genome to be the parental molecule, which differs from the reference genome by any variants that the simulated genome may have that are not present in the reference.
In other cases, it is safe to set :code:`parental_start` and :code:`parent_cigar` to the same as :code:`ref_start` and :code:`ref_cigar`.
Strand should be either "+" or "-" and the sequence is a string of ACGT values (N is not allowed and should be replaced with an A, C, G, or T prior to running BEERS).

Output
==================

BEERS2 produces output through the entire process, so that individual aspects of library preparation and sequencing can be investigated.
However, for most uses, the final outputs will be all that is needed.
These are placed in the :code:`results/` subdirectory of whatever output directory was set with :code:`--directory`.
Output types can be set in the config but are either:

 - FASTQ file: these contain the reads much like a real sequencing run output.
 - SAM/BAM: these contain the reads and quality scores as well as a ideal alignment, giving the true simulated origin of the molecule, even if such an origin can no longer be recreated from its sequence.
   This causes some unusual alignments, such as "100S" alignments that skip all bases but still have a position,
   indicating that the read is entirely adapter sequence or PolyA or other things that do not originate in the genome, but that the original source fragment originated at the specified location.

For all types of output, the read ids contain additional useful information.
An example below is annotated::

    @BEERS:1:1:1103:32736:10931:ENSMUST00000055032_1.74.941.cdna.cdna.7	1:N:0:AGCGCTAG+AACCGCGG
    @BEERS:{flowcell id}:{lane id}:{x}:{y}:{input_molecule ID}.{BEERS additional fields, separated by .}  {read number}:N:0:{i5 barcode}+{i7 barcode}

Read numbers are either 1 or 2 to denote forward or reverse reads if paired end.
Here the additional fields are added by most steps of the pipeline when producing a new molecule, so that they can be distinguished from other molecules.
For example, the 'cdna' markers are added by the FirstStrandSynthesisStep and SecondStrandSynthesisStep steps.
The final number (7 in the example), is added by PCR amplification, meaning that one can identify all PCR duplicates by looking for ids that differ only in that final number.
If using input from CAMPAREE, the transcript name will have a :code:`_1` or :code:`_2`  at the end, indicating which allele it derived from.

