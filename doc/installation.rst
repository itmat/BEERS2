Installation
============

First set up a virtual environment using Python 3.11::

    python -m venv venv_beers
    source venv_beers/bin/activate

You'll know the virtual environment is activated because the virtual environment path with precede
your terminal prompt.
Now we need to install the BEERS2 package.
Note that we do not recommend doing this outside of the newly created virtual environment::

    pip install git+https://github.com/itmat/BEERS2


Now BEERS2 should be installed.
We can verify correct installation by running an example dataset.

Example Dataset
===============

We will run a simplified 'baby' example, with a reduced mouse genome and provided output to test the installation and Snakemake configuration.
We first download the example data and config file::

    wget -c https://s3.amazonaws.com/itmat.data/BEERS2/examples/baby_mouse_example.tar.gz -O - | tar -xz
    cd baby_mouse_example/

Now we are ready to run BEERS2 using this dataset::

    beers --configfile baby.config.yaml --jobs 1

Verify that the run has been successful by examing :code:`results/` which should contain output FASTQ and SAM files.
To confirm that the expected results were produced, compare md5 hashes to this reference::

    $md5sum results/*sam

    c02704029f6e16a77ea6890ebf92375a  results/S1_L1.sam
    35dca3ee057e7a8bac8b8a5d790b054f  results/S1_L2.sam
    7c963445bab5153c63c87d507ba589aa  results/S1_unidentifiedL1.sam
    00dcaa96464e2e207aa2989650d0b8f3  results/S1_unidentifiedL2.sam
    e1b0b04a7e3ca9274b5fa59b03479f97  results/S2_L1.sam
    4a7e526ea6b333b19ff40689f659e680  results/S2_L2.sam
    19937368990e8bff3ed6b86a8a6718fc  results/S2_unidentifiedL1.sam
    bef563a2384e92961bd9701e5f5e2720  results/S2_unidentifiedL2.sam


The :code:`--jobs 1` option sets to run this on a single-core locally.
Increasing this number will allow Snakemake to run multiple processes simultaneously on the machine your execute the :code:`beers` command from.
