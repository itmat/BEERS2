Installation
============

First set up a virtual environment using Python 3.10::

    python -m venv ./venv_beers
    source ./venv_beers/bin/activate

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

    $ md5sum reuslts/*sam
    beaa4988ffa1cc50fda5d14b0dfef7df  results/S1_L1.sam
    6656813664c7b91db480c5dd5a3ab6d0  results/S1_L2.sam
    a5db0b010d9407a76b3da4452073600b  results/S1_unidentified_L1.sam
    ed257189a7d5915046e5327b70afa1d5  results/S1_unidentified_L2.sam
    fd73b1c094e9256713c33e065a287ca6  results/S2_L1.sam
    1f7d6df7bb21161245ced8486d5fb487  results/S2_L2.sam
    2ffa4652a19ab01436408c1da1314398  results/S2_unidentified_L1.sam
    d1124ed6d0149d977e52aafd862a0f6a  results/S2_unidentified_L2.sam

The :code:`--jobs 1` option sets to run this on a single-core locally.
Increasing this number will allow Snakemake to run multiple processes simultaneously on the machine your execute the :code:`beers` command from.
