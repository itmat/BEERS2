# BEERS2

## Release status

This is publicly funded research, therefore for the sake of full transparency we are developing BEERS2 as a public repository. Please stay tuned until our first beta release, sometime later this summer.

In the meantime please check out CAMPAREE, a member of the BEERS2 suite of tools that simulates the collection of RNA molecules in a biological sample, prior to RNA-seq library preparation.
https://github.com/itmat/CAMPAREE

## Installation on LSF cluster for developers

Somewhere under your home directory, clone the develop branch of the BEERS2 repository:

```bash
git clone -b develop git@github.com:itmat/BEERS2.git
```

Now set up a virtual environment using python 3.6.

```bash
cd BEERS2
python3 -m venv ./venv_beers
```

I put venv* in .gitignore so you can use any name you want if you start it with venv and not have to worry about accidentally committing it.

Now activate the environment thus:

```bash
source ./venv_beers/bin/activate
```

You'll know the virtual environment is activated because the virtual environment path with precede
your terminal prompt.  Now you need to add the python packages/modules upon which BEERS2 depends.  You
do that by installing the packages/modules listed in the requirements_dev.txt file like so:

```bash
pip install -r requirements_dev.txt
```

The requirements_dev.txt file is supposed to be a superset of the requirements.txt file and in fact,
pulls in the requirements.txt file.  Any packages/modules needed exclusively for development should
be listed in the requirements_dev.txt file.  Requirements needed for a user to run the code should
live in the requirements.txt file.

Next, we need to put the BEERS2 package where python can find it.  And this is where the setup.py
file on the top level comes in.  From the top level directory once again, do the following:

```bash
pip install -e .
```

This takes the current directory, packages it and creates a link to the packaged version in
<virtualenv>/lib/python3.6/site-packages.  The file name is beers.egg-link.  This allows python
to find the beer package and subpackages while we can continue to edit them in place.

## Funding

Work on the BEERS2 project is supported by R21-LM012763-01A1: “The Next Generation of RNA-Seq Simulators for Benchmarking Analyses” (PI: Gregory R. Grant).
