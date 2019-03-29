# CAMPAREE

## Installation on PMACS Cluster for Developers

Somewhere under your home directory, clone the develop branch of the CAMPAREE repository:

```bash
git clone -b develop git@github.com:itmat/CAMPAREE.git
```

Now set up a virtual environment using python 3.6.  I use conda on laptops but on PMACS I stick
with python's venv module and I place the virual environment inside ny project:

```bash
cd CAMPAREE
python3 -m venv ./venv_camparee
```

I put venv* in .gitignore so you can use any name you want if you start it with venv and not have
to worry about accidentally committing it.

Now activate the environment thus:

```bash
source ./venv_camparee/bin/activate
```

You'll know the virtual environement is activated because the virtual environment path with precede
your terminal prompt.  Now you need to add the python packages/modules upon which CAMPAREE depends.
You do that by installing the packages/modules listed in the requirements_dev.txt file like so:

```bash
pip install -r requirements_dev.txt
```

The requirements_dev.txt file is supposed to be a superset of the requirements.txt file and in fact,
pulls in the requirements.txt file.  Any packages/modules needed exclusively for development should
be listed in the requirements_dev.txt file.  Requirements needed for a user to run the code should
live in the requirements.txt file.

Next, we need to put the CAMPAREE package where python can find it.  And this is where the setup.py
file on the top level comes in.  From the top level directory once again, do the following:

```bash
pip install -e .
```

This takes the current directory, packages it and creates a link to the packaged version in
<virtualenv>/lib/python3.6/site-packages.  The file name is camparee.egg-link.  This allows
python to find the CAMPAREE package and subpackages while we can continue to edit them in place.

Similarly, we need to point the python toward the BEERS_UTILS package as well. Again, from the
top-level CAMPAREE directory, run the following:

```bash
pip install -e /path/to/BEERS_UTILS
```

Next go to the configuration directory and cp config.json to personal config file
(_e.g._, my_config.json).  You can put it anywhere you like.  You will have to reference it
when running CAMPAREE.  Open your version and modify all the absolute pathnames to conform to
your directory structure.  Modify any parameters you wish to alter and save it.

There is 1 command that you can find in the bin directory under the top level, called run_camparee.py.
Calling help on it will show you what is currently possible.

```bash
 ./run_camparee.py -h
usage: run_camparee.py [-h] -c CONFIG [-r RUN_ID] [-d]


CAMPAREE - RNA molecule simulator

optional arguments:
  -h, --help            show this help message and exit

required named arguments:
  -c CONFIG, --config CONFIG
                        Full path to configuration file.

optional named arguments - these override configuration file arguments.:
  -r RUN_ID, --run_id RUN_ID
                        Integer used to specify run id.
  -d, --debug           Indicates whether additional diagnostics are printed.

```

The run id and the path to the configuration file are both required.  The -d is a debug switch.
Without it, exception tracebacks will not appear.

I have been using 100 as a seed to get reproducible results, I would suggest others use other
seeds to avoid us getting tunnel vision.

The expression pipeline is more difficult to use as it requires the reference genome and the pair
of alignment files presently (bam and bai) and really only runs the variants finder portion of
the pipeline.  I threw in a BeagleStep that will eventually call the Beagle process.  For now, I
put my own Java program as a parameter to that step so I'd have something to run.  You can
add your own external process as a placemarker for now, if you like.

## Requirements for Users

If the user chooses to supply his/her own reference genome, it should be edited so that a
sequence contains no line breaks.

If the user declines to provide gender for each sample, the sample will not have X,Y, MT
data.  If the user neglects to provide gender for just some of the samples, X,Y,MT data
will be generated for those samples that have gender and a warning will be issued to
the user.
