# BEERS2.0

## Installation on PMACS Cluster for Developers

Somewhere under your home directory, clone the develop branch of the BEERS2.0 respository:

```bash
git clone -b develop git@github.com:itmat/BEERS2.0.git
```

Now set up a virtual environment using python 3.6.  I use conda on laptops but on PMACS I stick
with python's venv module and I place the virual environment inside ny project:

```bash
cd BEERS2.0
python3 -m venv ./venv_beers
```

I put venv* in .gitignore so you can use any name you want if you start it with venv and not have
to worry about accidentally committing it.

Now activate the environment thus:

```bash
source ./venv_beers/bin/activate
```

You'll know the virtual environement is activated because the virtual environment path with precede
your terminal prompt.  Now you need to add the python packages/modules upon which BEERS depends.  You
do that by installing the packages/modules listed in the requirements_dev.txt file like so:

```bash
pip install -r requirements_dev.txt
```

The requirements_dev.txt file is supposed to be a superset of the requirements.txt file and in fact,
pulls in the requirements.txt file.  Any packages/modules needed exclusively for development should
be listed in the requirements_dev.txt file.  Requirements needed for a user to run the code should
live in the requirements.txt file.

Next, we need to put the beers package where python can find it.  And this is where the setup.py
file on the top level comes in.  From the top level directory once again, do the following:

```bash
pip install -e .
```

This takes the current directory, packages it and creates a link to the packaged version in
<virtualenv>/lib/python3.6/site-packages.  The file name is beers.egg-link.  This allows python
to find the beer package and subpackages while we can continue to edit them in place.

Next go to the configuration directory and cp config.json to personal config file
(_e.g._, my_config.json).  You can put it anywhere you like.  You will have to reference it
when running beers.  Open your version and modify all the absolute pathnames to conform to your
directory structure.  Modify any parameters you wish to alter and save it.

There is 1 command that you can find in the bin directory under the top level, called run_beers.py.
Calling help on it will show you what is currently possible.

```bash
 ./run_beers.py -h
usage: run_beers.py [-h] -c CONFIG [-r RUN_ID] [-d]
                    {expression_pipeline,library_prep_pipeline,sequence_pipeline}
                    ...

BEERS Simulator, Version 2.0

positional arguments:
  {expression_pipeline,library_prep_pipeline,sequence_pipeline}
                        pipeline subcommand
    expression_pipeline
                        Run the expression pipeline only
    library_prep_pipeline
                        Run the library prep pipeline only
    sequence_pipeline   Run the sequence pipeline only

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

Of the three subcommands, expression_pipeline, library_prep_pipeline, and sequence_pipeline, the
library_prep_pipeline is probably the easiest to run currently.  You would run it from the bin 
directory thus:

```bash
./run_beers -r123 -d -c ../config/my_config.json library_prep_pipeline
```

The run id and the path to the configuration file are both required.  The -d is a debug switch.
Without it, exception tracebacks will not appear.  The library_prep_pipeline currently accepts just
one molecule packet which it locates via the configuration file.  For example:

```json
     "input": {
      "directory_path": "/home/crislawrence/Documents/beers_project/BEERS2.0/data/library_prep",
      "molecule_packet_filename": "molecule_packet_plus_source.pickle"
    }
```

We only have the one packet so it is kind of precious right now.  A copy of
molecule_packet_plus_source.pickle is available under /projects/itmatlab/for_cris  feel free to
grab it.  It has 10K molecules (all polyadenylated) derived from Test_data.1002_baseline.sorted.bam.

I have been using 100 as a seed to get reproducible results, I would suggest others use other
seeds to avoid us getting tunnel vision.

One can use the molecule_packet output from the library prep pipeline as input for the sequence
pipeline but again, you will need to tell the sequence pipeline where to find it via the
configuration file,  For example:

```bash
 "input": {
      "directory_path": "/home/crislawrence/Documents/beers_project/BEERS2.0/data/library_prep/output",
      "molecule_packet_filename": "final_output.pickle"
    }
```

Running the pipeline one stage at a time is a bit inconvenient presently.  We have yet to write
the stages together into a complete pipeline.

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
