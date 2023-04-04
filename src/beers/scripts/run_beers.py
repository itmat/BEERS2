def main():
    """
    Run the BEERS2 simulation pipeline.
    Arguments will be passed down to the Snakemake orchestrator.
    """
    import argparse
    import beers
    import importlib.resources
    import snakemake
    import sys

    parser = argparse.ArgumentParser(
        description="BEERS2 RNA-Seq simulator",
        add_help=False,
    )

    parser.add_argument(
        "--configfile",
        help="Location of the simulation configuration YAML file",
        required=True,
    )

    args, rest = parser.parse_known_args()

    snakefile = importlib.resources.files(beers) / "Snakefile"

    arguments = ["--configfile", args.configfile]
    arguments += ["--snakefile", str(snakefile)]
    arguments += rest

    sys.argv = ["snakemake"] + arguments
    snakemake.main()
