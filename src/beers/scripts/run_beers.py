def main():
    '''
    This script runs the BEERS2 using Snakemake.
    All arguments passed to this will be given to the `snakemake` command.
    '''
    import argparse

    parser = argparse.ArgumentParser(
        description="BEERS2 RNA-seq Simulator",
        add_help=False,
    )

    args, other_args = parser.parse_known_args()

    import beers
    import importlib.resources
    # Find the Snakefile relative to the beers module
    snakemake_file = importlib.resources.path(beers, "Snakefile")

    # Start the snakemake
    cmd = [
        "snakemake",
        "--snakefile",
        snakemake_file,
        *other_args
    ]
    import sys
    if sys.stdout.isatty():
        import pty
        # We use pty.spawn not subprocess.run in order to
        # get the nice colored output from Snakemake
        pty.spawn(cmd)
    else:
        # Run Snakemake without tty
        import subprocess
        try:
            subprocess.run(
                cmd,
                check=True
            )
        except subprocess.CalledProcessError:
            exit(1) # Had an error, want to signal that the subprocess failed
