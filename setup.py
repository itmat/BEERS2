# Use pip install -e .

from setuptools import setup, find_packages

setup(
    name='beers',
    version='2.0',
    packages=find_packages(),
    scripts=["bin/run_beers.py"],
)
