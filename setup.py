# Use pip install -e .

from setuptools import setup, find_packages
from Cython.Build import cythonize

setup(
    name='beers',
    version='2.0',
    packages=find_packages(),
    scripts=["bin/run_beers.py"],
    ext_modules = cythonize(["beers/sequence/sequence_by_synthesis_helper.pyx"]),
)
