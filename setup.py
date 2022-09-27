# Use pip install -e .

from setuptools import setup, find_packages
from Cython.Build import cythonize
import numpy

setup(
    name='beers',
    version='2.0',
    packages=find_packages(),
    scripts=["bin/run_beers"],
    ext_modules = cythonize(["src/beers/sequence/sequence_by_synthesis_helper.pyx"]),
    include_dirs=[numpy.get_include()],
)
