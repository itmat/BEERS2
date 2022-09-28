from setuptools import setup, find_packages
from Cython.Build import cythonize
import numpy

setup(
    name='beers',
    version='2.0',
    packages=find_packages("src"),
    package_dir={"": "src"},
    scripts=["bin/run_beers"],
    ext_modules = cythonize(["src/beers/sequence/sequence_by_synthesis_helper.pyx"], language_level=3),
    package_data = {"beers": ["Snakefile"]},
    include_dirs=[numpy.get_include()],
)
