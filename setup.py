from setuptools import setup, find_packages
from Cython.Build import cythonize
import numpy

setup(
    name='beers',
    version='2.0',
    packages=find_packages("src"),
    python_requires=">=3.11",
    package_dir={"": "src"},
    entry_points={
        'console_scripts': [
            'run_beers=beers.scripts.run_beers:main',
        ],
    },
    ext_modules = cythonize(["src/beers/sequence/sequence_by_synthesis_helper.pyx"], language_level=3),
    package_data = {"beers": ["Snakefile"]},
    include_dirs=[numpy.get_include()],
)
