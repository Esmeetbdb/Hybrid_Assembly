from setuptools import setup
#from distutils.core import setup
from Cython.Build import cythonize

import numpy

setup(
    name = 'hybrid assembler',
    version = '1.1.0',

    url = "https://github.com/Esmeetbdb/Hybrid_Assembly",
    author = "Esmee ten Berk de Boer",

    ext_modules = cythonize(["process_fasta.py", "find_overlaps_np_parallel.pyx"]),
    include_dirs=[numpy.get_include()]
)
