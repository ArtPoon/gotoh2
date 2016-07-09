from setuptools import setup, Extension
import numpy as np

setup (name = 'gotoh2',
       version = '0.1',
       description = "C implementation of Gotoh pairwise alignment algorithm to be wrapped in Python",
       py_modules = ['gotoh2.aligner'],
       ext_modules = [Extension('gotoh2.Cgotoh2', sources = ['gotoh2/src/_gotoh2.c'])],
       include_dirs = [np.get_include()],
       zip_safe = False  #avoid headache with permissions on ~/.python-eggs
)
