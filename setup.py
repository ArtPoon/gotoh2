from setuptools import setup, Extension
import numpy as np

gotoh2 = Extension('Cgotoh2',
                    sources = ['_gotoh2.c'],
                    define_macros=[('__PYTHON__', None)])

setup (name = 'Cgotoh2',
       version = '0.1',
       description = "C implementation of Gotoh pairwise alignment algorithm to be wrapped in Python",
       ext_modules = [gotoh2],
       include_dirs = [np.get_include()],
       zip_safe = False  #avoid headache with permissions on ~/.python-eggs
)
