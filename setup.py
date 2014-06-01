import os
import sys
from distutils.core import setup, Extension

INC_DIR = ['src']
LIB_DIR = []
LIBRARIES = ['m']

extra_compile_args = []

#### is there some way to check this compatibility version on OSX ?
##if sys.platform == 'darwin':
##    extra_compile_args.extend(['-Wl','-compatibility_version,0.0.0','-current_version,0.0.0'])
    
moduleerfa = Extension('_erfa',
                       include_dirs = INC_DIR,
                       libraries = LIBRARIES,
                       library_dirs = LIB_DIR,
                       extra_compile_args = extra_compile_args,
                       sources = ['src/erfa.c',
                                  'src/_erfamodule.c'])

setup (name = 'erfa_python',
       version = '2014.06.02',
       description = 'Python wrapper for ERFA library',
       url = 'https://github.com/nirinA/erfa_python',
       author = 'nirinA raseliarison',
       author_email = 'nirina.raseliarison@gmail.com',
       py_modules=['erfa'],
       ext_modules = [moduleerfa, ],
       license="Public Domain")
