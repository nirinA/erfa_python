import os, sys
from distutils.core import setup, Extension

INC_DIR = ['$HOME/include',
           '/usr/include',
           '/usr/local/include']
LIB_DIR = ['$HOME/lib',
           '/usr/lib64',
           '/usr/lib',
           '/usr/local/lib',]

extra_compile_args = []

#### is there some way to check this compatibility version on OSX ?
##if sys.platform == 'darwin':
##    extra_compile_args.append(['-Wl','-compatibility_version,0.0.0','-current_version,0.0.0'])
    
moduleerfa = Extension('_erfa',
                       include_dirs = INC_DIR,
                       libraries = ['erfa', 'm'],
                       library_dirs = LIB_DIR,
                       extra_compile_args = extra_compile_args,
                       sources = ['src/_erfamodule.c'])

setup (name = 'erfa_python',
       version = '2014.05.30',
       description = 'Python wrapper for ERFA library',
       url = 'https://github.com/nirinA/erfa_python',
       author = 'nirinA raseliarison',
       author_email = 'nirina.raseliarison@gmail.com',
       py_modules=['erfa'],
       ext_modules = [moduleerfa, ],
       license="Public Domain")
