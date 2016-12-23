import os
import sys
import subprocess
from distutils.core import setup, Extension

INC_DIR = ['src']
LIB_DIR = []
LIBRARIES = ['m']
sources = ['src/erfa.c', 'src/_erfamodule.c']
extra_compile_args = []

if sys.platform == 'darwin':
#### is there some way to check this compatibility version on OSX ?
##    extra_compile_args.extend(['-Wl','-compatibility_version,0.0.0','-current_version,0.0.0'])
    sources.pop(0)
    subprocess.Popen(('gcc', '-Wall', '-O', 'fPIC', '-shared', '-I./src', 
                      '-Wl,-compatibility_version,2.0.0,-current_version,2.0.0',
                      '-o', 'liberfa.dylib', './src/erfa.c', '-lm'))
    LIBRARIES.append('erfa')
    LIB_DIR.append('.')

moduleerfa = Extension('_erfa',
                       include_dirs = INC_DIR,
                       libraries = LIBRARIES,
                       library_dirs = LIB_DIR,
                       extra_compile_args = extra_compile_args,
                       sources = sources)

setup (name = 'erfa_python',
       version = '2016.12.24',
       description = 'Python wrapper for ERFA library',
       url = 'https://github.com/nirinA/erfa_python',
       author = 'nirinA raseliarison',
       author_email = 'nirina.raseliarison@gmail.com',
       py_modules=['erfa'],
       ext_modules = [moduleerfa, ],
       license="Public Domain")
