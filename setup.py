from distutils.core import setup, Extension

INC_DIR = ['$HOME/include',
           '/usr/include',
           '/usr/local/include']
LIB_DIR = ['$HOME/lib',
           '/usr/lib64',
           '/usr/lib',
           '/usr/local/lib',]

moduleerfa = Extension('erfa',
                    include_dirs = INC_DIR,
                    libraries = ['erfa', 'm'],
                    library_dirs = LIB_DIR,
                    sources = ['src/erfamodule.c'])

setup (name = 'erfa_python',
       version = '2013.08.15',
       description = 'Python wrapper for ERFA library',
       ext_modules = [moduleerfa, ],
       license="Public Domain")
