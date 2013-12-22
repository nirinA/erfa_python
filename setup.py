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
       version = '2013.12.03',
       description = 'Python wrapper for ERFA library',
       url = 'https://github.com/nirinA/erfa_python',
       author = 'nirinA raseliarison',
       author_email = 'nirina.raseliarison@gmail.com',
       ext_modules = [moduleerfa, ],
       classifiers = [
           'Development Status :: alpha',
           'Programming Language :: Python :: 2.7',
           'Programming Language :: Python :: 3.x',
           'Topic :: Astronomy',
           ],
       license="Public Domain")
