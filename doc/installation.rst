============
installation
============

-----------------------
installing ERFA library
-----------------------

erfa_python requires the erfa_c library to be installed.

    https://github.com/liberfa/erfa

you may need to build a shared library if Python is built with
shared library enable.

erfa headers (erfa.h and erfam.h) are expected to
reside at $HOME/include, /usr/include or /usr/local/include.

libraries (liberfa.a and liberfa.so) are expected to reside
at $HOME/lib, /usr/lib64, /usr/lib, or /usr/local/lib.

if your erfa_c library is in other location, you need to adjust
INC_DIR and LIB_DIR in ``setup.py`` to point out their location.

Note for building erfa_python on OSX:
------------------------------------

on OSX system, to be compatible with erfa_python,
liberfa needs to be compiled with some extra flags:

```
gcc -Wall -O -fPIC -shared -Wl,-compatibility_version,2.0.0,-current_version,2.0.0 -o liberfa.dylib erfa.c
```

see:

    https://github.com/nirinA/erfa_python/issues/4

thanks to @vhaasteren for reporting and solving this issue.

*********************************************
building ERFA against the latest SOFA release
*********************************************

to build the ERFA library that includes the latest
routines released within the SOFA library, you need to
run this script:

    https://github.com/liberfa/erfa-fetch/blob/master/sofa_deriver.py

and run the usual:

    make

to build and install ERFA.

-----------------------------------
building and installing erfa_python
-----------------------------------

now, do the usual:

    python setup.py build

to build and (may need root privilege):

    python setup.py install

to install.

Note: if previous version is already installed in ``site-packages``,
be sure to remove all related files before installing the new version


-------
testing
-------

a test-suit adapted from t_erfa_c.c is provided.
do:

    python erfa_test.py

to validate erfa_python module.

--------------------------
building the documentation
--------------------------

you need sphinx and its dependencies to build the doc.

    cd doc
    
    mkdir _build
    
    sphinx-build . _build

also, usual on-line documentation is available with:

    pydoc erfa

and

   >>> help(erfa)
