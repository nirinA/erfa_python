===========
erfa_python
===========

This is a Python wrapper for ERFA, Essential Routines for 
Fundamental Astronomy, a library derived from SOFA,
Standards of Fundamental Astronomy.

ERFA : https://github.com/liberfa/erfa

SOFA : http://www.iausofa.org

-----------
requirements
-----------

to use this wrapper, you need:

* Python 2.7 - 3.4

* a C compiler to build the module.

--------
download
--------

erfa_python is available from Pypi:

    https://pypi.python.org/pypi/erfa_python

and from github:

    https://github.com/nirinA/erfa_python/
    
------------
installation
------------

unpack the tarball, cd to ``erfa_python`` directory and do the usual:

```
    python setup.py build
```

to build and (may need root privilege):

```
    python setup.py install
```

to install.

Note: if previous version is already installed in ``site-packages``,
be sure to remove all related files before installing the new version

-------
testing
-------

a test-suit adapted from t_erfa_c.c is provided.
do:

```
    python erfa_test.py
```

to validate erfa_python module.

-------------
documentation
-------------

on-line documentation is available with:

```
    pydoc erfa
```

to generate a html version:

```
    pydoc -w erfa
```
    
and, from the interpreter:

```
   >>> help(erfa)
   >>> help(erfa.<function>)
```

C-api reference for liberfa is available in ``doc/api.rst``.
a html version can be built:

```
    cd doc   
    mkdir _build  
    sphinx-build . _build
```

you need ``sphinx`` and its dependencies to build the doc.

--------
examples
--------

some examples adapted from SOFA documentation
can be found in ``examples/`` directory.

***********
. ex_ast.py
***********

  example of Astrometry Tools

***********
. ex_pn.py
***********

  example of tools for Earth Attitude 

***********
. ex_ts.py
***********

  example of Time Scale and Calendar Tools

***********
. ex_ephem.py
***********

  basic example of Ephemerids for Sun and Major planets

----------
what's new
----------

***********
. 20140601
**********

- ``erfa`` source is included, so a pre-installed ``liberfa`` is no longer required.

***********
. 20140530
**********

- update ex_ast.py
- fix _erfamodule.c docstring

***********
. 20140127
**********

- support for Python-3.4.

***********
. 20131226
**********

- completed routines from the latest SOFA release.
- rename and slightly reworked the C extension to _erfamodule.c

***********
. 20131222
**********

- added a few routines from the latest SOFA release.

***********
. 20130830
**********

- registered to PyPi

***********
. 20130815
**********

- rewrote the wrapper from SOFA to ERFA

