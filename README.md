===========
erfa_python
===========

This is a Python wrapper for ERFA, Essential Routines for 
Fundamental Astronomy, a library derived from SOFA,
Standards of Fundamental Astronomy.

-----------
requirements
-----------

To use this wrapper, you need:

    * a working liberfa:
        
    https://github.com/liberfa/erfa

    * Python 2.7 to 3.4

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

unpack the tarball, cd to ``erfa_python`` directory and:
    
```
python setup.py install
```

for a complete installation procedure, including
the installation of ERFA, see ``doc/installation.rst``

--------
examples
--------

some examples adapted from SOFA documentation
can be found in ``examples/`` directory.

. ex_ast.py
***********

  example of Astrometry Tools

. ex_pn.py
***********

  example of tools for Earth Attitude 

. ex_ts.py
***********

  example of Time Scale and Calendar Tools

. ex_ephem.py
***********

  basic example of Ephemerids for Sun and Major planets

----------
what's new
----------

. 20140530
**********

- update ex_ast.py
- fix _erfamodule.c docstring

. 20140127
**********

- support for Python-3.4.

. 20131226
**********

- completed routines from the latest SOFA release.
- rename and slightly reworked the C extension to _erfamodule.c

. 20131222
**********

- added a few routines from the latest SOFA release.

. 20130830
**********

- registered to PyPi

. 20130815
**********

- rewrote the wrapper from SOFA to ERFA

