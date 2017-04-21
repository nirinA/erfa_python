===========
erfa_python
===========

This is a Python wrapper for ERFA, Essential Routines for 
Fundamental Astronomy, a library derived and 'freed' from SOFA,
Standards of Fundamental Astronomy.

ERFA : https://github.com/liberfa/erfa

SOFA : http://www.iausofa.org

-----------
requirements
-----------

to use this wrapper, you need:

* Python 2.7 - 3.6

* C compiler to build the module.

  . newer version of gcc (4.7+) is prefered.

----------
obtaining
----------

***********
download tarball
***********

latest version is available from pypi:

    https://pypi.python.org/pypi/erfa_python

or from github:

    https://github.com/nirinA/erfa_python/releases

**********
clone from github
**********

you can clone the project repository at:
    
    https://github.com/nirinA/erfa_python/


------------
building
------------

* to build the module, unpack the tarball,
cd to ``erfa_python`` directory and do the usual:

```
    python setup.py build
```

note for OS X:
**************

* on OSX platform, it seems that ``erfa`` needs to be build as a library,
and some 'compatibility version' needs to be adjusted
thanks to @vhaasteren for reporting and solving this issue, see:

https://github.com/nirinA/erfa_python/issues/4

------------
installing
------------

* if previous version of erfa_python is already installed in ``site-packages``,
be sure to remove all related files before installing the new version.

* install (may need root privilege) with:

```
    python setup.py install
```

-------
testing
-------

a test-suite adapted from t_erfa_c.c is provided.
do:

```
    python erfa_test.py
```

to validate erfa_python module.

note:
*****
some test fail with older version of gcc (4.5).
thanks to @yanwang2012 for this report, see:
    
https://github.com/nirinA/erfa_python/issues/6

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

C-api reference of liberfa is available in ``doc/api.rst``.
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
. 20170421
***********

- update from SOFA release 13.
- this release replaces the IAU 1976 value of astronomical unit with the IAU 2012 value.
- codes are all backward compatible but results are not, especially those using ```erfa.ASTROM```
- ```erfa_test.py``` has been updated accordingly


***********
. 2016123
***********

- update from SOFA release 12, revised on 2016-12-23.
- tests for fk52h and starpm are relaxed, as they failed with Python3.6 built with all optimization enabled.

***********
. 20160503
***********

- add new routines from SOFA release 12.


***********
. 20150211
***********

- add new functions ``icrs2g`` and ``g2icrs`` from latest SOFA release

***********
. 20141016
**********

- upload a windows binary to PyPi.

***********
. 20140921
**********

- update ``erfa.c`` from latest ``sofa`` release.

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
- renamed and slightly reworked the C extension to _erfamodule.c

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

