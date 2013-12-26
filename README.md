===========
erfa_python
===========

this is a Python wrapper for ERFA (Essential Routines for 
Fundamental Astronomy).

for installation, see ``doc/installation.rst``

----------
what's new
----------

. 20130815
**********

- registered to PyPi

. 20131222
**********

- added a few routines from the latest SOFA release.

. 20131226
**********

- completed routines from the latest SOFA release
- rename and sligthly reworked the C extension to _erfamodule.c
- a module erfa.py was created especially to handle routines that use
LDBODY. some missing routines were also added into this module.
this change *is compatible* with previous version, and users do not
need to change their code.
