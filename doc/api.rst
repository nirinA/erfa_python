===============
ERFA: notes and references
===============

------
a2af
------

Given:
lution (Note 1)
e in radians

or '-'
ees, arcminutes, arcseconds, fraction

se days to hms

interpreted as follows:
on
0
0
0
0
0
0
0
0
1
0.1
0.01
0.001
0.000...
 useful value for ndp is determined by the
ormat of doubles on the target platform, and
ing idmsf[3].  On a typical platform, for
 available floating-point precision might
.  However, the practical limit is typically
pacity of a 32-bit int, or ndp=4 if int is

f angle may exceed 2pi.  In cases where it
o the caller to test for and handle the
very nearly 2pi and rounds up to 360 degrees,
[0]=360 and setting idmsf[0-3] to zero.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
a2tf
------

Given:
lution (Note 1)
e in radians

or '-'
s, minutes, seconds, fraction

se days to hms

interpreted as follows:
on
0
0
0
0
0
0
0
0
1
0.1
0.01
0.001
0.000...
 useful value for ndp is determined by the
ormat of doubles on the target platform, and
ing ihmsf[3].  On a typical platform, for
 available floating-point precision might
.  However, the practical limit is typically
pacity of a 32-bit int, or ndp=4 if int is

f angle may exceed 2pi.  In cases where it
o the caller to test for and handle the
very nearly 2pi and rounds up to 24 hours,
[0]=24 and setting ihmsf[0-3] to zero.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
ab
------

Given:
atural direction to the source (unit vector)
bserver barycentric velocity in units of c
istance between the Sun and the observer (au)
qrt(1-|v|^2): reciprocal of Lorenz factor

roper direction to source (unit vector)

ed on Expr. (7.40) in the Explanatory
Seidelmann 2013), but with the following

han approximate normalization is applied.
 potential term from Expr. (7) in
 added, taking into account only the Sun's
is has a maximum effect of about
d.
 the maximum accuracy will be limited by the
For example, if the ERFA eraEpv00 function is
o 5 microarcseconds could occur.

nn, P. K. (eds), Explanatory Supplement to
anac, 3rd ed., University Science Books

"A practical relativistic model for micro-
 in space", Astr. J. 125, 1580-1597 (2003).

product of two p-vectors
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The algorithm is based on Expr. (7.40) in the Explanatory
     Supplement (Urban & Seidelmann 2013), but with the following
     changes:

     o  Rigorous rather than approximate normalization is applied.

     o  The gravitational potential term from Expr. (7) in
        Klioner (2003) is added, taking into account only the Sun's
        contribution.  This has a maximum effect of about
        0.4 microarcsecond.

  2) In almost all cases, the maximum accuracy will be limited by the
     supplied velocity.  For example, if the ERFA eraEpv00 function is
     used, errors of up to 5 microarcseconds could occur.

  References:

     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
     the Astronomical Almanac, 3rd ed., University Science Books
     (2013).

     Klioner, Sergei A., "A practical relativistic model for micro-
     arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).


References:

     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
     the Astronomical Almanac, 3rd ed., University Science Books
     (2013).

     Klioner, Sergei A., "A practical relativistic model for micro-
     arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).

  Called:
     eraPdp       scalar product of two p-vectors

------
af2a
------

Given:
gn:  '-' = negative, otherwise positive
grees
cminutes
cseconds

gle in radians
e):
atus:  0 = OK
       1 = ideg outside range 0-359
       2 = iamin outside range 0-59
       3 = asec outside range 0-59.999...

ted even if any of the range checks fail.
n and/or asec produce a warning status, but
is used in the conversion.
le errors, the status value reflects only the
 taking precedence.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
anp
------

Given:
angle (radians)
e):
angle in range 0-2pi
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
anpm
------

Given:
angle (radians)
e):
angle in range +/-pi
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
apcg
------

Given:
TDB as a 2-part...
...Julian Date (Note 1)
Earth barycentric pos/vel (au, au/day)
Earth heliocentric position (au)

star-independent astrometry parameters:
 PM time interval (SSB, Julian years)
 SSB to observer (vector, au)
 Sun to observer (unit vector)
 distance from Sun to observer (au)
 barycentric observer velocity (vector, c)
 sqrt(1-|v|^2): reciprocal of Lorenz factor
 bias-precession-nutation matrix
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged

ate2 is a Julian Date, apportioned in any
en the two arguments.  For example,
uld be expressed in any of these ways, among

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 method is best matched to the way the
internally and will deliver the optimum
 method and the date & time methods are both
ween resolution and convenience.  For most
 function the choice will not be at all

ad of TDB without any significant impact on

with respect to BCRS axes.
al functions that inserts into the astrom
endent parameters needed for the chain of
mations ICRS <-> GCRS <-> CIRS <-> observed.
s support different classes of observer and
sformation chain:
   observer        transformation
   geocentric      ICRS <-> GCRS
   terrestrial     ICRS <-> CIRS
   terrestrial     ICRS <-> observed
   space           ICRS <-> GCRS
   terrestrial     update Earth rotation
   terrestrial     CIRS <-> observed
ing in "13" use contemporary ERFA models to
ephemerides.  The others accept ephemerides
er.
rom ICRS to GCRS covers space motion,
ection, and aberration.  From GCRS to CIRS
 and precession-nutation.  From CIRS to
nt of Earth rotation, polar motion, diurnal
lax (unless subsumed into the ICRS <-> GCRS
 atmospheric refraction.
e astrom produced by this function is used by
cq*.

try parameters, ICRS-GCRS, space observer
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TDB date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TDB)=2450123.7 could be expressed in any of these ways, among
     others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in cases
     where the loss of several decimal digits of resolution is
     acceptable.  The J2000 method is best matched to the way the
     argument is handled internally and will deliver the optimum
     resolution.  The MJD method and the date & time methods are both
     good compromises between resolution and convenience.  For most
     applications of this function the choice will not be at all
     critical.

     TT can be used instead of TDB without any significant impact on
     accuracy.

  2) All the vectors are with respect to BCRS axes.

  3) This is one of several functions that inserts into the astrom
     structure star-independent parameters needed for the chain of
     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

     The various functions support different classes of observer and
     portions of the transformation chain:

          functions         observer        transformation

       eraApcg eraApcg13    geocentric      ICRS <-> GCRS
       eraApci eraApci13    terrestrial     ICRS <-> CIRS
       eraApco eraApco13    terrestrial     ICRS <-> observed
       eraApcs eraApcs13    space           ICRS <-> GCRS
       eraAper eraAper13    terrestrial     update Earth rotation
       eraApio eraApio13    terrestrial     CIRS <-> observed

     Those with names ending in "13" use contemporary ERFA models to
     compute the various ephemerides.  The others accept ephemerides
     supplied by the caller.

     The transformation from ICRS to GCRS covers space motion,
     parallax, light deflection, and aberration.  From GCRS to CIRS
     comprises frame bias and precession-nutation.  From CIRS to
     observed takes account of Earth rotation, polar motion, diurnal
     aberration and parallax (unless subsumed into the ICRS <-> GCRS
     transformation), and atmospheric refraction.

  4) The context structure astrom produced by this function is used by
     eraAtciq* and eraAticq*.


------
apcg13
------

Given:
B as a 2-part...
.Julian Date (Note 1)

ar-independent astrometry parameters:
 PM time interval (SSB, Julian years)
 SSB to observer (vector, au)
 Sun to observer (unit vector)
 distance from Sun to observer (au)
 barycentric observer velocity (vector, c)
 sqrt(1-|v|^2): reciprocal of Lorenz factor
 bias-precession-nutation matrix
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged

ate2 is a Julian Date, apportioned in any
en the two arguments.  For example,
uld be expressed in any of these ways, among

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 method is best matched to the way the
internally and will deliver the optimum
 method and the date & time methods are both
ween resolution and convenience.  For most
 function the choice will not be at all

ad of TDB without any significant impact on

with respect to BCRS axes.
aller wishes to supply his own Earth
ion eraApcg can be used instead of the present

al functions that inserts into the astrom
endent parameters needed for the chain of
mations ICRS <-> GCRS <-> CIRS <-> observed.
s support different classes of observer and
sformation chain:
   observer        transformation
   geocentric      ICRS <-> GCRS
   terrestrial     ICRS <-> CIRS
   terrestrial     ICRS <-> observed
   space           ICRS <-> GCRS
   terrestrial     update Earth rotation
   terrestrial     CIRS <-> observed
ing in "13" use contemporary ERFA models to
ephemerides.  The others accept ephemerides
er.
rom ICRS to GCRS covers space motion,
ection, and aberration.  From GCRS to CIRS
 and precession-nutation.  From CIRS to
nt of Earth rotation, polar motion, diurnal
lax (unless subsumed into the ICRS <-> GCRS
 atmospheric refraction.
e astrom produced by this function is used by
cq*.

osition and velocity
try parameters, ICRS-GCRS, geocenter
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TDB date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TDB)=2450123.7 could be expressed in any of these ways, among
     others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in cases
     where the loss of several decimal digits of resolution is
     acceptable.  The J2000 method is best matched to the way the
     argument is handled internally and will deliver the optimum
     resolution.  The MJD method and the date & time methods are both
     good compromises between resolution and convenience.  For most
     applications of this function the choice will not be at all
     critical.

     TT can be used instead of TDB without any significant impact on
     accuracy.

  2) All the vectors are with respect to BCRS axes.

  3) In cases where the caller wishes to supply his own Earth
     ephemeris, the function eraApcg can be used instead of the present
     function.

  4) This is one of several functions that inserts into the astrom
     structure star-independent parameters needed for the chain of
     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

     The various functions support different classes of observer and
     portions of the transformation chain:

          functions         observer        transformation

       eraApcg eraApcg13    geocentric      ICRS <-> GCRS
       eraApci eraApci13    terrestrial     ICRS <-> CIRS
       eraApco eraApco13    terrestrial     ICRS <-> observed
       eraApcs eraApcs13    space           ICRS <-> GCRS
       eraAper eraAper13    terrestrial     update Earth rotation
       eraApio eraApio13    terrestrial     CIRS <-> observed

     Those with names ending in "13" use contemporary ERFA models to
     compute the various ephemerides.  The others accept ephemerides
     supplied by the caller.

     The transformation from ICRS to GCRS covers space motion,
     parallax, light deflection, and aberration.  From GCRS to CIRS
     comprises frame bias and precession-nutation.  From CIRS to
     observed takes account of Earth rotation, polar motion, diurnal
     aberration and parallax (unless subsumed into the ICRS <-> GCRS
     transformation), and atmospheric refraction.

  5) The context structure astrom produced by this function is used by
     eraAtciq* and eraAticq*.


------
apci
------

Given:
TDB as a 2-part...
...Julian Date (Note 1)
Earth barycentric position/velocity (au, au/day)
Earth heliocentric position (au)
CIP X,Y (components of unit vector)
the CIO locator s (radians)

star-independent astrometry parameters:
 PM time interval (SSB, Julian years)
 SSB to observer (vector, au)
 Sun to observer (unit vector)
 distance from Sun to observer (au)
 barycentric observer velocity (vector, c)
 sqrt(1-|v|^2): reciprocal of Lorenz factor
 bias-precession-nutation matrix
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged

ate2 is a Julian Date, apportioned in any
en the two arguments.  For example,
uld be expressed in any of these ways, among

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 method is best matched to the way the
internally and will deliver the optimum
 method and the date & time methods are both
ween resolution and convenience.  For most
 function the choice will not be at all

ad of TDB without any significant impact on

with respect to BCRS axes.
aller does not wish to provide the Earth
O, the function eraApci13 can be used instead
ion.  This computes the required quantities
ctions.
al functions that inserts into the astrom
endent parameters needed for the chain of
mations ICRS <-> GCRS <-> CIRS <-> observed.
s support different classes of observer and
sformation chain:
   observer        transformation
   geocentric      ICRS <-> GCRS
   terrestrial     ICRS <-> CIRS
   terrestrial     ICRS <-> observed
   space           ICRS <-> GCRS
   terrestrial     update Earth rotation
   terrestrial     CIRS <-> observed
ing in "13" use contemporary ERFA models to
ephemerides.  The others accept ephemerides
er.
rom ICRS to GCRS covers space motion,
ection, and aberration.  From GCRS to CIRS
 and precession-nutation.  From CIRS to
nt of Earth rotation, polar motion, diurnal
lax (unless subsumed into the ICRS <-> GCRS
 atmospheric refraction.
e astrom produced by this function is used by
cq*.

try parameters, ICRS-GCRS, geocenter
al-to-intermediate matrix, given X,Y and s
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TDB date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TDB)=2450123.7 could be expressed in any of these ways, among
     others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in cases
     where the loss of several decimal digits of resolution is
     acceptable.  The J2000 method is best matched to the way the
     argument is handled internally and will deliver the optimum
     resolution.  The MJD method and the date & time methods are both
     good compromises between resolution and convenience.  For most
     applications of this function the choice will not be at all
     critical.

     TT can be used instead of TDB without any significant impact on
     accuracy.

  2) All the vectors are with respect to BCRS axes.

  3) In cases where the caller does not wish to provide the Earth
     ephemeris and CIP/CIO, the function eraApci13 can be used instead
     of the present function.  This computes the required quantities
     using other ERFA functions.

  4) This is one of several functions that inserts into the astrom
     structure star-independent parameters needed for the chain of
     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

     The various functions support different classes of observer and
     portions of the transformation chain:

          functions         observer        transformation

       eraApcg eraApcg13    geocentric      ICRS <-> GCRS
       eraApci eraApci13    terrestrial     ICRS <-> CIRS
       eraApco eraApco13    terrestrial     ICRS <-> observed
       eraApcs eraApcs13    space           ICRS <-> GCRS
       eraAper eraAper13    terrestrial     update Earth rotation
       eraApio eraApio13    terrestrial     CIRS <-> observed

     Those with names ending in "13" use contemporary ERFA models to
     compute the various ephemerides.  The others accept ephemerides
     supplied by the caller.

     The transformation from ICRS to GCRS covers space motion,
     parallax, light deflection, and aberration.  From GCRS to CIRS
     comprises frame bias and precession-nutation.  From CIRS to
     observed takes account of Earth rotation, polar motion, diurnal
     aberration and parallax (unless subsumed into the ICRS <-> GCRS
     transformation), and atmospheric refraction.

  5) The context structure astrom produced by this function is used by
     eraAtciq* and eraAticq*.


------
apci13
------

Given:
DB as a 2-part...
..Julian Date (Note 1)

tar-independent astrometry parameters:
 PM time interval (SSB, Julian years)
 SSB to observer (vector, au)
 Sun to observer (unit vector)
 distance from Sun to observer (au)
 barycentric observer velocity (vector, c)
 sqrt(1-|v|^2): reciprocal of Lorenz factor
 bias-precession-nutation matrix
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
quation of the origins (ERA-GST)

ate2 is a Julian Date, apportioned in any
en the two arguments.  For example,
uld be expressed in any of these ways, among

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 method is best matched to the way the
internally and will deliver the optimum
 method and the date & time methods are both
ween resolution and convenience.  For most
 function the choice will not be at all

ad of TDB without any significant impact on

with respect to BCRS axes.
aller wishes to supply his own Earth
O, the function eraApci can be used instead
ion.
al functions that inserts into the astrom
endent parameters needed for the chain of
mations ICRS <-> GCRS <-> CIRS <-> observed.
s support different classes of observer and
sformation chain:
   observer        transformation
   geocentric      ICRS <-> GCRS
   terrestrial     ICRS <-> CIRS
   terrestrial     ICRS <-> observed
   space           ICRS <-> GCRS
   terrestrial     update Earth rotation
   terrestrial     CIRS <-> observed
ing in "13" use contemporary ERFA models to
ephemerides.  The others accept ephemerides
er.
rom ICRS to GCRS covers space motion,
ection, and aberration.  From GCRS to CIRS
 and precession-nutation.  From CIRS to
nt of Earth rotation, polar motion, diurnal
lax (unless subsumed into the ICRS <-> GCRS
 atmospheric refraction.
e astrom produced by this function is used by
cq*.

osition and velocity
al NPB matrix, IAU 2006/2000A
 CIP X,Y coordinates from NPB matrix
 locator s, given X,Y, IAU 2006
try parameters, ICRS-CIRS
n of the origins, given NPB matrix and s
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TDB date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TDB)=2450123.7 could be expressed in any of these ways, among
     others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in cases
     where the loss of several decimal digits of resolution is
     acceptable.  The J2000 method is best matched to the way the
     argument is handled internally and will deliver the optimum
     resolution.  The MJD method and the date & time methods are both
     good compromises between resolution and convenience.  For most
     applications of this function the choice will not be at all
     critical.

     TT can be used instead of TDB without any significant impact on
     accuracy.

  2) All the vectors are with respect to BCRS axes.

  3) In cases where the caller wishes to supply his own Earth
     ephemeris and CIP/CIO, the function eraApci can be used instead
     of the present function.

  4) This is one of several functions that inserts into the astrom
     structure star-independent parameters needed for the chain of
     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

     The various functions support different classes of observer and
     portions of the transformation chain:

          functions         observer        transformation

       eraApcg eraApcg13    geocentric      ICRS <-> GCRS
       eraApci eraApci13    terrestrial     ICRS <-> CIRS
       eraApco eraApco13    terrestrial     ICRS <-> observed
       eraApcs eraApcs13    space           ICRS <-> GCRS
       eraAper eraAper13    terrestrial     update Earth rotation
       eraApio eraApio13    terrestrial     CIRS <-> observed

     Those with names ending in "13" use contemporary ERFA models to
     compute the various ephemerides.  The others accept ephemerides
     supplied by the caller.

     The transformation from ICRS to GCRS covers space motion,
     parallax, light deflection, and aberration.  From GCRS to CIRS
     comprises frame bias and precession-nutation.  From CIRS to
     observed takes account of Earth rotation, polar motion, diurnal
     aberration and parallax (unless subsumed into the ICRS <-> GCRS
     transformation), and atmospheric refraction.

  5) The context structure astrom produced by this function is used by
     eraAtciq* and eraAticq*.


------
apco
------

Given:
TDB as a 2-part...
...Julian Date (Note 1)
Earth barycentric PV (au, au/day, Note 2)
Earth heliocentric P (au, Note 2)
CIP X,Y (components of unit vector)
the CIO locator s (radians)
Earth rotation angle (radians)
longitude (radians, east +ve, Note 3)
latitude (geodetic, radians, Note 3)
height above ellipsoid (m, geodetic, Note 3)
polar motion coordinates (radians, Note 4)
the TIO locator s' (radians, Note 4)
refraction constant A (radians, Note 5)
refraction constant B (radians, Note 5)

star-independent astrometry parameters:
 PM time interval (SSB, Julian years)
 SSB to observer (vector, au)
 Sun to observer (unit vector)
 distance from Sun to observer (au)
 barycentric observer velocity (vector, c)
 sqrt(1-|v|^2): reciprocal of Lorenz factor
 bias-precession-nutation matrix
 longitude + s' (radians)
 polar motion xp wrt local meridian (radians)
 polar motion yp wrt local meridian (radians)
 sine of geodetic latitude
 cosine of geodetic latitude
 magnitude of diurnal aberration vector
 "local" Earth rotation angle (radians)
 refraction constant A (radians)
 refraction constant B (radians)

ate2 is a Julian Date, apportioned in any
en the two arguments.  For example,
uld be expressed in any of these ways, among

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 method is best matched to the way the
internally and will deliver the optimum
 method and the date & time methods are both
ween resolution and convenience.  For most
 function the choice will not be at all

ad of TDB without any significant impact on

and all the astrom vectors, are with respect

rdinates are with respect to the ERFA_WGS84
  TAKE CARE WITH THE LONGITUDE SIGN
gitude required by the present function is
ast-positive, in accordance with geographical

ordinates (in radians) of the Celestial
th respect to the International Terrestrial
e IERS Conventions), measured along the
eg west respectively.  sp is the TIO locator
h positions the Terrestrial Intermediate
r.  For many applications, xp, yp and
be set to zero.
r motion is stored in a form rotated onto the

ants refa and refb are for use in a
3(Z) model, where Z is the observed
ith distance and dZ is the amount of

ake great care with units, as even unlikely
parameters are accepted and processed in
models used.
aller does not wish to provide the Earth
 rotation information and refraction
ion eraApco13 can be used instead of the
his starts from UTC and weather readings etc.
e values using other ERFA functions.
al functions that inserts into the astrom
endent parameters needed for the chain of
mations ICRS <-> GCRS <-> CIRS <-> observed.
s support different classes of observer and
sformation chain:
   observer        transformation
   geocentric      ICRS <-> GCRS
   terrestrial     ICRS <-> CIRS
   terrestrial     ICRS <-> observed
   space           ICRS <-> GCRS
   terrestrial     update Earth rotation
   terrestrial     CIRS <-> observed
ing in "13" use contemporary ERFA models to
ephemerides.  The others accept ephemerides
er.
rom ICRS to GCRS covers space motion,
ection, and aberration.  From GCRS to CIRS
 and precession-nutation.  From CIRS to
nt of Earth rotation, polar motion, diurnal
lax (unless subsumed into the ICRS <-> GCRS
 atmospheric refraction.
e astrom produced by this function is used by
eraAtciq* and eraAticq*.

try parameters: update ERA
al-to-intermediate matrix, given X,Y and s
n/velocity of terrestrial station
 of transpose of r-matrix and pv-vector
try parameters, ICRS-GCRS, space observer
matrix
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TDB date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TDB)=2450123.7 could be expressed in any of these ways, among
     others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in cases
     where the loss of several decimal digits of resolution is
     acceptable.  The J2000 method is best matched to the way the
     argument is handled internally and will deliver the optimum
     resolution.  The MJD method and the date & time methods are both
     good compromises between resolution and convenience.  For most
     applications of this function the choice will not be at all
     critical.

     TT can be used instead of TDB without any significant impact on
     accuracy.

  2) The vectors eb, eh, and all the astrom vectors, are with respect
     to BCRS axes.

  3) The geographical coordinates are with respect to the ERFA_WGS84
     reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN
     CONVENTION:  the longitude required by the present function is
     right-handed, i.e. east-positive, in accordance with geographical
     convention.

  4) xp and yp are the coordinates (in radians) of the Celestial
     Intermediate Pole with respect to the International Terrestrial
     Reference System (see IERS Conventions), measured along the
     meridians 0 and 90 deg west respectively.  sp is the TIO locator
     s', in radians, which positions the Terrestrial Intermediate
     Origin on the equator.  For many applications, xp, yp and
     (especially) sp can be set to zero.

     Internally, the polar motion is stored in a form rotated onto the
     local meridian.

  5) The refraction constants refa and refb are for use in a
     dZ = A*tan(Z)+B*tan^3(Z) model, where Z is the observed
     (i.e. refracted) zenith distance and dZ is the amount of
     refraction.

  6) It is advisable to take great care with units, as even unlikely
     values of the input parameters are accepted and processed in
     accordance with the models used.

  7) In cases where the caller does not wish to provide the Earth
     Ephemeris, the Earth rotation information and refraction
     constants, the function eraApco13 can be used instead of the
     present function.  This starts from UTC and weather readings etc.
     and computes suitable values using other ERFA functions.

  8) This is one of several functions that inserts into the astrom
     structure star-independent parameters needed for the chain of
     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

     The various functions support different classes of observer and
     portions of the transformation chain:

          functions         observer        transformation

       eraApcg eraApcg13    geocentric      ICRS <-> GCRS
       eraApci eraApci13    terrestrial     ICRS <-> CIRS
       eraApco eraApco13    terrestrial     ICRS <-> observed
       eraApcs eraApcs13    space           ICRS <-> GCRS
       eraAper eraAper13    terrestrial     update Earth rotation
       eraApio eraApio13    terrestrial     CIRS <-> observed

     Those with names ending in "13" use contemporary ERFA models to
     compute the various ephemerides.  The others accept ephemerides
     supplied by the caller.

     The transformation from ICRS to GCRS covers space motion,
     parallax, light deflection, and aberration.  From GCRS to CIRS
     comprises frame bias and precession-nutation.  From CIRS to
     observed takes account of Earth rotation, polar motion, diurnal
     aberration and parallax (unless subsumed into the ICRS <-> GCRS
     transformation), and atmospheric refraction.

  9) The context structure astrom produced by this function is used by
     eraAtioq, eraAtoiq, eraAtciq* and eraAticq*.


------
apco13
------

Given:
C as a 2-part...
.quasi Julian Date (Notes 1,2)
1-UTC (seconds, Note 3)
ngitude (radians, east +ve, Note 4)
titude (geodetic, radians, Note 4)
ight above ellipsoid (m, geodetic, Notes 4,6)
lar motion coordinates (radians, Note 5)
essure at the observer (hPa = mB, Note 6)
bient temperature at the observer (deg C)
lative humidity at the observer (range 0-1)
velength (micrometers, Note 7)

ar-independent astrometry parameters:
 PM time interval (SSB, Julian years)
 SSB to observer (vector, au)
 Sun to observer (unit vector)
 distance from Sun to observer (au)
 barycentric observer velocity (vector, c)
 sqrt(1-|v|^2): reciprocal of Lorenz factor
 bias-precession-nutation matrix
 longitude + s' (radians)
 polar motion xp wrt local meridian (radians)
 polar motion yp wrt local meridian (radians)
 sine of geodetic latitude
 cosine of geodetic latitude
 magnitude of diurnal aberration vector
 "local" Earth rotation angle (radians)
 refraction constant A (radians)
 refraction constant B (radians)
uation of the origins (ERA-GST)
e):
atus: +1 = dubious year (Note 2)
       0 = OK
      -1 = unacceptable date

Julian Date (see Note 2), apportioned in any
een the two arguments, for example where utc1
umber and utc2 is the fraction of a day.
unambiguously represent UTC during a leap
al measures are taken.  The convention in the
 that the JD day represents UTC days whether
, 86400 or 86401 SI seconds.
 use the function eraDtf2d to convert from
ime of day into 2-part quasi Julian Date, as
eap-second-ambiguity convention just

"dubious year" flags UTCs that predate the
 time scale or that are too far in the
d.  See eraDat for further details.
d in IERS bulletins.  It increases by exactly
nd of each positive UTC leap second,
 to keep UT1-UTC within +/- 0.9s.  n.b. This
eview, and in the future UT1-UTC may grow
 limit.
ordinates are with respect to the ERFA_WGS84
.  TAKE CARE WITH THE LONGITUDE SIGN:  the
by the present function is east-positive
, in accordance with geographical convention.
,yp can be obtained from IERS bulletins.  The
dinates (in radians) of the Celestial
ith respect to the International Terrestrial
ee IERS Conventions 2003), measured along the
deg west respectively.  For many
d yp can be set to zero.
ar motion is stored in a form rotated onto

bove the ellipsoid of the observing station
nown but phpa, the pressure in hPa (=mB), is
ate estimate of hm can be obtained from the

tsl * log ( phpa / 1013.25 );
proximate sea-level air temperature in K
Quantities, C.W.Allen, 3rd edition, section
 the pressure phpa is not known, it can be
height of the observing station, hm, as

5 * exp ( -hm / ( 29.3 * tsl ) );
 the refraction is nearly proportional to
at an accurate phpa value is important for

cifies the observing wavelength in
ransition from optical to radio is assumed to
eters (about 3000 GHz).
take great care with units, as even unlikely
 parameters are accepted and processed in
 models used.
caller wishes to supply his own Earth
tation information and refraction constants,
o can be used instead of the present function.
ral functions that inserts into the astrom
pendent parameters needed for the chain of
rmations ICRS <-> GCRS <-> CIRS <-> observed.
ns support different classes of observer and
nsformation chain:
   observer        transformation
   geocentric      ICRS <-> GCRS
   terrestrial     ICRS <-> CIRS
   terrestrial     ICRS <-> observed
   space           ICRS <-> GCRS
   terrestrial     update Earth rotation
   terrestrial     CIRS <-> observed
ding in "13" use contemporary ERFA models to
 ephemerides.  The others accept ephemerides
ler.
from ICRS to GCRS covers space motion,
lection, and aberration.  From GCRS to CIRS
s and precession-nutation.  From CIRS to
unt of Earth rotation, polar motion, diurnal
llax (unless subsumed into the ICRS <-> GCRS
d atmospheric refraction.
re astrom produced by this function is used
iq, eraAtciq* and eraAticq*.

TAI
TT
UT1
osition and velocity
al NPB matrix, IAU 2006/2000A
 CIP X,Y coordinates from NPB matrix
 locator s, given X,Y, IAU 2006
otation angle, IAU 2000
 locator s', IERS 2000
ion constants for given ambient conditions
try parameters, ICRS-observed
n of the origins, given NPB matrix and s
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1)  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
      convenient way between the two arguments, for example where utc1
      is the Julian Day Number and utc2 is the fraction of a day.

      However, JD cannot unambiguously represent UTC during a leap
      second unless special measures are taken.  The convention in the
      present function is that the JD day represents UTC days whether
      the length is 86399, 86400 or 86401 SI seconds.

      Applications should use the function eraDtf2d to convert from
      calendar date and time of day into 2-part quasi Julian Date, as
      it implements the leap-second-ambiguity convention just
      described.

  2)  The warning status "dubious year" flags UTCs that predate the
      introduction of the time scale or that are too far in the
      future to be trusted.  See eraDat for further details.

  3)  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
      one second at the end of each positive UTC leap second,
      introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
      practice is under review, and in the future UT1-UTC may grow
      essentially without limit.

  4)  The geographical coordinates are with respect to the ERFA_WGS84
      reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
      longitude required by the present function is east-positive
      (i.e. right-handed), in accordance with geographical convention.

  5)  The polar motion xp,yp can be obtained from IERS bulletins.  The
      values are the coordinates (in radians) of the Celestial
      Intermediate Pole with respect to the International Terrestrial
      Reference System (see IERS Conventions 2003), measured along the
      meridians 0 and 90 deg west respectively.  For many
      applications, xp and yp can be set to zero.

      Internally, the polar motion is stored in a form rotated onto
      the local meridian.

  6)  If hm, the height above the ellipsoid of the observing station
      in meters, is not known but phpa, the pressure in hPa (=mB), is
      available, an adequate estimate of hm can be obtained from the
      expression

            hm = -29.3 * tsl * log ( phpa / 1013.25 );

      where tsl is the approximate sea-level air temperature in K
      (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
      52).  Similarly, if the pressure phpa is not known, it can be
      estimated from the height of the observing station, hm, as
      follows:

            phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

      Note, however, that the refraction is nearly proportional to
      the pressure and that an accurate phpa value is important for
      precise work.

  7)  The argument wl specifies the observing wavelength in
      micrometers.  The transition from optical to radio is assumed to
      occur at 100 micrometers (about 3000 GHz).

  8)  It is advisable to take great care with units, as even unlikely
      values of the input parameters are accepted and processed in
      accordance with the models used.

  9)  In cases where the caller wishes to supply his own Earth
      ephemeris, Earth rotation information and refraction constants,
      the function eraApco can be used instead of the present function.

  10) This is one of several functions that inserts into the astrom
      structure star-independent parameters needed for the chain of
      astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

      The various functions support different classes of observer and
      portions of the transformation chain:

          functions         observer        transformation

       eraApcg eraApcg13    geocentric      ICRS <-> GCRS
       eraApci eraApci13    terrestrial     ICRS <-> CIRS
       eraApco eraApco13    terrestrial     ICRS <-> observed
       eraApcs eraApcs13    space           ICRS <-> GCRS
       eraAper eraAper13    terrestrial     update Earth rotation
       eraApio eraApio13    terrestrial     CIRS <-> observed

      Those with names ending in "13" use contemporary ERFA models to
      compute the various ephemerides.  The others accept ephemerides
      supplied by the caller.

      The transformation from ICRS to GCRS covers space motion,
      parallax, light deflection, and aberration.  From GCRS to CIRS
      comprises frame bias and precession-nutation.  From CIRS to
      observed takes account of Earth rotation, polar motion, diurnal
      aberration and parallax (unless subsumed into the ICRS <-> GCRS
      transformation), and atmospheric refraction.

  11) The context structure astrom produced by this function is used
      by eraAtioq, eraAtoiq, eraAtciq* and eraAticq*.


------
apcs
------

Given:
TDB as a 2-part...
...Julian Date (Note 1)
observer's geocentric pos/vel (m, m/s)
Earth barycentric PV (au, au/day)
Earth heliocentric P (au)

star-independent astrometry parameters:
 PM time interval (SSB, Julian years)
 SSB to observer (vector, au)
 Sun to observer (unit vector)
 distance from Sun to observer (au)
 barycentric observer velocity (vector, c)
 sqrt(1-|v|^2): reciprocal of Lorenz factor
 bias-precession-nutation matrix
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged

ate2 is a Julian Date, apportioned in any
en the two arguments.  For example,
uld be expressed in any of these ways, among

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 method is best matched to the way the
internally and will deliver the optimum
 method and the date & time methods are both
ween resolution and convenience.  For most
 function the choice will not be at all

ad of TDB without any significant impact on

with respect to BCRS axes.
rguments for (i) the observer's geocentric
y and (ii) the Earth ephemeris is done for
eocentric, terrestrial and Earth orbit cases.
cations it maybe more convenient to specify
tion and velocity and to supply the
and velocity information directly instead of
Earth.  However, note the different units:
ocentric vectors, au and au/day for the
ycentric vectors.
aller does not wish to provide the Earth
ion eraApcs13 can be used instead of the
his computes the Earth ephemeris using the
00.
al functions that inserts into the astrom
endent parameters needed for the chain of
mations ICRS <-> GCRS <-> CIRS <-> observed.
s support different classes of observer and
sformation chain:
   observer        transformation
   geocentric      ICRS <-> GCRS
   terrestrial     ICRS <-> CIRS
   terrestrial     ICRS <-> observed
   space           ICRS <-> GCRS
   terrestrial     update Earth rotation
   terrestrial     CIRS <-> observed
ing in "13" use contemporary ERFA models to
ephemerides.  The others accept ephemerides
er.
rom ICRS to GCRS covers space motion,
ection, and aberration.  From GCRS to CIRS
 and precession-nutation.  From CIRS to
nt of Earth rotation, polar motion, diurnal
lax (unless subsumed into the ICRS <-> GCRS
 atmospheric refraction.
e astrom produced by this function is used by
cq*.

vector
 of p-vector
se p-vector into modulus and direction
ize r-matrix to identity
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TDB date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TDB)=2450123.7 could be expressed in any of these ways, among
     others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in cases
     where the loss of several decimal digits of resolution is
     acceptable.  The J2000 method is best matched to the way the
     argument is handled internally and will deliver the optimum
     resolution.  The MJD method and the date & time methods are both
     good compromises between resolution and convenience.  For most
     applications of this function the choice will not be at all
     critical.

     TT can be used instead of TDB without any significant impact on
     accuracy.

  2) All the vectors are with respect to BCRS axes.

  3) Providing separate arguments for (i) the observer's geocentric
     position and velocity and (ii) the Earth ephemeris is done for
     convenience in the geocentric, terrestrial and Earth orbit cases.
     For deep space applications it maybe more convenient to specify
     zero geocentric position and velocity and to supply the
     observer's position and velocity information directly instead of
     with respect to the Earth.  However, note the different units:
     m and m/s for the geocentric vectors, au and au/day for the
     heliocentric and barycentric vectors.

  4) In cases where the caller does not wish to provide the Earth
     ephemeris, the function eraApcs13 can be used instead of the
     present function.  This computes the Earth ephemeris using the
     ERFA function eraEpv00.

  5) This is one of several functions that inserts into the astrom
     structure star-independent parameters needed for the chain of
     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

     The various functions support different classes of observer and
     portions of the transformation chain:

          functions         observer        transformation

       eraApcg eraApcg13    geocentric      ICRS <-> GCRS
       eraApci eraApci13    terrestrial     ICRS <-> CIRS
       eraApco eraApco13    terrestrial     ICRS <-> observed
       eraApcs eraApcs13    space           ICRS <-> GCRS
       eraAper eraAper13    terrestrial     update Earth rotation
       eraApio eraApio13    terrestrial     CIRS <-> observed

     Those with names ending in "13" use contemporary ERFA models to
     compute the various ephemerides.  The others accept ephemerides
     supplied by the caller.

     The transformation from ICRS to GCRS covers space motion,
     parallax, light deflection, and aberration.  From GCRS to CIRS
     comprises frame bias and precession-nutation.  From CIRS to
     observed takes account of Earth rotation, polar motion, diurnal
     aberration and parallax (unless subsumed into the ICRS <-> GCRS
     transformation), and atmospheric refraction.

  6) The context structure astrom produced by this function is used by
     eraAtciq* and eraAticq*.


------
apcs13
------

Given:
TDB as a 2-part...
...Julian Date (Note 1)
observer's geocentric pos/vel (Note 3)

star-independent astrometry parameters:
 PM time interval (SSB, Julian years)
 SSB to observer (vector, au)
 Sun to observer (unit vector)
 distance from Sun to observer (au)
 barycentric observer velocity (vector, c)
 sqrt(1-|v|^2): reciprocal of Lorenz factor
 bias-precession-nutation matrix
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged

ate2 is a Julian Date, apportioned in any
en the two arguments.  For example,
uld be expressed in any of these ways, among

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 method is best matched to the way the
internally and will deliver the optimum
 method and the date & time methods are both
ween resolution and convenience.  For most
 function the choice will not be at all

ad of TDB without any significant impact on

with respect to BCRS axes.
ion and velocity pv are geocentric but with
, and in units of m and m/s.  No assumptions
mity to the Earth, and the function can be
applications as well as Earth orbit and

aller wishes to supply his own Earth
ion eraApcs can be used instead of the present

al functions that inserts into the astrom
endent parameters needed for the chain of
mations ICRS <-> GCRS <-> CIRS <-> observed.
s support different classes of observer and
sformation chain:
   observer        transformation
   geocentric      ICRS <-> GCRS
   terrestrial     ICRS <-> CIRS
   terrestrial     ICRS <-> observed
   space           ICRS <-> GCRS
   terrestrial     update Earth rotation
   terrestrial     CIRS <-> observed
ing in "13" use contemporary ERFA models to
ephemerides.  The others accept ephemerides
er.
rom ICRS to GCRS covers space motion,
ection, and aberration.  From GCRS to CIRS
 and precession-nutation.  From CIRS to
nt of Earth rotation, polar motion, diurnal
lax (unless subsumed into the ICRS <-> GCRS
 atmospheric refraction.
e astrom produced by this function is used by
cq*.

osition and velocity
try parameters, ICRS-GCRS, space observer
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TDB date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TDB)=2450123.7 could be expressed in any of these ways, among
     others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in cases
     where the loss of several decimal digits of resolution is
     acceptable.  The J2000 method is best matched to the way the
     argument is handled internally and will deliver the optimum
     resolution.  The MJD method and the date & time methods are both
     good compromises between resolution and convenience.  For most
     applications of this function the choice will not be at all
     critical.

     TT can be used instead of TDB without any significant impact on
     accuracy.

  2) All the vectors are with respect to BCRS axes.

  3) The observer's position and velocity pv are geocentric but with
     respect to BCRS axes, and in units of m and m/s.  No assumptions
     are made about proximity to the Earth, and the function can be
     used for deep space applications as well as Earth orbit and
     terrestrial.

  4) In cases where the caller wishes to supply his own Earth
     ephemeris, the function eraApcs can be used instead of the present
     function.

  5) This is one of several functions that inserts into the astrom
     structure star-independent parameters needed for the chain of
     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

     The various functions support different classes of observer and
     portions of the transformation chain:

          functions         observer        transformation

       eraApcg eraApcg13    geocentric      ICRS <-> GCRS
       eraApci eraApci13    terrestrial     ICRS <-> CIRS
       eraApco eraApco13    terrestrial     ICRS <-> observed
       eraApcs eraApcs13    space           ICRS <-> GCRS
       eraAper eraAper13    terrestrial     update Earth rotation
       eraApio eraApio13    terrestrial     CIRS <-> observed

     Those with names ending in "13" use contemporary ERFA models to
     compute the various ephemerides.  The others accept ephemerides
     supplied by the caller.

     The transformation from ICRS to GCRS covers space motion,
     parallax, light deflection, and aberration.  From GCRS to CIRS
     comprises frame bias and precession-nutation.  From CIRS to
     observed takes account of Earth rotation, polar motion, diurnal
     aberration and parallax (unless subsumed into the ICRS <-> GCRS
     transformation), and atmospheric refraction.

  6) The context structure astrom produced by this function is used by
     eraAtciq* and eraAticq*.


------
aper
------

Given:
Earth rotation angle (radians, Note 2)
star-independent astrometry parameters:
 not used
 not used
 not used
 not used
 not used
 not used
 not used
 longitude + s' (radians)
 not used
 not used
 not used
 not used
 not used
 not used
 not used
 not used

star-independent astrometry parameters:
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 "local" Earth rotation angle (radians)
 unchanged
 unchanged

 to enable sidereal-tracking applications to
putation of the bulk of the astrometry
e Earth rotation is updated.
d as equinox based positions, such as
 apparent (RA,Dec), the supplied theta can be
idereal time rather than Earth rotation

13 can be used instead of the present
 from UT1 rather than ERA itself.
al functions that inserts into the astrom
endent parameters needed for the chain of
mations ICRS <-> GCRS <-> CIRS <-> observed.
s support different classes of observer and
sformation chain:
   observer        transformation
   geocentric      ICRS <-> GCRS
   terrestrial     ICRS <-> CIRS
   terrestrial     ICRS <-> observed
   space           ICRS <-> GCRS
   terrestrial     update Earth rotation
   terrestrial     CIRS <-> observed
ing in "13" use contemporary ERFA models to
ephemerides.  The others accept ephemerides
er.
rom ICRS to GCRS covers space motion,
ection, and aberration.  From GCRS to CIRS
 and precession-nutation.  From CIRS to
nt of Earth rotation, polar motion, diurnal
lax (unless subsumed into the ICRS <-> GCRS
 atmospheric refraction.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
aper13
------

Given:
UT1 as a 2-part...
...Julian Date (Note 1)
star-independent astrometry parameters:
 not used
 not used
 not used
 not used
 not used
 not used
 not used
 longitude + s' (radians)
 not used
 not used
 not used
 not used
 not used
 not used
 not used
 not used

star-independent astrometry parameters:
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 "local" Earth rotation angle (radians)
 unchanged
 unchanged

ot UTC) ut11+ut12 is a Julian Date,
onvenient way between the arguments ut11 and
JD(UT1)=2450123.7 could be expressed in any
 others:
  ut12
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 and MJD methods are good compromises
nd convenience.  The date & time method is
algorithm used:  maximum precision is
t11 argument is for 0hrs UT1 on the day in
2 argument lies in the range 0 to 1, or vice

 to provide the Earth rotation angle itself,
 can be used instead.  One use of this
titute Greenwich apparent sidereal time and
quinox based transformations directly.
al functions that inserts into the astrom
endent parameters needed for the chain of
mations ICRS <-> GCRS <-> CIRS <-> observed.
s support different classes of observer and
sformation chain:
   observer        transformation
   geocentric      ICRS <-> GCRS
   terrestrial     ICRS <-> CIRS
   terrestrial     ICRS <-> observed
   space           ICRS <-> GCRS
   terrestrial     update Earth rotation
   terrestrial     CIRS <-> observed
ing in "13" use contemporary ERFA models to
ephemerides.  The others accept ephemerides
er.
rom ICRS to GCRS covers space motion,
ection, and aberration.  From GCRS to CIRS
 and precession-nutation.  From CIRS to
nt of Earth rotation, polar motion, diurnal
lax (unless subsumed into the ICRS <-> GCRS
 atmospheric refraction.

try parameters: update ERA
otation angle, IAU 2000
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The UT1 date (n.b. not UTC) ut11+ut12 is a Julian Date,
     apportioned in any convenient way between the arguments ut11 and
     ut12.  For example, JD(UT1)=2450123.7 could be expressed in any
     of these ways, among others:

            ut11           ut12

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in cases
     where the loss of several decimal digits of resolution is
     acceptable.  The J2000 and MJD methods are good compromises
     between resolution and convenience.  The date & time method is
     best matched to the algorithm used:  maximum precision is
     delivered when the ut11 argument is for 0hrs UT1 on the day in
     question and the ut12 argument lies in the range 0 to 1, or vice
     versa.

  2) If the caller wishes to provide the Earth rotation angle itself,
     the function eraAper can be used instead.  One use of this
     technique is to substitute Greenwich apparent sidereal time and
     thereby to support equinox based transformations directly.

  3) This is one of several functions that inserts into the astrom
     structure star-independent parameters needed for the chain of
     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

     The various functions support different classes of observer and
     portions of the transformation chain:

          functions         observer        transformation

       eraApcg eraApcg13    geocentric      ICRS <-> GCRS
       eraApci eraApci13    terrestrial     ICRS <-> CIRS
       eraApco eraApco13    terrestrial     ICRS <-> observed
       eraApcs eraApcs13    space           ICRS <-> GCRS
       eraAper eraAper13    terrestrial     update Earth rotation
       eraApio eraApio13    terrestrial     CIRS <-> observed

     Those with names ending in "13" use contemporary ERFA models to
     compute the various ephemerides.  The others accept ephemerides
     supplied by the caller.

     The transformation from ICRS to GCRS covers space motion,
     parallax, light deflection, and aberration.  From GCRS to CIRS
     comprises frame bias and precession-nutation.  From CIRS to
     observed takes account of Earth rotation, polar motion, diurnal
     aberration and parallax (unless subsumed into the ICRS <-> GCRS
     transformation), and atmospheric refraction.


------
apio
------

Given:
he TIO locator s' (radians, Note 1)
arth rotation angle (radians)
ongitude (radians, east +ve, Note 2)
eodetic latitude (radians, Note 2)
eight above ellipsoid (m, geodetic Note 2)
olar motion coordinates (radians, Note 3)
efraction constant A (radians, Note 4)
efraction constant B (radians, Note 4)

tar-independent astrometry parameters:
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 longitude + s' (radians)
 polar motion xp wrt local meridian (radians)
 polar motion yp wrt local meridian (radians)
 sine of geodetic latitude
 cosine of geodetic latitude
 magnitude of diurnal aberration vector
 "local" Earth rotation angle (radians)
 refraction constant A (radians)
 refraction constant B (radians)

s', is a tiny quantity needed only by the
tions.  It can either be set to zero or
ERFA function eraSp00.
rdinates are with respect to the ERFA_WGS84
  TAKE CARE WITH THE LONGITUDE SIGN:  the
y the present function is east-positive
 in accordance with geographical convention.
yp can be obtained from IERS bulletins.  The
inates (in radians) of the Celestial
th respect to the International Terrestrial
e IERS Conventions 2003), measured along the
eg west respectively.  For many applications,
 to zero.
r motion is stored in a form rotated onto the

ants refa and refb are for use in a
3(Z) model, where Z is the observed
ith distance and dZ is the amount of

ake great care with units, as even unlikely
parameters are accepted and processed in
models used.
aller does not wish to provide the Earth
 and refraction constants, the function
d instead of the present function.  This
weather readings etc. and computes suitable
RFA functions.
al functions that inserts into the astrom
endent parameters needed for the chain of
mations ICRS <-> GCRS <-> CIRS <-> observed.
s support different classes of observer and
sformation chain:
   observer        transformation
   geocentric      ICRS <-> GCRS
   terrestrial     ICRS <-> CIRS
   terrestrial     ICRS <-> observed
   space           ICRS <-> GCRS
   terrestrial     update Earth rotation
   terrestrial     CIRS <-> observed
ing in "13" use contemporary ERFA models to
ephemerides.  The others accept ephemerides
er.
rom ICRS to GCRS covers space motion,
ection, and aberration.  From GCRS to CIRS
 and precession-nutation.  From CIRS to
nt of Earth rotation, polar motion, diurnal
lax (unless subsumed into the ICRS <-> GCRS
 atmospheric refraction.
e astrom produced by this function is used by
q.

n/velocity of terrestrial station
try parameters: update ERA
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) sp, the TIO locator s', is a tiny quantity needed only by the
     most precise applications.  It can either be set to zero or
     predicted using the ERFA function eraSp00.

  2) The geographical coordinates are with respect to the ERFA_WGS84
     reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
     longitude required by the present function is east-positive
     (i.e. right-handed), in accordance with geographical convention.

  3) The polar motion xp,yp can be obtained from IERS bulletins.  The
     values are the coordinates (in radians) of the Celestial
     Intermediate Pole with respect to the International Terrestrial
     Reference System (see IERS Conventions 2003), measured along the
     meridians 0 and 90 deg west respectively.  For many applications,
     xp and yp can be set to zero.

     Internally, the polar motion is stored in a form rotated onto the
     local meridian.

  4) The refraction constants refa and refb are for use in a
     dZ = A*tan(Z)+B*tan^3(Z) model, where Z is the observed
     (i.e. refracted) zenith distance and dZ is the amount of
     refraction.

  5) It is advisable to take great care with units, as even unlikely
     values of the input parameters are accepted and processed in
     accordance with the models used.

  6) In cases where the caller does not wish to provide the Earth
     rotation information and refraction constants, the function
     eraApio13 can be used instead of the present function.  This
     starts from UTC and weather readings etc. and computes suitable
     values using other ERFA functions.

  7) This is one of several functions that inserts into the astrom
     structure star-independent parameters needed for the chain of
     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

     The various functions support different classes of observer and
     portions of the transformation chain:

          functions         observer        transformation

       eraApcg eraApcg13    geocentric      ICRS <-> GCRS
       eraApci eraApci13    terrestrial     ICRS <-> CIRS
       eraApco eraApco13    terrestrial     ICRS <-> observed
       eraApcs eraApcs13    space           ICRS <-> GCRS
       eraAper eraAper13    terrestrial     update Earth rotation
       eraApio eraApio13    terrestrial     CIRS <-> observed

     Those with names ending in "13" use contemporary ERFA models to
     compute the various ephemerides.  The others accept ephemerides
     supplied by the caller.

     The transformation from ICRS to GCRS covers space motion,
     parallax, light deflection, and aberration.  From GCRS to CIRS
     comprises frame bias and precession-nutation.  From CIRS to
     observed takes account of Earth rotation, polar motion, diurnal
     aberration and parallax (unless subsumed into the ICRS <-> GCRS
     transformation), and atmospheric refraction.

  8) The context structure astrom produced by this function is used by
     eraAtioq and eraAtoiq.


------
apio13
------

Given:
TC as a 2-part...
..quasi Julian Date (Notes 1,2)
T1-UTC (seconds)
ongitude (radians, east +ve, Note 3)
eodetic latitude (radians, Note 3)
eight above ellipsoid (m, geodetic Notes 4,6)
olar motion coordinates (radians, Note 5)
ressure at the observer (hPa = mB, Note 6)
mbient temperature at the observer (deg C)
elative humidity at the observer (range 0-1)
avelength (micrometers, Note 7)

tar-independent astrometry parameters:
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 unchanged
 longitude + s' (radians)
 polar motion xp wrt local meridian (radians)
 polar motion yp wrt local meridian (radians)
 sine of geodetic latitude
 cosine of geodetic latitude
 magnitude of diurnal aberration vector
 "local" Earth rotation angle (radians)
 refraction constant A (radians)
 refraction constant B (radians)
e):
tatus: +1 = dubious year (Note 2)
        0 = OK
       -1 = unacceptable date

Julian Date (see Note 2), apportioned in any
een the two arguments, for example where utc1
umber and utc2 is the fraction of a day.
unambiguously represent UTC during a leap
al measures are taken.  The convention in the
 that the JD day represents UTC days whether
, 86400 or 86401 SI seconds.
 use the function eraDtf2d to convert from
ime of day into 2-part quasi Julian Date, as
eap-second-ambiguity convention just

"dubious year" flags UTCs that predate the
 time scale or that are too far in the future
 eraDat for further details.
d in IERS bulletins.  It increases by exactly
nd of each positive UTC leap second,
 to keep UT1-UTC within +/- 0.9s.  n.b. This
eview, and in the future UT1-UTC may grow
 limit.
ordinates are with respect to the ERFA_WGS84
.  TAKE CARE WITH THE LONGITUDE SIGN:  the
by the present function is east-positive
, in accordance with geographical convention.
,yp can be obtained from IERS bulletins.  The
dinates (in radians) of the Celestial
ith respect to the International Terrestrial
ee IERS Conventions 2003), measured along the
deg west respectively.  For many applications,
t to zero.
ar motion is stored in a form rotated onto

bove the ellipsoid of the observing station
nown but phpa, the pressure in hPa (=mB), is
ate estimate of hm can be obtained from the

tsl * log ( phpa / 1013.25 );
proximate sea-level air temperature in K
Quantities, C.W.Allen, 3rd edition, section
 the pressure phpa is not known, it can be
height of the observing station, hm, as

5 * exp ( -hm / ( 29.3 * tsl ) );
 the refraction is nearly proportional to the
n accurate phpa value is important for

cifies the observing wavelength in
ransition from optical to radio is assumed to
eters (about 3000 GHz).
take great care with units, as even unlikely
 parameters are accepted and processed in
 models used.
caller wishes to supply his own Earth
n and refraction constants, the function
instead of the present function.
ral functions that inserts into the astrom
pendent parameters needed for the chain of
rmations ICRS <-> GCRS <-> CIRS <-> observed.
ns support different classes of observer and
nsformation chain:
   observer        transformation
   geocentric      ICRS <-> GCRS
   terrestrial     ICRS <-> CIRS
   terrestrial     ICRS <-> observed
   space           ICRS <-> GCRS
   terrestrial     update Earth rotation
   terrestrial     CIRS <-> observed
ding in "13" use contemporary ERFA models to
 ephemerides.  The others accept ephemerides
ler.
from ICRS to GCRS covers space motion,
lection, and aberration.  From GCRS to CIRS
s and precession-nutation.  From CIRS to
unt of Earth rotation, polar motion, diurnal
llax (unless subsumed into the ICRS <-> GCRS
d atmospheric refraction.
re astrom produced by this function is used
Atoiq.

TAI
TT
UT1
 locator s', IERS 2000
otation angle, IAU 2000
ion constants for given ambient conditions
try parameters, CIRS-observed
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1)  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
      convenient way between the two arguments, for example where utc1
      is the Julian Day Number and utc2 is the fraction of a day.

      However, JD cannot unambiguously represent UTC during a leap
      second unless special measures are taken.  The convention in the
      present function is that the JD day represents UTC days whether
      the length is 86399, 86400 or 86401 SI seconds.

      Applications should use the function eraDtf2d to convert from
      calendar date and time of day into 2-part quasi Julian Date, as
      it implements the leap-second-ambiguity convention just
      described.

  2)  The warning status "dubious year" flags UTCs that predate the
      introduction of the time scale or that are too far in the future
      to be trusted.  See eraDat for further details.

  3)  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
      one second at the end of each positive UTC leap second,
      introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
      practice is under review, and in the future UT1-UTC may grow
      essentially without limit.

  4)  The geographical coordinates are with respect to the ERFA_WGS84
      reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
      longitude required by the present function is east-positive
      (i.e. right-handed), in accordance with geographical convention.

  5)  The polar motion xp,yp can be obtained from IERS bulletins.  The
      values are the coordinates (in radians) of the Celestial
      Intermediate Pole with respect to the International Terrestrial
      Reference System (see IERS Conventions 2003), measured along the
      meridians 0 and 90 deg west respectively.  For many applications,
      xp and yp can be set to zero.

      Internally, the polar motion is stored in a form rotated onto
      the local meridian.

  6)  If hm, the height above the ellipsoid of the observing station
      in meters, is not known but phpa, the pressure in hPa (=mB), is
      available, an adequate estimate of hm can be obtained from the
      expression

            hm = -29.3 * tsl * log ( phpa / 1013.25 );

      where tsl is the approximate sea-level air temperature in K
      (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
      52).  Similarly, if the pressure phpa is not known, it can be
      estimated from the height of the observing station, hm, as
      follows:

            phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

      Note, however, that the refraction is nearly proportional to the
      pressure and that an accurate phpa value is important for
      precise work.

  7)  The argument wl specifies the observing wavelength in
      micrometers.  The transition from optical to radio is assumed to
      occur at 100 micrometers (about 3000 GHz).

  8)  It is advisable to take great care with units, as even unlikely
      values of the input parameters are accepted and processed in
      accordance with the models used.

  9)  In cases where the caller wishes to supply his own Earth
      rotation information and refraction constants, the function
      eraApc can be used instead of the present function.

  10) This is one of several functions that inserts into the astrom
      structure star-independent parameters needed for the chain of
      astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

      The various functions support different classes of observer and
      portions of the transformation chain:

          functions         observer        transformation

       eraApcg eraApcg13    geocentric      ICRS <-> GCRS
       eraApci eraApci13    terrestrial     ICRS <-> CIRS
       eraApco eraApco13    terrestrial     ICRS <-> observed
       eraApcs eraApcs13    space           ICRS <-> GCRS
       eraAper eraAper13    terrestrial     update Earth rotation
       eraApio eraApio13    terrestrial     CIRS <-> observed

      Those with names ending in "13" use contemporary ERFA models to
      compute the various ephemerides.  The others accept ephemerides
      supplied by the caller.

      The transformation from ICRS to GCRS covers space motion,
      parallax, light deflection, and aberration.  From GCRS to CIRS
      comprises frame bias and precession-nutation.  From CIRS to
      observed takes account of Earth rotation, polar motion, diurnal
      aberration and parallax (unless subsumed into the ICRS <-> GCRS
      transformation), and atmospheric refraction.

  11) The context structure astrom produced by this function is used
      by eraAtioq and eraAtoiq.


------
atci13
------

Given:
 right ascension at J2000.0 (radians, Note 1)
 declination at J2000.0 (radians, Note 1)
roper motion (radians/year; Note 2)
proper motion (radians/year)
llax (arcsec)
al velocity (km/s, +ve if receding)
as a 2-part...
ulian Date (Note 3)

 geocentric RA,Dec (radians)
tion of the origins (ERA-GST, Note 5)

ch other than J2000.0 (for example from the
hich has an epoch of J1991.25) will require a
eraPmsafe before use.
 RA is dRA/dt rather than cos(Dec)*dRA/dt.
ate2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ould be expressed in any of these ways, among

  date2
     0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 method is best matched to the way the
internally and will deliver the optimum
 method and the date & time methods are both
ween resolution and convenience.  For most
 function the choice will not be at all

ad of TDB without any significant impact on

cy is better than 1 milliarcsecond, limited
sion-nutation model that is used, namely
y close to solar system bodies, additional
ral milliarcseconds can occur because of
ection;  however, the Sun's contribution is
to first order.  The accuracy limitations of
aEpv00 (used to compute Earth position and
bute aberration errors of up to
Light deflection at the Sun's limb is
 mas level.
ation to (equinox based) apparent place be
 (CIO based) intermediate place, subtract the
ins from the returned right ascension:
raAnp function can then be applied, as
e result in the conventional 0-2pi range.)

try parameters, ICRS-CIRS, 2013
CRS to CIRS
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Star data for an epoch other than J2000.0 (for example from the
     Hipparcos catalog, which has an epoch of J1991.25) will require a
     preliminary call to eraPmsafe before use.

  2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

  3) The TDB date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TDB)=2450123.8g could be expressed in any of these ways, among
     others:

            date1          date2

         2450123.8g           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in cases
     where the loss of several decimal digits of resolution is
     acceptable.  The J2000 method is best matched to the way the
     argument is handled internally and will deliver the optimum
     resolution.  The MJD method and the date & time methods are both
     good compromises between resolution and convenience.  For most
     applications of this function the choice will not be at all
     critical.

     TT can be used instead of TDB without any significant impact on
     accuracy.

  4) The available accuracy is better than 1 milliarcsecond, limited
     mainly by the precession-nutation model that is used, namely
     IAU 2000A/2006.  Very close to solar system bodies, additional
     errors of up to several milliarcseconds can occur because of
     unmodeled light deflection;  however, the Sun's contribution is
     taken into account, to first order.  The accuracy limitations of
     the ERFA function eraEpv00 (used to compute Earth position and
     velocity) can contribute aberration errors of up to
     5 microarcseconds.  Light deflection at the Sun's limb is
     uncertain at the 0.4 mas level.

  5) Should the transformation to (equinox based) apparent place be
     required rather than (CIO based) intermediate place, subtract the
     equation of the origins from the returned right ascension:
     RA = RI - EO. (The eraAnp function can then be applied, as
     required, to keep the result in the conventional 0-2pi range.)


------
atciq
------

Given:
RS RA,Dec at J2000.0 (radians)
 proper motion (radians/year; Note 3)
c proper motion (radians/year)
rallax (arcsec)
dial velocity (km/s, +ve if receding)
ar-independent astrometry parameters:
 PM time interval (SSB, Julian years)
 SSB to observer (vector, au)
 Sun to observer (unit vector)
 distance from Sun to observer (au)
 barycentric observer velocity (vector, c)
 sqrt(1-|v|^2): reciprocal of Lorenz factor
 bias-precession-nutation matrix
 longitude + s' (radians)
 polar motion xp wrt local meridian (radians)
 polar motion yp wrt local meridian (radians)
 sine of geodetic latitude
 cosine of geodetic latitude
 magnitude of diurnal aberration vector
 "local" Earth rotation angle (radians)
 refraction constant A (radians)
 refraction constant B (radians)

RS RA,Dec (radians)

with respect to BCRS axes.
ch other than J2000.0 (for example from the
hich has an epoch of J1991.25) will require a
eraPmsafe before use.
 RA is dRA/dt rather than cos(Dec)*dRA/dt.

motion and parallax
eflection by the Sun
 aberration
 of r-matrix and pv-vector
r to spherical
ze angle into range 0 to 2pi
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) All the vectors are with respect to BCRS axes.

  2) Star data for an epoch other than J2000.0 (for example from the
     Hipparcos catalog, which has an epoch of J1991.25) will require a
     preliminary call to eraPmsafe before use.

  3) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.


------
atciqn
------

Given:
ICRS RA,Dec at J2000.0 (radians)
RA proper motion (radians/year; Note 3)
Dec proper motion (radians/year)
parallax (arcsec)
radial velocity (km/s, +ve if receding)
star-independent astrometry parameters:
 PM time interval (SSB, Julian years)
 SSB to observer (vector, au)
 Sun to observer (unit vector)
 distance from Sun to observer (au)
 barycentric observer velocity (vector, c)
 sqrt(1-|v|^2): reciprocal of Lorenz factor
 bias-precession-nutation matrix
 longitude + s' (radians)
 polar motion xp wrt local meridian (radians)
 polar motion yp wrt local meridian (radians)
 sine of geodetic latitude
 cosine of geodetic latitude
 magnitude of diurnal aberration vector
 "local" Earth rotation angle (radians)
 refraction constant A (radians)
 refraction constant B (radians)
 number of bodies (Note 3)
data for each of the n bodies (Notes 3,4):
  mass of the body (solar masses, Note 5)
  deflection limiter (Note 6)
  barycentric PV of the body (au, au/day)

RS RA,Dec (radians)

ch other than J2000.0 (for example from the
hich has an epoch of J1991.25) will require a
eraPmsafe before use.
 RA is dRA/dt rather than cos(Dec)*dRA/dt.
s n entries, one for each body to be
0, no gravitational light deflection will be
r the Sun.
include an entry for the Sun as well as for
body to be taken into account.  The entries
er in which the light passes the body.
b struct for body i, the mass parameter
ired, be adjusted in order to allow for such
e field.
er parameter b[i].dl is phi^2/2, where phi is
on (in radians) between star and body at
plied.  As phi shrinks below the chosen
ction is artificially reduced, reaching zero
le values suitable for a terrestrial
ith masses, are as follows:
m        b[i].dl
         6e-6
5435     3e-9
8574     3e-10
dation of the contents of the b array is
ed masses must be greater than zero, the
y vectors must be right, and the deflection
 zero.

motion and parallax
eflection by n bodies
 aberration
 of r-matrix and pv-vector
r to spherical
ze angle into range 0 to 2pi
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Star data for an epoch other than J2000.0 (for example from the
     Hipparcos catalog, which has an epoch of J1991.25) will require a
     preliminary call to eraPmsafe before use.

  2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

  3) The struct b contains n entries, one for each body to be
     considered.  If n = 0, no gravitational light deflection will be
     applied, not even for the Sun.

  4) The struct b should include an entry for the Sun as well as for
     any planet or other body to be taken into account.  The entries
     should be in the order in which the light passes the body.

  5) In the entry in the b struct for body i, the mass parameter
     b[i].bm can, as required, be adjusted in order to allow for such
     effects as quadrupole field.

  6) The deflection limiter parameter b[i].dl is phi^2/2, where phi is
     the angular separation (in radians) between star and body at
     which limiting is applied.  As phi shrinks below the chosen
     threshold, the deflection is artificially reduced, reaching zero
     for phi = 0.   Example values suitable for a terrestrial
     observer, together with masses, are as follows:

        body i     b[i].bm        b[i].dl

        Sun        1.0            6e-6
        Jupiter    0.00095435     3e-9
        Saturn     0.00028574     3e-10

  7) For efficiency, validation of the contents of the b array is
     omitted.  The supplied masses must be greater than zero, the
     position and velocity vectors must be right, and the deflection
     limiter greater than zero.


------
atciqz
------

Given:
RS astrometric RA,Dec (radians)
ar-independent astrometry parameters:
 PM time interval (SSB, Julian years)
 SSB to observer (vector, au)
 Sun to observer (unit vector)
 distance from Sun to observer (au)
 barycentric observer velocity (vector, c)
 sqrt(1-|v|^2): reciprocal of Lorenz factor
 bias-precession-nutation matrix
 longitude + s' (radians)
 polar motion xp wrt local meridian (radians)
 polar motion yp wrt local meridian (radians)
 sine of geodetic latitude
 cosine of geodetic latitude
 magnitude of diurnal aberration vector
 "local" Earth rotation angle (radians)
 refraction constant A (radians)
 refraction constant B (radians)

RS RA,Dec (radians)

with respect to BCRS axes.

nn, P. K. (eds), Explanatory Supplement to
anac, 3rd ed., University Science Books

"A practical relativistic model for micro-
 in space", Astr. J. 125, 1580-1597 (2003).

al coordinates to unit vector
eflection due to Sun
 aberration
 of r-matrix and p-vector
r to spherical
ze angle into range +/- pi
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

References:

     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
     the Astronomical Almanac, 3rd ed., University Science Books
     (2013).

     Klioner, Sergei A., "A practical relativistic model for micro-
     arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).

  Called:
     eraS2c       spherical coordinates to unit vector
     eraLdsun     light deflection due to Sun
     eraAb        stellar aberration
     eraRxp       product of r-matrix and p-vector
     eraC2s       p-vector to spherical
     eraAnp       normalize angle into range +/- pi

------
atco13
------

Given:
 right ascension at J2000.0 (radians, Note 1)
roper motion (radians/year; Note 2)
proper motion (radians/year)
llax (arcsec)
al velocity (km/s, +ve if receding)
as a 2-part...
uasi Julian Date (Notes 3-4)
UTC (seconds, Note 5)
itude (radians, east +ve, Note 6)
tude (geodetic, radians, Note 6)
ht above ellipsoid (m, geodetic, Notes 6,8)
r motion coordinates (radians, Note 7)
sure at the observer (hPa = mB, Note 8)
ent temperature at the observer (deg C)
tive humidity at the observer (range 0-1)
length (micrometers, Note 9)

rved azimuth (radians: N=0,E=90)
rved zenith distance (radians)
rved hour angle (radians)
rved declination (radians)
rved right ascension (CIO-based, radians)
tion of the origins (ERA-GST)
e):
us: +1 = dubious year (Note 4)
     0 = OK
    -1 = unacceptable date

och other than J2000.0 (for example from the
which has an epoch of J1991.25) will require
to eraPmsafe before use.
n RA is dRA/dt rather than cos(Dec)*dRA/dt.
Julian Date (see Note 2), apportioned in any
een the two arguments, for example where utc1
umber and utc2 is the fraction of a day.
unambiguously represent UTC during a leap
al measures are taken.  The convention in the
 that the JD day represents UTC days whether
, 86400 or 86401 SI seconds.
 use the function eraDtf2d to convert from
ime of day into 2-part quasi Julian Date, as
eap-second-ambiguity convention just

"dubious year" flags UTCs that predate the
 time scale or that are too far in the
d.  See eraDat for further details.
d in IERS bulletins.  It increases by exactly
nd of each positive UTC leap second,
 to keep UT1-UTC within +/- 0.9s.  n.b. This
eview, and in the future UT1-UTC may grow
 limit.
ordinates are with respect to the ERFA_WGS84
.  TAKE CARE WITH THE LONGITUDE SIGN:  the
by the present function is east-positive
, in accordance with geographical convention.
,yp can be obtained from IERS bulletins.  The
dinates (in radians) of the Celestial
ith respect to the International Terrestrial
ee IERS Conventions 2003), measured along the
deg west respectively.  For many
d yp can be set to zero.
bove the ellipsoid of the observing station
nown but phpa, the pressure in hPa (=mB),
equate estimate of hm can be obtained from

tsl * log ( phpa / 1013.25 );
proximate sea-level air temperature in K
Quantities, C.W.Allen, 3rd edition, section
 the pressure phpa is not known, it can be
height of the observing station, hm, as

5 * exp ( -hm / ( 29.3 * tsl ) );
 the refraction is nearly proportional to
at an accurate phpa value is important for

cifies the observing wavelength in
ransition from optical to radio is assumed to
eters (about 3000 GHz).
 result is limited by the corrections for
se a simple A*tan(z) + B*tan^3(z) model.
rological parameters are known accurately and
local effects, the predicted observed
be within 0.05 arcsec (optical) or 1 arcsec
h distance of less than 70 degrees, better
ical or radio) at 85 degrees and better
ical) or 30 arcmin (radio) at the horizon.
 the complementary functions eraAtco13 and
consistent to better than 1 microarcsecond
ial sphere.  With refraction included,
ff at high zenith distances, but is still
csec at 85 degrees.
ans the position that would be seen by a
y aligned theodolite.  (Zenith distance is
titude in order to reflect the fact that no
or depression of the horizon.)  This is
rved HA,Dec via the standard rotation, using
de (corrected for polar motion), while the
are related simply through the Earth rotation
longitude.  "Observed" RA,Dec or HA,Dec thus
that would be seen by a perfect equatorial
 aligned to the Earth's axis of rotation.
take great care with units, as even unlikely
 parameters are accepted and processed in
 models used.

try parameters, ICRS-observed, 2013
CRS to CIRS
IRS to observed
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1)  Star data for an epoch other than J2000.0 (for example from the
      Hipparcos catalog, which has an epoch of J1991.25) will require
      a preliminary call to eraPmsafe before use.

  2)  The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

  3)  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
      convenient way between the two arguments, for example where utc1
      is the Julian Day Number and utc2 is the fraction of a day.

      However, JD cannot unambiguously represent UTC during a leap
      second unless special measures are taken.  The convention in the
      present function is that the JD day represents UTC days whether
      the length is 86399, 86400 or 86401 SI seconds.

      Applications should use the function eraDtf2d to convert from
      calendar date and time of day into 2-part quasi Julian Date, as
      it implements the leap-second-ambiguity convention just
      described.

  4)  The warning status "dubious year" flags UTCs that predate the
      introduction of the time scale or that are too far in the
      future to be trusted.  See eraDat for further details.

  5)  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
      one second at the end of each positive UTC leap second,
      introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
      practice is under review, and in the future UT1-UTC may grow
      essentially without limit.

  6)  The geographical coordinates are with respect to the ERFA_WGS84
      reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
      longitude required by the present function is east-positive
      (i.e. right-handed), in accordance with geographical convention.

  7)  The polar motion xp,yp can be obtained from IERS bulletins.  The
      values are the coordinates (in radians) of the Celestial
      Intermediate Pole with respect to the International Terrestrial
      Reference System (see IERS Conventions 2003), measured along the
      meridians 0 and 90 deg west respectively.  For many
      applications, xp and yp can be set to zero.

  8)  If hm, the height above the ellipsoid of the observing station
      in meters, is not known but phpa, the pressure in hPa (=mB),
      is available, an adequate estimate of hm can be obtained from
      the expression

            hm = -29.3 * tsl * log ( phpa / 1013.25 );

      where tsl is the approximate sea-level air temperature in K
      (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
      52).  Similarly, if the pressure phpa is not known, it can be
      estimated from the height of the observing station, hm, as
      follows:

            phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

      Note, however, that the refraction is nearly proportional to
      the pressure and that an accurate phpa value is important for
      precise work.

  9)  The argument wl specifies the observing wavelength in
      micrometers.  The transition from optical to radio is assumed to
      occur at 100 micrometers (about 3000 GHz).

  10) The accuracy of the result is limited by the corrections for
      refraction, which use a simple A*tan(z) + B*tan^3(z) model.
      Providing the meteorological parameters are known accurately and
      there are no gross local effects, the predicted observed
      coordinates should be within 0.05 arcsec (optical) or 1 arcsec
      (radio) for a zenith distance of less than 70 degrees, better
      than 30 arcsec (optical or radio) at 85 degrees and better
      than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.

      Without refraction, the complementary functions eraAtco13 and
      eraAtoc13 are self-consistent to better than 1 microarcsecond
      all over the celestial sphere.  With refraction included,
      consistency falls off at high zenith distances, but is still
      better than 0.05 arcsec at 85 degrees.

  11) "Observed" Az,ZD means the position that would be seen by a
      perfect geodetically aligned theodolite.  (Zenith distance is
      used rather than altitude in order to reflect the fact that no
      allowance is made for depression of the horizon.)  This is
      related to the observed HA,Dec via the standard rotation, using
      the geodetic latitude (corrected for polar motion), while the
      observed HA and RA are related simply through the Earth rotation
      angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
      means the position that would be seen by a perfect equatorial
      with its polar axis aligned to the Earth's axis of rotation.

  12) It is advisable to take great care with units, as even unlikely
      values of the input parameters are accepted and processed in
      accordance with the models used.


------
atic13
------

Given:
geocentric RA,Dec (radians)
s a 2-part...
lian Date (Note 1)

astrometric RA,Dec (radians)
ion of the origins (ERA-GST, Note 4)

ate2 is a Julian Date, apportioned in any
en the two arguments.  For example,
uld be expressed in any of these ways, among

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 method is best matched to the way the
internally and will deliver the optimum
 method and the date & time methods are both
ween resolution and convenience.  For most
 function the choice will not be at all

ad of TDB without any significant impact on

 are used for the aberration and light
ns so that the functions eraAtic13 (or
i13 (or eraAtciq) are accurate inverses;
the Sun's disk the discrepancy is only about

cy is better than 1 milliarcsecond, limited
sion-nutation model that is used, namely
y close to solar system bodies, additional
ral milliarcseconds can occur because of
ection;  however, the Sun's contribution is
to first order.  The accuracy limitations of
aEpv00 (used to compute Earth position and
bute aberration errors of up to
Light deflection at the Sun's limb is
 mas level.
ation to (equinox based) J2000.0 mean place
han (CIO based) ICRS coordinates, subtract the
ins from the returned right ascension:
eraAnp function can then be applied, as
e result in the conventional 0-2pi range.)

try parameters, ICRS-CIRS, 2013
IRS to ICRS astrometric
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TDB date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TDB)=2450123.7 could be expressed in any of these ways, among
     others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in cases
     where the loss of several decimal digits of resolution is
     acceptable.  The J2000 method is best matched to the way the
     argument is handled internally and will deliver the optimum
     resolution.  The MJD method and the date & time methods are both
     good compromises between resolution and convenience.  For most
     applications of this function the choice will not be at all
     critical.

     TT can be used instead of TDB without any significant impact on
     accuracy.

  2) Iterative techniques are used for the aberration and light
     deflection corrections so that the functions eraAtic13 (or
     eraAticq) and eraAtci13 (or eraAtciq) are accurate inverses;
     even at the edge of the Sun's disk the discrepancy is only about
     1 nanoarcsecond.

  3) The available accuracy is better than 1 milliarcsecond, limited
     mainly by the precession-nutation model that is used, namely
     IAU 2000A/2006.  Very close to solar system bodies, additional
     errors of up to several milliarcseconds can occur because of
     unmodeled light deflection;  however, the Sun's contribution is
     taken into account, to first order.  The accuracy limitations of
     the ERFA function eraEpv00 (used to compute Earth position and
     velocity) can contribute aberration errors of up to
     5 microarcseconds.  Light deflection at the Sun's limb is
     uncertain at the 0.4 mas level.

  4) Should the transformation to (equinox based) J2000.0 mean place
     be required rather than (CIO based) ICRS coordinates, subtract the
     equation of the origins from the returned right ascension:
     RA = RI - EO.  (The eraAnp function can then be applied, as
     required, to keep the result in the conventional 0-2pi range.)


------
aticq
------

Given:
RS RA,Dec (radians)
ar-independent astrometry parameters:
 PM time interval (SSB, Julian years)
 SSB to observer (vector, au)
 Sun to observer (unit vector)
 distance from Sun to observer (au)
 barycentric observer velocity (vector, c)
 sqrt(1-|v|^2): reciprocal of Lorenz factor
 bias-precession-nutation matrix
 longitude + s' (radians)
 polar motion xp wrt local meridian (radians)
 polar motion yp wrt local meridian (radians)
 sine of geodetic latitude
 cosine of geodetic latitude
 magnitude of diurnal aberration vector
 "local" Earth rotation angle (radians)
 refraction constant A (radians)
 refraction constant B (radians)

RS astrometric RA,Dec (radians)

n into account in the light deflection

 are used for the aberration and light
ns so that the functions eraAtic13 (or
i13 (or eraAtciq) are accurate inverses;
the Sun's disk the discrepancy is only about


al coordinates to unit vector
 of transpose of r-matrix and p-vector
vector
 aberration
eflection by the Sun
r to spherical
ze angle into range +/- pi
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Only the Sun is taken into account in the light deflection
     correction.

  2) Iterative techniques are used for the aberration and light
     deflection corrections so that the functions eraAtic13 (or
     eraAticq) and eraAtci13 (or eraAtciq) are accurate inverses;
     even at the edge of the Sun's disk the discrepancy is only about
     1 nanoarcsecond.


------
aticqn
------

Given:
IRS RA,Dec (radians)
tar-independent astrometry parameters:
 PM time interval (SSB, Julian years)
 SSB to observer (vector, au)
 Sun to observer (unit vector)
 distance from Sun to observer (au)
 barycentric observer velocity (vector, c)
 sqrt(1-|v|^2): reciprocal of Lorenz factor
 bias-precession-nutation matrix
 longitude + s' (radians)
 polar motion xp wrt local meridian (radians)
 polar motion yp wrt local meridian (radians)
 sine of geodetic latitude
 cosine of geodetic latitude
 magnitude of diurnal aberration vector
 "local" Earth rotation angle (radians)
 refraction constant A (radians)
 refraction constant B (radians)
 number of bodies (Note 3)
data for each of the n bodies (Notes 3,4):
 mass of the body (solar masses, Note 5)
 deflection limiter (Note 6)
 barycentric PV of the body (au, au/day)

RS astrometric RA,Dec (radians)

 are used for the aberration and light
ns so that the functions eraAticqn and
te inverses; even at the edge of the Sun's
 is only about 1 nanoarcsecond.
flecting body to be taken into account is the
nction can be used instead.
s n entries, one for each body to be
0, no gravitational light deflection will be
r the Sun.
include an entry for the Sun as well as for
body to be taken into account.  The entries
er in which the light passes the body.
b struct for body i, the mass parameter
ired, be adjusted in order to allow for such
e field.
er parameter b[i].dl is phi^2/2, where phi is
on (in radians) between star and body at
plied.  As phi shrinks below the chosen
ction is artificially reduced, reaching zero
le values suitable for a terrestrial
ith masses, are as follows:
m        b[i].dl
         6e-6
5435     3e-9
8574     3e-10
dation of the contents of the b array is
ed masses must be greater than zero, the
y vectors must be right, and the deflection
 zero.

al coordinates to unit vector
 of transpose of r-matrix and p-vector
vector
 aberration
eflection by n bodies
r to spherical
ze angle into range +/- pi
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Iterative techniques are used for the aberration and light
     deflection corrections so that the functions eraAticqn and
     eraAtciqn are accurate inverses; even at the edge of the Sun's
     disk the discrepancy is only about 1 nanoarcsecond.

  2) If the only light-deflecting body to be taken into account is the
     Sun, the eraAticq function can be used instead.

  3) The struct b contains n entries, one for each body to be
     considered.  If n = 0, no gravitational light deflection will be
     applied, not even for the Sun.

  4) The struct b should include an entry for the Sun as well as for
     any planet or other body to be taken into account.  The entries
     should be in the order in which the light passes the body.

  5) In the entry in the b struct for body i, the mass parameter
     b[i].bm can, as required, be adjusted in order to allow for such
     effects as quadrupole field.

  6) The deflection limiter parameter b[i].dl is phi^2/2, where phi is
     the angular separation (in radians) between star and body at
     which limiting is applied.  As phi shrinks below the chosen
     threshold, the deflection is artificially reduced, reaching zero
     for phi = 0.   Example values suitable for a terrestrial
     observer, together with masses, are as follows:

        body i     b[i].bm        b[i].dl

        Sun        1.0            6e-6
        Jupiter    0.00095435     3e-9
        Saturn     0.00028574     3e-10

  7) For efficiency, validation of the contents of the b array is
     omitted.  The supplied masses must be greater than zero, the
     position and velocity vectors must be right, and the deflection
     limiter greater than zero.


------
atio13
------

Given:
 right ascension (CIO-based, radians)
 declination (radians)
as a 2-part...
uasi Julian Date (Notes 1,2)
UTC (seconds, Note 3)
itude (radians, east +ve, Note 4)
etic latitude (radians, Note 4)
ht above ellipsoid (m, geodetic Notes 4,6)
r motion coordinates (radians, Note 5)
sure at the observer (hPa = mB, Note 6)
ent temperature at the observer (deg C)
tive humidity at the observer (range 0-1)
length (micrometers, Note 7)

rved azimuth (radians: N=0,E=90)
rved zenith distance (radians)
rved hour angle (radians)
rved declination (radians)
rved right ascension (CIO-based, radians)
e):
us: +1 = dubious year (Note 2)
     0 = OK
    -1 = unacceptable date

Julian Date (see Note 2), apportioned in any
een the two arguments, for example where utc1
umber and utc2 is the fraction of a day.
unambiguously represent UTC during a leap
al measures are taken.  The convention in the
 that the JD day represents UTC days whether
, 86400 or 86401 SI seconds.
 use the function eraDtf2d to convert from
ime of day into 2-part quasi Julian Date, as
eap-second-ambiguity convention just

"dubious year" flags UTCs that predate the
 time scale or that are too far in the
d.  See eraDat for further details.
d in IERS bulletins.  It increases by exactly
nd of each positive UTC leap second,
 to keep UT1-UTC within +/- 0.9s.  n.b. This
eview, and in the future UT1-UTC may grow
 limit.
ordinates are with respect to the ERFA_WGS84
.  TAKE CARE WITH THE LONGITUDE SIGN:  the
by the present function is east-positive
, in accordance with geographical convention.
,yp can be obtained from IERS bulletins.  The
dinates (in radians) of the Celestial
ith respect to the International Terrestrial
ee IERS Conventions 2003), measured along the
deg west respectively.  For many
d yp can be set to zero.
bove the ellipsoid of the observing station
nown but phpa, the pressure in hPa (=mB), is
ate estimate of hm can be obtained from the

tsl * log ( phpa / 1013.25 );
proximate sea-level air temperature in K
Quantities, C.W.Allen, 3rd edition, section
 the pressure phpa is not known, it can be
height of the observing station, hm, as

5 * exp ( -hm / ( 29.3 * tsl ) );
 the refraction is nearly proportional to
at an accurate phpa value is important for

cifies the observing wavelength in
ransition from optical to radio is assumed to
eters (about 3000 GHz).
ans the position that would be seen by a
y aligned theodolite.  (Zenith distance is
titude in order to reflect the fact that no
or depression of the horizon.)  This is
rved HA,Dec via the standard rotation, using
de (corrected for polar motion), while the
are related simply through the Earth rotation
longitude.  "Observed" RA,Dec or HA,Dec thus
that would be seen by a perfect equatorial
 aligned to the Earth's axis of rotation.
 result is limited by the corrections for
se a simple A*tan(z) + B*tan^3(z) model.
rological parameters are known accurately and
local effects, the predicted astrometric
be within 0.05 arcsec (optical) or 1 arcsec
h distance of less than 70 degrees, better
ical or radio) at 85 degrees and better
ical) or 30 arcmin (radio) at the horizon.
unctions eraAtio13 and eraAtoi13 are self-
r than 1 microarcsecond all over the

take great care with units, as even unlikely
 parameters are accepted and processed in
 models used.

try parameters, CIRS-observed, 2013
IRS to observed
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1)  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
      convenient way between the two arguments, for example where utc1
      is the Julian Day Number and utc2 is the fraction of a day.

      However, JD cannot unambiguously represent UTC during a leap
      second unless special measures are taken.  The convention in the
      present function is that the JD day represents UTC days whether
      the length is 86399, 86400 or 86401 SI seconds.

      Applications should use the function eraDtf2d to convert from
      calendar date and time of day into 2-part quasi Julian Date, as
      it implements the leap-second-ambiguity convention just
      described.

  2)  The warning status "dubious year" flags UTCs that predate the
      introduction of the time scale or that are too far in the
      future to be trusted.  See eraDat for further details.

  3)  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
      one second at the end of each positive UTC leap second,
      introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
      practice is under review, and in the future UT1-UTC may grow
      essentially without limit.

  4)  The geographical coordinates are with respect to the ERFA_WGS84
      reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
      longitude required by the present function is east-positive
      (i.e. right-handed), in accordance with geographical convention.

  5)  The polar motion xp,yp can be obtained from IERS bulletins.  The
      values are the coordinates (in radians) of the Celestial
      Intermediate Pole with respect to the International Terrestrial
      Reference System (see IERS Conventions 2003), measured along the
      meridians 0 and 90 deg west respectively.  For many
      applications, xp and yp can be set to zero.

  6)  If hm, the height above the ellipsoid of the observing station
      in meters, is not known but phpa, the pressure in hPa (=mB), is
      available, an adequate estimate of hm can be obtained from the
      expression

            hm = -29.3 * tsl * log ( phpa / 1013.25 );

      where tsl is the approximate sea-level air temperature in K
      (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
      52).  Similarly, if the pressure phpa is not known, it can be
      estimated from the height of the observing station, hm, as
      follows:

            phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

      Note, however, that the refraction is nearly proportional to
      the pressure and that an accurate phpa value is important for
      precise work.

  7)  The argument wl specifies the observing wavelength in
      micrometers.  The transition from optical to radio is assumed to
      occur at 100 micrometers (about 3000 GHz).

  8)  "Observed" Az,ZD means the position that would be seen by a
      perfect geodetically aligned theodolite.  (Zenith distance is
      used rather than altitude in order to reflect the fact that no
      allowance is made for depression of the horizon.)  This is
      related to the observed HA,Dec via the standard rotation, using
      the geodetic latitude (corrected for polar motion), while the
      observed HA and RA are related simply through the Earth rotation
      angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
      means the position that would be seen by a perfect equatorial
      with its polar axis aligned to the Earth's axis of rotation.

  9)  The accuracy of the result is limited by the corrections for
      refraction, which use a simple A*tan(z) + B*tan^3(z) model.
      Providing the meteorological parameters are known accurately and
      there are no gross local effects, the predicted astrometric
      coordinates should be within 0.05 arcsec (optical) or 1 arcsec
      (radio) for a zenith distance of less than 70 degrees, better
      than 30 arcsec (optical or radio) at 85 degrees and better
      than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.

  10) The complementary functions eraAtio13 and eraAtoi13 are self-
      consistent to better than 1 microarcsecond all over the
      celestial sphere.

  11) It is advisable to take great care with units, as even unlikely
      values of the input parameters are accepted and processed in
      accordance with the models used.


------
atioq
------

Given:
RS right ascension
RS declination
ar-independent astrometry parameters:
 PM time interval (SSB, Julian years)
 SSB to observer (vector, au)
 Sun to observer (unit vector)
 distance from Sun to observer (au)
 barycentric observer velocity (vector, c)
 sqrt(1-|v|^2): reciprocal of Lorenz factor
 bias-precession-nutation matrix
 longitude + s' (radians)
 polar motion xp wrt local meridian (radians)
 polar motion yp wrt local meridian (radians)
 sine of geodetic latitude
 cosine of geodetic latitude
 magnitude of diurnal aberration vector
 "local" Earth rotation angle (radians)
 refraction constant A (radians)
 refraction constant B (radians)

served azimuth (radians: N=0,E=90)
served zenith distance (radians)
served hour angle (radians)
served declination (radians)
served right ascension (CIO-based, radians)

s zenith distance rather than altitude in
 fact that no allowance is made for
rizon.
result is limited by the corrections for
e a simple A*tan(z) + B*tan^3(z) model.
ological parameters are known accurately and
ocal effects, the predicted observed
e within 0.05 arcsec (optical) or 1 arcsec
 distance of less than 70 degrees, better
cal or radio) at 85 degrees and better
cal) or 30 arcmin (radio) at the horizon.
the complementary functions eraAtioq and
nsistent to better than 1 microarcsecond all
phere.  With refraction included, consistency
nith distances, but is still better than
grees.
ake great care with units, as even unlikely
parameters are accepted and processed in
models used.
btained from a star catalog mean place by
otion, parallax, the Sun's gravitational lens
ation and precession-nutation.  For star
S, these effects can be applied by means of
 functions.  Starting from classical "mean
tional transformations will be needed first.
ns the position that would be seen by a
 aligned theodolite.  This is obtained from
llowing for Earth orientation and diurnal
 from equator to horizon coordinates, and
efraction.  The HA,Dec is obtained by
quatorial coordinates, and is the position
y a perfect equatorial with its polar axis
's axis of rotation.  Finally, the RA is
ing the HA from the local ERA.
 CIRS-to-observed-place parameters in ASTROM
 eraApio[13] or eraApco[13].  If nothing has
y except the time, eraAper[13] may be used to
e adjustment to the astrom structure.

al coordinates to unit vector
r to spherical
ze angle into range 0 to 2pi
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) This function returns zenith distance rather than altitude in
     order to reflect the fact that no allowance is made for
     depression of the horizon.

  2) The accuracy of the result is limited by the corrections for
     refraction, which use a simple A*tan(z) + B*tan^3(z) model.
     Providing the meteorological parameters are known accurately and
     there are no gross local effects, the predicted observed
     coordinates should be within 0.05 arcsec (optical) or 1 arcsec
     (radio) for a zenith distance of less than 70 degrees, better
     than 30 arcsec (optical or radio) at 85 degrees and better
     than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.

     Without refraction, the complementary functions eraAtioq and
     eraAtoiq are self-consistent to better than 1 microarcsecond all
     over the celestial sphere.  With refraction included, consistency
     falls off at high zenith distances, but is still better than
     0.05 arcsec at 85 degrees.

  3) It is advisable to take great care with units, as even unlikely
     values of the input parameters are accepted and processed in
     accordance with the models used.

  4) The CIRS RA,Dec is obtained from a star catalog mean place by
     allowing for space motion, parallax, the Sun's gravitational lens
     effect, annual aberration and precession-nutation.  For star
     positions in the ICRS, these effects can be applied by means of
     the eraAtci13 (etc.) functions.  Starting from classical "mean
     place" systems, additional transformations will be needed first.

  5) "Observed" Az,El means the position that would be seen by a
     perfect geodetically aligned theodolite.  This is obtained from
     the CIRS RA,Dec by allowing for Earth orientation and diurnal
     aberration, rotating from equator to horizon coordinates, and
     then adjusting for refraction.  The HA,Dec is obtained by
     rotating back into equatorial coordinates, and is the position
     that would be seen by a perfect equatorial with its polar axis
     aligned to the Earth's axis of rotation.  Finally, the RA is
     obtained by subtracting the HA from the local ERA.

  6) The star-independent CIRS-to-observed-place parameters in ASTROM
     may be computed with eraApio[13] or eraApco[13].  If nothing has
     changed significantly except the time, eraAper[13] may be used to
     perform the requisite adjustment to the astrom structure.


------
atoc13
------

Given:
 of coordinates - "R", "H" or "A" (Notes 1,2)
rved Az, HA or RA (radians; Az is N=0,E=90)
rved ZD or Dec (radians)
as a 2-part...
uasi Julian Date (Notes 3,4)
UTC (seconds, Note 5)
itude (radians, east +ve, Note 6)
etic latitude (radians, Note 6)
ht above ellipsoid (m, geodetic Notes 6,8)
r motion coordinates (radians, Note 7)
sure at the observer (hPa = mB, Note 8)
ent temperature at the observer (deg C)
tive humidity at the observer (range 0-1)
length (micrometers, Note 9)

 astrometric RA,Dec (radians)
e):
us: +1 = dubious year (Note 4)
     0 = OK
    -1 = unacceptable date

ans the position that would be seen by a
y aligned theodolite.  (Zenith distance is
titude in order to reflect the fact that no
or depression of the horizon.)  This is
rved HA,Dec via the standard rotation, using
de (corrected for polar motion), while the
are related simply through the Earth rotation
longitude.  "Observed" RA,Dec or HA,Dec thus
that would be seen by a perfect equatorial
 aligned to the Earth's axis of rotation.
acter of the type argument is significant.
s that ob1 and ob2 are the observed right
nation;  "H" or "h" indicates that they are
e) and declination;  anything else ("A" or
 indicates that ob1 and ob2 are azimuth
0 deg) and zenith distance.
Julian Date (see Note 2), apportioned in any
een the two arguments, for example where utc1
umber and utc2 is the fraction of a day.
unambiguously represent UTC during a leap
al measures are taken.  The convention in the
 that the JD day represents UTC days whether
, 86400 or 86401 SI seconds.
 use the function eraDtf2d to convert from
ime of day into 2-part quasi Julian Date, as
eap-second-ambiguity convention just

"dubious year" flags UTCs that predate the
 time scale or that are too far in the
d.  See eraDat for further details.
d in IERS bulletins.  It increases by exactly
nd of each positive UTC leap second,
 to keep UT1-UTC within +/- 0.9s.  n.b. This
eview, and in the future UT1-UTC may grow
 limit.
ordinates are with respect to the ERFA_WGS84
.  TAKE CARE WITH THE LONGITUDE SIGN:  the
by the present function is east-positive
, in accordance with geographical convention.
,yp can be obtained from IERS bulletins.  The
dinates (in radians) of the Celestial
ith respect to the International Terrestrial
ee IERS Conventions 2003), measured along the
deg west respectively.  For many
d yp can be set to zero.
bove the ellipsoid of the observing station
nown but phpa, the pressure in hPa (=mB), is
ate estimate of hm can be obtained from the

tsl * log ( phpa / 1013.25 );
proximate sea-level air temperature in K
Quantities, C.W.Allen, 3rd edition, section
 the pressure phpa is not known, it can be
height of the observing station, hm, as

5 * exp ( -hm / ( 29.3 * tsl ) );
 the refraction is nearly proportional to
at an accurate phpa value is important for

cifies the observing wavelength in
ransition from optical to radio is assumed to
eters (about 3000 GHz).
 result is limited by the corrections for
se a simple A*tan(z) + B*tan^3(z) model.
rological parameters are known accurately and
local effects, the predicted astrometric
be within 0.05 arcsec (optical) or 1 arcsec
h distance of less than 70 degrees, better
ical or radio) at 85 degrees and better
ical) or 30 arcmin (radio) at the horizon.
 the complementary functions eraAtco13 and
consistent to better than 1 microarcsecond
ial sphere.  With refraction included,
ff at high zenith distances, but is still
csec at 85 degrees.
take great care with units, as even unlikely
 parameters are accepted and processed in
 models used.

try parameters, ICRS-observed
bserved to CIRS
IRS to ICRS
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1)  "Observed" Az,ZD means the position that would be seen by a
      perfect geodetically aligned theodolite.  (Zenith distance is
      used rather than altitude in order to reflect the fact that no
      allowance is made for depression of the horizon.)  This is
      related to the observed HA,Dec via the standard rotation, using
      the geodetic latitude (corrected for polar motion), while the
      observed HA and RA are related simply through the Earth rotation
      angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
      means the position that would be seen by a perfect equatorial
      with its polar axis aligned to the Earth's axis of rotation.

  2)  Only the first character of the type argument is significant.
      "R" or "r" indicates that ob1 and ob2 are the observed right
      ascension and declination;  "H" or "h" indicates that they are
      hour angle (west +ve) and declination;  anything else ("A" or
      "a" is recommended) indicates that ob1 and ob2 are azimuth
      (north zero, east 90 deg) and zenith distance.

  3)  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
      convenient way between the two arguments, for example where utc1
      is the Julian Day Number and utc2 is the fraction of a day.

      However, JD cannot unambiguously represent UTC during a leap
      second unless special measures are taken.  The convention in the
      present function is that the JD day represents UTC days whether
      the length is 86399, 86400 or 86401 SI seconds.

      Applications should use the function eraDtf2d to convert from
      calendar date and time of day into 2-part quasi Julian Date, as
      it implements the leap-second-ambiguity convention just
      described.

  4)  The warning status "dubious year" flags UTCs that predate the
      introduction of the time scale or that are too far in the
      future to be trusted.  See eraDat for further details.

  5)  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
      one second at the end of each positive UTC leap second,
      introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
      practice is under review, and in the future UT1-UTC may grow
      essentially without limit.

  6)  The geographical coordinates are with respect to the ERFA_WGS84
      reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
      longitude required by the present function is east-positive
      (i.e. right-handed), in accordance with geographical convention.

  7)  The polar motion xp,yp can be obtained from IERS bulletins.  The
      values are the coordinates (in radians) of the Celestial
      Intermediate Pole with respect to the International Terrestrial
      Reference System (see IERS Conventions 2003), measured along the
      meridians 0 and 90 deg west respectively.  For many
      applications, xp and yp can be set to zero.

  8)  If hm, the height above the ellipsoid of the observing station
      in meters, is not known but phpa, the pressure in hPa (=mB), is
      available, an adequate estimate of hm can be obtained from the
      expression

            hm = -29.3 * tsl * log ( phpa / 1013.25 );

      where tsl is the approximate sea-level air temperature in K
      (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
      52).  Similarly, if the pressure phpa is not known, it can be
      estimated from the height of the observing station, hm, as
      follows:

            phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

      Note, however, that the refraction is nearly proportional to
      the pressure and that an accurate phpa value is important for
      precise work.

  9)  The argument wl specifies the observing wavelength in
      micrometers.  The transition from optical to radio is assumed to
      occur at 100 micrometers (about 3000 GHz).

  10) The accuracy of the result is limited by the corrections for
      refraction, which use a simple A*tan(z) + B*tan^3(z) model.
      Providing the meteorological parameters are known accurately and
      there are no gross local effects, the predicted astrometric
      coordinates should be within 0.05 arcsec (optical) or 1 arcsec
      (radio) for a zenith distance of less than 70 degrees, better
      than 30 arcsec (optical or radio) at 85 degrees and better
      than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.

      Without refraction, the complementary functions eraAtco13 and
      eraAtoc13 are self-consistent to better than 1 microarcsecond
      all over the celestial sphere.  With refraction included,
      consistency falls off at high zenith distances, but is still
      better than 0.05 arcsec at 85 degrees.

  11) It is advisable to take great care with units, as even unlikely
      values of the input parameters are accepted and processed in
      accordance with the models used.


------
atoi13
------

Given:
 of coordinates - "R", "H" or "A" (Notes 1,2)
rved Az, HA or RA (radians; Az is N=0,E=90)
rved ZD or Dec (radians)
as a 2-part...
uasi Julian Date (Notes 3,4)
UTC (seconds, Note 5)
itude (radians, east +ve, Note 6)
etic latitude (radians, Note 6)
ht above the ellipsoid (meters, Notes 6,8)
r motion coordinates (radians, Note 7)
sure at the observer (hPa = mB, Note 8)
ent temperature at the observer (deg C)
tive humidity at the observer (range 0-1)
length (micrometers, Note 9)

 right ascension (CIO-based, radians)
 declination (radians)
e):
us: +1 = dubious year (Note 2)
     0 = OK
    -1 = unacceptable date

ans the position that would be seen by a
y aligned theodolite.  (Zenith distance is
titude in order to reflect the fact that no
or depression of the horizon.)  This is
rved HA,Dec via the standard rotation, using
de (corrected for polar motion), while the
are related simply through the Earth rotation
longitude.  "Observed" RA,Dec or HA,Dec thus
that would be seen by a perfect equatorial
 aligned to the Earth's axis of rotation.
acter of the type argument is significant.
s that ob1 and ob2 are the observed right
nation;  "H" or "h" indicates that they are
e) and declination;  anything else ("A" or
 indicates that ob1 and ob2 are azimuth
0 deg) and zenith distance.
Julian Date (see Note 2), apportioned in any
een the two arguments, for example where utc1
umber and utc2 is the fraction of a day.
unambiguously represent UTC during a leap
al measures are taken.  The convention in the
 that the JD day represents UTC days whether
, 86400 or 86401 SI seconds.
 use the function eraDtf2d to convert from
ime of day into 2-part quasi Julian Date, as
eap-second-ambiguity convention just

"dubious year" flags UTCs that predate the
 time scale or that are too far in the
d.  See eraDat for further details.
d in IERS bulletins.  It increases by exactly
nd of each positive UTC leap second,
 to keep UT1-UTC within +/- 0.9s.  n.b. This
eview, and in the future UT1-UTC may grow
 limit.
ordinates are with respect to the ERFA_WGS84
.  TAKE CARE WITH THE LONGITUDE SIGN:  the
by the present function is east-positive
, in accordance with geographical convention.
,yp can be obtained from IERS bulletins.  The
dinates (in radians) of the Celestial
ith respect to the International Terrestrial
ee IERS Conventions 2003), measured along the
deg west respectively.  For many
d yp can be set to zero.
bove the ellipsoid of the observing station
nown but phpa, the pressure in hPa (=mB), is
ate estimate of hm can be obtained from the

tsl * log ( phpa / 1013.25 );
proximate sea-level air temperature in K
Quantities, C.W.Allen, 3rd edition, section
 the pressure phpa is not known, it can be
height of the observing station, hm, as

5 * exp ( -hm / ( 29.3 * tsl ) );
 the refraction is nearly proportional to
at an accurate phpa value is important for

cifies the observing wavelength in
ransition from optical to radio is assumed to
eters (about 3000 GHz).
 result is limited by the corrections for
se a simple A*tan(z) + B*tan^3(z) model.
rological parameters are known accurately and
local effects, the predicted astrometric
be within 0.05 arcsec (optical) or 1 arcsec
h distance of less than 70 degrees, better
ical or radio) at 85 degrees and better
ical) or 30 arcmin (radio) at the horizon.
 the complementary functions eraAtio13 and
consistent to better than 1 microarcsecond
ial sphere.  With refraction included,
ff at high zenith distances, but is still
csec at 85 degrees.
take great care with units, as even unlikely
 parameters are accepted and processed in
 models used.

try parameters, CIRS-observed, 2013
bserved to CIRS
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1)  "Observed" Az,ZD means the position that would be seen by a
      perfect geodetically aligned theodolite.  (Zenith distance is
      used rather than altitude in order to reflect the fact that no
      allowance is made for depression of the horizon.)  This is
      related to the observed HA,Dec via the standard rotation, using
      the geodetic latitude (corrected for polar motion), while the
      observed HA and RA are related simply through the Earth rotation
      angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
      means the position that would be seen by a perfect equatorial
      with its polar axis aligned to the Earth's axis of rotation.

  2)  Only the first character of the type argument is significant.
      "R" or "r" indicates that ob1 and ob2 are the observed right
      ascension and declination;  "H" or "h" indicates that they are
      hour angle (west +ve) and declination;  anything else ("A" or
      "a" is recommended) indicates that ob1 and ob2 are azimuth
      (north zero, east 90 deg) and zenith distance.

  3)  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
      convenient way between the two arguments, for example where utc1
      is the Julian Day Number and utc2 is the fraction of a day.

      However, JD cannot unambiguously represent UTC during a leap
      second unless special measures are taken.  The convention in the
      present function is that the JD day represents UTC days whether
      the length is 86399, 86400 or 86401 SI seconds.

      Applications should use the function eraDtf2d to convert from
      calendar date and time of day into 2-part quasi Julian Date, as
      it implements the leap-second-ambiguity convention just
      described.

  4)  The warning status "dubious year" flags UTCs that predate the
      introduction of the time scale or that are too far in the
      future to be trusted.  See eraDat for further details.

  5)  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
      one second at the end of each positive UTC leap second,
      introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
      practice is under review, and in the future UT1-UTC may grow
      essentially without limit.

  6)  The geographical coordinates are with respect to the ERFA_WGS84
      reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
      longitude required by the present function is east-positive
      (i.e. right-handed), in accordance with geographical convention.

  7)  The polar motion xp,yp can be obtained from IERS bulletins.  The
      values are the coordinates (in radians) of the Celestial
      Intermediate Pole with respect to the International Terrestrial
      Reference System (see IERS Conventions 2003), measured along the
      meridians 0 and 90 deg west respectively.  For many
      applications, xp and yp can be set to zero.

  8)  If hm, the height above the ellipsoid of the observing station
      in meters, is not known but phpa, the pressure in hPa (=mB), is
      available, an adequate estimate of hm can be obtained from the
      expression

            hm = -29.3 * tsl * log ( phpa / 1013.25 );

      where tsl is the approximate sea-level air temperature in K
      (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
      52).  Similarly, if the pressure phpa is not known, it can be
      estimated from the height of the observing station, hm, as
      follows:

            phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

      Note, however, that the refraction is nearly proportional to
      the pressure and that an accurate phpa value is important for
      precise work.

  9)  The argument wl specifies the observing wavelength in
      micrometers.  The transition from optical to radio is assumed to
      occur at 100 micrometers (about 3000 GHz).

  10) The accuracy of the result is limited by the corrections for
      refraction, which use a simple A*tan(z) + B*tan^3(z) model.
      Providing the meteorological parameters are known accurately and
      there are no gross local effects, the predicted astrometric
      coordinates should be within 0.05 arcsec (optical) or 1 arcsec
      (radio) for a zenith distance of less than 70 degrees, better
      than 30 arcsec (optical or radio) at 85 degrees and better
      than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.

      Without refraction, the complementary functions eraAtio13 and
      eraAtoi13 are self-consistent to better than 1 microarcsecond
      all over the celestial sphere.  With refraction included,
      consistency falls off at high zenith distances, but is still
      better than 0.05 arcsec at 85 degrees.

  12) It is advisable to take great care with units, as even unlikely
      values of the input parameters are accepted and processed in
      accordance with the models used.


------
atoiq
------

Given:
pe of coordinates: "R", "H" or "A" (Note 1)
served Az, HA or RA (radians; Az is N=0,E=90)
served ZD or Dec (radians)
ar-independent astrometry parameters:
 PM time interval (SSB, Julian years)
 SSB to observer (vector, au)
 Sun to observer (unit vector)
 distance from Sun to observer (au)
 barycentric observer velocity (vector, c)
 sqrt(1-|v|^2): reciprocal of Lorenz factor
 bias-precession-nutation matrix
 longitude + s' (radians)
 polar motion xp wrt local meridian (radians)
 polar motion yp wrt local meridian (radians)
 sine of geodetic latitude
 cosine of geodetic latitude
 magnitude of diurnal aberration vector
 "local" Earth rotation angle (radians)
 refraction constant A (radians)
 refraction constant B (radians)

RS right ascension (CIO-based, radians)
RS declination (radians)

ns the position that would be seen by a
 aligned theodolite.  This is related to
via the standard rotation, using the geodetic
for polar motion), while the observed HA and
y through the Earth rotation angle and the
served" RA,Dec or HA,Dec thus means the
be seen by a perfect equatorial with its
o the Earth's axis of rotation.  By removing
ace the effects of atmospheric refraction and
the CIRS RA,Dec is obtained.
cter of the type argument is significant.
 that ob1 and ob2 are the observed right
ation;  "H" or "h" indicates that they are
) and declination;  anything else ("A" or
indicates that ob1 and ob2 are azimuth (north
nd zenith distance.  (Zenith distance is used
 in order to reflect the fact that no
r depression of the horizon.)
result is limited by the corrections for
e a simple A*tan(z) + B*tan^3(z) model.
ological parameters are known accurately and
ocal effects, the predicted observed
e within 0.05 arcsec (optical) or 1 arcsec
 distance of less than 70 degrees, better
cal or radio) at 85 degrees and better than
or 30 arcmin (radio) at the horizon.
the complementary functions eraAtioq and
nsistent to better than 1 microarcsecond all
phere.  With refraction included, consistency
nith distances, but is still better than
grees.
ake great care with units, as even unlikely
parameters are accepted and processed in
models used.

al coordinates to unit vector
r to spherical
ze angle into range 0 to 2pi
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) "Observed" Az,El means the position that would be seen by a
     perfect geodetically aligned theodolite.  This is related to
     the observed HA,Dec via the standard rotation, using the geodetic
     latitude (corrected for polar motion), while the observed HA and
     RA are related simply through the Earth rotation angle and the
     site longitude.  "Observed" RA,Dec or HA,Dec thus means the
     position that would be seen by a perfect equatorial with its
     polar axis aligned to the Earth's axis of rotation.  By removing
     from the observed place the effects of atmospheric refraction and
     diurnal aberration, the CIRS RA,Dec is obtained.

  2) Only the first character of the type argument is significant.
     "R" or "r" indicates that ob1 and ob2 are the observed right
     ascension and declination;  "H" or "h" indicates that they are
     hour angle (west +ve) and declination;  anything else ("A" or
     "a" is recommended) indicates that ob1 and ob2 are azimuth (north
     zero, east 90 deg) and zenith distance.  (Zenith distance is used
     rather than altitude in order to reflect the fact that no
     allowance is made for depression of the horizon.)

  3) The accuracy of the result is limited by the corrections for
     refraction, which use a simple A*tan(z) + B*tan^3(z) model.
     Providing the meteorological parameters are known accurately and
     there are no gross local effects, the predicted observed
     coordinates should be within 0.05 arcsec (optical) or 1 arcsec
     (radio) for a zenith distance of less than 70 degrees, better
     than 30 arcsec (optical or radio) at 85 degrees and better than
     20 arcmin (optical) or 30 arcmin (radio) at the horizon.

     Without refraction, the complementary functions eraAtioq and
     eraAtoiq are self-consistent to better than 1 microarcsecond all
     over the celestial sphere.  With refraction included, consistency
     falls off at high zenith distances, but is still better than
     0.05 arcsec at 85 degrees.

  4) It is advisable to take great care with units, as even unlikely
     values of the input parameters are accepted and processed in
     accordance with the models used.


------
bi00
------


Notes:

  1) The frame bias corrections in longitude and obliquity (radians)
     are required in order to correct for the offset between the GCRS
     pole and the mean J2000.0 pole.  They define, with respect to the
     GCRS frame, a J2000.0 mean pole that is consistent with the rest
     of the IAU 2000A precession-nutation model.

  2) In addition to the displacement of the pole, the complete
     description of the frame bias requires also an offset in right
     ascension.  This is not part of the IAU 2000A model, and is from
     Chapront et al. (2002).  It is returned in radians.

  3) This is a supplemented implementation of one aspect of the IAU
     2000A nutation model, formally adopted by the IAU General
     Assembly in 2000, namely MHB2000 (Mathews et al. 2002).


References:

     Chapront, J., Chapront-Touze, M. & Francou, G., Astron.
     Astrophys., 387, 700, 2002.

     Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation
     and precession   New nutation series for nonrigid Earth and
     insights into the Earth's interior", J.Geophys.Res., 107, B4,
     2002.  The MHB2000 code itself was obtained on 9th September 2002
     from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.

------
bp00
------

Given:
        TT as a 2-part Julian Date (Note 1)

3][3]   frame bias matrix (Note 2)
3][3]   precession matrix (Note 3)
3][3]   bias-precession matrix (Note 4)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
orms vectors from GCRS to mean J2000.0 by

orms vectors from J2000.0 mean equator and
tor and equinox of date by applying

forms vectors from GCRS to mean equator and
pplying frame bias then precession.  It is

 re-use the same array in the returned
ys are filled in the order given.

ias components, IAU 2000
0 precession adjustments
ize r-matrix to identity
around X-axis
around Y-axis
around Z-axis
matrix
 of two r-matrices

 Celestial Intermediate Pole and Celestial
sistent with the IAU 2000A precession-
ron.Astrophys. 400, 1145-1154 (2003)
phemeris origin (CEO) was renamed "celestial
igin" (CIO) by IAU 2006 Resolution 2.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

             date1         date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix rb transforms vectors from GCRS to mean J2000.0 by
     applying frame bias.

  3) The matrix rp transforms vectors from J2000.0 mean equator and
     equinox to mean equator and equinox of date by applying
     precession.

  4) The matrix rbp transforms vectors from GCRS to mean equator and
     equinox of date by applying frame bias then precession.  It is
     the product rp x rb.

  5) It is permissible to re-use the same array in the returned
     arguments.  The arrays are filled in the order given.


------
bp06
------

Given:
        TT as a 2-part Julian Date (Note 1)

3][3]   frame bias matrix (Note 2)
3][3]   precession matrix (Note 3)
3][3]   bias-precession matrix (Note 4)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
orms vectors from GCRS to mean J2000.0 by

orms vectors from mean J2000.0 to mean of
cession.
forms vectors from GCRS to mean of date by
then precession.  It is the product rp x rb.
 re-use the same array in the returned
ys are filled in the order given.

ecession F-W angles, IAU 2006
les to r-matrix
ix, IAU 2006
se r-matrix
 of two r-matrices
matrix

ace, P.T., 2006, Astron.Astrophys. 450, 855
taine, N., 2006, Astron.Astrophys. 459, 981
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

             date1         date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix rb transforms vectors from GCRS to mean J2000.0 by
     applying frame bias.

  3) The matrix rp transforms vectors from mean J2000.0 to mean of
     date by applying precession.

  4) The matrix rbp transforms vectors from GCRS to mean of date by
     applying frame bias then precession.  It is the product rp x rb.

  5) It is permissible to re-use the same array in the returned
     arguments.  The arrays are filled in the order given.


References:

     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981

------
bpn2xy
------

Given:
3]  celestial-to-true matrix (Note 1)

    Celestial Intermediate Pole (Note 2)

sforms vectors from GCRS to true equator (and
ate, and therefore the Celestial Intermediate
the bottom row of the matrix.
e components of the Celestial Intermediate
the Geocentric Celestial Reference System.

 Celestial Intermediate Pole and Celestial
sistent with the IAU 2000A precession-
ron.Astrophys. 400, 1145-1154

phemeris origin (CEO) was renamed "celestial
igin" (CIO) by IAU 2006 Resolution 2.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
c2i00a
------

Given:
     TT as a 2-part Julian Date (Note 1)

][3] celestial-to-intermediate matrix (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
he first stage in the transformation from
rial coordinates:
R_3(ERA) * rc2i * [CRS]
[CRS]
tor in the Geocentric Celestial Reference
a vector in the International Terrestrial
e IERS Conventions 2003), ERA is the Earth
POM is the polar motion matrix.
ly less accurate result (about 1 mas), can be
stead the eraC2i00b function.

al NPB matrix, IAU 2000A
al-to-intermediate matrix, given NPB matrix

 Celestial Intermediate Pole and Celestial
sistent with the IAU 2000A precession-
ron.Astrophys. 400, 1145-1154

phemeris origin (CEO) was renamed "celestial
igin" (CIO) by IAU 2006 Resolution 2.
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix rc2i is the first stage in the transformation from
     celestial to terrestrial coordinates:

        [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]

               =  rc2t * [CRS]

     where [CRS] is a vector in the Geocentric Celestial Reference
     System and [TRS] is a vector in the International Terrestrial
     Reference System (see IERS Conventions 2003), ERA is the Earth
     Rotation Angle and RPOM is the polar motion matrix.

  3) A faster, but slightly less accurate result (about 1 mas), can be
     obtained by using instead the eraC2i00b function.


References:

     "Expressions for the Celestial Intermediate Pole and Celestial
     Ephemeris Origin consistent with the IAU 2000A precession-
     nutation model", Astron.Astrophys. 400, 1145-1154
     (2003)

     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
          intermediate origin" (CIO) by IAU 2006 Resolution 2.

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

------
c2i00b
------

Given:
     TT as a 2-part Julian Date (Note 1)

][3] celestial-to-intermediate matrix (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
he first stage in the transformation from
rial coordinates:
R_3(ERA) * rc2i * [CRS]
[CRS]
tor in the Geocentric Celestial Reference
a vector in the International Terrestrial
e IERS Conventions 2003), ERA is the Earth
POM is the polar motion matrix.
 is faster, but slightly less accurate (about
C2i00a function.

al NPB matrix, IAU 2000B
al-to-intermediate matrix, given NPB matrix

 Celestial Intermediate Pole and Celestial
sistent with the IAU 2000A precession-
ron.Astrophys. 400, 1145-1154

phemeris origin (CEO) was renamed "celestial
igin" (CIO) by IAU 2006 Resolution 2.
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix rc2i is the first stage in the transformation from
     celestial to terrestrial coordinates:

        [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]

               =  rc2t * [CRS]

     where [CRS] is a vector in the Geocentric Celestial Reference
     System and [TRS] is a vector in the International Terrestrial
     Reference System (see IERS Conventions 2003), ERA is the Earth
     Rotation Angle and RPOM is the polar motion matrix.

  3) The present function is faster, but slightly less accurate (about
     1 mas), than the eraC2i00a function.


References:

     "Expressions for the Celestial Intermediate Pole and Celestial
     Ephemeris Origin consistent with the IAU 2000A precession-
     nutation model", Astron.Astrophys. 400, 1145-1154
     (2003)

     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
          intermediate origin" (CIO) by IAU 2006 Resolution 2.

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

------
c2i06a
------

Given:
     TT as a 2-part Julian Date (Note 1)

][3] celestial-to-intermediate matrix (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
he first stage in the transformation from
rial coordinates:
R_3(ERA) * rc2i * [CRS]
[CRS]
tor in the Geocentric Celestial Reference
a vector in the International Terrestrial
e IERS Conventions 2003), ERA is the Earth
POM is the polar motion matrix.

al NPB matrix, IAU 2006/2000A
 CIP X,Y coordinates from NPB matrix
 locator s, given X,Y, IAU 2006
al-to-intermediate matrix, given X,Y and s

it, G. (eds.), 2004, IERS Conventions (2003),
No. 32, BKG
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix rc2i is the first stage in the transformation from
     celestial to terrestrial coordinates:

        [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]

               =  RC2T * [CRS]

     where [CRS] is a vector in the Geocentric Celestial Reference
     System and [TRS] is a vector in the International Terrestrial
     Reference System (see IERS Conventions 2003), ERA is the Earth
     Rotation Angle and RPOM is the polar motion matrix.


References:

     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
     IERS Technical Note No. 32, BKG

------
c2ibpn
------

Given:
     TT as a 2-part Julian Date (Note 1)
][3] celestial-to-true matrix (Note 2)

][3] celestial-to-intermediate matrix (Note 3)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
sforms vectors from GCRS to true equator (and
ate.  Only the CIP (bottom row) is used.
he first stage in the transformation from
rial coordinates:
3(ERA) * rc2i * [CRS]
RS]
tor in the Geocentric Celestial Reference
a vector in the International Terrestrial
e IERS Conventions 2003), ERA is the Earth
POM is the polar motion matrix.
es not include "00", This function is in fact
2000 models.

 CIP X,Y coordinates from NPB matrix
al-to-intermediate matrix, given X,Y

 Celestial Intermediate Pole and Celestial
sistent with the IAU 2000A precession-
ron.Astrophys. 400, 1145-1154 (2003)
phemeris origin (CEO) was renamed "celestial
igin" (CIO) by IAU 2006 Resolution 2.
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix rbpn transforms vectors from GCRS to true equator (and
     CIO or equinox) of date.  Only the CIP (bottom row) is used.

  3) The matrix rc2i is the first stage in the transformation from
     celestial to terrestrial coordinates:

        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]

              = RC2T * [CRS]

     where [CRS] is a vector in the Geocentric Celestial Reference
     System and [TRS] is a vector in the International Terrestrial
     Reference System (see IERS Conventions 2003), ERA is the Earth
     Rotation Angle and RPOM is the polar motion matrix.

  4) Although its name does not include "00", This function is in fact
     specific to the IAU 2000 models.


References:
     "Expressions for the Celestial Intermediate Pole and Celestial
     Ephemeris Origin consistent with the IAU 2000A precession-
     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
          intermediate origin" (CIO) by IAU 2006 Resolution 2.

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

------
c2ixy
------

Given:
     TT as a 2-part Julian Date (Note 1)
     Celestial Intermediate Pole (Note 2)

][3] celestial-to-intermediate matrix (Note 3)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
ediate Pole coordinates are the x,y components
n the Geocentric Celestial Reference System.
he first stage in the transformation from
rial coordinates:
3(ERA) * rc2i * [CRS]
RS]
tor in the Geocentric Celestial Reference
a vector in the International Terrestrial
e IERS Conventions 2003), ERA is the Earth
POM is the polar motion matrix.
es not include "00", This function is in fact
2000 models.

al-to-intermediate matrix, given X,Y and s
 locator s, given X,Y, IAU 2000A

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The Celestial Intermediate Pole coordinates are the x,y components
     of the unit vector in the Geocentric Celestial Reference System.

  3) The matrix rc2i is the first stage in the transformation from
     celestial to terrestrial coordinates:

        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]

              = RC2T * [CRS]

     where [CRS] is a vector in the Geocentric Celestial Reference
     System and [TRS] is a vector in the International Terrestrial
     Reference System (see IERS Conventions 2003), ERA is the Earth
     Rotation Angle and RPOM is the polar motion matrix.

  4) Although its name does not include "00", This function is in fact
     specific to the IAU 2000 models.


------
c2ixys
------

Given:
    Celestial Intermediate Pole (Note 1)
    the CIO locator s (Note 2)

]   celestial-to-intermediate matrix (Note 3)

ediate Pole coordinates are the x,y
it vector in the Geocentric Celestial

n radians) positions the Celestial
on the equator of the CIP.
he first stage in the transformation from
rial coordinates:
3(ERA) * rc2i * [CRS]
RS]
tor in the Geocentric Celestial Reference
a vector in the International Terrestrial
e IERS Conventions 2003), ERA is the Earth
POM is the polar motion matrix.

ize r-matrix to identity
around Z-axis
around Y-axis

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The Celestial Intermediate Pole coordinates are the x,y
     components of the unit vector in the Geocentric Celestial
     Reference System.

  2) The CIO locator s (in radians) positions the Celestial
     Intermediate Origin on the equator of the CIP.

  3) The matrix rc2i is the first stage in the transformation from
     celestial to terrestrial coordinates:

        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]

              = RC2T * [CRS]

     where [CRS] is a vector in the Geocentric Celestial Reference
     System and [TRS] is a vector in the International Terrestrial
     Reference System (see IERS Conventions 2003), ERA is the Earth
     Rotation Angle and RPOM is the polar motion matrix.


------
c2s
------

Given:
p-vector

longitude angle (radians)
latitude angle (radians)

e any magnitude; only its direction is used.
heta and phi are returned.
 theta is returned.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
c2t00a
------

Given:
    TT as a 2-part Julian Date (Note 1)
    UT1 as a 2-part Julian Date (Note 1)
    coordinates of the pole (radians, Note 2)

]   celestial-to-terrestrial matrix (Note 3)

 tta+ttb and uta+utb are Julian Dates,
onvenient way between the arguments uta and
D(UT1)=2450123.7 could be expressed in any of
hers:
   utb
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution is
00 and MJD methods are good compromises
nd convenience.  In the case of uta,utb, the
s best matched to the Earth rotation angle
imum precision is delivered when the uta
 UT1 on the day in question and the utb
 range 0 to 1, or vice versa.
 yp are the coordinates (in radians) of the
te Pole with respect to the International
e System (see IERS Conventions 2003),
eridians to 0 and 90 deg west respectively.
sforms from celestial to terrestrial

3(ERA) * RC2I * [CRS]
RS]
tor in the Geocentric Celestial Reference
a vector in the International Terrestrial
e IERS Conventions 2003), RC2I is the
diate matrix, ERA is the Earth rotation
e polar motion matrix.
ly less accurate result (about 1 mas), can
 instead the eraC2t00b function.

al-to-intermediate matrix, IAU 2000A
otation angle, IAU 2000
 locator s', IERS 2000
otion matrix
O-based celestial-to-terrestrial matrix

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
     apportioned in any convenient way between the arguments uta and
     utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
     these ways, among others:

             uta            utb

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution is
     acceptable.  The J2000 and MJD methods are good compromises
     between resolution and convenience.  In the case of uta,utb, the
     date & time method is best matched to the Earth rotation angle
     algorithm used:  maximum precision is delivered when the uta
     argument is for 0hrs UT1 on the day in question and the utb
     argument lies in the range 0 to 1, or vice versa.

  2) The arguments xp and yp are the coordinates (in radians) of the
     Celestial Intermediate Pole with respect to the International
     Terrestrial Reference System (see IERS Conventions 2003),
     measured along the meridians to 0 and 90 deg west respectively.

  3) The matrix rc2t transforms from celestial to terrestrial
     coordinates:

        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]

              = rc2t * [CRS]

     where [CRS] is a vector in the Geocentric Celestial Reference
     System and [TRS] is a vector in the International Terrestrial
     Reference System (see IERS Conventions 2003), RC2I is the
     celestial-to-intermediate matrix, ERA is the Earth rotation
     angle and RPOM is the polar motion matrix.

  4) A faster, but slightly less accurate result (about 1 mas), can
     be obtained by using instead the eraC2t00b function.


------
c2t00b
------

Given:
    TT as a 2-part Julian Date (Note 1)
    UT1 as a 2-part Julian Date (Note 1)
    coordinates of the pole (radians, Note 2)

]   celestial-to-terrestrial matrix (Note 3)

 tta+ttb and uta+utb are Julian Dates,
onvenient way between the arguments uta and
D(UT1)=2450123.7 could be expressed in any of
hers:
   utb
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution is
00 and MJD methods are good compromises
nd convenience.  In the case of uta,utb, the
s best matched to the Earth rotation angle
imum precision is delivered when the uta
 UT1 on the day in question and the utb
 range 0 to 1, or vice versa.
 yp are the coordinates (in radians) of the
te Pole with respect to the International
e System (see IERS Conventions 2003),
eridians to 0 and 90 deg west respectively.
sforms from celestial to terrestrial

3(ERA) * RC2I * [CRS]
RS]
tor in the Geocentric Celestial Reference
a vector in the International Terrestrial
e IERS Conventions 2003), RC2I is the
diate matrix, ERA is the Earth rotation
e polar motion matrix.
 is faster, but slightly less accurate (about
C2t00a function.

al-to-intermediate matrix, IAU 2000B
otation angle, IAU 2000
otion matrix
O-based celestial-to-terrestrial matrix

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
     apportioned in any convenient way between the arguments uta and
     utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
     these ways, among others:

             uta            utb

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution is
     acceptable.  The J2000 and MJD methods are good compromises
     between resolution and convenience.  In the case of uta,utb, the
     date & time method is best matched to the Earth rotation angle
     algorithm used:  maximum precision is delivered when the uta
     argument is for 0hrs UT1 on the day in question and the utb
     argument lies in the range 0 to 1, or vice versa.

  2) The arguments xp and yp are the coordinates (in radians) of the
     Celestial Intermediate Pole with respect to the International
     Terrestrial Reference System (see IERS Conventions 2003),
     measured along the meridians to 0 and 90 deg west respectively.

  3) The matrix rc2t transforms from celestial to terrestrial
     coordinates:

        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]

              = rc2t * [CRS]

     where [CRS] is a vector in the Geocentric Celestial Reference
     System and [TRS] is a vector in the International Terrestrial
     Reference System (see IERS Conventions 2003), RC2I is the
     celestial-to-intermediate matrix, ERA is the Earth rotation
     angle and RPOM is the polar motion matrix.

  4) The present function is faster, but slightly less accurate (about
     1 mas), than the eraC2t00a function.


------
c2t06a
------

Given:
    TT as a 2-part Julian Date (Note 1)
    UT1 as a 2-part Julian Date (Note 1)
    coordinates of the pole (radians, Note 2)

]   celestial-to-terrestrial matrix (Note 3)

 tta+ttb and uta+utb are Julian Dates,
onvenient way between the arguments uta and
D(UT1)=2450123.7 could be expressed in any of
hers:
   utb
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution is
00 and MJD methods are good compromises
nd convenience.  In the case of uta,utb, the
s best matched to the Earth rotation angle
imum precision is delivered when the uta
 UT1 on the day in question and the utb
 range 0 to 1, or vice versa.
 yp are the coordinates (in radians) of the
te Pole with respect to the International
e System (see IERS Conventions 2003),
eridians to 0 and 90 deg west respectively.
sforms from celestial to terrestrial

3(ERA) * RC2I * [CRS]
RS]
tor in the Geocentric Celestial Reference
a vector in the International Terrestrial
e IERS Conventions 2003), RC2I is the
diate matrix, ERA is the Earth rotation
e polar motion matrix.

al-to-intermediate matrix, IAU 2006/2000A
otation angle, IAU 2000
 locator s', IERS 2000
otion matrix
O-based celestial-to-terrestrial matrix

it, G. (eds.), 2004, IERS Conventions (2003),
No. 32, BKG
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
     apportioned in any convenient way between the arguments uta and
     utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
     these ways, among others:

             uta            utb

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution is
     acceptable.  The J2000 and MJD methods are good compromises
     between resolution and convenience.  In the case of uta,utb, the
     date & time method is best matched to the Earth rotation angle
     algorithm used:  maximum precision is delivered when the uta
     argument is for 0hrs UT1 on the day in question and the utb
     argument lies in the range 0 to 1, or vice versa.

  2) The arguments xp and yp are the coordinates (in radians) of the
     Celestial Intermediate Pole with respect to the International
     Terrestrial Reference System (see IERS Conventions 2003),
     measured along the meridians to 0 and 90 deg west respectively.

  3) The matrix rc2t transforms from celestial to terrestrial
     coordinates:

        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]

              = rc2t * [CRS]

     where [CRS] is a vector in the Geocentric Celestial Reference
     System and [TRS] is a vector in the International Terrestrial
     Reference System (see IERS Conventions 2003), RC2I is the
     celestial-to-intermediate matrix, ERA is the Earth rotation
     angle and RPOM is the polar motion matrix.


------
c2tcio
------

Given:
]    celestial-to-intermediate matrix
     Earth rotation angle (radians)
]    polar-motion matrix

]    celestial-to-terrestrial matrix

ucts the rotation matrix that transforms
tial system into vectors in the terrestrial
starting from precomputed components, namely
ates from celestial coordinates to the
the Earth rotation angle and the polar motion
the present function is when generating a
to-terrestrial matrices where only the Earth
es, avoiding the considerable overhead of
ession-nutation more often than necessary to
cy objectives.
ween the arguments is as follows:
3(ERA) * rc2i * [CRS]
RS]
tor in the Geocentric Celestial Reference
a vector in the International Terrestrial
e IERS Conventions 2003).

matrix
around Z-axis
 of two r-matrices

it, G. (eds.), 2004, IERS Conventions (2003),
No. 32, BKG
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) This function constructs the rotation matrix that transforms
     vectors in the celestial system into vectors in the terrestrial
     system.  It does so starting from precomputed components, namely
     the matrix which rotates from celestial coordinates to the
     intermediate frame, the Earth rotation angle and the polar motion
     matrix.  One use of the present function is when generating a
     series of celestial-to-terrestrial matrices where only the Earth
     Rotation Angle changes, avoiding the considerable overhead of
     recomputing the precession-nutation more often than necessary to
     achieve given accuracy objectives.

  2) The relationship between the arguments is as follows:

        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]

              = rc2t * [CRS]

     where [CRS] is a vector in the Geocentric Celestial Reference
     System and [TRS] is a vector in the International Terrestrial
     Reference System (see IERS Conventions 2003).


------
c2teqx
------

Given:
 celestial-to-true matrix
 Greenwich (apparent) Sidereal Time (radians)
 polar-motion matrix

 celestial-to-terrestrial matrix (Note 2)

ucts the rotation matrix that transforms
tial system into vectors in the terrestrial
starting from precomputed components, namely
ates from celestial coordinates to the
inox of date, the Greenwich Apparent Sidereal
otion matrix.  One use of the present function
 series of celestial-to-terrestrial matrices
eal Time changes, avoiding the considerable
ing the precession-nutation more often than
 given accuracy objectives.
ween the arguments is as follows:
3(gst) * rbpn * [CRS]
RS]
tor in the Geocentric Celestial Reference
a vector in the International Terrestrial
e IERS Conventions 2003).

matrix
around Z-axis
 of two r-matrices

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) This function constructs the rotation matrix that transforms
     vectors in the celestial system into vectors in the terrestrial
     system.  It does so starting from precomputed components, namely
     the matrix which rotates from celestial coordinates to the
     true equator and equinox of date, the Greenwich Apparent Sidereal
     Time and the polar motion matrix.  One use of the present function
     is when generating a series of celestial-to-terrestrial matrices
     where only the Sidereal Time changes, avoiding the considerable
     overhead of recomputing the precession-nutation more often than
     necessary to achieve given accuracy objectives.

  2) The relationship between the arguments is as follows:

        [TRS] = rpom * R_3(gst) * rbpn * [CRS]

              = rc2t * [CRS]

     where [CRS] is a vector in the Geocentric Celestial Reference
     System and [TRS] is a vector in the International Terrestrial
     Reference System (see IERS Conventions 2003).


------
c2tpe
------

Given:
     TT as a 2-part Julian Date (Note 1)
     UT1 as a 2-part Julian Date (Note 1)
     nutation (Note 2)
     coordinates of the pole (radians, Note 3)

[3]  celestial-to-terrestrial matrix (Note 4)

 tta+ttb and uta+utb are Julian Dates,
onvenient way between the arguments uta and
D(UT1)=2450123.7 could be expressed in any of
hers:
   utb
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution is
00 and MJD methods are good compromises
nd convenience.  In the case of uta,utb, the
s best matched to the Earth rotation angle
imum precision is delivered when the uta
 UT1 on the day in question and the utb
 range 0 to 1, or vice versa.
sible for providing the nutation components;
e and obliquity, in radians and are with
ox and ecliptic of date.  For high-accuracy
ore nutation should be included as well as
orrections to the position of the CIP.
 yp are the coordinates (in radians) of the
te Pole with respect to the International
e System (see IERS Conventions 2003),
eridians to 0 and 90 deg west respectively.
sforms from celestial to terrestrial

3(GST) * RBPN * [CRS]
RS]
tor in the Geocentric Celestial Reference
a vector in the International Terrestrial
e IERS Conventions 2003), RBPN is the
tion matrix, GST is the Greenwich (apparent)
OM is the polar motion matrix.
es not include "00", This function is in fact
2000 models.

ecession/nutation results, IAU 2000
ch mean sidereal time, IAU 2000
 locator s', IERS 2000
n of the equinoxes, IAU 2000
otion matrix
uinox-based celestial-to-terrestrial matrix

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
     apportioned in any convenient way between the arguments uta and
     utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
     these ways, among others:

             uta            utb

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution is
     acceptable.  The J2000 and MJD methods are good compromises
     between resolution and convenience.  In the case of uta,utb, the
     date & time method is best matched to the Earth rotation angle
     algorithm used:  maximum precision is delivered when the uta
     argument is for 0hrs UT1 on the day in question and the utb
     argument lies in the range 0 to 1, or vice versa.

  2) The caller is responsible for providing the nutation components;
     they are in longitude and obliquity, in radians and are with
     respect to the equinox and ecliptic of date.  For high-accuracy
     applications, free core nutation should be included as well as
     any other relevant corrections to the position of the CIP.

  3) The arguments xp and yp are the coordinates (in radians) of the
     Celestial Intermediate Pole with respect to the International
     Terrestrial Reference System (see IERS Conventions 2003),
     measured along the meridians to 0 and 90 deg west respectively.

  4) The matrix rc2t transforms from celestial to terrestrial
     coordinates:

        [TRS] = RPOM * R_3(GST) * RBPN * [CRS]

              = rc2t * [CRS]

     where [CRS] is a vector in the Geocentric Celestial Reference
     System and [TRS] is a vector in the International Terrestrial
     Reference System (see IERS Conventions 2003), RBPN is the
     bias-precession-nutation matrix, GST is the Greenwich (apparent)
     Sidereal Time and RPOM is the polar motion matrix.

  5) Although its name does not include "00", This function is in fact
     specific to the IAU 2000 models.


------
c2txy
------

Given:
    TT as a 2-part Julian Date (Note 1)
    UT1 as a 2-part Julian Date (Note 1)
    Celestial Intermediate Pole (Note 2)
    coordinates of the pole (radians, Note 3)

]   celestial-to-terrestrial matrix (Note 4)

 tta+ttb and uta+utb are Julian Dates,
onvenient way between the arguments uta and
D(UT1)=2450123.7 could be expressed in any o
hers:
   utb
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution is
00 and MJD methods are good compromises
nd convenience.  In the case of uta,utb, the
s best matched to the Earth rotation angle
imum precision is delivered when the uta
 UT1 on the day in question and the utb
 range 0 to 1, or vice versa.
ediate Pole coordinates are the x,y
it vector in the Geocentric Celestial

 yp are the coordinates (in radians) of the
te Pole with respect to the International
e System (see IERS Conventions 2003),
eridians to 0 and 90 deg west respectively.
sforms from celestial to terrestrial

3(ERA) * RC2I * [CRS]
RS]
tor in the Geocentric Celestial Reference
a vector in the International Terrestrial
e IERS Conventions 2003), ERA is the Earth
POM is the polar motion matrix.
es not include "00", This function is in fact
2000 models.

al-to-intermediate matrix, given X,Y
otation angle, IAU 2000
 locator s', IERS 2000
otion matrix
O-based celestial-to-terrestrial matrix

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
     apportioned in any convenient way between the arguments uta and
     utb.  For example, JD(UT1)=2450123.7 could be expressed in any o
     these ways, among others:

             uta            utb

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution is
     acceptable.  The J2000 and MJD methods are good compromises
     between resolution and convenience.  In the case of uta,utb, the
     date & time method is best matched to the Earth rotation angle
     algorithm used:  maximum precision is delivered when the uta
     argument is for 0hrs UT1 on the day in question and the utb
     argument lies in the range 0 to 1, or vice versa.

  2) The Celestial Intermediate Pole coordinates are the x,y
     components of the unit vector in the Geocentric Celestial
     Reference System.

  3) The arguments xp and yp are the coordinates (in radians) of the
     Celestial Intermediate Pole with respect to the International
     Terrestrial Reference System (see IERS Conventions 2003),
     measured along the meridians to 0 and 90 deg west respectively.

  4) The matrix rc2t transforms from celestial to terrestrial
     coordinates:

        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]

              = rc2t * [CRS]

     where [CRS] is a vector in the Geocentric Celestial Reference
     System and [TRS] is a vector in the International Terrestrial
     Reference System (see IERS Conventions 2003), ERA is the Earth
     Rotation Angle and RPOM is the polar motion matrix.

  5) Although its name does not include "00", This function is in fact
     specific to the IAU 2000 models.


------
cal2jd
------

Given:
ar, month, day in Gregorian calendar (Note 1)

D zero-point: always 2400000.5
dified Julian Date for 0 hrs
e):
atus:
  0 = OK
 -1 = bad year   (Note 3: JD not computed)
 -2 = bad month  (JD not computed)
 -3 = bad day    (JD computed)

s valid from -4800 March 1, but this
ts dates before -4799 January 1.
eturned in two pieces, in the usual ERFA
igned to preserve time resolution.  The
able as a single number by adding djm0 and

nversion is from the "Proleptic Gregorian
nt is taken of the date(s) of adoption of
ar, nor is the AD/BC numbering convention


nt to the Astronomical Almanac,
n (ed), University Science Books (1992),
.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
cp
------

Given:
   p-vector to be copied

   copy
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
cpv
------

Given:
   position/velocity vector to be copied

   copy

vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
cr
------

Given:
]    r-matrix to be copied

]    copy

vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
d2dtf
------

Given:
me scale ID (Note 1)
solution (Note 2)
me as a 2-part Julian Date (Notes 3,4)

ar, month, day in Gregorian calendar (Note 5)
urs, minutes, seconds, fraction (Note 1)
e):
atus: +1 = dubious year (Note 5)
       0 = OK
      -1 = unacceptable date (Note 6)

 time scale.  Only the value "UTC" (in upper
, and enables handling of leap seconds (see

 decimal places in the seconds field, and can
l as positive values, such as:
on
0
0
0
0
1
0.1
0.01
0.001
orm dependent, but a safe range is -5 to +9.
, apportioned in any convenient way between
or example where d1 is the Julian Day Number
on of a day.  In the case of UTC, where the
atical, special conventions apply:  see the

sly represent UTC during a leap second unless
 taken.  The ERFA internal convention is that
resents UTC days whether the length is 86399,
conds.  In the 1960-1972 era there were
ther direction) each time the linear UTC(TAI)
ed, and these "mini-leaps" are also included
on.
dubious year" flags UTCs that predate the
time scale or that are too far in the future
eraDat for further details.
ions and limitations, see eraCal2jd.

regorian calendar
se days to hms
T) = TAI-UTC
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) scale identifies the time scale.  Only the value "UTC" (in upper
     case) is significant, and enables handling of leap seconds (see
     Note 4).

  2) ndp is the number of decimal places in the seconds field, and can
     have negative as well as positive values, such as:

     ndp         resolution
     -4            1 00 00
     -3            0 10 00
     -2            0 01 00
     -1            0 00 10
      0            0 00 01
      1            0 00 00.1
      2            0 00 00.01
      3            0 00 00.001

     The limits are platform dependent, but a safe range is -5 to +9.

  3) d1+d2 is Julian Date, apportioned in any convenient way between
     the two arguments, for example where d1 is the Julian Day Number
     and d2 is the fraction of a day.  In the case of UTC, where the
     use of JD is problematical, special conventions apply:  see the
     next note.

  4) JD cannot unambiguously represent UTC during a leap second unless
     special measures are taken.  The ERFA internal convention is that
     the quasi-JD day represents UTC days whether the length is 86399,
     86400 or 86401 SI seconds.  In the 1960-1972 era there were
     smaller jumps (in either direction) each time the linear UTC(TAI)
     expression was changed, and these "mini-leaps" are also included
     in the ERFA convention.

  5) The warning status "dubious year" flags UTCs that predate the
     introduction of the time scale or that are too far in the future
     to be trusted.  See eraDat for further details.

  6) For calendar conventions and limitations, see eraCal2jd.


------
d2tf
------

Given:
lution (Note 1)
rval in days

or '-'
s, minutes, seconds, fraction

interpreted as follows:
on
0
0
0
0
0
0
0
0
1
0.1
0.01
0.001
0.000...
 useful value for ndp is determined by the
rmat of double on the target platform, and
ing ihmsf[3].  On a typical platform, for
available floating-point precision might
.  However, the practical limit is typically
pacity of a 32-bit int, or ndp=4 if int is

f days may exceed 1.0.  In cases where it
o the caller to test for and handle the
ery nearly 1.0 and rounds up to 24 hours,
[0]=24 and setting ihmsf[0-3] to zero.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
dat
------

Given:
  year (Notes 1 and 2)
  month (Note 2)
  day (Notes 2 and 3)
  fraction of day (Note 4)

minus UTC, seconds
e):
us (Note 5):
= dubious year (Note 1)
= OK
= bad year
= bad month
= bad day (Note 3)
= bad fraction (Note 4)
= internal error (Note 5)

nuary 1.0 (JD 2436934.5) and it is improper
 with an earlier date.  If this is attempted,
ether with a warning status.
 cannot, in principle, be predicted in
check for dates beyond the valid range is
d against gross errors, a year five or more
ar of the present function (see the constant
ubious.  In this case a warning status is
ult is computed in the normal way.
nd too-late years, the warning status is +1.
m the error status -1, which signifies a year
ld not be computed.
e is for a day which ends with a leap second,
turned is for the period leading up to the
 date is for a day which begins as a leap
-TAI returned is for the period following the

be in the normal calendar range, for example
il.  The "almanac" convention of allowing
y 0 and December 32 is not supported in this
o avoid confusion near leap seconds.
is used only for dates before the
 seconds, the first of which occurred at the
tested for validity (0 to 1 is the valid
sed;  if invalid, zero is used and status -4
ny applications, setting fd to zero is
ulting error is always less than 3 ms (and
).
urned in the case where there are multiple
 first error detected.  For example, if the
 and 32 respectively, status -2 (bad month)
he "internal error" status refers to a
ble but causes some compilers to issue a

id result is not available, zero is returned.

January 1 onwards, the expressions from the
.navy.mil/ser7/tai-utc.dat are used.
1961 January 1 is taken from 2.58.1 (p87) of
 Supplement.

an calendar to JD
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
     to call the function with an earlier date.  If this is attempted,
     zero is returned together with a warning status.

     Because leap seconds cannot, in principle, be predicted in
     advance, a reliable check for dates beyond the valid range is
     impossible.  To guard against gross errors, a year five or more
     after the release year of the present function (see the constant
     IYV) is considered dubious.  In this case a warning status is
     returned but the result is computed in the normal way.

     For both too-early and too-late years, the warning status is +1.
     This is distinct from the error status -1, which signifies a year
     so early that JD could not be computed.

  2) If the specified date is for a day which ends with a leap second,
     the UTC-TAI value returned is for the period leading up to the
     leap second.  If the date is for a day which begins as a leap
     second ends, the UTC-TAI returned is for the period following the
     leap second.

  3) The day number must be in the normal calendar range, for example
     1 through 30 for April.  The "almanac" convention of allowing
     such dates as January 0 and December 32 is not supported in this
     function, in order to avoid confusion near leap seconds.

  4) The fraction of day is used only for dates before the
     introduction of leap seconds, the first of which occurred at the
     end of 1971.  It is tested for validity (0 to 1 is the valid
     range) even if not used;  if invalid, zero is used and status -4
     is returned.  For many applications, setting fd to zero is
     acceptable;  the resulting error is always less than 3 ms (and
     occurs only pre-1972).

  5) The status value returned in the case where there are multiple
     errors refers to the first error detected.  For example, if the
     month and day are 13 and 32 respectively, status -2 (bad month)
     will be returned.  The "internal error" status refers to a
     case that is impossible but causes some compilers to issue a
     warning.

  6) In cases where a valid result is not available, zero is returned.

  References:

  1) For dates from 1961 January 1 onwards, the expressions from the
     file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.

  2) The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
     the 1992 Explanatory Supplement.


References:

  1) For dates from 1961 January 1 onwards, the expressions from the
     file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.

  2) The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
     the 1992 Explanatory Supplement.

  Called:
     eraCal2jd    Gregorian calendar to JD

------
dtdb
------

Given:
  date, TDB (Notes 1-3)
  universal time (UT1, fraction of one day)
  longitude (east positive, radians)
  distance from Earth spin axis (km)
  distance north of equatorial plane (km)
e):
  TDB-TT (seconds)

 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
, formally, barycentric dynamical time (TDB),
mical time (TT) can be used with no practical
cy of the prediction.
s a coordinate time that is realized as an
om International Atomic Time, TAI.  TT is a
sformation of geocentric coordinate time TCG,
ale for the Geocentric Celestial Reference

time, and is a specific linear transformation
inate time TCB, which is the time scale for
stial Reference System, BCRS.
CB depends on the masses and positions of the
system and the velocity of the Earth.  It is
difference, the residual being of a periodic
er, which is modeled by the present function,
nual) sinusoidal term of amplitude
6 seconds, plus planetary terms up to about
 lunar and diurnal terms up to 2 microseconds.
rom the changing transverse Doppler effect
d-shift as the observer (on the Earth's
 variations in speed (with respect to the
nal potential.
as the same as TCB but with a rate adjustment
TT, which is convenient for many applications.
ssive attempts to define TDB is set out in
 by the IAU General Assembly in 2006, which
TCB) transformation that is consistent with
ystem ephemerides.  Future ephemerides will
ed transformations between TCG and TCB, which
near drift between TDB and TT;  however, any
ly to exceed 1 nanosecond per century.
T model used in the present function is that of
 (1990), in its full form.  It was originally
 (private communications with P.T.Wallace,
ubroutine.  The present C function contains an
irhead code.  The numerical results are
ed by the changes, the differences with
ead & Bretagnon original being at the 1e-20 s

 of the model is from Moyer (1981) and
fundamental arguments adapted from
It is an approximation to the expression
), where v is the barycentric velocity of
geocentric position of the observer and
ght.
for u and v, the topocentric part of the
ed, and the function will return the Fairhead
lone.
1950-2050, the absolute accuracy is better
ds relative to time ephemerides obtained by
egrations based on the JPL DE405 solar system

that the present function is merely a model,
ntegration of solar-system ephemerides is the
r predicting the relationship between TCG and
n TT and TDB.

agnon, P., Astron.Astrophys., 229, 240-247

3.
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
h., 23, 33 (1981).
ial Astrometry, Adam Hilger (1983).
al., Explanatory Supplement to the
, Chapter 2, University Science Books (1992).
on, P., Chapront, J., Chapront-Touze, M.,
, J., Astron.Astrophys., 282, 663-683 (1994).
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

     Although the date is, formally, barycentric dynamical time (TDB),
     the terrestrial dynamical time (TT) can be used with no practical
     effect on the accuracy of the prediction.

  2) TT can be regarded as a coordinate time that is realized as an
     offset of 32.184s from International Atomic Time, TAI.  TT is a
     specific linear transformation of geocentric coordinate time TCG,
     which is the time scale for the Geocentric Celestial Reference
     System, GCRS.

  3) TDB is a coordinate time, and is a specific linear transformation
     of barycentric coordinate time TCB, which is the time scale for
     the Barycentric Celestial Reference System, BCRS.

  4) The difference TCG-TCB depends on the masses and positions of the
     bodies of the solar system and the velocity of the Earth.  It is
     dominated by a rate difference, the residual being of a periodic
     character.  The latter, which is modeled by the present function,
     comprises a main (annual) sinusoidal term of amplitude
     approximately 0.00166 seconds, plus planetary terms up to about
     20 microseconds, and lunar and diurnal terms up to 2 microseconds.
     These effects come from the changing transverse Doppler effect
     and gravitational red-shift as the observer (on the Earth's
     surface) experiences variations in speed (with respect to the
     BCRS) and gravitational potential.

  5) TDB can be regarded as the same as TCB but with a rate adjustment
     to keep it close to TT, which is convenient for many applications.
     The history of successive attempts to define TDB is set out in
     Resolution 3 adopted by the IAU General Assembly in 2006, which
     defines a fixed TDB(TCB) transformation that is consistent with
     contemporary solar-system ephemerides.  Future ephemerides will
     imply slightly changed transformations between TCG and TCB, which
     could introduce a linear drift between TDB and TT;  however, any
     such drift is unlikely to exceed 1 nanosecond per century.

  6) The geocentric TDB-TT model used in the present function is that of
     Fairhead & Bretagnon (1990), in its full form.  It was originally
     supplied by Fairhead (private communications with P.T.Wallace,
     1990) as a Fortran subroutine.  The present C function contains an
     adaptation of the Fairhead code.  The numerical results are
     essentially unaffected by the changes, the differences with
     respect to the Fairhead & Bretagnon original being at the 1e-20 s
     level.

     The topocentric part of the model is from Moyer (1981) and
     Murray (1983), with fundamental arguments adapted from
     Simon et al. 1994.  It is an approximation to the expression
     ( v / c ) . ( r / c ), where v is the barycentric velocity of
     the Earth, r is the geocentric position of the observer and
     c is the speed of light.

     By supplying zeroes for u and v, the topocentric part of the
     model can be nullified, and the function will return the Fairhead
     & Bretagnon result alone.

  7) During the interval 1950-2050, the absolute accuracy is better
     than +/- 3 nanoseconds relative to time ephemerides obtained by
     direct numerical integrations based on the JPL DE405 solar system
     ephemeris.

  8) It must be stressed that the present function is merely a model,
     and that numerical integration of solar-system ephemerides is the
     definitive method for predicting the relationship between TCG and
     TCB and hence between TT and TDB.


References:

     Fairhead, L., & Bretagnon, P., Astron.Astrophys., 229, 240-247
     (1990).

     IAU 2006 Resolution 3.

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Moyer, T.D., Cel.Mech., 23, 33 (1981).

     Murray, C.A., Vectorial Astrometry, Adam Hilger (1983).

     Seidelmann, P.K. et al., Explanatory Supplement to the
     Astronomical Almanac, Chapter 2, University Science Books (1992).

     Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G. & Laskar, J., Astron.Astrophys., 282, 663-683 (1994).

------
dtf2d
------

Given:
me scale ID (Note 1)
ar, month, day in Gregorian calendar (Note 2)
ur, minute
conds

part Julian Date (Notes 3,4)
e):
atus: +3 = both of next two
      +2 = time is after end of day (Note 5)
      +1 = dubious year (Note 6)
       0 = OK
      -1 = bad year
      -2 = bad month
      -3 = bad day
      -4 = bad hour
      -5 = bad minute
      -6 = bad second (<0)

 time scale.  Only the value "UTC" (in upper
, and enables handling of leap seconds (see

ions and limitations, see eraCal2jd.
ts, d1+d2, is Julian Date, where normally d1
mber and d2 is the fraction of a day.  In the
he use of JD is problematical, special
see the next note.
sly represent UTC during a leap second unless
 taken.  The ERFA internal convention is that
resents UTC days whether the length is 86399,
conds.  In the 1960-1972 era there were
ther direction) each time the linear UTC(TAI)
ed, and these "mini-leaps" are also included
on.
time is after end of day" usually means that
greater than 60.0.  However, in a day ending
 limit changes to 61.0 (or 59.0 in the case
econd).
dubious year" flags UTCs that predate the
time scale or that are too far in the future
eraDat for further details.
continuous and regular time scales (TAI, TT,
 the result d1+d2 a Julian Date, strictly
her cases (UT1 and UTC) the result must be
tion;  in particular the difference between
not be interpreted as a precise time


an calendar to JD
T) = TAI-UTC
regorian calendar
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) scale identifies the time scale.  Only the value "UTC" (in upper
     case) is significant, and enables handling of leap seconds (see
     Note 4).

  2) For calendar conventions and limitations, see eraCal2jd.

  3) The sum of the results, d1+d2, is Julian Date, where normally d1
     is the Julian Day Number and d2 is the fraction of a day.  In the
     case of UTC, where the use of JD is problematical, special
     conventions apply:  see the next note.

  4) JD cannot unambiguously represent UTC during a leap second unless
     special measures are taken.  The ERFA internal convention is that
     the quasi-JD day represents UTC days whether the length is 86399,
     86400 or 86401 SI seconds.  In the 1960-1972 era there were
     smaller jumps (in either direction) each time the linear UTC(TAI)
     expression was changed, and these "mini-leaps" are also included
     in the ERFA convention.

  5) The warning status "time is after end of day" usually means that
     the sec argument is greater than 60.0.  However, in a day ending
     in a leap second the limit changes to 61.0 (or 59.0 in the case
     of a negative leap second).

  6) The warning status "dubious year" flags UTCs that predate the
     introduction of the time scale or that are too far in the future
     to be trusted.  See eraDat for further details.

  7) Only in the case of continuous and regular time scales (TAI, TT,
     TCG, TCB and TDB) is the result d1+d2 a Julian Date, strictly
     speaking.  In the other cases (UT1 and UTC) the result must be
     used with circumspection;  in particular the difference between
     two such results cannot be interpreted as a precise time
     interval.


------
eceq06
------

Given:
T as a 2-part Julian date (Note 1)
cliptic longitude and latitude (radians)

CRS right ascension and declination (radians)
te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
ade about whether the coordinates represent
 astrometric effects such as parallax or

s approximately that from ecliptic longitude
quinox and ecliptic of date) to mean J2000.0
declination, with only frame bias (always
 disturb this classical picture.

al coordinates to unit vector
 to ecliptic rotation matrix, IAU 2006
 of transpose of r-matrix and p-vector
ctor to spherical coordinates
ze angle into range 0 to 2pi
ze angle into range +/- pi
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
ecm06
------

Given:
        TT as a 2-part Julian date (Note 1)

3][3]   ICRS to ecliptic rotation matrix

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
 sense
S,
ctor with respect to ICRS right ascension
 and E_ep is the same vector with respect to
tic and equinox of date.
tor, merely a direction, typically of unit
ound to any particular spatial origin, such
 SSB.  No assumptions are made about whether
ght and embodies astrometric effects such as
on.  The transformation is approximately that
 right ascension and declination and ecliptic
de, with only frame bias (always less than
his classical picture.

liquity, IAU 2006
ix, IAU 2006
ize r-matrix to identity
around X-axis
 of two r-matrices
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  1) The matrix is in the sense

        E_ep = rm x P_ICRS,

     where P_ICRS is a vector with respect to ICRS right ascension
     and declination axes and E_ep is the same vector with respect to
     the (inertial) ecliptic and equinox of date.

  2) P_ICRS is a free vector, merely a direction, typically of unit
     magnitude, and not bound to any particular spatial origin, such
     as the Earth, Sun or SSB.  No assumptions are made about whether
     it represents starlight and embodies astrometric effects such as
     parallax or aberration.  The transformation is approximately that
     between mean J2000.0 right ascension and declination and ecliptic
     longitude and latitude, with only frame bias (always less than
     25 mas) to disturb this classical picture.


------
ee00
------

Given:
   TT as a 2-part Julian Date (Note 1)
   mean obliquity (Note 2)
   nutation in longitude (Note 3)
e):
   equation of the equinoxes (Note 4)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
dians, is mean of date.
 in radians, operates in the following sense:
t ST = GMST + equation of the equinoxes
ible with the IAU 2000 resolutions.  For
 IERS Conventions 2003 and Capitaine et al.


n of the equinoxes complementary terms

ce, P.T. and McCarthy, D.D., "Expressions to
00 definition of UT1", Astronomy &
135-1149 (2003)
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The obliquity, in radians, is mean of date.

  3) The result, which is in radians, operates in the following sense:

        Greenwich apparent ST = GMST + equation of the equinoxes

  4) The result is compatible with the IAU 2000 resolutions.  For
     further details, see IERS Conventions 2003 and Capitaine et al.
     (2002).


References:

     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
     implement the IAU 2000 definition of UT1", Astronomy &
     Astrophysics, 406, 1135-1149 (2003)

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

------
ee00a
------

Given:
   TT as a 2-part Julian Date (Note 1)
e):
   equation of the equinoxes (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
 in radians, operates in the following sense:
t ST = GMST + equation of the equinoxes
ible with the IAU 2000 resolutions.  For
 IERS Conventions 2003 and Capitaine et al.


0 precession adjustments
liquity, IAU 1980
n, IAU 2000A
n of the equinoxes, IAU 2000

ce, P.T. and McCarthy, D.D., "Expressions to
00 definition of UT1", Astronomy &
135-1149 (2003).
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004).
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The result, which is in radians, operates in the following sense:

        Greenwich apparent ST = GMST + equation of the equinoxes

  3) The result is compatible with the IAU 2000 resolutions.  For
     further details, see IERS Conventions 2003 and Capitaine et al.
     (2002).


References:

     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
     implement the IAU 2000 definition of UT1", Astronomy &
     Astrophysics, 406, 1135-1149 (2003).

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004).

------
ee00b
------

Given:
   TT as a 2-part Julian Date (Note 1)
e):
   equation of the equinoxes (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
 in radians, operates in the following sense:
t ST = GMST + equation of the equinoxes
ible with the IAU 2000 resolutions except
en compromised for the sake of speed.  For
 McCarthy & Luzum (2001), IERS Conventions
t al. (2003).

0 precession adjustments
liquity, IAU 1980
n, IAU 2000B
n of the equinoxes, IAU 2000

ce, P.T. and McCarthy, D.D., "Expressions to
00 definition of UT1", Astronomy &
135-1149 (2003)
um, B.J., "An abridged model of the
of the celestial pole", Celestial Mechanics &
 85, 37-49 (2003)
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The result, which is in radians, operates in the following sense:

        Greenwich apparent ST = GMST + equation of the equinoxes

  3) The result is compatible with the IAU 2000 resolutions except
     that accuracy has been compromised for the sake of speed.  For
     further details, see McCarthy & Luzum (2001), IERS Conventions
     2003 and Capitaine et al. (2003).


References:

     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
     implement the IAU 2000 definition of UT1", Astronomy &
     Astrophysics, 406, 1135-1149 (2003)

     McCarthy, D.D. & Luzum, B.J., "An abridged model of the
     precession-nutation of the celestial pole", Celestial Mechanics &
     Dynamical Astronomy, 85, 37-49 (2003)

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

------
ee06a
------

Given:
   TT as a 2-part Julian Date (Note 1)
e):
   equation of the equinoxes (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
 in radians, operates in the following sense:
t ST = GMST + equation of the equinoxes

ze angle into range +/- pi
ch apparent sidereal time, IAU 2006/2000A
ch mean sidereal time, IAU 2006

it, G. (eds.), 2004, IERS Conventions (2003),
No. 32, BKG
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The result, which is in radians, operates in the following sense:

        Greenwich apparent ST = GMST + equation of the equinoxes


------
eect00
------

Given:
  TT as a 2-part Julian Date (Note 1)
e):
  complementary terms (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
erms" are part of the equation of the
sically the difference between apparent and



ps)
tation in longitude and eps is the obliquity
f the rotation of the Earth were constant in
e classical formulation would lead to
ies in the UT1 timescale traceable to side-
n-nutation.  In order to eliminate these
omplementary terms" were introduced in 1994
 effect from 1997 (Capitaine and Gontier,

+ EE
omplementary terms are included as part of
equinoxes rather than as part of the mean
 slightly compromises the "geometrical"
an sidereal time but is otherwise

 computes CT in the above expression,
2000 resolutions (Capitaine et al., 2002, and
3).

omaly of the Moon
omaly of the Sun
gument of the latitude of the Moon
ongation of the Moon from the Sun
ngitude of the Moon's ascending node
ngitude of Venus
ngitude of Earth
 accumulated precession in longitude

ier, A.-M., Astron. Astrophys., 275,

ce, P.T. and McCarthy, D.D., "Expressions to
00 definition of UT1", Astronomy &
135-1149 (2003)
ecommendation 3 (1994)
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The "complementary terms" are part of the equation of the
     equinoxes (EE), classically the difference between apparent and
     mean Sidereal Time:

        GAST = GMST + EE

     with:

        EE = dpsi * cos(eps)

     where dpsi is the nutation in longitude and eps is the obliquity
     of date.  However, if the rotation of the Earth were constant in
     an inertial frame the classical formulation would lead to
     apparent irregularities in the UT1 timescale traceable to side-
     effects of precession-nutation.  In order to eliminate these
     effects from UT1, "complementary terms" were introduced in 1994
     (IAU, 1994) and took effect from 1997 (Capitaine and Gontier,
     1993):

        GAST = GMST + CT + EE

     By convention, the complementary terms are included as part of
     the equation of the equinoxes rather than as part of the mean
     Sidereal Time.  This slightly compromises the "geometrical"
     interpretation of mean sidereal time but is otherwise
     inconsequential.

     The present function computes CT in the above expression,
     compatible with IAU 2000 resolutions (Capitaine et al., 2002, and
     IERS Conventions 2003).


References:

     Capitaine, N. & Gontier, A.-M., Astron. Astrophys., 275,
     645-650 (1993)

     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
     implement the IAU 2000 definition of UT1", Astronomy &
     Astrophysics, 406, 1135-1149 (2003)

     IAU Resolution C7, Recommendation 3 (1994)

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

------
eform
------

Given:
ipsoid identifier (Note 1)

atorial radius (meters, Note 2)
ttening (Note 2)
e):
tus:  0 = OK
     -1 = illegal identifier (Note 3)

a number that specifies the choice of
  The following are supported:




ignificance outside the ERFA software.  For
 ERFA_WGS84 etc. are defined in erfam.h.
ters are returned in the form of equatorial
 and flattening (f).  The latter is a number
 around 1/298.
n unsupported n value is supplied, zero a and
ell as error status.

e World Geodetic System 1984, National
Agency Technical Report 8350.2, Third

odesique 66-2, 187 (1992).
fense World Geodetic System 1972, World
ittee, May 1974.
nt to the Astronomical Almanac,
n (ed), University Science Books (1992),

, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The identifier n is a number that specifies the choice of
     reference ellipsoid.  The following are supported:

        n    ellipsoid

        1     ERFA_WGS84
        2     ERFA_GRS80
        3     ERFA_WGS72

     The n value has no significance outside the ERFA software.  For
     convenience, symbols ERFA_WGS84 etc. are defined in erfam.h.

  2) The ellipsoid parameters are returned in the form of equatorial
     radius in meters (a) and flattening (f).  The latter is a number
     around 0.00335, i.e. around 1/298.

  3) For the case where an unsupported n value is supplied, zero a and
     f are returned, as well as error status.


References:

     Department of Defense World Geodetic System 1984, National
     Imagery and Mapping Agency Technical Report 8350.2, Third
     Edition, p3-2.

     Moritz, H., Bull. Geodesique 66-2, 187 (1992).

     The Department of Defense World Geodetic System 1972, World
     Geodetic System Committee, May 1974.

     Explanatory Supplement to the Astronomical Almanac,
     P. Kenneth Seidelmann (ed), University Science Books (1992),
     p220.

------
eo06a
------

Given:
   TT as a 2-part Julian Date (Note 1)
e):
   equation of the origins in radians

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
origins is the distance between the true
stial intermediate origin and, equivalently,
en Earth rotation angle and Greenwich
me (ERA-GST).  It comprises the precession
ight ascension plus the equation of the
 the small correction terms).

al NPB matrix, IAU 2006/2000A
 CIP X,Y coordinates from NPB matrix
 locator s, given X,Y, IAU 2006
n of the origins, given NPB matrix and s

ace, P.T., 2006, Astron.Astrophys. 450, 855
taine, N., 2006, Astron.Astrophys. 459, 981
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The equation of the origins is the distance between the true
     equinox and the celestial intermediate origin and, equivalently,
     the difference between Earth rotation angle and Greenwich
     apparent sidereal time (ERA-GST).  It comprises the precession
     (since J2000.0) in right ascension plus the equation of the
     equinoxes (including the small correction terms).


References:

     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981

------
eors
------

Given:
classical nutation x precession x bias matrix
the quantity s (the CIO locator)
e):
the equation of the origins in radians.

 origins is the distance between the true
estial intermediate origin and, equivalently,
een Earth rotation angle and Greenwich
ime (ERA-GST).  It comprises the precession
right ascension plus the equation of the
g the small correction terms).
om Wallace & Capitaine (2006).

ace, P.T., 2006, Astron.Astrophys. 450, 855
ine, N., 2006, Astron.Astrophys. 459, 981
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
epb
------

Given:
  Julian Date (see note)
e):
  Besselian Epoch.

upplied in two pieces, in the usual ERFA
igned to preserve time resolution.  The
able as a single number by adding dj1 and
solution is achieved if dj1 is 2451545.0


Astron.Astrophys., 73, 282.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
epb2jd
------

Given:
esselian Epoch (e.g. 1957.3)

JD zero-point: always 2400000.5
odified Julian Date

eturned in two pieces, in the usual ERFA
igned to preserve time resolution.  The
able as a single number by adding djm0 and


Astron.Astrophys. 73, 282.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
epj
------

Given:
  Julian Date (see note)
e):
  Julian Epoch

upplied in two pieces, in the usual ERFA
igned to preserve time resolution.  The
able as a single number by adding dj1 and
solution is achieved if dj1 is 2451545.0


Astron.Astrophys. 73, 282.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
epj2jd
------

Given:
ulian Epoch (e.g. 1996.8)

JD zero-point: always 2400000.5
odified Julian Date

eturned in two pieces, in the usual ERFA
igned to preserve time resolution.  The
able as a single number by adding djm0 and


Astron.Astrophys. 73, 282.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
epv00
------

Given:
       TDB date (Note 1)

2][3]  heliocentric Earth position/velocity
2][3]  barycentric Earth position/velocity
e):
       status: 0 = OK
              +1 = warning: date outside
                   the range 1900-2100 AD

ate2 is a Julian Date, apportioned in any
en the two arguments.  For example,
uld be expressed in any of these ways, among

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 method is best matched to the way the
internally and will deliver the optimum
 method and the date & time methods are both
ween resolution and convenience.  However,
result is more likely to be limited by the
n the way the date has been expressed.
instead of TDB in most applications.
s pvh and pvb contain the following:
  }
  } heliocentric position, AU
  }
  }
  } heliocentric velocity, AU/d
  }
  }
  } barycentric position, AU
  }
  }
  } barycentric velocity, AU/d
  }
 respect to the Barycentric Celestial
he time unit is one day in TDB.
MPLIFIED SOLUTION from the planetary theory
, P. Bretagnon, 2001, Celes. Mechanics &
4, 205-213) and is an adaptation of original
d by P. Bretagnon (private comm., 2000).
 time span 1900-2100 with this simplified
 DE405 ephemeris give the following results:
       RMS    max

ror    3.7   11.2   km
ror    1.4    5.0   mm/s

ror    4.6   13.4   km
ror    1.4    4.9   mm/s
 JPL DE406 ephemeris show that by 1800 and
rors are approximately double their 1900-2100
500 the deterioration is a factor of 10 and
actor of 60.  The velocity accuracy falls off
ate.
 use the same array for pvh and pvb, which
ycentric values.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
eqec06
------

Given:
T as a 2-part Julian date (Note 1)
CRS right ascension and declination (radians)

cliptic longitude and latitude (radians)
te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
ade about whether the coordinates represent
 astrometric effects such as parallax or

s approximately that from mean J2000.0 right
ation to ecliptic longitude and latitude
liptic of date), with only frame bias (always
 disturb this classical picture.

al coordinates to unit vector
 to ecliptic rotation matrix, IAU 2006
 of r-matrix and p-vector
ctor to spherical coordinates
ze angle into range 0 to 2pi
ze angle into range +/- pi
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
eqeq94
------

Given:
     TDB date (Note 1)
e):
     equation of the equinoxes (Note 2)

 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
 in radians, operates in the following sense:
t ST = GMST + equation of the equinoxes

ze angle into range +/- pi
n, IAU 1980
liquity, IAU 1980

ecommendation 3 (1994).
ier, A.-M., 1993, Astron. Astrophys., 275,

, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The result, which is in radians, operates in the following sense:

        Greenwich apparent ST = GMST + equation of the equinoxes


References:

     IAU Resolution C7, Recommendation 3 (1994).

     Capitaine, N. & Gontier, A.-M., 1993, Astron. Astrophys., 275,
     645-650.

------
era00
------

Given:
UT1 as a 2-part Julian Date (see note)
e):
Earth rotation angle (radians), range 0-2pi

 is a Julian Date, apportioned in any
en the arguments dj1 and dj2.  For example,
uld be expressed in any of these ways,

   dj2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 and MJD methods are good compromises
nd convenience.  The date & time method is
algorithm used:  maximum precision is
j1 argument is for 0hrs UT1 on the day in
 argument lies in the range 0 to 1, or vice

pted from Expression 22 of Capitaine et al.
ment has been expressed in days directly,
sion, integer contributions have been
e formulation is given in IERS Conventions
 14.

ze angle into range 0 to 2pi

 B. and McCarthy D.D, 2000, Astron.
-405.
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The UT1 date dj1+dj2 is a Julian Date, apportioned in any
     convenient way between the arguments dj1 and dj2.  For example,
     JD(UT1)=2450123.7 could be expressed in any of these ways,
     among others:

             dj1            dj2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 and MJD methods are good compromises
     between resolution and convenience.  The date & time method is
     best matched to the algorithm used:  maximum precision is
     delivered when the dj1 argument is for 0hrs UT1 on the day in
     question and the dj2 argument lies in the range 0 to 1, or vice
     versa.

  2) The algorithm is adapted from Expression 22 of Capitaine et al.
     2000.  The time argument has been expressed in days directly,
     and, to retain precision, integer contributions have been
     eliminated.  The same formulation is given in IERS Conventions
     (2003), Chap. 5, Eq. 14.


References:

     Capitaine N., Guinot B. and McCarthy D.D, 2000, Astron.
     Astrophys., 355, 398-405.

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

------
fad03
------

Given:
 Julian centuries since J2000.0 (Note 1)
e):
adians (Note 2)

 TDB, it is usually more convenient to use
ignificant difference.
is as adopted in IERS Conventions (2003) and
 (1994).

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
non, P., Chapront, J., Chapront-Touze, M.,
 J. 1994, Astron.Astrophys. 282, 663-683
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Though t is strictly TDB, it is usually more convenient to use
     TT, which makes no significant difference.

  2) The expression used is as adopted in IERS Conventions (2003) and
     is from Simon et al. (1994).


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

------
fae03
------

Given:
 Julian centuries since J2000.0 (Note 1)
e):
 longitude of Earth, radians (Note 2)

 TDB, it is usually more convenient to use
ignificant difference.
is as adopted in IERS Conventions (2003) and
t al. (1999) after Simon et al. (1994).

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
non, P., Chapront, J., Chapront-Touze, M.,
 J. 1994, Astron.Astrophys. 282, 663-683
 B., Kinoshita, H., Folgueira, M. 1999,
p.Ser. 135, 111
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Though t is strictly TDB, it is usually more convenient to use
     TT, which makes no significant difference.

  2) The expression used is as adopted in IERS Conventions (2003) and
     comes from Souchay et al. (1999) after Simon et al. (1994).


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     Astron.Astrophys.Supp.Ser. 135, 111

------
faf03
------

Given:
 Julian centuries since J2000.0 (Note 1)
e):
adians (Note 2)

 TDB, it is usually more convenient to use
ignificant difference.
is as adopted in IERS Conventions (2003) and
 (1994).

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
non, P., Chapront, J., Chapront-Touze, M.,
 J. 1994, Astron.Astrophys. 282, 663-683
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Though t is strictly TDB, it is usually more convenient to use
     TT, which makes no significant difference.

  2) The expression used is as adopted in IERS Conventions (2003) and
     is from Simon et al. (1994).


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

------
faju03
------

Given:
 Julian centuries since J2000.0 (Note 1)
e):
 longitude of Jupiter, radians (Note 2)

 TDB, it is usually more convenient to use
ignificant difference.
is as adopted in IERS Conventions (2003) and
t al. (1999) after Simon et al. (1994).

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
non, P., Chapront, J., Chapront-Touze, M.,
 J. 1994, Astron.Astrophys. 282, 663-683
 B., Kinoshita, H., Folgueira, M. 1999,
p.Ser. 135, 111
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Though t is strictly TDB, it is usually more convenient to use
     TT, which makes no significant difference.

  2) The expression used is as adopted in IERS Conventions (2003) and
     comes from Souchay et al. (1999) after Simon et al. (1994).


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     Astron.Astrophys.Supp.Ser. 135, 111

------
fal03
------

Given:
 Julian centuries since J2000.0 (Note 1)
e):
adians (Note 2)

 TDB, it is usually more convenient to use
ignificant difference.
is as adopted in IERS Conventions (2003) and
 (1994).

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
non, P., Chapront, J., Chapront-Touze, M.,
 J. 1994, Astron.Astrophys. 282, 663-683
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Though t is strictly TDB, it is usually more convenient to use
     TT, which makes no significant difference.

  2) The expression used is as adopted in IERS Conventions (2003) and
     is from Simon et al. (1994).


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

------
falp03
------

Given:
 Julian centuries since J2000.0 (Note 1)
e):
radians (Note 2)

 TDB, it is usually more convenient to use
ignificant difference.
is as adopted in IERS Conventions (2003) and
 (1994).

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
non, P., Chapront, J., Chapront-Touze, M.,
 J. 1994, Astron.Astrophys. 282, 663-683
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Though t is strictly TDB, it is usually more convenient to use
     TT, which makes no significant difference.

  2) The expression used is as adopted in IERS Conventions (2003) and
     is from Simon et al. (1994).


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

------
fama03
------

Given:
 Julian centuries since J2000.0 (Note 1)
e):
 longitude of Mars, radians (Note 2)

 TDB, it is usually more convenient to use
ignificant difference.
is as adopted in IERS Conventions (2003) and
t al. (1999) after Simon et al. (1994).

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
non, P., Chapront, J., Chapront-Touze, M.,
 J. 1994, Astron.Astrophys. 282, 663-683
 B., Kinoshita, H., Folgueira, M. 1999,
p.Ser. 135, 111
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Though t is strictly TDB, it is usually more convenient to use
     TT, which makes no significant difference.

  2) The expression used is as adopted in IERS Conventions (2003) and
     comes from Souchay et al. (1999) after Simon et al. (1994).


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     Astron.Astrophys.Supp.Ser. 135, 111

------
fame03
------

Given:
 Julian centuries since J2000.0 (Note 1)
e):
 longitude of Mercury, radians (Note 2)

 TDB, it is usually more convenient to use
ignificant difference.
is as adopted in IERS Conventions (2003) and
t al. (1999) after Simon et al. (1994).

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
non, P., Chapront, J., Chapront-Touze, M.,
 J. 1994, Astron.Astrophys. 282, 663-683
 B., Kinoshita, H., Folgueira, M. 1999,
p.Ser. 135, 111
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Though t is strictly TDB, it is usually more convenient to use
     TT, which makes no significant difference.

  2) The expression used is as adopted in IERS Conventions (2003) and
     comes from Souchay et al. (1999) after Simon et al. (1994).


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     Astron.Astrophys.Supp.Ser. 135, 111

------
fane03
------

Given:
 Julian centuries since J2000.0 (Note 1)
e):
 longitude of Neptune, radians (Note 2)

 TDB, it is usually more convenient to use
ignificant difference.
is as adopted in IERS Conventions (2003) and
n et al. (1994).

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
non, P., Chapront, J., Chapront-Touze, M.,
 J. 1994, Astron.Astrophys. 282, 663-683
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Though t is strictly TDB, it is usually more convenient to use
     TT, which makes no significant difference.

  2) The expression used is as adopted in IERS Conventions (2003) and
     is adapted from Simon et al. (1994).


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

------
faom03
------

Given:
 Julian centuries since J2000.0 (Note 1)
e):
a, radians (Note 2)

 TDB, it is usually more convenient to use
ignificant difference.
is as adopted in IERS Conventions (2003) and
 (1994).

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
non, P., Chapront, J., Chapront-Touze, M.,
 J. 1994, Astron.Astrophys. 282, 663-683
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Though t is strictly TDB, it is usually more convenient to use
     TT, which makes no significant difference.

  2) The expression used is as adopted in IERS Conventions (2003) and
     is from Simon et al. (1994).


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

------
fapa03
------

Given:
 Julian centuries since J2000.0 (Note 1)
e):
ral precession in longitude, radians (Note 2)

 TDB, it is usually more convenient to use
ignificant difference.
is as adopted in IERS Conventions (2003).  It
ita & Souchay (1990) and comes originally
1977).

uchay J. 1990, Celest.Mech. and Dyn.Astron.

e, T., Fricke, W. & Morando, B. 1977,
, 1-16
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Though t is strictly TDB, it is usually more convenient to use
     TT, which makes no significant difference.

  2) The expression used is as adopted in IERS Conventions (2003).  It
     is taken from Kinoshita & Souchay (1990) and comes originally
     from Lieske et al. (1977).


References:

     Kinoshita, H. and Souchay J. 1990, Celest.Mech. and Dyn.Astron.
     48, 187

     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
     Astron.Astrophys. 58, 1-16

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

------
fasa03
------

Given:
 Julian centuries since J2000.0 (Note 1)
e):
 longitude of Saturn, radians (Note 2)

 TDB, it is usually more convenient to use
ignificant difference.
is as adopted in IERS Conventions (2003) and
t al. (1999) after Simon et al. (1994).

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
non, P., Chapront, J., Chapront-Touze, M.,
 J. 1994, Astron.Astrophys. 282, 663-683
 B., Kinoshita, H., Folgueira, M. 1999,
p.Ser. 135, 111
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Though t is strictly TDB, it is usually more convenient to use
     TT, which makes no significant difference.

  2) The expression used is as adopted in IERS Conventions (2003) and
     comes from Souchay et al. (1999) after Simon et al. (1994).


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     Astron.Astrophys.Supp.Ser. 135, 111

------
faur03
------

Given:
 Julian centuries since J2000.0 (Note 1)
ue):
 longitude of Uranus, radians (Note 2)

 TDB, it is usually more convenient to use
ignificant difference.
is as adopted in IERS Conventions (2003) and
n et al. (1994).

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
non, P., Chapront, J., Chapront-Touze, M.,
 J. 1994, Astron.Astrophys. 282, 663-683
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Though t is strictly TDB, it is usually more convenient to use
     TT, which makes no significant difference.

  2) The expression used is as adopted in IERS Conventions (2003) and
     is adapted from Simon et al. (1994).


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

------
fave03
------

Given:
 Julian centuries since J2000.0 (Note 1)
e):
 longitude of Venus, radians (Note 2)

 TDB, it is usually more convenient to use
ignificant difference.
is as adopted in IERS Conventions (2003) and
t al. (1999) after Simon et al. (1994).

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
non, P., Chapront, J., Chapront-Touze, M.,
 J. 1994, Astron.Astrophys. 282, 663-683
 B., Kinoshita, H., Folgueira, M. 1999,
p.Ser. 135, 111
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Though t is strictly TDB, it is usually more convenient to use
     TT, which makes no significant difference.

  2) The expression used is as adopted in IERS Conventions (2003) and
     comes from Souchay et al. (1999) after Simon et al. (1994).


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     Astron.Astrophys.Supp.Ser. 135, 111

------
fk52h
------


Notes:

  1) This function transforms FK5 star positions and proper motions
     into the system of the Hipparcos catalog.

  2) The proper motions in RA are dRA/dt rather than
     cos(Dec)*dRA/dt, and are per year rather than per century.

  3) The FK5 to Hipparcos transformation is modeled as a pure
     rotation and spin;  zonal errors in the FK5 catalog are not
     taken into account.

  4) See also eraH2fk5, eraFk5hz, eraHfk5z.


------
fk5hip
------


Notes:

  1) This function models the FK5 to Hipparcos transformation as a
     pure rotation and spin;  zonal errors in the FK5 catalogue are
     not taken into account.

  2) The r-matrix r5h operates in the sense:

           P_Hipparcos = r5h x P_FK5

     where P_FK5 is a p-vector in the FK5 frame, and P_Hipparcos is
     the equivalent Hipparcos p-vector.

  3) The r-vector s5h represents the time derivative of the FK5 to
     Hipparcos rotation.  The units are radians per year (Julian,
     TDB).


------
fk5hz
------

Given:
  FK5 RA (radians), equinox J2000.0, at date
  FK5 Dec (radians), equinox J2000.0, at date
  TDB date (Notes 1,2)

  Hipparcos RA (radians)
  Hipparcos Dec (radians)

ts a star position from the FK5 system to
, in such a way that the Hipparcos proper
ause such a star has, in general, a non-zero
 FK5 system, the function requires the date
n in the FK5 system was determined.
te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
 transformation is modeled as a pure
zonal errors in the FK5 catalogue are not

d by this function is in the Hipparcos
 at date date1+date2.
raH2fk5, eraHfk5z.

al coordinates to unit vector
Hipparcos rotation and spin
y p-vector by scalar
r to r-matrix
 of transpose of r-matrix and p-vector
product of two p-vectors
r to spherical
ze angle into range 0 to 2pi

hle, 2000, Astron.Astrophys. 354, 732-739.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) This function converts a star position from the FK5 system to
     the Hipparcos system, in such a way that the Hipparcos proper
     motion is zero.  Because such a star has, in general, a non-zero
     proper motion in the FK5 system, the function requires the date
     at which the position in the FK5 system was determined.

  2) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  3) The FK5 to Hipparcos transformation is modeled as a pure
     rotation and spin;  zonal errors in the FK5 catalogue are not
     taken into account.

  4) The position returned by this function is in the Hipparcos
     reference system but at date date1+date2.

  5) See also eraFk52h, eraH2fk5, eraHfk5z.


------
fw2m
------

Given:
    F-W angle gamma_bar (radians)
    F-W angle phi_bar (radians)
    F-W angle psi (radians)
    F-W angle epsilon (radians)

]   rotation matrix

 points:
liptic pole,

ole of date,

illiams angles are as follows:
E


P
ing the combined effects of frame bias,
ion is:
.R_3(-psi).R_1(phib).R_3(gamb)
ices can be constructed, depending on the

ation x precession x frame bias matrix,
 precession angles, generate the nutation
d them to the psi_bar and epsilon_A angles,
ent function.
cession x frame bias matrix, generate the
ngles and call the present function.
me bias matrix, generate the four precession
2000.0 and call the present function.
d precession-only matrices can if necessary
ning these three appropriately.

ize r-matrix to identity
around Z-axis
around X-axis

006, Celest.Mech.Dyn.Astron. 94, 351
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Naming the following points:

           e = J2000.0 ecliptic pole,
           p = GCRS pole,
           E = ecliptic pole of date,
     and   P = CIP,

     the four Fukushima-Williams angles are as follows:

        gamb = gamma = epE
        phib = phi = pE
        psi = psi = pEP
        eps = epsilon = EP

  2) The matrix representing the combined effects of frame bias,
     precession and nutation is:

        NxPxB = R_1(-eps).R_3(-psi).R_1(phib).R_3(gamb)

  3) Three different matrices can be constructed, depending on the
     supplied angles:

     o  To obtain the nutation x precession x frame bias matrix,
        generate the four precession angles, generate the nutation
        components and add them to the psi_bar and epsilon_A angles,
        and call the present function.

     o  To obtain the precession x frame bias matrix, generate the
        four precession angles and call the present function.

     o  To obtain the frame bias matrix, generate the four precession
        angles for date J2000.0 and call the present function.

     The nutation-only and precession-only matrices can if necessary
     be obtained by combining these three appropriately.


------
fw2xy
------

Given:
-W angle gamma_bar (radians)
-W angle phi_bar (radians)
-W angle psi (radians)
-W angle epsilon (radians)

IP unit vector X,Y

 points:
liptic pole,

ole of date,

illiams angles are as follows:
E


P
ing the combined effects of frame bias,
ion is:
).R_3(-psi).R_1(phib).R_3(gamb)
x,y are elements [2][0] and [2][1] of the
0, they are essentially angles in radians.

les to r-matrix
 CIP X,Y coordinates from NPB matrix

006, Celest.Mech.Dyn.Astron. 94, 351
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) Naming the following points:

           e = J2000.0 ecliptic pole,
           p = GCRS pole
           E = ecliptic pole of date,
     and   P = CIP,

     the four Fukushima-Williams angles are as follows:

        gamb = gamma = epE
        phib = phi = pE
        psi = psi = pEP
        eps = epsilon = EP

  2) The matrix representing the combined effects of frame bias,
     precession and nutation is:

        NxPxB = R_1(-epsA).R_3(-psi).R_1(phib).R_3(gamb)

     The returned values x,y are elements [2][0] and [2][1] of the
     matrix.  Near J2000.0, they are essentially angles in radians.


------
g2icrs
------

Given:
alactic longitude (radians)
alactic latitude (radians)

CRS right ascension (radians)
CRS declination (radians)

of Galactic coordinates was defined with
bsolete reference system FK4 B1950.0.  When
tem in a modern context, several factors have
ount:
K4 positions of the E-terms of aberration.
the FK4 proper motion system by differential

50.0 equinox rather than the now-standard

ween ICRS and the J2000.0 mean place system.
gue (Perryman & ESA 1997) provides a rotation
ms directly between ICRS and Galactic
 above factors taken into account.  The
om three angles, namely the ICRS coordinates
 and the longitude of the ascending node of
 on the ICRS equator.  They are given in
mal places and for canonical purposes are
In the Hipparcos Catalogue the matrix
o 10 decimal places (about 20 microarcsec).
function the matrix elements have been
canonical three angles and are given to 30

mation is performed by the function eraIcrs2g.

ze angle into range 0 to 2pi
ze angle into range +/- pi
al coordinates to unit vector
 of transpose of r-matrix and p-vector
r to spherical

A, 1997, ESA SP-1200, The Hipparcos and Tycho
tric and photometric star catalogues
 Hipparcos Space Astrometry Mission.  ESA
n, Noordwijk, Netherlands.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The IAU 1958 system of Galactic coordinates was defined with
     respect to the now obsolete reference system FK4 B1950.0.  When
     interpreting the system in a modern context, several factors have
     to be taken into account:

     . The inclusion in FK4 positions of the E-terms of aberration.

     . The distortion of the FK4 proper motion system by differential
       Galactic rotation.

     . The use of the B1950.0 equinox rather than the now-standard
       J2000.0.

     . The frame bias between ICRS and the J2000.0 mean place system.

     The Hipparcos Catalogue (Perryman & ESA 1997) provides a rotation
     matrix that transforms directly between ICRS and Galactic
     coordinates with the above factors taken into account.  The
     matrix is derived from three angles, namely the ICRS coordinates
     of the Galactic pole and the longitude of the ascending node of
     the galactic equator on the ICRS equator.  They are given in
     degrees to five decimal places and for canonical purposes are
     regarded as exact.  In the Hipparcos Catalogue the matrix
     elements are given to 10 decimal places (about 20 microarcsec).
     In the present ERFA function the matrix elements have been
     recomputed from the canonical three angles and are given to 30
     decimal places.

  2) The inverse transformation is performed by the function eraIcrs2g.


------
gc2gd
------

Given:
llipsoid identifier (Note 1)
eocentric vector (Note 2)

ongitude (radians, east +ve, Note 3)
atitude (geodetic, radians, Note 3)
eight above ellipsoid (geodetic, Notes 2,3)
e):
tatus:  0 = OK
       -1 = illegal identifier (Note 3)
       -2 = internal error (Note 3)

a number that specifies the choice of
  The following are supported:




ignificance outside the ERFA software.  For
 ERFA_WGS84 etc. are defined in erfam.h.
r (xyz, given) and height (height, returned)

eans that the identifier n is illegal.  An
heoretically impossible.  In all error cases,
e set to -1e9.
mation is performed in the function eraGd2gc.

eference ellipsoids
ric to geodetic transformation, general
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The identifier n is a number that specifies the choice of
     reference ellipsoid.  The following are supported:

        n    ellipsoid

        1     ERFA_WGS84
        2     ERFA_GRS80
        3     ERFA_WGS72

     The n value has no significance outside the ERFA software.  For
     convenience, symbols ERFA_WGS84 etc. are defined in erfam.h.

  2) The geocentric vector (xyz, given) and height (height, returned)
     are in meters.

  3) An error status -1 means that the identifier n is illegal.  An
     error status -2 is theoretically impossible.  In all error cases,
     all three results are set to -1e9.

  4) The inverse transformation is performed in the function eraGd2gc.


------
gc2gde
------

Given:
quatorial radius (Notes 2,4)
lattening (Note 3)
eocentric vector (Note 4)

ongitude (radians, east +ve)
atitude (geodetic, radians)
eight above ellipsoid (geodetic, Note 4)
e):
tatus:  0 = OK
       -1 = illegal f
       -2 = illegal a

ed on the GCONV2H Fortran subroutine by
e reference).
s, a, can be in any units, but meters is
ice.
s (for the Earth) a value around 0.00335,

s, a, and the geocentric vector, xyz,
 same units, and determine the units of
 height.
status < 0), elong, phi and height are

mation is performed in the function

or a standard ellipsoid (such as ERFA_WGS84) can
 performed by calling eraGc2gd, which uses a
entify the required A and F values.

sformation from Cartesian to geodetic
ted by Halley's method", J.Geodesy (2006)

, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
gd2gc
------

Given:
llipsoid identifier (Note 1)
ongitude (radians, east +ve)
atitude (geodetic, radians, Note 3)
eight above ellipsoid (geodetic, Notes 2,3)

eocentric vector (Note 2)
e):
tatus:  0 = OK
       -1 = illegal identifier (Note 3)
       -2 = illegal case (Note 3)

a number that specifies the choice of
  The following are supported:




ignificance outside the ERFA software.  For
 ERFA_WGS84 etc. are defined in erfam.h.
given) and the geocentric vector (xyz,
ers.
formed on the arguments elong, phi and
atus -1 means that the identifier n is
tatus -2 protects against cases that would
xceptions.  In all error cases, xyz is set

mation is performed in the function eraGc2gd.

eference ellipsoids
c to geocentric transformation, general
vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The identifier n is a number that specifies the choice of
     reference ellipsoid.  The following are supported:

        n    ellipsoid

        1     ERFA_WGS84
        2     ERFA_GRS80
        3     ERFA_WGS72

     The n value has no significance outside the ERFA software.  For
     convenience, symbols ERFA_WGS84 etc. are defined in erfam.h.

  2) The height (height, given) and the geocentric vector (xyz,
     returned) are in meters.

  3) No validation is performed on the arguments elong, phi and
     height.  An error status -1 means that the identifier n is
     illegal.  An error status -2 protects against cases that would
     lead to arithmetic exceptions.  In all error cases, xyz is set
     to zeros.

  4) The inverse transformation is performed in the function eraGc2gd.


------
gd2gce
------

Given:
quatorial radius (Notes 1,4)
lattening (Notes 2,4)
ongitude (radians, east +ve)
atitude (geodetic, radians, Note 4)
eight above ellipsoid (geodetic, Notes 3,4)

eocentric vector (Note 3)
e):
tatus:  0 = OK
       -1 = illegal case (Note 4)

s, a, can be in any units, but meters is
ice.
s (for the Earth) a value around 0.00335,

s, a, and the height, height, must be
its, and determine the units of the
vector, xyz.
formed on individual arguments.  The error
gainst (unrealistic) cases that would lead
ions.  If an error occurs, xyz is unchanged.
mation is performed in the function

or a standard ellipsoid (such as ERFA_WGS84) can
 performed by calling eraGd2gc,  which uses a
entify the required a and f values.

al Astronomy, Cambridge University Press,
p96.
nt to the Astronomical Almanac,
n (ed), University Science Books (1992),

, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The equatorial radius, a, can be in any units, but meters is
     the conventional choice.

  2) The flattening, f, is (for the Earth) a value around 0.00335,
     i.e. around 1/298.

  3) The equatorial radius, a, and the height, height, must be
     given in the same units, and determine the units of the
     returned geocentric vector, xyz.

  4) No validation is performed on individual arguments.  The error
     status -1 protects against (unrealistic) cases that would lead
     to arithmetic exceptions.  If an error occurs, xyz is unchanged.

  5) The inverse transformation is performed in the function
     eraGc2gde.

  6) The transformation for a standard ellipsoid (such as ERFA_WGS84) can
     more conveniently be performed by calling eraGd2gc,  which uses a
     numerical code to identify the required a and f values.


References:

     Green, R.M., Spherical Astronomy, Cambridge University Press,
     (1985) Section 4.5, p96.

     Explanatory Supplement to the Astronomical Almanac,
     P. Kenneth Seidelmann (ed), University Science Books (1992),
     Section 4.22, p202.

------
gmst00
------

Given:
 UT1 as a 2-part Julian Date (Notes 1,2)
 TT as a 2-part Julian Date (Notes 1,2)
e):
 Greenwich mean sidereal time (radians)

 uta+utb and tta+ttb respectively, are both
ioned in any convenient way between the
 example, JD=2450123.7 could be expressed in
mong others:
  Part B
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
e case of UT;  the TT is not at all critical
he J2000 and MJD methods are good compromises
nd convenience.  For UT, the date & time
ed to the algorithm that is used by the Earth
ion, called internally:  maximum precision is
ta argument is for 0hrs UT1 on the day in
 argument lies in the range 0 to 1, or vice

required, UT1 to predict the Earth rotation
e effects of precession.  If UT1 is used for
s of order 100 microarcseconds result.
ble with the IAU 2000 resolutions and must be
tion with other IAU 2000 compatible
recession-nutation and equation of the

ed in the range 0 to 2pi.
m Capitaine et al. (2003) and IERS


otation angle, IAU 2000
ze angle into range 0 to 2pi

ce, P.T. and McCarthy, D.D., "Expressions to
00 definition of UT1", Astronomy &
135-1149 (2003)
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
     Julian Dates, apportioned in any convenient way between the
     argument pairs.  For example, JD=2450123.7 could be expressed in
     any of these ways, among others:

            Part A         Part B

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable (in the case of UT;  the TT is not at all critical
     in this respect).  The J2000 and MJD methods are good compromises
     between resolution and convenience.  For UT, the date & time
     method is best matched to the algorithm that is used by the Earth
     Rotation Angle function, called internally:  maximum precision is
     delivered when the uta argument is for 0hrs UT1 on the day in
     question and the utb argument lies in the range 0 to 1, or vice
     versa.

  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
     and TT to predict the effects of precession.  If UT1 is used for
     both purposes, errors of order 100 microarcseconds result.

  3) This GMST is compatible with the IAU 2000 resolutions and must be
     used only in conjunction with other IAU 2000 compatible
     components such as precession-nutation and equation of the
     equinoxes.

  4) The result is returned in the range 0 to 2pi.

  5) The algorithm is from Capitaine et al. (2003) and IERS
     Conventions 2003.


References:

     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
     implement the IAU 2000 definition of UT1", Astronomy &
     Astrophysics, 406, 1135-1149 (2003)

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

------
gmst06
------

Given:
 UT1 as a 2-part Julian Date (Notes 1,2)
 TT as a 2-part Julian Date (Notes 1,2)
e):
 Greenwich mean sidereal time (radians)

 uta+utb and tta+ttb respectively, are both
ioned in any convenient way between the
 example, JD=2450123.7 could be expressed in
mong others:
 Part B
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
e case of UT;  the TT is not at all critical
he J2000 and MJD methods are good compromises
nd convenience.  For UT, the date & time
ed to the algorithm that is used by the Earth
ion, called internally:  maximum precision is
ta argument is for 0hrs UT1 on the day in
 argument lies in the range 0 to 1, or vice

required, UT1 to predict the Earth rotation
e effects of precession.  If UT1 is used for
s of order 100 microarcseconds result.
ble with the IAU 2006 precession and must not
recession models.
ed in the range 0 to 2pi.

otation angle, IAU 2000
ze angle into range 0 to 2pi

ce, P.T. & Chapront, J., 2005,
2, 355
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
     Julian Dates, apportioned in any convenient way between the
     argument pairs.  For example, JD=2450123.7 could be expressed in
     any of these ways, among others:

            Part A        Part B

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable (in the case of UT;  the TT is not at all critical
     in this respect).  The J2000 and MJD methods are good compromises
     between resolution and convenience.  For UT, the date & time
     method is best matched to the algorithm that is used by the Earth
     rotation angle function, called internally:  maximum precision is
     delivered when the uta argument is for 0hrs UT1 on the day in
     question and the utb argument lies in the range 0 to 1, or vice
     versa.

  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
     and TT to predict the effects of precession.  If UT1 is used for
     both purposes, errors of order 100 microarcseconds result.

  3) This GMST is compatible with the IAU 2006 precession and must not
     be used with other precession models.

  4) The result is returned in the range 0 to 2pi.


------
gmst82
------

Given:
 UT1 Julian Date (see note)
e):
 Greenwich mean sidereal time (radians)

 is a Julian Date, apportioned in any
en the arguments dj1 and dj2.  For example,
uld be expressed in any of these ways,

   dj2
   0          (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
  0.2         (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 and MJD methods are good compromises
nd convenience.  The date & time method is
algorithm used:  maximum accuracy (or, at
) is delivered when the dj1 argument is for
in question and the dj2 argument lies in the
e versa.
ed on the IAU 1982 expression.  This is
giving the GMST at 0 hours UT1.  In fact, it
 between the GMST and the UT, the steady
awing-ahead of ST with respect to UT.  When
ed, the expression happens to equal the GMST
day.
e entire UT1 (the sum of the two arguments
 directly as the argument for the standard
t term of which is adjusted by 12 hours to
noon phasing of Julian Date.  The UT1 is then
whole days to conserve accuracy.

ze angle into range 0 to 2pi

International Astronomical Union,

 Astrophys. 105, 359-361 (1982).
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The UT1 date dj1+dj2 is a Julian Date, apportioned in any
     convenient way between the arguments dj1 and dj2.  For example,
     JD(UT1)=2450123.7 could be expressed in any of these ways,
     among others:

             dj1            dj2

         2450123.7          0          (JD method)
          2451545        -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5         0.2         (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 and MJD methods are good compromises
     between resolution and convenience.  The date & time method is
     best matched to the algorithm used:  maximum accuracy (or, at
     least, minimum noise) is delivered when the dj1 argument is for
     0hrs UT1 on the day in question and the dj2 argument lies in the
     range 0 to 1, or vice versa.

  2) The algorithm is based on the IAU 1982 expression.  This is
     always described as giving the GMST at 0 hours UT1.  In fact, it
     gives the difference between the GMST and the UT, the steady
     4-minutes-per-day drawing-ahead of ST with respect to UT.  When
     whole days are ignored, the expression happens to equal the GMST
     at 0 hours UT1 each day.

  3) In this function, the entire UT1 (the sum of the two arguments
     dj1 and dj2) is used directly as the argument for the standard
     formula, the constant term of which is adjusted by 12 hours to
     take account of the noon phasing of Julian Date.  The UT1 is then
     added, but omitting whole days to conserve accuracy.


References:

     Transactions of the International Astronomical Union,
     XVIII B, 67 (1983).

     Aoki et al., Astron. Astrophys. 105, 359-361 (1982).

------
gst00a
------

Given:
 UT1 as a 2-part Julian Date (Notes 1,2)
 TT as a 2-part Julian Date (Notes 1,2)
e):
 Greenwich apparent sidereal time (radians)

 uta+utb and tta+ttb respectively, are both
ioned in any convenient way between the
 example, JD=2450123.7 could be expressed in
mong others:
 Part B
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
e case of UT;  the TT is not at all critical
he J2000 and MJD methods are good compromises
nd convenience.  For UT, the date & time
ed to the algorithm that is used by the Earth
ion, called internally:  maximum precision is
ta argument is for 0hrs UT1 on the day in
 argument lies in the range 0 to 1, or vice

required, UT1 to predict the Earth rotation
e effects of precession-nutation.  If UT1 is
es, errors of order 100 microarcseconds

ble with the IAU 2000 resolutions and must be
tion with other IAU 2000 compatible
recession-nutation.
ed in the range 0 to 2pi.
m Capitaine et al. (2003) and IERS


ch mean sidereal time, IAU 2000
n of the equinoxes, IAU 2000A
ze angle into range 0 to 2pi

ce, P.T. and McCarthy, D.D., "Expressions to
00 definition of UT1", Astronomy &
135-1149 (2003)
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
     Julian Dates, apportioned in any convenient way between the
     argument pairs.  For example, JD=2450123.7 could be expressed in
     any of these ways, among others:

            Part A        Part B

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable (in the case of UT;  the TT is not at all critical
     in this respect).  The J2000 and MJD methods are good compromises
     between resolution and convenience.  For UT, the date & time
     method is best matched to the algorithm that is used by the Earth
     Rotation Angle function, called internally:  maximum precision is
     delivered when the uta argument is for 0hrs UT1 on the day in
     question and the utb argument lies in the range 0 to 1, or vice
     versa.

  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
     and TT to predict the effects of precession-nutation.  If UT1 is
     used for both purposes, errors of order 100 microarcseconds
     result.

  3) This GAST is compatible with the IAU 2000 resolutions and must be
     used only in conjunction with other IAU 2000 compatible
     components such as precession-nutation.

  4) The result is returned in the range 0 to 2pi.

  5) The algorithm is from Capitaine et al. (2003) and IERS
     Conventions 2003.


References:

     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
     implement the IAU 2000 definition of UT1", Astronomy &
     Astrophysics, 406, 1135-1149 (2003)

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

------
gst00b
------

Given:
 UT1 as a 2-part Julian Date (Notes 1,2)
e):
 Greenwich apparent sidereal time (radians)

 is a Julian Date, apportioned in any
en the argument pair.  For example,
e expressed in any of these ways, among

   utb
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 and MJD methods are good compromises
nd convenience.  For UT, the date & time
ed to the algorithm that is used by the Earth
ion, called internally:  maximum precision is
ta argument is for 0hrs UT1 on the day in
 argument lies in the range 0 to 1, or vice

ible with the IAU 2000 resolutions, except
en compromised for the sake of speed and
espects:
 of TDB (or TT) to compute the precession
and the equation of the equinoxes.  This
of order 0.1 mas at present.
dged nutation model (McCarthy & Luzum, 2001)
ng errors of up to 1 mas.
ble with the IAU 2000 resolutions and must be
tion with other IAU 2000 compatible
recession-nutation.
ed in the range 0 to 2pi.
m Capitaine et al. (2003) and IERS


ch mean sidereal time, IAU 2000
n of the equinoxes, IAU 2000B
ze angle into range 0 to 2pi

ce, P.T. and McCarthy, D.D., "Expressions to
00 definition of UT1", Astronomy &
135-1149 (2003)
um, B.J., "An abridged model of the
of the celestial pole", Celestial Mechanics &
 85, 37-49 (2003)
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The UT1 date uta+utb is a Julian Date, apportioned in any
     convenient way between the argument pair.  For example,
     JD=2450123.7 could be expressed in any of these ways, among
     others:

             uta            utb

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in cases
     where the loss of several decimal digits of resolution is
     acceptable.  The J2000 and MJD methods are good compromises
     between resolution and convenience.  For UT, the date & time
     method is best matched to the algorithm that is used by the Earth
     Rotation Angle function, called internally:  maximum precision is
     delivered when the uta argument is for 0hrs UT1 on the day in
     question and the utb argument lies in the range 0 to 1, or vice
     versa.

  2) The result is compatible with the IAU 2000 resolutions, except
     that accuracy has been compromised for the sake of speed and
     convenience in two respects:

     . UT is used instead of TDB (or TT) to compute the precession
       component of GMST and the equation of the equinoxes.  This
       results in errors of order 0.1 mas at present.

     . The IAU 2000B abridged nutation model (McCarthy & Luzum, 2001)
       is used, introducing errors of up to 1 mas.

  3) This GAST is compatible with the IAU 2000 resolutions and must be
     used only in conjunction with other IAU 2000 compatible
     components such as precession-nutation.

  4) The result is returned in the range 0 to 2pi.

  5) The algorithm is from Capitaine et al. (2003) and IERS
     Conventions 2003.


References:

     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
     implement the IAU 2000 definition of UT1", Astronomy &
     Astrophysics, 406, 1135-1149 (2003)

     McCarthy, D.D. & Luzum, B.J., "An abridged model of the
     precession-nutation of the celestial pole", Celestial Mechanics &
     Dynamical Astronomy, 85, 37-49 (2003)

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

------
gst06
------

Given:
   UT1 as a 2-part Julian Date (Notes 1,2)
   TT as a 2-part Julian Date (Notes 1,2)
]  nutation x precession x bias matrix
e):
   Greenwich apparent sidereal time (radians)

 uta+utb and tta+ttb respectively, are both
ioned in any convenient way between the
 example, JD=2450123.7 could be expressed in
mong others:
 Part B
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
e case of UT;  the TT is not at all critical
he J2000 and MJD methods are good compromises
nd convenience.  For UT, the date & time
ed to the algorithm that is used by the Earth
ion, called internally:  maximum precision is
ta argument is for 0hrs UT1 on the day in
 argument lies in the range 0 to 1, or vice

required, UT1 to predict the Earth rotation
e effects of precession-nutation.  If UT1 is
es, errors of order 100 microarcseconds

n uses the IAU 2006 series for s+XY/2, it is
t of the precession-nutation model and can in
h any equinox-based NPB matrix.
ed in the range 0 to 2pi.

 CIP X,Y coordinates from NPB matrix
 locator s, given X,Y, IAU 2006
ze angle into range 0 to 2pi
otation angle, IAU 2000
n of the origins, given NPB matrix and s

taine, N., 2006, Astron.Astrophys. 459, 981
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
     Julian Dates, apportioned in any convenient way between the
     argument pairs.  For example, JD=2450123.7 could be expressed in
     any of these ways, among others:

            Part A        Part B

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable (in the case of UT;  the TT is not at all critical
     in this respect).  The J2000 and MJD methods are good compromises
     between resolution and convenience.  For UT, the date & time
     method is best matched to the algorithm that is used by the Earth
     rotation angle function, called internally:  maximum precision is
     delivered when the uta argument is for 0hrs UT1 on the day in
     question and the utb argument lies in the range 0 to 1, or vice
     versa.

  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
     and TT to predict the effects of precession-nutation.  If UT1 is
     used for both purposes, errors of order 100 microarcseconds
     result.

  3) Although the function uses the IAU 2006 series for s+XY/2, it is
     otherwise independent of the precession-nutation model and can in
     practice be used with any equinox-based NPB matrix.

  4) The result is returned in the range 0 to 2pi.


------
gst06a
------

Given:
 UT1 as a 2-part Julian Date (Notes 1,2)
 TT as a 2-part Julian Date (Notes 1,2)
e):
 Greenwich apparent sidereal time (radians)

 uta+utb and tta+ttb respectively, are both
ioned in any convenient way between the
 example, JD=2450123.7 could be expressed in
mong others:
 Part B
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
e case of UT;  the TT is not at all critical
he J2000 and MJD methods are good compromises
nd convenience.  For UT, the date & time
ed to the algorithm that is used by the Earth
ion, called internally:  maximum precision is
ta argument is for 0hrs UT1 on the day in
 argument lies in the range 0 to 1, or vice

required, UT1 to predict the Earth rotation
e effects of precession-nutation.  If UT1 is
es, errors of order 100 microarcseconds

ble with the IAU 2000/2006 resolutions and
 conjunction with IAU 2006 precession and

ed in the range 0 to 2pi.

al NPB matrix, IAU 2006/2000A
ch apparent ST, IAU 2006, given NPB matrix

taine, N., 2006, Astron.Astrophys. 459, 981
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
     Julian Dates, apportioned in any convenient way between the
     argument pairs.  For example, JD=2450123.7 could be expressed in
     any of these ways, among others:

            Part A        Part B

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable (in the case of UT;  the TT is not at all critical
     in this respect).  The J2000 and MJD methods are good compromises
     between resolution and convenience.  For UT, the date & time
     method is best matched to the algorithm that is used by the Earth
     rotation angle function, called internally:  maximum precision is
     delivered when the uta argument is for 0hrs UT1 on the day in
     question and the utb argument lies in the range 0 to 1, or vice
     versa.

  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
     and TT to predict the effects of precession-nutation.  If UT1 is
     used for both purposes, errors of order 100 microarcseconds
     result.

  3) This GAST is compatible with the IAU 2000/2006 resolutions and
     must be used only in conjunction with IAU 2006 precession and
     IAU 2000A nutation.

  4) The result is returned in the range 0 to 2pi.


------
gst94
------

Given:
 UT1 as a 2-part Julian Date (Notes 1,2)
e):
 Greenwich apparent sidereal time (radians)

 is a Julian Date, apportioned in any
en the argument pair.  For example,
e expressed in any of these ways, among

   utb
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 and MJD methods are good compromises
nd convenience.  For UT, the date & time
ed to the algorithm that is used by the Earth
ion, called internally:  maximum precision is
ta argument is for 0hrs UT1 on the day in
 argument lies in the range 0 to 1, or vice

ible with the IAU 1982 and 1994 resolutions,
 has been compromised for the sake of
UT is used instead of TDB (or TT) to compute
equinoxes.
ed only in conjunction with contemporaneous
s 1976 precession, 1980 obliquity and 1982
 compatible with the IAU 2000 resolutions.
ed in the range 0 to 2pi.

ch mean sidereal time, IAU 1982
n of the equinoxes, IAU 1994
ze angle into range 0 to 2pi

nt to the Astronomical Almanac,
n (ed), University Science Books (1992)
ecommendation 3 (1994)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The UT1 date uta+utb is a Julian Date, apportioned in any
     convenient way between the argument pair.  For example,
     JD=2450123.7 could be expressed in any of these ways, among
     others:

             uta            utb

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in cases
     where the loss of several decimal digits of resolution is
     acceptable.  The J2000 and MJD methods are good compromises
     between resolution and convenience.  For UT, the date & time
     method is best matched to the algorithm that is used by the Earth
     Rotation Angle function, called internally:  maximum precision is
     delivered when the uta argument is for 0hrs UT1 on the day in
     question and the utb argument lies in the range 0 to 1, or vice
     versa.

  2) The result is compatible with the IAU 1982 and 1994 resolutions,
     except that accuracy has been compromised for the sake of
     convenience in that UT is used instead of TDB (or TT) to compute
     the equation of the equinoxes.

  3) This GAST must be used only in conjunction with contemporaneous
     IAU standards such as 1976 precession, 1980 obliquity and 1982
     nutation.  It is not compatible with the IAU 2000 resolutions.

  4) The result is returned in the range 0 to 2pi.


References:

     Explanatory Supplement to the Astronomical Almanac,
     P. Kenneth Seidelmann (ed), University Science Books (1992)

     IAU Resolution C7, Recommendation 3 (1994)

------
h2fk5
------


Notes:

  1) This function transforms Hipparcos star positions and proper
     motions into FK5 J2000.0.

  2) The proper motions in RA are dRA/dt rather than
     cos(Dec)*dRA/dt, and are per year rather than per century.

  3) The FK5 to Hipparcos transformation is modeled as a pure
     rotation and spin;  zonal errors in the FK5 catalog are not
     taken into account.

  4) See also eraFk52h, eraFk5hz, eraHfk5z.


------
hfk5z
------

Given:
    Hipparcos RA (radians)
    Hipparcos Dec (radians)
    TDB date (Note 1)
nox J2000.0, date date1+date2):
    RA (radians)
    Dec (radians)
    FK5 RA proper motion (rad/year, Note 4)
    Dec proper motion (rad/year, Note 4)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
 RA is dRA/dt rather than cos(Dec)*dRA/dt.
 transformation is modeled as a pure rotation
ors in the FK5 catalogue are not taken into

 that Hipparcos should be a close
inertial frame, so that distant objects have
 such objects have (in general) non-zero
, and this function returns those fictitious

d by this function is in the FK5 J2000.0
 at date date1+date2.
raH2fk5, eraFk5zhz.

al coordinates to unit vector
Hipparcos rotation and spin
 of r-matrix and p-vector
y p-vector by scalar
 of two r-matrices
 of transpose of r-matrix and p-vector
product of two p-vectors
or to spherical
ze angle into range 0 to 2pi

hle, 2000, Astron.Astrophys. 354, 732-739.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

  3) The FK5 to Hipparcos transformation is modeled as a pure rotation
     and spin;  zonal errors in the FK5 catalogue are not taken into
     account.

  4) It was the intention that Hipparcos should be a close
     approximation to an inertial frame, so that distant objects have
     zero proper motion;  such objects have (in general) non-zero
     proper motion in FK5, and this function returns those fictitious
     proper motions.

  5) The position returned by this function is in the FK5 J2000.0
     reference system but at date date1+date2.

  6) See also eraFk52h, eraH2fk5, eraFk5zhz.


------
icrs2g
------

Given:
CRS right ascension (radians)
CRS declination (radians)

alactic longitude (radians)
alactic latitude (radians)

of Galactic coordinates was defined with
bsolete reference system FK4 B1950.0.  When
tem in a modern context, several factors have
ount:
K4 positions of the E-terms of aberration.
the FK4 proper motion system by differential

50.0 equinox rather than the now-standard

ween ICRS and the J2000.0 mean place system.
gue (Perryman & ESA 1997) provides a rotation
ms directly between ICRS and Galactic
 above factors taken into account.  The
om three angles, namely the ICRS coordinates
 and the longitude of the ascending node of
 on the ICRS equator.  They are given in
mal places and for canonical purposes are
In the Hipparcos Catalogue the matrix
o 10 decimal places (about 20 microarcsec).
function the matrix elements have been
canonical three angles and are given to 30

mation is performed by the function eraG2icrs.

ze angle into range 0 to 2pi
ze angle into range +/- pi
al coordinates to unit vector
 of r-matrix and p-vector
r to spherical

A, 1997, ESA SP-1200, The Hipparcos and Tycho
tric and photometric star catalogues
 Hipparcos Space Astrometry Mission.  ESA
n, Noordwijk, Netherlands.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The IAU 1958 system of Galactic coordinates was defined with
     respect to the now obsolete reference system FK4 B1950.0.  When
     interpreting the system in a modern context, several factors have
     to be taken into account:

     . The inclusion in FK4 positions of the E-terms of aberration.

     . The distortion of the FK4 proper motion system by differential
       Galactic rotation.

     . The use of the B1950.0 equinox rather than the now-standard
       J2000.0.

     . The frame bias between ICRS and the J2000.0 mean place system.

     The Hipparcos Catalogue (Perryman & ESA 1997) provides a rotation
     matrix that transforms directly between ICRS and Galactic
     coordinates with the above factors taken into account.  The
     matrix is derived from three angles, namely the ICRS coordinates
     of the Galactic pole and the longitude of the ascending node of
     the galactic equator on the ICRS equator.  They are given in
     degrees to five decimal places and for canonical purposes are
     regarded as exact.  In the Hipparcos Catalogue the matrix
     elements are given to 10 decimal places (about 20 microarcsec).
     In the present ERFA function the matrix elements have been
     recomputed from the canonical three angles and are given to 30
     decimal places.

  2) The inverse transformation is performed by the function eraG2icrs.


------
ir
------



------
jd2cal
------

Given:
ulian Date (Notes 1, 2)

ear
onth
ay
raction of day
e):
tatus:
  0 = OK
 -1 = unacceptable date (Note 3)

ate is -68569.5 (-4900 March 1).  The
ed is 1e9.
pportioned in any convenient way between
d dj2.  For example, JD=2450123.7 could
of these ways, among others:
   dj2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
nversion is from the "proleptic Gregorian
nt is taken of the date(s) of adoption of
ar, nor is the AD/BC numbering convention


nt to the Astronomical Almanac,
n (ed), University Science Books (1992),
.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
jdcalf
------

Given:
umber of decimal places of days in fraction
j1+dj2 = Julian Date (Note 1)

ear, month, day, fraction in Gregorian
alendar
e):
tatus:
 -1 = date out of range
  0 = OK
 +1 = NDP not 0-9 (interpreted as 0)

pportioned in any convenient way between
d dj2.  For example, JD=2450123.7 could
of these ways, among others:
   dj2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
nversion is from the "Proleptic Gregorian
nt is taken of the date(s) of adoption of
ar, nor is the AD/BC numbering convention

n eraJd2cal.
ess if internal overflows are to be
which use 16-bit integers.

regorian calendar

nt to the Astronomical Almanac,
n (ed), University Science Books (1992),
.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The Julian Date is apportioned in any convenient way between
     the arguments dj1 and dj2.  For example, JD=2450123.7 could
     be expressed in any of these ways, among others:

             dj1            dj2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

  2) In early eras the conversion is from the "Proleptic Gregorian
     Calendar";  no account is taken of the date(s) of adoption of
     the Gregorian Calendar, nor is the AD/BC numbering convention
     observed.

  3) Refer to the function eraJd2cal.

  4) NDP should be 4 or less if internal overflows are to be
     avoided on machines which use 16-bit integers.


------
ld
------

Given:
ss of the gravitating body (solar masses)
rection from observer to source (unit vector)
rection from body to source (unit vector)
rection from body to observer (unit vector)
stance from body to observer (au)
flection limiter (Note 4)

server to deflected source (unit vector)

ed on Expr. (70) in Klioner (2003) and
Explanatory Supplement (Urban & Seidelmann
rrangement to minimize the effects of machine

m can, as required, be adjusted in order to
ts as quadrupole field.
tion of the deflecting body should ideally
me of closest approach of the light ray to

er parameter dlim is phi^2/2, where phi is
on (in radians) between source and body at
plied.  As phi shrinks below the chosen
ction is artificially reduced, reaching zero

p1 is not normalized, but the consequential
magnitude is always negligible.
p1 can be the same array.
light deflection taking into account the
everal bodies, call the present function for
ion, in decreasing order of distance from the

dation is omitted.  The supplied vectors must
, and the deflection limiter non-zero and


nn, P. K. (eds), Explanatory Supplement to
anac, 3rd ed., University Science Books

"A practical relativistic model for micro-
 in space", Astr. J. 125, 1580-1597 (2003).

product of two p-vectors
product of two p-vectors
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The algorithm is based on Expr. (70) in Klioner (2003) and
     Expr. (7.63) in the Explanatory Supplement (Urban & Seidelmann
     2013), with some rearrangement to minimize the effects of machine
     precision.

  2) The mass parameter bm can, as required, be adjusted in order to
     allow for such effects as quadrupole field.

  3) The barycentric position of the deflecting body should ideally
     correspond to the time of closest approach of the light ray to
     the body.

  4) The deflection limiter parameter dlim is phi^2/2, where phi is
     the angular separation (in radians) between source and body at
     which limiting is applied.  As phi shrinks below the chosen
     threshold, the deflection is artificially reduced, reaching zero
     for phi = 0.

  5) The returned vector p1 is not normalized, but the consequential
     departure from unit magnitude is always negligible.

  6) The arguments p and p1 can be the same array.

  7) To accumulate total light deflection taking into account the
     contributions from several bodies, call the present function for
     each body in succession, in decreasing order of distance from the
     observer.

  8) For efficiency, validation is omitted.  The supplied vectors must
     be of unit magnitude, and the deflection limiter non-zero and
     positive.

  References:

     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
     the Astronomical Almanac, 3rd ed., University Science Books
     (2013).

     Klioner, Sergei A., "A practical relativistic model for micro-
     arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).


References:

     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
     the Astronomical Almanac, 3rd ed., University Science Books
     (2013).

     Klioner, Sergei A., "A practical relativistic model for micro-
     arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).

  Called:
     eraPdp       scalar product of two p-vectors
     eraPxp       vector product of two p-vectors

------
ldn
------

Given:
umber of bodies (note 1)
ata for each of the n bodies (Notes 1,2):
 mass of the body (solar masses, Note 3)
 deflection limiter (Note 4)
 barycentric PV of the body (au, au/day)
arycentric position of the observer (au)
bserver to star coord direction (unit vector)

 observer to deflected star (unit vector)
 n entries, one for each body to be
0, no gravitational light deflection will be
r the Sun.
nclude an entry for the Sun as well as for
body to be taken into account.  The entries
er in which the light passes the body.
b array for body i, the mass parameter
ired, be adjusted in order to allow for such
e field.
er parameter b[i].dl is phi^2/2, where phi is
on (in radians) between star and body at
plied.  As phi shrinks below the chosen
ction is artificially reduced, reaching zero
le values suitable for a terrestrial
ith masses, are as follows:
m        b[i].dl
         6e-6
5435     3e-9
8574     3e-10
starlight passes the body before reaching the
s placed back along its barycentric track by
that point to the observer.  For cases where
 the observer no such shift is applied.  If
t is preferred, the user has the option of
aLd function.  Similarly, eraLd can be used
source is nearby, not a star.
sn is not normalized, but the consequential
magnitude is always negligible.
 sn can be the same array.
dation is omitted.  The supplied masses must
, the position and velocity vectors must be
ction limiter greater than zero.

nn, P. K. (eds), Explanatory Supplement to
anac, 3rd ed., University Science Books
4.

vector
product of two p-vectors
r minus p-vector
r plus scaled p-vector
se p-vector into modulus and direction
eflection by a solar-system body
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
ldsun
------

Given:
rection from observer to star (unit vector)
rection from Sun to observer (unit vector)
stance from Sun to observer (au)

server to deflected star (unit vector)

ed to be sufficiently distant that its
 the Sun and the observer are essentially

strained when the angle between the star and
n is less than about 9 arcsec, falling to
tion. (The chosen threshold is within the
olar-system applications.)
p1 can be the same array.

eflection by a solar-system body
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The source is presumed to be sufficiently distant that its
     directions seen from the Sun and the observer are essentially
     the same.

  2) The deflection is restrained when the angle between the star and
     the center of the Sun is less than about 9 arcsec, falling to
     zero for zero separation. (The chosen threshold is within the
     solar limb for all solar-system applications.)

  3) The arguments p and p1 can be the same array.


------
lteceq
------

Given:
ulian epoch (TT)
cliptic longitude and latitude (radians)

CRS right ascension and declination (radians)
ade about whether the coordinates represent
 astrometric effects such as parallax or

s approximately that from ecliptic longitude
quinox and ecliptic of date) to mean J2000.0
declination, with only frame bias (always
 disturb this classical picture.
2011, 2012) 400 millennia precession model
2006 precession at J2000.0 and stays within
during the 20th and 21st centuries.  It is
cseconds throughout the historical period,
enths of a degree at the end of the
e span.

al coordinates to unit vector
 to ecliptic rotation matrix, long term
 of transpose of r-matrix and p-vector
ctor to spherical coordinates
ze angle into range 0 to 2pi
ze angle into range +/- pi

e, N. and Wallace, P., 2011, New precession
r long time intervals, Astron.Astrophys. 534,

e, N. and Wallace, P., 2012, New precession
r long time intervals (Corrigendum),
, C1
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

References:

    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
    expressions, valid for long time intervals, Astron.Astrophys. 534,
    A22

    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
    expressions, valid for long time intervals (Corrigendum),
    Astron.Astrophys. 541, C1

------
ltecm
------

Given:
   Julian epoch (TT)

   ICRS to ecliptic rotation matrix

 sense
S,
ctor with respect to ICRS right ascension
 and E_ep is the same vector with respect to
tic and equinox of epoch epj.
tor, merely a direction, typically of unit
ound to any particular spatial origin, such
 SSB.  No assumptions are made about whether
ght and embodies astrometric effects such as
on.  The transformation is approximately that
 right ascension and declination and ecliptic
de, with only frame bias (always less than
his classical picture.
2011, 2012) 400 millennia precession model
2006 precession at J2000.0 and stays within
during the 20th and 21st centuries.  It is
cseconds throughout the historical period,
enths of a degree at the end of the
e span.

 pole, long term
c pole, long term
product
ze vector

e, N. and Wallace, P., 2011, New precession
r long time intervals, Astron.Astrophys. 534,

e, N. and Wallace, P., 2012, New precession
r long time intervals (Corrigendum),
, C1
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The matrix is in the sense

        E_ep = rm x P_ICRS,

     where P_ICRS is a vector with respect to ICRS right ascension
     and declination axes and E_ep is the same vector with respect to
     the (inertial) ecliptic and equinox of epoch epj.

  2) P_ICRS is a free vector, merely a direction, typically of unit
     magnitude, and not bound to any particular spatial origin, such
     as the Earth, Sun or SSB.  No assumptions are made about whether
     it represents starlight and embodies astrometric effects such as
     parallax or aberration.  The transformation is approximately that
     between mean J2000.0 right ascension and declination and ecliptic
     longitude and latitude, with only frame bias (always less than
     25 mas) to disturb this classical picture.

  3) The Vondrak et al. (2011, 2012) 400 millennia precession model
     agrees with the IAU 2006 precession at J2000.0 and stays within
     100 microarcseconds during the 20th and 21st centuries.  It is
     accurate to a few arcseconds throughout the historical period,
     worsening to a few tenths of a degree at the end of the
     +/- 200,000 year time span.


References:

    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
    expressions, valid for long time intervals, Astron.Astrophys. 534,
    A22

    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
    expressions, valid for long time intervals (Corrigendum),
    Astron.Astrophys. 541, C1

------
lteqec
------

Given:
ulian epoch (TT)
CRS right ascension and declination (radians)

cliptic longitude and latitude (radians)
ade about whether the coordinates represent
 astrometric effects such as parallax or

s approximately that from mean J2000.0 right
ation to ecliptic longitude and latitude
liptic of date), with only frame bias (always
 disturb this classical picture.
2011, 2012) 400 millennia precession model
2006 precession at J2000.0 and stays within
during the 20th and 21st centuries.  It is
cseconds throughout the historical period,
enths of a degree at the end of the
e span.

al coordinates to unit vector
 to ecliptic rotation matrix, long term
 of r-matrix and p-vector
ctor to spherical coordinates
ze angle into range 0 to 2pi
ze angle into range +/- pi

e, N. and Wallace, P., 2011, New precession
r long time intervals, Astron.Astrophys. 534,

e, N. and Wallace, P., 2012, New precession
r long time intervals (Corrigendum),
, C1
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

References:

    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
    expressions, valid for long time intervals, Astron.Astrophys. 534,
    A22

    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
    expressions, valid for long time intervals (Corrigendum),
    Astron.Astrophys. 541, C1

------
ltp
------

Given:
   Julian epoch (TT)

   precession matrix, J2000.0 to date

 sense
2000,
ector with respect to the J2000.0 mean
and P_date is the same vector with respect to
nox of epoch epj.
2011, 2012) 400 millennia precession model
2006 precession at J2000.0 and stays within
during the 20th and 21st centuries.  It is
cseconds throughout the historical period,
enths of a degree at the end of the
e span.

 pole, long term
c pole, long term
product
ze vector

e, N. and Wallace, P., 2011, New precession
r long time intervals, Astron.Astrophys. 534,

e, N. and Wallace, P., 2012, New precession
r long time intervals (Corrigendum),
, C1
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The matrix is in the sense

        P_date = rp x P_J2000,

     where P_J2000 is a vector with respect to the J2000.0 mean
     equator and equinox and P_date is the same vector with respect to
     the equator and equinox of epoch epj.

  2) The Vondrak et al. (2011, 2012) 400 millennia precession model
     agrees with the IAU 2006 precession at J2000.0 and stays within
     100 microarcseconds during the 20th and 21st centuries.  It is
     accurate to a few arcseconds throughout the historical period,
     worsening to a few tenths of a degree at the end of the
     +/- 200,000 year time span.


References:

    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
    expressions, valid for long time intervals, Astron.Astrophys. 534,
    A22

    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
    expressions, valid for long time intervals (Corrigendum),
    Astron.Astrophys. 541, C1

------
ltpb
------

Given:
   Julian epoch (TT)

   precession-bias matrix, J2000.0 to date

 sense
ICRS,
ctor in the Geocentric Celestial Reference
s the vector with respect to the Celestial
ce System at that date but with nutation

bias formulation is used, of sub-
acy compared with a full 3D rotation.
2011, 2012) 400 millennia precession model
2006 precession at J2000.0 and stays within
during the 20th and 21st centuries.  It is
cseconds throughout the historical period,
enths of a degree at the end of the
e span.

e, N. and Wallace, P., 2011, New precession
r long time intervals, Astron.Astrophys. 534,

e, N. and Wallace, P., 2012, New precession
r long time intervals (Corrigendum),
, C1
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The matrix is in the sense

        P_date = rpb x P_ICRS,

     where P_ICRS is a vector in the Geocentric Celestial Reference
     System, and P_date is the vector with respect to the Celestial
     Intermediate Reference System at that date but with nutation
     neglected.

  2) A first order frame bias formulation is used, of sub-
     microarcsecond accuracy compared with a full 3D rotation.

  3) The Vondrak et al. (2011, 2012) 400 millennia precession model
     agrees with the IAU 2006 precession at J2000.0 and stays within
     100 microarcseconds during the 20th and 21st centuries.  It is
     accurate to a few arcseconds throughout the historical period,
     worsening to a few tenths of a degree at the end of the
     +/- 200,000 year time span.


References:

    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
    expressions, valid for long time intervals, Astron.Astrophys. 534,
    A22

    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
    expressions, valid for long time intervals (Corrigendum),
    Astron.Astrophys. 541, C1

------
ltpecl
------

Given:
   Julian epoch (TT)

   ecliptic pole unit vector

is with respect to the J2000.0 mean equator

2011, 2012) 400 millennia precession model
2006 precession at J2000.0 and stays within
during the 20th and 21st centuries.  It is
cseconds throughout the historical period,
enths of a degree at the end of the
e span.

e, N. and Wallace, P., 2011, New precession
r long time intervals, Astron.Astrophys. 534,

e, N. and Wallace, P., 2012, New precession
r long time intervals (Corrigendum),
, C1
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The returned vector is with respect to the J2000.0 mean equator
     and equinox.

  2) The Vondrak et al. (2011, 2012) 400 millennia precession model
     agrees with the IAU 2006 precession at J2000.0 and stays within
     100 microarcseconds during the 20th and 21st centuries.  It is
     accurate to a few arcseconds throughout the historical period,
     worsening to a few tenths of a degree at the end of the
     +/- 200,000 year time span.


References:

    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
    expressions, valid for long time intervals, Astron.Astrophys. 534,
    A22

    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
    expressions, valid for long time intervals (Corrigendum),
    Astron.Astrophys. 541, C1

------
ltpequ
------

Given:
   Julian epoch (TT)

   equator pole unit vector

is with respect to the J2000.0 mean equator

2011, 2012) 400 millennia precession model
2006 precession at J2000.0 and stays within
during the 20th and 21st centuries.  It is
cseconds throughout the historical period,
enths of a degree at the end of the
e span.

e, N. and Wallace, P., 2011, New precession
r long time intervals, Astron.Astrophys. 534,

e, N. and Wallace, P., 2012, New precession
r long time intervals (Corrigendum),
, C1
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The returned vector is with respect to the J2000.0 mean equator
     and equinox.

  2) The Vondrak et al. (2011, 2012) 400 millennia precession model
     agrees with the IAU 2006 precession at J2000.0 and stays within
     100 microarcseconds during the 20th and 21st centuries.  It is
     accurate to a few arcseconds throughout the historical period,
     worsening to a few tenths of a degree at the end of the
     +/- 200,000 year time span.


References:

    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
    expressions, valid for long time intervals, Astron.Astrophys. 534,
    A22

    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
    expressions, valid for long time intervals (Corrigendum),
    Astron.Astrophys. 541, C1

------
num00a
------

Given:
         TT as a 2-part Julian Date (Note 1)

3][3]    nutation matrix

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
in the sense V(true) = rmatn * V(mean), where
 is with respect to the true equatorial triad
ctor V(mean) is with respect to the mean
date.
ly less accurate result (about 1 mas), can be
stead the eraNum00b function.

ecession/nutation, IAU 2000A

nt to the Astronomical Almanac,
n (ed), University Science Books (1992),
4).
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix operates in the sense V(true) = rmatn * V(mean), where
     the p-vector V(true) is with respect to the true equatorial triad
     of date and the p-vector V(mean) is with respect to the mean
     equatorial triad of date.

  3) A faster, but slightly less accurate result (about 1 mas), can be
     obtained by using instead the eraNum00b function.


------
num00b
------

Given:
        TT as a 2-part Julian Date (Note 1)

3][3]   nutation matrix

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
in the sense V(true) = rmatn * V(mean), where
 is with respect to the true equatorial triad
ctor V(mean) is with respect to the mean
date.
 is faster, but slightly less accurate (about
Num00a function.

ecession/nutation, IAU 2000B

nt to the Astronomical Almanac,
n (ed), University Science Books (1992),
4).
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix operates in the sense V(true) = rmatn * V(mean), where
     the p-vector V(true) is with respect to the true equatorial triad
     of date and the p-vector V(mean) is with respect to the mean
     equatorial triad of date.

  3) The present function is faster, but slightly less accurate (about
     1 mas), than the eraNum00a function.


------
num06a
------

Given:
          TT as a 2-part Julian Date (Note 1)

[3][3]    nutation matrix

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
in the sense V(true) = rmatn * V(mean), where
 is with respect to the true equatorial triad
ctor V(mean) is with respect to the mean
date.

liquity, IAU 2006
n, IAU 2006/2000A
tation matrix

nt to the Astronomical Almanac,
n (ed), University Science Books (1992),
4).
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix operates in the sense V(true) = rmatn * V(mean), where
     the p-vector V(true) is with respect to the true equatorial triad
     of date and the p-vector V(mean) is with respect to the mean
     equatorial triad of date.


------
numat
------

Given:
       mean obliquity of date (Note 1)
       nutation (Note 2)

][3]   nutation matrix (Note 3)

liquity epsa, must be consistent with the
models from which dpsi and deps were obtained.
sible for providing the nutation components;
e and obliquity, in radians and are with
ox and ecliptic of date.
in the sense V(true) = rmatn * V(mean),
(true) is with respect to the true
date and the p-vector V(mean) is with
equatorial triad of date.

ize r-matrix to identity
around X-axis
around Z-axis

nt to the Astronomical Almanac,
n (ed), University Science Books (1992),
4).
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:


  1) The supplied mean obliquity epsa, must be consistent with the
     precession-nutation models from which dpsi and deps were obtained.

  2) The caller is responsible for providing the nutation components;
     they are in longitude and obliquity, in radians and are with
     respect to the equinox and ecliptic of date.

  3) The matrix operates in the sense V(true) = rmatn * V(mean),
     where the p-vector V(true) is with respect to the true
     equatorial triad of date and the p-vector V(mean) is with
     respect to the mean equatorial triad of date.


------
nut00a
------

Given:
   TT as a 2-part Julian Date (Note 1)

   nutation, luni-solar + planetary (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
nts in longitude and obliquity are in radians
the equinox and ecliptic of date.  The
 is assumed to be the Lieske et al. (1977)
rcsec.
and planetary nutations are included.  The
rect planetary nutations and the
 lunar and terrestrial orbits.
s the MHB2000 nutation series with the
ns for planetary nutations.  It is an
e nutation part of the IAU 2000A precession-
ally adopted by the IAU General Assembly in
 (Mathews et al. 2002), but with the free
 see Note 4) omitted.
el also contains contributions to the
de and obliquity due to the free-excitation
ation during the period 1979-2000.  These FCN
e-dependent and unpredictable, are NOT
ent function and, if required, must be
ed.  With the FCN corrections included, the
ivers a pole which is at current epochs
ndred microarcseconds.  The omission of FCN
rrors of about that size.
 provides classical nutation.  The MHB2000
h it is adapted, deals also with (i) the
GCRS and mean poles and (ii) the adjustments
iquity due to the changed precession rates.
ctions, namely frame bias and precession
ported by the ERFA functions eraBi00  and

m also provides "total" nutations, comprising
f the frame bias, precession adjustments,
and planetary nutation.  These total
d in combination with an existing IAU 1976
ation, such as eraPmat76,  to deliver GCRS-
of sub-mas accuracy at current dates.
hree shortcomings in the MHB2000 model that
ccount if more accurate or definitive results
llace 2002):
tal nutations are simply arithmetic sums,
 the various components are successive Euler
is slight lack of rigor leads to cross terms
mas after a century.  The rigorous procedure
 GCRS-to-true rotation matrix by applying the
on and nutation in that order.
recession adjustments are stated to be with
ske et al. (1977), the MHB2000 model does
ich set of Euler angles are to be used and
ments are to be applied.  The most literal
rward procedure is to adopt the 4-rotation
_A, omega_A, xi_A option, and to add DPSIPR
EPSPR to both omega_A and eps_A.
del predates the determination by Chapront
of a 14.6 mas displacement between the
quinox and the origin of the ICRS frame.  It
r, be noted that neglecting this displacement
ng star coordinates does not lead to a
e in right ascension, only a small second-
on in the pattern of the precession-nutation

he ERFA functions do not generate the "total
 though they can of course easily be
 eraBi00, eraPr00 and the present function
ts.
ntains 41 instances where the same frequency
es, of which 38 are duplicates and three are
p the present code close to the original MHB
l inefficiency has not been corrected.

omaly of the Moon
gument of the latitude of the Moon
ngitude of the Moon's ascending node
ngitude of Mercury
ngitude of Venus
ngitude of Earth
ngitude of Mars
ngitude of Jupiter
ngitude of Saturn
ngitude of Uranus
 accumulated precession in longitude

nt-Touze, M. & Francou, G. 2002,
7, 700
e, T., Fricke, W. & Morando, B. 1977,
, 1-16
ng, T.A., Buffet, B.A. 2002, J.Geophys.Res.
00 code itself was obtained on 9th September
usno.navy.mil/conv2000/chapter5/IAU2000A.
non, P., Chapront, J., Chapront-Touze, M.,
 J. 1994, Astron.Astrophys. 282, 663-683
 B., Kinoshita, H., Folgueira, M. 1999,
p.Ser. 135, 111
ware for Implementing the IAU 2000
S Workshop 5.1 (2002)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The nutation components in longitude and obliquity are in radians
     and with respect to the equinox and ecliptic of date.  The
     obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
     value of 84381.448 arcsec.

     Both the luni-solar and planetary nutations are included.  The
     latter are due to direct planetary nutations and the
     perturbations of the lunar and terrestrial orbits.

  3) The function computes the MHB2000 nutation series with the
     associated corrections for planetary nutations.  It is an
     implementation of the nutation part of the IAU 2000A precession-
     nutation model, formally adopted by the IAU General Assembly in
     2000, namely MHB2000 (Mathews et al. 2002), but with the free
     core nutation (FCN - see Note 4) omitted.

  4) The full MHB2000 model also contains contributions to the
     nutations in longitude and obliquity due to the free-excitation
     of the free-core-nutation during the period 1979-2000.  These FCN
     terms, which are time-dependent and unpredictable, are NOT
     included in the present function and, if required, must be
     independently computed.  With the FCN corrections included, the
     present function delivers a pole which is at current epochs
     accurate to a few hundred microarcseconds.  The omission of FCN
     introduces further errors of about that size.

  5) The present function provides classical nutation.  The MHB2000
     algorithm, from which it is adapted, deals also with (i) the
     offsets between the GCRS and mean poles and (ii) the adjustments
     in longitude and obliquity due to the changed precession rates.
     These additional functions, namely frame bias and precession
     adjustments, are supported by the ERFA functions eraBi00  and
     eraPr00.

  6) The MHB2000 algorithm also provides "total" nutations, comprising
     the arithmetic sum of the frame bias, precession adjustments,
     luni-solar nutation and planetary nutation.  These total
     nutations can be used in combination with an existing IAU 1976
     precession implementation, such as eraPmat76,  to deliver GCRS-
     to-true predictions of sub-mas accuracy at current dates.
     However, there are three shortcomings in the MHB2000 model that
     must be taken into account if more accurate or definitive results
     are required (see Wallace 2002):

       (i) The MHB2000 total nutations are simply arithmetic sums,
           yet in reality the various components are successive Euler
           rotations.  This slight lack of rigor leads to cross terms
           that exceed 1 mas after a century.  The rigorous procedure
           is to form the GCRS-to-true rotation matrix by applying the
           bias, precession and nutation in that order.

      (ii) Although the precession adjustments are stated to be with
           respect to Lieske et al. (1977), the MHB2000 model does
           not specify which set of Euler angles are to be used and
           how the adjustments are to be applied.  The most literal
           and straightforward procedure is to adopt the 4-rotation
           epsilon_0, psi_A, omega_A, xi_A option, and to add DPSIPR
           to psi_A and DEPSPR to both omega_A and eps_A.

     (iii) The MHB2000 model predates the determination by Chapront
           et al. (2002) of a 14.6 mas displacement between the
           J2000.0 mean equinox and the origin of the ICRS frame.  It
           should, however, be noted that neglecting this displacement
           when calculating star coordinates does not lead to a
           14.6 mas change in right ascension, only a small second-
           order distortion in the pattern of the precession-nutation
           effect.

     For these reasons, the ERFA functions do not generate the "total
     nutations" directly, though they can of course easily be
     generated by calling eraBi00, eraPr00 and the present function
     and adding the results.

  7) The MHB2000 model contains 41 instances where the same frequency
     appears multiple times, of which 38 are duplicates and three are
     triplicates.  To keep the present code close to the original MHB
     algorithm, this small inefficiency has not been corrected.


References:

     Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
     Astron.Astrophys. 387, 700

     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
     Astron.Astrophys. 58, 1-16

     Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
     107, B4.  The MHB_2000 code itself was obtained on 9th September
     2002 from ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     Astron.Astrophys.Supp.Ser. 135, 111

     Wallace, P.T., "Software for Implementing the IAU 2000
     Resolutions", in IERS Workshop 5.1 (2002)

------
nut00b
------

Given:
    TT as a 2-part Julian Date (Note 1)

    nutation, luni-solar + planetary (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
nts in longitude and obliquity are in radians
the equinox and ecliptic of date.  The
 is assumed to be the Lieske et al. (1977)
rcsec.  (The errors that result from using
he IAU 2006 value of 84381.406 arcsec can be

onsists only of luni-solar terms, but
d offset which compensates for certain long-
ms (Note 7).
implementation of the IAU 2000B abridged
lly adopted by the IAU General Assembly in
computes the MHB_2000_SHORT luni-solar
um 2001), but without the associated
precession rate adjustments and the offset
 J2000.0 mean poles.
MHB2000) nutation model contains nearly 1400
B model (McCarthy & Luzum 2003) contains only
ional simplifications, yet still delivers
uracy at present epochs.  This combination of
kes the IAU 2000B abridged nutation model
actical applications.
s a pole accurate to 1 mas from 1900 to 2100
 1 mas, very occasionally just outside
U 2000A model, which is implemented in the
q.v.), delivers considerably greater accuracy
owever, to realize this improved accuracy,
essentially unpredictable free-core-nutation
ncluded.
 provides classical nutation.  The
ithm, from which it is adapted, deals also
 between the GCRS and mean poles and (ii) the
tude and obliquity due to the changed
hese additional functions, namely frame bias
tments, are supported by the ERFA functions
.
lgorithm also provides "total" nutations,
metic sum of the frame bias, precession
ation (luni-solar + planetary).  These total
d in combination with an existing IAU 1976
ation, such as eraPmat76,  to deliver GCRS-
of mas accuracy at current epochs.  However,
e eraNut00a  function (q.v. for the reasons),
o not generate the "total nutations"
ey be required, they could of course easily
ing eraBi00, eraPr00 and the present function
ts.
includes "planetary bias" terms that are
mpensate for long-period nutations.  The
 McCarthy & Luzum (2003), namely
nd Depsilon = +1.6339 mas, are optimized for
" method described in Note 6.  The Luzum
n this ERFA implementation, namely -0.135 mas
optimized for the "rigorous" method, where
on and nutation are applied separately and in
the interval 1995-2050, the ERFA
ers a maximum error of 1.001 mas (not


e, T., Fricke, W., Morando, B., "Expressions
uantities based upon the IAU /1976/ system of
ts", Astron.Astrophys. 58, 1-2, 1-16. (1977)
ommunication, 2001 (Fortran code

um, B.J., "An abridged model of the
of the celestial pole", Cel.Mech.Dyn.Astron.

non, P., Chapront, J., Chapront-Touze, M.,
 J., Astron.Astrophys. 282, 663-683 (1994)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The nutation components in longitude and obliquity are in radians
     and with respect to the equinox and ecliptic of date.  The
     obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
     value of 84381.448 arcsec.  (The errors that result from using
     this function with the IAU 2006 value of 84381.406 arcsec can be
     neglected.)

     The nutation model consists only of luni-solar terms, but
     includes also a fixed offset which compensates for certain long-
     period planetary terms (Note 7).

  3) This function is an implementation of the IAU 2000B abridged
     nutation model formally adopted by the IAU General Assembly in
     2000.  The function computes the MHB_2000_SHORT luni-solar
     nutation series (Luzum 2001), but without the associated
     corrections for the precession rate adjustments and the offset
     between the GCRS and J2000.0 mean poles.

  4) The full IAU 2000A (MHB2000) nutation model contains nearly 1400
     terms.  The IAU 2000B model (McCarthy & Luzum 2003) contains only
     77 terms, plus additional simplifications, yet still delivers
     results of 1 mas accuracy at present epochs.  This combination of
     accuracy and size makes the IAU 2000B abridged nutation model
     suitable for most practical applications.

     The function delivers a pole accurate to 1 mas from 1900 to 2100
     (usually better than 1 mas, very occasionally just outside
     1 mas).  The full IAU 2000A model, which is implemented in the
     function eraNut00a (q.v.), delivers considerably greater accuracy
     at current dates;  however, to realize this improved accuracy,
     corrections for the essentially unpredictable free-core-nutation
     (FCN) must also be included.

  5) The present function provides classical nutation.  The
     MHB_2000_SHORT algorithm, from which it is adapted, deals also
     with (i) the offsets between the GCRS and mean poles and (ii) the
     adjustments in longitude and obliquity due to the changed
     precession rates.  These additional functions, namely frame bias
     and precession adjustments, are supported by the ERFA functions
     eraBi00  and eraPr00.

  6) The MHB_2000_SHORT algorithm also provides "total" nutations,
     comprising the arithmetic sum of the frame bias, precession
     adjustments, and nutation (luni-solar + planetary).  These total
     nutations can be used in combination with an existing IAU 1976
     precession implementation, such as eraPmat76,  to deliver GCRS-
     to-true predictions of mas accuracy at current epochs.  However,
     for symmetry with the eraNut00a  function (q.v. for the reasons),
     the ERFA functions do not generate the "total nutations"
     directly.  Should they be required, they could of course easily
     be generated by calling eraBi00, eraPr00 and the present function
     and adding the results.

  7) The IAU 2000B model includes "planetary bias" terms that are
     fixed in size but compensate for long-period nutations.  The
     amplitudes quoted in McCarthy & Luzum (2003), namely
     Dpsi = -1.5835 mas and Depsilon = +1.6339 mas, are optimized for
     the "total nutations" method described in Note 6.  The Luzum
     (2001) values used in this ERFA implementation, namely -0.135 mas
     and +0.388 mas, are optimized for the "rigorous" method, where
     frame bias, precession and nutation are applied separately and in
     that order.  During the interval 1995-2050, the ERFA
     implementation delivers a maximum error of 1.001 mas (not
     including FCN).


References:

     Lieske, J.H., Lederle, T., Fricke, W., Morando, B., "Expressions
     for the precession quantities based upon the IAU /1976/ system of
     astronomical constants", Astron.Astrophys. 58, 1-2, 1-16. (1977)

     Luzum, B., private communication, 2001 (Fortran code
     MHB_2000_SHORT)

     McCarthy, D.D. & Luzum, B.J., "An abridged model of the
     precession-nutation of the celestial pole", Cel.Mech.Dyn.Astron.
     85, 37-49 (2003)

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J., Astron.Astrophys. 282, 663-683 (1994)

------
nut06a
------

Given:
   TT as a 2-part Julian Date (Note 1)

   nutation, luni-solar + planetary (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
nts in longitude and obliquity are in radians
the mean equinox and ecliptic of date,
model (Hilton et al. 2006, Capitaine et al.

omputes the IAU 2000A nutation, then applies
the consequences of the change in obliquity
liptic to the IAU 2006 ecliptic and (ii) the
 the Earth's dynamical form factor J2.
 provides classical nutation, complementing
ias and IAU 2006 precession.  It delivers a
rent epochs accurate to a few tens of
rt from the free core nutation.

n, IAU 2000A

nt-Touze, M. & Francou, G. 2002,
7, 700
e, T., Fricke, W. & Morando, B. 1977,
, 1-16
ng, T.A., Buffet, B.A. 2002, J.Geophys.Res.
00 code itself was obtained on 9th September
usno.navy.mil/conv2000/chapter5/IAU2000A.
non, P., Chapront, J., Chapront-Touze, M.,
 J. 1994, Astron.Astrophys. 282, 663-683
 B., Kinoshita, H., Folgueira, M. 1999,
p.Ser. 135, 111
ware for Implementing the IAU 2000
S Workshop 5.1 (2002)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The nutation components in longitude and obliquity are in radians
     and with respect to the mean equinox and ecliptic of date,
     IAU 2006 precession model (Hilton et al. 2006, Capitaine et al.
     2005).

  3) The function first computes the IAU 2000A nutation, then applies
     adjustments for (i) the consequences of the change in obliquity
     from the IAU 1980 ecliptic to the IAU 2006 ecliptic and (ii) the
     secular variation in the Earth's dynamical form factor J2.

  4) The present function provides classical nutation, complementing
     the IAU 2000 frame bias and IAU 2006 precession.  It delivers a
     pole which is at current epochs accurate to a few tens of
     microarcseconds, apart from the free core nutation.


References:

     Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
     Astron.Astrophys. 387, 700

     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
     Astron.Astrophys. 58, 1-16

     Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
     107, B4.  The MHB_2000 code itself was obtained on 9th September
     2002 from ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.

     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
     Astron.Astrophys.Supp.Ser. 135, 111

     Wallace, P.T., "Software for Implementing the IAU 2000
     Resolutions", in IERS Workshop 5.1 (2002)

------
nut80
------

Given:
    TT as a 2-part Julian Date (Note 1)

    nutation in longitude (radians)
    nutation in obliquity (radians)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
nts are with respect to the ecliptic of


ze angle into range +/- pi

nt to the Astronomical Almanac,
n (ed), University Science Books (1992),
.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The nutation components are with respect to the ecliptic of
     date.


------
nutm80
------

Given:
e          TDB date (Note 1)

e[3][3]    nutation matrix

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
in the sense V(true) = rmatn * V(mean),
(true) is with respect to the true
date and the p-vector V(mean) is with
equatorial triad of date.

n, IAU 1980
liquity, IAU 1980
tation matrix
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix operates in the sense V(true) = rmatn * V(mean),
     where the p-vector V(true) is with respect to the true
     equatorial triad of date and the p-vector V(mean) is with
     respect to the mean equatorial triad of date.


------
obl06
------

Given:
  TT as a 2-part Julian Date (Note 1)
e):
  obliquity of the ecliptic (radians, Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
gle between the ecliptic and mean equator of


006, Celest.Mech.Dyn.Astron. 94, 351
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
obl80
------

Given:
    TT as a 2-part Julian Date (Note 1)
e):
    obliquity of the ecliptic (radians, Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
gle between the ecliptic and mean equator of


nt to the Astronomical Almanac,
n (ed), University Science Books (1992),
p114).
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
p06e
------

Given:
   TT as a 2-part Julian Date (Note 1)

   epsilon_0
   psi_A
   omega_A
   P_A
   Q_A
   pi_A
   Pi_A
   obliquity epsilon_A
   chi_A
   z_A
   zeta_A
   theta_A
   p_A
   F-W angle gamma_J2000
   F-W angle phi_J2000
   F-W angle psi_J2000

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
s the set of equinox based angles for the
3" precession theory, adopted by the IAU in
e set out in Table 1 of Hilton et al. (2006):
bliquity at J2000.0
uni-solar precession
nclination of equator wrt J2000.0 ecliptic
cliptic pole x, J2000.0 ecliptic triad
cliptic pole -y, J2000.0 ecliptic triad
ngle between moving and J2000.0 ecliptics
ongitude of ascending node of the ecliptic
bliquity of the ecliptic
lanetary precession
quatorial precession: -3rd 323 Euler angle
quatorial precession: -1st 323 Euler angle
quatorial precession: 2nd 323 Euler angle
eneral precession
2000.0 RA difference of ecliptic poles
2000.0 codeclination of ecliptic pole
ongitude difference of equator poles, J2000.0
are all radians.
 Table 1 also contains angles that depend on
 the P03 precession theory itself, namely the
 and nutation.  The quoted polynomials are
unctions:
the polynomial parts of the X and Y series.
he polynomial part of the s+XY/2 series.
ts the series for the Fukushima-Williams
th respect to the GCRS pole (i.e. the variants
 bias).
tipulated that the choice of parameterization
, and so an IAU compliant precession
e constructed using various combinations of
by the present function.
 used by ERFA is the version of the Fukushima-
 refers directly to the GCRS pole.  These
ated by calling the function eraPfw06.  ERFA
rect computation of the CIP GCRS X,Y by
 calling eraXy06.
n the different parameterizations is at the
el in the present era.
precession formulation that refers to the GCRS
 dynamical pole, it may (depending on the
 necessary to introduce the frame bias

 re-use the same variable in the returned
tities are stored in the stated order.

006, Celest.Mech.Dyn.Astron. 94, 351

liquity, IAU 2006
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) This function returns the set of equinox based angles for the
     Capitaine et al. "P03" precession theory, adopted by the IAU in
     2006.  The angles are set out in Table 1 of Hilton et al. (2006):

     eps0   epsilon_0   obliquity at J2000.0
     psia   psi_A       luni-solar precession
     oma    omega_A     inclination of equator wrt J2000.0 ecliptic
     bpa    P_A         ecliptic pole x, J2000.0 ecliptic triad
     bqa    Q_A         ecliptic pole -y, J2000.0 ecliptic triad
     pia    pi_A        angle between moving and J2000.0 ecliptics
     bpia   Pi_A        longitude of ascending node of the ecliptic
     epsa   epsilon_A   obliquity of the ecliptic
     chia   chi_A       planetary precession
     za     z_A         equatorial precession: -3rd 323 Euler angle
     zetaa  zeta_A      equatorial precession: -1st 323 Euler angle
     thetaa theta_A     equatorial precession: 2nd 323 Euler angle
     pa     p_A         general precession
     gam    gamma_J2000 J2000.0 RA difference of ecliptic poles
     phi    phi_J2000   J2000.0 codeclination of ecliptic pole
     psi    psi_J2000   longitude difference of equator poles, J2000.0

     The returned values are all radians.

  3) Hilton et al. (2006) Table 1 also contains angles that depend on
     models distinct from the P03 precession theory itself, namely the
     IAU 2000A frame bias and nutation.  The quoted polynomials are
     used in other ERFA functions:

     . eraXy06  contains the polynomial parts of the X and Y series.

     . eraS06  contains the polynomial part of the s+XY/2 series.

     . eraPfw06  implements the series for the Fukushima-Williams
       angles that are with respect to the GCRS pole (i.e. the variants
       that include frame bias).

  4) The IAU resolution stipulated that the choice of parameterization
     was left to the user, and so an IAU compliant precession
     implementation can be constructed using various combinations of
     the angles returned by the present function.

  5) The parameterization used by ERFA is the version of the Fukushima-
     Williams angles that refers directly to the GCRS pole.  These
     angles may be calculated by calling the function eraPfw06.  ERFA
     also supports the direct computation of the CIP GCRS X,Y by
     series, available by calling eraXy06.

  6) The agreement between the different parameterizations is at the
     1 microarcsecond level in the present era.

  7) When constructing a precession formulation that refers to the GCRS
     pole rather than the dynamical pole, it may (depending on the
     choice of angles) be necessary to introduce the frame bias
     explicitly.

  8) It is permissible to re-use the same variable in the returned
     arguments.  The quantities are stored in the stated order.

  Reference:

     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351


------
p2pv
------

Given:
     p-vector

]    pv-vector

vector
vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
p2s
------

Given:
  p-vector

  longitude angle (radians)
  latitude angle (radians)
  radial distance

heta, phi and r are returned.
 theta is returned.

r to spherical
 of p-vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) If P is null, zero theta, phi and r are returned.

  2) At either pole, zero theta is returned.


------
pap
------

Given:
rection of reference point
rection of point whose PA is required
e):
sition angle of b with respect to a (radians)

sition angle, in radians, of direction b with
 a.  It is in the range -pi to +pi.  The
f b is a small distance "north" of a the
proximately zero, and if b is a small
 the position angle is approximately +pi/2.
need not be of unit length.
the two directions are the same or if either

pole, the result is ill-defined.

se p-vector into modulus and direction
 of p-vector
product of two p-vectors
r minus p-vector
product of two p-vectors
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The result is the position angle, in radians, of direction b with
     respect to direction a.  It is in the range -pi to +pi.  The
     sense is such that if b is a small distance "north" of a the
     position angle is approximately zero, and if b is a small
     distance "east" of a the position angle is approximately +pi/2.

  2) The vectors a and b need not be of unit length.

  3) Zero is returned if the two directions are the same or if either
     vector is null.

  4) If vector a is at a pole, the result is ill-defined.


------
pas
------

Given:
ngitude of point A (e.g. RA) in radians
titude of point A (e.g. Dec) in radians
ngitude of point B
titude of point B
e):
sition angle of B with respect to A

aring (position angle), in radians, of point
int A.  It is in the range -pi to +pi.  The
f B is a small distance "east" of point A,
ximately +pi/2.
the two points are coincident.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
pb06
------

Given:
  TT as a 2-part Julian Date (Note 1)

  1st rotation: radians cw around z
  3rd rotation: radians cw around z
  2nd rotation: radians ccw around y

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
mulated precession angles zeta_A, z_A,
tained in the usual way, namely through
ns, because of the frame bias.  The latter
e angles undergo rapid changes near this
ead the results of decomposing the
ix obtained by using the Fukushima-Williams
ot suffer from the problem.  The
s values which can be used in the
tion and which include frame bias.
 returned in the conventional order, which
he order of the corresponding Euler
ession-bias matrix is
) x R_3(-zeta).
theta_A angles be required that do not
they are available by calling the ERFA


ix, IAU 2006
around Z-axis
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The traditional accumulated precession angles zeta_A, z_A,
     theta_A cannot be obtained in the usual way, namely through
     polynomial expressions, because of the frame bias.  The latter
     means that two of the angles undergo rapid changes near this
     date.  They are instead the results of decomposing the
     precession-bias matrix obtained by using the Fukushima-Williams
     method, which does not suffer from the problem.  The
     decomposition returns values which can be used in the
     conventional formulation and which include frame bias.

  3) The three angles are returned in the conventional order, which
     is not the same as the order of the corresponding Euler
     rotations.  The precession-bias matrix is
     R_3(-z) x R_2(+theta) x R_3(-zeta).

  4) Should zeta_A, z_A, theta_A angles be required that do not
     contain frame bias, they are available by calling the ERFA
     function eraP06e.


------
pdp
------

Given:
 first p-vector
 second p-vector
e):
 a . b
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
pfw06
------

Given:
  TT as a 2-part Julian Date (Note 1)

  F-W angle gamma_bar (radians)
  F-W angle phi_bar (radians)
  F-W angle psi_bar (radians)
  F-W angle epsilon_A (radians)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
 points:
liptic pole,

tic pole of date,
of date,
illiams angles are as follows:
= epE
pE
pEP
= EP
ing the combined effects of frame bias and

R_3(-psib).R_1(phib).R_3(gamb)
ing the combined effects of frame bias,
ion is simply:
-dE).R_3(-psib-dP).R_1(phib).R_3(gamb)
the nutation components with respect to the


006, Celest.Mech.Dyn.Astron. 94, 351

liquity, IAU 2006
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) Naming the following points:

           e = J2000.0 ecliptic pole,
           p = GCRS pole,
           E = mean ecliptic pole of date,
     and   P = mean pole of date,

     the four Fukushima-Williams angles are as follows:

        gamb = gamma_bar = epE
        phib = phi_bar = pE
        psib = psi_bar = pEP
        epsa = epsilon_A = EP

  3) The matrix representing the combined effects of frame bias and
     precession is:

        PxB = R_1(-epsa).R_3(-psib).R_1(phib).R_3(gamb)

  4) The matrix representing the combined effects of frame bias,
     precession and nutation is simply:

        NxPxB = R_1(-epsa-dE).R_3(-psib-dP).R_1(phib).R_3(gamb)

     where dP and dE are the nutation components with respect to the
     ecliptic of date.

  Reference:

     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351


------
plan94
------

Given:
TDB date part A (Note 1)
TDB date part B (Note 1)
planet (1=Mercury, 2=Venus, 3=EMB, 4=Mars,
    5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune)

planet p,v (heliocentric, J2000.0, AU,AU/d)
e):
status: -1 = illegal NP (outside 1-8)
         0 = OK
        +1 = warning: year outside 1000-3000
        +2 = warning: failed to converge

 is in the TDB time scale (in practice TT can
lian Date, apportioned in any convenient way
ments.  For example, JD(TDB)=2450123.7 could
of these ways, among others:
  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 method is best matched to the way the
internally and will deliver the optimum
 method and the date & time methods are both
ween resolution and convenience.  The limited
ent algorithm is such that any of the methods

de the range 1-8 is supplied, an error status
is returned and the pv vector set to zeroes.
is for the Earth-Moon Barycenter.  To obtain
ition and velocity of the Earth, use instead
aEpv00.
, the array pv contains the following:
 }
 } heliocentric position, AU
 }
 }
 } heliocentric velocity, AU/d
 }
is equatorial and is with respect to the
inox of epoch J2000.0.
 to J.L. Simon, P. Bretagnon, J. Chapront,
. Francou and J. Laskar (Bureau des
rance).  From comparisons with JPL
y quote the following maximum errors
00-2050:
rcsec)    B (arcsec)      R (km)
4             1             300
5             1             800
6             1            1000
7             1            7700
1             5           76000
1            13          267000
6             7          712000
1             1          253000
00-3000, they report that the accuracy is no
 that over 1800-2050.  Outside 1000-3000 the

resent function with the JPL DE200 ephemeris
MS errors over the interval 1960-2025:
ition (km)     velocity (m/s)
  334               0.437
 1060               0.855
 2010               0.815
 7690               1.98
71700               7.70
99000              19.4
64000              16.4
58000              14.4
DE200 over the interval 1800-2100 gave the
solute differences.  (The results using
ly the same.)
sec)   B (arcsec)     R (km)   Rdot (m/s)
           1            500       0.7
           1           1100       0.9
           1           1300       1.0
           1           9000       2.5
           6          82000       8.2
          14         263000      24.6
           7         661000      27.4
           2         248000      21.4
implementation of the original Simon et al.
 from the original in the following respects:
rtran.
plied in two parts.
eturned only in equatorial Cartesian form;
ngitude, latitude and radius vector are not

n the J2000.0 equatorial frame, not ecliptic.
-line: there are fewer calls to subroutines.
/warning status values are used.
ler's-equation-solver is used (avoiding
recision complex).
t are nested to minimize rounding errors.
 constants are used to avoid mixed-mode

anges affects the result significantly.
indicates the most serious condition
xecution of the function.  Illegal np is
serious, overriding failure to converge,
precedence over the remote date warning.

ze angle into range 0 to 2pi
 Bretagnon, P., Chapront, J.,
uze, M., Francou, G., and Laskar, J.,
rophys. 282, 663 (1994).
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The date date1+date2 is in the TDB time scale (in practice TT can
     be used) and is a Julian Date, apportioned in any convenient way
     between the two arguments.  For example, JD(TDB)=2450123.7 could
     be expressed in any of these ways, among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in cases
     where the loss of several decimal digits of resolution is
     acceptable.  The J2000 method is best matched to the way the
     argument is handled internally and will deliver the optimum
     resolution.  The MJD method and the date & time methods are both
     good compromises between resolution and convenience.  The limited
     accuracy of the present algorithm is such that any of the methods
     is satisfactory.

  2) If an np value outside the range 1-8 is supplied, an error status
     (function value -1) is returned and the pv vector set to zeroes.

  3) For np=3 the result is for the Earth-Moon Barycenter.  To obtain
     the heliocentric position and velocity of the Earth, use instead
     the ERFA function eraEpv00.

  4) On successful return, the array pv contains the following:

        pv[0][0]   x      }
        pv[0][1]   y      } heliocentric position, AU
        pv[0][2]   z      }

        pv[1][0]   xdot   }
        pv[1][1]   ydot   } heliocentric velocity, AU/d
        pv[1][2]   zdot   }

     The reference frame is equatorial and is with respect to the
     mean equator and equinox of epoch J2000.0.

  5) The algorithm is due to J.L. Simon, P. Bretagnon, J. Chapront,
     M. Chapront-Touze, G. Francou and J. Laskar (Bureau des
     Longitudes, Paris, France).  From comparisons with JPL
     ephemeris DE102, they quote the following maximum errors
     over the interval 1800-2050:

                     L (arcsec)    B (arcsec)      R (km)

        Mercury          4             1             300
        Venus            5             1             800
        EMB              6             1            1000
        Mars            17             1            7700
        Jupiter         71             5           76000
        Saturn          81            13          267000
        Uranus          86             7          712000
        Neptune         11             1          253000

     Over the interval 1000-3000, they report that the accuracy is no
     worse than 1.5 times that over 1800-2050.  Outside 1000-3000 the
     accuracy declines.

     Comparisons of the present function with the JPL DE200 ephemeris
     give the following RMS errors over the interval 1960-2025:

                      position (km)     velocity (m/s)

        Mercury            334               0.437
        Venus             1060               0.855
        EMB               2010               0.815
        Mars              7690               1.98
        Jupiter          71700               7.70
        Saturn          199000              19.4
        Uranus          564000              16.4
        Neptune         158000              14.4

     Comparisons against DE200 over the interval 1800-2100 gave the
     following maximum absolute differences.  (The results using
     DE406 were essentially the same.)

                   L (arcsec)   B (arcsec)     R (km)   Rdot (m/s)

        Mercury        7            1            500       0.7
        Venus          7            1           1100       0.9
        EMB            9            1           1300       1.0
        Mars          26            1           9000       2.5
        Jupiter       78            6          82000       8.2
        Saturn        87           14         263000      24.6
        Uranus        86            7         661000      27.4
        Neptune       11            2         248000      21.4

  6) The present ERFA re-implementation of the original Simon et al.
     Fortran code differs from the original in the following respects:

       *  C instead of Fortran.

       *  The date is supplied in two parts.

       *  The result is returned only in equatorial Cartesian form;
          the ecliptic longitude, latitude and radius vector are not
          returned.

       *  The result is in the J2000.0 equatorial frame, not ecliptic.

       *  More is done in-line: there are fewer calls to subroutines.

       *  Different error/warning status values are used.

       *  A different Kepler's-equation-solver is used (avoiding
          use of double precision complex).

       *  Polynomials in t are nested to minimize rounding errors.

       *  Explicit double constants are used to avoid mixed-mode
          expressions.

     None of the above changes affects the result significantly.

  7) The returned status indicates the most serious condition
     encountered during execution of the function.  Illegal np is
     considered the most serious, overriding failure to converge,
     which in turn takes precedence over the remote date warning.


------
pm
------

Given:
 p-vector
e):
 modulus
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
pmat00
------

Given:
         TT as a 2-part Julian Date (Note 1)

3][3]    bias-precession matrix (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
in the sense V(date) = rbp * V(GCRS), where
 is with respect to the Geocentric Celestial
U, 2000) and the p-vector V(date) is with
equatorial triad of the given date.

ias and precession matrices, IAU 2000

ional Astronomical Union, Vol. XXIVB;  Proc.
y, Manchester, UK.  Resolutions B1.3, B1.6.

, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix operates in the sense V(date) = rbp * V(GCRS), where
     the p-vector V(GCRS) is with respect to the Geocentric Celestial
     Reference System (IAU, 2000) and the p-vector V(date) is with
     respect to the mean equatorial triad of the given date.


------
pmat06
------

Given:
         TT as a 2-part Julian Date (Note 1)

3][3]    bias-precession matrix (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
in the sense V(date) = rbp * V(GCRS), where
 is with respect to the Geocentric Celestial
U, 2000) and the p-vector V(date) is with
equatorial triad of the given date.

ecession F-W angles, IAU 2006
les to r-matrix

ace, P.T., 2006, Astron.Astrophys. 450, 855
taine, N., 2006, Astron.Astrophys. 459, 981
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix operates in the sense V(date) = rbp * V(GCRS), where
     the p-vector V(GCRS) is with respect to the Geocentric Celestial
     Reference System (IAU, 2000) and the p-vector V(date) is with
     respect to the mean equatorial triad of the given date.


References:

     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981

------
pmat76
------

Given:
     ending date, TT (Note 1)

][3] precession matrix, J2000.0 -> date1+date2

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
in the sense V(date) = RMATP * V(J2000),
(J2000) is with respect to the mean
epoch J2000.0 and the p-vector V(date)
he mean equatorial triad of the given

thod itself is rigorous, the precession
 through canonical polynomials which are
ited time span.  In addition, the IAU 1976
nown to be imperfect.  The absolute accuracy
lation is better than 0.1 arcsec from
tter than 1 arcsec from 1640AD to 2360AD,
arcsec for the whole of the period
e errors exceed 10 arcsec outside the
AD, exceed 100 arcsec outside 4200BC to
00 arcsec outside 6800BC to 8200AD.

ated precession angles, IAU 1976
ize r-matrix to identity
around Z-axis
around Y-axis
matrix

Astron.Astrophys. 73, 282.
, p283.
SNO circular no. 163, pA2.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix operates in the sense V(date) = RMATP * V(J2000),
     where the p-vector V(J2000) is with respect to the mean
     equatorial triad of epoch J2000.0 and the p-vector V(date)
     is with respect to the mean equatorial triad of the given
     date.

  3) Though the matrix method itself is rigorous, the precession
     angles are expressed through canonical polynomials which are
     valid only for a limited time span.  In addition, the IAU 1976
     precession rate is known to be imperfect.  The absolute accuracy
     of the present formulation is better than 0.1 arcsec from
     1960AD to 2040AD, better than 1 arcsec from 1640AD to 2360AD,
     and remains below 3 arcsec for the whole of the period
     500BC to 3000AD.  The errors exceed 10 arcsec outside the
     range 1200BC to 3900AD, exceed 100 arcsec outside 4200BC to
     5600AD and exceed 1000 arcsec outside 6800BC to 8200AD.


References:

     Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
      equations (6) & (7), p283.

     Kaplan,G.H., 1981. USNO circular no. 163, pA2.

------
pmp
------

Given:
    first p-vector
    second p-vector

    a - b

 re-use the same array for any of the

, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
pmpx
------

Given:
RS RA,Dec at catalog epoch (radians)
 proper motion (radians/year; Note 1)
c proper motion (radians/year)
rallax (arcsec)
dial velocity (km/s, +ve if receding)
oper motion time interval (SSB, Julian years)
B to observer vector (au)

ordinate direction (BCRS unit vector)

 RA is dRA/dt rather than cos(Dec)*dRA/dt.
me interval is for when the starlight
stem barycenter.
r iteration, the Roemer effect (i.e. the
ion of the proper motion coming from the
 is applied approximately, using the
r at the catalog epoch.

manac, pp B39-B41.
nn, P. K. (eds), Explanatory Supplement to
anac, 3rd ed., University Science Books


product of two p-vectors
se p-vector into modulus and direction
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

  2) The proper motion time interval is for when the starlight
     reaches the solar system barycenter.

  3) To avoid the need for iteration, the Roemer effect (i.e. the
     small annual modulation of the proper motion coming from the
     changing light time) is applied approximately, using the
     direction of the star at the catalog epoch.

  References:

     1984 Astronomical Almanac, pp B39-B41.

     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
     the Astronomical Almanac, 3rd ed., University Science Books
     (2013), Section 7.2.


References:

     1984 Astronomical Almanac, pp B39-B41.

     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
     the Astronomical Almanac, 3rd ed., University Science Books
     (2013), Section 7.2.

  Called:
     eraPdp       scalar product of two p-vectors
     eraPn        decompose p-vector into modulus and direction

------
pmsafe
------

Given:
ight ascension (radians), before
eclination (radians), before
A proper motion (radians/year), before
ec proper motion (radians/year), before
arallax (arcseconds), before
adial velocity (km/s, +ve = receding), before
before" epoch, part A (Note 1)
before" epoch, part B (Note 1)
after" epoch, part A (Note 1)
after" epoch, part B (Note 1)

ight ascension (radians), after
eclination (radians), after
A proper motion (radians/year), after
ec proper motion (radians/year), after
arallax (arcseconds), after
adial velocity (km/s, +ve = receding), after
e):
tatus:
-1 = system error (should not occur)
 0 = no warnings or errors
 1 = distance overridden (Note 6)
 2 = excessive velocity (Note 7)
 4 = solution didn't converge (Note 8)
se = binary logical OR of the above warnings

ing TDB epochs ep1a+ep1b and ep2a+ep2b are
ioned in any convenient way between the two
r example, JD(TDB)=2450123.7 could be
these ways, among others:
   epNb
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 method is best matched to the way the
internally and will deliver the optimum
 method and the date & time methods are both
ween resolution and convenience.
ormal star-catalog conventions, the object's
declination are freed from the effects of
 The frame, which is aligned to the catalog
 is Lorentzian and centered on the SSB.
re the rate of change of the right ascension
he catalog epoch and are in radians per TDB

ial velocity are in the same frame.
units.  The star coordinates are in radians
ns in radians per Julian year, but the
conds.
 is in terms of coordinate angle, not true
og uses arcseconds for both RA and Dec proper
er motion will need to be divided by cos(Dec)

 at constant speed, in the inertial frame, is

or zero or negative) parallax is overridden
bject is at a finite but very large distance,
t the proper motion is equivalent to a large
t 0.1c using the chosen constant).  A warning
 to the status if this action has been taken.
y is a significant fraction of c (see the
 function eraStarpv), it is arbitrarily set
action occurs, 2 is added to the status.
ustment carried out in the eraStarpv function
e calculation.  If the process fails to
t number of iterations, 4 is added to the


etween two points
star catalog data for space motion
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The starting and ending TDB epochs ep1a+ep1b and ep2a+ep2b are
     Julian Dates, apportioned in any convenient way between the two
     parts (A and B).  For example, JD(TDB)=2450123.7 could be
     expressed in any of these ways, among others:

            epNa            epNb

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in cases
     where the loss of several decimal digits of resolution is
     acceptable.  The J2000 method is best matched to the way the
     argument is handled internally and will deliver the optimum
     resolution.  The MJD method and the date & time methods are both
     good compromises between resolution and convenience.

  2) In accordance with normal star-catalog conventions, the object's
     right ascension and declination are freed from the effects of
     secular aberration.  The frame, which is aligned to the catalog
     equator and equinox, is Lorentzian and centered on the SSB.

     The proper motions are the rate of change of the right ascension
     and declination at the catalog epoch and are in radians per TDB
     Julian year.

     The parallax and radial velocity are in the same frame.

  3) Care is needed with units.  The star coordinates are in radians
     and the proper motions in radians per Julian year, but the
     parallax is in arcseconds.

  4) The RA proper motion is in terms of coordinate angle, not true
     angle.  If the catalog uses arcseconds for both RA and Dec proper
     motions, the RA proper motion will need to be divided by cos(Dec)
     before use.

  5) Straight-line motion at constant speed, in the inertial frame, is
     assumed.

  6) An extremely small (or zero or negative) parallax is overridden
     to ensure that the object is at a finite but very large distance,
     but not so large that the proper motion is equivalent to a large
     but safe speed (about 0.1c using the chosen constant).  A warning
     status of 1 is added to the status if this action has been taken.

  7) If the space velocity is a significant fraction of c (see the
     constant VMAX in the function eraStarpv), it is arbitrarily set
     to zero.  When this action occurs, 2 is added to the status.

  8) The relativistic adjustment carried out in the eraStarpv function
     involves an iterative calculation.  If the process fails to
     converge within a set number of iterations, 4 is added to the
     status.


------
pn
------

Given:
    p-vector

    modulus
    unit vector

sult is null.  Otherwise the result is a unit

 re-use the same array for any of the


 of p-vector
vector
y p-vector by scalar
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) If p is null, the result is null.  Otherwise the result is a unit
     vector.

  2) It is permissible to re-use the same array for any of the
     arguments.


------
pn00
------

Given:
         TT as a 2-part Julian Date (Note 1)
         nutation (Note 2)

         mean obliquity (Note 3)
3][3]    frame bias matrix (Note 4)
3][3]    precession matrix (Note 5)
3][3]    bias-precession matrix (Note 6)
3][3]    nutation matrix (Note 7)
3][3]    GCRS-to-true matrix (Note 8)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
sible for providing the nutation components;
e and obliquity, in radians and are with
ox and ecliptic of date.  For high-accuracy
ore nutation should be included as well as
orrections to the position of the CIP.
liquity is consistent with the IAU 2000
models.
orms vectors from GCRS to J2000.0 mean
by applying frame bias.
orms vectors from J2000.0 mean equator and
tor and equinox of date by applying

forms vectors from GCRS to mean equator and
pplying frame bias then precession.  It is

orms vectors from mean equator and equinox of
 and equinox of date by applying the nutation
ary).
sforms vectors from GCRS to true equator and
 is the product rn x rbp, applying frame
 nutation in that order.
 re-use the same array in the returned
ys are filled in the order given.

0 precession adjustments
liquity, IAU 1980
ias and precession matrices, IAU 2000
matrix
tation matrix
 of two r-matrices

ont, J., Lambert, S. and Wallace, P.,
 Celestial Intermediate Pole and Celestial
sistent with the IAU 2000A precession-
ron.Astrophys. 400, 1145-1154 (2003)
phemeris origin (CEO) was renamed "celestial
igin" (CIO) by IAU 2006 Resolution 2.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The caller is responsible for providing the nutation components;
     they are in longitude and obliquity, in radians and are with
     respect to the equinox and ecliptic of date.  For high-accuracy
     applications, free core nutation should be included as well as
     any other relevant corrections to the position of the CIP.

  3) The returned mean obliquity is consistent with the IAU 2000
     precession-nutation models.

  4) The matrix rb transforms vectors from GCRS to J2000.0 mean
     equator and equinox by applying frame bias.

  5) The matrix rp transforms vectors from J2000.0 mean equator and
     equinox to mean equator and equinox of date by applying
     precession.

  6) The matrix rbp transforms vectors from GCRS to mean equator and
     equinox of date by applying frame bias then precession.  It is
     the product rp x rb.

  7) The matrix rn transforms vectors from mean equator and equinox of
     date to true equator and equinox of date by applying the nutation
     (luni-solar + planetary).

  8) The matrix rbpn transforms vectors from GCRS to true equator and
     equinox of date.  It is the product rn x rbp, applying frame
     bias, precession and nutation in that order.

  9) It is permissible to re-use the same array in the returned
     arguments.  The arrays are filled in the order given.


------
pn00a
------

Given:
         TT as a 2-part Julian Date (Note 1)

         nutation (Note 2)
         mean obliquity (Note 3)
3][3]    frame bias matrix (Note 4)
3][3]    precession matrix (Note 5)
3][3]    bias-precession matrix (Note 6)
3][3]    nutation matrix (Note 7)
3][3]    GCRS-to-true matrix (Notes 8,9)

ate2 is a Julian Date, apportioned in any
een the two arguments.  For example,
uld be expressed in any of these ways,

   date2
     0.0       (JD method)
 -1421.3       (J2000 method)
 50123.2       (MJD method)
     0.2       (date & time method)
e most natural and convenient to use in
s of several decimal digits of resolution
 J2000 method is best matched to the way
dled internally and will deliver the
  The MJD method and the date & time methods
omises between resolution and convenience.
ents (luni-solar + planetary, IAU 2000A) in
uity are in radians and with respect to the
c of date.  Free core nutation is omitted;
racy, use the eraPn00  function, where the
 are caller-specified.  For faster but
ate results, use the eraPn00b function.
is consistent with the IAU 2000 precession.
forms vectors from GCRS to J2000.0 mean
 by applying frame bias.
forms vectors from J2000.0 mean equator and
ator and equinox of date by applying

sforms vectors from GCRS to mean equator and
applying frame bias then precession.  It is
.
forms vectors from mean equator and equinox
ator and equinox of date by applying the
r + planetary).
nsforms vectors from GCRS to true equator and
t is the product rn x rbp, applying frame
d nutation in that order.
es of the IAU 2000A Celestial Intermediate
3,1-3) of the GCRS-to-true matrix,

o re-use the same array in the returned
ays are filled in the order given.

n, IAU 2000A
ecession/nutation results, IAU 2000

ont, J., Lambert, S. and Wallace, P.,
 Celestial Intermediate Pole and Celestial
sistent with the IAU 2000A precession-
ron.Astrophys. 400, 1145-1154 (2003)
phemeris origin (CEO) was renamed "celestial
igin" (CIO) by IAU 2006 Resolution 2.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1)  The TT date date1+date2 is a Julian Date, apportioned in any
      convenient way between the two arguments.  For example,
      JD(TT)=2450123.7 could be expressed in any of these ways,
      among others:

             date1          date2

          2450123.7           0.0       (JD method)
          2451545.0       -1421.3       (J2000 method)
          2400000.5       50123.2       (MJD method)
          2450123.5           0.2       (date & time method)

      The JD method is the most natural and convenient to use in
      cases where the loss of several decimal digits of resolution
      is acceptable.  The J2000 method is best matched to the way
      the argument is handled internally and will deliver the
      optimum resolution.  The MJD method and the date & time methods
      are both good compromises between resolution and convenience.

  2)  The nutation components (luni-solar + planetary, IAU 2000A) in
      longitude and obliquity are in radians and with respect to the
      equinox and ecliptic of date.  Free core nutation is omitted;
      for the utmost accuracy, use the eraPn00  function, where the
      nutation components are caller-specified.  For faster but
      slightly less accurate results, use the eraPn00b function.

  3)  The mean obliquity is consistent with the IAU 2000 precession.

  4)  The matrix rb transforms vectors from GCRS to J2000.0 mean
      equator and equinox by applying frame bias.

  5)  The matrix rp transforms vectors from J2000.0 mean equator and
      equinox to mean equator and equinox of date by applying
      precession.

  6)  The matrix rbp transforms vectors from GCRS to mean equator and
      equinox of date by applying frame bias then precession.  It is
      the product rp x rb.

  7)  The matrix rn transforms vectors from mean equator and equinox
      of date to true equator and equinox of date by applying the
      nutation (luni-solar + planetary).

  8)  The matrix rbpn transforms vectors from GCRS to true equator and
      equinox of date.  It is the product rn x rbp, applying frame
      bias, precession and nutation in that order.

  9)  The X,Y,Z coordinates of the IAU 2000A Celestial Intermediate
      Pole are elements (3,1-3) of the GCRS-to-true matrix,
      i.e. rbpn[2][0-2].

  10) It is permissible to re-use the same array in the returned
      arguments.  The arrays are filled in the order given.


------
pn00b
------

Given:
         TT as a 2-part Julian Date (Note 1)

         nutation (Note 2)
         mean obliquity (Note 3)
3][3]    frame bias matrix (Note 4)
3][3]    precession matrix (Note 5)
3][3]    bias-precession matrix (Note 6)
3][3]    nutation matrix (Note 7)
3][3]    GCRS-to-true matrix (Notes 8,9)

ate2 is a Julian Date, apportioned in any
een the two arguments.  For example,
uld be expressed in any of these ways,

   date2
     0.0       (JD method)
 -1421.3       (J2000 method)
 50123.2       (MJD method)
     0.2       (date & time method)
e most natural and convenient to use in
s of several decimal digits of resolution
 J2000 method is best matched to the way
dled internally and will deliver the
  The MJD method and the date & time methods
omises between resolution and convenience.
ents (luni-solar + planetary, IAU 2000B) in
uity are in radians and with respect to the
c of date.  For more accurate results, but
eased computation, use the eraPn00a function.
racy, use the eraPn00  function, where the
 are caller-specified.
is consistent with the IAU 2000 precession.
forms vectors from GCRS to J2000.0 mean
 by applying frame bias.
forms vectors from J2000.0 mean equator and
ator and equinox of date by applying

sforms vectors from GCRS to mean equator and
applying frame bias then precession.  It is
.
forms vectors from mean equator and equinox
ator and equinox of date by applying the
r + planetary).
nsforms vectors from GCRS to true equator and
t is the product rn x rbp, applying frame
d nutation in that order.
es of the IAU 2000B Celestial Intermediate
3,1-3) of the GCRS-to-true matrix,

o re-use the same array in the returned
ays are filled in the stated order.

n, IAU 2000B
ecession/nutation results, IAU 2000

ont, J., Lambert, S. and Wallace, P.,
 Celestial Intermediate Pole and Celestial
sistent with the IAU 2000A precession-
ron.Astrophys. 400, 1145-1154 (2003).
phemeris origin (CEO) was renamed "celestial
igin" (CIO) by IAU 2006 Resolution 2.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1)  The TT date date1+date2 is a Julian Date, apportioned in any
      convenient way between the two arguments.  For example,
      JD(TT)=2450123.7 could be expressed in any of these ways,
      among others:

             date1          date2

          2450123.7           0.0       (JD method)
          2451545.0       -1421.3       (J2000 method)
          2400000.5       50123.2       (MJD method)
          2450123.5           0.2       (date & time method)

      The JD method is the most natural and convenient to use in
      cases where the loss of several decimal digits of resolution
      is acceptable.  The J2000 method is best matched to the way
      the argument is handled internally and will deliver the
      optimum resolution.  The MJD method and the date & time methods
      are both good compromises between resolution and convenience.

  2)  The nutation components (luni-solar + planetary, IAU 2000B) in
      longitude and obliquity are in radians and with respect to the
      equinox and ecliptic of date.  For more accurate results, but
      at the cost of increased computation, use the eraPn00a function.
      For the utmost accuracy, use the eraPn00  function, where the
      nutation components are caller-specified.

  3)  The mean obliquity is consistent with the IAU 2000 precession.

  4)  The matrix rb transforms vectors from GCRS to J2000.0 mean
      equator and equinox by applying frame bias.

  5)  The matrix rp transforms vectors from J2000.0 mean equator and
      equinox to mean equator and equinox of date by applying
      precession.

  6)  The matrix rbp transforms vectors from GCRS to mean equator and
      equinox of date by applying frame bias then precession.  It is
      the product rp x rb.

  7)  The matrix rn transforms vectors from mean equator and equinox
      of date to true equator and equinox of date by applying the
      nutation (luni-solar + planetary).

  8)  The matrix rbpn transforms vectors from GCRS to true equator and
      equinox of date.  It is the product rn x rbp, applying frame
      bias, precession and nutation in that order.

  9)  The X,Y,Z coordinates of the IAU 2000B Celestial Intermediate
      Pole are elements (3,1-3) of the GCRS-to-true matrix,
      i.e. rbpn[2][0-2].

  10) It is permissible to re-use the same array in the returned
      arguments.  The arrays are filled in the stated order.


------
pn06
------

Given:
         TT as a 2-part Julian Date (Note 1)
         nutation (Note 2)

         mean obliquity (Note 3)
3][3]    frame bias matrix (Note 4)
3][3]    precession matrix (Note 5)
3][3]    bias-precession matrix (Note 6)
3][3]    nutation matrix (Note 7)
3][3]    GCRS-to-true matrix (Note 8)

ate2 is a Julian Date, apportioned in any
een the two arguments.  For example,
uld be expressed in any of these ways,

   date2
     0.0       (JD method)
 -1421.3       (J2000 method)
 50123.2       (MJD method)
     0.2       (date & time method)
e most natural and convenient to use in
s of several decimal digits of resolution
 J2000 method is best matched to the way
dled internally and will deliver the
  The MJD method and the date & time methods
omises between resolution and convenience.
nsible for providing the nutation components;
de and obliquity, in radians and are with
nox and ecliptic of date.  For high-accuracy
core nutation should be included as well as
corrections to the position of the CIP.
bliquity is consistent with the IAU 2006

forms vectors from GCRS to J2000.0 mean
 by applying frame bias.
forms vectors from J2000.0 mean equator and
ator and equinox of date by applying

sforms vectors from GCRS to mean equator and
applying frame bias then precession.  It is
.
forms vectors from mean equator and equinox
ator and equinox of date by applying the
r + planetary).
nsforms vectors from GCRS to true equator and
t is the product rn x rbp, applying frame
d nutation in that order.
es of the Celestial Intermediate Pole are
 the GCRS-to-true matrix, i.e. rbpn[2][0-2].
o re-use the same array in the returned
ays are filled in the stated order.

ecession F-W angles, IAU 2006
les to r-matrix
matrix
se r-matrix
 of two r-matrices

ace, P.T., 2006, Astron.Astrophys. 450, 855
taine, N., 2006, Astron.Astrophys. 459, 981
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1)  The TT date date1+date2 is a Julian Date, apportioned in any
      convenient way between the two arguments.  For example,
      JD(TT)=2450123.7 could be expressed in any of these ways,
      among others:

             date1          date2

          2450123.7           0.0       (JD method)
          2451545.0       -1421.3       (J2000 method)
          2400000.5       50123.2       (MJD method)
          2450123.5           0.2       (date & time method)

      The JD method is the most natural and convenient to use in
      cases where the loss of several decimal digits of resolution
      is acceptable.  The J2000 method is best matched to the way
      the argument is handled internally and will deliver the
      optimum resolution.  The MJD method and the date & time methods
      are both good compromises between resolution and convenience.

  2)  The caller is responsible for providing the nutation components;
      they are in longitude and obliquity, in radians and are with
      respect to the equinox and ecliptic of date.  For high-accuracy
      applications, free core nutation should be included as well as
      any other relevant corrections to the position of the CIP.

  3)  The returned mean obliquity is consistent with the IAU 2006
      precession.

  4)  The matrix rb transforms vectors from GCRS to J2000.0 mean
      equator and equinox by applying frame bias.

  5)  The matrix rp transforms vectors from J2000.0 mean equator and
      equinox to mean equator and equinox of date by applying
      precession.

  6)  The matrix rbp transforms vectors from GCRS to mean equator and
      equinox of date by applying frame bias then precession.  It is
      the product rp x rb.

  7)  The matrix rn transforms vectors from mean equator and equinox
      of date to true equator and equinox of date by applying the
      nutation (luni-solar + planetary).

  8)  The matrix rbpn transforms vectors from GCRS to true equator and
      equinox of date.  It is the product rn x rbp, applying frame
      bias, precession and nutation in that order.

  9)  The X,Y,Z coordinates of the Celestial Intermediate Pole are
      elements (3,1-3) of the GCRS-to-true matrix, i.e. rbpn[2][0-2].

  10) It is permissible to re-use the same array in the returned
      arguments.  The arrays are filled in the stated order.


References:

     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981

------
pn06a
------

Given:
         TT as a 2-part Julian Date (Note 1)

         nutation (Note 2)
         mean obliquity (Note 3)
3][3]    frame bias matrix (Note 4)
3][3]    precession matrix (Note 5)
3][3]    bias-precession matrix (Note 6)
3][3]    nutation matrix (Note 7)
3][3]    GCRS-to-true matrix (Notes 8,9)

ate2 is a Julian Date, apportioned in any
een the two arguments.  For example,
uld be expressed in any of these ways,

   date2
     0.0       (JD method)
 -1421.3       (J2000 method)
 50123.2       (MJD method)
     0.2       (date & time method)
e most natural and convenient to use in
s of several decimal digits of resolution
 J2000 method is best matched to the way
dled internally and will deliver the
  The MJD method and the date & time methods
omises between resolution and convenience.
ents (luni-solar + planetary, IAU 2000A) in
uity are in radians and with respect to the
c of date.  Free core nutation is omitted;
racy, use the eraPn06 function, where the
 are caller-specified.
is consistent with the IAU 2006 precession.
forms vectors from GCRS to mean J2000.0 by
.
forms vectors from mean J2000.0 to mean of
ecession.
sforms vectors from GCRS to mean of date by
 then precession.  It is the product rp x rb.
forms vectors from mean of date to true of
e nutation (luni-solar + planetary).
nsforms vectors from GCRS to true of date
is the product rn x rbp, applying frame bias,
tion in that order.
es of the IAU 2006/2000A Celestial
re elements (3,1-3) of the GCRS-to-true
][0-2].
o re-use the same array in the returned
ays are filled in the stated order.

n, IAU 2006/2000A
ecession/nutation results, IAU 2006

ace, P.T., 2006, Astron.Astrophys. 450, 855
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1)  The TT date date1+date2 is a Julian Date, apportioned in any
      convenient way between the two arguments.  For example,
      JD(TT)=2450123.7 could be expressed in any of these ways,
      among others:

             date1          date2

          2450123.7           0.0       (JD method)
          2451545.0       -1421.3       (J2000 method)
          2400000.5       50123.2       (MJD method)
          2450123.5           0.2       (date & time method)

      The JD method is the most natural and convenient to use in
      cases where the loss of several decimal digits of resolution
      is acceptable.  The J2000 method is best matched to the way
      the argument is handled internally and will deliver the
      optimum resolution.  The MJD method and the date & time methods
      are both good compromises between resolution and convenience.

  2)  The nutation components (luni-solar + planetary, IAU 2000A) in
      longitude and obliquity are in radians and with respect to the
      equinox and ecliptic of date.  Free core nutation is omitted;
      for the utmost accuracy, use the eraPn06 function, where the
      nutation components are caller-specified.

  3)  The mean obliquity is consistent with the IAU 2006 precession.

  4)  The matrix rb transforms vectors from GCRS to mean J2000.0 by
      applying frame bias.

  5)  The matrix rp transforms vectors from mean J2000.0 to mean of
      date by applying precession.

  6)  The matrix rbp transforms vectors from GCRS to mean of date by
      applying frame bias then precession.  It is the product rp x rb.

  7)  The matrix rn transforms vectors from mean of date to true of
      date by applying the nutation (luni-solar + planetary).

  8)  The matrix rbpn transforms vectors from GCRS to true of date
      (CIP/equinox).  It is the product rn x rbp, applying frame bias,
      precession and nutation in that order.

  9)  The X,Y,Z coordinates of the IAU 2006/2000A Celestial
      Intermediate Pole are elements (3,1-3) of the GCRS-to-true
      matrix, i.e. rbpn[2][0-2].

  10) It is permissible to re-use the same array in the returned
      arguments.  The arrays are filled in the stated order.


------
pnm00a
------

Given:
    TT as a 2-part Julian Date (Note 1)

3][3]    classical NPB matrix (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
in the sense V(date) = rbpn * V(GCRS), where
 is with respect to the true equatorial triad
and the p-vector V(GCRS) is with respect to
tial Reference System (IAU, 2000).
ly less accurate result (about 1 mas), can be
stead the eraPnm00b function.

ecession/nutation, IAU 2000A

ional Astronomical Union, Vol. XXIVB;  Proc.
y, Manchester, UK.  Resolutions B1.3, B1.6.

, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix operates in the sense V(date) = rbpn * V(GCRS), where
     the p-vector V(date) is with respect to the true equatorial triad
     of date date1+date2 and the p-vector V(GCRS) is with respect to
     the Geocentric Celestial Reference System (IAU, 2000).

  3) A faster, but slightly less accurate result (about 1 mas), can be
     obtained by using instead the eraPnm00b function.


------
pnm00b
------

Given:
     TT as a 2-part Julian Date (Note 1)

][3] bias-precession-nutation matrix (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
in the sense V(date) = rbpn * V(GCRS), where
 is with respect to the true equatorial triad
and the p-vector V(GCRS) is with respect to
tial Reference System (IAU, 2000).
 is faster, but slightly less accurate (about
Pnm00a function.

ecession/nutation, IAU 2000B

ional Astronomical Union, Vol. XXIVB;  Proc.
y, Manchester, UK.  Resolutions B1.3, B1.6.

, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix operates in the sense V(date) = rbpn * V(GCRS), where
     the p-vector V(date) is with respect to the true equatorial triad
     of date date1+date2 and the p-vector V(GCRS) is with respect to
     the Geocentric Celestial Reference System (IAU, 2000).

  3) The present function is faster, but slightly less accurate (about
     1 mas), than the eraPnm00a function.


------
pnm06a
------

Given:
     TT as a 2-part Julian Date (Note 1)

][3] bias-precession-nutation matrix (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
in the sense V(date) = rnpb * V(GCRS), where
 is with respect to the true equatorial triad
and the p-vector V(GCRS) is with respect to
tial Reference System (IAU, 2000).

ecession F-W angles, IAU 2006
n, IAU 2006/2000A
les to r-matrix

ace, P.T., 2006, Astron.Astrophys. 450, 855.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix operates in the sense V(date) = rnpb * V(GCRS), where
     the p-vector V(date) is with respect to the true equatorial triad
     of date date1+date2 and the p-vector V(GCRS) is with respect to
     the Geocentric Celestial Reference System (IAU, 2000).


------
pnm80
------

Given:
e         TDB date (Note 1)

e[3][3]   combined precession/nutation matrix

ate2 is a Julian Date, apportioned in any
en the two arguments.  For example,
uld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
in the sense V(date) = rmatpn * V(J2000),
(date) is with respect to the true equatorial
date2 and the p-vector V(J2000) is with
equatorial triad of epoch J2000.0.

ion matrix, IAU 1976
n matrix, IAU 1980
 of two r-matrices

nt to the Astronomical Almanac,
n (ed), University Science Books (1992),

, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TDB date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TDB)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The matrix operates in the sense V(date) = rmatpn * V(J2000),
     where the p-vector V(date) is with respect to the true equatorial
     triad of date date1+date2 and the p-vector V(J2000) is with
     respect to the mean equatorial triad of epoch J2000.0.


------
pom00
------

Given:
oordinates of the pole (radians, Note 1)
he TIO locator s' (radians, Note 2)

]   polar-motion matrix (Note 3)

 yp are the coordinates (in radians) of the
te Pole with respect to the International
e System (see IERS Conventions 2003),
eridians to 0 and 90 deg west respectively.
he TIO locator s', in radians, which
trial Intermediate Origin on the equator.  It
ar motion observations by numerical
is in essence unpredictable.  However, it is
ar drift of about 47 microarcseconds per
be taken into account by using s' = -47*t,
 since J2000.0.  The function eraSp00
oximation.
in the sense V(TRS) = rpom * V(CIP), meaning
 rotation when computing the pointing
tial source.

ize r-matrix to identity
around Z-axis
around Y-axis
around X-axis

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The arguments xp and yp are the coordinates (in radians) of the
     Celestial Intermediate Pole with respect to the International
     Terrestrial Reference System (see IERS Conventions 2003),
     measured along the meridians to 0 and 90 deg west respectively.

  2) The argument sp is the TIO locator s', in radians, which
     positions the Terrestrial Intermediate Origin on the equator.  It
     is obtained from polar motion observations by numerical
     integration, and so is in essence unpredictable.  However, it is
     dominated by a secular drift of about 47 microarcseconds per
     century, and so can be taken into account by using s' = -47*t,
     where t is centuries since J2000.0.  The function eraSp00
     implements this approximation.

  3) The matrix operates in the sense V(TRS) = rpom * V(CIP), meaning
     that it is the final rotation when computing the pointing
     direction to a celestial source.


------
ppp
------

Given:
    first p-vector
    second p-vector

    a + b

 re-use the same array for any of the

, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
ppsp
------

Given:
 first p-vector
 scalar (multiplier for b)
 second p-vector

 a + s*b

r any of a, b and apsb to be the same array.

y p-vector by scalar
r plus p-vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
pr00
------

Given:
e  TT as a 2-part Julian Date (Note 1)

e  precession corrections (Notes 2,3)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
tments are expressed as "nutation
ions in longitude and obliquity with respect
ox and ecliptic.
ion adjustments are stated to be with respect
77), the MHB2000 model does not specify which
are to be used and how the adjustments are to
t literal and straightforward procedure is to
 epsilon_0, psi_A, omega_A, xi_A option, and
_A and depspr to both omega_A and eps_A.
ation of one aspect of the IAU 2000A nutation
ted by the IAU General Assembly in 2000,
ews et al. 2002).

e, T., Fricke, W. & Morando, B., "Expressions
uantities based upon the IAU (1976) System of
ts", Astron.Astrophys., 58, 1-16 (1977)
ng, T.A., Buffet, B.A., "Modeling of nutation
 nutation series for nonrigid Earth and
rth's interior", J.Geophys.Res., 107, B4,
ode itself was obtained on 9th September 2002
.navy.mil/conv2000/chapter5/IAU2000A.
ware for Implementing the IAU 2000
S Workshop 5.1 (2002).
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The precession adjustments are expressed as "nutation
     components", corrections in longitude and obliquity with respect
     to the J2000.0 equinox and ecliptic.

  3) Although the precession adjustments are stated to be with respect
     to Lieske et al. (1977), the MHB2000 model does not specify which
     set of Euler angles are to be used and how the adjustments are to
     be applied.  The most literal and straightforward procedure is to
     adopt the 4-rotation epsilon_0, psi_A, omega_A, xi_A option, and
     to add dpsipr to psi_A and depspr to both omega_A and eps_A.

  4) This is an implementation of one aspect of the IAU 2000A nutation
     model, formally adopted by the IAU General Assembly in 2000,
     namely MHB2000 (Mathews et al. 2002).


References:

     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B., "Expressions
     for the precession quantities based upon the IAU (1976) System of
     Astronomical Constants", Astron.Astrophys., 58, 1-16 (1977)

     Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation
     and precession   New nutation series for nonrigid Earth and
     insights into the Earth's interior", J.Geophys.Res., 107, B4,
     2002.  The MHB2000 code itself was obtained on 9th September 2002
     from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.

     Wallace, P.T., "Software for Implementing the IAU 2000
     Resolutions", in IERS Workshop 5.1 (2002).

------
prec76
------

Given:
le    TDB starting date (Note 1)
le    TDB ending date (Note 1)

le    1st rotation: radians cw around z
le    3rd rotation: radians cw around z
le    2nd rotation: radians ccw around y

e02 and date11+date12 are Julian Dates,
onvenient way between the arguments daten1
mple, JD(TDB)=2450123.7 could be expressed in
mong others:
daten2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in cases
veral decimal digits of resolution is
00 method is best matched to the way the
internally and will deliver the optimum
 The MJD method and the date & time methods
mises between resolution and convenience.
 expressed using different methods, but at
ome resolution.
ession angles zeta, z, theta are expressed
lynomials which are valid only for a limited
ion, the IAU 1976 precession rate is known to
bsolute accuracy of the present formulation
rcsec from 1960AD to 2040AD, better than
 to 2360AD, and remains below 3 arcsec for
iod 500BC to 3000AD.  The errors exceed
e range 1200BC to 3900AD, exceed 100 arcsec
00AD and exceed 1000 arcsec outside 6800BC to

 returned in the conventional order, which
he order of the corresponding Euler
ession matrix is
) x R_3(-zeta).

Astron.Astrophys. 73, 282, equations

, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
pv2p
------

Given:
     pv-vector

     p-vector

vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
pv2s
------

Given:
]  pv-vector

   longitude angle (radians)
   latitude angle (radians)
   radial distance
   rate of change of theta
   rate of change of phi
   rate of change of r

 of pv is null, theta, phi, td and pd
This is handled by extrapolating the
t time by using the velocity part of
origin without changing the direction
onent.  If the position and velocity
 both null, zeroes are returned for all

 pole, theta, td and pd are indeterminate.
 are returned for all three.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
pvdpv
------

Given:
]      first pv-vector
]      second pv-vector

       a . b (see note)

velocity components of the two pv-vectors are
 bv ), the result, a . b, is the pair of
ap . bv + av . bp ).  The two numbers are the
wo p-vectors and its derivative.

product of two p-vectors
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
pvm
------

Given:
  pv-vector

  modulus of position component
  modulus of velocity component

 of p-vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
pvmpv
------

Given:
      first pv-vector
      second pv-vector

      a - b

 re-use the same array for any of the


r minus p-vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
pvppv
------

Given:
]      first pv-vector
]      second pv-vector

]      a + b

 re-use the same array for any of the


r plus p-vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
pvstar
------


Notes:

  1) The specified pv-vector is the coordinate direction (and its rate
     of change) for the date at which the light leaving the star
     reached the solar-system barycenter.

  2) The star data returned by this function are "observables" for an
     imaginary observer at the solar-system barycenter.  Proper motion
     and radial velocity are, strictly, in terms of barycentric
     coordinate time, TCB.  For most practical applications, it is
     permissible to neglect the distinction between TCB and ordinary
     "proper" time on Earth (TT/TAI).  The result will, as a rule, be
     limited by the intrinsic accuracy of the proper-motion and
     radial-velocity data;  moreover, the supplied pv-vector is likely
     to be merely an intermediate result (for example generated by the
     function eraStarpv), so that a change of time unit will cancel
     out overall.

     In accordance with normal star-catalog conventions, the object's
     right ascension and declination are freed from the effects of
     secular aberration.  The frame, which is aligned to the catalog
     equator and equinox, is Lorentzian and centered on the SSB.

     Summarizing, the specified pv-vector is for most stars almost
     identical to the result of applying the standard geometrical
     "space motion" transformation to the catalog data.  The
     differences, which are the subject of the Stumpff paper cited
     below, are:

     (i) In stars with significant radial velocity and proper motion,
     the constantly changing light-time distorts the apparent proper
     motion.  Note that this is a classical, not a relativistic,
     effect.

     (ii) The transformation complies with special relativity.

  3) Care is needed with units.  The star coordinates are in radians
     and the proper motions in radians per Julian year, but the
     parallax is in arcseconds; the radial velocity is in km/s, but
     the pv-vector result is in AU and AU/day.

  4) The proper motions are the rate of change of the right ascension
     and declination at the catalog epoch and are in radians per Julian
     year.  The RA proper motion is in terms of coordinate angle, not
     true angle, and will thus be numerically larger at high
     declinations.

  5) Straight-line motion at constant speed in the inertial frame is
     assumed.  If the speed is greater than or equal to the speed of
     light, the function aborts with an error status.

  6) The inverse transformation is performed by the function eraStarpv.


------
pvtob
------

Given:
 longitude (radians, east +ve, Note 1)
 latitude (geodetic, radians, Note 1)
 height above ref. ellipsoid (geodetic, m)
 coordinates of the pole (radians, Note 2)
 the TIO locator s' (radians, Note 2)
 Earth rotation angle (radians, Note 3)

 position/velocity vector (m, m/s, CIRS)

dinates are with respect to the ERFA_WGS84

ordinates (in radians) of the Celestial
th respect to the International Terrestrial
e IERS Conventions), measured along the
eg west respectively.  sp is the TIO locator
h positions the Terrestrial Intermediate
r.  For many applications, xp, yp and
be set to zero.
h apparent sidereal time instead of Earth
result is with respect to the true equator
 i.e. with the x-axis at the equinox rather
ntermediate origin.
re meters per UT1 second, not per SI second.
have any practical consequences in the modern

formed on the arguments.  Error cases that
etic exceptions are trapped by the eraGd2gc
sult set to zeros.

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
nn, P. K. (eds), Explanatory Supplement to
anac, 3rd ed., University Science Books
3.3.

c to geocentric transformation
otion matrix
 of transpose of r-matrix and p-vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The terrestrial coordinates are with respect to the ERFA_WGS84
     reference ellipsoid.

  2) xp and yp are the coordinates (in radians) of the Celestial
     Intermediate Pole with respect to the International Terrestrial
     Reference System (see IERS Conventions), measured along the
     meridians 0 and 90 deg west respectively.  sp is the TIO locator
     s', in radians, which positions the Terrestrial Intermediate
     Origin on the equator.  For many applications, xp, yp and
     (especially) sp can be set to zero.

  3) If theta is Greenwich apparent sidereal time instead of Earth
     rotation angle, the result is with respect to the true equator
     and equinox of date, i.e. with the x-axis at the equinox rather
     than the celestial intermediate origin.

  4) The velocity units are meters per UT1 second, not per SI second.
     This is unlikely to have any practical consequences in the modern
     era.

  5) No validation is performed on the arguments.  Error cases that
     could lead to arithmetic exceptions are trapped by the eraGd2gc
     function, and the result set to zeros.

  References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
     the Astronomical Almanac, 3rd ed., University Science Books
     (2013), Section 7.4.3.3.


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
     the Astronomical Almanac, 3rd ed., University Science Books
     (2013), Section 7.4.3.3.

  Called:
     eraGd2gc     geodetic to geocentric transformation
     eraPom00     polar motion matrix
     eraTrxp      product of transpose of r-matrix and p-vector

------
pvu
------

Given:
      time interval
]     pv-vector

]     p updated, v unchanged

r the position component of the vector
e units from the existing date".
 must match those of the velocity.
r pv and upv to be the same array.

r plus scaled p-vector
vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) "Update" means "refer the position component of the vector
     to a new date dt time units from the existing date".

  2) The time units of dt must match those of the velocity.

  3) It is permissible for pv and upv to be the same array.


------
pvup
------

Given:
       time interval
]      pv-vector

       p-vector

r the position component of the vector to a
ts from the existing date".
 must match those of the velocity.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
pvxpv
------

Given:
]      first pv-vector
]      second pv-vector

]      a x b

velocity components of the two pv-vectors are
 bv ), the result, a x b, is the pair of
p x bv + av x bp ).  The two vectors are the
 two p-vectors and its derivative.
 re-use the same array for any of the


-vector
product of two p-vectors
r plus p-vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) If the position and velocity components of the two pv-vectors are
     ( ap, av ) and ( bp, bv ), the result, a x b, is the pair of
     vectors ( ap x bp, ap x bv + av x bp ).  The two vectors are the
     cross-product of the two p-vectors and its derivative.

  2) It is permissible to re-use the same array for any of the
     arguments.


------
pxp
------

Given:
    first p-vector
    second p-vector

    a x b

 re-use the same array for any of the

, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
refco
------

Given:
sure at the observer (hPa = millibar)
ent temperature at the observer (deg C)
tive humidity at the observer (range 0-1)
length (micrometers)

Z coefficient (radians)
3 Z coefficient (radians)

peed and accuracy to give good results in
erformance at low altitudes is not paramount.
ained across a range of conditions, and
cal/IR and radio.
effects of (i) height above sea level (apart
ssure itself), (ii) latitude (i.e. the
rth), (iii) variations in tropospheric lapse
sive effects in the radio.
 using the following range of conditions:
, 0.0065, 0.0075 deg/meter
0, 75 degrees
000 meters ASL
 height -10% to +5% in steps of 5%
eg to +20 deg with respect to 280 deg at SL
0, 0.5, 1
.6, ... 2 micron, + radio
5, 45, 75 degrees
spect to raytracing through a model
llows:
   worst         RMS
   62 mas       8 mas
  319 mas      49 mas
set of conditions:
K/meter
s


 K

gstroms
follows:
    eraRefco   Saastamoinen
      10.27        10.27
      21.20        21.19
      33.61        33.60
      48.83        48.81
      58.18        58.16
      69.30        69.27
      82.99        82.95
     100.54       100.50
     124.26       124.20
     158.68       158.61
     177.37       177.31
     200.38       200.32
     229.43       229.42
     267.29       267.41
     318.55       319.10
     arcsec       arcsec
amoinen's formula (which includes terms
en from Hohenkerk and Sinclair (1985).
nge 0-100 selects the optical/IR case and is
eters.  Any value outside this range selects

ameters are silently limited to
values.  Zero pressure is permissible, and
returned.
on several sources, as follows:
he saturation vapour pressure of water as
perature and temperature is taken from
4.7) of Gill (1982).
he water vapour pressure, given the
re and the relative humidity, is from
ation (2.5.5).
of air is a function of temperature,
ater-vapour pressure and, in the case
velength.  The formulae for the two cases are
henkerk & Sinclair (1985) and Rueger (2002).
eta, the ratio of the scale height of the
 geocentric distance of the observer, is
uation (9) from Stone (1996).  The
ved at empirically, consist of (i) a small
 coefficient and (ii) a humidity term for the

the refraction constants as a function of
from Green (1987), Equation (4.31).

M.L. (ed), "Refraction Effects in the Neutral
 of Experimental Physics: Astrophysics 12B,
.
mosphere-Ocean Dynamics", Academic Press,

cal Astronomy", Cambridge University Press,

inclair, A.T., NAO Technical Note No. 63,

ctive Index Formulae for Electronic Distance
io and Millimetre Waves", in Unisurv Report
eying and Spatial Information Systems,
uth Wales, Sydney, Australia, 2002.
A.S.P. 108, 1051-1058, 1996.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The model balances speed and accuracy to give good results in
     applications where performance at low altitudes is not paramount.
     Performance is maintained across a range of conditions, and
     applies to both optical/IR and radio.

  2) The model omits the effects of (i) height above sea level (apart
     from the reduced pressure itself), (ii) latitude (i.e. the
     flattening of the Earth), (iii) variations in tropospheric lapse
     rate and (iv) dispersive effects in the radio.

     The model was tested using the following range of conditions:

       lapse rates 0.0055, 0.0065, 0.0075 deg/meter
       latitudes 0, 25, 50, 75 degrees
       heights 0, 2500, 5000 meters ASL
       pressures mean for height -10% to +5% in steps of 5%
       temperatures -10 deg to +20 deg with respect to 280 deg at SL
       relative humidity 0, 0.5, 1
       wavelengths 0.4, 0.6, ... 2 micron, + radio
       zenith distances 15, 45, 75 degrees

     The accuracy with respect to raytracing through a model
     atmosphere was as follows:

                            worst         RMS

       optical/IR           62 mas       8 mas
       radio               319 mas      49 mas

     For this particular set of conditions:

       lapse rate 0.0065 K/meter
       latitude 50 degrees
       sea level
       pressure 1005 mb
       temperature 280.15 K
       humidity 80%
       wavelength 5740 Angstroms

     the results were as follows:

       ZD       raytrace     eraRefco   Saastamoinen

       10         10.27        10.27        10.27
       20         21.19        21.20        21.19
       30         33.61        33.61        33.60
       40         48.82        48.83        48.81
       45         58.16        58.18        58.16
       50         69.28        69.30        69.27
       55         82.97        82.99        82.95
       60        100.51       100.54       100.50
       65        124.23       124.26       124.20
       70        158.63       158.68       158.61
       72        177.32       177.37       177.31
       74        200.35       200.38       200.32
       76        229.45       229.43       229.42
       78        267.44       267.29       267.41
       80        319.13       318.55       319.10

      deg        arcsec       arcsec       arcsec

     The values for Saastamoinen's formula (which includes terms
     up to tan^5) are taken from Hohenkerk and Sinclair (1985).

  3) A wl value in the range 0-100 selects the optical/IR case and is
     wavelength in micrometers.  Any value outside this range selects
     the radio case.

  4) Outlandish input parameters are silently limited to
     mathematically safe values.  Zero pressure is permissible, and
     causes zeroes to be returned.

  5) The algorithm draws on several sources, as follows:

     a) The formula for the saturation vapour pressure of water as
        a function of temperature and temperature is taken from
        Equations (A4.5-A4.7) of Gill (1982).

     b) The formula for the water vapour pressure, given the
        saturation pressure and the relative humidity, is from
        Crane (1976), Equation (2.5.5).

     c) The refractivity of air is a function of temperature,
        total pressure, water-vapour pressure and, in the case
        of optical/IR, wavelength.  The formulae for the two cases are
        developed from Hohenkerk & Sinclair (1985) and Rueger (2002).

     d) The formula for beta, the ratio of the scale height of the
        atmosphere to the geocentric distance of the observer, is
        an adaption of Equation (9) from Stone (1996).  The
        adaptations, arrived at empirically, consist of (i) a small
        adjustment to the coefficient and (ii) a humidity term for the
        radio case only.

     e) The formulae for the refraction constants as a function of
        n-1 and beta are from Green (1987), Equation (4.31).


References:

     Crane, R.K., Meeks, M.L. (ed), "Refraction Effects in the Neutral
     Atmosphere", Methods of Experimental Physics: Astrophysics 12B,
     Academic Press, 1976.

     Gill, Adrian E., "Atmosphere-Ocean Dynamics", Academic Press,
     1982.

     Green, R.M., "Spherical Astronomy", Cambridge University Press,
     1987.

     Hohenkerk, C.Y., & Sinclair, A.T., NAO Technical Note No. 63,
     1985.

     Rueger, J.M., "Refractive Index Formulae for Electronic Distance
     Measurement with Radio and Millimetre Waves", in Unisurv Report
     S-68, School of Surveying and Spatial Information Systems,
     University of New South Wales, Sydney, Australia, 2002.

     Stone, Ronald C., P.A.S.P. 108, 1051-1058, 1996.

------
rm2v
------

Given:
]    rotation matrix

     rotation vector (Note 1)

scribes a rotation through some angle about
called the Euler axis.  The "rotation vector"
ction has the same direction as the Euler axis,
 the angle in radians.  (The magnitude and
arated by means of the function eraPn.)
the result.  If r is not a rotation matrix
ned;  r must be proper (i.e. have a positive
l orthogonal (inverse = transpose).
rotates clockwise as seen looking along
from the origin.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
rv2m
------

Given:
    rotation vector (Note 1)

]    rotation matrix

scribes a rotation through some angle about
called the Euler axis.  The "rotation vector"
ction has the same direction as the Euler
ude is the angle in radians.
it matrix is returned.
rotates clockwise as seen looking along the
 the origin.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
rx
------

Given:
   angle (radians)

   r-matrix, rotated

n with positive phi incorporates in the
an additional rotation, about the x-axis,
n looking towards the origin from positive x.
ion can be represented by this matrix:
         0      )
                )
)   + sin(phi)  )
                )
)   + cos(phi)  )
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
rxp
------

Given:
]    r-matrix
     p-vector

     r * p

r p and rp to be the same array.

vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
rxpv
------

Given:
]    r-matrix
]    pv-vector

]    r * pv

r pv and rpv to be the same array.

 of r-matrix and p-vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
rxr
------

Given:
]    first r-matrix
]    second r-matrix

]    a * b

 re-use the same array for any of the


matrix
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
ry
------

Given:
   angle (radians)

   r-matrix, rotated

n with positive theta incorporates in the
an additional rotation, about the y-axis,
n looking towards the origin from positive y.
ion can be represented by this matrix:
    0      - sin(theta)  )
                         )
    1           0        )
                         )
    0      + cos(theta)  )
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
rz
------

Given:
   angle (radians)

   r-matrix, rotated

n with positive psi incorporates in the
an additional rotation, about the z-axis,
n looking towards the origin from positive z.
ion can be represented by this matrix:
+ sin(psi)     0  )
                  )
+ cos(psi)     0  )
                  )
     0         1  )
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
s00
------

Given:
    TT as a 2-part Julian Date (Note 1)
    CIP coordinates (Note 3)
e):
    the CIO locator s in radians (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
 the difference between the right ascensions
 two systems:  the two systems are the GCRS
 the point is the ascending node of the
antity s remains below 0.1 arcsecond
.
ompute s is in fact for s+XY/2, where X and Y
onents of the CIP unit vector;  this series
 a direct series for s would be.  This
Y to be supplied by the caller, who is
iding values that are consistent with the

ent with the IAU 2000A precession-nutation.

omaly of the Moon
omaly of the Sun
gument of the latitude of the Moon
ongation of the Moon from the Sun
ngitude of the Moon's ascending node
ngitude of Venus
ngitude of Earth
 accumulated precession in longitude

ont, J., Lambert, S. and Wallace, P.,
 Celestial Intermediate Pole and Celestial
sistent with the IAU 2000A precession-
ron.Astrophys. 400, 1145-1154 (2003)
phemeris origin (CEO) was renamed "celestial
igin" (CIO) by IAU 2006 Resolution 2.
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The CIO locator s is the difference between the right ascensions
     of the same point in two systems:  the two systems are the GCRS
     and the CIP,CIO, and the point is the ascending node of the
     CIP equator.  The quantity s remains below 0.1 arcsecond
     throughout 1900-2100.

  3) The series used to compute s is in fact for s+XY/2, where X and Y
     are the x and y components of the CIP unit vector;  this series
     is more compact than a direct series for s would be.  This
     function requires X,Y to be supplied by the caller, who is
     responsible for providing values that are consistent with the
     supplied date.

  4) The model is consistent with the IAU 2000A precession-nutation.


References:

     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
     "Expressions for the Celestial Intermediate Pole and Celestial
     Ephemeris Origin consistent with the IAU 2000A precession-
     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
          intermediate origin" (CIO) by IAU 2006 Resolution 2.

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

------
s00a
------

Given:
   TT as a 2-part Julian Date (Note 1)
e):
   the CIO locator s in radians (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
 the difference between the right ascensions
 two systems.  The two systems are the GCRS
 the point is the ascending node of the
O locator s remains a small fraction of
ut 1900-2100.
ompute s is in fact for s+XY/2, where X and Y
onents of the CIP unit vector;  this series
 a direct series for s would be.  The present
ll IAU 2000A nutation model when predicting
aster results, with no significant loss of
ained via the function eraS00b, which uses
B truncated model.

al NPB matrix, IAU 2000A
 CIP X,Y from the BPN matrix
 locator s, given X,Y, IAU 2000A

ont, J., Lambert, S. and Wallace, P.,
 Celestial Intermediate Pole and Celestial
sistent with the IAU 2000A precession-
ron.Astrophys. 400, 1145-1154 (2003)
phemeris origin (CEO) was renamed "celestial
igin" (CIO) by IAU 2006 Resolution 2.
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The CIO locator s is the difference between the right ascensions
     of the same point in two systems.  The two systems are the GCRS
     and the CIP,CIO, and the point is the ascending node of the
     CIP equator.  The CIO locator s remains a small fraction of
     1 arcsecond throughout 1900-2100.

  3) The series used to compute s is in fact for s+XY/2, where X and Y
     are the x and y components of the CIP unit vector;  this series
     is more compact than a direct series for s would be.  The present
     function uses the full IAU 2000A nutation model when predicting
     the CIP position.  Faster results, with no significant loss of
     accuracy, can be obtained via the function eraS00b, which uses
     instead the IAU 2000B truncated model.


References:

     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
     "Expressions for the Celestial Intermediate Pole and Celestial
     Ephemeris Origin consistent with the IAU 2000A precession-
     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
          intermediate origin" (CIO) by IAU 2006 Resolution 2.

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

------
s00b
------

Given:
   TT as a 2-part Julian Date (Note 1)
e):
   the CIO locator s in radians (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
 the difference between the right ascensions
 two systems.  The two systems are the GCRS
 the point is the ascending node of the
O locator s remains a small fraction of
ut 1900-2100.
ompute s is in fact for s+XY/2, where X and Y
onents of the CIP unit vector;  this series
 a direct series for s would be.  The present
U 2000B truncated nutation model when
osition.  The function eraS00a uses instead
odel, but with no significant increase in
 cost in speed.

al NPB matrix, IAU 2000B
 CIP X,Y from the BPN matrix
 locator s, given X,Y, IAU 2000A

ont, J., Lambert, S. and Wallace, P.,
 Celestial Intermediate Pole and Celestial
sistent with the IAU 2000A precession-
ron.Astrophys. 400, 1145-1154 (2003)
phemeris origin (CEO) was renamed "celestial
igin" (CIO) by IAU 2006 Resolution 2.
it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The CIO locator s is the difference between the right ascensions
     of the same point in two systems.  The two systems are the GCRS
     and the CIP,CIO, and the point is the ascending node of the
     CIP equator.  The CIO locator s remains a small fraction of
     1 arcsecond throughout 1900-2100.

  3) The series used to compute s is in fact for s+XY/2, where X and Y
     are the x and y components of the CIP unit vector;  this series
     is more compact than a direct series for s would be.  The present
     function uses the IAU 2000B truncated nutation model when
     predicting the CIP position.  The function eraS00a uses instead
     the full IAU 2000A model, but with no significant increase in
     accuracy and at some cost in speed.


References:

     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
     "Expressions for the Celestial Intermediate Pole and Celestial
     Ephemeris Origin consistent with the IAU 2000A precession-
     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
          intermediate origin" (CIO) by IAU 2006 Resolution 2.

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

------
s06
------

Given:
    TT as a 2-part Julian Date (Note 1)
    CIP coordinates (Note 3)
e):
    the CIO locator s in radians (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
 the difference between the right ascensions
 two systems:  the two systems are the GCRS
 the point is the ascending node of the
antity s remains below 0.1 arcsecond
.
ompute s is in fact for s+XY/2, where X and Y
onents of the CIP unit vector;  this series
 a direct series for s would be.  This
Y to be supplied by the caller, who is
iding values that are consistent with the

ent with the "P03" precession (Capitaine et
y IAU 2006 Resolution 1, 2006, and the
with P03 adjustments).

omaly of the Moon
omaly of the Sun
gument of the latitude of the Moon
ongation of the Moon from the Sun
ngitude of the Moon's ascending node
ngitude of Venus
ngitude of Earth
 accumulated precession in longitude

ce, P.T. & Chapront, J., 2003, Astron.

t, G. (eds.) 2004, IERS Conventions (2003),
No. 32, BKG
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The CIO locator s is the difference between the right ascensions
     of the same point in two systems:  the two systems are the GCRS
     and the CIP,CIO, and the point is the ascending node of the
     CIP equator.  The quantity s remains below 0.1 arcsecond
     throughout 1900-2100.

  3) The series used to compute s is in fact for s+XY/2, where X and Y
     are the x and y components of the CIP unit vector;  this series
     is more compact than a direct series for s would be.  This
     function requires X,Y to be supplied by the caller, who is
     responsible for providing values that are consistent with the
     supplied date.

  4) The model is consistent with the "P03" precession (Capitaine et
     al. 2003), adopted by IAU 2006 Resolution 1, 2006, and the
     IAU 2000A nutation (with P03 adjustments).


References:

     Capitaine, N., Wallace, P.T. & Chapront, J., 2003, Astron.
     Astrophys. 432, 355

     McCarthy, D.D., Petit, G. (eds.) 2004, IERS Conventions (2003),
     IERS Technical Note No. 32, BKG

------
s06a
------

Given:
   TT as a 2-part Julian Date (Note 1)
e):
   the CIO locator s in radians (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
 the difference between the right ascensions
 two systems.  The two systems are the GCRS
 the point is the ascending node of the
O locator s remains a small fraction of
ut 1900-2100.
ompute s is in fact for s+XY/2, where X and Y
onents of the CIP unit vector;  this series is
direct series for s would be.  The present
ll IAU 2000A nutation model when predicting


al NPB matrix, IAU 2006/2000A
 CIP X,Y coordinates from NPB matrix
 locator s, given X,Y, IAU 2006

ont, J., Lambert, S. and Wallace, P.,
 Celestial Intermediate Pole and Celestial
sistent with the IAU 2000A precession-
ron.Astrophys. 400, 1145-1154 (2003)
phemeris origin (CEO) was renamed "celestial
igin" (CIO) by IAU 2006 Resolution 2.
ace, P.T., 2006, Astron.Astrophys. 450, 855
it, G. (eds.), 2004, IERS Conventions (2003),
No. 32, BKG
taine, N., 2006, Astron.Astrophys. 459, 981
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The CIO locator s is the difference between the right ascensions
     of the same point in two systems.  The two systems are the GCRS
     and the CIP,CIO, and the point is the ascending node of the
     CIP equator.  The CIO locator s remains a small fraction of
     1 arcsecond throughout 1900-2100.

  3) The series used to compute s is in fact for s+XY/2, where X and Y
     are the x and y components of the CIP unit vector;  this series is
     more compact than a direct series for s would be.  The present
     function uses the full IAU 2000A nutation model when predicting
     the CIP position.


References:

     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
     "Expressions for the Celestial Intermediate Pole and Celestial
     Ephemeris Origin consistent with the IAU 2000A precession-
     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
          intermediate origin" (CIO) by IAU 2006 Resolution 2.

     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
     IERS Technical Note No. 32, BKG

     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981

------
s2c
------

Given:
  longitude angle (radians)
  latitude angle (radians)

  direction cosines
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
s2p
------

Given:
 longitude angle (radians)
 latitude angle (radians)
 radial distance

 Cartesian coordinates

al coordinates to unit vector
y p-vector by scalar
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
s2pv
------

Given:
     longitude angle (radians)
     latitude angle (radians)
     radial distance
     rate of change of theta
     rate of change of phi
     rate of change of r

]    pv-vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
s2xpv
------

Given:
  scalar to multiply position component by
  scalar to multiply velocity component by
  pv-vector

  pv-vector: p scaled by s1, v scaled by s2

r pv and spv to be the same array.

y p-vector by scalar
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
sepp
------

Given:
first p-vector (not necessarily unit length)
second p-vector (not necessarily unit length)
e):
angular separation (radians, always positive)

null, a zero result is returned.
on is most simply formulated in terms of
ever, this gives poor accuracy for angles
he present algorithm uses both cross product
deliver full accuracy whatever the size of


product of two p-vectors
 of p-vector
product of two p-vectors
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) If either vector is null, a zero result is returned.

  2) The angular separation is most simply formulated in terms of
     scalar product.  However, this gives poor accuracy for angles
     near zero and pi.  The present algorithm uses both cross product
     and dot product, to deliver full accuracy whatever the size of
     the angle.


------
seps
------

Given:
first longitude (radians)
first latitude (radians)
second longitude (radians)
second latitude (radians)
e):
angular separation (radians)

al coordinates to unit vector
 separation between two p-vectors
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
sp00
------

Given:
   TT as a 2-part Julian Date (Note 1)
e):
   the TIO locator s' in radians (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
s obtained from polar motion observations by
n, and so is in essence unpredictable.
ated by a secular drift of about
er century, which is the approximation
sent function.

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
starpm
------

Given:
ght ascension (radians), before
clination (radians), before
 proper motion (radians/year), before
c proper motion (radians/year), before
rallax (arcseconds), before
dial velocity (km/s, +ve = receding), before
efore" epoch, part A (Note 1)
efore" epoch, part B (Note 1)
fter" epoch, part A (Note 1)
fter" epoch, part B (Note 1)

ght ascension (radians), after
clination (radians), after
 proper motion (radians/year), after
c proper motion (radians/year), after
rallax (arcseconds), after
dial velocity (km/s, +ve = receding), after
e):
atus:
 -1 = system error (should not occur)
  0 = no warnings or errors
  1 = distance overridden (Note 6)
  2 = excessive velocity (Note 7)
  4 = solution didn't converge (Note 8)
lse = binary logical OR of the above warnings

ing TDB dates ep1a+ep1b and ep2a+ep2b are
ioned in any convenient way between the two
r example, JD(TDB)=2450123.7 could be
these ways, among others:
  epnb
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
ormal star-catalog conventions, the object's
declination are freed from the effects of
 The frame, which is aligned to the catalog
 is Lorentzian and centered on the SSB.
re the rate of change of the right ascension
he catalog epoch and are in radians per TDB

ial velocity are in the same frame.
units.  The star coordinates are in radians
ns in radians per Julian year, but the
conds.
 is in terms of coordinate angle, not true
og uses arcseconds for both RA and Dec proper
er motion will need to be divided by cos(Dec)

 at constant speed, in the inertial frame,

or zero or negative) parallax is interpreted
ect is on the "celestial sphere", the radius
rary (large) value (see the eraStarpv
ue used).  When the distance is overridden in
, initially zero, has 1 added to it.
y is a significant fraction of c (see the
 function eraStarpv), it is arbitrarily set
action occurs, 2 is added to the status.
ustment carried out in the eraStarpv function
e calculation.  If the process fails to
t number of iterations, 4 is added to the


talog data to space motion pv-vector
a pv-vector
product of two p-vectors
otion pv-vector to star catalog data
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The starting and ending TDB dates ep1a+ep1b and ep2a+ep2b are
     Julian Dates, apportioned in any convenient way between the two
     parts (A and B).  For example, JD(TDB)=2450123.7 could be
     expressed in any of these ways, among others:

             epna          epnb

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) In accordance with normal star-catalog conventions, the object's
     right ascension and declination are freed from the effects of
     secular aberration.  The frame, which is aligned to the catalog
     equator and equinox, is Lorentzian and centered on the SSB.

     The proper motions are the rate of change of the right ascension
     and declination at the catalog epoch and are in radians per TDB
     Julian year.

     The parallax and radial velocity are in the same frame.

  3) Care is needed with units.  The star coordinates are in radians
     and the proper motions in radians per Julian year, but the
     parallax is in arcseconds.

  4) The RA proper motion is in terms of coordinate angle, not true
     angle.  If the catalog uses arcseconds for both RA and Dec proper
     motions, the RA proper motion will need to be divided by cos(Dec)
     before use.

  5) Straight-line motion at constant speed, in the inertial frame,
     is assumed.

  6) An extremely small (or zero or negative) parallax is interpreted
     to mean that the object is on the "celestial sphere", the radius
     of which is an arbitrary (large) value (see the eraStarpv
     function for the value used).  When the distance is overridden in
     this way, the status, initially zero, has 1 added to it.

  7) If the space velocity is a significant fraction of c (see the
     constant VMAX in the function eraStarpv), it is arbitrarily set
     to zero.  When this action occurs, 2 is added to the status.

  8) The relativistic adjustment carried out in the eraStarpv function
     involves an iterative calculation.  If the process fails to
     converge within a set number of iterations, 4 is added to the
     status.


------
starpv
------


Notes:

  1) The star data accepted by this function are "observables" for an
     imaginary observer at the solar-system barycenter.  Proper motion
     and radial velocity are, strictly, in terms of barycentric
     coordinate time, TCB.  For most practical applications, it is
     permissible to neglect the distinction between TCB and ordinary
     "proper" time on Earth (TT/TAI).  The result will, as a rule, be
     limited by the intrinsic accuracy of the proper-motion and
     radial-velocity data;  moreover, the pv-vector is likely to be
     merely an intermediate result, so that a change of time unit
     would cancel out overall.

     In accordance with normal star-catalog conventions, the object's
     right ascension and declination are freed from the effects of
     secular aberration.  The frame, which is aligned to the catalog
     equator and equinox, is Lorentzian and centered on the SSB.

  2) The resulting position and velocity pv-vector is with respect to
     the same frame and, like the catalog coordinates, is freed from
     the effects of secular aberration.  Should the "coordinate
     direction", where the object was located at the catalog epoch, be
     required, it may be obtained by calculating the magnitude of the
     position vector pv[0][0-2] dividing by the speed of light in
     AU/day to give the light-time, and then multiplying the space
     velocity pv[1][0-2] by this light-time and adding the result to
     pv[0][0-2].

     Summarizing, the pv-vector returned is for most stars almost
     identical to the result of applying the standard geometrical
     "space motion" transformation.  The differences, which are the
     subject of the Stumpff paper referenced below, are:

     (i) In stars with significant radial velocity and proper motion,
     the constantly changing light-time distorts the apparent proper
     motion.  Note that this is a classical, not a relativistic,
     effect.

     (ii) The transformation complies with special relativity.

  3) Care is needed with units.  The star coordinates are in radians
     and the proper motions in radians per Julian year, but the
     parallax is in arcseconds; the radial velocity is in km/s, but
     the pv-vector result is in AU and AU/day.

  4) The RA proper motion is in terms of coordinate angle, not true
     angle.  If the catalog uses arcseconds for both RA and Dec proper
     motions, the RA proper motion will need to be divided by cos(Dec)
     before use.

  5) Straight-line motion at constant speed, in the inertial frame,
     is assumed.

  6) An extremely small (or zero or negative) parallax is interpreted
     to mean that the object is on the "celestial sphere", the radius
     of which is an arbitrary (large) value (see the constant PXMIN).
     When the distance is overridden in this way, the status,
     initially zero, has 1 added to it.

  7) If the space velocity is a significant fraction of c (see the
     constant VMAX), it is arbitrarily set to zero.  When this action
     occurs, 2 is added to the status.

  8) The relativistic adjustment involves an iterative calculation.
     If the process fails to converge within a set number (IMAX) of
     iterations, 4 is added to the status.

  9) The inverse transformation is performed by the function
     eraPvstar.


------
sxp
------

Given:
 scalar
 p-vector

 s * p

r p and sp to be the same array.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
sxpv
------

Given:
    scalar
    pv-vector

    s * pv

r pv and spv to be the same array

y pv-vector by two scalars
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
taitt
------

Given:
 TAI as a 2-part Julian Date

 TT as a 2-part Julian Date
e):
 status:  0 = OK

Date, apportioned in any convenient way
ments, for example where tai1 is the Julian
is the fraction of a day.  The returned


it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
nt to the Astronomical Almanac,
n (ed), University Science Books (1992)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Explanatory Supplement to the Astronomical Almanac,
     P. Kenneth Seidelmann (ed), University Science Books (1992)

------
taiut1
------

Given:
 TAI as a 2-part Julian Date
 UT1-TAI in seconds

 UT1 as a 2-part Julian Date
e):
 status:  0 = OK

Date, apportioned in any convenient way
ments, for example where tai1 is the Julian
is the fraction of a day.  The returned
t.
e. UT1-TAI, is an observed quantity, and is
tabulations.

nt to the Astronomical Almanac,
n (ed), University Science Books (1992)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
taiutc
------

Given:
TAI as a 2-part Julian Date (Note 1)

UTC as a 2-part quasi Julian Date (Notes 1-3)
e):
status: +1 = dubious year (Note 4)
         0 = OK
        -1 = unacceptable date

Date, apportioned in any convenient way
ments, for example where tai1 is the Julian
is the fraction of a day.  The returned utc1
logous pair, except that a special convention
h the problem of leap seconds - see the next

sly represent UTC during a leap second unless
 taken.  The convention in the present
 JD day represents UTC days whether the
00 or 86401 SI seconds.  In the 1960-1972 era
umps (in either direction) each time the
ession was changed, and these "mini-leaps"
 the ERFA convention.
f can be used to transform the UTC quasi-JD
nd clock time, including UTC leap second

dubious year" flags UTCs that predate the
time scale or that are too far in the future
eraDat for further details.

TAI

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
nt to the Astronomical Almanac,
n (ed), University Science Books (1992)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) tai1+tai2 is Julian Date, apportioned in any convenient way
     between the two arguments, for example where tai1 is the Julian
     Day Number and tai2 is the fraction of a day.  The returned utc1
     and utc2 form an analogous pair, except that a special convention
     is used, to deal with the problem of leap seconds - see the next
     note.

  2) JD cannot unambiguously represent UTC during a leap second unless
     special measures are taken.  The convention in the present
     function is that the JD day represents UTC days whether the
     length is 86399, 86400 or 86401 SI seconds.  In the 1960-1972 era
     there were smaller jumps (in either direction) each time the
     linear UTC(TAI) expression was changed, and these "mini-leaps"
     are also included in the ERFA convention.

  3) The function eraD2dtf can be used to transform the UTC quasi-JD
     into calendar date and clock time, including UTC leap second
     handling.

  4) The warning status "dubious year" flags UTCs that predate the
     introduction of the time scale or that are too far in the future
     to be trusted.  See eraDat for further details.


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Explanatory Supplement to the Astronomical Almanac,
     P. Kenneth Seidelmann (ed), University Science Books (1992)

------
tcbtdb
------

Given:
 TCB as a 2-part Julian Date

 TDB as a 2-part Julian Date
e):
 status:  0 = OK

Date, apportioned in any convenient way
ments, for example where tcb1 is the Julian
is the fraction of a day.  The returned
t.
 Assembly introduced a conventional linear
en TDB and TCB.  This transformation
drift between TCB and terrestrial time TT,
imately centered on TT.  Because the
 TT and TCB depends on the adopted solar
e degree of alignment between TDB and TT over
vary according to which ephemeris is used.
f TDB attempted to avoid this problem by
 and TT should differ only by periodic
good description of the nature of the
ded precise mathematical formulation.  The
relationship adopted in 2006 sidestepped
hilst delivering a TDB that in practice was
es before that date.
he same as Teph, the time argument for the
emerides.

B3
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
tcgtt
------

Given:
 TCG as a 2-part Julian Date

 TT as a 2-part Julian Date
e):
 status:  0 = OK

Date, apportioned in any convenient way
ments, for example where tcg1 is the Julian
 is the fraction of a day.  The returned


it, G. (eds.), IERS Conventions (2003),.
No. 32, BKG (2004)
B1.9
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),.
     IERS Technical Note No. 32, BKG (2004)

     IAU 2000 Resolution B1.9

------
tdbtcb
------

Given:
 TDB as a 2-part Julian Date

 TCB as a 2-part Julian Date
e):
 status:  0 = OK

Date, apportioned in any convenient way
ments, for example where tdb1 is the Julian
is the fraction of a day.  The returned
t.
 Assembly introduced a conventional linear
en TDB and TCB.  This transformation
drift between TCB and terrestrial time TT,
imately centered on TT.  Because the
 TT and TCB depends on the adopted solar
e degree of alignment between TDB and TT over
vary according to which ephemeris is used.
f TDB attempted to avoid this problem by
 and TT should differ only by periodic
good description of the nature of the
ded precise mathematical formulation.  The
relationship adopted in 2006 sidestepped
hilst delivering a TDB that in practice was
es before that date.
he same as Teph, the time argument for the
emerides.

B3
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
tdbtt
------

Given:
 TDB as a 2-part Julian Date
 TDB-TT in seconds

 TT as a 2-part Julian Date
e):
 status:  0 = OK

Date, apportioned in any convenient way
ments, for example where tdb1 is the Julian
is the fraction of a day.  The returned

resents the quasi-periodic component of the
tween TT and TCB.  It is dependent upon the
 ephemeris, and can be obtained by numerical
rrogating a precomputed time ephemeris or by
uch as that implemented in the ERFA function
ity is dominated by an annual term of 1.7 ms

he same as Teph, the time argument for the
emerides.

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
3
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) tdb1+tdb2 is Julian Date, apportioned in any convenient way
     between the two arguments, for example where tdb1 is the Julian
     Day Number and tdb2 is the fraction of a day.  The returned
     tt1,tt2 follow suit.

  2) The argument dtr represents the quasi-periodic component of the
     GR transformation between TT and TCB.  It is dependent upon the
     adopted solar-system ephemeris, and can be obtained by numerical
     integration, by interrogating a precomputed time ephemeris or by
     evaluating a model such as that implemented in the ERFA function
     eraDtdb.   The quantity is dominated by an annual term of 1.7 ms
     amplitude.

  3) TDB is essentially the same as Teph, the time argument for the
     JPL solar system ephemerides.


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     IAU 2006 Resolution 3

------
test_erfa
------



------
tf2a
------

Given:
gn:  '-' = negative, otherwise positive
urs
nutes
conds

gle in radians
e):
atus:  0 = OK
       1 = ihour outside range 0-23
       2 = imin outside range 0-59
       3 = sec outside range 0-59.999...

ted even if any of the range checks fail.
n and/or sec produce a warning status, but
is used in the conversion.
le errors, the status value reflects only the
 taking precedence.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
tf2d
------

Given:
gn:  '-' = negative, otherwise positive
urs
nutes
conds

terval in days
e):
atus:  0 = OK
       1 = ihour outside range 0-23
       2 = imin outside range 0-59
       3 = sec outside range 0-59.999...

ted even if any of the range checks fail.
n and/or sec produce a warning status, but
is used in the conversion.
le errors, the status value reflects only the
 taking precedence.
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
tr
------

Given:
]    r-matrix

]    transpose

r r and rt to be the same array.

matrix
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
trxp
------

Given:
]   r-matrix
    p-vector

    r * p

r p and trp to be the same array.

se r-matrix
 of r-matrix and p-vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
trxpv
------

Given:
]    r-matrix
]    pv-vector

]    r * pv

r pv and trpv to be the same array.

se r-matrix
 of r-matrix and pv-vector
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
tttai
------

Given:
 TT as a 2-part Julian Date

 TAI as a 2-part Julian Date
e):
 status:  0 = OK

te, apportioned in any convenient way between
or example where tt1 is the Julian Day Number
ion of a day.  The returned tai1,tai2 follow


it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
nt to the Astronomical Almanac,
n (ed), University Science Books (1992)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Explanatory Supplement to the Astronomical Almanac,
     P. Kenneth Seidelmann (ed), University Science Books (1992)

------
tttcg
------

Given:
 TT as a 2-part Julian Date

 TCG as a 2-part Julian Date
e):
 status:  0 = OK

te, apportioned in any convenient way between
or example where tt1 is the Julian Day Number
ion of a day.  The returned tcg1,tcg2 follow


it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
B1.9
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     IAU 2000 Resolution B1.9

------
tttdb
------

Given:
 TT as a 2-part Julian Date
 TDB-TT in seconds

 TDB as a 2-part Julian Date
e):
 status:  0 = OK

te, apportioned in any convenient way between
or example where tt1 is the Julian Day Number
ion of a day.  The returned tdb1,tdb2 follow

resents the quasi-periodic component of the
tween TT and TCB.  It is dependent upon the
 ephemeris, and can be obtained by numerical
rrogating a precomputed time ephemeris or by
uch as that implemented in the ERFA function
ity is dominated by an annual term of 1.7 ms

he same as Teph, the time argument for the JPL
ides.

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
3
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) tt1+tt2 is Julian Date, apportioned in any convenient way between
     the two arguments, for example where tt1 is the Julian Day Number
     and tt2 is the fraction of a day.  The returned tdb1,tdb2 follow
     suit.

  2) The argument dtr represents the quasi-periodic component of the
     GR transformation between TT and TCB.  It is dependent upon the
     adopted solar-system ephemeris, and can be obtained by numerical
     integration, by interrogating a precomputed time ephemeris or by
     evaluating a model such as that implemented in the ERFA function
     eraDtdb.   The quantity is dominated by an annual term of 1.7 ms
     amplitude.

  3) TDB is essentially the same as Teph, the time argument for the JPL
     solar system ephemerides.


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     IAU 2006 Resolution 3

------
ttut1
------

Given:
 TT as a 2-part Julian Date
 TT-UT1 in seconds

 UT1 as a 2-part Julian Date
e):
 status:  0 = OK

te, apportioned in any convenient way between
or example where tt1 is the Julian Day Number
ion of a day.  The returned ut11,ut12 follow

lassical Delta T.

nt to the Astronomical Almanac,
n (ed), University Science Books (1992)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
ut1tai
------

Given:
 UT1 as a 2-part Julian Date
 UT1-TAI in seconds

 TAI as a 2-part Julian Date
e):
 status:  0 = OK

Date, apportioned in any convenient way
ments, for example where ut11 is the Julian
is the fraction of a day.  The returned
t.
e. UT1-TAI, is an observed quantity, and is
tabulations.

nt to the Astronomical Almanac,
n (ed), University Science Books (1992)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
ut1tt
------

Given:
 UT1 as a 2-part Julian Date
 TT-UT1 in seconds

 TT as a 2-part Julian Date
e):
 status:  0 = OK

Date, apportioned in any convenient way
ments, for example where ut11 is the Julian
is the fraction of a day.  The returned

lassical Delta T.

nt to the Astronomical Almanac,
n (ed), University Science Books (1992)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.

------
ut1utc
------

Given:
UT1 as a 2-part Julian Date (Note 1)
Delta UT1: UT1-UTC in seconds (Note 2)

UTC as a 2-part quasi Julian Date (Notes 3,4)
e):
status: +1 = dubious year (Note 5)
         0 = OK
        -1 = unacceptable date

Date, apportioned in any convenient way
ments, for example where ut11 is the Julian
is the fraction of a day.  The returned utc1
logous pair, except that a special convention
h the problem of leap seconds - see Note 3.
ained from tabulations provided by the
Rotation and Reference Systems Service.  The
ly by 1s at a leap second;  however, close to
gorithm used here is tolerant of the "wrong"
g made.
sly represent UTC during a leap second unless
 taken.  The convention in the present
 returned quasi JD day UTC1+UTC2 represents
 length is 86399, 86400 or 86401 SI seconds.
f can be used to transform the UTC quasi-JD
nd clock time, including UTC leap second

dubious year" flags UTCs that predate the
time scale or that are too far in the future
eraDat for further details.

regorian calendar
T) = TAI-UTC
an calendar to JD

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
nt to the Astronomical Almanac,
n (ed), University Science Books (1992)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) ut11+ut12 is Julian Date, apportioned in any convenient way
     between the two arguments, for example where ut11 is the Julian
     Day Number and ut12 is the fraction of a day.  The returned utc1
     and utc2 form an analogous pair, except that a special convention
     is used, to deal with the problem of leap seconds - see Note 3.

  2) Delta UT1 can be obtained from tabulations provided by the
     International Earth Rotation and Reference Systems Service.  The
     value changes abruptly by 1s at a leap second;  however, close to
     a leap second the algorithm used here is tolerant of the "wrong"
     choice of value being made.

  3) JD cannot unambiguously represent UTC during a leap second unless
     special measures are taken.  The convention in the present
     function is that the returned quasi JD day UTC1+UTC2 represents
     UTC days whether the length is 86399, 86400 or 86401 SI seconds.

  4) The function eraD2dtf can be used to transform the UTC quasi-JD
     into calendar date and clock time, including UTC leap second
     handling.

  5) The warning status "dubious year" flags UTCs that predate the
     introduction of the time scale or that are too far in the future
     to be trusted.  See eraDat for further details.


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Explanatory Supplement to the Astronomical Almanac,
     P. Kenneth Seidelmann (ed), University Science Books (1992)

------
utctai
------

Given:
UTC as a 2-part quasi Julian Date (Notes 1-4)

TAI as a 2-part Julian Date (Note 5)
e):
status: +1 = dubious year (Note 3)
         0 = OK
        -1 = unacceptable date

ulian Date (see Note 2), apportioned in any
en the two arguments, for example where utc1
mber and utc2 is the fraction of a day.
sly represent UTC during a leap second unless
 taken.  The convention in the present
 JD day represents UTC days whether the
00 or 86401 SI seconds.  In the 1960-1972 era
umps (in either direction) each time the
ession was changed, and these "mini-leaps"
 the ERFA convention.
dubious year" flags UTCs that predate the
time scale or that are too far in the future
eraDat for further details.
d converts from calendar date and time of day
ate, and in the case of UTC implements the
y convention described above.
I2 are such that their sum is the TAI Julian


regorian calendar
T) = TAI-UTC
an calendar to JD

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
nt to the Astronomical Almanac,
n (ed), University Science Books (1992)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
     convenient way between the two arguments, for example where utc1
     is the Julian Day Number and utc2 is the fraction of a day.

  2) JD cannot unambiguously represent UTC during a leap second unless
     special measures are taken.  The convention in the present
     function is that the JD day represents UTC days whether the
     length is 86399, 86400 or 86401 SI seconds.  In the 1960-1972 era
     there were smaller jumps (in either direction) each time the
     linear UTC(TAI) expression was changed, and these "mini-leaps"
     are also included in the ERFA convention.

  3) The warning status "dubious year" flags UTCs that predate the
     introduction of the time scale or that are too far in the future
     to be trusted.  See eraDat for further details.

  4) The function eraDtf2d converts from calendar date and time of day
     into 2-part Julian Date, and in the case of UTC implements the
     leap-second-ambiguity convention described above.

  5) The returned TAI1,TAI2 are such that their sum is the TAI Julian
     Date.


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Explanatory Supplement to the Astronomical Almanac,
     P. Kenneth Seidelmann (ed), University Science Books (1992)

------
utcut1
------

Given:
UTC as a 2-part quasi Julian Date (Notes 1-4)
Delta UT1 = UT1-UTC in seconds (Note 5)

UT1 as a 2-part Julian Date (Note 6)
e):
status: +1 = dubious year (Note 3)
         0 = OK
        -1 = unacceptable date

ulian Date (see Note 2), apportioned in any
en the two arguments, for example where utc1
mber and utc2 is the fraction of a day.
sly represent UTC during a leap second unless
 taken.  The convention in the present
 JD day represents UTC days whether the
00 or 86401 SI seconds.
dubious year" flags UTCs that predate the
time scale or that are too far in the future
eraDat for further details.
d converts from calendar date and time of
an Date, and in the case of UTC implements
guity convention described above.
ained from tabulations provided by the
Rotation and Reference Systems Service.
esponsibility to supply a dut1 argument
TC value that matches the given UTC.
12 are such that their sum is the UT1 Julian


it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
nt to the Astronomical Almanac,
n (ed), University Science Books (1992)

regorian calendar
T) = TAI-UTC
TAI
UT1
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
     convenient way between the two arguments, for example where utc1
     is the Julian Day Number and utc2 is the fraction of a day.

  2) JD cannot unambiguously represent UTC during a leap second unless
     special measures are taken.  The convention in the present
     function is that the JD day represents UTC days whether the
     length is 86399, 86400 or 86401 SI seconds.

  3) The warning status "dubious year" flags UTCs that predate the
     introduction of the time scale or that are too far in the future
     to be trusted.  See eraDat for further details.

  4) The function eraDtf2d converts from calendar date and time of
     day into 2-part Julian Date, and in the case of UTC implements
     the leap-second-ambiguity convention described above.

  5) Delta UT1 can be obtained from tabulations provided by the
     International Earth Rotation and Reference Systems Service.
     It is the caller's responsibility to supply a dut1 argument
     containing the UT1-UTC value that matches the given UTC.

  6) The returned ut11,ut12 are such that their sum is the UT1 Julian
     Date.

  References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Explanatory Supplement to the Astronomical Almanac,
     P. Kenneth Seidelmann (ed), University Science Books (1992)


References:

     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
     IERS Technical Note No. 32, BKG (2004)

     Explanatory Supplement to the Astronomical Almanac,
     P. Kenneth Seidelmann (ed), University Science Books (1992)

  Called:
     eraJd2cal    JD to Gregorian calendar
     eraDat       delta(AT) = TAI-UTC
     eraUtctai    UTC to TAI
     eraTaiut1    TAI to UT1

------
xy06
------

Given:
    TT as a 2-part Julian Date (Note 1)

    CIP X,Y coordinates (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
are those of the unit vector towards the
te pole.  They represent the combined effects
ssion and nutation.
ments used are as adopted in IERS Conventions
Simon et al. (1994) and Souchay et al.

ve to the angles-based method, via the ERFA
d as used in eraXys06a for example.  The two
 1 microarcsecond level (at present), a
mpared with the intrinsic accuracy of the
 would be unwise to mix the two methods
ries-based) in a single application.

omaly of the Moon
omaly of the Sun
gument of the latitude of the Moon
ongation of the Moon from the Sun
ngitude of the Moon's ascending node
ngitude of Mercury
ngitude of Venus
ngitude of Earth
ngitude of Mars
ngitude of Jupiter
ngitude of Saturn
ngitude of Uranus
ngitude of Neptune
 accumulated precession in longitude

ce, P.T. & Chapront, J., 2003,
12, 567
ace, P.T., 2006, Astron.Astrophys. 450, 855
it, G. (eds.), 2004, IERS Conventions (2003),
No. 32, BKG
on, P., Chapront, J., Chapront-Touze, M.,
, J., Astron.Astrophys., 1994, 282, 663
 B., Kinoshita, H., Folgueira, M., 1999,
p.Ser. 135, 111
taine, N., 2006, Astron.Astrophys. 459, 981
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The X,Y coordinates are those of the unit vector towards the
     celestial intermediate pole.  They represent the combined effects
     of frame bias, precession and nutation.

  3) The fundamental arguments used are as adopted in IERS Conventions
     (2003) and are from Simon et al. (1994) and Souchay et al.
     (1999).

  4) This is an alternative to the angles-based method, via the ERFA
     function eraFw2xy and as used in eraXys06a for example.  The two
     methods agree at the 1 microarcsecond level (at present), a
     negligible amount compared with the intrinsic accuracy of the
     models.  However, it would be unwise to mix the two methods
     (angles-based and series-based) in a single application.


References:

     Capitaine, N., Wallace, P.T. & Chapront, J., 2003,
     Astron.Astrophys., 412, 567

     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
     IERS Technical Note No. 32, BKG

     Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
     Francou, G. & Laskar, J., Astron.Astrophys., 1994, 282, 663

     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M., 1999,
     Astron.Astrophys.Supp.Ser. 135, 111

     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981

------
xys00a
------

Given:
  TT as a 2-part Julian Date (Note 1)

  Celestial Intermediate Pole (Note 2)
  the CIO locator s (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
ediate Pole coordinates are the x,y
it vector in the Geocentric Celestial

n radians) positions the Celestial
on the equator of the CIP.
ly less accurate result (about 1 mas for
d by using instead the eraXys00b function.

al NPB matrix, IAU 2000A
 CIP X,Y coordinates from NPB matrix
 locator s, given X,Y, IAU 2000A

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The Celestial Intermediate Pole coordinates are the x,y
     components of the unit vector in the Geocentric Celestial
     Reference System.

  3) The CIO locator s (in radians) positions the Celestial
     Intermediate Origin on the equator of the CIP.

  4) A faster, but slightly less accurate result (about 1 mas for
     X,Y), can be obtained by using instead the eraXys00b function.


------
xys00b
------

Given:
  TT as a 2-part Julian Date (Note 1)

  Celestial Intermediate Pole (Note 2)
  the CIO locator s (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
ediate Pole coordinates are the x,y
it vector in the Geocentric Celestial

n radians) positions the Celestial
on the equator of the CIP.
 is faster, but slightly less accurate (about
the eraXys00a function.

al NPB matrix, IAU 2000B
 CIP X,Y coordinates from NPB matrix
 locator s, given X,Y, IAU 2000A

it, G. (eds.), IERS Conventions (2003),
No. 32, BKG (2004)
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The Celestial Intermediate Pole coordinates are the x,y
     components of the unit vector in the Geocentric Celestial
     Reference System.

  3) The CIO locator s (in radians) positions the Celestial
     Intermediate Origin on the equator of the CIP.

  4) The present function is faster, but slightly less accurate (about
     1 mas in X,Y), than the eraXys00a function.


------
xys06a
------

Given:
 TT as a 2-part Julian Date (Note 1)

 Celestial Intermediate Pole (Note 2)
 the CIO locator s (Note 2)

te2 is a Julian Date, apportioned in any
en the two arguments.  For example,
ld be expressed in any of these ways,

  date2
    0.0       (JD method)
-1421.3       (J2000 method)
50123.2       (MJD method)
    0.2       (date & time method)
 most natural and convenient to use in
 of several decimal digits of resolution
J2000 method is best matched to the way
led internally and will deliver the
 The MJD method and the date & time methods
mises between resolution and convenience.
ediate Pole coordinates are the x,y components
n the Geocentric Celestial Reference System.
n radians) positions the Celestial
on the equator of the CIP.
ns for generating X and Y are also available:
ace (2006) and eraXy06.

al NPB matrix, IAU 2006/2000A
 CIP X,Y coordinates from NPB matrix
 locator s, given X,Y, IAU 2006

ace, P.T., 2006, Astron.Astrophys. 450, 855
taine, N., 2006, Astron.Astrophys. 459, 981
, NumFOCUS Foundation.
n, from the SOFA library.  See notes at end of file.
Notes:

  1) The TT date date1+date2 is a Julian Date, apportioned in any
     convenient way between the two arguments.  For example,
     JD(TT)=2450123.7 could be expressed in any of these ways,
     among others:

            date1          date2

         2450123.7           0.0       (JD method)
         2451545.0       -1421.3       (J2000 method)
         2400000.5       50123.2       (MJD method)
         2450123.5           0.2       (date & time method)

     The JD method is the most natural and convenient to use in
     cases where the loss of several decimal digits of resolution
     is acceptable.  The J2000 method is best matched to the way
     the argument is handled internally and will deliver the
     optimum resolution.  The MJD method and the date & time methods
     are both good compromises between resolution and convenience.

  2) The Celestial Intermediate Pole coordinates are the x,y components
     of the unit vector in the Geocentric Celestial Reference System.

  3) The CIO locator s (in radians) positions the Celestial
     Intermediate Origin on the equator of the CIP.

  4) Series-based solutions for generating X and Y are also available:
     see Capitaine & Wallace (2006) and eraXy06.


References:

     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981

------
zp
------



------
zpv
------



------
zr
------



