==========
erfa_c API
==========

-------
A2af
-------

Decompose radians into degrees, arcminutes, arcseconds, fraction.

Given:
   ndp     int     resolution (Note 1)
   angle   double  angle in radians

Returned:
   sign    char    '+' or '-'
   idmsf   int[4]  degrees, arcminutes, arcseconds, fraction

Called:
   eraD2tf      decompose days to hms

Notes:

1) The argument ndp is interpreted as follows:

   ndp         resolution
    :      ...0000 00 00
   -7         1000 00 00
   -6          100 00 00
   -5           10 00 00
   -4            1 00 00
   -3            0 10 00
   -2            0 01 00
   -1            0 00 10
    0            0 00 01
    1            0 00 00.1
    2            0 00 00.01
    3            0 00 00.001
    :            0 00 00.000...

2) The largest positive useful value for ndp is determined by the
   size of angle, the format of doubles on the target platform, and
   the risk of overflowing idmsf[3].  On a typical platform, for
   angle up to 2pi, the available floating-point precision might
   correspond to ndp=12.  However, the practical limit is typically
   ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
   only 16 bits.

3) The absolute value of angle may exceed 2pi.  In cases where it
   does not, it is up to the caller to test for and handle the
   case where angle is very nearly 2pi and rounds up to 360 degrees,
   by testing for idmsf[0]=360 and setting idmsf[0-3] to zero.



-------
A2tf
-------

Decompose radians into hours, minutes, seconds, fraction.

Given:
   ndp     int     resolution (Note 1)
   angle   double  angle in radians

Returned:
   sign    char    '+' or '-'
   ihmsf   int[4]  hours, minutes, seconds, fraction

Called:
   eraD2tf      decompose days to hms

Notes:

1) The argument ndp is interpreted as follows:

   ndp         resolution
    :      ...0000 00 00
   -7         1000 00 00
   -6          100 00 00
   -5           10 00 00
   -4            1 00 00
   -3            0 10 00
   -2            0 01 00
   -1            0 00 10
    0            0 00 01
    1            0 00 00.1
    2            0 00 00.01
    3            0 00 00.001
    :            0 00 00.000...

2) The largest positive useful value for ndp is determined by the
   size of angle, the format of doubles on the target platform, and
   the risk of overflowing ihmsf[3].  On a typical platform, for
   angle up to 2pi, the available floating-point precision might
   correspond to ndp=12.  However, the practical limit is typically
   ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
   only 16 bits.

3) The absolute value of angle may exceed 2pi.  In cases where it
   does not, it is up to the caller to test for and handle the
   case where angle is very nearly 2pi and rounds up to 24 hours,
   by testing for ihmsf[0]=24 and setting ihmsf(0-3) to zero.



-------
Af2a
-------

Convert degrees, arcminutes, arcseconds to radians.

Given:
   s         char    sign:  '-' = negative, otherwise positive
   ideg      int     degrees
   iamin     int     arcminutes
   asec      double  arcseconds

Returned:
   rad       double  angle in radians

Returned (function value):
             int     status:  0 = OK
                              1 = ideg outside range 0-359
                              2 = iamin outside range 0-59
                              3 = asec outside range 0-59.999...

Notes:

1)  The result is computed even if any of the range checks fail.

2)  Negative ideg, iamin and/or asec produce a warning status, but
    the absolute value is used in the conversion.

3)  If there are multiple errors, the status value reflects only the
    first, the smallest taking precedence.



------
Anp
------

Normalize angle into the range 0 <= a < 2pi.

Given:
   a        double     angle (radians)

Returned (function value):
            double     angle in range 0-2pi



-------
Anpm
-------

Normalize angle into the range -pi <= a < +pi.

Given:
   a        double     angle (radians)

Returned (function value):
            double     angle in range +/-pi



-------
Bi00
-------

Frame bias components of IAU 2000 precession-nutation models (part
of MHB2000 with additions).

Returned:
   dpsibi,depsbi  double  longitude and obliquity corrections
   dra            double  the ICRS RA of the J2000.0 mean equinox

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



-------
Bp00
-------

Frame bias and precession, IAU 2000.

Given:
   date1,date2  double         TT as a 2-part Julian Date (Note 1)

Returned:
   rb           double[3][3]   frame bias matrix (Note 2)
   rp           double[3][3]   precession matrix (Note 3)
   rbp          double[3][3]   bias-precession matrix (Note 4)

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

Called:
   eraBi00      frame bias components, IAU 2000
   eraPr00      IAU 2000 precession adjustments
   eraIr        initialize r-matrix to identity
   eraRx        rotate around X-axis
   eraRy        rotate around Y-axis
   eraRz        rotate around Z-axis
   eraCr        copy r-matrix
   eraRxr       product of two r-matrices

Reference:
   "Expressions for the Celestial Intermediate Pole and Celestial
   Ephemeris Origin consistent with the IAU 2000A precession-
   nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

   n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.



-------
Bp06
-------

Frame bias and precession, IAU 2006.

Given:
   date1,date2  double         TT as a 2-part Julian Date (Note 1)

Returned:
   rb           double[3][3]   frame bias matrix (Note 2)
   rp           double[3][3]   precession matrix (Note 3)
   rbp          double[3][3]   bias-precession matrix (Note 4)

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

Called:
   eraPfw06     bias-precession F-W angles, IAU 2006
   eraFw2m      F-W angles to r-matrix
   eraPmat06    PB matrix, IAU 2006
   eraTr        transpose r-matrix
   eraRxr       product of two r-matrices

References:

   Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

   Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981



---------
Bpn2xy
---------

Extract from the bias-precession-nutation matrix the X,Y coordinates
of the Celestial Intermediate Pole.

Given:
   rbpn      double[3][3]  celestial-to-true matrix (Note 1)

Returned:
   x,y       double        Celestial Intermediate Pole (Note 2)

Notes:

1) The matrix rbpn transforms vectors from GCRS to true equator (and
   CIO or equinox) of date, and therefore the Celestial Intermediate
   Pole unit vector is the bottom row of the matrix.

2) The arguments x,y are components of the Celestial Intermediate
   Pole unit vector in the Geocentric Celestial Reference System.

Reference:

   "Expressions for the Celestial Intermediate Pole and Celestial
   Ephemeris Origin consistent with the IAU 2000A precession-
   nutation model", Astron.Astrophys. 400, 1145-1154
   (2003)

   n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.



---------
C2i00a
---------

Form the celestial-to-intermediate matrix for a given date using the
IAU 2000A precession-nutation model.

Given:
   date1,date2 double       TT as a 2-part Julian Date (Note 1)

Returned:
   rc2i        double[3][3] celestial-to-intermediate matrix (Note 2)

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

Called:
   eraPnm00a    classical NPB matrix, IAU 2000A
   eraC2ibpn    celestial-to-intermediate matrix, given NPB matrix

References:

   "Expressions for the Celestial Intermediate Pole and Celestial
   Ephemeris Origin consistent with the IAU 2000A precession-
   nutation model", Astron.Astrophys. 400, 1145-1154
   (2003)

   n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



---------
C2i00b
---------

Form the celestial-to-intermediate matrix for a given date using the
IAU 2000B precession-nutation model.

Given:
   date1,date2 double       TT as a 2-part Julian Date (Note 1)

Returned:
   rc2i        double[3][3] celestial-to-intermediate matrix (Note 2)

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

Called:
   eraPnm00b    classical NPB matrix, IAU 2000B
   eraC2ibpn    celestial-to-intermediate matrix, given NPB matrix

References:

   "Expressions for the Celestial Intermediate Pole and Celestial
   Ephemeris Origin consistent with the IAU 2000A precession-
   nutation model", Astron.Astrophys. 400, 1145-1154
   (2003)

   n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



---------
c2i06a
---------

Form the celestial-to-intermediate matrix for a given date using the
IAU 2006 precession and IAU 2000A nutation models.

Given:
   date1,date2 double       TT as a 2-part Julian Date (Note 1)

Returned:
   rc2i        double[3][3] celestial-to-intermediate matrix (Note 2)

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

Called:
   eraPnm06a    classical NPB matrix, IAU 2006/2000A
   eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
   eraS06       the CIO locator s, Given X,Y, IAU 2006
   eraC2ixys    celestial-to-intermediate matrix, Given X,Y and s

References:

   McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
   IERS Technical Note No. 32, BKG



---------
C2ibpn
---------

Form the celestial-to-intermediate matrix for a given date given
the bias-precession-nutation matrix.  IAU 2000.

Given:
   date1,date2 double       TT as a 2-part Julian Date (Note 1)
   rbpn        double[3][3] celestial-to-true matrix (Note 2)

Returned:
   rc2i        double[3][3] celestial-to-intermediate matrix (Note 3)

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

Called:
   eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
   eraC2ixy     celestial-to-intermediate matrix, given X,Y

References:
   "Expressions for the Celestial Intermediate Pole and Celestial
   Ephemeris Origin consistent with the IAU 2000A precession-
   nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

   n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



--------
C2ixy
--------

Form the celestial to intermediate-frame-of-date matrix for a given
date when the CIP X,Y coordinates are known.  IAU 2000.

Given:
   date1,date2 double       TT as a 2-part Julian Date (Note 1)
   x,y         double       Celestial Intermediate Pole (Note 2)

Returned:
   rc2i        double[3][3] celestial-to-intermediate matrix (Note 3)

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

Called:
   eraC2ixys    celestial-to-intermediate matrix, given X,Y and s
   eraS00       the CIO locator s, given X,Y, IAU 2000A

Reference:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



---------
C2ixys
---------

Form the celestial to intermediate-frame-of-date matrix given the CIP
X,Y and the CIO locator s.

Given:
   x,y      double         Celestial Intermediate Pole (Note 1)
   s        double         the CIO locator s (Note 2)

Returned:
   rc2i     double[3][3]   celestial-to-intermediate matrix (Note 3)

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

Called:
   eraIr        initialize r-matrix to identity
   eraRz        rotate around Z-axis
   eraRy        rotate around Y-axis

Reference:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



------
C2s
------

P-vector to spherical coordinates.

Given:
   p      double[3]    p-vector

Returned:
   theta  double       longitude angle (radians)
   phi    double       latitude angle (radians)

Notes:

1) The vector p can have any magnitude; only its direction is used.

2) If p is null, zero theta and phi are returned.

3) At either pole, zero theta is returned.



---------
C2t00a
---------

Form the celestial to terrestrial matrix given the date, the UT1 and
the polar motion, using the IAU 2000A nutation model.

Given:
   tta,ttb  double         TT as a 2-part Julian Date (Note 1)
   uta,utb  double         UT1 as a 2-part Julian Date (Note 1)
   xp,yp    double         coordinates of the pole (radians, Note 2)

Returned:
   rc2t     double[3][3]   celestial-to-terrestrial matrix (Note 3)

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

Called:
   eraC2i00a    celestial-to-intermediate matrix, IAU 2000A
   eraEra00     Earth rotation angle, IAU 2000
   eraSp00      the TIO locator s', IERS 2000
   eraPom00     polar motion matrix
   eraC2tcio    form CIO-based celestial-to-terrestrial matrix

Reference:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



---------
C2t00b
---------

Form the celestial to terrestrial matrix given the date, the UT1 and
the polar motion, using the IAU 2000B nutation model.

Given:
   tta,ttb  double         TT as a 2-part Julian Date (Note 1)
   uta,utb  double         UT1 as a 2-part Julian Date (Note 1)
   xp,yp    double         coordinates of the pole (radians, Note 2)

Returned:
   rc2t     double[3][3]   celestial-to-terrestrial matrix (Note 3)

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

Called:
   eraC2i00b    celestial-to-intermediate matrix, IAU 2000B
   eraEra00     Earth rotation angle, IAU 2000
   eraPom00     polar motion matrix
   eraC2tcio    form CIO-based celestial-to-terrestrial matrix

Reference:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



---------
C2t06a
---------

Form the celestial to terrestrial matrix given the date, the UT1 and
the polar motion, using the IAU 2006 precession and IAU 2000A
nutation models.

Given:
   tta,ttb  double         TT as a 2-part Julian Date (Note 1)
   uta,utb  double         UT1 as a 2-part Julian Date (Note 1)
   xp,yp    double         coordinates of the pole (radians, Note 2)

Returned:
   rc2t     double[3][3]   celestial-to-terrestrial matrix (Note 3)

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

Called:
   eraC2i06a    celestial-to-intermediate matrix, IAU 2006/2000A
   eraEra00     Earth rotation angle, IAU 2000
   eraSp00      the TIO locator s', IERS 2000
   eraPom00     polar motion matrix
   eraC2tcio    form CIO-based celestial-to-terrestrial matrix

Reference:

   McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
   IERS Technical Note No. 32, BKG



---------
C2tcio
---------

Assemble the celestial to terrestrial matrix from CIO-based
components (the celestial-to-intermediate matrix, the Earth Rotation
Angle and the polar motion matrix).

Given:
   rc2i     double[3][3]    celestial-to-intermediate matrix
   era      double          Earth rotation angle
   rpom     double[3][3]    polar-motion matrix

Returned:
   rc2t     double[3][3]    celestial-to-terrestrial matrix

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

Called:
   eraCr        copy r-matrix
   eraRz        rotate around Z-axis
   eraRxr       product of two r-matrices

Reference:

   McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
   IERS Technical Note No. 32, BKG



---------
C2teqx
---------

Assemble the celestial to terrestrial matrix from equinox-based
components (the celestial-to-true matrix, the Greenwich Apparent
Sidereal Time and the polar motion matrix).

Given:
   rbpn     double[3][3]    celestial-to-true matrix
   gst      double          Greenwich (apparent) Sidereal Time
   rpom     double[3][3]    polar-motion matrix

Returned:
   rc2t     double[3][3]    celestial-to-terrestrial matrix (Note 2)

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

Called:
   eraCr        copy r-matrix
   eraRz        rotate around Z-axis
   eraRxr       product of two r-matrices

Reference:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



--------
C2tpe
--------

Form the celestial to terrestrial matrix given the date, the UT1,
the nutation and the polar motion.  IAU 2000.

Given:
   tta,ttb    double        TT as a 2-part Julian Date (Note 1)
   uta,utb    double        UT1 as a 2-part Julian Date (Note 1)
   dpsi,deps  double        nutation (Note 2)
   xp,yp      double        coordinates of the pole (radians, Note 3)

Returned:
   rc2t       double[3][3]  celestial-to-terrestrial matrix (Note 4)

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

Called:
   eraPn00      bias/precession/nutation results, IAU 2000
   eraGmst00    Greenwich mean sidereal time, IAU 2000
   eraSp00      the TIO locator s', IERS 2000
   eraEe00      equation of the equinoxes, IAU 2000
   eraPom00     polar motion matrix
   eraC2teqx    form equinox-based celestial-to-terrestrial matrix

Reference:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



--------
C2txy
--------

Form the celestial to terrestrial matrix given the date, the UT1,
the CIP coordinates and the polar motion.  IAU 2000.

Given:
   tta,ttb  double         TT as a 2-part Julian Date (Note 1)
   uta,utb  double         UT1 as a 2-part Julian Date (Note 1)
   x,y      double         Celestial Intermediate Pole (Note 2)
   xp,yp    double         coordinates of the pole (radians, Note 3)

Returned:
   rc2t     double[3][3]   celestial-to-terrestrial matrix (Note 4)

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

Called:
   eraC2ixy     celestial-to-intermediate matrix, given X,Y
   eraEra00     Earth rotation angle, IAU 2000
   eraSp00      the TIO locator s', IERS 2000
   eraPom00     polar motion matrix
   eraC2tcio    form CIO-based celestial-to-terrestrial matrix

eference:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



---------
Cal2jd
---------

Gregorian Calendar to Julian Date.

Given:
   iy,im,id  int     year, month, day in Gregorian calendar (Note 1)

Returned:
   djm0      double  MJD zero-point: always 2400000.5
   djm       double  Modified Julian Date for 0 hrs

Returned (function value):
             int     status:
                         0 = OK
                        -1 = bad year   (Note 3: JD not computed)
                        -2 = bad month  (JD not computed)
                        -3 = bad day    (JD computed)

Notes:

1) The algorithm used is valid from -4800 March 1, but this
   implementation rejects dates before -4799 January 1.

2) The Julian Date is returned in two pieces, in the usual ERFA
   manner, which is designed to preserve time resolution.  The
   Julian Date is available as a single number by adding djm0 and
   djm.

3) In early eras the conversion is from the "Proleptic Gregorian
   Calendar";  no account is taken of the date(s) of adoption of
   the Gregorian Calendar, nor is the AD/BC numbering convention
   observed.

Reference:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992),
   Section 12.92 (p604).



-----
Cp
-----

Copy a p-vector.

Given:
   p        double[3]     p-vector to be copied

Returned:
   c        double[3]     copy



------
Cpv
------

Copy a position/velocity vector.

Given:
   pv     double[2][3]    position/velocity vector to be copied

Returned:
   c      double[2][3]    copy

Called:
   eraCp        copy p-vector



-----
Cr
-----

Copy an r-matrix.

Given:
   r        double[3][3]    r-matrix to be copied

Returned:
 char[]     double[3][3]    copy

Called:
   eraCp        copy p-vector



--------
D2dtf
--------

Format for output a 2-part Julian Date (or in the case of UTC a
quasi-JD form that includes special provision for leap seconds).

Given:
   scale     char[]  time scale ID (Note 1)
   ndp       int     resolution (Note 2)
   d1,d2     double  time as a 2-part Julian Date (Notes 3,4)

Returned:
   iy,im,id  int     year, month, day in Gregorian calendar (Note 5)
   ihmsf     int[4]  hours, minutes, seconds, fraction (Note 1)

Returned (function value):
             int     status: +1 = dubious year (Note 5)
                              0 = OK
                             -1 = unacceptable date (Note 6)

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
   86400 or 86401 SI seconds.

5) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale and that are too far in the future
   to be trusted.  See eraDat for further details.

6) For calendar conventions and limitations, see eraCal2jd.

Called:
   eraJd2cal    JD to Gregorian calendar
   eraD2tf      decompose days to hms
   eraDat       delta(AT) = TAI-UTC



-------
D2tf
-------

Decompose days to hours, minutes, seconds, fraction.

Given:
   ndp     int     resolution (Note 1)
   days    double  interval in days

Returned:
   sign    char    '+' or '-'
   ihmsf   int[4]  hours, minutes, seconds, fraction

Notes:

1) The argument ndp is interpreted as follows:

   ndp         resolution
    :      ...0000 00 00
   -7         1000 00 00
   -6          100 00 00
   -5           10 00 00
   -4            1 00 00
   -3            0 10 00
   -2            0 01 00
   -1            0 00 10
    0            0 00 01
    1            0 00 00.1
    2            0 00 00.01
    3            0 00 00.001
    :            0 00 00.000...

2) The largest positive useful value for ndp is determined by the
   size of days, the format of double on the target platform, and
   the risk of overflowing ihmsf[3].  On a typical platform, for
   days up to 1.0, the available floating-point precision might
   correspond to ndp=12.  However, the practical limit is typically
   ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
   only 16 bits.

3) The absolute value of days may exceed 1.0.  In cases where it
   does not, it is up to the caller to test for and handle the
   case where days is very nearly 1.0 and rounds up to 24 hours,
   by testing for ihmsf[0]=24 and setting ihmsf[0-3] to zero.



------
Dat
------

For a given UTC date, calculate delta(AT) = TAI-UTC.

   :------------------------------------------:
   :                                          :
   :                 IMPORTANT                :
   :                                          :
   :  A new version of this function must be  :
   :  produced whenever a new leap second is  :
   :  announced.  There are four items to     :
   :  change on each such occasion:           :
   :                                          :
   :  1) A new line must be added to the set  :
   :     of statements that initialize the    :
   :     array "changes".                     :
   :                                          :
   :  2) The parameter IYV must be set to     :
   :     the current year.                    :
   :                                          :
   :  3) The "Latest leap second" comment     :
   :     below must be set to the new leap    :
   :     second date.                         :
   :                                          :
   :  4) The "This revision" comment, later,  :
   :     must be set to the current date.     :
   :                                          :
   :  Change (2) must also be carried out     :
   :  whenever the function is re-issued,     :
   :  even if no leap seconds have been       :
   :  added.                                  :
   :                                          :
   :  Latest leap second:  2012 June 30       :
   :                                          :
   :__________________________________________:

Given:
   iy     int      UTC:  year (Notes 1 and 2)
   im     int            month (Note 2)
   id     int            day (Notes 2 and 3)
   fd     double         fraction of day (Note 4)

Returned:
   deltat double   TAI minus UTC, seconds

Returned (function value):
          int      status (Note 5):
                     1 = dubious year (Note 1)
                     0 = OK
                    -1 = bad year
                    -2 = bad month
                    -3 = bad day (Note 3)
                    -4 = bad fraction (Note 4)

Notes:

1) UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
   to call the function with an earlier date.  If this is attempted,
   zero is returned together with a warning status.

   Because leap seconds cannot, in principle, be predicted in
   advance, a reliable check for dates beyond the valid range is
   impossible.  To guard against gross errors, a year five or more
   after the release year of the present function (see parameter
   IYV) is considered dubious.  In this case a warning status is
   returned but the result is computed in the normal way.

   For both too-early and too-late years, the warning status is
   j=+1.  This is distinct from the error status j=-1, which
   signifies a year so early that JD could not be computed.

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
   range) even if not used;  if invalid, zero is used and status
   j=-4 is returned.  For many applications, setting fd to zero is
   acceptable;  the resulting error is always less than 3 ms (and
   occurs only pre-1972).

5) The status value returned in the case where there are multiple
   errors refers to the first error detected.  For example, if the
   month and day are 13 and 32 respectively, j=-2 (bad month)
   will be returned.

6) In cases where a valid result is not available, zero is returned.

References:

1) For dates from 1961 January 1 onwards, the expressions from the
   file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.

2) The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
   the 1992 Explanatory Supplement.

Called:
   eraCal2jd    Gregorian calendar to Julian Day number



-------
Dtdb
-------

An approximation to TDB-TT, the difference between barycentric
dynamical time and terrestrial time, for an observer on the Earth.

The different time scales - proper, coordinate and realized - are
related to each other:

          TAI             <-  physically realized
           :
        offset            <-  observed (nominally +32.184s)
           :
          TT              <-  terrestrial time
           :
  rate adjustment (L_G)   <-  definition of TT
           :
          TCG             <-  time scale for GCRS
           :
    "periodic" terms      <-  eraDtdb  is an implementation
           :
  rate adjustment (L_C)   <-  function of solar-system ephemeris
           :
          TCB             <-  time scale for BCRS
           :
  rate adjustment (-L_B)  <-  definition of TDB
           :
          TDB             <-  TCB scaled to track TT
           :
    "periodic" terms      <-  -eraDtdb is an approximation
           :
          TT              <-  terrestrial time

Adopted values for the various constants can be found in the IERS
Conventions (McCarthy & Petit 2003).

Given:
   date1,date2   double  date, TDB (Notes 1-3)
   ut            double  universal time (UT1, fraction of one day)
   elong         double  longitude (east positive, radians)
   u             double  distance from Earth spin axis (km)
   v             double  distance north of equatorial plane (km)

Returned (function value):
                 double  TDB-TT (seconds)

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



--------
Dtf2d
--------

Encode date and time fields into 2-part Julian Date (or in the case
of UTC a quasi-JD form that includes special provision for leap
seconds).

Given:
   scale     char[]  time scale ID (Note 1)
   iy,im,id  int     year, month, day in Gregorian calendar (Note 2)
   ihr,imn   int     hour, minute
   sec       double  seconds

Returned:
   d1,d2     double  2-part Julian Date (Notes 3,4)

Returned (function value):
             int     status: +3 = both of next two
                             +2 = time is after end of day (Note 5)
                             +1 = dubious year (Note 6)
                              0 = OK
                             -1 = bad year
                             -2 = bad month
                             -3 = bad day
                             -4 = bad hour
                             -5 = bad minute
                             -6 = bad second (<0)

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
   86400 or 86401 SI seconds.

5) The warning status "time is after end of day" usually means that
   the sec argument is greater than 60.0.  However, in a day ending
   in a leap second the limit changes to 61.0 (or 59.0 in the case
   of a negative leap second).

6) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale and that are too far in the future
   to be trusted.  See eraDat for further details.

7) Only in the case of continuous and regular time scales (TAI, TT,
   TCG, TCB and TDB) is the result d1+d2 a Julian Date, strictly
   speaking.  In the other cases (UT1 and UTC) the result must be
   used with circumspection;  in particular the difference between
   two such results cannot be interpreted as a precise time
   interval.

Called:
   eraCal2jd    Gregorian calendar to JD
   eraDat       delta(AT) = TAI-UTC
   eraJd2cal    JD to Gregorian calendar



-------
Ee00
-------

The equation of the equinoxes, compatible with IAU 2000 resolutions,
given the nutation in longitude and the mean obliquity.

Given:
   date1,date2  double    TT as a 2-part Julian Date (Note 1)
   epsa         double    mean obliquity (Note 2)
   dpsi         double    nutation in longitude (Note 3)

Returned (function value):
                double    equation of the equinoxes (Note 4)

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

Called:
   eraEect00    equation of the equinoxes complementary terms

References:

   Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
   implement the IAU 2000 definition of UT1", Astronomy &
   Astrophysics, 406, 1135-1149 (2003)

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



--------
Ee00a
--------

Equation of the equinoxes, compatible with IAU 2000 resolutions.

Given:
   date1,date2  double    TT as a 2-part Julian Date (Note 1)

Returned (function value):
                double    equation of the equinoxes (Note 2)

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

Called:
   eraPr00      IAU 2000 precession adjustments
   eraObl80     mean obliquity, IAU 1980
   eraNut00a    nutation, IAU 2000A
   eraEe00      equation of the equinoxes, IAU 2000

References:

   Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
   implement the IAU 2000 definition of UT1", Astronomy &
   Astrophysics, 406, 1135-1149 (2003).

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004).



--------
Ee00b
--------

Equation of the equinoxes, compatible with IAU 2000 resolutions but
using the truncated nutation model IAU 2000B.

Given:
   date1,date2  double    TT as a 2-part Julian Date (Note 1)

Returned (function value):
                double    equation of the equinoxes (Note 2)

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

Called:
   eraPr00      IAU 2000 precession adjustments
   eraObl80     mean obliquity, IAU 1980
   eraNut00b    nutation, IAU 2000B
   eraEe00      equation of the equinoxes, IAU 2000

References:

   Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
   implement the IAU 2000 definition of UT1", Astronomy &
   Astrophysics, 406, 1135-1149 (2003)

   McCarthy, D.D. & Luzum, B.J., "An abridged model of the
   precession-nutation of the celestial pole", Celestial Mechanics &
   Dynamical Astronomy, 85, 37-49 (2003)

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



--------
Ee06a
--------

Equation of the equinoxes, compatible with IAU 2000 resolutions and
IAU 2006/2000A precession-nutation.

Given:
   date1,date2  double    TT as a 2-part Julian Date (Note 1)

Returned (function value):
                double    equation of the equinoxes (Note 2)

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

Called:
   eraAnpm      normalize angle into range +/- pi
   eraGst06a    Greenwich apparent sidereal time, IAU 2006/2000A
   eraGmst06    Greenwich mean sidereal time, IAU 2006

Reference:

   McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
   IERS Technical Note No. 32, BKG



---------
Eect00
---------

Equation of the equinoxes complementary terms, consistent with
IAU 2000 resolutions.

Given:
   date1,date2  double   TT as a 2-part Julian Date (Note 1)

Returned (function value):
                double   complementary terms (Note 2)

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

Called:
   eraFal03     mean anomaly of the Moon
   eraFalp03    mean anomaly of the Sun
   eraFaf03     mean argument of the latitude of the Moon
   eraFad03     mean elongation of the Moon from the Sun
   eraFaom03    mean longitude of the Moon's ascending node
   eraFave03    mean longitude of Venus
   eraFae03     mean longitude of Earth
   eraFapa03    general accumulated precession in longitude

References:

   Capitaine, N. & Gontier, A.-M., Astron. Astrophys., 275,
   645-650 (1993)

   Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
   implement the IAU 2000 definition of UT1", Astronomy &
   Astrophysics, 406, 1135-1149 (2003)

   IAU Resolution C7, Recommendation 3 (1994)

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



--------
Eform
--------

Earth reference ellipsoids.

Given:
   n    int         ellipsoid identifier (Note 1)

Returned:
   a    double      equatorial radius (meters, Note 2)
   f    double      flattening (Note 2)

Returned (function value):
        int         status:  0 = OK
                            -1 = illegal identifier (Note 3)

Notes:

1) The identifier n is a number that specifies the choice of
   reference ellipsoid.  The following are supported:

      n    ellipsoid

      1     WGS84
      2     GRS80
      3     WGS72

   The n value has no significance outside the ERFA software.  For
   convenience, symbols WGS84 etc. are defined in erfam.h.

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



--------
Eo06a
--------

Equation of the origins, IAU 2006 precession and IAU 2000A nutation.

Given:
   date1,date2  double    TT as a 2-part Julian Date (Note 1)

Returned (function value):
                double    equation of the origins in radians

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

Called:
   eraPnm06a    classical NPB matrix, IAU 2006/2000A
   eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
   eraS06       the CIO locator s, given X,Y, IAU 2006
   eraEors      equation of the origins, Given NPB matrix and s

References:

   Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

   Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981



-------
Eors
-------

Equation of the origins, given the classical NPB matrix and the
quantity s.

Given:
   rnpb  double[3][3]  classical nutation x precession x bias matrix
   s     double        the quantity s (the CIO locator)

Returned (function value):
         double        the equation of the origins in radians.

Notes:

1)  The equation of the origins is the distance between the true
    equinox and the celestial intermediate origin and, equivalently,
    the difference between Earth rotation angle and Greenwich
    apparent sidereal time (ERA-GST).  It comprises the precession
    (since J2000.0) in right ascension plus the equation of the
    equinoxes (including the small correction terms).

2)  The algorithm is from Wallace & Capitaine (2006).

eferences:

   Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

   Wallace, P. & Capitaine, N., 2006, Astron.Astrophys. 459, 981



------
Epb
------

Julian Date to Besselian Epoch.

Given:
   dj1,dj2    double     Julian Date (see note)

Returned (function value):
              double     Besselian Epoch.

Note:

   The Julian Date is supplied in two pieces, in the usual ERFA
   manner, which is designed to preserve time resolution.  The
   Julian Date is available as a single number by adding dj1 and
   dj2.  The maximum resolution is achieved if dj1 is 2451545D0
   (J2000.0).

Reference:

   Lieske,J.H., 1979. Astron.Astrophys.,73,282.



---------
Epb2jd
---------

Besselian Epoch to Julian Date.

Given:
   epb      double    Besselian Epoch (e.g. 1957.3D0)

Returned:
   djm0     double    MJD zero-point: always 2400000.5
   djm      double    Modified Julian Date

Note:

   The Julian Date is returned in two pieces, in the usual ERFA
   manner, which is designed to preserve time resolution.  The
   Julian Date is available as a single number by adding djm0 and
   djm.

Reference:

   Lieske, J.H., 1979, Astron.Astrophys. 73, 282.



------
Epj
------

Julian Date to Julian Epoch.

Given:
   dj1,dj2    double     Julian Date (see note)

Returned (function value):
              double     Julian Epoch

Note:

   The Julian Date is supplied in two pieces, in the usual ERFA
   manner, which is designed to preserve time resolution.  The
   Julian Date is available as a single number by adding dj1 and
   dj2.  The maximum resolution is achieved if dj1 is 2451545D0
   (J2000.0).

Reference:

   Lieske, J.H., 1979, Astron.Astrophys. 73, 282.



---------
Epj2jd
---------

Julian Epoch to Julian Date.

Given:
   epj      double    Julian Epoch (e.g. 1996.8D0)

Returned:
   djm0     double    MJD zero-point: always 2400000.5
   djm      double    Modified Julian Date

Note:

   The Julian Date is returned in two pieces, in the usual ERFA
   manner, which is designed to preserve time resolution.  The
   Julian Date is available as a single number by adding djm0 and
   djm.

Reference:

   Lieske, J.H., 1979, Astron.Astrophys. 73, 282.



--------
Epv00
--------

Earth position and velocity, heliocentric and barycentric, with
respect to the Barycentric Celestial Reference System.

Given:
   date1,date2  double        TDB date (Note 1)

Returned:
   pvh          double[2][3]  heliocentric Earth position/velocity
   pvb          double[2][3]  barycentric Earth position/velocity

Returned (function value):
                int           status: 0 = OK
                                     +1 = warning: date outside
                                          the range 1900-2100 AD

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
   good compromises between resolution and convenience.  However,
   the accuracy of the result is more likely to be limited by the
   algorithm itself than the way the date has been expressed.

   n.b. TT can be used instead of TDB in most applications.

2) On return, the arrays pvh and pvb contain the following:

      pvh[0][0]  x       }
      pvh[0][1]  y       } heliocentric position, AU
      pvh[0][2]  z       }

      pvh[1][0]  xdot    }
      pvh[1][1]  ydot    } heliocentric velocity, AU/d
      pvh[1][2]  zdot    }

      pvb[0][0]  x       }
      pvb[0][1]  y       } barycentric position, AU
      pvb[0][2]  z       }

      pvb[1][0]  xdot    }
      pvb[1][1]  ydot    } barycentric velocity, AU/d
      pvb[1][2]  zdot    }

   The vectors are with respect to the Barycentric Celestial
   Reference System.  The time unit is one day in TDB.

3) The function is a SIMPLIFIED SOLUTION from the planetary theory
   VSOP2000 (X. Moisson, P. Bretagnon, 2001, Celes. Mechanics &
   Dyn. Astron., 80, 3/4, 205-213) and is an adaptation of original
   Fortran code supplied by P. Bretagnon (private comm., 2000).

4) Comparisons over the time span 1900-2100 with this simplified
   solution and the JPL DE405 ephemeris give the following results:

                              RMS    max
         Heliocentric:
            position error    3.7   11.2   km
            velocity error    1.4    5.0   mm/s

         Barycentric:
            position error    4.6   13.4   km
            velocity error    1.4    4.9   mm/s

   Comparisons with the JPL DE406 ephemeris show that by 1800 and
   2200 the position errors are approximately double their 1900-2100
   size.  By 1500 and 2500 the deterioration is a factor of 10 and
   by 1000 and 3000 a factor of 60.  The velocity accuracy falls off
   at about half that rate.

5) It is permissible to use the same array for pvh and pvb, which
   will receive the barycentric values.



---------
Eqeq94
---------

Equation of the equinoxes, IAU 1994 model.

Given:
   date1,date2   double     TDB date (Note 1)

Returned (function value):
                 double     equation of the equinoxes (Note 2)

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

Called:
   eraNut80     nutation, IAU 1980
   eraObl80     mean obliquity, IAU 1980

References:

   IAU Resolution C7, Recommendation 3 (1994).

   Capitaine, N. & Gontier, A.-M., 1993, Astron. Astrophys., 275,
   645-650.



--------
Era00
--------

Earth rotation angle (IAU 2000 model).

Given:
   dj1,dj2   double    UT1 as a 2-part Julian Date (see note)

Returned (function value):
             double    Earth rotation angle (radians), range 0-2pi

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

Called:
   eraAnp       normalize angle into range 0 to 2pi

References:

   Capitaine N., Guinot B. and McCarthy D.D, 2000, Astron.
   Astrophys., 355, 398-405.

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



--------
Fad03
--------

Fundamental argument, IERS Conventions (2003):
mean elongation of the Moon from the Sun.

Given:
   t     double    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
         double    D, radians (Note 2)

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



--------
Fae03
--------

Fundamental argument, IERS Conventions (2003):
mean longitude of Earth.

Given:
   t     double    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
         double    mean longitude of Earth, radians (Note 2)

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



--------
Faf03
--------

Fundamental argument, IERS Conventions (2003):
mean longitude of the Moon minus mean longitude of the ascending
node.

Given:
   t     double    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
         double    F, radians (Note 2)

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



---------
Faju03
---------

Fundamental argument, IERS Conventions (2003):
mean longitude of Jupiter.

Given:
   t     double    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
         double    mean longitude of Jupiter, radians (Note 2)

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



--------
Fal03
--------

Fundamental argument, IERS Conventions (2003):
mean anomaly of the Moon.

Given:
   t     double    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
         double    l, radians (Note 2)

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



---------
Falp03
---------

Fundamental argument, IERS Conventions (2003):
mean anomaly of the Sun.

Given:
   t     double    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
         double    l', radians (Note 2)

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



---------
Fama03
---------

Fundamental argument, IERS Conventions (2003):
mean longitude of Mars.

Given:
   t     double    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
         double    mean longitude of Mars, radians (Note 2)

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



---------
Fame03
---------

Fundamental argument, IERS Conventions (2003):
mean longitude of Mercury.

Given:
   t     double    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
         double    mean longitude of Mercury, radians (Note 2)

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



---------
Fane03
---------

Fundamental argument, IERS Conventions (2003):
mean longitude of Neptune.

Given:
   t     double    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
         double    mean longitude of Neptune, radians (Note 2)

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



---------
Faom03
---------

Fundamental argument, IERS Conventions (2003):
mean longitude of the Moon's ascending node.

Given:
   t     double    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
         double    Omega, radians (Note 2)

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



---------
Fapa03
---------

Fundamental argument, IERS Conventions (2003):
general accumulated precession in longitude.

Given:
   t     double    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
         double    general precession in longitude, radians (Note 2)

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



---------
Fasa03
---------

Fundamental argument, IERS Conventions (2003):
mean longitude of Saturn.

Given:
   t     double    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
         double    mean longitude of Saturn, radians (Note 2)

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



---------
Faur03
---------

Fundamental argument, IERS Conventions (2003):
mean longitude of Uranus.

Given:
   t     double    TDB, Julian centuries since J2000.0 (Note 1)

Returned  (function value):
         double    mean longitude of Uranus, radians (Note 2)

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



---------
Fave03
---------

Fundamental argument, IERS Conventions (2003):
mean longitude of Venus.

Given:
   t     double    TDB, Julian centuries since J2000.0 (Note 1)

Returned (function value):
         double    mean longitude of Venus, radians (Note 2)

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



--------
Fk52h
--------

Transform FK5 (J2000.0) star data into the Hipparcos system.

Given (all FK5, equinox J2000.0, epoch J2000.0):
   r5      double    RA (radians)
   d5      double    Dec (radians)
   dr5     double    proper motion in RA (dRA/dt, rad/Jyear)
   dd5     double    proper motion in Dec (dDec/dt, rad/Jyear)
   px5     double    parallax (arcsec)
   rv5     double    radial velocity (km/s, positive = receding)

Returned (all Hipparcos, epoch J2000.0):
   rh      double    RA (radians)
   dh      double    Dec (radians)
   drh     double    proper motion in RA (dRA/dt, rad/Jyear)
   ddh     double    proper motion in Dec (dDec/dt, rad/Jyear)
   pxh     double    parallax (arcsec)
   rvh     double    radial velocity (km/s, positive = receding)

Notes:

1) This function transforms FK5 star positions and proper motions
   into the system of the Hipparcos catalog.

2) The proper motions in RA are dRA/dt rather than
   cos(Dec)*dRA/dt, and are per year rather than per century.

3) The FK5 to Hipparcos transformation is modeled as a pure
   rotation and spin;  zonal errors in the FK5 catalog are not
   taken into account.

4) See also eraH2fk5, eraFk5hz, eraHfk5z.

Called:
   eraStarpv    star catalog data to space motion pv-vector
   eraFk5hip    FK5 to Hipparcos rotation and spin
   eraRxp       product of r-matrix and p-vector
   eraPxp       vector product of two p-vectors
   eraPpp       p-vector plus p-vector
   eraPvstar    space motion pv-vector to star catalog data

Reference:

   F.Mignard & M.Froeschle, Astron. Astrophys. 354, 732-739 (2000).



---------
Fk5hip
---------

FK5 to Hipparcos rotation and spin.

Returned:
   r5h   double[3][3]  r-matrix: FK5 rotation wrt Hipparcos (Note 2)
   s5h   double[3]     r-vector: FK5 spin wrt Hipparcos (Note 3)

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

Called:
   eraRv2m      r-vector to r-matrix

Reference:

   F.Mignard & M.Froeschle, Astron. Astrophys. 354, 732-739 (2000).



--------
Fk5hz
--------

Transform an FK5 (J2000.0) star position into the system of the
Hipparcos catalogue, assuming zero Hipparcos proper motion.

Given:
   r5           double   FK5 RA (radians), equinox J2000.0, at date
   d5           double   FK5 Dec (radians), equinox J2000.0, at date
   date1,date2  double   TDB date (Notes 1,2)

Returned:
   rh           double   Hipparcos RA (radians)
   dh           double   Hipparcos Dec (radians)

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

Called:
   eraS2c       spherical coordinates to unit vector
   eraFk5hip    FK5 to Hipparcos rotation and spin
   eraSxp       multiply p-vector by scalar
   eraRv2m      r-vector to r-matrix
   eraTrxp      product of transpose of r-matrix and p-vector
   eraPxp       vector product of two p-vectors
   eraC2s       p-vector to spherical
   eraAnp       normalize angle into range 0 to 2pi

Reference:

   F.Mignard & M.Froeschle, 2000, Astron.Astrophys. 354, 732-739.



-------
Fw2m
-------

Form rotation matrix given the Fukushima-Williams angles.

Given:
   gamb     double         F-W angle gamma_bar (radians)
   phib     double         F-W angle phi_bar (radians)
   psi      double         F-W angle psi (radians)
   eps      double         F-W angle epsilon (radians)

Returned:
   r        double[3][3]   rotation matrix

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

Called:
   eraIr        initialize r-matrix to identity
   eraRz        rotate around Z-axis
   eraRx        rotate around X-axis

Reference:

   Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351



--------
Fw2xy
--------

CIP X,Y given Fukushima-Williams bias-precession-nutation angles.

Given:
   gamb     double    F-W angle gamma_bar (radians)
   phib     double    F-W angle phi_bar (radians)
   psi      double    F-W angle psi (radians)
   eps      double    F-W angle epsilon (radians)

Returned:
   x,y      double    CIP X,Y ("radians")

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

   X,Y are elements (3,1) and (3,2) of the matrix.

Called:
   eraFw2m      F-W angles to r-matrix
   eraBpn2xy    extract CIP X,Y coordinates from NPB matrix

Reference:

   Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351



--------
Gc2gd
--------

Transform geocentric coordinates to geodetic using the specified
reference ellipsoid.

Given:
   n       int        ellipsoid identifier (Note 1)
   xyz     double[3]  geocentric vector (Note 2)

Returned:
   elong   double     longitude (radians, east +ve)
   phi     double     latitude (geodetic, radians, Note 3)
   height  double     height above ellipsoid (geodetic, Notes 2,3)

Returned (function value):
          int         status:  0 = OK
                              -1 = illegal identifier (Note 3)
                              -2 = internal error (Note 3)

Notes:

1) The identifier n is a number that specifies the choice of
   reference ellipsoid.  The following are supported:

      n    ellipsoid

      1     WGS84
      2     GRS80
      3     WGS72

   The n value has no significance outside the ERFA software.  For
   convenience, symbols WGS84 etc. are defined in erfam.h.

2) The geocentric vector (xyz, given) and height (height, returned)
   are in meters.

3) An error status -1 means that the identifier n is illegal.  An
   error status -2 is theoretically impossible.  In all error cases,
   phi and height are both set to -1e9.

4) The inverse transformation is performed in the function eraGd2gc.

Called:
   eraEform     Earth reference ellipsoids
   eraGc2gde    geocentric to geodetic transformation, general



---------
Gc2gde
---------

Transform geocentric coordinates to geodetic for a reference
ellipsoid of specified form.

Given:
   a       double     equatorial radius (Notes 2,4)
   f       double     flattening (Note 3)
   xyz     double[3]  geocentric vector (Note 4)

Returned:
   elong   double     longitude (radians, east +ve)
   phi     double     latitude (geodetic, radians)
   height  double     height above ellipsoid (geodetic, Note 4)

Returned (function value):
           int        status:  0 = OK
                              -1 = illegal f
                              -2 = illegal a

Notes:

1) This function is based on the GCONV2H Fortran subroutine by
   Toshio Fukushima (see reference).

2) The equatorial radius, a, can be in any units, but meters is
   the conventional choice.

3) The flattening, f, is (for the Earth) a value around 0.00335,
   i.e. around 1/298.

4) The equatorial radius, a, and the geocentric vector, xyz,
   must be given in the same units, and determine the units of
   the returned height, height.

5) If an error occurs (status < 0), elong, phi and height are
   unchanged.

6) The inverse transformation is performed in the function
   eraGd2gce.

7) The transformation for a standard ellipsoid (such as WGS84) can
   more conveniently be performed by calling eraGc2gd, which uses a
   numerical code to identify the required A and F values.

Reference:

   Fukushima, T., "Transformation from Cartesian to geodetic
   coordinates accelerated by Halley's method", J.Geodesy (2006)
   79: 689-693



--------
Gd2gc
--------

Transform geodetic coordinates to geocentric using the specified
reference ellipsoid.

Given:
   n       int        ellipsoid identifier (Note 1)
   elong   double     longitude (radians, east +ve)
   phi     double     latitude (geodetic, radians, Note 3)
   height  double     height above ellipsoid (geodetic, Notes 2,3)

Returned:
   xyz     double[3]  geocentric vector (Note 2)

Returned (function value):
           int        status:  0 = OK
                              -1 = illegal identifier (Note 3)
                              -2 = illegal case (Note 3)

Notes:

1) The identifier n is a number that specifies the choice of
   reference ellipsoid.  The following are supported:

      n    ellipsoid

      1     WGS84
      2     GRS80
      3     WGS72

   The n value has no significance outside the ERFA software.  For
   convenience, symbols WGS84 etc. are defined in erfam.h.

2) The height (height, given) and the geocentric vector (xyz,
   returned) are in meters.

3) No validation is performed on the arguments elong, phi and
   height.  An error status -1 means that the identifier n is
   illegal.  An error status -2 protects against cases that would
   lead to arithmetic exceptions.  In all error cases, xyz is set
   to zeros.

4) The inverse transformation is performed in the function eraGc2gd.

Called:
   eraEform     Earth reference ellipsoids
   eraGd2gce    geodetic to geocentric transformation, general
   eraZp        zero p-vector



---------
Gd2gce
---------

Transform geodetic coordinates to geocentric for a reference
ellipsoid of specified form.

Given:
   a       double     equatorial radius (Notes 1,4)
   f       double     flattening (Notes 2,4)
   elong   double     longitude (radians, east +ve)
   phi     double     latitude (geodetic, radians, Note 4)
   height  double     height above ellipsoid (geodetic, Notes 3,4)

Returned:
   xyz     double[3]  geocentric vector (Note 3)

Returned (function value):
           int        status:  0 = OK
                              -1 = illegal case (Note 4)
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

6) The transformation for a standard ellipsoid (such as WGS84) can
   more conveniently be performed by calling eraGd2gc,  which uses a
   numerical code to identify the required a and f values.

References:

   Green, R.M., Spherical Astronomy, Cambridge University Press,
   (1985) Section 4.5, p96.

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992),
   Section 4.22, p202.



---------
Gmst00
---------

Greenwich mean sidereal time (model consistent with IAU 2000
resolutions).

Given:
   uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
   tta,ttb    double    TT as a 2-part Julian Date (Notes 1,2)

Returned (function value):
              double    Greenwich mean sidereal time (radians)

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

Called:
   eraEra00     Earth rotation angle, IAU 2000
   eraAnp       normalize angle into range 0 to 2pi

References:

   Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
   implement the IAU 2000 definition of UT1", Astronomy &
   Astrophysics, 406, 1135-1149 (2003)

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



---------
Gmst06
---------

Greenwich mean sidereal time (consistent with IAU 2006 precession).

Given:
   uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
   tta,ttb    double    TT as a 2-part Julian Date (Notes 1,2)

Returned (function value):
              double    Greenwich mean sidereal time (radians)

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

Called:
   eraEra00     Earth rotation angle, IAU 2000
   eraAnp       normalize angle into range 0 to 2pi

Reference:

   Capitaine, N., Wallace, P.T. & Chapront, J., 2005,
   Astron.Astrophys. 432, 355



---------
Gmst82
---------

Universal Time to Greenwich mean sidereal time (IAU 1982 model).

Given:
   dj1,dj2    double    UT1 Julian Date (see note)

Returned (function value):
              double    Greenwich mean sidereal time (radians)

Notes:

1) The UT1 date dj1+dj2 is a Julian Date, apportioned in any
   convenient way between the arguments dj1 and dj2.  For example,
   JD(UT1)=2450123.7 could be expressed in any of these ways,
   among others:

           dj1            dj2

       2450123.7D0        0D0        (JD method)
        2451545D0      -1421.3D0     (J2000 method)
       2400000.5D0     50123.2D0     (MJD method)
       2450123.5D0       0.2D0       (date & time method)

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

Called:
   eraAnp       normalize angle into range 0 to 2pi

References:

   Transactions of the International Astronomical Union,
   XVIII B, 67 (1983).

   Aoki et al., Astron. Astrophys. 105, 359-361 (1982).



---------
Gst00a
---------

Greenwich apparent sidereal time (consistent with IAU 2000
resolutions).

Given:
   uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
   tta,ttb    double    TT as a 2-part Julian Date (Notes 1,2)

Returned (function value):
              double    Greenwich apparent sidereal time (radians)

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

Called:
   eraGmst00    Greenwich mean sidereal time, IAU 2000
   eraEe00a     equation of the equinoxes, IAU 2000A
   eraAnp       normalize angle into range 0 to 2pi

References:

   Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
   implement the IAU 2000 definition of UT1", Astronomy &
   Astrophysics, 406, 1135-1149 (2003)

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



---------
Gst00b
---------

Greenwich apparent sidereal time (consistent with IAU 2000
resolutions but using the truncated nutation model IAU 2000B).

Given:
   uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)

Returned (function value):
              double    Greenwich apparent sidereal time (radians)

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

Called:
   eraGmst00    Greenwich mean sidereal time, IAU 2000
   eraEe00b     equation of the equinoxes, IAU 2000B
   eraAnp       normalize angle into range 0 to 2pi

References:

   Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
   implement the IAU 2000 definition of UT1", Astronomy &
   Astrophysics, 406, 1135-1149 (2003)

   McCarthy, D.D. & Luzum, B.J., "An abridged model of the
   precession-nutation of the celestial pole", Celestial Mechanics &
   Dynamical Astronomy, 85, 37-49 (2003)

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



--------
Gst06
--------

Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.

Given:
   uta,utb  double        UT1 as a 2-part Julian Date (Notes 1,2)
   tta,ttb  double        TT as a 2-part Julian Date (Notes 1,2)
   rnpb     double[3][3]  nutation x precession x bias matrix

Returned (function value):
            double        Greenwich apparent sidereal time (radians)

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

Called:
   eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
   eraS06       the CIO locator s, given X,Y, IAU 2006
   eraAnp       normalize angle into range 0 to 2pi
   eraEra00     Earth rotation angle, IAU 2000
   eraEors      equation of the origins, given NPB matrix and s

Reference:

   Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981



---------
Gst06a
---------

Greenwich apparent sidereal time (consistent with IAU 2000 and 2006
resolutions).

Given:
   uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
   tta,ttb    double    TT as a 2-part Julian Date (Notes 1,2)

Returned (function value):
              double    Greenwich apparent sidereal time (radians)

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

Called:
   eraPnm06a    classical NPB matrix, IAU 2006/2000A
   eraGst06     Greenwich apparent ST, IAU 2006, given NPB matrix

Reference:

   Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981



--------
Gst94
--------

Greenwich apparent sidereal time (consistent with IAU 1982/94
resolutions).

Given:
   uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)

Returned (function value):
              double    Greenwich apparent sidereal time (radians)

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

Called:
   eraGmst82    Greenwich mean sidereal time, IAU 1982
   eraEqeq94    equation of the equinoxes, IAU 1994
   eraAnp       normalize angle into range 0 to 2pi

References:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992)

   IAU Resolution C7, Recommendation 3 (1994)



--------
H2fk5
--------

Transform Hipparcos star data into the FK5 (J2000.0) system.

Given (all Hipparcos, epoch J2000.0):
   rh      double    RA (radians)
   dh      double    Dec (radians)
   drh     double    proper motion in RA (dRA/dt, rad/Jyear)
   ddh     double    proper motion in Dec (dDec/dt, rad/Jyear)
   pxh     double    parallax (arcsec)
   rvh     double    radial velocity (km/s, positive = receding)

Returned (all FK5, equinox J2000.0, epoch J2000.0):
   r5      double    RA (radians)
   d5      double    Dec (radians)
   dr5     double    proper motion in RA (dRA/dt, rad/Jyear)
   dd5     double    proper motion in Dec (dDec/dt, rad/Jyear)
   px5     double    parallax (arcsec)
   rv5     double    radial velocity (km/s, positive = receding)

Notes:

1) This function transforms Hipparcos star positions and proper
   motions into FK5 J2000.0.

2) The proper motions in RA are dRA/dt rather than
   cos(Dec)*dRA/dt, and are per year rather than per century.

3) The FK5 to Hipparcos transformation is modeled as a pure
   rotation and spin;  zonal errors in the FK5 catalog are not
   taken into account.

4) See also eraFk52h, eraFk5hz, eraHfk5z.

Called:
   eraStarpv    star catalog data to space motion pv-vector
   eraFk5hip    FK5 to Hipparcos rotation and spin
   eraRv2m      r-vector to r-matrix
   eraRxp       product of r-matrix and p-vector
   eraTrxp      product of transpose of r-matrix and p-vector
   eraPxp       vector product of two p-vectors
   eraPmp       p-vector minus p-vector
   eraPvstar    space motion pv-vector to star catalog data

Reference:

   F.Mignard & M.Froeschle, Astron. Astrophys. 354, 732-739 (2000).



--------
Hfk5z
--------

Transform a Hipparcos star position into FK5 J2000.0, assuming
zero Hipparcos proper motion.

Given:
   rh            double    Hipparcos RA (radians)
   dh            double    Hipparcos Dec (radians)
   date1,date2   double    TDB date (Note 1)

Returned (all FK5, equinox J2000.0, date date1+date2):
   r5            double    RA (radians)
   d5            double    Dec (radians)
   dr5           double    FK5 RA proper motion (rad/year, Note 4)
   dd5           double    Dec proper motion (rad/year, Note 4)

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

Called:
   eraS2c       spherical coordinates to unit vector
   eraFk5hip    FK5 to Hipparcos rotation and spin
   eraRxp       product of r-matrix and p-vector
   eraSxp       multiply p-vector by scalar
   eraRxr       product of two r-matrices
   eraTrxp      product of transpose of r-matrix and p-vector
   eraPxp       vector product of two p-vectors
   eraPv2s      pv-vector to spherical
   eraAnp       normalize angle into range 0 to 2pi

Reference:

   F.Mignard & M.Froeschle, 2000, Astron.Astrophys. 354, 732-739.



-----
Ir
-----

Initialize an r-matrix to the identity matrix.

Returned:
   r       double[3][3]    r-matrix



---------
Jd2cal
---------

Julian Date to Gregorian year, month, day, and fraction of a day.

Given:
   dj1,dj2   double   Julian Date (Notes 1, 2)

Returned (arguments):
   iy        int      year
   im        int      month
   id        int      day
   fd        double   fraction of day

Returned (function value):
             int      status:
                         0 = OK
                        -1 = unacceptable date (Note 3)

Notes:

1) The earliest valid date is -68569.5 (-4900 March 1).  The
   largest value accepted is 10^9.

2) The Julian Date is apportioned in any convenient way between
   the arguments dj1 and dj2.  For example, JD=2450123.7 could
   be expressed in any of these ways, among others:

          dj1             dj2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

3) In early eras the conversion is from the "proleptic Gregorian
   calendar";  no account is taken of the date(s) of adoption of
   the Gregorian calendar, nor is the AD/BC numbering convention
   observed.

Reference:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992),
   Section 12.92 (p604).



---------
Jdcalf
---------

Julian Date to Gregorian Calendar, expressed in a form convenient
for formatting messages:  rounded to a specified precision.

Given:
   ndp       int      number of decimal places of days in fraction
   dj1,dj2   double   dj1+dj2 = Julian Date (Note 1)

Returned:
   iymdf     int[4]   year, month, day, fraction in Gregorian
                      calendar

Returned (function value):
             int      status:
                        -1 = date out of range
                         0 = OK
                        +1 = NDP not 0-9 (interpreted as 0)

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

Called:
   eraJd2cal    JD to Gregorian calendar

Reference:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992),
   Section 12.92 (p604).



---------
Num00a
---------

Form the matrix of nutation for a given date, IAU 2000A model.

Given:
   date1,date2  double          TT as a 2-part Julian Date (Note 1)

Returned:
   rmatn        double[3][3]    nutation matrix

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

Called:
   eraPn00a     bias/precession/nutation, IAU 2000A

Reference:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992),
   Section 3.222-3 (p114).



---------
Num00b
---------

Form the matrix of nutation for a given date, IAU 2000B model.

Given:
   date1,date2  double         TT as a 2-part Julian Date (Note 1)

Returned:
   rmatn        double[3][3]   nutation matrix

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

Called:
   eraPn00b     bias/precession/nutation, IAU 2000B

Reference:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992),
   Section 3.222-3 (p114).



---------
Num06a
---------

Form the matrix of nutation for a given date, IAU 2006/2000A model.

Given:
   date1,date2   double          TT as a 2-part Julian Date (Note 1)

Returned:
   rmatn         double[3][3]    nutation matrix

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

Called:
   eraObl06     mean obliquity, IAU 2006
   eraNut06a    nutation, IAU 2006/2000A
   eraNumat     form nutation matrix

Reference:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992),
   Section 3.222-3 (p114).



--------
Numat
--------

Form the matrix of nutation.

Given:
   epsa        double         mean obliquity of date (Note 1)
   dpsi,deps   double         nutation (Note 2)

Returned:
   rmatn       double[3][3]   nutation matrix (Note 3)

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

Called:
   eraIr        initialize r-matrix to identity
   eraRx        rotate around X-axis
   eraRz        rotate around Z-axis

Reference:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992),
   Section 3.222-3 (p114).



---------
Nut00a
---------

Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation
with free core nutation omitted).

Given:
   date1,date2   double   TT as a 2-part Julian Date (Note 1)

Returned:
   dpsi,deps     double   nutation, luni-solar + planetary (Note 2)

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

Called:
   eraFal03     mean anomaly of the Moon
   eraFaf03     mean argument of the latitude of the Moon
   eraFaom03    mean longitude of the Moon's ascending node
   eraFame03    mean longitude of Mercury
   eraFave03    mean longitude of Venus
   eraFae03     mean longitude of Earth
   eraFama03    mean longitude of Mars
   eraFaju03    mean longitude of Jupiter
   eraFasa03    mean longitude of Saturn
   eraFaur03    mean longitude of Uranus
   eraFapa03    general accumulated precession in longitude

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



---------
Nut00b
---------

Nutation, IAU 2000B model.

Given:
   date1,date2   double    TT as a 2-part Julian Date (Note 1)

Returned:
   dpsi,deps     double    nutation, luni-solar + planetary (Note 2)

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



---------
Nut06a
---------

IAU 2000A nutation with adjustments to match the IAU 2006
precession.

Given:
   date1,date2   double   TT as a 2-part Julian Date (Note 1)

Returned:
   dpsi,deps     double   nutation, luni-solar + planetary (Note 2)

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

Called:
   eraNut00a    nutation, IAU 2000A

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



--------
Nut80
--------

Nutation, IAU 1980 model.

Given:
   date1,date2   double    TT as a 2-part Julian Date (Note 1)

Returned:
   dpsi          double    nutation in longitude (radians)
   deps          double    nutation in obliquity (radians)

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

Called:
   eraAnpm      normalize angle into range +/- pi

Reference:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992),
   Section 3.222 (p111).



---------
Nutm80
---------

Form the matrix of nutation for a given date, IAU 1980 model.

Given:
   date1,date2    double          TDB date (Note 1)

Returned:
   rmatn          double[3][3]    nutation matrix

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

Called:
   eraNut80     nutation, IAU 1980
   eraObl80     mean obliquity, IAU 1980
   eraNumat     form nutation matrix



--------
Obl06
--------

Mean obliquity of the ecliptic, IAU 2006 precession model.

Given:
   date1,date2  double   TT as a 2-part Julian Date (Note 1)

Returned (function value):
                double   obliquity of the ecliptic (radians, Note 2)

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

2) The result is the angle between the ecliptic and mean equator of
   date date1+date2.

Reference:

   Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351



--------
Obl80
--------

Mean obliquity of the ecliptic, IAU 1980 model.

Given:
   date1,date2   double    TT as a 2-part Julian Date (Note 1)

Returned (function value):
                 double    obliquity of the ecliptic (radians, Note 2)

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

2) The result is the angle between the ecliptic and mean equator of
   date date1+date2.

Reference:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992),
   Expression 3.222-1 (p114).



-------
P06e
-------

Precession angles, IAU 2006, equinox based.

Given:
   date1,date2   double   TT as a 2-part Julian Date (Note 1)

Returned (see Note 2):
   eps0          double   epsilon_0
   psia          double   psi_A
   oma           double   omega_A
   bpa           double   P_A
   bqa           double   Q_A
   pia           double   pi_A
   bpia          double   Pi_A
   epsa          double   obliquity epsilon_A
   chia          double   chi_A
   za            double   z_A
   zetaa         double   zeta_A
   thetaa        double   theta_A
   pa            double   p_A
   gam           double   F-W angle gamma_J2000
   phi           double   F-W angle phi_J2000
   psi           double   F-W angle psi_J2000

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

Called:
   eraObl06     mean obliquity, IAU 2006



-------
P2pv
-------

Extend a p-vector to a pv-vector by appending a zero velocity.

Given:
   p        double[3]       p-vector

Returned:
   pv       double[2][3]    pv-vector

Called:
   eraCp        copy p-vector
   eraZp        zero p-vector



------
P2s
------

P-vector to spherical polar coordinates.

Given:
   p        double[3]    p-vector

Returned:
   theta    double       longitude angle (radians)
   phi      double       latitude angle (radians)
   r        double       radial distance

Notes:

1) If P is null, zero theta, phi and r are returned.

2) At either pole, zero theta is returned.

Called:
   eraC2s       p-vector to spherical
   eraPm        modulus of p-vector



------
Pap
------

Position-angle from two p-vectors.

Given:
   a      double[3]  direction of reference point
   b      double[3]  direction of point whose PA is required

Returned (function value):
          double     position angle of b with respect to a (radians)

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

Called:
   eraPn        decompose p-vector into modulus and direction
   eraPm        modulus of p-vector
   eraPxp       vector product of two p-vectors
   eraPmp       p-vector minus p-vector
   eraPdp       scalar product of two p-vectors



------
Pas
------

Position-angle from spherical coordinates.

Given:
   al     double     longitude of point A (e.g. RA) in radians
   ap     double     latitude of point A (e.g. Dec) in radians
   bl     double     longitude of point B
   bp     double     latitude of point B

Returned (function value):
          double     position angle of B with respect to A

Notes:

1) The result is the bearing (position angle), in radians, of point
   B with respect to point A.  It is in the range -pi to +pi.  The
   sense is such that if B is a small distance "east" of point A,
   the bearing is approximately +pi/2.

2) Zero is returned if the two points are coincident.



-------
Pb06
-------

This function forms three Euler angles which implement general
precession from epoch J2000.0, using the IAU 2006 model.  Frame
bias (the offset between ICRS and mean J2000.0) is included.

Given:
   date1,date2  double   TT as a 2-part Julian Date (Note 1)

Returned:
   bzeta        double   1st rotation: radians cw around z
   bz           double   3rd rotation: radians cw around z
   btheta       double   2nd rotation: radians ccw around y

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

Called:
   eraPmat06    PB matrix, IAU 2006
   eraRz        rotate around Z-axis



------
Pdp
------

p-vector inner (=scalar=dot) product.

Given:
   a      double[3]     first p-vector
   b      double[3]     second p-vector

Returned (function value):
          double        a . b



--------
Pfw06
--------

Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).

Given:
   date1,date2  double   TT as a 2-part Julian Date (Note 1)

Returned:
   gamb         double   F-W angle gamma_bar (radians)
   phib         double   F-W angle phi_bar (radians)
   psib         double   F-W angle psi_bar (radians)
   epsa         double   F-W angle epsilon_A (radians)

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

Called:
   eraObl06     mean obliquity, IAU 2006



---------
Plan94
---------

Approximate heliocentric position and velocity of a nominated major
planet:  Mercury, Venus, EMB, Mars, Jupiter, Saturn, Uranus or
Neptune (but not the Earth itself).

Given:
   date1  double       TDB date part A (Note 1)
   date2  double       TDB date part B (Note 1)
   np     int          planet (1=Mercury, 2=Venus, 3=EMB, 4=Mars,
                           5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune)

Returned (argument):
   pv     double[2][3] planet p,v (heliocentric, J2000.0, AU,AU/d)

Returned (function value):
          int          status: -1 = illegal NP (outside 1-8)
                                0 = OK
                               +1 = warning: year outside 1000-3000
                               +2 = warning: failed to converge

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

Called:
   eraAnp       normalize angle into range 0 to 2pi

Reference:  Simon, J.L, Bretagnon, P., Chapront, J.,
            Chapront-Touze, M., Francou, G., and Laskar, J.,
            Astron. Astrophys. 282, 663 (1994).



-----
Pm
-----

Modulus of p-vector.

Given:
   p      double[3]     p-vector

Returned (function value):
          double        modulus



---------
Pmat00
---------

Precession matrix (including frame bias) from GCRS to a specified
date, IAU 2000 model.

Given:
   date1,date2  double          TT as a 2-part Julian Date (Note 1)

Returned:
   rbp          double[3][3]    bias-precession matrix (Note 2)

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

Called:
   eraBp00      frame bias and precession matrices, IAU 2000

Reference:

   IAU: Trans. International Astronomical Union, Vol. XXIVB;  Proc.
   24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
   (2000)



---------
Pmat06
---------

Precession matrix (including frame bias) from GCRS to a specified
date, IAU 2006 model.

Given:
   date1,date2  double          TT as a 2-part Julian Date (Note 1)

Returned:
   rbp          double[3][3]    bias-precession matrix (Note 2)

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

Called:
   eraPfw06     bias-precession F-W angles, IAU 2006
   eraFw2m      F-W angles to r-matrix

References:

   Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

   Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981



---------
Pmat76
---------

Precession matrix from J2000.0 to a specified date, IAU 1976 model.

Given:
   date1,date2 double       ending date, TT (Note 1)

Returned:
   rmatp       double[3][3] precession matrix, J2000.0 -> date1+date2

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

Called:
   eraPrec76    accumulated precession angles, IAU 1976
   eraIr        initialize r-matrix to identity
   eraRz        rotate around Z-axis
   eraRy        rotate around Y-axis
   eraCr        copy r-matrix

References:

   Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
    equations (6) & (7), p283.

   Kaplan,G.H., 1981. USNO circular no. 163, pA2.



------
Pmp
------

P-vector subtraction.

Given:
   a        double[3]      first p-vector
   b        double[3]      second p-vector

Returned:
   amb      double[3]      a - b

Note:
   It is permissible to re-use the same array for any of the
   arguments.



-----
Pn
-----

Convert a p-vector into modulus and unit vector.

Given:
   p        double[3]      p-vector

Returned:
   r        double         modulus
   u        double[3]      unit vector

Notes:

1) If p is null, the result is null.  Otherwise the result is a unit
   vector.

2) It is permissible to re-use the same array for any of the
   arguments.

Called:
   eraPm        modulus of p-vector
   eraZp        zero p-vector
   eraSxp       multiply p-vector by scalar



-------
Pn00
-------

Precession-nutation, IAU 2000 model:  a multi-purpose function,
supporting classical (equinox-based) use directly and CIO-based
use indirectly.

Given:
   date1,date2  double          TT as a 2-part Julian Date (Note 1)
   dpsi,deps    double          nutation (Note 2)

Returned:
   epsa         double          mean obliquity (Note 3)
   rb           double[3][3]    frame bias matrix (Note 4)
   rp           double[3][3]    precession matrix (Note 5)
   rbp          double[3][3]    bias-precession matrix (Note 6)
   rn           double[3][3]    nutation matrix (Note 7)
   rbpn         double[3][3]    GCRS-to-true matrix (Note 8)

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

Called:
   eraPr00      IAU 2000 precession adjustments
   eraObl80     mean obliquity, IAU 1980
   eraBp00      frame bias and precession matrices, IAU 2000
   eraCr        copy r-matrix
   eraNumat     form nutation matrix
   eraRxr       product of two r-matrices

Reference:

   Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
   "Expressions for the Celestial Intermediate Pole and Celestial
   Ephemeris Origin consistent with the IAU 2000A precession-
   nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

   n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.



--------
Pn00a
--------

Precession-nutation, IAU 2000A model:  a multi-purpose function,
supporting classical (equinox-based) use directly and CIO-based
use indirectly.

Given:
   date1,date2  double          TT as a 2-part Julian Date (Note 1)

Returned:
   dpsi,deps    double          nutation (Note 2)
   epsa         double          mean obliquity (Note 3)
   rb           double[3][3]    frame bias matrix (Note 4)
   rp           double[3][3]    precession matrix (Note 5)
   rbp          double[3][3]    bias-precession matrix (Note 6)
   rn           double[3][3]    nutation matrix (Note 7)
   rbpn         double[3][3]    GCRS-to-true matrix (Notes 8,9)

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

9)  The X,Y,Z coordinates of the IAU 2000B Celestial Intermediate
    Pole are elements (3,1-3) of the matrix rbpn.

10) It is permissible to re-use the same array in the returned
    arguments.  The arrays are filled in the order given.

Called:
   eraNut00a    nutation, IAU 2000A
   eraPn00      bias/precession/nutation results, IAU 2000

Reference:

   Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
   "Expressions for the Celestial Intermediate Pole and Celestial
   Ephemeris Origin consistent with the IAU 2000A precession-
   nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

   n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.



--------
Pn00b
--------

Precession-nutation, IAU 2000B model:  a multi-purpose function,
supporting classical (equinox-based) use directly and CIO-based
use indirectly.

Given:
   date1,date2  double          TT as a 2-part Julian Date (Note 1)

Returned:
   dpsi,deps    double          nutation (Note 2)
   epsa         double          mean obliquity (Note 3)
   rb           double[3][3]    frame bias matrix (Note 4)
   rp           double[3][3]    precession matrix (Note 5)
   rbp          double[3][3]    bias-precession matrix (Note 6)
   rn           double[3][3]    nutation matrix (Note 7)
   rbpn         double[3][3]    GCRS-to-true matrix (Notes 8,9)

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
    Pole are elements (3,1-3) of the matrix rbpn.

10) It is permissible to re-use the same array in the returned
    arguments.  The arrays are filled in the stated order.

Called:
   eraNut00b    nutation, IAU 2000B
   eraPn00      bias/precession/nutation results, IAU 2000

Reference:

   Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
   "Expressions for the Celestial Intermediate Pole and Celestial
   Ephemeris Origin consistent with the IAU 2000A precession-
   nutation model", Astron.Astrophys. 400, 1145-1154 (2003).

   n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.



-------
Pn06
-------

Precession-nutation, IAU 2006 model:  a multi-purpose function,
supporting classical (equinox-based) use directly and CIO-based use
indirectly.

Given:
   date1,date2  double          TT as a 2-part Julian Date (Note 1)
   dpsi,deps    double          nutation (Note 2)

Returned:
   epsa         double          mean obliquity (Note 3)
   rb           double[3][3]    frame bias matrix (Note 4)
   rp           double[3][3]    precession matrix (Note 5)
   rbp          double[3][3]    bias-precession matrix (Note 6)
   rn           double[3][3]    nutation matrix (Note 7)
   rbpn         double[3][3]    GCRS-to-true matrix (Note 8)

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

9)  The X,Y,Z coordinates of the IAU 2000B Celestial Intermediate
    Pole are elements (3,1-3) of the matrix rbpn.

10) It is permissible to re-use the same array in the returned
    arguments.  The arrays are filled in the stated order.

Called:
   eraPfw06     bias-precession F-W angles, IAU 2006
   eraFw2m      F-W angles to r-matrix
   eraCr        copy r-matrix
   eraTr        transpose r-matrix
   eraRxr       product of two r-matrices

References:

   Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

   Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981



--------
Pn06a
--------

Precession-nutation, IAU 2006/2000A models:  a multi-purpose function,
supporting classical (equinox-based) use directly and CIO-based use
indirectly.

Given:
   date1,date2  double          TT as a 2-part Julian Date (Note 1)

Returned:
   dpsi,deps    double          nutation (Note 2)
   epsa         double          mean obliquity (Note 3)
   rb           double[3][3]    frame bias matrix (Note 4)
   rp           double[3][3]    precession matrix (Note 5)
   rbp          double[3][3]    bias-precession matrix (Note 6)
   rn           double[3][3]    nutation matrix (Note 7)
   rbpn         double[3][3]    GCRS-to-true matrix (Notes 8,9)

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
    Intermediate Pole are elements (1,1-3) of the matrix rbpn.

10) It is permissible to re-use the same array in the returned
    arguments.  The arrays are filled in the stated order.

Called:
   eraNut06a    nutation, IAU 2006/2000A
   eraPn06      bias/precession/nutation results, IAU 2006

Reference:

   Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855



---------
Pnm00a
---------

Form the matrix of precession-nutation for a given date (including
frame bias), equinox-based, IAU 2000A model.

Given:
   date1,date2  double     TT as a 2-part Julian Date (Note 1)

Returned:
   rbpn         double[3][3]    classical NPB matrix (Note 2)

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

Called:
   eraPn00a     bias/precession/nutation, IAU 2000A

Reference:

   IAU: Trans. International Astronomical Union, Vol. XXIVB;  Proc.
   24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
   (2000)



---------
Pnm00b
---------

Form the matrix of precession-nutation for a given date (including
frame bias), equinox-based, IAU 2000B model.

Given:
   date1,date2 double       TT as a 2-part Julian Date (Note 1)

Returned:
   rbpn        double[3][3] bias-precession-nutation matrix (Note 2)

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

Called:
   eraPn00b     bias/precession/nutation, IAU 2000B

Reference:

   IAU: Trans. International Astronomical Union, Vol. XXIVB;  Proc.
   24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
   (2000)



---------
Pnm06a
---------

Form the matrix of precession-nutation for a given date (including
frame bias), IAU 2006 precession and IAU 2000A nutation models.

Given:
   date1,date2 double       TT as a 2-part Julian Date (Note 1)

Returned:
   rnpb        double[3][3] bias-precession-nutation matrix (Note 2)

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

Called:
   eraPfw06     bias-precession F-W angles, IAU 2006
   eraNut06a    nutation, IAU 2006/2000A
   eraFw2m      F-W angles to r-matrix

Reference:

   Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855.



--------
Pnm80
--------

Form the matrix of precession/nutation for a given date, IAU 1976
precession model, IAU 1980 nutation model.

Given:
   date1,date2    double         TDB date (Note 1)

Returned:
   rmatpn         double[3][3]   combined precession/nutation matrix

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

Called:
   eraPmat76    precession matrix, IAU 1976
   eraNutm80    nutation matrix, IAU 1980
   eraRxr       product of two r-matrices

Reference:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992),
   Section 3.3 (p145).



---------
Pom00
---------

Form the matrix of polar motion for a given date, IAU 2000.

Given:
   xp,yp    double    coordinates of the pole (radians, Note 1)
   sp       double    the TIO locator s' (radians, Note 2)

Returned:
   rpom     double[3][3]   polar-motion matrix (Note 3)

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

Called:
   eraIr        initialize r-matrix to identity
   eraRz        rotate around Z-axis
   eraRy        rotate around Y-axis
   eraRx        rotate around X-axis

Reference:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



------
Ppp
------

P-vector addition.

Given:
   a        double[3]      first p-vector
   b        double[3]      second p-vector

Returned:
   apb      double[3]      a + b

Note:
   It is permissible to re-use the same array for any of the
   arguments.



-------
Ppsp
-------

P-vector plus scaled p-vector.

Given:
   a      double[3]     first p-vector
   s      double        scalar (multiplier for b)
   b      double[3]     second p-vector

Returned:
   apsb   double[3]     a + s*b

Note:
   It is permissible for any of a, b and apsb to be the same array.

Called:
   eraSxp       multiply p-vector by scalar
   eraPpp       p-vector plus p-vector



-------
Pr00
-------

Precession-rate part of the IAU 2000 precession-nutation models
(part of MHB2000).

Given:
   date1,date2    double  TT as a 2-part Julian Date (Note 1)

Returned:
   dpsipr,depspr  double  precession corrections (Notes 2,3)

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



---------
Prec76
---------

IAU 1976 precession model.

This function forms the three Euler angles which implement general
precession between two epochs, using the IAU 1976 model (as for
the FK5 catalog).

Given:
   ep01,ep02   double    TDB starting epoch (Note 1)
   ep11,ep12   double    TDB ending epoch (Note 1)

Returned:
   zeta        double    1st rotation: radians cw around z
   z           double    3rd rotation: radians cw around z
   theta       double    2nd rotation: radians ccw around y

Notes:

1) The epochs ep01+ep02 and ep11+ep12 are Julian Dates, apportioned
   in any convenient way between the arguments epn1 and epn2.  For
   example, JD(TDB)=2450123.7 could be expressed in any of these
   ways, among others:

           epn1          epn2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 method is best matched to the way the
   argument is handled internally and will deliver the optimum
   optimum resolution.  The MJD method and the date & time methods
   are both good compromises between resolution and convenience.
   The two epochs may be expressed using different methods, but at
   the risk of losing some resolution.

2) The accumulated precession angles zeta, z, theta are expressed
   through canonical polynomials which are valid only for a limited
   time span.  In addition, the IAU 1976 precession rate is known to
   be imperfect.  The absolute accuracy of the present formulation
   is better than 0.1 arcsec from 1960AD to 2040AD, better than
   1 arcsec from 1640AD to 2360AD, and remains below 3 arcsec for
   the whole of the period 500BC to 3000AD.  The errors exceed
   10 arcsec outside the range 1200BC to 3900AD, exceed 100 arcsec
   outside 4200BC to 5600AD and exceed 1000 arcsec outside 6800BC to
   8200AD.

3) The three angles are returned in the conventional order, which
   is not the same as the order of the corresponding Euler
   rotations.  The precession matrix is
   R_3(-z) x R_2(+theta) x R_3(-zeta).

Reference:

   Lieske, J.H., 1979, Astron.Astrophys. 73, 282, equations
   (6) & (7), p283.



-------
Pv2p
-------

Discard velocity component of a pv-vector.

Given:
   pv      double[2][3]     pv-vector

Returned:
   p       double[3]        p-vector

Called:
   eraCp        copy p-vector



-------
Pv2s
-------

Convert position/velocity from Cartesian to spherical coordinates.

Given:
   pv       double[2][3]  pv-vector

Returned:
   theta    double        longitude angle (radians)
   phi      double        latitude angle (radians)
   r        double        radial distance
   td       double        rate of change of theta
   pd       double        rate of change of phi
   rd       double        rate of change of r

Notes:

1) If the position part of pv is null, theta, phi, td and pd
   are indeterminate.  This is handled by extrapolating the
   position through unit time by using the velocity part of
   pv.  This moves the origin without changing the direction
   of the velocity component.  If the position and velocity
   components of pv are both null, zeroes are returned for all
   six results.

2) If the position is a pole, theta, td and pd are indeterminate.
   In such cases zeroes are returned for all three.



--------
Pvdpv
--------

Inner (=scalar=dot) product of two pv-vectors.

Given:
   a        double[2][3]      first pv-vector
   b        double[2][3]      second pv-vector

Returned:
   adb      double[2]         a . b (see note)

Note:

   If the position and velocity components of the two pv-vectors are
   ( ap, av ) and ( bp, bv ), the result, a . b, is the pair of
   numbers ( ap . bp , ap . bv + av . bp ).  The two numbers are the
   dot-product of the two p-vectors and its derivative.

Called:
   eraPdp       scalar product of two p-vectors



------
Pvm
------

Modulus of pv-vector.

Given:
   pv     double[2][3]   pv-vector

Returned:
   r      double         modulus of position component
   s      double         modulus of velocity component

Called:
   eraPm        modulus of p-vector



--------
Pvmpv
--------

Subtract one pv-vector from another.

Given:
   a       double[2][3]      first pv-vector
   b       double[2][3]      second pv-vector

Returned:
   amb     double[2][3]      a - b

Note:
   It is permissible to re-use the same array for any of the
   arguments.

Called:
   eraPmp       p-vector minus p-vector



--------
Pvppv
--------

Add one pv-vector to another.

Given:
   a        double[2][3]      first pv-vector
   b        double[2][3]      second pv-vector

Returned:
   apb      double[2][3]      a + b

Note:
   It is permissible to re-use the same array for any of the
   arguments.

Called:
   eraPpp       p-vector plus p-vector



---------
Pvstar
---------

Convert star position+velocity vector to catalog coordinates.

Given (Note 1):
   pv     double[2][3]   pv-vector (AU, AU/day)

Returned (Note 2):
   ra     double         right ascension (radians)
   dec    double         declination (radians)
   pmr    double         RA proper motion (radians/year)
   pmd    double         Dec proper motion (radians/year)
   px     double         parallax (arcsec)
   rv     double         radial velocity (km/s, positive = receding)

Returned (function value):
          int            status:
                            0 = OK
                           -1 = superluminal speed (Note 5)
                           -2 = null position vector

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

Called:
   eraPn        decompose p-vector into modulus and direction
   eraPdp       scalar product of two p-vectors
   eraSxp       multiply p-vector by scalar
   eraPmp       p-vector minus p-vector
   eraPm        modulus of p-vector
   eraPpp       p-vector plus p-vector
   eraPv2s      pv-vector to spherical
   eraAnp       normalize angle into range 0 to 2pi

Reference:

   Stumpff, P., 1985, Astron.Astrophys. 144, 232-240.



------
Pvu
------

Update a pv-vector.

Given:
   dt       double           time interval
   pv       double[2][3]     pv-vector

Returned:
   upv      double[2][3]     p updated, v unchanged

Notes:

1) "Update" means "refer the position component of the vector
   to a new date dt time units from the existing date".

2) The time units of dt must match those of the velocity.

3) It is permissible for pv and upv to be the same array.

Called:
   eraPpsp      p-vector plus scaled p-vector
   eraCp        copy p-vector



-------
Pvup
-------

Update a pv-vector, discarding the velocity component.

Given:
   dt       double            time interval
   pv       double[2][3]      pv-vector

Returned:
   p        double[3]         p-vector

Notes:

1) "Update" means "refer the position component of the vector to a
   new date dt time units from the existing date".

2) The time units of dt must match those of the velocity.



--------
Pvxpv
--------

Outer (=vector=cross) product of two pv-vectors.

Given:
   a        double[2][3]      first pv-vector
   b        double[2][3]      second pv-vector

Returned:
   axb      double[2][3]      a x b

Notes:

1) If the position and velocity components of the two pv-vectors are
   ( ap, av ) and ( bp, bv ), the result, a x b, is the pair of
   vectors ( ap x bp, ap x bv + av x bp ).  The two vectors are the
   cross-product of the two p-vectors and its derivative.

2) It is permissible to re-use the same array for any of the
   arguments.

Called:
   eraCpv       copy pv-vector
   eraPxp       vector product of two p-vectors
   eraPpp       p-vector plus p-vector



------
Pxp
------

p-vector outer (=vector=cross) product.

Given:
   a        double[3]      first p-vector
   b        double[3]      second p-vector

Returned:
   axb      double[3]      a x b

Note:
   It is permissible to re-use the same array for any of the
   arguments.



-------
Rm2v
-------

Express an r-matrix as an r-vector.

Given:
   r        double[3][3]    rotation matrix

Returned:
   w        double[3]       rotation vector (Note 1)

Notes:

1) A rotation matrix describes a rotation through some angle about
   some arbitrary axis called the Euler axis.  The "rotation vector"
   returned by this function has the same direction as the Euler axis,
   and its magnitude is the angle in radians.  (The magnitude and
   direction can be separated by means of the function eraPn.)

2) If r is null, so is the result.  If r is not a rotation matrix
   the result is undefined;  r must be proper (i.e. have a positive
   determinant) and real orthogonal (inverse = transpose).

3) The reference frame rotates clockwise as seen looking along
   the rotation vector from the origin.



-------
Rv2m
-------

Form the r-matrix corresponding to a given r-vector.

Given:
   w        double[3]      rotation vector (Note 1)

Returned:
   r        double[3][3]    rotation matrix

Notes:

1) A rotation matrix describes a rotation through some angle about
   some arbitrary axis called the Euler axis.  The "rotation vector"
   supplied to This function has the same direction as the Euler
   axis, and its magnitude is the angle in radians.

2) If w is null, the unit matrix is returned.

3) The reference frame rotates clockwise as seen looking along the
   rotation vector from the origin.



-----
Rx
-----

Rotate an r-matrix about the x-axis.

Given:
   phi    double          angle (radians)

Given and returned:
   r      double[3][3]    r-matrix, rotated

Notes:

1) Calling this function with positive phi incorporates in the
   supplied r-matrix r an additional rotation, about the x-axis,
   anticlockwise as seen looking towards the origin from positive x.

2) The additional rotation can be represented by this matrix:

       (  1        0            0      )
       (                               )
       (  0   + cos(phi)   + sin(phi)  )
       (                               )
       (  0   - sin(phi)   + cos(phi)  )



------
Rxp
------

Multiply a p-vector by an r-matrix.

Given:
   r        double[3][3]    r-matrix
   p        double[3]       p-vector

Returned:
   rp       double[3]       r * p

Note:
   It is permissible for p and rp to be the same array.

Called:
   eraCp        copy p-vector



-------
Rxpv
-------

Multiply a pv-vector by an r-matrix.

Given:
   r        double[3][3]    r-matrix
   pv       double[2][3]    pv-vector

Returned:
   rpv      double[2][3]    r * pv

Note:
   It is permissible for pv and rpv to be the same array.

Called:
   eraRxp       product of r-matrix and p-vector



------
Rxr
------

Multiply two r-matrices.

Given:
   a        double[3][3]    first r-matrix
   b        double[3][3]    second r-matrix

Returned:
   atb      double[3][3]    a * b

Note:
   It is permissible to re-use the same array for any of the
   arguments.

Called:
   eraCr        copy r-matrix



-----
Ry
-----

Rotate an r-matrix about the y-axis.

Given:
   theta  double          angle (radians)

Given and returned:
   r      double[3][3]    r-matrix, rotated

Notes:

1) Calling this function with positive theta incorporates in the
   supplied r-matrix r an additional rotation, about the y-axis,
   anticlockwise as seen looking towards the origin from positive y.

2) The additional rotation can be represented by this matrix:

       (  + cos(theta)     0      - sin(theta)  )
       (                                        )
       (       0           1           0        )
       (                                        )
       (  + sin(theta)     0      + cos(theta)  )



-----
Rz
-----

Rotate an r-matrix about the z-axis.

Given:
   psi    double          angle (radians)

Given and returned:
   r      double[3][3]    r-matrix, rotated

Notes:

1) Calling this function with positive psi incorporates in the
   supplied r-matrix r an additional rotation, about the z-axis,
   anticlockwise as seen looking towards the origin from positive z.

2) The additional rotation can be represented by this matrix:

       (  + cos(psi)   + sin(psi)     0  )
       (                                 )
       (  - sin(psi)   + cos(psi)     0  )
       (                                 )
       (       0            0         1  )



------
S00
------

The CIO locator s, positioning the Celestial Intermediate Origin on
the equator of the Celestial Intermediate Pole, given the CIP's X,Y
coordinates.  Compatible with IAU 2000A precession-nutation.

Given:
   date1,date2   double    TT as a 2-part Julian Date (Note 1)
   x,y           double    CIP coordinates (Note 3)

Returned (function value):
                 double    the CIO locator s in radians (Note 2)

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

Called:
   eraFal03     mean anomaly of the Moon
   eraFalp03    mean anomaly of the Sun
   eraFaf03     mean argument of the latitude of the Moon
   eraFad03     mean elongation of the Moon from the Sun
   eraFaom03    mean longitude of the Moon's ascending node
   eraFave03    mean longitude of Venus
   eraFae03     mean longitude of Earth
   eraFapa03    general accumulated precession in longitude

References:

   Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
   "Expressions for the Celestial Intermediate Pole and Celestial
   Ephemeris Origin consistent with the IAU 2000A precession-
   nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

   n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



-------
S00a
-------

The CIO locator s, positioning the Celestial Intermediate Origin on
the equator of the Celestial Intermediate Pole, using the IAU 2000A
precession-nutation model.

Given:
   date1,date2  double    TT as a 2-part Julian Date (Note 1)

Returned (function value):
                double    the CIO locator s in radians (Note 2)

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

Called:
   eraPnm00a     classical NPB matrix, IAU 2000A
   eraBnp2xy     extract CIP X,Y from the BPN matrix
   eraS00        the CIO locator s, given X,Y, IAU 2000A

References:

   Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
   "Expressions for the Celestial Intermediate Pole and Celestial
   Ephemeris Origin consistent with the IAU 2000A precession-
   nutation model", Astron.Astrophys. 400, 1145-1154 (2003)

   n.b. The celestial ephemeris origin (CEO) was renamed "celestial
        intermediate origin" (CIO) by IAU 2006 Resolution 2.

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



-------
S00b
-------

The CIO locator s, positioning the Celestial Intermediate Origin on
the equator of the Celestial Intermediate Pole, using the IAU 2000B
precession-nutation model.

Given:
   date1,date2  double    TT as a 2-part Julian Date (Note 1)

Returned (function value):
                double    the CIO locator s in radians (Note 2)

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

Called:
   eraPnm00b     classical NPB matrix, IAU 2000B
   eraBnp2xy     extract CIP X,Y from the BPN matrix
   eraS00        the CIO locator s, given X,Y, IAU 2000A

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
S06
------

The CIO locator s, positioning the Celestial Intermediate Origin on
the equator of the Celestial Intermediate Pole, given the CIP's X,Y
coordinates.  Compatible with IAU 2006/2000A precession-nutation.

Given:
   date1,date2   double    TT as a 2-part Julian Date (Note 1)
   x,y           double    CIP coordinates (Note 3)

Returned (function value):
                 double    the CIO locator s in radians (Note 2)

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

Called:
   eraFal03     mean anomaly of the Moon
   eraFalp03    mean anomaly of the Sun
   eraFaf03     mean argument of the latitude of the Moon
   eraFad03     mean elongation of the Moon from the Sun
   eraFaom03    mean longitude of the Moon's ascending node
   eraFave03    mean longitude of Venus
   eraFae03     mean longitude of Earth
   eraFapa03    general accumulated precession in longitude

References:

   Capitaine, N., Wallace, P.T. & Chapront, J., 2003, Astron.
   Astrophys. 432, 355

   McCarthy, D.D., Petit, G. (eds.) 2004, IERS Conventions (2003),
   IERS Technical Note No. 32, BKG



-------
S06a
-------

The CIO locator s, positioning the Celestial Intermediate Origin on
the equator of the Celestial Intermediate Pole, using the IAU 2006
precession and IAU 2000A nutation models.

Given:
   date1,date2  double    TT as a 2-part Julian Date (Note 1)

Returned (function value):
                double    the CIO locator s in radians (Note 2)

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

Called:
   eraPnm06a    classical NPB matrix, IAU 2006/2000A
   eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
   eraS06       the CIO locator s, given X,Y, IAU 2006

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
S2c
------

Convert spherical coordinates to Cartesian.

Given:
   theta    double       longitude angle (radians)
   phi      double       latitude angle (radians)

Returned:
   c        double[3]    direction cosines



------
S2p
------

Convert spherical polar coordinates to p-vector.

Given:
   theta   double       longitude angle (radians)
   phi     double       latitude angle (radians)
   r       double       radial distance

Returned:
   p       double[3]    Cartesian coordinates

Called:
   eraS2c       spherical coordinates to unit vector
   eraSxp       multiply p-vector by scalar



-------
S2pv
-------

Convert position/velocity from spherical to Cartesian coordinates.

Given:
   theta    double          longitude angle (radians)
   phi      double          latitude angle (radians)
   r        double          radial distance
   td       double          rate of change of theta
   pd       double          rate of change of phi
   rd       double          rate of change of r

Returned:
   pv       double[2][3]    pv-vector



--------
S2xpv
--------

Multiply a pv-vector by two scalars.

Given:
   s1     double         scalar to multiply position component by
   s2     double         scalar to multiply velocity component by
   pv     double[2][3]   pv-vector

Returned:
   spv    double[2][3]   pv-vector: p scaled by s1, v scaled by s2

Note:
   It is permissible for pv and spv to be the same array.

Called:
   eraSxp       multiply p-vector by scalar



-------
Sepp
-------

Angular separation between two p-vectors.

Given:
   a      double[3]    first p-vector (not necessarily unit length)
   b      double[3]    second p-vector (not necessarily unit length)

Returned (function value):
          double       angular separation (radians, always positive)

Notes:

1) If either vector is null, a zero result is returned.

2) The angular separation is most simply formulated in terms of
   scalar product.  However, this gives poor accuracy for angles
   near zero and pi.  The present algorithm uses both cross product
   and dot product, to deliver full accuracy whatever the size of
   the angle.

Called:
   eraPxp       vector product of two p-vectors
   eraPm        modulus of p-vector
   eraPdp       scalar product of two p-vectors



-------
Seps
-------

Angular separation between two sets of spherical coordinates.

Given:
   al     double       first longitude (radians)
   ap     double       first latitude (radians)
   bl     double       second longitude (radians)
   bp     double       second latitude (radians)

Returned (function value):
          double       angular separation (radians)

Called:
   eraS2c       spherical coordinates to unit vector
   eraSepp      angular separation between two p-vectors



-------
Sp00
-------

The TIO locator s', positioning the Terrestrial Intermediate Origin
on the equator of the Celestial Intermediate Pole.

Given:
   date1,date2  double    TT as a 2-part Julian Date (Note 1)

Returned (function value):
                double    the TIO locator s' in radians (Note 2)

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

2) The TIO locator s' is obtained from polar motion observations by
   numerical integration, and so is in essence unpredictable.
   However, it is dominated by a secular drift of about
   47 microarcseconds per century, which is the approximation
   evaluated by the present function.

Reference:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



---------
Starpm
---------

Star proper motion:  update star catalog data for space motion.

Given:
   ra1    double     right ascension (radians), before
   dec1   double     declination (radians), before
   pmr1   double     RA proper motion (radians/year), before
   pmd1   double     Dec proper motion (radians/year), before
   px1    double     parallax (arcseconds), before
   rv1    double     radial velocity (km/s, +ve = receding), before
   ep1a   double     "before" epoch, part A (Note 1)
   ep1b   double     "before" epoch, part B (Note 1)
   ep2a   double     "after" epoch, part A (Note 1)
   ep2b   double     "after" epoch, part B (Note 1)

Returned:
   ra2    double     right ascension (radians), after
   dec2   double     declination (radians), after
   pmr2   double     RA proper motion (radians/year), after
   pmd2   double     Dec proper motion (radians/year), after
   px2    double     parallax (arcseconds), after
   rv2    double     radial velocity (km/s, +ve = receding), after

Returned (function value):
          int        status:
                        -1 = system error (should not occur)
                         0 = no warnings or errors
                         1 = distance overridden (Note 6)
                         2 = excessive velocity (Note 7)
                         4 = solution didn't converge (Note 8)
                      else = binary logical OR of the above warnings

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
   constant VMAX in the function eraStarpv),  it is arbitrarily set
   to zero.  When this action occurs, 2 is added to the status.

8) The relativistic adjustment carried out in the eraStarpv function
   involves an iterative calculation.  If the process fails to
   converge within a set number of iterations, 4 is added to the
   status.

Called:
   eraStarpv    star catalog data to space motion pv-vector
   eraPvu       update a pv-vector
   eraPdp       scalar product of two p-vectors
   eraPvstar    space motion pv-vector to star catalog data



---------
Starpv
---------

Convert star catalog coordinates to position+velocity vector.

Given (Note 1):
   ra     double        right ascension (radians)
   dec    double        declination (radians)
   pmr    double        RA proper motion (radians/year)
   pmd    double        Dec proper motion (radians/year)
   px     double        parallax (arcseconds)
   rv     double        radial velocity (km/s, positive = receding)

Returned (Note 2):
   pv     double[2][3]  pv-vector (AU, AU/day)

Returned (function value):
          int           status:
                            0 = no warnings
                            1 = distance overridden (Note 6)
                            2 = excessive speed (Note 7)
                            4 = solution didn't converge (Note 8)
                         else = binary logical OR of the above

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

Called:
   eraS2pv      spherical coordinates to pv-vector
   eraPm        modulus of p-vector
   eraZp        zero p-vector
   eraPn        decompose p-vector into modulus and direction
   eraPdp       scalar product of two p-vectors
   eraSxp       multiply p-vector by scalar
   eraPmp       p-vector minus p-vector
   eraPpp       p-vector plus p-vector

Reference:

   Stumpff, P., 1985, Astron.Astrophys. 144, 232-240.



------
Sxp
------

Multiply a p-vector by a scalar.

Given:
   s      double        scalar
   p      double[3]     p-vector

Returned:
   sp     double[3]     s * p

Note:
   It is permissible for p and sp to be the same array.



-------
Sxpv
-------

Multiply a pv-vector by a scalar.

Given:
   s       double          scalar
   pv      double[2][3]    pv-vector

Returned:
   spv     double[2][3]    s * pv

Note:
   It is permissible for pv and psv to be the same array

Called:
   eraS2xpv     multiply pv-vector by two scalars



--------
Taitt
--------

Time scale transformation:  International Atomic Time, TAI, to
Terrestrial Time, TT.

Given:
   tai1,tai2  double    TAI as a 2-part Julian Date

Returned:
   tt1,tt2    double    TT as a 2-part Julian Date

Returned (function value):
              int       status:  0 = OK

Note:

   tai1+tai2 is Julian Date, apportioned in any convenient way
   between the two arguments, for example where tai1 is the Julian
   Day Number and tai2 is the fraction of a day.  The returned
   tt1,tt2 follow suit.

References:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992)



---------
Taiut1
---------

Time scale transformation:  International Atomic Time, TAI, to
Universal Time, UT1.

Given:
   tai1,tai2  double    TAI as a 2-part Julian Date
   dta        double    UT1-TAI in seconds

Returned:
   ut11,ut12  double    UT1 as a 2-part Julian Date

Returned (function value):
              int       status:  0 = OK

Notes:

1) tai1+tai2 is Julian Date, apportioned in any convenient way
   between the two arguments, for example where tai1 is the Julian
   Day Number and tai2 is the fraction of a day.  The returned
   UT11,UT12 follow suit.

2) The argument dta, i.e. UT1-TAI, is an observed quantity, and is
   available from IERS tabulations.

Reference:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992)



---------
Taiutc
---------

Time scale transformation:  International Atomic Time, TAI, to
Coordinated Universal Time, UTC.

Given:
   tai1,tai2  double   TAI as a 2-part Julian Date (Note 1)

Returned:
   utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-3)

Returned (function value):
              int      status: +1 = dubious year (Note 4)
                                0 = OK
                               -1 = unacceptable date

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
   length is 86399, 86400 or 86401 SI seconds.

3) The function eraD2dtf can be used to transform the UTC quasi-JD
   into calendar date and clock time, including UTC leap second
   handling.

4) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale and that are too far in the future
   to be trusted.  See eraDat for further details.

Called:
   eraJd2cal    JD to Gregorian calendar
   eraDat       delta(AT) = TAI-UTC
   eraCal2jd    Gregorian calendar to JD

References:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992)



---------
Tcbtdb
---------

Time scale transformation:  Barycentric Coordinate Time, TCB, to
Barycentric Dynamical Time, TDB.

Given:
   tcb1,tcb2  double    TCB as a 2-part Julian Date

Returned:
   tdb1,tdb2  double    TDB as a 2-part Julian Date

Returned (function value):
              int       status:  0 = OK

Notes:

1) tcb1+tcb2 is Julian Date, apportioned in any convenient way
   between the two arguments, for example where tcb1 is the Julian
   Day Number and tcb2 is the fraction of a day.  The returned
   tdb1,tdb2 follow suit.

2) The 2006 IAU General Assembly introduced a conventional linear
   transformation between TDB and TCB.  This transformation
   compensates for the drift between TCB and terrestrial time TT,
   and keeps TDB approximately centered on TT.  Because the
   relationship between TT and TCB depends on the adopted solar
   system ephemeris, the degree of alignment between TDB and TT over
   long intervals will vary according to which ephemeris is used.
   Former definitions of TDB attempted to avoid this problem by
   stipulating that TDB and TT should differ only by periodic
   effects.  This is a good description of the nature of the
   relationship but eluded precise mathematical formulation.  The
   conventional linear relationship adopted in 2006 sidestepped
   these difficulties whilst delivering a TDB that in practice was
   consistent with values before that date.

3) TDB is essentially the same as Teph, the time argument for the
   JPL solar system ephemerides.

Reference:

   IAU 2006 Resolution B3



--------
Tcgtt
--------

Time scale transformation:  Geocentric Coordinate Time, TCG, to
Terrestrial Time, TT.

Given:
   tcg1,tcg2  double    TCG as a 2-part Julian Date

Returned:
   tt1,tt2    double    TT as a 2-part Julian Date

Returned (function value):
              int       status:  0 = OK

Note:

   tcg1+tcg2 is Julian Date, apportioned in any convenient way
   between the two arguments, for example where tcg1 is the Julian
   Day Number and tcg22 is the fraction of a day.  The returned
   tt1,tt2 follow suit.

References:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),.
   IERS Technical Note No. 32, BKG (2004)

   IAU 2000 Resolution B1.9



---------
Tdbtcb
---------

Time scale transformation:  Barycentric Dynamical Time, TDB, to
Barycentric Coordinate Time, TCB.

Given:
   tdb1,tdb2  double    TDB as a 2-part Julian Date

Returned:
   tcb1,tcb2  double    TCB as a 2-part Julian Date

Returned (function value):
              int       status:  0 = OK

Notes:

1) tdb1+tdb2 is Julian Date, apportioned in any convenient way
   between the two arguments, for example where tdb1 is the Julian
   Day Number and tdb2 is the fraction of a day.  The returned
   tcb1,tcb2 follow suit.

2) The 2006 IAU General Assembly introduced a conventional linear
   transformation between TDB and TCB.  This transformation
   compensates for the drift between TCB and terrestrial time TT,
   and keeps TDB approximately centered on TT.  Because the
   relationship between TT and TCB depends on the adopted solar
   system ephemeris, the degree of alignment between TDB and TT over
   long intervals will vary according to which ephemeris is used.
   Former definitions of TDB attempted to avoid this problem by
   stipulating that TDB and TT should differ only by periodic
   effects.  This is a good description of the nature of the
   relationship but eluded precise mathematical formulation.  The
   conventional linear relationship adopted in 2006 sidestepped
   these difficulties whilst delivering a TDB that in practice was
   consistent with values before that date.

3) TDB is essentially the same as Teph, the time argument for the
   JPL solar system ephemerides.

Reference:

   IAU 2006 Resolution B3



--------
Tdbtt
--------

Time scale transformation:  Barycentric Dynamical Time, TDB, to
Terrestrial Time, TT.

Given:
   tdb1,tdb2  double    TDB as a 2-part Julian Date
   dtr        double    TDB-TT in seconds

Returned:
   tt1,tt2    double    TT as a 2-part Julian Date

Returned (function value):
              int       status:  0 = OK

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



-------
Tf2a
-------

Convert hours, minutes, seconds to radians.

Given:
   s         char    sign:  '-' = negative, otherwise positive
   ihour     int     hours
   imin      int     minutes
   sec       double  seconds

Returned:
   rad       double  angle in radians

Returned (function value):
             int     status:  0 = OK
                              1 = ihour outside range 0-23
                              2 = imin outside range 0-59
                              3 = sec outside range 0-59.999...

Notes:

1)  The result is computed even if any of the range checks fail.

2)  Negative ihour, imin and/or sec produce a warning status, but
    the absolute value is used in the conversion.

3)  If there are multiple errors, the status value reflects only the
    first, the smallest taking precedence.



-------
Tf2d
-------

Convert hours, minutes, seconds to days.

Given:
   s         char    sign:  '-' = negative, otherwise positive
   ihour     int     hours
   imin      int     minutes
   sec       double  seconds

Returned:
   days      double  interval in days

Returned (function value):
             int     status:  0 = OK
                              1 = ihour outside range 0-23
                              2 = imin outside range 0-59
                              3 = sec outside range 0-59.999...

Notes:

1)  The result is computed even if any of the range checks fail.

2)  Negative ihour, imin and/or sec produce a warning status, but
    the absolute value is used in the conversion.

3)  If there are multiple errors, the status value reflects only the
    first, the smallest taking precedence.



-----
Tr
-----

Transpose an r-matrix.

Given:
   r        double[3][3]    r-matrix

Returned:
   rt       double[3][3]    transpose

Note:
   It is permissible for r and rt to be the same array.

Called:
   eraCr        copy r-matrix



-------
Trxp
-------

Multiply a p-vector by the transpose of an r-matrix.

Given:
   r        double[3][3]   r-matrix
   p        double[3]      p-vector

Returned:
   trp      double[3]      r * p

Note:
   It is permissible for p and trp to be the same array.

Called:
   eraTr        transpose r-matrix
   eraRxp       product of r-matrix and p-vector



--------
Trxpv
--------

Multiply a pv-vector by the transpose of an r-matrix.

Given:
   r        double[3][3]    r-matrix
   pv       double[2][3]    pv-vector

Returned:
   trpv     double[2][3]    r * pv

Note:
   It is permissible for pv and trpv to be the same array.

Called:
   eraTr        transpose r-matrix
   eraRxpv      product of r-matrix and pv-vector



--------
Tttai
--------

Time scale transformation:  Terrestrial Time, TT, to International
Atomic Time, TAI.

Given:
   tt1,tt2    double    TT as a 2-part Julian Date

Returned:
   tai1,tai2  double    TAI as a 2-part Julian Date

Returned (function value):
              int       status:  0 = OK

Note:

   tt1+tt2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where tt1 is the Julian Day Number
   and tt2 is the fraction of a day.  The returned tai1,tai2 follow
   suit.

References:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992)



--------
Tttcg
--------

Time scale transformation:  Terrestrial Time, TT, to Geocentric
Coordinate Time, TCG.

Given:
   tt1,tt2    double    TT as a 2-part Julian Date

Returned:
   tcg1,tcg2  double    TCG as a 2-part Julian Date

Returned (function value):
              int       status:  0 = OK

Note:

   tt1+tt2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where tt1 is the Julian Day Number
   and tt2 is the fraction of a day.  The returned tcg1,tcg2 follow
   suit.

References:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)

   IAU 2000 Resolution B1.9



--------
Tttdb
--------

Time scale transformation:  Terrestrial Time, TT, to Barycentric
Dynamical Time, TDB.

Given:
   tt1,tt2    double    TT as a 2-part Julian Date
   dtr        double    TDB-TT in seconds

Returned:
   tdb1,tdb2  double    TDB as a 2-part Julian Date

Returned (function value):
              int       status:  0 = OK

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



--------
Ttut1
--------

Time scale transformation:  Terrestrial Time, TT, to Universal Time,
UT1.

Given:
   tt1,tt2    double    TT as a 2-part Julian Date
   dt         double    TT-UT1 in seconds

Returned:
   ut11,ut12  double    UT1 as a 2-part Julian Date

Returned (function value):
              int       status:  0 = OK

Notes:

1) tt1+tt2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where tt1 is the Julian Day Number
   and tt2 is the fraction of a day.  The returned ut11,ut12 follow
   suit.

2) The argument dt is classical Delta T.

Reference:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992)



---------
Ut1tai
---------

Time scale transformation:  Universal Time, UT1, to International
Atomic Time, TAI.

Given:
   ut11,ut12  double    UT1 as a 2-part Julian Date
   dta        double    UT1-TAI in seconds

Returned:
   tai1,tai2  double    TAI as a 2-part Julian Date

Returned (function value):
              int       status:  0 = OK

Notes:

1) ut11+ut12 is Julian Date, apportioned in any convenient way
   between the two arguments, for example where ut11 is the Julian
   Day Number and ut12 is the fraction of a day.  The returned
   tai1,tai2 follow suit.

2) The argument dta, i.e. UT1-TAI, is an observed quantity, and is
   available from IERS tabulations.

Reference:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992)



--------
Ut1tt
--------

Time scale transformation:  Universal Time, UT1, to Terrestrial
Time, TT.

Given:
   ut11,ut12  double    UT1 as a 2-part Julian Date
   dt         double    TT-UT1 in seconds

Returned:
   tt1,tt2    double    TT as a 2-part Julian Date

Returned (function value):
              int       status:  0 = OK

Notes:

1) ut11+ut12 is Julian Date, apportioned in any convenient way
   between the two arguments, for example where ut11 is the Julian
   Day Number and ut12 is the fraction of a day.  The returned
   tt1,tt2 follow suit.

2) The argument dt is classical Delta T.

Reference:

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992)



---------
Ut1utc
---------

Time scale transformation:  Universal Time, UT1, to Coordinated
Universal Time, UTC.

Given:
   ut11,ut12  double   UT1 as a 2-part Julian Date (Note 1)
   dut1       double   Delta UT1: UT1-UTC in seconds (Note 2)

Returned:
   utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 3,4)

Returned (function value):
              int      status: +1 = dubious year (Note 5)
                                0 = OK
                               -1 = unacceptable date

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
   introduction of the time scale and that are too far in the future
   to be trusted.  See eraDat for further details.

Called:
   eraJd2cal    JD to Gregorian calendar
   eraDat       delta(AT) = TAI-UTC
   eraCal2jd    Gregorian calendar to JD

References:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992)



---------
Utctai
---------

Time scale transformation:  Coordinated Universal Time, UTC, to
International Atomic Time, TAI.

Given:
   utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-4)

Returned:
   tai1,tai2  double   TAI as a 2-part Julian Date (Note 5)

Returned (function value):
              int      status: +1 = dubious year (Note 3)
                                0 = OK
                               -1 = unacceptable date

Notes:

1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
   convenient way between the two arguments, for example where utc1
   is the Julian Day Number and utc2 is the fraction of a day.

2) JD cannot unambiguously represent UTC during a leap second unless
   special measures are taken.  The convention in the present
   function is that the JD day represents UTC days whether the
   length is 86399, 86400 or 86401 SI seconds.

3) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale and that are too far in the future
   to be trusted.  See eraDat  for further details.

4) The function eraDtf2d converts from calendar date and time of day
   into 2-part Julian Date, and in the case of UTC implements the
   leap-second-ambiguity convention described above.

5) The returned TAI1,TAI2 are such that their sum is the TAI Julian
   Date.

Called:
   eraJd2cal    JD to Gregorian calendar
   eraDat       delta(AT) = TAI-UTC
   eraCal2jd    Gregorian calendar to JD

References:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992)



---------
Utcut1
---------

Time scale transformation:  Coordinated Universal Time, UTC, to
Universal Time, UT1.

Given:
   utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-4)
   dut1       double   Delta UT1 = UT1-UTC in seconds (Note 5)

Returned:
   ut11,ut12  double   UT1 as a 2-part Julian Date (Note 6)

Returned (function value):
              int      status: +1 = dubious year (Note 7)
                                0 = OK
                               -1 = unacceptable date

Notes:

1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
   convenient way between the two arguments, for example where utc1
   is the Julian Day Number and utc2 is the fraction of a day.

2) JD cannot unambiguously represent UTC during a leap second unless
   special measures are taken.  The convention in the present
   function is that the JD day represents UTC days whether the
   length is 86399, 86400 or 86401 SI seconds.

3) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale and that are too far in the future
   to be trusted.  See eraDat  for further details.

4) The function eraDtf2d  converts from calendar date and time of
   day into 2-part Julian Date, and in the case of UTC implements
   the leap-second-ambiguity convention described above.

5) Delta UT1 can be obtained from tabulations provided by the
   International Earth Rotation and Reference Systems Service.  It
   It is the caller's responsibility to supply a DUT argument
   containing the UT1-UTC value that matches the given UTC.

6) The returned ut11,ut12 are such that their sum is the UT1 Julian
   Date.

7) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale and that are too far in the future
   to be trusted.  See eraDat for further details.

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



-------
Xy06
-------

X,Y coordinates of celestial intermediate pole from series based
on IAU 2006 precession and IAU 2000A nutation.

Given:
   date1,date2  double     TT as a 2-part Julian Date (Note 1)

Returned:
   x,y          double     CIP X,Y coordinates (Note 2)

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

Called:
   eraFal03     mean anomaly of the Moon
   eraFalp03    mean anomaly of the Sun
   eraFaf03     mean argument of the latitude of the Moon
   eraFad03     mean elongation of the Moon from the Sun
   eraFaom03    mean longitude of the Moon's ascending node
   eraFame03    mean longitude of Mercury
   eraFave03    mean longitude of Venus
   eraFae03     mean longitude of Earth
   eraFama03    mean longitude of Mars
   eraFaju03    mean longitude of Jupiter
   eraFasa03    mean longitude of Saturn
   eraFaur03    mean longitude of Uranus
   eraFane03    mean longitude of Neptune
   eraFapa03    general accumulated precession in longitude

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



---------
Xys00a
---------

For a given TT date, compute the X,Y coordinates of the Celestial
Intermediate Pole and the CIO locator s, using the IAU 2000A
precession-nutation model.

Given:
   date1,date2  double   TT as a 2-part Julian Date (Note 1)

Returned:
   x,y          double   Celestial Intermediate Pole (Note 2)
   s            double   the CIO locator s (Note 2)

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

Called:
   eraPnm00a    classical NPB matrix, IAU 2000A
   eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
   eraS00       the CIO locator s, given X,Y, IAU 2000A

Reference:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



---------
Xys00b
---------

For a given TT date, compute the X,Y coordinates of the Celestial
Intermediate Pole and the CIO locator s, using the IAU 2000B
precession-nutation model.

Given:
   date1,date2  double   TT as a 2-part Julian Date (Note 1)

Returned:
   x,y          double   Celestial Intermediate Pole (Note 2)
   s            double   the CIO locator s (Note 2)

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

Called:
   eraPnm00b    classical NPB matrix, IAU 2000B
   eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
   eraS00       the CIO locator s, given X,Y, IAU 2000A

Reference:

   McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
   IERS Technical Note No. 32, BKG (2004)



---------
Xys06a
---------

For a given TT date, compute the X,Y coordinates of the Celestial
Intermediate Pole and the CIO locator s, using the IAU 2006
precession and IAU 2000A nutation models.

Given:
   date1,date2  double  TT as a 2-part Julian Date (Note 1)

Returned:
   x,y          double  Celestial Intermediate Pole (Note 2)
   s            double  the CIO locator s (Note 2)

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

Called:
   eraPnm06a    classical NPB matrix, IAU 2006/2000A
   eraBpn2xy    extract CIP X,Y coordinates from NPB matrix
   eraS06       the CIO locator s, given X,Y, IAU 2006

References:

   Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

   Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981



-----
Zp
-----

Zero a p-vector.

Returned:
   p        double[3]      p-vector



------
Zpv
------

Zero a pv-vector.

Returned:
   pv       double[2][3]      pv-vector

Called:
   eraZp        zero p-vector



-----
Zr
-----

Initialize an r-matrix to the null matrix.

Returned:
   r        double[3][3]    r-matrix



