'''example from sofa_ast_c.pdf'''
import math
import erfa

def reprd(title, r, d):
    print("%25s"%title, 
          " %s%2.2d %2.2d %2.2d.%7.7d"%erfa.a2tf(7, r),
          " %s%2.2d %2.2d %2.2d.%6.6d"%erfa.a2af(6, d))
    
# site longitude, latitude (radians) and height above the geoid (m).
phi = erfa.af2a(-15,57,42.8)
elong = erfa.af2a(-5,41,54.2)
hm = 625.0

# Ambient pressure (HPa), temperature (C) and rel. humidity (frac).
phpa = 952.0
tc = 18.5
rh = 0.83

# Effective color (microns).
wl = 0.55

# UTC date
utc1, utc2 = erfa.dtf2d(2013, 4, 2, 23, 15, 43.55, "UTC")

# TT date
tai1, tai2 = erfa.utctai(utc1, utc2)
tt1, tt2 = erfa.taitt(tai1, tai2)

# EOPs: polar motion in radians, UT1-UTC in seconds. 
xp = 50.995e-3 * erfa.DAS2R
yp = 376.723e-3 * erfa.DAS2R
dut1 = 155.0675e-3
##print('xp, yp', xp, yp)
# Corrections to IAU 2000A CIP (radians). 
dx = 0.269e-3 * erfa.DAS2R
dy = -0.274e-3 * erfa.DAS2R

# Star ICRS RA,Dec (radians).
rc = erfa.tf2a(14,34,16.81183)
dc = erfa.af2a(-12,31,10.3965)
##print('rc, dc', rc,dc)
#
reprd("ICRS, epoch J2000.0:", rc, dc )

# Proper motion: RA/Dec derivatives, epoch J2000.0.
pr = math.atan2(-354.45e-3 * erfa.DAS2R, math.cos(dc))
pd = 595.35e-3 * erfa.DAS2R
#print('pr, pd: ', pr, pd)

# Parallax (arcsec) and recession speed (km/s).
px = 164.99e-3
rv = 0.0

# ICRS to CIRS (geocentric observer).
ri, di, eo = erfa.atci13(rc, dc, pr, pd, px, rv, tt1, tt2)
#print('ri, di', ri, di)
#
reprd ( "catalog -> CIRS:", ri, di )

# CIRS to ICRS (astrometric).
rca, dca, eo = erfa.atic13 ( ri, di, tt1, tt2)
#
reprd ( "CIRS -> astrometric:", rca, dca )

#ICRS (astrometric) to CIRS (geocentric observer).
ri, di, eo = erfa.atci13 ( rca, dca, 0.0, 0.0, 0.0, 0.0, tt1, tt2)
reprd ( "astrometric -> CIRS:", ri, di );

# Apparent place.
ra = erfa.anp ( ri - eo )
da = di
reprd ( "geocentric apparent:", ra, da )

# CIRS to topocentric.
aot, zot, hot, dot, rot = erfa.atio13(ri, di, utc1, utc2, dut1,
                                      elong, phi, hm, xp, yp,
                                      0.0, 0.0, 0.0, 0.0)
reprd( "CIRS -> topocentric:", rot, dot )

# CIRS to topocentric.
aob, zob, hob, dob, rob = erfa.atio13(ri, di, utc1, utc2, dut1,
                                      elong, phi, hm, xp, yp,
                                      phpa, tc, rh, wl)
reprd( "CIRS -> observed:", rob, dob )

# ICRS to observed.
aob, zob, hob, dob, rob, eo = erfa.atco13 ( rc, dc, pr, pd, px, rv,
                                            utc1, utc2, dut1,
                                            elong, phi, hm, xp, yp,
                                            phpa, tc, rh, wl)
reprd ( "ICRS -> observed:", rob, dob )

# ICRS to CIRS using some user-supplied parameters.
# SOFA heliocentric Earth ephemeris.
pvh, pvb = erfa.epv00(tt1, tt2)

# JPL DE405 barycentric Earth ephemeris.
pvb = ((-0.9741704366519668, -0.2115201000882231, -0.0917583114068277),
       (0.0036436589347388, -0.0154287318503146, -0.0066892203821059))

# era 2000A CIP.
r = erfa.pnm00a(tt1, tt2)
x, y = erfa.bpn2xy(r)

# Apply IERS corrections.
x += dx
y += dy
# SOFA CIO locator. */
s = erfa.s06(tt1, tt2, x, y)

# Populate the context.
astrom = erfa.apci(tt1, tt2, pvb, pvh[0], x, y, s)

# Carry out the transformation and report the results.
ri, di = erfa.atciq(rc, dc, pr, pd, px, rv, *astrom)
reprd ( "ICRS -> CIRS (JPL, IERS):", ri, di )

# The same but with Saturn then Jupiter then Sun light deflection.
b0 = erfa.LDBODY((0.00028574, 3e-10,
                  ((-7.8101442680818964, -5.6095668114887358, -1.9807981923749924),
                   (0.0030723248971152, -0.0040699547707598, -0.0018133584165345))))
b1 = erfa.LDBODY((0.00095435, 3e-9,
                  ((0.7380987962351833, .6365869247538951, 1.9693136030111202),
                   (-0.0075581692172088, 0.0012691372216750, 0.0007279990012801))))
b2 = erfa.LDBODY((1.0, 6e-6,
                  ((-0.0007121743770509, -0.0023047830339257, -0.0010586596574639),
                   (0.0000062923521264, -0.0000003308883872, -0.0000002964866231))))
b= [b0, b1, b2]
ri, di = erfa.atciqn(rc, dc, pr, pd, px, rv, b, *astrom)
reprd ( "ICRS -> CIRS (+ planets):", ri, di )

# CIRS to ICRS (astrometric).
rca, dca = erfa.aticqn(ri, di, b, *astrom)
reprd ( "CIRS -> astrometric:", rca, dca )
