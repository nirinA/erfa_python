# -*- coding: utf-8 -*-
'''ex from sofa_pn.pdf
'''

import math
import erfa
##import scipy

from pprint import pprint

# some init
AS2R = erfa.DAS2R ##4.848136811095359935899141E-6
R2AS = erfa.DR2AS ##206264.8062470963551564734
D2PI = 6.283185307179586476925287E0

# UTC
IY = 2007
IM = 4
ID = 5
IH = 12
MIN = 0
SEC = 0.0

# Polar motion (arcsec->radians)
XP = 0.0349282 * AS2R
YP = 0.4833163 * AS2R

# UT1-UTC (s)
DUT1 = -0.072073685

# Nutation corrections wrt IAU 1976/1980 (mas->radians)
DDP80 = -55.0655 * AS2R/1000.
DDE80 = -6.3580 * AS2R/1000.

# CIP offsets wrt IAU 2000A (mas->radians)
DX00 = 0.1725 * AS2R/1000.
DY00 = -0.2650 * AS2R/1000.

# CIP offsets wrt IAU 2006/2000A (mas->radians)
DX06 = 0.1750 * AS2R/1000.
DY06 = -0.2259 * AS2R/1000.

# TT (MJD)
DJMJD0, DATE = erfa.cal2jd(IY, IM, ID)
TIME = ( 60*(60*IH + MIN) + SEC ) / 86400.
UTC = DATE + TIME
DAT = erfa.dat(IY, IM, ID, TIME)
TAI = UTC + DAT/86400.
TT = TAI + 32.184/86400.

# UT1
TUT = TIME + DUT1/86400.
UT1 = DATE + TUT

print("UTC :%4d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%3.3d"%erfa.d2dtf(3, *erfa.dtf2d(IY, IM, ID, IH, MIN, SEC)))
print("TT  = 2400000.5 + %.17f "%TT)
print("UT1 = 2400000.5 + %.17f "%UT1)

print('''
======================================
IAU 1976/1980/1982/1994, equinox based
======================================
''')

# IAU 1976 precession matrix, J2000.0 to date.
RP = erfa.pmat76(DJMJD0, TT)

# IAU 1980 nutation
DP80, DE80 = erfa.nut80(DJMJD0, TT)

# Add adjustments: frame bias, precession-rates, geophysical
DPSI = DP80 + DDP80
DEPS = DE80 + DDE80

# Mean obliquity
EPSA = erfa.obl80(DJMJD0, TT)

# Build the rotation matrix
RN = erfa.numat(EPSA, DPSI, DEPS)

# Combine the matrices: PN = N x P
RNPB = erfa.rxr(RN, RP)
print("NPB matrix, equinox based:")
pprint(RNPB)

# Equation of the equinoxes, including nutation correction
EE = erfa.eqeq94(DJMJD0, TT) + DDP80 * math.cos(EPSA)

# Greenwich apparent sidereal time (IAU 1982/1994).
GST = erfa.anp(erfa.gmst82(DJMJD0+DATE, TUT) + EE)
print("GST = %.17f radians"%GST)
print("    = %.17f degrees"%math.degrees(GST))
print("    = %s%dd%dm%d.%ds"%erfa.a2af(6, GST))
print("    = %s%dh%dm%d.%ds"%erfa.a2tf(6, GST))

# Form celestial-terrestrial matrix (no polar motion yet).
##RC2TI = erfa.cr(RNPB)
##RC2TI = erfa.rz(GST, RC2TI)
RC2TI = erfa.rz(GST, RNPB)
print("celestial to terrestrial matrix (no polar motion)")
pprint(RC2TI)

# Polar motion matrix (TIRS->ITRS, IERS 1996).
RPOM = erfa.ir()
RPOM = erfa.rx(-YP, RPOM)
RPOM = erfa.ry(-XP, RPOM)

# Form celestial-terrestrial matrix (including polar motion).
RC2IT = erfa.rxr(RPOM, RC2TI)
print("celestial to terrestrial matrix (including polar motion)")
pprint(RC2IT)

##A = scipy.matrix(RC2IT)

print('''
============================================
IAU 2000A, CIO based, using classical angles
============================================
''')
# CIP and CIO, IAU 2000A.
X, Y, S = erfa.xys00a(DJMJD0, TT)

# Add CIP corrections.
X = X + DX00
Y = Y + DY00

print("CIP corrections")
print('X = %.17f\nY = %.17f\nS = %.17f'%(X, Y, S*R2AS))

# GCRS to CIRS matrix.
RC2I = erfa.c2ixys(X, Y, S)
print("NPB matrix, CIO based")
pprint(RC2I)

# Earth rotation angle.
ERA = erfa.era00(DJMJD0+DATE, TUT)
print("Earth rotation angle")
print('ERA = %.17f radians'%ERA)
print("    = %.17f degrees"%math.degrees(ERA))
print("    = %s%dd%dm%d.%ds"%erfa.a2af(6, ERA))
print("    = %s%dh%dm%d.%ds"%erfa.a2tf(6, ERA))

# Form celestial-terrestrial matrix (no polar motion yet).
##RC2TI = erfa.cr(RC2I)
##RC2TI = erfa.rz(ERA, RC2TI)
RC2TI = erfa.rz(ERA, RC2I)
print("celestial to terrestrial matrix (no polar motion)")
pprint(RC2TI)

# Polar motion matrix (TIRS->ITRS, IERS 2003).
RPOM = erfa.pom00(XP, YP, erfa.sp00(DJMJD0,TT))

# Form celestial-terrestrial matrix (including polar motion).
RC2IT = erfa.rxr(RPOM, RC2TI)
print("celestial to terrestrial matrix (including polar motion)")
pprint(RC2IT)

##B = scipy.matrix(RC2IT)

print('''
================================================
IAU 2000A, equinox based, using classical angles
================================================
''')

# Nutation, IAU 2000A.
DP00, DE00 = erfa.nut00a(DJMJD0, TT)

# Precession-nutation quantities, IAU 2000.
EPSA, RB, RP, RPB, RN, RNPB = erfa.pn00(DJMJD0, TT, DP00, DE00)

# Transform dX,dY corrections from GCRS to mean of date.
V1 = (DX00, DY00, 0.)
V2 = erfa.rxp(RNPB, V1)
DDP00 = V2[1] / math.sin(EPSA)
DDE00 = V2[2]

# Corrected nutation.
DPSI = DP00 + DDP00
DEPS = DE00 + DDE00

# Build the rotation matrix.
RN = erfa.numat(EPSA, DPSI, DEPS)

# Combine the matrices: PN = N x P.
RNPB = erfa.rxr(RN, RPB)
print("NPB matrix, equinox based:")
pprint(RNPB)

# Greenwich apparent sidereal time (IAU 1982/1994).
GST = erfa.anp(erfa.gmst00(DJMJD0+DATE, TUT, DJMJD0, TT)+erfa.ee00(DJMJD0, TT, EPSA, DPSI))
print("GST = %.17f radians"%GST)
print("    = %.17f degrees"%math.degrees(GST))
print("    = %s%dd%dm%d.%ds"%erfa.a2af(6, GST))
print("    = %s%dh%dm%d.%ds"%erfa.a2tf(6, GST))

# Form celestial-terrestrial matrix (no polar motion yet).
RC2TI = erfa.rz(GST, RNPB)
print("celestial to terrestrial matrix (no polar motion)")
pprint(RC2TI)

# Polar motion matrix (TIRS->ITRS, IERS 2003).
RPOM = erfa.pom00(XP, YP, erfa.sp00(DJMJD0,TT))

# Form celestial-terrestrial matrix (including polar motion).
RC2IT = erfa.rxr(RPOM, RC2TI)
print("celestial to terrestrial matrix (including polar motion)")
pprint(RC2IT)

print('''
=================================================
IAU 2006/2000A, CIO based, using classical angles
=================================================
''')

# CIP and CIO , IAU 2006/2000A
X, Y, S = erfa.xys06a(DJMJD0, TT)

# Add CIP corrections.
X = X + DX06
Y = Y + DY06

print("CIP corrections")
print('X = %.17f\nY = %.17f\nS = %.17f'%(X, Y, S*R2AS))

# GCRS to CIRS matrix
RC2I = erfa.c2ixys(X, Y, S)
print("NPB matrix, CIO based")
pprint(RC2I)

# Earth rotation angle
ERA = erfa.era00(DJMJD0+DATE, TUT)
print("ERA = %.17f radians"%ERA)
print("    = %.17f degrees"%math.degrees(ERA))
print("    = %s%dd%dm%d.%ds"%erfa.a2af(6, ERA))
print("    = %s%dh%dm%d.%ds"%erfa.a2tf(6, ERA))

# Form celestial-terrestrial matrix (no polar motion yet).
RC2TI = erfa.rz(ERA, RC2I)
print("celestial to terrestrial matrix (no polar motion)")
pprint(RC2TI)

# Polar motion matrix (TIRS->ITRS, IERS 2003).
RPOM = erfa.pom00(XP, YP, erfa.sp00(DJMJD0,TT))

# Form celestial-terrestrial matrix (including polar motion).
RC2IT = erfa.rxr(RPOM, RC2TI)
print("celestial to terrestrial matrix (including polar motion)")
pprint(RC2IT)

print('''
===========================================
IAU 2006/2000A, CIO based, using X,Y series
===========================================
''')

# CIP and CIO , IAU 2006/2000A
X, Y = erfa.xy06(DJMJD0, TT)
S = erfa.s06(DJMJD0, TT, X, Y)

# Add CIP corrections.
X = X + DX06
Y = Y + DY06

print("CIP corrections")
print('X = %.17f\nY = %.17f\nS = %.17f'%(X, Y, S*R2AS))

# GCRS to CIRS matrix
RC2I = erfa.c2ixys(X, Y, S)
print("NPB matrix, CIO based")
pprint(RC2I)

# Earth rotation angle
ERA = erfa.era00(DJMJD0+DATE, TUT)
print("Earth rotation angle")
print("ERA = %.17f radians"%ERA)
print("    = %.17f degrees"%math.degrees(ERA))
print("    = %s%dd%dm%d.%ds"%erfa.a2af(6, ERA))
print("    = %s%dh%dm%d.%ds"%erfa.a2tf(6, ERA))

# Form celestial-terrestrial matrix (no polar motion yet).
RC2TI = erfa.rz(ERA, RC2I)
print("celestial to terrestrial matrix (no polar motion)")
pprint(RC2TI)

# Polar motion matrix (TIRS->ITRS, IERS 2003).
RPOM = erfa.pom00(XP, YP, erfa.sp00(DJMJD0,TT))

# Form celestial-terrestrial matrix (including polar motion).
RC2IT = erfa.rxr(RPOM, RC2TI)
print("celestial to terrestrial matrix (including polar motion)")
pprint(RC2IT)
