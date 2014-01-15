# -*- coding: utf-8 -*-
'''SOFA example, from sofa_ts_c.pdf
'''
from __future__ import print_function
import math
import erfa

print('''UTC to TT
transform 2010 July 24, 11:18:07.318 (UTC) into Terrestrial Time (TT)
and report it rounded to 1ms precision''')
# encode UTC date and time into internal format
u1, u2 = erfa.dtf2d(2010,7,24,11,18,7.318)

# transform UTC to TAI, then TAI to TT
a1, a2 = erfa.utctai(u1, u2)
t1, t2 = erfa.taitt(a1, a2)

# decode and report the TT
y, m, d, h, mn, sc, f = erfa.d2dtf(3, t1, t2)
print("UTC: 2010 July 24, 11:18:07.318")
print("TT: %4d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%3.3d"%(y, m, d, h, mn, sc, f))

print('=====')

print('''
TAI to UTC
take a time expressed as TAI, encode it into
the internal format and transform it into UTC''')

# encode TAI date and time into internal format
a1, a2 = erfa.dtf2d(2009,1,1,0,0,33.7, "TAI")

# decode and report TAI
print("TAI :%4d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%3.3d"%erfa.d2dtf(3, a1, a2))

# transform TAI to UTC
u1, u2 = erfa.taiutc(a1, a2)

# decode and report UTC
print("UTC :%4d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%3.3d"%erfa.d2dtf(3, u1, u2))
print('=====')
print('convert time interval')
# The time
ih = 23
im = 5
s = 11.63079
print( "%2d:%2.2d:%9.6f\n"%(ih, im, s))

print('Express as a fraction of 1 day:')
f = erfa.tf2d(ih, im, s)
print ( "%14.12f\n"%f )

print('Six hours earlier:')
f -= 0.25

#Report to 1 ms precision.
ihmsf = erfa.d2tf(3, f)
print ( "%s%2d:%2.2d:%2.2d.%3.3d\n"%ihmsf)

print('=====')
print('''
Date and time.''')
iy = 2008; im = 2; id = 29;
ihour = 23; imin = 59; sec = 59.9;
print("%4d/%2.2d/%2.2d%3d:%2.2d:%4.1f\n"%
         (iy, im, id, ihour, imin, sec))
print('Express as 2-part JD.')
d1, d2 = erfa.cal2jd ( iy, im, id)
d = erfa.tf2d (ihour, imin, sec)
d2 += d
print("%9.1f +%13.6f =%15.6f\n"%( d1, d2, d1 + d))
print('Express as calendar date and fraction of a day.')
iy, im, id, fd = erfa.jd2cal(d1, d2)
d = id + fd
print( "%4d/%2.2d/%9.6f\n"%(iy, im, d))
print('Round to 0.001 day.')
iymdf = erfa.jdcalf ( 3, d1, d2)
print( "%4d/%2.2d/%2.2d.%3.3d\n"%iymdf)
print('=====')
print('''
Besselian and Julian epochs''')
print('Julian Date.')
d = 2457073.05631
print ( "%13.5f"%d )
print('Transform into Besselian epoch:')
e = erfa.epb ( 0.0, d)
print( "B%15.10f"%e )
print('Transform back.')
d1, d2 = erfa.epb2jd ( e )
print ( "%17.9f"%(d1+d2) )
print('The same for Julian epoch.')
e = erfa.epj ( 0.0, d )
print ( "J%15.10f"%e )
d1, d2 = erfa.epj2jd ( e )
print ( "%17.9f"%(d1+d2) )

print('=====')
print('''
transform UTC to other times.
an observer at north latitude +19°28'52''.5,
west longitude 155°55'59''.6,
at sea level, on 2006 January 15 at 21:24:37.5 UTC
requires the time in all other supported time scales''')

# site terrestrial coordinates (WGS84)
latnd = +19
latnm = 28
slatn = 52.5
lonwd = -155
lonwm = 55
slonw = 59.6
hm = 0.

# transform to geocentric
phi = erfa.af2a(latnd, latnm, slatn)
elon = erfa.af2a(lonwd, lonwm, slonw)
xyz = erfa.gd2gc(1, elon, phi, hm)
u = math.hypot(xyz[0], xyz[1])
v = xyz[2]

# UTC date and time
iy = 2006
mo = 1
d = 15
ih = 21
im = 24
sec = 37.5

# transform into intenal format
utc1, utc2 = erfa.dtf2d(iy,mo,d,ih,im,sec, "UTC")

# UT1-UTC from IERS
dut = .3341

# UTC -> UT1
ut11, ut12 = erfa.utcut1(utc1, utc2, dut)

# Extract fraction for TDB-TT calculation, later.
ut = math.fmod(math.fmod(ut11,1.0)+math.fmod(ut12,1.0),1.0)

# UTC -> TAI -> TT -> TCG
tai1, tai2 = erfa.utctai(utc1, utc2)
tt1, tt2 = erfa.taitt(tai1, tai2)
tcg1, tcg2 = erfa.tttcg(tt1, tt2)

# TDB-TT (using TT as a substitute for TDB).
dtr = erfa.dtdb(tt1, tt2, ut, elon, u, v)

# TT -> TDB -> TCB.
tdb1, tdb2 = erfa.tttdb(tt1, tt2, dtr)
tcb1, tcb2 = erfa.tdbtcb(tdb1, tdb2)

# report
print("UTC %5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d"%erfa.d2dtf(6, utc1, utc2))
print("UT1 %5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d"%erfa.d2dtf(6, ut11, ut12, "ut1"))
print("TAI %5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d"%erfa.d2dtf(6, tai1, tai2, "tai"))
print("TT  %5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d"%erfa.d2dtf(6, tt1, tt2, "tt"))
print("TCG %5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d"%erfa.d2dtf(6, tcg1, tcg2, "tcg"))
print("TDB %5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d"%erfa.d2dtf(6, tdb1, tdb2, "tdb"))
print("TCB %5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d"%erfa.d2dtf(6, tcb1, tcb2, "tcb"))
