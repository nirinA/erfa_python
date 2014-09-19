'''display Sun and major planets ephemerides.
'''
##   usage:
##   python ex_ephem.py -L 47:38:25 -l -18:57:03 -a 1258.0 2014/02/14T23:58:53 
##
##   required argument: date
##   optional argument: -L --longitude
##                      -l --latitude
##                      -a --altitude
##
##TODO:
##   * add refraction,
##   * latitude,
##   * diurnal aberration,
##   * parallax
##   * Moon

from __future__ import print_function
import sys
import math
import time
import argparse
try:
    import enum
    Planet = enum.Enum('Planet', 'Mercury Venus EMB  Mars Jupiter Saturn Uranus Neptune')
except ImportError:
    class Planet(object):
        def __init__(self, n):
            p = {1:'Mercury',
                 2:'Venus',
                 3:'EMB',
                 4:'Mars',
                 5:'Jupiter',
                 6:'Saturn',
                 7:'Uranus',
                 8:'Neptune'}
            self.name = p[n]
    
import erfa

# observer location
latitude = '-18:57:03' #South
longitude = '47:38:25' #East
altitude = 1258.0

# Time
##utc = '2013/08/26T09:52:05'
utc= '2015/02/01T23:58:50'
tt = '2015/02/02T00:00:00'

def lmst(gmst, long):
    '''local mean sideral time
Given:
    gmst
    long in degrees, positive for place east of Greenwich'''
    return gmst + (3600.*long/15.)

# gast

def print_ln(title, a, d):
        print(title, end=': ')
        print("%s%dh:%dm:%d.%ds"%erfa.a2tf(4, a), end=', ')
        print("%s%dd:%dm:%d.%ds"%erfa.a2af(4, d))

    
def strputc(string):
    '''parse a utc string and return utc1, utc2'''
    Y,m,d,H,M,S,w,y,dst=time.strptime(string, '%Y/%m/%dT%H:%M:%S')
    return erfa.dtf2d(Y,m,d,H,M,S)

class Observer(object):
    def __init__(self, longitude, latitude, altitude):
        # site terrestrial coordinates (WGS84)
        self.latnd, self.latnm, self.slatn = map(int, latitude.split(':'))
        self.lonwd, self.lonwm, self.slonw = map(int, longitude.split(':'))
        self.hm = altitude
        # transform to geocentric
        self.phi = erfa.af2a(self.latnd, self.latnm, self.slatn)
        self.elon = erfa.af2a(self.lonwd, self.lonwm, self.slonw)
        self.xyz = erfa.gd2gc(1, self.elon, self.phi, self.hm)
        self.u = math.hypot(self.xyz[0], self.xyz[1])
        self.v = self.xyz[2]

class TimeScale(Observer):
    def __init__(self, utc, longitude, latitude, altitude):
        try:
            super().__init__(longitude, latitude, altitude)
        except TypeError:
            super(TimeScale, self).__init__(longitude, latitude, altitude)
        self.utc = utc
        self.utc1, self.utc2 = strputc(self.utc)
        self.tai1, self.tai2 = erfa.utctai(self.utc1, self.utc2)
        self.tt1, self.tt2 = erfa.taitt(self.tai1, self.tai2)
        self.tcg1, self.tcg2 = erfa.tttcg(self.tt1, self.tt2)
        self.dut1 = -0.1
        if (self.utc1+self.utc2) >= 2456708.5:
            self.dut1 = -0.2
        elif (self.utc1+self.utc2) >= 2456785.5:
            self.dut1 = -0.3
        elif (self.utc1+self.utc2) >= 2456925.5:
            self.dut1 = -0.4
        self.dt = 35+32.184-self.dut1
        self.ut11, self.ut12 = erfa.utcut1(self.utc1, self.utc2, self.dut1)
        # Extract fraction for TDB-TT calculation, later.
        self.ut = math.fmod(math.fmod(self.ut11,1.0)+math.fmod(self.ut12,1.0),1.0)
        # TDB-TT (using TT as a substitute for TDB).
        self.dtr = erfa.dtdb(self.tt1, self.tt2, self.ut, self.elon, self.u, self.v)
        self.tdb1, self.tdb2 = erfa.tttdb(self.tt1, self.tt2, self.dtr)
        self.tcb1, self.tcb2 = erfa.tdbtcb(self.tdb1, self.tdb2)

class Ephem(TimeScale, Observer):
    '''approximate RA and DEC of the Sun, the Moon and major planets
use plan94, then transform the heliocentric coordinates to equatorial'''
    def __init__(self, utc, longitude, latitude, altitude):
        try:
            super().__init__(utc, longitude, latitude, altitude)
        except TypeError:
            super(Ephem, self).__init__(utc, longitude, latitude, altitude)
        # julian centuries at the given UTC
        self.jc = (self.tdb1+self.tdb2-erfa.DJ00)/erfa.DJC
        # heliocentric position of major bodies
        self.p94 = [erfa.plan94(self.tdb1, self.tdb2, i) for i in range(1,9)]
        # Earth heliocentric position and velocity
        self.eph, self.epb = erfa.epv00(self.tdb1, self.tdb2)

        self.distances_from_earth = {}
        self.RA_DEC = {}

##        # precession angles
##        eps0,psia,oma,bpa,bqa,pia,bpia,epsa,chia,za,zetaa,thetaa,pa,gam,phi,psi = erfa.p06e(self.tt1, self.tt2)
##
##        # Fundamental argument of major bodies
##        mer = erfa.fame03(self.jc))
        self.sun()
        self.planets()

    def sun(self):
        '''RA and DEC of the Sun'''
        x, y, z = self.eph[0][0], self.eph[0][1], self.eph[0][2]
        self.distances_from_earth.update(
            {'Sun':math.sqrt(x*x + y*y +z*z)}
            )
        a, d = erfa.c2s((-x, -y, -z))
        a = erfa.anp(a)
        self.RA_DEC.update(
            {'Sun':(a, d)}
            )

    def planets(self):
        '''RA and DEC of major planets and EMB'''
        for i in range(8):
            x = self.p94[i][0][0] - self.eph[0][0]
            y = self.p94[i][0][1] - self.eph[0][1]
            z = self.p94[i][0][2] - self.eph[0][2]
            self.distances_from_earth.update(
                {Planet(i+1).name:math.sqrt(x*x + y*y +z*z)}
                )
            a, d = erfa.c2s((x,y,z))
            a = erfa.anp(a)
            self.RA_DEC.update(
                {Planet(i+1).name:(a, d)}
                )

## need more coeffs from periodic term, see ELP
    def moon(self):
        ## mean elongation of the Moon from the Sun
        m_elongation = erfa.fad03(self.jc)
        ## mean longitude of Earth
        e_longitude = erfa.fae03(self.jc)
        ## mean argument of the latitude of the Moon,
        ## mean longitude of the Moon minus mean longitude of the ascending node
        m_latitude = erfa.faf03(self.jc)
        ## mean anomaly of the Moon
        m_anomaly = erfa.fal03(self.jc)
        ## mean anomaly of the Sun
        s_anomaly = erfa.falp03(self.jc)
        ## mean longitude of the Moon's ascending node
        m_longitude_asc = erfa.faom03(self.jc)
        return m_elongation, m_latitude, m_anomaly, s_anomaly, m_longitude_asc

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        epilog = 'Ephemerides')
    parser.add_argument(
        'date',
        help = 'utc date')
    parser.add_argument(
        '-L', '--longitude',
        help = 'observer longitude')
    parser.add_argument(
        '-l', '--latitude',
        help = 'observer latitude')
    parser.add_argument(
        '-a', '--altitude',
        help = 'observer altitude')
    
    res = parser.parse_args(sys.argv[1:])
    if res.longitude:
        longitude = res.longitude
    else:
        longitude = '00:00:00'
    if res.latitude:
        latitude = res.latitude
    else:
        latitude = '00:00:00'
    if res.altitude:
        altitude = float(res.altitude)
    else:
        altitude = 0.0
    e = Ephem(res.date, longitude, latitude, altitude)
    print('Body : RA, DEC')
    print_ln('Sun', *e.RA_DEC['Sun'])
    for i in range(1,8):
        print_ln(Planet(i).name, *e.RA_DEC[Planet(i).name])
