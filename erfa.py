import math
from _erfa import *
from _erfa import __doc__

def ir():
    '''Initialize an r-matrix to the identity matrix.'''
    return [[1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.]]

def zp():
    '''Zero a p-vector.'''
    return [0., 0., 0.]

def zpv():
    '''Zero a pv-vector.'''
    return [[0., 0., 0.],
            [0., 0., 0.]]

def zr():
    '''Initialize an r-matrix to the null matrix.'''
    return [[0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.]]

def aper(theta, astrom):
    '''aper(theta, astrom) -> astrom
In the star-independent astrometry parameters, update only the
Earth rotation angle, supplied by the caller explicitly.'''
    pmt, eb, eh, em, v, bm1, bpn, along, phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb = astrom 
    return ASTROM((pmt, eb, eh, em, v, bm1, bpn, along, phi, xpl, ypl, sphi, cphi, diurab, theta+along, refa, refb))

def aper13(ut11, ut12, astrom):
    '''aper13(ut11, ut12, astrom) -> astrom
In the star-independent astrometry parameters, update only the
Earth rotation angle. The caller provides UT1, (n.b. not UTC).'''
    era = era00(ut11, ut12)
    return aper(era, astrom)

def ldn(ldbody, ob, sc):
    '''ldn(ldbody[], ob[3], sc[3])) -> sn[3]
For a star, apply light deflection by multiple solar-system bodies,
as part of transforming coordinate direction into natural direction.
Given:
    b    data for each of the n bodies
    ob   barycentric position of the observer (au)
    sc   observer to star coord direction (unit vector)
Returned:
    sn   observer to deflected star (unit vector)'''
    # Light time for 1 AU (days)
    CR = AULT/DAYSEC
    # Star direction prior to deflection.
    sn = list(cp(sc))
    # Body by body.
    for l in ldbody:
        # Body to observer vector at epoch of observation (au).
        v = pmp(ob, l.pv[0])
        # Minus the time since the light passed the body (days). 
        dt = pdp(sn, v) * CR
        # Neutralize if the star is "behind" the observer. 
        dt = min(dt, 0.)
        # Backtrack the body to the time the light was passing the body.
        ev = ppsp(v, -dt, l.pv[1])
        # Body to observer vector as magnitude and direction. 
        em, e = pn(ev)
        # Apply light deflection for this body. 
        sn = ld(l.bm, sn, sn, e, em, l.dl)
    return sn

def atciqn(rc, dc, pr, pd, px,rv, astrom, ldbody):
    '''atciqn(rc, dc, pr, pd, px, rv, astrom, ldbody) -> ri,di
Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
star-independent astrometry parameters plus a list of light-deflecting bodies.

Use of this function is appropriate when efficiency is important and
where many star positions are to be transformed for one date.  The
star-independent parameters can be obtained by calling one of the
functions apci[13], apcg[13], apco[13] or apcs[13].

If the only light-deflecting body to be taken into account is the
Sun, the atciq function can be used instead.  If in addition the
parallax and proper motions are zero, the aciqz function can be used.
Given:
    rc,dc      ICRS RA,Dec at J2000.0 (radians)
    pr         RA proper motion (radians/year)
    pd         Dec proper motion (radians/year)
    px         parallax (arcsec)
    rv         radial velocity (km/s, +ve if receding)
    astrom     star-independent astrometry parameters
    ldbody     data for each of the n bodies
Returned:
    ri,di      CIRS RA,Dec (radians)'''
    # Proper motion and parallax, giving BCRS coordinate direction.
    pco = pmpx(rc, dc, pr, pd, px, rv, astrom.pmt, astrom.eb)
    # Light deflection, giving BCRS natural direction.
    pnat = ldn(ldbody, astrom.eb, pco)
    # Aberration, giving GCRS proper direction.
    ppr = ab(pnat, astrom.v, astrom.em, astrom.bm1)
    # Bias-precession-nutation, giving CIRS proper direction.
    pi = rxp(astrom.bpn, ppr)
    # CIRS RA,Dec.
    w, di = c2s(pi)
    ri = anp(w)
    return ri, di

def aticqn(ri, di, astrom, ldbody):
    '''aticqn(ri, di, astrom, ldbody) -> rc, dc
Quick CIRS to ICRS astrometric place transformation, given the star-
independent astrometry parameters plus a list of light-deflecting
bodies.

Use of this function is appropriate when efficiency is important and
where many star positions are all to be transformed for one date.
The star-independent astrometry parameters can be obtained by
calling one of the functions apci[13], apcg[13], apco[13]
or eraApcs[13].

If the only light-deflecting body to be taken into account is the
Sun, the eraAticq function can be used instead.
Given:
    ri,di  CIRS RA,Dec (radians)
    astrom star-independent astrometry parameters:
    b      data for each of the n bodies
Returned:
    rc,dc  ICRS astrometric RA,Dec (radians)'''
    before = [0,0,0]
    after = [0,0,0]
    pnat = [0,0,0]
    pco = [0,0,0]
    # CIRS RA,Dec to Cartesian.
    pi = s2c(ri, di)
    # Bias-precession-nutation, giving GCRS proper direction.
    ppr = trxp(astrom.bpn, pi)
    # Aberration, giving GCRS natural direction.
    d = zp()
    for j in range(2):
        r2 = 0.
        for i in range(3):
            w = ppr[i] - d[i]
            before[i] = w
            r2 += w*w
        r = math.sqrt(r2)
        for i in range(3):
            before[i] /= r
        after = ab(before, astrom.v, astrom.em, astrom.bm1)
        r2 = 0
        for i in range(3):
            d[i] = after[i] - before[i]
            w = ppr[i] - d[i]
            pnat[i] = w
            r2 += w*w
        r = math.sqrt(r2)
        for i in range(3):
            pnat[i] /= r
    #Light deflection, giving BCRS coordinate direction.
    d = zp()
    for j in range(5):
        r2 = 0.
        for i in range(3):
            w = pnat[i] - d[i]
            before[i] = w
            r2 += w*w
        r = math.sqrt(r2)
        for i in range(3):
            before[i] /= r
        after = ldn(ldbody, astrom.eb, before)
        r2 = 0.0
        for i in range(3):
            d[i] = after[i] - before[i]
            w = pnat[i] - d[i]
            pco[i] = w
            r2 += w*w
        r = math.sqrt(r2)
        for i in range(3):
            pco[i] /= r
    # ICRS astrometric RA,Dec.
    w ,dc = c2s(pco)
    rc = anp(w)
    return rc, dc

__all__ = ['ASTROM', 'AULT', 'CMPS', 'D2PI', 'DAS2R', 'DAU', 'DAYSEC',
           'DC', 'DD2R', 'DJ00', 'DJC', 'DJM', 'DJM0', 'DJM00', 'DJM77',
           'DJY', 'DMAS2R', 'DPI', 'DR2AS', 'DR2D', 'DS2R', 'DTY',
           'ELB', 'ELG', 'GRS80', 'LDBODY', 'SRS', 'TDB0', 'TTMTAI',
           'TURNAS', 'WGS72', 'WGS84',
           'a2af', 'a2tf', 'ab', 'af2a', 'anp', 'anpm', 'apcg', 'apcg13',
           'apci', 'apci13', 'apco', 'apco13', 'apcs', 'apcs13', 'aper',
           'aper13', 'apio', 'apio13', 'atci13', 'atciq', 'atciqn',
           'atciqz', 'atco13', 'atic13', 'aticq', 'aticqn', 'atio13',
           'atioq', 'atoc13', 'atoi13', 'atoiq', 'bi00', 'bp00', 'bp06',
           'bpn2xy', 'c2i00a', 'c2i00b', 'c2i06a', 'c2ibpn', 'c2ixy',
           'c2ixys', 'c2s', 'c2t00a', 'c2t00b', 'c2t06a', 'c2tcio',
           'c2teqx', 'c2tpe', 'c2txy', 'cal2jd', 'cp', 'cpv', 'cr',
           'd2dtf', 'd2tf', 'dat', 'dtdb', 'dtf2d', 'ee00', 'ee00a',
           'ee00b', 'ee06a', 'eect00', 'eform', 'eo06a', 'eors', 'epb',
           'epb2jd', 'epj', 'epj2jd', 'epv00', 'eqeq94', 'era00', 'error',
           'fad03', 'fae03', 'faf03', 'faju03', 'fal03', 'falp03', 'fama03',
           'fame03', 'fane03', 'faom03', 'fapa03', 'fasa03', 'faur03', 'fave03',
           'fk52h', 'fk5hip', 'fk5hz', 'fw2m', 'fw2xy', 'gc2gd', 'gc2gde', 'gd2gc',
           'gd2gce', 'gmst00', 'gmst06', 'gmst82', 'gst00a', 'gst00b', 'gst06',
           'gst06a', 'gst94', 'h2fk5', 'hfk5z', 'ir', 'jd2cal', 'jdcalf',
           'ld', 'ldn', 'ldsun', 'math', 'num00a', 'num00b', 'num06a', 'numat',
           'nut00a', 'nut00b', 'nut06a', 'nut80', 'nutm80', 'obl06', 'obl80',
           'p06e', 'p2pv', 'p2s', 'pap', 'pas', 'pb06', 'pdp', 'pfw06', 'plan94',
           'pm', 'pmat00', 'pmat06', 'pmat76', 'pmp', 'pmpx', 'pmsafe', 'pn',
           'pn00', 'pn00a', 'pn00b', 'pn06', 'pn06a', 'pnm00a', 'pnm00b',
           'pnm06a', 'pnm80', 'pom00', 'ppp', 'ppsp', 'pr00', 'prec76',
           'pv2p', 'pv2s', 'pvdpv', 'pvm', 'pvmpv', 'pvppv', 'pvstar',
           'pvtob', 'pvu', 'pvup', 'pvxpv', 'pxp', 'refco', 'rm2v',
           'rv2m', 'rx', 'rxp', 'rxpv', 'rxr', 'ry', 'rz', 's00',
           's00a', 's00b', 's06', 's06a', 's2c', 's2p', 's2pv',
           's2xpv', 'sepp', 'seps', 'sp00', 'starpm', 'starpv',
           'sxp', 'sxpv', 'taitt', 'taiut1', 'taiutc', 'tcbtdb',
           'tcgtt', 'tdbtcb', 'tdbtt', 'tf2a', 'tf2d', 'tr', 'trxp',
           'trxpv', 'tttai', 'tttcg', 'tttdb', 'ttut1', 'ut1tai',
           'ut1tt', 'ut1utc', 'utctai', 'utcut1', 'xy06', 'xys00a',
           'xys00b', 'xys06a', 'zp', 'zpv', 'zr']
