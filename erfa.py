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
           'xys00b', 'xys06a', 'zp', 'zpv', 'zr',
           'icrs2g', 'g2icrs',
           'ltp', 'ltpb', 'ltpecl', 'ltpequ',
           'eceq06', 'ecm06', 'eqec06',
           'lteceq', 'ltecm', 'lteqec']
