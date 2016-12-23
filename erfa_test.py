import unittest
import sys
import math
import erfa

class Validate(unittest.TestCase):
    '''
     float:   self.assertAlmostEqual(value, expected, places=diff)       
     int:     self.assertEqual(value, expected)
     char:    self.assertEqual(value, expected)
    '''
## astrometry tools
    def test_ab(self):
        pnat = (-0.76321968546737951,-0.60869453983060384,-0.21676408580639883)
        v = (2.1044018893653786e-5,-8.9108923304429319e-5,-3.8633714797716569e-5)
        s = 0.99980921395708788
        bm1 = 0.99999999506209258
        ppr = erfa.ab(pnat, v, s, bm1)
        self.assertAlmostEqual(ppr[0], -0.7631631094219556269, places=12)
        self.assertAlmostEqual(ppr[1], -0.6087553082505590832, places=12)
        self.assertAlmostEqual(ppr[2], -0.2167926269368471279, places=12)

    def test_apcg(self):
        date1 = 2456165.5
        date2 = 0.401182685
        ebpv=((0.901310875,-0.417402664,-0.180982288),
              ( 0.00742727954,0.0140507459,0.00609045792))
        ehp=(0.903358544,-0.415395237,-0.180084014)
        astrom = erfa.apcg(date1,date2,ebpv,ehp)
        self.assertAlmostEqual(astrom.pmt, 12.65133794027378508, places=11)
        self.assertAlmostEqual(astrom.eb[0], 0.901310875, places=12)
        self.assertAlmostEqual(astrom.eb[1], -0.417402664, places=12,)
        self.assertAlmostEqual(astrom.eb[2], -0.180982288, places=12)
        self.assertAlmostEqual(astrom.eh[0], 0.8940025429324143045, places=12,)
        self.assertAlmostEqual(astrom.eh[1], -0.4110930268679817955, places=12)
        self.assertAlmostEqual(astrom.eh[2], -0.1782189004872870264, places=12)
        self.assertAlmostEqual(astrom.em, 1.010465295811013146, places=12)
        self.assertAlmostEqual(astrom.v[0], 0.4289638897813379954e-4, places=16)
        self.assertAlmostEqual(astrom.v[1], 0.8115034021720941898e-4, places=16)
        self.assertAlmostEqual(astrom.v[2], 0.3517555123437237778e-4, places=16)
        self.assertAlmostEqual(astrom.bm1, 0.9999999951686013336, places=15) #failed at 16
        self.assertAlmostEqual(astrom.bpn[0][0], 1.0, 10)
        self.assertAlmostEqual(astrom.bpn[1][0], 0.0, 10)
        self.assertAlmostEqual(astrom.bpn[2][0], 0.0, 10)
        self.assertAlmostEqual(astrom.bpn[0][1], 0.0, 10)
        self.assertAlmostEqual(astrom.bpn[1][1], 1.0, 10)
        self.assertAlmostEqual(astrom.bpn[2][1], 0.0, 10)
        self.assertAlmostEqual(astrom.bpn[0][2], 0.0, 10)
        self.assertAlmostEqual(astrom.bpn[1][2], 0.0, 10)
        self.assertAlmostEqual(astrom.bpn[2][2], 1.0, 10)

    def test_apcg13(self):
        date1 = 2456165.5
        date2 = 0.401182685
        astrom=erfa.apcg13(date1,date2)
        self.assertAlmostEqual(astrom.pmt, 12.65133794027378508, places=12)
        self.assertAlmostEqual(astrom.eb[0], 0.9013108747340644755, places=12)
        self.assertAlmostEqual(astrom.eb[1], -0.4174026640406119957, places=12)
        self.assertAlmostEqual(astrom.eb[2], -0.1809822877867817771, places=12)
        self.assertAlmostEqual(astrom.eh[0], 0.8940025429255499549, places=12)
        self.assertAlmostEqual(astrom.eh[1], -0.4110930268331896318, places=12)
        self.assertAlmostEqual(astrom.eh[2], -0.1782189006019749850, places=12)
        self.assertAlmostEqual(astrom.em, 1.010465295964664178, places=12)
        self.assertAlmostEqual(astrom.v[0], 0.4289638897157027528e-4, places=16)
        self.assertAlmostEqual(astrom.v[1], 0.8115034002544663526e-4, places=16)
        self.assertAlmostEqual(astrom.v[2], 0.3517555122593144633e-4, places=16)
        self.assertAlmostEqual(astrom.bm1, 0.9999999951686013498, places=12)
        self.assertAlmostEqual(astrom.bpn[0][0], 1.0, 10)
        self.assertAlmostEqual(astrom.bpn[1][0], 0.0, 10)
        self.assertAlmostEqual(astrom.bpn[2][0], 0.0, 10)
        self.assertAlmostEqual(astrom.bpn[0][1], 0.0, 10)
        self.assertAlmostEqual(astrom.bpn[1][1], 1.0, 10)
        self.assertAlmostEqual(astrom.bpn[2][1], 0.0, 10)
        self.assertAlmostEqual(astrom.bpn[0][2], 0.0, 10)
        self.assertAlmostEqual(astrom.bpn[1][2], 0.0, 10)
        self.assertAlmostEqual(astrom.bpn[2][2], 1.0, 10)

    def test_apci(self):
        date1 = 2456165.5
        date2 = 0.401182685
        ebpv = ((0.901310875,-0.417402664,-0.180982288),
                 (0.00742727954,0.0140507459,0.00609045792))
        ehp = (0.903358544,-0.415395237,-0.180084014)
        x =  0.0013122272
        y = -2.92808623e-5
        s =  3.05749468e-8
        astrom = erfa.apci(date1, date2, ebpv, ehp, x, y, s)
        self.assertAlmostEqual(astrom.pmt, 12.65133794027378508, places=11)
        self.assertAlmostEqual(astrom.eb[0], 0.901310875, places=12)
        self.assertAlmostEqual(astrom.eb[1], -0.417402664, places=12)
        self.assertAlmostEqual(astrom.eb[2], -0.180982288, places=12)
        self.assertAlmostEqual(astrom.eh[0], 0.8940025429324143045, places=12)
        self.assertAlmostEqual(astrom.eh[1], -0.4110930268679817955, places=12)
        self.assertAlmostEqual(astrom.eh[2], -0.1782189004872870264, places=12)
        self.assertAlmostEqual(astrom.em, 1.010465295811013146, places=12)
        self.assertAlmostEqual(astrom.v[0], 0.4289638897813379954e-4, places=16)
        self.assertAlmostEqual(astrom.v[1], 0.8115034021720941898e-4, places=16)
        self.assertAlmostEqual(astrom.v[2], 0.3517555123437237778e-4, places=16)
        self.assertAlmostEqual(astrom.bm1, 0.9999999951686013336, places=12)
        self.assertAlmostEqual(astrom.bpn[0][0], 0.9999991390295159156, places=12)
        self.assertAlmostEqual(astrom.bpn[1][0], 0.4978650072505016932e-7, places=12)
        self.assertAlmostEqual(astrom.bpn[2][0], 0.1312227200000000000e-2, places=12)
        self.assertAlmostEqual(astrom.bpn[0][1], -0.1136336653771609630e-7, places=12)
        self.assertAlmostEqual(astrom.bpn[1][1], 0.9999999995713154868, places=12)
        self.assertAlmostEqual(astrom.bpn[2][1], -0.2928086230000000000e-4, places=12)
        self.assertAlmostEqual(astrom.bpn[0][2], -0.1312227200895260194e-2, places=12)
        self.assertAlmostEqual(astrom.bpn[1][2], 0.2928082217872315680e-4, places=12)
        self.assertAlmostEqual(astrom.bpn[2][2], 0.9999991386008323373, places=12)

    def test_apci13(self):
        date1 = 2456165.5
        date2 = 0.401182685
        astrom, eo = erfa.apci13(date1, date2)
        self.assertAlmostEqual(astrom.pmt, 12.65133794027378508, places=11)
        self.assertAlmostEqual(astrom.eb[0], 0.9013108747340644755, places=12)
        self.assertAlmostEqual(astrom.eb[1], -0.4174026640406119957, places=12)
        self.assertAlmostEqual(astrom.eb[2], -0.1809822877867817771, places=12)
        self.assertAlmostEqual(astrom.eh[0], 0.8940025429255499549, places=12)
        self.assertAlmostEqual(astrom.eh[1], -0.4110930268331896318, places=12)
        self.assertAlmostEqual(astrom.eh[2], -0.1782189006019749850, places=12)
        self.assertAlmostEqual(astrom.em, 1.010465295964664178, places=12)
        self.assertAlmostEqual(astrom.v[0], 0.4289638897157027528e-4, places=16)
        self.assertAlmostEqual(astrom.v[1], 0.8115034002544663526e-4, places=16)
        self.assertAlmostEqual(astrom.v[2], 0.3517555122593144633e-4, places=16)
        self.assertAlmostEqual(astrom.bm1, 0.9999999951686013498, places=12)
        self.assertAlmostEqual(astrom.bpn[0][0], 0.9999992060376761710, places=12)
        self.assertAlmostEqual(astrom.bpn[1][0], 0.4124244860106037157e-7, places=12)
        self.assertAlmostEqual(astrom.bpn[2][0], 0.1260128571051709670e-2, places=12)
        self.assertAlmostEqual(astrom.bpn[0][1], -0.1282291987222130690e-7, places=12)
        self.assertAlmostEqual(astrom.bpn[1][1], 0.9999999997456835325, places=12)
        self.assertAlmostEqual(astrom.bpn[2][1], -0.2255288829420524935e-4, places=12)
        self.assertAlmostEqual(astrom.bpn[0][2], -0.1260128571661374559e-2, places=12)
        self.assertAlmostEqual(astrom.bpn[1][2], 0.2255285422953395494e-4, places=12)
        self.assertAlmostEqual(astrom.bpn[2][2], 0.9999992057833604343, places=12)
        self.assertAlmostEqual(eo, -0.2900618712657375647e-2, places=12)

    def test_apco(self):
        date1 = 2456384.5
        date2 = 0.970031644
        ebpv = ((-0.974170438,-0.211520082,-0.0917583024),
                (0.00364365824,-0.0154287319,-0.00668922024))
        ehp = (-0.973458265,-0.209215307,-0.0906996477)
        x = 0.0013122272
        y = -2.92808623e-5
        s = 3.05749468e-8
        theta = 3.14540971
        elong = -0.527800806
        phi = -1.2345856
        hm = 2738.0
        xp = 2.47230737e-7
        yp = 1.82640464e-6
        sp = -3.01974337e-11
        refa = 0.000201418779
        refb = -2.36140831e-7
        astrom = erfa.apco(date1, date2, ebpv, ehp, x, y, s,
                           theta, elong, phi, hm, xp, yp, sp,
                           refa, refb)
        self.assertAlmostEqual(astrom.pmt, 13.25248468622587269, places=11)
        self.assertAlmostEqual(astrom.eb[0], -0.9741827110630897003, places=12)
        self.assertAlmostEqual(astrom.eb[1], -0.2115130190135014340, places=12)
        self.assertAlmostEqual(astrom.eb[2], -0.09179840186968295686, places=12)
        self.assertAlmostEqual(astrom.eh[0], -0.9736425571689670428, places=12)
        self.assertAlmostEqual(astrom.eh[1], -0.2092452125848862201, places=12)
        self.assertAlmostEqual(astrom.eh[2], -0.09075578152261439954, places=12)
        self.assertAlmostEqual(astrom.em, 0.9998233241710617934, places=12)
        self.assertAlmostEqual(astrom.v[0], 0.2078704985147609823e-4, places=16)
        self.assertAlmostEqual(astrom.v[1], -0.8955360074407552709e-4, places=16)
        self.assertAlmostEqual(astrom.v[2], -0.3863338980073114703e-4, places=16)
        self.assertAlmostEqual(astrom.bm1, 0.9999999950277561600, places=12)
        self.assertAlmostEqual(astrom.bpn[0][0], 0.9999991390295159156, places=12)
        self.assertAlmostEqual(astrom.bpn[1][0], 0.4978650072505016932e-7, places=12)
        self.assertAlmostEqual(astrom.bpn[2][0], 0.1312227200000000000e-2, places=12)
        self.assertAlmostEqual(astrom.bpn[0][1], -0.1136336653771609630e-7, places=12)
        self.assertAlmostEqual(astrom.bpn[1][1], 0.9999999995713154868, places=12)
        self.assertAlmostEqual(astrom.bpn[2][1], -0.2928086230000000000e-4, places=12)
        self.assertAlmostEqual(astrom.bpn[0][2], -0.1312227200895260194e-2, places=12)
        self.assertAlmostEqual(astrom.bpn[1][2], 0.2928082217872315680e-4, places=12)
        self.assertAlmostEqual(astrom.bpn[2][2], 0.9999991386008323373, places=12)
        self.assertAlmostEqual(astrom.along, -0.5278008060301974337, places=12)
        self.assertAlmostEqual(astrom.xpl, 0.1133427418174939329e-5, places=17)
        self.assertAlmostEqual(astrom.ypl, 0.1453347595745898629e-5, places=17)
        self.assertAlmostEqual(astrom.sphi, -0.9440115679003211329, places=12)
        self.assertAlmostEqual(astrom.cphi, 0.3299123514971474711, places=12)
        self.assertAlmostEqual(astrom.diurab, 0, 10)
        self.assertAlmostEqual(astrom.eral, 2.617608903969802566, places=12)
        self.assertAlmostEqual(astrom.refa, 0.2014187790000000000e-3, places=15)
        self.assertAlmostEqual(astrom.refb, -0.2361408310000000000e-6, places=18)

    def test_apco13(self):
        utc1 = 2456384.5
        utc2 = 0.969254051
        dut1 = 0.1550675
        elong = -0.527800806
        phi = -1.2345856
        hm = 2738.0
        p = 2.47230737e-7
        xp = 2.47230737e-7
        yp = 1.82640464e-6
        phpa = 731.0
        tc = 12.8
        rh = 0.59
        wl = 0.55
        astrom, eo = erfa.apco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
        self.assertAlmostEqual(astrom.pmt, 13.25248468622475727, places=11)
        self.assertAlmostEqual(astrom.eb[0], -0.9741827107321449445, places=12)
        self.assertAlmostEqual(astrom.eb[1], -0.2115130190489386190, places=12)
        self.assertAlmostEqual(astrom.eb[2], -0.09179840189515518726, places=12)
        self.assertAlmostEqual(astrom.eh[0], -0.9736425572586866640, places=12)
        self.assertAlmostEqual(astrom.eh[1], -0.2092452121602867431, places=12)
        self.assertAlmostEqual(astrom.eh[2], -0.09075578153903832650, places=12)
        self.assertAlmostEqual(astrom.em, 0.9998233240914558422, places=12)
        self.assertAlmostEqual(astrom.v[0], 0.2078704986751370303e-4, places=16)
        self.assertAlmostEqual(astrom.v[1], -0.8955360100494469232e-4, places=16)
        self.assertAlmostEqual(astrom.v[2], -0.3863338978840051024e-4, places=16)
        self.assertAlmostEqual(astrom.bm1, 0.9999999950277561368, places=12)
        self.assertAlmostEqual(astrom.bpn[0][0], 0.9999991390295147999, places=12)
        self.assertAlmostEqual(astrom.bpn[1][0], 0.4978650075315529277e-7, places=12)
        self.assertAlmostEqual(astrom.bpn[2][0], 0.001312227200850293372, places=12)
        self.assertAlmostEqual(astrom.bpn[0][1], -0.1136336652812486604e-7, places=12)
        self.assertAlmostEqual(astrom.bpn[1][1], 0.9999999995713154865, places=12)
        self.assertAlmostEqual(astrom.bpn[2][1], -0.2928086230975367296e-4, places=12)
        self.assertAlmostEqual(astrom.bpn[0][2], -0.001312227201745553566, places=12)
        self.assertAlmostEqual(astrom.bpn[1][2], 0.2928082218847679162e-4, places=12)
        self.assertAlmostEqual(astrom.bpn[2][2], 0.9999991386008312212, places=12)
        self.assertAlmostEqual(astrom.along, -0.5278008060301974337, places=12)
        self.assertAlmostEqual(astrom.xpl, 0.1133427418174939329e-5, places=17)
        self.assertAlmostEqual(astrom.ypl, 0.1453347595745898629e-5, places=17)
        self.assertAlmostEqual(astrom.sphi, -0.9440115679003211329, places=12)
        self.assertAlmostEqual(astrom.cphi, 0.3299123514971474711, places=12)
        self.assertAlmostEqual(astrom.diurab, 0, 10)
        self.assertAlmostEqual(astrom.eral, 2.617608909189066140, places=12)
        self.assertAlmostEqual(astrom.refa, 0.2014187785940396921e-3, places=15)
        self.assertAlmostEqual(astrom.refb, -0.2361408314943696227e-6, places=18)
        self.assertAlmostEqual(eo, -0.003020548354802412839, places=14)

    def test_apcs(self):
        date1 = 2456384.5
        date2 = 0.970031644
        pv = ((-1836024.09,1056607.72,-5998795.26),
              (-77.0361767,-133.310856,0.0971855934))
        ebpv = ((-0.974170438,-0.211520082,-0.0917583024),
                (0.00364365824,-0.0154287319,-0.00668922024))
        ehp = (-0.973458265,-0.209215307,-0.0906996477)
        astrom = erfa.apcs(date1, date2, pv, ebpv, ehp)
        self.assertAlmostEqual(astrom.pmt, 13.25248468622587269, places=11)
        self.assertAlmostEqual(astrom.eb[0], -0.9741827110630456169, places=12)
        self.assertAlmostEqual(astrom.eb[1], -0.2115130190136085494, places=12)
        self.assertAlmostEqual(astrom.eb[2], -0.09179840186973175487, places=12)
        self.assertAlmostEqual(astrom.eh[0], -0.9736425571689386099, places=12)
        self.assertAlmostEqual(astrom.eh[1], -0.2092452125849967195, places=12)
        self.assertAlmostEqual(astrom.eh[2], -0.09075578152266466572, places=12)
        self.assertAlmostEqual(astrom.em, 0.9998233241710457140, places=12)
        self.assertAlmostEqual(astrom.v[0], 0.2078704985513566571e-4, places=16)
        self.assertAlmostEqual(astrom.v[1], -0.8955360074245006073e-4, places=16)
        self.assertAlmostEqual(astrom.v[2], -0.3863338980073572719e-4, places=16)
        self.assertAlmostEqual(astrom.bm1, 0.9999999950277561601, places=12)
        self.assertAlmostEqual(astrom.bpn[0][0], 1, 10)
        self.assertAlmostEqual(astrom.bpn[1][0], 0, 10)
        self.assertAlmostEqual(astrom.bpn[2][0], 0, 10)
        self.assertAlmostEqual(astrom.bpn[0][1], 0, 10)
        self.assertAlmostEqual(astrom.bpn[1][1], 1, 10)
        self.assertAlmostEqual(astrom.bpn[2][1], 0, 10)
        self.assertAlmostEqual(astrom.bpn[0][2], 0, 10)
        self.assertAlmostEqual(astrom.bpn[1][2], 0, 10)
        self.assertAlmostEqual(astrom.bpn[2][2], 1, 10)

    def test_apcs13(self):
        date1 = 2456165.5
        date2 = 0.401182685
        pv = ((-6241497.16,401346.896,-1251136.04),
              (-29.264597,-455.021831,0.0266151194))
        astrom = erfa.apcs13(date1, date2, pv)
        self.assertAlmostEqual(astrom.pmt, 12.65133794027378508, places=11)
        self.assertAlmostEqual(astrom.eb[0], 0.9012691529023298391, places=12)
        self.assertAlmostEqual(astrom.eb[1], -0.4173999812023068781, places=12)
        self.assertAlmostEqual(astrom.eb[2], -0.1809906511146821008, places=12)
        self.assertAlmostEqual(astrom.eh[0], 0.8939939101759726824, places=12)
        self.assertAlmostEqual(astrom.eh[1], -0.4111053891734599955, places=12)
        self.assertAlmostEqual(astrom.eh[2], -0.1782336880637689334, places=12)
        self.assertAlmostEqual(astrom.em, 1.010428384373318379, places=12)
        self.assertAlmostEqual(astrom.v[0], 0.4279877278327626511e-4, places=16)
        self.assertAlmostEqual(astrom.v[1], 0.7963255057040027770e-4, places=16)
        self.assertAlmostEqual(astrom.v[2], 0.3517564000441374759e-4, places=16)
        self.assertAlmostEqual(astrom.bm1, 0.9999999952947981330, places=12)
        self.assertAlmostEqual(astrom.bpn[0][0], 1, 10)
        self.assertAlmostEqual(astrom.bpn[1][0], 0, 10)
        self.assertAlmostEqual(astrom.bpn[2][0], 0, 10)
        self.assertAlmostEqual(astrom.bpn[0][1], 0, 10)
        self.assertAlmostEqual(astrom.bpn[1][1], 1, 10)
        self.assertAlmostEqual(astrom.bpn[2][1], 0, 10)
        self.assertAlmostEqual(astrom.bpn[0][2], 0, 10)
        self.assertAlmostEqual(astrom.bpn[1][2], 0, 10)
        self.assertAlmostEqual(astrom.bpn[2][2], 1, 10)

    def test_aper(self):
        theta = 5.678
        pmt = 0
        eb = (0,0,0)
        eh = (0,0,0)
        em = 0
        v = (0,0,0)
        bm1 = 0
        bpn = ((0,0,0),(0,0,0),(0,0,0))
        along = 1.234
        phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb = 0, 0, 0, 0, 0, 0, 0, 0, 0
        astrom = erfa.ASTROM((pmt, eb, eh, em, v, bm1, bpn, along,
                              phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb))
        astrom = erfa.aper(theta, astrom)
        self.assertAlmostEqual(astrom.eral, 6.912000000000000000, places=12)

    def test_aper13(self):
        ut11 = 2456165.5
        ut12 = 0.401182685
        pmt = 0
        eb = (0,0,0)
        eh = (0,0,0)
        em = 0
        v = (0,0,0)
        bm1 = 0
        bpn = ((0,0,0),(0,0,0),(0,0,0))
        along = 1.234
        phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb = 0, 0, 0, 0, 0, 0, 0, 0, 0
        astrom = erfa.ASTROM((pmt, eb, eh, em, v, bm1, bpn, along,
                              phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb))
        astrom = erfa.aper13(ut11, ut12, astrom)
        self.assertAlmostEqual(astrom.eral, 3.316236661789694933, places=12)

    def test_apio(self):
        sp = -3.01974337e-11
        theta = 3.14540971
        elong = -0.527800806
        phi = -1.2345856
        hm = 2738.0
        xp = 2.47230737e-7
        yp = 1.82640464e-6
        refa = 0.000201418779
        refb = -2.36140831e-7
        astrom = erfa.apio(sp, theta, elong, phi, hm, xp, yp, refa, refb)
        self.assertAlmostEqual(astrom.along, -0.5278008060301974337, places=12)
        self.assertAlmostEqual(astrom.xpl, 0.1133427418174939329e-5, places=17)
        self.assertAlmostEqual(astrom.ypl, 0.1453347595745898629e-5, places=17)
        self.assertAlmostEqual(astrom.sphi, -0.9440115679003211329, places=12)
        self.assertAlmostEqual(astrom.cphi, 0.3299123514971474711, places=12)
        self.assertAlmostEqual(astrom.diurab, 0.5135843661699913529e-6, places=12)
        self.assertAlmostEqual(astrom.eral, 2.617608903969802566, places=12)
        self.assertAlmostEqual(astrom.refa, 0.2014187790000000000e-3, places=15)
        self.assertAlmostEqual(astrom.refb, -0.2361408310000000000e-6, places=18)

    def test_apio13(self):
        utc1 = 2456384.5
        utc2 = 0.969254051
        dut1 = 0.1550675
        elong = -0.527800806
        phi = -1.2345856
        hm = 2738.0
        xp = 2.47230737e-7
        yp = 1.82640464e-6
        phpa = 731.0
        tc = 12.8
        rh = 0.59
        wl = 0.55
        astrom = erfa.apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
        self.assertAlmostEqual(astrom.along, -0.5278008060301974337, places=12)
        self.assertAlmostEqual(astrom.xpl, 0.1133427418174939329e-5, places=17)
        self.assertAlmostEqual(astrom.ypl, 0.1453347595745898629e-5, places=17)
        self.assertAlmostEqual(astrom.sphi, -0.9440115679003211329, places=12)
        self.assertAlmostEqual(astrom.cphi, 0.3299123514971474711, places=12)
        self.assertAlmostEqual(astrom.diurab, 0.5135843661699913529e-6, places=12)
        self.assertAlmostEqual(astrom.eral, 2.617608909189066140, places=12)
        self.assertAlmostEqual(astrom.refa, 0.2014187785940396921e-3, places=15)
        self.assertAlmostEqual(astrom.refb, -0.2361408314943696227e-6, places=18)

    def test_atci13(self):
        rc = 2.71
        dc = 0.174
        pr = 1e-5
        pd = 5e-6
        px = 0.1
        rv = 55.0
        date1 = 2456165.5
        date2 = 0.401182685
        ri, di, eo = erfa.atci13(rc, dc, pr, pd, px, rv, date1, date2)
        self.assertAlmostEqual(ri, 2.710121572969038991, places=12)
        self.assertAlmostEqual(di, 0.1729371367218230438, places=12)
        self.assertAlmostEqual(eo, -0.002900618712657375647, places=14)

    def test_atciq(self):
        date1 = 2456165.5
        date2 = 0.401182685
        astrom, eo = erfa.apci13(date1, date2)
        rc = 2.71
        dc = 0.174
        pr = 1e-5
        pd = 5e-6
        px = 0.1
        rv = 55.0
        ri, di = erfa.atciq(rc, dc, pr, pd, px, rv, *astrom)
        self.assertAlmostEqual(ri, 2.710121572969038991, places=12)
        self.assertAlmostEqual(di, 0.1729371367218230438, places=12)

    def test_atciqn(self):
        date1 = 2456165.5
        date2 = 0.401182685
        astrom, eo = erfa.apci13(date1, date2)
        rc = 2.71
        dc = 0.174
        pr = 1e-5
        pd = 5e-6
        px = 0.1
        rv = 55.0
        b0 = erfa.LDBODY((0.00028574,
                          3e-10,
                          ((-7.81014427,-5.60956681,-1.98079819),
                           (0.0030723249,-0.00406995477,-0.00181335842))))
        b1 = erfa.LDBODY((0.00095435,
                          3e-9,
                          ((0.738098796, 4.63658692,1.9693136),
                           (-0.00755816922, 0.00126913722, 0.000727999001))))
        b2 = erfa.LDBODY((1.0,
                          6e-6,
                          ((-0.000712174377, -0.00230478303, -0.00105865966),
                           (6.29235213e-6, -3.30888387e-7, -2.96486623e-7))))
        b = [b0, b1, b2]
        ri, di = erfa.atciqn(rc, dc, pr, pd, px, rv, b, *astrom)
        self.assertAlmostEqual(ri, 2.710122008105325582, places=12)
        self.assertAlmostEqual(di, 0.1729371916491459122, places=12)

    def test_atciqz(self):
        date1 = 2456165.5
        date2 = 0.401182685
        astrom, eo = erfa.apci13(date1, date2)
        rc = 2.71
        dc = 0.174
        ri, di = erfa.atciqz(rc, dc, *astrom)
        self.assertAlmostEqual(ri, 2.709994899247599271, places=12)
        self.assertAlmostEqual(di, 0.1728740720983623469, places=12)

    def test_atco13(self):
        rc = 2.71
        dc = 0.174
        pr = 1e-5
        pd = 5e-6
        px = 0.1
        rv = 55.0
        utc1 = 2456384.5
        utc2 = 0.969254051
        dut1 = 0.1550675
        elong = -0.527800806
        phi = -1.2345856
        hm = 2738.0
        xp = 2.47230737e-7
        yp = 1.82640464e-6
        phpa = 731.0
        tc = 12.8
        rh = 0.59
        wl = 0.55
        aob, zob, hob, dob, rob, eo = erfa.atco13(rc, dc, pr, pd, px, rv,
                                                  utc1, utc2, dut1, elong, phi, hm, xp, yp,
                                                  phpa, tc, rh, wl)
        self.assertAlmostEqual(aob, 0.09251774485358230653, places=12)
        self.assertAlmostEqual(zob, 1.407661405256767021, places=12)
        self.assertAlmostEqual(hob, -0.09265154431403157925, places=12)
        self.assertAlmostEqual(dob, 0.1716626560075591655, places=12)
        self.assertAlmostEqual(rob, 2.710260453503097719, places=12)
        self.assertAlmostEqual(eo, -0.003020548354802412839, places=14)

    def test_atic13(self):
        ri = 2.710121572969038991
        di = 0.1729371367218230438
        date1 = 2456165.5
        date2 = 0.401182685
        rc, dc, eo = erfa.atic13(ri, di, date1, date2)
        self.assertAlmostEqual(rc, 2.710126504531374930, places=12)
        self.assertAlmostEqual(dc, 0.1740632537628342320, places=12)
        self.assertAlmostEqual(eo, -0.002900618712657375647, places=14)

    def test_aticq(self):
        date1 = 2456165.5
        date2 = 0.401182685
        astrom, eo = erfa.apci13(date1, date2)
        ri = 2.710121572969038991
        di = 0.1729371367218230438
        rc, dc = erfa.aticq(ri, di, *astrom)
        self.assertAlmostEqual(rc, 2.710126504531374930, places=12)
        self.assertAlmostEqual(dc, 0.1740632537628342320, places=12)

    def test_aticqn(self):
        date1 = 2456165.5
        date2 = 0.401182685
        astrom, eo = erfa.apci13(date1, date2)
        ri = 2.709994899247599271
        di = 0.1728740720983623469
        b0 = erfa.LDBODY((0.00028574,
                          3e-10,
                          ((-7.81014427,-5.60956681,-1.98079819),
                           (0.0030723249,-0.00406995477,-0.00181335842))))
        b1 = erfa.LDBODY((0.00095435,
                          3e-9,
                          ((0.738098796, 4.63658692,1.9693136),
                           (-0.00755816922, 0.00126913722, 0.000727999001))))
        b2 = erfa.LDBODY((1.0,
                          6e-6,
                          ((-0.000712174377, -0.00230478303, -0.00105865966),
                           (6.29235213e-6, -3.30888387e-7, -2.96486623e-7))))
        l = [b0, b1, b2]
        rc, dc = erfa.aticqn(ri, di, l, *astrom)
        self.assertAlmostEqual(rc, 2.709999575032685412, places=12)
        self.assertAlmostEqual(dc, 0.1739999656317778034, places=12)

    def test_atio13(self):
        ri = 2.710121572969038991
        di = 0.1729371367218230438
        utc1 = 2456384.5
        utc2 = 0.969254051
        dut1 = 0.1550675
        elong = -0.527800806
        phi = -1.2345856
        hm = 2738.0
        xp = 2.47230737e-7
        yp = 1.82640464e-6
        phpa = 731.0
        tc = 12.8
        rh = 0.59
        wl = 0.55
        aob, zob, hob, dob, rob = erfa.atio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
        self.assertAlmostEqual(aob, 0.09233952224794989993, places=12)
        self.assertAlmostEqual(zob, 1.407758704513722461, places=12)
        self.assertAlmostEqual(hob, -0.09247619879782006106, places=12)
        self.assertAlmostEqual(dob, 0.1717653435758265198, places=12)
        self.assertAlmostEqual(rob, 2.710085107986886201, places=12)

    def test_atioq(self):
        utc1 = 2456384.5
        utc2 = 0.969254051
        dut1 = 0.1550675
        elong = -0.527800806
        phi = -1.2345856
        hm = 2738.0
        xp = 2.47230737e-7
        yp = 1.82640464e-6
        phpa = 731.0
        tc = 12.8
        rh = 0.59
        wl = 0.55
        astrom = erfa.apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
                             phpa, tc, rh, wl)
        ri = 2.710121572969038991
        di = 0.1729371367218230438
        aob, zob, hob, dob, rob = erfa.atioq(ri, di, *astrom)
        self.assertAlmostEqual(aob, 0.09233952224794989993, places=12)
        self.assertAlmostEqual(zob, 1.407758704513722461, places=12)
        self.assertAlmostEqual(hob, -0.09247619879782006106, places=12)
        self.assertAlmostEqual(dob, 0.1717653435758265198, places=12)
        self.assertAlmostEqual(rob, 2.710085107986886201, places=12)

    def test_atoc13(self):
        utc1 = 2456384.5
        utc2 = 0.969254051
        dut1 = 0.1550675
        elong = -0.527800806
        phi = -1.2345856
        hm = 2738.0
        xp = 2.47230737e-7
        yp = 1.82640464e-6
        phpa = 731.0
        tc = 12.8
        rh = 0.59
        wl = 0.55
        ob1 = 2.710085107986886201
        ob2 = 0.1717653435758265198
        rc, dc = erfa.atoc13("R", ob1, ob2, utc1, utc2, dut1,
                             elong, phi, hm, xp, yp, phpa, tc, rh, wl)
        self.assertAlmostEqual(rc, 2.709956744661000609, places=12)
        self.assertAlmostEqual(dc, 0.1741696500895398562, places=12)
        ob1 = -0.09247619879782006106
        ob2 = 0.1717653435758265198
        rc, dc = erfa.atoc13("H", ob1, ob2, utc1, utc2, dut1,
                             elong, phi, hm, xp, yp, phpa, tc, rh, wl)
        self.assertAlmostEqual(rc, 2.709956744661000609, places=12)
        self.assertAlmostEqual(dc, 0.1741696500895398562, places=12)
        ob1 = 0.09233952224794989993
        ob2 = 1.407758704513722461
        rc, dc = erfa.atoc13("A", ob1, ob2, utc1, utc2, dut1,
                             elong, phi, hm, xp, yp, phpa, tc, rh, wl)
        self.assertAlmostEqual(rc, 2.709956744661000609, places=12)
        self.assertAlmostEqual(dc, 0.1741696500895398562, places=12)

    def test_atoi13(self):
        utc1 = 2456384.5
        utc2 = 0.969254051
        dut1 = 0.1550675
        elong = -0.527800806
        phi = -1.2345856
        hm = 2738.0
        xp = 2.47230737e-7
        yp = 1.82640464e-6
        phpa = 731.0
        tc = 12.8
        rh = 0.59
        wl = 0.55
        ob1 = 2.710085107986886201
        ob2 = 0.1717653435758265198
        ri, di = erfa.atoi13("R", ob1, ob2, utc1, utc2, dut1,
                             elong, phi, hm, xp, yp, phpa, tc, rh, wl)       
        self.assertAlmostEqual(ri, 2.710121574449135955, places=12)
        self.assertAlmostEqual(di, 0.1729371839114567725, places=12)
        ob1 = -0.09247619879782006106
        ob2 = 0.1717653435758265198
        ri, di = erfa.atoi13("H", ob1, ob2, utc1, utc2, dut1,
                             elong, phi, hm, xp, yp, phpa, tc, rh, wl)       
        self.assertAlmostEqual(ri, 2.710121574449135955, places=12)
        self.assertAlmostEqual(di, 0.1729371839114567725, places=12)
        ob1 = 0.09233952224794989993
        ob2 = 1.407758704513722461
        ri, di = erfa.atoi13("A", ob1, ob2, utc1, utc2, dut1,
                             elong, phi, hm, xp, yp, phpa, tc, rh, wl)       
        self.assertAlmostEqual(ri, 2.710121574449135955, places=12)
        self.assertAlmostEqual(di, 0.1729371839114567728, places=12)

    def test_atoiq(self):
        utc1 = 2456384.5
        utc2 = 0.969254051
        dut1 = 0.1550675
        elong = -0.527800806
        phi = -1.2345856
        hm = 2738.0
        xp = 2.47230737e-7
        yp = 1.82640464e-6
        phpa = 731.0
        tc = 12.8
        rh = 0.59
        wl = 0.55
        astrom = erfa.apio13(utc1, utc2, dut1,
                             elong, phi, hm, xp, yp, phpa, tc, rh, wl)
        ob1 = 2.710085107986886201
        ob2 = 0.1717653435758265198
        ri, di = erfa.atoiq("R", ob1, ob2, *astrom)
        self.assertAlmostEqual(ri, 2.710121574449135955, places=12)
        self.assertAlmostEqual(di, 0.1729371839114567725, places=12)
        ob1 = -0.09247619879782006106
        ob2 = 0.1717653435758265198
        ri, di = erfa.atoiq("H", ob1, ob2, *astrom)
        self.assertAlmostEqual(ri, 2.710121574449135955, places=12)
        self.assertAlmostEqual(di, 0.1729371839114567725, places=12)
        ob1 = 0.09233952224794989993
        ob2 = 1.407758704513722461
        ri, di = erfa.atoiq("A", ob1, ob2, *astrom)
        self.assertAlmostEqual(ri, 2.710121574449135955, places=12)
        self.assertAlmostEqual(di, 0.1729371839114567728, places=12)

    def test_ld(self):
        bm = 0.00028574
        p = (-0.763276255, -0.608633767, -0.216735543)
        q = (-0.763276255, -0.608633767, -0.216735543)
        e = (0.76700421, 0.605629598, 0.211937094)
        em = 8.91276983
        dlim = 3e-10
        p1 = erfa.ld(bm, p, q, e, em, dlim)
        self.assertAlmostEqual(p1[0], -0.7632762548968159627, places=12)
        self.assertAlmostEqual(p1[1], -0.6086337670823762701, places=12)
        self.assertAlmostEqual(p1[2], -0.2167355431320546947, places=12)

    def test_ldn(self):
        ob = (-0.974170437, -0.2115201, -0.0917583114)
        sc = (-0.763276255, -0.608633767, -0.216735543)
        b0 = erfa.LDBODY((0.00028574,
                          3e-10,
                          ((-7.81014427,-5.60956681,-1.98079819),
                           (0.0030723249,-0.00406995477,-0.00181335842))))
        b1 = erfa.LDBODY((0.00095435,
                          3e-9,
                          ((0.738098796, 4.63658692,1.9693136),
                           (-0.00755816922, 0.00126913722, 0.000727999001))))
        b2 = erfa.LDBODY((1.0,
                          6e-6,
                          ((-0.000712174377, -0.00230478303, -0.00105865966),
                           (6.29235213e-6, -3.30888387e-7, -2.96486623e-7))))
        l = [b0, b1, b2]
        sn = erfa.ldn(l, ob, sc)
        self.assertAlmostEqual(sn[0], -0.7632762579693333866, places=12)
        self.assertAlmostEqual(sn[1], -0.6086337636093002660, places=12)
        self.assertAlmostEqual(sn[2], -0.2167355420646328159, places=12)

    def test_ldsun(self):
        p = (-0.763276255, -0.608633767, -0.216735543)
        e = (-0.973644023, -0.20925523, -0.0907169552)
        em = 0.999809214
        p1 = erfa.ldsun(p, e, em)
        self.assertAlmostEqual(p1[0], -0.7632762580731413169, places=12)
        self.assertAlmostEqual(p1[1], -0.6086337635262647900, places=12)
        self.assertAlmostEqual(p1[2], -0.2167355419322321302, places=12)

    def test_pmpx(self):
        rc = 1.234
        dc = 0.789
        pr = 1e-5
        pd = -2e-5
        px = 1e-2
        rv = 10.0
        pmt = 8.75
        pob = (0.9, 0.4, 0.1)
        pco = erfa.pmpx(rc, dc, pr, pd, px, rv, pmt, pob)
        self.assertAlmostEqual(pco[0], 0.2328137623960308440, places=12)
        self.assertAlmostEqual(pco[1], 0.6651097085397855317, places=12)
        self.assertAlmostEqual(pco[2], 0.7095257765896359847, places=12)

    def test_pmsafe(self):
        ra1 = 1.234
        dec1 = 0.789
        pmr1 = 1e-5
        pmd1 = -2e-5
        px1 = 1e-2
        rv1 = 10.0
        ep1a = 2400000.5
        ep1b = 48348.5625
        ep2a = 2400000.5
        ep2b = 51544.5
        ra2, dec2, pmr2, pmd2, px2, rv2 = erfa.pmsafe(ra1, dec1, pmr1, pmd1, px1,
                                                      rv1, ep1a, ep1b, ep2a, ep2b)
        self.assertAlmostEqual(ra2, 1.234087484501017061, places=12)
        self.assertAlmostEqual(dec2, 0.7888249982450468574, places=12)
        self.assertAlmostEqual(pmr2, 0.9996457663586073988e-5, places=12)
        self.assertAlmostEqual(pmd2, -0.2000040085106737816e-4, places=16)
        self.assertAlmostEqual(px2, 0.9999997295356765185e-2, places=12)
        self.assertAlmostEqual(rv2, 10.38468380113917014, places=10)

    def test_pvtob(self):
        elong = 2.0
        phi = 0.5
        hm = 3000.0
        xp = 1e-6
        yp = -0.5e-6
        sp = 1e-8
        theta = 5.0
        pv = erfa.pvtob(elong, phi, hm, xp, yp, sp, theta)
        self.assertAlmostEqual(pv[0][0], 4225081.367071159207, places=5)
        self.assertAlmostEqual(pv[0][1], 3681943.215856198144, places=5)
        self.assertAlmostEqual(pv[0][2], 3041149.399241260785, places=5)
        self.assertAlmostEqual(pv[1][0], -268.4915389365998787, places=9)
        self.assertAlmostEqual(pv[1][1], 308.0977983288903123, places=9)
        self.assertAlmostEqual(pv[1][2], 0, 0)

    def test_refco(self):
        phpa = 800.0
        tc = 10.0
        rh = 0.9
        wl = 0.4
        refa, refb = erfa.refco(phpa, tc, rh, wl)
        self.assertAlmostEqual(refa, 0.2264949956241415009e-3, places=15)
        self.assertAlmostEqual(refb, -0.2598658261729343970e-6, places=18)

## astronomy library
    def test_bi00(self):
        dp, de, dr = erfa.bi00()
        self.assertAlmostEqual(dp, -0.2025309152835086613e-6, places=15)
        self.assertAlmostEqual(de, -0.3306041454222147847e-7, places=15)
        self.assertAlmostEqual(dr, -0.7078279744199225506e-7, places=15)

    def test_bp00(self):
        rb, rp,rbp = erfa.bp00(2400000.5, 50123.9999)
        self.assertAlmostEqual(rb[0][0], 0.9999999999999942498, places=12)
        self.assertAlmostEqual(rb[0][1], -0.7078279744199196626e-7, places=16)
        self.assertAlmostEqual(rb[0][2], 0.8056217146976134152e-7, places=16)
        self.assertAlmostEqual(rb[1][0], 0.7078279477857337206e-7, places=16)
        self.assertAlmostEqual(rb[1][1], 0.9999999999999969484, places=12)
        self.assertAlmostEqual(rb[1][2], 0.3306041454222136517e-7, places=16)
        self.assertAlmostEqual(rb[2][0], -0.8056217380986972157e-7, places=16)
        self.assertAlmostEqual(rb[2][1], -0.3306040883980552500e-7, places=16)
        self.assertAlmostEqual(rb[2][2], 0.9999999999999962084, places=12)
        self.assertAlmostEqual(rp[0][0], 0.9999995504864048241, places=12)
        self.assertAlmostEqual(rp[0][1], 0.8696113836207084411e-3, places=14)
        self.assertAlmostEqual(rp[0][2], 0.3778928813389333402e-3, places=14)
        self.assertAlmostEqual(rp[1][0], -0.8696113818227265968e-3, places=14)
        self.assertAlmostEqual(rp[1][1], 0.9999996218879365258, places=12)
        self.assertAlmostEqual(rp[1][2], -0.1690679263009242066e-6, places=14)
        self.assertAlmostEqual(rp[2][0], -0.3778928854764695214e-3, places=14)
        self.assertAlmostEqual(rp[2][1], -0.1595521004195286491e-6, places=14)
        self.assertAlmostEqual(rp[2][2], 0.9999999285984682756, places=12)
        self.assertAlmostEqual(rbp[0][0], 0.9999995505175087260, places=12)
        self.assertAlmostEqual(rbp[0][1], 0.8695405883617884705e-3, places=14)
        self.assertAlmostEqual(rbp[0][2], 0.3779734722239007105e-3, places=14)
        self.assertAlmostEqual(rbp[1][0], -0.8695405990410863719e-3, places=14)
        self.assertAlmostEqual(rbp[1][1], 0.9999996219494925900, places=12)
        self.assertAlmostEqual(rbp[1][2], -0.1360775820404982209e-6, places=14)
        self.assertAlmostEqual(rbp[2][0], -0.3779734476558184991e-3, places=14)
        self.assertAlmostEqual(rbp[2][1], -0.1925857585832024058e-6, places=14)
        self.assertAlmostEqual(rbp[2][2], 0.9999999285680153377, places=12)

    def test_bp06(self):
        rb, rp,rbp = erfa.bp06(2400000.5, 50123.9999)
        self.assertAlmostEqual(rb[0][0], 0.9999999999999942497, places=12)
        self.assertAlmostEqual(rb[0][1], -0.7078368960971557145e-7, places=14)
        self.assertAlmostEqual(rb[0][2], 0.8056213977613185606e-7, places=14)
        self.assertAlmostEqual(rb[1][0], 0.7078368694637674333e-7, places=14)
        self.assertAlmostEqual(rb[1][1], 0.9999999999999969484, places=12)
        self.assertAlmostEqual(rb[1][2], 0.3305943742989134124e-7, places=14)
        self.assertAlmostEqual(rb[2][0], -0.8056214211620056792e-7, places=14)
        self.assertAlmostEqual(rb[2][1], -0.3305943172740586950e-7, places=14)
        self.assertAlmostEqual(rb[2][2], 0.9999999999999962084, places=12)

        self.assertAlmostEqual(rp[0][0], 0.9999995504864960278, places=12)
        self.assertAlmostEqual(rp[0][1], 0.8696112578855404832e-3, places=14)
        self.assertAlmostEqual(rp[0][2], 0.3778929293341390127e-3, places=14)
        self.assertAlmostEqual(rp[1][0], -0.8696112560510186244e-3, places=14)
        self.assertAlmostEqual(rp[1][1], 0.9999996218880458820, places=12)
        self.assertAlmostEqual(rp[1][2], -0.1691646168941896285e-6, places=14)
        self.assertAlmostEqual(rp[2][0], -0.3778929335557603418e-3, places=14)
        self.assertAlmostEqual(rp[2][1], -0.1594554040786495076e-6, places=14)
        self.assertAlmostEqual(rp[2][2], 0.9999999285984501222, places=12)

        self.assertAlmostEqual(rbp[0][0], 0.9999995505176007047, places=12)
        self.assertAlmostEqual(rbp[0][1], 0.8695404617348208406e-3, places=14)
        self.assertAlmostEqual(rbp[0][2], 0.3779735201865589104e-3, places=14)
        self.assertAlmostEqual(rbp[1][0], -0.8695404723772031414e-3, places=14)
        self.assertAlmostEqual(rbp[1][1], 0.9999996219496027161, places=12)
        self.assertAlmostEqual(rbp[1][2], -0.1361752497080270143e-6, places=14)
        self.assertAlmostEqual(rbp[2][0], -0.3779734957034089490e-3, places=14)
        self.assertAlmostEqual(rbp[2][1], -0.1924880847894457113e-6, places=14)
        self.assertAlmostEqual(rbp[2][2], 0.9999999285679971958, places=12)

    def test_bpn2xy(self):
        rbpn = ((9.999962358680738e-1,-2.516417057665452e-3,-1.093569785342370e-3),
                (2.516462370370876e-3,9.999968329010883e-1,4.006159587358310e-5),
                (1.093465510215479e-3,-4.281337229063151e-5,9.999994012499173e-1))
        x, y = erfa.bpn2xy(rbpn)
        self.assertAlmostEqual(x, 1.093465510215479e-3, places=12)
        self.assertAlmostEqual(y, -4.281337229063151e-5, places=12)

    def test_c2i00a(self):
        rc2i = erfa.c2i00a(2400000.5, 53736.0)
        self.assertAlmostEqual(rc2i[0][0], 0.9999998323037165557, places=12)
        self.assertAlmostEqual(rc2i[0][1], 0.5581526348992140183e-9, places=12)
        self.assertAlmostEqual(rc2i[0][2], -0.5791308477073443415e-3, places=12)

        self.assertAlmostEqual(rc2i[1][0], -0.2384266227870752452e-7, places=12)
        self.assertAlmostEqual(rc2i[1][1], 0.9999999991917405258, places=12)
        self.assertAlmostEqual(rc2i[1][2], -0.4020594955028209745e-4, places=12)

        self.assertAlmostEqual(rc2i[2][0], 0.5791308472168152904e-3, places=12)
        self.assertAlmostEqual(rc2i[2][1], 0.4020595661591500259e-4, places=12)
        self.assertAlmostEqual(rc2i[2][2], 0.9999998314954572304, places=12)

    def test_c2i00b(self):
        rc2i = erfa.c2i00b(2400000.5, 53736.0)
        self.assertAlmostEqual(rc2i[0][0], 0.9999998323040954356, places=12)
        self.assertAlmostEqual(rc2i[0][1], 0.5581526349131823372e-9, places=12)
        self.assertAlmostEqual(rc2i[0][2], -0.5791301934855394005e-3, places=12)

        self.assertAlmostEqual(rc2i[1][0], -0.2384239285499175543e-7, places=12)
        self.assertAlmostEqual(rc2i[1][1], 0.9999999991917574043, places=12)
        self.assertAlmostEqual(rc2i[1][2], -0.4020552974819030066e-4, places=12)

        self.assertAlmostEqual(rc2i[2][0], 0.5791301929950208873e-3, places=12)
        self.assertAlmostEqual(rc2i[2][1], 0.4020553681373720832e-4, places=12)
        self.assertAlmostEqual(rc2i[2][2], 0.9999998314958529887, places=12)

    def test_c2i06a(self):
        rc2i = erfa.c2i06a(2400000.5, 53736.0)
        self.assertAlmostEqual(rc2i[0][0], 0.9999998323037159379, places=12)
        self.assertAlmostEqual(rc2i[0][1], 0.5581121329587613787e-9, places=12)
        self.assertAlmostEqual(rc2i[0][2], -0.5791308487740529749e-3, places=12)

        self.assertAlmostEqual(rc2i[1][0], -0.2384253169452306581e-7, places=12)
        self.assertAlmostEqual(rc2i[1][1], 0.9999999991917467827, places=12)
        self.assertAlmostEqual(rc2i[1][2], -0.4020579392895682558e-4, places=12)

        self.assertAlmostEqual(rc2i[2][0], 0.5791308482835292617e-3, places=12)
        self.assertAlmostEqual(rc2i[2][1], 0.4020580099454020310e-4, places=12)
        self.assertAlmostEqual(rc2i[2][2], 0.9999998314954628695, places=12)

    def test_c2ibpn(self):
        rbpn = ((9.999962358680738e-1,-2.516417057665452e-3,-1.093569785342370e-3),
                (2.516462370370876e-3,9.999968329010883e-1,4.006159587358310e-5),
                (1.093465510215479e-3,-4.281337229063151e-5,9.999994012499173e-1))
        rc2i=erfa.c2ibpn(2400000.5, 50123.9999, rbpn)
        self.assertAlmostEqual(rc2i[0][0], 0.9999994021664089977, places=12)
        self.assertAlmostEqual(rc2i[0][1], -0.3869195948017503664e-8, places=12)
        self.assertAlmostEqual(rc2i[0][2], -0.1093465511383285076e-2, places=12)

        self.assertAlmostEqual(rc2i[1][0], 0.5068413965715446111e-7, places=12)
        self.assertAlmostEqual(rc2i[1][1], 0.9999999990835075686, places=12)
        self.assertAlmostEqual(rc2i[1][2], 0.4281334246452708915e-4, places=12)

        self.assertAlmostEqual(rc2i[2][0], 0.1093465510215479000e-2, places=12)
        self.assertAlmostEqual(rc2i[2][1], -0.4281337229063151000e-4, places=12)
        self.assertAlmostEqual(rc2i[2][2], 0.9999994012499173103, places=12)

    def test_c2ixy(self):
        x = 0.5791308486706011000e-3
        y = 0.4020579816732961219e-4
        rc2i = erfa.c2ixy(2400000.5, 53736, x, y)
        self.assertAlmostEqual(rc2i[0][0], 0.9999998323037157138, places=12)
        self.assertAlmostEqual(rc2i[0][1], 0.5581526349032241205e-9, places=12)
        self.assertAlmostEqual(rc2i[0][2], -0.5791308491611263745e-3, places=12)

        self.assertAlmostEqual(rc2i[1][0], -0.2384257057469842953e-7, places=12)
        self.assertAlmostEqual(rc2i[1][1], 0.9999999991917468964, places=12)
        self.assertAlmostEqual(rc2i[1][2], -0.4020579110172324363e-4, places=12)

        self.assertAlmostEqual(rc2i[2][0], 0.5791308486706011000e-3, places=12)
        self.assertAlmostEqual(rc2i[2][1], 0.4020579816732961219e-4, places=12)
        self.assertAlmostEqual(rc2i[2][2], 0.9999998314954627590, places=12)

    def test_c2ixys(self):
        x =  0.5791308486706011000e-3
        y =  0.4020579816732961219e-4
        s = -0.1220040848472271978e-7
        rc2i = erfa.c2ixys(x, y, s)
        self.assertAlmostEqual(rc2i[0][0], 0.9999998323037157138, places=12)
        self.assertAlmostEqual(rc2i[0][1], 0.5581984869168499149e-9, places=12)
        self.assertAlmostEqual(rc2i[0][2], -0.5791308491611282180e-3, places=12)

        self.assertAlmostEqual(rc2i[1][0], -0.2384261642670440317e-7, places=12)
        self.assertAlmostEqual(rc2i[1][1], 0.9999999991917468964, places=12)
        self.assertAlmostEqual(rc2i[1][2], -0.4020579110169668931e-4, places=12)

        self.assertAlmostEqual(rc2i[2][0], 0.5791308486706011000e-3, places=12)
        self.assertAlmostEqual(rc2i[2][1], 0.4020579816732961219e-4, places=12)
        self.assertAlmostEqual(rc2i[2][2], 0.9999998314954627590, places=12)

    def test_c2t00a(self):
        tta = 2400000.5
        uta = 2400000.5
        ttb = 53736.0
        utb = 53736.0
        xp = 2.55060238e-7
        yp = 1.860359247e-6
        rc2t = erfa.c2t00a(tta, ttb, uta, utb, xp, yp)
        self.assertAlmostEqual(rc2t[0][0], -0.1810332128307182668, places=12)
        self.assertAlmostEqual(rc2t[0][1], 0.9834769806938457836, places=12)
        self.assertAlmostEqual(rc2t[0][2], 0.6555535638688341725e-4, places=12)

        self.assertAlmostEqual(rc2t[1][0], -0.9834768134135984552, places=12)
        self.assertAlmostEqual(rc2t[1][1], -0.1810332203649520727, places=12)
        self.assertAlmostEqual(rc2t[1][2], 0.5749801116141056317e-3, places=12)

        self.assertAlmostEqual(rc2t[2][0], 0.5773474014081406921e-3, places=12)
        self.assertAlmostEqual(rc2t[2][1], 0.3961832391770163647e-4, places=12)
        self.assertAlmostEqual(rc2t[2][2], 0.9999998325501692289, places=12)
        
    def test_c2t00b(self):
        tta = 2400000.5
        uta = 2400000.5
        ttb = 53736.0
        utb = 53736.0
        xp = 2.55060238e-7
        yp = 1.860359247e-6
        rc2t = erfa.c2t00b(tta, ttb, uta, utb, xp, yp)
        self.assertAlmostEqual(rc2t[0][0], -0.1810332128439678965, places=12)
        self.assertAlmostEqual(rc2t[0][1], 0.9834769806913872359, places=12)
        self.assertAlmostEqual(rc2t[0][2], 0.6555565082458415611e-4, places=12)

        self.assertAlmostEqual(rc2t[1][0], -0.9834768134115435923, places=12)
        self.assertAlmostEqual(rc2t[1][1], -0.1810332203784001946, places=12)
        self.assertAlmostEqual(rc2t[1][2], 0.5749793922030017230e-3, places=12)

        self.assertAlmostEqual(rc2t[2][0], 0.5773467471863534901e-3, places=12)
        self.assertAlmostEqual(rc2t[2][1], 0.3961790411549945020e-4, places=12)
        self.assertAlmostEqual(rc2t[2][2], 0.9999998325505635738, places=12)

    def test_c2t06a(self):
        tta = 2400000.5
        uta = 2400000.5
        ttb = 53736.0
        utb = 53736.0
        xp = 2.55060238e-7
        yp = 1.860359247e-6
        rc2t = erfa.c2t06a(tta, ttb, uta, utb, xp, yp)
        self.assertAlmostEqual(rc2t[0][0], -0.1810332128305897282, places=12)
        self.assertAlmostEqual(rc2t[0][1], 0.9834769806938592296, places=12)
        self.assertAlmostEqual(rc2t[0][2], 0.6555550962998436505e-4, places=12)

        self.assertAlmostEqual(rc2t[1][0], -0.9834768134136214897, places=12)
        self.assertAlmostEqual(rc2t[1][1], -0.1810332203649130832, places=12)
        self.assertAlmostEqual(rc2t[1][2], 0.5749800844905594110e-3, places=12)

        self.assertAlmostEqual(rc2t[2][0], 0.5773474024748545878e-3, places=12)
        self.assertAlmostEqual(rc2t[2][1], 0.3961816829632690581e-4, places=12)
        self.assertAlmostEqual(rc2t[2][2], 0.9999998325501747785, places=12)

    def test_c2tcio(self):
        c = ((0.9999998323037164738,0.5581526271714303683e-9,-0.5791308477073443903e-3),
             (-0.2384266227524722273e-7,0.9999999991917404296,-0.4020594955030704125e-4),
             (0.5791308472168153320e-3,.4020595661593994396e-4,0.9999998314954572365))
        era = 1.75283325530307
        p = ((0.9999999999999674705,-0.1367174580728847031e-10,0.2550602379999972723e-6),
             (0.1414624947957029721e-10,0.9999999999982694954,-0.1860359246998866338e-5),
             (-0.2550602379741215275e-6,0.1860359247002413923e-5,0.9999999999982369658))
        rc2t = erfa.c2tcio(c,era,p)
        self.assertAlmostEqual(rc2t[0][0], -0.1810332128307110439, places=12)
        self.assertAlmostEqual(rc2t[0][1], 0.9834769806938470149, places=12)
        self.assertAlmostEqual(rc2t[0][2], 0.6555535638685466874e-4, places=12)

        self.assertAlmostEqual(rc2t[1][0], -0.9834768134135996657, places=12)
        self.assertAlmostEqual(rc2t[1][1], -0.1810332203649448367, places=12)
        self.assertAlmostEqual(rc2t[1][2], 0.5749801116141106528e-3, places=12)

        self.assertAlmostEqual(rc2t[2][0], 0.5773474014081407076e-3, places=12)
        self.assertAlmostEqual(rc2t[2][1], 0.3961832391772658944e-4, places=12)
        self.assertAlmostEqual(rc2t[2][2], 0.9999998325501691969, places=12)

    def test_c2teqx(self):
        c = ((0.9999989440476103608,-0.1332881761240011518e-2,-0.5790767434730085097e-3),
             (0.1332858254308954453e-2,0.9999991109044505944,-0.4097782710401555759e-4),
             (0.5791308472168153320e-3,0.4020595661593994396e-4,0.9999998314954572365))
        gst = 1.754166138040730516
        p = ((0.9999999999999674705,-0.1367174580728847031e-10,0.2550602379999972723e-6),
             (0.1414624947957029721e-10,0.9999999999982694954,-0.1860359246998866338e-5),
             (-0.2550602379741215275e-6,0.1860359247002413923e-5,0.9999999999982369658))
        rc2t = erfa.c2teqx(c,gst,p)
        self.assertAlmostEqual(rc2t[0][0], -0.1810332128528685730, places=12)
        self.assertAlmostEqual(rc2t[0][1], 0.9834769806897685071, places=12)
        self.assertAlmostEqual(rc2t[0][2], 0.6555535639982634449e-4, places=12)

        self.assertAlmostEqual(rc2t[1][0], -0.9834768134095211257, places=12)
        self.assertAlmostEqual(rc2t[1][1], -0.1810332203871023800, places=12)
        self.assertAlmostEqual(rc2t[1][2], 0.5749801116126438962e-3, places=12)

        self.assertAlmostEqual(rc2t[2][0], 0.5773474014081539467e-3, places=12)
        self.assertAlmostEqual(rc2t[2][1], 0.3961832391768640871e-4, places=12)
        self.assertAlmostEqual(rc2t[2][2], 0.9999998325501691969, places=12)

    def test_c2tpe(self):
        tta = 2400000.5
        uta = 2400000.5
        ttb = 53736.0
        utb = 53736.0
        deps =  0.4090789763356509900
        dpsi = -0.9630909107115582393e-5
        xp = 2.55060238e-7
        yp = 1.860359247e-6
        rc2t = erfa.c2tpe(tta, ttb, uta, utb, dpsi, deps, xp, yp)
        self.assertAlmostEqual(rc2t[0][0], -0.1813677995763029394, places=12)
        self.assertAlmostEqual(rc2t[0][1], 0.9023482206891683275, places=12)
        self.assertAlmostEqual(rc2t[0][2], -0.3909902938641085751, places=12)

        self.assertAlmostEqual(rc2t[1][0], -0.9834147641476804807, places=12)
        self.assertAlmostEqual(rc2t[1][1], -0.1659883635434995121, places=12)
        self.assertAlmostEqual(rc2t[1][2], 0.7309763898042819705e-1, places=12)

        self.assertAlmostEqual(rc2t[2][0], 0.1059685430673215247e-2, places=12)
        self.assertAlmostEqual(rc2t[2][1], 0.3977631855605078674, places=12)
        self.assertAlmostEqual(rc2t[2][2], 0.9174875068792735362, places=12)

    def test_c2txy(self):
        tta = 2400000.5
        uta = 2400000.5
        ttb = 53736.0
        utb = 53736.0
        x = 0.5791308486706011000e-3
        y = 0.4020579816732961219e-4
        xp = 2.55060238e-7
        yp = 1.860359247e-6;
        rc2t = erfa.c2txy(tta, ttb, uta, utb, x, y, xp, yp)
        self.assertAlmostEqual(rc2t[0][0], -0.1810332128306279253, places=12)
        self.assertAlmostEqual(rc2t[0][1], 0.9834769806938520084, places=12)
        self.assertAlmostEqual(rc2t[0][2], 0.6555551248057665829e-4, places=12)

        self.assertAlmostEqual(rc2t[1][0], -0.9834768134136142314, places=12)
        self.assertAlmostEqual(rc2t[1][1], -0.1810332203649529312, places=12)
        self.assertAlmostEqual(rc2t[1][2], 0.5749800843594139912e-3, places=12)

        self.assertAlmostEqual(rc2t[2][0], 0.5773474028619264494e-3, places=12)
        self.assertAlmostEqual(rc2t[2][1], 0.3961816546911624260e-4, places=12)
        self.assertAlmostEqual(rc2t[2][2], 0.9999998325501746670, places=12)

    def test_cal2jd(self):
        dmj0, dmj = erfa.cal2jd(2003, 6, 1)
        self.assertAlmostEqual(dmj0, 2400000.5, 9)
        self.assertAlmostEqual(dmj, 52791.0, 9)

    def test_d2dtf(self):
        y,m,d,H,M,S,F = erfa.d2dtf(5, 2400000.5, 49533.99999)
        self.assertEqual(y, 1994)
        self.assertEqual(m, 6)
        self.assertEqual(d, 30)
        self.assertEqual(H, 23)
        self.assertEqual(M, 59)
        self.assertEqual(S, 60)
        self.assertEqual(F, 13599)

    def test_d2tf(self):
        c,H,M,S,F = erfa.d2tf(4, -0.987654321)
        self.assertEqual(c, '-')
        self.assertEqual(H, 23)
        self.assertEqual(M, 42)
        self.assertEqual(S, 13)
        self.assertEqual(F, 3333)

    def test_dat(self):
        d = erfa.dat(2003, 6, 1, 0.0)
        self.assertAlmostEqual(d, 32.0, 9)
        d = erfa.dat(2008, 1, 17, 0.0)
        self.assertAlmostEqual(d, 33.0, 9)
        d = erfa.dat(2015, 9, 1, 0.0)
        self.assertAlmostEqual(d, 36.0, 9)
        d = erfa.dat(2017, 9, 1, 0.0)
        self.assertAlmostEqual(d, 37.0, 9)
        
    def test_dtdb(self):
        d = erfa.dtdb(2448939.5, 0.123, 0.76543, 5.0123, 5525.242, 3190.0)
        self.assertAlmostEqual(d, -0.1280368005936998991e-2, 17)

    def test_dtf2d(self):
        jd1, jd2 = erfa.dtf2d(1994, 6, 30, 23, 59, 60.13599)
        self.assertAlmostEqual(jd1+jd2, 2449534.49999, 6)

    def test_ee00(self):
        epsa =  0.4090789763356509900
        dpsi = -0.9630909107115582393e-5
        ee = erfa.ee00(2400000.5, 53736.0, epsa, dpsi)
        self.assertAlmostEqual(ee, -0.8834193235367965479e-5, 18)

    def test_ee00a(self):
        ee = erfa.ee00a(2400000.5, 53736.0)
        self.assertAlmostEqual(ee, -0.8834192459222588227e-5, 17) ##fail at 18

    def test_ee00b(self):
        ee = erfa.ee00b(2400000.5, 53736.0)
        self.assertAlmostEqual(ee, -0.8835700060003032831e-5, 18)

    def test_ee06a(self):
        ee = erfa.ee06a(2400000.5, 53736.0)
        self.assertAlmostEqual(ee, -0.8834195072043790156e-5, 15) ##fail at 16

    def test_eect00(self):
        ct = erfa.eect00(2400000.5, 53736.0)
        self.assertAlmostEqual(ct, 0.2046085004885125264e-8, 20)

    def test_eform(self):
        #a, f = erfa.eform(0)
        self.assertRaises(erfa.error, erfa.eform, 0)
        a, f = erfa.eform(1)
        self.assertAlmostEqual(a, 6378137.0, 10)
        self.assertAlmostEqual(f,  0.3352810664747480720e-2, 18)
        a, f = erfa.eform(2)
        self.assertAlmostEqual(a, 6378137.0, 10)
        self.assertAlmostEqual(f, 0.3352810681182318935e-2, 18)
        a, f = erfa.eform(3)
        self.assertAlmostEqual(a, 6378135.0, 10)
        self.assertAlmostEqual(f, 0.3352779454167504862e-2, 18)
        #a, f = erfa.eform(4)
        self.assertRaises(erfa.error, erfa.eform, 4)

    def test_eo06a(self):
        eo = erfa.eo06a(2400000.5, 53736.0)
        self.assertAlmostEqual(eo, -0.1332882371941833644e-2, 15)
        
    def test_eors(self):
        r = ((0.9999989440476103608,-0.1332881761240011518e-2,-0.5790767434730085097e-3),
             (0.1332858254308954453e-2,0.9999991109044505944,-0.4097782710401555759e-4),
             (0.5791308472168153320e-3,0.4020595661593994396e-4,0.9999998314954572365))
        s = -0.1220040848472271978e-7
        eo = erfa.eors(r,s)
        self.assertAlmostEqual(eo, -0.1332882715130744606e-2, 15)

    def test_epb(self):
        b = erfa.epb(2415019.8135, 30103.18648)
        self.assertAlmostEqual(b, 1982.418424159278580, 12)

    def test_epb2jd(self):
        dj0,dj1 = erfa.epb2jd(1957.3)
        self.assertAlmostEqual(dj0, 2400000.5, 9)
        self.assertAlmostEqual(dj1, 35948.1915101513, 9)

    def test_epj(self):
        j = erfa.epj(2451545, -7392.5)
        self.assertAlmostEqual(j, 1979.760438056125941, 12)

    def test_epj2jd(self):
        dj0,dj1 = erfa.epj2jd(1996.8)
        self.assertAlmostEqual(dj0, 2400000.5, 9)
        self.assertAlmostEqual(dj1, 50375.7, 9)

    def test_epv00(self):
        pvh,pvb = erfa.epv00(2400000.5, 53411.52501161)
        self.assertAlmostEqual(pvh[0][0], -0.7757238809297706813, places=14)
        self.assertAlmostEqual(pvh[0][1], 0.5598052241363340596, places=13) #fail at 14
        self.assertAlmostEqual(pvh[0][2], 0.2426998466481686993, places=14)

        self.assertAlmostEqual(pvh[1][0], -0.1091891824147313846e-1, places=15)
        self.assertAlmostEqual(pvh[1][1], -0.1247187268440845008e-1, places=15)
        self.assertAlmostEqual(pvh[1][2], -0.5407569418065039061e-2, places=15)
        
        self.assertAlmostEqual(pvb[0][0], -0.7714104440491111971, places=14)
        self.assertAlmostEqual(pvb[0][1], 0.5598412061824171323, places=13)  #fail at 14
        self.assertAlmostEqual(pvb[0][2], 0.2425996277722452400, places=14)

        self.assertAlmostEqual(pvb[1][0], -0.1091874268116823295e-1, places=15)
        self.assertAlmostEqual(pvb[1][1], -0.1246525461732861538e-1, places=15)
        self.assertAlmostEqual(pvb[1][2], -0.5404773180966231279e-2, places=15)

    def test_eqeq94(self):
        ee = erfa.eqeq94(2400000.5, 41234.0)
        self.assertAlmostEqual(ee, 0.5357758254609256894e-4, 17)

    def test_era00(self):
        era = erfa.era00(2400000.5, 54388.0)
        self.assertAlmostEqual(era, 0.4022837240028158102, 12)

    def test_fad03(self):
        d = erfa.fad03(0.80)
        self.assertAlmostEqual(d, 1.946709205396925672, 12)

    def test_fae03(self):
        e = erfa.fae03(0.80)
        self.assertAlmostEqual(e, 1.744713738913081846, 12)
        
    def test_faf03(self):
        f = erfa.faf03(0.80)
        self.assertAlmostEqual(f, 0.2597711366745499518, 11) # failed at 12
        
    def test_faju03(self):
        l = erfa.faju03(0.80)
        self.assertAlmostEqual(l, 5.275711665202481138, 12)
        
    def test_fal03(self):
        l = erfa.fal03(0.80)
        self.assertAlmostEqual(l, 5.132369751108684150, 12)
        
    def test_falp03(self):
        lp = erfa.falp03(0.80)
        self.assertAlmostEqual(lp, 6.226797973505507345, 12)
        
    def test_fama03(self):
        l = erfa.fama03(0.80)
        self.assertAlmostEqual(l, 3.275506840277781492, 12)
        
    def test_fame03(self):
        l = erfa.fame03(0.80)
        self.assertAlmostEqual(l, 5.417338184297289661, 12)
        
    def test_fane03(self):
        l = erfa.fane03(0.80)
        self.assertAlmostEqual(l, 2.079343830860413523, 12)
        
    def test_faom03(self):
        l = erfa.faom03(0.80)
        self.assertAlmostEqual(l, -5.973618440951302183, 12)
        
    def test_fapa03(self):
        l = erfa.fapa03(0.80)
        self.assertAlmostEqual(l, 0.1950884762240000000e-1, 12)
        
    def test_fasa03(self):
        l = erfa.fasa03(0.80)
        self.assertAlmostEqual(l, 5.371574539440827046, 12)
        
    def test_faur03(self):
        l = erfa.faur03(0.80)
        self.assertAlmostEqual(l, 5.180636450180413523, 12)
        
    def test_fave03(self):
        l = erfa.fave03(0.80)
        self.assertAlmostEqual(l, 3.424900460533758000, 12)

    def test_fk52h(self):
        r5  =  1.76779433
        d5  = -0.2917517103
        dr5 = -1.91851572e-7
        dd5 = -5.8468475e-6
        px5 =  0.379210
        rv5 = -7.6
        rh, dh, drh, ddh, pxh, rvh = erfa.fk52h(r5, d5, dr5, dd5, px5, rv5)
        self.assertAlmostEqual(rh, 1.767794226299947632, places=14)
        self.assertAlmostEqual(dh,  -0.2917516070530391757, places=14)
        self.assertAlmostEqual(drh, -0.19618741256057224e-6, places=19)
        self.assertAlmostEqual(ddh, -0.58459905176693911e-5, places=19)
        self.assertAlmostEqual(pxh,  0.37921, places=14)
        self.assertAlmostEqual(rvh, -7.6000000940000254, places=10)
        
    def test_fk5hip(self):
        r5h, s5h = erfa.fk5hip()
        self.assertAlmostEqual(r5h[0][0], 0.9999999999999928638, places=14)
        self.assertAlmostEqual(r5h[0][1], 0.1110223351022919694e-6, places=17)
        self.assertAlmostEqual(r5h[0][2], 0.4411803962536558154e-7, places=16) # failed at 17
        self.assertAlmostEqual(r5h[1][0], -0.1110223308458746430e-6, places=17)
        self.assertAlmostEqual(r5h[1][1], 0.9999999999999891830, places=14)
        self.assertAlmostEqual(r5h[1][2], -0.9647792498984142358e-7, places=17)
        self.assertAlmostEqual(r5h[2][0], -0.4411805033656962252e-7, places=16) # failed at 17
        self.assertAlmostEqual(r5h[2][1], 0.9647792009175314354e-7, places=17)
        self.assertAlmostEqual(r5h[2][2], 0.9999999999999943728, places=14)
        self.assertAlmostEqual(s5h[0], -0.1454441043328607981e-8, places=17)
        self.assertAlmostEqual(s5h[1], 0.2908882086657215962e-8, places=17)
        self.assertAlmostEqual(s5h[2], 0.3393695767766751955e-8, places=17)

    def test_fk5hz(self):
        r5 =  1.76779433
        d5 = -0.2917517103
        rh,dh = erfa.fk5hz(r5,d5,2400000.5, 54479.0)
        self.assertAlmostEqual(rh, 1.767794191464423978, 12)
        self.assertAlmostEqual(dh, -0.2917516001679884419, 12)

    def test_fw2m(self):
        gamb = -0.2243387670997992368e-5
        phib =  0.4091014602391312982
        psi  = -0.9501954178013015092e-3
        eps  =  0.4091014316587367472
        r = erfa.fw2m(gamb, phib, psi, eps)
        self.assertAlmostEqual(r[0][0], 0.9999995505176007047, places=12)
        self.assertAlmostEqual(r[0][1], 0.8695404617348192957e-3, places=12)
        self.assertAlmostEqual(r[0][2], 0.3779735201865582571e-3, places=12)

        self.assertAlmostEqual(r[1][0], -0.8695404723772016038e-3, places=12)
        self.assertAlmostEqual(r[1][1], 0.9999996219496027161, places=12)
        self.assertAlmostEqual(r[1][2], -0.1361752496887100026e-6, places=12)

        self.assertAlmostEqual(r[2][0], -0.3779734957034082790e-3, places=12)
        self.assertAlmostEqual(r[2][1], -0.1924880848087615651e-6, places=12)
        self.assertAlmostEqual(r[2][2], 0.9999999285679971958, places=12)
        
    def test_fw2xy(self):
        gamb = -0.2243387670997992368e-5
        phib =  0.4091014602391312982
        psi  = -0.9501954178013015092e-3
        eps  =  0.4091014316587367472
        x, y = erfa.fw2xy(gamb, phib, psi, eps)
        self.assertAlmostEqual(x,-0.3779734957034082790e-3, 14)
        self.assertAlmostEqual(y,-0.1924880848087615651e-6, 14)

    def test_gc2gd(self):
        xyz = (2e6, 3e6, 5.244e6)
        #e,p,h = erfa.gc2gd(0, xyz)
        self.assertRaises(erfa.error, erfa.gc2gd, 0, xyz)
        e,p,h = erfa.gc2gd(1, xyz)
        self.assertAlmostEqual(e, 0.9827937232473290680, places=14)
        self.assertAlmostEqual(p, 0.97160184819075459, places=14)
        self.assertAlmostEqual(h, 331.4172461426059892, places=8)
        e,p,h = erfa.gc2gd(2, xyz)
        self.assertAlmostEqual(e, 0.9827937232473290680, places=14)
        self.assertAlmostEqual(p, 0.97160184820607853, places=14)
        self.assertAlmostEqual(h, 331.41731754844348, places=8)
        e,p,h = erfa.gc2gd(3, xyz)
        self.assertAlmostEqual(e, 0.9827937232473290680, places=14)
        self.assertAlmostEqual(p, 0.9716018181101511937, places=14)
        self.assertAlmostEqual(h, 333.2770726130318123, places=8)
        #e,p,h = erfa.gc2gd(4, xyz)
        self.assertRaises(erfa.error, erfa.gc2gd, 4, xyz)

    def test_gc2gde(self):
        a = 6378136.0
        f = 0.0033528
        xyz = (2e6, 3e6, 5.244e6)
        e, p, h = erfa.gc2gde(a, f, xyz)
        self.assertAlmostEqual(e, 0.9827937232473290680, places=14)
        self.assertAlmostEqual(p, 0.9716018377570411532, places=14)
        self.assertAlmostEqual(h, 332.36862495764397, places=8)

    def test_gd2gc(self):
        e = 3.1
        p = -0.5
        h = 2500.0
        self.assertRaises(erfa.error, erfa.gd2gc, 0, e, p, h)
        xyz = erfa.gd2gc(1,e,p,h)
        self.assertAlmostEqual(xyz[0], -5599000.5577049947, places=7)
        self.assertAlmostEqual(xyz[1], 233011.67223479203, places=7)
        self.assertAlmostEqual(xyz[2], -3040909.4706983363, places=7)
        xyz = erfa.gd2gc(2,e,p,h)
        self.assertAlmostEqual(xyz[0], -5599000.5577260984, places=7)
        self.assertAlmostEqual(xyz[1], 233011.6722356703, places=7)
        self.assertAlmostEqual(xyz[2], -3040909.4706095476, places=7)
        xyz = erfa.gd2gc(3,e,p,h)
        self.assertAlmostEqual(xyz[0], -5598998.7626301490, places=7)
        self.assertAlmostEqual(xyz[1], 233011.5975297822, places=7)
        self.assertAlmostEqual(xyz[2], -3040908.6861467111, places=7)        
        self.assertRaises(erfa.error, erfa.gd2gc, 4, e, p, h)

    def test_gd2gce(self):
        a = 6378136.0
        f = 0.0033528
        e = 3.1
        p = -0.5
        h = 2500.0
        xyz = erfa.gd2gce(a, f, e, p, h)
        self.assertAlmostEqual(xyz[0], -5598999.6665116328, places=7)
        self.assertAlmostEqual(xyz[1], 233011.63514630572, places=7)
        self.assertAlmostEqual(xyz[2], -3040909.0517314132, places=7)

    def test_gmst00(self):
        g = erfa.gmst00(2400000.5, 53736.0, 2400000.5, 53736.0)
        self.assertAlmostEqual(g, 1.754174972210740592, 14)
        
    def test_gmst06(self):
        g = erfa.gmst06(2400000.5, 53736.0, 2400000.5, 53736.0)
        self.assertAlmostEqual(g, 1.754174971870091203, 14)
        
    def test_gmst82(self):
        g = erfa.gmst82(2400000.5, 53736.0)
        self.assertAlmostEqual(g, 1.754174981860675096, 14)
        
    def test_gst00a(self):
        g = erfa.gst00a(2400000.5, 53736.0, 2400000.5, 53736.0)
        self.assertAlmostEqual(g, 1.754166138018281369, 14)
        
    def test_gst00b(self):
        g = erfa.gst00b(2400000.5, 53736.0)
        self.assertAlmostEqual(g, 1.754166136510680589, 14)

    def test_gst06(self):
        rnpb = ((0.9999989440476103608,-0.1332881761240011518e-2,-0.5790767434730085097e-3),
                (0.1332858254308954453e-2,0.9999991109044505944,-0.4097782710401555759e-4),
                (0.5791308472168153320e-3,0.4020595661593994396e-4,0.9999998314954572365))
        g = erfa.gst06(2400000.5, 53736.0, 2400000.5, 53736.0, rnpb)
        self.assertAlmostEqual(g, 1.754166138018167568, 14)

    def test_gst06a(self):
        g = erfa.gst06a(2400000.5, 53736.0, 2400000.5, 53736.0)
        self.assertAlmostEqual(g, 1.754166137675019159, 14)
        
    def test_gst94(self):
        g = erfa.gst94(2400000.5, 53736.0)
        self.assertAlmostEqual(g, 1.754166136020645203, 14)

    def test_h2fk5(self):
        rh  =  1.767794352
        dh  = -0.2917512594
        drh = -2.76413026e-6
        ddh = -5.92994449e-6
        pxh =  0.379210
        rvh = -7.6;
        r5, d5, dr5, dd5, px5, rv5 = erfa.h2fk5(rh, dh, drh, ddh, pxh, rvh)
        self.assertAlmostEqual(r5, 1.767794455700065506, places=13)
        self.assertAlmostEqual(d5, -0.2917513626469638890, places=13)
        self.assertAlmostEqual(dr5, -0.27597945024511204e-5, places=18)
        self.assertAlmostEqual(dd5, -0.59308014093262838e-5, places=18)
        self.assertAlmostEqual(px5, 0.37921, places=13)
        self.assertAlmostEqual(rv5, -7.6000001309071126, places=10)

    def test_hfk5z(self):
        rh =  1.767794352
        dh = -0.2917512594
        r5, d5, dr5, dd5 = erfa.hfk5z(rh, dh, 2400000.5, 54479.0)
        self.assertAlmostEqual(r5, 1.767794490535581026, places=13)
        self.assertAlmostEqual(d5, -0.2917513695320114258, places=14)
        self.assertAlmostEqual(dr5, 0.4335890983539243029e-8, places=22)
        self.assertAlmostEqual(dd5, -0.8569648841237745902e-9, places=23)

    def test_jd2cal(self):
        y, m, d, fd = erfa.jd2cal(2400000.5, 50123.9999)
        self.assertEqual(y, 1996)
        self.assertEqual(m, 2)
        self.assertEqual(d, 10)
        self.assertAlmostEqual(fd, 0.9999, 7)
        
    def test_jdcalf(self):
        y, m, d, fd = erfa.jdcalf(4, 2400000.5, 50123.9999)
        self.assertEqual(y, 1996)
        self.assertEqual(m, 2)
        self.assertEqual(d, 10)
        self.assertEqual(fd, 9999)

    def test_num00a(self):
        rmatn = erfa.num00a(2400000.5, 53736.0)
        self.assertAlmostEqual(rmatn[0][0], 0.9999999999536227949, places=12)
        self.assertAlmostEqual(rmatn[0][1], 0.8836238544090873336e-5, places=12)
        self.assertAlmostEqual(rmatn[0][2], 0.3830835237722400669e-5, places=12)

        self.assertAlmostEqual(rmatn[1][0], -0.8836082880798569274e-5, places=12)
        self.assertAlmostEqual(rmatn[1][1], 0.9999999991354655028, places=12)
        self.assertAlmostEqual(rmatn[1][2], -0.4063240865362499850e-4, places=12)

        self.assertAlmostEqual(rmatn[2][0], -0.3831194272065995866e-5, places=12)
        self.assertAlmostEqual(rmatn[2][1], 0.4063237480216291775e-4, places=12)
        self.assertAlmostEqual(rmatn[2][2], 0.9999999991671660338, places=12)

    def test_num00b(self):
        rmatn = erfa.num00b(2400000.5, 53736.0)
        self.assertAlmostEqual(rmatn[0][0], 0.9999999999536069682, places=12)
        self.assertAlmostEqual(rmatn[0][1], 0.8837746144871248011e-5, places=12)
        self.assertAlmostEqual(rmatn[0][2], 0.3831488838252202945e-5, places=12)

        self.assertAlmostEqual(rmatn[1][0], -0.8837590456632304720e-5, places=12)
        self.assertAlmostEqual(rmatn[1][1], 0.9999999991354692733, places=12)
        self.assertAlmostEqual(rmatn[1][2], -0.4063198798559591654e-4, places=12)

        self.assertAlmostEqual(rmatn[2][0], -0.3831847930134941271e-5, places=12)
        self.assertAlmostEqual(rmatn[2][1], 0.4063195412258168380e-4, places=12)
        self.assertAlmostEqual(rmatn[2][2], 0.9999999991671806225, places=12)

    def test_num06a(self):
        rmatn = erfa.num06a(2400000.5, 53736)
        self.assertAlmostEqual(rmatn[0][0], 0.9999999999536227668, places=12)
        self.assertAlmostEqual(rmatn[0][1], 0.8836241998111535233e-5, places=12)
        self.assertAlmostEqual(rmatn[0][2], 0.3830834608415287707e-5, places=12)

        self.assertAlmostEqual(rmatn[1][0], -0.8836086334870740138e-5, places=12)
        self.assertAlmostEqual(rmatn[1][1], 0.9999999991354657474, places=12)
        self.assertAlmostEqual(rmatn[1][2], -0.4063240188248455065e-4, places=12)

        self.assertAlmostEqual(rmatn[2][0], -0.3831193642839398128e-5, places=12)
        self.assertAlmostEqual(rmatn[2][1], 0.4063236803101479770e-4, places=12)
        self.assertAlmostEqual(rmatn[2][2], 0.9999999991671663114, places=12)

    def test_numat(self):
        epsa =  0.4090789763356509900
        dpsi = -0.9630909107115582393e-5
        deps =  0.4063239174001678826e-4
        rmatn = erfa.numat(epsa, dpsi, deps)
        self.assertAlmostEqual(rmatn[0][0], 0.9999999999536227949, places=12)
        self.assertAlmostEqual(rmatn[0][1], 0.8836239320236250577e-5, places=12)
        self.assertAlmostEqual(rmatn[0][2], 0.3830833447458251908e-5, places=12)

        self.assertAlmostEqual(rmatn[1][0], -0.8836083657016688588e-5, places=12)
        self.assertAlmostEqual(rmatn[1][1], 0.9999999991354654959, places=12)
        self.assertAlmostEqual(rmatn[1][2], -0.4063240865361857698e-4, places=12)

        self.assertAlmostEqual(rmatn[2][0], -0.3831192481833385226e-5, places=12)
        self.assertAlmostEqual(rmatn[2][1], 0.4063237480216934159e-4, places=12)
        self.assertAlmostEqual(rmatn[2][2], 0.9999999991671660407, places=12)

    def test_nut00a(self):
        dpsi, deps = erfa.nut00a(2400000.5, 53736.0)
        self.assertAlmostEqual(dpsi, -0.9630909107115518431e-5, 13)
        self.assertAlmostEqual(deps, 0.4063239174001678710e-4, 13)
        
    def test_nut00b(self):
        dpsi, deps = erfa.nut00b(2400000.5, 53736.0)
        self.assertAlmostEqual(dpsi, -0.9632552291148362783e-5, 13)
        self.assertAlmostEqual(deps, 0.4063197106621159367e-4, 13)
        
    def test_nut06a(self):
        dpsi, deps = erfa.nut06a(2400000.5, 53736.0)
        self.assertAlmostEqual(dpsi, -0.9630912025820308797e-5, 13)
        self.assertAlmostEqual(deps, 0.4063238496887249798e-4, 13)
        
    def test_nut80(self):
        dpsi, deps = erfa.nut80(2400000.5, 53736.0)
        self.assertAlmostEqual(dpsi, -0.9643658353226563966e-5, 13)
        self.assertAlmostEqual(deps, 0.4060051006879713322e-4, 13)

    def test_nutm80(self):
        rmatn = erfa.nutm80(2400000.5, 53736)
        self.assertAlmostEqual(rmatn[0][0], 0.9999999999534999268, places=12)
        self.assertAlmostEqual(rmatn[0][1], 0.8847935789636432161e-5, places=12)
        self.assertAlmostEqual(rmatn[0][2], 0.3835906502164019142e-5, places=12)

        self.assertAlmostEqual(rmatn[1][0], -0.8847780042583435924e-5, places=12)
        self.assertAlmostEqual(rmatn[1][1], 0.9999999991366569963, places=12)
        self.assertAlmostEqual(rmatn[1][2], -0.4060052702727130809e-4, places=12)

        self.assertAlmostEqual(rmatn[2][0], -0.3836265729708478796e-5, places=12)
        self.assertAlmostEqual(rmatn[2][1], 0.4060049308612638555e-4, places=12)
        self.assertAlmostEqual(rmatn[2][2], 0.9999999991684415129, places=12)

    def test_obl06(self):
        obl = erfa.obl06(2400000.5, 54388.0)
        self.assertAlmostEqual(obl, 0.4090749229387258204, 14)
        
    def test_obl80(self):
        obl = erfa.obl80(2400000.5, 54388.0)
        self.assertAlmostEqual(obl, 0.4090751347643816218, 14)

    def test_p06e(self):
        eps0,psia,oma,bpa,bqa,pia,bpia,epsa,chia,za,zetaa,thetaa,pa,gam,phi,psi = erfa.p06e(2400000.5, 52541.0)
        self.assertAlmostEqual(eps0, 0.4090926006005828715, places=14)
        self.assertAlmostEqual(psia, 0.6664369630191613431e-3, places=14)
        self.assertAlmostEqual(oma , 0.4090925973783255982, places=14)
        self.assertAlmostEqual(bpa, 0.5561149371265209445e-6, places=14)
        self.assertAlmostEqual(bqa, -0.6191517193290621270e-5, places=14)
        self.assertAlmostEqual(pia, 0.6216441751884382923e-5, places=14)
        self.assertAlmostEqual(bpia, 3.052014180023779882, places=14)
        self.assertAlmostEqual(epsa, 0.4090864054922431688, places=14)
        self.assertAlmostEqual(chia, 0.1387703379530915364e-5, places=14)
        self.assertAlmostEqual(za, 0.2921789846651790546e-3, places=14)
        self.assertAlmostEqual(zetaa, 0.3178773290332009310e-3, places=14)
        self.assertAlmostEqual(thetaa, 0.2650932701657497181e-3, places=14)
        self.assertAlmostEqual(pa, 0.6651637681381016344e-3, places=14)
        self.assertAlmostEqual(gam, 0.1398077115963754987e-5, places=14)
        self.assertAlmostEqual(phi, 0.4090864090837462602, places=14)
        self.assertAlmostEqual(psi, 0.6664464807480920325e-3, places=14)

    def test_pb06(self):
        bzeta, bz, btheta = erfa.pb06(2400000.5, 50123.9999)
        self.assertAlmostEqual(bzeta, -0.5092634016326478238e-3, places=12)
        self.assertAlmostEqual(bz, -0.3602772060566044413e-3, places=12)
        self.assertAlmostEqual(btheta, -0.3779735537167811177e-3, places=12)

    def test_pfw06(self):
        gamb, phib, psib, epsa = erfa.pfw06(2400000.5, 50123.9999)
        self.assertAlmostEqual(gamb, -0.2243387670997995690e-5, places=16)
        self.assertAlmostEqual(phib,  0.4091014602391312808, places=12)
        self.assertAlmostEqual(psib, -0.9501954178013031895e-3, places=14)
        self.assertAlmostEqual(epsa,  0.4091014316587367491, places=12)

    def test_plan94(self):
        self.assertRaises(erfa.error, erfa.plan94, 2400000.5, -320000, 0)
        self.assertRaises(erfa.error, erfa.plan94, 2400000.5, -320000, 10)
        ##self.assertRaises(erfa.error, erfa.plan94, 2400000.5, -320000, 3)
        pv = erfa.plan94(2400000.5, -320000, 3)
        self.assertAlmostEqual(pv[0][0], 0.9308038666832975759, places=11)
        self.assertAlmostEqual(pv[0][1], 0.3258319040261346000, places=11)
        self.assertAlmostEqual(pv[0][2], 0.1422794544481140560, places=11)

        self.assertAlmostEqual(pv[1][0], -0.6429458958255170006e-2, places=11)
        self.assertAlmostEqual(pv[1][1], 0.1468570657704237764e-1, places=11)
        self.assertAlmostEqual(pv[1][2], 0.6406996426270981189e-2, places=11)

        pv = erfa.plan94(2400000.5, 43999.9, 1)
        self.assertAlmostEqual(pv[0][0], 0.2945293959257430832, places=11)
        self.assertAlmostEqual(pv[0][1], -0.2452204176601049596, places=11)
        self.assertAlmostEqual(pv[0][2], -0.1615427700571978153, places=11)

        self.assertAlmostEqual(pv[1][0], 0.1413867871404614441e-1, places=11)
        self.assertAlmostEqual(pv[1][1], 0.1946548301104706582e-1, places=11)
        self.assertAlmostEqual(pv[1][2], 0.8929809783898904786e-2, places=11)

    def test_pmat00(self):
        rbp = erfa.pmat00(2400000.5, 50123.9999)
        self.assertAlmostEqual(rbp[0][0], 0.9999995505175087260, places=12)
        self.assertAlmostEqual(rbp[0][1], 0.8695405883617884705e-3, places=14)
        self.assertAlmostEqual(rbp[0][2], 0.3779734722239007105e-3, places=14)

        self.assertAlmostEqual(rbp[1][0], -0.8695405990410863719e-3, places=14)
        self.assertAlmostEqual(rbp[1][1], 0.9999996219494925900, places=12)
        self.assertAlmostEqual(rbp[1][2], -0.1360775820404982209e-6, places=14)

        self.assertAlmostEqual(rbp[2][0], -0.3779734476558184991e-3, places=14)
        self.assertAlmostEqual(rbp[2][1], -0.1925857585832024058e-6, places=14)
        self.assertAlmostEqual(rbp[2][2], 0.9999999285680153377, places=12)

    def test_pmat06(self):
        rbp = erfa.pmat06(2400000.5, 50123.9999)
        self.assertAlmostEqual(rbp[0][0], 0.9999995505176007047, places=12)
        self.assertAlmostEqual(rbp[0][1], 0.8695404617348208406e-3, places=14)
        self.assertAlmostEqual(rbp[0][2], 0.3779735201865589104e-3, places=14)
        self.assertAlmostEqual(rbp[1][0], -0.8695404723772031414e-3, places=14)
        self.assertAlmostEqual(rbp[1][1], 0.9999996219496027161, places=12)
        self.assertAlmostEqual(rbp[1][2], -0.1361752497080270143e-6, places=14)
        self.assertAlmostEqual(rbp[2][0], -0.3779734957034089490e-3, places=14)
        self.assertAlmostEqual(rbp[2][1], -0.1924880847894457113e-6, places=14)
        self.assertAlmostEqual(rbp[2][2], 0.9999999285679971958, places=12)

    def test_pmat76(self):
        rmatp = erfa.pmat76(2400000.5, 50123.9999)
        self.assertAlmostEqual(rmatp[0][0], 0.9999995504328350733, places=12)
        self.assertAlmostEqual(rmatp[0][1], 0.8696632209480960785e-3, places=14)
        self.assertAlmostEqual(rmatp[0][2], 0.3779153474959888345e-3, places=14)
        self.assertAlmostEqual(rmatp[1][0], -0.8696632209485112192e-3, places=14)
        self.assertAlmostEqual(rmatp[1][1], 0.9999996218428560614, places=12)
        self.assertAlmostEqual(rmatp[1][2], -0.1643284776111886407e-6, places=14)
        self.assertAlmostEqual(rmatp[2][0], -0.3779153474950335077e-3, places=14)
        self.assertAlmostEqual(rmatp[2][1], -0.1643306746147366896e-6, places=14)
        self.assertAlmostEqual(rmatp[2][2], 0.9999999285899790119, places=12)

    def test_pn00(self):
        dpsi = -0.9632552291149335877e-5
        deps =  0.4063197106621141414e-4
        epsa, rb, rp, rbp, rn, rbpn = erfa.pn00(2400000.5, 53736.0, dpsi, deps)
        self.assertAlmostEqual(epsa, 0.4090791789404229916, 12)
        self.assertAlmostEqual(rb[0][0], 0.9999999999999942498, places=12)
        self.assertAlmostEqual(rb[0][1], -0.7078279744199196626e-7, places=18)
        self.assertAlmostEqual(rb[0][2], 0.8056217146976134152e-7, places=18)
        self.assertAlmostEqual(rb[1][0], 0.7078279477857337206e-7, places=18)
        self.assertAlmostEqual(rb[1][1], 0.9999999999999969484, places=12)
        self.assertAlmostEqual(rb[1][2], 0.3306041454222136517e-7, places=18)
        self.assertAlmostEqual(rb[2][0], -0.8056217380986972157e-7, places=18)
        self.assertAlmostEqual(rb[2][1], -0.3306040883980552500e-7, places=18)
        self.assertAlmostEqual(rb[2][2], 0.9999999999999962084, places=12)
        self.assertAlmostEqual(rp[0][0], 0.9999989300532289018, places=12)
        self.assertAlmostEqual(rp[0][1], -0.1341647226791824349e-2, places=14)
        self.assertAlmostEqual(rp[0][2], -0.5829880927190296547e-3, places=14)
        self.assertAlmostEqual(rp[1][0], 0.1341647231069759008e-2, places=14)
        self.assertAlmostEqual(rp[1][1], 0.9999990999908750433, places=12)
        self.assertAlmostEqual(rp[1][2], -0.3837444441583715468e-6, places=14)
        self.assertAlmostEqual(rp[2][0], 0.5829880828740957684e-3, places=14)
        self.assertAlmostEqual(rp[2][1], -0.3984203267708834759e-6, places=14)
        self.assertAlmostEqual(rp[2][2], 0.9999998300623538046, places=12)
        self.assertAlmostEqual(rbp[0][0], 0.9999989300052243993, places=12)
        self.assertAlmostEqual(rbp[0][1], -0.1341717990239703727e-2, places=14)
        self.assertAlmostEqual(rbp[0][2], -0.5829075749891684053e-3, places=14)
        self.assertAlmostEqual(rbp[1][0], 0.1341718013831739992e-2, places=14)
        self.assertAlmostEqual(rbp[1][1], 0.9999990998959191343, places=12)
        self.assertAlmostEqual(rbp[1][2], -0.3505759733565421170e-6, places=14)
        self.assertAlmostEqual(rbp[2][0], 0.5829075206857717883e-3, places=14)
        self.assertAlmostEqual(rbp[2][1], -0.4315219955198608970e-6, places=14)
        self.assertAlmostEqual(rbp[2][2], 0.9999998301093036269, places=12)
        self.assertAlmostEqual(rn[0][0], 0.9999999999536069682, places=12)
        self.assertAlmostEqual(rn[0][1], 0.8837746144872140812e-5, places=16)
        self.assertAlmostEqual(rn[0][2], 0.3831488838252590008e-5, places=16)
        self.assertAlmostEqual(rn[1][0], -0.8837590456633197506e-5, places=16)
        self.assertAlmostEqual(rn[1][1], 0.9999999991354692733, places=12)
        self.assertAlmostEqual(rn[1][2], -0.4063198798559573702e-4, places=15) #failed at 16
        self.assertAlmostEqual(rn[2][0], -0.3831847930135328368e-5, places=16)
        self.assertAlmostEqual(rn[2][1], 0.4063195412258150427e-4, places=16)
        self.assertAlmostEqual(rn[2][2], 0.9999999991671806225, places=12)
        self.assertAlmostEqual(rbpn[0][0], 0.9999989440499982806, places=12)
        self.assertAlmostEqual(rbpn[0][1], -0.1332880253640848301e-2, places=14)
        self.assertAlmostEqual(rbpn[0][2], -0.5790760898731087295e-3, places=14)
        self.assertAlmostEqual(rbpn[1][0], 0.1332856746979948745e-2, places=14)
        self.assertAlmostEqual(rbpn[1][1], 0.9999991109064768883, places=12)
        self.assertAlmostEqual(rbpn[1][2], -0.4097740555723063806e-4, places=14)
        self.assertAlmostEqual(rbpn[2][0], 0.5791301929950205000e-3, places=14)
        self.assertAlmostEqual(rbpn[2][1], 0.4020553681373702931e-4, places=14)
        self.assertAlmostEqual(rbpn[2][2], 0.9999998314958529887, places=12)

    def test_pn00a(self):
        dpsi, deps, epsa, rb, rp, rbp, rn, rbpn = erfa.pn00a(2400000.5, 53736.0)
        self.assertAlmostEqual(dpsi, -0.9630909107115518431e-5, places=12)
        self.assertAlmostEqual(deps,  0.4063239174001678710e-4, places=12)
        self.assertAlmostEqual(epsa,  0.4090791789404229916, places=12)
        self.assertAlmostEqual(rb[0][0], 0.9999999999999942498, places=12)
        self.assertAlmostEqual(rb[0][1], -0.7078279744199196626e-7, places=16)
        self.assertAlmostEqual(rb[0][2], 0.8056217146976134152e-7, places=16)
        self.assertAlmostEqual(rb[1][0], 0.7078279477857337206e-7, places=16)
        self.assertAlmostEqual(rb[1][1], 0.9999999999999969484, places=12)
        self.assertAlmostEqual(rb[1][2], 0.3306041454222136517e-7, places=16)
        self.assertAlmostEqual(rb[2][0], -0.8056217380986972157e-7, places=16)
        self.assertAlmostEqual(rb[2][1], -0.3306040883980552500e-7, places=16)
        self.assertAlmostEqual(rb[2][2], 0.9999999999999962084, places=12)
        self.assertAlmostEqual(rp[0][0], 0.9999989300532289018, places=12)
        self.assertAlmostEqual(rp[0][1], -0.1341647226791824349e-2, places=14)
        self.assertAlmostEqual(rp[0][2], -0.5829880927190296547e-3, places=14)
        self.assertAlmostEqual(rp[1][0], 0.1341647231069759008e-2, places=14)
        self.assertAlmostEqual(rp[1][1], 0.9999990999908750433, places=12)
        self.assertAlmostEqual(rp[1][2], -0.3837444441583715468e-6, places=14)
        self.assertAlmostEqual(rp[2][0], 0.5829880828740957684e-3, places=14)
        self.assertAlmostEqual(rp[2][1], -0.3984203267708834759e-6, places=14)
        self.assertAlmostEqual(rp[2][2], 0.9999998300623538046, places=12)
        self.assertAlmostEqual(rbp[0][0], 0.9999989300052243993, places=12)
        self.assertAlmostEqual(rbp[0][1], -0.1341717990239703727e-2, places=14)
        self.assertAlmostEqual(rbp[0][2], -0.5829075749891684053e-3, places=14)
        self.assertAlmostEqual(rbp[1][0], 0.1341718013831739992e-2, places=14)
        self.assertAlmostEqual(rbp[1][1], 0.9999990998959191343, places=12)
        self.assertAlmostEqual(rbp[1][2], -0.3505759733565421170e-6, places=14)
        self.assertAlmostEqual(rbp[2][0], 0.5829075206857717883e-3, places=14)
        self.assertAlmostEqual(rbp[2][1], -0.4315219955198608970e-6, places=14)
        self.assertAlmostEqual(rbp[2][2], 0.9999998301093036269, places=12)
        self.assertAlmostEqual(rn[0][0], 0.9999999999536227949, places=12)
        self.assertAlmostEqual(rn[0][1], 0.8836238544090873336e-5, places=14)
        self.assertAlmostEqual(rn[0][2], 0.3830835237722400669e-5, places=14)
        self.assertAlmostEqual(rn[1][0], -0.8836082880798569274e-5, places=14)
        self.assertAlmostEqual(rn[1][1], 0.9999999991354655028, places=12)
        self.assertAlmostEqual(rn[1][2], -0.4063240865362499850e-4, places=14)
        self.assertAlmostEqual(rn[2][0], -0.3831194272065995866e-5, places=14)
        self.assertAlmostEqual(rn[2][1], 0.4063237480216291775e-4, places=14)
        self.assertAlmostEqual(rn[2][2], 0.9999999991671660338, places=12)
        self.assertAlmostEqual(rbpn[0][0], 0.9999989440476103435, places=12)
        self.assertAlmostEqual(rbpn[0][1], -0.1332881761240011763e-2, places=14)
        self.assertAlmostEqual(rbpn[0][2], -0.5790767434730085751e-3, places=14)
        self.assertAlmostEqual(rbpn[1][0], 0.1332858254308954658e-2, places=14)
        self.assertAlmostEqual(rbpn[1][1], 0.9999991109044505577, places=12)
        self.assertAlmostEqual(rbpn[1][2], -0.4097782710396580452e-4, places=14)
        self.assertAlmostEqual(rbpn[2][0], 0.5791308472168152904e-3, places=14)
        self.assertAlmostEqual(rbpn[2][1], 0.4020595661591500259e-4, places=14)
        self.assertAlmostEqual(rbpn[2][2], 0.9999998314954572304, places=12)

    def test_pn00b(self):
        dpsi, deps, epsa, rb, rp, rbp, rn, rbpn = erfa.pn00b(2400000.5, 53736.0)
        self.assertAlmostEqual(dpsi, -0.9632552291148362783e-5, places=12)
        self.assertAlmostEqual(deps,  0.4063197106621159367e-4, places=12)
        self.assertAlmostEqual(epsa,  0.4090791789404229916, places=12)
        self.assertAlmostEqual(rb[0][0], 0.9999999999999942498, places=12)
        self.assertAlmostEqual(rb[0][1], -0.7078279744199196626e-7, places=16)
        self.assertAlmostEqual(rb[0][2], 0.8056217146976134152e-7, places=16)
        self.assertAlmostEqual(rb[1][0], 0.7078279477857337206e-7, places=16)
        self.assertAlmostEqual(rb[1][1], 0.9999999999999969484, places=12)
        self.assertAlmostEqual(rb[1][2], 0.3306041454222136517e-7, places=16)
        self.assertAlmostEqual(rb[2][0], -0.8056217380986972157e-7, places=16)
        self.assertAlmostEqual(rb[2][1], -0.3306040883980552500e-7, places=16)
        self.assertAlmostEqual(rb[2][2], 0.9999999999999962084, places=12)
        self.assertAlmostEqual(rp[0][0], 0.9999989300532289018, places=12)
        self.assertAlmostEqual(rp[0][1], -0.1341647226791824349e-2, places=14)
        self.assertAlmostEqual(rp[0][2], -0.5829880927190296547e-3, places=14)
        self.assertAlmostEqual(rp[1][0], 0.1341647231069759008e-2, places=14)
        self.assertAlmostEqual(rp[1][1], 0.9999990999908750433, places=12)
        self.assertAlmostEqual(rp[1][2], -0.3837444441583715468e-6, places=14)
        self.assertAlmostEqual(rp[2][0], 0.5829880828740957684e-3, places=14)
        self.assertAlmostEqual(rp[2][1], -0.3984203267708834759e-6, places=14)
        self.assertAlmostEqual(rp[2][2], 0.9999998300623538046, places=12)
        self.assertAlmostEqual(rbp[0][0], 0.9999989300052243993, places=12)
        self.assertAlmostEqual(rbp[0][1], -0.1341717990239703727e-2, places=14)
        self.assertAlmostEqual(rbp[0][2], -0.5829075749891684053e-3, places=14)
        self.assertAlmostEqual(rbp[1][0], 0.1341718013831739992e-2, places=14)
        self.assertAlmostEqual(rbp[1][1], 0.9999990998959191343, places=12)
        self.assertAlmostEqual(rbp[1][2], -0.3505759733565421170e-6, places=14)
        self.assertAlmostEqual(rbp[2][0], 0.5829075206857717883e-3, places=14)
        self.assertAlmostEqual(rbp[2][1], -0.4315219955198608970e-6, places=14)
        self.assertAlmostEqual(rbp[2][2], 0.9999998301093036269, places=12)
        self.assertAlmostEqual(rn[0][0], 0.9999999999536069682, places=12)
        self.assertAlmostEqual(rn[0][1], 0.8837746144871248011e-5, places=14)
        self.assertAlmostEqual(rn[0][2], 0.3831488838252202945e-5, places=14)
        self.assertAlmostEqual(rn[1][0], -0.8837590456632304720e-5, places=14)
        self.assertAlmostEqual(rn[1][1], 0.9999999991354692733, places=12)
        self.assertAlmostEqual(rn[1][2], -0.4063198798559591654e-4, places=14)
        self.assertAlmostEqual(rn[2][0], -0.3831847930134941271e-5, places=14)
        self.assertAlmostEqual(rn[2][1], 0.4063195412258168380e-4, places=14)
        self.assertAlmostEqual(rn[2][2], 0.9999999991671806225, places=12)
        self.assertAlmostEqual(rbpn[0][0], 0.9999989440499982806, places=12)
        self.assertAlmostEqual(rbpn[0][1], -0.1332880253640849194e-2, places=14)
        self.assertAlmostEqual(rbpn[0][2], -0.5790760898731091166e-3, places=14)
        self.assertAlmostEqual(rbpn[1][0], 0.1332856746979949638e-2, places=14)
        self.assertAlmostEqual(rbpn[1][1], 0.9999991109064768883, places=12)
        self.assertAlmostEqual(rbpn[1][2], -0.4097740555723081811e-4, places=14)
        self.assertAlmostEqual(rbpn[2][0], 0.5791301929950208873e-3, places=14)
        self.assertAlmostEqual(rbpn[2][1], 0.4020553681373720832e-4, places=14)
        self.assertAlmostEqual(rbpn[2][2], 0.9999998314958529887, places=12)

    def test_pn06(self):
        dpsi = -0.9632552291149335877e-5
        deps =  0.4063197106621141414e-4
        epsa, rb, rp, rbp, rn, rbpn = erfa.pn06(2400000.5, 53736.0, dpsi, deps)
        self.assertAlmostEqual(epsa, 0.4090789763356509926, places=12)
        self.assertAlmostEqual(rb[0][0], 0.9999999999999942497, places=12)
        self.assertAlmostEqual(rb[0][1], -0.7078368960971557145e-7, places=14)
        self.assertAlmostEqual(rb[0][2], 0.8056213977613185606e-7, places=14)
        self.assertAlmostEqual(rb[1][0], 0.7078368694637674333e-7, places=14)
        self.assertAlmostEqual(rb[1][1], 0.9999999999999969484, places=12)
        self.assertAlmostEqual(rb[1][2], 0.3305943742989134124e-7, places=14)
        self.assertAlmostEqual(rb[2][0], -0.8056214211620056792e-7, places=14)
        self.assertAlmostEqual(rb[2][1], -0.3305943172740586950e-7, places=14)
        self.assertAlmostEqual(rb[2][2], 0.9999999999999962084, places=12)
        self.assertAlmostEqual(rp[0][0], 0.9999989300536854831, places=12)
        self.assertAlmostEqual(rp[0][1], -0.1341646886204443795e-2, places=14)
        self.assertAlmostEqual(rp[0][2], -0.5829880933488627759e-3, places=14)
        self.assertAlmostEqual(rp[1][0], 0.1341646890569782183e-2, places=14)
        self.assertAlmostEqual(rp[1][1], 0.9999990999913319321, places=12)
        self.assertAlmostEqual(rp[1][2], -0.3835944216374477457e-6, places=14)
        self.assertAlmostEqual(rp[2][0], 0.5829880833027867368e-3, places=14)
        self.assertAlmostEqual(rp[2][1], -0.3985701514686976112e-6, places=14)
        self.assertAlmostEqual(rp[2][2], 0.9999998300623534950, places=12)
        self.assertAlmostEqual(rbp[0][0], 0.9999989300056797893, places=12)
        self.assertAlmostEqual(rbp[0][1], -0.1341717650545059598e-2, places=14)
        self.assertAlmostEqual(rbp[0][2], -0.5829075756493728856e-3, places=14)
        self.assertAlmostEqual(rbp[1][0], 0.1341717674223918101e-2, places=14)
        self.assertAlmostEqual(rbp[1][1], 0.9999990998963748448, places=12)
        self.assertAlmostEqual(rbp[1][2], -0.3504269280170069029e-6, places=14)
        self.assertAlmostEqual(rbp[2][0], 0.5829075211461454599e-3, places=14)
        self.assertAlmostEqual(rbp[2][1], -0.4316708436255949093e-6, places=14)
        self.assertAlmostEqual(rbp[2][2], 0.9999998301093032943, places=12)
        self.assertAlmostEqual(rn[0][0], 0.9999999999536069682, places=12)
        self.assertAlmostEqual(rn[0][1], 0.8837746921149881914e-5, places=14)
        self.assertAlmostEqual(rn[0][2], 0.3831487047682968703e-5, places=14)
        self.assertAlmostEqual(rn[1][0], -0.8837591232983692340e-5, places=14)
        self.assertAlmostEqual(rn[1][1], 0.9999999991354692664, places=12)
        self.assertAlmostEqual(rn[1][2], -0.4063198798558931215e-4, places=14)
        self.assertAlmostEqual(rn[2][0], -0.3831846139597250235e-5, places=14)
        self.assertAlmostEqual(rn[2][1], 0.4063195412258792914e-4, places=14)
        self.assertAlmostEqual(rn[2][2], 0.9999999991671806293, places=12)
        self.assertAlmostEqual(rbpn[0][0], 0.9999989440504506688, places=12)
        self.assertAlmostEqual(rbpn[0][1], -0.1332879913170492655e-2, places=14)
        self.assertAlmostEqual(rbpn[0][2], -0.5790760923225655753e-3, places=14)
        self.assertAlmostEqual(rbpn[1][0], 0.1332856406595754748e-2, places=14)
        self.assertAlmostEqual(rbpn[1][1], 0.9999991109069366795, places=12)
        self.assertAlmostEqual(rbpn[1][2], -0.4097725651142641812e-4, places=14)
        self.assertAlmostEqual(rbpn[2][0], 0.5791301952321296716e-3, places=14)
        self.assertAlmostEqual(rbpn[2][1], 0.4020538796195230577e-4, places=14)
        self.assertAlmostEqual(rbpn[2][2], 0.9999998314958576778, places=12)

    def test_pn06a(self):
        dpsi, deps, epsa, rb, rp, rbp, rn, rbpn = erfa.pn06a(2400000.5, 53736.0)
        self.assertAlmostEqual(dpsi, -0.9630912025820308797e-5, places=12)
        self.assertAlmostEqual(deps,  0.4063238496887249798e-4, places=12)
        self.assertAlmostEqual(epsa,  0.4090789763356509926, places=12)
        self.assertAlmostEqual(rb[0][0], 0.9999999999999942497, places=12)
        self.assertAlmostEqual(rb[0][1], -0.7078368960971557145e-7, places=14)
        self.assertAlmostEqual(rb[0][2], 0.8056213977613185606e-7, places=14)
        self.assertAlmostEqual(rb[1][0], 0.7078368694637674333e-7, places=14)
        self.assertAlmostEqual(rb[1][1], 0.9999999999999969484, places=12)
        self.assertAlmostEqual(rb[1][2], 0.3305943742989134124e-7, places=14)
        self.assertAlmostEqual(rb[2][0], -0.8056214211620056792e-7, places=14)
        self.assertAlmostEqual(rb[2][1], -0.3305943172740586950e-7, places=14)
        self.assertAlmostEqual(rb[2][2], 0.9999999999999962084, places=12)
        self.assertAlmostEqual(rp[0][0], 0.9999989300536854831, places=12)
        self.assertAlmostEqual(rp[0][1], -0.1341646886204443795e-2, places=14)
        self.assertAlmostEqual(rp[0][2], -0.5829880933488627759e-3, places=14)
        self.assertAlmostEqual(rp[1][0], 0.1341646890569782183e-2, places=14)
        self.assertAlmostEqual(rp[1][1], 0.9999990999913319321, places=12)
        self.assertAlmostEqual(rp[1][2], -0.3835944216374477457e-6, places=14)
        self.assertAlmostEqual(rp[2][0], 0.5829880833027867368e-3, places=14)
        self.assertAlmostEqual(rp[2][1], -0.3985701514686976112e-6, places=14)
        self.assertAlmostEqual(rp[2][2], 0.9999998300623534950, places=12)
        self.assertAlmostEqual(rbp[0][0], 0.9999989300056797893, places=12)
        self.assertAlmostEqual(rbp[0][1], -0.1341717650545059598e-2, places=14)
        self.assertAlmostEqual(rbp[0][2], -0.5829075756493728856e-3, places=14)
        self.assertAlmostEqual(rbp[1][0], 0.1341717674223918101e-2, places=14)
        self.assertAlmostEqual(rbp[1][1], 0.9999990998963748448, places=12)
        self.assertAlmostEqual(rbp[1][2], -0.3504269280170069029e-6, places=14)
        self.assertAlmostEqual(rbp[2][0], 0.5829075211461454599e-3, places=14)
        self.assertAlmostEqual(rbp[2][1], -0.4316708436255949093e-6, places=14)
        self.assertAlmostEqual(rbp[2][2], 0.9999998301093032943, places=12)
        self.assertAlmostEqual(rn[0][0], 0.9999999999536227668, places=12)
        self.assertAlmostEqual(rn[0][1], 0.8836241998111535233e-5, places=14)
        self.assertAlmostEqual(rn[0][2], 0.3830834608415287707e-5, places=14)
        self.assertAlmostEqual(rn[1][0], -0.8836086334870740138e-5, places=14)
        self.assertAlmostEqual(rn[1][1], 0.9999999991354657474, places=12)
        self.assertAlmostEqual(rn[1][2], -0.4063240188248455065e-4, places=14)
        self.assertAlmostEqual(rn[2][0], -0.3831193642839398128e-5, places=14)
        self.assertAlmostEqual(rn[2][1], 0.4063236803101479770e-4, places=14)
        self.assertAlmostEqual(rn[2][2], 0.9999999991671663114, places=12)
        self.assertAlmostEqual(rbpn[0][0], 0.9999989440480669738, places=12)
        self.assertAlmostEqual(rbpn[0][1], -0.1332881418091915973e-2, places=14)
        self.assertAlmostEqual(rbpn[0][2], -0.5790767447612042565e-3, places=14)
        self.assertAlmostEqual(rbpn[1][0], 0.1332857911250989133e-2, places=14)
        self.assertAlmostEqual(rbpn[1][1], 0.9999991109049141908, places=12)
        self.assertAlmostEqual(rbpn[1][2], -0.4097767128546784878e-4, places=14)
        self.assertAlmostEqual(rbpn[2][0], 0.5791308482835292617e-3, places=14)
        self.assertAlmostEqual(rbpn[2][1], 0.4020580099454020310e-4, places=14)
        self.assertAlmostEqual(rbpn[2][2], 0.9999998314954628695, places=12)

    def test_pnm00a(self):
        rbpn = erfa.pnm00a(2400000.5, 50123.9999)
        self.assertAlmostEqual(rbpn[0][0], 0.9999995832793134257, places=12)
        self.assertAlmostEqual(rbpn[0][1], 0.8372384254137809439e-3, places=14)
        self.assertAlmostEqual(rbpn[0][2], 0.3639684306407150645e-3, places=14)
        self.assertAlmostEqual(rbpn[1][0], -0.8372535226570394543e-3, places=14)
        self.assertAlmostEqual(rbpn[1][1], 0.9999996486491582471, places=12)
        self.assertAlmostEqual(rbpn[1][2], 0.4132915262664072381e-4, places=14)
        self.assertAlmostEqual(rbpn[2][0], -0.3639337004054317729e-3, places=14)
        self.assertAlmostEqual(rbpn[2][1], -0.4163386925461775873e-4, places=14)
        self.assertAlmostEqual(rbpn[2][2], 0.9999999329094390695, places=12)

    def test_pnm00b(self):
        rbpn = erfa.pnm00b(2400000.5, 50123.9999)
        self.assertAlmostEqual(rbpn[0][0], 0.9999995832776208280, places=12)
        self.assertAlmostEqual(rbpn[0][1], 0.8372401264429654837e-3, places=14)
        self.assertAlmostEqual(rbpn[0][2], 0.3639691681450271771e-3, places=14)
        self.assertAlmostEqual(rbpn[1][0], -0.8372552234147137424e-3, places=14)
        self.assertAlmostEqual(rbpn[1][1], 0.9999996486477686123, places=12)
        self.assertAlmostEqual(rbpn[1][2], 0.4132832190946052890e-4, places=14)
        self.assertAlmostEqual(rbpn[2][0], -0.3639344385341866407e-3, places=14)
        self.assertAlmostEqual(rbpn[2][1], -0.4163303977421522785e-4, places=14)
        self.assertAlmostEqual(rbpn[2][2], 0.9999999329092049734, places=12)

    def test_pnm06a(self):
        rbpn = erfa.pnm06a(2400000.5, 50123.9999)
        self.assertAlmostEqual(rbpn[0][0], 0.9999995832794205484, places=12)
        self.assertAlmostEqual(rbpn[0][1], 0.8372382772630962111e-3, places=14)
        self.assertAlmostEqual(rbpn[0][2], 0.3639684771140623099e-3, places=14)
        self.assertAlmostEqual(rbpn[1][0], -0.8372533744743683605e-3, places=14)
        self.assertAlmostEqual(rbpn[1][1], 0.9999996486492861646, places=12)
        self.assertAlmostEqual(rbpn[1][2], 0.4132905944611019498e-4, places=14)
        self.assertAlmostEqual(rbpn[2][0], -0.3639337469629464969e-3, places=14)
        self.assertAlmostEqual(rbpn[2][1], -0.4163377605910663999e-4, places=14)
        self.assertAlmostEqual(rbpn[2][2], 0.9999999329094260057, places=12)

    def test_pnm80(self):
        rmatpn = erfa.pnm80(2400000.5, 50123.9999)
        self.assertAlmostEqual(rmatpn[0][0], 0.9999995831934611169, places=12)
        self.assertAlmostEqual(rmatpn[0][1], 0.8373654045728124011e-3, places=14)
        self.assertAlmostEqual(rmatpn[0][2], 0.3639121916933106191e-3, places=14)
        self.assertAlmostEqual(rmatpn[1][0], -0.8373804896118301316e-3, places=14)
        self.assertAlmostEqual(rmatpn[1][1], 0.9999996485439674092, places=12)
        self.assertAlmostEqual(rmatpn[1][2], 0.4130202510421549752e-4, places=14)
        self.assertAlmostEqual(rmatpn[2][0], -0.3638774789072144473e-3, places=14)
        self.assertAlmostEqual(rmatpn[2][1], -0.4160674085851722359e-4, places=14)
        self.assertAlmostEqual(rmatpn[2][2], 0.9999999329310274805, places=12)

    def test_pom00(self):
        xp =  2.55060238e-7
        yp =  1.860359247e-6
        sp = -0.1367174580728891460e-10
        rpom = erfa.pom00(xp, yp, sp)
        self.assertAlmostEqual(rpom[0][0], 0.9999999999999674721, places=12)
        self.assertAlmostEqual(rpom[0][1], -0.1367174580728846989e-10, places=16)
        self.assertAlmostEqual(rpom[0][2], 0.2550602379999972345e-6, places=16)
        self.assertAlmostEqual(rpom[1][0], 0.1414624947957029801e-10, places=16)
        self.assertAlmostEqual(rpom[1][1], 0.9999999999982695317, places=12)
        self.assertAlmostEqual(rpom[1][2], -0.1860359246998866389e-5, places=16)
        self.assertAlmostEqual(rpom[2][0], -0.2550602379741215021e-6, places=16)
        self.assertAlmostEqual(rpom[2][1], 0.1860359247002414021e-5, places=16)
        self.assertAlmostEqual(rpom[2][2], 0.9999999999982370039, places=12)

    def test_pr00(self):
        dpsipr, depspr = erfa.pr00(2400000.5, 53736)
        self.assertAlmostEqual(dpsipr, -0.8716465172668347629e-7, places=22)
        self.assertAlmostEqual(depspr, -0.7342018386722813087e-8, places=22)

    def test_prec76(self):
        ep01 = 2400000.5
        ep02 = 33282.0
        ep11 = 2400000.5
        ep12 = 51544.0
        zeta, z, theta = erfa.prec76(ep01, ep02, ep11, ep12)
        self.assertAlmostEqual(zeta,  0.5588961642000161243e-2, places=12)
        self.assertAlmostEqual(z,     0.5589922365870680624e-2, places=12)
        self.assertAlmostEqual(theta, 0.4858945471687296760e-2, places=12)

    def test_pvstar(self):
        pv = ((126668.5912743160601,2136.792716839935195,-245251.2339876830091),
              (-0.4051854035740712739e-2,-0.6253919754866173866e-2,0.1189353719774107189e-1))
        ra, dec, pmr, pmd, px, rv = erfa.pvstar(pv)
        self.assertAlmostEqual(ra, 0.1686756e-1, places=12)
        self.assertAlmostEqual(dec, -1.093989828, places=12)
        self.assertAlmostEqual(pmr, -0.178323516e-4, places=16)
        self.assertAlmostEqual(pmd, 0.2336024047e-5, places=16)
        self.assertAlmostEqual(px, 0.74723, places=12)
        self.assertAlmostEqual(rv, -21.6, places=10) ## failed at 11

    def test_s00(self):
        x = 0.5791308486706011000e-3
        y = 0.4020579816732961219e-4
        s = erfa.s00(2400000.5, 53736.0, x, y)
        self.assertAlmostEqual(s, -0.1220036263270905693e-7, 18)
        
    def test_s00a(self):
        s = erfa.s00a(2400000.5, 52541.0)
        self.assertAlmostEqual(s, -0.1340684448919163584e-7, 18)
        
    def test_s00b(self):
        s = erfa.s00b(2400000.5, 52541.0)
        self.assertAlmostEqual(s, -0.1340695782951026584e-7, 18)
        
    def test_s06(self):
        x = 0.5791308486706011000e-3
        y = 0.4020579816732961219e-4
        s = erfa.s06(2400000.5, 53736.0, x, y)
        self.assertAlmostEqual(s, -0.1220032213076463117e-7, 18)
        
    def test_s06a(self):
        s = erfa.s06a(2400000.5, 52541.0)
        self.assertAlmostEqual(s, -0.1340680437291812383e-7, 18)

    def test_sp00(self):
        s = erfa.sp00(2400000.5, 52541.0)
        self.assertAlmostEqual(s, -0.6216698469981019309e-11, 12)

    def test_starpm(self):
        ra1 =   0.01686756
        dec1 = -1.093989828
        pmr1 = -1.78323516e-5
        pmd1 =  2.336024047e-6
        px1 =   0.74723
        rv1 = -21.6
        ra2, dec2, pmr2, pmd2, px2, rv2 = erfa.starpm(ra1, dec1, pmr1, pmd1, px1, rv1,
                 2400000.5, 50083.0, 2400000.5, 53736.0)
        self.assertAlmostEqual(ra2, 0.01668919069414242368, places=13)
        self.assertAlmostEqual(dec2, -1.093966454217127879, places=13)
        self.assertAlmostEqual(pmr2, -0.1783662682155932702e-4, places=17)
        self.assertAlmostEqual(pmd2, 0.2338092915987603664e-5, places=17)
        self.assertAlmostEqual(px2, 0.7473533835323493644, places=13)
        self.assertAlmostEqual(rv2, -21.59905170476860786, places=10)

    def test_starpv(self):
        ra =   0.01686756
        dec = -1.093989828
        pmr = -1.78323516e-5
        pmd =  2.336024047e-6
        px =   0.74723
        rv = -21.6
        pv = erfa.starpv(ra, dec, pmr, pmd, px, rv)
        self.assertAlmostEqual(pv[0][0], 126668.5912743160601, places=10)
        self.assertAlmostEqual(pv[0][1], 2136.792716839935195, places=12)
        self.assertAlmostEqual(pv[0][2], -245251.2339876830091, places=10)
        self.assertAlmostEqual(pv[1][0], -0.4051854035740712739e-2, places=13)
        self.assertAlmostEqual(pv[1][1], -0.6253919754866173866e-2, places=15)
        self.assertAlmostEqual(pv[1][2], 0.1189353719774107189e-1, places=13)

    def test_taitt(self):
        t1, t2 = erfa.taitt(2453750.5, 0.892482639)
        self.assertAlmostEqual(t1, 2453750.5, places=6)
        self.assertAlmostEqual(t2, 0.892855139, places=12)

    def test_taiut1(self):
        u1, u2 = erfa.taiut1(2453750.5, 0.892482639, -32.6659)
        self.assertAlmostEqual(u1, 2453750.5, places=6)
        self.assertAlmostEqual(u2, 0.8921045614537037037, places=12)

    def test_taiutc(self):
        u1, u2 = erfa.taiutc(2453750.5, 0.892482639)
        self.assertAlmostEqual(u1, 2453750.5, places=6)
        self.assertAlmostEqual(u2, 0.8921006945555555556, places=12)

    def test_tcbtdb(self):
        b1, b2 = erfa.tcbtdb(2453750.5, 0.893019599)
        self.assertAlmostEqual(b1, 2453750.5, places=6)
        self.assertAlmostEqual(b2, 0.8928551362746343397, places=12)
        
    def test_tcgtt(self):
        t1, t2 = erfa.tcgtt(2453750.5,  0.892862531)
        self.assertAlmostEqual(t1, 2453750.5, places=6)
        self.assertAlmostEqual(t2, 0.8928551387488816828, places=12)

    def test_tdbtcb(self):
        b1, b2 = erfa.tdbtcb(2453750.5, 0.892855137)
        self.assertAlmostEqual(b1, 2453750.5, places=6)
        self.assertAlmostEqual(b2, 0.8930195997253656716, places=12)

    def test_tdbtt(self):
        t1, t2 = erfa.tdbtt(2453750.5,  0.892855137, -0.000201)
        self.assertAlmostEqual(t1, 2453750.5, places=6)
        self.assertAlmostEqual(t2, 0.8928551393263888889, places=12)

    def test_tttai(self):
        t1, t2 = erfa.tttai(2453750.5, 0.892482639)
        self.assertAlmostEqual(t1, 2453750.5, places=6)
        self.assertAlmostEqual(t2, 0.892110139, places=12)

    def test_tttcg(self):
        t1, t2 = erfa.tttcg(2453750.5, 0.892482639)
        self.assertAlmostEqual(t1, 2453750.5, places=6)
        self.assertAlmostEqual(t2, 0.8924900312508587113, places=12)

    def test_tttdb(self):
        t1, t2 = erfa.tttdb(2453750.5, 0.892855139, -0.000201)
        self.assertAlmostEqual(t1, 2453750.5, places=6)
        self.assertAlmostEqual(t2, 0.8928551366736111111, places=12)

    def test_ttut1(self):
        t1, t2 = erfa.ttut1(2453750.5, 0.892855139, 64.8499)
        self.assertAlmostEqual(t1, 2453750.5, places=6)
        self.assertAlmostEqual(t2, 0.8921045614537037037, places=12)

    def test_ut1tai(self):
        a1, a2 = erfa.ut1tai(2453750.5, 0.892104561, -32.6659)
        self.assertAlmostEqual(a1, 2453750.5, places=6)
        self.assertAlmostEqual(a2, 0.8924826385462962963, places=12)
        
    def test_ut1tt(self):
        a1, a2 = erfa.ut1tt(2453750.5, 0.892104561, 64.8499)
        self.assertAlmostEqual(a1, 2453750.5, places=6)
        self.assertAlmostEqual(a2, 0.8928551385462962963, places=15)
        
    def test_ut1utc(self):
        a1, a2 = erfa.ut1utc(2453750.5, 0.892104561, 0.3341)
        self.assertAlmostEqual(a1, 2453750.5, places=6)
        self.assertAlmostEqual(a2, 0.8921006941018518519, places=13)
        
    def test_utctai(self):
        u1, u2 = erfa.utctai(2453750.5, 0.892100694)
        self.assertAlmostEqual(u1, 2453750.5, places=6)
        self.assertAlmostEqual(u2, 0.8924826384444444444, places=13)
        
    def test_utcut1(self):
        u1, u2 = erfa.utcut1(2453750.5, 0.892100694, 0.3341)
        self.assertAlmostEqual(u1, 2453750.5, places=6)
        self.assertAlmostEqual(u2, 0.8921045608981481481, places=13)

    def test_xy06(self):
        x, y = erfa.xy06(2400000.5, 53736.0)
        self.assertAlmostEqual(x, 0.5791308486706010975e-3, places=16)
        self.assertAlmostEqual(y, 0.4020579816732958141e-4, places=17)

    def test_xys00a(self):
        x, y, s = erfa.xys00a(2400000.5, 53736.0)
        self.assertAlmostEqual(x, 0.5791308472168152904e-3, places=16)
        self.assertAlmostEqual(y, 0.4020595661591500259e-4, places=17)
        self.assertAlmostEqual(s, -0.1220040848471549623e-7, places=20)
        
    def test_xys00b(self):
        x, y, s = erfa.xys00b(2400000.5, 53736.0)
        self.assertAlmostEqual(x, 0.5791301929950208873e-3, places=16)
        self.assertAlmostEqual(y, 0.4020553681373720832e-4, places=16)
        self.assertAlmostEqual(s, -0.1220027377285083189e-7, places=19)
        
    def test_xys06a(self):
        x, y, s = erfa.xys06a(2400000.5, 53736.0)
        self.assertAlmostEqual(x, 0.5791308482835292617e-3, places=16)
        self.assertAlmostEqual(y, 0.4020580099454020310e-4, places=15)
        self.assertAlmostEqual(s, -0.1220032294164579896e-7, places=19)

## vector matrix library

    def test_a2af(self):
        r = erfa.a2af(4, 2.345)
        self.assertEqual(r[0], '+')
        self.assertEqual(r[1], 134)
        self.assertEqual(r[2], 21)
        self.assertEqual(r[3], 30)
        self.assertEqual(r[4], 9706)

    def test_a2tf(self):
        r = erfa.a2tf(4, -3.01234)
        self.assertEqual(r[0], '-')
        self.assertEqual(r[1], 11)
        self.assertEqual(r[2], 30)
        self.assertEqual(r[3], 22)
        self.assertEqual(r[4], 6484)

    def test_af2a(self):
        r = erfa.af2a(-45, 13, 27.2)
        self.assertAlmostEqual(r, -0.7893115794313644842, 15)
        r = erfa.af2a(45, 13, 27.2)
        self.assertAlmostEqual(r, 0.7893115794313644842, 15)

    def test_anp(self):
        r = erfa.anp(-0.1)
        self.assertAlmostEqual(r, 6.183185307179586477, 15)
        
    def test_anpm(self):
        r = erfa.anpm(-4.0)
        self.assertAlmostEqual(r, 2.283185307179586477, 15)

    def test_c2s(self):
        t,p = erfa.c2s((100.,-50.,25.))
        self.assertAlmostEqual(t, -0.4636476090008061162, 15)
        self.assertAlmostEqual(p, 0.2199879773954594463, 15)

    def test_cp(self):
        p = (0.3,1.2,-2.5)
        c = erfa.cp(p)
        self.assertAlmostEqual(c[0], 0.3, places=12)        
        self.assertAlmostEqual(c[1], 1.2, places=12)        
        self.assertAlmostEqual(c[2], -2.5, places=12)        

    def test_cpv(self):
        pv = ((0.3,1.2,-2.5),
              (-0.5,3.1,0.9))
        c = erfa.cpv(pv)
        self.assertAlmostEqual(c[0][0],  0.3, 0.0,)
        self.assertAlmostEqual(c[0][1],  1.2, 0.0,)
        self.assertAlmostEqual(c[0][2], -2.5, 0.0,)
        self.assertAlmostEqual(c[1][0], -0.5, 0.0,)
        self.assertAlmostEqual(c[1][1],  3.1, 0.0,)
        self.assertAlmostEqual(c[1][2],  0.9, 0.0,)

    def test_cr(self):
        a = ((2.0,3.0,2.0),
             (3.0,2.0,3.0),
             (3.0,4.0,5.0))
        c = erfa.cr(a)
        self.assertAlmostEqual(c[0][0], 2.0, places=12)
        self.assertAlmostEqual(c[0][1], 3.0, places=12)
        self.assertAlmostEqual(c[0][2], 2.0, places=12)
        self.assertAlmostEqual(c[1][0], 3.0, places=12)
        self.assertAlmostEqual(c[1][1], 2.0, places=12)
        self.assertAlmostEqual(c[1][2], 3.0, places=12)
        self.assertAlmostEqual(c[2][0], 3.0, places=12)
        self.assertAlmostEqual(c[2][1], 4.0, places=12)
        self.assertAlmostEqual(c[2][2], 5.0, places=12)

    def test_p2pv(self):
        pv = erfa.p2pv((0.25,1.2,3.0))
        self.assertAlmostEqual(pv[0][0], 0.25, places=12)
        self.assertAlmostEqual(pv[0][1], 1.2, places=12)
        self.assertAlmostEqual(pv[0][2], 3.0, places=12)
        self.assertAlmostEqual(pv[1][0], 0.0, places=12)
        self.assertAlmostEqual(pv[1][1], 0.0, places=12)
        self.assertAlmostEqual(pv[1][2], 0.0, places=12)
        
    def test_p2s(self):
        theta, phi, r = erfa.p2s((100,-50,25))
        self.assertAlmostEqual(theta, -0.4636476090008061162, 12)
        self.assertAlmostEqual(phi, 0.2199879773954594463, 12)
        self.assertAlmostEqual(r, 114.5643923738960002, 12)

    def test_pap(self):
        a = (1.,0.1,0.2)
        b= (-3.,1e-3,0.2)
        theta = erfa.pap(a,b)
        self.assertAlmostEqual(theta, 0.3671514267841113674, 12)

    def test_pas(self):
        p = erfa.pas(1.0,0.1,0.2,-1.0)
        self.assertAlmostEqual(p, -2.724544922932270424, 12)

    def test_pdp(self):
        ab =erfa.pdp((2.,2.,3.),(1.,3.,4.))
        self.assertAlmostEqual(ab, 20, 12)

    def test_pm(self):
##        m = erfa.pm((2.,2.,3.))
##        r = math.sqrt(2.*2.+2.*2.+3.*3.)
##        self.assertAlmostEqual(m, r, 14)
        m = erfa.pm((0.3,1.2,-2.5))
        self.assertAlmostEqual(m, 2.789265136196270604, 14)
        
    def test_pmp(self):
        a = (2.0,2.0,3.0)
        b = (1.0,3.0,4.0)
        amb= erfa.pmp(a,b)
        self.assertAlmostEqual(amb[0], 1.0, 12)
        self.assertAlmostEqual(amb[1], -1.0, 12)
        self.assertAlmostEqual(amb[2], -1.0, 12)

    def test_pn(self):
        r, u = erfa.pn((0.3,1.2,-2.5))
        self.assertAlmostEqual(r, 2.789265136196270604, 12)
        self.assertAlmostEqual(u[0], 0.1075552109073112058, 12)
        self.assertAlmostEqual(u[1], 0.4302208436292448232, 12)
        self.assertAlmostEqual(u[2], -0.8962934242275933816, 12)

    def test_ppp(self):
        apb = erfa.ppp((2.0,2.0,3.0),(1.0,3.0,4.0))
        self.assertAlmostEqual(apb[0], 3.0, places=12)
        self.assertAlmostEqual(apb[1], 5.0, places=12)
        self.assertAlmostEqual(apb[2], 7.0, places=12)

    def test_ppsp(self):
        apsb = erfa.ppsp((2.0,2.0,3.0),5.0,(1.0,3.0,4.0))
        self.assertAlmostEqual(apsb[0], 7.0, places=12)
        self.assertAlmostEqual(apsb[1], 17.0, places=12)
        self.assertAlmostEqual(apsb[2], 23.0, places=12)

    def test_pv2p(self):
        p = erfa.pv2p(((0.3,1.2,-2.5),(-0.5,3.1,0.9)))
        self.assertAlmostEqual(p[0],  0.3, 0.0)
        self.assertAlmostEqual(p[1],  1.2, 0.0)
        self.assertAlmostEqual(p[2], -2.5, 0.0)

    def test_pv2s(self):
        pv = ((-0.4514964673880165,0.03093394277342585,0.05594668105108779),
              (1.292270850663260e-5,2.652814182060692e-6,2.568431853930293e-6))
        theta, phi, r, td, pd, rd = erfa.pv2s(pv)
        self.assertAlmostEqual(theta, 3.073185307179586515, places=12)
        self.assertAlmostEqual(phi, 0.1229999999999999992, places=12)
        self.assertAlmostEqual(r, 0.4559999999999999757, places=12)
        self.assertAlmostEqual(td, -0.7800000000000000364e-5, places=16)
        self.assertAlmostEqual(pd, 0.9010000000000001639e-5, places=16)
        self.assertAlmostEqual(rd, -0.1229999999999999832e-4, places=16)

    def test_pvdpv(self):
        a = ((2.,2.,3.),(6.,0.,4.))
        b = ((1.,3.,4.),(0.,2.,8.))
        adb = erfa.pvdpv(a, b)
        self.assertAlmostEqual(adb[0], 20.0, places=12)
        self.assertAlmostEqual(adb[1], 50.0, places=12)

    def test_pvm(self):
        pv = ((0.3,1.2,-2.5),(0.45,-0.25,1.1))
        r, s = erfa.pvm(pv)
        self.assertAlmostEqual(r, 2.789265136196270604, 12)
        self.assertAlmostEqual(s, 1.214495780149111922, 12)

    def test_pvmpv(self):
        a = ((2.0,2.0,3.0),(5.0,6.0,3.0))
        b = ((1.0,3.0,4.0),(3.0,2.0,1.0))
        amb = erfa.pvmpv(a,b)
        self.assertAlmostEqual(amb[0][0],  1.0, places=12)
        self.assertAlmostEqual(amb[0][1], -1.0, places=12)
        self.assertAlmostEqual(amb[0][2], -1.0, places=12)
        self.assertAlmostEqual(amb[1][0],  2.0, places=12)
        self.assertAlmostEqual(amb[1][1],  4.0, places=12)
        self.assertAlmostEqual(amb[1][2],  2.0, places=12)

    def test_pvppv(self):
        a = ((2.0,2.0,3.0),(5.0,6.0,3.0))
        b = ((1.0,3.0,4.0),(3.0,2.0,1.0))
        apb = erfa.pvppv(a,b)
        self.assertAlmostEqual(apb[0][0], 3.0, places=12)
        self.assertAlmostEqual(apb[0][1], 5.0, places=12)
        self.assertAlmostEqual(apb[0][2], 7.0, places=12)
        self.assertAlmostEqual(apb[1][0], 8.0, places=12)
        self.assertAlmostEqual(apb[1][1], 8.0, places=12)
        self.assertAlmostEqual(apb[1][2], 4.0, places=12)

    def test_pvu(self):
        dt = 2920.0
        pv = ((126668.5912743160734,2136.792716839935565,-245251.2339876830229),
              (-0.4051854035740713039e-2,-0.6253919754866175788e-2,0.1189353719774107615e-1))
        upv = erfa.pvu(dt, pv)
        self.assertAlmostEqual(upv[0][0], 126656.7598605317105, places=12)
        self.assertAlmostEqual(upv[0][1], 2118.531271155726332, places=12)
        self.assertAlmostEqual(upv[0][2], -245216.5048590656190, places=12)
        self.assertAlmostEqual(upv[1][0], -0.4051854035740713039e-2, places=12)
        self.assertAlmostEqual(upv[1][1], -0.6253919754866175788e-2, places=12)
        self.assertAlmostEqual(upv[1][2], 0.1189353719774107615e-1, places=12)

    def test_pvup(self):
        dt = 2920.0
        pv = ((126668.5912743160734,2136.792716839935565,-245251.2339876830229),
              (-0.4051854035740713039e-2,-0.6253919754866175788e-2,0.1189353719774107615e-1))
        p = erfa.pvup(dt, pv)
        self.assertAlmostEqual(p[0],  126656.7598605317105,   1e-12)
        self.assertAlmostEqual(p[1],    2118.531271155726332, places=12)
        self.assertAlmostEqual(p[2], -245216.5048590656190,   1e-12)

    def test_pvxpv(self):
        a = ((2.0,2.0,3.0),(6.0,0.0,4.0))
        b = ((1.0,3.0,4.0),(0.0,2.0,8.0))
        axb =erfa.pvxpv(a, b)
        self.assertAlmostEqual(axb[0][0],  -1.0, places=12)
        self.assertAlmostEqual(axb[0][1],  -5.0, places=12)
        self.assertAlmostEqual(axb[0][2],   4.0, places=12)
        self.assertAlmostEqual(axb[1][0],  -2.0, places=12,)
        self.assertAlmostEqual(axb[1][1], -36.0, places=12,)
        self.assertAlmostEqual(axb[1][2],  22.0, places=12,)

    def test_pxp(self):
        axb = erfa.pxp((2.0,2.0,3.0),(1.0,3.0,4.0))
        self.assertAlmostEqual(axb[0], -1.0, places=12)
        self.assertAlmostEqual(axb[1], -5.0, places=12)
        self.assertAlmostEqual(axb[2],  4.0, places=12)

    def test_rm2v(self):
        w = erfa.rm2v(((0.0,-0.8,-0.6),
                      (0.8,-0.36,0.48),
                      (0.6,0.48,-0.64)))
        self.assertAlmostEqual(w[0],  0.0, places=12)
        self.assertAlmostEqual(w[1],  1.413716694115406957, places=12)
        self.assertAlmostEqual(w[2], -1.884955592153875943, places=12)

    def test_rv2m(self):
        r = erfa.rv2m((0.0, 1.41371669, -1.88495559))
        self.assertAlmostEqual(r[0][0], -0.7071067782221119905, places=14)
        self.assertAlmostEqual(r[0][1], -0.5656854276809129651, places=14)
        self.assertAlmostEqual(r[0][2], -0.4242640700104211225, places=14)
        self.assertAlmostEqual(r[1][0],  0.5656854276809129651, places=14)
        self.assertAlmostEqual(r[1][1], -0.0925483394532274246, places=14)
        self.assertAlmostEqual(r[1][2], -0.8194112531408833269, places=14)
        self.assertAlmostEqual(r[2][0],  0.4242640700104211225, places=14)
        self.assertAlmostEqual(r[2][1], -0.8194112531408833269, places=14)
        self.assertAlmostEqual(r[2][2],  0.3854415612311154341, places=14)

    def test_rx(self):
        phi = 0.3456789
        r = ((2.0,3.0,2.0),
             (3.0,2.0,3.0),
             (3.0,4.0,5.0))
        r = erfa.rx(phi, r)
        self.assertAlmostEqual(r[0][0], 2.0, 0.0)
        self.assertAlmostEqual(r[0][1], 3.0, 0.0)
        self.assertAlmostEqual(r[0][2], 2.0, 0.0)
        self.assertAlmostEqual(r[1][0], 3.839043388235612460, places=12)
        self.assertAlmostEqual(r[1][1], 3.237033249594111899, places=12)
        self.assertAlmostEqual(r[1][2], 4.516714379005982719, places=12)
        self.assertAlmostEqual(r[2][0], 1.806030415924501684, places=12)
        self.assertAlmostEqual(r[2][1], 3.085711545336372503, places=12)
        self.assertAlmostEqual(r[2][2], 3.687721683977873065, places=12)

    def test_rxp(self):
        r = ((2.0,3.0,2.0),
             (3.0,2.0,3.0),
             (3.0,4.0,5.0))
        p = (0.2,1.5,0.1)
        rp = erfa.rxp(r, p)
        self.assertAlmostEqual(rp[0], 5.1, places=12)
        self.assertAlmostEqual(rp[1], 3.9, places=12)
        self.assertAlmostEqual(rp[2], 7.1, places=12)

    def test_rxpv(self):
        r = ((2.0,3.0,2.0),
             (3.0,2.0,3.0),
             (3.0,4.0,5.0))
        pv = ((0.2,1.5,0.1),(1.5,0.2,0.1))
        rpv = erfa.rxpv(r, pv)
        self.assertAlmostEqual(rpv[0][0], 5.1, places=12)
        self.assertAlmostEqual(rpv[1][0], 3.8, places=12)
        self.assertAlmostEqual(rpv[0][1], 3.9, places=12)
        self.assertAlmostEqual(rpv[1][1], 5.2, places=12)
        self.assertAlmostEqual(rpv[0][2], 7.1, places=12)
        self.assertAlmostEqual(rpv[1][2], 5.8, places=12)

    def test_rxr(self):
        a = ((2.0,3.0,2.0),
             (3.0,2.0,3.0),
             (3.0,4.0,5.0))
        b = ((1.0,2.0,2.0),
             (4.0,1.0,1.0),
             (3.0,0.0,1.0))
        atb = erfa.rxr(a, b)
        self.assertAlmostEqual(atb[0][0], 20.0, places=12)
        self.assertAlmostEqual(atb[0][1],  7.0, places=12)
        self.assertAlmostEqual(atb[0][2],  9.0, places=12)
        self.assertAlmostEqual(atb[1][0], 20.0, places=12)
        self.assertAlmostEqual(atb[1][1],  8.0, places=12)
        self.assertAlmostEqual(atb[1][2], 11.0, places=12)
        self.assertAlmostEqual(atb[2][0], 34.0, places=12)
        self.assertAlmostEqual(atb[2][1], 10.0, places=12)
        self.assertAlmostEqual(atb[2][2], 15.0, places=12)

    def test_ry(self):
        theta = 0.3456789
        r = ((2.0,3.0,2.0),
             (3.0,2.0,3.0),
             (3.0,4.0,5.0))
        r = erfa.ry(theta, r)
        self.assertAlmostEqual(r[0][0], 0.8651847818978159930, places=12)
        self.assertAlmostEqual(r[0][1], 1.467194920539316554, places=12)
        self.assertAlmostEqual(r[0][2], 0.1875137911274457342, places=12)
        self.assertAlmostEqual(r[1][0], 3, places=12)
        self.assertAlmostEqual(r[1][1], 2, places=12)
        self.assertAlmostEqual(r[1][2], 3, places=12)
        self.assertAlmostEqual(r[2][0], 3.500207892850427330, places=12)
        self.assertAlmostEqual(r[2][1], 4.779889022262298150, places=12)
        self.assertAlmostEqual(r[2][2], 5.381899160903798712, places=12)

    def test_rz(self):
        psi = 0.3456789
        r = ((2.0,3.0,2.0),
             (3.0,2.0,3.0),
             (3.0,4.0,5.0))
        r = erfa.rz(psi, r)
        self.assertAlmostEqual(r[0][0], 2.898197754208926769, places=12)
        self.assertAlmostEqual(r[0][1], 3.500207892850427330, places=12)
        self.assertAlmostEqual(r[0][2], 2.898197754208926769, places=12)
        self.assertAlmostEqual(r[1][0], 2.144865911309686813, places=12)
        self.assertAlmostEqual(r[1][1], 0.865184781897815993, places=12)
        self.assertAlmostEqual(r[1][2], 2.144865911309686813, places=12)
        self.assertAlmostEqual(r[2][0], 3.0, places=12)
        self.assertAlmostEqual(r[2][1], 4.0, places=12)
        self.assertAlmostEqual(r[2][2], 5.0, places=12)

    def test_s2c(self):
        c = erfa.s2c(3.0123, -0.999)
        self.assertAlmostEqual(c[0], -0.5366267667260523906, places=12)
        self.assertAlmostEqual(c[1],  0.0697711109765145365, places=12)
        self.assertAlmostEqual(c[2], -0.8409302618566214041, places=12)

    def test_s2p(self):
        p = erfa.s2p(-3.21, 0.123, 0.456)
        self.assertAlmostEqual(p[0], -0.4514964673880165228, places=12)
        self.assertAlmostEqual(p[1],  0.0309339427734258688, places=12)
        self.assertAlmostEqual(p[2],  0.0559466810510877933, places=12)

    def test_s2pv(self):
        pv = erfa.s2pv(-3.21, 0.123, 0.456, -7.8e-6, 9.01e-6, -1.23e-5)
        self.assertAlmostEqual(pv[0][0], -0.4514964673880165228, places=12)
        self.assertAlmostEqual(pv[0][1],  0.0309339427734258688, places=12)
        self.assertAlmostEqual(pv[0][2],  0.0559466810510877933, places=12)
        self.assertAlmostEqual(pv[1][0],  0.1292270850663260170e-4, places=16)
        self.assertAlmostEqual(pv[1][1],  0.2652814182060691422e-5, places=16)
        self.assertAlmostEqual(pv[1][2],  0.2568431853930292259e-5, places=16)

    def test_s2xpv(self):
        s1 = 2.0
        s2 = 3.0
        pv = ((0.3,1.2,-2.5),
              (0.5,2.3,-0.4))
        spv = erfa.s2xpv(s1, s2, pv)
        self.assertAlmostEqual(spv[0][0],  0.6, places=12)
        self.assertAlmostEqual(spv[0][1],  2.4, places=12)
        self.assertAlmostEqual(spv[0][2], -5.0, places=12)
        self.assertAlmostEqual(spv[1][0],  1.5, places=12)
        self.assertAlmostEqual(spv[1][1],  6.9, places=12)
        self.assertAlmostEqual(spv[1][2], -1.2, places=12)

    def test_sepp(self):
        a = (1.0,0.1,0.2)
        b = (-3.0,1e-3,0.2)
        s = erfa.sepp(a,b)
        self.assertAlmostEqual(s, 2.860391919024660768, 12)

    def test_seps(self):
        al =  1.0
        ap =  0.1
        bl =  0.2
        bp = -3.0
        s = erfa.seps(al, ap, bl, bp)
        self.assertAlmostEqual(s, 2.346722016996998842, 14)

    def test_sxp(self):
        s = 2.0
        p = (0.3,1.2,-2.5)
        sp = erfa.sxp(s, p)
        self.assertAlmostEqual(sp[0],  0.6, 0.0)
        self.assertAlmostEqual(sp[1],  2.4, 0.0)
        self.assertAlmostEqual(sp[2], -5.0, 0.0)

    def test_sxpv(self):
        s = 2.0
        pv = ((0.3,1.2,-2.5),(0.5,3.2,-0.7))
        spv = erfa.sxpv(s, pv)
        self.assertAlmostEqual(spv[0][0],  0.6, 0.0)
        self.assertAlmostEqual(spv[0][1],  2.4, 0.0)
        self.assertAlmostEqual(spv[0][2], -5.0, 0.0)
        self.assertAlmostEqual(spv[1][0],  1.0, 0.0)
        self.assertAlmostEqual(spv[1][1],  6.4, 0.0)
        self.assertAlmostEqual(spv[1][2], -1.4, 0.0)

    def test_tf2a(self):
        a = erfa.tf2a(4,58,20.2)
        self.assertAlmostEqual(a, 1.301739278189537429, 12)
        
    def test_tf2d(self):
        d = erfa.tf2d(23,55,10.9)
        self.assertAlmostEqual(d, 0.9966539351851851852, 12)

    def test_tr(self):
        r = ((2.0,3.0,2.0),
             (3.0,2.0,3.0),
             (3.0,4.0,5.0))
        rt = erfa.tr(r)
        self.assertAlmostEqual(rt[0][0], 2.0, 0.0)
        self.assertAlmostEqual(rt[0][1], 3.0, 0.0)
        self.assertAlmostEqual(rt[0][2], 3.0, 0.0)
        self.assertAlmostEqual(rt[1][0], 3.0, 0.0)
        self.assertAlmostEqual(rt[1][1], 2.0, 0.0)
        self.assertAlmostEqual(rt[1][2], 4.0, 0.0)
        self.assertAlmostEqual(rt[2][0], 2.0, 0.0)
        self.assertAlmostEqual(rt[2][1], 3.0, 0.0)
        self.assertAlmostEqual(rt[2][2], 5.0, 0.0)        
        
    def test_trxp(self):
        r = ((2.0,3.0,2.0),
             (3.0,2.0,3.0),
             (3.0,4.0,5.0))
        p = (0.2,1.5,0.1)
        trp = erfa.trxp(r, p)
        self.assertAlmostEqual(trp[0], 5.2, places=12)
        self.assertAlmostEqual(trp[1], 4.0, places=12)
        self.assertAlmostEqual(trp[2], 5.4, places=12)

    def test_trxpv(self):
        r = ((2.0,3.0,2.0),
             (3.0,2.0,3.0),
             (3.0,4.0,5.0))
        pv = ((0.2,1.5,0.1),(1.5,0.2,0.1))
        trpv = erfa.trxpv(r, pv)
        self.assertAlmostEqual(trpv[0][0], 5.2, places=12)
        self.assertAlmostEqual(trpv[0][1], 4.0, places=12)
        self.assertAlmostEqual(trpv[0][2], 5.4, places=12)
        self.assertAlmostEqual(trpv[1][0], 3.9, places=12)
        self.assertAlmostEqual(trpv[1][1], 5.3, places=12)
        self.assertAlmostEqual(trpv[1][2], 4.1, places=12)

    def test_icrs2g(self):
        dr =  5.9338074302227188048671087
        dd = -1.1784870613579944551540570
        dl, db = erfa.icrs2g(dr,dd)
        self.assertAlmostEqual(dl,  5.5850536063818546461558, places=14)
        self.assertAlmostEqual(db, -0.7853981633974483096157, places=14)

    def test_g2icrs(self):
        dl =  5.5850536063818546461558105
        db = -0.7853981633974483096156608
        dr, dd = erfa.g2icrs(dl, db)
        self.assertAlmostEqual(dr,  5.9338074302227188048671, places=14)
        self.assertAlmostEqual(dd, -1.1784870613579944551541, places=14)

    def test_ltp(self):
        epj = 1666.666
        rp = erfa.ltp(epj)
        self.assertAlmostEqual(rp[0][0], 0.9967044141159213819, places=14)
        self.assertAlmostEqual(rp[0][1], 0.7437801893193210840e-1, places=14)
        self.assertAlmostEqual(rp[0][2], 0.3237624409345603401e-1, places=14)
        self.assertAlmostEqual(rp[1][0], -0.7437802731819618167e-1, places=14)
        self.assertAlmostEqual(rp[1][1], 0.9972293894454533070, places=14)
        self.assertAlmostEqual(rp[1][2], -0.1205768842723593346e-2, places=14)
        self.assertAlmostEqual(rp[2][0], -0.3237622482766575399e-1, places=14)
        self.assertAlmostEqual(rp[2][1], -0.1206286039697609008e-2, places=14)
        self.assertAlmostEqual(rp[2][2], 0.9994750246704010914, places=14)

    def test_ltpb(self):
        epj = 1666.666
        rpb = erfa.ltpb(epj)
        self.assertAlmostEqual(rpb[0][0], 0.9967044167723271851, places=14)
        self.assertAlmostEqual(rpb[0][1], 0.7437794731203340345e-1, places=14)
        self.assertAlmostEqual(rpb[0][2], 0.3237632684841625547e-1, places=14)
        self.assertAlmostEqual(rpb[1][0], -0.7437795663437177152e-1, places=14)
        self.assertAlmostEqual(rpb[1][1], 0.9972293947500013666, places=14)
        self.assertAlmostEqual(rpb[1][2], -0.1205741865911243235e-2, places=14)
        self.assertAlmostEqual(rpb[2][0], -0.3237630543224664992e-1, places=14)
        self.assertAlmostEqual(rpb[2][1], -0.1206316791076485295e-2, places=14)
        self.assertAlmostEqual(rpb[2][2], 0.9994750220222438819, places=14)

    def test_ltpecl(self):
        epj = -1500.0
        vec = erfa.ltpecl(epj)
        self.assertAlmostEqual(vec[0], 0.4768625676477096525e-3, places=14)
        self.assertAlmostEqual(vec[1], -0.4052259533091875112, places=14)
        self.assertAlmostEqual(vec[2], 0.9142164401096448012, places=14)

    def test_ltpequ(self):
        epj = -2500.0
        veq = erfa.ltpequ(epj)
        self.assertAlmostEqual(veq[0], -0.3586652560237326659, places=14)
        self.assertAlmostEqual(veq[1], -0.1996978910771128475, places=14)
        self.assertAlmostEqual(veq[2], 0.9118552442250819624, places=14)

    def test_eceq06(self):
        date1 = 2456165.5
        date2 = 0.401182685
        dl = 5.1
        db = -0.9
        dr, dd = erfa.eceq06(date1, date2, dl, db)
        self.assertAlmostEqual(dr, 5.533459733613627767, places=14)
        self.assertAlmostEqual(dd, -1.246542932554480576, places=14)

    def test_ecm06(self):
        date1 = 2456165.5
        date2 = 0.401182685
        rm = erfa.ecm06(date1, date2)
        self.assertAlmostEqual(rm[0][0], 0.9999952427708701137, places=14)
        self.assertAlmostEqual(rm[0][1], -0.2829062057663042347e-2, places=14)
        self.assertAlmostEqual(rm[0][2], -0.1229163741100017629e-2, places=14)
        self.assertAlmostEqual(rm[1][0], 0.3084546876908653562e-2, places=14)
        self.assertAlmostEqual(rm[1][1], 0.9174891871550392514, places=14)
        self.assertAlmostEqual(rm[1][2], 0.3977487611849338124, places=14)
        self.assertAlmostEqual(rm[2][0], 0.2488512951527405928e-5, places=14)
        self.assertAlmostEqual(rm[2][1], -0.3977506604161195467, places=14)
        self.assertAlmostEqual(rm[2][2], 0.9174935488232863071, places=14)

    def test_eqec06(self):
        date1 = 1234.5
        date2 = 2440000.5
        dr = 1.234
        dd = 0.987
        dl, db = erfa.eqec06(date1, date2, dr, dd)
        self.assertAlmostEqual(dl, 1.342509918994654619, places=14)
        self.assertAlmostEqual(db, 0.5926215259704608132, places=14)

    def test_lteceq(self):
        epj = 2500.0
        dl = 1.5
        db = 0.6
        dr, dd = erfa.lteceq(epj, dl, db)
        self.assertAlmostEqual(dr, 1.275156021861921167, places=14)
        self.assertAlmostEqual(dd, 0.9966573543519204791, places=14)

    def test_ltecm(self):
        epj = -3000.0
        rm = erfa.ltecm(epj)
        self.assertAlmostEqual(rm[0][0], 0.3564105644859788825, places=14)
        self.assertAlmostEqual(rm[0][1], 0.8530575738617682284, places=14)
        self.assertAlmostEqual(rm[0][2], 0.3811355207795060435, places=14)
        self.assertAlmostEqual(rm[1][0], -0.9343283469640709942, places=14)
        self.assertAlmostEqual(rm[1][1], 0.3247830597681745976, places=14)
        self.assertAlmostEqual(rm[1][2], 0.1467872751535940865, places=14)
        self.assertAlmostEqual(rm[2][0], 0.1431636191201167793e-2, places=14)
        self.assertAlmostEqual(rm[2][1], -0.4084222566960599342, places=14)
        self.assertAlmostEqual(rm[2][2], 0.9127919865189030899, places=14)

    def test_lteqec(self):
        epj = -1500.0
        dr = 1.234
        dd = 0.987
        dl, db = erfa.lteqec(epj, dr, dd)
        self.assertAlmostEqual(dl, 0.5039483649047114859, places=14)
        self.assertAlmostEqual(db, 0.5848534459726224882, places=14)
        
## some older version of Python 2.7
## failed to use 'support'
if __name__ == '__main__':
    if sys.version_info.major >= 3:
        try:
            from test import support
        except ImportError:
            from test import test_support as support
        support.run_unittest(Validate)
    else:
        unittest.main()
