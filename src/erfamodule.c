#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <math.h>
#include <time.h>
#include "erfa.h"
#include "erfam.h"

static PyObject *erfaError;

static PyObject *
erfa_bi00(PyObject *self)
{
    double dpsibi, depsbi, dra;
    eraBi00(&dpsibi, &depsbi, &dra);
    if (!PyErr_Occurred()) {
        return Py_BuildValue("ddd", dpsibi, depsbi, dra);
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_bi00_doc,
"bi00() -> dpsibi,depsbi,dra\n\n"
"Frame bias components of IAU 2000 precession-nutation models (part of MHB2000 with additions).\n"
"Returned:\n"
"    dpsibi,depsbi    obliquity and correction\n"
"    dra              the ICRS RA of the J2000.0 mean equinox");

static PyObject *
erfa_bp00(PyObject *self, PyObject *args)
{
    double d1, d2, rb[3][3], rp[3][3], rbp[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1, &d2);
    if (ok) {
        eraBp00(d1, d2, rb, rp, rbp);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd)),((ddd)(ddd)(ddd)),((ddd)(ddd)(ddd))",
        rb[0][0],rb[0][1],rb[0][2],rb[1][0],rb[1][1],rb[1][2],rb[2][0],rb[2][1],rb[2][2],
        rp[0][0],rp[0][1],rp[0][2],rp[1][0],rp[1][1],rp[1][2],rp[2][0],rp[2][1],rp[2][2],
        rbp[0][0],rbp[0][1],rbp[0][2],rbp[1][0],rbp[1][1],rbp[1][2],rbp[2][0],rbp[2][1],rbp[2][2]
        );
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_bp00_doc,
"bp00(d1, d2) -> rb, rp, rbp\n\n"
"Frame bias and precession, IAU 2000.\n"
"Given:\n"
"    d1, d2     TT as a 2-part Julian Date\n"
"Returned:\n"
"    rb         frame bias matrix\n"
"    rp         precession matrix\n"
"    rbp        bias-precession matrix");

static PyObject *
erfa_bp06(PyObject *self, PyObject *args)
{
    double d1, d2, rb[3][3], rp[3][3], rbp[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1, &d2);
    if (ok) {
        eraBp06(d1, d2, rb, rp, rbp);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd)),((ddd)(ddd)(ddd)),((ddd)(ddd)(ddd))",
        rb[0][0],rb[0][1],rb[0][2],rb[1][0],rb[1][1],rb[1][2],rb[2][0],rb[2][1],rb[2][2],
        rp[0][0],rp[0][1],rp[0][2],rp[1][0],rp[1][1],rp[1][2],rp[2][0],rp[2][1],rp[2][2],
        rbp[0][0],rbp[0][1],rbp[0][2],rbp[1][0],rbp[1][1],rbp[1][2],rbp[2][0],rbp[2][1],rbp[2][2]
        );
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_bp06_doc,
"bp06(d1, d2) -> rb, rp, rbp\n\n"
"Frame bias and precession, IAU 2006.\n"
"Given:\n"
"    d1, d2     TT as 2-part Julian Date\n"
"Returned:\n"
"    rb         frame bias matrix\n"
"    p          precession matrix)\n"
"    rbp        bias-precession matrix");

static PyObject *
erfa_bpn2xy(PyObject *self, PyObject *args)
{
    double r00,r01,r02,r10,r11,r12,r20,r21,r22;
    double x, y, rbpn[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))",
    &r00,&r01,&r02,&r10,&r11,&r12,&r20,&r21,&r22);
    if (ok) {
        rbpn[0][0] = r00;
        rbpn[0][1] = r01;
        rbpn[0][2] = r02;
        rbpn[1][0] = r10;
        rbpn[1][1] = r11;
        rbpn[1][2] = r12;
        rbpn[2][0] = r20;
        rbpn[2][1] = r21;
        rbpn[2][2] = r22;
        eraBpn2xy(rbpn, &x, &y);
        return Py_BuildValue("dd", x, y);
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_bpn2xy_doc,
"bpn2xy(rbpn[3][3] -> x, y\n\n"
"Extract from the bias-precession-nutation matrix\n"
"the X,Y coordinates of the Celestial Intermediate Pole.\n"
"Given:\n"
"    rbpn       celestial-to-true matrix\n"
"Returned:\n"
"    x,y        celestial Intermediate Pole");

static PyObject *
erfa_c2i00a(PyObject *self, PyObject *args)
{
    double d1, d2, rc2i[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1, &d2);
    if (ok) {
        eraC2i00a(d1, d2, rc2i);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2i[0][0],rc2i[0][1],rc2i[0][2],
        rc2i[1][0],rc2i[1][1],rc2i[1][2],
        rc2i[2][0],rc2i[2][1],rc2i[2][2]
        );
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_c2i00a_doc,
"c2i00a(d1, d2) -> rc2i\n\n"
"Form the celestial-to-intermediate matrix\n"
"for a given date using the IAU 2000A precession-nutation model.\n"
"Given:\n"
"    d1, d2     TT as 2-part Julian date\n"
"Returned:\n"
"    rc2i       celestial-to-intermediate matrix\n");

static PyObject *
erfa_c2i00b(PyObject *self, PyObject *args)
{
    double d1, d2, rc2i[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1, &d2);
    if (ok) {
        eraC2i00b(d1, d2, rc2i);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2i[0][0],rc2i[0][1],rc2i[0][2],
        rc2i[1][0],rc2i[1][1],rc2i[1][2],
        rc2i[2][0],rc2i[2][1],rc2i[2][2]
        );
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_c2i00b_doc,
"c2i00b(d1, d2) -> rc2i\n\n"
"Form the celestial-to-intermediate matrix for a given\n"
"date using the IAU 2000B precession-nutation model.\n"
"Given:\n"
"    d1,d2      TT as 2-part Julian date\n"
"Returned:\n"
"    rc2i       celestial-to-intermediate matrix");

static PyObject *
erfa_c2i06a(PyObject *self, PyObject *args)
{
    double d1, d2, rc2i[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1, &d2);
    if (ok) {
        eraC2i06a(d1, d2, rc2i);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2i[0][0],rc2i[0][1],rc2i[0][2],
        rc2i[1][0],rc2i[1][1],rc2i[1][2],
        rc2i[2][0],rc2i[2][1],rc2i[2][2]
        );
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_c2i06a_doc,
"c2i06a(d1, d2) -> rc2i\n\n"
"Form the celestial-to-intermediate matrix for a given date \n"
"using the IAU 2006 precession and IAU 200A nutation model.\n"
"Given:\n"
"    d1,d2      TT as 2-part Julian date\n"
"Returned:\n"
"    rc2i       celestial-to-intermediate matrix");

static PyObject *
erfa_c2ibpn(PyObject *self, PyObject *args)
{
    double r00,r01,r02,r10,r11,r12,r20,r21,r22;
    double d1, d2, rbpn[3][3], rc2i[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dd((ddd)(ddd)(ddd))", &d1, &d2,
    &r00,&r01,&r02,&r10,&r11,&r12,&r20,&r21,&r22);
    if (ok) {
        rbpn[0][0] = r00;
        rbpn[0][1] = r01;
        rbpn[0][2] = r02;
        rbpn[1][0] = r10;
        rbpn[1][1] = r11;
        rbpn[1][2] = r12;
        rbpn[2][0] = r20;
        rbpn[2][1] = r21;
        rbpn[2][2] = r22;
        eraC2ibpn(d1, d2, rbpn, rc2i);
        return Py_BuildValue("((ddd)(ddd)(ddd))",
                                rc2i[0][0],rc2i[0][1],rc2i[0][2],
                                rc2i[1][0],rc2i[1][1],rc2i[1][2],
                                rc2i[2][0],rc2i[2][1],rc2i[2][2]
        );
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_c2ibpn_doc,
"c2ibpn(d1, d2, rbpn) -> rc2i\n\n"
"Form the celestial-to-intermediate matrix for a given date\n"
"and given the bias-precession-nutation matrix. IAU 2000.\n"
"Given:\n"
"    d1,d2      TT as 2-part Julian date\n"
"    rbpn       celestial-to-true matrix\n"
"Returned:\n"
"    rc2i       celestial-to-intermediate matrix");

static PyObject *
erfa_c2ixy(PyObject *self, PyObject *args)
{
    double d1, d2, x, y, rc2i[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dddd", &d1, &d2, &x, &y);
    if (ok) {
        eraC2ixy(d1,d2,x,y,rc2i);
        return Py_BuildValue("((ddd)(ddd)(ddd))",
                                rc2i[0][0],rc2i[0][1],rc2i[0][2],
                                rc2i[1][0],rc2i[1][1],rc2i[1][2],
                                rc2i[2][0],rc2i[2][1],rc2i[2][2]
        );    
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_c2ixy_doc,
"c2ixy(d1, d2, x, y) -> rc2i\n\n"
"Form the celestial to intermediate-frame-of-date matrix for a given\n"
"date when the CIP X,Y coordinates are known. IAU 2000.\n"
"Given:\n"
"    d1,d2      TT as 2-part Julian Date\n"
"    x, y       the Celestial Intermediate Pole\n"
"Returned:\n"
"    rc2i       the celestial-to-intermediate matrix");

static PyObject *
erfa_c2ixys(PyObject *self, PyObject *args)
{
    double x, y, s, rc2i[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "ddd", &x, &y, &s);
    if (ok) {
        eraC2ixys(x,y,s,rc2i);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2i[0][0],rc2i[0][1],rc2i[0][2],
        rc2i[1][0],rc2i[1][1],rc2i[1][2],
        rc2i[2][0],rc2i[2][1],rc2i[2][2]
        );    
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_c2ixys_doc,
"c2ixys(x, y, s) -> rc2i\n\n"
"Form the celestial to intermediate-frame-of-date matrix\n"
" given the CIP X,Y and the CIO locator s.\n"
"Given:\n"
"    x, y       Celestial Intermediate Pole\n"
"    s          the CIO locator \n"
"Returned:\n"
"   rc2i        celestial-to-intermediate matrix");

static PyObject *
erfa_c2t00a(PyObject *self, PyObject *args)
{
    double tta,ttb,uta,utb,xp,yp,rc2t[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dddddd", &tta,&ttb,&uta,&utb,&xp,&yp);
    if (ok) {
        eraC2t00a(tta,ttb,uta,utb,xp,yp,rc2t);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2t[0][0],rc2t[0][1],rc2t[0][2],
        rc2t[1][0],rc2t[1][1],rc2t[1][2],
        rc2t[2][0],rc2t[2][1],rc2t[2][2]
        ); 
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_c2t00a_doc,
"c2t00a(tta,ttb,uta,utb,xp,yp) -> rc2t\n\n"
"Form the celestial to terrestrial matrix given the date, the UT1 and\n"
"the polar motion, using the IAU 2000A nutation model.\n"
"Given:\n"
"    tta,ttb    TT as 2-part Julian Date\n"
"    uta,utb    UT1 as 2-part Julian Date\n"
"    xp, yp     the coordinates of the pole (radians)\n"
"Returned:\n"
"    rc2t       the celestial-to-terrestrial matrix");

static PyObject *
erfa_c2t00b(PyObject *self, PyObject *args)
{
    double tta,ttb,uta,utb,xp,yp,rc2t[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dddddd", &tta,&ttb,&uta,&utb,&xp,&yp);
    if (ok) {
        eraC2t00b(tta,ttb,uta,utb,xp,yp,rc2t);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2t[0][0],rc2t[0][1],rc2t[0][2],
        rc2t[1][0],rc2t[1][1],rc2t[1][2],
        rc2t[2][0],rc2t[2][1],rc2t[2][2]
        ); 
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_c2t00b_doc,
"c2t00b(tta,ttb,uta,utb,xp,yp) -> rc2t\n\n"
"Form the celestial to terrestrial matrix given the date, the UT1 and\n"
"the polar motion, using the IAU 2000B nutation model.\n"
"Given:\n"
"    tta,ttb    TT as 2-part Julian Date\n"
"    uta,utb    UT1 as 2-part Julian Date\n"
"    xp, yp     the coordinates of the pole (radians)\n"
"Returned:\n"
"    rc2t       the celestial-to-terrestrial matrix");

static PyObject *
erfa_c2t06a(PyObject *self, PyObject *args)
{
    double tta,ttb,uta,utb,xp,yp,rc2t[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dddddd", &tta,&ttb,&uta,&utb,&xp,&yp);
    if (ok) {
        eraC2t06a(tta,ttb,uta,utb,xp,yp,rc2t);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2t[0][0],rc2t[0][1],rc2t[0][2],
        rc2t[1][0],rc2t[1][1],rc2t[1][2],
        rc2t[2][0],rc2t[2][1],rc2t[2][2]
        ); 
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_c2t06a_doc,
"c2t06a(tta,ttb,uta,utb,xp,yp) -> rc2t\n\n"
"Form the celestial to terrestrial matrix given the date, the UT1 and\n"
"the polar motion, using the IAU 2006 precession and IAU 2000A nutation models.\n"
"Given:\n"
"    tta,ttb    TT as 2-part Julian Date\n"
"    uta,utb    UT1 as 2-part Julian Date\n"
"    xp, yp     the coordinates of the pole (radians)\n"
"Returned:\n"
"    rc2t       the celestial-to-terrestrial matrix"
);

static PyObject *
erfa_c2tcio(PyObject *self, PyObject *args)
{
    double c00,c01,c02,c10,c11,c12,c20,c21,c22,rc2i[3][3];
    double era, rc2t[3][3];
    double p00,p01,p02,p10,p11,p12,p20,p21,p22,rpom[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))d((ddd)(ddd)(ddd))",
    &c00,&c01,&c02,&c10,&c11,&c12,&c20,&c21,&c22,
    &era,
    &p00,&p01,&p02,&p10,&p11,&p12,&p20,&p21,&p22);
    if (ok) {
        rc2i[0][0] = c00;
        rc2i[0][1] = c01;
        rc2i[0][2] = c02;
        rc2i[1][0] = c10;
        rc2i[1][1] = c11;
        rc2i[1][2] = c12;
        rc2i[2][0] = c20;
        rc2i[2][1] = c21;
        rc2i[2][2] = c22;
        rpom[0][0] = p00;
        rpom[0][1] = p01;
        rpom[0][2] = p02;
        rpom[1][0] = p10;
        rpom[1][1] = p11;
        rpom[1][2] = p12;
        rpom[2][0] = p20;
        rpom[2][1] = p21;
        rpom[2][2] = p22;
        eraC2tcio(rc2i, era, rpom, rc2t);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2t[0][0],rc2t[0][1],rc2t[0][2],
        rc2t[1][0],rc2t[1][1],rc2t[1][2],
        rc2t[2][0],rc2t[2][1],rc2t[2][2]
        );
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_c2tcio_doc,
"c2tcio(rc2i, era, rpom) -> rc2t\n\n"
"Assemble the celestial to terrestrial matrix from CIO-based components\n"
"(the celestial-to-intermediate matrix, the Earth Rotation Angle and the polar motion matrix)\n"
"Given:\n"
"    rc2        celestial-to-intermediate matrix\n"
"    era        Earth rotation angle\n"
"    rpom       polar-motion matrix\n"
"Returned:t\n"
"    rc2t       the celestial-to-terrestrial matrix");

static PyObject *
erfa_c2teqx(PyObject *self, PyObject *args)
{
    double c00,c01,c02,c10,c11,c12,c20,c21,c22,rc2i[3][3];
    double gst, rc2t[3][3];
    double p00,p01,p02,p10,p11,p12,p20,p21,p22,rpom[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))d((ddd)(ddd)(ddd))",
    &c00,&c01,&c02,&c10,&c11,&c12,&c20,&c21,&c22,
    &gst,
    &p00,&p01,&p02,&p10,&p11,&p12,&p20,&p21,&p22);
    if (ok) {
        rc2i[0][0] = c00;
        rc2i[0][1] = c01;
        rc2i[0][2] = c02;
        rc2i[1][0] = c10;
        rc2i[1][1] = c11;
        rc2i[1][2] = c12;
        rc2i[2][0] = c20;
        rc2i[2][1] = c21;
        rc2i[2][2] = c22;
        rpom[0][0] = p00;
        rpom[0][1] = p01;
        rpom[0][2] = p02;
        rpom[1][0] = p10;
        rpom[1][1] = p11;
        rpom[1][2] = p12;
        rpom[2][0] = p20;
        rpom[2][1] = p21;
        rpom[2][2] = p22;
        eraC2teqx(rc2i, gst, rpom, rc2t);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2t[0][0],rc2t[0][1],rc2t[0][2],
        rc2t[1][0],rc2t[1][1],rc2t[1][2],
        rc2t[2][0],rc2t[2][1],rc2t[2][2]
        );
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_c2teqx_doc,
"c2teqx(rbpn, gst, rpom -> rc2t\n\n"
"Assemble the celestial to terrestrial matrix from equinox-based\n"
"components (the celestial-to-true matrix, the Greenwich Apparent\n"
"Sidereal Time and the polar motion matrix).\n"
"Given:\n"
"    rbpn       celestial-to-true matrix\n"
"    gst        Greenwich (apparent) Sidereal Time\n"
"    rpom       polar-motion matrix\n"
"Returned:\n"
"    rc2t       the celestial-to-terrestrial matrix");

static PyObject *
erfa_c2tpe(PyObject *self, PyObject *args)
{
    double tta,ttb,uta,utb,dpsi,deps,xp,yp,rc2t[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dddddddd", &tta,&ttb,&uta,&utb,&dpsi,&deps,&xp,&yp);
    if (ok) {
        eraC2tpe(tta,ttb,uta,utb,dpsi,deps,xp,yp,rc2t);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2t[0][0],rc2t[0][1],rc2t[0][2],
        rc2t[1][0],rc2t[1][1],rc2t[1][2],
        rc2t[2][0],rc2t[2][1],rc2t[2][2]
        ); 
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_c2tpe_doc,
"c2tpe(tta,ttb,uta,utb,dpsi,deps,xp,yp) -> rc2t\n\n"
"Form the celestial to terrestrial matrix given the date, the UT1,\n"
"the nutation and the polar motion.  IAU 2000.\n"
"Given:\n"
"    tta,ttb    TT as 2-part Julian Date\n"
"    uta,utb    UT1 as 2-part Julian Date\n"
"    dpsi,deps  nutation\n"
"    xp,yp      coordinates of the pole\n"
"Returned:\n"
"    rc2t       celestial-to-terrestrial matrix");

static PyObject *
erfa_c2txy(PyObject *self, PyObject *args)
{
    double tta,ttb,uta,utb,x,y,xp,yp,rc2t[3][3];
    int ok;
    ok = PyArg_ParseTuple(args, "dddddddd", &tta,&ttb,&uta,&utb,&x,&y,&xp,&yp);
    if (ok) {
        eraC2txy(tta,ttb,uta,utb,x,y,xp,yp,rc2t);
        return Py_BuildValue(
        "((ddd)(ddd)(ddd))",
        rc2t[0][0],rc2t[0][1],rc2t[0][2],
        rc2t[1][0],rc2t[1][1],rc2t[1][2],
        rc2t[2][0],rc2t[2][1],rc2t[2][2]
        ); 
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_c2txy_doc,
"c2txy(tta,ttb,uta,utb,x,y,xp,yp) -> rc2t\n\n"
"Form the celestial to terrestrial matrix given the date, the UT1,\n"
"the CIP coordinates and the polar motion.  IAU 2000."
"Given:\n"
"    tta,ttb    TT as 2-part Julian Date\n"
"    uta,utb    UT1 as 2-part Julian Date\n"
"    x,y        Celestial Intermediate Pole\n"
"    xp,yp      coordinates of the pole\n"
"Returned:\n"
"    rc2t       celestial-to-terrestrial matrix");

static PyObject *
erfa_cal2jd(PyObject *self, PyObject *args)
{
    int iy, im, id;
    double dmj0, dmj;
    int ok, status;
    ok = PyArg_ParseTuple(args, "iii", &iy, &im, &id);
    if (ok) {
        status = eraCal2jd(iy, im, id, &dmj0, &dmj);
    }
    else {
        return NULL;
    }
    if (status < 0){
        if (status == -1){
            PyErr_SetString(erfaError, "bad year");
            return NULL;
        }
        else if (status == -2){
            PyErr_SetString(erfaError, "bad month");
            return NULL;
        }
        else if (status == -3){
            PyErr_SetString(erfaError, "bad day");
            return NULL;
        }
    }
    return Py_BuildValue("dd", dmj0, dmj);
}

PyDoc_STRVAR(erfa_cal2jd_doc,
"cal2jd(year, month, day) -> 2400000.5,djm\n\n"
"Gregorian Calendar to Julian Date.\n"
"Given:\n"
"    year, month      day in Gregorian calendar\n"
"Returned:\n"
"    2400000.5,djm    MJD zero-point and Modified Julian Date for 0 hrs");

static PyObject *
erfa_d2dtf(PyObject *self, PyObject *args)
{
    double d1, d2;
    int n, ok, status;
    int iy, im, id, ihmsf[4];
    char *scale = "UTC";
    ok = PyArg_ParseTuple(args, "idd|s", &n, &d1, &d2, &scale);
    if (ok) {
        status = eraD2dtf(scale, n, d1, d2, &iy, &im, &id, ihmsf);
    }
    else {
        return NULL;
    }
    if (status > 0){
        PyErr_SetString(erfaError, "doubious year: date predates UTC or too far in the future");
        return NULL;
    }
    if (status < 0){
        PyErr_SetString(erfaError, "unaceptable date");
        return NULL;
    }
    return Py_BuildValue("iiiiiii", iy,im,id,ihmsf[0],ihmsf[1],ihmsf[2],ihmsf[3]);
}

PyDoc_STRVAR(erfa_d2dtf_doc,
"d2dtf(n, d1, d2 | scale) -> year, month, day, hours, minutes, seconds, fraction\n\n"
"Format for output a 2-part Julian Date (or in the case of UTC a\n"
"quasi-JD form that includes special provision for leap seconds).\n"
"Given:\n"
"    ndp        resolution\n"
"    d1,d2      time as a 2-part Julian Date\n"
"    scale      optional time scale ID, default to ''UTC''\n"
"Returned:\n"
"   iy,im,id    year, month, day in Gregorian calendar\n"
"   ihmsf       hours, minutes, seconds, fraction\n"

);

static PyObject *
erfa_dat(PyObject *self, PyObject *args)
{
    int y, m, d, status;
    double fd, deltat;
    if (! PyArg_ParseTuple(args, "iiid", &y, &m, &d, &fd))
        return NULL;
    status = eraDat(y,m,d,fd,&deltat);
    if (status > 0){
        PyErr_SetString(erfaError, "doubious year: date before UTC:1960 January 1.0.");
        return NULL;
    }
    else if (status < 0){
        if (status == -1){
            PyErr_SetString(erfaError, "unaceptable date, bad year");
            return NULL;
        }
        else if (status == -2){
            PyErr_SetString(erfaError, "unaceptable date, bad month");
            return NULL;
        }
        else if (status == -3){
            PyErr_SetString(erfaError, "unaceptable date, bad day");
            return NULL;
        }      
        else if (status == -4){
            PyErr_SetString(erfaError, "bad fraction day, should be < 1.");
            return NULL;
        }      
    }
    return PyFloat_FromDouble(deltat);
}

PyDoc_STRVAR(erfa_dat_doc,
"dat(y,m,d,fd) -> delta t\n\n"
"For a given UTC date, calculate delta(AT) = TAI-UTC.\n"
"Given:\n"
"    y          year\n"
"    m          month\n"
"    d          day\n"
"    fd         fraction of day\n"
"Returned:\n"
"    deltat     TAI minus UTC, seconds");

static PyObject *
erfa_dtdb(PyObject *self, PyObject *args)
{
    double tdbtt, jd1, jd2, ut1, elon, u, v;
    if (! PyArg_ParseTuple(args, "dddddd", &jd1, &jd2, &ut1, &elon, &u, &v))
        return NULL;
    tdbtt = eraDtdb(jd1, jd2, ut1, elon, u, v);
    return PyFloat_FromDouble(tdbtt);
}

PyDoc_STRVAR(erfa_dtdb_doc,
"dtdb(d1, d2, ut, elon, u, v) -> TDB-TT\n\n"
"An approximation to TDB-TT, the difference between barycentric\n"
"dynamical time and terrestrial time, for an observer on the Earth.\n"
"Given:\n"
"    d1,d2      TDB date as 2-part Julian Date\n"
"    ut         UT1 universal time as fraction of one day\n"
"    elong      longitude (east positive, radians)\n"
"    u          distance from Earth spin axis (km)\n"
"    v          distance north of equatorial plane (km)");

static PyObject *
erfa_dtf2d(PyObject *self, PyObject *args)
{
    int iy, im, id, ihr, imn, ok, status;
    double sec, dj1, dj2;
    char *scale = "UTC";
    ok = PyArg_ParseTuple(args, "iiiiid|s", &iy,&im,&id,&ihr,&imn,&sec,&scale);
    if (ok) {
        status = eraDtf2d(scale, iy, im, id, ihr, imn, sec, &dj1, &dj2);
        if (status == 3) {
            PyErr_SetString(erfaError, "both of next two");
            return NULL;
        }
        else if (status == 2) {
            PyErr_SetString(erfaError, "time is after end of day");
            return NULL;
        }
        else if (status == 1) {
            PyErr_SetString(erfaError, "dubious year");
            return NULL;
        }
        else if (status == -1) {
            PyErr_SetString(erfaError, "bad year");
            return NULL;
        }
        else if (status == -2) {
            PyErr_SetString(erfaError, "bad month");
            return NULL;
        }
        else if (status == -3) {
            PyErr_SetString(erfaError, "bad day");
            return NULL;
        }
        else if (status == -4) {
            PyErr_SetString(erfaError, "bad hour");
            return NULL;
        }
        else if (status == -5) {
            PyErr_SetString(erfaError, "bad minute");
            return NULL;
        }
        else if (status == -6) {
            PyErr_SetString(erfaError, "bad second (<0)");
            return NULL;
        }
        else {
            return  Py_BuildValue("dd",dj1, dj2);
        }
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_dtf2d_doc,
"dtf2d(y,m,d,hr,mn,sec |scale) -> d1,d2\n\n"
"Encode date and time fields into 2-part Julian Date (or in the case\n"
"of UTC a quasi-JD form that includes special provision for leap seconds).\n"
"Given:\n"
"    y,m,d   year, month, day in Gregorian calendar\n"
"    hr,mn   hour, minute\n"
"    sec     seconds\n"
"    scale   optional time scale ID, default ''UTC''\n"
"Returned:\n"
"   d1,d2     2-part Julian Date");

static PyObject *
erfa_ee00(PyObject *self, PyObject *args)
{
    double d1,d2,epsa,dpsi,ee;
    int ok;
    ok = PyArg_ParseTuple(args, "dddd", &d1,&d2,&epsa,&dpsi);
    if (ok) {
        ee = eraEe00(d1,d2,epsa,dpsi);
        return PyFloat_FromDouble(ee);
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_ee00_doc,
"ee00(d1,d2,epsa,dpsi) -> ee\n\n"
"The equation of the equinoxes compatible with IAU 2000 resolutions,\n"
"given the nutation in longitude and the mean obliquity.\n"
"Given:\n"
"    d1,d2      TT as 2-part Julian Date\n"
"    epsa       mean obliquity\n"
"    dpsi       nutation in longitude\n"
"Returned:\n"
"    ee         equation of the equinoxes");

static PyObject *
erfa_ee00a(PyObject *self, PyObject *args)
{
    double d1,d2,ee;
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1,&d2);
    if (ok) {
        ee = eraEe00a(d1,d2);
        return Py_BuildValue("d",ee); //PyFloat_FromDouble(ee);
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_ee00a_doc,
"ee00a(d1,d2) -> ee\n\n"
"equation of the equinoxes, compatible with IAU 2000 resolutions.\n"
"Given:\n"
"    d1,d2      TT as 2-part Julian Date\n"
"Returned:\n"
"    ee         equation of the equinoxes");

static PyObject *
erfa_ee00b(PyObject *self, PyObject *args)
{
    double d1,d2,ee;
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1,&d2);
    if (ok) {
        ee = eraEe00b(d1,d2);
        return Py_BuildValue("d",ee); //PyFloat_FromDouble(ee);
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_ee00b_doc,
"ee00b(d1,d2) -> ee\n\n"
"Equation of the equinoxes, compatible with IAU 2000 resolutions\n"
"but using the truncated nutation model IAU 2000B.\n"
"Given:\n"
"    d1,d2      TT as 2-part Julian Date\n"
"Returned:\n"
"    ee         equation of the equinoxes");

static PyObject *
erfa_ee06a(PyObject *self, PyObject *args)
{
    double d1,d2,ee;
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1,&d2);
    if (ok) {
        ee = eraEe06a(d1,d2);
        return Py_BuildValue("d",ee); //PyFloat_FromDouble(ee);
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_ee06a_doc,
"ee06a(d1,d2) -> ee\n\n"
"Equation of the equinoxes, compatible with IAU 2000 resolutions and \n"
"IAU 2006/2000A precession-nutation.\n"
"Given:\n"
"    d1,d2      TT as 2-part Julian Date\n"
"Returned:\n"
"    ee         equation of the equinoxes");

static PyObject *
erfa_eect00(PyObject *self, PyObject *args)
{
    double d1,d2,ct;
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1,&d2);
    if (ok) {
        ct = eraEect00(d1,d2);
        return Py_BuildValue("d",ct);
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_eect00_doc,
"eect00(d1,d2) -> ct\n\n"
"Equation of the equinoxes complementary terms, consistent with IAU 2000 resolutions.\n"
"Given:\n"
"    d1,d2      TT as 2-part Julian Date\n"
"Returned:\n"
"    ct        complementary terms");

static PyObject *
erfa_eform(PyObject *self, PyObject *args)
{
    int n, status, ok;
    double a, f;
    ok = PyArg_ParseTuple(args, "i", &n);
    if (ok) {
        status = eraEform(n, &a, &f);
        if (status) {
            PyErr_SetString(erfaError, "illegal identifier; n should be 1,2 or 3");
            return NULL;
        }            
        return Py_BuildValue("dd",a,f);
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_eform_doc,
"eform(n) -> a, f\n\n"
"Earth reference ellipsoids.\n"
"Given:\n"
"    n          ellipsoid identifier\n\n"
"    n=1 WGS84\n"
"    n=2 GRS80\n"
"    n=3 WGS72\n"
"Returned:\n"
"    a          equatorial radius (meters)\n"
"    f          flattening");

static PyObject *
erfa_eo06a(PyObject *self, PyObject *args)
{
    double d1, d2, eo;
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1, &d2);
    if (ok) {
        eo = eraEo06a(d1, d2);
        return Py_BuildValue("d",eo);
    }
    else {
        return NULL;
    }    
}

PyDoc_STRVAR(erfa_eo06a_doc,
"eo06a(d1, d2) -> eo\n"
"equation of the origins, IAU 2006 precession and IAU 2000A nutation.\n"
"Given:\n"
"    d1,d2      TT as 2-part Julian Date\n"
"Returned:\n"
"    eo         equation of the origins (radians)");

static PyObject *
erfa_eors(PyObject *self, PyObject *args)
{
    double rnpb[3][3], s, eo;
    double r00,r01,r02,r10,r11,r12,r20,r21,r22;
    int ok;
    ok = PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))d",
        &r00,&r01,&r02,&r10,&r11,&r12,&r20,&r21,&r22,&s);
    if (ok) {
        rnpb[0][0] = r00;
        rnpb[0][1] = r01;
        rnpb[0][2] = r02;
        rnpb[1][0] = r10;
        rnpb[1][1] = r11;
        rnpb[1][2] = r12;
        rnpb[2][0] = r20;
        rnpb[2][1] = r21;
        rnpb[2][2] = r22;
        eo = eraEors(rnpb, s);
        return Py_BuildValue("d",eo);
    }
    else {
        return NULL;
    }        
}

PyDoc_STRVAR(erfa_eors_doc,
"eors(rnpb, s) -> eo\n\n"
"Equation of the origins, given the classical NPB matrix and the quantity s\n"
"Given:\n"
"    rnpb       classical nutation x precession x bias matrix\n"
"    s          CIO locator\n"
"Returned:\n"
"    eo         equation of the origins in radians");

static PyObject *
erfa_epb(PyObject *self, PyObject *args)
{
    double d1, d2, b;
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1, &d2);
    if (ok) {
        b = eraEpb(d1,d2);
        return Py_BuildValue("d",b);
    }
    else {
        return NULL;
    }                
}

PyDoc_STRVAR(erfa_epb_doc,
"epb(d1, d2) -> b\n\n"
"Julian Date to Besselian Epoch\n"
"Given:\n"
"    d1,d2      2-part Julian Date\n"
"Returned:\n"
"    b          Besselian Epoch.");

static PyObject *
erfa_epb2jd(PyObject *self, PyObject *args)
{
    double epb, jd0, jd1;
    if (!PyArg_ParseTuple(args, "d", &epb)) {
        return NULL;
    }
    eraEpb2jd(epb, &jd0, &jd1);
    return Py_BuildValue("dd",jd0,jd1);
}

PyDoc_STRVAR(erfa_epb2jd_doc,
"epb2jd(epb) -> 2400000.5 djm\n\n"
"Given:\n"
"    epb        Besselian Epoch,\n"
"Returned:\n"
"    djm        Modified Julian Date");

static PyObject *
erfa_epj(PyObject *self, PyObject *args)
{
    double d1, d2, j;
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1, &d2);
    if (ok) {
        j = eraEpj(d1,d2);
        return Py_BuildValue("d",j);
    }
    else {
        return NULL;
    }                
}

PyDoc_STRVAR(erfa_epj_doc,
"epj(d1, d2) -> b\n\n"
"Julian Date to Julian Epoch.\n"
"Given:\n"
"    d1,d2      2-part Julian Date\n"
"Returned:\n"
"    b          Julian Epoch");

static PyObject *
erfa_epj2jd(PyObject *self, PyObject *args)
{
    double epj, jd0, jd1;
    if (!PyArg_ParseTuple(args, "d", &epj)) {
        return NULL;
    }
    eraEpj2jd(epj, &jd0, &jd1);
    return Py_BuildValue("dd",jd0,jd1);
}

PyDoc_STRVAR(erfa_epj2jd_doc,
"epj2jd(epj) -> 2400000.5 djm\n\n"
"Julian Epoch to Julian Date\n"
"Given:\n"
"    epj        Julian Epoch\n"
"Returned:\n"
"    djm        Modified Julian Date");

static PyObject *
erfa_epv00(PyObject *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"dj1", "dj2", "prec", NULL};
    double dj1, dj2;
    int prec=0;
    double pvh[2][3] = {{[0]=0.},{[0]=0.}};
    double pvb[2][3] = {{[0]=0.},{[0]=0.}};
    if (! PyArg_ParseTupleAndKeywords(args, kwds, "dd|i", kwlist, &dj1, &dj2, &prec))
        return NULL;
    int status = eraEpv00(dj1, dj2, pvh, pvb);
    if (status == 1) {
        if (prec) {
           return Py_BuildValue("((ddd)(ddd))((ddd)(ddd))",
           pvh[0][0], pvh[0][1], pvh[0][2], pvh[1][0], pvh[1][1], pvh[1][2],
           pvb[0][0], pvb[0][1], pvb[0][2], pvb[1][0], pvb[1][1], pvb[1][2]);
        }
        else {
            PyErr_SetString(erfaError, "date outside the range 1900-2100 AD, set prec!=0 to force computation");
            return NULL;
        }
    }
    else {
       return Py_BuildValue("((ddd)(ddd))((ddd)(ddd))",
       pvh[0][0], pvh[0][1], pvh[0][2], pvh[1][0], pvh[1][1], pvh[1][2],
       pvb[0][0], pvb[0][1], pvb[0][2], pvb[1][0], pvb[1][1], pvb[1][2]);
    }
}

PyDoc_STRVAR(erfa_epv00_doc,
"epv00(d1,d2 |prec=0) -> pvh, pvb\n\n"
"Earth position and velocity, heliocentric and barycentric,\n"
"with respect to the Barycentric Celestial Reference System.\n"
"Given:\n"
"    d1,d2      TDB as 2-part Julian Date\n"
"    prec!=0    optional, to compute outside the range 1900-2100 AD\n"
"Returned:\n"
"    pvh        heliocentric Earth position/velocity\n"
"    pvb        barycentric Earth position/velocity");

static PyObject *
erfa_eqeq94(PyObject *self, PyObject *args)
{
    double d1, d2, ee;
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1, &d2);
    if (ok) {
        ee = eraEqeq94(d1,d2);
        return Py_BuildValue("d",ee);
    }
    else {
        return NULL;
    }   
}

PyDoc_STRVAR(erfa_eqeq94_doc,
"eqeq94(d1,d2) -> ee\n\n"
"Equation of the equinoxes, IAU 1994 model.\n"
"Given:\n"
"    d1,d2      TDB as 2-part Julian Date\n"
"Returned:\n"
"    ee         equation of the equinoxes");

static PyObject *
erfa_era00(PyObject *self, PyObject *args)
{
    double d1, d2, era;
    int ok;
    ok = PyArg_ParseTuple(args, "dd", &d1, &d2);
    if (ok) {
        era = eraEra00(d1,d2);
        return Py_BuildValue("d",era);
    }
    else {
        return NULL;
    }   
}

PyDoc_STRVAR(erfa_era00_doc,
"era00(d1,d2) -> era\n\n"
"Earth rotation angle (IAU 2000 model).\n"
"Given:\n"
"    d1,d2      UT1 as 2-part Julian Date (d1,d2)\n"
"Returned:\n"
"    era         Earth rotation angle (radians), range 0-2pi");

static PyObject *
erfa_fad03(PyObject *self, PyObject *args)
{
    double t,d;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    d = eraFad03(t);
    return Py_BuildValue("d",d);
}

PyDoc_STRVAR(erfa_fad03_doc,
"fad03(t) -> d\n\n"
"Fundamental argument, IERS Conventions (2003):\n"
"mean elongation of the Moon from the Sun.\n"
"Given:\n"
"    t          TDB as Julian centuries since J2000.0\n"
"Returned:\n"
"    d          mean elongation of the Moon from the Sun, radians.");

static PyObject *
erfa_fae03(PyObject *self, PyObject *args)
{
    double t,e;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    e = eraFae03(t);
    return Py_BuildValue("d",e);
}

PyDoc_STRVAR(erfa_fae03_doc,
"fae03(t) -> e\n\n"
"Fundamental argument, IERS Conventions (2003):\n"
"mean longitude of Earth.\n"
"Given:\n"
"    t          TDB as Julian centuries since J2000.0\n"
"Returned:\n"
"    e          mean longitude of Earth, radians.");

static PyObject *
erfa_faf03(PyObject *self, PyObject *args)
{
    double t,f;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    f = eraFaf03(t);
    return Py_BuildValue("d",f);
}

PyDoc_STRVAR(erfa_faf03_doc,
"faf03(t) -> f\n\n"
"Fundamental argument, IERS Conventions (2003):\n"
"mean longitude of the Moon minus mean longitude of the ascending node.\n"
"Given:\n"
"    t          TDB as Julian centuries since J2000.0\n"
"Returned:\n"
"    f          mean longitude of the Moon minus mean longitude of the ascending node, radians.");

static PyObject *
erfa_faju03(PyObject *self, PyObject *args)
{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFaju03(t);
    return Py_BuildValue("d",l);
}

PyDoc_STRVAR(erfa_faju03_doc,
"faju03(t) -> l\n\n"
"Fundamental argument, IERS Conventions (2003):\n"
"mean longitude of Jupiter.\n"
"Given:\n"
"    t          TDB as Julian centuries since J2000.0\n"
"Returned:\n"
"    l          mean longitude of Jupiter., in radians.");

static PyObject *
erfa_fal03(PyObject *self, PyObject *args)
{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFal03(t);
    return Py_BuildValue("d",l);
}

PyDoc_STRVAR(erfa_fal03_doc,
"fal03(t) -> l\n\n"
"Fundamental argument, IERS Conventions (2003):\n"
"mean anomaly of the Moon.\n"
"Given:\n"
"    t          TDB as Julian centuries since J2000.0\n"
"Returned:\n"
"    l          mean anomaly of the Moon, in radians.");

static PyObject *
erfa_falp03(PyObject *self, PyObject *args)
{
    double t,lp;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    lp = eraFalp03(t);
    return Py_BuildValue("d",lp);
}

PyDoc_STRVAR(erfa_falp03_doc,
"fal03(t) -> lp\n"
"Fundamental argument, IERS Conventions (2003):\n"
"mean anomaly of the Sun.\n"
"Given:\n"
"    t          TDB as Julian centuries since J2000.0\n"
"Returned\n"
"    lp         mean anomaly of the Sun, in radians.");

static PyObject *
erfa_fama03(PyObject *self, PyObject *args)
{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFama03(t);
    return Py_BuildValue("d",l);
}

PyDoc_STRVAR(erfa_fama03_doc,
"fama03(t) -> l\n\n"
"Fundamental argument, IERS Conventions (2003):\n"
"mean longitude of Mars.\n"
"Given:\n"
"    t          TDB as Julian centuries since J2000.0\n"
"Returned:\n"
"    l          mean longitude of Mars, in radians.");

static PyObject *
erfa_fame03(PyObject *self, PyObject *args)
{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFame03(t);
    return Py_BuildValue("d",l);
}

PyDoc_STRVAR(erfa_fame03_doc,
"fame03(t) -> l\n\n"
"Fundamental argument, IERS Conventions (2003):\n"
"mean longitude of Mercury.\n"
"Given:\n"
"    t          TDB as Julian centuries since J2000.0\n"
"Returned:\n"
"    l          mean longitude of Mercury, in radians.");

static PyObject *
erfa_fane03(PyObject *self, PyObject *args)
{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFane03(t);
    return Py_BuildValue("d",l);
}

PyDoc_STRVAR(erfa_fane03_doc,
"fane03(t) -> l\n\n"
"Fundamental argument, IERS Conventions (2003):\n"
"mean longitude of Neptune.\n"
"Given:\n"
"    t          TDB as Julian centuries since J2000.0\n"
"Returned:\n"
"    l          mean longitude of Neptune, in radians.");

static PyObject *
erfa_faom03(PyObject *self, PyObject *args)
{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFaom03(t);
    return Py_BuildValue("d",l);
}

PyDoc_STRVAR(erfa_faom03_doc,
"faom03(t) -> l\n\n"
"Fundamental argument, IERS Conventions (2003):\n"
"mean longitude of Moon's ascending node.\n"
"Given:\n"
"    t          TDB as Julian centuries since J2000.0\n"
"Returned:\n"
"    l          mean longitude of Moon's ascending node, in radians.");

static PyObject *
erfa_fapa03(PyObject *self, PyObject *args)
{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFapa03(t);
    return Py_BuildValue("d",l);
}

PyDoc_STRVAR(erfa_fapa03_doc,
"fapa03(t) -> l\n\n"
"Fundamental argument, IERS Conventions (2003):\n"
"general accumulated precession in longitude.\n"
"Given:\n"
"    t          TDB as Julian centuries since J2000.0\n"
"Returned:\n"
"    l          general accumulated precession in longitude, in radians.");

static PyObject *
erfa_fasa03(PyObject *self, PyObject *args)
{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFasa03(t);
    return Py_BuildValue("d",l);
}

PyDoc_STRVAR(erfa_fasa03_doc,
"fasa03(t) -> l\n"
"Fundamental argument, IERS Conventions (2003):\n"
"mean longitude of Saturn.\n"
"Given:\n"
"    t          TDB as Julian centuries since J2000.0\n"
"Returned:\n"
"    l          , in radians, the mean longitude of Saturn.");

static PyObject *
erfa_faur03(PyObject *self, PyObject *args)
{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFaur03(t);
    return Py_BuildValue("d",l);
}

PyDoc_STRVAR(erfa_faur03_doc,
"faur03(t) -> l\n\n"
"Fundamental argument, IERS Conventions (2003):\n"
"mean longitude of Uranus.\n"
"Given:\n"
"    t          TDB as Julian centuries since J2000.,\n"
"Returned:\n"
"    l          mean longitude of Uranus, in radians.");

static PyObject *
erfa_fave03(PyObject *self, PyObject *args)
{
    double t,l;
    if (!PyArg_ParseTuple(args, "d", &t)) {
        return NULL;
    }
    l = eraFave03(t);
    return Py_BuildValue("d",l);
}

PyDoc_STRVAR(erfa_fave03_doc,
"faver03(t) -> l\n\n"
"Fundamental argument, IERS Conventions (2003)\n"
"mean longitude of Venus."
"Given:\n"
"    t          TDB as Julian centuries since J2000.0\n"
"Returned:\n"
"    l          mean longitude of Venus, in radians.");

static PyObject *
erfa_fk52h(PyObject *self, PyObject *args)
{
    double r5, d5, dr5, dd5, px5, rv5;
    double rh, dh, drh, ddh, pxh, rvh;
    if (!PyArg_ParseTuple(args, "dddddd", &r5, &d5, &dr5, &dd5, &px5, &rv5)) {
        return NULL;
    }
    eraFk52h(r5, d5, dr5, dd5, px5, rv5,
             &rh, &dh, &drh, &ddh, &pxh, &rvh);
    return Py_BuildValue("dddddd", rh, dh, drh, ddh, pxh, rvh);
}

PyDoc_STRVAR(erfa_fk52h_doc,
"fk52h(r5, d5, dr5, dd5, px5,rv5) -> rh, dh, drh, ddh, pxh, rvh)\n\n"
"Transform FK5 (J2000.0) star data into the Hipparcos system.\n"
"Given (all FK5, equinox J2000.0, epoch J2000.0):\n"
"    r5         RA (radians)\n"
"    d5         Dec (radians)\n"
"    dr5        proper motion in RA (dRA/dt, rad/Jyear)\n"
"    dd5        proper motion in Dec (dDec/dt, rad/Jyear)\n"
"    px5        parallax (arcsec)\n"
"    rv5        radial velocity (km/s, positive = receding)\n"
"Returned (all Hipparcos, epoch J2000.0):\n"
"    rh         RA (radians)\n"
"    dh         Dec (radians)\n"
"    drh        proper motion in RA (dRA/dt, rad/Jyear)\n"
"    ddh        proper motion in Dec (dDec/dt, rad/Jyear)\n"
"    pxh        parallax (arcsec)\n"
"    rvh        radial velocity (km/s, positive = receding)");

static PyObject *
erfa_fk5hip(PyObject *self)
{
    double r5h[3][3], s5h[3];
    eraFk5hip(r5h, s5h);
    return Py_BuildValue("((ddd)(ddd)(ddd))(ddd)",
        r5h[0][0],r5h[0][1],r5h[0][2],
        r5h[1][0],r5h[1][1],r5h[1][2],
        r5h[2][0],r5h[2][1],r5h[2][2],
        s5h[0],s5h[1],s5h[2]);
}

PyDoc_STRVAR(erfa_fk5hip_doc,
"fk5hip() -> r5h, s5h\n\n"
"FK5 to Hipparcos rotation and spin.\n"
"Returned:\n"
"    r5h        r-matrix: FK5 rotation wrt Hipparcos \n"
"    s5         r-vector: FK5 spin wrt Hipparcos.");

static PyObject *
erfa_fk5hz(PyObject *self, PyObject *args)
{
    double r5,d5,d1,d2,rh,dh;
    if (!PyArg_ParseTuple(args, "dddd", &r5, &d5, &d1, &d2)) {
        return NULL;
    }
    eraFk5hz(r5,d5,d1,d2,&rh,&dh);
    return Py_BuildValue("dd", rh, dh);
}

PyDoc_STRVAR(erfa_fk5hz_doc,
"fk5hz(r5, d5, d1, d2) -> rh, dh\n\n"
"Transform an FK5 (J2000.0) star position into the system of the\n"
"Hipparcos catalogue, assuming zero Hipparcos proper motion.\n"
"Given:\n"
"    r5         FK5 RA (radians), equinox J2000.0, at date d1,d2\n"
"    d5         FK5 Dec (radians), equinox J2000.0, at date d1,d2\n"
"    d1,d2      TDB date \n"
"Returned:\n"
"    rh         Hipparcos RA (radians)\n"
"    dh         Hipparcos Dec (radians)");

static PyObject *
erfa_fw2m(PyObject *self, PyObject *args)
{
    double gamb, phib, psi, eps, r[3][3];
    if (!PyArg_ParseTuple(args, "dddd", &gamb, &phib, &psi, &eps)) {
        return NULL;
    }
    eraFw2m(gamb, phib, psi, eps, r);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
        r[0][0],r[0][1],r[0][2],
        r[1][0],r[1][1],r[1][2],
        r[2][0],r[2][1],r[2][2]);
}

PyDoc_STRVAR(erfa_fw2m_doc,
"fw2m(gamb, phib, psi, eps) -> r\n\n"
"Form rotation matrix given the Fukushima-Williams angles.\n"
"Given:\n"
"    gamb       F-W angle gamma_bar (radians)\n"
"    phib       F-W angle phi_bar (radians)\n"
"    si         F-W angle psi (radians)\n"
"    eps        F-W angle epsilon (radians)\n"
"  Returned:\n"
"    r          rotation matrix");

static PyObject *
erfa_fw2xy(PyObject *self, PyObject *args)
{
    double gamb, phib, psi, eps, x, y;
    if (!PyArg_ParseTuple(args, "dddd", &gamb, &phib, &psi, &eps)) {
        return NULL;
    }
    eraFw2xy(gamb, phib, psi, eps, &x, &y);
    return Py_BuildValue("dd", x, y);
}

PyDoc_STRVAR(erfa_fw2xy_doc,
"fw2xy(gamb, phib, psi, eps) -> x, y\n\n"
"Form CIP X,Y given Fukushima-Williams bias-precession-nutation angles.\n"
"Given:\n"
"    gamb       F-W angle gamma_bar (radians)\n"
"    phib       F-W angle phi_bar (radians)\n"
"    psi        F-W angle psi (radians)\n"
"    eps        F-W angle epsilon (radians)\n"
"  Returned:\n"
"    x,y        CIP X,Y (radians)");

static PyObject *
erfa_gc2gd(PyObject *self, PyObject *args)
{
    double xyz[3], elong, phi, height;
    double x, y, z;
    int ok, n, status;
    ok = PyArg_ParseTuple(args, "n(ddd)", &n, &x, &y, &z);
    if (ok) {
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
        status = eraGc2gd(n, xyz, &elong, &phi, &height);
        if (status == -1) {
            PyErr_SetString(erfaError, "illegal identifier; n should be 1,2 or 3");
            return NULL;
        }            
        else if (status == -2) {
            PyErr_SetString(erfaError, "internal error");
            return NULL;
        }  
        else {
            return Py_BuildValue("ddd", elong, phi, height);
        }
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_gc2gd_doc,
"gc2gd(n, xyz[3]) -> elong, phi, height\n\n"
"Transform geocentric coordinates to geodetic\n"
"using the specified reference ellipsoid.\n"
"Given:\n"
"    n          ellipsoid identifier\n"
"    xyz        geocentric vector\n"
"Returned:\n"
"    elong      longitude (radians, east +ve)\n"
"    phi        latitude (geodetic, radians)\n"
"    height     height above ellipsoid (geodetic)");

static PyObject *
erfa_gc2gde(PyObject *self, PyObject *args)
{
    double xyz[3], a, f, elong, phi, height;
    double x, y, z;
    int ok, status;
    ok = PyArg_ParseTuple(args, "dd(ddd)", &a, &f, &x, &y, &z);
    if (ok) {
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
        status = eraGc2gde(a, f, xyz, &elong, &phi, &height);
        if (status == -1) {
            PyErr_SetString(erfaError, "illegal f");
            return NULL;
        }            
        else if (status == -2) {
            PyErr_SetString(erfaError, "illegal a");
            return NULL;
        }  
        else {
            return Py_BuildValue("ddd", elong, phi, height);
        }
    }
    else {
        return NULL;
    }    
}

PyDoc_STRVAR(erfa_gc2gde_doc,
"gc2gde(a, f, xyz) -> elong, phi, height\n\n"
"Transform geocentric coordinates to geodetic\n"
" for a reference ellipsoid of specified form.\n"
"Given:\n"
"    a          equatorial radius\n"
"    f          flattening\n"
"    xyz        geocentric vector\n"
"Returned:\n"
"    elong      longitude (radians, east +ve)\n"
"    phi        latitude (geodetic, radians)\n"
"    height     height above ellipsoid (geodetic)");

static PyObject *
erfa_gd2gc(PyObject *self, PyObject *args)
{
    double elong, phi, height, xyz[3];
    int n, status, ok;
    ok = PyArg_ParseTuple(args, "iddd", &n, &elong, &phi, &height);
    if (ok) {
        status = eraGd2gc(n, elong, phi, height, xyz);
        if (status == -1) {
            PyErr_SetString(erfaError, "illegal identifier; n should be 1,2 or 3");
            return NULL;
        }            
        else if (status == -2) {
            PyErr_SetString(erfaError, "illegal case");
            return NULL;
        }  
        else {
            return Py_BuildValue("ddd", xyz[0], xyz[1], xyz[2]);
        }
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_gd2gc_doc,
"gd2gc(n, elong, phi, height) -> xyz\n\n"
"Transform geodetic coordinates to geocentric\n"
" using the specified reference ellipsoid.\n"
"Given:\n"
"    n          ellipsoid identifier\n"
"    elong      longitude (radians, east +ve)\n"
"    phi        latitude (geodetic, radians)\n"
"    height     height above ellipsoid (geodetic in meters)\n"
"  Returned:\n"
"    xyz        geocentric vector in meters");

static PyObject *
erfa_gd2gce(PyObject *self, PyObject *args)
{
    double a, f, elong, phi, height, xyz[3];
    int status, ok;
    ok = PyArg_ParseTuple(args, "ddddd", &a, &f, &elong, &phi, &height);
    if (ok) {
        status = eraGd2gce(a, f, elong, phi, height, xyz);
        if (status == -1) {
            PyErr_SetString(erfaError, "illegal case");
            return NULL;
        }  
        else {
            return Py_BuildValue("ddd", xyz[0], xyz[1], xyz[2]);
        }
    }
    else {
        return NULL;
    }
}

PyDoc_STRVAR(erfa_gd2gce_doc,
"gd2gce(a, f, elong, phi, height) -> xyz\n\n"
"Transform geodetic coordinates to geocentric for a reference\n"
" for a reference ellipsoid of specified form\n"
"Given:\n"
"    a          equatorial radius in meters\n"
"    f          flattening\n"
"    elong      longitude (radians, east +ve)\n"
"    phi        latitude (geodetic, radians)\n"
"    height     height above ellipsoid (geodetic in meters)\n"
"  Returned:\n"
"    xyz        geocentric vector in meters");

static PyObject *
erfa_gmst00(PyObject *self, PyObject *args)
{
    double uta, utb, tta, ttb, g;
    if (!PyArg_ParseTuple(args, "dddd", &uta, &utb, &tta, &ttb)) {
        return NULL;
    }
    g = eraGmst00(uta, utb, tta, ttb);
    return Py_BuildValue("d", g);
}

PyDoc_STRVAR(erfa_gmst00_doc,
"gmst00(uta, utb, tta, ttb) -> gmst\n\n"
"Greenwich mean sidereal time\n"
"(model consistent with IAU 2000 resolutions).\n"
"Given:\n"
"    uta,utb    UT1 as a 2-part Julian Date\n"
"    tta,ttb    TT as a 2-part Julian Date\n"
"Returned:\n"
"    g          Greenwich mean sidereal time (radians)");

static PyObject *
erfa_gmst06(PyObject *self, PyObject *args)
{
    double uta, utb, tta, ttb, g;
    if (!PyArg_ParseTuple(args, "dddd", &uta, &utb, &tta, &ttb)) {
        return NULL;
    }
    g = eraGmst06(uta, utb, tta, ttb);
    return Py_BuildValue("d", g);
}

PyDoc_STRVAR(erfa_gmst06_doc,
"gmst06(uta, utb, tta, ttb) -> gmst\n\n"
"Greenwich mean sidereal time\n"
"(model consistent with IAU 2006resolutions).\n"
"Given:\n"
"    uta,utb    UT1 as a 2-part Julian Date\n"
"    tta,ttb    TT as a 2-part Julian Date\n"
"Returned:\n"
"    g          Greenwich mean sidereal time (radians)");

static PyObject *
erfa_gmst82(PyObject *self, PyObject *args)
{
    double dj1, dj2, g;
    if (!PyArg_ParseTuple(args, "dd", &dj1, &dj2)) {
        return NULL;
    }
    g = eraGmst82(dj1, dj2);
    return Py_BuildValue("d", g);
}

PyDoc_STRVAR(erfa_gmst82_doc,
"gmst82(d1, d2) -> gmst\n\n"
"Universal Time to Greenwich mean sidereal time (IAU 1982 model)\n"
"Given:\n"
"    d1,d2      UT1 as a 2-part Julian Date\n"
"Returned:\n"
"    g          Greenwich mean sidereal time (radians)");

static PyObject *
erfa_gst00a(PyObject *self, PyObject *args)
{
    double uta, utb, tta, ttb, g;
    if (!PyArg_ParseTuple(args, "dddd", &uta, &utb, &tta, &ttb)) {
        return NULL;
    }
    g = eraGst00a(uta, utb, tta, ttb);
    return Py_BuildValue("d", g);
}

PyDoc_STRVAR(erfa_gst00a_doc,
"gst00a(uta, utb, tta, ttb) -> gast\n\n"
"Greenwich apparent sidereal time\n"
"(model consistent with IAU 2000resolutions).\n"
"Given:\n"
"    uta,utb    UT1 as a 2-part Julian Date\n"
"    tta,ttb    TT as a 2-part Julian Date\n"
"Returned:\n"
"    g          Greenwich apparent sidereal time (radians)");

static PyObject *
erfa_gst00b(PyObject *self, PyObject *args)
{
    double uta, utb, g;
    if (!PyArg_ParseTuple(args, "dd", &uta, &utb)) {
        return NULL;
    }
    g = eraGst00b(uta, utb);
    return Py_BuildValue("d", g);
}

PyDoc_STRVAR(erfa_gst00b_doc,
"gst00b(uta, utb) -> gast\n\n"
"Greenwich apparent sidereal time (model consistent with IAU 2000\n"
"resolutions but using the truncated nutation model IAU 2000B).\n"
"Given:\n"
"    uta,utb    UT1 as a 2-part Julian Date\n"
"Returned:\n"
"    g          Greenwich apparent sidereal time (radians)");

static PyObject *
erfa_gst06(PyObject *self, PyObject *args)
{
    double uta, utb, tta, ttb, rnpb[3][3];
    double r00,r01,r02,r10,r11,r12,r20,r21,r22, g;
    int ok;
    ok = PyArg_ParseTuple(args, "dddd((ddd)(ddd)(ddd))",
        &uta, &utb, &tta, &ttb, 
        &r00,&r01,&r02,&r10,&r11,&r12,&r20,&r21,&r22);
    if (ok) {
        rnpb[0][0] = r00;
        rnpb[0][1] = r01;
        rnpb[0][2] = r02;
        rnpb[1][0] = r10;
        rnpb[1][1] = r11;
        rnpb[1][2] = r12;
        rnpb[2][0] = r20;
        rnpb[2][1] = r21;
        rnpb[2][2] = r22;
        g = eraGst06(uta, utb, tta, ttb, rnpb);
        return Py_BuildValue("d",g);
    }
    else {
        return NULL;
    }        
}

PyDoc_STRVAR(erfa_gst06_doc,
"gst06(uta, utb, tta, ttb, rnpb[3][3]) -> gast\n\n"
"Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.\n"
"Given:\n"
"    uta,utb    UT1 as a 2-part Julian Date\n"
"    tta,ttb    TT as a 2-part Julian Date\n"
"    rnpb       nutation x precession x bias matrix\n"
"Returned:\n"
"    g          Greenwich apparent sidereal time");

static PyObject *
erfa_gst06a(PyObject *self, PyObject *args)
{
    double uta, utb, tta, ttb, g;
    if (!PyArg_ParseTuple(args, "dddd", &uta, &utb, &tta, &ttb)) {
        return NULL;
    }
    g = eraGst06a(uta, utb, tta, ttb);
    return Py_BuildValue("d", g);
}

PyDoc_STRVAR(erfa_gst06a_doc,
"gst06a(uta, utb, tta, ttb) -> gast\n\n"
"Greenwich apparent sidereal time\n"
"(model consistent with IAU 2000and 2006 resolutions).\n"
"Given:\n"
"    uta,utb    UT1 as a 2-part Julian Date\n"
"    tta,ttb    TT as a 2-part Julian Date\n"
"Returned:\n"
"    g          Greenwich apparent sidereal time (radians)");

static PyObject *
erfa_gst94(PyObject *self, PyObject *args)
{
    double uta, utb, g;
    if (!PyArg_ParseTuple(args, "dd", &uta, &utb)) {
        return NULL;
    }
    g = eraGst94(uta, utb);
    return Py_BuildValue("d", g);
}

PyDoc_STRVAR(erfa_gst94_doc,
"gst94(uta, utb) -> gast\n\n"
"Greenwich apparent sidereal time\n"
"(consistent with IAU 1982/94 resolutions)\n"
"Given:\n"
"    uta,utb    UT1 as a 2-part Julian Date\n"
"Returned:\n"
"    g          Greenwich mean sidereal time (radians)");

static PyObject *
erfa_h2fk5(PyObject *self, PyObject *args)
{
    double rh, dh, drh, ddh, pxh, rvh;
    double r5, d5, dr5, dd5, px5, rv5;
    if (!PyArg_ParseTuple(args, "dddddd", &rh, &dh, &drh, &ddh, &pxh, &rvh)) {
        return NULL;
    }
    eraH2fk5(rh, dh, drh, ddh, pxh, rvh, &r5, &d5, &dr5, &dd5, &px5, &rv5);
    return Py_BuildValue("dddddd", r5, d5, dr5, dd5, px5, rv5);
}

PyDoc_STRVAR(erfa_h2fk5_doc,
"h2fk5(rh, dh, drh, ddh, pxh, rvh) -> r5, d5, dr5, dd5, px5, rv5\n\n"
"Transform Hipparcos star data into the FK5 (J2000.0) system.\n"
"Given (all Hipparcos, epoch J2000.0):\n"
"     rh    RA (radians)\n"
"     dh    Dec (radians)\n"
"     drh   proper motion in RA (dRA/dt, rad/Jyear)\n"
"     ddh   proper motion in Dec (dDec/dt, rad/Jyear)\n"
"     pxh   parallax (arcsec)\n"
"     rvh   radial velocity (km/s, positive = receding)\n"
"Returned (all FK5, equinox J2000.0, epoch J2000.0):\n"
"     r5    RA (radians)\n"
"     d5    Dec (radians)\n"
"     dr5   proper motion in RA (dRA/dt, rad/Jyear)\n"
"     dd5   proper motion in Dec (dDec/dt, rad/Jyear)\n"
"     px5   parallax (arcsec)\n"
"     rv5   radial velocity (km/s, positive = receding)");

static PyObject *
erfa_hfk5z(PyObject *self, PyObject *args)
{
    double rh, dh, d1, d2;
    double r5, d5, dr5, dd5;
    if (!PyArg_ParseTuple(args, "dddd", &rh, &dh, &d1, &d2)) {
        return NULL;
    }
    eraHfk5z(rh, dh, d1, d2, &r5, &d5, &dr5, &dd5);
    return Py_BuildValue("dddd", r5, d5, dr5, dd5);    
}

PyDoc_STRVAR(erfa_hfk5z_doc,
"hf5kz(rh, dh, d1, d2) -> r5, d5, dr5, dd5\n\n"
"Transform a Hipparcos star position into FK5 J2000.0, assuming\n"
"zero Hipparcos proper motion.\n"
"Given:\n"
"     rh    Hipparcos RA (radians)\n"
"     dh    Hipparcos Dec (radians)\n"
"     d1,d2 TDB date\n"
"Returned (all FK5, equinox J2000.0, date date1+date2):\n"
"     r5    RA (radians)\n"
"     d5    Dec (radians)\n"
"     dr5   FK5 RA proper motion (rad/year)\n"
"     dd5   Dec proper motion (rad/year)");

static PyObject *
erfa_jd2cal(PyObject *self, PyObject *args)
{
    double dj1, dj2, fd;
    int iy, im, id, status;
    if (!PyArg_ParseTuple(args, "dd", &dj1, &dj2)) {
        return NULL;
    }
    status = eraJd2cal(dj1, dj2, &iy, &im, &id, &fd);
    if (status) {
        PyErr_SetString(erfaError, "unacceptable date");
        return NULL;
    }
    return Py_BuildValue("iiid", iy, im, id, fd);    
}

PyDoc_STRVAR(erfa_jd2cal_doc,
"jd2cal(dj1, dj2) -> year, month, day, fraction of day\n\n"
"Julian Date to Gregorian year, month, day, and fraction of a day.\n"
"Given:\n"
"     dj1,dj2   2-part Julian Date\n"
"Returned:\n"
"     iy    year\n"
"     im    month\n"
"     id    day\n"
"     fd    fraction of day");

static PyObject *
erfa_jdcalf(PyObject *self, PyObject *args)
{
    double dj1, dj2;
    int ndp, status, iymdf[4];
    if (!PyArg_ParseTuple(args, "idd", &ndp, &dj1, &dj2)) {
        return NULL;
    }
    status = eraJdcalf(ndp, dj1, dj2, iymdf);
    if (status == -1) {
        PyErr_SetString(erfaError, "date out of range");
        return NULL;
    }
    if (status == 1) {
        PyErr_SetString(erfaError, "n not in 0-9 (interpreted as 0)");
        return NULL;
    }
    else {
        return Py_BuildValue("iiii", iymdf[0],iymdf[1],iymdf[2],iymdf[3]);
    }
}

PyDoc_STRVAR(erfa_jdcalf_doc,
"jdcalf(n, d1, d2) -> y, m, d, fd\n\n"
"Julian Date to Gregorian Calendar, expressed in a form convenient\n"
"for formatting messages:  rounded to a specified precision.\n"
"Given:\n"
"     n         number of decimal places of days in fraction\n"
"     d1,d2     2-part Julian Date\n"
"Returned:\n"
"     y         year\n"
"     m         month\n"
"     d         day\n"
"     fd        fraction of day");

static PyObject *
erfa_num00a(PyObject *self, PyObject *args)
{
    double d1, d2, rmatn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNum00a(d1, d2, rmatn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                            rmatn[0][0],rmatn[0][1],rmatn[0][2],
                            rmatn[1][0],rmatn[1][1],rmatn[1][2],
                            rmatn[2][0],rmatn[2][1],rmatn[2][2]);    
}

PyDoc_STRVAR(erfa_num00a_doc,
"num00a(d1, d2) -> rmatn\n\n"
"Form the matrix of nutation for a given date, IAU 2000A model.\n"
"Given:\n"
"     d1,d2     TT as a 2-part Julian Date\n"
"Returned:\n"
"     rmatn     nutation matrix");

static PyObject *
erfa_num00b(PyObject *self, PyObject *args)
{
    double d1, d2, rmatn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNum00b(d1, d2, rmatn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rmatn[0][0],rmatn[0][1],rmatn[0][2],
    rmatn[1][0],rmatn[1][1],rmatn[1][2],
    rmatn[2][0],rmatn[2][1],rmatn[2][2]);    
}

PyDoc_STRVAR(erfa_num00b_doc,
"num00b(d1, d2) -> rmatn\n\n"
"Form the matrix of nutation for a given date, IAU 2000B model.\n"
"Given:\n"
"     d1,d2     TT as a 2-part Julian Date\n"
"Returned:\n"
"     rmatn     nutation matrix");

static PyObject *
erfa_num06a(PyObject *self, PyObject *args)
{
    double d1, d2, rmatn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNum06a(d1, d2, rmatn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rmatn[0][0],rmatn[0][1],rmatn[0][2],
    rmatn[1][0],rmatn[1][1],rmatn[1][2],
    rmatn[2][0],rmatn[2][1],rmatn[2][2]);    
}

PyDoc_STRVAR(erfa_num06a_doc,
"num06a(d1, d2) -> rmatn\n\n"
"Form the matrix of nutation for a given date, IAU 2006/2000A model.\n"
"Given:\n"
"     d1,d2     TT as a 2-part Julian Date\n"
"Returned:\n"
"     rmatn     nutation matrix");

static PyObject *
erfa_numat(PyObject *self, PyObject *args)
{
    double epsa, dpsi, deps, rmatn[3][3];
    if (!PyArg_ParseTuple(args, "ddd", &epsa, &dpsi, &deps)) {
        return NULL;
    }
    eraNumat(epsa, dpsi, deps, rmatn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rmatn[0][0],rmatn[0][1],rmatn[0][2],
    rmatn[1][0],rmatn[1][1],rmatn[1][2],
    rmatn[2][0],rmatn[2][1],rmatn[2][2]);    
}

PyDoc_STRVAR(erfa_numat_doc,
"numat(epsa, dpsi, deps) -> rmatn\n\n"
"Form the matrix of nutation.\n"
"Given:\n"
"     epsa          mean obliquity of date\n"
"     dpsi,deps     nutation\n"
"Returned:\n"
"     rmatn         nutation matrix");

static PyObject *
erfa_nut00a(PyObject *self, PyObject *args)
{
    double d1, d2, dpsi, deps;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNut00a(d1, d2, &dpsi, &deps);
    return Py_BuildValue("dd", dpsi, deps);
}

PyDoc_STRVAR(erfa_nut00a_doc,
"nut00a(d1, d2) -> dpsi, deps\n\n"
"Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation\n"
"with free core nutation omitted).\n"
"Given:\n"
"    d1,d2      TT as a 2-part Julian Date\n"
"Returned:\n"
"   dpsi,deps   nutation, luni-solar + planetary\n");

static PyObject *
erfa_nut00b(PyObject *self, PyObject *args)
{
    double d1, d2, dpsi, deps;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNut00b(d1, d2, &dpsi, &deps);
    return Py_BuildValue("dd", dpsi, deps);
}

PyDoc_STRVAR(erfa_nut00b_doc,
"nut00b(d1, d2) -> dpsi, deps\n\n"
"Nutation, IAU 2000B.\n"
"Given:\n"
"    d1,d2       TT as a 2-part Julian Date\n"
"Returned:\n"
"   dpsi,deps   nutation, luni-solar + planetary\n");

static PyObject *
erfa_nut06a(PyObject *self, PyObject *args)
{
    double d1, d2, dpsi, deps;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNut06a(d1, d2, &dpsi, &deps);
    return Py_BuildValue("dd", dpsi, deps);
}

PyDoc_STRVAR(erfa_nut06a_doc,
"nut06a(d1, d2) -> dpsi, deps\n\n"
"IAU 2000A nutation with adjustments to match the IAU 2006 precession.\n"
"Given:\n"
"    d1,d2      TT as a 2-part Julian Date\n"
"Returned:\n"
"   dpsi,deps   nutation, luni-solar + planetary\n");

static PyObject *
erfa_nut80(PyObject *self, PyObject *args)
{
    double d1, d2, dpsi, deps;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNut80(d1, d2, &dpsi, &deps);
    return Py_BuildValue("dd", dpsi, deps);
}

PyDoc_STRVAR(erfa_nut80_doc,
"nut80(d1, d2) -> dpsi, deps\n\n"
"Nutation, IAU 1980 model.\n"
"Given:\n"
"    d1,d2      TT as a 2-part Julian Date\n"
"Returned:\n"
"    dpsi        nutation in longitude (radians)\n"
"    deps        nutation in obliquity (radians)\n");

static PyObject *
erfa_nutm80(PyObject *self, PyObject *args)
{
    double d1, d2, rmatn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraNutm80(d1, d2, rmatn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rmatn[0][0],rmatn[0][1],rmatn[0][2],
    rmatn[1][0],rmatn[1][1],rmatn[1][2],
    rmatn[2][0],rmatn[2][1],rmatn[2][2]);    
}

PyDoc_STRVAR(erfa_nutm80_doc,
"nutm80(d1, d2) -> rmatn\n\n"
"Form the matrix of nutation for a given date, IAU 1980 model.\n"
"Given:\n"
"    d1,d2      TDB as a 2-part Julian Date\n"
"Returned:\n"
"    rmatn[     nutation matrix");

static PyObject *
erfa_obl06(PyObject *self, PyObject *args)
{
    double d1, d2, obl;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    obl = eraObl06(d1, d2);
        return Py_BuildValue("d", obl);
}

PyDoc_STRVAR(erfa_obl06_doc,
"obl06(d1, d2) -> obl\n\n"
"Mean obliquity of the ecliptic, IAU 2006 precession model.\n"
"Given:\n"
"    d1,d2      TT as a 2-part Julian Date\n"
"Returned:\n"
"    obl        obliquity of the ecliptic (radians)");

static PyObject *
erfa_obl80(PyObject *self, PyObject *args)
{
    double d1, d2, obl;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    obl = eraObl80(d1, d2);
    return Py_BuildValue("d", obl);
}

PyDoc_STRVAR(erfa_obl80_doc,
"obl80(d1, d2) -> obl\n\n"
"Mean obliquity of the ecliptic, IAU 1980 model.\n"
"Given:\n"
"    d1,d2      TT as a 2-part Julian Date\n"
"Returned:\n"
"    obl        obliquity of the ecliptic (radians)");

static PyObject *
erfa_p06e(PyObject *self, PyObject *args)
{
    double d1,d2, eps0,psia,oma,bpa,bqa,pia,bpia,epsa,chia,za,zetaa,thetaa,pa,gam,phi,psi;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraP06e(d1,d2,
    &eps0,&psia,&oma,&bpa,&bqa,&pia,&bpia,&epsa,
    &chia,&za,&zetaa,&thetaa,&pa,&gam,&phi,&psi);
    return Py_BuildValue("dddddddddddddddd",
    eps0,psia,oma,bpa,bqa,pia,bpia,epsa,chia,za,zetaa,thetaa,pa,gam,phi,psi);
}

PyDoc_STRVAR(erfa_p06e_doc,
"p06e(d1, d2) -> eps0,psia,oma,bpa,bqa,pia,bpia,epsa,chia,za,zetaa,thetaa,pa,gam,phi,psi\n\n"
"Precession angles, IAU 2006, equinox based.\n"
"Given:\n"
"   d1,d2   TT as a 2-part Julian Date\n"
"Returned:\n"
"   eps0   epsilon_0   obliquity at J2000.0\n"
"   psia   psi_A       luni-solar precession\n"
"   oma    omega_A     inclination of equator wrt J2000.0 ecliptic\n"
"   bpa    P_A         ecliptic pole x, J2000.0 ecliptic triad\n"
"   bqa    Q_A         ecliptic pole -y, J2000.0 ecliptic triad\n"
"   pia    pi_A        angle between moving and J2000.0 ecliptics\n"
"   bpia   Pi_A        longitude of ascending node of the ecliptic\n"
"   epsa   epsilon_A   obliquity of the ecliptic\n"
"   chia   chi_A       planetary precession\n"
"   za     z_A         equatorial precession: -3rd 323 Euler angle\n"
"   zetaa  zeta_A      equatorial precession: -1st 323 Euler angle\n"
"   thetaa theta_A     equatorial precession: 2nd 323 Euler angle\n"
"   pa     p_A         general precession\n"
"   gam    gamma_J2000 J2000.0 RA difference of ecliptic poles\n"
"   phi    phi_J2000   J2000.0 codeclination of ecliptic pole\n"
"   psi    psi_J2000   longitude difference of equator poles, J2000.0");

static PyObject *
erfa_pb06(PyObject *self, PyObject *args)
{
    double d1, d2, bzeta, bz, btheta;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPb06(d1, d2, &bzeta, &bz, &btheta);
    return Py_BuildValue("ddd", bzeta, bz, btheta);
}

PyDoc_STRVAR(erfa_pb06_doc,
"pb06(d1, d2) -> bzeta, bz, btheta\n\n"
"Forms three Euler angles which implement general\n"
"precession from epoch J2000.0, using the IAU 2006 model.\n"
"Framebias (the offset between ICRS and mean J2000.0) is included.\n"
"Given:\n"
"   d1,d2   TT as a 2-part Julian Date\n"
"Returned:\n"
"   bzeta   1st rotation: radians cw around z\n"
"   bz      3rd rotation: radians cw around z\n"
"   btheta  2nd rotation: radians ccw around y");

static PyObject *
erfa_pfw06(PyObject *self, PyObject *args)
{
    double d1, d2, gamb, phib, psib, epsa;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPfw06(d1,d2, &gamb, &phib, &psib, &epsa);
    return Py_BuildValue("dddd", gamb, phib, psib, epsa);
}

PyDoc_STRVAR(erfa_pfw06_doc,
"pfw06(d1, d2) -> gamb, phib, psib, epsa\n\n"
"Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).\n"
"Given:\n"
"   d1,d2   TT as a 2-part Julian Date\n"
"Returned:\n"
"   gamb    F-W angle gamma_bar (radians)\n"
"   phib    F-W angle phi_bar (radians)\n"
"   psib    F-W angle psi_bar (radians)\n"
"   epsa    F-W angle epsilon_A (radians)");

static PyObject *
erfa_plan94(PyObject *self, PyObject *args)
{
    double d1, d2, pv[2][3];
    int np, status;
    if (!PyArg_ParseTuple(args, "ddi", &d1, &d2, &np)) {
        return NULL;
    }
    status = eraPlan94(d1, d2, np, pv);
    if (status == -1){
        PyErr_SetString(erfaError, "illegal np,  not in range(1,8) for planet");
        return NULL;
    }
/*    else if (status == 1){
        PyErr_SetString(erfaError, "year outside range(1000:3000)");
        return NULL;
    }*/
    else if (status == 2){
        PyErr_SetString(erfaError, "computation failed to converge");
        return NULL;
    }
    else {
        return Py_BuildValue("(ddd)(ddd)",
        pv[0][0], pv[0][1], pv[0][2],
        pv[1][0], pv[1][1], pv[1][2]);
    }
}

PyDoc_STRVAR(erfa_plan94_doc,
"plan94(d1, d2, np) -> pv\n\n"
"Approximate heliocentric position and velocity of a nominated major\n"
"planet:  Mercury, Venus, EMB, Mars, Jupiter, Saturn, Uranus or\n"
"Neptune (but not the Earth itself).\n"
"Given:\n"
"    d1         TDB date part A\n"
"    d2         TDB date part B\n"
"    np         planet (1=Mercury, 2=Venus, 3=EMB, 4=Mars,\n"
"                       5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune)\n"
"Returned:\n"
"    pv         planet p,v (heliocentric, J2000.0, AU,AU/d)");

static PyObject *
erfa_pmat00(PyObject *self, PyObject *args)
{
    double d1, d2, rbp[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPmat00(d1, d2, rbp);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rbp[0][0],rbp[0][1],rbp[0][2],
    rbp[1][0],rbp[1][1],rbp[1][2],
    rbp[2][0],rbp[2][1],rbp[2][2]);    
}

PyDoc_STRVAR(erfa_pmat00_doc,
"pmat00(d1, d2) -> rbp\n\n"
"Precession matrix (including frame bias) from GCRS to a specified\n"
"date, IAU 2000 model.\n"
"Given:\n"
"   d1,d2       TT as a 2-part Julian Date\n"
"Returned:\n"
"   rbp         bias-precession matrix");

static PyObject *
erfa_pmat06(PyObject *self, PyObject *args)
{
    double d1, d2, rbp[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPmat06(d1, d2, rbp);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rbp[0][0],rbp[0][1],rbp[0][2],
    rbp[1][0],rbp[1][1],rbp[1][2],
    rbp[2][0],rbp[2][1],rbp[2][2]);    
}

PyDoc_STRVAR(erfa_pmat06_doc,
"pmat06(d1, d2) -> rbp\n\n"
"Precession matrix (including frame bias) from GCRS to a specified\n"
"date, IAU 2006 model.\n"
"Given:\n"
"   d1,d2       TT as a 2-part Julian Date\n"
"Returned:\n"
"   rbp         bias-precession matrix");

static PyObject *
erfa_pmat76(PyObject *self, PyObject *args)
{
    double d1, d2, rmatp[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPmat76(d1, d2, rmatp);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rmatp[0][0],rmatp[0][1],rmatp[0][2],
    rmatp[1][0],rmatp[1][1],rmatp[1][2],
    rmatp[2][0],rmatp[2][1],rmatp[2][2]);    
}

PyDoc_STRVAR(erfa_pmat76_doc,
"pmat76(d1, d2) -> rmatp\n\n"
"Precession matrix from J2000.0 to a specified date, IAU 1976 model.\n"
"Given:\n"
"   d1,d2       TT ending date as a 2-part Julian Date\n"
"Returned:\n"
"   rmatp       precession matrix, J2000.0 -> d1+d2");

static PyObject *
erfa_pn00(PyObject *self, PyObject *args)
{
    double d1,d2,dpsi,deps,epsa;
    double rb[3][3],rp[3][3],rbp[3][3],rn[3][3],rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dddd", &d1, &d2, &dpsi, &deps)) {
        return NULL;
    }
    eraPn00(d1,d2,dpsi,deps,&epsa,rb,rp,rbp,rn,rbpn);
    return Py_BuildValue("d((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))",
    epsa,
    rb[0][0],rb[0][1],rb[0][2],rb[1][0],rb[1][1],rb[1][2],rb[2][0],rb[2][1],rb[2][2],
    rp[0][0],rp[0][1],rp[0][2],rp[1][0],rp[1][1],rp[1][2],rp[2][0],rp[2][1],rp[2][2],
    rbp[0][0],rbp[0][1],rbp[0][2],rbp[1][0],rbp[1][1],rbp[1][2],rbp[2][0],rbp[2][1],rbp[2][2],
    rn[0][0],rn[0][1],rn[0][2],rn[1][0],rn[1][1],rn[1][2],rn[2][0],rn[2][1],rn[2][2],
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}

PyDoc_STRVAR(erfa_pn00_doc,
"pn00(d1,d2,dpsi,deps) -> epsa,rb,rp,rbp,rn,rbpn\n\n"
"Precession-nutation, IAU 2000 model:  a multi-purpose function,\n"
"supporting classical (equinox-based) use directly and CIO-based\n"
"use indirectly.\n"
"Given:\n"
"   d1,d2       TT as a 2-part Julian Date\n"
"   dpsi,deps   nutation\n"
"Returned:\n"
"   epsa        mean obliquity\n"
"   rb          frame bias matrix\n"
"   rp          precession matrix\n"
"   rbp         bias-precession matrix\n"
"   rn          nutation matrix\n"
"   rbpn        GCRS-to-true matrix");

static PyObject *
erfa_pn00a(PyObject *self, PyObject *args)
{
    double d1,d2,dpsi,deps,epsa;
    double rb[3][3],rp[3][3],rbp[3][3],rn[3][3],rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPn00a(d1,d2,&dpsi,&deps,&epsa,rb,rp,rbp,rn,rbpn);
    return Py_BuildValue("ddd((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))",
    dpsi,deps,epsa,
    rb[0][0],rb[0][1],rb[0][2],rb[1][0],rb[1][1],rb[1][2],rb[2][0],rb[2][1],rb[2][2],
    rp[0][0],rp[0][1],rp[0][2],rp[1][0],rp[1][1],rp[1][2],rp[2][0],rp[2][1],rp[2][2],
    rbp[0][0],rbp[0][1],rbp[0][2],rbp[1][0],rbp[1][1],rbp[1][2],rbp[2][0],rbp[2][1],rbp[2][2],
    rn[0][0],rn[0][1],rn[0][2],rn[1][0],rn[1][1],rn[1][2],rn[2][0],rn[2][1],rn[2][2],
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}

PyDoc_STRVAR(erfa_pn00a_doc,
"pn00a(d1,d2) -> dpsi,deps,epsa,rb,rp,rbp,rn,rbpn\n\n"
"Precession-nutation, IAU 2000A model:  a multi-purpose function,\n"
"supporting classical (equinox-based) use directly and CIO-based\n"
"use indirectly.\n"
"Given:\n"
"   d1,d2       TT as a 2-part Julian Date\n"
"Returned:\n"
"   dpsi,deps   nutation\n"
"   epsa        mean obliquity\n"
"   rb          frame bias matrix\n"
"   rp          precession matrix\n"
"   rbp         bias-precession matrix\n"
"   rn          nutation matrix\n"
"   rbpn        GCRS-to-true matrix");

static PyObject *
erfa_pn00b(PyObject *self, PyObject *args)
{
    double d1,d2,dpsi,deps,epsa;
    double rb[3][3],rp[3][3],rbp[3][3],rn[3][3],rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPn00b(d1,d2,&dpsi,&deps,&epsa,rb,rp,rbp,rn,rbpn);
    return Py_BuildValue("ddd((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))",
    dpsi,deps,epsa,
    rb[0][0],rb[0][1],rb[0][2],rb[1][0],rb[1][1],rb[1][2],rb[2][0],rb[2][1],rb[2][2],
    rp[0][0],rp[0][1],rp[0][2],rp[1][0],rp[1][1],rp[1][2],rp[2][0],rp[2][1],rp[2][2],
    rbp[0][0],rbp[0][1],rbp[0][2],rbp[1][0],rbp[1][1],rbp[1][2],rbp[2][0],rbp[2][1],rbp[2][2],
    rn[0][0],rn[0][1],rn[0][2],rn[1][0],rn[1][1],rn[1][2],rn[2][0],rn[2][1],rn[2][2],
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}

PyDoc_STRVAR(erfa_pn00b_doc,
"pn00b(d1,d2) -> dpsi,deps,epsa,rb,rp,rbp,rn,rbpn\n\n"
"Precession-nutation, IAU 2000B model:  a multi-purpose function,\n"
"supporting classical (equinox-based) use directly and CIO-based\n"
"use indirectly.\n"
"Given:\n"
"   d1,d2       TT as a 2-part Julian Date\n"
"Returned:\n"
"   dpsi,deps   nutation\n"
"   epsa        mean obliquity\n"
"   rb          frame bias matrix\n"
"   rp          precession matrix\n"
"   rbp         bias-precession matrix\n"
"   rn[         nutation matrix\n"
"   rbpn        GCRS-to-true matrix");

static PyObject *
erfa_pn06(PyObject *self, PyObject *args)
{
    double d1,d2,dpsi,deps,epsa;
    double rb[3][3],rp[3][3],rbp[3][3],rn[3][3],rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dddd", &d1, &d2, &dpsi, &deps)) {
        return NULL;
    }
    eraPn06(d1,d2,dpsi,deps,&epsa,rb,rp,rbp,rn,rbpn);
    return Py_BuildValue("d((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))",
    epsa,
    rb[0][0],rb[0][1],rb[0][2],rb[1][0],rb[1][1],rb[1][2],rb[2][0],rb[2][1],rb[2][2],
    rp[0][0],rp[0][1],rp[0][2],rp[1][0],rp[1][1],rp[1][2],rp[2][0],rp[2][1],rp[2][2],
    rbp[0][0],rbp[0][1],rbp[0][2],rbp[1][0],rbp[1][1],rbp[1][2],rbp[2][0],rbp[2][1],rbp[2][2],
    rn[0][0],rn[0][1],rn[0][2],rn[1][0],rn[1][1],rn[1][2],rn[2][0],rn[2][1],rn[2][2],
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}

PyDoc_STRVAR(erfa_pn06_doc,
"pn06(d1,d2,dpsi,deps) -> epsa,rb,rp,rbp,rn,rbpn\n\n"
"Precession-nutation, IAU 2006 model:  a multi-purpose function,\n"
"supporting classical (equinox-based) use directly and CIO-based\n"
"use indirectly.\n"
"Given:\n"
"   d1,d2       TT as a 2-part Julian Date\n"
"   dpsi,deps   nutation\n"
"Returned:\n"
"   epsa        mean obliquity\n"
"   rb          frame bias matrix\n"
"   rp          precession matrix\n"
"   rbp         bias-precession matrix\n"
"   rn          nutation matrix\n"
"   rbpn        GCRS-to-true matrix");

static PyObject *
erfa_pn06a(PyObject *self, PyObject *args)
{
    double d1,d2,dpsi,deps,epsa;
    double rb[3][3],rp[3][3],rbp[3][3],rn[3][3],rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPn06a(d1,d2,&dpsi,&deps,&epsa,rb,rp,rbp,rn,rbpn);
    return Py_BuildValue("ddd((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))",
    dpsi,deps,epsa,
    rb[0][0],rb[0][1],rb[0][2],rb[1][0],rb[1][1],rb[1][2],rb[2][0],rb[2][1],rb[2][2],
    rp[0][0],rp[0][1],rp[0][2],rp[1][0],rp[1][1],rp[1][2],rp[2][0],rp[2][1],rp[2][2],
    rbp[0][0],rbp[0][1],rbp[0][2],rbp[1][0],rbp[1][1],rbp[1][2],rbp[2][0],rbp[2][1],rbp[2][2],
    rn[0][0],rn[0][1],rn[0][2],rn[1][0],rn[1][1],rn[1][2],rn[2][0],rn[2][1],rn[2][2],
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}

PyDoc_STRVAR(erfa_pn06a_doc,
"pn06a(d1,d2) -> dpsi,deps,epsa,rb,rp,rbp,rn,rbpn\n\n"
"Precession-nutation, IAU 2006/2000A model:  a multi-purpose function,\n"
"supporting classical (equinox-based) use directly and CIO-based\n"
"use indirectly.\n"
"Given:\n"
"   d1,d2       TT as a 2-part Julian Date\n"
"Returned:\n"
"   dpsi,deps   nutation\n"
"   epsa        mean obliquity\n"
"   rb          frame bias matrix\n"
"   rp          precession matrix\n"
"   rbp         bias-precession matrix\n"
"   rn          nutation matrix\n"
"   rbpn        GCRS-to-true matrix");

static PyObject *
erfa_pnm00a(PyObject *self, PyObject *args)
{
    double d1,d2,rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPnm00a(d1, d2, rbpn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}

PyDoc_STRVAR(erfa_pnm00a_doc,
"pnm00a(d1, d2) -> rbpn\n\n"
"Form the matrix of precession-nutation for a given date (including\n"
"frame bias), equinox-based, IAU 2000A model.\n"
"Given:\n"
"   d1,d2       TT as a 2-part Julian Date\n"
"Returned:\n"
"   rbpn        classical NPB matrix");

static PyObject *
erfa_pnm00b(PyObject *self, PyObject *args)
{
    double d1,d2,rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPnm00b(d1, d2, rbpn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}

PyDoc_STRVAR(erfa_pnm00b_doc,
"pnm00b(d1, d2) -> rbpn\n\n"
"Form the matrix of precession-nutation for a given date (including\n"
"frame bias), equinox-based, IAU 2000B model.\n"
"Given:\n"
"   d1,d2       TT as a 2-part Julian Date\n"
"Returned:\n"
"   rbpn        classical NPB matrix");

static PyObject *
erfa_pnm06a(PyObject *self, PyObject *args)
{
    double d1,d2,rbpn[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPnm06a(d1, d2, rbpn);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rbpn[0][0],rbpn[0][1],rbpn[0][2],rbpn[1][0],rbpn[1][1],rbpn[1][2],rbpn[2][0],rbpn[2][1],rbpn[2][2]);
}

PyDoc_STRVAR(erfa_pnm06a_doc,
"pnm06a(d1, d2) -> rbpn\n\n"
"Form the matrix of precession-nutation for a given date (including\n"
"frame bias), IAU 2006 precession and IAU 2000A nutation models.\n"
"Given:\n"
"   d1,d2       TT as a 2-part Julian Date\n"
"Returned:\n"
"   rbpn        classical NPB matrix");

static PyObject *
erfa_pnm80(PyObject *self, PyObject *args)
{
    double d1, d2, rmatp[3][3];
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPnm80(d1, d2, rmatp);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rmatp[0][0],rmatp[0][1],rmatp[0][2],
    rmatp[1][0],rmatp[1][1],rmatp[1][2],
    rmatp[2][0],rmatp[2][1],rmatp[2][2]);    
}

PyDoc_STRVAR(erfa_pnm80_doc,
"pnm80(d1, d2) -> rmatp\n\n"
"Form the matrix of precession/nutation for a given date, IAU 1976\n"
"precession model, IAU 1980 nutation model.\n"
"Given:\n"
"   d1,d2           TDB date as a 2-part Julian Date\n"
"Returned:\n"
"   rmatp           combined precession/nutation matrix");

static PyObject *
erfa_pom00(PyObject *self, PyObject *args)
{    
    double xp, yp, sp, rpom[3][3];
    if (!PyArg_ParseTuple(args, "ddd", &xp, &yp, &sp)) {
        return NULL;
    }
    eraPom00(xp, yp, sp, rpom);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    rpom[0][0],rpom[0][1],rpom[0][2],
    rpom[1][0],rpom[1][1],rpom[1][2],
    rpom[2][0],rpom[2][1],rpom[2][2]);    
}

PyDoc_STRVAR(erfa_pom00_doc,
"pom00(xp, yp, sp) -> rpom\n\n"
"Form the matrix of polar motion for a given date, IAU 2000.\n"
"Given:\n"
"   xp,yp       coordinates of the pole (radians)\n"
"   sp          the TIO locator s' (radians)\n"
"Returned:\n"
"   rpom        polar-motion matrix");

static PyObject *
erfa_pr00(PyObject *self, PyObject *args)
{
    double d1, d2, dpsipr, depspr;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraPr00(d1, d2, &dpsipr, &depspr);
    return Py_BuildValue("dd", dpsipr, depspr);
}

PyDoc_STRVAR(erfa_pr00_doc,
"pr00(d1,d2) -> dpsipr,depspr\n\n"
"Precession-rate part of the IAU 2000 precession-nutation models\n"
"(part of MHB2000).\n"
"Given:\n"
"   d1,d2           TT as a 2-part Julian Date\n"
"Returned:\n"
"   dpsipr,depspr   precession corrections");

static PyObject *
erfa_prec76(PyObject *self, PyObject *args)
{
    double ep01, ep02, ep11, ep12, zeta, z, theta;
    if (!PyArg_ParseTuple(args, "dddd", &ep01, &ep02, &ep11, &ep12)) {
        return NULL;
    }
    eraPrec76(ep01, ep02, ep11, ep12, &zeta, &z, &theta);
    return Py_BuildValue("ddd", zeta, z, theta);
}

PyDoc_STRVAR(erfa_prec76_doc,
"prec76(ep01, ep02, ep11, ep12) -> zeta, z, theta\n\n"
"IAU 1976 precession model.\n"
"Given:\n"
"   ep01,ep02   TDB starting epoch as 2-part Julian Date\n"
"   ep11,ep12   TDB ending epoch as 2-part Julian Date\n"
"Returned:\n"
"   zeta        1st rotation: radians cw around z\n"
"   z           3rd rotation: radians cw around z\n"
"   theta       2nd rotation: radians ccw around y");

static PyObject *
erfa_pvstar(PyObject *self, PyObject *args)
{
    double pv[2][3], ra, dec, pmr, pmd, px, rv;
    double pv00, pv01, pv02, pv10, pv11, pv12;
    int status;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))",
        &pv00, &pv01, &pv02, &pv10, &pv11, &pv12)) {
        return NULL;
    }
    pv[0][0] = pv00;
    pv[0][1] = pv01;
    pv[0][2] = pv02;
    pv[1][0] = pv10;
    pv[1][1] = pv11;
    pv[1][2] = pv12;
    status = eraPvstar(pv, &ra, &dec, &pmr, &pmd, &px, &rv);
    if (status == -1) {
        PyErr_SetString(erfaError, "superluminal speed");
        return NULL;
    }
    else if (status == -2) {
        PyErr_SetString(erfaError, "null position vector");
        return NULL;
    }
    else {
    return Py_BuildValue("dddddd", ra, dec, pmr, pmd, px, rv);
    }
}
PyDoc_STRVAR(erfa_pvstar_doc,
"pvstar(pv[2][3]) -> ra, dec, pmr, pmd, px, rv\n\n"
"Convert star position+velocity vector to catalog coordinates.\n"
"Given:\n"
"   pv[2][3]    pv-vector (AU, AU/day)\n"
"Returned:\n"
"   ra          right ascension (radians)\n"
"   dec         declination (radians)\n"
"   pmr         RA proper motion (radians/year)\n"
"   pmd         Dec proper motion (radians/year)\n"
"   px          parallax (arcsec)\n"
"   rv          radial velocity (km/s, positive = receding)");

static PyObject *
erfa_s00(PyObject *self, PyObject *args)
{
    double d1, d2, x, y, s;
    if (!PyArg_ParseTuple(args, "dddd", &d1, &d2, &x, &y)) {
        return NULL;
    }
    s = eraS00(d1, d2, x, y);
    return Py_BuildValue("d", s);
}

PyDoc_STRVAR(erfa_s00_doc,
"s00(d1, d2, x, y) -> s\n\n"
"The CIO locator s, positioning the Celestial Intermediate Origin on\n"
"the equator of the Celestial Intermediate Pole, given the CIP's X,Y\n"
"coordinates.  Compatible with IAU 2000A precession-nutation.\n"
"Given:\n"
"   d1,d2   TT as a 2-part Julian Date\n"
"   x,y     CIP coordinates\n"
"Returned:\n"
"   s       the CIO locator s in radians");

static PyObject *
erfa_s00a(PyObject *self, PyObject *args)
{
    double d1, d2,s;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    s = eraS00a(d1, d2);
    return Py_BuildValue("d", s);
}

PyDoc_STRVAR(erfa_s00a_doc,
"s00a(d1, d2) -> s\n\n"
"The CIO locator s, positioning the Celestial Intermediate Origin on\n"
"the equator of the Celestial Intermediate Pole, using the IAU 2000A\n"
"precession-nutation model.\n"
"Given:\n"
"   d1,d2   TT as a 2-part Julian Date\n"
"Returned:\n"
"   s       the CIO locator s in radians");

static PyObject *
erfa_s00b(PyObject *self, PyObject *args)
{
    double d1, d2,s;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    s = eraS00b(d1, d2);
    return Py_BuildValue("d", s);
}

PyDoc_STRVAR(erfa_s00b_doc,
"s00b(d1, d2) -> s\n\n"
"The CIO locator s, positioning the Celestial Intermediate Origin on\n"
"the equator of the Celestial Intermediate Pole, using the IAU 2000B\n"
"precession-nutation model.\n"
"Given:\n"
"   d1,d2   TT as a 2-part Julian Date\n"
"Returned:\n"
"   s       the CIO locator s in radians");

static PyObject *
erfa_s06(PyObject *self, PyObject *args)
{
    double d1, d2, x, y, s;
    if (!PyArg_ParseTuple(args, "dddd", &d1, &d2, &x, &y)) {
        return NULL;
    }
    s = eraS06(d1, d2, x, y);
    return Py_BuildValue("d", s);
}

PyDoc_STRVAR(erfa_s06_doc,
"s06(d1, d2, x, y) -> s\n\n"
"The CIO locator s, positioning the Celestial Intermediate Origin on\n"
"the equator of the Celestial Intermediate Pole, given the CIP's X,Y\n"
"coordinates.  Compatible with IAU 2006/2000A precession-nutation.\n"
"Given:\n"
"   d1,d2   TT as a 2-part Julian Date\n"
"   x,y     CIP coordinates\n"
"Returned:\n"
"   s       the CIO locator s in radians");

static PyObject *
erfa_s06a(PyObject *self, PyObject *args)
{
    double d1, d2,s;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    s = eraS06a(d1, d2);
    return Py_BuildValue("d", s);
}

PyDoc_STRVAR(erfa_s06a_doc,
"s06a(d1, d2) -> s\n\n"
"The CIO locator s, positioning the Celestial Intermediate Origin on\n"
"the equator of the Celestial Intermediate Pole, using the IAU 2006\n"
"precession and IAU 2000A nutation model.\n"
"Given:\n"
"   d1,d2   TT as a 2-part Julian Date\n"
"Returned:\n"
"   s       the CIO locator s in radians");

static PyObject *
erfa_sp00(PyObject *self, PyObject *args)
{
    double d1, d2, s;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    s = eraSp00(d1, d2);  
    return Py_BuildValue("d", s);
}

PyDoc_STRVAR(erfa_sp00_doc,
"sp00(d1, d2) -> s\n\n"
"The TIO locator s', positioning the Terrestrial Intermediate Origin\n"
"on the equator of the Celestial Intermediate Pole.\n"
"Given:\n"
"   d1,d2   TT as a 2-part Julian Date\n"
"Returned:\n"
"   s       the TIO locator s' in radians");

static PyObject *
erfa_starpm(PyObject *self, PyObject *args)
{
   double ra1, dec1, pmr1, pmd1, px1, rv1;
   double ep1a, ep1b, ep2a, ep2b;
   double ra2, dec2, pmr2, pmd2, px2, rv2;
   int status;
   if (!PyArg_ParseTuple(args, "dddddddddd",
        &ra1, &dec1, &pmr1, &pmd1, &px1, &rv1,
        &ep1a, &ep1b, &ep2a, &ep2b)) {
        return NULL;
    }
    status = eraStarpm(ra1, dec1, pmr1, pmd1, px1, rv1,
                 ep1a, ep1b, ep2a, ep2b,
                 &ra2, &dec2, &pmr2, &pmd2, &px2, &rv2);
    if (status) {
        if (status == -1) {
            PyErr_SetString(erfaError, "system error, normally should NOT occur");
            return NULL;
        }
        else if (status == 1) {
            PyErr_SetString(erfaError, "distance overriden, extremely small (or zero or negative) parallax");
            return NULL;
        }
        else if (status == 2) {
            PyErr_SetString(erfaError, "excessive velocity");
            return NULL;
        }
        else if (status == 4) {
            PyErr_SetString(erfaError, "solution didn't converge");
            return NULL;
        }
        else {
            PyErr_SetString(erfaError, "binary logical OR of other error above");
            return NULL;
        }
    }
    else {
        return Py_BuildValue("dddddd", ra2, dec2, pmr2, pmd2, px2, rv2);
    }
}

PyDoc_STRVAR(erfa_starpm_doc,
"starpm(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b) -> ra2, dec2, pmr2, pmd2, px2, rv2\n\n"
"Star proper motion:  update star catalog data for space motion.\n"
"Given:\n"
"   ra1     right ascension (radians), before\n"
"   dec1    declination (radians), before\n"
"   pmr1    RA proper motion (radians/year), before\n"
"   pmd1    Dec proper motion (radians/year), before\n"
"   px1     parallax (arcseconds), before\n"
"   rv1     radial velocity (km/s, +ve = receding), before\n"
"   ep1a    ''before'' epoch, part A \n"
"   ep1b    ''before'' epoch, part B\n"
"   ep2a    ''after'' epoch, part A\n"
"   ep2b    ''after'' epoch, part B\n"
"Returned:\n"
"   ra2     right ascension (radians), after\n"
"   dec2    declination (radians), after\n"
"   pmr2    RA proper motion (radians/year), after\n"
"   pmd2    Dec proper motion (radians/year), after\n"
"   px2     parallax (arcseconds), after\n"
"   rv2     radial velocity (km/s, +ve = receding), after");

static PyObject *
erfa_starpv(PyObject *self, PyObject *args)
{
   double ra, dec, pmr, pmd, px, rv;
   double pv[2][3];
   int status;
   if (!PyArg_ParseTuple(args, "dddddd",
        &ra, &dec, &pmr, &pmd, &px, &rv)) {
        return NULL;
    }
    status = eraStarpv(ra, dec, pmr, pmd, px, rv, pv);
    if (status) {
        if (status == 1) {
            PyErr_SetString(erfaError, "distance overriden, extremely small (or zero or negative) parallax");
            return NULL;
        }
        else if (status == 2) {
            PyErr_SetString(erfaError, "excessive velocity");
            return NULL;
        }
        else if (status == 4) {
            PyErr_SetString(erfaError, "solution didn't converge");
            return NULL;
        }
        else {
            PyErr_SetString(erfaError, "binary logical OR of other error above");
            return NULL;
        }
    }
    else {
        return Py_BuildValue("(ddd)(ddd)",
            pv[0][0],pv[0][1],pv[0][2],
            pv[1][0],pv[1][1],pv[1][2]);
    }
}

PyDoc_STRVAR(erfa_starpv_doc,
"starpv(ra1, dec1, pmr1, pmd1, px1, rv1) -> pv\n\n"
"Convert star catalog coordinates to position+velocity vector.\n"
"Given:\n"
"   ra      right ascension (radians)\n"
"   dec     declination (radians)\n"
"   pmr     RA proper motion (radians/year)\n"
"   pmd     Dec proper motion (radians/year)\n"
"   px      parallax (arcseconds)\n"
"   rv  radial velocity (km/s, positive = receding)\n"
"Returned:\n"
"   pv      pv-vector (AU, AU/day)");

static PyObject *
erfa_taitt(PyObject *self, PyObject *args)
{
    double tai1, tai2, tt1, tt2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &tai1, &tai2)) {
        return NULL;
    }
    status = eraTaitt(tai1, tai2, &tt1, &tt2);
    if (status) {
        PyErr_SetString(erfaError, "internal error...");
        return NULL;
    }
    return Py_BuildValue("dd", tt1, tt2);
}

PyDoc_STRVAR(erfa_taitt_doc,
"taitt(tai1, tai2) -> tt1, tt2\n\n"
"Time scale transformation:  International Atomic Time, TAI, to\n"
"Terrestrial Time, TT.\n"
"Given:\n"
"   tai1,tai2   TAI as a 2-part Julian Date\n"
"Returned:\n"
"   tt1,tt2     TT as a 2-part Julian Date");

static PyObject *
erfa_taiut1(PyObject *self, PyObject *args)
{
    double tai1, tai2, dta, ut11, ut12;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &tai1, &tai2, &dta)) {
        return NULL;
    }
    status = eraTaiut1(tai1, tai2, dta, &ut11, &ut12);
    if (status) {
        PyErr_SetString(erfaError, "internal error...");
        return NULL;
    }
    return Py_BuildValue("dd", ut11, ut12);
}

PyDoc_STRVAR(erfa_taiut1_doc,
"taiut1(tai1, tai2, dta) -> ut11, ut12\n\n"
"Time scale transformation:  International Atomic Time, TAI, to\n"
"Universal Time, UT1..\n"
"The argument dta, i.e. UT1-TAI, is an observed quantity, and is\n"
"available from IERS tabulations.\n"
"Given:\n"
"   tai1,tai2   TAI as a 2-part Julian Date\n"
"   dta         UT1-TAI in seconds\n"
"Returned:\n"
"   ut11,ut12     TT as a 2-part Julian Date");

static PyObject *
erfa_taiutc(PyObject *self, PyObject *args)
{
    double tai1, tai2, utc1, utc2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &tai1, &tai2)) {
        return NULL;
    }
    status = eraTaiutc(tai1, tai2, &utc1, &utc2);
    if (status) {
        if (status == 1) {
            PyErr_SetString(erfaError, "dubious year");
            return NULL;
        }
        else if (status == -1) {
            PyErr_SetString(erfaError, "unacceptable date");
            return NULL;
        }
    }
    return Py_BuildValue("dd", utc1, utc2);
}

PyDoc_STRVAR(erfa_taiutc_doc,
"taiutc(tai1, tai2) -> utc1, utc2\n\n"
"Time scale transformation:  International Atomic Time, TAI, to\n"
"Coordinated Universal Time, UTC.\n"
"Given:\n"
"   tai1,tai2   TAI as a 2-part Julian Date\n"
"Returned:\n"
"   utc1,utc2   TT as a 2-part Julian Date");

static PyObject *
erfa_tcbtdb(PyObject *self, PyObject *args)
{
    double tcb1, tcb2, tdb1, tdb2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &tcb1, &tcb2)) {
        return NULL;
    }
    status = eraTcbtdb(tcb1, tcb2, &tdb1, &tdb2);
    if (status) {
        PyErr_SetString(erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", tdb1, tdb2);
}

PyDoc_STRVAR(erfa_tcbtdb_doc,
"tcbtdb(tcb1, tcb2) -> tdb1, tdb2\n\n"
"Time scale transformation:  Barycentric Coordinate Time, TCB, to\n"
"Barycentric Dynamical Time, TDB.\n"
"Given:\n"
"   tcb1,tcb2   TCB as a 2-part Julian Date\n"
"Returned:\n"
"   tdb1,tdb2   TDB as a 2-part Julian Date");

static PyObject *
erfa_tcgtt(PyObject *self, PyObject *args)
{
    double tcg1, tcg2, tt1, tt2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &tcg1, &tcg2)) {
        return NULL;
    }
    status = eraTcgtt(tcg1, tcg2, &tt1, &tt2);
    if (status) {
        PyErr_SetString(erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", tt1, tt2);
}

PyDoc_STRVAR(erfa_tcgtt_doc,
"tcgtt(tcg1, tcg2) -> tt1, tt2\n\n"
"Time scale transformation:  Geocentric Coordinate Time, TCG, to\n"
"Terrestrial Time, TT.\n"
"Given:\n"
"   tcg1,tcg2   TCG as a 2-part Julian Date\n"
"Returned:\n"
"   tt1,tt2   TT as a 2-part Julian Date");

static PyObject *
erfa_tdbtcb(PyObject *self, PyObject *args)
{
    double tdb1, tdb2, tcb1, tcb2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &tdb1, &tdb2)) {
        return NULL;
    }
    status = eraTdbtcb(tdb1, tdb2, &tcb1, &tcb2);
    if (status) {
        PyErr_SetString(erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", tcb1, tcb2);
}

PyDoc_STRVAR(erfa_tdbtcb_doc,
"tdbtcb(tdb1, tdb2) -> tcb1, tcb2\n\n"
"Time scale transformation:  Barycentric Dynamical Time, TDB, to\n"
"Barycentric Coordinate Time, TCB.\n"
"Given:\n"
"   tdb1,tdb2   TDB as a 2-part Julian Date\n"
"Returned:\n"
"   tcb1,tcb2   TCB as a 2-part Julian Date");

static PyObject *
erfa_tdbtt(PyObject *self, PyObject *args)
{
    double tdb1, tdb2, dtr, tt1, tt2;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &tdb1, &tdb2, &dtr)) {
        return NULL;
    }
    status = eraTdbtt(tdb1, tdb2, dtr, &tt1, &tt2);
    if (status) {
        PyErr_SetString(erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", tt1, tt2);
}

PyDoc_STRVAR(erfa_tdbtt_doc,
"tdbtt(tdb1, tdb2, dtr) -> tt1, tt2\n\n"
"Time scale transformation: Barycentric Dynamical Time, TDB, to\n"
"Terrestrial Time, TT.\n"
"Given:\n"
"   tdb1,tdb2   TDB as a 2-part Julian Date\n"
"   dtr         TDB-TT in seconds\n"
"Returned:\n"
"   tt1,tt2   TT as a 2-part Julian Date");

static PyObject *
erfa_tttai(PyObject *self, PyObject *args)
{
    double tai1, tai2, tt1, tt2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &tt1, &tt2)) {
        return NULL;
    }
    status = eraTttai(tt1, tt2, &tai1, &tai2);
    if (status) {
        PyErr_SetString(erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", tai1, tai2);
}

PyDoc_STRVAR(erfa_tttai_doc,
"tttai(tt1, tt2) -> tai1, tai2\n\n"
"Time scale transformation: Terrestrial Time, TT, to\n"
"International Atomic Time, TAI.\n"
"Given:\n"
"   tt1,tt2     TT as a 2-part Julian Date\n"
"Returned:\n"
"   tai1,tai2   TAI as a 2-part Julian Date");

static PyObject *
erfa_tttcg(PyObject *self, PyObject *args)
{
    double tcg1, tcg2, tt1, tt2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &tt1, &tt2)) {
        return NULL;
    }
    status = eraTttcg(tt1, tt2, &tcg1, &tcg2);
    if (status) {
        PyErr_SetString(erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", tcg1, tcg2);
}

PyDoc_STRVAR(erfa_tttcg_doc,
"tttcg(tt1, tt2) -> tcg1, tcg2\n\n"
"Time scale transformation: Terrestrial Time, TT, to Geocentric\n"
"Coordinate Time, TCG.\n"
"Given:\n"
"   tt1,tt2     TT as a 2-part Julian Date\n"
"Returned:\n"
"   tcg1,tcg2   TCG as a 2-part Julian Date");

static PyObject *
erfa_tttdb(PyObject *self, PyObject *args)
{
    double tdb1, tdb2, dtr, tt1, tt2;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &tt1, &tt2, &dtr)) {
        return NULL;
    }
    status = eraTttdb(tt1, tt2, dtr, &tdb1, &tdb2);
    if (status) {
        PyErr_SetString(erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", tdb1, tdb2);
}

PyDoc_STRVAR(erfa_tttdb_doc,
"tttdb(tt1, tt2, dtr) -> tdb1, tdb2\n\n"
"Time scale transformation: Terrestrial Time, TT, to\n"
"Barycentric Dynamical Time, TDB.\n"
"Given:\n"
"   tt1,tt2     TT as a 2-part Julian Date\n"
"   dtr         TDB-TT in seconds\n"
"Returned:\n"
"   tdb1,tdb2   TDB as a 2-part Julian Date");

static PyObject *
erfa_ttut1(PyObject *self, PyObject *args)
{
    double ut11, ut12, dt, tt1, tt2;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &tt1, &tt2, &dt)) {
        return NULL;
    }
    status = eraTtut1(tt1, tt2, dt, &ut11, &ut12);
    if (status) {
        PyErr_SetString(erfaError, "should NOT occur");
        return NULL;
    }
    return Py_BuildValue("dd", ut11, ut12);
}

PyDoc_STRVAR(erfa_ttut1_doc,
"ttut1(tt1, tt2, dt) -> ut11, ut12\n\n"
"Time scale transformation: Terrestrial Time, TT, to\n"
"Universal Time UT1.\n"
"Given:\n"
"   tt1,tt2     TT as a 2-part Julian Date\n"
"   dt          TT-UT1 in seconds\n"
"Returned:\n"
"   ut11,ut12   UT1 as a 2-part Julian Date");

static PyObject *
erfa_ut1tai(PyObject *self, PyObject *args)
{
    double tai1, tai2, dta, ut11, ut12;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &ut11, &ut12, &dta)) {
        return NULL;
    }
    status = eraUt1tai(ut11, ut12, dta, &tai1, &tai2);
    if (status) {
        PyErr_SetString(erfaError, "internal error...");
        return NULL;
    }
    return Py_BuildValue("dd", tai1, tai2);
}

PyDoc_STRVAR(erfa_ut1tai_doc,
"ut1tai(ut11, ut12, dta) -> tai1, tai2\n\n"
"Time scale transformation: Universal Time, UT1, to\n"
"International Atomic Time, TAI.\n"
"The argument dta, i.e. UT1-TAI, is an observed quantity, and is\n"
"available from IERS tabulations.\n"
"Given:\n"
"   ut11,ut12     TT as a 2-part Julian Date\n"
"   dta         UT1-TAI in seconds\n"
"Returned:\n"
"   tai1,tai2   TAI as a 2-part Julian Date");

static PyObject *
erfa_ut1tt(PyObject *self, PyObject *args)
{
    double tt1, tt2, dt, ut11, ut12;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &ut11, &ut12, &dt)) {
        return NULL;
    }
    status = eraUt1tt(ut11, ut12, dt, &tt1, &tt2);
    if (status) {
        PyErr_SetString(erfaError, "internal error...");
        return NULL;
    }
    return Py_BuildValue("dd", tt1, tt2);
}

PyDoc_STRVAR(erfa_ut1tt_doc,
"ut1tt(ut11, ut12, dt) -> tt1, tt2\n\n"
"Time scale transformation: Universal Time, UT1, to\n"
"Terrestrial Time, TT.\n"
"Given:\n"
"   ut11,ut12   UT1 as a 2-part Julian Date\n"
"   dt          TT-UT1 in seconds, dt is classical Delta T\n"
"Returned:\n"
"   tt1,tt2     TT as a 2-part Julian Date");

static PyObject *
erfa_ut1utc(PyObject *self, PyObject *args)
{
    double utc1, utc2, dut1, ut11, ut12;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &ut11, &ut12, &dut1)) {
        return NULL;
    }
    status = eraUt1utc(ut11, ut12, dut1, &utc1, &utc2);
    if (status) {
        if (status == 1) {
            PyErr_SetString(erfaError, "dubious year");
            return NULL;
        }
        else if (status == -1) {
            PyErr_SetString(erfaError, "unacceptable date");
            return NULL;
        }
    }
    return Py_BuildValue("dd", utc1, utc2);
}

PyDoc_STRVAR(erfa_ut1utc_doc,
"ut1utc(ut11, ut12, dut1) -> utc1, utc2\n\n"
"Time scale transformation: Universal Time, UT1, to\n"
"Coordinated Universal Time, UTC.\n"
"Given:\n"
"   ut11,ut12   UT1 as a 2-part Julian Date\n"
"   dut1        UT1-UTC in seconds, Delta UT1\n"
"Returned:\n"
"   utc1,utc2   UTC as a 2-part Julian Date");

static PyObject *
erfa_utctai(PyObject *self, PyObject *args)
{
    double utc1, utc2, tai1, tai2;
    int status;
    if (!PyArg_ParseTuple(args, "dd", &utc1, &utc2)) {
        return NULL;
    }
    status = eraUtctai(utc1, utc2, &tai1, &tai2);
    if (status) {
        if (status == 1) {
            PyErr_SetString(erfaError, "dubious year");
            return NULL;
        }
        else if (status == -1) {
            PyErr_SetString(erfaError, "unacceptable date");
            return NULL;
        }
    }
    return Py_BuildValue("dd", tai1, tai2);
}

PyDoc_STRVAR(erfa_utctai_doc,
"utctai(utc1, utc2) -> tai1, tai2\n\n"
"Time scale transformation: Coordinated Universal Time, UTC, to\n"
"International Atomic Time, TAI.\n"
"Given:\n"
"   utc1,uc12   UTC as a 2-part Julian Date\n"
"Returned:\n"
"   tai1,tai2   TAI as a 2-part Julian Date");

static PyObject *
erfa_utcut1(PyObject *self, PyObject *args)
{
    double utc1, utc2, dut1, ut11, ut12;
    int status;
    if (!PyArg_ParseTuple(args, "ddd", &utc1, &utc2, &dut1)) {
        return NULL;
    }
    status = eraUtcut1(utc1, utc2, dut1, &ut11, &ut12);
    if (status) {
        if (status == 1) {
            PyErr_SetString(erfaError, "dubious year");
            return NULL;
        }
        else if (status == -1) {
            PyErr_SetString(erfaError, "unacceptable date");
            return NULL;
        }
    }
    return Py_BuildValue("dd", ut11, ut12);
}

PyDoc_STRVAR(erfa_utcut1_doc,
"utcut1(utc1, utc2, dut1) -> ut11, ut12\n\n"
"Time scale transformation: Coordinated Universal Time, UTC, to\n"
"Universal Time, UT1.\n"
"Given:\n"
"   utc1,utc2   UTC as a 2-part Julian Date\n"
"   dut1        UT1-UTC in seconds, Delta UT1\n"
"Returned:\n"
"   ut11,ut12   UT1 as a 2-part Julian Date");

static PyObject *
erfa_xy06(PyObject *self, PyObject *args)
{
    double d1, d2, x, y;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraXy06(d1, d2, &x, &y);
    return Py_BuildValue("dd", x, y);
}

PyDoc_STRVAR(erfa_xy06_doc,
"xy06(d1, d2) -> x, y\n\n"
"X,Y coordinates of celestial intermediate pole from series based\n"
"on IAU 2006 precession and IAU 2000A nutation.\n"
"Given:\n"
"   d1,d2   TT as a 2-part Julian Date\n"
"Returned:\n"
"   x,y     CIP X,Y coordinates");

static PyObject *
erfa_xys00a(PyObject *self, PyObject *args)
{
    double d1, d2, x, y,  s;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraXys00a(d1, d2, &x, &y, &s);
    return Py_BuildValue("ddd", x, y, s);
}

PyDoc_STRVAR(erfa_xys00a_doc,
"xys00a(d1, d2) -> x, y, s\n\n"
"For a given TT date, compute the X,Y coordinates of the Celestial\n"
"Intermediate Pole and the CIO locator s, using the IAU 2000A\n"
"precession-nutation model.\n"
"Given:\n"
"   d1,d2   TT as a 2-part Julian Date\n"
"Returned:\n"
"   x,y     Celestial Intermediate Pole\n"
"   s       the CIO locator s");

static PyObject *
erfa_xys00b(PyObject *self, PyObject *args)
{
    double d1, d2, x, y,  s;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraXys00b(d1, d2, &x, &y, &s);
    return Py_BuildValue("ddd", x, y, s);
}

PyDoc_STRVAR(erfa_xys00b_doc,
"xys00b(d1, d2) -> x, y, s\n\n"
"For a given TT date, compute the X,Y coordinates of the Celestial\n"
"Intermediate Pole and the CIO locator s, using the IAU 2000B\n"
"precession-nutation model.\n"
"Given:\n"
"   d1,d2   TT as a 2-part Julian Date\n"
"Returned:\n"
"   x,y     Celestial Intermediate Pole\n"
"   s       the CIO locator s");

static PyObject *
erfa_xys06a(PyObject *self, PyObject *args)
{
    double d1, d2, x, y,  s;
    if (!PyArg_ParseTuple(args, "dd", &d1, &d2)) {
        return NULL;
    }
    eraXys06a(d1, d2, &x, &y, &s);
    return Py_BuildValue("ddd", x, y, s);
}

PyDoc_STRVAR(erfa_xys06a_doc,
"xys06a(d1, d2) -> x, y, s\n\n"
"For a given TT date, compute the X,Y coordinates of the Celestial\n"
"Intermediate Pole and the CIO locator s, using the IAU 2006\n"
"precession and IAU 2000A nutation model.\n"
"Given:\n"
"   d1,d2   TT as a 2-part Julian Date\n"
"Returned:\n"
"   x,y     Celestial Intermediate Pole\n"
"   s       the CIO locator s");

static PyObject *
erfa_a2af(PyObject *self, PyObject *args)
{
    int ndp, idmsf[4];
    char sign;
    double a;
    if (!PyArg_ParseTuple(args, "id", &ndp, &a)) {
        return NULL;
    }
    eraA2af(ndp, a, &sign, idmsf);
    return Py_BuildValue("Ciiii",sign,idmsf[0],idmsf[1],idmsf[2],idmsf[3]);
}

PyDoc_STRVAR(erfa_a2af_doc,
"a2af(ndp, angle) -> sign, degrees, arcminutes, arcseconds, fraction\n\n"
"Decompose radians into degrees, arcminutes, arcseconds, fraction.\n"
"Given:\n"
"   ndp     resolution\n"
"   angle   angle in radians\n"
"Returned:\n"
"   sign    char    '+' or '-'\n"
"   idmsf[4]  degrees, arcminutes, arcseconds, fraction");

static PyObject *
erfa_a2tf(PyObject *self, PyObject *args)
{
    int ndp, ihmsf[4];
    char sign;
    double a;
    if (!PyArg_ParseTuple(args, "id", &ndp, &a)) {
        return NULL;
    }
    eraA2tf(ndp, a, &sign, ihmsf);
    return Py_BuildValue("Ciiii",sign,ihmsf[0],ihmsf[1],ihmsf[2],ihmsf[3]);
}

PyDoc_STRVAR(erfa_a2tf_doc,
"a2tf(ndp, angle) -> sign, hours, minutes, seconds, fraction\n\n"
"Decompose radians into hours, minutes, seconds, fraction.\n"
"Given:\n"
"   ndp     resolution\n"
"   angle   angle in radians\n"
"Returned:\n"
"   sign    '+' or '-'\n"
"   hours, minutes, seconds, fraction");

static PyObject *
erfa_af2a(PyObject *self, PyObject *args)
{
    int ideg, iamin;
    double asec, rad;
    int sign = '+';
    if (!PyArg_ParseTuple(args, "iid", &ideg, &iamin, &asec)) {
        return NULL;
    }
    if (ideg < 0) {
        sign = '-';
    }
    eraAf2a(sign, ideg, iamin, asec, &rad);
    return Py_BuildValue("d", rad);
}

PyDoc_STRVAR(erfa_af2a_doc,
"af2a(degrees, arcminutes, arcseconds) -> radians\n"
"Convert degrees, arcminutes, arcseconds to radians.\n"
"Given:\n"
"   s       sign:  '-' = negative, otherwise positive\n"
"   ideg    degrees\n"
"   iamin   arcminutes\n"
"   asec    arcseconds\n"
"Returned:\n"
"   rad     angle in radians");

static PyObject *
erfa_anp(PyObject *self, PyObject *arg)
{
    double a = PyFloat_AsDouble(arg);
    if ((a == -1.0) && PyErr_Occurred()) {
        PyErr_SetString(erfaError, "cannot convert angle to float!");
        return NULL;
    }
    return Py_BuildValue("d", eraAnp(a));
}

PyDoc_STRVAR(erfa_anp_doc,
"anp(a) -> 0 <= a < 2pi\n\n"
"Normalize angle into the range 0 <= a < 2pi.\n"
"Given:\n"
"   a   angle (radians)\n"
"Returned:\n"
"   a   angle in range 0-2pi");

static PyObject *
erfa_anpm(PyObject *self, PyObject *arg)
{
    double a = PyFloat_AsDouble(arg);
    if ((a == -1.0) && PyErr_Occurred()) {
        PyErr_SetString(erfaError, "cannot convert angle to float!");
        return NULL;
    }
    return Py_BuildValue("d", eraAnpm(a));
}

PyDoc_STRVAR(erfa_anpm_doc,
"anpm(a) -> -pi <= a < +pi\n\n"
"Normalize angle into the range -pi <= a < +pi.\n"
"Given:\n"
"   a   angle (radians)\n"
"Returned:\n"
"   a   angle in range 0-2pi");

static PyObject *
erfa_c2s(PyObject *self, PyObject *args)
{
    double p0, p1, p2, theta, phi, p[3];
    int ok;
    ok = PyArg_ParseTuple(args, "(ddd)", &p0, &p1, &p2);
    if (ok) {
        p[0]=p0;
        p[1]=p1;
        p[2]=p2;
    }
    eraC2s(p, &theta, &phi);
    return Py_BuildValue("dd", theta, phi);
}

PyDoc_STRVAR(erfa_c2s_doc,
"c2s((p0,p1,p2)) -> theta, phi\n\n"
"P-vector to spherical coordinates.\n"
"Given:\n"
"   p[3]    p-vector\n"
"Returned:\n"
"   theta   longitude angle (radians)\n"
"   phi     latitude angle (radians)");

static PyObject *
erfa_cp(PyObject *self, PyObject *args)
{
    double p[3], c[3];
    double p0,p1,p2;
    if (!PyArg_ParseTuple(args, "(ddd)",
                        &p0,&p1,&p2)) {
        return NULL;
    }
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    eraCp(p, c);
    return Py_BuildValue("ddd",
    c[0],c[1],c[2]);
}

PyDoc_STRVAR(erfa_cp_doc,
"cp(p) -> c\n\n"
"Copy a p-vector.\n"
"Given:\n"
"   p           p-vector to be copied\n"
"  Returned:\n"
"   c           copy");

static PyObject *
erfa_cpv(PyObject *self, PyObject *args)
{
    double pv[2][3], c[2][3];
    double p00,p01,p02,p10,p11,p12;
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)",
                        &p00,&p01,&p02,&p10,&p11,&p12)) {
        return NULL;
    }
    pv[0][0] = p00;
    pv[0][1] = p01;
    pv[0][2] = p02;
    pv[1][0] = p10;
    pv[1][1] = p11;
    pv[1][2] = p12;
    eraCpv(pv, c);
    return Py_BuildValue("(ddd)(ddd)",
    c[0][0],c[0][1],c[0][2],c[1][0],c[1][1],c[1][2]);    
}

PyDoc_STRVAR(erfa_cpv_doc,
"cp(pv) -> c\n\n"
"Copy a position/velocity vector.\n"
"Given:\n"
"   pv          position/velocity vector to be copied\n"
"  Returned:\n"
"   c           copy");

static PyObject *
erfa_cr(PyObject *self, PyObject *args)
{
    double r[3][3], c[3][3];
    double r00,r01,r02,r10,r11,r12,r20,r21,r22;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))",
                        &r00,&r01,&r02,&r10,&r11,&r12,&r20,&r21,&r22)) {
        return NULL;
    }
    r[0][0] = r00;
    r[0][1] = r01;
    r[0][2] = r02;
    r[1][0] = r10;
    r[1][1] = r11;
    r[1][2] = r12;
    r[2][0] = r20;
    r[2][1] = r21;
    r[2][2] = r22;
    eraCr(r, c);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
    c[0][0],c[0][1],c[0][2],
    c[1][0],c[1][1],c[1][2],
    c[2][0],c[2][1],c[2][2]);    
}

PyDoc_STRVAR(erfa_cr_doc,
"cr(r) -> c\n\n"
"Copy an r-matrix.\n"
"Given:\n"
"   r           r-matrix to be copied\n"
"  Returned:\n"
"   c           copy");

static PyObject *
erfa_d2tf(PyObject *self, PyObject *args)
{
    double d;
    char sign = '+';
    int ndp, df[4];
    if (! PyArg_ParseTuple(args, "id",&ndp, &d)) {
        return NULL;
    }
    eraD2tf(ndp, d, &sign, df);
    return Py_BuildValue("Ciiii", sign,df[0],df[1],df[2],df[3]);
}

PyDoc_STRVAR(erfa_d2tf_doc,
"d2tf(n, d) -> +/-,h, m, s, f\n\n"
"Decompose days to hours, minutes, seconds, fraction.\n"
"Given:\n"
"    n          resolution\n"
"    d          interval in days\n"
"Returned:\n"
"    sign    '+' or '-'\n"
"    h          hours\n"
"    m          minutes\n"
"    s          seconds\n"
"    f          fraction");

static PyObject *
erfa_p2pv(PyObject *self, PyObject *args)
{
    double p[3], pv[2][3];
    double p0,p1,p2;
    if (!PyArg_ParseTuple(args, "(ddd)",
                        &p0,&p1,&p2)) {
        return NULL;
    }
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    eraP2pv(p, pv);
    return Py_BuildValue("(ddd)(ddd)",
    pv[0][0],pv[0][1],pv[0][2],pv[1][0],pv[1][1],pv[1][2]);
}

PyDoc_STRVAR(erfa_p2pv_doc,
"p2pv(p) -> pv\n\n"
"Extend a p-vector to a pv-vector by appending a zero velocity.\n"
"Given:\n"
"    p          p-vector\n"
"Returned:\n"
"    pv         pv-vector");

static PyObject *
erfa_p2s(PyObject *self, PyObject *args)
{
    double p[3], p0, p1, p2, theta, phi, r;
    if (!PyArg_ParseTuple(args, "(ddd)",&p0, &p1, &p2)) {
        return NULL;
    }
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    eraP2s(p, &theta, &phi, &r);
    return Py_BuildValue("ddd", theta, phi, r);
}

PyDoc_STRVAR(erfa_p2s_doc,
"p2s(p) -> theta, phi, r\n\n"
"P-vector to spherical polar coordinates.\n"
"Given:\n"
"    p          p-vector\n"
"Returned:\n"
"    theta      longitude angle (radians)\n"
"    phi        latitude angle (radians)\n"
"    r          radial distance");

static PyObject *
erfa_pap(PyObject *self, PyObject *args)
{
    double a[3], b[3], theta;
    double a0,a1,a2,b0,b1,b2;
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)",&a0,&a1,&a2,&b0,&b1,&b2)) {
        return NULL;
    }
    a[0] = a0;
    a[1] = a1;
    a[2] = a2;
    b[0] = b0;
    b[1] = b1;
    b[2] = b2;
    theta = eraPap(a, b);
    return Py_BuildValue("d", theta);
}

PyDoc_STRVAR(erfa_pap_doc,
"pap(a, b) -> theta\n\n"
"Position-angle from two p-vectors.\n"
"Given:\n"
"   a       direction of reference point\n"
"   b       direction of point whose PA is required\n"
"Returned:\n"
"   theta   position angle of b with respect to a (radians)");

static PyObject *
erfa_pas(PyObject *self, PyObject *args)
{
    double al, ap, bl, bp, p;
    if (!PyArg_ParseTuple(args, "dddd",&al, &ap, &bl, &bp)) {
        return NULL;
    }
    p = eraPas(al, ap, bl, bp);
    return Py_BuildValue("d", p);
}

PyDoc_STRVAR(erfa_pas_doc,
"pas(al, ap, bl, bp) -> p\n\n"
"Position-angle from spherical coordinates.\n"
"Given:\n"
"   al  longitude of point A (e.g. RA) in radians\n"
"   ap  latitude of point A (e.g. Dec) in radians\n"
"   bl  longitude of point B\n"
"   bp  latitude of point B\n"
"Returned:\n"
"   p   position angle of B with respect to A");

static PyObject *
erfa_pdp(PyObject *self, PyObject *args)
{
    double a[3], b[3], ab;
    double a0,a1,a2, b0,b1,b2;
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)",&a0,&a1,&a2,&b0,&b1,&b2)) {
        return NULL;
    }
    a[0] = a0;
    a[1] = a1;
    a[2] = a2;
    b[0] = b0;
    b[1] = b1;
    b[2] = b2;
    ab = eraPdp(a,b);
    return Py_BuildValue("d", ab);
}

PyDoc_STRVAR(erfa_pdp_doc,
"pdp(a, b -> a.b\n\n"
"p-vector inner (=scalar=dot) product.\n"
"Given:\n"
"   a       first p-vector\n"
"   b       second p-vector\n"
"Returned:\n"
"   ab      a . b");

static PyObject *
erfa_pm(PyObject *self, PyObject *args)
{
    double p[3], p0, p1, p2, m;
    if (!PyArg_ParseTuple(args, "(ddd)",&p0, &p1, &p2)) {
        return NULL;
    }
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    m = eraPm(p);
    return Py_BuildValue("d", m);
}

PyDoc_STRVAR(erfa_pm_doc,
"pm(p) -> modulus\n\n"
"Modulus of p-vector.\n"
"Given:\n"
"   p       p-vector\n"
"Returned:\n"
"   m       modulus");

static PyObject *
erfa_pmp(PyObject *self, PyObject *args)
{
    double a[3], b[3], amb[3];
    double a0,a1,a2,b0,b1,b2;
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)",
        &a0, &a1, &a2, &b0, &b1, &b2)) {
        return NULL;
    }
    a[0] = a0;
    a[1] = a1;
    a[2] = a2;
    b[0] = b0;
    b[1] = b1;
    b[2] = b2;
    eraPmp(a, b, amb);
    return Py_BuildValue("ddd", amb[0], amb[1], amb[2]);
}

PyDoc_STRVAR(erfa_pmp_doc,
"pmp(a, b) -> amb = a-b\n\n"
"P-vector subtraction.\n"
"Given:\n"
"   a           first p-vector\n"
"   b           second p-vector\n"
"Returned:\n"
"   amb         a - b");

static PyObject *
erfa_pn(PyObject *self, PyObject *args)
{
    double p[3], p0, p1, p2, r, u[3];
    if (!PyArg_ParseTuple(args, "(ddd)",&p0, &p1, &p2)) {
        return NULL;
    }
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    eraPn(p, &r, u);
    return Py_BuildValue("d(ddd)", r, u[0], u[1], u[2]);    
}

PyDoc_STRVAR(erfa_pn_doc,
"pn(p) -> r,u\n\n"
"Convert a p-vector into modulus and unit vector.\n"
"Given:\n"
"   p           p-vector\n"
"Returned:\n"
"   r           modulus\n"
"   u           unit vector");

static PyObject *
erfa_ppp(PyObject *self, PyObject *args)
{
    double a[3], b[3], apb[3];
    double a0,a1,a2,b0,b1,b2;
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)",
        &a0, &a1, &a2, &b0, &b1, &b2)) {
        return NULL;
    }
    a[0] = a0;
    a[1] = a1;
    a[2] = a2;
    b[0] = b0;
    b[1] = b1;
    b[2] = b2;
    eraPpp(a, b, apb);
    return Py_BuildValue("ddd", apb[0], apb[1], apb[2]);
}

PyDoc_STRVAR(erfa_ppp_doc,
"ppp(a, b) -> apb = a+b\n\n"
"P-vector addition.\n"
"Given:\n"
"   a           first p-vector\n"
"   b           second p-vector\n"
"Returned:\n"
"   apb         a + b");

static PyObject *
erfa_ppsp(PyObject *self, PyObject *args)
{
    double a[3], b[3], apsb[3];
    double a0,a1,a2,s,b0,b1,b2;
    if (!PyArg_ParseTuple(args, "(ddd)d(ddd)",
        &a0, &a1, &a2, &s, &b0, &b1, &b2)) {
        return NULL;
    }
    a[0] = a0;
    a[1] = a1;
    a[2] = a2;
    b[0] = b0;
    b[1] = b1;
    b[2] = b2;
    eraPpsp(a, s, b, apsb);
    return Py_BuildValue("ddd", apsb[0], apsb[1], apsb[2]);
}

PyDoc_STRVAR(erfa_ppsp_doc,
"ppsp(a, s, b) -> apsb = a + s*b\n\n"
"P-vector plus scaled p-vector.\n"
"Given:\n"
"   a           first p-vector\n"
"   s           scalar (multiplier for b)\n"
"   b           second p-vector\n"
"Returned:\n"
"   apsb        a + s*b");

static PyObject *
erfa_pv2p(PyObject *self, PyObject *args)
{
    double pv[2][3], p[3];
    double p00,p01,p02,p10,p11,p12;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))",
                        &p00,&p01,&p02,&p10,&p11,&p12)) {
        return NULL;
    }
    pv[0][0] = p00;
    pv[0][1] = p01;
    pv[0][2] = p02;
    pv[1][0] = p10;
    pv[1][1] = p11;
    pv[1][2] = p12;
    eraPv2p(pv, p);
    return Py_BuildValue("ddd", p[0],p[1],p[2]);    
}

PyDoc_STRVAR(erfa_pv2p_doc,
"pv2p(pv) -> p\n\n"
"Discard velocity component of a pv-vector.\n"
"Given:\n"
"   pv          pv-vector\n"
"Returned:\n"
"   p           p-vector");

static PyObject *
erfa_pv2s(PyObject *self, PyObject *args)
{
    double pv[2][3], theta, phi, r, td, pd, rd;
    double pv00, pv01, pv02, pv10, pv11, pv12;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))",
        &pv00, &pv01, &pv02, &pv10, &pv11, &pv12)) {
        return NULL;
    }
    pv[0][0] = pv00;
    pv[0][1] = pv01;
    pv[0][2] = pv02;
    pv[1][0] = pv10;
    pv[1][1] = pv11;
    pv[1][2] = pv12;
    eraPv2s(pv, &theta, &phi, &r, &td, &pd, &rd);
    return Py_BuildValue("dddddd", theta, phi, r, td, pd, rd);
}

PyDoc_STRVAR(erfa_pv2s_doc,
"pv2s(pv) -> theta, phi, r, td, pd, rd\n\n"
"Convert position/velocity from Cartesian to spherical coordinates.\n"
"Given:\n"
"   pv          pv-vector\n"
"Returned:\n"
"   theta       longitude angle (radians)\n"
"   phi         latitude angle (radians)\n"
"   r           radial distance\n"
"   td          rate of change of theta\n"
"   pd          rate of change of phi\n"
"   rd          rate of change of r");

static PyObject *
erfa_pvdpv(PyObject *self, PyObject *args)
{
    double pva[2][3], pvb[2][3], adb[2];
    double pva00, pva01, pva02, pva10, pva11, pva12;
    double pvb00, pvb01, pvb02, pvb10, pvb11, pvb12;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))((ddd)(ddd))",
        &pva00, &pva01, &pva02, &pva10, &pva11, &pva12,
        &pvb00, &pvb01, &pvb02, &pvb10, &pvb11, &pvb12)) {
        return NULL;
    }
    pva[0][0] = pva00;
    pva[0][1] = pva01;
    pva[0][2] = pva02;
    pva[1][0] = pva10;
    pva[1][1] = pva11;
    pva[1][2] = pva12;
    pvb[0][0] = pvb00;
    pvb[0][1] = pvb01;
    pvb[0][2] = pvb02;
    pvb[1][0] = pvb10;
    pvb[1][1] = pvb11;
    pvb[1][2] = pvb12;
    eraPvdpv(pva, pvb, adb);
    return Py_BuildValue("dd", adb[0], adb[1]);
}

PyDoc_STRVAR(erfa_pvdpv_doc,
"pvdpv(a, b) -> adb = a.b\n\n"
"Inner (=scalar=dot) product of two pv-vectors.\n"
"Given:\n"
"   a           first pv-vector\n"
"   b           second pv-vector\n"
"Returned:\n"
"   adb         a . b");

static PyObject *
erfa_pvm(PyObject *self, PyObject *args)
{
    double pv[2][3], r, s;
    double pv00, pv01, pv02, pv10, pv11, pv12;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))",
        &pv00, &pv01, &pv02, &pv10, &pv11, &pv12)) {
        return NULL;
    }
    pv[0][0] = pv00;
    pv[0][1] = pv01;
    pv[0][2] = pv02;
    pv[1][0] = pv10;
    pv[1][1] = pv11;
    pv[1][2] = pv12;
    eraPvm(pv, &r, &s);
    return Py_BuildValue("dd", r, s);
}

PyDoc_STRVAR(erfa_pvm_doc,
"pvm(pv) -> r,s\n\n"
"Modulus of pv-vector.\n"
"Given:\n"
"   pv[2][3]    pv-vector\n"
"Returned:\n"
"   r           modulus of position component\n"
"   s           modulus of velocity component");

static PyObject *
erfa_pvmpv(PyObject *self, PyObject *args)
{
    double pva[2][3], pvb[2][3], amb[2][3];
    double pva00, pva01, pva02, pva10, pva11, pva12;
    double pvb00, pvb01, pvb02, pvb10, pvb11, pvb12;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))((ddd)(ddd))",
        &pva00, &pva01, &pva02, &pva10, &pva11, &pva12,
        &pvb00, &pvb01, &pvb02, &pvb10, &pvb11, &pvb12)) {
        return NULL;
    }
    pva[0][0] = pva00;
    pva[0][1] = pva01;
    pva[0][2] = pva02;
    pva[1][0] = pva10;
    pva[1][1] = pva11;
    pva[1][2] = pva12;
    pvb[0][0] = pvb00;
    pvb[0][1] = pvb01;
    pvb[0][2] = pvb02;
    pvb[1][0] = pvb10;
    pvb[1][1] = pvb11;
    pvb[1][2] = pvb12;
    eraPvmpv(pva, pvb, amb);
    return Py_BuildValue("(ddd)(ddd)",
                           amb[0][0], amb[0][1], amb[0][2],
                           amb[1][0], amb[1][1], amb[1][2]);
}

PyDoc_STRVAR(erfa_pvmpv_doc,
"pvmpv(a, b) -> amb = a-b\n\n"
"Subtract one pv-vector from another.\n"
"Given:\n"
"   a           first pv-vector\n"
"   b           second pv-vector\n"
"Returned:\n"
"   amb         a - b");

static PyObject *
erfa_pvppv(PyObject *self, PyObject *args)
{
    double pva[2][3], pvb[2][3], apb[2][3];
    double pva00, pva01, pva02, pva10, pva11, pva12;
    double pvb00, pvb01, pvb02, pvb10, pvb11, pvb12;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))((ddd)(ddd))",
        &pva00, &pva01, &pva02, &pva10, &pva11, &pva12,
        &pvb00, &pvb01, &pvb02, &pvb10, &pvb11, &pvb12)) {
        return NULL;
    }
    pva[0][0] = pva00;
    pva[0][1] = pva01;
    pva[0][2] = pva02;
    pva[1][0] = pva10;
    pva[1][1] = pva11;
    pva[1][2] = pva12;
    pvb[0][0] = pvb00;
    pvb[0][1] = pvb01;
    pvb[0][2] = pvb02;
    pvb[1][0] = pvb10;
    pvb[1][1] = pvb11;
    pvb[1][2] = pvb12;
    eraPvppv(pva, pvb, apb);
    return Py_BuildValue("(ddd)(ddd)",
                           apb[0][0], apb[0][1], apb[0][2],
                           apb[1][0], apb[1][1], apb[1][2]);
}

PyDoc_STRVAR(erfa_pvppv_doc,
"pvppv(a, b) -> apb = a+b\n\n"
"Add one pv-vector to another.\n"
"Given:\n"
"   a           first pv-vector\n"
"   b           second pv-vector\n"
"Returned:\n"
"   apb         a + b");

static PyObject *
erfa_pvu(PyObject *self, PyObject *args)
{
    double pv[2][3], dt, upv[2][3];
    double pv00, pv01, pv02, pv10, pv11, pv12;
    if (!PyArg_ParseTuple(args, "d((ddd)(ddd))",
        &dt, &pv00, &pv01, &pv02, &pv10, &pv11, &pv12)) {
        return NULL;
    }
    pv[0][0] = pv00;
    pv[0][1] = pv01;
    pv[0][2] = pv02;
    pv[1][0] = pv10;
    pv[1][1] = pv11;
    pv[1][2] = pv12;
    eraPvu(dt, pv, upv);
    return Py_BuildValue("(ddd)(ddd)",
                           upv[0][0], upv[0][1], upv[0][2],
                           upv[1][0], upv[1][1], upv[1][2]);
}

PyDoc_STRVAR(erfa_pvu_doc,
"pvu(dt, pv) -> upv\n\n"
"Update a pv-vector.\n"
"Given:\n"
"   dt          time interval\n"
"   pv          pv-vector\n"
"Returned:\n"
"   upv         p updated, v unchanged");

static PyObject *
erfa_pvup(PyObject *self, PyObject *args)
{
    double pv[2][3], dt, p[3];
    double pv00, pv01, pv02, pv10, pv11, pv12;
    if (!PyArg_ParseTuple(args, "d((ddd)(ddd))",
        &dt, &pv00, &pv01, &pv02, &pv10, &pv11, &pv12)) {
        return NULL;
    }
    pv[0][0] = pv00;
    pv[0][1] = pv01;
    pv[0][2] = pv02;
    pv[1][0] = pv10;
    pv[1][1] = pv11;
    pv[1][2] = pv12;
    eraPvup(dt, pv, p);
    return Py_BuildValue("ddd", p[0], p[1], p[2]);
}

PyDoc_STRVAR(erfa_pvup_doc,
"pvup(dt, pv) -> p\n\n"
"Update a pv-vector, discarding the velocity component.\n"
"Given:\n"
"   dt          time interval\n"
"   pv          pv-vector\n"
"Returned:\n"
"   p           p-vector");

static PyObject *
erfa_pvxpv(PyObject *self, PyObject *args)
{
    double pva[2][3], pvb[2][3], axb[2][3];
    double pva00, pva01, pva02, pva10, pva11, pva12;
    double pvb00, pvb01, pvb02, pvb10, pvb11, pvb12;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd))((ddd)(ddd))",
        &pva00, &pva01, &pva02, &pva10, &pva11, &pva12,
        &pvb00, &pvb01, &pvb02, &pvb10, &pvb11, &pvb12)) {
        return NULL;
    }
    pva[0][0] = pva00;
    pva[0][1] = pva01;
    pva[0][2] = pva02;
    pva[1][0] = pva10;
    pva[1][1] = pva11;
    pva[1][2] = pva12;
    pvb[0][0] = pvb00;
    pvb[0][1] = pvb01;
    pvb[0][2] = pvb02;
    pvb[1][0] = pvb10;
    pvb[1][1] = pvb11;
    pvb[1][2] = pvb12;
    eraPvxpv(pva, pvb, axb);
    return Py_BuildValue("(ddd)(ddd)",
                           axb[0][0], axb[0][1], axb[0][2],
                           axb[1][0], axb[1][1], axb[1][2]);
}

PyDoc_STRVAR(erfa_pvxpv_doc,
"pvxpv(a, b) -> axb = a x b\n\n"
"Outer (=vector=cross) product of two pv-vectors.\n"
"Given:\n"
"   a           first pv-vector\n"
"   b           second pv-vector\n"
"Returned:\n"
"   axb         a x b");

static PyObject *
erfa_pxp(PyObject *self, PyObject *args)
{
    double a[3], b[3], axb[3];
    double a0,a1,a2,b0,b1,b2;
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)",
        &a0, &a1, &a2, &b0, &b1, &b2)) {
        return NULL;
    }
    a[0] = a0;
    a[1] = a1;
    a[2] = a2;
    b[0] = b0;
    b[1] = b1;
    b[2] = b2;
    eraPxp(a, b, axb);
    return Py_BuildValue("ddd", axb[0], axb[1], axb[2]);
}

PyDoc_STRVAR(erfa_pxp_doc,
"pxp(a, b) -> axb = a x b\n\n"
"p-vector outer (=vector=cross) product.\n"
"Given:\n"
"   a           first p-vector\n"
"   b           second p-vector\n"
"Returned:\n"
"   axb         a x b");

static PyObject *
erfa_rm2v(PyObject *self, PyObject *args)
{
    double r[3][3], w[3];
    double r00,r01,r02,r10,r11,r12,r20,r21,r22;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))",
                        &r00,&r01,&r02,&r10,&r11,&r12,&r20,&r21,&r22)) {
        return NULL;
    }
    r[0][0] = r00;
    r[0][1] = r01;
    r[0][2] = r02;
    r[1][0] = r10;
    r[1][1] = r11;
    r[1][2] = r12;
    r[2][0] = r20;
    r[2][1] = r21;
    r[2][2] = r22;
    eraRm2v(r, w);
    return Py_BuildValue("ddd", w[0], w[1], w[2]);    
}

PyDoc_STRVAR(erfa_rm2v_doc,
"rm2v(r) -> w\n\n"
"Express an r-matrix as an r-vector.\n"
"Given:\n"
"   r          rotation matrix\n"
"Returned:\n"
"   w          rotation vector");

static PyObject *
erfa_rv2m(PyObject *self, PyObject *args)
{    
    double r[3][3], w[3];
    double w0, w1, w2;
    if (!PyArg_ParseTuple(args, "(ddd)", &w0,&w1,&w2)) {
        return NULL;
    }
    w[0] = w0;
    w[1] = w1;
    w[2] = w2;
    eraRv2m(w, r);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                            r[0][0],r[0][1],r[0][2],
                            r[1][0],r[1][1],r[1][2],
                            r[2][0],r[2][1],r[2][2]);        
}

PyDoc_STRVAR(erfa_rv2m_doc,
"rv2m(w) -> r\n\n"
"Form the r-matrix corresponding to a given r-vector.\n"
"Given:\n"
"   w           rotation vector\n"
"Returned:\n"
"   r           rotation matrix");

static PyObject *
erfa_rx(PyObject *self, PyObject *args)
{
    double r[3][3], phi;
    double r00,r01,r02,r10,r11,r12,r20,r21,r22;
    if (!PyArg_ParseTuple(args, "d((ddd)(ddd)(ddd))", &phi,
                        &r00,&r01,&r02,&r10,&r11,&r12,&r20,&r21,&r22)) {
        return NULL;
    }
    r[0][0] = r00;
    r[0][1] = r01;
    r[0][2] = r02;
    r[1][0] = r10;
    r[1][1] = r11;
    r[1][2] = r12;
    r[2][0] = r20;
    r[2][1] = r21;
    r[2][2] = r22;
    eraRx(phi, r);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                            r[0][0],r[0][1],r[0][2],
                            r[1][0],r[1][1],r[1][2],
                            r[2][0],r[2][1],r[2][2]);    
}

PyDoc_STRVAR(erfa_rx_doc,
"rx(phi, r) -> r\n\n"
"Rotate an r-matrix about the x-axis.\n"
"Given:\n"
"   phi         angle (radians)\n"
"Given and returned:\n"
"   r           r-matrix, rotated");

static PyObject *
erfa_rxp(PyObject *self, PyObject *args)
{
    double r[3][3], p[3], rp[3];
    double r00,r01,r02,r10,r11,r12,r20,r21,r22, p0, p1, p2;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))(ddd)", 
                                &r00,&r01,&r02,
                                &r10,&r11,&r12,
                                &r20,&r21,&r22,
                                &p0, &p1, &p2)) {
        return NULL;
    }
    r[0][0] = r00;
    r[0][1] = r01;
    r[0][2] = r02;
    r[1][0] = r10;
    r[1][1] = r11;
    r[1][2] = r12;
    r[2][0] = r20;
    r[2][1] = r21;
    r[2][2] = r22;
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    eraRxp(r, p, rp);
    return Py_BuildValue("ddd", rp[0], rp[1], rp[2]);
}

PyDoc_STRVAR(erfa_rxp_doc,
"rxp(r, p) -> rp\n\n"
"Multiply a p-vector by an r-matrix.\n"
"Given:\n"
"   r           r-matrix\n"
"   p           p-vector\n"
"Returned:\n"
"   rp          r * p");

static PyObject *
erfa_rxpv(PyObject *self, PyObject *args)
{
    double r[3][3], pv[2][3], rpv[2][3];
    double r00,r01,r02,r10,r11,r12,r20,r21,r22;
    double pv00, pv01, pv02, pv10, pv11, pv12;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))((ddd)(ddd))", 
                                &r00,&r01,&r02,
                                &r10,&r11,&r12,
                                &r20,&r21,&r22,
                                &pv00, &pv01, &pv02, &pv10, &pv11, &pv12)) {
        return NULL;
    }
    r[0][0] = r00;
    r[0][1] = r01;
    r[0][2] = r02;
    r[1][0] = r10;
    r[1][1] = r11;
    r[1][2] = r12;
    r[2][0] = r20;
    r[2][1] = r21;
    r[2][2] = r22;
    pv[0][0] = pv00;
    pv[0][1] = pv01;
    pv[0][2] = pv02;
    pv[1][0] = pv10;
    pv[1][1] = pv11;
    pv[1][2] = pv12;
    eraRxpv(r, pv, rpv);
    return Py_BuildValue("(ddd)(ddd)",
                          rpv[0][0], rpv[0][1], rpv[0][2],
                          rpv[1][0], rpv[1][1], rpv[1][2]);
}

PyDoc_STRVAR(erfa_rxpv_doc,
"rxpv(r, pv) -> rpv\n\n"
"Multiply a pv-vector by an r-matrix.\n"
"Given:\n"
"   r           r-matrix\n"
"   pv          pv-vector\n"
"Returned:\n"
"   rpv         r * pv");

static PyObject *
erfa_rxr(PyObject *self, PyObject *args)
{
    double a[3][3], b[3][3], atb[3][3];
    double a00,a01,a02,a10,a11,a12,a20,a21,a22;
    double b00,b01,b02,b10,b11,b12,b20,b21,b22;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))((ddd)(ddd)(ddd))",
                        &a00,&a01,&a02,&a10,&a11,&a12,&a20,&a21,&a22,
                        &b00,&b01,&b02,&b10,&b11,&b12,&b20,&b21,&b22)) {
        return NULL;
    }
    a[0][0] = a00;
    a[0][1] = a01;
    a[0][2] = a02;
    a[1][0] = a10;
    a[1][1] = a11;
    a[1][2] = a12;
    a[2][0] = a20;
    a[2][1] = a21;
    a[2][2] = a22;

    b[0][0] = b00;
    b[0][1] = b01;
    b[0][2] = b02;
    b[1][0] = b10;
    b[1][1] = b11;
    b[1][2] = b12;
    b[2][0] = b20;
    b[2][1] = b21;
    b[2][2] = b22;
    eraRxr(a, b, atb);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                           atb[0][0],atb[0][1],atb[0][2],
                           atb[1][0],atb[1][1],atb[1][2],
                           atb[2][0],atb[2][1],atb[2][2]);    
}

PyDoc_STRVAR(erfa_rxr_doc,
"rxr(a, b -> atb\n\n"
"Multiply two r-matrices.\n"
"Given:\n"
"   a           first r-matrix\n"
"   b           second r-matrix\n"
"Returned:\n"
"   atb         a * b");

static PyObject *
erfa_ry(PyObject *self, PyObject *args)
{
    double r[3][3], theta;
    double r00,r01,r02,r10,r11,r12,r20,r21,r22;
    if (!PyArg_ParseTuple(args, "d((ddd)(ddd)(ddd))", &theta,
                        &r00,&r01,&r02,&r10,&r11,&r12,&r20,&r21,&r22)) {
        return NULL;
    }
    r[0][0] = r00;
    r[0][1] = r01;
    r[0][2] = r02;
    r[1][0] = r10;
    r[1][1] = r11;
    r[1][2] = r12;
    r[2][0] = r20;
    r[2][1] = r21;
    r[2][2] = r22;
    eraRy(theta, r);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                            r[0][0],r[0][1],r[0][2],
                            r[1][0],r[1][1],r[1][2],
                            r[2][0],r[2][1],r[2][2]);    
}

PyDoc_STRVAR(erfa_ry_doc,
"ry(theta, r) -> r\n\n"
"Rotate an r-matrix about the y-axis.\n"
"Given:\n"
"   theta       angle (radians)\n"
"Given and returned:\n"
"   r           r-matrix, rotated");

static PyObject *
erfa_rz(PyObject *self, PyObject *args)
{
    double r[3][3], psi;
    double r00,r01,r02,r10,r11,r12,r20,r21,r22;
    if (!PyArg_ParseTuple(args, "d((ddd)(ddd)(ddd))", &psi,
                        &r00,&r01,&r02,&r10,&r11,&r12,&r20,&r21,&r22)) {
        return NULL;
    }
    r[0][0] = r00;
    r[0][1] = r01;
    r[0][2] = r02;
    r[1][0] = r10;
    r[1][1] = r11;
    r[1][2] = r12;
    r[2][0] = r20;
    r[2][1] = r21;
    r[2][2] = r22;
    eraRz(psi, r);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                            r[0][0],r[0][1],r[0][2],
                            r[1][0],r[1][1],r[1][2],
                            r[2][0],r[2][1],r[2][2]);    
}

PyDoc_STRVAR(erfa_rz_doc,
"rz(psi, r) -> r\n\n"
"Rotate an r-matrix about the z-axis.\n"
"Given:\n"
"   psi         angle (radians)\n"
"Given and returned:\n"
"   r           r-matrix, rotated");

static PyObject *
erfa_s2c(PyObject *self, PyObject *args)
{
    double theta, phi, c[3];
    if (!PyArg_ParseTuple(args, "dd", &theta, &phi)) {
        return NULL;
    }
    eraS2c(theta, phi, c);
    return Py_BuildValue("ddd", c[0], c[1], c[2]);
}

PyDoc_STRVAR(erfa_s2c_doc,
"s2c(theta, phi) -> c\n\n"
"Convert spherical coordinates to Cartesian.\n"
"Given:\n"
"    theta   longitude angle (radians)\n"
"    phi     latitude angle (radians)\n"
"Returned:\n"
"    c       direction cosines");

static PyObject *
erfa_s2p(PyObject *self, PyObject *args)
{
    double theta, phi, r, p[3];
    if (!PyArg_ParseTuple(args, "ddd", &theta, &phi, &r)) {
        return NULL;
    }
    eraS2p(theta, phi, r, p);
    return Py_BuildValue("ddd", p[0], p[1], p[2]);
}

PyDoc_STRVAR(erfa_s2p_doc,
"s2p(theta, phi, r) -> p\n\n"
"Convert spherical polar coordinates to p-vector.\n"
"Given:\n"
"   theta   longitude angle (radians)\n"
"   phi     latitude angle (radians)\n"
"   r       radial distance\n"
"Returned:\n"
"   p       direction cosines");

static PyObject *
erfa_s2pv(PyObject *self, PyObject *args)
{
    double theta, phi, r, td, pd, rd, pv[2][3];
    if (!PyArg_ParseTuple(args, "dddddd",
        &theta, &phi, &r, &td, &pd, &rd)) {
        return NULL;
    }
    eraS2pv(theta, phi, r, td, pd, rd, pv);
    return Py_BuildValue("(ddd)(ddd)",
                        pv[0][0], pv[0][1], pv[0][2],
                        pv[1][0], pv[1][1], pv[1][2]);
}

PyDoc_STRVAR(erfa_s2pv_doc,
"s2pv(theta, phi, r, td, pd, rd) -> pv\n\n"
"Convert position/velocity from spherical to Cartesian coordinates.\n"
"Given:\n"
"   theta   longitude angle (radians)\n"
"   phi     latitude angle (radians)\n"
"   r       radial distance\n"
"   td      rate of change of theta\n"
"   pd      rate of change of phi\n"
"   rd      rate of change of r\n"
"Returned:\n"
"   pv      pv-vector");

static PyObject *
erfa_s2xpv(PyObject *self, PyObject *args)
{
    double s1, s2, pv[2][3], spv[2][3];
    double pv00, pv01, pv02, pv10, pv11, pv12;
    if (!PyArg_ParseTuple(args, "dd((ddd)(ddd))",
        &s1, &s2, &pv00, &pv01, &pv02, &pv10, &pv11, &pv12)) {
        return NULL;
    }
    pv[0][0] = pv00;
    pv[0][1] = pv01;
    pv[0][2] = pv02;
    pv[1][0] = pv10;
    pv[1][1] = pv11;
    pv[1][2] = pv12;
    eraS2xpv(s1, s2, pv, spv);
    return Py_BuildValue("(ddd)(ddd)",
                        spv[0][0], spv[0][1], spv[0][2],
                        spv[1][0], spv[1][1], spv[1][2]);
}

PyDoc_STRVAR(erfa_s2xpv_doc,
"s2xpv(s1, s2, pv) -> spv\n\n"
"Multiply a pv-vector by two scalars.\n"
"Given:\n"
"   s1          scalar to multiply position component by\n"
"   s2          scalar to multiply velocity component by\n"
"   pv          pv-vector\n"
"Returned:\n"
"   spv         pv-vector: p scaled by s1, v scaled by s2");

static PyObject *
erfa_sepp(PyObject *self, PyObject *args)
{
    double a[3], b[3], s;
    double a0, a1, a2, b0, b1, b2;
    if (!PyArg_ParseTuple(args, "(ddd)(ddd)",
        &a0, &a1, &a2, &b0, &b1, &b2)) {
        return NULL;
    }
    a[0] = a0;
    a[1] = a1;
    a[2] = a2;
    b[0] = b0;
    b[1] = b1;
    b[2] = b2;
    s = eraSepp(a, b);
    return Py_BuildValue("d", s);
}

PyDoc_STRVAR(erfa_sepp_doc,
"sepp(a, b) -> s\n\n"
"Angular separation between two p-vectors.\n"
"Given:\n"
"   a       first p-vector (not necessarily unit length)\n"
"   b       second p-vector (not necessarily unit length)\n"
"Returned:\n"
"   s       angular separation (radians, always positive)");

static PyObject *
erfa_seps(PyObject *self, PyObject *args)
{
    double al, ap, bl, bp, s;
    if (!PyArg_ParseTuple(args, "dddd",
        &al, &ap, &bl, &bp)) {
        return NULL;
    }
    s = eraSeps(al, ap, bl, bp);
    return Py_BuildValue("d", s);
}

PyDoc_STRVAR(erfa_seps_doc,
"seps(al, ap, bl, bp) -> s\n\n"
"Angular separation between two sets of spherical coordinates.\n"
"Given:\n"
"   al  first longitude (radians)\n"
"   ap  first latitude (radians)\n"
"   bl  second longitude (radians)\n"
"   bp  second latitude (radians)\n"
"Returned:\n"
"   s   angular separation (radians)");

static PyObject *
erfa_sxp(PyObject *self, PyObject *args)
{
    double s, p[3], sp[3], p0, p1, p2;
    if (!PyArg_ParseTuple(args, "d(ddd)", &s,
                                &p0, &p1, &p2)) {
        return NULL;
    }
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    eraSxp(s, p, sp);
    return Py_BuildValue("ddd", sp[0], sp[1], sp[2]);
}

PyDoc_STRVAR(erfa_sxp_doc,
"sxp(s, p) -> sp\n\n"
"Multiply a p-vector by a scalar.\n"
"Given:\n"
"   s           scalar\n"
"   p           p-vector\n"
"Returned:\n"
"   sp          s * p");

static PyObject *
erfa_sxpv(PyObject *self, PyObject *args)
{
    double s, pv[2][3], spv[2][3];
    double pv00, pv01, pv02, pv10, pv11, pv12;
    if (!PyArg_ParseTuple(args, "d((ddd)(ddd))", &s,
                                &pv00, &pv01, &pv02, &pv10, &pv11, &pv12)) {
        return NULL;
    }
    pv[0][0] = pv00;
    pv[0][1] = pv01;
    pv[0][2] = pv02;
    pv[1][0] = pv10;
    pv[1][1] = pv11;
    pv[1][2] = pv12;
    eraSxpv(s, pv, spv);
    return Py_BuildValue("(ddd)(ddd)",
                        spv[0][0], spv[0][1], spv[0][2],
                        spv[1][0], spv[1][1], spv[1][2]);
}

PyDoc_STRVAR(erfa_sxpv_doc,
"sxpv(s, pv) -> spv\n\n"
"Multiply a pv-vector by a scalar.\n"
"Given:\n"
"   s           scalar\n"
"   pv          pv-vector\n"
"Returned:\n"
"   spv         s * pv");

static PyObject *
erfa_tf2a(PyObject *self, PyObject *args)
{
    int ihour, imin, status;
    double sec, rad;
    int s = '+';
    if (!PyArg_ParseTuple(args, "iid", &ihour, &imin, &sec)) {
        return NULL;
    }
    if (ihour < 0) s = '-';
    status = eraTf2a(s, ihour, imin, sec, &rad);
    if (status == 1) {
        PyErr_WarnEx(PyExc_Warning, "hour outside range 0-23", 1);
    }
    else if (status == 2) {
        PyErr_WarnEx(PyExc_Warning, "min outside range 0-59", 1);
    }
    else if (status == 3) {
        PyErr_WarnEx(PyExc_Warning, "sec outside range 0-59", 1);
    }
    return Py_BuildValue("d", rad);
}

PyDoc_STRVAR(erfa_tf2a_doc,
"tf2a(hour, min, sec) -> rad\n\n"
"Convert hours, minutes, seconds to radians.\n"
"Given:\n"
"   hour    hours\n"
"   min     minutes\n"
"   sec     seconds\n"
"Returned:\n"
"   rad     angle in radians");

static PyObject *
erfa_tf2d(PyObject *self, PyObject *args)
{
    int ihour, imin, status;
    double sec, days;
    int s = '+';
    if (!PyArg_ParseTuple(args, "iid", &ihour, &imin, &sec)) {
        return NULL;
    }
    if (ihour < 0) s = '-';
    status = eraTf2d(s, ihour, imin, sec, &days);
    if (status == 1) {
        PyErr_WarnEx(PyExc_Warning, "hour outside range 0-23", 1);
    }
    else if (status == 2) {
        PyErr_WarnEx(PyExc_Warning, "min outside range 0-59", 1);
    }
    else if (status == 3) {
        PyErr_WarnEx(PyExc_Warning, "sec outside range 0-59", 1);
    }
    return Py_BuildValue("d", days);
}

PyDoc_STRVAR(erfa_tf2d_doc,
"tf2d(hour, min, sec) -> days\n\n"
"Convert hours, minutes, seconds to days.\n"
"Given:\n"
"   hour    hours\n"
"   min     minutes\n"
"   sec     seconds\n"
"Returned:\n"
"   days    interval in days");

static PyObject *
erfa_tr(PyObject *self, PyObject *args)
{
    double r[3][3], rt[3][3];
    double r00,r01,r02,r10,r11,r12,r20,r21,r22;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))", 
                                &r00,&r01,&r02,&r10,&r11,&r12,&r20,&r21,&r22)) {
        return NULL;
    }
    r[0][0] = r00;
    r[0][1] = r01;
    r[0][2] = r02;
    r[1][0] = r10;
    r[1][1] = r11;
    r[1][2] = r12;
    r[2][0] = r20;
    r[2][1] = r21;
    r[2][2] = r22;
    eraTr(r, rt);
    return Py_BuildValue("((ddd)(ddd)(ddd))",
                            rt[0][0],rt[0][1],rt[0][2],
                            rt[1][0],rt[1][1],rt[1][2],
                            rt[2][0],rt[2][1],rt[2][2]);    
}

PyDoc_STRVAR(erfa_tr_doc,
"tr(r) -> rt\n\n"
"Transpose an r-matrix.\n"
"Given:\n"
"   r           r-matrix\n"
"Given and returned:\n"
"   rt          transpose");

static PyObject *
erfa_trxp(PyObject *self, PyObject *args)
{
    double r[3][3], p[3], trp[3];
    double r00,r01,r02,r10,r11,r12,r20,r21,r22,p0,p1,p2;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))(ddd)", 
                                &r00,&r01,&r02,&r10,&r11,&r12,&r20,&r21,&r22,
                                &p0,&p1,&p2)) {
        return NULL;
    }
    r[0][0] = r00;
    r[0][1] = r01;
    r[0][2] = r02;
    r[1][0] = r10;
    r[1][1] = r11;
    r[1][2] = r12;
    r[2][0] = r20;
    r[2][1] = r21;
    r[2][2] = r22;
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    eraTrxp(r, p, trp);
    return Py_BuildValue("ddd", trp[0], trp[1], trp[2]);    
}

PyDoc_STRVAR(erfa_trxp_doc,
"trxp(r, p) -> trp\n\n"
"Multiply a p-vector by the transpose of an r-matrix.\n"
"Given:\n"
"   r           r-matrix\n"
"   p           p-vector\n"
"Given and returned:\n"
"   trp         r * p");

static PyObject *
erfa_trxpv(PyObject *self, PyObject *args)
{
    double r[3][3], pv[2][3], trpv[2][3];
    double r00,r01,r02,r10,r11,r12,r20,r21,r22;
    double pv00, pv01, pv02, pv10, pv11, pv12;
    if (!PyArg_ParseTuple(args, "((ddd)(ddd)(ddd))((ddd)(ddd))", 
                                &r00,&r01,&r02,&r10,&r11,&r12,&r20,&r21,&r22,
                                &pv00, &pv01, &pv02, &pv10, &pv11, &pv12)) {
        return NULL;
    }
    r[0][0] = r00;
    r[0][1] = r01;
    r[0][2] = r02;
    r[1][0] = r10;
    r[1][1] = r11;
    r[1][2] = r12;
    r[2][0] = r20;
    r[2][1] = r21;
    r[2][2] = r22;
    pv[0][0] = pv00;
    pv[0][1] = pv01;
    pv[0][2] = pv02;
    pv[1][0] = pv10;
    pv[1][1] = pv11;
    pv[1][2] = pv12;
    eraTrxpv(r, pv, trpv);
    return Py_BuildValue("(ddd)(ddd)",
                        trpv[0][0], trpv[0][1], trpv[0][2],
                        trpv[1][0], trpv[1][1], trpv[1][2]);
}

PyDoc_STRVAR(erfa_trxpv_doc,
"trxpv(r, pv) -> trpv\n\n"
"Multiply a pv-vector by the transpose of an r-matrix.\n"
"Given:\n"
"   r           r-matrix\n"
"   pv          pv-vector\n"
"Given and returned:\n"
"   trpv        r * pv");


static PyMethodDef erfa_methods[] = {
    {"bi00", (PyCFunction)erfa_bi00, METH_NOARGS, erfa_bi00_doc},
    {"bp00", erfa_bp00, METH_VARARGS, erfa_bp00_doc},
    {"bp06", erfa_bp06, METH_VARARGS, erfa_bp06_doc},
    {"bpn2xy", erfa_bpn2xy, METH_VARARGS, erfa_bpn2xy_doc},
    {"c2i00a", erfa_c2i00a, METH_VARARGS, erfa_c2i00a_doc},
    {"c2i00b", erfa_c2i00b, METH_VARARGS, erfa_c2i00b_doc},
    {"c2i06a", erfa_c2i06a, METH_VARARGS, erfa_c2i06a_doc},
    {"c2ibpn", erfa_c2ibpn, METH_VARARGS, erfa_c2ibpn_doc},
    {"c2ixy", erfa_c2ixy, METH_VARARGS, erfa_c2ixy_doc},
    {"c2ixys", erfa_c2ixys, METH_VARARGS, erfa_c2ixys_doc},
    {"c2t00a", erfa_c2t00a, METH_VARARGS, erfa_c2t00a_doc},
    {"c2t00b", erfa_c2t00b, METH_VARARGS, erfa_c2t00b_doc},
    {"c2t06a", erfa_c2t06a, METH_VARARGS, erfa_c2t06a_doc},
    {"c2tcio", erfa_c2tcio, METH_VARARGS, erfa_c2tcio_doc},
    {"c2teqx", erfa_c2teqx, METH_VARARGS, erfa_c2teqx_doc},
    {"c2tpe", erfa_c2tpe, METH_VARARGS, erfa_c2tpe_doc},
    {"c2txy", erfa_c2txy, METH_VARARGS, erfa_c2txy_doc},
    {"cal2jd", erfa_cal2jd, METH_VARARGS, erfa_cal2jd_doc},
    {"d2dtf", erfa_d2dtf, METH_VARARGS, erfa_d2dtf_doc},
    {"dat", erfa_dat, METH_VARARGS, erfa_dat_doc},
    {"dtdb", erfa_dtdb, METH_VARARGS, erfa_dtdb_doc},
    {"dtf2d", erfa_dtf2d, METH_VARARGS, erfa_dtf2d_doc},
    {"ee00", erfa_ee00, METH_VARARGS, erfa_ee00_doc},
    {"ee00a", erfa_ee00a, METH_VARARGS, erfa_ee00a_doc},
    {"ee00b", erfa_ee00b, METH_VARARGS, erfa_ee00b_doc},
    {"ee06a", erfa_ee06a, METH_VARARGS, erfa_ee06a_doc},
    {"eect00", erfa_eect00, METH_VARARGS, erfa_eect00_doc},
    {"eform", erfa_eform, METH_VARARGS, erfa_eform_doc},
    {"eo06a", erfa_eo06a, METH_VARARGS, erfa_eo06a_doc},
    {"eors", erfa_eors, METH_VARARGS, erfa_eors_doc},
    {"epb", erfa_epb, METH_VARARGS, erfa_epb_doc},
    {"epb2jd", erfa_epb2jd, METH_VARARGS, erfa_epb2jd_doc},
    {"epj", erfa_epj, METH_VARARGS, erfa_epj_doc},
    {"epj2jd", erfa_epj2jd, METH_VARARGS, erfa_epj2jd_doc},
    {"epv00", (PyCFunction)erfa_epv00, METH_VARARGS|METH_KEYWORDS, erfa_epv00_doc},
    {"eqeq94", erfa_eqeq94, METH_VARARGS, erfa_eqeq94_doc},
    {"era00", erfa_era00, METH_VARARGS, erfa_era00_doc},
    {"fad03", erfa_fad03, METH_VARARGS, erfa_fad03_doc},
    {"fae03", erfa_fae03, METH_VARARGS, erfa_fae03_doc},
    {"faf03", erfa_faf03, METH_VARARGS, erfa_faf03_doc},
    {"faju03", erfa_faju03, METH_VARARGS, erfa_faju03_doc},
    {"fal03", erfa_fal03, METH_VARARGS, erfa_fal03_doc},
    {"falp03", erfa_falp03, METH_VARARGS, erfa_falp03_doc},
    {"fama03", erfa_fama03, METH_VARARGS, erfa_fama03_doc},
    {"fame03", erfa_fame03, METH_VARARGS, erfa_fame03_doc},
    {"fane03", erfa_fane03, METH_VARARGS, erfa_fane03_doc},
    {"faom03", erfa_faom03, METH_VARARGS, erfa_faom03_doc},
    {"fapa03", erfa_fapa03, METH_VARARGS, erfa_fapa03_doc},
    {"fasa03", erfa_fasa03, METH_VARARGS, erfa_fasa03_doc},
    {"faur03", erfa_faur03, METH_VARARGS, erfa_faur03_doc},
    {"fave03", erfa_fave03, METH_VARARGS, erfa_fave03_doc},
    {"fk52h", erfa_fk52h, METH_VARARGS, erfa_fk52h_doc},
    {"fk5hip", (PyCFunction)erfa_fk5hip, METH_NOARGS, erfa_fk5hip_doc},
    {"fk5hz", erfa_fk5hz, METH_VARARGS, erfa_fk5hz_doc},
    {"fw2m", erfa_fw2m, METH_VARARGS, erfa_fw2m_doc},
    {"fw2xy", erfa_fw2xy, METH_VARARGS, erfa_fw2xy_doc},
    {"gc2gd", erfa_gc2gd, METH_VARARGS, erfa_gc2gd_doc},
    {"gc2gde", erfa_gc2gde, METH_VARARGS, erfa_gc2gde_doc},
    {"gd2gc", erfa_gd2gc, METH_VARARGS, erfa_gd2gc_doc},
    {"gd2gce", erfa_gd2gce, METH_VARARGS, erfa_gd2gce_doc},
    {"gmst00", erfa_gmst00, METH_VARARGS, erfa_gmst00_doc},
    {"gmst06", erfa_gmst06, METH_VARARGS, erfa_gmst06_doc},
    {"gmst82", erfa_gmst82, METH_VARARGS, erfa_gmst82_doc},
    {"gst00a", erfa_gst00a, METH_VARARGS, erfa_gst00a_doc},
    {"gst00b", erfa_gst00b, METH_VARARGS, erfa_gst00b_doc},
    {"gst06", erfa_gst06, METH_VARARGS, erfa_gst06_doc},
    {"gst06a", erfa_gst06a, METH_VARARGS, erfa_gst06a_doc},
    {"gst94", erfa_gst94, METH_VARARGS, erfa_gst94_doc},
    {"h2fk5", erfa_h2fk5, METH_VARARGS, erfa_h2fk5_doc},
    {"hfk5z", erfa_hfk5z, METH_VARARGS, erfa_hfk5z_doc},
    {"jd2cal", erfa_jd2cal, METH_VARARGS, erfa_jd2cal_doc},
    {"jdcalf", erfa_jdcalf, METH_VARARGS, erfa_jdcalf_doc},
    {"num00a", erfa_num00a, METH_VARARGS, erfa_num00a_doc},
    {"num00b", erfa_num00b, METH_VARARGS, erfa_num00b_doc},
    {"num06a", erfa_num06a, METH_VARARGS, erfa_num06a_doc},
    {"numat", erfa_numat, METH_VARARGS, erfa_numat_doc},
    {"nut00a", erfa_nut00a, METH_VARARGS, erfa_nut00a_doc},
    {"nut00b", erfa_nut00b, METH_VARARGS, erfa_nut00b_doc},
    {"nut06a", erfa_nut06a, METH_VARARGS, erfa_nut06a_doc},
    {"nut80", erfa_nut80, METH_VARARGS, erfa_nut80_doc},
    {"nutm80", erfa_nutm80, METH_VARARGS, erfa_nutm80_doc},
    {"obl06", erfa_obl06, METH_VARARGS, erfa_obl06_doc},
    {"obl80", erfa_obl80, METH_VARARGS, erfa_obl80_doc},
    {"p06e", erfa_p06e, METH_VARARGS, erfa_p06e_doc},
    {"pb06", erfa_pb06, METH_VARARGS, erfa_pb06_doc},
    {"pfw06", erfa_pfw06, METH_VARARGS, erfa_pfw06_doc},
    {"plan94", erfa_plan94, METH_VARARGS, erfa_plan94_doc},
    {"pmat00", erfa_pmat00, METH_VARARGS, erfa_pmat00_doc},
    {"pmat06", erfa_pmat06, METH_VARARGS, erfa_pmat06_doc},
    {"pmat76", erfa_pmat76, METH_VARARGS, erfa_pmat76_doc},
    {"pn00", erfa_pn00, METH_VARARGS, erfa_pn00_doc},
    {"pn00a", erfa_pn00a, METH_VARARGS, erfa_pn00a_doc},
    {"pn00b", erfa_pn00b, METH_VARARGS, erfa_pn00b_doc},
    {"pn06", erfa_pn06, METH_VARARGS, erfa_pn06_doc},
    {"pn06a", erfa_pn06a, METH_VARARGS, erfa_pn06a_doc},
    {"pnm00a", erfa_pnm00a, METH_VARARGS, erfa_pnm00a_doc},
    {"pnm00b", erfa_pnm00b, METH_VARARGS, erfa_pnm00b_doc},
    {"pnm06a", erfa_pnm06a, METH_VARARGS, erfa_pnm06a_doc},
    {"pnm80", erfa_pnm80, METH_VARARGS, erfa_pnm80_doc},
    {"pom00", erfa_pom00, METH_VARARGS, erfa_pom00_doc},
    {"pr00", erfa_pr00, METH_VARARGS, erfa_pr00_doc},
    {"prec76", erfa_prec76, METH_VARARGS, erfa_prec76_doc},
    {"pvstar", erfa_pvstar, METH_VARARGS, erfa_pvstar_doc},
    {"s00", erfa_s00, METH_VARARGS, erfa_s00_doc},
    {"s00a", erfa_s00a, METH_VARARGS, erfa_s00a_doc},
    {"s00b", erfa_s00b, METH_VARARGS, erfa_s00b_doc},
    {"s06", erfa_s06, METH_VARARGS, erfa_s06_doc},
    {"s06a", erfa_s06a, METH_VARARGS, erfa_s06a_doc},
    {"sp00", erfa_sp00, METH_VARARGS, erfa_sp00_doc},
    {"starpm", erfa_starpm, METH_VARARGS, erfa_starpm_doc},
    {"starpv", erfa_starpv, METH_VARARGS, erfa_starpv_doc},
    {"taitt", erfa_taitt, METH_VARARGS, erfa_taitt_doc},
    {"taiut1", erfa_taiut1, METH_VARARGS, erfa_taiut1_doc},
    {"taiutc", erfa_taiutc, METH_VARARGS, erfa_taiutc_doc},
    {"tcbtdb", erfa_tcbtdb, METH_VARARGS, erfa_tcbtdb_doc},
    {"tcgtt", erfa_tcgtt, METH_VARARGS, erfa_tcgtt_doc},
    {"tdbtcb", erfa_tdbtcb, METH_VARARGS, erfa_tdbtcb_doc},
    {"tdbtt", erfa_tdbtt, METH_VARARGS, erfa_tdbtt_doc},
    {"tttai", erfa_tttai, METH_VARARGS, erfa_tttai_doc},
    {"tttcg", erfa_tttcg, METH_VARARGS, erfa_tttcg_doc},
    {"tttdb", erfa_tttdb, METH_VARARGS, erfa_tttdb_doc},
    {"ttut1", erfa_ttut1, METH_VARARGS, erfa_ttut1_doc},
    {"ut1tai", erfa_ut1tai, METH_VARARGS, erfa_ut1tai_doc},
    {"ut1tt", erfa_ut1tt, METH_VARARGS, erfa_ut1tt_doc},
    {"ut1utc", erfa_ut1utc, METH_VARARGS, erfa_ut1utc_doc},
    {"utctai", erfa_utctai, METH_VARARGS, erfa_utctai_doc},
    {"utcut1", erfa_utcut1, METH_VARARGS, erfa_utcut1_doc},
    {"xy06", erfa_xy06, METH_VARARGS, erfa_xy06_doc},
    {"xys00a", erfa_xys00a, METH_VARARGS, erfa_xys00a_doc},
    {"xys00b", erfa_xys00b, METH_VARARGS, erfa_xys00b_doc},
    {"xys06a", erfa_xys06a, METH_VARARGS, erfa_xys06a_doc},
    {"a2af", erfa_a2af, METH_VARARGS, erfa_a2af_doc},
    {"a2tf", erfa_a2tf, METH_VARARGS, erfa_a2tf_doc},
    {"af2a", erfa_af2a, METH_VARARGS, erfa_af2a_doc},
    {"anp", erfa_anp, METH_O, erfa_anp_doc},
    {"anpm", erfa_anpm, METH_O, erfa_anpm_doc},
    {"c2s", erfa_c2s, METH_VARARGS, erfa_c2s_doc},
    {"cp", erfa_cp, METH_VARARGS, erfa_cp_doc},
    {"cpv", erfa_cpv, METH_VARARGS, erfa_cpv_doc},
    {"cr", erfa_cr, METH_VARARGS, erfa_cr_doc},
    {"d2tf", erfa_d2tf, METH_VARARGS, erfa_d2tf_doc},
    {"p2pv", erfa_p2pv, METH_VARARGS, erfa_p2pv_doc},
    {"p2s", erfa_p2s, METH_VARARGS, erfa_p2s_doc},
    {"pap", erfa_pap, METH_VARARGS, erfa_pap_doc},
    {"pas", erfa_pas, METH_VARARGS, erfa_pas_doc},
    {"pdp", erfa_pdp, METH_VARARGS, erfa_pdp_doc},
    {"pm", erfa_pm, METH_VARARGS, erfa_pm_doc},
    {"pmp", erfa_pmp, METH_VARARGS, erfa_pmp_doc},
    {"pn", erfa_pn, METH_VARARGS, erfa_pn_doc},
    {"ppp", erfa_ppp, METH_VARARGS, erfa_ppp_doc},
    {"ppsp", erfa_ppsp, METH_VARARGS, erfa_ppsp_doc},
    {"pv2p", erfa_pv2p, METH_VARARGS, erfa_pv2p_doc},
    {"pv2s", erfa_pv2s, METH_VARARGS, erfa_pv2s_doc},
    {"pvdpv", erfa_pvdpv, METH_VARARGS, erfa_pvdpv_doc},
    {"pvm", erfa_pvm, METH_VARARGS, erfa_pvm_doc},
    {"pvmpv", erfa_pvmpv, METH_VARARGS, erfa_pvmpv_doc},
    {"pvppv", erfa_pvppv, METH_VARARGS, erfa_pvppv_doc},
    {"pvu", erfa_pvu, METH_VARARGS, erfa_pvu_doc},
    {"pvup", erfa_pvup, METH_VARARGS, erfa_pvup_doc},
    {"pvxpv", erfa_pvxpv, METH_VARARGS, erfa_pvxpv_doc},
    {"pxp", erfa_pxp, METH_VARARGS, erfa_pxp_doc},
    {"rm2v", erfa_rm2v, METH_VARARGS, erfa_rm2v_doc},
    {"rv2m", erfa_rv2m, METH_VARARGS, erfa_rv2m_doc},
    {"rx", erfa_rx, METH_VARARGS, erfa_rx_doc},
    {"rxp", erfa_rxp, METH_VARARGS, erfa_rxp_doc},
    {"rxpv", erfa_rxpv, METH_VARARGS, erfa_rxpv_doc},
    {"rxr", erfa_rxr, METH_VARARGS, erfa_rxr_doc},
    {"ry", erfa_ry, METH_VARARGS, erfa_ry_doc},
    {"rz", erfa_rz, METH_VARARGS, erfa_rz_doc},
    {"s2c", erfa_s2c, METH_VARARGS, erfa_s2c_doc},
    {"s2p", erfa_s2p, METH_VARARGS, erfa_s2p_doc},
    {"s2pv", erfa_s2pv, METH_VARARGS, erfa_s2pv_doc},
    {"s2xpv", erfa_s2xpv, METH_VARARGS, erfa_s2xpv_doc},
    {"sepp", erfa_sepp, METH_VARARGS, erfa_sepp_doc},
    {"seps", erfa_seps, METH_VARARGS, erfa_seps_doc},
    {"sxp", erfa_sxp, METH_VARARGS, erfa_sxp_doc},
    {"sxpv", erfa_sxpv, METH_VARARGS, erfa_sxpv_doc},
    {"tf2a", erfa_tf2a, METH_VARARGS, erfa_tf2a_doc},
    {"tf2d", erfa_tf2d, METH_VARARGS, erfa_tf2d_doc},
    {"tr", erfa_tr, METH_VARARGS, erfa_tr_doc},
    {"trxp", erfa_trxp, METH_VARARGS, erfa_trxp_doc},
    {"trxpv", erfa_trxpv, METH_VARARGS, erfa_trxpv_doc},
    {NULL,		NULL}		/* sentinel */
};

PyDoc_STRVAR(module_doc,
"This module provides ERFA,\n\
the Essential Routine for Fundamental Astronomy,\n\
interface to Python\n\
");

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef erfamodule = {
	PyModuleDef_HEAD_INIT,
	"erfa",
	module_doc,
	-1,
	erfa_methods,
	NULL,
	NULL,
	NULL,
	NULL
};

PyMODINIT_FUNC
PyInit_erfa(void)
{
	PyObject *m;
	m = PyModule_Create(&erfamodule);
	if (m == NULL)
            return NULL;
        erfaError = PyErr_NewException("erfa.error", NULL, NULL);
        Py_INCREF(erfaError);
        PyModule_AddObject(m, "error", erfaError);
        PyModule_AddObject(m, "S2R", PyFloat_FromDouble(DS2R));
        PyModule_AddObject(m, "D2R", PyFloat_FromDouble(DD2R));
        PyModule_AddObject(m, "R2AS", PyFloat_FromDouble(DR2AS));
        PyModule_AddObject(m, "AS2R", PyFloat_FromDouble(DAS2R));
        PyModule_AddObject(m, "TURNAS", PyFloat_FromDouble(TURNAS));
        PyModule_AddObject(m, "MAS2R", PyFloat_FromDouble(DMAS2R));
        PyModule_AddObject(m, "DAYSEC", PyFloat_FromDouble(DAYSEC));
        PyModule_AddObject(m, "JY", PyFloat_FromDouble(DJY));
        PyModule_AddObject(m, "JC", PyFloat_FromDouble(DJC));
        PyModule_AddObject(m, "JM", PyFloat_FromDouble(DJM));
        PyModule_AddObject(m, "J00", PyFloat_FromDouble(DJ00));
        PyModule_AddObject(m, "JM0", PyFloat_FromDouble(DJM0));
        PyModule_AddObject(m, "JM00", PyFloat_FromDouble(DJM00));
        PyModule_AddObject(m, "JM77", PyFloat_FromDouble(DJM77));
        PyModule_AddObject(m, "TTMTAI", PyFloat_FromDouble(TTMTAI));
        PyModule_AddObject(m, "AU", PyFloat_FromDouble(DAU));
        PyModule_AddObject(m, "C", PyFloat_FromDouble(DC));
        PyModule_AddObject(m, "ELG", PyFloat_FromDouble(ELG));
        PyModule_AddObject(m, "ELB", PyFloat_FromDouble(ELB));
        PyModule_AddObject(m, "TDB0", PyFloat_FromDouble(TDB0));
        return m;
}
#else
PyMODINIT_FUNC
initerfa(void)
{
	PyObject *m;
	m = Py_InitModule3("erfa", erfa_methods, module_doc);
	if (m == NULL)
            goto finally;
        erfaError = PyErr_NewException("erfa.error", NULL, NULL);
        Py_INCREF(erfaError);
        PyModule_AddObject(m, "error", erfaError);
        PyModule_AddObject(m, "S2R", PyFloat_FromDouble(DS2R));
        PyModule_AddObject(m, "D2R", PyFloat_FromDouble(DD2R));
        PyModule_AddObject(m, "R2AS", PyFloat_FromDouble(DR2AS));
        PyModule_AddObject(m, "AS2R", PyFloat_FromDouble(DAS2R));
        PyModule_AddObject(m, "TURNAS", PyFloat_FromDouble(TURNAS));
        PyModule_AddObject(m, "MAS2R", PyFloat_FromDouble(DMAS2R));
        PyModule_AddObject(m, "DAYSEC", PyFloat_FromDouble(DAYSEC));
        PyModule_AddObject(m, "JY", PyFloat_FromDouble(DJY));
        PyModule_AddObject(m, "JC", PyFloat_FromDouble(DJC));
        PyModule_AddObject(m, "JM", PyFloat_FromDouble(DJM));
        PyModule_AddObject(m, "J00", PyFloat_FromDouble(DJ00));
        PyModule_AddObject(m, "JM0", PyFloat_FromDouble(DJM0));
        PyModule_AddObject(m, "JM00", PyFloat_FromDouble(DJM00));
        PyModule_AddObject(m, "JM77", PyFloat_FromDouble(DJM77));
        PyModule_AddObject(m, "TTMTAI", PyFloat_FromDouble(TTMTAI));
        PyModule_AddObject(m, "AU", PyFloat_FromDouble(DAU));
        PyModule_AddObject(m, "C", PyFloat_FromDouble(DC));
        PyModule_AddObject(m, "ELG", PyFloat_FromDouble(ELG));
        PyModule_AddObject(m, "ELB", PyFloat_FromDouble(ELB));
        PyModule_AddObject(m, "TDB0", PyFloat_FromDouble(TDB0));

        finally:
        return;
 
}
#endif
