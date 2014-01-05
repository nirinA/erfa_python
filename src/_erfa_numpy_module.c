#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "structseq.h"
#include <math.h>
#include <time.h>
#include "erfa.h"
#include "erfam.h"
#include "numpy/npy_3kcompat.h"
#include "numpy/arrayobject.h"

static PyObject *_erfaError;

static int initialized;
/* local prototype */
double *pyvector_to_Carrayptrs(PyArrayObject *arrayin);
double **pymatrix_to_Carrayptrs(PyArrayObject *arrayin);
double **ptrvector(long n);
void free_Carrayptr(double **v);

static PyTypeObject AstromType;

static PyTypeObject LdbodyType;
/* local function */
double **pymatrix_to_Carrayptrs(PyArrayObject *arrayin)
{
    double **c, *a;
    int i, n, m;
    n = arrayin->dimensions[0];
    m = arrayin->dimensions[1];
    c = ptrvector(n);
    /* pointer to arrayin data as double */
    a = (double *)arrayin->data;
    for (i=0;i<n;i++) {
        c[i] = a + i*m;
    }
    return c;
}

double **ptrvector(long n)
{
    double **v;
    v = (double **)malloc((size_t)(n *sizeof(double)));
    if (!v) {
        PyErr_SetString(_erfaError, "malloc failed in **ptrvector");
        exit(0);
    }
    return v;
}

double *pyvector_to_Carrayptrs(PyArrayObject *arrayin)
{
    /* pointer to arrayin data as double */
    return (double *) arrayin->data; 
}

void free_Carrayptr(double **v){
    free((char *)v);
}

static PyObject *
_to_py_vector(double v[3])
{
    double *cv;
    PyArrayObject *pyout;
    PyArray_Descr * dsc;
    dsc = PyArray_DescrFromType(NPY_DOUBLE);
    npy_intp dims[] = {3};
    int j;
    pyout=(PyArrayObject *) PyArray_NewFromDescr(&PyArray_Type,
                                                 dsc,
                                                 1,
                                                 dims,
                                                 NULL,
                                                 NULL,
                                                 0,
                                                 NULL);
    if (NULL == pyout)  return NULL;
    cv = pyvector_to_Carrayptrs(pyout);
    for(j=0;j<3;j++) cv[j] = v[j];
    Py_INCREF(pyout); 
    return PyArray_Return(pyout);                     
}

static PyObject *
_to_py_posvel(double pv[2][3])
{
    double **cpv;
    PyArrayObject *pyout;
    PyArray_Descr * dsc;
    dsc = PyArray_DescrFromType(NPY_DOUBLE);
    npy_intp dims[] = {2,3};
    int i, j;
    pyout=(PyArrayObject *) PyArray_NewFromDescr(&PyArray_Type,
                                                 dsc,
                                                 2,
                                                 dims,
                                                 NULL,
                                                 NULL,
                                                 0,
                                                 NULL);
    if (NULL == pyout)  return NULL;
    cpv = pymatrix_to_Carrayptrs(pyout);
    for (i=0;i<2;i++) {
        for(j=0;j<3;j++) {
            cpv[i][j] = pv[i][j];
        }
    }
    //free?
    Py_INCREF(pyout); 
    return PyArray_Return(pyout);                     
}

static PyObject *
_to_py_matrix(double m[3][3])
{
    double **cm;
    PyArrayObject *pyout;
    PyArray_Descr * dsc;
    dsc = PyArray_DescrFromType(NPY_DOUBLE);
    npy_intp dims[] = {3,3};
    int i, j;
    pyout=(PyArrayObject *) PyArray_NewFromDescr(&PyArray_Type,
                                                 dsc,
                                                 2,
                                                 dims,
                                                 NULL,
                                                 NULL,
                                                 0,
                                                 NULL);
    if (NULL == pyout)  return NULL;
    cm = pymatrix_to_Carrayptrs(pyout);
    for (i=0;i<3;i++) {
        for(j=0;j<3;j++) {
            cm[i][j] = m[i][j];
        }
    }
    //free?
    Py_INCREF(pyout); 
    return PyArray_Return(pyout);                     
}

static PyObject *
_to_py_astrom(eraASTROM *a)
{
    PyObject *v = PyStructSequence_New(&AstromType);
    if (v == NULL)
        return NULL;
#define SET(i,val) PyStructSequence_SET_ITEM(v, i, PyFloat_FromDouble((double) val))
    SET(0, a->pmt);
    PyStructSequence_SET_ITEM(v, 1, _to_py_vector(a->eb));
    PyStructSequence_SET_ITEM(v, 2, _to_py_vector(a->eh));
    SET(3, a->em);
    PyStructSequence_SET_ITEM(v, 4, _to_py_vector(a->v));
    SET(5, a->bm1);
    PyStructSequence_SET_ITEM(v, 6, _to_py_matrix(a->bpn));
    SET(7, a->along);
    SET(8, a->phi);
    SET(9, a->xpl);
    SET(10, a->ypl);
    SET(11, a->sphi);
    SET(12, a->cphi);
    SET(13, a->diurab);
    SET(14, a->eral);
    SET(15, a->refa);
    SET(16, a->refb);
    if (PyErr_Occurred()) {
        Py_XDECREF(v);
        return NULL;
    }
#undef SET
    return v;
}

static PyStructSequence_Field ASTROM_type_fields[] = {
    {"pmt", "PM time interval (SSB, Julian years)"},
    {"eb", "SSB to observer (vector, au)"},
    {"eh", "Sun to observer (unit vector)"},
    {"em", "distance from Sun to observer (au)"},
    {"v", "barycentric observer velocity (vector, c"},
    {"bm1", "sqrt(1-|v|^2): reciprocal of Lorenz factor"},
    {"bpn", "bias-precession-nutation matrix"},
    {"along", "longitude + s' + dERA(DUT) (radians)"},
    {"phi", "geodetic latitude (radians)"},
    {"xpl", "polar motion xp wrt local meridian (radians)"},
    {"ypl", "polar motion yp wrt local meridian (radians)"},
    {"sphi", "sine of geodetic latitude"},
    {"cphi", "cosine of geodetic latitude"},
    {"diurab", "magnitude of diurnal aberration vector"},
    {"eral", "local Earth rotation angle (radians)"},
    {"refa", "refraction constant A (radians)"},
    {"refb", "refraction constant B (radians)"},
    {0}
};

static PyStructSequence_Desc ASTROM_type_desc = {
    "_erfa.ASTROM",
    "Star-independent astrometry parameters\n"
"(Vectors eb, eh, em and v are all with respect to BCRS axes.)",
    ASTROM_type_fields,
    17,
};

static PyStructSequence_Field LDBODY_type_fields[] = {
    {"bm", "mass of the body (solar masses)"},
    {"dl", "deflection limiter (radians^2/2)"},
    {"pv", "barycentric PV of the body (au, au/day)"},
    {0}
};

static PyStructSequence_Desc LDBODY_type_desc = {
    "_erfa.LDBODY",
    "Body parameters for light deflection",
    LDBODY_type_fields,
    3,
};
static PyObject *
_erfa_ab(PyObject *self, PyObject *args)
{
    double pnat[3], v[3], s, bm1, ppr[3];
    double *cpnat, *cv, *cppr;
    PyArrayObject *pypnat, *pyv, *pyout;
    PyArray_Descr * dsc;
    dsc = PyArray_DescrFromType(NPY_DOUBLE);
    npy_intp dims[] = {3};
    int j;
    if (!PyArg_ParseTuple(args, "O!O!dd",
                                 &PyArray_Type, &pypnat,
                                 &PyArray_Type, &pyv,
                                 &s,&bm1))
        return NULL;
    if ((NULL == pypnat) || (NULL == pyv)) return NULL;
    pyout=(PyArrayObject *) PyArray_NewFromDescr(&PyArray_Type,
                                                 dsc,
                                                 1,
                                                 dims,
                                                 NULL,
                                                 NULL,
                                                 0,
                                                 NULL);
    if (NULL == pyout)  return NULL;
    cpnat=pyvector_to_Carrayptrs(pypnat);
    cv=pyvector_to_Carrayptrs(pyv);
    cppr=pyvector_to_Carrayptrs(pyout);
    for(j=0;j<3;j++) {
        pnat[j] = cpnat[j];
        v[j] = cv[j];
    }
    eraAb(pnat, v, s, bm1, ppr);
    for(j=0;j<3;j++)   cppr[j] = ppr[j];
    Py_INCREF(pyout); 
    return PyArray_Return(pyout);                     
}

PyDoc_STRVAR(_erfa_ab_doc,
"\nab(pnat[3], v[3], s, bm1) -> ppr\n"
"Apply aberration to transform natural direction into proper direction.\n"
"Given:\n"
"    pnat       natural direction to the source (unit vector)\n"
"    v          observer barycentric velocity in units of c\n"
"    s          distance between the Sun and the observer (au)\n"
"    bm1        sqrt(1-|v|^2): reciprocal of Lorenz factor\n"
"Returned:\n"
"    ppr        proper direction to source (unit vector)");

static PyObject *
_erfa_apcs(PyObject *self, PyObject *args)
{
    double date1, date2, pv[2][3], ebpv[2][3], ehp[3];
    double **cpv, **cebpv, *cehp;
    PyArrayObject *pypv, *pyebpv, *pyehp;
    eraASTROM astrom;
    int i, j;
    if (!PyArg_ParseTuple(args, "ddO!O!O!",
                          &date1, &date2,
                          &PyArray_Type, &pypv,
                          &PyArray_Type, &pyebpv,
                          &PyArray_Type, &pyehp))      
        return NULL;
    if ((NULL == pypv) || (NULL == pyebpv) || (NULL == pyehp)) return NULL;
    cehp = pyvector_to_Carrayptrs(pyehp);
    cebpv = pymatrix_to_Carrayptrs(pyebpv);
    cpv = pymatrix_to_Carrayptrs(pypv);
    for (i=0;i<3;i++) ehp[i] = cehp[i];
    for (i=0;i<2;i++) {
        for (j=0;j<3;j++) {
            ebpv[i][j] = cebpv[i][j];
            pv[i][j] = cpv[i][j];
        }
    }
    eraApcs(date1, date2, pv, ebpv, ehp, &astrom);
    free_Carrayptr(cebpv);
    free_Carrayptr(cpv);
    return _to_py_astrom(&astrom);
}

PyDoc_STRVAR(_erfa_apcs_doc,
"\napcs(date1, date2, pv[2][3], ebpv[2][3], ehp[3]) -> astrom\n"
"For an observer whose geocentric position and velocity are known,\n"
"prepare star-independent astrometry parameters for transformations\n"
"between ICRS and GCRS.  The Earth ephemeris is supplied by the caller.\n"
"\n"
"The parameters produced by this function are required in the\n"
"parallax, light deflection and aberration parts of the astrometric\n"
"transformation chain.\n"
"Given:\n"
"    date1      TDB as a 2-part...\n"
"    date2      ...Julian Date\n"
"    pv         observer's geocentric pos/vel (m, m/s)\n"
"    ebpv       Earth barycentric pos/vel (au, au/day)\n"
"    ehp        Earth heliocentric position (au)\n"
"Returned:\n"
"    astrom     star-independent astrometry parameters");

static PyObject *
_erfa_aticq(PyObject *self, PyObject *args)
{
    double rc, dc, ri, di;
    double pmt, eb[3], eh[3], em, v[3], bm1, bpn[3][3];
    double along, phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb;
    eraASTROM astrom;
    if (!PyArg_ParseTuple(args, "ddd(ddd)(ddd)d(ddd)d((ddd)(ddd)(ddd))dddddddddd",
                          &ri, &di,
                          &pmt, &eb[0], &eb[1], &eb[2],
                          &eh[0], &eh[1], &eh[2], &em,
                          &v[0], &v[1], &v[2], &bm1,
                          &bpn[0][0],&bpn[0][1],&bpn[0][2],
                          &bpn[1][0],&bpn[1][1],&bpn[1][2],      
                          &bpn[2][0],&bpn[2][1],&bpn[2][2],
                          &along, &phi, &xpl, &ypl, &sphi, &cphi, &diurab,
                          &eral, &refa, &refb))      
        return NULL;
    astrom.pmt = pmt;
    astrom.eb[0] = eb[0];
    astrom.eb[1] = eb[1];
    astrom.eb[2] = eb[2];
    astrom.eh[0] = eh[0];
    astrom.eh[1] = eh[1];
    astrom.eh[2] = eh[2];
    astrom.em = em;
    astrom.v[0] = v[0];
    astrom.v[1] = v[1];
    astrom.v[2] = v[2];
    astrom.bm1 = bm1;
    astrom.bpn[0][0] = bpn[0][0];
    astrom.bpn[0][1] = bpn[0][1];
    astrom.bpn[0][2] = bpn[0][2];
    astrom.bpn[1][0] = bpn[1][0];
    astrom.bpn[1][1] = bpn[1][1];
    astrom.bpn[1][2] = bpn[1][2];
    astrom.bpn[2][0] = bpn[2][0];
    astrom.bpn[2][1] = bpn[2][1];
    astrom.bpn[2][2] = bpn[2][2];
    astrom.along = along;
    astrom.phi = phi;
    astrom.xpl = xpl;
    astrom.ypl = ypl;
    astrom.sphi = sphi;
    astrom.cphi = cphi;
    astrom.diurab = diurab;
    astrom.eral = eral;
    astrom.refa = refa;
    astrom.refb = refb;
    eraAticq(ri, di, &astrom, &rc, &dc);
    return Py_BuildValue("dd", rc, dc);
}

PyDoc_STRVAR(_erfa_aticq_doc,
"\naticq(ri,di, astrom) -> rc, dc\n"
"Quick CIRS RA,Dec to ICRS astrometric place, given the star-\n"
"independent astrometry parameters.\n"
"\n"
"Use of this function is appropriate when efficiency is important and\n"
"where many star positions are all to be transformed for one date.  The\n"
"star-independent parameters can be obtained by calling one of the\n"
"functions apci[13], apcg[13], apco[13] or apcs[13].\n"
"Given:\n"
"    ri,di      CIRS RA,Dec (radians)\n"
"    astrom     star-independent astrometry parameters\n"
"Returned:\n"
"    rc,dc      ICRS astrometric RA,Dec (radians)\n");

static PyObject *
_erfa_plan94(PyObject *self, PyObject *args)
{
    double d1, d2, pv[2][3];
    int np, status;
    if (!PyArg_ParseTuple(args, "ddi", &d1, &d2, &np)) {
        return NULL;
    }
    status = eraPlan94(d1, d2, np, pv);
    if (status == -1){
        PyErr_SetString(_erfaError, "illegal np,  not in range(1,8) for planet");
        return NULL;
    }
    else if (status == 1){
        PyErr_WarnEx(PyExc_Warning, "year outside range(1000:3000)", 1);
    }
    else if (status == 2){
        PyErr_WarnEx(PyExc_Warning,  "computation failed to converge", 1);
    }
    return _to_py_posvel(pv);
}

PyDoc_STRVAR(_erfa_plan94_doc,
"\nplan94(d1, d2, np) -> pv\n\n"
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


static PyMethodDef _erfa_methods[] = {
    {"ab", _erfa_ab, METH_VARARGS, _erfa_ab_doc},
    {"apcs", _erfa_apcs, METH_VARARGS, _erfa_apcs_doc},
    {"aticq", _erfa_aticq, METH_VARARGS, _erfa_aticq_doc},
    {"plan94", _erfa_plan94, METH_VARARGS, _erfa_plan94_doc},
    {NULL,		NULL}		/* sentinel */
};

PyDoc_STRVAR(module_doc,
"This module provides ERFA,\n\
the Essential Routine for Fundamental Astronomy,\n\
interface to Python\n\
");

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef _erfamodule = {
	PyModuleDef_HEAD_INIT,
	"_erfa",
	module_doc,
	-1,
	_erfa_methods,
	NULL,
	NULL,
	NULL,
	NULL
};

PyMODINIT_FUNC
PyInit__erfa(void)
{
	PyObject *m;
	m = PyModule_Create(&_erfamodule);
	import_array();
	if (m == NULL)
            return NULL;
#else
PyMODINIT_FUNC
init_erfa(void)
{
	PyObject *m;
	m = Py_InitModule3("_erfa", _erfa_methods, module_doc);
	import_array();
	if (m == NULL)
            goto finally;
#endif
        _erfaError = PyErr_NewException("_erfa.error", NULL, NULL);
        Py_INCREF(_erfaError);
        PyModule_AddObject(m, "error", _erfaError);
        PyModule_AddObject(m, "DPI", PyFloat_FromDouble(ERFA_DPI));
        PyModule_AddObject(m, "D2PI", PyFloat_FromDouble(ERFA_D2PI));
        PyModule_AddObject(m, "DR2D", PyFloat_FromDouble(ERFA_DR2D));
        PyModule_AddObject(m, "DD2R", PyFloat_FromDouble(ERFA_DD2R));
        PyModule_AddObject(m, "DR2AS", PyFloat_FromDouble(ERFA_DR2AS));
        PyModule_AddObject(m, "DAS2R", PyFloat_FromDouble(ERFA_DAS2R));
        PyModule_AddObject(m, "DS2R", PyFloat_FromDouble(ERFA_DS2R));
        PyModule_AddObject(m, "TURNAS", PyFloat_FromDouble(ERFA_TURNAS));
        PyModule_AddObject(m, "DMAS2R", PyFloat_FromDouble(ERFA_DMAS2R));
        PyModule_AddObject(m, "DTY", PyFloat_FromDouble(ERFA_DTY));
        PyModule_AddObject(m, "DAYSEC", PyFloat_FromDouble(ERFA_DAYSEC));
        PyModule_AddObject(m, "DJY", PyFloat_FromDouble(ERFA_DJY));
        PyModule_AddObject(m, "DJC", PyFloat_FromDouble(ERFA_DJC));
        PyModule_AddObject(m, "DJM", PyFloat_FromDouble(ERFA_DJM));
        PyModule_AddObject(m, "DJ00", PyFloat_FromDouble(ERFA_DJ00));
        PyModule_AddObject(m, "DJM0", PyFloat_FromDouble(ERFA_DJM0));
        PyModule_AddObject(m, "DJM00", PyFloat_FromDouble(ERFA_DJM00));
        PyModule_AddObject(m, "DJM77", PyFloat_FromDouble(ERFA_DJM77));
        PyModule_AddObject(m, "TTMTAI", PyFloat_FromDouble(ERFA_TTMTAI));
        PyModule_AddObject(m, "DAU", PyFloat_FromDouble(ERFA_DAU));
        PyModule_AddObject(m, "CMPS", PyFloat_FromDouble(ERFA_CMPS));
        PyModule_AddObject(m, "AULT", PyFloat_FromDouble(ERFA_AULT));
        PyModule_AddObject(m, "DC", PyFloat_FromDouble(ERFA_DC));
        PyModule_AddObject(m, "ELG", PyFloat_FromDouble(ERFA_ELG));
        PyModule_AddObject(m, "ELB", PyFloat_FromDouble(ERFA_ELB));
        PyModule_AddObject(m, "TDB0", PyFloat_FromDouble(ERFA_TDB0));
        PyModule_AddObject(m, "SRS", PyFloat_FromDouble(ERFA_SRS));
        PyModule_AddObject(m, "WGS84", PyLong_FromLong(ERFA_WGS84));
        PyModule_AddObject(m, "GRS80", PyLong_FromLong(ERFA_GRS80));
        PyModule_AddObject(m, "WGS72", PyLong_FromLong(ERFA_WGS72));

        if (!initialized) {
            PyStructSequence_InitType(&AstromType, &ASTROM_type_desc);
            PyStructSequence_InitType(&LdbodyType, &LDBODY_type_desc);

        }
        Py_INCREF(&AstromType);
        PyModule_AddObject(m, "ASTROM", (PyObject*) &AstromType);
        Py_INCREF(&LdbodyType);
        PyModule_AddObject(m, "LDBODY", (PyObject*) &LdbodyType);

        initialized = 1;

#if PY_VERSION_HEX >= 0x03000000
        return m;
#else
        finally:
        return;
#endif
}
