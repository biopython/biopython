/* Copyright 2000 by Jeffrey Chang.  All rights reserved.
 * This code is part of the Biopython distribution and governed by its
 * license.  Please see the LICENSE file that should have been included
 * as part of this package.
 *
 * cmathfnsmodule.c
 * Created 3 Jun 2000
 */

#include "Python.h"
#include <math.h>

#include "csupport.h"



/************************************** Exported Functions ***********/

static char cmathfns_intd__doc__[] = 
"intd(x[, digits_after_decimal]) -> int x, rounded\n\
\n\
Represent a floating point number with some digits after the\n\
decimal point as an integer.  This is useful when floating point\n\
comparisons are failing due to precision problems.  e.g.\n\
intd(5.35, 1) -> 54.\n\
\n\
";

static PyObject *cmathfns_intd(
    PyObject *self, PyObject *args, PyObject *keywds)
{
    PyObject *digits_after_decimal = Py_None;
    double x, digits;
    double precision;

    static char *kwlist[] = {"x", "digits_after_decimal", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, keywds, "d|O", kwlist, 
				    &x, &digits_after_decimal))
	return NULL;

    if(digits_after_decimal == Py_None)
	digits = 0;
    else {
	digits = PyNumber_AsDouble(digits_after_decimal);
	if(PyErr_Occurred()) {
	    return NULL;
	}
    }
    precision = pow(10, digits);
    if(x >= 0)
	x = (int)(x * precision + 0.5);
    else
	x = (int)(x * precision - 0.5);
    return PyFloat_FromDouble(x);
}




static char cmathfns_fcmp__doc__[] = 
"fcmp(x, y, precision) -> -1, 0, or 1";

static PyObject *cmathfns_fcmp(
    PyObject *self, PyObject *args, PyObject *keywds)
{
    double x, y, precision;
    int result;

    static char *kwlist[] = {"x", "y", "precision", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, keywds, "ddd", kwlist, 
				    &x, &y, &precision))
	return NULL;

    if(fabs(x-y) < precision)
	result = 0;
    else if(x < y)
	result = -1;
    else result = 1;
    return PyInt_FromLong(result);
}



static char cmathfns_safe_log__doc__[] = 
"safe_log(n, zero=None, neg=None) -> log(n)\n\
\n\
Calculate the log of n.  If n is 0, returns the value of zero.  If n is\n\
negative, returns the value of neg.\n\
\n\
";

static PyObject *cmathfns_safe_log(
    PyObject *self, PyObject *args, PyObject *keywds)
{
    PyObject *zero = Py_None,
	*neg = Py_None;
    double n;

    static char *kwlist[] = {"n", "zero", "neg", NULL};

    if(!PyArg_ParseTupleAndKeywords(args, keywds, "d|OO", kwlist, 
				    &n, &zero, &neg))
	return NULL;
    
    if(n < 0) {
	Py_INCREF(neg);
	return neg;
    } else if(n < 1E-100) {
	Py_INCREF(zero);
	return zero;
    }

    return PyFloat_FromDouble(log(n));
}




/************************************** Module definition stuff ******/

static PyMethodDef cmathfnsMethods[] = {
    {"fcmp", (PyCFunction)cmathfns_fcmp, METH_VARARGS|METH_KEYWORDS, 
     cmathfns_fcmp__doc__},
    {"intd", (PyCFunction)cmathfns_intd, METH_VARARGS|METH_KEYWORDS, 
     cmathfns_intd__doc__},
    {"safe_log", (PyCFunction)cmathfns_safe_log, METH_VARARGS|METH_KEYWORDS, 
     cmathfns_safe_log__doc__},
    {NULL, NULL}
};

static char cmathfns__doc__[] =
"This provides helper functions for the mathfns module.\n\
You should never import this module on its own.\n\
\n\
";

void initcmathfns(void)
{
    (void) Py_InitModule3("cmathfns", cmathfnsMethods, cmathfns__doc__);
}
