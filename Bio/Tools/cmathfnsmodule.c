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


/* Return a PyNumber as a double.
 * Raises a TypeError if I can't do it.
 */
static double PyNumber_AsDouble(PyObject *py_num)
{
    double val;
    PyObject *floatobj;

    if(!PyNumber_Check(py_num)) {
	PyErr_SetString(PyExc_TypeError, "I received a non-number");
	return 0.0;
    }
    if((floatobj = PyNumber_Float(py_num)) == NULL)
	return 0.0;
    val = PyFloat_AsDouble(floatobj);
    Py_DECREF(floatobj);
    return val;
}


/************************************** Exported Functions ***********/

static char cmathfns_safe_log__doc__[] = 
"safe_log(n, zero=None, neg=None) -> log(n)\n\
\n\
Calculate the log of n.  If n is 0, returns the value of zero.  If n is\n\
negative, returns the value of neg.\n\
\n\
";

static PyObject *
cmathfns_safe_log(self, args, keywds)
     PyObject *self;
     PyObject *args;
     PyObject *keywds;
{
    int i;
    PyObject *nobj,
	*zero = Py_None,
	*neg = Py_None;
    double n, logn;

    static char *kwlist[] = {"n", "zero", "neg", NULL};

    if(!PyArg_ParseTupleAndKeywords(args, keywds, "O|OO", kwlist, 
				    &nobj, &zero, &neg))
	return NULL;
    
    n = PyNumber_AsDouble(nobj);
    if(PyErr_Occurred()) {
	return NULL;
    }

    if(n < 0) {
	Py_INCREF(neg);
	return neg;
    } else if(n < 1E-100) {
	Py_INCREF(zero);
	return zero;
    }

    logn = log(n);
    return PyFloat_FromDouble(logn);
}




/************************************** Module definition stuff ******/

static PyMethodDef cmathfnsMethods[] = {
    {"safe_log", (PyCFunction)cmathfns_safe_log, METH_VARARGS|METH_KEYWORDS, 
     cmathfns_safe_log__doc__},
    {NULL, NULL}
};

static char cmathfns__doc__[] =
"This provides helper functions for the mathfns module.\n\
You should never import this module on its own.\n\
\n\
";

void initcmathfns()
{
    (void) Py_InitModule3("cmathfns", cmathfnsMethods, cmathfns__doc__);
}
