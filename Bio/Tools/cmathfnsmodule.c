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
