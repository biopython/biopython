/* Copyright 2001 by Jeffrey Chang.  All rights reserved.
 * This code is part of the Biopython distribution and governed by its
 * license.  Please see the LICENSE file that should have been included
 * as part of this package.
 *
 * csupportmodule.c
 * Created 30 Sep 2001
 *
 * Optimized C routines that complement support.py.
 */

#include "Python.h"


#define _PRECISION 1000

/* Functions in this module. */

static PyObject *
csupport_rint(self, args, keywds)
     PyObject *self;
     PyObject *args;
     PyObject *keywds;
{
    double x;
    int precision = _PRECISION;
    int rint_x;

    static char *kwlist[] = {"precision", NULL};

    if(!PyArg_ParseTupleAndKeywords(args, keywds, "d|l", kwlist,
				    &x, &precision))
	return NULL;
    rint_x = (int)(x * precision + 0.5);
    return PyInt_FromLong(rint_x);
}
/* Module definition stuff */

static PyMethodDef csupportMethods[] = {
    {"rint", csupport_rint, METH_VARARGS|METH_KEYWORDS},
    {NULL, NULL}
};

static char csupport__doc__[] =
"XXX document here\n\
\n\
";

void initcsupport()
{
    (void) Py_InitModule3("csupport", csupportMethods, csupport__doc__);
}
