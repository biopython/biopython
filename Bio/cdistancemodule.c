/* Copyright 2002 by Jeffrey Chang.  All rights reserved.
 * This code is part of the Biopython distribution and governed by its
 * license.  Please see the LICENSE file that should have been included
 * as part of this package.
 *
 * cdistancemodule.c
 * Created 12 Jan 2002
 *
 * Optimized C routines that complement distance.py.
 */

#include "Python.h"
#include <math.h>
#include "csupport.h"



/* Functions in this module. */

static char cdistance_euclidean__doc__[] = 
"euclidean(x, y) -> euclidean distance between x and y";

static PyObject *
cdistance_euclidean(PyObject *self, PyObject *args)
{
    int i;
    int length;
    PyObject *py_x, *py_y;
    PyObject *py_fastx, *py_fasty;
    PyObject *py_xi, *py_yi;
    double xi, yi;
    double sum;

    /* This could be optimized a lot more if I operated on the Numeric
     * Array directly.  However, that would seriously compromise the
     * module's portability...
     */
    if(!PyArg_ParseTuple(args, "OO", &py_x, &py_y))
	return NULL;
    length = PySequence_Size(py_x);
    if(length < 0)
	return NULL;
    if(length != PySequence_Size(py_y)) {
	PyErr_SetString(PyExc_ValueError, "vectors must be same length");
	return NULL;
    }
    if(!(py_fastx = PySequence_Fast(py_x, "could not get fast sequence")))
	return NULL;
    if(!(py_fasty = PySequence_Fast(py_y, "could not get fast sequence"))) {
	Py_DECREF(py_fastx);
	return NULL;
    }

    sum = 0.0;
    for(i=0; i<length; i++) {
	py_xi = PySequence_Fast_GET_ITEM(py_fastx, i);
	py_yi = PySequence_Fast_GET_ITEM(py_fasty, i);
	xi = PyNumber_AsDouble(py_xi);
	yi = PyNumber_AsDouble(py_yi);
	if(PyErr_Occurred())
	    break;
	sum += pow(xi-yi, 2.0);
    }
    Py_DECREF(py_fastx);
    Py_DECREF(py_fasty);
    if(PyErr_Occurred())
	return NULL;
    return PyFloat_FromDouble(sqrt(sum));
}



/* Module definition stuff */

static PyMethodDef cdistanceMethods[] = {
    {"euclidean", cdistance_euclidean, METH_VARARGS, 
     cdistance_euclidean__doc__},
    {NULL, NULL}
};

static char cdistance__doc__[] =
"XXX document here\n\
\n\
";

void initcdistance(void)
{
    (void) Py_InitModule3("cdistance", cdistanceMethods, cdistance__doc__);
}
