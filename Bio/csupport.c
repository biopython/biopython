/* Copyright 2002 by Jeffrey Chang.  All rights reserved.
 * This code is part of the Biopython distribution and governed by its
 * license.  Please see the LICENSE file that should have been included
 * as part of this package.
 *
 * csupport.c
 * Created 27 January 2002
 *
 * Miscellaneous useful C functions not to be exported as a python
 * module.
 *
 */

#include "Python.h"


/* Return a PyNumber as a double.
 * Raises a TypeError if I can't do it.
 */
double PyNumber_AsDouble(PyObject *py_num)
{
    double val;
    PyObject *floatobj;

    if((floatobj = PyNumber_Float(py_num)) == NULL)
	return(0.0);
    val = PyFloat_AsDouble(floatobj);
    Py_DECREF(floatobj);
    return val;
}
