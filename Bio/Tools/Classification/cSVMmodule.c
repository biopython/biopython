/* Copyright 2000 by Jeffrey Chang.  All rights reserved.
 * This code is part of the Biopython distribution and governed by its
 * license.  Please see the LICENSE file that should have been included
 * as part of this package.
 *
 * cSVMmodule.c
 * Created 9 Apr 2000
 *
 * Optimized C routines that complement SVM.py.
 */

#include "Python.h"


/* Return a PyNumber as a double.
 * Raises a TypeError if I can't do it.
 */
static double PyNumber_AsDouble(PyObject *py_num)
{
    double val;
    PyObject *floatobj;

    if(!PyNumber_Check(py_num)) {
	PyErr_SetString(PyExc_TypeError, "I received a non-number");
	return(0.0);
    }
    if((floatobj = PyNumber_Float(py_num)) == NULL)
	return(0.0);
    val = PyFloat_AsDouble(floatobj);
    Py_DECREF(floatobj);
    return val;
}



/* Functions in this module. */

static char cSVM_classify__doc__[] = 
"classify(svm, x) -> num";

static PyObject *
cSVM_classify(self, args)
     PyObject *self;
     PyObject *args;
{
    int i;
    PyObject *svm, *x;
    double sum;
    PyObject *ys, *alphas, *xs, *kernel_fn, *b;
    PyObject *ys_i, *alphas_i, *xs_i;
    PyObject *arglist, *retval;
    double y_value, alpha_value, kernel_value, b_value;
    
    if(!PyArg_ParseTuple(args, "OO", &svm, &x))
	return NULL;
    
    xs = PyObject_GetAttrString(svm, "xs");
    if(xs == NULL) {
	return NULL;
    }
    ys = PyObject_GetAttrString(svm, "ys");
    if(ys == NULL) {
	Py_DECREF(xs);
	return NULL;
    }
    alphas = PyObject_GetAttrString(svm, "alphas");
    if(alphas == NULL) {
	Py_DECREF(xs); Py_DECREF(ys);
	return NULL;
    }
    kernel_fn = PyObject_GetAttrString(svm, "kernel_fn");
    if(kernel_fn == NULL) {
	Py_DECREF(xs); Py_DECREF(ys); Py_DECREF(alphas);
	return NULL;
    }
    if(!PyCallable_Check(kernel_fn)) {
	PyErr_SetString(PyExc_TypeError, "kernel_fn is not callable");
	Py_DECREF(xs); Py_DECREF(ys); Py_DECREF(alphas); Py_DECREF(kernel_fn);
	return NULL;
    }
    b = PyObject_GetAttrString(svm, "b");
    if(b == NULL) {
	Py_DECREF(xs); Py_DECREF(ys); Py_DECREF(alphas); Py_DECREF(kernel_fn);
	return NULL;
    }
    
    sum = 0.0;
    for(i=0; ; i++) {
	xs_i = PySequence_GetItem(xs, i);
	if(xs_i == NULL) {
	    PyErr_Clear();
	    break;
	}
	
	alphas_i = PySequence_GetItem(alphas, i);
	if(alphas_i == NULL)
	    break;
	alpha_value = PyNumber_AsDouble(alphas_i);
	if(PyErr_Occurred())
	    break;
	if((alpha_value < 1E-5) && (alpha_value > -1E-5))
	    continue;
	
	ys_i = PySequence_GetItem(ys, i);
	if(ys_i == NULL)
	    break;
	y_value = PyNumber_AsDouble(ys_i);
	if(PyErr_Occurred())
	    break;
	
	arglist = Py_BuildValue("(OO)", xs_i, x);
	retval = PyObject_CallObject(kernel_fn, arglist);
	Py_DECREF(arglist);
	if(retval == NULL)
	    break;
	kernel_value = PyNumber_AsDouble(retval);
	Py_DECREF(retval);
	if(PyErr_Occurred())
	    break;
	sum = sum + y_value*alpha_value*kernel_value;
    }
    
    if(!PyErr_Occurred()) {
	b_value = PyNumber_AsDouble(b);
	if(!PyErr_Occurred())
	    sum = sum - b_value;
    }
    
    Py_DECREF(xs);
    Py_DECREF(ys);
    Py_DECREF(alphas);
    Py_DECREF(kernel_fn);
    Py_DECREF(b);
    
    if(PyErr_Occurred())
	return NULL;
    
    return PyFloat_FromDouble(sum);
}



static char cSVM__sparse_dot__doc__[] = 
"_sparse_dot(x, y) -> num";

static PyObject *
cSVM__sparse_dot(self, args)
     PyObject *self;
     PyObject *args;
{
    int i;
    int size;
    PyObject *x, *y;
    PyObject *xkeys, *key, *xvalue, *yvalue;
    double sum, xdouble, ydouble;

    if(!PyArg_ParseTuple(args, "O!O!", &PyDict_Type, &x, &PyDict_Type, &y))
	return NULL;

    sum = 0.0;
    xkeys = PyDict_Keys(x);
    if(xkeys == NULL) {
	PyErr_SetString(PyExc_SystemError, "Could not get keys for x");
	return NULL;
    }
    size = PyList_Size(xkeys);
    for(i=0; i<size; i++) {
	key = PyList_GetItem(xkeys, i);
	if(key == NULL) {
	    /* Checking for safety.  This will probably only occur if someone
	       were doing something nasty with threads. */
	    PyErr_SetString(PyExc_SystemError, "x dictionary has shrunk");
	    break;
	}
	
	yvalue = PyDict_GetItem(y, key);
	if(yvalue == NULL)
	    continue;
	xvalue = PyDict_GetItem(x, key);
	
	xdouble = PyNumber_AsDouble(xvalue);
	if(PyErr_Occurred())
	    break;
	ydouble = PyNumber_AsDouble(yvalue);
	if(PyErr_Occurred())
	    break;
	
	sum = sum + xdouble * ydouble;
    }
    
    Py_DECREF(xkeys);
    
    if(PyErr_Occurred()) {
	return NULL;
    }
    return PyFloat_FromDouble(sum);
}



static char cSVM__dot__doc__[] = 
"_dot(x, y) -> num";

static PyObject *
cSVM__dot(self, args)
     PyObject *self;
     PyObject *args;
{
    int i;
    double sum;
    PyObject *x, *y, *xobj, *yobj;

    if(!PyArg_ParseTuple(args, "OO", &x, &y))
	return NULL;

    if(!PySequence_Check(x) || !PySequence_Check(y)) {
	PyErr_SetString(PyExc_ValueError, "x and y should be sequences");
	return NULL;
    }
    if(PySequence_Length(x) != PySequence_Length(y)) {
	PyErr_SetString(PyExc_ValueError, "x and y should be the same length");
	return NULL;
    }

    sum = 0.0;
    for(i=0; ; i++) {
	xobj = PySequence_GetItem(x, i);
	if(!xobj) {
	    PyErr_Clear();
	    break;
	}
	yobj = PySequence_GetItem(y, i);

	sum = sum + PyNumber_AsDouble(xobj) * PyNumber_AsDouble(yobj);
	Py_DECREF(xobj);
	Py_DECREF(yobj);
	if(PyErr_Occurred()) { /* Check if PyNumber_AsDouble failed */
	    break;
	}
    }

    if(PyErr_Occurred()) {
	return NULL;
    }
    return PyFloat_FromDouble(sum);
}




/* Module definition stuff */

static PyMethodDef cSVMMethods[] = {
    {"_sparse_dot", cSVM__sparse_dot, METH_VARARGS, cSVM__sparse_dot__doc__},
    {"_dot", cSVM__dot, METH_VARARGS, cSVM__dot__doc__},
    {"classify", cSVM_classify, METH_VARARGS, cSVM_classify__doc__},
    {NULL, NULL}
};

static char cSVM__doc__[] =
"XXX document here\n\
\n\
";

void initcSVM()
{
    (void) Py_InitModule3("cSVM", cSVMMethods, cSVM__doc__);
}
