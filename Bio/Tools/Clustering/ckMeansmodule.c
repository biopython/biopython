/* Copyright 2002 by Jeffrey Chang.  All rights reserved.
 * This code is part of the Biopython distribution and governed by its
 * license.  Please see the LICENSE file that should have been included
 * as part of this package.
 *
 * ckMeansmodule.c
 * Created 12 Jan 2002
 *
 * Optimized C routines that complement kKmeans.py.
 */

#include "Python.h"
#include "csupport.h"

/* Functions in this module. */

static char ckMeans__find_closest_centroid__doc__[] = 
"_find_closest_centroid(vector, centroids, distance_fn) -> index of closest centroid";

static PyObject *
ckMeans__find_closest_centroid(self, args)
     PyObject *self;
     PyObject *args;
{
    int i;
    int length;
    PyObject *py_vector, *py_centroids, *py_distance_fn;
    PyObject *py_fast_centroids;
    PyObject *py_centroids_i;
    PyObject *py_arglist, *py_result;
    double dist;
    int closest_index;
    double closest_dist;

    if(!PyArg_ParseTuple(args, "OOO", 
			 &py_vector, &py_centroids, &py_distance_fn))
	return NULL;
    if(!PySequence_Check(py_centroids)) {
	PyErr_SetString(PyExc_ValueError, "centroids should be a list");
	return NULL;
    }
    if(!PyCallable_Check(py_distance_fn)) {
	PyErr_SetString(PyExc_ValueError, "distance_fn is not callable");
	return NULL;
    }
    if((length = PySequence_Length(py_centroids)) < 0)
	return NULL;
    if(!(py_fast_centroids = PySequence_Fast(py_centroids, "fast failed")))
	return NULL;

    closest_index = 0;   /* prevent compiler warnings */
    closest_dist = 0.0;
    for(i=0; i<length; i++) {
	py_centroids_i = PySequence_Fast_GET_ITEM(py_fast_centroids, i);
	if(!(py_arglist = Py_BuildValue("(OO)", py_vector, py_centroids_i)))
	    break;
	if(!(py_result = PyObject_CallObject(py_distance_fn, py_arglist)))
	    break;
	dist = PyNumber_AsDouble(py_result);
	if(PyErr_Occurred())
	    break;
	if((i == 0) || (dist < closest_dist)) {
	    closest_dist = dist;
	    closest_index = i;
	}
    }
    if(PyErr_Occurred())
	return NULL;
    return PyInt_FromLong((long)closest_index);
}


/* Module definition stuff */

static PyMethodDef ckMeansMethods[] = {
    {"_find_closest_centroid", ckMeans__find_closest_centroid, METH_VARARGS, 
     ckMeans__find_closest_centroid__doc__},
    {NULL, NULL}
};

static char ckMeans__doc__[] =
"XXX document here\n\
\n\
";

void initckMeans()
{
    (void) Py_InitModule3("ckMeans", ckMeansMethods, ckMeans__doc__);
}
