/* Copyright 2000 by Jeffrey Chang.  All rights reserved.
 * This code is part of the Biopython distribution and governed by its
 * license.  Please see the LICENSE file that should have been included
 * as part of this package.
 *
 * clistfnsmodule.c
 * Created 3 Jun 2000
 */

#include "Python.h"
#include <math.h>




/************************************** Exported Functions ***********/

static char clistfns_count__doc__[] = 
"count(items) -> dict of counts of each item\n\
\n\
Count the number of times each item appears in a list of data.\n\
\n\
";

static PyObject *
clistfns_count(self, args)
     PyObject *self;
     PyObject *args;
{
    int i;
    PyObject *items, *counts;
    PyObject *item, *count;
    long int current;

    if(!PyArg_ParseTuple(args, "O", &items))
	return NULL;
    if(!PySequence_Check(items)) {
	PyErr_SetString(PyExc_TypeError, "expected mapping type");
	return NULL;
    }
    
    if(!(counts = PyDict_New()))
	return NULL;
    
    /* Go through the loop, counting how often each item appears. */
    i = 0;
    while(1) {
	if(!(item = PySequence_GetItem(items, i))) 
	    break; /* no more numbers */

	if(!(count = PyDict_GetItem(counts, item))) {
	    /* This is a new item.  Set count to 1. */
	    count = PyInt_FromLong(1);
	    PyDict_SetItem(counts, item, count);
	    Py_DECREF(count);
	    if(PyErr_Occurred()) 
		/* item not hashable? */
		return NULL;
	} else {
	    if(!PyInt_Check(count)) {
		PyErr_SetString(PyExc_SystemError, "non-int in dictionary");
		return NULL;
	    }
	    current = PyInt_AsLong(count);
	    count = PyInt_FromLong(current+1);
	    if(PyErr_Occurred()) 
		return NULL;
	    PyDict_SetItem(counts, item, count);
	    Py_DECREF(count);
	    if(PyErr_Occurred()) 
		return NULL;
	}

	i++;
    }
    PyErr_Clear(); /* clear the exception set by PySequence_GetItem */
    
    return counts;
}


static char clistfns_contents__doc__[] = 
"contents(items) -> dict of item -> percentage\n\
\n\
Summarize the contents of the list in terms of the percentages of each\n\
item.  For example, if an item appears 3 times in a list with 10 items,\n\
it is in 0.3 of the list\n\
\n\
";

static PyObject *
clistfns_contents(self, args)
     PyObject *self;
     PyObject *args;
{
    int i;
    PyObject *items, *counts, *count;
    PyObject *percentages, *perc;
    PyObject *keys, *key;
    long c;
    double total;
    int real_error;

    if(!PyArg_ParseTuple(args, "O", &items))
	return NULL;
    if(!PySequence_Check(items)) {
	PyErr_SetString(PyExc_TypeError, "expected mapping type");
	return NULL;
    }
    if((total = PySequence_Length(items)) == -1) {
	PyErr_SetString(PyExc_ValueError, "I couldn't get length of item.");
	return NULL;
    }
    
    counts = clistfns_count(self, args);
    if(PyErr_Occurred())
	return NULL;

    if(!(percentages = PyDict_New())) {
	Py_DECREF(counts);
	return NULL;
    }
    
    /* Loop through every element in counts, calculating the probabilities. */
    if(!(keys = PyMapping_Keys(counts))) {
	Py_DECREF(counts);
	return NULL;
    }

    /* Go through the loop, counting how often each item appears. */
    real_error = 0;
    i = 0;
    while(1) {
	if(!(key = PySequence_GetItem(keys, i))) 
	    break; /* no more numbers */
	if(!(count = PyDict_GetItem(counts, key))) {
	    real_error = 1;
	    break;
	}
	if(!PyInt_Check(count)) {
	    real_error = 1;
	    PyErr_SetString(PyExc_SystemError, "non-int in dictionary");
	}
	c = PyInt_AsLong(count);
	perc = PyFloat_FromDouble((double)c / total);
	PyDict_SetItem(percentages, key, perc);
	Py_DECREF(perc);
	if(PyErr_Occurred()) {
	    real_error = 1;
	    break;
	}
	i++;
    }
    if(real_error) {
	Py_DECREF(percentages);
	percentages = NULL;
    }
    else PyErr_Clear(); /* clear the exception set by PySequence_GetItem */
    
    Py_DECREF(counts);
    Py_DECREF(keys);
    
    return percentages;
}


/************************************** Module definition stuff ******/

static PyMethodDef clistfnsMethods[] = {
    {"count", clistfns_count, METH_VARARGS, clistfns_count__doc__},
    {"contents", clistfns_contents, METH_VARARGS, clistfns_contents__doc__},
    {NULL, NULL}
};

static char clistfns__doc__[] =
"This provides helper functions for the listfns module.\n\
You should never import this module on its own.\n\
\n\
";

void initclistfns()
{
    (void) Py_InitModule3("clistfns", clistfnsMethods, clistfns__doc__);
}
