/* cMarkovModelmodule.c
 * jchang
 * Created: 1/13/01
 * Last modified: 11/26/02
 *
 * This optimizes some of the functions in MarkovModel.py.
 */

#include "Python.h"
#include "csupport.h"


/* Functions in this module. */

static char cMarkovModel__logadd__doc__[] = 
"_logadd(logx, logy) -> log(x+y)\n";

static PyObject *cMarkovModel__logadd(PyObject *self, PyObject *args)
{
    PyObject *py_logx, *py_logy;
    double logx, logy, minxy;
    double sum;

    if(!PyArg_ParseTuple(args, "OO", &py_logx, &py_logy))
	return NULL;
    logx = PyNumber_AsDouble(py_logx);
    logy = PyNumber_AsDouble(py_logy);
    if(PyErr_Occurred())
	return NULL;

    if(logy-logx > 100.0) {
	Py_INCREF(py_logy);
	return py_logy;
    } else if (logx-logy > 100.0) {
	Py_INCREF(py_logx);
	return py_logx;
    }
    minxy = (logx < logy) ? logx : logy;
    sum = minxy + log(exp(logx-minxy) + exp(logy-minxy));
    return PyFloat_FromDouble(sum);
}
	

/* Module definition stuff */

static PyMethodDef CMarkovModelMethods[] = {
  {"_logadd", cMarkovModel__logadd, METH_VARARGS, cMarkovModel__logadd__doc__},
  {NULL, NULL}
};

static char cMarkovModel__doc__[] =
"This module provides optimized replacement functions for MarkovModel.\n\
";

void initcMarkovModel(void)
{
  Py_InitModule3("cMarkovModel", CMarkovModelMethods, cMarkovModel__doc__);
}



