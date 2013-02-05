#include <Python.h>
#include "numpy/arrayobject.h"



static PyObject*
calculate(const char sequence[], int s, PyObject* matrix, npy_intp m)
{
    npy_intp n = s - m + 1;
    npy_intp i, j;
    char c;
    double score;
    int ok;
    PyObject* result;
    float* p;
    npy_intp shape = (npy_intp)n;
    float nan = 0.0;
    nan /= nan;
    if ((int)shape!=n)
    {
        PyErr_SetString(PyExc_ValueError, "integer overflow");
        return NULL;
    }
    result = PyArray_SimpleNew(1, &shape, NPY_FLOAT32);
    if (!result)
    {
        PyErr_SetString(PyExc_MemoryError, "failed to create output data");
        return NULL;
    }
    p = PyArray_DATA(result);
    for (i = 0; i < n; i++)
    {
        score = 0.0;
        ok = 1;
        for (j = 0; j < m; j++)
        {
            c = sequence[i+j];
            switch (c)
            {
                case 'A':
                case 'a':
                    score += *((double*)PyArray_GETPTR2(matrix, j, 0)); break;
                case 'C':
                case 'c':
                    score += *((double*)PyArray_GETPTR2(matrix, j, 1)); break;
                case 'G':
                case 'g':
                    score += *((double*)PyArray_GETPTR2(matrix, j, 2)); break;
                case 'T':
                case 't':
                    score += *((double*)PyArray_GETPTR2(matrix, j, 3)); break;
                default:
                    ok = 0;
            }
        }
        if (ok) *p = (float)score;
        else *p = nan;
        p++;
    }
    return result;
}

static char calculate__doc__[] = 
"    calculate(sequence, pwm) -> array of score values\n"
"\n"
"This function calculates the position-weight matrix scores for all\n"
"positions along the sequence, and returns them as a Numerical Python\n"
"array.\n";

static PyObject*
py_calculate(PyObject* self, PyObject* args, PyObject* keywords)
{
    const char* sequence;
    PyObject* matrix = NULL;
    static char* kwlist[] = {"sequence", "matrix", NULL};
    npy_intp m;
    int s;
    PyObject* result;
    if(!PyArg_ParseTupleAndKeywords(args, keywords, "s#O&", kwlist,
                                    &sequence,
                                    &s,
                                    PyArray_Converter,
                                    &matrix)) return NULL;

    if (PyArray_TYPE(matrix) != NPY_DOUBLE)
    {
        PyErr_SetString(PyExc_ValueError,
            "position-weight matrix should contain floating-point values");
        result = NULL;
    }
    else if (PyArray_NDIM(matrix) != 2) /* Checking number of dimensions */
    {
        result = PyErr_Format(PyExc_ValueError, 
            "position-weight matrix has incorrect rank (%d expected 2)",
            PyArray_NDIM(matrix));
    }
    else if(PyArray_DIM(matrix, 1) != 4)
    {
        result = PyErr_Format(PyExc_ValueError,
            "position-weight matrix should have four columns (%" NPY_INTP_FMT
            " columns found)", PyArray_DIM(matrix, 1));
    }
    else
    {
        m = PyArray_DIM(matrix, 0);
        result = calculate(sequence, s, matrix, m);
    }
    Py_DECREF(matrix);
    return result;
}


static struct PyMethodDef methods[] = {
   {"calculate", (PyCFunction)py_calculate, METH_VARARGS | METH_KEYWORDS, calculate__doc__},
   {NULL,          NULL, 0, NULL} /* sentinel */
};


#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_pwm",
        "Fast calculations involving position-weight matrices",
        -1,
        methods,
        NULL,
        NULL,
        NULL,
        NULL
};

PyObject*
PyInit__pwm(void)

#else

void init_pwm(void)
#endif
{
  PyObject *m;

  import_array();

#if PY_MAJOR_VERSION >= 3
  m = PyModule_Create(&moduledef);
  if (m==NULL) return NULL;
#else
  m = Py_InitModule4("_pwm",
                     methods,
                     "Fast calculations involving position-weight matrices",
                     NULL,
                     PYTHON_API_VERSION);
  if (m==NULL) return;
#endif

  if (PyErr_Occurred()) Py_FatalError("can't initialize module _pwm");
#if PY_MAJOR_VERSION >= 3
    return m;
#endif
}
