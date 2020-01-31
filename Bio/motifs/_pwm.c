/* Copyright 2009-2020 by Michiel de Hoon.  All rights reserved.
 *
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */

#include <Python.h>
#include <math.h>


static void
calculate(const char sequence[], int s, Py_ssize_t m, double* matrix,
          Py_ssize_t n, float* scores)
{
    Py_ssize_t i, j;
    char c;
    double score;
    int ok;
    float* p = scores;
#ifndef NAN
    float NAN = 0.0;
    NAN /= NAN;
#endif

    for (i = 0; i < n; i++)
    {
        score = 0.0;
        ok = 1;
        for (j = 0; j < m; j++)
        {
            c = sequence[i+j];
            switch (c)
            {
              /* Handling mixed case input here rather than converting it to
                 uppercase in Python code first, since doing so could use too
                 much memory if sequence is too long (e.g. chromosome or
                 plasmid). */
                case 'A':
                case 'a':
                    score += matrix[j*4+0]; break;
                case 'C':
                case 'c':
                    score += matrix[j*4+1]; break;
                case 'G':
                case 'g':
                    score += matrix[j*4+2]; break;
                case 'T':
                case 't':
                    score += matrix[j*4+3]; break;
                default:
                    ok = 0;
            }
        }
        if (ok) *p = (float)score;
        else *p = NAN;
        p++;
    }
}

static int
matrix_converter(PyObject* object, void* address)
{
    const int flags = PyBUF_C_CONTIGUOUS | PyBUF_FORMAT;
    char datatype;
    Py_buffer* view = address;

    if (object == NULL) goto exit;
    if (PyObject_GetBuffer(object, view, flags) == -1) {
        PyErr_SetString(PyExc_RuntimeError,
                        "position-weight matrix is not an array");
        return 0;
    }
    datatype = view->format[0];
    switch (datatype) {
        case '@':
        case '=':
        case '<':
        case '>':
        case '!': datatype = view->format[1]; break;
        default: break;
    }
    if (datatype != 'd') {
        PyErr_Format(PyExc_RuntimeError,
            "position-weight matrix data format incorrect "
            "('%c', expected 'd')", datatype);
        goto exit;
    }
    if (view->ndim != 2) {
        PyErr_Format(PyExc_RuntimeError,
            "position-weight matrix has incorrect rank (%d expected 2)",
            view->ndim);
        goto exit;
    }
    if (view->shape[1] != 4) {
        PyErr_Format(PyExc_RuntimeError,
            "position-weight matrix should have four columns "
            "(%zd columns found)", view->shape[1]);
        goto exit;
    }
    return Py_CLEANUP_SUPPORTED;

exit:
    PyBuffer_Release(view);
    return 0;
}

static int
scores_converter(PyObject* object, void* address)
{
    const int flags = PyBUF_C_CONTIGUOUS | PyBUF_FORMAT;
    char datatype;
    Py_buffer* view = address;

    if (object == NULL) goto exit;
    if (PyObject_GetBuffer(object, view, flags) == -1) return 0;
    datatype = view->format[0];
    switch (datatype) {
        case '@':
        case '=':
        case '<':
        case '>':
        case '!': datatype = view->format[1]; break;
        default: break;
    }
    if (datatype != 'f') {
        PyErr_Format(PyExc_RuntimeError,
            "scores array has incorrect data format ('%c', expected 'f')",
            datatype);
        goto exit;
    }
    if (view->ndim != 1) {
        PyErr_Format(PyExc_ValueError,
            "scores array has incorrect rank (%d expected 1)",
            view->ndim);
        goto exit;
    }
    return Py_CLEANUP_SUPPORTED;

exit:
    PyBuffer_Release(view);
    return 0;
}

static char calculate__doc__[] =
"    calculate(sequence, pwm, scores)\n"
"\n"
"This function calculates the position-weight matrix scores for all\n"
"positions along the sequence for position-weight matrix pwm, and stores\n"
"them in the provided numpy array scores.\n";

static PyObject*
py_calculate(PyObject* self, PyObject* args, PyObject* keywords)
{
    const char* sequence;
    static char* kwlist[] = {"sequence", "matrix", "scores", NULL};
    Py_ssize_t m;
    Py_ssize_t n;
    int s;
    PyObject* result = NULL;
    Py_buffer scores;
    Py_buffer matrix;

    matrix.obj = NULL;
    scores.obj = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, keywords, "s#O&O&", kwlist,
                                     &sequence,
                                     &s,
                                     matrix_converter, &matrix,
                                     scores_converter, &scores)) return NULL;
    m = matrix.shape[0];
    n = scores.shape[0];
    if (n == s - m + 1) {
        calculate(sequence, s, m, matrix.buf, n, scores.buf);
        Py_INCREF(Py_None);
        result = Py_None;
    }
    else {
        PyErr_SetString(PyExc_RuntimeError,
                        "size of scores array is inconsistent");
    }

    matrix_converter(NULL, &matrix);
    scores_converter(NULL, &scores);
    return result;
}

static struct PyMethodDef methods[] = {
   {"calculate",
    (PyCFunction)py_calculate,
    METH_VARARGS | METH_KEYWORDS,
    PyDoc_STR(calculate__doc__),
   },
   {NULL, NULL, 0, NULL} // sentinel
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_pwm",
    PyDoc_STR("Fast calculations involving position-weight matrices"),
    -1,
    methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject*
PyInit__pwm(void)
{
    return PyModule_Create(&moduledef);
}
