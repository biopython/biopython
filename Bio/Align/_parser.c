/* Copyright 2024 by Michiel de Hoon.  All rights reserved.
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */



#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <float.h>


typedef struct {
    Py_buffer* buffers;
    Py_ssize_t n;
    Py_ssize_t m;
} Alignment;


static int
alignment_converter(PyObject* argument, void* pointer)
{
    Py_ssize_t i, n, m;
    Py_buffer* view;
    Alignment* alignment = pointer;
    PyObject* sequence;
    if (argument == NULL) {
        n = alignment->n;
        for (i = 0; i < n; i++) PyBuffer_Release(&alignment->buffers[i]);
        PyMem_Free(alignment->buffers);
        return 1;
    }
    if (!PyList_Check(argument)) {
       PyErr_SetString(PyExc_ValueError, "expected a list");
       alignment->n = 0;
       return 0;
    }
    n = PyList_GET_SIZE(argument);
    alignment->buffers = PyMem_Calloc(n+1, sizeof(Py_buffer));
    if (!alignment->buffers) {
        PyErr_SetNone(PyExc_MemoryError);
        return 0;
    }
    alignment->n = n;
    for (i = 0; i < n; i++) {
        sequence = PyList_GET_ITEM(argument, i);
        view = &alignment->buffers[i];
        if (PyObject_GetBuffer(sequence, view,
                               PyBUF_FORMAT | PyBUF_C_CONTIGUOUS) == 0) {
            if (view->ndim != 1) {
                PyErr_Format(PyExc_ValueError,
                             "sequence %zi has incorrect rank %d (expected 1)",
                             i, view->ndim);
                goto error;
            }
            if (view->itemsize != 1) {
                PyErr_Format(PyExc_ValueError,
                    "sequence %zi has incorrect item size %d (expected 1)",
                    i, view->itemsize);
                goto error;
            }
            if (i == 0) {
                m = view->len;
                if (m == 0) {
                    PyErr_SetString(PyExc_ValueError,
                                    "sequence 0 has zero length");
                    goto error;
                }
            }
            else {
                if (view->len != m) {
                    PyErr_Format(PyExc_ValueError,
                        "sequence %zi has inconsistent length %zi "
                        "(expected %zi)", i, view->len, m);
                    goto error;
                }
            }
        }
    }
    alignment->m = m;
    return Py_CLEANUP_SUPPORTED;

error:
    for (i = 0; i < n; i++) {
        view = &alignment->buffers[i];
        PyBuffer_Release(view);
    }
    PyMem_Free(alignment->buffers);
    return 0;
}

static int
sizes_converter(PyObject* argument, void* pointer)
{
    Py_buffer* view = pointer;
    if (argument == NULL) {
        PyBuffer_Release(view);
        return 1;
    }
    if (PyObject_GetBuffer(argument,
                           view,
                           PyBUF_FORMAT | PyBUF_C_CONTIGUOUS) != 0) return 0;
    if (view->itemsize != sizeof(long)) {
        PyErr_Format(PyExc_ValueError,
            "sizes has incorrect item size %zi (expected %zi)",
            view->itemsize, sizeof(long));
        PyBuffer_Release(view);
        return 0;
    }
    if (view->ndim != 1) {
        PyErr_Format(PyExc_ValueError,
                     "sizes array has incorrect rank %d (expected 2)",
                     view->ndim);
        PyBuffer_Release(view);
        return 0;
    }
    return Py_CLEANUP_SUPPORTED;
}

static int
strands_converter(PyObject* argument, void* pointer)
{
    Py_buffer* view = pointer;
    if (!argument) {
        PyBuffer_Release(view);
        return 1;
    }
    if (PyObject_GetBuffer(argument, view,
                           PyBUF_FORMAT | PyBUF_C_CONTIGUOUS) != 0) return 0;
    if (view->ndim != 1) {
        PyErr_Format(PyExc_ValueError,
                     "strands has incorrect rank %d (expected 1)",
                     view->ndim);
        PyBuffer_Release(view);
        return 0;
    }
    if (view->itemsize != sizeof(long)) {
        PyErr_Format(PyExc_ValueError,
            "strands has incorrect item size %zi (expected %zi)",
            view->itemsize, sizeof(long));
        PyBuffer_Release(view);
        return 0;
    }
    return Py_CLEANUP_SUPPORTED;
}

static int
coordinates_converter(PyObject* argument, void* pointer)
{
    Py_buffer* view = pointer;
    if (argument == NULL) {
        PyBuffer_Release(view);
        return 1;
    }
    if (PyObject_GetBuffer(argument,
                           view,
                           PyBUF_FORMAT | PyBUF_C_CONTIGUOUS) != 0) return 0;
    if (view->itemsize != sizeof(long)) {
        PyErr_Format(PyExc_ValueError,
            "coordinates has incorrect item size %zi (expected %zi)",
            view->itemsize, sizeof(long));
        PyBuffer_Release(view);
        return 0;
    }
    if (view->ndim != 2) {
        PyErr_Format(PyExc_ValueError,
                     "coordinates array has incorrect rank %d (expected 2)",
                     view->ndim);
        PyBuffer_Release(view);
        return 0;
    }
    return Py_CLEANUP_SUPPORTED;
}

/* Module definition */

static char _parser__doc__[] =
"C extension module implementing fast functions for alignment file parsing";

static PyObject*
calculate_alignment_coordinates_columns(PyObject* self, PyObject* args)
{
    Py_ssize_t i, n, m;
    Alignment alignment;
    Py_ssize_t k;
    char* s;
    Py_ssize_t p;
    Py_ssize_t position;
    Py_ssize_t* positions;
    PyObject* result = NULL;

    if (!PyArg_ParseTuple(args, "O&",
                                alignment_converter,  &alignment))
        return NULL;

    n = alignment.n;
    m = alignment.m;

    positions = PyMem_Calloc(n, sizeof(Py_ssize_t));
    if (positions == NULL) goto exit;
    k = 0;
    position = 0;
    while (1) {
        p = position;
        position = m;
        for (i = 0; i < n; i++) {
            if (positions[i] < position) position = positions[i];
        }
        for (i = 0; i < n; i++) {
            s = alignment.buffers[i].buf + p;
        }
        k++;
        if (position == m) break;
        for (i = 0; i < n; i++) {
            if (position < positions[i]) continue;
            s = alignment.buffers[i].buf + position;
            if (*s == '-') {
                for (p = position + 1; p < m; p++) {
                    if (*(++s) != '-') break;
                }
            }
            else {
                for (p = position + 1; p < m; p++) {
                    if (*(++s) == '-') break;
                }
            }
            positions[i] = p;
        }
    }
    PyMem_Free(positions);
    result = PyLong_FromSsize_t(k);
exit:
    for (i = 0; i < n; i++) {
        PyBuffer_Release(&alignment.buffers[i]);
    }
    PyMem_Free(alignment.buffers);
    return result;
}

static PyObject*
fill_alignment_coordinates(PyObject* self, PyObject* args)
{
    Py_ssize_t i, n, m;
    Alignment alignment;
    Py_buffer strands;
    Py_buffer coordinates;
    Py_ssize_t k;
    char* s;
    Py_ssize_t p;
    Py_ssize_t previous;
    Py_ssize_t position;
    Py_ssize_t* positions = NULL;
    Py_ssize_t* indices = NULL;
    long* c;
    long* sign;
    Py_ssize_t q;
    Py_ssize_t step;
    PyObject* sequences = NULL;
    Py_buffer sizes;
    PyObject* b;

    if (!PyArg_ParseTuple(args, "O&O&O&O&",
                                alignment_converter,  &alignment,
                                coordinates_converter,  &coordinates,
                                strands_converter,  &strands,
                                sizes_converter,  &sizes))
        return NULL;

    n = alignment.n;
    m = alignment.m;
    if (coordinates.shape[0] != n) { 
        PyErr_Format(PyExc_ValueError,
            "coordinates has incorrect number of rows %zi (expected %zi)",
            coordinates.shape[0], n);
        goto exit;
    }
    q = coordinates.shape[1];
    c = coordinates.buf;
    sign = strands.buf;

    sequences = PyTuple_New(n);
    if (sequences == NULL) goto exit;
    for (i = 0; i < n; i++) {
        long size = ((long*)sizes.buf)[i];
        b = PyBytes_FromStringAndSize(NULL, size);
        if (!b) {
            Py_DECREF(sequences);
            sequences = NULL;
            goto exit;
        }
        PyTuple_SET_ITEM(sequences, i, b);
    }

    positions = PyMem_Calloc(n, sizeof(Py_ssize_t));
    if (positions == NULL) goto exit;
    k = 1;
    position = 0;
    for (i = 0; i < n; i++) {
        s = alignment.buffers[i].buf;
        if (*s == '-') {
            for (p = position + 1; p < m; p++) if (*(++s) != '-') break;
        }
        else {
            for (p = position + 1; p < m; p++) if (*(++s) == '-') break;
        }
        positions[i] = p;
    }
    indices = PyMem_Calloc(n, sizeof(Py_ssize_t));
    if (indices == NULL) goto exit;
    while (1) {
        previous = position;
        position = m;
        for (i = 0; i < n; i++) {
            if (positions[i] < position) position = positions[i];
        }
        step = position - previous;
        if (position == m) break;
        if (k >= q - 1) {
            PyErr_SetString(PyExc_ValueError,
                            "coordinates array size is insufficient.");
            Py_DECREF(sequences);
            sequences = NULL;
            goto exit;
        }
        for (i = 0; i < n; i++) {
            b = PyTuple_GET_ITEM(sequences, i);
            s = alignment.buffers[i].buf + position;
            if (positions[i] == position) {
                if (*s == '-') {
                    memcpy(PyBytes_AS_STRING(b) + indices[i], alignment.buffers[i].buf + previous, step);
                    for (p = position + 1; p < m; p++) {
                        if (*(++s) != '-') break;
                    }
                    c[i * q + k] = c[i * q + k - 1] + sign[i] * step;
                    indices[i] += step;
                }
                else {
                    for (p = position + 1; p < m; p++) {
                        if (*(++s) == '-') break;
                    }
                    c[i * q + k] = c[i * q + k - 1];
                }
                positions[i] = p;
            }
            else {
                if (*s == '-') {
                    c[i * q + k] = c[i * q + k - 1];
                }
                else {
                    c[i * q + k] = c[i * q + k - 1] + sign[i] * step;
                    memcpy(PyBytes_AS_STRING(b) + indices[i], alignment.buffers[i].buf + previous, step);
                    indices[i] += step;
                }
            }
        }
        k++;
    }
    for (i = 0; i < n; i++) {
        b = PyTuple_GET_ITEM(sequences, i);
        s = alignment.buffers[i].buf + previous;
        if (*s == '-') {
            c[i * q + k] = c[i * q + k - 1];
        }
        else {
            c[i * q + k] = c[i * q + k - 1] + sign[i] * step;
            memcpy(PyBytes_AS_STRING(b) + indices[i], s, step);
            indices[i] += step;
        }
    }
exit:
    if (positions) PyMem_Free(positions);
    if (indices) PyMem_Free(indices);
    for (i = 0; i < n; i++) {
        PyBuffer_Release(&alignment.buffers[i]);
    }
    PyMem_Free(alignment.buffers);
    strands_converter(NULL, &strands);
    coordinates_converter(NULL, &coordinates);
    if (sequences == NULL) return NULL;
    if (k + 1 != q) {
        PyErr_SetString(PyExc_ValueError,
                        "unexpected size of coordinates array");
        Py_DECREF(sequences);
        return NULL;
    }
    return sequences;
}


static struct PyMethodDef methods[] = {
    {"calculate_alignment_coordinates_columns",
     (PyCFunction) calculate_alignment_coordinates_columns,
     METH_VARARGS,
     "Calculate the number of columns needed for the coordinates array to represent the sequence alignment",
    },
    {"fill_alignment_coordinates",
     (PyCFunction) fill_alignment_coordinates,
     METH_VARARGS,
     "Fill in the coordinates array based on the gapped sequence alignment",
    },
    {NULL, NULL, 0, NULL} /* sentinel */
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_parser",
        _parser__doc__,
        -1,
        methods,
        NULL,
        NULL,
        NULL,
        NULL
};

PyObject *
PyInit__parser(void)
{
    return PyModule_Create(&moduledef);
}
