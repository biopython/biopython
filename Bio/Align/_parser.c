/* Copyright 2024 by Michiel de Hoon.  All rights reserved.
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */


#include <Python.h>
#include <string.h>
#include <stdbool.h>

static PyTypeObject CoordinatesType;

typedef struct {
    PyObject_HEAD
    long** data;
    Py_ssize_t n;
    Py_ssize_t m;
    Py_ssize_t k;
    char eol;
} Coordinates;

/* Destructor function */
static void
Coordinates_dealloc(Coordinates *self)
{
    long** data = self->data;
    ssize_t k;
    const ssize_t n = self->n;
    if (data) {
        for (k = 0; k < n; k++) {
            if (data[k] == NULL) break;
            PyMem_Free(data[k]);
        }
        PyMem_Free(data);
    }
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static int
array_converter(PyObject* argument, void* pointer)
{
    Py_buffer* view = pointer;
    Coordinates* self;

    if (!PyObject_TypeCheck(view->obj, &CoordinatesType)) {
        PyErr_SetString(PyExc_RuntimeError,
                        "expected an object of the Coordinates class");
        return 0;
    }

    self = (Coordinates*) view->obj;

    if (PyObject_GetBuffer(argument, view, PyBUF_CONTIG) != 0) {
        PyErr_SetString(PyExc_RuntimeError,
                        "argument does not implement the buffer protocol");
        return 0;
    }
    if (view->ndim != 2) {
        PyErr_Format(PyExc_RuntimeError,
                     "buffer has incorrect rank %d (expected 2)",
                      view->ndim);
    }
    else if (view->shape[0] != self->n) {
        PyErr_Format(PyExc_RuntimeError,
                     "buffer has incorrect number of rows %zd (expected %zd)",
                      view->shape[0], self->n);
    }
    else if (view->shape[1] != self->k) {
        PyErr_Format(PyExc_RuntimeError,
                     "buffer has incorrect number of columns %zd (expected %zd)",
                      view->shape[1], self->k);
    }
    else if (view->itemsize != sizeof(long)) {
        PyErr_Format(PyExc_RuntimeError,
                    "buffer has unexpected item byte size "
                    "(%ld, expected %ld)", view->itemsize, sizeof(long));
    }
    else return Py_CLEANUP_SUPPORTED;

    PyBuffer_Release(view);
    return 0;
}

static PyObject*
Coordinates_feed(Coordinates* self, PyObject* args)
{
    PyObject* line;
    const char* buffer;
    const char* s;
    const char eol = self->eol;
    Py_ssize_t n = self->n;
    Py_ssize_t m = self->m;
    Py_ssize_t size = 2;
    Py_ssize_t i = 0;
    Py_ssize_t p = 0;
    Py_ssize_t offset = 0;
    long** data;
    long* row;
    char c;

    if (!PyArg_ParseTuple(args, "S|n:feed", &line, &offset)) return NULL;

    buffer = PyBytes_AS_STRING(line) + offset;

    s = buffer;
    row = PyMem_Malloc(size*sizeof(long));
    if (!row) return NULL;
    if (*s == '-') row[i++] = 0;

    data = PyMem_Realloc(self->data, (n+1)*size*sizeof(long*));
    if (!data) {
        PyMem_Free(row);
        return NULL;
    }
    self->data = data;
    data[n] = row;

    while (*s != '\0' && *s != eol) {
        if (*s == '-') {
            do s++; while (*s == '-');
        }
        else {
            p -= (s - buffer);
            do c = *(++s); while (c != '-' && c != eol && c != '\0');
            p += (s - buffer);
        }

        if (i == size) {
            size *= 2;
            row = PyMem_Realloc(row, size*sizeof(long));
            if (!row) {
                PyMem_Free(data[n]);
                return NULL;
            }
            data[n] = row;
        }
        row[i++] = s - buffer;
    }
    row = PyMem_Realloc(row, i*sizeof(long));
    if (!row) {
        PyMem_Free(data[n]);
        return NULL;
    }
    data[n] = row;
    m = s - buffer;
    if (n == 0) self->m = m;
    else if (buffer + m != s) {
        PyErr_Format(PyExc_ValueError,
                     "line has length %zi (expected %zi)", m, self->m);
        return NULL;
    }

    n++;
    self->n = n;

    PyObject* sequence = PyBytes_FromStringAndSize(NULL, p);
    if (!sequence) return NULL;
    char* d;
    d = PyBytes_AS_STRING(sequence);
    int ii;
    Py_ssize_t end = 0;
    Py_ssize_t start;
    Py_ssize_t step;
    s = buffer;
    bool gap = false;
    ii = 0;
    if (row[ii] == 0) {
        gap = true;
        ii++;
    }
    for ( ; ii < i; ii++) {
        start = end;
        end = row[ii];
        step = end - start;
        gap = !gap;
        if (gap) {
            s = buffer + start;
            memcpy(d, s, step);
            d += step;
        }
    }
    *d = '\0';

    PyObject* result = Py_BuildValue("lO", m, sequence);
    Py_DECREF(sequence);
    return result;
}

static PyObject*
Coordinates_fill(Coordinates* self, PyObject* args)
{
    Py_buffer view;
    Py_ssize_t i, j, k, n, m;
    Py_ssize_t index;
    long step;
    Py_ssize_t* indices = NULL;
    Py_ssize_t** pointer = NULL;
    bool* gaps = NULL;
    long* data;

    n = self->n;
    if (n == 0) Py_RETURN_NONE;

    view.obj = (PyObject*) self;

    if (!PyArg_ParseTuple(args, "O&:fill", array_converter, &view))
        return NULL;

    data = view.buf;
    k = view.shape[1];
    if (n != view.shape[0]) {
        PyErr_Format(PyExc_ValueError,
                     "expected an array with %zd rows (found %zd rows)",
                     n, view.shape[0]);
        return 0;
    }
    for (i = 0; i < n; i++) data[i*k] = 0;

    m = self->m;

    indices = PyMem_Calloc(n, sizeof(Py_ssize_t));
    if (!indices) goto exit;

    gaps = PyMem_Malloc(n * sizeof(bool));
    if (!gaps) goto exit;

    pointer = PyMem_Calloc(n, sizeof(Py_ssize_t*));
    if (!pointer) goto exit;

    for (i = 0; i < n; i++) {
        pointer[i] = self->data[i];
        if (pointer[i][0] == 0) {
            gaps[i] = true;
            pointer[i]++;
        }
        else {
            gaps[i] = false;
        }
    }

    j = 0;
    index = 0;
    do {
        j++;
        for (i = 0; i < n; i++) {
            if (index == indices[i]) indices[i] = *pointer[i];
        }
        step = m;
	for (i = 0; i < n; i++) if (indices[i] < step) step = indices[i];
        step -= index;
        for (i = 0; i < n; i++) {
            if (gaps[i] == true) data[i*k+j] = data[i*k+j-1];
            else data[i*k+j] = data[i*k+j-1] + step;
            if (step + index == indices[i]) {
                pointer[i]++;
                gaps[i] = !gaps[i];
            }
        }
        index += step;
    }
    while (index < m);

exit:
    if (indices) PyMem_Free(indices);
    if (pointer) PyMem_Free(pointer);
    if (gaps) PyMem_Free(gaps);
    Py_RETURN_NONE;
}

static PyObject*
Coordinates_get_shape(Coordinates* self, void* closure)
{
    Py_ssize_t i;
    Py_ssize_t index;
    Py_ssize_t min_index;
    long** data = self->data;
    const Py_ssize_t n = self->n;
    const Py_ssize_t m = self->m;
    Py_ssize_t k = 1;

    if (n > 0) {
        data = PyMem_Calloc(n, sizeof(long*));
        if (!data) return NULL;
        memcpy(data, self->data, n*sizeof(long*));
        for (i = 0; i < n; i++) {
            index = *(data[i]);
            if (index == 0) {
                k = 0;
                break;
            }
        }
        while (1) {
            k++;
            min_index = m;
            for (i = 0; i < n; i++) {
                index = *(data[i]);
                if (index < min_index) min_index = index;
            }
            if (min_index == m) break;
            for (i = 0; i < n; i++) {
                index = *(data[i]);
	        if (index == min_index) data[i]++;
            }
        }
    }

    self->k = k;
    return Py_BuildValue("ll", n, k);
}

static PyGetSetDef Coordinates_getset[] = {
    {"shape", (getter)Coordinates_get_shape,
     NULL, "return the shape of the coordinates matrix", NULL,
    },
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

static PyMethodDef Coordinates_methods[] = {
    {"feed", (PyCFunction)Coordinates_feed, METH_VARARGS, "Feed a bytes object to the parser. The parser will read from the buffer until it finds the end-of-line character, and return the number of bytes read."},
    {"fill",
     (PyCFunction)Coordinates_fill,
     METH_VARARGS,
    "Fill in the coordinates array based on the alignment lines fed to the parser so far.",
    },
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

static PyObject*
Coordinates_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    char eol = '\n';  /* end-of-line character */
    static char *kwlist[] = {"eol", NULL};

    Coordinates *self;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|c", kwlist, &eol))
        return NULL;

    self= (Coordinates *)type->tp_alloc(type, 0);
    if (!self) return NULL;
    self->eol = eol;
    self->n = 0;
    return (PyObject *)self;
}

static PyTypeObject CoordinatesType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_parser.Coordinates",
    .tp_doc = "Coordinates objects",
    .tp_basicsize = sizeof(Coordinates),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_dealloc = (destructor)Coordinates_dealloc,
    .tp_methods = Coordinates_methods,
    .tp_getset = Coordinates_getset,
    .tp_new = (newfunc)Coordinates_new,
};

/* Module initialization function */
static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_parser",
    .m_doc = "fast C implementation of utility functions for alignment parsing",
    .m_size = -1,
};

PyMODINIT_FUNC
PyInit__parser(void)
{
    PyObject *module;

    if (PyType_Ready(&CoordinatesType) < 0)
        return NULL;

    module = PyModule_Create(&moduledef);
    if (module == NULL)
        return NULL;

    Py_INCREF(&CoordinatesType);
    PyModule_AddObject(module, "Coordinates", (PyObject *)&CoordinatesType);
    return module;
}
