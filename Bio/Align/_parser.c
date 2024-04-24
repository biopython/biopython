/* Copyright 2024 by Michiel de Hoon.  All rights reserved.
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */


#include <Python.h>
#include <string.h>


static PyTypeObject CoordinatesType;

typedef struct {
    PyObject_HEAD
    Py_ssize_t shape[2];
    long** positions;
    size_t kmax;
} Coordinates;

/* Destructor function */
static void
Coordinates_dealloc(Coordinates *self)
{
    long** positions = self->positions;
    ssize_t k;
    const ssize_t kmax = self->kmax;
    for (k = 0; k < kmax; k++) {
        if (positions[k] == NULL) break;
        PyMem_Free(positions[k]);
    }
    PyMem_Free(self->positions);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static int
converter(PyObject* argument, void* pointer)
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
    else if (view->shape[0] != self->shape[0]) {
        PyErr_Format(PyExc_RuntimeError,
                     "buffer has incorrect number of rows %zd (expected %zd)",
                      view->shape[0], self->shape[0]);
    }
    else if (view->shape[1] != self->shape[1]) {
        PyErr_Format(PyExc_RuntimeError,
                     "buffer has incorrect number of columns %zd (expected %zd)",
                      view->shape[1], self->shape[1]);
    }
    else if (view->itemsize != sizeof(long)) {
        PyErr_Format(PyExc_RuntimeError,
                    "buffer has unexpected item byte size "
                    "(%ld, expected %ld)", view->itemsize, sizeof(long));
    }
    else return Py_CLEANUP_SUPPORTED;

    PyBuffer_Release(view);
    return 1;
}

static PyObject*
Coordinates_fill(Coordinates* self, PyObject* args)
{
    Py_buffer view;
    ssize_t n = self->shape[0];
    ssize_t m = self->shape[1];
    ssize_t i, j;
    long** positions = self->positions;
    long* data;

    view.obj = (PyObject*) self;

    if (!PyArg_ParseTuple(args, "O&:fill", converter, &view))
        return NULL;

    data  = view.buf;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            data[i*m+j] = positions[j][i];
        }
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject*
Coordinates_get_shape(Coordinates* self, void* closure)
{
    return Py_BuildValue("ll", self->shape[0], self->shape[1]);
}

static PyGetSetDef Coordinates_getset[] = {
    {"shape", (getter)Coordinates_get_shape,
     NULL, "return the shape of the coordinates matrix", NULL,
    },
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

static PyMethodDef Coordinates_methods[] = {
    {"fill", (PyCFunction)Coordinates_fill, METH_VARARGS, "Fill a buffer object with the calculated indices"},
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

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
};

static PyObject*
parse_printed_alignment(PyObject* module, PyObject* args)
{
    PyObject* lines;
    PyObject* line;
    Py_ssize_t i, k, n, m;
    const char* s;
    Py_ssize_t index;
    long step;
    Py_ssize_t* indices = NULL;
    Coordinates *coordinates = NULL;
    PyObject* sequences = NULL;
    PyObject* sequence;
    char** destinations = NULL;
    long** positions;
    long* steps = NULL;
    Py_ssize_t kmax = 2;

    PyObject* result = NULL;

    if (!PyArg_ParseTuple(args, "O!:parse_printed_alignment",
                                &PyList_Type, &lines))
        return NULL;

    n = PyList_GET_SIZE(lines);
    for (i = 0; i < n; i++) {
        line = PyList_GET_ITEM(lines, i);
        if (!PyBytes_Check(line)) {
            PyErr_Format(PyExc_ValueError,
                         "line [%zi] is not a bytes object", i);
            return 0;
        }
        if (i == 0) m = PyBytes_GET_SIZE(line);
        else if (PyBytes_GET_SIZE(line) != m) {
            PyErr_Format(PyExc_ValueError,
                "all lines must have the same length "
                "(line [0] has length %zi, line [%zi] has length %zi).",
                m, i, PyBytes_GET_SIZE(line));
            return 0;
        }
    }

    coordinates = (Coordinates *)PyType_GenericAlloc(&CoordinatesType, 0);
    if (coordinates == NULL) return NULL;
    if (n == 0) {
        coordinates->shape[0] = 0;
        coordinates->shape[1] = 0;
        result = Py_BuildValue("[]O", coordinates);
        Py_DECREF(coordinates);
        return result;
    }

    indices = PyMem_Calloc(n, sizeof(Py_ssize_t));
    if (!indices) goto exit;

    positions = PyMem_Calloc(kmax, sizeof(long*));
    if (!positions) goto exit;
    positions[0] = PyMem_Calloc(n, sizeof(long));
    if (!positions[0]) {
        PyMem_Free(positions);
        goto exit;
    }
    coordinates->positions = positions;
    coordinates->kmax = kmax;

    steps = PyMem_Calloc(kmax, sizeof(long));
    if (!steps) goto exit;

    k = 0;
    index = 0;
    do {
        k++;
        for (i = 0; i < n; i++) {
            if (index == indices[i]) {
                line = PyList_GET_ITEM(lines, i);
                s = PyBytes_AS_STRING(line) + index;
                if (*s == '-') {
                    do s++; while (*s == '-');
                }
                else {
                    do {
                        s++;
                        if (*s == '\0') {
                            if (s - PyBytes_AS_STRING(line) != m) break;
                            /* otherwise, it is an internal null character */
                        }
                        else if (*s == '-') break;
                    } while (1);
                }
                indices[i] = s - PyBytes_AS_STRING(line);
            }
        }
        if (k == kmax) {
            long* new_steps;
            kmax *= 2;
            new_steps = PyMem_Realloc(steps, kmax*sizeof(long));
            if (new_steps == NULL) goto exit;
            steps = new_steps;
            positions = PyMem_Realloc(positions, kmax*sizeof(long*));
            if (positions == NULL) {
                positions = coordinates->positions;
                goto exit;
            }
            coordinates->positions = positions;
            coordinates->kmax = kmax;
            memset(positions + k, '\0', k * sizeof(long*));
        }
        step = m;
        for (i = 0; i < n; i++) if (indices[i] < step) step = indices[i];
        step -= index;
        positions[k] = PyMem_Malloc(n*sizeof(long));
        if (!positions[k]) goto exit;
        for (i = 0; i < n; i++) {
            s = PyBytes_AS_STRING(PyList_GET_ITEM(lines, i));
            if (s[index] == '-') {
                positions[k][i] = positions[k-1][i];
            }
            else {
                positions[k][i] = positions[k-1][i] + step;
            }
        }
        index += step;
        steps[k] = step;
    }
    while (index < m);

    sequences = PyList_New(n);
    if (sequences == NULL) goto exit;
    destinations = PyMem_Calloc(n, sizeof(char*));
    if (destinations == NULL) goto exit;

    for (i = 0; i < n; i++) {
        sequence = PyBytes_FromStringAndSize(NULL, positions[k][i]);
        if (!sequence) goto exit;
        PyList_SET_ITEM(sequences, i, sequence);
        destinations[i] = PyBytes_AS_STRING(sequence);
    }

    index = 0;
    for (m = 1; m <= k; m++) {
        step = steps[m];
        for (i = 0; i < n; i++) {
            if (positions[m][i] == positions[m-1][i]) continue;
            s = PyBytes_AS_STRING(PyList_GET_ITEM(lines, i)) + index;
            memcpy(destinations[i], s, step);
            destinations[i] += step;
        }
        index += step;
    }

    coordinates->shape[0] = n;
    coordinates->shape[1] = m;

    result = Py_BuildValue("OO", sequences, coordinates);
exit:
    Py_XDECREF(sequences);
    Py_XDECREF(coordinates);
    if (indices) PyMem_Free(indices);
    if (steps) PyMem_Free(steps);
    if (destinations) PyMem_Free(destinations);
    return result;
}

static char parse_printed_alignment__doc__[] =
"Infer the coordinates from a printed alignment."
"\n"
"This method is primarily employed by the alignment parsers in Bio.Align,\n"
"though it may be useful for other purposes.\n"
"\n"
"For an alignment consisting of N sequences, printed as N lines with the same\n"
"number of columns, where gaps are represented by dashes, this method\n"
"calculates the sequence coordinates that define the alignment.\n"
"This function returns a tuple of the sequences and the coordinates.\n"
"The sequences object is a tuple of bytes objects containing the sequence\n"
"data after removing the dashes. The coordinates object returned supports\n"
"the buffer protocol, allowing the coordinates to be converted into a NumPy\n"
"array of integers without copying.\n";


static struct PyMethodDef methods[] = {
    {"parse_printed_alignment",
     (PyCFunction)parse_printed_alignment,
     METH_VARARGS,
     parse_printed_alignment__doc__
    },                             
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

/* Module initialization function */
static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_parser",
    .m_doc = "fast C implementation of utility functions for alignment parsing",
    .m_size = -1,
    .m_methods = methods,
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
