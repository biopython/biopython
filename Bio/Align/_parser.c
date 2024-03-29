/* Copyright 2024 by Michiel de Hoon.  All rights reserved.
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */


#include <Python.h>
#include <string.h>


typedef struct {
    PyObject_HEAD
    long *data;
    Py_ssize_t shape[2];
} Coordinates;

/* Destructor function */
static void
Coordinates_dealloc(Coordinates *self)
{
    if (self->data != NULL) {
        PyMem_Free(self->data);
    }
    Py_TYPE(self)->tp_free((PyObject *)self);
}

/* Buffer protocol getbuffer function */
static int
Coordinates_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
    Coordinates *self = (Coordinates *)obj;

    if (view == NULL) {
        PyErr_SetString(PyExc_BufferError, "NULL view in getbuffer");
        return -1;
    }

    view->obj = obj;
    view->buf = self->data;
    view->len = self->shape[0] * self->shape[1] * sizeof(long);
    view->readonly = 0;
    view->itemsize = sizeof(long);
    view->format = "l";
    view->ndim = 2;
    view->shape = self->shape;
    view->strides = NULL;
    view->suboffsets = NULL;
    view->internal = NULL;

    Py_INCREF(obj);

    return 0;
}

static PyBufferProcs Coordinates_as_buffer = {
    .bf_getbuffer = Coordinates_getbuffer,
};

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


static PyTypeObject CoordinatesType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_parser.Coordinates",
    .tp_doc = "Coordinates objects",
    .tp_basicsize = sizeof(Coordinates),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_dealloc = (destructor)Coordinates_dealloc,
    .tp_as_buffer = &Coordinates_as_buffer,
    .tp_getset = Coordinates_getset,
};

static int
converter(PyObject* argument, void* pointer)
{
    Py_buffer* lines;
    Py_ssize_t i, n;
    PyObject* item;
    PyObject** items;

    argument = PySequence_Fast(argument,
                               "argument must be a tuple, list, or any other "
                               "object supporting the sequence protocol");
    if (!argument) return 0;
    n = PySequence_Fast_GET_SIZE(argument);
    items = PySequence_Fast_ITEMS(argument);
    lines = PyMem_Calloc(n+1, sizeof(Py_buffer));  /* one more as sentinel */
    if (!lines) {
        Py_DECREF(argument);
        return 0;
    }
    *((Py_buffer**)pointer) = lines;
    for (i = 0; i < n; i++) {
        item = items[i];
        if (PyObject_GetBuffer(item, &lines[i], PyBUF_CONTIG_RO) == 0) {
            if (lines[i].itemsize == 1 && lines[i].ndim == 1) continue;
            PyErr_Format(PyExc_ValueError,
                         "line [%zi] does not contain a one-dimensional array "
                         "of single bytes.", i);
        }
        else {
            PyErr_Clear();
            if (PyUnicode_Check(item)
             && PyUnicode_READY(item) == 0
             && PyUnicode_IS_ASCII(item)) {
                PyBuffer_FillInfo(&lines[i],
                                  item,
                                  PyUnicode_DATA(item),
                                  PyUnicode_GET_LENGTH(item),
                                  1,
                                  PyBUF_CONTIG_RO);
                continue;
            }
            else {
                PyErr_Format(PyExc_TypeError,
                             "line [%zi] is neither a bytes-like object "
                             "nor an ASCII string", i);
            }
        }
        n = i;
        for (i = 0; i < n; i++) PyBuffer_Release(&lines[i]);
        Py_DECREF(argument);
        return 0;
    }
    Py_DECREF(argument);
    return 1;
}


static PyObject*
parse_printed_alignment(PyObject* module, PyObject* args)
{
    Py_buffer* lines;
    Py_ssize_t i, j, k, n, m, p;
    const char* s;
    char** destination = NULL;
    Py_ssize_t index, previous;
    long step;
    Py_ssize_t* indices = NULL;
    Py_ssize_t length;
    Py_ssize_t* lengths = NULL;
    Coordinates *coordinates = NULL;
    PyObject* sequences = NULL;
    PyObject* sequence;

    PyObject* result = NULL;

    if (!PyArg_ParseTuple(args, "O&:parse_printed_alignment", converter, &lines))
        return NULL;

    coordinates = (Coordinates *)PyType_GenericAlloc(&CoordinatesType, 0);
    if (coordinates == NULL) {
        PyMem_Free(lines);
        return NULL;
    }
    if (lines[0].buf == NULL) {
        coordinates->shape[0] = 0;
        coordinates->shape[1] = 0;
        PyMem_Free(lines);
        return Py_BuildValue("[]O", coordinates);
    }
    m = lines[0].shape[0];
    n = 1;
    for (i = 1; lines[i].buf; i++) {
        if (lines[i].shape[0] != m) {
            PyErr_Format(PyExc_ValueError,
                "all lines must have the same length "
                "(line [0] has length %zi, line [%zi] has length %zi).",
                m, i, lines[i].shape[0]);
            goto exit;
        }
        n++;
    }
    indices = PyMem_Calloc(n, sizeof(Py_ssize_t));
    if (!indices) goto exit;
    lengths = PyMem_Calloc(n, sizeof(Py_ssize_t));
    if (!lengths) goto exit;
    destination = PyMem_Malloc(n*sizeof(char*));
    if (!destination) goto exit;
    p = 0;
    while (1) {
        index = m;
        for (i = 0; i < n; i++) if (indices[i] < index) index = indices[i];
        p++;
        if (index == m) break;
        for (i = 0; i < n; i++) {
            if (indices[i] == index) {
                s = lines[i].buf + index;
                if (*s == '-') {
                    for (j = index+1; j < m; j++) if (*(++s) != '-') break;
                }
                else {
                    for (j = index+1; j < m; j++) if (*(++s) == '-') break;
                    lengths[i] += j - index;
                }
                indices[i] = j;
            }
        }
    }

    sequences = PyList_New(n);
    if (sequences == NULL) goto exit;
    for (i = 0; i < n ; i++) {
        if (PyUnicode_Check(lines[i].obj)) {
            sequence = PyUnicode_New(lengths[i], 127);
            if (!sequence) goto exit;
            destination[i] = PyUnicode_DATA(sequence);
        }
        else {
            sequence = PyBytes_FromStringAndSize(NULL, lengths[i]);
            if (!sequence) goto exit;
            destination[i] = PyBytes_AS_STRING(sequence);
        }
        PyList_SET_ITEM(sequences, i, sequence);
    }
    coordinates = (Coordinates *)PyType_GenericAlloc(&CoordinatesType, 0);
    if (coordinates == NULL) goto exit;
    coordinates->data = PyMem_Malloc(n*p*sizeof(long));
    if (!coordinates->data) goto exit;
    for (i = 0; i < n; i++) coordinates->data[i*p] = 0;
    coordinates->shape[0] = n;
    coordinates->shape[1] = p;

    memset(indices, '\0', n * sizeof(Py_ssize_t));
    memset(lengths, '\0', n * sizeof(Py_ssize_t));
    k = 1;
    index = 0;
    do {
        previous = index;
        for (i = 0; i < n; i++) {
            if (index == indices[i]) {
                length = lengths[i];
                s = lines[i].buf + index;
                if (*s == '-') {
                    for (j = index+1; j < m; j++) if (*(++s) != '-') break;
                }
                else {
                    for (j = index+1; j < m; j++) if (*(++s) == '-') break;
                    s = lines[i].buf + index;
                    memcpy(destination[i]+length, s, j - index);
                    lengths[i] += j - index;
                }
                indices[i] = j;
            }
        }
        index = m;
        for (i = 0; i < n; i++) if (indices[i] < index) index = indices[i];
        step = index - previous;
        for (i = 0; i < n; i++) {
            s = lines[i].buf;
            if (s[previous] == '-') {
                coordinates->data[i*p+k] = coordinates->data[i*p+k-1];
            }
            else {
                coordinates->data[i*p+k] = coordinates->data[i*p+k-1] + step;
            }
        }
        k++;
    }
    while (index < m);

    result = Py_BuildValue("OO", sequences, coordinates);
exit:
    Py_XDECREF(sequences);
    Py_XDECREF(coordinates);
    if (indices) PyMem_Free(indices);
    if (lengths) PyMem_Free(lengths);
    if (destination) PyMem_Free(destination);
    for (i = 0; i < n; i++) PyBuffer_Release(&lines[i]);
    PyMem_Free(lines);
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
