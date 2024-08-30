/* Copyright 2024 by Michiel de Hoon.  All rights reserved.
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */


#include <Python.h>
#include <stdbool.h>


static PyTypeObject ParserType;

typedef struct {
    PyObject_HEAD
    Py_uintptr_t** data;
    /* Array of n Py_uintptr_t* pointers; each pointer points to an array of
       Py_uintptr_t of variable length.  This array contains, for each sequence,
       the positions in the printed alignment at which a letter is followed by
       a gap, or vice-versa.  If the first character is a gap, then the first
       column in data is 0.
     */
    Py_ssize_t n;  /* number of sequences in the alignment */
    Py_ssize_t m;  /* number of columns in the printed alignment */
    Py_ssize_t k;  /* number of columns of the coordinates array */
    char eol;      /* optional end-of-line character (set to '\n' by default) */
} Parser;

static void
Parser_dealloc(Parser *self)
{
    Py_ssize_t i;
    const Py_ssize_t n = self->n;
    Py_uintptr_t** data = self->data;
    if (data) {
        for (i = 0; i < n; i++) {
            if (data[i] == NULL) break;
            PyMem_Free(data[i]);
        }
        PyMem_Free(data);
    }
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static int
array_converter(PyObject* argument, void* pointer)
/* Check if the NumPy array received in argument has the correct shape and
 * data type, and fill the Py_buffer structure referred to by pointer.
 */
{
    Py_buffer* view = pointer;
    Parser* self;

    if (!PyObject_TypeCheck(view->obj, &ParserType)) {
        PyErr_SetString(PyExc_RuntimeError,
                        "expected an object of the PrintedAlignmentParser class");
        return 0;
    }

    self = (Parser*) view->obj;

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
    else if (view->itemsize != sizeof(int64_t)) {
        PyErr_Format(PyExc_RuntimeError,
                    "buffer has unexpected item byte size "
                    "(%ld, expected %ld)", view->itemsize, sizeof(int64_t));
    }
    else return 1;  /* return status 1 to indicate a successful converstion */

    PyBuffer_Release(view);
    return 0;  /* return status 0 to indicate that converstion failed */
}

PyDoc_STRVAR(
    Parser_feed__doc__,
    "feed(self, line, offset=0)\n"
    "--\n"
    "\n"
    "Feed one line of the printed alignment into the parser.\n"
    "\n"
    "The line must be a bytes object. The parser will read from line\n"
    "until it finds the end-of-line character (defined by self->eol)\n"
    "or a null character.\n"
    "\n"
    "The parser skips the first offset bytes.\n"
    "\n"
    "Any dashes in line are interpreted as gaps.\n"
    "This method finds the gap locations and stores them in self.\n"
    "\n"
    "The return value is the tuple (nbytes, sequence), in which\n"
    " - nbytes is the number of bytes read from line (not counting\n"
    "   the end-of-line character); this is equal to the number of\n"
    "   columns in the printed alignment.\n"
    " - sequence is a bytes object with the contents of line after\n"
    "   removal of the dashes; sequence can be used to create the\n"
    "   ungapped sequence object stored in the sequences attribute\n"
    "   of the Alignment object.");

static PyObject*
Parser_feed(Parser* self, PyObject* args, PyObject *kwds)
{
    PyObject* line = NULL;
    PyObject* sequence;
    PyObject* result;
    const char* buffer;
    const char* s;
    char* d;
    const char eol = self->eol;
    Py_ssize_t n = self->n;
    Py_ssize_t m = self->m;
    Py_ssize_t size = 2;
    Py_ssize_t i = 0;
    Py_ssize_t p = 0;
    Py_ssize_t offset = 0;
    Py_ssize_t start, end, step;
    Py_uintptr_t** data;
    Py_uintptr_t* row;
    char c;
    bool gap = false;

    if (!PyArg_ParseTuple(args, "S|n:feed", &line, &offset)) return NULL;

    buffer = PyBytes_AS_STRING(line) + offset;

    s = buffer;
    row = PyMem_Malloc(size*sizeof(Py_uintptr_t));
    if (!row) return NULL;
    if (*s == '-') row[i++] = 0;

    data = PyMem_Realloc(self->data, (n+1)*size*sizeof(Py_uintptr_t*));
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
            start = s - buffer;
            do c = *(++s); while (c != '-' && c != eol && c != '\0');
            end = s - buffer;
            step = end - start;
            p += step;
        }

        if (i == size) {
            size *= 2;
            row = PyMem_Realloc(row, size*sizeof(Py_uintptr_t));
            if (!row) {
                PyMem_Free(data[n]);
                return NULL;
            }
            data[n] = row;
        }
        row[i++] = s - buffer;
    }
    row = PyMem_Realloc(row, i*sizeof(Py_uintptr_t));
    if (!row) {
        PyMem_Free(data[n]);
        return NULL;
    }
    data[n] = row;
    m = s - buffer;
    if (n == 0) self->m = m;
    else if (buffer + m != s) {
        PyErr_Format(PyExc_ValueError,
                     "line has length %zd (expected %zd)", m, self->m);
        PyMem_Free(row);
        return NULL;
    }

    n++;
    self->n = n;

    sequence = PyBytes_FromStringAndSize(NULL, p);
    if (!sequence) return NULL;
    d = PyBytes_AS_STRING(sequence);
    end = 0;
    s = buffer;
    p = 0;
    if (row[p] == 0) {
        gap = true;
        p++;
    }
    for ( ; p < i; p++) {
        start = end;
        end = row[p];
        step = end - start;
        gap = !gap;
        if (gap) {
            s = buffer + start;
            memcpy(d, s, step);
            d += step;
        }
    }
    *d = '\0';

    result = Py_BuildValue("nN", m, sequence);
    if (result == NULL) Py_DECREF(sequence);
    return result;
}

PyDoc_STRVAR(
    Parser_fill__doc__,
    "fill(self, arr)\n"
    "--\n"
    "\n"
    "Fill in the coordinates array based on the alignment lines fed\n"
    "to the parser so far.\n"
    "\n"
    "The argument arr must be a 2D numpy array of data type int,\n"
    "with the number of rows equal to the number of lines fed into\n"
    "the parser so far, and the number of columns equal to the number\n"
    "of columns needed to store the coordinates array.\n"
    "The appropriate number of rows and columns can be obtained in\n"
    "advance using the self.shape attribute.\n"
    "\n"
    "This method stores the alignment coordinates in arr, and returns\n"
    "None.\n"
);

static PyObject*
Parser_fill(Parser* self, PyObject* args)
{
    Py_buffer view;
    Py_ssize_t i, j, k, n, m, p;
    Py_ssize_t start;
    Py_ssize_t step;
    Py_ssize_t end;
    Py_ssize_t* starts = NULL;
    Py_uintptr_t** data = NULL;
    bool* gaps = NULL;
    int64_t* buffer;

    n = self->n;
    if (n == 0) Py_RETURN_NONE;

    view.obj = (PyObject*) self;

    if (!PyArg_ParseTuple(args, "O&:fill", array_converter, &view))
        return NULL;

    buffer = view.buf;
    k = view.shape[1];
    if (n != view.shape[0]) {
        PyErr_Format(PyExc_ValueError,
                     "expected an array with %zd rows (found %zd rows)",
                     n, view.shape[0]);
        goto exit;
    }
    for (i = 0; i < n; i++) buffer[i*k] = 0;

    m = self->m;

    starts = PyMem_Calloc(n, sizeof(Py_ssize_t));
    if (!starts) goto exit;

    gaps = PyMem_Malloc(n * sizeof(bool));
    if (!gaps) goto exit;

    data = PyMem_Malloc(n * sizeof(Py_ssize_t*));
    if (!data) goto exit;

    for (i = 0; i < n; i++) {
        data[i] = self->data[i];
        if (data[i][0] == 0) {
            gaps[i] = true;
            data[i]++;
        }
        else {
            gaps[i] = false;
        }
    }

    j = 0;
    start = 0;
    do {
        j++;
        for (i = 0; i < n; i++) {
            if (start == starts[i]) starts[i] = *data[i];
        }
        end = m;
	for (i = 0; i < n; i++) if (starts[i] < end) end = starts[i];
        step = end - start;
        for (i = 0; i < n; i++) {
            p = i*k+j;
            if (gaps[i] == true) buffer[p] = buffer[p-1];
            else buffer[i*k+j] = (int64_t) (buffer[p-1] + step);
            if (end == starts[i]) {
                data[i]++;
                gaps[i] = !gaps[i];
            }
        }
        start = end;
    }
    while (start < m);

exit:
    PyBuffer_Release(&view);
    if (starts) PyMem_Free(starts);
    if (data) PyMem_Free(data);
    if (gaps) PyMem_Free(gaps);
    Py_RETURN_NONE;
}

PyDoc_STRVAR(
    Parser_shape__doc__,
    "Return the required shape of the coordinates array.\n"
    "Typically, this attribute is used before calling the fill method.\n"
);

static PyObject*
Parser_get_shape(Parser* self, void* closure)
{
    Py_ssize_t i;
    Py_ssize_t index;
    Py_ssize_t min_index;
    Py_uintptr_t** data = self->data;
    const Py_ssize_t n = self->n;
    const Py_ssize_t m = self->m;
    Py_ssize_t k = 1;

    if (n > 0) {
        data = PyMem_Malloc(n * sizeof(Py_uintptr_t*));
        if (!data) return NULL;
        memcpy(data, self->data, n*sizeof(Py_uintptr_t*));
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
    return Py_BuildValue("nn", n, k);
}

static PyGetSetDef Parser_getset[] = {
    {"shape", (getter)Parser_get_shape,
     NULL,
     Parser_shape__doc__,
     NULL,
    },
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

static PyMethodDef Parser_methods[] = {
    {"feed",
     (PyCFunction)Parser_feed,
     METH_VARARGS,
     Parser_feed__doc__,
    },
    {"fill",
     (PyCFunction)Parser_fill,
     METH_VARARGS,
     Parser_fill__doc__,
    },
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

static PyObject*
Parser_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    char eol = '\n';  /* end-of-line character */
    static char *kwlist[] = {"eol", NULL};

    Parser *self;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|c", kwlist, &eol))
        return NULL;

    self = (Parser *)type->tp_alloc(type, 0);
    if (!self) return NULL;
    self->eol = eol;
    self->n = 0;
    return (PyObject *)self;
}

PyDoc_STRVAR(
    Parser__doc__,
    "PrintedAlignmentParser(eol=b'\n')\n"
    "--\n"
    "\n"
    "Create a fast parser for printed alignments.\n"
    "\n"
    "The argument eol must be a bytes object of length 1, and specifies\n"
    "the end-of-line character, defaulting to '\\n'"
    "\n"
    "As an example, to parse this printed alignment:\n"
    "\n"
    "ACCGGGTTTT\n"
    "AC-GAG--TT\n"
    "AC--AG--TT\n"
    "\n"
    "use\n"
    "\n"
    ">>> parser = PrintedAlignmentParser()\n"
    ">>> nbytes1, seq1 = parser.feed(b'ACCGGGTTTT')\n"
    ">>> nbytes2, seq2 = parser.feed(b'AC-GAG--TT')\n"
    ">>> nbytes3, seq3 = parser.feed(b'AC--AG--TT')\n"
    ">>> nbytes1\n"
    "10\n"
    ">>> nbytes2\n"
    "10\n"
    ">>> nbytes3\n"
    "10\n"
    ">>> seq1\n"
    "b'ACCGGGTTTT'\n"
    ">>> seq2\n"
    "b'ACGAGTT'\n"
    ">>> seq3\n"
    "b'ACAGTT'\n"
    ">>> parser.shape\n"
    "(3, 7)\n"
    ">>> import numpy as np\n"
    ">>> coordinates = np.zeros((3, 7), int)\n"
    ">>> parser.fill(coordinates)\n"
    ">>> coordinates\n"
    "array([[ 0,  2,  3,  4,  6,  8, 10],\n"
    "       [ 0,  2,  2,  3,  5,  5,  7],\n"
    "       [ 0,  2,  2,  2,  4,  4,  6]])\n"
    ">>> from Bio.Align import Alignment\n"
    ">>> from Bio.Seq import Seq\n"
    ">>> sequences = (Seq(seq1), Seq(seq2), Seq(seq3))\n"
    ">>> alignment = Alignment(sequences, coordinates)\n"
    ">>> print(alignment)\n"
    "                  0 ACCGGGTTTT 10\n"
    "                  0 AC-GAG--TT  7\n"
    "                  0 AC--AG--TT  6\n"
    "\n"
);

static PyTypeObject ParserType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_aligncore.PrintedAlignmentParser",
    .tp_doc = Parser__doc__,
    .tp_basicsize = sizeof(Parser),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_dealloc = (destructor)Parser_dealloc,
    .tp_methods = Parser_methods,
    .tp_getset = Parser_getset,
    .tp_new = (newfunc)Parser_new,
};

/* Module initialization function */
static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_aligncore",
    .m_doc = "fast C implementation of a parser for printed alignments; for internal use.",
    .m_size = -1,
};

PyMODINIT_FUNC
PyInit__aligncore(void)
{
    PyObject *module;

    if (PyType_Ready(&ParserType) < 0)
        return NULL;

    module = PyModule_Create(&moduledef);
    if (module == NULL)
        return NULL;

    Py_INCREF(&ParserType);
    PyModule_AddObject(module, "PrintedAlignmentParser", (PyObject *)&ParserType);
    return module;
}
