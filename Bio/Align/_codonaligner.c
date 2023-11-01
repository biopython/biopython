/* Copyright 2023 by Michiel de Hoon.  All rights reserved.
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */


#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "float.h"


#define FRAMESHIFT_MINUS_TWO 0x1
#define FRAMESHIFT_MINUS_ONE 0x2
#define FRAMESHIFT_NONE 0x4
#define FRAMESHIFT_PLUS_ONE 0x8
#define FRAMESHIFT_PLUS_TWO 0x10

#define DONE 6
#define NONE 7

#define OVERFLOW_ERROR -1
#define MEMORY_ERROR -2

#define SAFE_ADD(t, s) \
{   if (s != OVERFLOW_ERROR) { \
        term = t; \
        if (term > PY_SSIZE_T_MAX - s) s = OVERFLOW_ERROR; \
        else s += term; \
    } \
}


typedef struct {
    unsigned char trace : 5;
    unsigned char path : 3;
} Trace;

typedef struct {
    PyObject_HEAD
    Trace** M;
    int nA;
    int nB;
    Py_ssize_t length;
} PathGenerator;

static PyObject*
PathGenerator_create_path(PathGenerator* self, int j) {
    PyObject* tuple;
    PyObject* target_row;
    PyObject* query_row;
    PyObject* value;
    int path;
    int i = 0;
    int k, l;
    int n = 1;
    int direction = 0;
    Trace** M = self->M;

    k = i;
    l = j;
    while (1) {
        path = M[k][l].path;
        if (!path) break;
        if (path % 3 != 0) {
            n += 2;
            direction = 3;
        }
        else if (path != direction) {
            n++;
            direction = path;
        }
        l += path;
        k++;
    }

    direction = 0;
    tuple = PyTuple_New(2);
    if (!tuple) return NULL;
    target_row = PyTuple_New(n);
    query_row = PyTuple_New(n);
    PyTuple_SET_ITEM(tuple, 0, target_row);
    PyTuple_SET_ITEM(tuple, 1, query_row);

    if (target_row && query_row) {
        k = 0;
        while (1) {
            path = M[i][j].path;
            if (path % 3 != 0) {
                value = PyLong_FromLong(i);
                if (!value) break;
                PyTuple_SET_ITEM(target_row, k, value);
                value = PyLong_FromLong(j);
                if (!value) break;
                PyTuple_SET_ITEM(query_row, k, value);
                k++;
                j += path - 3;
                value = PyLong_FromLong(i);
                if (!value) break;
                PyTuple_SET_ITEM(target_row, k, value);
                value = PyLong_FromLong(j);
                if (!value) break;
                PyTuple_SET_ITEM(query_row, k, value);
                k++;
                direction = 3;
            }
            else if (path != direction) {
                value = PyLong_FromLong(i);
                if (!value) break;
                PyTuple_SET_ITEM(target_row, k, value);
                value = PyLong_FromLong(j);
                if (!value) break;
                PyTuple_SET_ITEM(query_row, k, value);
                k++;
                direction = path;
                if (path == 0) return tuple;
            }
            i++;
            j += 3;
        }
    }
    Py_DECREF(tuple); /* all references were stolen */
    return PyErr_NoMemory();
}

static Py_ssize_t PathGenerator_length(PathGenerator* self) {
    Py_ssize_t count = self->length;
    if (count == 0) {
        int i;
        int j;
        int trace;
        const int nA = self->nA;
        const int nB = self->nB;
        Trace** M = self->M;
        Py_ssize_t term;
        Py_ssize_t* counts1 = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
        Py_ssize_t* counts2 = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
        if (counts1 == NULL || counts2 == NULL) {
            PyErr_NoMemory();
            count = MEMORY_ERROR;
            goto exit;
        }
        for (j = 0; j <= nB; j++) counts2[j] = 1;
        for (i = 1; i <= nA; i++) {
            memcpy(counts1, counts2, (nB+1)*sizeof(Py_ssize_t));
            for (j = 0; j <= nB; j++) {
                trace = M[i][j].trace;
                count = 0;
                if (trace & FRAMESHIFT_MINUS_TWO) SAFE_ADD(counts1[j-1], count);
                if (trace & FRAMESHIFT_MINUS_ONE) SAFE_ADD(counts1[j-2], count);
                if (trace & FRAMESHIFT_NONE) SAFE_ADD(counts1[j-3], count);
                if (trace & FRAMESHIFT_PLUS_ONE) SAFE_ADD(counts1[j-4], count);
                if (trace & FRAMESHIFT_PLUS_TWO) SAFE_ADD(counts1[j-5], count);
                counts2[j] = count;
            }
        }
        count = 0;
        for (j = 0; j <= nB; j++) count += counts2[j];
        self->length = count;
exit:
        PyMem_Free(counts1);
        PyMem_Free(counts2);
    }
    if (count == OVERFLOW_ERROR)
        PyErr_Format(PyExc_OverflowError,
                     "number of optimal alignments is larger than %zd",
                     PY_SSIZE_T_MAX);
    return count;
}

static void
PathGenerator_dealloc(PathGenerator* self)
{
    int i;
    const int nA = self->nA;
    Trace** M = self->M;
    if (M) {
        for (i = 0; i <= nA; i++) {
            if (!M[i]) break;
            PyMem_Free(M[i]);
        }
        PyMem_Free(M);
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
PathGenerator_next(PathGenerator* self)
{
    int i = 0;
    int j;
    int path;
    int trace = 0;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M;

    if (M[0][0].path == DONE) return NULL;
    for (j = 0; j <= nB; j++) {
        path = M[0][j].path;
        if (path) {
            /* We already have a path. Prune the path to see if there are
             * any alternative paths. */
            M[0][j].path = 0;
            while (1) {
                j += path;
                trace = M[i+1][j].trace;
                if (path == 1 && trace & FRAMESHIFT_MINUS_ONE) {
                    path = 2;
                    break;
                }
                else if (path <= 2 && trace & FRAMESHIFT_NONE) {
                    path = 3;
                    break;
                }
                else if (path <= 3 && trace & FRAMESHIFT_PLUS_ONE) {
                    path = 4;
                    break;
                }
                else if (path <= 4 &&  trace & FRAMESHIFT_PLUS_TWO) {
                    path = 5;
                    break;
                }
                i++;
                path = M[i][j].path;
                if (!path)
                    /* we reached the end of the alignment without finding
                     * an alternative path */
                    break;
            }
            if (path) {
                j -= path;
                M[i][j].path = path;
            }
            break;
        }
    }
    if (path == 0) {
        /* Find the next end point. */
        if (i == 0) {
            j = 0;
            i = nA;
        }
        else j++;
        for ( ; j <= nB; j++) {
            if (M[nA][j].trace) break;
        }
        if (j > nB) {
            /* No further end points, so we are done. */
            M[0][0].path = DONE;
            return NULL;
        }
    }
    /* Follow the traceback until we reach the origin. */
    while (1) {
        trace = M[i][j].trace;
        if (trace & FRAMESHIFT_MINUS_TWO) path = 1;
        else if (trace & FRAMESHIFT_MINUS_ONE) path = 2;
        else if (trace & FRAMESHIFT_NONE) path = 3;
        else if (trace & FRAMESHIFT_PLUS_ONE) path = 4;
        else if (trace & FRAMESHIFT_PLUS_TWO) path = 5;
        else break;
        j -= path;
        i--;
        M[i][j].path = path;
    }
    return PathGenerator_create_path(self, j);
}

static const char PathGenerator_reset__doc__[] = "reset the iterator";

static PyObject*
PathGenerator_reset(PathGenerator* self)
{
    Trace** M = self->M;
    if (M[0][0].path != NONE) M[0][0].path = 0;
    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef PathGenerator_methods[] = {
    {"reset",
     (PyCFunction)PathGenerator_reset,
     METH_NOARGS,
     PathGenerator_reset__doc__
    },
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

static PySequenceMethods PathGenerator_as_sequence = {
    (lenfunc)PathGenerator_length,  /* sq_length */
    NULL,                           /* sq_concat */
    NULL,                           /* sq_repeat */
    NULL,                           /* sq_item */
    NULL,                           /* sq_ass_item */
    NULL,                           /* sq_contains */
    NULL,                           /* sq_inplace_concat */
    NULL,                           /* sq_inplace_repeat */
};

static PyTypeObject PathGenerator_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "Path generator",               /* tp_name */
    sizeof(PathGenerator),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)PathGenerator_dealloc,  /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    0,                              /* tp_repr */
    0,                              /* tp_as_number */
    &PathGenerator_as_sequence,     /* tp_as_sequence */
    0,                              /* tp_as_mapping */
    0,                              /* tp_hash */
    0,                              /* tp_call */
    0,                              /* tp_str */
    0,                              /* tp_getattro */
    0,                              /* tp_setattro */
    0,                              /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,             /* tp_flags */
    0,                              /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    PyObject_SelfIter,              /* tp_iter */
    (iternextfunc)PathGenerator_next,      /* tp_iternext */
    PathGenerator_methods,          /* tp_methods */
};

typedef struct {
    PyObject_HEAD
    double match;
    double mismatch;
    double epsilon;
    char wildcard;
    double frameshift_minus_two_score;
    double frameshift_minus_one_score;
    double frameshift_plus_one_score;
    double frameshift_plus_two_score;
} Aligner;


static int
Aligner_init(Aligner *self, PyObject *args, PyObject *kwds)
{
    self->match = 1.0;
    self->mismatch = 0.0;
    self->epsilon = 1.e-6;
    self->wildcard = 'X';
    self->frameshift_minus_two_score = -3.0;
    self->frameshift_minus_one_score = -3.0;
    self->frameshift_plus_one_score = -3.0;
    self->frameshift_plus_two_score = -3.0;
    return 0;
}


static void
Aligner_dealloc(Aligner* self)
{
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject*
Aligner_repr(Aligner* self)
{
  const char text[] = "Codon aligner, implementing a dynamic programming algorithm to align a nucleotide sequence to an amino acid sequence";
  return PyUnicode_FromString(text);
}

static PyObject*
Aligner_str(Aligner* self)
{
    PyObject* match;
    PyObject* mismatch = NULL;
    PyObject* text = NULL;
    PyObject* frameshift_minus_two_score = NULL;
    PyObject* frameshift_minus_one_score = NULL;
    PyObject* frameshift_plus_one_score = NULL;
    PyObject* frameshift_plus_two_score = NULL;
    match = PyFloat_FromDouble(self->match);
    if (match == NULL) goto exit;
    mismatch = PyFloat_FromDouble(self->mismatch);
    if (mismatch == NULL) goto exit;
    frameshift_minus_two_score = PyFloat_FromDouble(self->frameshift_minus_two_score);
    if (frameshift_minus_two_score == NULL) goto exit;
    frameshift_minus_one_score = PyFloat_FromDouble(self->frameshift_minus_one_score);
    if (frameshift_minus_one_score == NULL) goto exit;
    frameshift_plus_one_score = PyFloat_FromDouble(self->frameshift_plus_one_score);
    if (frameshift_plus_one_score == NULL) goto exit;
    frameshift_plus_two_score = PyFloat_FromDouble(self->frameshift_plus_two_score);
    if (frameshift_plus_two_score == NULL) goto exit;
    text = PyUnicode_FromFormat("Codon aligner with parameters\n"
                                "  wildcard: '%c'\n"
                                "  match_score: %S\n"
                                "  mismatch_score: %S\n"
                                "  frameshift_minus_two_score: %S\n"
                                "  frameshift_minus_one_score: %S\n"
                                "  frameshift_plus_one_score: %S\n"
                                "  frameshift_plus_two_score: %S\n",
                                self->wildcard,
                                match,
                                mismatch,
                                frameshift_minus_two_score,
                                frameshift_minus_one_score,
                                frameshift_plus_one_score,
                                frameshift_plus_two_score);
exit:
    Py_XDECREF(match);
    Py_XDECREF(mismatch);
    Py_XDECREF(frameshift_minus_two_score);
    Py_XDECREF(frameshift_minus_one_score);
    Py_XDECREF(frameshift_plus_one_score);
    Py_XDECREF(frameshift_plus_two_score);
    return text;
}

static char Aligner_match_score__doc__[] = "match score";

static PyObject*
Aligner_get_match_score(Aligner* self, void* closure)
{
    return PyFloat_FromDouble(self->match);
}

static int
Aligner_set_match_score(Aligner* self, PyObject* value, void* closure)
{
    const double match = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_ValueError, "invalid match score");
        return -1;
    }
    self->match = match;
    return 0;
}

static char Aligner_mismatch_score__doc__[] = "mismatch score";

static PyObject*
Aligner_get_mismatch_score(Aligner* self, void* closure)
{
    return PyFloat_FromDouble(self->mismatch);
}

static int
Aligner_set_mismatch_score(Aligner* self, PyObject* value, void* closure)
{
    const double mismatch = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_ValueError, "invalid mismatch score");
        return -1;
    }
    self->mismatch = mismatch;
    return 0;
}

static char Aligner_epsilon__doc__[] = "roundoff epsilon";

static PyObject*
Aligner_get_epsilon(Aligner* self, void* closure)
{
    return PyFloat_FromDouble(self->epsilon);
}

static int
Aligner_set_epsilon(Aligner* self, PyObject* value, void* closure)
{   const double epsilon = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->epsilon = epsilon;
    return 0;
}

static PyObject*
Aligner_get_wildcard(Aligner* self, void* closure)
{
    return PyUnicode_FromFormat("%c", self->wildcard);
}

static int
Aligner_set_wildcard(Aligner* self, PyObject* value, void* closure)
{
    Py_UCS4 wildcard;
    if (!PyUnicode_Check(value)) {
        PyErr_SetString(PyExc_TypeError,
                        "wildcard should be a single ASCII character");
        return -1;
    }
    if (PyUnicode_READY(value) == -1) return -1;
    if (PyUnicode_GET_LENGTH(value) != 1) {
        PyErr_SetString(PyExc_ValueError,
                        "wildcard should be a single ASCII character");
        return -1;
    }
    wildcard = PyUnicode_READ_CHAR(value, 0);
    if (wildcard >= 256) {
        PyErr_SetString(PyExc_ValueError,
                        "wildcard should be a single ASCII character");
        return -1;
    }
    self->wildcard = (char)wildcard;
    return 0;
}

static char Aligner_wildcard__doc__[] = "wildcard character";

static PyObject*
Aligner_get_frameshift_minus_two_score(Aligner* self, void* closure)
{
    return PyFloat_FromDouble(self->frameshift_minus_two_score);
}

static int
Aligner_set_frameshift_minus_two_score(Aligner* self, PyObject* value, void* closure)
{
    const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->frameshift_minus_two_score = score;
    return 0;
}

static char Aligner_frameshift_minus_two_score__doc__[] = "score for a -2 frame shift";

static PyObject*
Aligner_get_frameshift_minus_one_score(Aligner* self, void* closure)
{
    return PyFloat_FromDouble(self->frameshift_minus_one_score);
}

static int
Aligner_set_frameshift_minus_one_score(Aligner* self, PyObject* value, void* closure)
{
    const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->frameshift_minus_one_score = score;
    return 0;
}

static char Aligner_frameshift_minus_one_score__doc__[] = "score for a -1 frame shift";

static PyObject*
Aligner_get_frameshift_plus_one_score(Aligner* self, void* closure)
{
    return PyFloat_FromDouble(self->frameshift_plus_one_score);
}

static int
Aligner_set_frameshift_plus_one_score(Aligner* self, PyObject* value, void* closure)
{
    const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->frameshift_plus_one_score = score;
    return 0;
}

static char Aligner_frameshift_plus_one_score__doc__[] = "score for a +1 frame shift";

static PyObject*
Aligner_get_frameshift_plus_two_score(Aligner* self, void* closure)
{
    return PyFloat_FromDouble(self->frameshift_plus_two_score);
}

static int
Aligner_set_frameshift_plus_two_score(Aligner* self, PyObject* value, void* closure)
{
    const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->frameshift_plus_two_score = score;
    return 0;
}

static char Aligner_frameshift_plus_two_score__doc__[] = "score for a +2 frame shift";


static PyObject*
Aligner_get_frameshift_minus_score(Aligner* self, void* closure)
{
    const double score = self->frameshift_minus_two_score;
    if (score != self->frameshift_minus_one_score) {
        PyErr_SetString(PyExc_ValueError, "-2 and -1 frame shift scores are different");
        return NULL;
    }
    return PyFloat_FromDouble(score);
}

static int
Aligner_set_frameshift_minus_score(Aligner* self, PyObject* value, void* closure)
{
    const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->frameshift_minus_one_score = score;
    self->frameshift_minus_two_score = score;
    return 0;
}

static char Aligner_frameshift_minus_score__doc__[] = "score for a negative frame shift";


static PyObject*
Aligner_get_frameshift_plus_score(Aligner* self, void* closure)
{
    const double score = self->frameshift_plus_two_score;
    if (score != self->frameshift_plus_one_score) {
        PyErr_SetString(PyExc_ValueError, "+2 and +1 frame shift scores are different");
        return NULL;
    }
    return PyFloat_FromDouble(score);
}

static int
Aligner_set_frameshift_plus_score(Aligner* self, PyObject* value, void* closure)
{
    const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->frameshift_plus_one_score = score;
    self->frameshift_plus_two_score = score;
    return 0;
}

static char Aligner_frameshift_plus_score__doc__[] = "score for a positive frame shift";


static PyObject*
Aligner_get_frameshift_two_score(Aligner* self, void* closure)
{
    const double score = self->frameshift_minus_two_score;
    if (score != self->frameshift_plus_two_score) {
        PyErr_SetString(PyExc_ValueError, "-2 and +2 frame shift scores are different");
        return NULL;
    }
    return PyFloat_FromDouble(score);
}

static int
Aligner_set_frameshift_two_score(Aligner* self, PyObject* value, void* closure)
{
    const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->frameshift_minus_two_score = score;
    self->frameshift_plus_two_score = score;
    return 0;
}

static char Aligner_frameshift_two_score__doc__[] = "score for a -2 or +2 frame shift";


static PyObject*
Aligner_get_frameshift_one_score(Aligner* self, void* closure)
{
    const double score = self->frameshift_minus_one_score;
    if (score != self->frameshift_plus_one_score) {
        PyErr_SetString(PyExc_ValueError, "-1 and +1 frame shift scores are different");
        return NULL;
    }
    return PyFloat_FromDouble(score);
}

static int
Aligner_set_frameshift_one_score(Aligner* self, PyObject* value, void* closure)
{
    const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->frameshift_minus_one_score = score;
    self->frameshift_plus_one_score = score;
    return 0;
}

static char Aligner_frameshift_one_score__doc__[] = "score for a -1 or +1 frame shift";


static PyObject*
Aligner_get_frameshift_score(Aligner* self, void* closure)
{
    const double score = self->frameshift_minus_one_score;
    if (score != self->frameshift_minus_two_score ||
        score != self->frameshift_plus_one_score ||
        score != self->frameshift_plus_two_score) {
        PyErr_SetString(PyExc_ValueError, "frame shift scores are different");
        return NULL;
    }
    return PyFloat_FromDouble(score);
}

static int
Aligner_set_frameshift_score(Aligner* self, PyObject* value, void* closure)
{
    const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->frameshift_minus_two_score = score;
    self->frameshift_minus_one_score = score;
    self->frameshift_plus_one_score = score;
    self->frameshift_plus_two_score = score;
    return 0;
}

static char Aligner_frameshift_score__doc__[] = "frame shift score";


static PyGetSetDef Aligner_getset[] = {
    {"match_score",
        (getter)Aligner_get_match_score,
        (setter)Aligner_set_match_score,
        Aligner_match_score__doc__, NULL},
    {"mismatch_score",
        (getter)Aligner_get_mismatch_score,
        (setter)Aligner_set_mismatch_score,
        Aligner_mismatch_score__doc__, NULL},
    {"match", /* synonym for match_score */
        (getter)Aligner_get_match_score,
        (setter)Aligner_set_match_score,
        Aligner_match_score__doc__, NULL},
    {"mismatch", /* synonym for mismatch_score */
        (getter)Aligner_get_mismatch_score,
        (setter)Aligner_set_mismatch_score,
        Aligner_mismatch_score__doc__, NULL},
    {"epsilon",
        (getter)Aligner_get_epsilon,
        (setter)Aligner_set_epsilon,
        Aligner_epsilon__doc__, NULL},
    {"wildcard",
        (getter)Aligner_get_wildcard,
        (setter)Aligner_set_wildcard,
        Aligner_wildcard__doc__, NULL},
    {"frameshift_minus_two_score",
        (getter)Aligner_get_frameshift_minus_two_score,
        (setter)Aligner_set_frameshift_minus_two_score,
        Aligner_frameshift_minus_two_score__doc__, NULL},
    {"frameshift_minus_one_score",
        (getter)Aligner_get_frameshift_minus_one_score,
        (setter)Aligner_set_frameshift_minus_one_score,
        Aligner_frameshift_minus_one_score__doc__, NULL},
    {"frameshift_plus_one_score",
        (getter)Aligner_get_frameshift_plus_one_score,
        (setter)Aligner_set_frameshift_plus_one_score,
        Aligner_frameshift_plus_one_score__doc__, NULL},
    {"frameshift_plus_two_score",
        (getter)Aligner_get_frameshift_plus_two_score,
        (setter)Aligner_set_frameshift_plus_two_score,
        Aligner_frameshift_plus_two_score__doc__, NULL},
    {"frameshift_minus_score",
        (getter)Aligner_get_frameshift_minus_score,
        (setter)Aligner_set_frameshift_minus_score,
        Aligner_frameshift_minus_score__doc__, NULL},
    {"frameshift_plus_score",
        (getter)Aligner_get_frameshift_plus_score,
        (setter)Aligner_set_frameshift_plus_score,
        Aligner_frameshift_plus_score__doc__, NULL},
    {"frameshift_one_score",
        (getter)Aligner_get_frameshift_one_score,
        (setter)Aligner_set_frameshift_one_score,
        Aligner_frameshift_one_score__doc__, NULL},
    {"frameshift_two_score",
        (getter)Aligner_get_frameshift_two_score,
        (setter)Aligner_set_frameshift_two_score,
        Aligner_frameshift_two_score__doc__, NULL},
    {"frameshift_score",
        (getter)Aligner_get_frameshift_score,
        (setter)Aligner_set_frameshift_score,
        Aligner_frameshift_score__doc__, NULL},
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

/* ----------------- alignment algorithms ----------------- */

#define COMPARE_SCORE (cA == wildcard || cB == wildcard) ? 0 : (cA == cB) ? match : mismatch


static const char Aligner_score__doc__[] = "calculates the alignment score";

static PyObject*
Aligner_score(Aligner* self, PyObject* args, PyObject* keywords)
{
    char* sA;
    char* sB[3];
    Py_ssize_t nA;
    Py_ssize_t nB;
    Py_ssize_t nB0;
    Py_ssize_t nB1;
    Py_ssize_t nB2;
    Py_buffer bA;
    Py_buffer bB0;
    Py_buffer bB1;
    Py_buffer bB2;
    PyObject* result = NULL;

    const double match = self->match;
    const double mismatch = self->mismatch;
    const char wildcard = self->wildcard;
    int i;
    int j;
    int div;
    int mod;
    char cA;
    char cB;
    const double frameshift_minus_two_score = self->frameshift_minus_two_score;
    const double frameshift_minus_one_score = self->frameshift_minus_one_score;
    const double frameshift_plus_one_score = self->frameshift_plus_one_score;
    const double frameshift_plus_two_score = self->frameshift_plus_two_score;
    double score;
    double temp;
    double* row = NULL;

    static char *kwlist[] = {"sA", "sB0", "sB1", "sB2", NULL};

    if(!PyArg_ParseTupleAndKeywords(args, keywords, "y*y*y*y*",
                                    kwlist, &bA, &bB0, &bB1, &bB2))
        return NULL;

    nA = bA.len;
    nB0 = bB0.len;
    nB1 = bB1.len;
    nB2 = bB2.len;
    if (nB1 == nB0 && nB2 == nB0) nB = 3 * nB0 + 2;
    else if (nB1 == nB0 && nB2 == nB0 - 1) nB = 3 * nB0 + 1;
    else if (nB1 == nB0 - 1 && nB2 == nB0 - 1) nB = 3 * nB0;
    else {
        PyErr_Format(PyExc_RuntimeError,
                     "unexpected length of buffers (%zd, %zd, %zd)",
                     nB0, nB1, nB2);
        PyBuffer_Release(&bA);
        PyBuffer_Release(&bB0);
        PyBuffer_Release(&bB1);
        PyBuffer_Release(&bB2);
        return NULL;
    }

    sA = bA.buf;
    sB[0] = bB0.buf;
    sB[1] = bB1.buf;
    sB[2] = bB2.buf;

    row = PyMem_Malloc((nB+1)*sizeof(double));
    if (!row) goto exit;
    memset(row, '\0', (nB+1)*sizeof(double));

    for (i = 1; i <= nA; i++) {
        cA = sA[i-1];
        for (j = (int)nB; j > 0; j--) {
            score = -DBL_MAX;
            if (j >= 3) {
                div = (j - 3) / 3;
                mod = (j - 3) % 3;
                cB = sB[mod][div];
                temp = row[j-1] + (COMPARE_SCORE) + frameshift_minus_two_score;
                if (temp > score) score = temp;
                temp = row[j-2] + (COMPARE_SCORE) + frameshift_minus_one_score;
                if (temp > score) score = temp;
                temp = row[j-3] + (COMPARE_SCORE);
                if (temp > score) score = temp;
            }
            if (j >= 4) {
                temp = row[j-4] + (COMPARE_SCORE) + frameshift_plus_one_score;
                if (temp > score) score = temp;
            }
            if (j >= 5) {
                temp = row[j-5] + (COMPARE_SCORE) + frameshift_plus_two_score;
                if (temp > score) score = temp;
            }
            row[j] = score;
        }
    }
    score = -DBL_MAX;
    for (j = 0; j <= nB; j++) {
        temp = row[j];
        if (temp > score) score = temp;
    }

    result = PyFloat_FromDouble(score);

exit:
    PyBuffer_Release(&bA);
    PyBuffer_Release(&bB0);
    PyBuffer_Release(&bB1);
    PyBuffer_Release(&bB2);
    PyMem_Free(row);
    if (result == NULL) {
        return PyErr_NoMemory();
    }
    return result;
}

static const char Aligner_align__doc__[] = "align two sequences";

static PyObject*
Aligner_align(Aligner* self, PyObject* args, PyObject* keywords)
{
    char* sA;
    char* sB[3];
    Py_ssize_t nA;
    Py_ssize_t nB;
    Py_ssize_t nB0;
    Py_ssize_t nB1;
    Py_ssize_t nB2;
    Py_buffer bA;
    Py_buffer bB0;
    Py_buffer bB1;
    Py_buffer bB2;
    PyObject* result = NULL;

    const double match = self->match;
    const double mismatch = self->mismatch;
    const char wildcard = self->wildcard;
    int i;
    int j;
    int div;
    int mod;
    char cA;
    char cB;
    const double epsilon = self->epsilon;
    const double frameshift_minus_two_score = self->frameshift_minus_two_score;
    const double frameshift_minus_one_score = self->frameshift_minus_one_score;
    const double frameshift_plus_one_score = self->frameshift_plus_one_score;
    const double frameshift_plus_two_score = self->frameshift_plus_two_score;
    Trace** M;
    double score;
    double temp;
    unsigned char trace;
    double* row = NULL;
    PathGenerator* paths;

    static char *kwlist[] = {"sA", "sB0", "sB1", "sB2", NULL};

    if(!PyArg_ParseTupleAndKeywords(args, keywords, "y*y*y*y*",
                                    kwlist, &bA, &bB0, &bB1, &bB2))
        return NULL;

    nA = bA.len;
    nB0 = bB0.len;
    nB1 = bB1.len;
    nB2 = bB2.len;
    if (nB1 == nB0 && nB2 == nB0) nB = 3 * nB0 + 2;
    else if (nB1 == nB0 && nB2 == nB0 - 1) nB = 3 * nB0 + 1;
    else if (nB1 == nB0 - 1 && nB2 == nB0 - 1) nB = 3 * nB0;
    else {
        PyErr_Format(PyExc_RuntimeError,
                     "unexpected length of buffers (%zd, %zd, %zd)",
                     nB0, nB1, nB2);
        PyBuffer_Release(&bA);
        PyBuffer_Release(&bB0);
        PyBuffer_Release(&bB1);
        PyBuffer_Release(&bB2);
        return NULL;
    }

    sA = bA.buf;
    sB[0] = bB0.buf;
    sB[1] = bB1.buf;
    sB[2] = bB2.buf;

    paths = (PathGenerator*)PyType_GenericAlloc(&PathGenerator_Type, 0);
    if (!paths) goto exit;

    paths->nA = (int)nA;
    paths->nB = (int)nB;
    paths->M = NULL;
    paths->length = 0;

    M = PyMem_Malloc((nA+1)*sizeof(Trace*));
    if (!M) goto exit;
    paths->M = M;
    for (i = 0; i <= nA; i++) {
        M[i] = PyMem_Malloc((nB+1)*sizeof(Trace));
        if (!M[i]) {
            Py_DECREF(paths);
            PyErr_NoMemory();
            goto exit;
        }
        M[i][0].trace = 0;
    }
    memset(M[0], '\0', (nB+1)*sizeof(Trace));

    row = PyMem_Malloc((nB+1)*sizeof(double));
    if (!row) goto exit;
    memset(row, '\0', (nB+1)*sizeof(double));

    M = paths->M;
    for (i = 1; i <= nA; i++) {
        cA = sA[i-1];
        for (j = (int)nB; j > 0; j--) {
            score = -DBL_MAX;
            trace = 0;
            if (j >= 3) {
                div = (j - 3) / 3;
                mod = (j - 3) % 3;
                cB = sB[mod][div];
                temp = row[j-1] + (COMPARE_SCORE) + frameshift_minus_two_score;
                if (temp > score + epsilon) {
                    score = temp;
                    trace = FRAMESHIFT_MINUS_TWO;
                }
                else if (temp > score - epsilon) {
                    trace |= FRAMESHIFT_MINUS_TWO;
                }
                temp = row[j-2] + (COMPARE_SCORE) + frameshift_minus_one_score;
                if (temp > score + epsilon) {
                    score = temp;
                    trace = FRAMESHIFT_MINUS_ONE;
                }
                else if (temp > score - epsilon) {
                    trace |= FRAMESHIFT_MINUS_ONE;
                }
                temp = row[j-3] + (COMPARE_SCORE);
                if (temp > score + epsilon) {
                    score = temp;
                    trace = FRAMESHIFT_NONE;
                }
                else if (temp > score - epsilon) {
                    trace |= FRAMESHIFT_NONE;
                }
            }
            if (j >= 4) {
                temp = row[j-4] + (COMPARE_SCORE) + frameshift_plus_one_score;
                if (temp > score + epsilon) {
                    score = temp;
                    trace = FRAMESHIFT_PLUS_ONE;
                }
                else if (temp > score - epsilon) {
                    trace |= FRAMESHIFT_PLUS_ONE;
                }
            }
            if (j >= 5) {
                temp = row[j-5] + (COMPARE_SCORE) + frameshift_plus_two_score;
                if (temp > score + epsilon) {
                    score = temp;
                    trace = FRAMESHIFT_PLUS_TWO;
                }
                else if (temp > score - epsilon) {
                    trace |= FRAMESHIFT_PLUS_TWO;
                }
            }
            M[i][j].trace = trace;
            row[j] = score;
        }
    }
    score = -DBL_MAX;
    for (j = 0; j <= nB; j++) {
        temp = row[j];
        if (temp > score) score = temp;
    }
    for (j = 0; j <= nB; j++) {
        temp = row[j];
        if (temp < score - epsilon) M[nA][j].trace = 0;
        else M[nA][j].path = 0;
    }
    result = Py_BuildValue("fN", score, paths);

exit:
    PyBuffer_Release(&bA);
    PyBuffer_Release(&bB0);
    PyBuffer_Release(&bB1);
    PyBuffer_Release(&bB2);
    PyMem_Free(row);
    if (result == NULL) {
        Py_XDECREF(paths);
        return PyErr_NoMemory();
    }
    return result;
}

static char Aligner_doc[] =
"The CodonAligner class implements a dynamic programming algorithm to\n"
"align a nucleotide sequence to an amino acid sequence.\n";

static PyMethodDef Aligner_methods[] = {
    {"score",
     (PyCFunction)Aligner_score,
     METH_VARARGS | METH_KEYWORDS,
     Aligner_score__doc__
    },
    {"align",
     (PyCFunction)Aligner_align,
     METH_VARARGS | METH_KEYWORDS,
     Aligner_align__doc__
    },
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

static PyTypeObject AlignerType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_codonaligner.CodonAligner",  /* tp_name */
    sizeof(Aligner),               /* tp_basicsize */
    0,                             /* tp_itemsize */
    (destructor)Aligner_dealloc,   /* tp_dealloc */
    0,                             /* tp_print */
    0,                             /* tp_getattr */
    0,                             /* tp_setattr */
    0,                             /* tp_compare */
    (reprfunc)Aligner_repr,        /* tp_repr */
    0,                             /* tp_as_number */
    0,                             /* tp_as_sequence */
    0,                             /* tp_as_mapping */
    0,                             /* tp_hash */
    0,                             /* tp_call */
    (reprfunc)Aligner_str,         /* tp_str */
    0,                             /* tp_getattro */
    0,                             /* tp_setattro */
    0,                             /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    Aligner_doc,                   /* tp_doc */
    0,                             /* tp_traverse */
    0,                             /* tp_clear */
    0,                             /* tp_richcompare */
    0,                             /* tp_weaklistoffset */
    0,                             /* tp_iter */
    0,                             /* tp_iternext */
    Aligner_methods,               /* tp_methods */
    0,                             /* tp_members */
    Aligner_getset,                /* tp_getset */
    0,                             /* tp_base */
    0,                             /* tp_dict */
    0,                             /* tp_descr_get */
    0,                             /* tp_descr_set */
    0,                             /* tp_dictoffset */
    (initproc)Aligner_init,        /* tp_init */
};


/* Module definition */

static char _codonaligner__doc__[] =
"C extension module implementing a dynamic programming algorithm to align a nucleotide sequence to an amino acid sequence";

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_codonaligner",
        _codonaligner__doc__,
        -1,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL
};

PyObject *
PyInit__codonaligner(void)
{
    PyObject* module;
    AlignerType.tp_new = PyType_GenericNew;

    if (PyType_Ready(&AlignerType) < 0 || PyType_Ready(&PathGenerator_Type) < 0)
        return NULL;

    module = PyModule_Create(&moduledef);
    if (!module) return NULL;

    Py_INCREF(&AlignerType);
    /* Reference to AlignerType will be stolen by PyModule_AddObject
     * only if it is successful. */
    if (PyModule_AddObject(module,
                           "CodonAligner", (PyObject*) &AlignerType) < 0) {
        Py_DECREF(&AlignerType);
        Py_DECREF(module);
        return NULL;
    }

    return module;
}
