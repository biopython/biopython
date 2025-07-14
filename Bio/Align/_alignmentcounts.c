/* Copyright 2018-2025 by Michiel de Hoon.  All rights reserved.
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */



#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <float.h>
#include <stdbool.h>
#include <inttypes.h>
#include "_pairwisealigner.h"
#include "substitution_matrices/_arraycore.h"


static PyTypeObject* Aligner_Type = NULL;
static PyTypeObject* Array_Type = NULL;
/* these will be set when initializing the module */

static PyTypeObject AlignmentCounts_Type;
/* defined in this module */


typedef struct {
    PyObject_HEAD
    Py_ssize_t open_left_insertions;
    Py_ssize_t extend_left_insertions;
    Py_ssize_t open_left_deletions;
    Py_ssize_t extend_left_deletions;
    Py_ssize_t open_internal_insertions;
    Py_ssize_t extend_internal_insertions;
    Py_ssize_t open_internal_deletions;
    Py_ssize_t extend_internal_deletions;
    Py_ssize_t open_right_insertions;
    Py_ssize_t extend_right_insertions;
    Py_ssize_t open_right_deletions;
    Py_ssize_t extend_right_deletions;
    Py_ssize_t aligned;
    Py_ssize_t identities;
    Py_ssize_t mismatches;
    Py_ssize_t positives;
    double gap_score;  // Py_NAN
    double substitution_score;  // Py_NAN
} AlignmentCounts;

static PyObject*
AlignmentCounts_str(AlignmentCounts* self)
{
    const double gap_score = self->gap_score;
    const double substitution_score = self->substitution_score;
    const double score = gap_score + substitution_score;
    const Py_ssize_t open_left_insertions = self->open_left_insertions;
    const Py_ssize_t extend_left_insertions = self->extend_left_insertions;
    const Py_ssize_t open_left_deletions = self->open_left_deletions;
    const Py_ssize_t extend_left_deletions = self->extend_left_deletions;
    const Py_ssize_t open_internal_insertions = self->open_internal_insertions;
    const Py_ssize_t extend_internal_insertions = self->extend_internal_insertions;
    const Py_ssize_t open_internal_deletions = self->open_internal_deletions;
    const Py_ssize_t extend_internal_deletions = self->extend_internal_deletions;
    const Py_ssize_t open_right_insertions = self->open_right_insertions;
    const Py_ssize_t extend_right_insertions = self->extend_right_insertions;
    const Py_ssize_t open_right_deletions = self->open_right_deletions;
    const Py_ssize_t extend_right_deletions = self->extend_right_deletions;
    const Py_ssize_t aligned = self->aligned;
    const Py_ssize_t identities = self->identities;
    const Py_ssize_t mismatches = self->mismatches;
    const Py_ssize_t positives = self->positives;
    const Py_ssize_t left_insertions = open_left_insertions
                                     + extend_left_insertions;
    const Py_ssize_t internal_insertions = open_internal_insertions
                                         + extend_internal_insertions;
    const Py_ssize_t right_insertions = open_right_insertions
                                      + extend_right_insertions;
    const Py_ssize_t left_deletions = open_left_deletions
                                    + extend_left_deletions;
    const Py_ssize_t internal_deletions = open_internal_deletions
                                        + extend_internal_deletions;
    const Py_ssize_t right_deletions = open_right_deletions
                                     + extend_right_deletions;
    const Py_ssize_t left_gaps = left_insertions + left_deletions;
    const Py_ssize_t internal_gaps = internal_insertions + internal_deletions;
    const Py_ssize_t right_gaps = right_insertions + right_deletions;
    const Py_ssize_t gaps = left_gaps + internal_gaps + right_gaps;
    char text[2048];
    char* p = text;
    p += sprintf(p, "AlignmentCounts object with\n");
    /* using arguments to PyOS_double_to_string as in
     * float_repr in the Python C source code.
     */
    if (!isnan(score)) {
        char* s = PyOS_double_to_string(score, 'r', 0, Py_DTSF_ADD_DOT_0, NULL);
        if (!s) return NULL;
        p += sprintf(p, "    score = %s:\n", s);
        PyMem_Free(s);
        if (!isnan(substitution_score)) {
            char* s = PyOS_double_to_string(substitution_score, 'r', 0, Py_DTSF_ADD_DOT_0, NULL);
            if (!s) return NULL;
            p += sprintf(p, "        substitution_score = %s,\n", s);
            PyMem_Free(s);
        }
        if (!isnan(gap_score)) {
            char* s = PyOS_double_to_string(gap_score, 'r', 0, Py_DTSF_ADD_DOT_0, NULL);
            if (!s) return NULL;
            p += sprintf(p, "        gap_score = %s.\n", s);
            PyMem_Free(s);
        }
    }
    else {
        if (!isnan(substitution_score)) {
            char* s = PyOS_double_to_string(substitution_score, 'r', 0, Py_DTSF_ADD_DOT_0, NULL);
            if (!s) return NULL;
            p += sprintf(p, "    substitution_score = %s,\n", s);
            PyMem_Free(s);
        }
        if (!isnan(gap_score)) {
            char* s = PyOS_double_to_string(gap_score, 'r', 0, Py_DTSF_ADD_DOT_0, NULL);
            if (!s) return NULL;
            p += sprintf(p, "    gap_score = %s.\n", s);
            PyMem_Free(s);
        }
    }
    p += sprintf(p, "    aligned = %zd:\n", aligned);
    p += sprintf(p, "        identities = %zd,\n", identities);
    if (positives != -1)
        p += sprintf(p, "        positives = %zd,\n", positives);
    p += sprintf(p, "        mismatches = %zd.\n", mismatches);
    p += sprintf(p, "    gaps = %zd:\n", gaps);
    p += sprintf(p, "        left_gaps = %zd:\n", left_gaps);
    p += sprintf(p, "            left_insertions = %zd:\n", left_insertions);
    p += sprintf(p, "                open_left_insertions = %zd,\n", open_left_insertions);
    p += sprintf(p, "                extend_left_insertions = %zd;\n", extend_left_insertions);
    p += sprintf(p, "            left_deletions = %zd:\n", left_deletions);
    p += sprintf(p, "                open_left_deletions = %zd,\n", open_left_deletions);
    p += sprintf(p, "                extend_left_deletions = %zd;\n", extend_left_deletions);
    p += sprintf(p, "        internal_gaps = %zd:\n", internal_gaps);
    p += sprintf(p, "            internal_insertions = %zd:\n", internal_insertions);
    p += sprintf(p, "                open_internal_insertions = %zd,\n", open_internal_insertions);
    p += sprintf(p, "                extend_internal_insertions = %zd;\n", extend_internal_insertions);
    p += sprintf(p, "            internal_deletions = %zd:\n", internal_deletions);
    p += sprintf(p, "                open_internal_deletions = %zd,\n", open_internal_deletions);
    p += sprintf(p, "                extend_internal_deletions = %zd;\n", extend_internal_deletions);
    p += sprintf(p, "        right_gaps = %zd:\n", right_gaps);
    p += sprintf(p, "            right_insertions = %zd:\n", right_insertions);
    p += sprintf(p, "                open_right_insertions = %zd,\n", open_right_insertions);
    p += sprintf(p, "                extend_right_insertions = %zd;\n", extend_right_insertions);
    p += sprintf(p, "            right_deletions = %zd:\n", right_deletions);
    p += sprintf(p, "                open_right_deletions = %zd,\n", open_right_deletions);
    p += sprintf(p, "                extend_right_deletions = %zd.\n", extend_right_deletions);
    return PyUnicode_FromString(text);
}

static PyObject*
AlignmentCounts_repr(AlignmentCounts* self)
{
    Py_uintptr_t id = 0;
    const double substitution_score = self->substitution_score;
    const double gap_score = self->gap_score;
    const double score = gap_score + substitution_score;
    const Py_ssize_t aligned = self->aligned;
    const Py_ssize_t identities = self->identities;
    const Py_ssize_t mismatches = self->mismatches;
    const Py_ssize_t positives = self->positives;
    const Py_ssize_t gaps = self->open_left_insertions
                          + self->extend_left_insertions
                          + self->open_internal_insertions
                          + self->extend_internal_insertions
                          + self->open_right_insertions
                          + self->extend_right_insertions
                          + self->open_left_deletions
                          + self->extend_left_deletions
                          + self->open_internal_deletions
                          + self->extend_internal_deletions
                          + self->open_right_deletions
                          + self->extend_right_deletions;
    char text[1024];
    char* p = text;
    p += sprintf(p, "<AlignmentCounts object (");
    /* using arguments to PyOS_double_to_string as in
     * float_repr in the Python C source code.
     */
    if (!isnan(score)) {
        char* s = PyOS_double_to_string(score, 'r', 0, Py_DTSF_ADD_DOT_0, NULL);
        if (!s) return NULL;
        p += sprintf(p, "score = %s; ", s);
        PyMem_Free(s);
    }
    if (!isnan(substitution_score)) {
        char* s = PyOS_double_to_string(substitution_score, 'r', 0, Py_DTSF_ADD_DOT_0, NULL);
        if (!s) return NULL;
        p += sprintf(p, "substitution score = %s; ", s);
        PyMem_Free(s);
    }
    if (!isnan(gap_score)) {
        char* s = PyOS_double_to_string(gap_score, 'r', 0, Py_DTSF_ADD_DOT_0, NULL);
        if (!s) return NULL;
        p += sprintf(p, "gap score = %s; ", s);
        PyMem_Free(s);
    }
    p += sprintf(p, "%zd aligned letters; ", aligned);
    p += sprintf(p, "%zd identities; ", identities);
    p += sprintf(p, "%zd mismatches; ", mismatches);
    if (positives != -1)
        p += sprintf(p, "%zd positives; ", positives);
#ifdef PYPY_VERSION
    // For PyPy, use PyObject_CallFunction to get id(self)
    PyObject* builtins = PyEval_GetBuiltins();
    PyObject* id_func = PyDict_GetItemString(builtins, "id");
    PyObject* id_result = PyObject_CallFunctionObjArgs(id_func, self, NULL);
    if (id_result) {
        if (PyLong_Check(id_result)) {
            id = (Py_uintptr_t)PyLong_AsUnsignedLongLong(id_result);
        }
        Py_DECREF(id_result);
    }
#else
    // In CPython, id(self) is just the address
    id = (Py_uintptr_t)self;
#endif
    sprintf(p, "%zd gaps) at 0x%" PRIxPTR ">", gaps, id);

    return PyUnicode_FromString(text);
}

static char AlignmentCounts_doc[] =
"AlignmentCounts objects store the number of aligned characters, identities,\n"
"mismatches, gaps, and score of an alignment.\n";

static char AlignmentCounts_aligned__doc__[] = "number of aligned characters in the alignment";

static PyObject*
AlignmentCounts_get_aligned(AlignmentCounts* self, void* closure)
{
    return PyLong_FromSsize_t(self->aligned);
}

static char AlignmentCounts_identities__doc__[] = "number of matched letters in the alignment";

static PyObject*
AlignmentCounts_get_identities(AlignmentCounts* self, void* closure)
{
    return PyLong_FromSsize_t(self->identities);
}

static char AlignmentCounts_mismatches__doc__[] = "number of mismatched letters in the alignment";

static PyObject*
AlignmentCounts_get_mismatches(AlignmentCounts* self, void* closure)
{
    return PyLong_FromSsize_t(self->mismatches);
}

static char AlignmentCounts_positives__doc__[] = "number of aligned letters with a positive substitution score";

static PyObject*
AlignmentCounts_get_positives(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t positives = self->positives;
    if (positives == -1) {
        Py_INCREF(Py_None);
        return Py_None;
    }
    else return PyLong_FromSsize_t(positives);
}

static char AlignmentCounts_score__doc__[] = "alignment score, or None if unknown";

static PyObject*
AlignmentCounts_get_score(AlignmentCounts* self, void* closure)
{
    const double score = self->gap_score + self->substitution_score;
    if (isnan(score)) {
        Py_INCREF(Py_None);
        return Py_None;
    }
    return PyFloat_FromDouble(score);
}

static char AlignmentCounts_gap_score__doc__[] = "total gap score, or None if unknown";

static PyObject*
AlignmentCounts_get_gap_score(AlignmentCounts* self, void* closure)
{
    const double gap_score = self->gap_score;
    if (isnan(gap_score)) {
        Py_INCREF(Py_None);
        return Py_None;
    }
    return PyFloat_FromDouble(gap_score);
}

static char AlignmentCounts_substitution_score__doc__[] = "total substitution score of letters aligned to each other, or None if unknown";

static PyObject*
AlignmentCounts_get_substitution_score(AlignmentCounts* self, void* closure)
{
    const double substitution_score = self->substitution_score;
    if (isnan(substitution_score)) {
        Py_INCREF(Py_None);
        return Py_None;
    }
    return PyFloat_FromDouble(substitution_score);
}

static char AlignmentCounts_open_left_insertions__doc__[] = "number of insertion gaps opened on the left side of the alignment";

static PyObject*
AlignmentCounts_get_open_left_insertions(AlignmentCounts* self, void* closure)
{
    return PyLong_FromSsize_t(self->open_left_insertions);
}

static char AlignmentCounts_extend_left_insertions__doc__[] = "number of insertion gap extensions on the left side of the alignment";

static PyObject*
AlignmentCounts_get_extend_left_insertions(AlignmentCounts* self, void* closure)
{
    return PyLong_FromSsize_t(self->extend_left_insertions);
}

static char AlignmentCounts_left_insertions__doc__[] = "number of letters inserted on the left side of the alignment";

static PyObject*
AlignmentCounts_get_left_insertions(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_left_insertions = self->open_left_insertions;
    const Py_ssize_t extend_left_insertions = self->extend_left_insertions;
    return PyLong_FromSsize_t(open_left_insertions + extend_left_insertions);
}

static char AlignmentCounts_open_left_deletions__doc__[] = "the number of deletion gaps opened on the left side of the alignment";

static PyObject*
AlignmentCounts_get_open_left_deletions(AlignmentCounts* self, void* closure)
{
    return PyLong_FromSsize_t(self->open_left_deletions);
}

static char AlignmentCounts_extend_left_deletions__doc__[] = "number of deletion gap extensions on the left side of the alignment";

static PyObject*
AlignmentCounts_get_extend_left_deletions(AlignmentCounts* self, void* closure)
{
    return PyLong_FromSsize_t(self->extend_left_deletions);
}

static char AlignmentCounts_left_deletions__doc__[] = "number of characters deleted on the left side of the alignment";

static PyObject*
AlignmentCounts_get_left_deletions(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_left_deletions = self->open_left_deletions;
    const Py_ssize_t extend_left_deletions = self->extend_left_deletions;
    return PyLong_FromSsize_t(open_left_deletions + extend_left_deletions);
}

static char AlignmentCounts_open_left_gaps__doc__[] = "number of gaps opened on the left side of the alignment";

static PyObject*
AlignmentCounts_get_open_left_gaps(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_left_insertions = self->open_left_insertions;
    const Py_ssize_t open_left_deletions = self->open_left_deletions;
    return PyLong_FromSsize_t(open_left_insertions + open_left_deletions);
}

static char AlignmentCounts_extend_left_gaps__doc__[] = "number of gap extensions on the left side of the alignment";

static PyObject*
AlignmentCounts_get_extend_left_gaps(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t extend_left_insertions = self->extend_left_insertions;
    const Py_ssize_t extend_left_deletions = self->extend_left_deletions;
    return PyLong_FromSsize_t(extend_left_insertions + extend_left_deletions);
}

static char AlignmentCounts_left_gaps__doc__[] = "total gap length on the left side of the alignment";

static PyObject*
AlignmentCounts_get_left_gaps(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_left_insertions = self->open_left_insertions;
    const Py_ssize_t extend_left_insertions = self->extend_left_insertions;
    const Py_ssize_t open_left_deletions = self->open_left_deletions;
    const Py_ssize_t extend_left_deletions = self->extend_left_deletions;
    return PyLong_FromSsize_t(open_left_insertions
                            + extend_left_insertions
                            + open_left_deletions
                            + extend_left_deletions);
}

static char AlignmentCounts_open_internal_insertions__doc__[] = "number of insertion gaps opened in the interior of the alignment";

static PyObject*
AlignmentCounts_get_open_internal_insertions(AlignmentCounts* self, void* closure)
{
    return PyLong_FromSsize_t(self->open_internal_insertions);
}

static char AlignmentCounts_extend_internal_insertions__doc__[] = "number of insertion gas extensions in the interior of the alignment";

static PyObject*
AlignmentCounts_get_extend_internal_insertions(AlignmentCounts* self, void* closure)
{
    return PyLong_FromSsize_t(self->extend_internal_insertions);
}

static char AlignmentCounts_internal_insertions__doc__[] = "number of letters inserted in the interior of the alignment";

static PyObject*
AlignmentCounts_get_internal_insertions(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_internal_insertions = self->open_internal_insertions;
    const Py_ssize_t extend_internal_insertions = self->extend_internal_insertions;
    return PyLong_FromSsize_t(open_internal_insertions + extend_internal_insertions);
}

static char AlignmentCounts_open_internal_deletions__doc__[] = "number of deletion gaps opened in the interior of the alignment";

static PyObject*
AlignmentCounts_get_open_internal_deletions(AlignmentCounts* self, void* closure)
{
    return PyLong_FromSsize_t(self->open_internal_deletions);
}

static char AlignmentCounts_extend_internal_deletions__doc__[] = "number of deletion gap exensions  in the interior of the alignment";

static PyObject*
AlignmentCounts_get_extend_internal_deletions(AlignmentCounts* self, void* closure)
{
    return PyLong_FromSsize_t(self->extend_internal_deletions);
}

static char AlignmentCounts_internal_deletions__doc__[] = "number of characters deleted from the alignment";

static PyObject*
AlignmentCounts_get_internal_deletions(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_internal_deletions = self->open_internal_deletions;
    const Py_ssize_t extend_internal_deletions = self->extend_internal_deletions;
    return PyLong_FromSsize_t(open_internal_deletions + extend_internal_deletions);
}

static PyObject*
AlignmentCounts_get_open_internal_gaps(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_internal_insertions = self->open_internal_insertions;
    const Py_ssize_t open_internal_deletions = self->open_internal_deletions;
    return PyLong_FromSsize_t(open_internal_insertions
                            + open_internal_deletions);
}

static char AlignmentCounts_open_internal_gaps__doc__[] = "number of gaps opened in the interior of the alignment";

static PyObject*
AlignmentCounts_get_extend_internal_gaps(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t extend_internal_insertions = self->extend_internal_insertions;
    const Py_ssize_t extend_internal_deletions = self->extend_internal_deletions;
    return PyLong_FromSsize_t(extend_internal_insertions
                            + extend_internal_deletions);
}

static char AlignmentCounts_extend_internal_gaps__doc__[] = "number of gap extensions in the interior of the alignment";

static PyObject*
AlignmentCounts_get_internal_gaps(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_internal_insertions = self->open_internal_insertions;
    const Py_ssize_t extend_internal_insertions = self->extend_internal_insertions;
    const Py_ssize_t open_internal_deletions = self->open_internal_deletions;
    const Py_ssize_t extend_internal_deletions = self->extend_internal_deletions;
    return PyLong_FromSsize_t(open_internal_insertions
                            + extend_internal_insertions
                            + open_internal_deletions
                            + extend_internal_deletions);
}

static char AlignmentCounts_internal_gaps__doc__[] = "total length of gaps within the alignment";

static char AlignmentCounts_open_right_insertions__doc__[] = "the number of insertion gaps opened on the right side of the alignment";

static PyObject*
AlignmentCounts_get_open_right_insertions(AlignmentCounts* self, void* closure)
{
    return PyLong_FromSsize_t(self->open_right_insertions);
}

static char AlignmentCounts_extend_right_insertions__doc__[] = "the number of insertion gap extensions on the right side of the alignment";

static PyObject*
AlignmentCounts_get_extend_right_insertions(AlignmentCounts* self, void* closure)
{
    return PyLong_FromSsize_t(self->extend_right_insertions);
}

static char AlignmentCounts_right_insertions__doc__[] = "number of letters inserted on the right side of the alignment";

static PyObject*
AlignmentCounts_get_right_insertions(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_right_insertions = self->open_right_insertions;
    const Py_ssize_t extend_right_insertions = self->extend_right_insertions;
    return PyLong_FromSsize_t(open_right_insertions + extend_right_insertions);
}

static char AlignmentCounts_open_right_deletions__doc__[] = "number of deletion gaps opened on the right side of the alignment";

static PyObject*
AlignmentCounts_get_open_right_deletions(AlignmentCounts* self, void* closure)
{
    return PyLong_FromSsize_t(self->open_right_deletions);
}

static char AlignmentCounts_extend_right_deletions__doc__[] = "number of deletion gap extensions on the right side of the alignment";

static PyObject*
AlignmentCounts_get_extend_right_deletions(AlignmentCounts* self, void* closure)
{
    return PyLong_FromSsize_t(self->extend_right_deletions);
}

static char AlignmentCounts_right_deletions__doc__[] = "number of letters deleted on the right side of the alignment";

static PyObject*
AlignmentCounts_get_right_deletions(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_right_deletions = self->open_right_deletions;
    const Py_ssize_t extend_right_deletions = self->extend_right_deletions;
    return PyLong_FromSsize_t(open_right_deletions + extend_right_deletions);
}

static char AlignmentCounts_open_right_gaps__doc__[] = "number of gaps opened on the right side of the alignment";

static PyObject*
AlignmentCounts_get_open_right_gaps(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_right_insertions = self->open_right_insertions;
    const Py_ssize_t open_right_deletions = self->open_right_deletions;
    return PyLong_FromSsize_t(open_right_insertions + open_right_deletions);
}

static char AlignmentCounts_extend_right_gaps__doc__[] = "number of gap extensions on the right side of the alignment";

static PyObject*
AlignmentCounts_get_extend_right_gaps(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t extend_right_insertions = self->extend_right_insertions;
    const Py_ssize_t extend_right_deletions = self->extend_right_deletions;
    return PyLong_FromSsize_t(extend_right_insertions + extend_right_deletions);
}

static char AlignmentCounts_right_gaps__doc__[] = "total gap length on the right side of the alignment";

static PyObject*
AlignmentCounts_get_right_gaps(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_right_insertions = self->open_right_insertions;
    const Py_ssize_t extend_right_insertions = self->extend_right_insertions;
    const Py_ssize_t open_right_deletions = self->open_right_deletions;
    const Py_ssize_t extend_right_deletions = self->extend_right_deletions;
    return PyLong_FromSsize_t(open_right_insertions
                            + extend_right_insertions
                            + open_right_deletions
                            + extend_right_deletions);
}

static char AlignmentCounts_open_gaps__doc__[] = "number of geps opened in the alignment";

static PyObject*
AlignmentCounts_get_open_gaps(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_left_insertions = self->open_left_insertions;
    const Py_ssize_t open_left_deletions = self->open_left_deletions;
    const Py_ssize_t open_internal_insertions = self->open_internal_insertions;
    const Py_ssize_t open_internal_deletions = self->open_internal_deletions;
    const Py_ssize_t open_right_insertions = self->open_right_insertions;
    const Py_ssize_t open_right_deletions = self->open_right_deletions;
    const Py_ssize_t open_gaps = open_left_insertions
                               + open_left_deletions
                               + open_internal_insertions
                               + open_internal_deletions
                               + open_right_insertions
                               + open_right_deletions;
    return PyLong_FromSsize_t(open_gaps);
}

static char AlignmentCounts_extend_gaps__doc__[] = "number of gep extensions in the alignment";

static PyObject*
AlignmentCounts_get_extend_gaps(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t extend_left_insertions = self->extend_left_insertions;
    const Py_ssize_t extend_left_deletions = self->extend_left_deletions;
    const Py_ssize_t extend_internal_insertions = self->extend_internal_insertions;
    const Py_ssize_t extend_internal_deletions = self->extend_internal_deletions;
    const Py_ssize_t extend_right_insertions = self->extend_right_insertions;
    const Py_ssize_t extend_right_deletions = self->extend_right_deletions;
    const Py_ssize_t extend_gaps = extend_left_insertions
                                 + extend_left_deletions
                                 + extend_internal_insertions
                                 + extend_internal_deletions
                                 + extend_right_insertions
                                 + extend_right_deletions;
    return PyLong_FromSsize_t(extend_gaps);
}

static char AlignmentCounts_insertions__doc__[] = "total number of letters inserted";

static PyObject*
AlignmentCounts_get_insertions(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_left_insertions = self->open_left_insertions;
    const Py_ssize_t extend_left_insertions = self->extend_left_insertions;
    const Py_ssize_t open_internal_insertions = self->open_internal_insertions;
    const Py_ssize_t extend_internal_insertions = self->extend_internal_insertions;
    const Py_ssize_t open_right_insertions = self->open_right_insertions;
    const Py_ssize_t extend_right_insertions = self->extend_right_insertions;
    return PyLong_FromSsize_t(open_left_insertions
                            + extend_left_insertions
                            + open_internal_insertions
                            + extend_internal_insertions
                            + open_right_insertions
                            + extend_right_insertions);
}

static char AlignmentCounts_open_insertions__doc__[] = "number of insertion geps opened in the alignment";

static PyObject*
AlignmentCounts_get_open_insertions(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_left_insertions = self->open_left_insertions;
    const Py_ssize_t open_internal_insertions = self->open_internal_insertions;
    const Py_ssize_t open_right_insertions = self->open_right_insertions;
    const Py_ssize_t open_insertions = open_left_insertions
                                     + open_internal_insertions
                                     + open_right_insertions;
    return PyLong_FromSsize_t(open_insertions);
}

static char AlignmentCounts_extend_insertions__doc__[] = "number of insertion gap extensions in the alignment";

static PyObject*
AlignmentCounts_get_extend_insertions(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t extend_left_insertions = self->extend_left_insertions;
    const Py_ssize_t extend_internal_insertions = self->extend_internal_insertions;
    const Py_ssize_t extend_right_insertions = self->extend_right_insertions;
    const Py_ssize_t extend_insertions = extend_left_insertions
                                       + extend_internal_insertions
                                       + extend_right_insertions;
    return PyLong_FromSsize_t(extend_insertions);
}

static char AlignmentCounts_deletions__doc__[] = "total number of letters deleted";

static PyObject*
AlignmentCounts_get_deletions(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_left_deletions = self->open_left_deletions;
    const Py_ssize_t extend_left_deletions = self->extend_left_deletions;
    const Py_ssize_t open_internal_deletions = self->open_internal_deletions;
    const Py_ssize_t extend_internal_deletions = self->extend_internal_deletions;
    const Py_ssize_t open_right_deletions = self->open_right_deletions;
    const Py_ssize_t extend_right_deletions = self->extend_right_deletions;
    return PyLong_FromSsize_t(open_left_deletions
                            + extend_left_deletions
                            + open_internal_deletions
                            + extend_internal_deletions
                            + open_right_deletions
                            + extend_right_deletions);
}

static char AlignmentCounts_open_deletions__doc__[] = "number of deletion geps opened in the alignment";

static PyObject*
AlignmentCounts_get_open_deletions(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_left_deletions = self->open_left_deletions;
    const Py_ssize_t open_internal_deletions = self->open_internal_deletions;
    const Py_ssize_t open_right_deletions = self->open_right_deletions;
    const Py_ssize_t open_gaps = open_left_deletions
                               + open_internal_deletions
                               + open_right_deletions;
    return PyLong_FromSsize_t(open_gaps);
}

static char AlignmentCounts_extend_deletions__doc__[] = "number of deletion gep extensions in the alignment";

static PyObject*
AlignmentCounts_get_extend_deletions(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t extend_left_deletions = self->extend_left_deletions;
    const Py_ssize_t extend_internal_deletions = self->extend_internal_deletions;
    const Py_ssize_t extend_right_deletions = self->extend_right_deletions;
    const Py_ssize_t extend_gaps = extend_left_deletions
                                 + extend_internal_deletions
                                 + extend_right_deletions;
    return PyLong_FromSsize_t(extend_gaps);
}

static char AlignmentCounts_gaps__doc__[] = "total gap length";

static PyObject*
AlignmentCounts_get_gaps(AlignmentCounts* self, void* closure)
{
    const Py_ssize_t open_left_insertions = self->open_left_insertions;
    const Py_ssize_t extend_left_insertions = self->extend_left_insertions;
    const Py_ssize_t open_left_deletions = self->open_left_deletions;
    const Py_ssize_t extend_left_deletions = self->extend_left_deletions;
    const Py_ssize_t open_internal_insertions = self->open_internal_insertions;
    const Py_ssize_t extend_internal_insertions = self->extend_internal_insertions;
    const Py_ssize_t open_internal_deletions = self->open_internal_deletions;
    const Py_ssize_t extend_internal_deletions = self->extend_internal_deletions;
    const Py_ssize_t open_right_insertions = self->open_right_insertions;
    const Py_ssize_t extend_right_insertions = self->extend_right_insertions;
    const Py_ssize_t open_right_deletions = self->open_right_deletions;
    const Py_ssize_t extend_right_deletions = self->extend_right_deletions;
    const Py_ssize_t gaps = open_left_insertions
                          + extend_left_insertions
                          + open_left_deletions
                          + extend_left_deletions
                          + open_internal_insertions
                          + extend_internal_insertions
                          + open_internal_deletions
                          + extend_internal_deletions
                          + open_right_insertions
                          + extend_right_insertions
                          + open_right_deletions
                          + extend_right_deletions;
    return PyLong_FromSsize_t(gaps);
}

static PyGetSetDef AlignmentCounts_getset[] = {
    {"score",
        (getter)AlignmentCounts_get_score, NULL,
        AlignmentCounts_score__doc__, NULL},
    {"aligned",
        (getter)AlignmentCounts_get_aligned, NULL,
        AlignmentCounts_aligned__doc__, NULL},
    {"substitution_score",
        (getter)AlignmentCounts_get_substitution_score, NULL,
        AlignmentCounts_substitution_score__doc__, NULL},
    {"identities",
        (getter)AlignmentCounts_get_identities, NULL,
        AlignmentCounts_identities__doc__, NULL},
    {"mismatches",
        (getter)AlignmentCounts_get_mismatches, NULL,
        AlignmentCounts_mismatches__doc__, NULL},
    {"positives",
        (getter)AlignmentCounts_get_positives, NULL,
        AlignmentCounts_positives__doc__, NULL},
    {"gap_score",
        (getter)AlignmentCounts_get_gap_score, NULL,
        AlignmentCounts_gap_score__doc__, NULL},
    {"gaps",
        (getter)AlignmentCounts_get_gaps, NULL,
        AlignmentCounts_gaps__doc__, NULL},
    {"open_gaps",
        (getter)AlignmentCounts_get_open_gaps, NULL,
        AlignmentCounts_open_gaps__doc__, NULL},
    {"extend_gaps",
        (getter)AlignmentCounts_get_extend_gaps, NULL,
        AlignmentCounts_extend_gaps__doc__, NULL},
    {"open_left_gaps",
        (getter)AlignmentCounts_get_open_left_gaps, NULL,
        AlignmentCounts_open_left_gaps__doc__, NULL},
    {"open_right_gaps",
        (getter)AlignmentCounts_get_open_right_gaps, NULL,
        AlignmentCounts_open_right_gaps__doc__, NULL},
    {"open_internal_gaps",
        (getter)AlignmentCounts_get_open_internal_gaps, NULL,
        AlignmentCounts_open_internal_gaps__doc__, NULL},
    {"extend_left_gaps",
        (getter)AlignmentCounts_get_extend_left_gaps, NULL,
        AlignmentCounts_extend_left_gaps__doc__, NULL},
    {"extend_right_gaps",
        (getter)AlignmentCounts_get_extend_right_gaps, NULL,
        AlignmentCounts_extend_right_gaps__doc__, NULL},
    {"extend_internal_gaps",
        (getter)AlignmentCounts_get_extend_internal_gaps, NULL,
        AlignmentCounts_extend_internal_gaps__doc__, NULL},
    {"open_left_insertions",
        (getter)AlignmentCounts_get_open_left_insertions, NULL,
        AlignmentCounts_open_left_insertions__doc__, NULL},
    {"open_left_deletions",
        (getter)AlignmentCounts_get_open_left_deletions, NULL,
        AlignmentCounts_open_left_deletions__doc__, NULL},
    {"open_right_insertions",
        (getter)AlignmentCounts_get_open_right_insertions, NULL,
        AlignmentCounts_open_right_insertions__doc__, NULL},
    {"open_right_deletions",
        (getter)AlignmentCounts_get_open_right_deletions, NULL,
        AlignmentCounts_open_right_deletions__doc__, NULL},
    {"open_internal_insertions",
        (getter)AlignmentCounts_get_open_internal_insertions, NULL,
        AlignmentCounts_open_internal_insertions__doc__, NULL},
    {"open_internal_deletions",
        (getter)AlignmentCounts_get_open_internal_deletions, NULL,
        AlignmentCounts_open_internal_deletions__doc__, NULL},
    {"extend_left_insertions",
        (getter)AlignmentCounts_get_extend_left_insertions, NULL,
        AlignmentCounts_extend_left_insertions__doc__, NULL},
    {"extend_left_deletions",
        (getter)AlignmentCounts_get_extend_left_deletions, NULL,
        AlignmentCounts_extend_left_deletions__doc__, NULL},
    {"extend_right_insertions",
        (getter)AlignmentCounts_get_extend_right_insertions, NULL,
        AlignmentCounts_extend_right_insertions__doc__, NULL},
    {"extend_right_deletions",
        (getter)AlignmentCounts_get_extend_right_deletions, NULL,
        AlignmentCounts_extend_right_deletions__doc__, NULL},
    {"extend_internal_insertions",
        (getter)AlignmentCounts_get_extend_internal_insertions, NULL,
        AlignmentCounts_extend_internal_insertions__doc__, NULL},
    {"extend_internal_deletions",
        (getter)AlignmentCounts_get_extend_internal_deletions, NULL,
        AlignmentCounts_extend_internal_deletions__doc__, NULL},
    {"left_insertions",
        (getter)AlignmentCounts_get_left_insertions, NULL,
        AlignmentCounts_left_insertions__doc__, NULL},
    {"left_deletions",
        (getter)AlignmentCounts_get_left_deletions, NULL,
        AlignmentCounts_left_deletions__doc__, NULL},
    {"right_insertions",
        (getter)AlignmentCounts_get_right_insertions, NULL,
        AlignmentCounts_right_insertions__doc__, NULL},
    {"right_deletions",
        (getter)AlignmentCounts_get_right_deletions, NULL,
        AlignmentCounts_right_deletions__doc__, NULL},
    {"internal_insertions",
        (getter)AlignmentCounts_get_internal_insertions, NULL,
        AlignmentCounts_internal_insertions__doc__, NULL},
    {"internal_deletions",
        (getter)AlignmentCounts_get_internal_deletions, NULL,
        AlignmentCounts_internal_deletions__doc__, NULL},
    {"insertions",
        (getter)AlignmentCounts_get_insertions, NULL,
        AlignmentCounts_insertions__doc__, NULL},
    {"open_insertions",
        (getter)AlignmentCounts_get_open_insertions, NULL,
        AlignmentCounts_open_insertions__doc__, NULL},
    {"extend_insertions",
        (getter)AlignmentCounts_get_extend_insertions, NULL,
        AlignmentCounts_extend_insertions__doc__, NULL},
    {"deletions",
        (getter)AlignmentCounts_get_deletions, NULL,
        AlignmentCounts_deletions__doc__, NULL},
    {"open_deletions",
        (getter)AlignmentCounts_get_open_deletions, NULL,
        AlignmentCounts_open_deletions__doc__, NULL},
    {"extend_deletions",
        (getter)AlignmentCounts_get_extend_deletions, NULL,
        AlignmentCounts_extend_deletions__doc__, NULL},
    {"left_gaps",
        (getter)AlignmentCounts_get_left_gaps, NULL,
        AlignmentCounts_left_gaps__doc__, NULL},
    {"right_gaps",
        (getter)AlignmentCounts_get_right_gaps, NULL,
        AlignmentCounts_right_gaps__doc__, NULL},
    {"internal_gaps",
        (getter)AlignmentCounts_get_internal_gaps, NULL,
        AlignmentCounts_internal_gaps__doc__, NULL},
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

static int
sequence_converter(PyObject* argument, void* pointer)
{
    Py_buffer* view = pointer;
    const int flag = PyBUF_FORMAT | PyBUF_C_CONTIGUOUS;

    if (PyObject_GetBuffer(argument, view, flag) == 0) {
        if (view->ndim == 1 && view->len > 0) {
            if ((strcmp(view->format, "i") == 0
              || strcmp(view->format, "l") == 0)
              && view->itemsize == sizeof(int)) {
                /* buffer contains int values */ return 1;
}
            if ((strcmp(view->format, "c") == 0
              || strcmp(view->format, "b") == 0
              || strcmp(view->format, "B") == 0)
              && view->itemsize == sizeof(char))
                /* buffer contains int values */ return 1;
        }
        PyBuffer_Release(view);
    } else PyErr_Clear();
    view->buf = NULL;
    if (PySequence_Check(argument)) {
        view->obj = argument;
        return 1;
    }
    else if (argument == Py_None) {
        view->obj = NULL;
        return 1;
    }
    else return 0;
}
 
static int
coordinates_converter(PyObject* argument, void* pointer)
{
    Py_buffer* view = pointer;
    if (argument == NULL) {
        PyBuffer_Release(view);
        return 1;
    }

    if (PyObject_GetBuffer(argument, view, PyBUF_STRIDED) != 0) return 0;
    if (view->ndim != 2 || view->itemsize != sizeof(Py_ssize_t)) {
        PyErr_SetString(PyExc_ValueError,
            "coordinates must be a 2D array of Py_ssize_t integers");
        PyBuffer_Release(view);
        return 0;
    }
    return Py_CLEANUP_SUPPORTED;
}

static int
strands_converter(PyObject* argument, void* pointer)
{
    Py_buffer* view = pointer;
    if (argument == NULL) {
        PyBuffer_Release(view);
        return 1;
    }

    if (PyObject_GetBuffer(argument, view, PyBUF_CONTIG_RO) != 0) return 0;
    if (view->ndim != 1 || view->itemsize != sizeof(bool)) {
        PyErr_SetString(PyExc_ValueError,
            "strands must be a 1D array of bool variables");
        PyBuffer_Release(view);
        return 0;
    }
    return Py_CLEANUP_SUPPORTED;
}

static int
substitution_matrix_converter(PyObject* argument, void* pointer)
{
    const int flag = PyBUF_FORMAT | PyBUF_ND;
    Py_buffer* view = pointer;
    if (argument == NULL) {
        PyBuffer_Release(view);
        return 1;
    }
    if (PyObject_GetBuffer(argument, view, flag) != 0) {
        PyErr_SetString(PyExc_ValueError, "expected a matrix");
        return 0;
    }
    if (view->ndim != 2) {
        PyErr_Format(PyExc_ValueError,
         "substitution matrix has incorrect rank (%d expected 2)",
          view->ndim);
        PyBuffer_Release(view);
        return 0;
    }
    if (view->len == 0) {
        PyErr_SetString(PyExc_ValueError, "substitution matrix has zero size");
        PyBuffer_Release(view);
        return 0;
    }
    if (strcmp(view->format, "d") != 0) {
        PyErr_SetString(PyExc_ValueError,
                "substitution matrix should contain float values");
        PyBuffer_Release(view);
        return 0;
    }
    if (view->itemsize != sizeof(double)) {
        PyErr_Format(PyExc_RuntimeError,
                    "substitution matrix has unexpected item byte size "
                    "(%zd, expected %zd)", view->itemsize, sizeof(double));
        PyBuffer_Release(view);
        return 0;
    }
    if (view->shape[0] != view->shape[1]) {
        PyErr_Format(PyExc_ValueError,
                    "substitution matrix should be square "
                    "(found a %zd x %zd matrix)",
                    view->shape[0], view->shape[1]);
        PyBuffer_Release(view);
        return 0;
    }
    return Py_CLEANUP_SUPPORTED;
}

static bool inline
check_indices(int c, Py_ssize_t j, Py_ssize_t l, Py_ssize_t m)
{
    if (c < 0) {
        PyErr_Format(PyExc_ValueError,
            "sequences[%zd][%zd] is negative (%d)", j, l, c);
        return false;
    }
    if (c >= m) { \
        PyErr_Format(PyExc_ValueError,
            "sequence[%zd][%zd] is out of bound"
            " (%d, should be < %zd)", j, l, c, m);
        return false;
    }
    return true;
}

static bool inline
map_indices(int* cA, int* cB, int* mapping)
{
    *cA = mapping[*cA];
    *cB = mapping[*cB];
    if (*cA == MISSING_LETTER || *cB == MISSING_LETTER) {
        PyErr_SetString(PyExc_ValueError,
            "sequence contains letters not in the alphabet");
        return false;
    }
    return true;
}

static inline PyObject* get_lazy_data(PyObject* sequence, Py_ssize_t start, Py_ssize_t end, Py_ssize_t j, char** b)
{
    PyObject* obj = PySequence_GetSlice(sequence, start, end);
    if (!obj) return NULL;
    if (PyBytes_Check(obj)) {
        if (PyBytes_GET_SIZE(obj) == end - start) {
            *b = PyBytes_AS_STRING(obj) - start;
            return obj;
        }
        PyErr_Format(PyExc_ValueError,
            "alignment.sequences[%zd][%zd:%zd] did not return a bytes object of size %zd",
            j, start, end, end - start);
        Py_DECREF(obj);
    }
    return obj;
}

static inline void reset_lazy_data(PyObject* obj, char** b) {
    if (obj) {
        Py_DECREF(obj);
        *b = NULL;
    }
}

static void inline
add_identities_mismatches(int cA, int cB, int wildcard,
                          Py_ssize_t* identities, Py_ssize_t* mismatches)
{
    if (cA == wildcard || cB == wildcard) return;
    else if (cA == cB) (*identities)++;
    else (*mismatches)++;
}

static void inline
add_identities_mismatches_score(int cA, int cB, int wildcard,
                                Py_buffer* substitution_matrix,
                                Py_ssize_t* identities,
                                Py_ssize_t* mismatches,
                                Py_ssize_t* positives,
                                double* substitution_score)
{
    const double value = *((double*)substitution_matrix->buf
                           + cA * substitution_matrix->shape[0]
                           + cB);
    if (cA == cB) (*identities)++;
    else (*mismatches)++;
    if (value > 0) (*positives)++;
    *substitution_score += value;
}

static inline bool add_gaps(int* path, const int direction, bool strand,
                            PyObject* score_function,
                            Py_ssize_t start, Py_ssize_t end,
                            Py_ssize_t left, Py_ssize_t right,
                            Py_ssize_t gapsize,
                            Py_ssize_t* open_left_gaps,
                            Py_ssize_t* extend_left_gaps,
                            Py_ssize_t* open_internal_gaps,
                            Py_ssize_t* extend_internal_gaps,
                            Py_ssize_t* open_right_gaps,
                            Py_ssize_t* extend_right_gaps,
                            double* gap_score)
{
    if (*path == direction) {
        if (start == left)
            *extend_left_gaps += gapsize;
        else if (end == right)
            *extend_right_gaps += gapsize;
        else
            *extend_internal_gaps += gapsize;
    }
    else {
        if (start == left) {
            (*open_left_gaps)++;
            *extend_left_gaps += gapsize - 1;
        }
        else if (end == right) {
            (*open_right_gaps)++;
            *extend_right_gaps += gapsize - 1;
        }
        else {
            (*open_internal_gaps)++;
            *extend_internal_gaps += gapsize - 1; \
        } \
        if (score_function) {
            double value;
            PyObject* result;
            if (strand)
                result = PyObject_CallFunction(score_function, "ii",
                                               right - start, gapsize);
            else
                result = PyObject_CallFunction(score_function,
                                               "ii", start, gapsize);
            if (result == NULL) return false;
            value = PyFloat_AsDouble(result);
            Py_DECREF(result);
            if (value == -1.0 && PyErr_Occurred()) return false;
            *gap_score += value;
        }
        *path = direction;
    }
    return true;
}

static inline Py_buffer* get_buffer(Py_buffer *sequence, int** i, char** b)
{
    *i = NULL;
    *b = NULL;
    if (sequence->buf) {
        switch (sequence->format[0]) {
            case 'i':
            case 'l':
                *i = sequence->buf;
                break;
            case 'c':
            case 'b':
            case 'B':
                *b = sequence->buf;
                break;
        }
    }
    return sequence;
}

static PyObject* 
AlignmentCounts_new(PyTypeObject *type, PyObject *args, PyObject *keywords)
{
    Py_ssize_t jA, jB;
    Py_ssize_t n = 0;
    PyObject* sequence;
    Aligner* aligner = NULL;
    PyObject* sequences;
    Py_buffer* sequence_buffers = NULL;
    Py_buffer coordinates = {0};
    Py_buffer strands = {0};
    AlignmentCounts* counts = NULL;
    int wildcard = -1;
    Py_buffer substitution_matrix = {0};
    PyObject* argument = NULL;

    Py_ssize_t k, lA, lB;
    int cA, cB;

    Py_buffer* sequenceA;
    Py_buffer* sequenceB;
    int* iA = NULL;
    int* iB = NULL;
    char* bA = NULL;
    char* bB = NULL;
    bool strandA;
    bool strandB;

    int* mapping = NULL;
    Py_ssize_t m = 0;

    PyObject* insertion_score_function = NULL;
    PyObject* deletion_score_function = NULL;

    Py_ssize_t open_left_insertions = 0, extend_left_insertions = 0;
    Py_ssize_t open_left_deletions = 0, extend_left_deletions = 0;
    Py_ssize_t open_internal_insertions = 0, extend_internal_insertions = 0;
    Py_ssize_t open_internal_deletions = 0, extend_internal_deletions = 0;
    Py_ssize_t open_right_insertions = 0, extend_right_insertions = 0;
    Py_ssize_t open_right_deletions = 0, extend_right_deletions = 0;

    Py_ssize_t aligned = 0;
    Py_ssize_t identities = 0;
    Py_ssize_t mismatches = 0;
    Py_ssize_t positives = -1;

    double gap_score = 0.0;
    double substitution_score = 0.0;

    Py_ssize_t n_columns;
    Py_ssize_t row_stride;
    Py_ssize_t column_stride;

    int path = 0;

    Py_ssize_t leftA, leftB;
    Py_ssize_t rightA, rightB;
    Py_ssize_t startA, startB;
    Py_ssize_t endA, endB;
    Py_ssize_t* buffer;

    PyObject* oA = NULL;
    PyObject* oB = NULL;

    static char *kwlist[] = {"sequences", "coordinates", "strands", "argument", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywords, "O!O&O&|O", kwlist,
                                     &PyList_Type, &sequences,
                                     coordinates_converter, &coordinates,
                                     strands_converter , &strands,
                                     &argument))
        return 0;

    if (argument == NULL) {
    }
    else if (PyObject_TypeCheck(argument, Aligner_Type)) {
        aligner = (Aligner*)argument;
        if (aligner->substitution_matrix.obj) {
            substitution_matrix = aligner->substitution_matrix;
            Py_INCREF(substitution_matrix.obj);
        }
        wildcard = aligner->wildcard;
    }
    else if (PyUnicode_Check(argument)) {
        if (PyUnicode_READY(argument) == -1) goto exit;
        if (PyUnicode_GET_LENGTH(argument) != 1) {
            PyErr_SetString(PyExc_ValueError,
                            "wildcard should be a single character, or None");
            goto exit;
        }
        wildcard = PyUnicode_READ_CHAR(argument, 0);
    }
    else {
        if (!substitution_matrix_converter(argument,
                                           &substitution_matrix)) {
            if (!Aligner_Type) 
                PyErr_SetString(PyExc_RuntimeError, "Aligner_Type is NULL");
            goto exit;
        }
    }

    n_columns = coordinates.shape[1];
    row_stride = coordinates.strides[0] / sizeof(Py_ssize_t);
    column_stride = coordinates.strides[1] / sizeof(Py_ssize_t);
    buffer = coordinates.buf;

    if (substitution_matrix.obj) positives = 0;

    n = PyList_GET_SIZE(sequences);
    if (n != coordinates.shape[0]) {
        PyErr_SetString(PyExc_ValueError,
            "number of rows in coordinates must equal the number of sequences");
        goto exit;
    }
    if (n != strands.shape[0]) {
        PyErr_SetString(PyExc_ValueError,
            "size of strands must equal the number of sequences");
        goto exit;
    }

    sequence_buffers = PyMem_Calloc(n, sizeof(Py_buffer));
    if (!sequence_buffers) goto exit;

    for (k = 0; k < n; k++) {
        sequence = PyList_GET_ITEM(sequences, k);
        if (!sequence_converter(sequence, &sequence_buffers[k])) {
            n = k;
            goto exit;
        }
    }

    if (aligner) {
        insertion_score_function = aligner->insertion_score_function;
        deletion_score_function = aligner->deletion_score_function;
    }

    counts = (AlignmentCounts*)PyType_GenericAlloc(&AlignmentCounts_Type, 0);
    if (!counts) goto exit;

    if (substitution_matrix.obj) {
        m = substitution_matrix.shape[0];
        if (PyObject_IsInstance(substitution_matrix.obj,
                               (PyObject*)Array_Type)) {
            PyTypeObject* basetype = Array_Type->tp_base;
            Fields* fields = (Fields*)((intptr_t)substitution_matrix.obj + basetype->tp_basicsize);
            Py_buffer* mapping_buffer = &fields->mapping;
            mapping = mapping_buffer->buf;
            if (mapping) m = mapping_buffer->len / mapping_buffer->itemsize;
        }
    }

    for (jA = 0; jA < n; jA++) {
        oA = NULL;
        sequenceA = get_buffer(&sequence_buffers[jA], &iA, &bA);
        strandA = ((bool*)(strands.buf))[jA];
        for (jB = jA + 1; jB < n; jB++) {
            oB = NULL;
            sequenceB = get_buffer(&sequence_buffers[jB], &iB, &bB);
            strandB = ((bool*)(strands.buf))[jB];
            leftA = buffer[jA * row_stride + 0];
            leftB = buffer[jB * row_stride + 0];
            rightA = buffer[jA * row_stride + (n_columns - 1) * column_stride];
            rightB = buffer[jB * row_stride + (n_columns - 1) * column_stride];
            startA = leftA;
            startB = leftB;
            for (k = 1; k < n_columns; k++) {
                endA = buffer[jA * row_stride + k * column_stride];
                endB = buffer[jB * row_stride + k * column_stride];
                if (startA == endA && startB == endB) {
                }
                else if (startA == endA) {
                    if (!add_gaps(&path, HORIZONTAL, strandA,
                                  insertion_score_function,
                                  startA, endA, leftA, rightA, endB - startB,
                                  &open_left_insertions,
                                  &extend_left_insertions,
                                  &open_internal_insertions,
                                  &extend_internal_insertions,
                                  &open_right_insertions,
                                  &extend_right_insertions,
                                  &gap_score)) goto error;
                }
                else if (startB == endB) {
                    if (!add_gaps(&path, VERTICAL, strandB,
                                  deletion_score_function,
                                  startB, endB, leftB, rightB, endA - startA,
                                  &open_left_deletions,
                                  &extend_left_deletions,
                                  &open_internal_deletions,
                                  &extend_internal_deletions,
                                  &open_right_deletions,
                                  &extend_right_deletions,
                                  &gap_score)) goto error;
                }
                else if (sequenceA->obj == NULL || sequenceB->obj == NULL) {
                    path = DIAGONAL;
                    aligned += endA - startA;
                }
                else {
                    path = DIAGONAL;
                    aligned += endA - startA;
                    if (iA == NULL && bA == NULL) {
                        oA = get_lazy_data(sequenceA->obj, startA, endA, jA, &bA);
                        if (!oA) goto error;
                    }
                    if (iB == NULL && bB == NULL) {
                        oB = get_lazy_data(sequenceB->obj, startB, endB, jB, &bB);
                        if (!oB) goto error;
                    }
                    if (substitution_matrix.obj == NULL) {
                        if (iA && iB) {
                            for (lA = startA, lB = startB;
                                 lA < endA && lB < endB;
                                 lA++, lB++) {
                                add_identities_mismatches(iA[lA],
                                                          iB[lB],
                                                          wildcard,
                                                          &identities,
                                                          &mismatches);
                            }
                        }
                        else if (iA && bB) {
                            for (lA = startA, lB = startB;
                                 lA < endA && lB < endB;
                                 lA++, lB++) {
                                add_identities_mismatches(iA[lA],
                                                          (int) bB[lB],
                                                          wildcard,
                                                          &identities,
                                                          &mismatches);
                            }
                        }
                        else if (iB && bA) {
                            for (lA = startA, lB = startB;
                                 lA < endA && lB < endB;
                                 lA++, lB++) {
                                add_identities_mismatches((int) bA[lA],
                                                          iB[lB],
                                                          wildcard,
                                                          &identities,
                                                          &mismatches);
                            }
                        }
                        else if (bA && bB) {
                            for (lA = startA, lB = startB;
                                 lA < endA && lB < endB;
                                 lA++, lB++) {
                                add_identities_mismatches((int) bA[lA],
                                                          (int) bB[lB],
                                                          wildcard,
                                                          &identities,
                                                          &mismatches);
                            }
                        }
                    }
                    else if (mapping == NULL) {
                        if (iA && iB) {
                            for (lA = startA, lB = startB;
                                 lA < endA && lB < endB;
                                 lA++, lB++) {
                                cA = iA[lA];
                                cB = iB[lB];
                                if (!check_indices(cA, jA, lA, m)) goto error;
                                if (!check_indices(cB, jB, lB, m)) goto error;
                                add_identities_mismatches_score(cA, cB,
                                    wildcard,
                                    &substitution_matrix,
                                    &identities,
                                    &mismatches,
                                    &positives,
                                    &substitution_score);
                            }
                        } else if (iA && bB) {
                            for (lA = startA, lB = startB;
                                 lA < endA && lB < endB;
                                 lA++, lB++) {
                                cA = iA[lA];
                                cB = (int)bB[lB];
                                if (!check_indices(cA, jA, lA, m)) goto error;
                                if (!check_indices(cB, jB, lB, m)) goto error;
                                add_identities_mismatches_score(cA, cB,
                                    wildcard,
                                    &substitution_matrix,
                                    &identities,
                                    &mismatches,
                                    &positives,
                                    &substitution_score);
                            }
                        } else if (iB && bA) {
                            for (lA = startA, lB = startB;
                                 lA < endA && lB < endB;
                                 lA++, lB++) {
                                cA = (int) bA[lA];
                                cB = iB[lB];
                                if (!check_indices(cA, jA, lA, m)) goto error;
                                if (!check_indices(cB, jB, lB, m)) goto error;
                                add_identities_mismatches_score(cA, cB,
                                    wildcard,
                                    &substitution_matrix,
                                    &identities,
                                    &mismatches,
                                    &positives,
                                    &substitution_score);
                            }
                        }
                        else if (bA && bB) {
                            for (lA = startA, lB = startB;
                                 lA < endA && lB < endB;
                                 lA++, lB++) {
                                cA = (int) bA[lA];
                                cB = (int) bB[lB];
                                if (!check_indices(cA, jA, lA, m)) goto error;
                                if (!check_indices(cB, jB, lB, m)) goto error;
                                add_identities_mismatches_score(cA, cB,
                                    wildcard,
                                    &substitution_matrix,
                                    &identities,
                                    &mismatches,
                                    &positives,
                                    &substitution_score);
                            }
                        }
                    }
                    else {
                        if (iA && iB) {
                            for (lA = startA, lB = startB;
                                 lA < endA && lB < endB;
                                 lA++, lB++) {
                                cA = iA[lA];
                                cB = iB[lB];
                                if (!check_indices(cA, jA, lA, m)) goto error;
                                if (!check_indices(cB, jB, lB, m)) goto error;
                                if (!map_indices(&cA, &cB, mapping)) goto error;
                                add_identities_mismatches_score(cA, cB,
                                    wildcard,
                                    &substitution_matrix,
                                    &identities,
                                    &mismatches,
                                    &positives,
                                    &substitution_score);
                            }
                        } else if (iA && bB) {
                            for (lA = startA, lB = startB;
                                 lA < endA && lB < endB;
                                 lA++, lB++) {
                                cA = iA[lA];
                                cB = (int)bB[lB];
                                if (!check_indices(cA, jA, lA, m)) goto error;
                                if (!check_indices(cB, jB, lB, m)) goto error;
                                if (!map_indices(&cA, &cB, mapping)) goto error;
                                add_identities_mismatches_score(cA, cB,
                                    wildcard,
                                    &substitution_matrix,
                                    &identities,
                                    &mismatches,
                                    &positives,
                                    &substitution_score);
                            }
                        } else if (iB && bA) {
                            for (lA = startA, lB = startB;
                                 lA < endA && lB < endB;
                                 lA++, lB++) {
                                cA = (int) bA[lA];
                                cB = iB[lB];
                                if (!check_indices(cA, jA, lA, m)) goto error;
                                if (!check_indices(cB, jB, lB, m)) goto error;
                                if (!map_indices(&cA, &cB, mapping)) goto error;
                                add_identities_mismatches_score(cA, cB,
                                    wildcard,
                                    &substitution_matrix,
                                    &identities,
                                    &mismatches,
                                    &positives,
                                    &substitution_score);
                            }
                        }
                        else if (bA && bB) {
                            for (lA = startA, lB = startB;
                                 lA < endA && lB < endB;
                                 lA++, lB++) {
                                cA = (int) bA[lA];
                                cB = (int) bB[lB];
                                if (!check_indices(cA, jA, lA, m)) goto error;
                                if (!check_indices(cB, jB, lB, m)) goto error;
                                if (!map_indices(&cA, &cB, mapping)) goto error;
                                add_identities_mismatches_score(cA, cB,
                                    wildcard,
                                    &substitution_matrix,
                                    &identities,
                                    &mismatches,
                                    &positives,
                                    &substitution_score);
                            }
                        }
                    }
                    reset_lazy_data(oA, &bA);
                    reset_lazy_data(oB, &bB);
                }
                startA = endA;
                startB = endB;
            }
        }
    }

    counts->open_left_insertions = open_left_insertions;
    counts->extend_left_insertions = extend_left_insertions;
    counts->open_left_deletions = open_left_deletions;
    counts->extend_left_deletions = extend_left_deletions;
    counts->open_internal_insertions = open_internal_insertions;
    counts->extend_internal_insertions = extend_internal_insertions;
    counts->open_internal_deletions = open_internal_deletions;
    counts->extend_internal_deletions = extend_internal_deletions;
    counts->open_right_insertions = open_right_insertions;
    counts->extend_right_insertions = extend_right_insertions;
    counts->open_right_deletions = open_right_deletions;
    counts->extend_right_deletions = extend_right_deletions;
    counts->aligned = aligned;
    counts->identities = identities;
    counts->mismatches = mismatches;
    counts->positives = positives;

    if (substitution_matrix.obj == NULL) {
        if (aligner) {
            substitution_score = aligner->match * counts->identities
                               + aligner->mismatch * counts->mismatches;
        }
        else substitution_score = Py_NAN;
    }
    counts->substitution_score = substitution_score;

    if (aligner) {
        if (aligner->insertion_score_function == NULL) {
            gap_score += counts->open_left_insertions
                       * aligner->open_left_insertion_score
                       + counts->extend_left_insertions
                       * aligner->extend_left_insertion_score
                       + counts->open_internal_insertions
                       * aligner->open_internal_insertion_score
                       + counts->extend_internal_insertions
                       * aligner->extend_internal_insertion_score
                       + counts->open_right_insertions
                       * aligner->open_right_insertion_score
                       + counts->extend_right_insertions
                       * aligner->extend_right_insertion_score;
        }
        if (aligner->deletion_score_function == NULL) {
            gap_score += counts->open_left_deletions
                       * aligner->open_left_deletion_score
                       + counts->extend_left_deletions
                       * aligner->extend_left_deletion_score
                       + counts->open_internal_deletions
                       * aligner->open_internal_deletion_score
                       + counts->extend_internal_deletions
                       * aligner->extend_internal_deletion_score
                       + counts->open_right_deletions
                       * aligner->open_right_deletion_score
                       + counts->extend_right_deletions
                       * aligner->extend_right_deletion_score;
        }
    }
    else gap_score = Py_NAN;
    counts->gap_score = gap_score;

    goto exit;

error:
    Py_XDECREF(oA);
    Py_XDECREF(oB);
    Py_DECREF(counts);
    counts = NULL;

exit:
    if (sequence_buffers) {
        for (k = 0; k < n; k++)
            if (sequence_buffers[k].buf) PyBuffer_Release(&sequence_buffers[k]);
        PyMem_Free(sequence_buffers);
    }
    coordinates_converter(NULL, &coordinates);
    strands_converter(NULL, &strands);
    substitution_matrix_converter(NULL, &substitution_matrix);

    return (PyObject*) counts;
}

static PyTypeObject AlignmentCounts_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "AlignmentCounts",
    .tp_basicsize = sizeof(AlignmentCounts),
    .tp_repr = (reprfunc)AlignmentCounts_repr,
    .tp_str = (reprfunc)AlignmentCounts_str,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_doc = AlignmentCounts_doc,
    .tp_getset = AlignmentCounts_getset,
    .tp_new = (newfunc)AlignmentCounts_new,
};


static char _alignmentcounts__doc__[] =
"C extension module implementing the AlignmentCounts class";

/* Module definition */

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_alignmentcounts",
    .m_doc = _alignmentcounts__doc__,
    .m_size = -1,
};

PyObject *
PyInit__alignmentcounts(void)
{
    PyObject* module;

    if (PyType_Ready(&AlignmentCounts_Type) < 0)
        return NULL;

    module = PyModule_Create(&moduledef);
    if (!module) return NULL;

    PyObject *mod;
    mod = PyImport_ImportModule("Bio.Align.substitution_matrices._arraycore");
    if (!mod) {
        Py_DECREF(module);
        return NULL;
    }
    Array_Type = (PyTypeObject*)PyObject_GetAttrString(mod, "Array");
    Py_DECREF(mod);

    mod = PyImport_ImportModule("Bio.Align._pairwisealigner");
    if (!mod) {
        Py_DECREF(module);
        return NULL;
    }
    Aligner_Type = (PyTypeObject*)PyObject_GetAttrString(mod, "PairwiseAligner");
    Py_DECREF(mod);
    if (!Aligner_Type) {
        Py_DECREF(module);
        return NULL;
    }

    Py_INCREF(&AlignmentCounts_Type);
    if (PyModule_AddObject(module,
                           "AlignmentCounts",
                           (PyObject*) &AlignmentCounts_Type) < 0) {
        Py_DECREF(module);
        Py_DECREF(&AlignmentCounts_Type);
        return NULL;
    }
    
    return module;
}
