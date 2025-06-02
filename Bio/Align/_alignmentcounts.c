/* Copyright 2018-2019 by Michiel de Hoon.  All rights reserved.
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */



#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <float.h>
#include <stdbool.h>
#include "_pairwisealigner.h"
#include "substitution_matrices/_arraycore.h"

#define HORIZONTAL 0x1
#define VERTICAL 0x2
#define DIAGONAL 0x4


#define MISSING_LETTER -1



static PyTypeObject* Aligner_Type = NULL;


static Array_get_mapping_buffer_signature Array_get_mapping_buffer;
/* this will be set when initializing the module */


static int _map_indices(Py_buffer* view, Py_buffer* substitution_matrix) {
    Py_ssize_t i;
    int index;
    int* indices = view->buf;
    const Py_ssize_t n = view->len / view->itemsize;
    Py_buffer buffer;
    Array_get_mapping_buffer(substitution_matrix->obj, &buffer);
    if (buffer.obj) {
        PyBuffer_Release(&buffer);
    }
    else {
        const Py_ssize_t m = substitution_matrix->shape[0];
        for (i = 0; i < n; i++) {
            index = indices[i];
            if (index >= m) {
                PyErr_Format(PyExc_ValueError,
                             "sequence item %zd is out of bound"
                             " (%d, should be < %zd)", i, index, m);
                return 0;
            }
        }
    }
    return 1;
}

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
AlignmentCounts_repr(AlignmentCounts* self)
{
    PyObject* representation = NULL;
    const double score = self->gap_score + self->substitution_score;
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
    if (isnan(score)) {
        if (positives == -1)
            representation = PyUnicode_FromFormat("AlignmentCounts("
"open_left_insertions=%zd, extend_left_insertions=%zd, "
"open_left_deletions=%zd, extend_left_deletions=%zd, "
"open_internal_insertions=%zd, extend_internal_insertions=%zd, "
"open_internal_deletions=%zd, extend_internal_deletions=%zd, "
"open_right_insertions=%zd, extend_right_insertions=%zd, "
"open_right_deletions=%zd, extend_right_deletions=%zd, "
"aligned=%zd, identities=%zd, mismatches=%zd)",
                open_left_insertions,
                extend_left_insertions,
                open_left_deletions,
                extend_left_deletions,
                open_internal_insertions,
                extend_internal_insertions,
                open_internal_deletions,
                extend_internal_deletions,
                open_right_insertions,
                extend_right_insertions,
                open_right_deletions,
                extend_right_deletions,
                aligned,
                identities,
                mismatches);
        else
            representation = PyUnicode_FromFormat("AlignmentCounts("
"open_left_insertions=%zd, extend_left_insertions=%zd, "
"open_left_deletions=%zd, extend_left_deletions=%zd, "
"open_internal_insertions=%zd, extend_internal_insertions=%zd, "
"open_internal_deletions=%zd, extend_internal_deletions=%zd, "
"open_right_insertions=%zd, extend_right_insertions=%zd, "
"open_right_deletions=%zd, extend_right_deletions=%zd, "
"aligned=%zd, identities=%zd, mismatches=%zd, positives=%zd)",
                open_left_insertions,
                extend_left_insertions,
                open_left_deletions,
                extend_left_deletions,
                open_internal_insertions,
                extend_internal_insertions,
                open_internal_deletions,
                extend_internal_deletions,
                open_right_insertions,
                extend_right_insertions,
                open_right_deletions,
                extend_right_deletions,
                aligned,
                identities,
                mismatches,
                positives);
    }
    else {
        char* score_text = PyOS_double_to_string(score, 'f', 6, 0, NULL);
        if (!score_text) return NULL;
        if (positives == -1)
            representation = PyUnicode_FromFormat("AlignmentCounts("
"open_left_insertions=%zd, extend_left_insertions=%zd, "
"open_left_deletions=%zd, extend_left_deletions=%zd, "
"open_internal_insertions=%zd, extend_internal_insertions=%zd, "
"open_internal_deletions=%zd, extend_internal_deletions=%zd, "
"open_right_insertions=%zd, extend_right_insertions=%zd, "
"open_right_deletions=%zd, extend_right_deletions=%zd, "
"aligned=%zd, identities=%zd, mismatches=%zd, score=%s)",
                open_left_insertions,
                extend_left_insertions,
                open_left_deletions,
                extend_left_deletions,
                open_internal_insertions,
                extend_internal_insertions,
                open_internal_deletions,
                extend_internal_deletions,
                open_right_insertions,
                extend_right_insertions,
                open_right_deletions,
                extend_right_deletions,
                aligned,
                identities,
                mismatches,
                score_text);
        else
            representation = PyUnicode_FromFormat("AlignmentCounts("
"open_left_insertions=%zd, extend_left_insertions=%zd, "
"open_left_deletions=%zd, extend_left_deletions=%zd, "
"open_internal_insertions=%zd, extend_internal_insertions=%zd, "
"open_internal_deletions=%zd, extend_internal_deletions=%zd, "
"open_right_insertions=%zd, extend_right_insertions=%zd, "
"open_right_deletions=%zd, extend_right_deletions=%zd, "
"aligned=%zd, identities=%zd, mismatches=%zd, positives=%zd, score=%s)",
                open_left_insertions,
                extend_left_insertions,
                open_left_deletions,
                extend_left_deletions,
                open_internal_insertions,
                extend_internal_insertions,
                open_internal_deletions,
                extend_internal_deletions,
                open_right_insertions,
                extend_right_insertions,
                open_right_deletions,
                extend_right_deletions,
                aligned,
                identities,
                mismatches,
                positives,
                score_text);
        PyMem_Free(score_text);
    }
    return representation;
}

static PyObject*
AlignmentCounts_str(AlignmentCounts* self)
{
    PyObject* str;
    const double score = self->gap_score + self->substitution_score;
    const Py_ssize_t aligned = self->aligned;
    const Py_ssize_t identities = self->identities;
    const Py_ssize_t mismatches = self->mismatches;
    const Py_ssize_t positives = self->positives;
    const Py_ssize_t insertions = self->open_left_insertions
                                + self->extend_left_insertions
                                + self->open_internal_insertions
                                + self->extend_internal_insertions
                                + self->open_right_insertions
                                + self->extend_right_insertions;
    const Py_ssize_t deletions = self->open_left_deletions
                               + self->extend_left_deletions
                               + self->open_internal_deletions
                               + self->extend_internal_deletions
                               + self->open_right_deletions
                               + self->extend_right_deletions;
    const Py_ssize_t gaps = insertions + deletions;
    if (isnan(score)) {
        if (positives == -1) {
            const char text[] = "Alignment with "
"%zd aligned letters (%zd identities, %zd mismatches) and "
"%zd gaps (%zd insertions, %zd deletions)";
            str = PyUnicode_FromFormat(text, aligned, identities, mismatches,
                                       gaps, insertions, deletions);
        }
        else {
            const char text[] = "Alignment with "
"%zd aligned letters (%zd identities, %zd mismatches, %zd positives) and "
"%zd gaps (%zd insertions, %zd deletions)";
            str = PyUnicode_FromFormat(text, aligned,
                                       identities, mismatches, positives,
                                       gaps, insertions, deletions);
        }
    }
    else {
        char* score_text = PyOS_double_to_string(score, 'f', 6, 0, NULL);
        if (!score_text) return NULL;
        if (positives == -1) {
            const char text[] = "Alignment with "
"%zd aligned letters (%zd identities, %zd mismatches) and "
"%zd gaps (%zd insertions, %zd deletions); alignment score = %s";
            str = PyUnicode_FromFormat(text, aligned, identities, mismatches,
                                       gaps, insertions, deletions, score_text);
        }
        else {
            const char text[] = "Alignment with "
"%zd aligned letters (%zd identities, %zd mismatches, %zd positives) and "
"%zd gaps (%zd insertions, %zd deletions); alignment score = %s";
            str = PyUnicode_FromFormat(text, aligned,
                                       identities, mismatches, positives,
                                       gaps, insertions, deletions, score_text);
        }
        PyMem_Free(score_text);
    }
    return str;
}

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

static char AlignmentCounts_open_gaps__doc__[] =  "number of geps opened in the alignment";

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

static char AlignmentCounts_extend_gaps__doc__[] =  "number of gep extensions in the alignment";

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

static char AlignmentCounts_gaps__doc__[] =  "total gap length";

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
    {"deletions",
        (getter)AlignmentCounts_get_deletions, NULL,
        AlignmentCounts_deletions__doc__, NULL},
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

static char AlignmentCounts_doc[] =
"AlignmentCounts objects store the number of aligned characters, identities,\n"
"mismatches, gaps, and score of an alignment.\n";

static PyTypeObject AlignmentCounts_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "AlignmentCounts",              /* tp_name */
    sizeof(AlignmentCounts),        /* tp_basicsize */
    0,                              /* tp_itemsize */
    0,                              /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_compare */
    (reprfunc)AlignmentCounts_repr, /* tp_repr */
    0,                              /* tp_as_number */
    0,                              /* tp_as_sequence */
    0,                              /* tp_as_mapping */
    0,                              /* tp_hash */
    0,                              /* tp_call */
    (reprfunc)AlignmentCounts_str,  /* tp_str */
    0,                              /* tp_getattro */
    0,                              /* tp_setattro */
    0,                              /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    AlignmentCounts_doc,            /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    0,                              /* tp_iter */
    0,                              /* tp_iternext */
    0,                              /* tp_methods */
    0,                              /* tp_members */
    AlignmentCounts_getset,         /* tp_getset */
};

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


/* ----------------- alignment algorithms ----------------- */


static int
sequence_converter(PyObject* argument, void* pointer)
{
    Py_buffer* view = pointer;
    const int flag = PyBUF_FORMAT | PyBUF_C_CONTIGUOUS;

    if (argument == NULL) {
        PyBuffer_Release(view);
        return 1;
    }

    if (PyObject_GetBuffer(argument, view, flag) != 0) {
        PyErr_SetString(PyExc_TypeError, "argument is not a sequence");
        return 0;
    }
    if (view->ndim != 1) {
        PyErr_Format(PyExc_ValueError,
                     "sequence has incorrect rank (%d expected 1)", view->ndim);
        return 0;
    }
    if (view->len == 0) {
        PyErr_SetString(PyExc_ValueError, "sequence has zero length");
        return 0;
    }
    if (strcmp(view->format, "i") != 0 && strcmp(view->format, "l") != 0) {
        PyErr_Format(PyExc_ValueError,
                     "sequence has incorrect data type '%s'", view->format);
        return 0;
    }
    if (view->itemsize != sizeof(int)) {
        PyErr_Format(PyExc_ValueError,
                    "sequence has unexpected item byte size "
                    "(%ld, expected %ld)", view->itemsize, sizeof(int));
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

/* Module definition */

static char _alignmentcounts__doc__[] =
"C extension module implementing the AlignmentCounts class";

static const char _calculate__doc__[] = "calculate the matches, mismatches, gaps, and score of the alignment";

static PyObject*
_calculate(PyObject* self, PyObject* args, PyObject* keywords)
{
    Py_ssize_t i;
    Py_ssize_t n = 0;
    PyObject* sequence;
    Aligner* aligner = NULL;
    PyObject* sequences;
    Py_buffer* buffers = NULL;
    Py_buffer coordinates = {0};
    Py_buffer strands = {0};
    AlignmentCounts* counts = NULL;
    int wildcard = -1;
    Py_buffer substitution_matrix = {0};
    PyObject* argument = NULL;

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
        const int ok = substitution_matrix_converter(argument,
                                                     &substitution_matrix);
        if (!ok) {
            if (!Aligner_Type) 
                PyErr_SetString(PyExc_RuntimeError, "Aligner_Type is NULL");
            goto exit;
        }
    }

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

    buffers = PyMem_Calloc(n, sizeof(Py_buffer));
    if (!buffers) goto exit;

    for (i = 0; i < n; i++) {
        sequence = PyList_GET_ITEM(sequences, i);
        if (sequence_converter(sequence, &buffers[i])) {
            if (substitution_matrix.obj) {
                if (!_map_indices(&buffers[i], &substitution_matrix)) goto exit;
            }
        }
        else {
            PyErr_Clear();  // to clear the exception raised by PyObject_GetBuffer
            buffers[i].buf = NULL;
            if (PySequence_Check(sequence)) buffers[i].obj = sequence;
            else if (sequence == Py_None) buffers[i].obj = NULL;
            else {
                n = i;
                goto exit;
            }
        }
    }

    counts = (AlignmentCounts*)PyType_GenericAlloc(&AlignmentCounts_Type, 0);
    if (!counts) goto exit;

    int* mapping = NULL;
    int m = 0;

    if (substitution_matrix.obj) {
        Py_buffer mapping_buffer;
        Array_get_mapping_buffer(substitution_matrix.obj, &mapping_buffer);
        if (mapping_buffer.obj) {
            mapping = mapping_buffer.buf;
            m = mapping_buffer.len / mapping_buffer.itemsize;
            PyBuffer_Release(&mapping_buffer);
        }
        else m = substitution_matrix.shape[0];
    }

    PyObject* insertion_score_function = NULL;
    PyObject* deletion_score_function = NULL;
    if (aligner) {
        insertion_score_function = aligner->insertion_score_function;
        deletion_score_function = aligner->deletion_score_function;
    }

    Py_ssize_t j, k, l1, l2;
    int cA, cB;

    Py_buffer* sequenceA;
    Py_buffer* sequenceB;
    int* sA;
    int* sB;
    bool strandA;
    bool strandB;

    Py_ssize_t open_left_insertions = 0, extend_left_insertions = 0;
    Py_ssize_t open_left_deletions = 0, extend_left_deletions = 0;
    Py_ssize_t open_internal_insertions = 0, extend_internal_insertions = 0;
    Py_ssize_t open_internal_deletions = 0, extend_internal_deletions = 0;
    Py_ssize_t open_right_insertions = 0, extend_right_insertions = 0;
    Py_ssize_t open_right_deletions = 0, extend_right_deletions = 0;
    Py_ssize_t aligned = 0;
    Py_ssize_t identities = 0;
    Py_ssize_t mismatches = 0;
    Py_ssize_t positives = substitution_matrix.obj ? 0 : -1;
    double gap_score = 0.0;
    double substitution_score = 0.0;

    const Py_ssize_t shape2 = coordinates.shape[1];
    const Py_ssize_t stride1 = coordinates.strides[0] / sizeof(Py_ssize_t);
    const Py_ssize_t stride2 = coordinates.strides[1] / sizeof(Py_ssize_t);
    Py_ssize_t left1, left2;
    Py_ssize_t right1, right2;
    Py_ssize_t start1, start2;
    Py_ssize_t end1, end2;
    Py_ssize_t* buffer = coordinates.buf;

    PyObject* oA = NULL;
    PyObject* oB = NULL;

    int path = 0;

    for (i = 0; i < n; i++) {
        sequenceA = &buffers[i];
        sA = sequenceA->buf;
        strandA = ((bool*)(strands.buf))[i];
        for (j = i + 1; j < n; j++) {
            sequenceB = &buffers[j];
            sB = sequenceB->buf;
            strandB = ((bool*)(strands.buf))[j];
            left1 = buffer[i * stride1 + 0];
            left2 = buffer[j * stride1 + 0];
            right1 = buffer[i * stride1 + (shape2 - 1) * stride2];
            right2 = buffer[j * stride1 + (shape2 - 1) * stride2];
            start1 = left1;
            start2 = left2;
            for (k = 1; k < shape2; k++) {
                end1 = buffer[i * stride1 + k * stride2];
                end2 = buffer[j * stride1 + k * stride2];
                if (start1 == end1 && start2 == end2) {
                }
                else if (start1 == end1) {
                    if (path == HORIZONTAL) {
                        if (start1 == left1)
                            extend_left_insertions += end2 - start2;
                        else if (end1 == right1)
                            extend_right_insertions += end2 - start2;
                        else
                            extend_internal_insertions += end2 - start2;
                    }
                    else {
                        if (start1 == left1) {
                            open_left_insertions++;
                            extend_left_insertions += end2 - start2 - 1;
                        }
                        else if (end1 == right1) {
                            open_right_insertions++;
                            extend_right_insertions += end2 - start2 - 1;
                        }
                        else {
                            open_internal_insertions++;
                            extend_internal_insertions += end2 - start2 - 1;
                        }
                        if (insertion_score_function) {
                            double value;
                            PyObject* result;
                            if (strandA)
                                result = PyObject_CallFunction(insertion_score_function,
                                                               "ii", right1 - start1, end2 - start2);
                            else
                                result = PyObject_CallFunction(insertion_score_function,
                                                               "ii", start1, end2 - start2);
                            if (result == NULL) goto error;
                            value = PyFloat_AsDouble(result);
                            Py_DECREF(result);
                            if (value == -1.0 && PyErr_Occurred()) goto error;
                            gap_score += value;
                        }
                        path = HORIZONTAL;
                    }
                }
                else if (start2 == end2) {
                    if (path == VERTICAL) {
                        if (start2 == left2)
                            extend_left_deletions += end1 - start1;
                        else if (end2 == right2)
                            extend_right_deletions += end1 - start1;
                        else
                            extend_internal_deletions += end1 - start1;
                    }
                    else {
                        if (start2 == left2) {
                            open_left_deletions++;
                            extend_left_deletions += end1 - start1 - 1;
                        }
                        else if (end2 == right2) {
                            open_right_deletions++;
                            extend_right_deletions += end1 - start1 - 1;
                        }
                        else {
                            open_internal_deletions++;
                            extend_internal_deletions += end1 - start1 - 1;
                        }
                        if (deletion_score_function) {
                            double value;
                            PyObject* result;
                            if (strandB)
                                result = PyObject_CallFunction(deletion_score_function,
                                                               "ii", right2 - start2, end1 - start1);
                            else
                                result = PyObject_CallFunction(deletion_score_function,
                                                               "ii", start2, end1 - start1);
                            if (result == NULL) return 0;
                            value = PyFloat_AsDouble(result);
                            Py_DECREF(result);
                            if (value == -1.0 && PyErr_Occurred()) goto error;
                            gap_score += value;
                        }
                        path = VERTICAL;
                    }
                }
                else if (sequenceA->obj == NULL || sequenceB->obj == NULL) {
                    path = DIAGONAL;
                    aligned += end1 - start1;
                }
                else {
                    char* bA = NULL;
                    char* bB = NULL;
                    path = DIAGONAL;
                    aligned += end1 - start1;
                    if (sA == NULL) {
                        if (PyBytes_Check(sequenceA->obj)) {
                            bA = PyBytes_AS_STRING(sequenceA->obj) + start1;
                            oA = NULL;
                        }
                        else {
                            oA = PySequence_GetSlice(sequenceA->obj, start1, end1);
                            if (!oA) goto error;
                            if (PyBytes_Check(oA)) {
                                if (PyBytes_GET_SIZE(oA) != end1 - start1) {
                                    PyErr_Format(PyExc_ValueError,
                                        "alignment.sequences[%d][%d:%d] did not return a bytes object of size %d",
                                        i, start1, end1, end1 - start1);
                                    goto error;
                                }
                                bA = PyBytes_AS_STRING(oA);
                            }
                        }
                    }
                    if (sB == NULL) {
                        if (PyBytes_Check(sequenceB->obj)) {
                            bB = PyBytes_AS_STRING(sequenceB->obj) + start2;
                            oB = NULL;
                        }
                        else {
                            oB = PySequence_GetSlice(sequenceB->obj, start2, end2);
                            if (!oB) goto error;
                            if (PyBytes_Check(oB)) {
                                if (PyBytes_GET_SIZE(oB) != end2 - start2) {
                                    PyErr_Format(PyExc_ValueError,
                                        "alignment.sequences[%d[%d:%d] did not return a bytes object of size %d",

                                        j, start2, end2, end2 - start2);
                                    goto error;
                                }
                                bB = PyBytes_AS_STRING(oB);
                            }
                        }
                    }
                    if (substitution_matrix.obj == NULL) {
                        if (sA && sB) {
                            for (l1 = start1, l2 = start2;
                                 l1 < end1 && l2 < end2;
                                 l1++, l2++) {
                                cA = sA[l1];
                                cB = sB[l2];
                                if (cA == wildcard || cB == wildcard) ;
                                else if (cA == cB) identities++;
                                else mismatches++;
                            }
                        }
                        else if (sA && bB) {
                            for (l1 = start1, l2 = start2;
                                 l1 < end1 && l2 < end2;
                                 l1++, l2++) {
                                cA = sA[l1];
                                cB = (int) bB[l2-start2];
                                if (cA == wildcard || cB == wildcard) ;
                                else if (cA == cB) identities++;
                                else mismatches++;
                            }
                            Py_XDECREF(oB);
                        }
                        else if (sB && bA) {
                            for (l1 = start1, l2 = start2;
                                 l1 < end1 && l2 < end2;
                                 l1++, l2++) {
                                cA = (int) bA[l1-start1];
                                cB = sB[l2];
                                if (cA == wildcard || cB == wildcard) ;
                                else if (cA == cB) identities++;
                                else mismatches++;
                            }
                            Py_XDECREF(oA);
                        }
                        else if (bA && bB) {
                            for (l1 = start1, l2 = start2;
                                 l1 < end1 && l2 < end2;
                                 l1++, l2++) {
                                cA = (int) bA[l1-start1];
                                cB = (int) bB[l2-start2];
                                if (cA == wildcard || cB == wildcard) ;
                                else if (cA == cB) identities++;
                                else mismatches++;
                            }
                            Py_XDECREF(oA);
                            Py_XDECREF(oB);
                        }
                    }
                    else {
                        double* ptr;
                        double value;
                        if (sA && sB) {
                            for (l1 = start1, l2 = start2;
                                 l1 < end1 && l2 < end2;
                                 l1++, l2++) {
                                cA = sA[l1];
                                cB = sB[l2];
                                if (cA < 0) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequences[%d][%zd] is negative (%d)",
                                        i, l1, cA);
                                    goto error;
                                }
                                if (cA >= m) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequence[%d][%zd] is out of bound"
                                        " (%d, should be < %zd)",
                                        i, l1, cA, m);
                                    goto error;
                                }
                                if (cB < 0) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequences[%d][%zd] is negative (%d)",
                                        j, l2, cB);
                                    goto error;
                                }
                                if (cB >= m) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequence[%d][%zd] is out of bound"
                                        " (%d, should be < %zd)",
                                        j, l2, cB, m);
                                    goto error;
                                }
                                if (mapping) {
                                    cA = mapping[cA];
                                    cB = mapping[cB];
                                    if (cA == MISSING_LETTER || cB == MISSING_LETTER) {
                                        PyErr_SetString(PyExc_ValueError,
                                            "sequence contains letters not in the alphabet");
                                        goto error;
                                    }
                                }
                                if (cA == cB) identities++;
                                else mismatches++;
                                ptr = (double*)substitution_matrix.buf
                                    + cA * substitution_matrix.shape[0] + cB;
                                value = *(double*)ptr;
                                if (value > 0) positives++;
                                substitution_score += value;
                            }
                        }
                        else if (sA && bB) {
                            for (l1 = start1, l2 = start2;
                                 l1 < end1 && l2 < end2;
                                 l1++, l2++) {
                                cA = sA[l1];
                                cB = (int) bB[l2-start2];
                                if (cA < 0) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequences[%d][%zd] is negative (%d)",
                                        i, l1, cA);
                                    Py_DECREF(oB);
                                    goto error;
                                }
                                if (cA >= m) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequence[%d][%zd] is out of bound"
                                        " (%d, should be < %zd)",
                                        i, l1, cA, m);
                                    Py_DECREF(oB);
                                    goto error;
                                }
                                if (cB < 0) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequences[%d][%zd] is negative (%d)",
                                        j, l2, cB);
                                    Py_DECREF(oB);
                                    goto error;
                                }
                                if (cB >= m) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequence[%d][%zd] is out of bound"
                                        " (%d, should be < %zd)",
                                        j, l2, cB, m);
                                    Py_DECREF(oB);
                                    goto error;
                                }
                                if (mapping) {
                                    cA = mapping[cA];
                                    cB = mapping[cB];
                                    if (cA == MISSING_LETTER || cB == MISSING_LETTER) {
                                        PyErr_SetString(PyExc_ValueError,
                                            "sequence contains letters not in the alphabet");
                                        goto error;
                                    }
                                }
                                if (cA == wildcard || cB == wildcard) ;
                                else if (cA == cB) identities++;
                                else mismatches++;
                                ptr = (double*)substitution_matrix.buf
                                    + cA * substitution_matrix.shape[0] + cB;
                                value = *(double*)ptr;
                                if (value > 0) positives++;
                                substitution_score += value;
                            }
                            Py_XDECREF(oB);
                        }
                        else if (sB && bA) {
                            for (l1 = start1, l2 = start2;
                                 l1 < end1 && l2 < end2;
                                 l1++, l2++) {
                                cA = (int) bA[l1-start1];
                                cB = sB[l2];
                                if (cA < 0) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequences[%d][%zd] is negative (%d)",
                                        i, l1, cA);
                                    Py_DECREF(oA);
                                    goto error;
                                }
                                if (cB < 0) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequences[%d][%zd] is negative (%d)",
                                        j, l2, cB);
                                    Py_DECREF(oA);
                                    goto error;
                                }
                                if (cA >= m) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequence[%d][%zd] is out of bound"
                                        " (%d, should be < %zd)",
                                        i, l1, cA, m);
                                    Py_DECREF(oA);
                                    goto error;
                                }
                                if (cB >= m) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequence[%d][%zd] is out of bound"
                                        " (%d, should be < %zd)",
                                        j, l2, cB, m);
                                    Py_DECREF(oA);
                                    goto error;
                                }
                                if (mapping) {
                                    cA = mapping[cA];
                                    cB = mapping[cB];
                                    if (cA == MISSING_LETTER || cB == MISSING_LETTER) {
                                        PyErr_SetString(PyExc_ValueError,
                                            "sequence contains letters not in the alphabet");
                                        goto error;
                                    }
                                }
                                if (cA == wildcard || cB == wildcard) ;
                                else if (cA == cB) identities++;
                                else mismatches++;
                                ptr = (double*)substitution_matrix.buf
                                    + cA * substitution_matrix.shape[0] + cB;
                                value = *(double*)ptr;
                                if (value > 0) positives++;
                                substitution_score += value;
                            }
                            Py_XDECREF(oA);
                        }
                        else if (bA && bB) {
                            for (l1 = start1, l2 = start2;
                                 l1 < end1 && l2 < end2;
                                 l1++, l2++) {
                                cA = (int) bA[l1-start1];
                                cB = (int) bB[l2-start2];
                                if (cA < 0) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequences[%d][%zd] is negative (%d)",
                                        i, l1, cA);
                                    Py_DECREF(oA);
                                    Py_DECREF(oB);
                                    goto error;
                                }
                                if (cA >= m) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequence[%d][%zd] is out of bound"
                                        " (%d, should be < %zd)",
                                        i, l1, cA, m);
                                    Py_DECREF(oA);
                                    Py_DECREF(oB);
                                    goto error;
                                }
                                if (cB < 0) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequences[%d][%zd] is negative (%d)",
                                        j, l2, cB);
                                    Py_DECREF(oA);
                                    Py_DECREF(oB);
                                    goto error;
                                }
                                if (cB >= m) {
                                    PyErr_Format(PyExc_ValueError,
                                        "sequence[%d][%zd] is out of bound"
                                        " (%d, should be < %zd)",
                                        j, l2, cB, m);
                                    Py_DECREF(oA);
                                    Py_DECREF(oB);
                                    goto error;
                                }
                                if (mapping) {
                                    cA = mapping[cA];
                                    cB = mapping[cB];
                                    if (cA == MISSING_LETTER || cB == MISSING_LETTER) {
                                        PyErr_SetString(PyExc_ValueError,
                                            "sequence contains letters not in the alphabet");
                                        goto error;
                                    }
                                }
                                if (cA == wildcard || cB == wildcard) ;
                                else if (cA == cB) identities++;
                                else mismatches++;
                                ptr = (double*)substitution_matrix.buf
                                    + cA * substitution_matrix.shape[0] + cB;
                                value = *(double*)ptr;
                                if (value > 0) positives++;
                                substitution_score += value;
                            }
                            Py_XDECREF(oA);
                            Py_XDECREF(oB);
                        }
                    }
                }
                start1 = end1;
                start2 = end2;
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

    if (aligner && counts->identities + counts->mismatches > 0) {
        if (substitution_matrix.obj == NULL) {
            substitution_score = aligner->match * counts->identities
                               + aligner->mismatch * counts->mismatches;
        }
    } else substitution_score = Py_NAN;
    counts->substitution_score = substitution_score;
    if (aligner && aligner->insertion_score_function == NULL) {
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
    if (aligner && aligner->deletion_score_function == NULL) {
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
    counts->gap_score = gap_score;
    goto exit;

error:
    Py_XDECREF(oA);
    Py_XDECREF(oB);
    Py_DECREF(counts);
    counts = NULL;

exit:
    if (buffers) {
        for (i = 0; i < n; i++)
            if (buffers[i].buf) PyBuffer_Release(&buffers[i]);
    }
    coordinates_converter(NULL, &coordinates);
    strands_converter(NULL, &strands);
    substitution_matrix_converter(NULL, &substitution_matrix);

    return (PyObject*) counts;
}

static PyMethodDef module_functions[] = {
    {"calculate",
     (PyCFunction)_calculate,
     METH_VARARGS | METH_KEYWORDS,
     _calculate__doc__
    },
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_alignmentcounts",
        _alignmentcounts__doc__,
        -1,
        module_functions,
        NULL,
        NULL,
        NULL,
        NULL
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

    PyObject *capsule = PyObject_GetAttrString(mod, "mapping_buffer_capsule");
    Py_DECREF(mod);
    if (!capsule) {
        Py_DECREF(module);
        return NULL;
    }         
    Array_get_mapping_buffer =
        (Array_get_mapping_buffer_signature)PyCapsule_GetPointer(capsule,
        "Bio.Align.substitution_matrices._arraycore.mapping_buffer_capsule");
    if (!Array_get_mapping_buffer) {
        Py_DECREF(module);
        return NULL;
    }

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

    return module;
}
