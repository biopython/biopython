/* Copyright 2018-2019 by Michiel de Hoon.  All rights reserved.
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */



#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "float.h"


#define HORIZONTAL 0x1
#define VERTICAL 0x2
#define DIAGONAL 0x4
#define STARTPOINT 0x8
#define ENDPOINT 0x10
#define M_MATRIX 0x1
#define Ix_MATRIX 0x2
#define Iy_MATRIX 0x4
#define DONE 0x3
#define NONE 0x7

#define OVERFLOW_ERROR -1
#define MEMORY_ERROR -2

#define ANY_LETTER -1
#define MISSING_LETTER -2
#define UNMAPPED -3

#define SAFE_ADD(t, s) \
{   if (s != OVERFLOW_ERROR) { \
        term = t; \
        if (term > PY_SSIZE_T_MAX - s) s = OVERFLOW_ERROR; \
        else s += term; \
    } \
}


typedef enum {NeedlemanWunschSmithWaterman,
              Gotoh,
              WatermanSmithBeyer,
              Unknown} Algorithm;

typedef enum {Global, Local} Mode;

typedef struct {
    unsigned char trace : 5;
    unsigned char path : 3;
} Trace;

typedef struct {
    unsigned char Ix : 4;
    unsigned char Iy : 4;
} TraceGapsGotoh;

typedef struct {
    int* MIx;
    int* IyIx;
    int* MIy;
    int* IxIy;
} TraceGapsWatermanSmithBeyer;

static PyObject*
_create_path(Trace** M, int i, int j) {
    PyObject* tuple;
    PyObject* row;
    PyObject* value;
    int path;
    const int ii = i;
    const int jj = j;
    int n = 1;
    int direction = 0;

    while (1) {
        path = M[i][j].path;
        if (!path) break;
        if (path != direction) {
            n++;
            direction = path;
        }
        switch (path) {
            case HORIZONTAL: j++; break;
            case VERTICAL: i++; break;
            case DIAGONAL: i++; j++; break;
        }
    }
    i = ii;
    j = jj;

    direction = 0;
    tuple = PyTuple_New(n);
    if (!tuple) return NULL;
    n = 0;
    while (1) {
        path = M[i][j].path;
        if (path != direction) {
            row = PyTuple_New(2);
            if (!row) break;
            value = PyLong_FromLong(i);
            if (!value) {
                Py_DECREF(row); /* all references were stolen */
                break;
            }
            PyTuple_SET_ITEM(row, 0, value);
            value = PyLong_FromLong(j);
            if (!value) {
                Py_DECREF(row); /* all references were stolen */
                break;
            }
            PyTuple_SET_ITEM(row, 1, value);
            PyTuple_SET_ITEM(tuple, n, row);
            n++;
            direction = path;
        }
        switch (path) {
            case HORIZONTAL: j++; break;
            case VERTICAL: i++; break;
            case DIAGONAL: i++; j++; break;
            default: return tuple;
        }
    }
    Py_DECREF(tuple); /* all references were stolen */
    return PyErr_NoMemory();
}

typedef struct {
    PyObject_HEAD
    Trace** M;
    union { TraceGapsGotoh** gotoh;
            TraceGapsWatermanSmithBeyer** waterman_smith_beyer; } gaps;
    int nA;
    int nB;
    int iA;
    int iB;
    Mode mode;
    Algorithm algorithm;
    Py_ssize_t length;
} PathGenerator;

static Py_ssize_t
PathGenerator_needlemanwunsch_length(PathGenerator* self)
{
    int i;
    int j;
    int trace;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M;
    Py_ssize_t term;
    Py_ssize_t count = MEMORY_ERROR;
    Py_ssize_t temp;
    Py_ssize_t* counts;
    counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!counts) goto exit;
    counts[0] = 1;
    for (j = 1; j <= nB; j++) {
        trace = M[0][j].trace;
        count = 0;
        if (trace & HORIZONTAL) SAFE_ADD(counts[j-1], count);
        counts[j] = count;
    }
    for (i = 1; i <= nA; i++) {
        trace = M[i][0].trace;
        count = 0;
        if (trace & VERTICAL) SAFE_ADD(counts[0], count);
        temp = counts[0];
        counts[0] = count;
        for (j = 1; j <= nB; j++) {
            trace = M[i][j].trace;
            count = 0;
            if (trace & HORIZONTAL) SAFE_ADD(counts[j-1], count);
            if (trace & VERTICAL) SAFE_ADD(counts[j], count);
            if (trace & DIAGONAL) SAFE_ADD(temp, count);
            temp = counts[j];
            counts[j] = count;
        }
    }
    PyMem_Free(counts);
exit:
    return count;
}

static Py_ssize_t
PathGenerator_smithwaterman_length(PathGenerator* self)
{
    int i;
    int j;
    int trace;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M;
    Py_ssize_t term;
    Py_ssize_t count = MEMORY_ERROR;
    Py_ssize_t total = 0;
    Py_ssize_t temp;
    Py_ssize_t* counts;
    counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!counts) goto exit;
    counts[0] = 1;
    for (j = 1; j <= nB; j++) counts[j] = 1;
    for (i = 1; i <= nA; i++) {
        temp = counts[0];
        counts[0] = 1;
        for (j = 1; j <= nB; j++) {
            trace = M[i][j].trace;
            count = 0;
            if (trace & DIAGONAL) SAFE_ADD(temp, count);
            if (M[i][j].trace & ENDPOINT) SAFE_ADD(count, total);
            if (trace & HORIZONTAL) SAFE_ADD(counts[j-1], count);
            if (trace & VERTICAL) SAFE_ADD(counts[j], count);
            temp = counts[j];
            if (count == 0 && (trace & STARTPOINT)) count = 1;
            counts[j] = count;
        }
    }
    count = total;
    PyMem_Free(counts);
exit:
    return count;
}

static Py_ssize_t
PathGenerator_gotoh_global_length(PathGenerator* self)
{
    int i;
    int j;
    int trace;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M;
    TraceGapsGotoh** gaps = self->gaps.gotoh;
    Py_ssize_t count = MEMORY_ERROR;
    Py_ssize_t term;
    Py_ssize_t M_temp;
    Py_ssize_t Ix_temp;
    Py_ssize_t Iy_temp;
    Py_ssize_t* M_counts = NULL;
    Py_ssize_t* Ix_counts = NULL;
    Py_ssize_t* Iy_counts = NULL;
    M_counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!M_counts) goto exit;
    Ix_counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!Ix_counts) goto exit;
    Iy_counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!Iy_counts) goto exit;
    M_counts[0] = 1;
    Ix_counts[0] = 0;
    Iy_counts[0] = 0;
    for (j = 1; j <= nB; j++) {
        M_counts[j] = 0;
        Ix_counts[j] = 0;
        Iy_counts[j] = 1;
    }
    for (i = 1; i <= nA; i++) {
        M_temp = M_counts[0];
        M_counts[0] = 0;
        Ix_temp = Ix_counts[0];
        Ix_counts[0] = 1;
        Iy_temp = Iy_counts[0];
        Iy_counts[0] = 0;
        for (j = 1; j <= nB; j++) {
            count = 0;
            trace = M[i][j].trace;
            if (trace & M_MATRIX) SAFE_ADD(M_temp, count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_temp, count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_temp, count);
            M_temp = M_counts[j];
            M_counts[j] = count;
            count = 0;
            trace = gaps[i][j].Ix;
            if (trace & M_MATRIX) SAFE_ADD(M_temp, count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_counts[j], count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_counts[j], count);
            Ix_temp = Ix_counts[j];
            Ix_counts[j] = count;
            count = 0;
            trace = gaps[i][j].Iy;
            if (trace & M_MATRIX) SAFE_ADD(M_counts[j-1], count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_counts[j-1], count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_counts[j-1], count);
            Iy_temp = Iy_counts[j];
            Iy_counts[j] = count;
        }
    }
    count = 0;
    if (M[nA][nB].trace) SAFE_ADD(M_counts[nB], count);
    if (gaps[nA][nB].Ix) SAFE_ADD(Ix_counts[nB], count);
    if (gaps[nA][nB].Iy) SAFE_ADD(Iy_counts[nB], count);
exit:
    if (M_counts) PyMem_Free(M_counts);
    if (Ix_counts) PyMem_Free(Ix_counts);
    if (Iy_counts) PyMem_Free(Iy_counts);
    return count;
}

static Py_ssize_t
PathGenerator_gotoh_local_length(PathGenerator* self)
{
    int i;
    int j;
    int trace;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M;
    TraceGapsGotoh** gaps = self->gaps.gotoh;
    Py_ssize_t term;
    Py_ssize_t count = MEMORY_ERROR;
    Py_ssize_t total = 0;
    Py_ssize_t M_temp;
    Py_ssize_t Ix_temp;
    Py_ssize_t Iy_temp;
    Py_ssize_t* M_counts = NULL;
    Py_ssize_t* Ix_counts = NULL;
    Py_ssize_t* Iy_counts = NULL;
    M_counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!M_counts) goto exit;
    Ix_counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!Ix_counts) goto exit;
    Iy_counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!Iy_counts) goto exit;
    M_counts[0] = 1;
    Ix_counts[0] = 0;
    Iy_counts[0] = 0;
    for (j = 1; j <= nB; j++) {
        M_counts[j] = 1;
        Ix_counts[j] = 0;
        Iy_counts[j] = 0;
    }
    for (i = 1; i <= nA; i++) {
        M_temp = M_counts[0];
        M_counts[0] = 1;
        Ix_temp = Ix_counts[0];
        Ix_counts[0] = 0;
        Iy_temp = Iy_counts[0];
        Iy_counts[0] = 0;
        for (j = 1; j <= nB; j++) {
            count = 0;
            trace = M[i][j].trace;
            if (trace & M_MATRIX) SAFE_ADD(M_temp, count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_temp, count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_temp, count);
            if (count == 0 && (trace & STARTPOINT)) count = 1;
            M_temp = M_counts[j];
            M_counts[j] = count;
            if (M[i][j].trace & ENDPOINT) SAFE_ADD(count, total);
            count = 0;
            trace = gaps[i][j].Ix;
            if (trace & M_MATRIX) SAFE_ADD(M_temp, count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_counts[j], count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_counts[j], count);
            Ix_temp = Ix_counts[j];
            Ix_counts[j] = count;
            count = 0;
            trace = gaps[i][j].Iy;
            if (trace & M_MATRIX) SAFE_ADD(M_counts[j-1], count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_counts[j-1], count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_counts[j-1], count);
            Iy_temp = Iy_counts[j];
            Iy_counts[j] = count;
        }
    }
    count = total;
exit:
    if (M_counts) PyMem_Free(M_counts);
    if (Ix_counts) PyMem_Free(Ix_counts);
    if (Iy_counts) PyMem_Free(Iy_counts);
    return count;
}

static Py_ssize_t
PathGenerator_waterman_smith_beyer_global_length(PathGenerator* self)
{
    int i;
    int j;
    int trace;
    int* p;
    int gap;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M;
    TraceGapsWatermanSmithBeyer** gaps = self->gaps.waterman_smith_beyer;
    Py_ssize_t count = MEMORY_ERROR;
    Py_ssize_t term;
    Py_ssize_t** M_count = NULL;
    Py_ssize_t** Ix_count = NULL;
    Py_ssize_t** Iy_count = NULL;
    M_count = PyMem_Malloc((nA+1)*sizeof(Py_ssize_t*));
    if (!M_count) goto exit;
    Ix_count = PyMem_Malloc((nA+1)*sizeof(Py_ssize_t*));
    if (!Ix_count) goto exit;
    Iy_count = PyMem_Malloc((nA+1)*sizeof(Py_ssize_t*));
    if (!Iy_count) goto exit;
    for (i = 0; i <= nA; i++) {
        M_count[i] = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
        if (!M_count[i]) goto exit;
        Ix_count[i] = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
        if (!Ix_count[i]) goto exit;
        Iy_count[i] = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
        if (!Iy_count[i]) goto exit;
    }
    for (i = 0; i <= nA; i++) {
        for (j = 0; j <= nB; j++) {
            count = 0;
            trace = M[i][j].trace;
            if (trace & M_MATRIX) SAFE_ADD(M_count[i-1][j-1], count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_count[i-1][j-1], count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_count[i-1][j-1], count);
            if (count == 0) count = 1; /* happens at M[0][0] only */
            M_count[i][j] = count;
            count = 0;
            p = gaps[i][j].MIx;
            if (p) {
                while (1) {
                    gap = *p;
                    if (!gap) break;
                    SAFE_ADD(M_count[i-gap][j], count);
                    p++;
                }
            }
            p = gaps[i][j].IyIx;
            if (p) {
                while (1) {
                    gap = *p;
                    if (!gap) break;
                    SAFE_ADD(Iy_count[i-gap][j], count);
                    p++;
                }
            }
            Ix_count[i][j] = count;
            count = 0;
            p = gaps[i][j].MIy;
            if (p) {
                while (1) {
                    gap = *p;
                    if (!gap) break;
                    SAFE_ADD(M_count[i][j-gap], count);
                    p++;
                }
            }
	    p = gaps[i][j].IxIy;
            if (p) {
                while (1) {
                    gap = *p;
                    if (!gap) break;
                    SAFE_ADD(Ix_count[i][j-gap], count);
                    p++;
                }
            }
            Iy_count[i][j] = count;
        }
    }
    count = 0;
    if (M[nA][nB].trace)
        SAFE_ADD(M_count[nA][nB], count);
    if (gaps[nA][nB].MIx[0] || gaps[nA][nB].IyIx[0])
        SAFE_ADD(Ix_count[nA][nB], count);
    if (gaps[nA][nB].MIy[0] || gaps[nA][nB].IxIy[0])
        SAFE_ADD(Iy_count[nA][nB], count);
exit:
    if (M_count) {
        if (Ix_count) {
            if (Iy_count) {
                for (i = 0; i <= nA; i++) {
                    if (!M_count[i]) break;
                    PyMem_Free(M_count[i]);
                    if (!Ix_count[i]) break;
                    PyMem_Free(Ix_count[i]);
                    if (!Iy_count[i]) break;
                    PyMem_Free(Iy_count[i]);
                }
                PyMem_Free(Iy_count);
            }
            PyMem_Free(Ix_count);
        }
        PyMem_Free(M_count);
    }
    return count;
}

static Py_ssize_t
PathGenerator_waterman_smith_beyer_local_length(PathGenerator* self)
{
    int i;
    int j;
    int trace;
    int* p;
    int gap;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M;
    TraceGapsWatermanSmithBeyer** gaps = self->gaps.waterman_smith_beyer;
    Py_ssize_t term;
    Py_ssize_t count = MEMORY_ERROR;
    Py_ssize_t total = 0;
    Py_ssize_t** M_count = NULL;
    Py_ssize_t** Ix_count = NULL;
    Py_ssize_t** Iy_count = NULL;
    M_count = PyMem_Malloc((nA+1)*sizeof(Py_ssize_t*));
    if (!M_count) goto exit;
    Ix_count = PyMem_Malloc((nA+1)*sizeof(Py_ssize_t*));
    if (!Ix_count) goto exit;
    Iy_count = PyMem_Malloc((nA+1)*sizeof(Py_ssize_t*));
    if (!Iy_count) goto exit;
    for (i = 0; i <= nA; i++) {
        M_count[i] = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
        if (!M_count[i]) goto exit;
        Ix_count[i] = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
        if (!Ix_count[i]) goto exit;
        Iy_count[i] = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
        if (!Iy_count[i]) goto exit;
    }
    for (i = 0; i <= nA; i++) {
        for (j = 0; j <= nB; j++) {
            count = 0;
            trace = M[i][j].trace;
            if (trace & M_MATRIX) SAFE_ADD(M_count[i-1][j-1], count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_count[i-1][j-1], count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_count[i-1][j-1], count);
            if (count == 0 && (trace & STARTPOINT)) count = 1;
            M_count[i][j] = count;
            if (M[i][j].trace & ENDPOINT) SAFE_ADD(count, total);
            count = 0;
            p = gaps[i][j].MIx;
            if (p) {
                while (1) {
                    gap = *p;
                    if (!gap) break;
                    SAFE_ADD(M_count[i-gap][j], count);
                    p++;
                }
            }
            p = gaps[i][j].IyIx;
            if (p) {
                while (1) {
                    gap = *p;
                    if (!gap) break;
                    SAFE_ADD(Iy_count[i-gap][j], count);
                    p++;
                }
            }
            Ix_count[i][j] = count;
            count = 0;
            p = gaps[i][j].MIy;
            if (p) {
                while (1) {
                    gap = *p;
                    if (!gap) break;
                    SAFE_ADD(M_count[i][j-gap], count);
                    p++;
                }
            }
            p = gaps[i][j].IxIy;
            if (p) {
                while (1) {
                    gap = *p;
                    if (!gap) break;
                    SAFE_ADD(Ix_count[i][j-gap], count);
                    p++;
                }
            }
            Iy_count[i][j] = count;
        }
    }
    count = total;
exit:
    if (M_count) {
        if (Ix_count) {
            if (Iy_count) {
                for (i = 0; i <= nA; i++) {
                    if (!M_count[i]) break;
                    PyMem_Free(M_count[i]);
                    if (!Ix_count[i]) break;
                    PyMem_Free(Ix_count[i]);
                    if (!Iy_count[i]) break;
                    PyMem_Free(Iy_count[i]);
                }
                PyMem_Free(Iy_count);
            }
            PyMem_Free(Ix_count);
        }
        PyMem_Free(M_count);
    }
    return count;
}

static Py_ssize_t PathGenerator_length(PathGenerator* self) {
    Py_ssize_t length = self->length;
    if (length == 0) {
        switch (self->algorithm) {
            case NeedlemanWunschSmithWaterman:
                switch (self->mode) {
                    case Global:
                        length = PathGenerator_needlemanwunsch_length(self);
                        break;
                    case Local:
                        length = PathGenerator_smithwaterman_length(self);
                        break;
                    default:
                        /* should not happen, but some compilers complain that
                         * that length can be used uninitialized.
                         */
                        PyErr_SetString(PyExc_RuntimeError, "Unknown mode");
                        return -1;
                }
                break;
            case Gotoh:
                switch (self->mode) {
                    case Global:
                        length = PathGenerator_gotoh_global_length(self);
                        break;
                    case Local:
                        length = PathGenerator_gotoh_local_length(self);
                        break;
                    default:
                        /* should not happen, but some compilers complain that
                         * that length can be used uninitialized.
                         */
                        PyErr_SetString(PyExc_RuntimeError, "Unknown mode");
                        return -1;
                }
                break;
            case WatermanSmithBeyer:
                switch (self->mode) {
                    case Global:
                        length = PathGenerator_waterman_smith_beyer_global_length(self);
                        break;
                    case Local:
                        length = PathGenerator_waterman_smith_beyer_local_length(self);
                        break;
                    default:
                        /* should not happen, but some compilers complain that
                         * that length can be used uninitialized.
                         */
                        PyErr_SetString(PyExc_RuntimeError, "Unknown mode");
                        return -1;
                }
                break;
            case Unknown:
            default:
                PyErr_SetString(PyExc_RuntimeError, "Unknown algorithm");
                return -1;
        }
        self->length = length;
    }
    switch (length) {
        case OVERFLOW_ERROR:
            PyErr_Format(PyExc_OverflowError,
                         "number of optimal alignments is larger than %zd",
                         PY_SSIZE_T_MAX);
            break;
        case MEMORY_ERROR:
            PyErr_SetNone(PyExc_MemoryError);
            break;
        default:
            break;
    }
    return length;
}

static void
PathGenerator_dealloc(PathGenerator* self)
{
    int i;
    const int nA = self->nA;
    const Algorithm algorithm = self->algorithm;
    Trace** M = self->M;
    if (M) {
        for (i = 0; i <= nA; i++) {
            if (!M[i]) break;
            PyMem_Free(M[i]);
        }
        PyMem_Free(M);
    }
    switch (algorithm) {
        case NeedlemanWunschSmithWaterman:
            break;
        case Gotoh: {
            TraceGapsGotoh** gaps = self->gaps.gotoh;
            if (gaps) {
                for (i = 0; i <= nA; i++) {
                    if (!gaps[i]) break;
                    PyMem_Free(gaps[i]);
                }
                PyMem_Free(gaps);
            }
            break;
        }
        case WatermanSmithBeyer: {
            TraceGapsWatermanSmithBeyer** gaps = self->gaps.waterman_smith_beyer;
            if (gaps) {
                int j;
                const int nB = self->nB;
                int* trace;
                for (i = 0; i <= nA; i++) {
                    if (!gaps[i]) break;
                    for (j = 0; j <= nB; j++) {
                        trace = gaps[i][j].MIx;
                        if (trace) PyMem_Free(trace);
                        trace = gaps[i][j].IyIx;
                        if (trace) PyMem_Free(trace);
                        trace = gaps[i][j].MIy;
                        if (trace) PyMem_Free(trace);
                        trace = gaps[i][j].IxIy;
                        if (trace) PyMem_Free(trace);
                    }
                    PyMem_Free(gaps[i]);
                }
                PyMem_Free(gaps);
            }
            break;
        }
        case Unknown:
        default:
            PyErr_WriteUnraisable((PyObject*)self);
            break;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* PathGenerator_next_needlemanwunsch(PathGenerator* self)
{
    int i = 0;
    int j = 0;
    int path;
    int trace = 0;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M;

    path = M[i][j].path;
    if (path == DONE) return NULL;
    if (path == 0) {
        /* Generate the first path. */
        i = nA;
        j = nB;
    }
    else {
        /* We already have a path. Prune the path to see if there are
         * any alternative paths. */
        while (1) {
            if (path == HORIZONTAL) {
                trace = M[i][++j].trace;
                if (trace & VERTICAL) {
                    M[--i][j].path = VERTICAL;
                    break;
                }
                if (trace & DIAGONAL) {
                    M[--i][--j].path = DIAGONAL;
                    break;
                }
            }
            else if (path == VERTICAL) {
                trace = M[++i][j].trace;
                if (trace & DIAGONAL) {
                    M[--i][--j].path = DIAGONAL;
                    break;
                }
            }
            else /* DIAGONAL */ {
                i++;
                j++;
            }
            path = M[i][j].path;
            if (!path) {
                /* we reached the end of the alignment without finding
                 * an alternative path */
                M[0][0].path = DONE;
                return NULL;
            }
        }
    }
    /* Follow the traceback until we reach the origin. */
    while (1) {
        trace = M[i][j].trace;
        if (trace & HORIZONTAL) M[i][--j].path = HORIZONTAL;
        else if (trace & VERTICAL) M[--i][j].path = VERTICAL;
        else if (trace & DIAGONAL) M[--i][--j].path = DIAGONAL;
        else break;
    }
    return _create_path(M, 0, 0);
}

static PyObject* PathGenerator_next_smithwaterman(PathGenerator* self)
{
    int trace = 0;
    int i = self->iA;
    int j = self->iB;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M;
    int path = M[0][0].path;

    if (path == DONE || path == NONE) return NULL;

    path = M[i][j].path;
    if (path) {
        /* We already have a path. Prune the path to see if there are
         * any alternative paths. */
        while (1) {
            if (path == HORIZONTAL) {
                trace = M[i][++j].trace;
                if (trace & VERTICAL) {
                    M[--i][j].path = VERTICAL;
                    break;
                }
                else if (trace & DIAGONAL) {
                    M[--i][--j].path = DIAGONAL;
                    break;
                }
            }
            else if (path == VERTICAL) {
                trace = M[++i][j].trace;
                if (trace & DIAGONAL) {
                    M[--i][--j].path = DIAGONAL;
                    break;
                }
            }
            else /* DIAGONAL */ {
                i++;
                j++;
            }
            path = M[i][j].path;
            if (!path) break;
        }
    }

    if (path) {
        trace = M[i][j].trace;
    } else {
        /* Find a suitable end point for a path.
         * Only allow end points ending at the M matrix. */
        while (1) {
            if (j < nB) j++;
            else if (i < nA) {
                i++;
                j = 0;
            }
            else {
                /* we reached the end of the sequences without finding
                 * an alternative path */
                M[0][0].path = DONE;
                return NULL;
            }
            trace = M[i][j].trace;
            if (trace & ENDPOINT) {
                trace &= DIAGONAL; /* exclude paths ending in a gap */
                break;
            }
        }
        M[i][j].path = 0;
    }

    /* Follow the traceback until we reach the origin. */
    while (1) {
        if (trace & HORIZONTAL) M[i][--j].path = HORIZONTAL;
        else if (trace & VERTICAL) M[--i][j].path = VERTICAL;
        else if (trace & DIAGONAL) M[--i][--j].path = DIAGONAL;
        else if (trace & STARTPOINT) {
            self->iA = i;
            self->iB = j;
            return _create_path(M, i, j);
        }
        else {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected trace in PathGenerator_next_smithwaterman");
            return NULL;
        }
        trace = M[i][j].trace;
    }
}

static PyObject* PathGenerator_next_gotoh_global(PathGenerator* self)
{
    int i = 0;
    int j = 0;
    int m;
    int path;
    int trace = 0;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M;
    TraceGapsGotoh** gaps = self->gaps.gotoh;

    m = M_MATRIX;
    path = M[i][j].path;
    if (path == DONE) return NULL;
    if (path == 0) {
        i = nA;
        j = nB;
    }
    else {
        /* We already have a path. Prune the path to see if there are
         * any alternative paths. */
        while (1) {
            path = M[i][j].path;
            if (path == 0) {
                switch (m) {
                    case M_MATRIX: m = Ix_MATRIX; break;
                    case Ix_MATRIX: m = Iy_MATRIX; break;
                    case Iy_MATRIX: m = 0; break;
                }
                break;
            }
            switch (path) {
                case HORIZONTAL: trace = gaps[i][++j].Iy; break;
                case VERTICAL: trace = gaps[++i][j].Ix; break;
                case DIAGONAL: trace = M[++i][++j].trace; break;
            }
            switch (m) {
                case M_MATRIX:
                    if (trace & Ix_MATRIX) {
                        m = Ix_MATRIX;
                        break;
                    }
                case Ix_MATRIX:
                    if (trace & Iy_MATRIX) {
                        m = Iy_MATRIX;
                        break;
                    }
                case Iy_MATRIX:
                default:
                    switch (path) {
                        case HORIZONTAL: m = Iy_MATRIX; break;
                        case VERTICAL: m = Ix_MATRIX; break;
                        case DIAGONAL: m = M_MATRIX; break;
                    }
                    continue;
            }
            switch (path) {
                case HORIZONTAL: j--; break;
                case VERTICAL: i--; break;
                case DIAGONAL: i--; j--; break;
            }
            M[i][j].path = path;
            break;
        }
    }

    if (path == 0) {
        /* Generate a new path. */
        switch (m) {
            case M_MATRIX:
                if (M[nA][nB].trace) {
                   /* m = M_MATRIX; */
                   break;
                }
            case Ix_MATRIX:
                if (gaps[nA][nB].Ix) {
                   m = Ix_MATRIX;
                   break;
                }
            case Iy_MATRIX:
                if (gaps[nA][nB].Iy) {
                   m = Iy_MATRIX;
                   break;
                }
            default:
                /* exhausted this generator */
                M[0][0].path = DONE;
                return NULL;
        }
    }

    switch (m) {
        case M_MATRIX:
            trace = M[i][j].trace;
            path = DIAGONAL;
            i--; j--;
            break;
        case Ix_MATRIX:
            trace = gaps[i][j].Ix;
            path = VERTICAL;
            i--;
            break;
        case Iy_MATRIX:
            trace = gaps[i][j].Iy;
            path = HORIZONTAL;
            j--;
            break;
    }

    while (1) {
        if (trace & M_MATRIX) {
            trace = M[i][j].trace;
            M[i][j].path = path;
            path = DIAGONAL;
            i--; j--;
        }
        else if (trace & Ix_MATRIX) {
            M[i][j].path = path;
            trace = gaps[i][j].Ix;
            path = VERTICAL;
            i--;
        }
        else if (trace & Iy_MATRIX) {
            M[i][j].path = path;
            trace = gaps[i][j].Iy;
            path = HORIZONTAL;
            j--;
        }
        else break;
    }
    return _create_path(M, 0, 0);
}

static PyObject* PathGenerator_next_gotoh_local(PathGenerator* self)
{
    int trace = 0;
    int i;
    int j;
    int m = M_MATRIX;
    int iA = self->iA;
    int iB = self->iB;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M;
    TraceGapsGotoh** gaps = self->gaps.gotoh;
    int path = M[0][0].path;

    if (path == DONE) return NULL;

    path = M[iA][iB].path;

    if (path) {
        i = iA;
        j = iB;
        while (1) {
            /* We already have a path. Prune the path to see if there are
             * any alternative paths. */
            path = M[i][j].path;
            if (path == 0) {
                m = M_MATRIX;
                iA = i;
                iB = j;
                break;
            }
            switch (path) {
                case HORIZONTAL: trace = gaps[i][++j].Iy; break;
                case VERTICAL: trace = gaps[++i][j].Ix; break;
                case DIAGONAL: trace = M[++i][++j].trace; break;
            }
            switch (m) {
                case M_MATRIX:
                    if (trace & Ix_MATRIX) {
                        m = Ix_MATRIX;
                        break;
                    }
                case Ix_MATRIX:
                    if (trace & Iy_MATRIX) {
                        m = Iy_MATRIX;
                        break;
                    }
                case Iy_MATRIX:
                default:
                    switch (path) {
                        case HORIZONTAL: m = Iy_MATRIX; break;
                        case VERTICAL: m = Ix_MATRIX; break;
                        case DIAGONAL: m = M_MATRIX; break;
                    }
                    continue;
            }
            switch (path) {
                case HORIZONTAL: j--; break;
                case VERTICAL: i--; break;
                case DIAGONAL: i--; j--; break;
            }
            M[i][j].path = path;
            break;
        }
    }

    if (path == 0) {
        /* Find the end point for a new path. */
        while (1) {
            if (iB < nB) iB++;
            else if (iA < nA) {
                iA++;
                iB = 0;
            }
            else {
                /* we reached the end of the alignment without finding
                 * an alternative path */
                M[0][0].path = DONE;
                return NULL;
            }
            if (M[iA][iB].trace & ENDPOINT) {
                M[iA][iB].path = 0;
                break;
            }
        }
        m = M_MATRIX;
        i = iA;
        j = iB;
    }

    while (1) {
        switch (m) {
            case M_MATRIX: trace = M[i][j].trace; break;
            case Ix_MATRIX: trace = gaps[i][j].Ix; break;
            case Iy_MATRIX: trace = gaps[i][j].Iy; break;
        }
        if (trace == STARTPOINT) {
            self->iA = i;
            self->iB = j;
            return _create_path(M, i, j);
        }
        switch (m) {
            case M_MATRIX:
                path = DIAGONAL;
                i--;
                j--;
                break;
            case Ix_MATRIX:
                path = VERTICAL;
                i--;
                break;
            case Iy_MATRIX:
                path = HORIZONTAL;
                j--;
                break;
        }
        if (trace & M_MATRIX) m = M_MATRIX;
        else if (trace & Ix_MATRIX) m = Ix_MATRIX;
        else if (trace & Iy_MATRIX) m = Iy_MATRIX;
        else {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected trace in PathGenerator_next_gotoh_local");
            return NULL;
        }
        M[i][j].path = path;
    }
    return NULL;
}

static PyObject*
PathGenerator_next_waterman_smith_beyer_global(PathGenerator* self)
{
    int i = 0, j = 0;
    int iA, iB;
    int trace;
    int* gapM;
    int* gapXY;

    int m = M_MATRIX;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M;
    TraceGapsWatermanSmithBeyer** gaps = self->gaps.waterman_smith_beyer;

    int gap;
    int path = M[0][0].path;

    if (path == DONE) return NULL;

    if (path) {
        /* We already have a path. Prune the path to see if there are
         * any alternative paths. */
        while (1) {
            if (!path) {
                m <<= 1;
                break;
            }
            switch (path) {
                case HORIZONTAL:
                    iA = i;
                    iB = j;
                    while (M[i][iB].path == HORIZONTAL) iB++;
                    break;
                case VERTICAL:
                    iA = i;
                    while (M[iA][j].path == VERTICAL) iA++;
                    iB = j;
                    break;
                case DIAGONAL:
                    iA = i + 1;
                    iB = j + 1;
                    break;
                default:
                    PyErr_SetString(PyExc_RuntimeError,
                        "Unexpected path in PathGenerator_next_waterman_smith_beyer_global");
                    return NULL;
            }
            if (i == iA) { /* HORIZONTAL */
                gapM = gaps[iA][iB].MIy;
                gapXY = gaps[iA][iB].IxIy;
                if (m == M_MATRIX) {
                    gap = iB - j;
                    while (*gapM != gap) gapM++;
                    gapM++;
                    gap = *gapM;
                    if (gap) {
                        j = iB - gap;
                        while (j < iB) M[i][--iB].path = HORIZONTAL;
                        break;
                    }
                } else if (m == Ix_MATRIX) {
                    gap = iB - j;
                    while (*gapXY != gap) gapXY++;
                    gapXY++;
                }
                gap = *gapXY;
                if (gap) {
                    m = Ix_MATRIX;
                    j = iB - gap;
                    while (j < iB) M[i][--iB].path = HORIZONTAL;
                    break;
                }
                /* no alternative found; continue pruning */
                m = Iy_MATRIX;
                j = iB;
            }
            else if (j == iB) { /* VERTICAL */
                gapM = gaps[iA][iB].MIx;
                gapXY = gaps[iA][iB].IyIx;
                if (m == M_MATRIX) {
                    gap = iA - i;
                    while (*gapM != gap) gapM++;
                    gapM++;
                    gap = *gapM;
                    if (gap) {
                        i = iA - gap;
                        while (i < iA) M[--iA][j].path = VERTICAL;
                        break;
                    }
                } else if (m == Iy_MATRIX) {
                    gap = iA - i;
                    while (*gapXY != gap) gapXY++;
                    gapXY++;
                }
                gap = *gapXY;
                if (gap) {
                    m = Iy_MATRIX;
                    i = iA - gap;
                    while (i < iA) M[--iA][j].path = VERTICAL;
                    break;
                }
                /* no alternative found; continue pruning */
                m = Ix_MATRIX;
                i = iA;
            }
            else { /* DIAGONAL */
                i = iA - 1;
                j = iB - 1;
                trace = M[iA][iB].trace;
                switch (m) {
                    case M_MATRIX:
                        if (trace & Ix_MATRIX) {
                            m = Ix_MATRIX;
                            M[i][j].path = DIAGONAL;
                            break;
                        }
                    case Ix_MATRIX:
                        if (trace & Iy_MATRIX) {
                            m = Iy_MATRIX;
                            M[i][j].path = DIAGONAL;
                            break;
                        }
                    case Iy_MATRIX:
                    default:
                        /* no alternative found; continue pruning */
                        m = M_MATRIX;
                        i = iA;
                        j = iB;
                        path = M[i][j].path;
                        continue;
                }
                /* alternative found; build path until starting point */
                break;
            }
            path = M[i][j].path;
        }
    }

    if (!path) {
        /* Find a suitable end point for a path. */
        switch (m) {
            case M_MATRIX:
                if (M[nA][nB].trace) {
                    /* m = M_MATRIX; */
                    break;
                }
            case Ix_MATRIX:
                if (gaps[nA][nB].MIx[0] || gaps[nA][nB].IyIx[0]) {
                    m = Ix_MATRIX;
                    break;
                }
            case Iy_MATRIX:
                if (gaps[nA][nB].MIy[0] || gaps[nA][nB].IxIy[0]) {
                    m = Iy_MATRIX;
                    break;
                }
            default:
                M[0][0].path = DONE;
                return NULL;
        }
        i = nA;
        j = nB;
    }

    /* Follow the traceback until we reach the origin. */
    while (1) {
        switch (m) {
            case M_MATRIX:
                trace = M[i][j].trace;
                if (trace & M_MATRIX) m = M_MATRIX;
                else if (trace & Ix_MATRIX) m = Ix_MATRIX;
                else if (trace & Iy_MATRIX) m = Iy_MATRIX;
                else return _create_path(M, i, j);
                i--;
                j--;
                M[i][j].path = DIAGONAL;
                break;
            case Ix_MATRIX:
                gap = gaps[i][j].MIx[0];
                if (gap) m = M_MATRIX;
                else {
                    gap = gaps[i][j].IyIx[0];
                    m = Iy_MATRIX;
                }
                iA = i - gap;
                while (iA < i) M[--i][j].path = VERTICAL;
                M[i][j].path = VERTICAL;
                break;
            case Iy_MATRIX:
                gap = gaps[i][j].MIy[0];
                if (gap) m = M_MATRIX;
                else {
                    gap = gaps[i][j].IxIy[0];
                    m = Ix_MATRIX;
                }
                iB = j - gap;
                while (iB < j) M[i][--j].path = HORIZONTAL;
                M[i][j].path = HORIZONTAL;
                break;
        }
    }
}

static PyObject*
PathGenerator_next_waterman_smith_beyer_local(PathGenerator* self)
{
    int i, j, m;
    int trace = 0;
    int* gapM;
    int* gapXY;

    int iA = self->iA;
    int iB = self->iB;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M;
    TraceGapsWatermanSmithBeyer** gaps = self->gaps.waterman_smith_beyer;

    int gap;
    int path = M[0][0].path;

    if (path == DONE) return NULL;
    m = 0;
    path = M[iA][iB].path;
    if (path) {
        /* We already have a path. Prune the path to see if there are
         * any alternative paths. */
        m = M_MATRIX;
        i = iA;
        j = iB;
        while (1) {
            path = M[i][j].path;
            switch (path) {
                case HORIZONTAL:
                    iA = i;
                    iB = j;
                    while (M[i][iB].path == HORIZONTAL) iB++;
                    break;
                case VERTICAL:
                    iA = i;
                    iB = j;
                    while (M[iA][j].path == VERTICAL) iA++;
                    break;
                case DIAGONAL:
                    iA = i + 1;
                    iB = j + 1;
                    break;
                default:
                    iA = -1;
                    break;
            }
            if (iA < 0) {
                m = 0;
                iA = i;
                iB = j;
                break;
            }
            if (i == iA) { /* HORIZONTAL */
                gapM = gaps[iA][iB].MIy;
                gapXY = gaps[iA][iB].IxIy;
                if (m == M_MATRIX) {
                    gap = iB - j;
                    while (*gapM != gap) gapM++;
                    gapM++;
                    gap = *gapM;
                    if (gap) {
                        j = iB - gap;
                        while (j < iB) M[i][--iB].path = HORIZONTAL;
                        break;
                    }
                } else if (m == Ix_MATRIX) {
                    gap = iB - j;
                    while (*gapXY != gap) gapXY++;
                    gapXY++;
                }
                gap = *gapXY;
                if (gap) {
                    m = Ix_MATRIX;
                    j = iB - gap;
                    M[i][j].path = HORIZONTAL;
                    while (iB > j) M[i][--iB].path = HORIZONTAL;
                    break;
                }
                /* no alternative found; continue pruning */
                m = Iy_MATRIX;
                j = iB;
            }
            else if (j == iB) { /* VERTICAL */
                gapM = gaps[iA][iB].MIx;
                gapXY = gaps[iA][iB].IyIx;
                if (m == M_MATRIX) {
                    gap = iA - i;
                    while (*gapM != gap) gapM++;
                    gapM++;
                    gap = *gapM;
                    if (gap) {
                        i = iA - gap;
                        while (i < iA) M[--iA][j].path = VERTICAL;
                        break;
                    }
                } else if (m == Iy_MATRIX) {
                    gap = iA - i;
                    while (*gapXY != gap) gapXY++;
                    gapXY++;
                }
                gap = *gapXY;
                if (gap) {
                    m = Iy_MATRIX;
                    i = iA - gap;
                    M[i][j].path = VERTICAL;
                    while (iA > i) M[--iA][j].path = VERTICAL;
                    break;
                }
                /* no alternative found; continue pruning */
                m = Ix_MATRIX;
                i = iA;
            }
            else { /* DIAGONAL */
                i = iA - 1;
                j = iB - 1;
                trace = M[iA][iB].trace;
                switch (m) {
                    case M_MATRIX:
                        if (trace & Ix_MATRIX) {
                            m = Ix_MATRIX;
                            M[i][j].path = DIAGONAL;
                            break;
                        }
                    case Ix_MATRIX:
                        if (trace & Iy_MATRIX) {
                            m = Iy_MATRIX;
                            M[i][j].path = DIAGONAL;
                            break;
                        }
                    case Iy_MATRIX:
                    default:
                        /* no alternative found; continue pruning */
                        m = M_MATRIX;
                        i = iA;
                        j = iB;
                        continue;
                }
                /* alternative found; build path until starting point */
                break;
            }
        }
    }
 
    if (m == 0) {
        /* We are at [nA][nB]. Find a suitable end point for a path. */
        while (1) {
            if (iB < nB) iB++;
            else if (iA < nA) {
                iA++;
                iB = 0;
            }
            else {
                /* exhausted this generator */
                M[0][0].path = DONE;
                return NULL;
            }
            if (M[iA][iB].trace & ENDPOINT) break;
        }
        M[iA][iB].path = 0;
        m = M_MATRIX;
        i = iA;
        j = iB;
    }

    /* Follow the traceback until we reach the origin. */
    while (1) {
        switch (m) {
            case Ix_MATRIX:
                gapM = gaps[i][j].MIx;
                gapXY = gaps[i][j].IyIx;
                iB = j;
                gap = *gapM;
                if (gap) m = M_MATRIX;
                else {
                    gap = *gapXY;
                    m = Iy_MATRIX;
                }
                iA = i - gap;
                while (i > iA) M[--i][iB].path = VERTICAL;
                break;
            case Iy_MATRIX:
                gapM = gaps[i][j].MIy;
                gapXY = gaps[i][j].IxIy;
                iA = i;
                gap = *gapM;
                if (gap) m = M_MATRIX;
                else {
                    gap = *gapXY;
                    m = Ix_MATRIX;
                }
                iB = j - gap;
                while (j > iB) M[iA][--j].path = HORIZONTAL;
                break;
            case M_MATRIX:
                iA = i-1;
                iB = j-1;
                trace = M[i][j].trace;
                if (trace & M_MATRIX) m = M_MATRIX;
                else if (trace & Ix_MATRIX) m = Ix_MATRIX;
                else if (trace & Iy_MATRIX) m = Iy_MATRIX;
                else if (trace == STARTPOINT) {
                    self->iA = i;
                    self->iB = j;
                    return _create_path(M, i, j);
                }
                else {
                    PyErr_SetString(PyExc_RuntimeError,
                        "Unexpected trace in PathGenerator_next_waterman_smith_beyer_local");
                    return NULL;
                }
                M[iA][iB].path = DIAGONAL;
                break;
        }
        i = iA;
        j = iB;
    }
}

static PyObject *
PathGenerator_next(PathGenerator* self)
{
    const Mode mode = self->mode;
    const Algorithm algorithm = self->algorithm;
    switch (algorithm) {
        case NeedlemanWunschSmithWaterman:
            switch (mode) {
                case Global:
                    return PathGenerator_next_needlemanwunsch(self);
                case Local:
                    return PathGenerator_next_smithwaterman(self);
            }
        case Gotoh:
            switch (mode) {
                case Global:
                    return PathGenerator_next_gotoh_global(self);
                case Local:
                    return PathGenerator_next_gotoh_local(self);
            }
        case WatermanSmithBeyer:
            switch (mode) {
                case Global:
                    return PathGenerator_next_waterman_smith_beyer_global(self);
                case Local:
                    return PathGenerator_next_waterman_smith_beyer_local(self);
            }
        case Unknown:
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unknown algorithm");
            return NULL;
    }
}

static const char PathGenerator_reset__doc__[] = "reset the iterator";

static PyObject*
PathGenerator_reset(PathGenerator* self)
{
    switch (self->mode) {
        case Local:
            self->iA = 0;
            self->iB = 0;
        case Global: {
            Trace** M = self->M;
            switch (self->algorithm) {
                case NeedlemanWunschSmithWaterman:
                case Gotoh: {
                    if (M[0][0].path != NONE) M[0][0].path = 0;
                    break;
                }
                case WatermanSmithBeyer: {
                    M[0][0].path = 0;
                    break;
                }
                case Unknown:
                default:
                    break;
            }
        }
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef PathGenerator_methods[] = {
    {"reset",
     (PyCFunction)PathGenerator_reset,
     METH_NOARGS,
     PathGenerator_reset__doc__
    },
    {NULL}  /* Sentinel */
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
    Mode mode;
    Algorithm algorithm;
    double match;
    double mismatch;
    double epsilon;
    double target_internal_open_gap_score;
    double target_internal_extend_gap_score;
    double target_left_open_gap_score;
    double target_left_extend_gap_score;
    double target_right_open_gap_score;
    double target_right_extend_gap_score;
    double query_internal_open_gap_score;
    double query_internal_extend_gap_score;
    double query_left_open_gap_score;
    double query_left_extend_gap_score;
    double query_right_open_gap_score;
    double query_right_extend_gap_score;
    PyObject* target_gap_function;
    PyObject* query_gap_function;
    Py_buffer substitution_matrix;
    PyObject* alphabet;
    signed char mapping[128];
} Aligner;


static Py_ssize_t
set_alphabet(Aligner* self, PyObject* alphabet)
{
    Py_ssize_t size;
    if (alphabet == NULL) {
        const char letters[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        size = strlen(letters);
        alphabet = PyUnicode_FromString(letters);
        if (!alphabet) {
            PyErr_SetString(PyExc_ValueError, "failed to initialize alphabet");
            return -1;
        }
    }
    else if (alphabet == Py_None) {
        if (self->alphabet) {
            Py_DECREF(self->alphabet);
            self->alphabet = NULL;
        }
        return 0;
    }
    else if (PyUnicode_Check(alphabet)) {
        if (PyUnicode_READY(alphabet) == -1) return 0;
        size = PyUnicode_GET_LENGTH(alphabet);
        if (size == 0) {
            PyErr_SetString(PyExc_ValueError, "alphabet has zero length");
            return 0;
        }
        if (PyUnicode_KIND(alphabet) == PyUnicode_1BYTE_KIND) {
            int i;
            /* don't set self->mapping until we are sure there are no
             * problems with the alphabet; define temporary storage: */
            signed char mapping[128];
            const char* characters = PyUnicode_DATA(alphabet);
            for (i = 0; i < 128; i++) mapping[i] = MISSING_LETTER;
            for (i = 0; i < size; i++) {
                int j = characters[i];
                if (j < 0 || j >= 128) break;
                if (mapping[j] != MISSING_LETTER) {
                    PyErr_Format(PyExc_ValueError,
                                 "alphabet contains '%c' more than once",
                                 characters[i]);
                    return 0;
                }
                mapping[j] = i;
            }
            if (i < size) {
                /* alphabet is not an ASCII string; cannot use mapping */
                *(self->mapping) = UNMAPPED;
            }
            else {
                /* alphabet is an ASCII string; use mapping for speed */
                memcpy(self->mapping, mapping, 128);
            }
        }
        Py_INCREF(alphabet);
    }
    else {
        /* alphabet is not a string; cannot use mapping */
        int i, j;
        PyObject* item;
        const char* character;
        Py_ssize_t k;
        PyObject* sequence;
        signed char mapping[128];
        sequence = PySequence_Fast(alphabet,
            "alphabet should support the sequence protocol (e.g.,\n"
            "strings, lists, and tuples can be valid alphabets).");
        if (!sequence) return -1;
        size = PySequence_Fast_GET_SIZE(sequence);
        for (i = 0; i < 128; i++) mapping[i] = MISSING_LETTER;
        for (i = 0; i < size; i++) {
            item = PySequence_Fast_GET_ITEM(sequence, i);
            if (!PyUnicode_Check(item)) break;
            if (PyUnicode_READY(item) == -1) {
                Py_DECREF(sequence);
                return -1;
            }
            k = PyUnicode_GET_LENGTH(item);
            if (k != 1) break;
            if (PyUnicode_KIND(item) != PyUnicode_1BYTE_KIND) break;
            character = PyUnicode_DATA(item);
            j = *character;
            if (j < 0 || j >= 128) break;
            if (mapping[j] != MISSING_LETTER) {
                PyErr_Format(PyExc_ValueError,
                             "alphabet contains '%c' more than once",
                             *character);
                return 0;
            }
            mapping[j] = i;
        }
        if (i < size) {
            /* alphabet is not an ASCII string; cannot use mapping */
            *(self->mapping) = UNMAPPED;
        }
        else {
            /* alphabet is an ASCII string; use mapping for speed */
            memcpy(self->mapping, mapping, 128);
        }
        Py_DECREF(sequence);
        Py_INCREF(alphabet);
    }
    Py_XDECREF(self->alphabet);
    self->alphabet = alphabet;
    return size;
}

static int
Aligner_init(Aligner *self, PyObject *args, PyObject *kwds)
{
    int i, j;
    Py_ssize_t n;
    self->alphabet = NULL;
    n = set_alphabet(self, NULL);
    if (n < 0) return -1;
    self->mode = Global;
    self->match = 1.0;
    self->mismatch = 0.0;
    self->epsilon = 1.e-6;
    self->target_internal_open_gap_score = 0;
    self->target_internal_extend_gap_score = 0;
    self->query_internal_open_gap_score = 0;
    self->query_internal_extend_gap_score = 0;
    self->target_left_open_gap_score = 0;
    self->target_left_extend_gap_score = 0;
    self->target_right_open_gap_score = 0;
    self->target_right_extend_gap_score = 0;
    self->query_left_open_gap_score = 0;
    self->query_left_extend_gap_score = 0;
    self->query_right_open_gap_score = 0;
    self->query_right_extend_gap_score = 0;
    self->target_gap_function = NULL;
    self->query_gap_function = NULL;
    self->substitution_matrix.obj = NULL;
    self->substitution_matrix.buf = NULL;
    self->algorithm = Unknown;
    for (i = 0; i < 128; i++) self->mapping[i] = MISSING_LETTER;
    i = (int)'A';
    for (j = 0; j < n; i++, j++) self->mapping[i] = j;
    i = (int)'a';
    for (j = 0; j < n; i++, j++) self->mapping[i] = j;
    /* X represents an undetermined letter. */
    i = (int)'x';
    self->mapping[i] = ANY_LETTER;
    i = (int)'X';
    self->mapping[i] = ANY_LETTER;
    return 0;
}

static void
Aligner_dealloc(Aligner* self)
{   Py_XDECREF(self->target_gap_function);
    Py_XDECREF(self->query_gap_function);
    if (self->substitution_matrix.obj) PyBuffer_Release(&self->substitution_matrix);
    Py_XDECREF(self->alphabet);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject*
Aligner_repr(Aligner* self)
{
  const char text[] = "Pairwise aligner, implementing the Needleman-Wunsch, Smith-Waterman, Gotoh, and Waterman-Smith-Beyer global and local alignment algorithms";
  return PyUnicode_FromString(text);
}

static PyObject*
Aligner_str(Aligner* self)
{
    int n;
    char text[1024];
    char* p = text;
    PyObject* substitution_matrix = self->substitution_matrix.obj;
    n = sprintf(text, "Pairwise sequence aligner with parameters\n");
    p += n;
    if (substitution_matrix) {
        n = sprintf(p, "  substitution_matrix: <%s object at %p>\n",
                       Py_TYPE(substitution_matrix)->tp_name, substitution_matrix);
        p += n;
    } else {
        n = sprintf(p, "  match_score: %f\n", self->match);
        p += n;
        n = sprintf(p, "  mismatch_score: %f\n", self->mismatch);
        p += n;
    }
    if (self->target_gap_function) {
        n = sprintf(p, "  target_gap_function: %%R\n");
        p += n;
    }
    else {
        n = sprintf(p, "  target_internal_open_gap_score: %f\n",
                       self->target_internal_open_gap_score);
        p += n;
        n = sprintf(p, "  target_internal_extend_gap_score: %f\n",
                       self->target_internal_extend_gap_score);
        p += n;
        n = sprintf(p, "  target_left_open_gap_score: %f\n",
                       self->target_left_open_gap_score);
        p += n;
        n = sprintf(p, "  target_left_extend_gap_score: %f\n",
                       self->target_left_extend_gap_score);
        p += n;
        n = sprintf(p, "  target_right_open_gap_score: %f\n",
                       self->target_right_open_gap_score);
        p += n;
        n = sprintf(p, "  target_right_extend_gap_score: %f\n",
                       self->target_right_extend_gap_score);
        p += n;
    }
    if (self->query_gap_function) {
        n = sprintf(p, "  query_gap_function: %%R\n");
        p += n;
    }
    else {
        n = sprintf(p, "  query_internal_open_gap_score: %f\n",
                       self->query_internal_open_gap_score);
        p += n;
        n = sprintf(p, "  query_internal_extend_gap_score: %f\n",
                       self->query_internal_extend_gap_score);
        p += n;
        n = sprintf(p, "  query_left_open_gap_score: %f\n",
                       self->query_left_open_gap_score);
        p += n;
        n = sprintf(p, "  query_left_extend_gap_score: %f\n",
                       self->query_left_extend_gap_score);
        p += n;
        n = sprintf(p, "  query_right_open_gap_score: %f\n",
                       self->query_right_open_gap_score);
        p += n;
        n = sprintf(p, "  query_right_extend_gap_score: %f\n",
                       self->query_right_extend_gap_score);
        p += n;
    }
    switch (self->mode) {
        case Global: n = sprintf(p, "  mode: global\n"); break;
        case Local: n = sprintf(p, "  mode: local\n"); break;
    }
    p += n;
    if (self->target_gap_function || self->query_gap_function)
        return PyUnicode_FromFormat(text, self->target_gap_function, self->query_gap_function);
    else if (self->target_gap_function)
        return PyUnicode_FromFormat(text, self->target_gap_function);
    else if (self->query_gap_function)
        return PyUnicode_FromFormat(text, self->query_gap_function);
    else
        return PyUnicode_FromString(text);
}

static char Aligner_mode__doc__[] = "alignment mode ('global' or 'local')";

static PyObject*
Aligner_get_mode(Aligner* self, void* closure)
{   const char* message = NULL;
    switch (self->mode) {
        case Global: message = "global"; break;
        case Local: message = "local"; break;
    }
    return PyUnicode_FromString(message);
}

static int
Aligner_set_mode(Aligner* self, PyObject* value, void* closure)
{
    if (PyUnicode_Check(value)) {
        if (PyUnicode_CompareWithASCIIString(value, "global") == 0) {
            self->mode = Global;
            return 0;
        }
        if (PyUnicode_CompareWithASCIIString(value, "local") == 0) {
            self->mode = Local;
            return 0;
        }
    }
    PyErr_SetString(PyExc_ValueError,
                    "invalid mode (expected 'global' or 'local'");
    return -1;
}

static char Aligner_match_score__doc__[] = "match score";

static PyObject*
Aligner_get_match_score(Aligner* self, void* closure)
{   if (self->substitution_matrix.obj) {
        Py_INCREF(Py_None);
        return Py_None;
    }
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
    if (self->substitution_matrix.obj) {
        if (set_alphabet(self, NULL) < 0) return -1;
        PyBuffer_Release(&self->substitution_matrix);
    }
    self->match = match;
    return 0;
}

static char Aligner_mismatch_score__doc__[] = "mismatch score";

static PyObject*
Aligner_get_mismatch_score(Aligner* self, void* closure)
{   if (self->substitution_matrix.obj) {
        Py_INCREF(Py_None);
        return Py_None;
    }
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
    if (self->substitution_matrix.obj) {
        if (set_alphabet(self, NULL) < 0) return -1;
        PyBuffer_Release(&self->substitution_matrix);
    }
    self->mismatch = mismatch;
    return 0;
}

static char Aligner_substitution_matrix__doc__[] = "substitution_matrix";

static PyObject*
Aligner_get_substitution_matrix(Aligner* self, void* closure)
{   PyObject* object = self->substitution_matrix.obj;
    if (!object) object = Py_None;
    Py_INCREF(object);
    return object;
}

static int
Aligner_set_substitution_matrix(Aligner* self, PyObject* values, void* closure)
{
    PyObject* alphabet;
    Py_ssize_t size = -1;
    Py_buffer view;
    const int flag = PyBUF_FORMAT | PyBUF_ND;
    if (values == Py_None) {
        if (self->substitution_matrix.obj)
            PyBuffer_Release(&self->substitution_matrix);
        return 0;
    }
    if (PyObject_GetBuffer(values, &view, flag) != 0) {
        PyErr_SetString(PyExc_ValueError, "expected a matrix");
        return -1;
    }
    if (view.ndim != 2) {
        PyErr_Format(PyExc_ValueError,
         "substitution matrix has incorrect rank (%d expected 2)",
          view.ndim);
        PyBuffer_Release(&view);
        return -1;
    }
    if (view.len == 0) {
        PyErr_SetString(PyExc_ValueError, "substitution matrix has zero size");
        PyBuffer_Release(&view);
        return -1;
    }
    if (strcmp(view.format, "d") != 0) {
        PyErr_SetString(PyExc_ValueError,
                "substitution matrix should contain float values");
        PyBuffer_Release(&view);
        return -1;
    }
    if (view.itemsize != sizeof(double)) {
        PyErr_Format(PyExc_RuntimeError,
                    "substitution matrix has unexpected item byte size "
                    "(%zd, expected %zd)", view.itemsize, sizeof(double));
        PyBuffer_Release(&view);
        return -1;
    }
    if (view.shape[0] != view.shape[1]) {
        PyErr_Format(PyExc_ValueError,
                    "substitution matrix should be square "
                    "(found a %zd x %zd matrix)",
                    view.shape[0], view.shape[1]);
        PyBuffer_Release(&view);
        return -1;
    }
    alphabet = PyObject_GetAttrString(values, "alphabet");
    if (alphabet) {
        size = set_alphabet(self, alphabet);
        Py_DECREF(alphabet);
    } else {
        /* Set a substitution matrix without setting an alphabet; useful
         * when aligning integers. */
        PyErr_Clear();
        size = set_alphabet(self, NULL);
    }
    if (size < 0) {
        PyBuffer_Release(&view);
        return -1;
    }
    if (self->substitution_matrix.obj) PyBuffer_Release(&self->substitution_matrix);
    self->substitution_matrix = view;
    return 0;
}

static char Aligner_alphabet__doc__[] = "alphabet";

static PyObject*
Aligner_get_alphabet(Aligner* self, void* closure)
{   PyObject* object = self->alphabet;
    if (!object) object = Py_None;
    Py_INCREF(object);
    return object;
}

static int
Aligner_set_alphabet(Aligner* self, PyObject* alphabet, void* closure)
{
    if (self->substitution_matrix.obj) {
        PyErr_SetString(PyExc_AttributeError,
            "can't set alphabet if a substitution matrix is used");
        return -1;
    }
    if (set_alphabet(self, alphabet) < 0) return -1;
    return 0;
}

static char Aligner_gap_score__doc__[] = "gap score";

static PyObject*
Aligner_get_gap_score(Aligner* self, void* closure)
{   
    if (self->target_gap_function || self->query_gap_function) {
        if (self->target_gap_function != self->query_gap_function) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        Py_INCREF(self->target_gap_function);
        return self->target_gap_function;
    }
    else {
        const double score = self->target_internal_open_gap_score;
        if (score != self->target_internal_extend_gap_score
         || score != self->target_left_open_gap_score
         || score != self->target_left_extend_gap_score
         || score != self->target_right_open_gap_score
         || score != self->target_right_extend_gap_score
         || score != self->query_internal_open_gap_score
         || score != self->query_internal_extend_gap_score
         || score != self->query_left_open_gap_score
         || score != self->query_left_extend_gap_score
         || score != self->query_right_open_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_gap_score(Aligner* self, PyObject* value, void* closure)
{   if (PyCallable_Check(value)) {
        Py_XDECREF(self->target_gap_function);
        Py_XDECREF(self->query_gap_function);
        Py_INCREF(value);
        Py_INCREF(value);
        self->target_gap_function = value;
        self->query_gap_function = value;
    }
    else {
        const double score = PyFloat_AsDouble(value);
        if (PyErr_Occurred()) return -1;
        if (self->target_gap_function) {
            Py_DECREF(self->target_gap_function);
            self->target_gap_function = NULL;
        }
        if (self->query_gap_function) {
            Py_DECREF(self->query_gap_function);
            self->query_gap_function = NULL;
        }
        self->target_internal_open_gap_score = score;
        self->target_internal_extend_gap_score = score;
        self->target_left_open_gap_score = score;
        self->target_left_extend_gap_score = score;
        self->target_right_open_gap_score = score;
        self->target_right_extend_gap_score = score;
        self->query_internal_open_gap_score = score;
        self->query_internal_extend_gap_score = score;
        self->query_left_open_gap_score = score;
        self->query_left_extend_gap_score = score;
        self->query_right_open_gap_score = score;
        self->query_right_extend_gap_score = score;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_open_gap_score__doc__[] = "internal and end open gap score";

static PyObject*
Aligner_get_open_gap_score(Aligner* self, void* closure)
{   
    if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_internal_open_gap_score;
        if (score != self->target_left_open_gap_score
         || score != self->target_right_open_gap_score
         || score != self->query_internal_open_gap_score
         || score != self->query_left_open_gap_score
         || score != self->query_right_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_internal_open_gap_score = score;
    self->target_left_open_gap_score = score;
    self->target_right_open_gap_score = score;
    self->query_internal_open_gap_score = score;
    self->query_left_open_gap_score = score;
    self->query_right_open_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_extend_gap_score__doc__[] = "extend gap score";

static PyObject*
Aligner_get_extend_gap_score(Aligner* self, void* closure)
{   
    if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_internal_extend_gap_score;
        if (score != self->target_left_extend_gap_score
         || score != self->target_right_extend_gap_score
         || score != self->query_internal_extend_gap_score
         || score != self->query_left_extend_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_internal_extend_gap_score = score;
    self->target_left_extend_gap_score = score;
    self->target_right_extend_gap_score = score;
    self->query_internal_extend_gap_score = score;
    self->query_left_extend_gap_score = score;
    self->query_right_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_internal_gap_score__doc__[] = "internal gap score";

static PyObject*
Aligner_get_internal_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_internal_open_gap_score;
        if (score != self->target_internal_extend_gap_score
         || score != self->query_internal_open_gap_score
         || score != self->query_internal_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_internal_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_internal_open_gap_score = score;
    self->target_internal_extend_gap_score = score;
    self->query_internal_open_gap_score = score;
    self->query_internal_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_internal_open_gap_score__doc__[] = "internal open gap score";

static PyObject*
Aligner_get_internal_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_internal_open_gap_score;
        if (score != self->query_internal_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_internal_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_internal_open_gap_score = score;
    self->query_internal_open_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_internal_extend_gap_score__doc__[] = "internal extend gap score";

static PyObject*
Aligner_get_internal_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_internal_extend_gap_score;
        if (score != self->query_internal_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_internal_extend_gap_score(Aligner* self, PyObject* value,
                                      void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_internal_extend_gap_score = score;
    self->query_internal_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_end_gap_score__doc__[] = "end gap score";

static PyObject*
Aligner_get_end_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_open_gap_score;
        if (score != self->target_left_extend_gap_score
         || score != self->target_right_open_gap_score
         || score != self->target_right_extend_gap_score
         || score != self->query_left_open_gap_score
         || score != self->query_left_extend_gap_score
         || score != self->query_right_open_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_end_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_left_open_gap_score = score;
    self->target_left_extend_gap_score = score;
    self->target_right_open_gap_score = score;
    self->target_right_extend_gap_score = score;
    self->query_left_open_gap_score = score;
    self->query_left_extend_gap_score = score;
    self->query_right_open_gap_score = score;
    self->query_right_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_end_open_gap_score__doc__[] = "end open gap score";

static PyObject*
Aligner_get_end_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_open_gap_score;
        if (score != self->target_right_open_gap_score
         || score != self->query_left_open_gap_score
         || score != self->query_right_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_end_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_left_open_gap_score = score;
    self->target_right_open_gap_score = score;
    self->query_left_open_gap_score = score;
    self->query_right_open_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_end_extend_gap_score__doc__[] = "end extend gap score";

static PyObject*
Aligner_get_end_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_extend_gap_score;
        if (score != self->target_right_extend_gap_score
         || score != self->query_left_extend_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_end_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_left_extend_gap_score = score;
    self->target_right_extend_gap_score = score;
    self->query_left_extend_gap_score = score;
    self->query_right_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_left_gap_score__doc__[] = "left gap score";

static PyObject*
Aligner_get_left_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_open_gap_score;
        if (score != self->target_left_extend_gap_score
         || score != self->query_left_open_gap_score
         || score != self->query_left_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_left_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_left_open_gap_score = score;
    self->target_left_extend_gap_score = score;
    self->query_left_open_gap_score = score;
    self->query_left_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_right_gap_score__doc__[] = "right gap score";

static PyObject*
Aligner_get_right_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_right_open_gap_score;
        if (score != self->target_right_extend_gap_score
         || score != self->query_right_open_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_right_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_right_open_gap_score = score;
    self->target_right_extend_gap_score = score;
    self->query_right_open_gap_score = score;
    self->query_right_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_left_open_gap_score__doc__[] = "left open gap score";

static PyObject*
Aligner_get_left_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_open_gap_score;
        if (score != self->query_left_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_left_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_left_open_gap_score = score;
    self->query_left_open_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_left_extend_gap_score__doc__[] = "left extend gap score";

static PyObject*
Aligner_get_left_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_extend_gap_score;
        if (score != self->query_left_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_left_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_left_extend_gap_score = score;
    self->query_left_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_right_open_gap_score__doc__[] = "right open gap score";

static PyObject*
Aligner_get_right_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_right_open_gap_score;
        if (score != self->query_right_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_right_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_right_open_gap_score = score;
    self->query_right_open_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_right_extend_gap_score__doc__[] = "right extend gap score";

static PyObject*
Aligner_get_right_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_right_extend_gap_score;
        if (score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_right_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_right_extend_gap_score = score;
    self->query_right_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_open_gap_score__doc__[] = "target open gap score";

static PyObject*
Aligner_get_target_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_internal_open_gap_score;
        if (score != self->target_left_open_gap_score
         || score != self->target_right_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_internal_open_gap_score = score;
    self->target_left_open_gap_score = score;
    self->target_right_open_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_extend_gap_score__doc__[] = "target extend gap score";

static PyObject*
Aligner_get_target_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_internal_extend_gap_score;
        if (score != self->target_left_extend_gap_score
         || score != self->target_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_internal_extend_gap_score = score;
    self->target_left_extend_gap_score = score;
    self->target_right_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_gap_score__doc__[] = "target gap score";

static PyObject*
Aligner_get_target_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        Py_INCREF(self->target_gap_function);
        return self->target_gap_function;
    }
    else {
        const double score = self->target_internal_open_gap_score;
        if (score != self->target_internal_extend_gap_score
         || score != self->target_left_open_gap_score
         || score != self->target_left_extend_gap_score
         || score != self->target_right_open_gap_score
         || score != self->target_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_gap_score(Aligner* self, PyObject* value, void* closure)
{
    if (PyCallable_Check(value)) {
        Py_XDECREF(self->target_gap_function);
        Py_INCREF(value);
        self->target_gap_function = value;
    }
    else {
        const double score = PyFloat_AsDouble(value);
        if (PyErr_Occurred()) {
            PyErr_SetString(PyExc_ValueError,
                            "gap score should be numerical or callable");
            return -1;
        }
        self->target_internal_open_gap_score = score;
        self->target_internal_extend_gap_score = score;
        self->target_left_open_gap_score = score;
        self->target_left_extend_gap_score = score;
        self->target_right_open_gap_score = score;
        self->target_right_extend_gap_score = score;
        if (self->target_gap_function) {
            Py_DECREF(self->target_gap_function);
            self->target_gap_function = NULL;
        }
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_open_gap_score__doc__[] = "query gap open score";

static PyObject*
Aligner_get_query_open_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_internal_open_gap_score;
        if (score != self->query_left_open_gap_score
         || score != self->query_right_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_internal_open_gap_score = score;
    self->query_left_open_gap_score = score;
    self->query_right_open_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_extend_gap_score__doc__[] = "query gap extend score";

static PyObject*
Aligner_get_query_extend_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_internal_extend_gap_score;
        if (score != self->query_left_extend_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_internal_extend_gap_score = score;
    self->query_left_extend_gap_score = score;
    self->query_right_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_gap_score__doc__[] = "query gap score";

static PyObject*
Aligner_get_query_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        Py_INCREF(self->query_gap_function);
        return self->query_gap_function;
    }
    else {
        const double score = self->query_internal_open_gap_score;
        if (score != self->query_left_open_gap_score
         || score != self->query_right_open_gap_score
         || score != self->query_internal_extend_gap_score
         || score != self->query_left_extend_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_gap_score(Aligner* self, PyObject* value, void* closure)
{   if (PyCallable_Check(value)) {
        Py_XDECREF(self->query_gap_function);
        Py_INCREF(value);
        self->query_gap_function = value;
    }
    else {
        const double score = PyFloat_AsDouble(value);
        if (PyErr_Occurred()) {
            PyErr_SetString(PyExc_ValueError,
                            "gap score should be numerical or callable");
            return -1;
        }
        self->query_internal_open_gap_score = score;
        self->query_internal_extend_gap_score = score;
        self->query_left_open_gap_score = score;
        self->query_left_extend_gap_score = score;
        self->query_right_open_gap_score = score;
        self->query_right_extend_gap_score = score;
        if (self->query_gap_function) {
            Py_DECREF(self->query_gap_function);
            self->query_gap_function = NULL;
        }
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_internal_open_gap_score__doc__[] = "target internal open gap score";

static PyObject*
Aligner_get_target_internal_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->target_internal_open_gap_score);
}

static int
Aligner_set_target_internal_open_gap_score(Aligner* self,
                                           PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_internal_open_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_internal_extend_gap_score__doc__[] = "target internal extend gap score";

static PyObject*
Aligner_get_target_internal_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->target_internal_extend_gap_score);
}

static int
Aligner_set_target_internal_extend_gap_score(Aligner* self,
                                             PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_internal_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_internal_gap_score__doc__[] = "target internal gap score";

static PyObject*
Aligner_get_target_internal_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_internal_open_gap_score;
        if (score != self->target_internal_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_internal_gap_score(Aligner* self, PyObject* value,
                                      void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_internal_open_gap_score = score;
    self->target_internal_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_end_gap_score__doc__[] = "target end gap score";

static PyObject*
Aligner_get_target_end_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_open_gap_score;
        if (score != self->target_left_extend_gap_score
         || score != self->target_right_open_gap_score
         || score != self->target_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_end_gap_score(Aligner* self, PyObject* value, void* closure) {
    const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_left_open_gap_score = score;
    self->target_left_extend_gap_score = score;
    self->target_right_open_gap_score = score;
    self->target_right_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_end_open_gap_score__doc__[] = "target end open gap score";

static PyObject*
Aligner_get_target_end_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_open_gap_score;
        if (score != self->target_right_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_end_open_gap_score(Aligner* self, PyObject* value,
                                      void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_left_open_gap_score = score;
    self->target_right_open_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_end_extend_gap_score__doc__[] = "target end extend gap score";

static PyObject*
Aligner_get_target_end_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_extend_gap_score;
        if (score != self->target_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_end_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_left_extend_gap_score = score;
    self->target_right_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_left_open_gap_score__doc__[] = "target left open score";

static PyObject*
Aligner_get_target_left_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->target_left_open_gap_score);
}

static int
Aligner_set_target_left_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_left_open_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_left_extend_gap_score__doc__[] = "target left extend score";

static PyObject*
Aligner_get_target_left_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->target_left_extend_gap_score);
}

static int
Aligner_set_target_left_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_left_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_left_gap_score__doc__[] = "target left score";

static PyObject*
Aligner_get_target_left_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_open_gap_score;
        if (score != self->target_left_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_left_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_left_open_gap_score = score;
    self->target_left_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_right_gap_score_open__doc__[] = "target right open score";

static PyObject*
Aligner_get_target_right_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->target_right_open_gap_score);
}

static int
Aligner_set_target_right_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_right_open_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_right_extend_gap_score__doc__[] = "target right extend score";

static PyObject*
Aligner_get_target_right_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->target_right_extend_gap_score);
}

static int
Aligner_set_target_right_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_right_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_right_gap_score__doc__[] = "target right score";

static PyObject*
Aligner_get_target_right_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_right_open_gap_score;
        if (score != self->target_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_right_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_right_open_gap_score = score;
    self->target_right_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_end_gap_score__doc__[] = "query end score";

static PyObject*
Aligner_get_query_end_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_left_open_gap_score;
        if (score != self->query_left_extend_gap_score
         || score != self->query_right_open_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_end_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_left_open_gap_score = score;
    self->query_left_extend_gap_score = score;
    self->query_right_open_gap_score = score;
    self->query_right_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_end_open_gap_score__doc__[] = "query end open score";

static PyObject*
Aligner_get_query_end_open_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_left_open_gap_score;
        if (score != self->query_right_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_end_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_left_open_gap_score = score;
    self->query_right_open_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_end_extend_gap_score__doc__[] = "query end extend score";

static PyObject*
Aligner_get_query_end_extend_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_left_extend_gap_score;
        if (score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_end_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_left_extend_gap_score = score;
    self->query_right_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_internal_open_gap_score__doc__[] = "query internal open gap score";

static PyObject*
Aligner_get_query_internal_open_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->query_internal_open_gap_score);
}

static int
Aligner_set_query_internal_open_gap_score(Aligner* self, PyObject* value,
                                          void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_internal_open_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_internal_extend_gap_score__doc__[] = "query internal extend gap score";

static PyObject*
Aligner_get_query_internal_extend_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->query_internal_extend_gap_score);
}

static int
Aligner_set_query_internal_extend_gap_score(Aligner* self, PyObject* value,
                                            void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_internal_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_internal_gap_score__doc__[] = "query internal gap score";

static PyObject*
Aligner_get_query_internal_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_internal_open_gap_score;
        if (score != self->query_internal_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_internal_gap_score(Aligner* self, PyObject* value,
                                     void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_internal_open_gap_score = score;
    self->query_internal_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_left_open_gap_score__doc__[] = "query left open score";

static PyObject*
Aligner_get_query_left_open_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->query_left_open_gap_score);
}

static int
Aligner_set_query_left_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_left_open_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_left_extend_gap_score__doc__[] = "query left extend score";

static PyObject*
Aligner_get_query_left_extend_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->query_left_extend_gap_score);
}

static int
Aligner_set_query_left_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_left_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_left_gap_score__doc__[] = "query left score";

static PyObject*
Aligner_get_query_left_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_left_open_gap_score;
        if (score != self->query_left_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_left_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_left_open_gap_score = score;
    self->query_left_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_right_open_gap_score__doc__[] = "query right open score";

static PyObject*
Aligner_get_query_right_open_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->query_right_open_gap_score);
}

static int
Aligner_set_query_right_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_right_open_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_right_extend_gap_score__doc__[] = "query right extend score";

static PyObject*
Aligner_get_query_right_extend_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->query_right_extend_gap_score);
}

static int
Aligner_set_query_right_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_right_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_right_gap_score__doc__[] = "query right score";

static PyObject*
Aligner_get_query_right_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_right_open_gap_score;
        if (score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_right_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_right_open_gap_score = score;
    self->query_right_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_epsilon__doc__[] = "roundoff epsilon";

static PyObject*
Aligner_get_epsilon(Aligner* self, void* closure)
{   return PyFloat_FromDouble(self->epsilon);
}

static int
Aligner_set_epsilon(Aligner* self, PyObject* value, void* closure)
{   const double epsilon = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->epsilon = epsilon;
    self->algorithm = Unknown;
    return 0;
}

static Algorithm _get_algorithm(Aligner* self)
{
    Algorithm algorithm = self->algorithm;
    if (algorithm == Unknown) {
        const double target_gap_open = self->target_internal_open_gap_score;
        const double query_gap_open = self->query_internal_open_gap_score;
        const double target_gap_extend = self->target_internal_extend_gap_score;
        const double query_gap_extend = self->query_internal_extend_gap_score;
        const double target_left_open = self->target_left_open_gap_score;
        const double target_left_extend = self->target_left_extend_gap_score;
        const double query_left_open = self->query_left_open_gap_score;
        const double target_right_open = self->target_right_open_gap_score;
        const double query_right_open = self->query_right_open_gap_score;
        const double target_right_extend = self->target_right_extend_gap_score;
        const double query_left_extend = self->query_left_extend_gap_score;
        const double query_right_extend = self->query_right_extend_gap_score;
        if (self->target_gap_function || self->query_gap_function)
            algorithm = WatermanSmithBeyer;
        else if (target_gap_open == target_gap_extend
              && query_gap_open == query_gap_extend
              && target_left_open == target_left_extend
              && target_right_open == target_right_extend
              && query_left_open == query_left_extend
              && query_right_open == query_right_extend)
            algorithm = NeedlemanWunschSmithWaterman;
        else
            algorithm = Gotoh;
        self->algorithm = algorithm;
    }
    return algorithm;
}


static char Aligner_algorithm__doc__[] = "alignment algorithm";

static PyObject*
Aligner_get_algorithm(Aligner* self, void* closure)
{
    const char* s = NULL;
    const Mode mode = self->mode;
    const Algorithm algorithm = _get_algorithm(self);
    switch (algorithm) {
        case NeedlemanWunschSmithWaterman:
            switch (mode) {
                case Global:
                    s = "Needleman-Wunsch";
                    break;
                case Local:
                    s = "Smith-Waterman";
                    break;
            }
            break;
        case Gotoh:
            switch (mode) {
                case Global:
                    s = "Gotoh global alignment algorithm";
                    break;
                case Local:
                    s = "Gotoh local alignment algorithm";
                    break;
            }
            break;
        case WatermanSmithBeyer:
            switch (mode) {
                case Global:
                    s = "Waterman-Smith-Beyer global alignment algorithm";
                    break;
                case Local:
                    s = "Waterman-Smith-Beyer local alignment algorithm";
                    break;
            }
            break;
        case Unknown:
        default:
            break;
    }
    return PyUnicode_FromString(s);
}

static PyGetSetDef Aligner_getset[] = {
    {"mode",
        (getter)Aligner_get_mode,
        (setter)Aligner_set_mode,
        Aligner_mode__doc__, NULL},
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
    {"substitution_matrix",
        (getter)Aligner_get_substitution_matrix,
        (setter)Aligner_set_substitution_matrix,
        Aligner_substitution_matrix__doc__, NULL},
    {"alphabet",
        (getter)Aligner_get_alphabet,
        (setter)Aligner_set_alphabet,
        Aligner_alphabet__doc__, NULL},
    {"gap_score",
        (getter)Aligner_get_gap_score,
        (setter)Aligner_set_gap_score,
        Aligner_gap_score__doc__, NULL},
    {"open_gap_score",
        (getter)Aligner_get_open_gap_score,
        (setter)Aligner_set_open_gap_score,
        Aligner_open_gap_score__doc__, NULL},
    {"extend_gap_score",
        (getter)Aligner_get_extend_gap_score,
        (setter)Aligner_set_extend_gap_score,
        Aligner_extend_gap_score__doc__, NULL},
    {"internal_gap_score",
        (getter)Aligner_get_internal_gap_score,
        (setter)Aligner_set_internal_gap_score,
        Aligner_internal_gap_score__doc__, NULL},
    {"internal_open_gap_score",
        (getter)Aligner_get_internal_open_gap_score,
        (setter)Aligner_set_internal_open_gap_score,
        Aligner_internal_open_gap_score__doc__, NULL},
    {"internal_extend_gap_score",
        (getter)Aligner_get_internal_extend_gap_score,
        (setter)Aligner_set_internal_extend_gap_score,
        Aligner_internal_extend_gap_score__doc__, NULL},
    {"end_gap_score",
        (getter)Aligner_get_end_gap_score,
        (setter)Aligner_set_end_gap_score,
        Aligner_end_gap_score__doc__, NULL},
    {"end_open_gap_score",
        (getter)Aligner_get_end_open_gap_score,
        (setter)Aligner_set_end_open_gap_score,
        Aligner_end_open_gap_score__doc__, NULL},
    {"end_extend_gap_score",
        (getter)Aligner_get_end_extend_gap_score,
        (setter)Aligner_set_end_extend_gap_score,
        Aligner_end_extend_gap_score__doc__, NULL},
    {"left_gap_score",
        (getter)Aligner_get_left_gap_score,
        (setter)Aligner_set_left_gap_score,
        Aligner_left_gap_score__doc__, NULL},
    {"left_open_gap_score",
        (getter)Aligner_get_left_open_gap_score,
        (setter)Aligner_set_left_open_gap_score,
        Aligner_left_open_gap_score__doc__, NULL},
    {"left_extend_gap_score",
        (getter)Aligner_get_left_extend_gap_score,
        (setter)Aligner_set_left_extend_gap_score,
        Aligner_left_extend_gap_score__doc__, NULL},
    {"right_gap_score",
        (getter)Aligner_get_right_gap_score,
        (setter)Aligner_set_right_gap_score,
        Aligner_right_gap_score__doc__, NULL},
    {"right_open_gap_score",
        (getter)Aligner_get_right_open_gap_score,
        (setter)Aligner_set_right_open_gap_score,
        Aligner_right_open_gap_score__doc__, NULL},
    {"right_extend_gap_score",
        (getter)Aligner_get_right_extend_gap_score,
        (setter)Aligner_set_right_extend_gap_score,
        Aligner_right_extend_gap_score__doc__, NULL},
    {"target_open_gap_score",
        (getter)Aligner_get_target_open_gap_score,
        (setter)Aligner_set_target_open_gap_score,
        Aligner_target_open_gap_score__doc__, NULL},
    {"target_extend_gap_score",
        (getter)Aligner_get_target_extend_gap_score,
        (setter)Aligner_set_target_extend_gap_score,
        Aligner_target_extend_gap_score__doc__, NULL},
    {"target_gap_score",
        (getter)Aligner_get_target_gap_score,
        (setter)Aligner_set_target_gap_score,
        Aligner_target_gap_score__doc__, NULL},
    {"query_open_gap_score",
        (getter)Aligner_get_query_open_gap_score,
        (setter)Aligner_set_query_open_gap_score,
        Aligner_query_open_gap_score__doc__, NULL},
    {"query_extend_gap_score",
        (getter)Aligner_get_query_extend_gap_score,
        (setter)Aligner_set_query_extend_gap_score,
        Aligner_query_extend_gap_score__doc__, NULL},
    {"query_gap_score",
        (getter)Aligner_get_query_gap_score,
        (setter)Aligner_set_query_gap_score,
        Aligner_query_gap_score__doc__, NULL},
    {"target_end_gap_score",
        (getter)Aligner_get_target_end_gap_score,
        (setter)Aligner_set_target_end_gap_score,
        Aligner_target_end_gap_score__doc__, NULL},
    {"target_end_open_gap_score",
        (getter)Aligner_get_target_end_open_gap_score,
        (setter)Aligner_set_target_end_open_gap_score,
        Aligner_target_end_open_gap_score__doc__, NULL},
    {"target_end_extend_gap_score",
        (getter)Aligner_get_target_end_extend_gap_score,
        (setter)Aligner_set_target_end_extend_gap_score,
        Aligner_target_end_extend_gap_score__doc__, NULL},
    {"target_internal_open_gap_score",
        (getter)Aligner_get_target_internal_open_gap_score,
        (setter)Aligner_set_target_internal_open_gap_score,
        Aligner_target_internal_open_gap_score__doc__, NULL},
    {"target_internal_extend_gap_score",
        (getter)Aligner_get_target_internal_extend_gap_score,
        (setter)Aligner_set_target_internal_extend_gap_score,
        Aligner_target_internal_extend_gap_score__doc__, NULL},
    {"target_internal_gap_score",
        (getter)Aligner_get_target_internal_gap_score,
        (setter)Aligner_set_target_internal_gap_score,
        Aligner_target_internal_gap_score__doc__, NULL},
    {"target_left_open_gap_score",
        (getter)Aligner_get_target_left_open_gap_score,
        (setter)Aligner_set_target_left_open_gap_score,
        Aligner_target_left_open_gap_score__doc__, NULL},
    {"target_left_extend_gap_score",
        (getter)Aligner_get_target_left_extend_gap_score,
        (setter)Aligner_set_target_left_extend_gap_score,
        Aligner_target_left_extend_gap_score__doc__, NULL},
    {"target_left_gap_score",
        (getter)Aligner_get_target_left_gap_score,
        (setter)Aligner_set_target_left_gap_score,
        Aligner_target_left_gap_score__doc__, NULL},
    {"target_right_open_gap_score",
        (getter)Aligner_get_target_right_open_gap_score,
        (setter)Aligner_set_target_right_open_gap_score,
        Aligner_target_right_gap_score_open__doc__, NULL},
    {"target_right_extend_gap_score",
        (getter)Aligner_get_target_right_extend_gap_score,
        (setter)Aligner_set_target_right_extend_gap_score,
        Aligner_target_right_extend_gap_score__doc__, NULL},
    {"target_right_gap_score",
        (getter)Aligner_get_target_right_gap_score,
        (setter)Aligner_set_target_right_gap_score,
        Aligner_target_right_gap_score__doc__, NULL},
    {"query_end_gap_score",
        (getter)Aligner_get_query_end_gap_score,
        (setter)Aligner_set_query_end_gap_score,
        Aligner_query_end_gap_score__doc__, NULL},
    {"query_end_open_gap_score",
        (getter)Aligner_get_query_end_open_gap_score,
        (setter)Aligner_set_query_end_open_gap_score,
        Aligner_query_end_open_gap_score__doc__, NULL},
    {"query_end_extend_gap_score",
        (getter)Aligner_get_query_end_extend_gap_score,
        (setter)Aligner_set_query_end_extend_gap_score,
        Aligner_query_end_extend_gap_score__doc__, NULL},
    {"query_internal_open_gap_score",
        (getter)Aligner_get_query_internal_open_gap_score,
        (setter)Aligner_set_query_internal_open_gap_score,
        Aligner_query_internal_open_gap_score__doc__, NULL},
    {"query_internal_extend_gap_score",
        (getter)Aligner_get_query_internal_extend_gap_score,
        (setter)Aligner_set_query_internal_extend_gap_score,
        Aligner_query_internal_extend_gap_score__doc__, NULL},
    {"query_internal_gap_score",
        (getter)Aligner_get_query_internal_gap_score,
        (setter)Aligner_set_query_internal_gap_score,
        Aligner_query_internal_gap_score__doc__, NULL},
    {"query_left_open_gap_score",
        (getter)Aligner_get_query_left_open_gap_score,
        (setter)Aligner_set_query_left_open_gap_score,
        Aligner_query_left_open_gap_score__doc__, NULL},
    {"query_left_extend_gap_score",
        (getter)Aligner_get_query_left_extend_gap_score,
        (setter)Aligner_set_query_left_extend_gap_score,
        Aligner_query_left_extend_gap_score__doc__, NULL},
    {"query_left_gap_score",
        (getter)Aligner_get_query_left_gap_score,
        (setter)Aligner_set_query_left_gap_score,
         Aligner_query_left_gap_score__doc__, NULL},
    {"query_right_open_gap_score",
        (getter)Aligner_get_query_right_open_gap_score,
        (setter)Aligner_set_query_right_open_gap_score,
        Aligner_query_right_open_gap_score__doc__, NULL},
    {"query_right_extend_gap_score",
        (getter)Aligner_get_query_right_extend_gap_score,
        (setter)Aligner_set_query_right_extend_gap_score,
        Aligner_query_right_extend_gap_score__doc__, NULL},
    {"query_right_gap_score",
        (getter)Aligner_get_query_right_gap_score,
        (setter)Aligner_set_query_right_gap_score,
        Aligner_query_right_gap_score__doc__, NULL},
    {"epsilon",
        (getter)Aligner_get_epsilon,
        (setter)Aligner_set_epsilon,
        Aligner_epsilon__doc__, NULL},
    {"algorithm",
        (getter)Aligner_get_algorithm,
        (setter)NULL,
        Aligner_algorithm__doc__, NULL},
    {NULL}  /* Sentinel */
};

#define SELECT_SCORE_GLOBAL(score1, score2, score3) \
    score = score1; \
    temp = score2; \
    if (temp > score) score = temp; \
    temp = score3; \
    if (temp > score) score = temp;

#define SELECT_SCORE_WATERMAN_SMITH_BEYER(score1, score2) \
    temp = score1 + gapscore; \
    if (temp > score) score = temp; \
    temp = score2 + gapscore; \
    if (temp > score) score = temp;

#define SELECT_SCORE_GOTOH_LOCAL_ALIGN(score1, score2, score3, score4) \
    score = score1; \
    temp = score2; \
    if (temp > score) score = temp; \
    temp = score3; \
    if (temp > score) score = temp; \
    score += score4; \
    if (score < 0) score = 0; \
    else if (score > maximum) maximum = score;

#define SELECT_SCORE_LOCAL3(score1, score2, score3) \
    score = score1; \
    temp = score2; \
    if (temp > score) score = temp; \
    temp = score3; \
    if (temp > score) score = temp; \
    if (score < 0) score = 0; \
    else if (score > maximum) maximum = score;

#define SELECT_SCORE_LOCAL1(score1) \
    score = score1; \
    if (score < 0) score = 0; \
    else if (score > maximum) maximum = score;

#define SELECT_TRACE_NEEDLEMAN_WUNSCH(hgap, vgap, align_score) \
    score = temp + (align_score); \
    trace = DIAGONAL; \
    temp = row[j-1] + hgap; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = HORIZONTAL; \
    } \
    else if (temp > score - epsilon) trace |= HORIZONTAL; \
    temp = row[j] + vgap; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = VERTICAL; \
    } \
    else if (temp > score - epsilon) trace |= VERTICAL; \
    temp = row[j]; \
    row[j] = score; \
    M[i][j].trace = trace;

#define SELECT_TRACE_SMITH_WATERMAN_HVD(align_score) \
    trace = DIAGONAL; \
    score = temp + (align_score); \
    temp = row[j-1] + gap_extend_A; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = HORIZONTAL; \
    } \
    else if (temp > score - epsilon) trace |= HORIZONTAL; \
    temp = row[j] + gap_extend_B; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = VERTICAL; \
    } \
    else if (temp > score - epsilon) trace |= VERTICAL; \
    if (score < epsilon) { \
        score = 0; \
        trace = STARTPOINT; \
    } \
    else if (trace & DIAGONAL && score > maximum - epsilon) { \
        if (score > maximum + epsilon) { \
            for ( ; im < i; im++, jm = 0) \
                for ( ; jm <= nB; jm++) M[im][jm].trace &= ~ENDPOINT; \
            for ( ; jm < j; jm++) M[im][jm].trace &= ~ENDPOINT; \
            im = i; \
            jm = j; \
        } \
        trace |= ENDPOINT; \
    } \
    M[i][j].trace = trace; \
    if (score > maximum) maximum = score; \
    temp = row[j]; \
    row[j] = score;

#define SELECT_TRACE_SMITH_WATERMAN_D(align_score) \
    score = temp + (align_score); \
    trace = DIAGONAL; \
    if (score < epsilon) { \
        score = 0; \
    } \
    else if (trace & DIAGONAL && score > maximum - epsilon) { \
        if (score > maximum + epsilon) { \
            for ( ; im < i; im++, jm = 0) \
                for ( ; jm <= nB; jm++) M[im][jm].trace &= ~ENDPOINT; \
            for ( ; jm < j; jm++) M[im][jm].trace &= ~ENDPOINT; \
            im = i; \
            jm = j; \
        } \
        trace |= ENDPOINT; \
    } \
    M[i][j].trace = trace; \
    if (score > maximum) maximum = score; \
    temp = row[j]; \
    row[j] = score

#define SELECT_TRACE_GOTOH_GLOBAL_GAP(matrix, score1, score2, score3) \
    trace = M_MATRIX; \
    score = score1; \
    temp = score2; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Ix_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Ix_MATRIX; \
    temp = score3; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Iy_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Iy_MATRIX; \
    gaps[i][j].matrix = trace;

#define SELECT_TRACE_GOTOH_GLOBAL_ALIGN \
    trace = M_MATRIX; \
    score = M_temp; \
    temp = Ix_temp; \
    if (temp > score + epsilon) { \
        score = Ix_temp; \
        trace = Ix_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Ix_MATRIX; \
    temp = Iy_temp; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Iy_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Iy_MATRIX; \
    M[i][j].trace = trace;

#define SELECT_TRACE_GOTOH_LOCAL_ALIGN(align_score) \
    trace = M_MATRIX; \
    score = M_temp; \
    if (Ix_temp > score + epsilon) { \
        score = Ix_temp; \
        trace = Ix_MATRIX; \
    } \
    else if (Ix_temp > score - epsilon) trace |= Ix_MATRIX; \
    if (Iy_temp > score + epsilon) { \
        score = Iy_temp; \
        trace = Iy_MATRIX; \
    } \
    else if (Iy_temp > score - epsilon) trace |= Iy_MATRIX; \
    score += (align_score); \
    if (score < epsilon) { \
        score = 0; \
        trace = STARTPOINT; \
    } \
    else if (score > maximum - epsilon) { \
        if (score > maximum + epsilon) { \
            maximum = score; \
            for ( ; im < i; im++, jm = 0) \
                for ( ; jm <= nB; jm++) M[im][jm].trace &= ~ENDPOINT; \
            for ( ; jm < j; jm++) M[im][jm].trace &= ~ENDPOINT; \
            im = i; \
            jm = j; \
        } \
        trace |= ENDPOINT; \
    } \
    M[i][j].trace = trace;

#define SELECT_TRACE_GOTOH_LOCAL_GAP(matrix, score1, score2, score3) \
    trace = M_MATRIX; \
    score = score1; \
    temp = score2; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Ix_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Ix_MATRIX; \
    temp = score3; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Iy_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Iy_MATRIX; \
    if (score < epsilon) { \
        score = -DBL_MAX; \
        trace = 0; \
    } \
    gaps[i][j].matrix = trace;

#define SELECT_TRACE_WATERMAN_SMITH_BEYER_GLOBAL_ALIGN(score4) \
    trace = M_MATRIX; \
    score = M_row[i-1][j-1]; \
    temp = Ix_row[i-1][j-1]; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Ix_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Ix_MATRIX; \
    temp = Iy_row[i-1][j-1]; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Iy_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Iy_MATRIX; \
    M_row[i][j] = score + score4; \
    M[i][j].trace = trace;

#define SELECT_TRACE_WATERMAN_SMITH_BEYER_GAP(score1, score2) \
    temp = score1 + gapscore; \
    if (temp > score - epsilon) { \
        if (temp > score + epsilon) { \
            score = temp; \
            nm = 0; \
            ng = 0; \
        } \
        gapM[nm] = gap; \
        nm++; \
    } \
    temp = score2 + gapscore; \
    if (temp > score - epsilon) { \
        if (temp > score + epsilon) { \
            score = temp; \
            nm = 0; \
            ng = 0; \
        } \
        gapXY[ng] = gap; \
        ng++; \
    }

#define SELECT_TRACE_WATERMAN_SMITH_BEYER_ALIGN(score1, score2, score3, score4) \
    trace = M_MATRIX; \
    score = score1; \
    if (score2 > score + epsilon) { \
        score = score2; \
        trace = Ix_MATRIX; \
    } \
    else if (score2 > score - epsilon) trace |= Ix_MATRIX; \
    if (score3 > score + epsilon) { \
        score = score3; \
        trace = Iy_MATRIX; \
    } \
    else if (score3 > score - epsilon) trace |= Iy_MATRIX; \
    score += score4; \
    if (score < epsilon) { \
        score = 0; \
        trace = STARTPOINT; \
    } \
    else if (score > maximum - epsilon) { \
        if (score > maximum + epsilon) { \
            maximum = score; \
            for ( ; im < i; im++, jm = 0) \
                for ( ; jm <= nB; jm++) M[im][jm].trace &= ~ENDPOINT; \
            for ( ; jm < j; jm++) M[im][jm].trace &= ~ENDPOINT; \
            im = i; \
            jm = j; \
        } \
        trace |= ENDPOINT; \
    } \
    M_row[i][j] = score; \
    M[i][j].trace = trace;

/* ----------------- alignment algorithms ----------------- */

#define NEEDLEMANWUNSCH_SCORE(align_score) \
    int i; \
    int j; \
    int kA; \
    int kB; \
    const double gap_extend_A = self->target_internal_extend_gap_score; \
    const double gap_extend_B = self->query_internal_extend_gap_score; \
    const double left_gap_extend_A = self->target_left_extend_gap_score; \
    const double right_gap_extend_A = self->target_right_extend_gap_score; \
    const double left_gap_extend_B = self->query_left_extend_gap_score; \
    const double right_gap_extend_B = self->query_right_extend_gap_score; \
    double score; \
    double temp; \
    double* row; \
\
    /* Needleman-Wunsch algorithm */ \
    row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!row) return PyErr_NoMemory(); \
\
    /* The top row of the score matrix is a special case, \
     * as there are no previously aligned characters. \
     */ \
    row[0] = 0.0; \
    for (j = 1; j <= nB; j++) row[j] = j * left_gap_extend_A; \
    for (i = 1; i < nA; i++) { \
        kA = sA[i-1]; \
        temp = row[0]; \
        row[0] = i * left_gap_extend_B; \
        for (j = 1; j < nB; j++) { \
            kB = sB[j-1]; \
            SELECT_SCORE_GLOBAL(temp + (align_score), \
                                row[j] + gap_extend_B, \
                                row[j-1] + gap_extend_A); \
            temp = row[j]; \
            row[j] = score; \
        } \
        kB = sB[nB-1]; \
        SELECT_SCORE_GLOBAL(temp + (align_score), \
                            row[nB] + right_gap_extend_B, \
                            row[nB-1] + gap_extend_A); \
        temp = row[nB]; \
        row[nB] = score; \
    } \
    kA = sA[nA-1]; \
    temp = row[0]; \
    row[0] = nA * right_gap_extend_B; \
    for (j = 1; j < nB; j++) { \
        kB = sB[j-1]; \
        SELECT_SCORE_GLOBAL(temp + (align_score), \
                            row[j] + gap_extend_B, \
                            row[j-1] + right_gap_extend_A); \
        temp = row[j]; \
        row[j] = score; \
    } \
    kB = sB[nB-1]; \
    SELECT_SCORE_GLOBAL(temp + (align_score), \
                        row[nB] + right_gap_extend_B, \
                        row[nB-1] + right_gap_extend_A); \
    PyMem_Free(row); \
    return PyFloat_FromDouble(score);


#define SMITHWATERMAN_SCORE(align_score) \
    int i; \
    int j; \
    int kA; \
    int kB; \
    const double gap_extend_A = self->target_internal_extend_gap_score; \
    const double gap_extend_B = self->query_internal_extend_gap_score; \
    double score; \
    double* row; \
    double temp; \
    double maximum = 0; \
\
    /* Smith-Waterman algorithm */ \
    row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!row) return PyErr_NoMemory(); \
\
    /* The top row of the score matrix is a special case, \
     * as there are no previously aligned characters. \
     */ \
    for (j = 0; j <= nB; j++) \
        row[j] = 0; \
    for (i = 1; i < nA; i++) { \
        kA = sA[i-1]; \
        temp = 0; \
        for (j = 1; j < nB; j++) { \
            kB = sB[j-1]; \
            SELECT_SCORE_LOCAL3(temp + (align_score), \
                                row[j] + gap_extend_B, \
                                row[j-1] + gap_extend_A); \
            temp = row[j]; \
            row[j] = score; \
        } \
        kB = sB[nB-1]; \
        SELECT_SCORE_LOCAL1(temp + (align_score)); \
        temp = row[nB]; \
        row[nB] = score; \
    } \
    kA = sA[nA-1]; \
    temp = 0; \
    for (j = 1; j < nB; j++) { \
        kB = sB[j-1]; \
        SELECT_SCORE_LOCAL1(temp + (align_score)); \
        temp = row[j]; \
        row[j] = score; \
    } \
    kB = sB[nB-1]; \
    SELECT_SCORE_LOCAL1(temp + (align_score)); \
    PyMem_Free(row); \
    return PyFloat_FromDouble(maximum);


#define NEEDLEMANWUNSCH_ALIGN(align_score) \
    int i; \
    int j; \
    int kA; \
    int kB; \
    const double gap_extend_A = self->target_internal_extend_gap_score; \
    const double gap_extend_B = self->query_internal_extend_gap_score; \
    const double left_gap_extend_A = self->target_left_extend_gap_score; \
    const double left_gap_extend_B = self->query_left_extend_gap_score; \
    const double right_gap_extend_A = self->target_right_extend_gap_score; \
    const double right_gap_extend_B = self->query_right_extend_gap_score; \
    const double epsilon = self->epsilon; \
    Trace** M; \
    double score; \
    int trace; \
    double temp; \
    double* row = NULL; \
    PathGenerator* paths; \
\
    /* Needleman-Wunsch algorithm */ \
    paths = PathGenerator_create_NWSW(nA, nB, Global); \
    if (!paths) return NULL; \
    row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!row) { \
        Py_DECREF(paths); \
        return PyErr_NoMemory(); \
    } \
    M = paths->M; \
    row[0] = 0; \
    for (j = 1; j <= nB; j++) row[j] = j * left_gap_extend_A; \
    for (i = 1; i < nA; i++) { \
        temp = row[0]; \
        row[0] = i * left_gap_extend_B; \
        kA = sA[i-1]; \
        for (j = 1; j < nB; j++) { \
            kB = sB[j-1]; \
            SELECT_TRACE_NEEDLEMAN_WUNSCH(gap_extend_A, gap_extend_B, align_score); \
        } \
        kB = sB[j-1]; \
        SELECT_TRACE_NEEDLEMAN_WUNSCH(gap_extend_A, right_gap_extend_B, align_score); \
    } \
    temp = row[0]; \
    row[0] = i * left_gap_extend_B; \
    kA = sA[nA-1]; \
    for (j = 1; j < nB; j++) { \
        kB = sB[j-1]; \
        SELECT_TRACE_NEEDLEMAN_WUNSCH(right_gap_extend_A, gap_extend_B, align_score); \
    } \
    kB = sB[j-1]; \
    SELECT_TRACE_NEEDLEMAN_WUNSCH(right_gap_extend_A, right_gap_extend_B, align_score); \
    PyMem_Free(row); \
    M[nA][nB].path = 0; \
    return Py_BuildValue("fN", score, paths);


#define SMITHWATERMAN_ALIGN(align_score) \
    int i; \
    int j; \
    int im = nA; \
    int jm = nB; \
    int kA; \
    int kB; \
    const double gap_extend_A = self->target_internal_extend_gap_score; \
    const double gap_extend_B = self->query_internal_extend_gap_score; \
    const double epsilon = self->epsilon; \
    Trace** M = NULL; \
    double maximum = 0; \
    double score = 0; \
    double* row = NULL; \
    double temp; \
    int trace; \
    PathGenerator* paths = NULL; \
\
    /* Smith-Waterman algorithm */ \
    paths = PathGenerator_create_NWSW(nA, nB, Local); \
    if (!paths) return NULL; \
    row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!row) { \
        Py_DECREF(paths); \
        return PyErr_NoMemory(); \
    } \
    M = paths->M; \
    for (j = 0; j <= nB; j++) row[j] = 0; \
    for (i = 1; i < nA; i++) { \
        temp = 0; \
        kA = sA[i-1]; \
        for (j = 1; j < nB; j++) { \
            kB = sB[j-1]; \
            SELECT_TRACE_SMITH_WATERMAN_HVD(align_score); \
        } \
        kB = sB[nB-1]; \
        SELECT_TRACE_SMITH_WATERMAN_D(align_score); \
    } \
    temp = 0; \
    kA = sA[nA-1]; \
    for (j = 1; j < nB; j++) { \
        kB = sB[j-1]; \
        SELECT_TRACE_SMITH_WATERMAN_D(align_score); \
    } \
    kB = sB[nB-1]; \
    SELECT_TRACE_SMITH_WATERMAN_D(align_score); \
    PyMem_Free(row); \
\
    /* As we don't allow zero-score extensions to alignments, \
     * we need to remove all traces towards an ENDPOINT. \
     * In addition, some points then won't have any path to a STARTPOINT. \
     * Here, use path as a temporary variable to indicate if the point \
     * is reachable from a STARTPOINT. If it is unreachable, remove all \
     * traces from it, and don't allow it to be an ENDPOINT. It may still \
     * be a valid STARTPOINT. */ \
    for (j = 0; j <= nB; j++) M[0][j].path = 1; \
    for (i = 1; i <= nA; i++) { \
        M[i][0].path = 1; \
        for (j = 1; j <= nB; j++) { \
            trace = M[i][j].trace; \
            /* Remove traces to unreachable points. */ \
            if (!M[i-1][j-1].path) trace &= ~DIAGONAL; \
            if (!M[i][j-1].path) trace &= ~HORIZONTAL; \
            if (!M[i-1][j].path) trace &= ~VERTICAL; \
            if (trace & (STARTPOINT | HORIZONTAL | VERTICAL | DIAGONAL)) { \
                /* The point is reachable. */ \
                if (trace & ENDPOINT) M[i][j].path = 0; /* no extensions after ENDPOINT */ \
                else M[i][j].path = 1; \
            } \
            else { \
                /* The point is not reachable. Then it is not a STARTPOINT, \
                 * all traces from it can be removed, and it cannot act as \
                 * an ENDPOINT. */ \
                M[i][j].path = 0; \
                trace = 0; \
            } \
            M[i][j].trace = trace; \
        } \
    } \
    if (maximum == 0) M[0][0].path = NONE; \
    else M[0][0].path = 0; \
    return Py_BuildValue("fN", maximum, paths);


#define GOTOH_GLOBAL_SCORE(align_score) \
    int i; \
    int j; \
    int kA; \
    int kB; \
    const double gap_open_A = self->target_internal_open_gap_score; \
    const double gap_open_B = self->query_internal_open_gap_score; \
    const double gap_extend_A = self->target_internal_extend_gap_score; \
    const double gap_extend_B = self->query_internal_extend_gap_score; \
    const double left_gap_open_A = self->target_left_open_gap_score; \
    const double left_gap_open_B = self->query_left_open_gap_score; \
    const double left_gap_extend_A = self->target_left_extend_gap_score; \
    const double left_gap_extend_B = self->query_left_extend_gap_score; \
    const double right_gap_open_A = self->target_right_open_gap_score; \
    const double right_gap_open_B = self->query_right_open_gap_score; \
    const double right_gap_extend_A = self->target_right_extend_gap_score; \
    const double right_gap_extend_B = self->query_right_extend_gap_score; \
    double* M_row = NULL; \
    double* Ix_row = NULL; \
    double* Iy_row = NULL; \
    double score; \
    double temp; \
    double M_temp; \
    double Ix_temp; \
    double Iy_temp; \
\
    /* Gotoh algorithm with three states */ \
    M_row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!M_row) goto exit; \
    Ix_row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!Ix_row) goto exit; \
    Iy_row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!Iy_row) goto exit; \
\
    /* The top row of the score matrix is a special case, \
     * as there are no previously aligned characters. \
     */ \
    M_row[0] = 0; \
    Ix_row[0] = -DBL_MAX; \
    Iy_row[0] = -DBL_MAX; \
    for (j = 1; j <= nB; j++) { \
        M_row[j] = -DBL_MAX; \
        Ix_row[j] = -DBL_MAX; \
        Iy_row[j] = left_gap_open_A + left_gap_extend_A * (j-1); \
    } \
\
    for (i = 1; i < nA; i++) { \
        M_temp = M_row[0]; \
        Ix_temp = Ix_row[0]; \
        Iy_temp = Iy_row[0]; \
        M_row[0] = -DBL_MAX; \
        Ix_row[0] = left_gap_open_B + left_gap_extend_B * (i-1); \
        Iy_row[0] = -DBL_MAX; \
        kA = sA[i-1]; \
        for (j = 1; j < nB; j++) { \
            kB = sB[j-1]; \
            SELECT_SCORE_GLOBAL(M_temp, \
                                Ix_temp, \
                                Iy_temp); \
            M_temp = M_row[j]; \
            M_row[j] = score + (align_score); \
            SELECT_SCORE_GLOBAL(M_temp + gap_open_B, \
                                Ix_row[j] + gap_extend_B, \
                                Iy_row[j] + gap_open_B); \
            Ix_temp = Ix_row[j]; \
            Ix_row[j] = score; \
            SELECT_SCORE_GLOBAL(M_row[j-1] + gap_open_A, \
                                Ix_row[j-1] + gap_open_A, \
                                Iy_row[j-1] + gap_extend_A); \
            Iy_temp = Iy_row[j]; \
            Iy_row[j] = score; \
        } \
        kB = sB[nB-1]; \
        SELECT_SCORE_GLOBAL(M_temp, \
                            Ix_temp, \
                            Iy_temp); \
        M_temp = M_row[nB]; \
        M_row[nB] = score + (align_score); \
        SELECT_SCORE_GLOBAL(M_temp + right_gap_open_B, \
                            Ix_row[nB] + right_gap_extend_B, \
                            Iy_row[nB] + right_gap_open_B); \
        Ix_row[nB] = score; \
        SELECT_SCORE_GLOBAL(M_row[nB-1] + gap_open_A, \
                            Iy_row[nB-1] + gap_extend_A, \
                            Ix_row[nB-1] + gap_open_A); \
        Iy_row[nB] = score; \
    } \
\
    M_temp = M_row[0]; \
    Ix_temp = Ix_row[0]; \
    Iy_temp = Iy_row[0]; \
    M_row[0] = -DBL_MAX; \
    Ix_row[0] = left_gap_open_B + left_gap_extend_B * (i-1); \
    Iy_row[0] = -DBL_MAX; \
    kA = sA[nA-1]; \
    for (j = 1; j < nB; j++) { \
        kB = sB[j-1]; \
        SELECT_SCORE_GLOBAL(M_temp, \
                            Ix_temp, \
                            Iy_temp); \
        M_temp = M_row[j]; \
        M_row[j] = score + (align_score); \
        SELECT_SCORE_GLOBAL(M_temp + gap_open_B, \
                            Ix_row[j] + gap_extend_B, \
                            Iy_row[j] + gap_open_B); \
        Ix_temp = Ix_row[j]; \
        Ix_row[j] = score; \
        SELECT_SCORE_GLOBAL(M_row[j-1] + right_gap_open_A, \
                            Iy_row[j-1] + right_gap_extend_A, \
                            Ix_row[j-1] + right_gap_open_A); \
        Iy_temp = Iy_row[j]; \
        Iy_row[j] = score; \
    } \
\
    kB = sB[nB-1]; \
    SELECT_SCORE_GLOBAL(M_temp, \
                        Ix_temp, \
                        Iy_temp); \
    M_temp = M_row[nB]; \
    M_row[nB] = score + (align_score); \
    SELECT_SCORE_GLOBAL(M_temp + right_gap_open_B, \
                        Ix_row[nB] + right_gap_extend_B, \
                        Iy_row[nB] + right_gap_open_B); \
    Ix_temp = Ix_row[nB]; \
    Ix_row[nB] = score; \
    SELECT_SCORE_GLOBAL(M_row[nB-1] + right_gap_open_A, \
                        Ix_row[nB-1] + right_gap_open_A, \
                        Iy_row[nB-1] + right_gap_extend_A); \
    Iy_temp = Iy_row[nB]; \
    Iy_row[nB] = score; \
\
    SELECT_SCORE_GLOBAL(M_row[nB], Ix_row[nB], Iy_row[nB]); \
    PyMem_Free(M_row); \
    PyMem_Free(Ix_row); \
    PyMem_Free(Iy_row); \
    return PyFloat_FromDouble(score); \
\
exit: \
    if (M_row) PyMem_Free(M_row); \
    if (Ix_row) PyMem_Free(Ix_row); \
    if (Iy_row) PyMem_Free(Iy_row); \
    return PyErr_NoMemory(); \


#define GOTOH_LOCAL_SCORE(align_score) \
    int i; \
    int j; \
    int kA; \
    int kB; \
    const double gap_open_A = self->target_internal_open_gap_score; \
    const double gap_open_B = self->query_internal_open_gap_score; \
    const double gap_extend_A = self->target_internal_extend_gap_score; \
    const double gap_extend_B = self->query_internal_extend_gap_score; \
    double* M_row = NULL; \
    double* Ix_row = NULL; \
    double* Iy_row = NULL; \
    double score; \
    double temp; \
    double M_temp; \
    double Ix_temp; \
    double Iy_temp; \
    double maximum = 0.0; \
\
    /* Gotoh algorithm with three states */ \
    M_row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!M_row) goto exit; \
    Ix_row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!Ix_row) goto exit; \
    Iy_row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!Iy_row) goto exit; \
 \
    /* The top row of the score matrix is a special case, \
     * as there are no previously aligned characters. \
     */ \
    M_row[0] = 0; \
    Ix_row[0] = -DBL_MAX; \
    Iy_row[0] = -DBL_MAX; \
    for (j = 1; j <= nB; j++) { \
        M_row[j] = -DBL_MAX; \
        Ix_row[j] = -DBL_MAX; \
        Iy_row[j] = 0; \
    } \
    for (i = 1; i < nA; i++) { \
        M_temp = M_row[0]; \
        Ix_temp = Ix_row[0]; \
        Iy_temp = Iy_row[0]; \
        M_row[0] = -DBL_MAX; \
        Ix_row[0] = 0; \
        Iy_row[0] = -DBL_MAX; \
        kA = sA[i-1]; \
        for (j = 1; j < nB; j++) { \
            kB = sB[j-1]; \
            SELECT_SCORE_GOTOH_LOCAL_ALIGN(M_temp, \
                                           Ix_temp, \
                                           Iy_temp, \
                                           (align_score)); \
            M_temp = M_row[j]; \
            M_row[j] = score; \
            SELECT_SCORE_LOCAL3(M_temp + gap_open_B, \
                                Ix_row[j] + gap_extend_B, \
                                Iy_row[j] + gap_open_B); \
            Ix_temp = Ix_row[j]; \
            Ix_row[j] = score; \
            SELECT_SCORE_LOCAL3(M_row[j-1] + gap_open_A, \
                                Ix_row[j-1] + gap_open_A, \
                                Iy_row[j-1] + gap_extend_A); \
            Iy_temp = Iy_row[j]; \
            Iy_row[j] = score; \
        } \
        kB = sB[nB-1]; \
        Ix_row[nB] = 0; \
        Iy_row[nB] = 0; \
        SELECT_SCORE_GOTOH_LOCAL_ALIGN(M_temp, \
                                       Ix_temp, \
                                       Iy_temp, \
                                       (align_score)); \
        M_temp = M_row[nB]; \
        M_row[nB] = score; \
    } \
    M_temp = M_row[0]; \
    Ix_temp = Ix_row[0]; \
    Iy_temp = Iy_row[0]; \
    M_row[0] = -DBL_MAX; \
    Ix_row[0] = 0; \
    Iy_row[0] = -DBL_MAX; \
    kA = sA[nA-1]; \
    for (j = 1; j < nB; j++) { \
        kB = sB[j-1]; \
        SELECT_SCORE_GOTOH_LOCAL_ALIGN(M_temp, \
                                       Ix_temp, \
                                       Iy_temp, \
                                       (align_score)); \
        M_temp = M_row[j]; \
        M_row[j] = score; \
        Ix_temp = Ix_row[j]; \
        Iy_temp = Iy_row[j]; \
        Ix_row[j] = 0; \
        Iy_row[j] = 0; \
    } \
    kB = sB[nB-1]; \
    SELECT_SCORE_GOTOH_LOCAL_ALIGN(M_temp, \
                                   Ix_temp, \
                                   Iy_temp, \
                                   (align_score)); \
    PyMem_Free(M_row); \
    PyMem_Free(Ix_row); \
    PyMem_Free(Iy_row); \
    return PyFloat_FromDouble(maximum); \
exit: \
    if (M_row) PyMem_Free(M_row); \
    if (Ix_row) PyMem_Free(Ix_row); \
    if (Iy_row) PyMem_Free(Iy_row); \
    return PyErr_NoMemory(); \


#define GOTOH_GLOBAL_ALIGN(align_score) \
    int i; \
    int j; \
    int kA; \
    int kB; \
    const double gap_open_A = self->target_internal_open_gap_score; \
    const double gap_open_B = self->query_internal_open_gap_score; \
    const double gap_extend_A = self->target_internal_extend_gap_score; \
    const double gap_extend_B = self->query_internal_extend_gap_score; \
    const double left_gap_open_A = self->target_left_open_gap_score; \
    const double left_gap_open_B = self->query_left_open_gap_score; \
    const double left_gap_extend_A = self->target_left_extend_gap_score; \
    const double left_gap_extend_B = self->query_left_extend_gap_score; \
    const double right_gap_open_A = self->target_right_open_gap_score; \
    const double right_gap_open_B = self->query_right_open_gap_score; \
    const double right_gap_extend_A = self->target_right_extend_gap_score; \
    const double right_gap_extend_B = self->query_right_extend_gap_score; \
    const double epsilon = self->epsilon; \
    TraceGapsGotoh** gaps = NULL; \
    Trace** M = NULL; \
    double* M_row = NULL; \
    double* Ix_row = NULL; \
    double* Iy_row = NULL; \
    double score; \
    int trace; \
    double temp; \
    double M_temp; \
    double Ix_temp; \
    double Iy_temp; \
    PathGenerator* paths; \
\
    /* Gotoh algorithm with three states */ \
    paths = PathGenerator_create_Gotoh(nA, nB, Global); \
    if (!paths) return NULL; \
    M_row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!M_row) goto exit; \
    Ix_row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!Ix_row) goto exit; \
    Iy_row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!Iy_row) goto exit; \
    M = paths->M; \
    gaps = paths->gaps.gotoh; \
 \
    /* Gotoh algorithm with three states */ \
    M_row[0] = 0; \
    Ix_row[0] = -DBL_MAX; \
    Iy_row[0] = -DBL_MAX; \
    for (j = 1; j <= nB; j++) { \
        M_row[j] = -DBL_MAX; \
        Ix_row[j] = -DBL_MAX; \
        Iy_row[j] = left_gap_open_A + left_gap_extend_A * (j-1); \
    } \
    for (i = 1; i < nA; i++) { \
        kA = sA[i-1]; \
        M_temp = M_row[0]; \
        Ix_temp = Ix_row[0]; \
        Iy_temp = Iy_row[0]; \
        M_row[0] = -DBL_MAX; \
        Ix_row[0] = left_gap_open_B + left_gap_extend_B * (i-1); \
        Iy_row[0] = -DBL_MAX; \
        for (j = 1; j < nB; j++) { \
            kB = sB[j-1]; \
            SELECT_TRACE_GOTOH_GLOBAL_ALIGN; \
            M_temp = M_row[j]; \
            M_row[j] = score + (align_score); \
            SELECT_TRACE_GOTOH_GLOBAL_GAP(Ix, \
                                          M_temp + gap_open_B, \
                                          Ix_row[j] + gap_extend_B, \
                                          Iy_row[j] + gap_open_B); \
            Ix_temp = Ix_row[j]; \
            Ix_row[j] = score; \
            SELECT_TRACE_GOTOH_GLOBAL_GAP(Iy, \
                                          M_row[j-1] + gap_open_A, \
                                          Ix_row[j-1] + gap_open_A, \
                                          Iy_row[j-1] + gap_extend_A); \
            Iy_temp = Iy_row[j]; \
            Iy_row[j] = score; \
        } \
        kB = sB[nB-1]; \
        SELECT_TRACE_GOTOH_GLOBAL_ALIGN; \
        M_temp = M_row[nB]; \
        M_row[nB] = score + (align_score); \
        SELECT_TRACE_GOTOH_GLOBAL_GAP(Ix, \
                                      M_temp + right_gap_open_B, \
                                      Ix_row[nB] + right_gap_extend_B, \
                                      Iy_row[nB] + right_gap_open_B); \
        Ix_temp = Ix_row[nB]; \
        Ix_row[nB] = score; \
        SELECT_TRACE_GOTOH_GLOBAL_GAP(Iy, \
                                      M_row[nB-1] + gap_open_A, \
                                      Ix_row[nB-1] + gap_open_A, \
                                      Iy_row[nB-1] + gap_extend_A); \
        Iy_temp = Iy_row[nB]; \
        Iy_row[nB] = score; \
    } \
    kA = sA[nA-1]; \
    M_temp = M_row[0]; \
    Ix_temp = Ix_row[0]; \
    Iy_temp = Iy_row[0]; \
    M_row[0] = -DBL_MAX; \
    Ix_row[0] = left_gap_open_B + left_gap_extend_B * (nA-1); \
    Iy_row[0] = -DBL_MAX; \
    for (j = 1; j < nB; j++) { \
        kB = sB[j-1]; \
        SELECT_TRACE_GOTOH_GLOBAL_ALIGN; \
        M_temp = M_row[j]; \
        M_row[j] = score + (align_score); \
        SELECT_TRACE_GOTOH_GLOBAL_GAP(Ix, \
                                      M_temp + gap_open_B, \
                                      Ix_row[j] + gap_extend_B, \
                                      Iy_row[j] + gap_open_B); \
        Ix_temp = Ix_row[j]; \
        Ix_row[j] = score; \
        SELECT_TRACE_GOTOH_GLOBAL_GAP(Iy, \
                                      M_row[j-1] + right_gap_open_A, \
                                      Ix_row[j-1] + right_gap_open_A, \
                                      Iy_row[j-1] + right_gap_extend_A); \
        Iy_temp = Iy_row[j]; \
        Iy_row[j] = score; \
    } \
    kB = sB[nB-1]; \
    SELECT_TRACE_GOTOH_GLOBAL_ALIGN; \
    M_temp = M_row[j]; \
    M_row[j] = score + (align_score); \
    SELECT_TRACE_GOTOH_GLOBAL_GAP(Ix, \
                                  M_temp + right_gap_open_B, \
                                  Ix_row[j] + right_gap_extend_B, \
                                  Iy_row[j] + right_gap_open_B); \
    Ix_row[nB] = score; \
    SELECT_TRACE_GOTOH_GLOBAL_GAP(Iy, \
                                  M_row[j-1] + right_gap_open_A, \
                                  Ix_row[j-1] + right_gap_open_A, \
                                  Iy_row[j-1] + right_gap_extend_A); \
    Iy_row[nB] = score; \
    M[nA][nB].path = 0; \
 \
    /* traceback */ \
    SELECT_SCORE_GLOBAL(M_row[nB], Ix_row[nB], Iy_row[nB]); \
    if (M_row[nB] < score - epsilon) M[nA][nB].trace = 0; \
    if (Ix_row[nB] < score - epsilon) gaps[nA][nB].Ix = 0; \
    if (Iy_row[nB] < score - epsilon) gaps[nA][nB].Iy = 0; \
    return Py_BuildValue("fN", score, paths); \
exit: \
    Py_DECREF(paths); \
    if (M_row) PyMem_Free(M_row); \
    if (Ix_row) PyMem_Free(Ix_row); \
    if (Iy_row) PyMem_Free(Iy_row); \
    return PyErr_NoMemory(); \


#define GOTOH_LOCAL_ALIGN(align_score) \
    int i; \
    int j; \
    int im = nA; \
    int jm = nB; \
    int kA; \
    int kB; \
    const double gap_open_A = self->target_internal_open_gap_score; \
    const double gap_open_B = self->query_internal_open_gap_score; \
    const double gap_extend_A = self->target_internal_extend_gap_score; \
    const double gap_extend_B = self->query_internal_extend_gap_score; \
    const double epsilon = self->epsilon; \
    Trace** M = NULL; \
    TraceGapsGotoh** gaps = NULL; \
    double* M_row = NULL; \
    double* Ix_row = NULL; \
    double* Iy_row = NULL; \
    double score; \
    int trace; \
    double temp; \
    double M_temp; \
    double Ix_temp; \
    double Iy_temp; \
    double maximum = 0.0; \
    PathGenerator* paths; \
 \
    /* Gotoh algorithm with three states */ \
    paths = PathGenerator_create_Gotoh(nA, nB, Local); \
    if (!paths) return NULL; \
    M = paths->M; \
    gaps = paths->gaps.gotoh; \
    M_row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!M_row) goto exit; \
    Ix_row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!Ix_row) goto exit; \
    Iy_row = PyMem_Malloc((nB+1)*sizeof(double)); \
    if (!Iy_row) goto exit; \
    M_row[0] = 0; \
    Ix_row[0] = -DBL_MAX; \
    Iy_row[0] = -DBL_MAX; \
    for (j = 1; j <= nB; j++) { \
        M_row[j] = 0; \
        Ix_row[j] = -DBL_MAX; \
        Iy_row[j] = -DBL_MAX; \
    } \
    for (i = 1; i < nA; i++) { \
        M_temp = M_row[0]; \
        Ix_temp = Ix_row[0]; \
        Iy_temp = Iy_row[0]; \
        M_row[0] = 0; \
        Ix_row[0] = -DBL_MAX; \
        Iy_row[0] = -DBL_MAX; \
        kA = sA[i-1]; \
        for (j = 1; j < nB; j++) { \
            kB = sB[j-1]; \
            SELECT_TRACE_GOTOH_LOCAL_ALIGN(align_score) \
            M_temp = M_row[j]; \
            M_row[j] = score; \
            SELECT_TRACE_GOTOH_LOCAL_GAP(Ix, \
                                     M_temp + gap_open_B, \
                                     Ix_row[j] + gap_extend_B, \
                                     Iy_row[j] + gap_open_B); \
            Ix_temp = Ix_row[j]; \
            Ix_row[j] = score; \
            SELECT_TRACE_GOTOH_LOCAL_GAP(Iy, \
                                     M_row[j-1] + gap_open_A, \
                                     Ix_row[j-1] + gap_open_A, \
                                     Iy_row[j-1] + gap_extend_A); \
            Iy_temp = Iy_row[j]; \
            Iy_row[j] = score; \
        } \
        kB = sB[nB-1]; \
        SELECT_TRACE_GOTOH_LOCAL_ALIGN(align_score) \
        M_temp = M_row[j]; \
        M_row[j] = score; \
        Ix_temp = Ix_row[nB]; \
        Ix_row[nB] = 0; \
        gaps[i][nB].Ix = 0; \
        Iy_temp = Iy_row[nB]; \
        Iy_row[nB] = 0; \
        gaps[i][nB].Iy = 0; \
    } \
    M_temp = M_row[0]; \
    M_row[0] = 0; \
    M[nA][0].trace = 0; \
    Ix_temp = Ix_row[0]; \
    Ix_row[0] = -DBL_MAX; \
    gaps[nA][0].Ix = 0; \
    gaps[nA][0].Iy = 0; \
    Iy_temp = Iy_row[0]; \
    Iy_row[0] = -DBL_MAX; \
    kA = sA[nA-1]; \
    for (j = 1; j < nB; j++) { \
        kB = sB[j-1]; \
        SELECT_TRACE_GOTOH_LOCAL_ALIGN(align_score) \
        M_temp = M_row[j]; \
        M_row[j] = score; \
        Ix_temp = Ix_row[j]; \
        Ix_row[j] = 0; \
        gaps[nA][j].Ix = 0; \
        Iy_temp = Iy_row[j]; \
        Iy_row[j] = 0; \
        gaps[nA][j].Iy = 0; \
    } \
    kB = sB[nB-1]; \
    SELECT_TRACE_GOTOH_LOCAL_ALIGN(align_score) \
    gaps[nA][nB].Ix = 0; \
    gaps[nA][nB].Iy = 0; \
\
    PyMem_Free(M_row); \
    PyMem_Free(Ix_row); \
    PyMem_Free(Iy_row); \
\
    /* As we don't allow zero-score extensions to alignments, \
     * we need to remove all traces towards an ENDPOINT. \
     * In addition, some points then won't have any path to a STARTPOINT. \
     * Here, use path as a temporary variable to indicate if the point \
     * is reachable from a STARTPOINT. If it is unreachable, remove all \
     * traces from it, and don't allow it to be an ENDPOINT. It may still \
     * be a valid STARTPOINT. */ \
    for (j = 0; j <= nB; j++) M[0][j].path = M_MATRIX; \
    for (i = 1; i <= nA; i++) { \
        M[i][0].path = M_MATRIX; \
        for (j = 1; j <= nB; j++) { \
            /* Remove traces to unreachable points. */ \
            trace = M[i][j].trace; \
            if (!(M[i-1][j-1].path & M_MATRIX)) trace &= ~M_MATRIX; \
            if (!(M[i-1][j-1].path & Ix_MATRIX)) trace &= ~Ix_MATRIX; \
            if (!(M[i-1][j-1].path & Iy_MATRIX)) trace &= ~Iy_MATRIX; \
            if (trace & (STARTPOINT | M_MATRIX | Ix_MATRIX | Iy_MATRIX)) { \
                /* The point is reachable. */ \
                if (trace & ENDPOINT) M[i][j].path = 0; /* no extensions after ENDPOINT */ \
                else M[i][j].path |= M_MATRIX; \
            } \
            else { \
                /* The point is not reachable. Then it is not a STARTPOINT, \
                 * all traces from it can be removed, and it cannot act as \
                 * an ENDPOINT. */ \
                M[i][j].path &= ~M_MATRIX; \
                trace = 0; \
            } \
            M[i][j].trace = trace; \
            trace = gaps[i][j].Ix; \
            if (!(M[i-1][j].path & M_MATRIX)) trace &= ~M_MATRIX; \
            if (!(M[i-1][j].path & Ix_MATRIX)) trace &= ~Ix_MATRIX; \
            if (!(M[i-1][j].path & Iy_MATRIX)) trace &= ~Iy_MATRIX; \
            if (trace & (M_MATRIX | Ix_MATRIX | Iy_MATRIX)) { \
                /* The point is reachable. */ \
                M[i][j].path |= Ix_MATRIX; \
            } \
            else { \
                /* The point is not reachable. Then \
                 * all traces from it can be removed. */ \
                M[i][j].path &= ~Ix_MATRIX; \
                trace = 0; \
            } \
            gaps[i][j].Ix = trace; \
            trace = gaps[i][j].Iy; \
            if (!(M[i][j-1].path & M_MATRIX)) trace &= ~M_MATRIX; \
            if (!(M[i][j-1].path & Ix_MATRIX)) trace &= ~Ix_MATRIX; \
            if (!(M[i][j-1].path & Iy_MATRIX)) trace &= ~Iy_MATRIX; \
            if (trace & (M_MATRIX | Ix_MATRIX | Iy_MATRIX)) { \
                /* The point is reachable. */ \
                M[i][j].path |= Iy_MATRIX; \
            } \
            else { \
                /* The point is not reachable. Then \
                 * all traces from it can be removed. */ \
                M[i][j].path &= ~Iy_MATRIX; \
                trace = 0; \
            } \
            gaps[i][j].Iy = trace; \
        } \
    } \
\
    /* traceback */ \
    if (maximum == 0) M[0][0].path = DONE; \
    else M[0][0].path = 0; \
    return Py_BuildValue("fN", maximum, paths); \
\
exit: \
    Py_DECREF(paths); \
    if (M_row) PyMem_Free(M_row); \
    if (Ix_row) PyMem_Free(Ix_row); \
    if (Iy_row) PyMem_Free(Iy_row); \
    return PyErr_NoMemory(); \


#define WATERMANSMITHBEYER_GLOBAL_SCORE(align_score) \
    int i; \
    int j; \
    int k; \
    int kA; \
    int kB; \
    double** M = NULL; \
    double** Ix = NULL; \
    double** Iy = NULL; \
    double score = 0.0; \
    double gapscore; \
    double temp; \
    int ok = 1; \
\
    /* Waterman-Smith-Beyer algorithm */ \
    M = PyMem_Malloc((nA+1)*sizeof(double*)); \
    if (!M) goto exit; \
    Ix = PyMem_Malloc((nA+1)*sizeof(double*)); \
    if (!Ix) goto exit; \
    Iy = PyMem_Malloc((nA+1)*sizeof(double*)); \
    if (!Iy) goto exit; \
    for (i = 0; i <= nA; i++) { \
        M[i] = PyMem_Malloc((nB+1)*sizeof(double)); \
        if (!M[i]) goto exit; \
        Ix[i] = PyMem_Malloc((nB+1)*sizeof(double)); \
        if (!Ix[i]) goto exit; \
        Iy[i] = PyMem_Malloc((nB+1)*sizeof(double)); \
        if (!Iy[i]) goto exit; \
    } \
\
    /* The top row of the score matrix is a special case, \
     *  as there are no previously aligned characters. \
     */ \
    M[0][0] = 0; \
    Ix[0][0] = -DBL_MAX; \
    Iy[0][0] = -DBL_MAX; \
    for (i = 1; i <= nA; i++) { \
        ok = _call_query_gap_function(self, 0, i, &score); \
        if (!ok) goto exit; \
        M[i][0] = -DBL_MAX; \
        Ix[i][0] = score; \
        Iy[i][0] = -DBL_MAX; \
    } \
    for (j = 1; j <= nB; j++) { \
        ok = _call_target_gap_function(self, 0, j, &score); \
        if (!ok) goto exit; \
        M[0][j] = -DBL_MAX; \
        Ix[0][j] = -DBL_MAX; \
        Iy[0][j] = score; \
    } \
    for (i = 1; i <= nA; i++) { \
        kA = sA[i-1]; \
        for (j = 1; j <= nB; j++) { \
            kB = sB[j-1]; \
            SELECT_SCORE_GLOBAL(M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1]); \
            M[i][j] = score + (align_score); \
            score = -DBL_MAX; \
            for (k = 1; k <= i; k++) { \
                ok = _call_query_gap_function(self, j, k, &gapscore); \
                if (!ok) goto exit; \
                SELECT_SCORE_WATERMAN_SMITH_BEYER(M[i-k][j], Iy[i-k][j]); \
            } \
            Ix[i][j] = score; \
            score = -DBL_MAX; \
            for (k = 1; k <= j; k++) { \
                ok = _call_target_gap_function(self, i, k, &gapscore); \
                if (!ok) goto exit; \
                SELECT_SCORE_WATERMAN_SMITH_BEYER(M[i][j-k], Ix[i][j-k]); \
            } \
            Iy[i][j] = score; \
        } \
    } \
    SELECT_SCORE_GLOBAL(M[nA][nB], Ix[nA][nB], Iy[nA][nB]); \
\
exit: \
    if (M) { \
        /* If M is NULL, then Ix is also NULL. */ \
        if (Ix) { \
            /* If Ix is NULL, then Iy is also NULL. */ \
            if (Iy) { \
                /* If Iy is NULL, then M[i], Ix[i], and Iy[i] are also NULL. */  \
                for (i = 0; i <= nA; i++) { \
                    if (!M[i]) break; \
                    PyMem_Free(M[i]);  \
                    if (!Ix[i]) break; \
                    PyMem_Free(Ix[i]); \
                    if (!Iy[i]) break; \
                    PyMem_Free(Iy[i]); \
                } \
                PyMem_Free(Iy); \
            } \
            PyMem_Free(Ix); \
        } \
        PyMem_Free(M); \
    } \
    if (!ok) return NULL; \
    return PyFloat_FromDouble(score); \


#define WATERMANSMITHBEYER_LOCAL_SCORE(align_score) \
    int i; \
    int j; \
    int gap; \
    int kA; \
    int kB; \
    double** M = NULL; \
    double** Ix = NULL; \
    double** Iy = NULL; \
    double score = 0.0; \
    double gapscore = 0.0; \
    double temp; \
    int ok = 1; \
    double maximum = 0.0; \
    PyObject* result = NULL; \
 \
    /* Waterman-Smith-Beyer algorithm */ \
    M = PyMem_Malloc((nA+1)*sizeof(double*)); \
    if (!M) goto exit; \
    Ix = PyMem_Malloc((nA+1)*sizeof(double*)); \
    if (!Ix) goto exit; \
    Iy = PyMem_Malloc((nA+1)*sizeof(double*)); \
    if (!Iy) goto exit; \
    for (i = 0; i <= nA; i++) { \
        M[i] = PyMem_Malloc((nB+1)*sizeof(double)); \
        if (!M[i]) goto exit; \
        Ix[i] = PyMem_Malloc((nB+1)*sizeof(double)); \
        if (!Ix[i]) goto exit; \
        Iy[i] = PyMem_Malloc((nB+1)*sizeof(double)); \
        if (!Iy[i]) goto exit; \
    } \
    /* The top row of the score matrix is a special case, \
     *  as there are no previously aligned characters. \
     */ \
    M[0][0] = 0; \
    Ix[0][0] = -DBL_MAX; \
    Iy[0][0] = -DBL_MAX; \
    for (i = 1; i <= nA; i++) { \
        M[i][0] = -DBL_MAX; \
        Ix[i][0] = 0; \
        Iy[i][0] = -DBL_MAX; \
    } \
    for (j = 1; j <= nB; j++) { \
        M[0][j] = -DBL_MAX; \
        Ix[0][j] = -DBL_MAX; \
        Iy[0][j] = 0; \
    } \
    for (i = 1; i <= nA; i++) { \
        kA = sA[i-1]; \
        for (j = 1; j <= nB; j++) { \
            kB = sB[j-1]; \
            SELECT_SCORE_GOTOH_LOCAL_ALIGN(M[i-1][j-1], \
                                           Ix[i-1][j-1], \
                                           Iy[i-1][j-1], \
                                           (align_score)); \
            M[i][j] = score; \
            if (i == nA || j == nB) { \
                Ix[i][j] = 0; \
                Iy[i][j] = 0; \
                continue; \
            } \
            score = 0.0; \
            for (gap = 1; gap <= i; gap++) { \
                ok = _call_query_gap_function(self, j, gap, &gapscore); \
                SELECT_SCORE_WATERMAN_SMITH_BEYER(M[i-gap][j], Iy[i-gap][j]); \
                if (!ok) goto exit; \
            } \
            if (score > maximum) maximum = score; \
            Ix[i][j] = score; \
            score = 0.0; \
            for (gap = 1; gap <= j; gap++) { \
                ok = _call_target_gap_function(self, i, gap, &gapscore); \
                if (!ok) goto exit; \
                SELECT_SCORE_WATERMAN_SMITH_BEYER(M[i][j-gap], Ix[i][j-gap]); \
            } \
            if (score > maximum) maximum = score; \
            Iy[i][j] = score; \
        } \
    } \
    SELECT_SCORE_GLOBAL(M[nA][nB], Ix[nA][nB], Iy[nA][nB]); \
    if (score > maximum) maximum = score; \
    result = PyFloat_FromDouble(maximum); \
exit: \
    if (M) { \
        /* If M is NULL, then Ix is also NULL. */ \
        if (Ix) { \
            /* If Ix is NULL, then Iy is also NULL. */ \
            if (Iy) { \
                /* If Iy is NULL, then M[i], Ix[i], and Iy[i] are \
                 * also NULL. */ \
                for (i = 0; i <= nA; i++) { \
                    if (!M[i]) break; \
                    PyMem_Free(M[i]); \
                    if (!Ix[i]) break; \
                    PyMem_Free(Ix[i]); \
                    if (!Iy[i]) break; \
                    PyMem_Free(Iy[i]); \
                } \
                PyMem_Free(Iy); \
            } \
            PyMem_Free(Ix); \
        } \
        PyMem_Free(M); \
    } \
    if (!ok) return NULL; \
    if (!result) return PyErr_NoMemory(); \
    return result; \


#define WATERMANSMITHBEYER_GLOBAL_ALIGN(align_score) \
    int i; \
    int j; \
    int gap; \
    int kA; \
    int kB; \
    const double epsilon = self->epsilon; \
    Trace** M; \
    TraceGapsWatermanSmithBeyer** gaps; \
    double** M_row = NULL; \
    double** Ix_row = NULL; \
    double** Iy_row = NULL; \
    int ng; \
    int nm; \
    double score; \
    double gapscore; \
    double temp; \
    int trace; \
    int* gapM; \
    int* gapXY; \
    int ok = 1; \
    PathGenerator* paths = NULL; \
 \
    /* Waterman-Smith-Beyer algorithm */ \
    paths = PathGenerator_create_WSB(nA, nB, Global); \
    if (!paths) return NULL; \
    M = paths->M; \
    gaps = paths->gaps.waterman_smith_beyer; \
    M_row = PyMem_Malloc((nA+1)*sizeof(double*)); \
    if (!M_row) goto exit; \
    Ix_row = PyMem_Malloc((nA+1)*sizeof(double*)); \
    if (!Ix_row) goto exit; \
    Iy_row = PyMem_Malloc((nA+1)*sizeof(double*)); \
    if (!Iy_row) goto exit; \
    for (i = 0; i <= nA; i++) { \
        M_row[i] = PyMem_Malloc((nB+1)*sizeof(double)); \
        if (!M_row[i]) goto exit; \
        M_row[i][0] = -DBL_MAX; \
        Ix_row[i] = PyMem_Malloc((nB+1)*sizeof(double)); \
        if (!Ix_row[i]) goto exit; \
        Ix_row[i][0] = 0; \
        Iy_row[i] = PyMem_Malloc((nB+1)*sizeof(double)); \
        if (!Iy_row[i]) goto exit; \
        Iy_row[i][0] = -DBL_MAX; \
    } \
    M_row[0][0] = 0; \
    Ix_row[0][0] = -DBL_MAX; \
    for (i = 1; i <= nB; i++) { \
        M_row[0][i] = -DBL_MAX; \
        Ix_row[0][i] = -DBL_MAX; \
        Iy_row[0][i] = 0; \
    } \
    for (i = 1; i <= nA; i++) { \
        ok = _call_query_gap_function(self, 0, i, &score); \
        if (!ok) goto exit; \
        Ix_row[i][0] = score; \
    } \
    for (j = 1; j <= nB; j++) { \
        ok = _call_target_gap_function(self, 0, j, &score); \
        if (!ok) goto exit; \
        Iy_row[0][j] = score; \
    } \
    for (i = 1; i <= nA; i++) { \
        kA = sA[i-1]; \
        for (j = 1; j <= nB; j++) { \
            kB = sB[j-1]; \
            SELECT_TRACE_WATERMAN_SMITH_BEYER_GLOBAL_ALIGN((align_score)); \
            gapM = PyMem_Malloc((i+1)*sizeof(int)); \
            if (!gapM) goto exit; \
            gaps[i][j].MIx = gapM; \
            gapXY = PyMem_Malloc((i+1)*sizeof(int)); \
            if (!gapXY) goto exit; \
            gaps[i][j].IyIx = gapXY; \
            nm = 0; \
            ng = 0; \
            score = -DBL_MAX; \
            for (gap = 1; gap <= i; gap++) { \
                ok = _call_query_gap_function(self, j, gap, &gapscore); \
                if (!ok) goto exit; \
                SELECT_TRACE_WATERMAN_SMITH_BEYER_GAP(M_row[i-gap][j], \
                                                      Iy_row[i-gap][j]); \
            } \
            gapM = PyMem_Realloc(gapM, (nm+1)*sizeof(int)); \
            if (!gapM) goto exit; \
            gaps[i][j].MIx = gapM; \
            gapM[nm] = 0; \
            gapXY = PyMem_Realloc(gapXY, (ng+1)*sizeof(int)); \
            if (!gapXY) goto exit; \
            gapXY[ng] = 0; \
            gaps[i][j].IyIx = gapXY; \
            Ix_row[i][j] = score; \
            gapM = PyMem_Malloc((j+1)*sizeof(int)); \
            if (!gapM) goto exit; \
            gaps[i][j].MIy = gapM; \
            gapXY = PyMem_Malloc((j+1)*sizeof(int)); \
            if (!gapXY) goto exit; \
            gaps[i][j].IxIy = gapXY; \
            nm = 0; \
            ng = 0; \
            score = -DBL_MAX; \
            for (gap = 1; gap <= j; gap++) { \
                ok = _call_target_gap_function(self, i, gap, &gapscore); \
                if (!ok) goto exit; \
                SELECT_TRACE_WATERMAN_SMITH_BEYER_GAP(M_row[i][j-gap], \
                                                      Ix_row[i][j-gap]); \
            } \
            Iy_row[i][j] = score; \
            gapM = PyMem_Realloc(gapM, (nm+1)*sizeof(int)); \
            if (!gapM) goto exit; \
            gaps[i][j].MIy = gapM; \
            gapM[nm] = 0; \
            gapXY = PyMem_Realloc(gapXY, (ng+1)*sizeof(int)); \
            if (!gapXY) goto exit; \
            gaps[i][j].IxIy = gapXY; \
            gapXY[ng] = 0; \
        } \
    } \
    /* traceback */ \
    SELECT_SCORE_GLOBAL(M_row[nA][nB], Ix_row[nA][nB], Iy_row[nA][nB]); \
    M[nA][nB].path = 0; \
    if (M_row[nA][nB] < score - epsilon) M[nA][nB].trace = 0; \
    if (Ix_row[nA][nB] < score - epsilon) { \
        gapM = PyMem_Realloc(gaps[nA][nB].MIx, sizeof(int)); \
        if (!gapM) goto exit; \
        gapM[0] = 0; \
        gaps[nA][nB].MIx = gapM; \
        gapXY = PyMem_Realloc(gaps[nA][nB].IyIx, sizeof(int)); \
        if (!gapXY) goto exit; \
        gapXY[0] = 0; \
        gaps[nA][nB].IyIx = gapXY; \
    } \
    if (Iy_row[nA][nB] < score - epsilon) { \
        gapM = PyMem_Realloc(gaps[nA][nB].MIy, sizeof(int)); \
        if (!gapM) goto exit; \
        gapM[0] = 0; \
        gaps[nA][nB].MIy = gapM; \
        gapXY = PyMem_Realloc(gaps[nA][nB].IxIy, sizeof(int)); \
        if (!gapXY) goto exit; \
        gapXY[0] = 0; \
        gaps[nA][nB].IxIy = gapXY; \
    } \
    for (i = 0; i <= nA; i++) { \
        PyMem_Free(M_row[i]); \
        PyMem_Free(Ix_row[i]); \
        PyMem_Free(Iy_row[i]); \
    } \
    PyMem_Free(M_row); \
    PyMem_Free(Ix_row); \
    PyMem_Free(Iy_row); \
    return Py_BuildValue("fN", score, paths); \
\
exit: \
    if (ok) /* otherwise, an exception was already set */ \
        PyErr_SetNone(PyExc_MemoryError); \
    Py_DECREF(paths); \
    if (M_row) { \
        /* If M is NULL, then Ix is also NULL. */ \
        if (Ix_row) { \
            /* If Ix is NULL, then Iy is also NULL. */ \
            if (Iy_row) { \
                /* If Iy is NULL, then M[i], Ix[i], and Iy[i] are also NULL. */ \
                for (i = 0; i <= nA; i++) { \
                    if (!M_row[i]) break; \
                    PyMem_Free(M_row[i]); \
                    if (!Ix_row[i]) break; \
                    PyMem_Free(Ix_row[i]); \
                    if (!Iy_row[i]) break; \
                    PyMem_Free(Iy_row[i]); \
                } \
                PyMem_Free(Iy_row); \
            } \
            PyMem_Free(Ix_row); \
        } \
        PyMem_Free(M_row); \
    } \
    return NULL; \


#define WATERMANSMITHBEYER_LOCAL_ALIGN(align_score) \
    int i; \
    int j; \
    int im = nA; \
    int jm = nB; \
    int gap; \
    int kA; \
    int kB; \
    const double epsilon = self->epsilon; \
    Trace** M = NULL; \
    TraceGapsWatermanSmithBeyer** gaps; \
    double** M_row; \
    double** Ix_row = NULL; \
    double** Iy_row = NULL; \
    double score; \
    double gapscore; \
    double temp; \
    int trace; \
    int* gapM; \
    int* gapXY; \
    int nm; \
    int ng; \
    int ok = 1; \
    double maximum = 0; \
    PathGenerator* paths = NULL; \
 \
    /* Waterman-Smith-Beyer algorithm */ \
    paths = PathGenerator_create_WSB(nA, nB, Local); \
    if (!paths) return NULL; \
    M = paths->M; \
    gaps = paths->gaps.waterman_smith_beyer; \
    M_row = PyMem_Malloc((nA+1)*sizeof(double*)); \
    if (!M_row) goto exit; \
    Ix_row = PyMem_Malloc((nA+1)*sizeof(double*)); \
    if (!Ix_row) goto exit; \
    Iy_row = PyMem_Malloc((nA+1)*sizeof(double*)); \
    if (!Iy_row) goto exit; \
    for (i = 0; i <= nA; i++) { \
        M_row[i] = PyMem_Malloc((nB+1)*sizeof(double)); \
        if (!M_row[i]) goto exit; \
        M_row[i][0] = 0; \
        Ix_row[i] = PyMem_Malloc((nB+1)*sizeof(double)); \
        if (!Ix_row[i]) goto exit; \
        Ix_row[i][0] = -DBL_MAX; \
        Iy_row[i] = PyMem_Malloc((nB+1)*sizeof(double)); \
        if (!Iy_row[i]) goto exit; \
        Iy_row[i][0] = -DBL_MAX; \
    } \
    for (i = 1; i <= nB; i++) { \
        M_row[0][i] = 0; \
        Ix_row[0][i] = -DBL_MAX; \
        Iy_row[0][i] = -DBL_MAX; \
    } \
    for (i = 1; i <= nA; i++) { \
        kA = sA[i-1]; \
        for (j = 1; j <= nB; j++) { \
            kB = sB[j-1]; \
            nm = 0; \
            ng = 0; \
            SELECT_TRACE_WATERMAN_SMITH_BEYER_ALIGN( \
                                           M_row[i-1][j-1], \
                                           Ix_row[i-1][j-1], \
                                           Iy_row[i-1][j-1], \
                                           (align_score)); \
            M[i][j].path = 0; \
            if (i == nA || j == nB) { \
                Ix_row[i][j] = score; \
                gaps[i][j].MIx = NULL; \
                gaps[i][j].IyIx = NULL; \
                gaps[i][j].MIy = NULL; \
                gaps[i][j].IxIy = NULL; \
                Iy_row[i][j] = score; \
                continue; \
            } \
            gapM = PyMem_Malloc((i+1)*sizeof(int)); \
            if (!gapM) goto exit; \
            gaps[i][j].MIx = gapM; \
            gapXY = PyMem_Malloc((i+1)*sizeof(int)); \
            if (!gapXY) goto exit; \
            gaps[i][j].IyIx = gapXY; \
            score = -DBL_MAX; \
            for (gap = 1; gap <= i; gap++) { \
                ok = _call_query_gap_function(self, j, gap, &gapscore); \
                if (!ok) goto exit; \
                SELECT_TRACE_WATERMAN_SMITH_BEYER_GAP(M_row[i-gap][j], \
                                                      Iy_row[i-gap][j]); \
            } \
            if (score < epsilon) { \
                score = -DBL_MAX; \
                nm = 0; \
                ng = 0; \
            } \
            else if (score > maximum) maximum = score; \
            gapM[nm] = 0; \
            gapXY[ng] = 0; \
            Ix_row[i][j] = score; \
            M[i][j].path = 0; \
            gapM = PyMem_Realloc(gapM, (nm+1)*sizeof(int)); \
            if (!gapM) goto exit; \
            gaps[i][j].MIx = gapM; \
            gapM[nm] = 0; \
            gapXY = PyMem_Realloc(gapXY, (ng+1)*sizeof(int)); \
            if (!gapXY) goto exit; \
            gaps[i][j].IyIx = gapXY; \
            gapXY[ng] = 0; \
            gapM = PyMem_Malloc((j+1)*sizeof(int)); \
            if (!gapM) goto exit; \
            gaps[i][j].MIy = gapM; \
            gapXY = PyMem_Malloc((j+1)*sizeof(int)); \
            if (!gapXY) goto exit; \
            gaps[i][j].IxIy = gapXY; \
            nm = 0; \
            ng = 0; \
            score = -DBL_MAX; \
            gapM[0] = 0; \
            for (gap = 1; gap <= j; gap++) { \
                ok = _call_target_gap_function(self, i, gap, &gapscore); \
                if (!ok) goto exit; \
                SELECT_TRACE_WATERMAN_SMITH_BEYER_GAP(M_row[i][j-gap], \
                                                      Ix_row[i][j-gap]); \
            } \
            if (score < epsilon) { \
                score = -DBL_MAX; \
                nm = 0; \
                ng = 0; \
            } \
            else if (score > maximum) maximum = score; \
            gapM = PyMem_Realloc(gapM, (nm+1)*sizeof(int)); \
            if (!gapM) goto exit; \
            gaps[i][j].MIy = gapM; \
            gapXY = PyMem_Realloc(gapXY, (ng+1)*sizeof(int)); \
            if (!gapXY) goto exit; \
            gaps[i][j].IxIy = gapXY; \
            gapM[nm] = 0; \
            gapXY[ng] = 0; \
            Iy_row[i][j] = score; \
            M[i][j].path = 0; \
        } \
    } \
    for (i = 0; i <= nA; i++) PyMem_Free(M_row[i]); \
    PyMem_Free(M_row); \
    for (i = 0; i <= nA; i++) PyMem_Free(Ix_row[i]); \
    PyMem_Free(Ix_row); \
    for (i = 0; i <= nA; i++) PyMem_Free(Iy_row[i]); \
    PyMem_Free(Iy_row); \
\
    /* As we don't allow zero-score extensions to alignments, \
     * we need to remove all traces towards an ENDPOINT. \
     * In addition, some points then won't have any path to a STARTPOINT. \
     * Here, use path as a temporary variable to indicate if the point \
     * is reachable from a STARTPOINT. If it is unreachable, remove all \
     * traces from it, and don't allow it to be an ENDPOINT. It may still \
     * be a valid STARTPOINT. */ \
    for (j = 0; j <= nB; j++) M[0][j].path = M_MATRIX; \
    for (i = 1; i <= nA; i++) { \
        M[i][0].path = M_MATRIX; \
        for (j = 1; j <= nB; j++) { \
            /* Remove traces to unreachable points. */ \
            trace = M[i][j].trace; \
            if (!(M[i-1][j-1].path & M_MATRIX)) trace &= ~M_MATRIX; \
            if (!(M[i-1][j-1].path & Ix_MATRIX)) trace &= ~Ix_MATRIX; \
            if (!(M[i-1][j-1].path & Iy_MATRIX)) trace &= ~Iy_MATRIX; \
            if (trace & (STARTPOINT | M_MATRIX | Ix_MATRIX | Iy_MATRIX)) { \
                /* The point is reachable. */ \
                if (trace & ENDPOINT) M[i][j].path = 0; /* no extensions after ENDPOINT */ \
                else M[i][j].path |= M_MATRIX; \
            } \
            else { \
                /* The point is not reachable. Then it is not a STARTPOINT, \
                 * all traces from it can be removed, and it cannot act as \
                 * an ENDPOINT. */ \
                M[i][j].path &= ~M_MATRIX; \
                trace = 0; \
            } \
            M[i][j].trace = trace; \
            if (i == nA || j == nB) continue; \
            gapM = gaps[i][j].MIx; \
            gapXY = gaps[i][j].IyIx; \
            nm = 0; \
            ng = 0; \
            for (im = 0; (gap = gapM[im]); im++) \
                if (M[i-gap][j].path & M_MATRIX) gapM[nm++] = gap; \
            gapM = PyMem_Realloc(gapM, (nm+1)*sizeof(int)); \
            if (!gapM) goto exit; \
            gapM[nm] = 0; \
            gaps[i][j].MIx = gapM; \
            for (im = 0; (gap = gapXY[im]); im++) \
                if (M[i-gap][j].path & Iy_MATRIX) gapXY[ng++] = gap; \
            gapXY = PyMem_Realloc(gapXY, (ng+1)*sizeof(int)); \
            if (!gapXY) goto exit; \
            gapXY[ng] = 0; \
            gaps[i][j].IyIx = gapXY; \
            if (nm==0 && ng==0) M[i][j].path &= ~Ix_MATRIX; /* not reachable */ \
            else M[i][j].path |= Ix_MATRIX; /* reachable */ \
            gapM = gaps[i][j].MIy; \
            gapXY = gaps[i][j].IxIy; \
            nm = 0; \
            ng = 0; \
            for (im = 0; (gap = gapM[im]); im++) \
                if (M[i][j-gap].path & M_MATRIX) gapM[nm++] = gap; \
            gapM = PyMem_Realloc(gapM, (nm+1)*sizeof(int)); \
            if (!gapM) goto exit; \
            gapM[nm] = 0; \
            gaps[i][j].MIy = gapM; \
            for (im = 0; (gap = gapXY[im]); im++) \
                if (M[i][j-gap].path & Ix_MATRIX) gapXY[ng++] = gap; \
            gapXY = PyMem_Realloc(gapXY, (ng+1)*sizeof(int)); \
            if (!gapXY) goto exit; \
            gapXY[ng] = 0; \
            gaps[i][j].IxIy = gapXY; \
            if (nm==0 && ng==0) M[i][j].path &= ~Iy_MATRIX; /* not reachable */ \
            else M[i][j].path |= Iy_MATRIX; /* reachable */ \
        } \
    } \
    /* traceback */ \
    if (maximum == 0) M[0][0].path = DONE; \
    else M[0][0].path = 0; \
    return Py_BuildValue("fN", maximum, paths); \
\
exit: \
    if (ok) /* otherwise, an exception was already set */ \
        PyErr_SetNone(PyExc_MemoryError); \
    Py_DECREF(paths); \
    if (M_row) { \
        /* If M is NULL, then Ix is also NULL. */ \
        if (Ix_row) { \
            /* If Ix is NULL, then Iy is also NULL. */ \
            if (Iy_row) { \
                /* If Iy is NULL, then M[i], Ix[i], and Iy[i] are also NULL. */ \
                for (i = 0; i <= nA; i++) { \
                    if (!M_row[i]) break; \
                    PyMem_Free(M_row[i]); \
                    if (!Ix_row[i]) break; \
                    PyMem_Free(Ix_row[i]); \
                    if (!Iy_row[i]) break; \
                    PyMem_Free(Iy_row[i]); \
                } \
                PyMem_Free(Iy_row); \
            } \
            PyMem_Free(Ix_row); \
        } \
        PyMem_Free(M_row); \
    } \
    return NULL; \


/* -------------- allocation & deallocation ------------- */

static PathGenerator*
PathGenerator_create_NWSW(Py_ssize_t nA, Py_ssize_t nB, Mode mode)
{
    int i;
    unsigned char trace = 0;
    Trace** M;
    PathGenerator* paths;

    paths = (PathGenerator*)PyType_GenericAlloc(&PathGenerator_Type, 0);
    if (!paths) return NULL;

    paths->iA = 0;
    paths->iB = 0;
    paths->nA = nA;
    paths->nB = nB;
    paths->M = NULL;
    paths->gaps.gotoh = NULL;
    paths->gaps.waterman_smith_beyer = NULL;
    paths->algorithm = NeedlemanWunschSmithWaterman;
    paths->mode = mode;
    paths->length = 0;

    M = PyMem_Malloc((nA+1)*sizeof(Trace*));
    paths->M = M;
    if (!M) goto exit;
    switch (mode) {
        case Global: trace = VERTICAL; break;
        case Local: trace = STARTPOINT; break;
    }
    for (i = 0; i <= nA; i++) {
        M[i] = PyMem_Malloc((nB+1)*sizeof(Trace));
        if (!M[i]) goto exit;
        M[i][0].trace = trace;
    }
    if (mode == Global) {
        M[0][0].trace = 0;
        trace = HORIZONTAL;
    }
    for (i = 1; i <= nB; i++) M[0][i].trace = trace;
    M[0][0].path = 0;
    return paths;
exit:
    Py_DECREF(paths);
    PyErr_SetNone(PyExc_MemoryError);
    return NULL;
}

static PathGenerator*
PathGenerator_create_Gotoh(Py_ssize_t nA, Py_ssize_t nB, Mode mode)
{
    int i;
    unsigned char trace;
    Trace** M;
    TraceGapsGotoh** gaps;
    PathGenerator* paths;

    switch (mode) {
        case Global: trace = 0; break;
        case Local: trace = STARTPOINT; break;
        default:
            /* Should not happen, but the compiler has no way of knowing that,
             * as the enum Mode does not restrict the value of mode, which can
             * be any integer. Include default: here to prevent compiler
             * warnings.
             */
            PyErr_Format(PyExc_RuntimeError,
                         "mode has unexpected value %d", mode);
            return NULL;
    }

    paths = (PathGenerator*)PyType_GenericAlloc(&PathGenerator_Type, 0);
    if (!paths) return NULL;

    paths->iA = 0;
    paths->iB = 0;
    paths->nA = nA;
    paths->nB = nB;
    paths->M = NULL;
    paths->gaps.gotoh = NULL;
    paths->algorithm = Gotoh;
    paths->mode = mode;
    paths->length = 0;

    M = PyMem_Malloc((nA+1)*sizeof(Trace*));
    if (!M) goto exit;
    paths->M = M;
    for (i = 0; i <= nA; i++) {
        M[i] = PyMem_Malloc((nB+1)*sizeof(Trace));
        if (!M[i]) goto exit;
        M[i][0].trace = trace;
    }
    gaps = PyMem_Malloc((nA+1)*sizeof(TraceGapsGotoh*));
    if (!gaps) goto exit;
    paths->gaps.gotoh = gaps;
    for (i = 0; i <= nA; i++) {
        gaps[i] = PyMem_Malloc((nB+1)*sizeof(TraceGapsGotoh));
        if (!gaps[i]) goto exit;
    }

    gaps[0][0].Ix = 0;
    gaps[0][0].Iy = 0;
    if (mode == Global) {
        for (i = 1; i <= nA; i++) {
            gaps[i][0].Ix = Ix_MATRIX;
            gaps[i][0].Iy = 0;
        }
        gaps[1][0].Ix = M_MATRIX;
        for (i = 1; i <= nB; i++) {
            M[0][i].trace = 0;
            gaps[0][i].Ix = 0;
            gaps[0][i].Iy = Iy_MATRIX;
        }
        gaps[0][1].Iy = M_MATRIX;
    }
    else if (mode == Local) {
        for (i = 1; i < nA; i++) {
            gaps[i][0].Ix = 0;
            gaps[i][0].Iy = 0;
        }
        for (i = 1; i <= nB; i++) {
            M[0][i].trace = trace;
            gaps[0][i].Ix = 0;
            gaps[0][i].Iy = 0;
        }
    }
    M[0][0].path = 0;

    return paths;
exit:
    Py_DECREF(paths);
    PyErr_SetNone(PyExc_MemoryError);
    return NULL;
}

static PathGenerator*
PathGenerator_create_WSB(Py_ssize_t nA, Py_ssize_t nB, Mode mode)
{
    int i, j;
    int* trace;
    Trace** M = NULL;
    TraceGapsWatermanSmithBeyer** gaps = NULL;
    PathGenerator* paths;

    paths = (PathGenerator*)PyType_GenericAlloc(&PathGenerator_Type, 0);
    if (!paths) return NULL;

    paths->iA = 0;
    paths->iB = 0;
    paths->nA = nA;
    paths->nB = nB;
    paths->M = NULL;
    paths->gaps.waterman_smith_beyer = NULL;
    paths->algorithm = WatermanSmithBeyer;
    paths->mode = mode;
    paths->length = 0;

    M = PyMem_Malloc((nA+1)*sizeof(Trace*));
    if (!M) goto exit;
    paths->M = M;
    for (i = 0; i <= nA; i++) {
        M[i] = PyMem_Malloc((nB+1)*sizeof(Trace));
        if (!M[i]) goto exit;
    }
    gaps = PyMem_Malloc((nA+1)*sizeof(TraceGapsWatermanSmithBeyer*));
    if (!gaps) goto exit;
    paths->gaps.waterman_smith_beyer = gaps;
    for (i = 0; i <= nA; i++) gaps[i] = NULL;
    for (i = 0; i <= nA; i++) {
        gaps[i] = PyMem_Malloc((nB+1)*sizeof(TraceGapsWatermanSmithBeyer));
        if (!gaps[i]) goto exit;
        for (j = 0; j <= nB; j++) {
            gaps[i][j].MIx = NULL;
            gaps[i][j].IyIx = NULL;
            gaps[i][j].MIy = NULL;
            gaps[i][j].IxIy = NULL;
        }
        M[i][0].path = 0;
        switch (mode) {
            case Global:
                M[i][0].trace = 0;
                trace = PyMem_Malloc(2*sizeof(int));
                if (!trace) goto exit;
                gaps[i][0].MIx = trace;
                trace[0] = i;
                trace[1] = 0;
                trace = PyMem_Malloc(sizeof(int));
                if (!trace) goto exit;
                gaps[i][0].IyIx = trace;
                trace[0] = 0;
                break;
            case Local:
                M[i][0].trace = STARTPOINT;
                break;
        }
    }
    for (i = 1; i <= nB; i++) {
        switch (mode) {
            case Global:
                M[0][i].trace = 0;
                trace = PyMem_Malloc(2*sizeof(int));
                if (!trace) goto exit;
                gaps[0][i].MIy = trace;
                trace[0] = i;
                trace[1] = 0;
                trace = PyMem_Malloc(sizeof(int));
                if (!trace) goto exit;
                gaps[0][i].IxIy = trace;
                trace[0] = 0;
                break;
            case Local:
                M[0][i].trace = STARTPOINT;
                break;
        }
    }
    M[0][0].path = 0;
    return paths;
exit:
    Py_DECREF(paths);
    PyErr_SetNone(PyExc_MemoryError);
    return NULL;
}

/* ----------------- alignment algorithms ----------------- */

#define MATRIX_SCORE scores[kA*n+kB]
#define COMPARE_SCORE (kA < 0 || kB < 0) ? 0 : (kA == kB) ? match : mismatch


static PyObject*
Aligner_needlemanwunsch_score_compare(Aligner* self,
                                      const int* sA, Py_ssize_t nA,
                                      const int* sB, Py_ssize_t nB)
{
    const double match = self->match;
    const double mismatch = self->mismatch;
    NEEDLEMANWUNSCH_SCORE(COMPARE_SCORE);
}

static PyObject*
Aligner_needlemanwunsch_score_matrix(Aligner* self,
                                     const int* sA, Py_ssize_t nA,
                                     const int* sB, Py_ssize_t nB)
{
    const Py_ssize_t n = self->substitution_matrix.shape[0];
    const double* scores = self->substitution_matrix.buf;
    NEEDLEMANWUNSCH_SCORE(MATRIX_SCORE);
}

static PyObject*
Aligner_smithwaterman_score_compare(Aligner* self,
                                    const int* sA, Py_ssize_t nA,
                                    const int* sB, Py_ssize_t nB)
{
    const double match = self->match;
    const double mismatch = self->mismatch;
    SMITHWATERMAN_SCORE(COMPARE_SCORE);
}

static PyObject*
Aligner_smithwaterman_score_matrix(Aligner* self,
                                   const int* sA, Py_ssize_t nA,
                                   const int* sB, Py_ssize_t nB)
{
    const Py_ssize_t n = self->substitution_matrix.shape[0];
    const double* scores = self->substitution_matrix.buf;
    SMITHWATERMAN_SCORE(MATRIX_SCORE);
}

static PyObject*
Aligner_needlemanwunsch_align_compare(Aligner* self,
                                      const int* sA, Py_ssize_t nA,
                                      const int* sB, Py_ssize_t nB)
{
    const double match = self->match;
    const double mismatch = self->mismatch;
    NEEDLEMANWUNSCH_ALIGN(COMPARE_SCORE);
}

static PyObject*
Aligner_needlemanwunsch_align_matrix(Aligner* self,
                                     const int* sA, Py_ssize_t nA,
                                     const int* sB, Py_ssize_t nB)
{
    const Py_ssize_t n = self->substitution_matrix.shape[0];
    const double* scores = self->substitution_matrix.buf;
    NEEDLEMANWUNSCH_ALIGN(MATRIX_SCORE);
}

static PyObject*
Aligner_smithwaterman_align_compare(Aligner* self,
                                    const int* sA, Py_ssize_t nA,
                                    const int* sB, Py_ssize_t nB)
{
    const double match = self->match;
    const double mismatch = self->mismatch;
    SMITHWATERMAN_ALIGN(COMPARE_SCORE);
}

static PyObject*
Aligner_smithwaterman_align_matrix(Aligner* self,
                                   const int* sA, Py_ssize_t nA,
                                   const int* sB, Py_ssize_t nB)
{
    const Py_ssize_t n = self->substitution_matrix.shape[0];
    const double* scores = self->substitution_matrix.buf;
    SMITHWATERMAN_ALIGN(MATRIX_SCORE);
}

static PyObject*
Aligner_gotoh_global_score_compare(Aligner* self,
                                   const int* sA, Py_ssize_t nA,
                                   const int* sB, Py_ssize_t nB)
{
    const double match = self->match;
    const double mismatch = self->mismatch;
    GOTOH_GLOBAL_SCORE(COMPARE_SCORE);
}

static PyObject*
Aligner_gotoh_global_score_matrix(Aligner* self,
                                  const int* sA, Py_ssize_t nA,
                                  const int* sB, Py_ssize_t nB)
{
    const Py_ssize_t n = self->substitution_matrix.shape[0];
    const double* scores = self->substitution_matrix.buf;
    GOTOH_GLOBAL_SCORE(MATRIX_SCORE);
}

static PyObject*
Aligner_gotoh_local_score_compare(Aligner* self,
                                  const int* sA, Py_ssize_t nA,
                                  const int* sB, Py_ssize_t nB)
{
    const double match = self->match;
    const double mismatch = self->mismatch;
    GOTOH_LOCAL_SCORE(COMPARE_SCORE);
}

static PyObject*
Aligner_gotoh_local_score_matrix(Aligner* self,
                                 const int* sA, Py_ssize_t nA,
                                 const int* sB, Py_ssize_t nB)
{
    const Py_ssize_t n = self->substitution_matrix.shape[0];
    const double* scores = self->substitution_matrix.buf;
    GOTOH_LOCAL_SCORE(MATRIX_SCORE);
}

static PyObject*
Aligner_gotoh_global_align_compare(Aligner* self,
                                   const int* sA, Py_ssize_t nA,
                                   const int* sB, Py_ssize_t nB)
{
    const double match = self->match;
    const double mismatch = self->mismatch;
    GOTOH_GLOBAL_ALIGN(COMPARE_SCORE);
}

static PyObject*
Aligner_gotoh_global_align_matrix(Aligner* self,
                                  const int* sA, Py_ssize_t nA,
                                  const int* sB, Py_ssize_t nB)
{
    const Py_ssize_t n = self->substitution_matrix.shape[0];
    const double* scores = self->substitution_matrix.buf;
    GOTOH_GLOBAL_ALIGN(MATRIX_SCORE);
}

static PyObject*
Aligner_gotoh_local_align_compare(Aligner* self,
                                  const int* sA, Py_ssize_t nA,
                                  const int* sB, Py_ssize_t nB)
{
    const double match = self->match;
    const double mismatch = self->mismatch;
    GOTOH_LOCAL_ALIGN(COMPARE_SCORE);
}

static PyObject*
Aligner_gotoh_local_align_matrix(Aligner* self,
                                 const int* sA, Py_ssize_t nA,
                                 const int* sB, Py_ssize_t nB)
{
    const Py_ssize_t n = self->substitution_matrix.shape[0];
    const double* scores = self->substitution_matrix.buf;
    GOTOH_LOCAL_ALIGN(MATRIX_SCORE);
}

static int
_call_query_gap_function(Aligner* aligner, int i, int j, double* score)
{
    double value;
    PyObject* result;
    PyObject* function = aligner->query_gap_function;
    if (!function)
        value = aligner->query_internal_open_gap_score
              + (j-1) * aligner->query_internal_extend_gap_score;
    else {
        result = PyObject_CallFunction(function, "ii", i, j);
        if (result == NULL) return 0;
        value = PyFloat_AsDouble(result);
        Py_DECREF(result);
        if (value == -1.0 && PyErr_Occurred()) return 0;
    }
    *score = value;
    return 1;
}

static int
_call_target_gap_function(Aligner* aligner, int i, int j, double* score)
{
    double value;
    PyObject* result;
    PyObject* function = aligner->target_gap_function;
    if (!function)
        value = aligner->target_internal_open_gap_score
              + (j-1) * aligner->target_internal_extend_gap_score;
    else {
        result = PyObject_CallFunction(function, "ii", i, j);
        if (result == NULL) return 0;
        value = PyFloat_AsDouble(result);
        Py_DECREF(result);
        if (value == -1.0 && PyErr_Occurred()) return 0;
    }
    *score = value;
    return 1;
}

static PyObject*
Aligner_watermansmithbeyer_global_score_compare(Aligner* self,
                                                const int* sA, Py_ssize_t nA,
                                                const int* sB, Py_ssize_t nB)
{
    const double match = self->match;
    const double mismatch = self->mismatch;
    WATERMANSMITHBEYER_GLOBAL_SCORE(COMPARE_SCORE);
}

static PyObject*
Aligner_watermansmithbeyer_global_score_matrix(Aligner* self,
                                               const int* sA, Py_ssize_t nA,
                                               const int* sB, Py_ssize_t nB)
{
    const Py_ssize_t n = self->substitution_matrix.shape[0];
    const double* scores = self->substitution_matrix.buf;
    WATERMANSMITHBEYER_GLOBAL_SCORE(MATRIX_SCORE);
}

static PyObject*
Aligner_watermansmithbeyer_local_score_compare(Aligner* self,
                                               const int* sA, Py_ssize_t nA,
                                               const int* sB, Py_ssize_t nB)
{
    const double match = self->match;
    const double mismatch = self->mismatch;
    WATERMANSMITHBEYER_LOCAL_SCORE(COMPARE_SCORE);
}

static PyObject*
Aligner_watermansmithbeyer_local_score_matrix(Aligner* self,
                                              const int* sA, Py_ssize_t nA,
                                              const int* sB, Py_ssize_t nB)
{
    const Py_ssize_t n = self->substitution_matrix.shape[0];
    const double* scores = self->substitution_matrix.buf;
    WATERMANSMITHBEYER_LOCAL_SCORE(MATRIX_SCORE);
}

static PyObject*
Aligner_watermansmithbeyer_global_align_compare(Aligner* self,
                                                const int* sA, Py_ssize_t nA,
                                                const int* sB, Py_ssize_t nB)
{
    const double match = self->match;
    const double mismatch = self->mismatch;
    WATERMANSMITHBEYER_GLOBAL_ALIGN(COMPARE_SCORE);
}

static PyObject*
Aligner_watermansmithbeyer_global_align_matrix(Aligner* self,
                                               const int* sA, Py_ssize_t nA,
                                               const int* sB, Py_ssize_t nB)
{
    const Py_ssize_t n = self->substitution_matrix.shape[0];
    const double* scores = self->substitution_matrix.buf;
    WATERMANSMITHBEYER_GLOBAL_ALIGN(MATRIX_SCORE);
}

static PyObject*
Aligner_watermansmithbeyer_local_align_compare(Aligner* self,
                                               const int* sA, Py_ssize_t nA,
                                               const int* sB, Py_ssize_t nB)
{
    const double match = self->match;
    const double mismatch = self->mismatch;
    WATERMANSMITHBEYER_LOCAL_ALIGN(COMPARE_SCORE);
}

static PyObject*
Aligner_watermansmithbeyer_local_align_matrix(Aligner* self,
                                              const int* sA, Py_ssize_t nA,
                                              const int* sB, Py_ssize_t nB)
{
    const Py_ssize_t n = self->substitution_matrix.shape[0];
    const double* scores = self->substitution_matrix.buf;
    WATERMANSMITHBEYER_LOCAL_ALIGN(MATRIX_SCORE);
}

static int*
convert_sequence_to_ints(const signed char mapping[], Py_ssize_t n,
                         const char s[])
{
    char c;
    Py_ssize_t i;
    int index;
    int* indices;
    if (n == 0) {
        PyErr_SetString(PyExc_ValueError, "sequence has zero length");
        return NULL;
    }
    indices = PyMem_Malloc(n*sizeof(int));
    if (!indices) {
        PyErr_NoMemory();
        return NULL;
    }
    for (i = 0; i < n; i++) {
        c = s[i];
        index = mapping[(int)c];
        if (index == MISSING_LETTER) {
            PyErr_SetString(PyExc_ValueError,
                "sequence contains letters not in the alphabet");
            PyMem_Free(indices);
            return NULL;
        }
        indices[i] = index;
    }
    return indices;
}

static int
convert_objects_to_ints(Py_buffer* view, PyObject* alphabet, PyObject* sequence)
{
    Py_ssize_t i, j;
    Py_ssize_t n;
    Py_ssize_t m;
    int* indices = NULL;
    PyObject *obj1, *obj2;
    int equal;
    view->buf = NULL;
    sequence = PySequence_Fast(sequence,
                               "argument should support the sequence protocol");
    if (!sequence) return 0;
    alphabet = PySequence_Fast(alphabet,
                               "alphabet should support the sequence protocol");
    if (!alphabet) goto exit;
    n = PySequence_Size(sequence);
    m = PySequence_Size(alphabet);
    indices = PyMem_Malloc(n*sizeof(int));
    if (!indices) {
        PyErr_NoMemory();
        goto exit;
    }
    for (i = 0; i < n; i++) {
        obj1 = PySequence_Fast_GET_ITEM(sequence, i);
        for (j = 0; j < m; j++) {
            obj2 = PySequence_Fast_GET_ITEM(alphabet, j);
            equal = PyObject_RichCompareBool(obj1, obj2, Py_EQ);
            if (equal == 1) /* obj1 == obj2 */ {
                indices[i] = j;
                break;
            }
            else if (equal == -1) /* error */ {
                PyMem_Del(indices);
                goto exit;
            }
            /* else (equal == 0) continue; */ /* not equal */
        }
        if (j == m) {
            PyErr_SetString(PyExc_ValueError, "failed to find object in alphabet");
            goto exit;
        }
    }
    view->buf = indices;
    view->itemsize = 1;
    view->len = n;
exit:
    if (sequence) {
        Py_DECREF(sequence);
        if (alphabet) Py_DECREF(alphabet);
    }
    if (view->buf) return 1;
    return 0;
}

static int
sequence_converter(PyObject* argument, void* pointer)
{
    Py_buffer* view = pointer;
    Py_ssize_t i;
    Py_ssize_t n;
    int index;
    int* indices;
    const char* s;
    const int flag = PyBUF_FORMAT | PyBUF_C_CONTIGUOUS;
    Aligner* aligner;
    signed char* mapping;

    if (argument == NULL) {
        if (view->obj) PyBuffer_Release(view);
        else {
            indices = view->buf;
            PyMem_Free(indices);
        }
        return 1;
    }

    aligner = (Aligner*)view->obj;
    view->obj = NULL;

    mapping = aligner->mapping;
    if (*mapping == UNMAPPED) {
        if (!convert_objects_to_ints(view, aligner->alphabet, argument)) return 0;
        return Py_CLEANUP_SUPPORTED;
    }
    
    s = PyUnicode_AsUTF8AndSize(argument, &n);
    if (s) {
        indices = convert_sequence_to_ints(aligner->mapping, n, s);
        if (!indices) return 0;
        view->buf = indices;
        view->itemsize = 1;
        view->len = n;
        return Py_CLEANUP_SUPPORTED;
    }
    PyErr_Clear();
    if (PyObject_GetBuffer(argument, view, flag) == -1) {
        PyErr_SetString(PyExc_ValueError, "sequence has unexpected format");
        return 0;
    }
    if (view->ndim != 1) {
        PyErr_Format(PyExc_ValueError,
                     "sequence has incorrect rank (%d expected 1)", view->ndim);
        return 0;
    }
    n = view->len / view->itemsize;
    if (n == 0) {
        PyErr_SetString(PyExc_ValueError, "sequence has zero length");
        return 0;
    }
    if (strcmp(view->format, "c") == 0
     || strcmp(view->format, "B") == 0) {
        if (view->itemsize != sizeof(char)) {
            PyErr_Format(PyExc_ValueError,
                        "sequence has unexpected item byte size "
                        "(%ld, expected %ld)", view->itemsize, sizeof(char));
            return 0;
        }
        indices = convert_sequence_to_ints(aligner->mapping, n, view->buf);
        if (!indices) return 0;
        view->itemsize = 1;
        view->len = n;
        view->buf = indices;
    }
    else if (strcmp(view->format, "i") == 0
          || strcmp(view->format, "l") == 0) {
        if (view->itemsize != sizeof(int)) {
            PyErr_Format(PyExc_ValueError,
                        "sequence has unexpected item byte size "
                        "(%ld, expected %ld)", view->itemsize, sizeof(int));
            return 0;
        }
        indices = view->buf;
        if (aligner->substitution_matrix.obj) {
            const Py_ssize_t m = aligner->substitution_matrix.shape[0];
            for (i = 0; i < n; i++) {
                index = indices[i];
                if (index < 0) {
                    PyErr_Format(PyExc_ValueError,
                                 "sequence item %zd is negative (%d)",
                                 i, index);
                    return 0;
                }
                if (index >= m) {
                    PyErr_Format(PyExc_ValueError,
                                 "sequence item %zd is out of bound"
                                 " (%d, should be < %zd)", i, index, m);
                    return 0;
                }
            }
        }
    }
    else {
        PyErr_Format(PyExc_ValueError,
                     "sequence has incorrect data type '%s'", view->format);
        return 0;
    }
    return Py_CLEANUP_SUPPORTED;
}
 
static const char Aligner_score__doc__[] = "calculates the alignment score";

static PyObject*
Aligner_score(Aligner* self, PyObject* args, PyObject* keywords)
{
    const int* sA;
    const int* sB;
    Py_ssize_t nA;
    Py_ssize_t nB;
    Py_buffer bA = {0};
    Py_buffer bB = {0};
    const Mode mode = self->mode;
    const Algorithm algorithm = _get_algorithm(self);
    PyObject* result = NULL;
    PyObject* substitution_matrix = self->substitution_matrix.obj;

    static char *kwlist[] = {"sequenceA", "sequenceB", NULL};

    bA.obj = (PyObject*)self;
    bB.obj = (PyObject*)self;
    if(!PyArg_ParseTupleAndKeywords(args, keywords, "O&O&", kwlist,
                                    sequence_converter, &bA,
                                    sequence_converter, &bB))
        return NULL;

    sA = bA.buf;
    nA = bA.len / bA.itemsize;
    sB = bB.buf;
    nB = bB.len / bB.itemsize;

    switch (algorithm) {
        case NeedlemanWunschSmithWaterman:
            switch (mode) {
                case Global:
                    if (substitution_matrix)
                        result = Aligner_needlemanwunsch_score_matrix(self, sA, nA, sB, nB);
                    else
                        result = Aligner_needlemanwunsch_score_compare(self, sA, nA, sB, nB);
                    break;
                case Local:
                    if (substitution_matrix)
                        result = Aligner_smithwaterman_score_matrix(self, sA, nA, sB, nB);
                    else
                        result = Aligner_smithwaterman_score_compare(self, sA, nA, sB, nB);
                    break;
            }
            break;
        case Gotoh:
            switch (mode) {
                case Global:
                    if (substitution_matrix)
                        result = Aligner_gotoh_global_score_matrix(self, sA, nA, sB, nB);
                    else
                        result = Aligner_gotoh_global_score_compare(self, sA, nA, sB, nB);
                    break;
                case Local:
                    if (substitution_matrix)
                        result = Aligner_gotoh_local_score_matrix(self, sA, nA, sB, nB);
                    else
                        result = Aligner_gotoh_local_score_compare(self, sA, nA, sB, nB);
                    break;
            }
            break;
        case WatermanSmithBeyer:
            switch (mode) {
                case Global:
                    if (substitution_matrix)
                        result = Aligner_watermansmithbeyer_global_score_matrix(self, sA, nA, sB, nB);
                    else
                        result = Aligner_watermansmithbeyer_global_score_compare(self, sA, nA, sB, nB);
                    break;
                case Local:
                    if (substitution_matrix)
                        result = Aligner_watermansmithbeyer_local_score_matrix(self, sA, nA, sB, nB);
                    else
                        result = Aligner_watermansmithbeyer_local_score_compare(self, sA, nA, sB, nB);
                    break;
            }
            break;
        case Unknown:
        default:
            PyErr_SetString(PyExc_RuntimeError, "unknown algorithm");
            break;
    }

    sequence_converter(NULL, &bA);
    sequence_converter(NULL, &bB);

    return result;
}

static const char Aligner_align__doc__[] = "align two sequences";

static PyObject*
Aligner_align(Aligner* self, PyObject* args, PyObject* keywords)
{
    const int* sA;
    const int* sB;
    Py_ssize_t nA;
    Py_ssize_t nB;
    Py_buffer bA = {0};
    Py_buffer bB = {0};
    const Mode mode = self->mode;
    const Algorithm algorithm = _get_algorithm(self);
    PyObject* result = NULL;
    PyObject* substitution_matrix = self->substitution_matrix.obj;

    static char *kwlist[] = {"sequenceA", "sequenceB", NULL};

    bA.obj = (PyObject*)self;
    bB.obj = (PyObject*)self;
    if(!PyArg_ParseTupleAndKeywords(args, keywords, "O&O&", kwlist,
                                    sequence_converter, &bA,
                                    sequence_converter, &bB))
        return NULL;

    sA = bA.buf;
    nA = bA.len / bA.itemsize;
    sB = bB.buf;
    nB = bB.len / bB.itemsize;

    switch (algorithm) {
        case NeedlemanWunschSmithWaterman:
            switch (mode) {
                case Global:
                    if (substitution_matrix)
                        result = Aligner_needlemanwunsch_align_matrix(self, sA, nA, sB, nB);
                    else
                        result = Aligner_needlemanwunsch_align_compare(self, sA, nA, sB, nB);
                    break;
                case Local:
                    if (substitution_matrix)
                        result = Aligner_smithwaterman_align_matrix(self, sA, nA, sB, nB);
                    else
                        result = Aligner_smithwaterman_align_compare(self, sA, nA, sB, nB);
                    break;
            }
            break;
        case Gotoh:
            switch (mode) {
                case Global:
                    if (substitution_matrix)
                        result = Aligner_gotoh_global_align_matrix(self, sA, nA, sB, nB);
                    else
                        result = Aligner_gotoh_global_align_compare(self, sA, nA, sB, nB);
                    break;
                case Local:
                    if (substitution_matrix)
                        result = Aligner_gotoh_local_align_matrix(self, sA, nA, sB, nB);
                    else
                        result = Aligner_gotoh_local_align_compare(self, sA, nA, sB, nB);
                    break;
            }
            break;
        case WatermanSmithBeyer:
            switch (mode) {
                case Global:
                    if (substitution_matrix)
                        result = Aligner_watermansmithbeyer_global_align_matrix(self, sA, nA, sB, nB);
                    else
                        result = Aligner_watermansmithbeyer_global_align_compare(self, sA, nA, sB, nB);
                    break;
                case Local:
                    if (substitution_matrix)
                        result = Aligner_watermansmithbeyer_local_align_matrix(self, sA, nA, sB, nB);
                    else
                        result = Aligner_watermansmithbeyer_local_align_compare(self, sA, nA, sB, nB);
                    break;
            }
            break;
        case Unknown:
        default:
            PyErr_SetString(PyExc_RuntimeError, "unknown algorithm");
            break;
    }

    sequence_converter(NULL, &bA);
    sequence_converter(NULL, &bB);

    return result;
}

static char Aligner_doc[] =
"Aligner.\n";

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
    {NULL}  /* Sentinel */
};

static PyTypeObject AlignerType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_algorithms.PairwiseAligner", /* tp_name */
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

static char _aligners__doc__[] =
"C extension module implementing pairwise alignment algorithms";

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_aligners",
        _aligners__doc__,
        -1,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL
};

PyObject *
PyInit__aligners(void)
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
                           "PairwiseAligner", (PyObject*) &AlignerType) < 0) {
        Py_DECREF(&AlignerType);
        Py_DECREF(module);
        return NULL;
    }

    return module;
}
