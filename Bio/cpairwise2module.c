/* Copyright 2002 by Jeffrey Chang.
 * Copyright 2016, 2019 by Markus Piotrowski.
 * All rights reserved.
 *
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 *
 * cpairwise2module.c
 * Created 30 Sep 2001
 *
 * Optimized C routines that complement pairwise2.py.
 */

#include "Python.h"


#define _PRECISION 1000
#define rint(x) (int)((x)*_PRECISION+0.5)

/* Functions in this module. */

static double calc_affine_penalty(int length, double open, double extend,
    int penalize_extend_when_opening)
{
    double penalty;

    if(length <= 0)
        return 0.0;
    penalty = open + extend * length;
    if(!penalize_extend_when_opening)
        penalty -= extend;
    return penalty;
}

static double _get_match_score(PyObject *py_sequenceA, PyObject *py_sequenceB,
                               PyObject *py_match_fn, int i, int j,
                               char *sequenceA, char *sequenceB,
                               int use_sequence_cstring,
                               double match, double mismatch,
                               int use_match_mismatch_scores)
{
    PyObject *py_A=NULL, *py_B=NULL;
    PyObject *py_arglist=NULL, *py_result=NULL;
    double score = 0;

    if(use_sequence_cstring && use_match_mismatch_scores) {
        score = (sequenceA[i] == sequenceB[j]) ? match : mismatch;
        return score;
    }
    /* Calculate the match score. */
    if(!(py_A = PySequence_GetItem(py_sequenceA, i)))
        goto _get_match_score_cleanup;
    if(!(py_B = PySequence_GetItem(py_sequenceB, j)))
        goto _get_match_score_cleanup;
    if(!(py_arglist = Py_BuildValue("(OO)", py_A, py_B)))
        goto _get_match_score_cleanup;

    if(!(py_result = PyEval_CallObject(py_match_fn, py_arglist)))
        goto _get_match_score_cleanup;
    score = PyFloat_AsDouble(py_result);

 _get_match_score_cleanup:
    if(py_A) {
        Py_DECREF(py_A);
    }
    if(py_B) {
        Py_DECREF(py_B);
    }
    if(py_arglist) {
        Py_DECREF(py_arglist);
    }
    if(py_result) {
        Py_DECREF(py_result);
    }
    return score;
}

#if PY_MAJOR_VERSION >= 3
static PyObject* _create_bytes_object(PyObject* o)
{
    PyObject* b;
    if (PyBytes_Check(o)) {
        return o;
    }
    if (!PyUnicode_Check(o)) {
        return NULL;
    }
    b = PyUnicode_AsASCIIString(o);
    if (!b) {
        PyErr_Clear();
        return NULL;
    }
    return b;
}
#endif

/* This function is a more-or-less straightforward port of the
 * equivalent function in pairwise2. Please see there for algorithm
 * documentation.
 */
static PyObject *cpairwise2__make_score_matrix_fast(PyObject *self,
                                                    PyObject *args)
{
    int i;
    int row, col;
    PyObject *py_sequenceA, *py_sequenceB, *py_match_fn;
#if PY_MAJOR_VERSION >= 3
    PyObject *py_bytesA, *py_bytesB;
#endif
    char *sequenceA=NULL, *sequenceB=NULL;
    int use_sequence_cstring;
    double open_A, extend_A, open_B, extend_B;
    int penalize_extend_when_opening, penalize_end_gaps_A, penalize_end_gaps_B;
    int align_globally, score_only;

    PyObject *py_match=NULL, *py_mismatch=NULL;
    double first_A_gap, first_B_gap;
    double match, mismatch;
    double score;
    double best_score = 0;
    double local_max_score = 0;
    int use_match_mismatch_scores;
    int lenA, lenB;
    double *score_matrix = NULL;
    unsigned char *trace_matrix = NULL;
    PyObject *py_score_matrix=NULL, *py_trace_matrix=NULL;

    double *col_cache_score = NULL;
    PyObject *py_retval = NULL;

    if(!PyArg_ParseTuple(args, "OOOddddi(ii)ii", &py_sequenceA, &py_sequenceB,
                         &py_match_fn, &open_A, &extend_A, &open_B, &extend_B,
                         &penalize_extend_when_opening,
                         &penalize_end_gaps_A, &penalize_end_gaps_B,
                         &align_globally, &score_only))
        return NULL;
    if(!PySequence_Check(py_sequenceA) || !PySequence_Check(py_sequenceB)) {
        PyErr_SetString(PyExc_TypeError,
                        "py_sequenceA and py_sequenceB should be sequences.");
        return NULL;
    }

    /* Optimize for the common case. Check to see if py_sequenceA and
       py_sequenceB are strings.  If they are, use the c string
       representation. */
#if PY_MAJOR_VERSION < 3
    use_sequence_cstring = 0;
    if(PyString_Check(py_sequenceA) && PyString_Check(py_sequenceB)) {
        sequenceA = PyString_AS_STRING(py_sequenceA);
        sequenceB = PyString_AS_STRING(py_sequenceB);
        use_sequence_cstring = 1;
    }
#else
    py_bytesA = _create_bytes_object(py_sequenceA);
    py_bytesB = _create_bytes_object(py_sequenceB);
    if (py_bytesA && py_bytesB) {
        sequenceA = PyBytes_AS_STRING(py_bytesA);
        sequenceB = PyBytes_AS_STRING(py_bytesB);
        use_sequence_cstring = 1;
    }
    else {
        Py_XDECREF(py_bytesA);
        Py_XDECREF(py_bytesB);
        use_sequence_cstring = 0;
    }
#endif

    if(!PyCallable_Check(py_match_fn)) {
        PyErr_SetString(PyExc_TypeError, "py_match_fn must be callable.");
        return NULL;
    }
    /* Optimize for the common case. Check to see if py_match_fn is
       an identity_match. If so, pull out the match and mismatch
       member variables and calculate the scores myself. */
    match = mismatch = 0;
    use_match_mismatch_scores = 0;
    if(!(py_match = PyObject_GetAttrString(py_match_fn, "match")))
        goto cleanup_after_py_match_fn;
    match = PyFloat_AsDouble(py_match);
    if(match==-1.0 && PyErr_Occurred())
        goto cleanup_after_py_match_fn;
    if(!(py_mismatch = PyObject_GetAttrString(py_match_fn, "mismatch")))
        goto cleanup_after_py_match_fn;
    mismatch = PyFloat_AsDouble(py_mismatch);
    if(mismatch==-1.0 && PyErr_Occurred())
        goto cleanup_after_py_match_fn;
    use_match_mismatch_scores = 1;

 cleanup_after_py_match_fn:
    if(PyErr_Occurred())
        PyErr_Clear();
    if(py_match) {
        Py_DECREF(py_match);
    }
    if(py_mismatch) {
        Py_DECREF(py_mismatch);
    }
    /* Cache some commonly used gap penalties */
    first_A_gap = calc_affine_penalty(1, open_A, extend_A,
                                      penalize_extend_when_opening);
    first_B_gap = calc_affine_penalty(1, open_B, extend_B,
                                      penalize_extend_when_opening);

    /* Allocate matrices for storing the results and initialize first row and col. */
    lenA = PySequence_Length(py_sequenceA);
    lenB = PySequence_Length(py_sequenceB);
    score_matrix = malloc((lenA+1)*(lenB+1)*sizeof(*score_matrix));
    if(!score_matrix) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory");
        goto _cleanup_make_score_matrix_fast;
    }
    for(i=0; i<(lenB+1); i++)
        score_matrix[i] = 0;
    for(i=0; i<(lenA+1)*(lenB+1); i += (lenB+1))
        score_matrix[i] = 0;
    /* If we only want the score, we don't need the trace matrix. */
    if (!score_only){
        trace_matrix = malloc((lenA+1)*(lenB+1)*sizeof(*trace_matrix));
        if(!trace_matrix) {
            PyErr_SetString(PyExc_MemoryError, "Out of memory");
            goto _cleanup_make_score_matrix_fast;
        }
        for(i=0; i<(lenB+1); i++)
            trace_matrix[i] = 0;
        for(i=0; i<(lenA+1)*(lenB+1); i += (lenB+1))
            trace_matrix[i] = 0;
        }
    else
        trace_matrix = malloc(1);

    /* Initialize the first row and col of the score matrix. */
    for(i=0; i<=lenA; i++) {
        if(penalize_end_gaps_B)
            score = calc_affine_penalty(i, open_B, extend_B,
                                        penalize_extend_when_opening);
        else
            score = 0;
        score_matrix[i*(lenB+1)] = score;
    }
    for(i=0; i<=lenB; i++) {
        if(penalize_end_gaps_A)
            score = calc_affine_penalty(i, open_A, extend_A,
                                        penalize_extend_when_opening);
        else
            score = 0;
        score_matrix[i] = score;
    }

    /* Now initialize the col cache. */
    col_cache_score = malloc((lenB+1)*sizeof(*col_cache_score));
    memset((void *)col_cache_score, 0, (lenB+1)*sizeof(*col_cache_score));
    for(i=0; i<=lenB; i++) {
        col_cache_score[i] = calc_affine_penalty(i, (2*open_B), extend_B,
                             penalize_extend_when_opening);
    }

    /* Fill in the score matrix. The row cache is calculated on the fly.*/
    for(row=1; row<=lenA; row++) {
        double row_cache_score = calc_affine_penalty(row, (2*open_A), extend_A,
                                 penalize_extend_when_opening);
        for(col=1; col<=lenB; col++) {
            double match_score, nogap_score;
            double row_open, row_extend, col_open, col_extend;
            int best_score_rint, row_score_rint, col_score_rint;
            unsigned char row_trace_score, col_trace_score, trace_score;

            /* Calculate the best score. */
            match_score = _get_match_score(py_sequenceA, py_sequenceB,
                                           py_match_fn, row-1, col-1,
                                           sequenceA, sequenceB,
                                           use_sequence_cstring,
                                           match, mismatch,
                                           use_match_mismatch_scores);
            if(match_score==-1.0 && PyErr_Occurred())
                goto _cleanup_make_score_matrix_fast;
            nogap_score = score_matrix[(row-1)*(lenB+1)+col-1] + match_score;

            if (!penalize_end_gaps_A && row==lenA) {
                row_open = score_matrix[(row)*(lenB+1)+col-1];
                row_extend = row_cache_score;
            }
            else {
                row_open = score_matrix[(row)*(lenB+1)+col-1] + first_A_gap;
                row_extend = row_cache_score + extend_A;
            }
            row_cache_score = (row_open > row_extend) ? row_open : row_extend;

            if (!penalize_end_gaps_B && col==lenB){
                col_open = score_matrix[(row-1)*(lenB+1)+col];
                col_extend = col_cache_score[col];
            }
            else {
                col_open = score_matrix[(row-1)*(lenB+1)+col] + first_B_gap;
                col_extend = col_cache_score[col] + extend_B;
            }
            col_cache_score[col] = (col_open > col_extend) ? col_open : col_extend;

            best_score = (row_cache_score > col_cache_score[col]) ? row_cache_score : col_cache_score[col];
            if(nogap_score > best_score)
                best_score = nogap_score;

            if (best_score > local_max_score)
                local_max_score = best_score;

            if(!align_globally && best_score < 0)
                score_matrix[row*(lenB+1)+col] = 0;
            else
                score_matrix[row*(lenB+1)+col] = best_score;

            if (!score_only) {
                row_score_rint = rint(row_cache_score);
                col_score_rint = rint(col_cache_score[col]);
                row_trace_score = 0;
                col_trace_score = 0;
                if (rint(row_open) == row_score_rint)
                    row_trace_score = row_trace_score|1;
                if (rint(row_extend) == row_score_rint)
                    row_trace_score = row_trace_score|8;
                if (rint(col_open) == col_score_rint)
                    col_trace_score = col_trace_score|4;
                if (rint(col_extend) == col_score_rint)
                    col_trace_score = col_trace_score|16;

                trace_score = 0;
                best_score_rint = rint(best_score);
                if (rint(nogap_score) == best_score_rint)
                    trace_score = trace_score|2;
                if (row_score_rint == best_score_rint)
                    trace_score += row_trace_score;
                if (col_score_rint == best_score_rint)
                    trace_score += col_trace_score;
                trace_matrix[row*(lenB+1)+col] = trace_score;
            }
        }
    }

    if (!align_globally)
        best_score = local_max_score;

    /* Save the score and traceback matrices into real python objects. */
	if(!score_only) {
		if(!(py_score_matrix = PyList_New(lenA+1)))
			goto _cleanup_make_score_matrix_fast;
		if(!(py_trace_matrix = PyList_New(lenA+1)))
			goto _cleanup_make_score_matrix_fast;

		for(row=0; row<=lenA; row++) {
			PyObject *py_score_row, *py_trace_row;
			if(!(py_score_row = PyList_New(lenB+1)))
				goto _cleanup_make_score_matrix_fast;
			PyList_SET_ITEM(py_score_matrix, row, py_score_row);
			if(!(py_trace_row = PyList_New(lenB+1)))
				goto _cleanup_make_score_matrix_fast;
			PyList_SET_ITEM(py_trace_matrix, row, py_trace_row);

			for(col=0; col<=lenB; col++) {
				PyObject *py_score, *py_trace;
				int offset = row*(lenB+1) + col;

				/* Set py_score_matrix[row][col] to the score. */
				if(!(py_score = PyFloat_FromDouble(score_matrix[offset])))
					goto _cleanup_make_score_matrix_fast;
				PyList_SET_ITEM(py_score_row, col, py_score);

				/* Set py_trace_matrix[row][col] to a list of indexes.  On
				   the edges of the matrix (row or column is 0), the
				   matrix should be [None]. */
				if(!row || !col) {
					if(!(py_trace = Py_BuildValue("B", 1)))
						goto _cleanup_make_score_matrix_fast;
					Py_INCREF(Py_None);
					PyList_SET_ITEM(py_trace_row, col, Py_None);
				}
				else {
					if(!(py_trace = Py_BuildValue("B", trace_matrix[offset])))
						goto _cleanup_make_score_matrix_fast;
					PyList_SET_ITEM(py_trace_row, col, py_trace);

				}
			}
		}
	}
	else {
		py_score_matrix = PyList_New(1);
		py_trace_matrix = PyList_New(1);
	}
    py_retval = Py_BuildValue("(OOd)", py_score_matrix, py_trace_matrix, best_score);

 _cleanup_make_score_matrix_fast:
    if(score_matrix)
        free(score_matrix);
    if(trace_matrix)
        free(trace_matrix);
    if(col_cache_score)
        free(col_cache_score);
    if(py_score_matrix){
        Py_DECREF(py_score_matrix);
    }
    if(py_trace_matrix){
        Py_DECREF(py_trace_matrix);
    }

#if PY_MAJOR_VERSION >= 3
    if (py_bytesA != NULL && py_bytesA != py_sequenceA) Py_DECREF(py_bytesA);
    if (py_bytesB != NULL && py_bytesB != py_sequenceB) Py_DECREF(py_bytesB);
#endif

    return py_retval;
}

static PyObject *cpairwise2_rint(PyObject *self, PyObject *args,
                                 PyObject *keywds)
{
    double x;
    int precision = _PRECISION;
    int rint_x;

    static char *kwlist[] = {"x", "precision", NULL};

    if(!PyArg_ParseTupleAndKeywords(args, keywds, "d|l", kwlist,
                                    &x, &precision))
        return NULL;
    rint_x = (int)(x * precision + 0.5);
#if PY_MAJOR_VERSION >= 3
    return PyLong_FromLong((long)rint_x);
#else
    return PyInt_FromLong((long)rint_x);
#endif
}

/* Module definition stuff */

static PyMethodDef cpairwise2Methods[] = {
    {"_make_score_matrix_fast",
     (PyCFunction)cpairwise2__make_score_matrix_fast, METH_VARARGS, ""},
    {"rint", (PyCFunction)cpairwise2_rint, METH_VARARGS|METH_KEYWORDS, ""},
    {NULL, NULL, 0, NULL}
};

static char cpairwise2__doc__[] =
"Optimized C routines that complement pairwise2.py. These are called from within pairwise2.py.\n\
\n\
";

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "cpairwise2",
        cpairwise2__doc__,
        -1,
        cpairwise2Methods,
        NULL,
        NULL,
        NULL,
        NULL
};

PyObject *
PyInit_cpairwise2(void)

#else

void
/* for Windows: _declspec(dllexport) initcpairwise2(void) */
initcpairwise2(void)
#endif

{
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
    if (module==NULL) return NULL;
    return module;
#else
    (void) Py_InitModule3("cpairwise2", cpairwise2Methods, cpairwise2__doc__);
#endif
}
