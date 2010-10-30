/* Copyright 2002 by Jeffrey Chang.  All rights reserved.
 * This code is part of the Biopython distribution and governed by its
 * license.  Please see the LICENSE file that should have been included
 * as part of this package.
 *
 * cpairwise2module.c
 * Created 30 Sep 2001
 *
 * Optimized C routines that complement pairwise2.py.
 */

#include "Python.h"


#define _PRECISION 1000
#define rint(x) (int)((x)*_PRECISION+0.5)


/* Return a PyNumber as a double.
 * Raises a TypeError if I can't do it.
 */
static double PyNumber_AsDouble(PyObject *py_num)
{
    double val;
    PyObject *floatobj;

    if((floatobj = PyNumber_Float(py_num)) == NULL)
	return(0.0);
    val = PyFloat_AsDouble(floatobj);
    Py_DECREF(floatobj);
    return val;
}

/* Functions in this module. */

double calc_affine_penalty(int length, double open, double extend, 
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

struct IndexList {
    int num_used;
    int num_allocated;
    int *indexes; /* Array of ints.  Even ints are rows, odd ones are cols. */
};

static void IndexList_init(struct IndexList *il) 
{
    /* il->num_used = 0;
    il->num_allocated = 0;
    il->indexes = NULL; */
    memset((void *)il, 0, sizeof(struct IndexList));
}

static void IndexList_free(struct IndexList *il) 
{
    if(il->indexes)
	free(il->indexes);
    IndexList_init(il);
}

static void IndexList_clear(struct IndexList *il)
{
    il->num_used = 0;
}

static int IndexList_contains(struct IndexList *il, 
			      const int row, const int col) 
{
    int i;
    int stop;

    stop = il->num_used*2;
    for(i=0; i<stop; i+=2) {
	if((il->indexes[i] == row) && (il->indexes[i+1] == col))
	    return 1;
    }
    return 0;
}

static int IndexList__verify_free_index(struct IndexList *il, int num_needed)
{
    int *indexes;
    int num_to_allocate;
    int num_bytes;

    if(il->num_allocated >= num_needed)
	return 1;

    /* Nearly all the cases are for list of length 1 or 2.
       Empirically, the code seems to run fastest when I allocate the
       exact amount for these. */
    if(num_needed <= 2)
	num_to_allocate = num_needed;
    else 
	num_to_allocate = num_needed*2;
    num_bytes = num_to_allocate*sizeof(int)*2;
    if(!(indexes = realloc((void *)il->indexes, num_bytes))) {
	PyErr_SetString(PyExc_MemoryError, "Out of memory");
	return 0;
    }
    il->indexes = indexes;
    il->num_allocated = num_to_allocate;
    return 1;
}

static void IndexList_append(struct IndexList *il, 
			     const int row, const int col) 
{
    int i;
    if(!IndexList__verify_free_index(il, il->num_used+1))
	return;
    i=il->num_used*2;
    il->indexes[i] = row;
    il->indexes[i+1] = col;
    il->num_used += 1;
}

static void IndexList_extend(struct IndexList *il1, struct IndexList *il2) 
{
    int i1, i2;
    int stop;

    if(!IndexList__verify_free_index(il1, il1->num_used+il2->num_used))
	return;
    stop=il2->num_used * 2;
    for(i1=il1->num_used*2, i2=0; i2<stop; i1+=2, i2+=2) {
	il1->indexes[i1] = il2->indexes[i2];
	il1->indexes[i1+1] = il2->indexes[i2+1];
    }
    il1->num_used += il2->num_used;
}

/* static void IndexList_copy(struct IndexList *il1, struct IndexList *il2) 
{
    IndexList_clear(il1);
    IndexList_extend(il1, il2);
    } */

/* static void IndexList_assign(struct IndexList *il1, struct IndexList *il2)
{
    if(il1->indexes)
	free(il1->indexes);
    il1->indexes = il2->indexes;
    il1->num_used = il2->num_used;
    il1->num_allocated = il2->num_allocated;
    } */


double _get_match_score(PyObject *py_sequenceA, PyObject *py_sequenceB, 
			PyObject *py_match_fn, int i, int j,
			char *sequenceA, char *sequenceB, 
			int use_sequence_cstring,
			double match, double mismatch, 
			int use_match_mismatch_scores)
{
    PyObject *py_A=NULL, 
	*py_B=NULL;
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
    score = PyNumber_AsDouble(py_result);
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


/* This function is a more-or-less straightforward port of the
 * equivalent function in pairwise2.  Please see there for algorithm
 * documentation.
 */
static PyObject *cpairwise2__make_score_matrix_fast(
    PyObject *self, PyObject *args)
{
    int i;
    int row, col;

    PyObject *py_sequenceA, *py_sequenceB, *py_match_fn;
    char *sequenceA=NULL, *sequenceB=NULL;
    int use_sequence_cstring;
    double open_A, extend_A, open_B, extend_B;
    int penalize_extend_when_opening, penalize_end_gaps;
    int align_globally, score_only;

    double first_A_gap, first_B_gap;
    double match, mismatch;
    int use_match_mismatch_scores;
    int lenA, lenB;
    double *score_matrix = (double *)NULL;
    struct IndexList *trace_matrix = (struct IndexList *)NULL;
    PyObject *py_score_matrix=NULL, *py_trace_matrix=NULL;

    double *row_cache_score = (double *)NULL, 
	*col_cache_score = (double *)NULL;
    struct IndexList *row_cache_index = (struct IndexList *)NULL, 
	*col_cache_index = (struct IndexList *)NULL;

    PyObject *py_retval = NULL;

    if(!PyArg_ParseTuple(args, "OOOddddiiii", &py_sequenceA, &py_sequenceB,
			 &py_match_fn, &open_A, &extend_A, &open_B, &extend_B,
			 &penalize_extend_when_opening, &penalize_end_gaps,
			 &align_globally, &score_only))
	return NULL;
    if(!PySequence_Check(py_sequenceA) || !PySequence_Check(py_sequenceB)) {
	PyErr_SetString(PyExc_TypeError, 
			"py_sequenceA and py_sequenceB should be sequences.");
	return NULL;
    }

    /* Optimize for the common case.  Check to see if py_sequenceA and
       py_sequenceB are strings.  If they are, use the c string
       representation. */
    use_sequence_cstring = 0;
    if(PyString_Check(py_sequenceA) && PyString_Check(py_sequenceB)) {
	sequenceA = PyString_AS_STRING(py_sequenceA);
	sequenceB = PyString_AS_STRING(py_sequenceB);
	use_sequence_cstring = 1;
    }

    if(!PyCallable_Check(py_match_fn)) {
	PyErr_SetString(PyExc_TypeError, "py_match_fn must be callable.");
	return NULL;
    }
    /* Optimize for the common case.  Check to see if py_match_fn is
       an identity_match.  If so, pull out the match and mismatch
       member variables and calculate the scores myself. */
    match = mismatch = 0;
    use_match_mismatch_scores = 0;
    if(PyInstance_Check(py_match_fn)) {
	PyObject *py_match=NULL, *py_mismatch=NULL;
	if(!(py_match = PyObject_GetAttrString(py_match_fn, "match")))
	    goto cleanup_after_py_match_fn;
	match = PyNumber_AsDouble(py_match);
	if(PyErr_Occurred())
	    goto cleanup_after_py_match_fn;
	if(!(py_mismatch = PyObject_GetAttrString(py_match_fn, "mismatch")))
	    goto cleanup_after_py_match_fn;
	mismatch = PyNumber_AsDouble(py_mismatch);
	if(PyErr_Occurred())
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
    }

    /* Cache some commonly used gap penalties */
    first_A_gap = calc_affine_penalty(1, open_A, extend_A, 
				      penalize_extend_when_opening);
    first_B_gap = calc_affine_penalty(1, open_B, extend_B,
				      penalize_extend_when_opening);

    /* Allocate matrices for storing the results and initalize them. */
    lenA = PySequence_Length(py_sequenceA);
    lenB = PySequence_Length(py_sequenceB);
    score_matrix = (double *)malloc(lenA*lenB*sizeof(*score_matrix));
    trace_matrix = (struct IndexList *)malloc(lenA*lenB*sizeof(*trace_matrix));
    if(!score_matrix || !trace_matrix) {
	PyErr_SetString(PyExc_MemoryError, "Out of memory");
	goto _cleanup_make_score_matrix_fast;
    }
    for(i=0; i<lenA*lenB; i++) {
	score_matrix[i] = 0;
	IndexList_init(&trace_matrix[i]);
    }

    /* Initialize the first row and col of the score matrix. */
    for(i=0; i<lenA; i++) {
	double score = _get_match_score(py_sequenceA, py_sequenceB, 
					py_match_fn, i, 0,
					sequenceA, sequenceB,
					use_sequence_cstring,
					match, mismatch,
					use_match_mismatch_scores);
	if(PyErr_Occurred())
	    goto _cleanup_make_score_matrix_fast;
	if(penalize_end_gaps)
	    score += calc_affine_penalty(i, open_B, extend_B, 
					 penalize_extend_when_opening);
	score_matrix[i*lenB] = score;
    }
    for(i=0; i<lenB; i++) {
	double score = _get_match_score(py_sequenceA, py_sequenceB, 
					py_match_fn, 0, i,
					sequenceA, sequenceB,
					use_sequence_cstring,
					match, mismatch,
					use_match_mismatch_scores);
	if(PyErr_Occurred())
	    goto _cleanup_make_score_matrix_fast;
	if(penalize_end_gaps)
	    score += calc_affine_penalty(i, open_A, extend_A, 
					 penalize_extend_when_opening);
	score_matrix[i] = score;
    }

    /* Now initialize the row and col cache. */
    row_cache_score = (double *)malloc((lenA-1)*sizeof(*row_cache_score));
    row_cache_index = (struct IndexList *)malloc((lenA-1)*
						 sizeof(*row_cache_index));
    col_cache_score = (double *)malloc((lenB-1)*sizeof(*col_cache_score));
    col_cache_index = (struct IndexList *)malloc((lenB-1)*
						 sizeof(*col_cache_index));
    if(!row_cache_score || !row_cache_index || 
       !col_cache_score || !col_cache_index) {
	PyErr_SetString(PyExc_MemoryError, "Out of memory");
	goto _cleanup_make_score_matrix_fast;
    }
    memset((void *)row_cache_score, 0, (lenA-1)*sizeof(*row_cache_score));
    memset((void *)row_cache_index, 0, (lenA-1)*sizeof(*row_cache_index));
    memset((void *)col_cache_score, 0, (lenB-1)*sizeof(*col_cache_score));
    memset((void *)col_cache_index, 0, (lenB-1)*sizeof(*col_cache_index));
    for(i=0; i<lenA-1; i++) {
	row_cache_score[i] = score_matrix[i*lenB] + first_A_gap;
	IndexList_append(&row_cache_index[i], i, 0);
    }
    for(i=0; i<lenB-1; i++) {
	col_cache_score[i] = score_matrix[i] + first_B_gap;
	IndexList_append(&col_cache_index[i], 0, i);
    }

    /* Fill in the score matrix. */
    for(row=1; row<lenA; row++) {
	for(col=1; col<lenB; col++) {
	    double nogap_score, row_score, col_score, best_score;
	    int best_score_rint;
	    struct IndexList *il;

	    double score, open_score, extend_score;
	    int open_score_rint, extend_score_rint;

	    /* Calculate the best score. */
	    nogap_score = score_matrix[(row-1)*lenB+col-1];
	    if(col > 1) {
		row_score = row_cache_score[row-1];
	    } else {
		row_score = nogap_score-1; /* Make sure it's not best score */
	    }
	    if(row > 1) {
		col_score = col_cache_score[col-1];
	    } else {
		col_score = nogap_score-1; /* Make sure it's not best score */
	    }

	    best_score = (row_score > col_score) ? row_score : col_score;
	    if(nogap_score > best_score)
		best_score = nogap_score;
	    best_score_rint = rint(best_score);

	    /* Set the score and traceback matrices. */
	    score = best_score + _get_match_score(py_sequenceA, py_sequenceB, 
						  py_match_fn, row, col,
						  sequenceA, sequenceB,
						  use_sequence_cstring,
						  match, mismatch,
						  use_match_mismatch_scores);
	    if(PyErr_Occurred())
		goto _cleanup_make_score_matrix_fast;
	    if(!align_globally && score < 0)
		score_matrix[row*lenB+col] = 0;
	    else
		score_matrix[row*lenB+col] = score;

	    il = &trace_matrix[row*lenB+col];
	    if(best_score_rint == rint(nogap_score)) {
		IndexList_append(il, row-1, col-1);
	    }
	    if(best_score_rint == rint(row_score)) {
		IndexList_extend(il, &row_cache_index[row-1]);
	    }
	    if(best_score_rint == rint(col_score)) {
		IndexList_extend(il, &col_cache_index[col-1]);
	    }

	    /* Update the cached column scores. */
	    open_score = score_matrix[(row-1)*lenB+col-1] + first_B_gap;
	    extend_score = col_cache_score[col-1] + extend_B;
	    open_score_rint = rint(open_score);
	    extend_score_rint = rint(extend_score);
	    if(open_score_rint > extend_score_rint) {
		col_cache_score[col-1] = open_score;
		IndexList_clear(&col_cache_index[col-1]);
		IndexList_append(&col_cache_index[col-1], row-1, col-1);
	    } else if(extend_score_rint > open_score_rint) {
		col_cache_score[col-1] = extend_score;
	    } else {
		col_cache_score[col-1] = open_score;
		if(!IndexList_contains(&col_cache_index[col-1], row-1, col-1))
		    IndexList_append(&col_cache_index[col-1], row-1, col-1);
	    }
	    
	    /* Update the cached row scores. */
	    open_score = score_matrix[(row-1)*lenB+col-1] + first_A_gap;
	    extend_score = row_cache_score[row-1] + extend_A;
	    open_score_rint = rint(open_score);
	    extend_score_rint = rint(extend_score);
	    if(open_score_rint > extend_score_rint) {
		row_cache_score[row-1] = open_score;
		IndexList_clear(&row_cache_index[row-1]);
		IndexList_append(&row_cache_index[row-1], row-1, col-1);
	    } else if(extend_score_rint > open_score_rint) {
		row_cache_score[row-1] = extend_score;
	    } else {
		row_cache_score[row-1] = open_score;
		if(!IndexList_contains(&row_cache_index[row-1], row-1, col-1))
		    IndexList_append(&row_cache_index[row-1], row-1, col-1);
	    }
	}
    }

    /* Save the score and traceback matrices into real python objects. */
    if(!(py_score_matrix = PyList_New(lenA)))
	goto _cleanup_make_score_matrix_fast;
    if(!(py_trace_matrix = PyList_New(lenA)))
	goto _cleanup_make_score_matrix_fast;
    for(row=0; row<lenA; row++) {
	PyObject *py_score_row, *py_trace_row;
	if(!(py_score_row = PyList_New(lenB)))
	    goto _cleanup_make_score_matrix_fast;
	PyList_SET_ITEM(py_score_matrix, row, py_score_row);
	if(!(py_trace_row = PyList_New(lenB)))
	    goto _cleanup_make_score_matrix_fast;
	PyList_SET_ITEM(py_trace_matrix, row, py_trace_row);

	for(col=0; col<lenB; col++) {
	    int i;
	    PyObject *py_score, *py_indexlist;
	    int offset = row*lenB + col;
	    struct IndexList *il = &trace_matrix[offset];

	    /* Set py_score_matrix[row][col] to the score. */
	    if(!(py_score = PyFloat_FromDouble(score_matrix[offset])))
		goto _cleanup_make_score_matrix_fast;
	    PyList_SET_ITEM(py_score_row, col, py_score);

	    if(score_only)
		continue;
	    /* Set py_trace_matrix[row][col] to a list of indexes.  On
	       the edges of the matrix (row or column is 0), the
	       matrix should be [None]. */
	    if(!row || !col) {
		if(!(py_indexlist = PyList_New(1)))
		    goto _cleanup_make_score_matrix_fast;
		Py_INCREF(Py_None);
		PyList_SET_ITEM(py_indexlist, 0, Py_None);
	    }
	    else {
		if(!(py_indexlist = PyList_New(il->num_used)))
		    goto _cleanup_make_score_matrix_fast;
		for(i=0; i<il->num_used; i++) {
		    PyObject *py_index=NULL;
		    int row = il->indexes[i*2],
			col = il->indexes[i*2+1];
		    if(!(py_index = Py_BuildValue("(ii)", row, col)))
			goto _cleanup_make_score_matrix_fast;
		    PyList_SET_ITEM(py_indexlist, i, py_index);
		}
	    }
	    PyList_SET_ITEM(py_trace_row, col, py_indexlist);
	}
    }

    py_retval = Py_BuildValue("(OO)", py_score_matrix, py_trace_matrix);


 _cleanup_make_score_matrix_fast:
    if(score_matrix)
	free(score_matrix);
    if(trace_matrix) {
	for(i=0; i<lenA*lenB; i++)
	    IndexList_free(&trace_matrix[i]);
	free(trace_matrix);
    }
    if(row_cache_score)
	free(row_cache_score);
    if(col_cache_score)
	free(col_cache_score);
    if(row_cache_index) {
	for(i=0; i<lenA-1; i++)
	    IndexList_free(&row_cache_index[i]);
	free(row_cache_index);
    }
    if(col_cache_index) {
	for(i=0; i<lenB-1; i++) {
	    IndexList_free(&col_cache_index[i]);
	}
	free(col_cache_index);
    }
    if(py_score_matrix) {
	Py_DECREF(py_score_matrix);
    }
    if(py_trace_matrix) {
	Py_DECREF(py_trace_matrix);
    }

    return py_retval;
}

static PyObject *cpairwise2_rint(
    PyObject *self, PyObject *args, PyObject *keywds)
{
    double x;
    int precision = _PRECISION;
    int rint_x;

    static char *kwlist[] = {"x", "precision", NULL};

    if(!PyArg_ParseTupleAndKeywords(args, keywds, "d|l", kwlist,
				    &x, &precision))
	return NULL;
    rint_x = (int)(x * precision + 0.5);
    return PyInt_FromLong(rint_x);
}

/* Module definition stuff */

static PyMethodDef cpairwise2Methods[] = {
    {"_make_score_matrix_fast", 
     (PyCFunction)cpairwise2__make_score_matrix_fast, METH_VARARGS, ""},
    {"rint", (PyCFunction)cpairwise2_rint, METH_VARARGS|METH_KEYWORDS, ""},
    {NULL, NULL, 0, NULL}
};

static char cpairwise2__doc__[] =
"XXX document here\n\
\n\
";

void initcpairwise2(void)
{
    (void) Py_InitModule3("cpairwise2", cpairwise2Methods, cpairwise2__doc__);
}
