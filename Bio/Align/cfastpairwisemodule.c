/* Copyright 2001 by Jeffrey Chang.  All rights reserved.
 * This code is part of the Biopython distribution and governed by its
 * license.  Please see the LICENSE file that should have been included
 * as part of this package.
 *
 * cfastpairwisemodule.c
 * Created 30 Sep 2001
 *
 * Optimized C routines that complement fastpairwise.py.
 */

#include "Python.h"
#include <string.h>


#define GO_D 1
#define GO_P 2
#define GO_Q 4

#define _PRECISION 1000
#define rint(x) ((int)((x) * _PRECISION + 0.5))



/* Functions in this module. */

static double PyNumber_AsDouble(PyObject *py_num)
{
    double val;
    PyObject *floatobj;

    if(!PyNumber_Check(py_num)) {
	PyErr_SetString(PyExc_TypeError, "I received a non-number");
	return(0.0);
    }
    if((floatobj = PyNumber_Float(py_num)) == NULL)
	return(0.0);
    val = PyFloat_AsDouble(floatobj);
    Py_DECREF(floatobj);
    return val;
}

double 
calc_affine_penalty(length, open, extend, count_first)
     int length;
     double open, extend;
     int count_first;
{
    double penalty;

    if(length <= 0)
	return 0.0;
    penalty = open + extend * length;
    if(!count_first)
	penalty -= extend;
    return penalty;
}


static PyObject *
cfastpairwise__make_score_matrix(self, args)
     PyObject *self;
     PyObject *args;
{
    int i, j;

    PyObject *py_sequenceA, *py_sequenceB, *py_match_fn;
    PyObject *py_A, *py_B;
    double open_A, extend_A, open_B, extend_B;
    int count_first, global_alignment, penalize_end_gaps;
    double forgiveness;
    int forgiveness_rint;

    PyObject *arglist, *result;
    PyObject *tmp;

    int M, N;   /* size of matrix */
    int dsize, csize;
    double *Dmatrix, *Pmatrix, *Qmatrix;
    unsigned char *direction_matrix;
    double first_A_gap;
    double first_B_gap;
    PyObject *py_Dmatrix, *py_direction_matrix,
	*py_Dmatrix_row, *py_direction_matrix_row;

    double Pscore, Qscore, Dscore, maxscore;
    int Pscore_rint, Qscore_rint, Dscore_rint, maxscore_rint;
    double p, q, d;
    int offset;
    char direction;

    if(!PyArg_ParseTuple(args, "OOddddiiiOd", &py_sequenceA, &py_sequenceB,
			 &open_A, &extend_A, &open_B, &extend_B,
			 &count_first, &global_alignment, &penalize_end_gaps,
			 &py_match_fn, &forgiveness))
	return NULL;
    if(!PySequence_Check(py_sequenceA) || !PySequence_Check(py_sequenceB)) {
	PyErr_SetString(PyExc_TypeError, 
			"py_sequenceA and py_sequenceB should be sequences.");
	return NULL;
    }
    if(!PyCallable_Check(py_match_fn)) {
	PyErr_SetString(PyExc_TypeError, "py_match_fn must be callable.");
	return NULL;
    }
    forgiveness_rint = rint(forgiveness);

    /* Allocate matrices for storing the results and initalize them to
       0. */
    M = PySequence_Length(py_sequenceA);
    N = PySequence_Length(py_sequenceB);
    Dmatrix = Pmatrix = Qmatrix = (double *)NULL;
    direction_matrix = (char *)NULL;
    dsize = M * N * sizeof(*Dmatrix);
    csize = M * N * sizeof(*direction_matrix);
    Dmatrix = (double *)malloc(dsize);
    Pmatrix = (double *)malloc(dsize);
    Qmatrix = (double *)malloc(dsize);
    direction_matrix = (char *)malloc(csize);
    if(!Dmatrix || !Pmatrix || !Qmatrix || !direction_matrix) {
	if(Dmatrix)
	    free(Dmatrix);
	if(Pmatrix)
	    free(Pmatrix);
	if(Qmatrix)
	    free(Qmatrix);
	if(direction_matrix)
	    free(direction_matrix);
	PyErr_SetString(PyExc_MemoryError, "Out of memory");
    }
    memset((void *)Dmatrix, 0, dsize);
    memset((void *)Pmatrix, 0, dsize);
    memset((void *)Qmatrix, 0, dsize);
    memset((void *)direction_matrix, 0, csize);

    /* Cache some commonly used gap penalties */
    first_A_gap = calc_affine_penalty(1, open_A, extend_A, count_first);
    first_B_gap = calc_affine_penalty(1, open_B, extend_B, count_first);

    for(i=0; i<M; i++) {
	for(j=0; j<N; j++) {
	    if(i>0 && i<M-1) {
		offset = (i-1)*N + j;
		p = Pmatrix[offset] + extend_B;
		d = Dmatrix[offset] + first_B_gap;
		Pscore = (p > d) ? p : d;
	    } else {
		/* Bad.  See comments in python code. */
		Pscore = -10000.;
	    }
	    if(j>0 && j<N-1) {
		offset = i*N + j;
		q = Qmatrix[offset] + extend_A;
		d = Dmatrix[offset] + first_A_gap;
		Qscore = (q > d) ? q : d;
	    } else {
		Qscore = -10000.;
	    }

	    /* Calculate the match score. */
	    py_A = PySequence_GetItem(py_sequenceA, i);
	    if(py_A == NULL) {
		i = M; j = N;
		break;
	    }
	    py_B = PySequence_GetItem(py_sequenceB, j);
	    if(py_B == NULL) {
		Py_DECREF(py_A);
		i = M; j = N;
		break;
	    }

	    arglist = Py_BuildValue("(OO)", py_A, py_B);
	    if(arglist == NULL) {
		Py_DECREF(py_A);
		Py_DECREF(py_B);
		i = M; j = N;
		break;
	    }
	    result = PyEval_CallObject(py_match_fn, arglist);
	    Py_DECREF(arglist);
	    Py_DECREF(py_A);
	    Py_DECREF(py_B);
	    if(result == NULL) {
		i = M; j = N;
		break;
	    }
	    Dscore = PyNumber_AsDouble(result);
	    Py_DECREF(result);
	    if(PyErr_Occurred()) {
		i = M; j = N;
		break;
	    }

	    if(i && j) {
		offset = (i-1)*N + j-1;
		Dscore += Dmatrix[offset];
	    } else if(penalize_end_gaps) {
		if(!i)
		    Dscore += calc_affine_penalty(j, open_B, extend_B, 
						  count_first);
		else
		    Dscore += calc_affine_penalty(i, open_A, extend_A,
						  count_first);
	    }

	    if(!global_alignment) {
		Pscore = (Pscore > 0) ? Pscore : 0;
		Qscore = (Qscore > 0) ? Qscore : 0;
		Dscore = (Dscore > 0) ? Dscore : 0;
	    }
	    
	    Pscore_rint = rint(Pscore);
	    Qscore_rint = rint(Qscore);
	    Dscore_rint = rint(Dscore);
	    maxscore = (Pscore > Qscore) ? Pscore : Qscore;
	    maxscore = (maxscore > Dscore) ? maxscore : Dscore;
	    maxscore_rint = rint(maxscore);

	    offset = i*N + j;
	    Pmatrix[offset] = Pscore;
	    Qmatrix[offset] = Qscore;
	    Dmatrix[offset] = maxscore;
	    
	    direction = 0;
	    if(maxscore_rint-Pscore_rint <= forgiveness_rint && i)
		direction |= GO_P;
	    if(maxscore_rint-Qscore_rint <= forgiveness_rint && j)
		direction |= GO_Q;
	    if(maxscore_rint-Dscore_rint <= forgiveness_rint)
		direction |= GO_D;
	    direction_matrix[offset] = direction;
	}
    }

    free(Pmatrix);
    free(Qmatrix);

    if(PyErr_Occurred()) {
	free(Dmatrix);
	free(direction_matrix);
	return NULL;
    }

    /* Now put Dmatrix and direction_matrix into proper Python lists. */
    if((py_Dmatrix = PyList_New(M)) == NULL) {
	free(Dmatrix);
	free(direction_matrix);
	return NULL;
    }
    if((py_direction_matrix = PyList_New(M)) == NULL) {
	Py_DECREF(py_Dmatrix);
	free(Dmatrix);
	free(direction_matrix);
	return NULL;
    }
    for(i=0; i<M; i++) {
	if((py_Dmatrix_row = PyList_New(N)) == NULL)
	    break;
	PyList_SetItem(py_Dmatrix, i, py_Dmatrix_row);
	if((py_direction_matrix_row = PyList_New(N)) == NULL)
	    break;
	PyList_SetItem(py_direction_matrix, i, py_direction_matrix_row);

	for(j=0; j<N; j++) {
	    offset = i*N + j;
	    if((tmp = PyFloat_FromDouble(Dmatrix[offset])) == NULL) {
		i=M; j=N;
		break;
	    }
	    PyList_SetItem(py_Dmatrix_row, j, tmp);
	    if((tmp = PyInt_FromLong((long)direction_matrix[offset])) == NULL){
		i=M; j=N;
		break;
	    }
	    PyList_SetItem(py_direction_matrix_row, j, tmp);
	}
    }
    free(Dmatrix);
    free(direction_matrix);

    if(PyErr_Occurred()) {
	Py_DECREF(py_Dmatrix);
	Py_DECREF(py_direction_matrix);
	return NULL;
    }

    tmp = Py_BuildValue("(OO)", py_Dmatrix, py_direction_matrix);
    Py_DECREF(py_Dmatrix);
    Py_DECREF(py_direction_matrix);
    return tmp;
}


static PyObject *
cfastpairwise__make_score_matrix_faster(self, args)
     PyObject *self;
     PyObject *args;
{
    int i, j;

    char *sequenceA, *sequenceB;
    PyObject *py_match_fn;
    double open_A, extend_A, open_B, extend_B;
    int count_first, global_alignment, penalize_end_gaps;
    double forgiveness;
    int forgiveness_rint;

    PyObject *tmp;

    int M, N;   /* size of matrix */
    int dsize, csize;
    double *Dmatrix, *Pmatrix, *Qmatrix;
    unsigned char *direction_matrix;
    double first_A_gap;
    double first_B_gap;
    PyObject *py_Dmatrix, *py_direction_matrix,
	*py_Dmatrix_row, *py_direction_matrix_row;
    double match, mismatch;

    double Pscore, Qscore, Dscore, maxscore;
    int Pscore_rint, Qscore_rint, Dscore_rint, maxscore_rint;
    double p, q, d;
    int offset;
    char direction;

    if(!PyArg_ParseTuple(args, "ssddddiiiOd", &sequenceA, &sequenceB,
			 &open_A, &extend_A, &open_B, &extend_B,
			 &count_first, &global_alignment, &penalize_end_gaps,
			 &py_match_fn, &forgiveness))
	return NULL;
    if((tmp = PyObject_GetAttrString(py_match_fn, "match")) == NULL)
	return NULL;
    match = PyNumber_AsDouble(tmp);
    Py_DECREF(tmp);
    if(PyErr_Occurred())
	return NULL;
    if((tmp = PyObject_GetAttrString(py_match_fn, "mismatch")) == NULL)
	return NULL;
    mismatch = PyNumber_AsDouble(tmp);
    Py_DECREF(tmp);
    if(PyErr_Occurred())
	return NULL;

    forgiveness_rint = rint(forgiveness);

    /* Allocate matrices for storing the results and initalize them to
       0. */
    M = strlen(sequenceA);
    N = strlen(sequenceB);
    Dmatrix = Pmatrix = Qmatrix = (double *)NULL;
    direction_matrix = (char *)NULL;
    dsize = M * N * sizeof(*Dmatrix);
    csize = M * N * sizeof(*direction_matrix);
    Dmatrix = (double *)malloc(dsize);
    Pmatrix = (double *)malloc(dsize);
    Qmatrix = (double *)malloc(dsize);
    direction_matrix = (char *)malloc(csize);
    if(!Dmatrix || !Pmatrix || !Qmatrix || !direction_matrix) {
	if(Dmatrix)
	    free(Dmatrix);
	if(Pmatrix)
	    free(Pmatrix);
	if(Qmatrix)
	    free(Qmatrix);
	if(direction_matrix)
	    free(direction_matrix);
	PyErr_SetString(PyExc_MemoryError, "Out of memory");
    }
    memset((void *)Dmatrix, 0, dsize);
    memset((void *)Pmatrix, 0, dsize);
    memset((void *)Qmatrix, 0, dsize);
    memset((void *)direction_matrix, 0, csize);

    /* Cache some commonly used gap penalties */
    first_A_gap = calc_affine_penalty(1, open_A, extend_A, count_first);
    first_B_gap = calc_affine_penalty(1, open_B, extend_B, count_first);

    for(i=0; i<M; i++) {
	for(j=0; j<N; j++) {
	    if(i>0 && i<M-1) {
		offset = (i-1)*N + j;
		p = Pmatrix[offset] + extend_B;
		d = Dmatrix[offset] + first_B_gap;
		Pscore = (p > d) ? p : d;
	    } else {
		/* Bad.  See comments in python code. */
		Pscore = -10000.;
	    }
	    if(j>0 && j<N-1) {
		offset = i*N + j;
		q = Qmatrix[offset] + extend_A;
		d = Dmatrix[offset] + first_A_gap;
		Qscore = (q > d) ? q : d;
	    } else {
		Qscore = -10000.;
	    }

	    /* Calculate the match score. */
	    Dscore = (sequenceA[i] == sequenceB[j]) ? match : mismatch;

	    if(i && j) {
		offset = (i-1)*N + j-1;
		Dscore += Dmatrix[offset];
	    } else if(penalize_end_gaps) {
		if(!i)
		    Dscore += calc_affine_penalty(j, open_B, extend_B, 
						  count_first);
		else
		    Dscore += calc_affine_penalty(i, open_A, extend_A,
						  count_first);
	    }

	    if(!global_alignment) {
		Pscore = (Pscore > 0) ? Pscore : 0;
		Qscore = (Qscore > 0) ? Qscore : 0;
		Dscore = (Dscore > 0) ? Dscore : 0;
	    }
	    
	    Pscore_rint = rint(Pscore);
	    Qscore_rint = rint(Qscore);
	    Dscore_rint = rint(Dscore);
	    maxscore = (Pscore > Qscore) ? Pscore : Qscore;
	    maxscore = (maxscore > Dscore) ? maxscore : Dscore;
	    maxscore_rint = rint(maxscore);

	    offset = i*N + j;
	    Pmatrix[offset] = Pscore;
	    Qmatrix[offset] = Qscore;
	    Dmatrix[offset] = maxscore;
	    
	    direction = 0;
	    if(maxscore_rint-Pscore_rint <= forgiveness_rint && i)
		direction |= GO_P;
	    if(maxscore_rint-Qscore_rint <= forgiveness_rint && j)
		direction |= GO_Q;
	    if(maxscore_rint-Dscore_rint <= forgiveness_rint)
		direction |= GO_D;
	    direction_matrix[offset] = direction;
	}
    }

    free(Pmatrix);
    free(Qmatrix);

    if(PyErr_Occurred()) {
	free(Dmatrix);
	free(direction_matrix);
	return NULL;
    }

    /* Now put Dmatrix and direction_matrix into proper Python lists. */
    if((py_Dmatrix = PyList_New(M)) == NULL) {
	free(Dmatrix);
	free(direction_matrix);
	return NULL;
    }
    if((py_direction_matrix = PyList_New(M)) == NULL) {
	Py_DECREF(py_Dmatrix);
	free(Dmatrix);
	free(direction_matrix);
	return NULL;
    }
    for(i=0; i<M; i++) {
	if((py_Dmatrix_row = PyList_New(N)) == NULL)
	    break;
	PyList_SetItem(py_Dmatrix, i, py_Dmatrix_row);
	if((py_direction_matrix_row = PyList_New(N)) == NULL)
	    break;
	PyList_SetItem(py_direction_matrix, i, py_direction_matrix_row);

	for(j=0; j<N; j++) {
	    offset = i*N + j;
	    if((tmp = PyFloat_FromDouble(Dmatrix[offset])) == NULL) {
		i=M; j=N;
		break;
	    }
	    PyList_SetItem(py_Dmatrix_row, j, tmp);
	    if((tmp = PyInt_FromLong((long)direction_matrix[offset])) == NULL){
		i=M; j=N;
		break;
	    }
	    PyList_SetItem(py_direction_matrix_row, j, tmp);
	}
    }
    free(Dmatrix);
    free(direction_matrix);

    if(PyErr_Occurred()) {
	Py_DECREF(py_Dmatrix);
	Py_DECREF(py_direction_matrix);
	return NULL;
    }

    tmp = Py_BuildValue("(OO)", py_Dmatrix, py_direction_matrix);
    Py_DECREF(py_Dmatrix);
    Py_DECREF(py_direction_matrix);
    return tmp;
}



/* Module definition stuff */

static PyMethodDef cfastpairwiseMethods[] = {
    {"_make_score_matrix", cfastpairwise__make_score_matrix, METH_VARARGS},
    {"_make_score_matrix_faster", cfastpairwise__make_score_matrix_faster, 
     METH_VARARGS},
    {NULL, NULL}
};

static char cfastpairwise__doc__[] =
"XXX document here\n\
\n\
";

void initcfastpairwise()
{
    (void) Py_InitModule3("cfastpairwise", cfastpairwiseMethods, 
			  cfastpairwise__doc__);
}
