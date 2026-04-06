/*  Copyright (C) 2025, Joaquim Calvera
 *
 *  This file is part of the Biopython distribution and governed by your
 *  choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 *  Please see the LICENSE file that should have been included as part of this
 *  package.
 *
 *  This file is based in part on 'pwm.c' (Copyright 2009–2020, Michiel de Hoon)
 *  and has been modified to support additional motif search algorithms and 
 *  outputs. 
 */
#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <math.h>
#include "search_algorithms.h"


/* Validate that the Position-Weight Matrix satisfies the necessary requirements*/
static int matrix_converter(PyObject* object, void* address)
{
    const int flags = PyBUF_C_CONTIGUOUS | PyBUF_FORMAT;
    char datatype;
    Py_buffer* view = address;
     
    if (object == NULL) goto exit;
    //PWM array must be C-contigous 
    if (PyObject_GetBuffer(object, view, flags) == -1) {
        PyErr_SetString(PyExc_RuntimeError,
                        "position-weight matrix is not an array");
        return 0;
    }
    datatype = view->format[0];
    switch (datatype) {
        case '@':
        case '=':
        case '<':
        case '>':
        case '!': datatype = view->format[1]; break;
        default: break;
    }
    //PWM must be a doubles array
    if (datatype != 'd') {
        PyErr_Format(PyExc_RuntimeError,
            "position-weight matrix data format incorrect "
            "('%c', expected 'd')", datatype);
        goto exit;
    }
    //PWM must be 2 dimensional
    if (view->ndim != 2) {
        PyErr_Format(PyExc_RuntimeError,
            "position-weight matrix has incorrect rank (%d expected 2)",
            view->ndim);
        goto exit;
    }
    //PWM must have 4 columns, one for nucleotide (A,C,G,T)
    if (view->shape[1] != 4) {
        PyErr_Format(PyExc_RuntimeError,
            "position-weight matrix should have four columns "
            "(%zd columns found)", view->shape[1]);
        goto exit;
    }
    return Py_CLEANUP_SUPPORTED;

exit:
    PyBuffer_Release(view);
    return 0;
}
 
  
static char search__doc__[] =
"search(sequence, matrix, threshold, algorithm, q)\n"
"\n"
"This function performs efficient motif search using either the 'lookahead',\n"
"'permuted_lookahead', or 'superalphabet' algorithms. It supports dynamic\n"
"thresholding and optimizations to discard unlikely matches early.\n"
"\n"
"Arguments:\n"
"  sequence  (str)    : DNA sequence to scan.\n"
"  matrix    (ndarray): Position-specific scoring matrix.\n"
"  threshold (float)  : Minimum score required for a match.\n"
"  algorithm (str)    : Search method: 'lookahead', 'permuted', 'superalphabet'.\n"
"  q         (int)    : Size of q-tuples for 'superalphabet' method (optional).\n"
"\n"
"Returns:\n"
"  2 Python lists: one for positions and another for scores.\n";
  
/*
    Python entry point
*/
static PyObject* py_search(PyObject* self, PyObject* args, PyObject* kwargs) 
{
    const char* sequence;          // DNA sequence input to search
    const char* algorithm;         // Algorithm selected
    Py_ssize_t s;                  // Sequence length
    float threshold;               // Threshold
    Py_buffer matrix;              // Position-weight matrix input
    DArray scores;                 // Dynamic array where matchs will be stored
    int status = -1;
    Py_ssize_t q;                  // Length of super-alphabet symbols (if super-alphabet algorithm selected)

    // Accepted parameter keywords
    static char* kwlist[] = {"sequence", "matrix", "threshold", "algorithm", "q", NULL};
    
    // Parse the parameters of the function into local variables
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "y#O&fyn", kwlist,
                                        &sequence, &s,
                                        matrix_converter, &matrix,
                                        &threshold,
                                        &algorithm,
                                        &q)) return NULL;
 
    Py_ssize_t m = matrix.shape[0];
    // Initialize scores (solution array)
    status = init_darray(&scores, 1000);
    if (status < 0) {
        free_darray(&scores);    
        PyBuffer_Release(&matrix);
        return PyErr_NoMemory();
    }

    // Lookup array to quickly get the value of a base of the sequence.
    int base_lookup[256];
    for (Py_ssize_t i = 0; i < 256; i++) base_lookup[i] = -1;
    base_lookup['A'] = base_lookup['a'] = 0;
    base_lookup['C'] = base_lookup['c'] = 1;
    base_lookup['G'] = base_lookup['g'] = 2;
    base_lookup['T'] = base_lookup['t'] = 3;

    AlgorithmType alg = parse_algorithm(algorithm);
    switch (alg) {
        case LOOKAHEAD:
            status = search_lookahead(sequence, s, matrix.buf, m, threshold, &scores, base_lookup); break;
        case PERMUTED_LOOKAHEAD:
            status = search_permuted_lookahead(sequence, s, matrix.buf, m, threshold, &scores, base_lookup); break;
        case SUPERALPHABET:
            status = search_superalphabet(sequence, s, matrix.buf, m, threshold, &scores, q, base_lookup); break;
        case UNKNOWN_ALGORITHM:
        default:
            status = -1;
            PyErr_Format(PyExc_ValueError, "Unknown algorithm: '%s'. Valid options are: 'lookahead', 'permuted', 'superalphabet'", algorithm);
    }

    if (status < 0) {
        free_darray(&scores);    
        PyBuffer_Release(&matrix);
        return NULL;
    }

   // Create and return 2 Python lists. One for positions and another for scores.
    PyObject *positions = PyList_New(scores.used);
    PyObject *scores_list = PyList_New(scores.used);

    if (!positions || !scores_list) {
        Py_XDECREF(positions);
        Py_XDECREF(scores_list);
        free_darray(&scores);
        PyBuffer_Release(&matrix);
        return PyErr_NoMemory();
    }

    for (Py_ssize_t i = 0; i < scores.used; i++) {
        PyObject *pos = PyLong_FromLong(scores.data[i].position);
        PyObject *score = PyFloat_FromDouble(scores.data[i].score);

        if (!pos || !score) {
            Py_XDECREF(pos);
            Py_XDECREF(score);
            Py_DECREF(positions);
            Py_DECREF(scores_list);
            free_darray(&scores);
            PyBuffer_Release(&matrix);
            return PyErr_NoMemory();
        }

        PyList_SET_ITEM(positions, i, pos);      
        PyList_SET_ITEM(scores_list, i, score);  
    }

    PyObject *result = PyTuple_New(2);
    if (!result) {
        Py_DECREF(positions);
        Py_DECREF(scores_list);
        free_darray(&scores);
        PyBuffer_Release(&matrix);
        return PyErr_NoMemory();
    }

    PyTuple_SET_ITEM(result, 0, positions);     
    PyTuple_SET_ITEM(result, 1, scores_list);   
    
    free_darray(&scores);    
    PyBuffer_Release(&matrix);
    return result;
}
  
/*Public methods in this module*/
static struct PyMethodDef methods[] = {
    {"search",              // Name of the function py_search visible from python
    (PyCFunction)py_search, 
    METH_VARARGS | METH_KEYWORDS,  
    PyDoc_STR(search__doc__),  
    },
    {NULL, NULL, 0, NULL} //sentinel
};
  
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,    
    "_searchmodule",   // Name of the module 
    PyDoc_STR("Motif search module."), 
    -1,                       
                               
    methods,                   
    NULL,               
    NULL,
    NULL,
    NULL
};
  
/* Initialization function*/
PyMODINIT_FUNC* PyInit__searchmodule(void)  
{
    return PyModule_Create(&moduledef);
}
  