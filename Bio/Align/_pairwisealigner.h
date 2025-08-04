/* Copyright 2018-2019 by Michiel de Hoon.  All rights reserved.
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */


#define PY_SSIZE_T_CLEAN
#include "Python.h"

#define HORIZONTAL 0x1
#define VERTICAL 0x2
#define DIAGONAL 0x4

typedef enum {NeedlemanWunschSmithWaterman,
              Gotoh,
              WatermanSmithBeyer,
              FOGSAA,
              Unknown} Algorithm;

typedef enum {Global, Local, FOGSAA_Mode} Mode;

typedef struct {
    PyObject_HEAD
    Mode mode;
    Algorithm algorithm;
    double match;
    double mismatch;
    double epsilon;
    double open_internal_insertion_score;
    double extend_internal_insertion_score;
    double open_left_insertion_score;
    double extend_left_insertion_score;
    double open_right_insertion_score;
    double extend_right_insertion_score;
    double open_internal_deletion_score;
    double extend_internal_deletion_score;
    double open_left_deletion_score;
    double extend_left_deletion_score;
    double open_right_deletion_score;
    double extend_right_deletion_score;
    bool open_internal_insertion_score_set;
    bool extend_internal_insertion_score_set;
    bool open_left_insertion_score_set;
    bool extend_left_insertion_score_set;
    bool open_right_insertion_score_set;
    bool extend_right_insertion_score_set;
    bool open_internal_deletion_score_set;
    bool extend_internal_deletion_score_set;
    bool open_left_deletion_score_set;
    bool extend_left_deletion_score_set;
    bool open_right_deletion_score_set;
    bool extend_right_deletion_score_set;
    PyObject* insertion_score_function;
    PyObject* deletion_score_function;
    Py_buffer substitution_matrix;
    PyObject* alphabet;
    int wildcard;
} Aligner;
