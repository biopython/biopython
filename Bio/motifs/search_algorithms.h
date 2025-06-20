#ifndef SEARCH_ALGORITHMS_H
#define SEARCH_ALGORITHMS_H

#include <Python.h>
#include "darray.h"

typedef enum {
    LOOKAHEAD,
    PERMUTED_LOOKAHEAD,
    SUPERALPHABET,
    UNKNOWN_ALGORITHM
} AlgorithmType;

typedef struct {
    int index;
    float priority;
} PermEntry;

AlgorithmType parse_algorithm(const char* name);

int search_lookahead(const char sequence[], Py_ssize_t s, double* matrix, 
                    Py_ssize_t m, float threshold, DArray* scores, const int* base_lookup);

int compare_permutation(const void* a, const void* b);
void compute_permutation(double* matrix, Py_ssize_t m, const float bkg[4], int* perm);
int search_permuted_lookahead(const char sequence[], Py_ssize_t s, double* matrix, 
                            Py_ssize_t m, float threshold, DArray* scores, const int* base_lookup);

int compute_mhat(const double* mat, Py_ssize_t m, Py_ssize_t q, double* mhat);
int search_superalphabet(const char sequence[], Py_ssize_t s, double* matrix, 
                        Py_ssize_t m, float threshold, DArray* scores, Py_ssize_t q, const int* base_lookup);

#endif 