#ifndef DARRAY_H
#define DARRAY_H

#include <Python.h>

typedef struct {
    Py_ssize_t position;
    float score;
} Hit;

typedef struct {
    Py_ssize_t used;
    Py_ssize_t size;
    Hit* data;
} DArray;

int init_darray(DArray *a, Py_ssize_t initial_size);
int insert_in_darray(DArray *a, Py_ssize_t pos, float score);
void free_darray(DArray *a);

#endif