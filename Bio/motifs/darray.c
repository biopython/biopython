#include "darray.h"

int init_darray(DArray *a, Py_ssize_t initial_size) {
    a->data = PyMem_Malloc(initial_size * sizeof(Hit));
    if (!a->data) return -1;
    a->used = 0;
    a->size = initial_size;
    return 0;
}
 
int insert_in_darray(DArray *a, Py_ssize_t pos, float score) {
  
    if (a->used == a->size) {
        //if the array is full
        Py_ssize_t newSize = a->size * 3 / 2 + 8; 
        
        if (newSize < a->size) {
            PyErr_SetString(PyExc_TypeError, "Results Array overflow");
            return -1;
        }
  
        Hit *tmp = realloc(a->data, newSize * sizeof(Hit));
        if (!tmp) {
            PyErr_NoMemory();// Sets the exception and returns NULL
            return -1;
        }
        a->data = tmp;
        a->size = newSize;
    }
  
    a->data[a->used].position = pos;
    a->data[a->used].score = score;
    a->used++;
    return 0;
}

void free_darray(DArray *a) {
    PyMem_Free(a->data);
    a->data = NULL;
    a->used = a->size = 0;
}
