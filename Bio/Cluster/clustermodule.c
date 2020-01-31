#include "Python.h"
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "cluster.h"


/* ========================================================================= */
/* -- Helper routines ------------------------------------------------------ */
/* ========================================================================= */

static char
extract_single_character(PyObject* object, const char variable[],
                         const char allowed[])
{
    Py_UCS4 ch;
    Py_ssize_t n;
    if (!PyUnicode_Check(object)) {
        PyErr_Format(PyExc_ValueError, "%s should be a string", variable);
        return 0;
    }
    if (PyUnicode_READY(object) == -1) return 0;
    n = PyUnicode_GET_LENGTH(object);
    if (n != 1) {
        PyErr_Format(PyExc_ValueError,
                     "%s should be a single character", variable);
        return 0;
    }
    ch = PyUnicode_READ_CHAR(object, 0);
    if (ch < 128) {
        const char c = ch;
        if (strchr(allowed, c)) return c;
    }
    PyErr_Format(PyExc_ValueError,
                 "unknown %s function specified (should be one of '%s')",
                 variable, allowed);
    return 0;
}

static int
distance_converter(PyObject* object, void* pointer)
{
    char c;

    c = extract_single_character(object, "dist", "ebcauxsk");
    if (c == 0) return 0;
    *((char*)pointer) = c;
    return 1;
}

static int
method_treecluster_converter(PyObject* object, void* pointer)
{
    char c;

    c = extract_single_character(object, "method", "csma");
    if (c == 0) return 0;
    *((char*)pointer) = c;
    return 1;
}

static int
method_kcluster_converter(PyObject* object, void* pointer)
{
    char c;

    c = extract_single_character(object, "method", "am");
    if (c == 0) return 0;
    *((char*)pointer) = c;
    return 1;
}

static int
method_clusterdistance_converter(PyObject* object, void* pointer)
{
    char c;

    c = extract_single_character(object, "method", "amsxv");
    if (c == 0) return 0;
    *((char*)pointer) = c;
    return 1;
}

/* -- data ----------------------------------------------------------------- */

typedef struct {
    int nrows;
    int ncols;
    double** values;
    Py_buffer view;
} Data;

static int
data_converter(PyObject* object, void* pointer)
{
    Data* data = pointer;
    int nrows;
    int ncols;
    int i;
    double** values = data->values;
    Py_buffer* view = &data->view;
    const char* p;
    Py_ssize_t stride;
    const int flag = PyBUF_ND | PyBUF_STRIDES;

    if (object == NULL) goto exit;
    if (object == Py_None) return 1;

    if (PyObject_GetBuffer(object, view, flag) == -1) {
        PyErr_SetString(PyExc_RuntimeError,
                        "data matrix has unexpected format.");
        return 0;
    }

    if (view->ndim != 2) {
        PyErr_Format(PyExc_RuntimeError,
                     "data matrix has incorrect rank %d (expected 2)",
                     view->ndim);
        goto exit;
    }
    if (view->itemsize != sizeof(double)) {
        PyErr_SetString(PyExc_RuntimeError,
                        "data matrix has incorrect data type");
        goto exit;
    }
    nrows = (int) view->shape[0];
    ncols = (int) view->shape[1];
    if (nrows != view->shape[0] || ncols != view->shape[1]) {
        PyErr_Format(PyExc_ValueError,
            "data matrix is too large (dimensions = %zd x %zd)",
            view->shape[0], view->shape[1]);
        goto exit;
    }
    if (nrows < 1 || ncols < 1) {
        PyErr_SetString(PyExc_ValueError, "data matrix is empty");
        goto exit;
    }
    stride = view->strides[0];
    if (view->strides[1] != view->itemsize) {
        PyErr_SetString(PyExc_RuntimeError, "data is not contiguous");
        goto exit;
    }
    values = PyMem_Malloc(nrows*sizeof(double*));
    if (!values) {
        PyErr_NoMemory();
        goto exit;
    }
    for (i = 0, p = view->buf; i < nrows; i++, p += stride)
        values[i] = (double*)p;
    data->values = values;
    data->nrows = nrows;
    data->ncols = ncols;
    return Py_CLEANUP_SUPPORTED;

exit:
    if (values) PyMem_Free(values);
    PyBuffer_Release(view);
    return 0;
}

/* -- mask ----------------------------------------------------------------- */

typedef struct {
    int** values;
    Py_buffer view;
} Mask;

static int
mask_converter(PyObject* object, void* pointer)
{
    Mask* mask = pointer;
    int nrows;
    int ncols;
    int i;
    int** values = mask->values;
    Py_buffer* view = &mask->view;
    const char* p;
    Py_ssize_t stride;
    const int flag = PyBUF_ND | PyBUF_STRIDES;

    if (object == NULL) goto exit;
    if (object == Py_None) return 1;

    if (PyObject_GetBuffer(object, view, flag) == -1) {
        PyErr_SetString(PyExc_RuntimeError, "mask has unexpected format.");
        return 0;
    }

    if (view->ndim != 2) {
        PyErr_Format(PyExc_ValueError,
                     "mask has incorrect rank %d (expected 2)", view->ndim);
        goto exit;
    }
    if (view->itemsize != sizeof(int)) {
        PyErr_SetString(PyExc_RuntimeError, "mask has incorrect data type");
        goto exit;
    }
    nrows = (int) view->shape[0];
    ncols = (int) view->shape[1];
    if (nrows != view->shape[0] || ncols != view->shape[1]) {
        PyErr_Format(PyExc_ValueError,
                     "mask is too large (dimensions = %zd x %zd)",
                     view->shape[0], view->shape[1]);
        goto exit;
    }
    stride = view->strides[0];
    if (view->strides[1] != view->itemsize) {
        PyErr_SetString(PyExc_RuntimeError, "mask is not contiguous");
        goto exit;
    }
    values = PyMem_Malloc(nrows*sizeof(int*));
    if (!values) {
        PyErr_NoMemory();
        goto exit;
    }
    for (i = 0, p = view->buf; i < nrows; i++, p += stride)
        values[i] = (int*)p;
    mask->values = values;
    return Py_CLEANUP_SUPPORTED;

exit:
    if (values) PyMem_Free(values);
    PyBuffer_Release(view);
    return 0;
}

/* -- 1d array ------------------------------------------------------------- */

static int
vector_converter(PyObject* object, void* pointer)
{
    Py_buffer* view = pointer;
    int ndata;
    const int flag = PyBUF_ND | PyBUF_C_CONTIGUOUS;

    if (object == NULL) goto exit;

    if (PyObject_GetBuffer(object, view, flag) == -1) {
        PyErr_SetString(PyExc_RuntimeError, "unexpected format.");
        return 0;
    }

    if (view->ndim != 1) {
        PyErr_Format(PyExc_ValueError, "incorrect rank %d (expected 1)",
                     view->ndim);
        goto exit;
    }
    if (view->itemsize != sizeof(double)) {
        PyErr_SetString(PyExc_RuntimeError, "array has incorrect data type");
        goto exit;
    }
    ndata = (int) view->shape[0];
    if (ndata != view->shape[0]) {
        PyErr_Format(PyExc_ValueError,
                     "array is too large (size = %zd)", view->shape[0]);
        goto exit;
    }
    return Py_CLEANUP_SUPPORTED;

exit:
    PyBuffer_Release(view);
    return 0;
}

static int
vector_none_converter(PyObject* object, void* pointer)
{
    if (object == Py_None) return 1;
    return vector_converter(object, pointer);
}

/* -- clusterid ------------------------------------------------------------ */

static int
check_clusterid(Py_buffer clusterid, int nitems) {
    int i, j;
    int *p = clusterid.buf;
    int nclusters = 0;
    int* number;

    if (nitems != clusterid.shape[0]) {
        PyErr_Format(PyExc_ValueError, "incorrect size (%zd, expected %d)",
                     clusterid.shape[0], nitems);
        return 0;
    }
    for (i = 0; i < nitems; i++) {
        j = p[i];
        if (j > nclusters) nclusters = j;
        if (j < 0) {
            PyErr_SetString(PyExc_ValueError, "negative cluster number found");
            return 0;
        }
    }
    nclusters++;
    /* -- Count the number of items in each cluster --------------------- */
    number = calloc(nclusters, sizeof(int));
    if (!number) {
        PyErr_NoMemory();
        return 0;
    }
    for (i = 0; i < nitems; i++) {
        j = p[i];
        number[j]++;
    }
    for (j = 0; j < nclusters; j++) if (number[j] == 0) break;
    PyMem_Free(number);
    if (j < nclusters) {
        PyErr_Format(PyExc_ValueError, "cluster %d is empty", j);
        return 0;
    }
    return nclusters;
}

/* -- distance ----------------------------------------------------------- */

typedef struct {
    int n;
    double** values;
    Py_buffer* views;
    Py_buffer view;
} Distancematrix;

static int
_convert_list_to_distancematrix(PyObject* list, Distancematrix* distances)
{
    int i;
    double** values;
    Py_buffer* view;
    Py_buffer* views;
    const int flag = PyBUF_ND | PyBUF_C_CONTIGUOUS;
    const int n = (int) PyList_GET_SIZE(list);

    if (n != PyList_GET_SIZE(list)) {
        PyErr_SetString(PyExc_ValueError, "distance matrix is too large");
        return 0;
    }
    values = PyMem_Malloc(n*sizeof(double*));
    if (!values) {
        PyErr_NoMemory();
        return 0;
    }
    distances->values = values;
    views = PyMem_Malloc(n*sizeof(Py_buffer));
    if (!views) {
        PyErr_NoMemory();
        return 0;
    }
    view = views;
    for (i = 0; i < n; i++, view++) {
        PyObject* item = PyList_GET_ITEM(list, i);
        view->len = -1;
        if (PyObject_GetBuffer(item, view, flag) == -1) {
            PyErr_Format(PyExc_RuntimeError, "failed to parse row %d.", i);
            view--;
            break;
        }
        if (view->ndim != 1) {
            PyErr_Format(PyExc_ValueError,
                         "row %d has incorrect rank (%d expected 1)",
                         i, view->ndim);
            break;
        }
        if (view->itemsize != sizeof(double)) {
            PyErr_Format(PyExc_RuntimeError,
                         "row %d has incorrect data type", i);
            break;
        }
        if (view->shape[0] != i) {
            PyErr_Format(PyExc_RuntimeError,
                         "row %d has incorrect size %zd (expected %d)",
                         i, view->shape[0], i);
            break;
        }
        values[i] = view->buf;
    }
    if (i < n) {
        for ( ; view >= views; view--) PyBuffer_Release(view);
        PyMem_Free(views);
        return 0;
    }
    distances->n = n;
    distances->view.len = 0;
    distances->views = views;
    distances->values = values;
    return 1;
}

static int
_convert_array_to_distancematrix(PyObject* array, Distancematrix* distances)
{
    int i;
    int n;
    double** values;
    double* p;
    Py_buffer* view = &distances->view;
    const int flag = PyBUF_ND | PyBUF_C_CONTIGUOUS;

    if (PyObject_GetBuffer(array, view, flag) == -1) {
        PyErr_SetString(PyExc_RuntimeError,
                        "distance matrix has unexpected format.");
        return 0;
    }

    if (view->len == 0) {
        PyBuffer_Release(view);
        PyErr_SetString(PyExc_ValueError, "distance matrix is empty");
        return 0;
    }
    if (view->itemsize != sizeof(double)) {
        PyErr_SetString(PyExc_RuntimeError,
                        "distance matrix has an incorrect data type");
        return 0;
    }
    if (view->ndim == 1) {
        int m = (int) view->shape[0];
        if (m != view->shape[0]) {
            PyErr_Format(PyExc_ValueError,
                         "distance matrix is too large (size = %zd)",
                         view->shape[0]);
            return 0;
        }
        n = (int)(1+sqrt(1+8*m)/2); /* rounds to (1+sqrt(1+8*m))/2 */
        if (n*n-n != 2 * m) {
            PyErr_SetString(PyExc_ValueError,
                            "distance matrix has unexpected size.");
            return 0;
        }
        distances->n = n;
        values = PyMem_Malloc(n*sizeof(double*));
        if (!values) {
            PyErr_NoMemory();
            return 0;
        }
        distances->values = values;
        for (p = view->buf, i = 0; i < n; p += i, i++) values[i] = p;
    }
    else if (view->ndim == 2) {
        n = (int) view->shape[0];
        if (n != view->shape[0]) {
            PyErr_Format(PyExc_ValueError,
                         "distance matrix is too large (size = %zd)",
                         view->shape[0]);
            return 0;
        }
        distances->n = n;
        if (view->shape[1] != n) {
            PyErr_SetString(PyExc_ValueError,
                            "distance matrix is not square.");
            return 0;
        }
        values = PyMem_Malloc(n*sizeof(double*));
        if (!values) {
            PyErr_NoMemory();
            return 0;
        }
        distances->values = values;
        for (p = view->buf, i = 0; i < n; p += n, i++) values[i] = p;
    }
    else {
        PyErr_Format(PyExc_ValueError,
                     "distance matrix has incorrect rank %d (expected 1 or 2)",
                     view->ndim);
        return 0;
    }
    return 1;
}

static int
distancematrix_converter(PyObject* argument, void* pointer)
{
    Distancematrix* distances = pointer;
    double** values;

    if (argument == NULL) goto exit;
    if (argument == Py_None) return 1;
    if (PyList_Check(argument)) {
        if (_convert_list_to_distancematrix(argument, distances))
            return Py_CLEANUP_SUPPORTED;
    }
    else {
        if (_convert_array_to_distancematrix(argument, distances))
            return Py_CLEANUP_SUPPORTED;
    }

exit:
    values = distances->values;
    if (values == NULL) return 0;
    else {
        int i;
        const int n = distances->n;
        Py_buffer* views = distances->views;
        if (views) {
            for (i = 0; i < n; i++) PyBuffer_Release(&views[i]);
            PyMem_Free(views);
        }
        else if (distances->view.len) {
            PyBuffer_Release(&distances->view);
        }
        PyMem_Free(values);
    }
    return 0;
}

/* -- celldata ------------------------------------------------------------- */

typedef struct {
    int nx;
    int ny;
    int nz;
    double*** values;
    Py_buffer view;
} Celldata;

static int
celldata_converter(PyObject* argument, void* pointer)
{
    int i, n;
    double* p;
    Celldata* celldata = pointer;
    double*** ppp = celldata->values;
    double** pp = ppp ? ppp[0] : NULL;
    int nx;
    int ny;
    int nz;
    Py_buffer* view = &celldata->view;
    const int flag = PyBUF_ND | PyBUF_C_CONTIGUOUS;

    if (argument == NULL) goto exit;

    if (PyObject_GetBuffer(argument, view, flag) == -1) {
        PyErr_SetString(PyExc_RuntimeError,
                        "celldata array has unexpected format.");
        return 0;
    }

    nx = (int) view->shape[0];
    ny = (int) view->shape[1];
    nz = (int) view->shape[2];
    if (nx != view->shape[0] || ny != view->shape[1] || nz != view->shape[2]) {
        PyErr_SetString(PyExc_RuntimeError, "celldata array too large");
        goto exit;
    }
    if (view->itemsize != sizeof(double)) {
        PyErr_SetString(PyExc_RuntimeError,
                        "celldata array has incorrect data type");
        goto exit;
    }
    pp = PyMem_Malloc(nx*ny*sizeof(double*));
    ppp = PyMem_Malloc(nx*sizeof(double**));
    if (!pp || !ppp) {
        PyErr_NoMemory();
        goto exit;
    }
    p = view->buf;
    n = nx * ny;
    for (i = 0; i < n; i++, p += nz) pp[i] = p;
    for (i = 0; i < nx; i++, pp += ny) ppp[i] = pp;
    celldata->values = ppp;
    celldata->nx = nx;
    celldata->ny = ny;
    celldata->nz = nz;
    return Py_CLEANUP_SUPPORTED;

exit:
    if (pp) PyMem_Free(pp);
    if (ppp) PyMem_Free(ppp);
    PyBuffer_Release(view);
    return 0;
}


/* -- index ---------------------------------------------------------------- */

static int
index_converter(PyObject* argument, void* pointer)
{
    Py_buffer* view = pointer;
    int n;
    const int flag = PyBUF_ND | PyBUF_C_CONTIGUOUS;

    if (argument == NULL) goto exit;

    if (PyObject_GetBuffer(argument, view, flag) == -1) {
        PyErr_SetString(PyExc_RuntimeError, "unexpected format.");
        return 0;
    }

    if (view->ndim != 1) {
        PyErr_Format(PyExc_ValueError, "incorrect rank %d (expected 1)",
                     view->ndim);
        goto exit;
    }
    if (view->itemsize != sizeof(int)) {
        PyErr_SetString(PyExc_RuntimeError,
                        "argument has incorrect data type");
        goto exit;
    }
    n = (int) view->shape[0];
    if (n != view->shape[0]) {
        PyErr_Format(PyExc_ValueError,
            "array size is too large (size = %zd)", view->shape[0]);
        goto exit;
    }
    return Py_CLEANUP_SUPPORTED;

exit:
    PyBuffer_Release(view);
    return 0;
}

/* -- index2d ------------------------------------------------------------- */

static int
index2d_converter(PyObject* argument, void* pointer)
{
    Py_buffer* view = pointer;
    int n;
    const int flag = PyBUF_ND | PyBUF_C_CONTIGUOUS;

    if (argument == NULL) goto exit;

    if (PyObject_GetBuffer(argument, view, flag) == -1) {
        PyErr_SetString(PyExc_RuntimeError, "unexpected format.");
        return 0;
    }

    if (view->ndim != 2) {
        PyErr_Format(PyExc_ValueError, "incorrect rank %d (expected 2)",
                     view->ndim);
        goto exit;
    }
    if (view->itemsize != sizeof(int)) {
        PyErr_SetString(PyExc_RuntimeError,
                        "argument has incorrect data type");
        goto exit;
    }
    n = (int) view->shape[0];
    if (n != view->shape[0]) {
        PyErr_Format(PyExc_ValueError,
            "array size is too large (size = %zd)", view->shape[0]);
        goto exit;
    }
    if (view->shape[1] != 2) {
        PyErr_Format(PyExc_ValueError,
            "array has %zd columns (expected 2)", view->shape[1]);
        goto exit;
    }
    return Py_CLEANUP_SUPPORTED;

exit:
    PyBuffer_Release(view);
    return 0;
}

/* ========================================================================= */
/* -- Classes -------------------------------------------------------------- */
/* ========================================================================= */

typedef struct {
    PyObject_HEAD
    Node node;
} PyNode;

static int
PyNode_init(PyNode *self, PyObject *args, PyObject *kwds)
{
    int left, right;
    double distance = 0.0;
    static char *kwlist[] = {"left", "right", "distance", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ii|d", kwlist,
                                      &left, &right, &distance))
        return -1;
    self->node.left = left;
    self->node.right = right;
    self->node.distance = distance;
    return 0;
}

static PyObject*
PyNode_repr(PyNode* self)
{
    char string[64];

    sprintf(string, "(%d, %d): %g",
                   self->node.left, self->node.right, self->node.distance);
    return PyUnicode_FromString(string);
}

static char PyNode_left__doc__[] =
"integer representing the first member of this node";

static PyObject*
PyNode_getleft(PyNode* self, void* closure)
{
    int left = self->node.left;

    return PyLong_FromLong((long)left);
}

static int
PyNode_setleft(PyNode* self, PyObject* value, void* closure)
{
    long left = PyLong_AsLong(value);

    if (PyErr_Occurred()) return -1;
    self->node.left = (int) left;
    return 0;
}

static char PyNode_right__doc__[] =
"integer representing the second member of this node";

static PyObject*
PyNode_getright(PyNode* self, void* closure)
{
    int right = self->node.right;

    return PyLong_FromLong((long)right);
}

static int
PyNode_setright(PyNode* self, PyObject* value, void* closure)
{
    long right = PyLong_AsLong(value);

    if (PyErr_Occurred()) return -1;
    self->node.right = (int) right;
    return 0;
}

static PyObject*
PyNode_getdistance(PyNode* self, void* closure)
{
    return PyFloat_FromDouble(self->node.distance);
}

static int
PyNode_setdistance(PyNode* self, PyObject* value, void* closure)
{
    const double distance = PyFloat_AsDouble(value);

    if (PyErr_Occurred()) return -1;
    self->node.distance = distance;
    return 0;
}

static char PyNode_distance__doc__[] =
"the distance between the two members of this node\n";

static PyGetSetDef PyNode_getset[] = {
    {"left",
     (getter)PyNode_getleft,
     (setter)PyNode_setleft,
     PyNode_left__doc__, NULL},
    {"right",
     (getter)PyNode_getright,
     (setter)PyNode_setright,
     PyNode_right__doc__, NULL},
    {"distance",
     (getter)PyNode_getdistance,
     (setter)PyNode_setdistance,
     PyNode_distance__doc__, NULL},
    {NULL}  /* Sentinel */
};

static char PyNode_doc[] =
"A Node object describes a single node in a hierarchical clustering tree.\n"
"The integer attributes 'left' and 'right' represent the two members that\n"
"make up this node; the floating point attribute 'distance' contains the\n"
"distance between the two members of this node.\n";

static PyTypeObject PyNodeType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_cluster.Node",           /* tp_name */
    sizeof(PyNode),            /* tp_basicsize */
    0,                         /* tp_itemsize */
    0,                         /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_compare */
    (reprfunc)PyNode_repr,     /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,          /*tp_flags*/
    PyNode_doc,                /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    0,                         /* tp_methods */
    0,                         /* tp_members */
    PyNode_getset,             /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)PyNode_init,     /* tp_init */
};

typedef struct {
    PyObject_HEAD
    Node* nodes;
    int n;
} PyTree;

static void
PyTree_dealloc(PyTree* self)
{
    if (self->n) PyMem_Free(self->nodes);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject*
PyTree_new(PyTypeObject *type, PyObject* args, PyObject* kwds)
{
    int i, j;
    int n;
    Node* nodes;
    PyObject* arg = NULL;
    int* flag;
    PyTree* self;

    self = (PyTree *)type->tp_alloc(type, 0);
    if (!self) return NULL;

    if (!PyArg_ParseTuple(args, "|O", &arg)) {
        Py_DECREF(self);
        return NULL;
    }

    if (arg == NULL) {
        self->n = 0;
        self->nodes = NULL;
        return (PyObject*)self;
    }

    if (!PyList_Check(arg)) {
        Py_DECREF(self);
        PyErr_SetString(PyExc_TypeError,
                        "Argument should be a list of Node objects");
        return NULL;
    }

    n = (int) PyList_GET_SIZE(arg);
    if (n != PyList_GET_SIZE(arg)) {
        Py_DECREF(self);
        PyErr_Format(PyExc_ValueError,
                     "List is too large (size = %zd)", PyList_GET_SIZE(arg));
        return NULL;
    }
    if (n < 1) {
        Py_DECREF(self);
        PyErr_SetString(PyExc_ValueError, "List is empty");
        return NULL;
    }
    nodes = PyMem_Malloc(n*sizeof(Node));
    if (!nodes) {
        Py_DECREF(self);
        return PyErr_NoMemory();
    }
    for (i = 0; i < n; i++) {
        PyNode* p;
        PyObject* row = PyList_GET_ITEM(arg, i);
        if (!PyType_IsSubtype(Py_TYPE(row), &PyNodeType)) {
            PyMem_Free(nodes);
            Py_DECREF(self);
            PyErr_Format(PyExc_TypeError,
                         "Row %d in list is not a Node object", i);
            return NULL;
        }
        p = (PyNode*)row;
        nodes[i] = p->node;
    }
    /* --- Check if this is a bona fide tree ------------------------------- */
    flag = PyMem_Malloc((2*n+1)*sizeof(int));
    if (!flag) {
        PyMem_Free(nodes);
        Py_DECREF(self);
        return PyErr_NoMemory();
    }
    for (i = 0; i < 2*n+1; i++) flag[i] = 0;
    for (i = 0; i < n; i++) {
        j = nodes[i].left;
        if (j < 0) {
            j = -j-1;
            if (j >= i) break;
        }
        else j += n;
        if (flag[j]) break;
        flag[j] = 1;
        j = nodes[i].right;
        if (j < 0) {
          j = -j-1;
          if (j >= i) break;
        }
        else j += n;
        if (flag[j]) break;
        flag[j] = 1;
    }
    PyMem_Free(flag);
    if (i < n) {
        /* break encountered */
        PyMem_Free(nodes);
        Py_DECREF(self);
        PyErr_SetString(PyExc_ValueError, "Inconsistent tree");
        return NULL;
    }
    self->n = n;
    self->nodes = nodes;
    return (PyObject*)self;
}

static PyObject*
PyTree_str(PyTree* self)
{
    int i;
    const int n = self->n;
    char string[128];
    Node node;
    PyObject* line;
    PyObject* output;
    PyObject* temp;

    output = PyUnicode_FromString("");
    for (i = 0; i < n; i++) {
        node = self->nodes[i];
        sprintf(string, "(%d, %d): %g", node.left, node.right, node.distance);
        if (i < n-1) strcat(string, "\n");
        line = PyUnicode_FromString(string);
        if (!line) {
            Py_DECREF(output);
            return NULL;
        }
        temp = PyUnicode_Concat(output, line);
        if (!temp) {
            Py_DECREF(output);
            Py_DECREF(line);
            return NULL;
        }
        output = temp;
    }
    return output;
}

static int
PyTree_length(PyTree *self)
{
    return self->n;
}

static PyObject*
PyTree_subscript(PyTree* self, PyObject* item)
{
    if (PyIndex_Check(item)) {
        PyNode* result;
        Py_ssize_t i;
        i = PyNumber_AsSsize_t(item, PyExc_IndexError);
        if (i == -1 && PyErr_Occurred())
            return NULL;
        if (i < 0)
            i += self->n;
        if (i < 0 || i >= self->n) {
            PyErr_SetString(PyExc_IndexError, "tree index out of range");
            return NULL;
        }
        result = (PyNode*) PyNodeType.tp_alloc(&PyNodeType, 0);
        if (!result) return PyErr_NoMemory();
        result->node = self->nodes[i];
        return (PyObject*) result;
    }
    else if (PySlice_Check(item)) {
        Py_ssize_t i, j;
        Py_ssize_t start, stop, step, slicelength;
        if (PySlice_GetIndicesEx(item, self->n, &start, &stop, &step,
                                 &slicelength) == -1) return NULL;
        if (slicelength == 0) return PyList_New(0);
        else {
            PyNode* node;
            PyObject* result = PyList_New(slicelength);
            if (!result) return PyErr_NoMemory();
            for (i = 0, j = start; i < slicelength; i++, j += step) {
                node = (PyNode*) PyNodeType.tp_alloc(&PyNodeType, 0);
                if (!node) {
                    Py_DECREF(result);
                    return PyErr_NoMemory();
                }
                node->node = self->nodes[j];
                PyList_SET_ITEM(result, i, (PyObject*)node);
            }
            return result;
        }
    }
    else {
        PyErr_Format(PyExc_TypeError,
                     "tree indices must be integers, not %.200s",
                     item->ob_type->tp_name);
        return NULL;
    }
}

static PyMappingMethods PyTree_mapping = {
    (lenfunc)PyTree_length,           /* mp_length */
    (binaryfunc)PyTree_subscript,     /* mp_subscript */
};

static char PyTree_scale__doc__[] =
"mytree.scale()\n"
"\n"
"Scale the node distances in the tree such that they are all between one\n"
"and zero.\n";

static PyObject*
PyTree_scale(PyTree* self)
{
    int i;
    const int n = self->n;
    Node* nodes = self->nodes;
    double maximum = DBL_MIN;

    for (i = 0; i < n; i++) {
        double distance = nodes[i].distance;
        if (distance > maximum) maximum = distance;
    }
    if (maximum != 0.0)
        for (i = 0; i < n; i++) nodes[i].distance /= maximum;
    Py_INCREF(Py_None);
    return Py_None;
}

static char PyTree_cut__doc__[] =
"mytree.cut(nclusters) -> array\n"
"\n"
"Divide the elements in a hierarchical clustering result mytree into\n"
"clusters, and return an array with the number of the cluster to which each\n"
"element was assigned. The number of clusters is given by nclusters.\n";

static PyObject*
PyTree_cut(PyTree* self, PyObject* args)
{
    int ok = -1;
    int nclusters;
    const int n = self->n + 1;
    Py_buffer indices = {0};

    if (!PyArg_ParseTuple(args, "O&i",
                          index_converter, &indices, &nclusters)) goto exit;
    if (nclusters < 1) {
        PyErr_SetString(PyExc_ValueError,
                        "requested number of clusters should be positive");
        goto exit;
    }
    if (nclusters > n) {
        PyErr_SetString(PyExc_ValueError,
                        "more clusters requested than items available");
        goto exit;
    }
    if (indices.shape[0] != n) {
        PyErr_SetString(PyExc_RuntimeError,
                        "indices array inconsistent with tree");
        goto exit;
    }
    ok = cuttree(n, self->nodes, nclusters, indices.buf);

exit:
    index_converter(NULL, &indices);
    if (ok == -1) return NULL;
    if (ok == 0) return PyErr_NoMemory();
    Py_INCREF(Py_None);
    return Py_None;
}

static char PyTree_sort__doc__[] =
"mytree.sort(order) -> array\n"
"\n"
"Sort a hierarchical clustering tree by switching the left and right\n"
"subnode of nodes such that the elements in the left-to-right order of the\n"
"tree tend to have increasing order values.\n"
"\n"
"Return the indices of the elements in the left-to-right order in the\n"
"hierarchical clustering tree, such that the element with index indices[i]\n"
"occurs at position i in the dendrogram.\n";

static PyObject*
PyTree_sort(PyTree* self, PyObject* args)
{
    int ok = -1;
    Py_buffer indices = {0};
    const int n = self->n;
    Py_buffer order = {0};

    if (n == 0) {
        PyErr_SetString(PyExc_ValueError, "tree is empty");
        return NULL;
    }
    if (!PyArg_ParseTuple(args, "O&O&",
                          index_converter, &indices,
                          vector_converter, &order)) goto exit;
    if (indices.shape[0] != n + 1) {
        PyErr_SetString(PyExc_RuntimeError,
                        "indices array inconsistent with tree");
        goto exit;
    }
    if (order.shape[0] != n + 1) {
        PyErr_Format(PyExc_ValueError,
            "order array has incorrect size %zd (expected %d)",
            order.shape[0], n + 1);
        goto exit;
    }
    ok = sorttree(n, self->nodes, order.buf, indices.buf);
exit:
    index_converter(NULL, &indices);
    vector_converter(NULL, &order);
    if (ok == -1) return NULL;
    if (ok == 0) return PyErr_NoMemory();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef PyTree_methods[] = {
    {"scale", (PyCFunction)PyTree_scale, METH_NOARGS, PyTree_scale__doc__},
    {"cut", (PyCFunction)PyTree_cut, METH_VARARGS, PyTree_cut__doc__},
    {"sort", (PyCFunction)PyTree_sort, METH_VARARGS, PyTree_sort__doc__},
    {NULL}  /* Sentinel */
};

static char PyTree_doc[] =
"Tree objects store a hierarchical clustering solution.\n"
"Individual nodes in the tree can be accessed with tree[i], where i is\n"
"an integer. Whereas the tree itself is a read-only object, tree[:]\n"
"returns a list of all the nodes, which can then be modified. To create\n"
"a new Tree from this list, use Tree(list).\n"
"See the description of the Node class for more information.";

static PyTypeObject PyTreeType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_cluster.Tree",             /* tp_name */
    sizeof(PyTree),              /* tp_basicsize */
    0,                           /* tp_itemsize */
    (destructor)PyTree_dealloc,  /* tp_dealloc */
    0,                           /* tp_print */
    0,                           /* tp_getattr */
    0,                           /* tp_setattr */
    0,                           /* tp_compare */
    0,                           /* tp_repr */
    0,                           /* tp_as_number */
    0,                           /* tp_as_sequence */
    &PyTree_mapping,             /* tp_as_mapping */
    0,                           /* tp_hash */
    0,                           /* tp_call */
    (reprfunc)PyTree_str,        /* tp_str */
    0,                           /* tp_getattro */
    0,                           /* tp_setattro */
    0,                           /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,          /*tp_flags*/
    PyTree_doc,                  /* tp_doc */
    0,                           /* tp_traverse */
    0,                           /* tp_clear */
    0,                           /* tp_richcompare */
    0,                           /* tp_weaklistoffset */
    0,                           /* tp_iter */
    0,                           /* tp_iternext */
    PyTree_methods,              /* tp_methods */
    NULL,                        /* tp_members */
    0,                           /* tp_getset */
    0,                           /* tp_base */
    0,                           /* tp_dict */
    0,                           /* tp_descr_get */
    0,                           /* tp_descr_set */
    0,                           /* tp_dictoffset */
    0,                           /* tp_init */
    0,                           /* tp_alloc */
    (newfunc)PyTree_new,         /* tp_new */
};

/* ========================================================================= */
/* -- Methods -------------------------------------------------------------- */
/* ========================================================================= */

/* version */
static char version__doc__[] =
"version() -> string\n"
"\n"
"Return the version number of the C Clustering Library as a string.\n";

static PyObject*
py_version(PyObject* self)
{
    return PyUnicode_FromString( CLUSTERVERSION );
}

/* kcluster */
static char kcluster__doc__[] =
"kcluster(data, nclusters, mask, weight, transpose, npass, method,\n"
"         dist, clusterid) -> None\n"
"\n"
"This function implements k-means clustering.\n"
"\n"
"Arguments:\n"
"\n"
" - data: nrows x ncols array containing the data to be clustered\n"
"\n"
" - nclusters: number of clusters (the 'k' in k-means)\n"
"\n"
" - mask: nrows x ncols array of integers, showing which data are\n"
"   missing. If mask[i,j] == 0, then data[i,j] is missing.\n"
"\n"
" - weight: the weights to be used when calculating distances\n"
" - transpose:\n"
"\n"
"   - if equal to 0, rows are clustered;\n"
"   - if equal to 1, columns are clustered.\n"
"\n"
" - npass: number of times the k-means clustering algorithm is\n"
"   performed, each time with a different (random) initial\n"
"   condition. If npass == 0, then the assignments in clusterid\n"
"   are used as the initial condition.\n"
"\n"
" - method: specifies how the center of a cluster is found:\n"
"\n"
"   - method == 'a': arithmetic mean\n"
"   - method == 'm': median\n"
"\n"
" - dist: specifies the distance function to be used:\n"
"\n"
"   - dist == 'e': Euclidean distance\n"
"   - dist == 'b': City Block distance\n"
"   - dist == 'c': Pearson correlation\n"
"   - dist == 'a': absolute value of the correlation\n"
"   - dist == 'u': uncentered correlation\n"
"   - dist == 'x': absolute uncentered correlation\n"
"   - dist == 's': Spearman's rank correlation\n"
"   - dist == 'k': Kendall's tau\n"
"\n"
" - clusterid: array in which the final clustering solution will be\n"
"   stored (output variable). If npass == 0, then clusterid is also used\n"
"   as an input variable, containing the initial condition from which\n"
"   the EM algorithm should start. In this case, the k-means algorithm\n"
"   is fully deterministic.\n"
"\n";

static PyObject*
py_kcluster(PyObject* self, PyObject* args, PyObject* keywords)
{
    int nclusters = 2;
    int nrows, ncols;
    int nitems;
    int ndata;
    Data data = {0};
    Mask mask = {0};
    Py_buffer weight = {0};
    int transpose = 0;
    int npass = 1;
    char method = 'a';
    char dist = 'e';
    Py_buffer clusterid = {0};
    double error;
    int ifound = 0;

    static char* kwlist[] = {"data",
                             "nclusters",
                             "mask",
                             "weight",
                             "transpose",
                             "npass",
                             "method",
                             "dist",
                             "clusterid",
                              NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywords, "O&iO&O&iiO&O&O&", kwlist,
                                     data_converter, &data,
                                     &nclusters,
                                     mask_converter, &mask,
                                     vector_converter, &weight,
                                     &transpose,
                                     &npass,
                                     method_kcluster_converter, &method,
                                     distance_converter, &dist,
                                     index_converter, &clusterid)) return NULL;
    if (!data.values) {
        PyErr_SetString(PyExc_RuntimeError, "data is None");
        goto exit;
    }
    if (!mask.values) {
        PyErr_SetString(PyExc_RuntimeError, "mask is None");
        goto exit;
    }
    if (data.nrows != mask.view.shape[0] ||
        data.ncols != mask.view.shape[1]) {
        PyErr_Format(PyExc_ValueError,
            "mask has incorrect dimensions %zd x %zd (expected %d x %d)",
            mask.view.shape[0], mask.view.shape[1], data.nrows, data.ncols);
        goto exit;
    }
    nrows = data.nrows;
    ncols = data.ncols;
    ndata = transpose ? nrows : ncols;
    nitems = transpose ? ncols : nrows;
    if (weight.shape[0] != ndata) {
        PyErr_Format(PyExc_ValueError,
                     "weight has incorrect size %zd (expected %d)",
                     weight.shape[0], ndata);
        goto exit;
    }
    if (nclusters < 1) {
        PyErr_SetString(PyExc_ValueError, "nclusters should be positive");
        goto exit;
    }
    if (nitems < nclusters) {
        PyErr_SetString(PyExc_ValueError,
                        "more clusters than items to be clustered");
        goto exit;
    }
    if (npass < 0) {
        PyErr_SetString(PyExc_RuntimeError, "expected a non-negative integer");
        goto exit;
    }
    else if (npass == 0) {
        int n = check_clusterid(clusterid, nitems);
        if (n == 0) goto exit;
        if (n != nclusters) {
            PyErr_SetString(PyExc_ValueError,
                            "more clusters requested than found in clusterid");
            goto exit;
        }
    }
    kcluster(nclusters,
             nrows,
             ncols,
             data.values,
             mask.values,
             weight.buf,
             transpose,
             npass,
             method,
             dist,
             clusterid.buf,
             &error,
             &ifound);
exit:
    data_converter(NULL, &data);
    mask_converter(NULL, &mask);
    vector_converter(NULL, &weight);
    index_converter(NULL, &clusterid);
    if (ifound) return Py_BuildValue("di", error, ifound);
    return NULL;
}
/* end of wrapper for kcluster */

/* kmedoids */
static char kmedoids__doc__[] =
"kmedoids(distance, nclusters, npass, clusterid) -> error, nfound\n"
"\n"
"This function implements k-medoids clustering.\n"
"\n"
"Arguments:\n"
" - distance: The distance matrix between the elements. There are three\n"
"   ways in which you can pass a distance matrix:\n"
"\n"
"   1. a 2D Numerical Python array (in which only the left-lower\n"
"      part of the array will be accessed);\n"
"   2. a 1D Numerical Python array containing the distances\n"
"      consecutively;\n"
"   3. a list of rows containing the lower-triangular part of\n"
"      the distance matrix.\n"
"\n"
"   Examples are:\n"
"\n"
"       >>> from numpy import array\n"
"       >>> distance = array([[0.0, 1.1, 2.3],\n"
"       ...                   [1.1, 0.0, 4.5],\n"
"       ...                   [2.3, 4.5, 0.0]])\n"
"       >>> # (option #1)\n"
"       >>> distance = array([1.1, 2.3, 4.5])\n"
"       >>> # (option #2)\n"
"       >>> distance = [array([]),\n"
"       ...             array([1.1]),\n"
"       ...             array([2.3, 4.5])]\n"
"       >>> # (option #3)\n"
"\n"
"   These three correspond to the same distance matrix.\n"
"\n"
" - nclusters: number of clusters (the 'k' in k-medoids)\n"
"\n"
" - npass: number of times the k-medoids clustering algorithm is\n"
"   performed, each time with a different (random) initial\n"
"   condition. If npass == 0, then the assignments in clusterid\n"
"   are used as the initial condition.\n"
"\n"
" - clusterid: array in which the final clustering solution will be\n"
"   stored (output variable). If npass == 0, then clusterid is also used\n"
"   as an input variable, containing the initial condition from which\n"
"   the EM algorithm should start. In this case, the k-medoids algorithm\n"
"   is fully deterministic.\n"
"\n"
"Return values:\n"
" - error: the within-cluster sum of distances for the returned k-means\n"
"   clustering solution;\n"
" - nfound: the number of times this solution was found.\n";

static PyObject*
py_kmedoids(PyObject* self, PyObject* args, PyObject* keywords)
{
    int nclusters = 2;
    Distancematrix distances = {0};
    Py_buffer clusterid = {0};
    int npass = 1;
    double error;
    int ifound = -2;

    static char* kwlist[] = {"distance",
                             "nclusters",
                             "npass",
                             "clusterid",
                              NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywords, "O&iiO&", kwlist,
                                     distancematrix_converter, &distances,
                                     &nclusters,
                                     &npass,
                                     index_converter, &clusterid)) return NULL;
    if (npass < 0) {
        PyErr_SetString(PyExc_RuntimeError, "expected a non-negative integer");
        goto exit;
    }
    else if (npass == 0) {
        int n = check_clusterid(clusterid, distances.n);
        if (n == 0) goto exit;
        if (n != nclusters) {
            PyErr_SetString(PyExc_RuntimeError,
                            "more clusters requested than found in clusterid");
            goto exit;
        }
    }
    if (nclusters <= 0) {
        PyErr_SetString(PyExc_ValueError,
                        "nclusters should be a positive integer");
        goto exit;
    }
    if (distances.n < nclusters) {
        PyErr_SetString(PyExc_ValueError,
                        "more clusters requested than items to be clustered");
        goto exit;
    }
    kmedoids(nclusters,
             distances.n,
             distances.values,
             npass,
             clusterid.buf,
             &error,
             &ifound);

exit:
    distancematrix_converter(NULL, &distances);
    index_converter(NULL, &clusterid);
    switch (ifound) {
        case -2:
            return NULL;
        case -1:
            return PyErr_NoMemory();
        case 0: /* should not occur */
            PyErr_SetString(PyExc_RuntimeError,
                        "error in kmedoids input arguments");
            return NULL;
        default:
            return Py_BuildValue("di", error, ifound);
    }
}
/* end of wrapper for kmedoids */

/* treecluster */
static char treecluster__doc__[] =
"treecluster(tree, data, mask, weight, transpose, dist, method,\n"
"            distancematrix) -> None\n"
"\n"
"This function implements the pairwise single, complete, centroid, and\n"
"average linkage hierarchical clustering methods.\n"
"\n"
"Arguments:\n"
" - tree: an empty Tree object; its nodes will be filled by treecluster\n"
"   to describe the hierarchical clustering result. See the description\n"
"   of the Tree class for more information.\n"
"\n"
" - data: nrows x ncols array containing the data to be clustered.\n"
"   Either data or distancematrix (see below) should be None.\n"
"\n"
" - mask: nrows x ncols array of integers, showing which data are\n"
"   missing. If mask[i,j]==0, then data[i,j] is missing.\n"
"\n"
" - weight: the weights to be used when calculating distances.\n"
"\n"
" - transpose:\n"
"\n"
"   - if equal to 0, rows are clustered;\n"
"   - if equal to 1, columns are clustered.\n"
"\n"
" - dist: specifies the distance function to be used:\n"
"\n"
"   - dist == 'e': Euclidean distance\n"
"   - dist == 'b': City Block distance\n"
"   - dist == 'c': Pearson correlation\n"
"   - dist == 'a': absolute value of the correlation\n"
"   - dist == 'u': uncentered correlation\n"
"   - dist == 'x': absolute uncentered correlation\n"
"   - dist == 's': Spearman's rank correlation\n"
"   - dist == 'k': Kendall's tau\n"
"\n"
" - method: specifies which linkage method is used:\n"
"\n"
"   - method == 's': Single pairwise linkage\n"
"   - method == 'm': Complete (maximum) pairwise linkage (default)\n"
"   - method == 'c': Centroid linkage\n"
"   - method == 'a': Average pairwise linkage\n"
"\n"
" - distancematrix:  The distance matrix between the elements.\n"
"   Either data (see above) or distancematrix should be None.\n"
"   There are three ways in which you can pass a distance matrix:\n"
"\n"
"   1. a 2D Numerical Python array (in which only the left-lower\n"
"      part of the array will be accessed);\n"
"   2. a 1D Numerical Python array containing the distances\n"
"      consecutively;\n"
"   3. a list of rows containing the lower-triangular part of\n"
"      the distance matrix.\n"
"\n"
"   Examples are:\n"
"\n"
"       >>> from numpy import array\n"
"       >>> distance = array([[0.0, 1.1, 2.3],\n"
"       ...                   [1.1, 0.0, 4.5],\n"
"       ...                   [2.3, 4.5, 0.0]])\n"
"       >>> # option 1.\n"
"       >>> distance = array([1.1, 2.3, 4.5])\n"
"       >>> # option 2.\n"
"       >>> distance = [array([]),\n"
"       ...             array([1.1]),\n"
"       ...             array([2.3, 4.5])]\n"
"       >>> # option 3.\n"
"\n"
"   These three correspond to the same distance matrix.\n"
"\n"
"   PLEASE NOTE:\n"
"   As the treecluster routine may shuffle the values in the\n"
"   distance matrix as part of the clustering algorithm, be sure\n"
"   to save this array in a different variable before calling\n"
"   treecluster if you need it later.\n"
"\n"
"Either data or distancematrix should be None. If distancematrix is None,\n"
"the hierarchical clustering solution is calculated from the values in\n"
"the argument data. Instead if data is None, the hierarchical clustering\n"
"solution is calculated from the distance matrix.\n"
"Pairwise centroid-linkage clustering can be calculated only from the data\n"
"and not from the distance matrix.\n"
"Pairwise single-, maximum-, and average-linkage clustering can be\n"
"calculated from either the data or from the distance matrix.\n";

static PyObject*
py_treecluster(PyObject* self, PyObject* args, PyObject* keywords)
{
    Data data = {0};
    Mask mask = {0};
    Py_buffer weight = {0};
    int transpose = 0;
    char dist = 'e';
    char method = 'm';
    Distancematrix distances = {0};
    PyTree* tree = NULL;
    Node* nodes;
    int nitems;

    static char* kwlist[] = {"tree",
                             "data",
                             "mask",
                             "weight",
                             "transpose",
                             "method",
                             "dist",
                             "distancematrix",
                              NULL };

    if (!PyArg_ParseTupleAndKeywords(args, keywords, "O!O&O&O&iO&O&O&", kwlist,
                                     &PyTreeType, &tree,
                                     data_converter, &data,
                                     mask_converter, &mask,
                                     vector_none_converter, &weight,
                                     &transpose,
                                     method_treecluster_converter, &method,
                                     distance_converter, &dist,
                                     distancematrix_converter, &distances))
        return NULL;

    if (tree->n != 0) {
        PyErr_SetString(PyExc_RuntimeError, "expected an empty tree");
        goto exit;
    }
    if (data.values != NULL && distances.values != NULL) {
        PyErr_SetString(PyExc_ValueError,
            "use either data or distancematrix, do not use both");
        goto exit;
    }
    if (data.values == NULL && distances.values == NULL) {
        PyErr_SetString(PyExc_ValueError,
                        "neither data nor distancematrix was given");
        goto exit;
    }

    if (data.values) /* use the values in data, not the distance matrix */ {
        int nrows;
        int ncols;
        int ndata;

        if (!mask.values) {
            PyErr_SetString(PyExc_RuntimeError, "mask is None");
            goto exit;
        }
        if (!weight.buf) {
            PyErr_SetString(PyExc_RuntimeError, "weight is None");
            goto exit;
        }
        nrows = data.nrows;
        ncols = data.ncols;
        if (nrows != mask.view.shape[0] || ncols != mask.view.shape[1]) {
            PyErr_Format(PyExc_ValueError,
                "mask has incorrect dimensions (%zd x %zd, expected %d x %d)",
                mask.view.shape[0], mask.view.shape[1],
                data.nrows, data.ncols);
            goto exit;
        }
        ndata = transpose ? nrows : ncols;
        nitems = transpose ? ncols : nrows;
        if (weight.shape[0] != ndata) {
            PyErr_Format(PyExc_RuntimeError,
                         "weight has incorrect size %zd (expected %d)",
                         weight.shape[0], ndata);
            goto exit;
        }

        nodes = treecluster(nrows,
                            ncols,
                            data.values,
                            mask.values,
                            weight.buf,
                            transpose,
                            dist,
                            method,
                            NULL);
    }
    else { /* use the distance matrix instead of the values in data */
        if (!strchr("sma", method)) {
            PyErr_SetString(PyExc_ValueError,
                            "argument method should be 's', 'm', or 'a' "
                            "when specifying the distance matrix");
            goto exit;
        }
        nitems = distances.n;
        nodes = treecluster(nitems,
                            nitems,
                            0,
                            0,
                            0,
                            transpose,
                            dist,
                            method,
                            distances.values);
    }

    if (!nodes) {
        PyErr_NoMemory();
        goto exit;
    }
    tree->nodes = nodes;
    tree->n = nitems-1;

exit:
    data_converter(NULL, &data);
    mask_converter(NULL, &mask);
    vector_none_converter(NULL, &weight);
    distancematrix_converter(NULL, &distances);
    if (tree == NULL || tree->n == 0) return NULL;
    Py_INCREF(Py_None);
    return Py_None;
}
/* end of wrapper for treecluster */

/* somcluster */
static char somcluster__doc__[] =
"somcluster(clusterid, celldata, data, mask, weight, transpose,\n"
"           inittau, niter, dist) -> None\n"
"\n"
"This function implements a self-organizing map on a rectangular grid.\n"
"\n"
"Arguments:\n"
" - clusterid: array with two columns, with the number of rows equal\n"
"   to the number of items being clustered. Upon return, each row\n"
"   in the array contains the x and y coordinates of the cell in the\n"
"   the rectangular SOM grid to which the item was assigned.\n"
"\n"
" - celldata: array with dimensions nxgrid x nygrid x number of columns\n"
"   if rows are being clustered, or nxgrid x nygrid x number of rows\n"
"   if columns are being clustered, where nxgrid is the horizontal\n"
"   dimension of the rectangular SOM map and nygrid is the vertical\n"
"   dimension of the rectangular SOM map.\n"
"   Upon return, each element [ix, iy] of this array contains the\n"
"   data for the centroid of the cluster in the SOM grid cell with\n"
"   coordinates [ix, iy].\n"
"\n"
" - data: nrows x ncols array containing the data to be clustered.\n"
"\n"
" - mask: nrows x ncols array of integers, showing which data are\n"
"   missing. If mask[i,j] == 0, then data[i,j] is missing.\n"
"\n"
" - weight: the weights to be used when calculating distances\n"
"\n"
" - transpose:\n"
"\n"
"   - if equal to 0, rows are clustered;\n"
"   - if equal to 1, columns are clustered.\n"
"\n"
" - inittau: the initial value of tau (the neighborbood function)\n"
"\n"
" - niter: the number of iterations\n"
"\n"
" - dist: specifies the distance function to be used:\n"
"\n"
"   - dist == 'e': Euclidean distance\n"
"   - dist == 'b': City Block distance\n"
"   - dist == 'c': Pearson correlation\n"
"   - dist == 'a': absolute value of the correlation\n"
"   - dist == 'u': uncentered correlation\n"
"   - dist == 'x': absolute uncentered correlation\n"
"   - dist == 's': Spearman's rank correlation\n"
"   - dist == 'k': Kendall's tau\n";

static PyObject*
py_somcluster(PyObject* self, PyObject* args, PyObject* keywords)
{
    int nrows;
    int ncols;
    int ndata;
    Data data = {0};
    Mask mask = {0};
    Py_buffer weight = {0};
    int transpose = 0;
    double inittau = 0.02;
    int niter = 1;
    char dist = 'e';
    Py_buffer indices = {0};
    Celldata celldata = {0};
    PyObject* result = NULL;

    static char* kwlist[] = {"clusterids",
                             "celldata",
                             "data",
                             "mask",
                             "weight",
                             "transpose",
                             "inittau",
                             "niter",
                             "dist",
                             NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywords, "O&O&O&O&O&idiO&", kwlist,
                                     index2d_converter, &indices,
                                     celldata_converter, &celldata,
                                     data_converter, &data,
                                     mask_converter, &mask,
                                     vector_converter, &weight,
                                     &transpose,
                                     &inittau,
                                     &niter,
                                     distance_converter, &dist)) return NULL;
    if (niter < 1) {
        PyErr_SetString(PyExc_ValueError,
                      "number of iterations (niter) should be positive");
        goto exit;
    }
    if (!data.values) {
        PyErr_SetString(PyExc_RuntimeError, "data is None");
        goto exit;
    }
    if (!mask.values) {
        PyErr_SetString(PyExc_RuntimeError, "mask is None");
        goto exit;
    }
    nrows = data.nrows;
    ncols = data.ncols;
    if (nrows != mask.view.shape[0] || ncols != mask.view.shape[1]) {
        PyErr_Format(PyExc_ValueError,
            "mask has incorrect dimensions (%zd x %zd, expected %d x %d)",
            mask.view.shape[0], mask.view.shape[1], data.nrows, data.ncols);
        goto exit;
    }
    ndata = transpose ? nrows : ncols;
    if (weight.shape[0] != ndata) {
        PyErr_Format(PyExc_RuntimeError,
                     "weight has incorrect size %zd (expected %d)",
                     weight.shape[0], ndata);
        goto exit;
    }
    if (celldata.nz != ndata) {
        PyErr_Format(PyExc_RuntimeError,
                    "the celldata array size is not consistent with the data "
                    "(last dimension is %d; expected %d)", celldata.nz, ndata);
        goto exit;
    }
    somcluster(nrows,
               ncols,
               data.values,
               mask.values,
               weight.buf,
               transpose,
               celldata.nx,
               celldata.ny,
               inittau,
               niter,
               dist,
               celldata.values,
               indices.buf);
    Py_INCREF(Py_None);
    result = Py_None;

exit:
    data_converter(NULL, &data);
    vector_converter(NULL, &weight);
    index2d_converter(NULL, &indices);
    celldata_converter(NULL, &celldata);
    return result;
}
/* end of wrapper for somcluster */

/* clusterdistance */
static char clusterdistance__doc__[] =
"clusterdistance(data, mask, weight, index1, index2, dist, method,\n"
"                transpose) -> distance between two clusters\n"
"\n"
"Arguments:\n"
"\n"
" - data: nrows x ncols array containing the data values.\n"
"\n"
" - mask: nrows x ncols array of integers, showing which data are\n"
"   missing. If mask[i,j] == 0, then data[i,j] is missing.\n"
"\n"
" - weight: the weights to be used when calculating distances\n"
"\n"
" - index1: 1D array identifying which items belong to the first\n"
"   cluster.\n"
"\n"
" - index2: 1D array identifying which items belong to the second\n"
"   cluster.\n"
"\n"
" - dist: specifies the distance function to be used:\n"
"\n"
"   - dist == 'e': Euclidean distance\n"
"   - dist == 'b': City Block distance\n"
"   - dist == 'c': Pearson correlation\n"
"   - dist == 'a': absolute value of the correlation\n"
"   - dist == 'u': uncentered correlation\n"
"   - dist == 'x': absolute uncentered correlation\n"
"   - dist == 's': Spearman's rank correlation\n"
"   - dist == 'k': Kendall's tau\n"
"\n"
" - method: specifies how the distance between two clusters is defined:\n"
"\n"
"   - method == 'a': the distance between the arithmetic means of the\n"
"     two clusters\n"
"   - method == 'm': the distance between the medians of the two\n"
"     clusters\n"
"   - method == 's': the smallest pairwise distance between members\n"
"     of the two clusters\n"
"   - method == 'x': the largest pairwise distance between members of\n"
"     the two clusters\n"
"   - method == 'v': average of the pairwise distances between\n"
"     members of the clusters\n"
"\n"
" - transpose:\n"
"\n"
"   - if equal to 0: clusters of rows are considered;\n"
"   - if equal to 1: clusters of columns are considered.\n"
"\n";

static PyObject*
py_clusterdistance(PyObject* self, PyObject* args, PyObject* keywords)
{
    double distance;
    int nrows;
    int ncols;
    int ndata;
    Data data = {0};
    Mask mask = {0};
    Py_buffer weight = {0};
    char dist = 'e';
    char method = 'a';
    int transpose = 0;
    Py_buffer index1 = {0};
    Py_buffer index2 = {0};
    PyObject* result = NULL;

    static char* kwlist[] = {"data",
                             "mask",
                             "weight",
                             "index1",
                             "index2",
                             "method",
                             "dist",
                             "transpose",
                              NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywords, "O&O&O&O&O&O&O&i", kwlist,
                                     data_converter, &data,
                                     mask_converter, &mask,
                                     vector_converter, &weight,
                                     index_converter, &index1,
                                     index_converter, &index2,
                                     method_clusterdistance_converter, &method,
                                     distance_converter, &dist,
                                     &transpose)) return NULL;
    if (!data.values) {
        PyErr_SetString(PyExc_RuntimeError, "data is None");
        goto exit;
    }
    if (!mask.values) {
        PyErr_SetString(PyExc_RuntimeError, "mask is None");
        goto exit;
    }
    nrows = data.nrows;
    ncols = data.ncols;
    ndata = transpose ? nrows : ncols;
    if (nrows != mask.view.shape[0] || ncols != mask.view.shape[1]) {
        PyErr_Format(PyExc_ValueError,
            "mask has incorrect dimensions (%zd x %zd, expected %d x %d)",
            mask.view.shape[0], mask.view.shape[1], data.nrows, data.ncols);
        goto exit;
    }
    if (weight.shape[0] != ndata) {
        PyErr_Format(PyExc_RuntimeError,
                     "weight has incorrect size %zd (expected %d)",
                     weight.shape[0], ndata);
        goto exit;
    }

    distance = clusterdistance(nrows,
                               ncols,
                               data.values,
                               mask.values,
                               weight.buf,
                               (int) index1.shape[0],
                               (int) index2.shape[0],
                               index1.buf,
                               index2.buf,
                               dist,
                               method,
                               transpose);

    if (distance < -0.5) /* Actually -1.0; avoiding roundoff errors */
        PyErr_SetString(PyExc_IndexError, "index out of range");
    else
        result = PyFloat_FromDouble(distance);
exit:
    data_converter(NULL, &data);
    mask_converter(NULL, &mask);
    vector_converter(NULL, &weight);
    index_converter(NULL, &index1);
    index_converter(NULL, &index2);
    return result;
}
/* end of wrapper for clusterdistance */

/* clustercentroids */
static char clustercentroids__doc__[] =
"clustercentroids(data, mask, clusterid, method, transpose) -> cdata, cmask\n"
"\n"
"The clustercentroids routine calculates the cluster centroids, given to\n"
"which cluster each element belongs. The centroid is defined as either\n"
"the mean or the median over all elements for each dimension.\n"
"\n"
"Arguments:\n"
" - data: nrows x ncols array containing the data values.\n"
"\n"
" - mask: nrows x ncols array of integers, showing which data are\n"
"   missing. If mask[i,j] == 0, then data[i,j] is missing.\n"
"\n"
" - clusterid: array containing the cluster number for each item.\n"
"   The cluster number should be non-negative.\n"
"\n"
" - method: specifies whether the centroid is calculated from the\n"
"   arithmetic mean (method == 'a', default) or the median\n"
"   (method == 'm') over each dimension.\n"
"\n"
" - transpose: if equal to 0, row clusters are considered;\n"
"   if equal to 1, column clusters are considered.\n"
"\n"
" - cdata: 2D array containing, upon return, the cluster centroids.\n"
"   If transpose == 0, then the dimensions of cdata should be\n"
"   nclusters x ncols.\n"
"   If transpose == 1, then the dimensions of cdata should be \n"
"   nrows x nclusters.\n"
"\n"
" - cmask: 2D array of integers describing, upon return,  which elements\n"
"   in cdata, if any, are missing.\n";

static PyObject*
py_clustercentroids(PyObject* self, PyObject* args, PyObject* keywords)
{
    int nrows;
    int ncols;
    int nclusters;
    Data data = {0};
    Mask mask = {0};
    Data cdata = {0};
    Mask cmask = {0};
    Py_buffer clusterid = {0};
    char method = 'a';
    int transpose = 0;
    int ok = -1;

    static char* kwlist[] = {"data",
                             "mask",
                             "clusterid",
                             "method",
                             "transpose",
                             "cdata",
                             "cmask",
                              NULL };

    if (!PyArg_ParseTupleAndKeywords(args, keywords, "O&O&O&O&iO&O&", kwlist,
                                     data_converter, &data,
                                     mask_converter, &mask,
                                     index_converter, &clusterid,
                                     method_kcluster_converter, &method,
                                     &transpose,
                                     data_converter, &cdata,
                                     mask_converter, &cmask)) return NULL;
    if (!data.values) {
        PyErr_SetString(PyExc_RuntimeError, "data is None");
        goto exit;
    }
    if (!mask.values) {
        PyErr_SetString(PyExc_RuntimeError, "mask is None");
        goto exit;
    }
    nrows = data.nrows;
    ncols = data.ncols;
    if (nrows != mask.view.shape[0] || ncols != mask.view.shape[1]) {
        PyErr_Format(PyExc_ValueError,
            "mask has incorrect dimensions (%zd x %zd, expected %d x %d)",
            mask.view.shape[0], mask.view.shape[1], data.nrows, data.ncols);
        goto exit;
    }
    if (transpose == 0) {
        nclusters = check_clusterid(clusterid, nrows);
        nrows = nclusters;
    }
    else {
        nclusters = check_clusterid(clusterid, ncols);
        ncols = nclusters;
    }
    if (nclusters == 0) goto exit;
    if (cdata.nrows != nrows) {
        PyErr_Format(PyExc_RuntimeError,
                     "cdata has incorrect number of rows (%d, expected %d)",
                     cdata.nrows, nrows);
        goto exit;
    }
    if (cdata.ncols != ncols) {
        PyErr_Format(PyExc_RuntimeError,
                     "cdata has incorrect number of columns (%d, expected %d)",
                     cdata.ncols, ncols);
        goto exit;
    }
    if (cmask.view.shape[0] != nrows) {
        PyErr_Format(PyExc_RuntimeError,
                     "cmask has incorrect number of rows (%zd, expected %d)",
                     cmask.view.shape[0], nrows);
        goto exit;
    }
    if (cmask.view.shape[1] != ncols) {
        PyErr_Format(PyExc_RuntimeError,
                     "cmask has incorrect number of columns "
                     "(%zd, expected %d)", cmask.view.shape[1], ncols);
        goto exit;
    }
    ok = getclustercentroids(nclusters,
                             data.nrows,
                             data.ncols,
                             data.values,
                             mask.values,
                             clusterid.buf,
                             cdata.values,
                             cmask.values,
                             transpose,
                             method);
exit:
    data_converter(NULL, &data);
    mask_converter(NULL, &mask);
    data_converter(NULL, &cdata);
    mask_converter(NULL, &cmask);
    index_converter(NULL, &clusterid);
    if (ok == -1) return NULL;
    if (ok == 0) return PyErr_NoMemory();
    Py_INCREF(Py_None);
    return Py_None;
}
/* end of wrapper for clustercentroids */

/* distancematrix */
static char distancematrix__doc__[] =
"distancematrix(data, mask, weight, transpose, dist, distancematrix)\n"
"              -> None\n"
"\n"
"This function calculuates the distance matrix between the data values.\n"
"\n"
"Arguments:\n"
"\n"
" - data: nrows x ncols array containing the data values.\n"
"\n"
" - mask: nrows x ncols array of integers, showing which data are\n"
"   missing. If mask[i,j] == 0, then data[i,j] is missing.\n"
"\n"
" - weight: the weights to be used when calculating distances.\n"
"\n"
" - transpose: if equal to 0: the distances between rows are\n"
"   calculated;\n"
"   if equal to 1, the distances beteeen columns are calculated.\n"
"\n"
" - dist: specifies the distance function to be used:\n"
"\n"
"   - dist == 'e': Euclidean distance\n"
"   - dist == 'b': City Block distance\n"
"   - dist == 'c': Pearson correlation\n"
"   - dist == 'a': absolute value of the correlation\n"
"   - dist == 'u': uncentered correlation\n"
"   - dist == 'x': absolute uncentered correlation\n"
"   - dist == 's': Spearman's rank correlation\n"
"   - dist == 'k': Kendall's tau\n"
"\n"
" - distancematrix: Upon return, the distance matrix as a list of 1D\n"
"   arrays. The number of columns in each row is equal to the row number\n"
"   (i.e., len(distancematrix[i]) == i).\n"
"   An example of the return value is:\n"
"\n"
"    matrix = [[],\n"
"              array([1.]),\n"
"              array([7., 3.]),\n"
"              array([4., 2., 6.])]\n"
"\n"
"This corresponds to the distance matrix:\n"
"\n"
"    [0.\t1.\t7.\t4.]\n"
"    [1.\t0.\t3.\t2.]\n"
"    [7.\t3.\t0.\t6.]\n"
"    [4.\t2.\t6.\t0.]\n";

static PyObject*
py_distancematrix(PyObject* self, PyObject* args, PyObject* keywords)
{
    PyObject* list;
    Distancematrix distances = {0};
    Data data = {0};
    Mask mask = {0};
    Py_buffer weight = {0};
    int transpose = 0;
    char dist = 'e';
    int nrows, ncols, ndata;
    PyObject* result = NULL;

    /* -- Read the input variables --------------------------------------- */
    static char* kwlist[] = {"data",
                             "mask",
                             "weight",
                             "transpose",
                             "dist",
                             "distancematrix",
                              NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywords, "O&O&O&iO&O!", kwlist,
                                     data_converter, &data,
                                     mask_converter, &mask,
                                     vector_converter, &weight,
                                     &transpose,
                                     distance_converter, &dist,
                                     &PyList_Type, &list)) return NULL;
    if (!data.values) {
        PyErr_SetString(PyExc_RuntimeError, "data is None");
        goto exit;
    }
    if (!mask.values) {
        PyErr_SetString(PyExc_RuntimeError, "mask is None");
        goto exit;
    }
    nrows = data.nrows;
    ncols = data.ncols;
    if (nrows != mask.view.shape[0] || ncols != mask.view.shape[1]) {
        PyErr_Format(PyExc_ValueError,
            "mask has incorrect dimensions (%zd x %zd, expected %d x %d)",
            mask.view.shape[0], mask.view.shape[1], data.nrows, data.ncols);
        goto exit;
    }
    ndata = (transpose == 0) ? ncols : nrows;
    if (weight.shape[0] != ndata) {
        PyErr_Format(PyExc_ValueError,
                     "weight has incorrect size %zd (expected %d)",
                     weight.shape[0], ndata);
        goto exit;
    }
    if (_convert_list_to_distancematrix(list, &distances) == 0) goto exit;

    distancematrix(nrows,
                   ncols,
                   data.values,
                   mask.values,
                   weight.buf,
                   dist,
                   transpose,
                   distances.values);

    Py_INCREF(Py_None);
    result = Py_None;
exit:
    data_converter(NULL, &data);
    mask_converter(NULL, &mask);
    vector_converter(NULL, &weight);
    distancematrix_converter(NULL, &distances);
    return result;
}
/* end of wrapper for distancematrix */

/* pca */
static char pca__doc__[] =
"pca(data, columnmean, coordinates, pc, eigenvalues) -> None\n"
"\n"
"This function calculates the principal component decomposition\n"
"of the values in data.\n"
"\n"
"Arguments:\n"
"\n"
" - data: nrows x ncols array containing the data values.\n"
"\n"
" - columnmean: array of size nrows) in which the mean of each column\n"
"               will be sorted.\n"
"\n"
" - coordinates: nrows x nmin array in which the coordinates of the\n"
"                data along the principal components will be stored;\n"
"                nmin is min(nrows, ncols).\n"
"\n"
" - pc : the principal components as an nmin x ncols array, where nmin\n"
"        is min(nrows, ncols).\n"
"\n"
" - eigenvalues: array of size min(nrows, ncols), in which the\n"
"                eigenvalues will be stored, sorted by the magnitude\n"
"                of the eigenvalues, with the largest eigenvalues\n"
"                appearing first.\n"
"\n"
"Adding the column means to the dot product of the coordinates and the\n"
"principal components, i.e.\n"
"\n"
"   columnmean + dot(coordinates, pc)\n"
"\n"
"recreates the data matrix.\n";

static PyObject*
py_pca(PyObject* self, PyObject* args)
{
    Py_buffer eigenvalues = {0};
    double** u;
    double** v;
    Data data = {0};
    Data pc = {0};
    Data coordinates = {0};
    Py_buffer mean = {0};
    int nrows, ncols;
    int nmin;
    int error = -2;
    double* p;
    double** values;
    int i, j;

    if (!PyArg_ParseTuple(args, "O&O&O&O&O&",
                          data_converter, &data,
                          vector_converter, &mean,
                          data_converter, &coordinates,
                          data_converter, &pc,
                          vector_converter, &eigenvalues)) return NULL;

    values = data.values;
    if (!values) {
        PyErr_SetString(PyExc_RuntimeError, "data is None");
        goto exit;
    }
    nrows = data.nrows;
    ncols = data.ncols;
    if (mean.shape[0] != ncols) {
        PyErr_Format(PyExc_RuntimeError,
                     "columnmean has inconsistent size %zd (expected %d)",
                     mean.shape[0], ncols);
        goto exit;
    }
    nmin = nrows < ncols ? nrows : ncols;
    if (pc.nrows != nmin || pc.ncols != ncols) {
        PyErr_Format(PyExc_RuntimeError,
                     "pc has inconsistent size %zd x %zd (expected %d x %d)",
                     mean.shape[0], mean.shape[1], nmin, ncols);
        goto exit;
    }
    if (coordinates.nrows != nrows || coordinates.ncols != nmin) {
        PyErr_Format(PyExc_RuntimeError,
            "coordinates has inconsistent size %zd x %zd (expected %d x %d)",
            mean.shape[0], mean.shape[1], nrows, nmin);
        goto exit;
    }
    if (nrows >= ncols) {
        u = coordinates.values;
        v = pc.values;
    }
    else { /* nrows < ncolums */
        u = pc.values;
        v = coordinates.values;
    }
    /* -- Calculate the mean of each column ------------------------------ */
    p = mean.buf;
    for (j = 0; j < ncols; j++) {
        p[j] = 0.0;
        for (i = 0; i < nrows; i++) p[j] += values[i][j];
        p[j] /= nrows;
    }
    /* --   Subtract the mean of each column ----------------------------- */
    for (i = 0; i < nrows; i++)
        for (j = 0; j < ncols; j++)
            u[i][j] = values[i][j] - p[j];
    /* -- Perform the principal component analysis ----------------------- */
    error = pca(nrows, ncols, u, v, eigenvalues.buf);
    /* ------------------------------------------------------------------- */
exit:
    data_converter(NULL, &data);
    vector_converter(NULL, &mean);
    data_converter(NULL, &pc);
    data_converter(NULL, &coordinates);
    vector_converter(NULL, &eigenvalues);
    if (error == 0) {
        Py_INCREF(Py_None);
        return Py_None;
    }
    if (error == -1) return PyErr_NoMemory();
    else if (error > 0)
        PyErr_SetString(PyExc_RuntimeError,
            "Singular value decomposition failed to converge");
    return NULL;
}
/* end of wrapper for pca */

/* ========================================================================= */
/* -- The methods table ---------------------------------------------------- */
/* ========================================================================= */


static struct PyMethodDef cluster_methods[] = {
    {"version", (PyCFunction) py_version, METH_NOARGS, version__doc__},
    {"kcluster",
     (PyCFunction) py_kcluster,
     METH_VARARGS | METH_KEYWORDS,
     kcluster__doc__
    },
    {"kmedoids",
     (PyCFunction) py_kmedoids,
     METH_VARARGS | METH_KEYWORDS,
     kmedoids__doc__
    },
    {"treecluster",
     (PyCFunction) py_treecluster,
     METH_VARARGS | METH_KEYWORDS,
     treecluster__doc__
    },
    {"somcluster",
     (PyCFunction) py_somcluster,
     METH_VARARGS | METH_KEYWORDS,
     somcluster__doc__
    },
    {"clusterdistance",
     (PyCFunction) py_clusterdistance,
     METH_VARARGS | METH_KEYWORDS,
     clusterdistance__doc__
    },
    {"clustercentroids",
     (PyCFunction) py_clustercentroids,
     METH_VARARGS | METH_KEYWORDS,
     clustercentroids__doc__
    },
    {"distancematrix",
     (PyCFunction) py_distancematrix,
     METH_VARARGS | METH_KEYWORDS,
     distancematrix__doc__
    },
    {"pca",
     (PyCFunction) py_pca,
     METH_VARARGS | METH_KEYWORDS,
     pca__doc__
    },
    {NULL, NULL, 0, NULL} /* sentinel */
};

/* ========================================================================= */
/* -- Initialization ------------------------------------------------------- */
/* ========================================================================= */

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_cluster",
    "C Clustering Library",
    -1,
    cluster_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *
PyInit__cluster(void)
{
    PyObject *module;

    PyNodeType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&PyNodeType) < 0)
        return NULL;
    if (PyType_Ready(&PyTreeType) < 0)
        return NULL;

    module = PyModule_Create(&moduledef);
    if (module == NULL) return NULL;

    Py_INCREF(&PyTreeType);
    if (PyModule_AddObject(module, "Tree", (PyObject*) &PyTreeType) < 0) {
        Py_DECREF(module);
        Py_DECREF(&PyTreeType);
        return NULL;
    }

    Py_INCREF(&PyNodeType);
    if (PyModule_AddObject(module, "Node", (PyObject*) &PyNodeType) < 0) {
        Py_DECREF(module);
        Py_DECREF(&PyNodeType);
        return NULL;
    }

    return module;
}
