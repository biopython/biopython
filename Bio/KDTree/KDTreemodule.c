#include <Python.h>
#include "KDTree.h"


typedef struct {
    PyObject_HEAD
    struct Neighbor neighbor;
} PyNeighbor;

static int
PyNeighbor_init(PyNeighbor *self, PyObject *args, PyObject *kwds)
{
    long int index1, index2;
    float radius = 0.0;
    static char *kwlist[] = {"index1", "index2", "radius", NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "ii|d", kwlist,
                                      &index1, &index2, &radius))
        return -1;
    self->neighbor.index1 = index1;
    self->neighbor.index2 = index2;
    self->neighbor.radius = radius;

    return 0;
}

static PyObject*
PyNeighbor_repr(PyNeighbor* self)
{
    char string[64];
    sprintf(string, "(%ld, %ld): %g",
            self->neighbor.index1, self->neighbor.index2, self->neighbor.radius);
#if PY_MAJOR_VERSION >= 3
    return PyUnicode_FromString(string);
#else
    return PyString_FromString(string);
#endif

}

static char PyNeighbor_index1__doc__[] =
"index of the first neighbor";

static PyObject*
PyNeighbor_getindex1(PyNeighbor* self, void* closure)
{
#if PY_MAJOR_VERSION >= 3
    return PyLong_FromLong(self->neighbor.index1);
#else
    return PyInt_FromLong(self->neighbor.index1);
#endif
}

static int
PyNeighbor_setindex1(PyNeighbor* self, PyObject* value, void* closure)
{
#if PY_MAJOR_VERSION >= 3
    long index1 = PyLong_AsLong(value);
#else
    long index1 = PyInt_AsLong(value);
#endif
    if (PyErr_Occurred()) return -1;
    self->neighbor.index1 = index1;
    return 0;
}

static char PyNeighbor_index2__doc__[] =
"index of the second neighbor";

static PyObject*
PyNeighbor_getindex2(PyNeighbor* self, void* closure)
{
#if PY_MAJOR_VERSION >= 3
    return PyLong_FromLong(self->neighbor.index2);
#else
    return PyInt_FromLong(self->neighbor.index2);
#endif
}

static int
PyNeighbor_setindex2(PyNeighbor* self, PyObject* value, void* closure)
{
#if PY_MAJOR_VERSION >= 3
    long index2 = PyLong_AsLong(value);
#else
    long index2 = PyInt_AsLong(value);
#endif
    if (PyErr_Occurred()) return -1;
    self->neighbor.index2 = index2;
    return 0;
}

static PyObject*
PyNeighbor_getradius(PyNeighbor* self, void* closure)
{
    float value = self->neighbor.radius;
    return PyFloat_FromDouble((double)value);
}

static int
PyNeighbor_setradius(PyNeighbor* self, PyObject* value, void* closure)
{
    const double radius = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->neighbor.radius = (float)radius;
    return 0;
}

static char PyNeighbor_radius__doc__[] =
"the radius\n";

static PyGetSetDef PyNeighbor_getset[] = {
    {"index1", (getter)PyNeighbor_getindex1, (setter)PyNeighbor_setindex1, PyNeighbor_index1__doc__, NULL},
    {"index2", (getter)PyNeighbor_getindex2, (setter)PyNeighbor_setindex2, PyNeighbor_index2__doc__, NULL},
    {"radius", (getter)PyNeighbor_getradius, (setter)PyNeighbor_setradius, PyNeighbor_radius__doc__, NULL},
    {NULL}  /* Sentinel */
};

static char PyNeighbor_doc[] =
"A neighbor pair; members are index1, index2, and radius.\n";

static PyTypeObject PyNeighborType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "KDTree.Neighbor",         /* tp_name*/
    sizeof(PyNeighbor),        /* tp_basicsize*/
    0,                         /* tp_itemsize*/
    0,                         /* tp_dealloc*/
    0,                         /* tp_print*/
    0,                         /* tp_getattr*/
    0,                         /* tp_setattr*/
    0,                         /* tp_compare*/
    (reprfunc)PyNeighbor_repr, /* tp_repr*/
    0,                         /* tp_as_number*/
    0,                         /* tp_as_sequence*/
    0,                         /* tp_as_mapping*/
    0,                         /* tp_hash */
    0,                         /* tp_call*/
    0,                         /* tp_str*/
    0,                         /* tp_getattro*/
    0,                         /* tp_setattro*/
    0,                         /* tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /* tp_flags*/
    PyNeighbor_doc,            /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    0,                         /* tp_methods */
    0,                         /* tp_members */
    PyNeighbor_getset,         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)PyNeighbor_init, /* tp_init */
};

typedef struct {
    PyObject_HEAD
    struct KDTree* tree;
} PyTree;

static void
PyTree_dealloc(PyTree* self)
{
    KDTree_destroy(self->tree);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
PyTree_init(PyTree* self, PyObject* args, PyObject* kwds)
{
    int dim;
    int bucket_size;
    struct KDTree* tree;

    if(!PyArg_ParseTuple(args, "ii:KDTree_init" ,&dim, &bucket_size)) return -1;

    if (dim <= 0 || bucket_size <= 0) {
        PyErr_SetString(PyExc_ValueError, "Both arguments should be positive");
        return -1;
    }

    tree = KDTree_init(dim, bucket_size);
    if (tree==NULL) {
        PyErr_SetString(PyExc_MemoryError, "Insufficient memory for tree");
        return -1;
    }

    self->tree = tree;
    return 0;
}

static PyObject*
PyTree_get_count(PyTree* self)
{
    long count;
    struct KDTree* tree = self->tree;
    PyObject* result;
    count = KDTree_get_count(tree);
#if PY_MAJOR_VERSION >= 3
    result = PyLong_FromLong(count);
#else
    result = PyInt_FromLong(count);
#endif
    if (!result) {
        PyErr_SetString (PyExc_MemoryError, "Failed to allocate memory for object.");
        return NULL;
    }
    return result;
}

static PyObject*
PyTree_neighbor_get_count(PyTree* self)
{
    long count;
    struct KDTree* tree = self->tree;
    PyObject* result;
    count = KDTree_neighbor_get_count(tree);
#if PY_MAJOR_VERSION >= 3
    result = PyLong_FromLong(count);
#else
    result = PyInt_FromLong(count);
#endif
    if (!result)
    {
        PyErr_SetString (PyExc_MemoryError, "Failed to allocate memory for object.");
        return NULL;
    }
    return result;
}

#define COPY2DARRAY(ctype) \
    for (i = 0; i < n; i++) { \
        for (j = 0; j < m; j++) { \
            coords[i*m+j] = *(ctype *) (p+i*rowstride+j*colstride); \
        } \
    }

static PyObject*
PyTree_set_data(PyTree* self, PyObject* args)
{
    float* coords;
    Py_ssize_t n, m, i, j;
    PyObject *obj;
    struct KDTree* tree = self->tree;
    int ok;
    Py_ssize_t rowstride, colstride;
    const char* p;
    const int flags = PyBUF_FORMAT | PyBUF_STRIDES;
    char datatype;
    Py_buffer view;

    if(!PyArg_ParseTuple(args, "O:KDTree_set_data",&obj)) return NULL;

    if (PyObject_GetBuffer(obj, &view, flags) == -1) return NULL;
    if (view.ndim != 2) {
        PyErr_SetString(PyExc_RuntimeError, "Array must be two-dimensional");
        return NULL;
    }
    n = view.shape[0];
    m = view.shape[1];
    rowstride = view.strides[0];
    colstride = view.strides[1];
    /* coord_data is deleted by the KDTree object */
    coords = malloc(m*n*sizeof(float));
    if (!coords) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory for coordinates.");
        goto exit;
    }
    p = view.buf;
    datatype = view.format[0];
    switch (datatype) {
        case '@':
        case '=':
        case '<':
        case '>':
        case '!': datatype = view.format[1]; break;
        default: break;
    }
    switch (datatype) {
        case 'd': COPY2DARRAY(double); break;
        case 'f': COPY2DARRAY(float); break;
        case 'i': COPY2DARRAY(int); break;
        case 'I': COPY2DARRAY(unsigned int); break;
        case 'l': COPY2DARRAY(long); break;
        case 'L': COPY2DARRAY(unsigned long); break;
        default:
            PyErr_Format(PyExc_RuntimeError,
                "array should contain numerical data (format character was %c).",
                datatype);
            goto exit;
    }
    ok = KDTree_set_data(tree, coords, n);
    if (!ok)
    {
        PyErr_SetString (PyExc_MemoryError, "Failed to allocate memory for nodes.");
        goto exit;
    }
    PyBuffer_Release(&view);
    Py_INCREF(Py_None);
    return Py_None;

exit:
    PyBuffer_Release(&view);
    if (coords) free(coords);
    return NULL;
}

#define COPY1DARRAY(ctype) \
    for (i = 0; i < n; i++) { \
        coords[i] = *(ctype *) (p+i*stride); \
    }

static PyObject*
PyTree_search_center_radius(PyTree* self, PyObject* args)
{
    PyObject *obj;
    double radius;
    long int n, i;
    float *coords;
    struct KDTree* tree = self->tree;
    int ok;
    const int flags = PyBUF_FORMAT | PyBUF_STRIDES;
    Py_ssize_t stride;
    Py_buffer view;
    char datatype;
    const char* p;

    if(!PyArg_ParseTuple(args, "Od:KDTree_search_center_radius", &obj ,&radius))
        return NULL;

    if(radius <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "Radius must be positive.");
        return NULL;
    }

    if (PyObject_GetBuffer(obj, &view, flags) == -1) return NULL;
    if (view.ndim != 1) {
        PyErr_SetString(PyExc_RuntimeError, "Array must be one-dimensional");
        return NULL;
    }
    n = view.shape[0];
    stride = view.strides[0];
    /* coord_data is deleted by the KDTree object */
    coords = malloc(n*sizeof(float));
    if (!coords) {
        PyErr_NoMemory();
        goto exit;
    }
    p = view.buf;
    datatype = view.format[0];
    switch (datatype) {
        case '@':
        case '=':
        case '<':
        case '>':
        case '!': datatype = view.format[1]; break;
        default: break;
    }
    switch (datatype) {
        case 'd': COPY1DARRAY(double); break;
        case 'f': COPY1DARRAY(float); break;
        case 'i': COPY1DARRAY(int); break;
        case 'I': COPY1DARRAY(unsigned int); break;
        case 'l': COPY1DARRAY(long); break;
        case 'L': COPY1DARRAY(unsigned long); break;
        default:
            PyErr_Format(PyExc_RuntimeError,
                "array should contain numerical data (format character was %c.",
                datatype);
            goto exit;
    }
    ok = KDTree_search_center_radius(tree, coords, radius);
    if (!ok) {
        PyErr_NoMemory();
        goto exit;
    }
    PyBuffer_Release(&view);
    Py_INCREF(Py_None);
    return Py_None;

exit:
    PyBuffer_Release(&view);
    if (coords) free(coords);
    return NULL;
}

static PyObject*
PyTree_neighbor_search(PyTree* self, PyObject* args)
{
    int ok;
    double radius;
    struct KDTree* tree = self->tree;
    struct Neighbor* neighbors;
    struct Neighbor* pp, *qq;
    PyObject* list;
    Py_ssize_t i, n;

    if(!PyArg_ParseTuple(args, "d:KDTree_neighbor_search", &radius))
        return NULL;

    if(radius <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "Radius must be positive.");
        return NULL;
    }

    ok = KDTree_neighbor_search(tree, radius, &neighbors);
    if (!ok)
    {
        PyErr_SetString(PyExc_MemoryError,
            "calculation failed due to lack of memory");
        return NULL;
    }

    pp = neighbors;
    n = 0;
    while (pp)
    {
        n+=1;
        pp = pp->next;
    }

    list = PyList_New(n);
    if (list)
    {
        PyNeighbor* p;
        pp = neighbors;
        for (i = 0; i < n; i++)
        {
            p = (PyNeighbor*) PyNeighborType.tp_alloc(&PyNeighborType, 0);
            if(!p)
            {
                PyErr_SetString(PyExc_MemoryError,
                    "could not create node for return value");
                Py_DECREF(list);
                return NULL;
            }
            p->neighbor = *pp;
            PyList_SET_ITEM(list, i, (PyObject*)p);
            qq = pp->next;
            free(pp);
            pp = qq;
        }
    }

    return list;
}

static PyObject*
PyTree_neighbor_simple_search(PyTree* self, PyObject* args)
{
    int ok;
    double radius;
    struct KDTree* tree = self->tree;
    struct Neighbor* neighbors;
    struct Neighbor* pp, *qq;
    PyObject* list;
    Py_ssize_t i, n;

    if(!PyArg_ParseTuple(args, "d:KDTree_neighbor_simple_search", &radius))
        return NULL;

    if(radius <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "Radius must be positive.");
        return NULL;
    }

    ok = KDTree_neighbor_simple_search(tree, radius, &neighbors);
    if (!ok)
    {
        PyErr_SetString(PyExc_MemoryError,
            "calculation failed due to lack of memory");
        return NULL;
    }

    pp = neighbors;
    n = 0;
    while (pp)
    {
        n+=1;
        pp = pp->next;
    }

    list = PyList_New(n);
    if (list)
    {
        PyNeighbor* p;
        pp = neighbors;
        for (i = 0; i < n; i++)
        {
            p = (PyNeighbor*) PyNeighborType.tp_alloc(&PyNeighborType, 0);
            if(!p)
            {
                PyErr_SetString(PyExc_MemoryError,
                    "could not create node for return value");
                Py_DECREF(list);
                return NULL;
            }
            p->neighbor = *pp;
            PyList_SET_ITEM(list, i, (PyObject*)p);
            qq = pp->next;
            free(pp);
            pp = qq;
        }
    }

    return list;
}

static char PyTree_get_indices__doc__[] =
"returns indices of coordinates within radius as a Numpy array\n";

static PyObject *PyTree_get_indices(PyTree *self, PyObject* args)
{
    struct KDTree* tree = self->tree;
    const int flags = PyBUF_C_CONTIGUOUS | PyBUF_FORMAT;
    char datatype;
    Py_buffer view;
    PyObject* object;

    if (!PyArg_ParseTuple(args, "O:KDTree_get_indices", &object)) return NULL;
    if (PyObject_GetBuffer(object, &view, flags) == -1)
        return NULL;
    datatype = view.format[0];
    switch (datatype) {
        case '@':
        case '=':
        case '<':
        case '>':
        case '!': datatype = view.format[1]; break;
        default: break;
    }
    if (datatype != 'l') {
        PyErr_Format(PyExc_RuntimeError,
            "array has incorrect data format ('%c', expected 'l')", datatype);
        PyBuffer_Release(&view);
        return NULL;
    }
    else if (view.ndim != 1) {
        PyErr_Format(PyExc_ValueError,
            "array has incorrect rank (%d expected 1)", view.ndim);
        PyBuffer_Release(&view);
        return NULL;
    }
    /* copy the data into the Numpy data pointer */
    KDTree_copy_indices(tree, (long int *) view.buf);
    PyBuffer_Release(&view);
    Py_INCREF(Py_None);
    return Py_None;
}

static char PyTree_get_radii__doc__[] =
"returns distances of coordinates within radius as a Numpy array.\n";

static PyObject *PyTree_get_radii(PyTree *self, PyObject* args)
{
    PyObject* object;
    const int flags = PyBUF_C_CONTIGUOUS | PyBUF_FORMAT;
    char datatype;
    Py_buffer view;
    struct KDTree* tree = self->tree;

    if (!PyArg_ParseTuple(args, "O:KDTree_get_radii", &object)) return NULL;
    if (PyObject_GetBuffer(object, &view, flags) == -1)
        return NULL;
    datatype = view.format[0];
    switch (datatype) {
        case '@':
        case '=':
        case '<':
        case '>':
        case '!': datatype = view.format[1]; break;
        default: break;
    }
    if (datatype != 'f') {
        PyErr_Format(PyExc_RuntimeError,
            "array has incorrect data format ('%c', expected 'f')", datatype);
        PyBuffer_Release(&view);
        return NULL;
    }
    else if (view.ndim != 1) {
        PyErr_Format(PyExc_ValueError,
            "array has incorrect rank (%d expected 1)", view.ndim);
        PyBuffer_Release(&view);
        return NULL;
    }

    /* copy the data into the Numpy data pointer */
    KDTree_copy_radii(tree, (float *) view.buf);

    PyBuffer_Release(&view);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef PyTree_methods[] = {
    {"get_count", (PyCFunction)PyTree_get_count, METH_NOARGS, NULL},
    {"set_data", (PyCFunction)PyTree_set_data, METH_VARARGS, NULL},
    {"search_center_radius", (PyCFunction)PyTree_search_center_radius, METH_VARARGS, NULL},
    {"neighbor_get_count", (PyCFunction)PyTree_neighbor_get_count, METH_NOARGS, NULL},
    {"neighbor_search", (PyCFunction)PyTree_neighbor_search, METH_VARARGS, NULL},
    {"neighbor_simple_search", (PyCFunction)PyTree_neighbor_simple_search, METH_VARARGS, NULL},
    {"get_indices", (PyCFunction)PyTree_get_indices, METH_VARARGS, PyTree_get_indices__doc__},
    {"get_radii", (PyCFunction)PyTree_get_radii, METH_VARARGS, PyTree_get_radii__doc__},
    {NULL}  /* Sentinel */
};

static char PyTree_doc[] = "C KDTree.\n";

static PyTypeObject PyTreeType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "C KDTree",                  /*tp_name*/
    sizeof(PyTree),              /*tp_basicsize*/
    0,                           /*tp_itemsize*/
    (destructor)PyTree_dealloc,  /*tp_dealloc*/
    0,                           /*tp_print*/
    0,                           /*tp_getattr*/
    0,                           /*tp_setattr*/
    0,                           /*tp_compare*/
    0,                           /*tp_repr*/
    0,                           /*tp_as_number*/
    0,                           /*tp_as_sequence*/
    0,                           /*tp_as_mapping*/
    0,                           /*tp_hash */
    0,                           /*tp_call*/
    0,                           /*tp_str*/
    0,                           /*tp_getattro*/
    0,                           /*tp_setattro*/
    0,                           /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,          /*tp_flags*/
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
    (initproc)PyTree_init,       /* tp_init */
};

/* ========================================================================== */
/* -- Initialization -------------------------------------------------------- */
/* ========================================================================== */

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_CKDTree",
        NULL,
        -1,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL
};

PyObject *
PyInit__CKDTree(void)

#else

void
init_CKDTree(void)
#endif
{
  PyObject *module;

  PyTreeType.tp_new = PyType_GenericNew;
  PyNeighborType.tp_new = PyType_GenericNew;
  if (PyType_Ready(&PyTreeType) < 0)
#if PY_MAJOR_VERSION >= 3
      return NULL;
#else
      return;
#endif
  if (PyType_Ready(&PyNeighborType) < 0)
#if PY_MAJOR_VERSION >= 3
      return NULL;
#else
      return;
#endif

#if PY_MAJOR_VERSION >= 3
  module = PyModule_Create(&moduledef);
  if (module==NULL) return NULL;
#else
  module = Py_InitModule("_CKDTree", NULL);
  if (module==NULL) return;
#endif

  Py_INCREF(&PyTreeType);
  Py_INCREF(&PyNeighborType);
  PyModule_AddObject(module, "KDTree", (PyObject*) &PyTreeType);
  PyModule_AddObject(module, "Neighbor", (PyObject*) &PyNeighborType);

  if (PyErr_Occurred()) Py_FatalError("can't initialize module _CKDTree");
#if PY_MAJOR_VERSION >= 3
  return module;
#endif
}
