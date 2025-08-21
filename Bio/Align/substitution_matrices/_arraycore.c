#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "_arraycore.h"


static PyTypeObject Array_Type;


static void
Array_dealloc(PyObject *self)
{
    PyTypeObject* basetype = Array_Type.tp_base;
    Fields* fields = (Fields*)((intptr_t)self + basetype->tp_basicsize);
    /* fields->alphabet may be NULL if this instance was created by numpy
     * and __array_finalize__ somehow failed.
     */
    Py_XDECREF(fields->alphabet);
    /* PyBuffer_Release won't do anything if fields->mapping.obj is NULL. */
    PyMem_Free(fields->mapping.buf);
    basetype->tp_dealloc(self);
}

static PyObject *
Array_finalize(PyObject *self, PyObject *obj)
{
    if (obj == Py_None) Py_RETURN_NONE;
    if (Py_TYPE(self) != Py_TYPE(obj)) {
        PyErr_SetString(PyExc_RuntimeError,
                        "__array_finalize__ argument is not an Array object");
        return NULL;
    }
    PyTypeObject* basetype = Array_Type.tp_base;
    Fields* self_fields = (Fields*)((intptr_t)self + basetype->tp_basicsize);
    const Fields* obj_fields = (Fields*)((intptr_t)obj + basetype->tp_basicsize);
    PyObject* alphabet = obj_fields->alphabet;
    if (alphabet) {
        Py_INCREF(alphabet);
        self_fields->alphabet = obj_fields->alphabet;
    }
    Py_RETURN_NONE;
}

static PyMethodDef Array_methods[] = {
    {"__array_finalize__", (PyCFunction)Array_finalize, METH_O,
     "Called by NumPy to finalize new views or copies"},
    {NULL}  /* Sentinel */
};

static PyObject *Array_get_alphabet(PyObject *self, void *closure) {
    PyTypeObject* basetype = Array_Type.tp_base;
    Fields* fields = (Fields*)((intptr_t)self + basetype->tp_basicsize);
    PyObject* alphabet = fields->alphabet;
    if (!alphabet) Py_RETURN_NONE;
    Py_INCREF(alphabet);
    return alphabet;
}

static int Array_set_alphabet(PyObject *self, PyObject *arg, void *closure) {
    Py_buffer view;
    const Py_ssize_t length = PySequence_Size(arg);
    PyTypeObject* basetype = Array_Type.tp_base;
    Fields* fields = (Fields*)((intptr_t)self + basetype->tp_basicsize);
    if (fields->alphabet) {
        PyErr_SetString(PyExc_ValueError, "the alphabet has already been set.");
        return -1;
    }
    if (!PySequence_Check(arg)) {
        PyErr_SetString(PyExc_TypeError,
            "alphabet must support the sequence protocol (e.g.,\n"
            "strings, lists, and tuples can be valid alphabets).");
        return -1;
    }
    if (PyObject_GetBuffer(self, &view, PyBUF_STRIDES) != 0) {
        PyErr_SetString(PyExc_RuntimeError, "failed to access matrix buffer");
        return -1;
    }
    switch (view.ndim) {
        case 1:
            if (view.shape[0] == length) break;
            PyErr_Format(PyExc_ValueError,
                "alphabet length %zd is inconsistent with array size "
                "%zd", length, view.shape[0]);
            PyBuffer_Release(&view);
            return -1;
        case 2:
            if ((view.shape[0] == length && view.shape[1] == length)
             || (view.shape[0] == length && view.shape[1] == 1)
             || (view.shape[0] == 1 && view.shape[1] == length)) break;
            PyErr_Format(PyExc_ValueError,
                "alphabet length %zd is inconsistent with array size "
                "(%zd, %zd)", length, view.shape[0], view.shape[1]);
            PyBuffer_Release(&view);
            return -1;
        default:
            PyErr_Format(PyExc_ValueError,
                         "substitution matrix has incorrect rank %d "
                         "(expected 1 or 2)", view.ndim);
            PyBuffer_Release(&view);
            return -1;
    }
    PyBuffer_Release(&view);
    if (PyUnicode_Check(arg)) {
        /* initialize mapping if alphabet is a string: */
        Py_ssize_t mapping_size;
        void* characters = PyUnicode_DATA(arg);
        int kind = PyUnicode_KIND(arg);
        int* mapping;
        Py_ssize_t i;
        switch (kind) {
            case PyUnicode_1BYTE_KIND: {
                mapping_size = 1 << 8 * sizeof(Py_UCS1);
                break;
            }
            case PyUnicode_2BYTE_KIND: {
                mapping_size = 1 << 8 * sizeof(Py_UCS2);
                break;
            }
            case PyUnicode_4BYTE_KIND: {
                mapping_size = 0x110000;  /* Maximum code point in Unicode 6.0
                                           * is 0x10ffff = 1114111 */
                break;
            }
            default:
                PyErr_SetString(PyExc_ValueError, "could not interpret alphabet");
                return -1;
        }
        mapping = PyMem_Malloc(mapping_size*sizeof(int));
        if (!mapping) return -1;
        for (i = 0; i < mapping_size; i++) mapping[i] = MISSING_LETTER;
        for (i = 0; i < length; i++) {
            Py_UCS4 character = PyUnicode_READ(kind, characters, i);
            if (mapping[character] != MISSING_LETTER) {
                PyObject* c = PyUnicode_FromKindAndData(kind, &character, 1);
                PyErr_Format(PyExc_ValueError,
                             "alphabet contains '%S' more than once", c);
                Py_XDECREF(c);
                PyMem_Free(mapping);
                return -1;
            }
            mapping[character] = i;
        }
        if (PyBuffer_FillInfo(&fields->mapping,
                              NULL,
                              mapping,
                              mapping_size * sizeof(int),
                              0,
                              PyBUF_SIMPLE) == -1) {
            PyMem_Free(mapping);
            return -1;
        }
        fields->mapping.itemsize = sizeof(int);
    }
    Py_INCREF(arg);
    fields->alphabet = arg;
    return 0;
}

static PyGetSetDef Array_getset[] = {
    {"alphabet", (getter)Array_get_alphabet, (setter)Array_set_alphabet, "alphabet", NULL},
    {NULL}
};

static PyTypeObject Array_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_arraycore.Array",
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_dealloc = (destructor)Array_dealloc,
    .tp_methods = Array_methods,
    .tp_getset = Array_getset,
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_arraycore",
    .m_doc = "Base module defining the Array base class",
    .m_size = -1,
};

PyMODINIT_FUNC
PyInit__arraycore(void)
{
    int result;
    PyObject *basemodule;
    PyObject *baseclass;
    PyTypeObject* basetype;

    PyObject *mod = PyModule_Create(&module);
    if (!mod) return NULL;

    // Import the module containing the base class
    basemodule = PyImport_ImportModule("numpy");
    if (!basemodule)
        return NULL;

    // Get the base class
    baseclass = PyObject_GetAttrString(basemodule, "ndarray");
    Py_DECREF(basemodule);
    if (!baseclass || !PyType_Check(baseclass)) {
        Py_XDECREF(basemodule);
        PyErr_SetString(PyExc_RuntimeError, "Failed to get numpy.ndarray");
        return NULL;
    }

    basetype = (PyTypeObject *)baseclass;
    if (!(basetype->tp_flags & Py_TPFLAGS_BASETYPE)) {
        PyErr_SetString(PyExc_RuntimeError,
                        "numpy ndarray class is not subclassable");
        return NULL;
    }
    if (basetype->tp_itemsize != 0) {
        PyErr_Format(PyExc_RuntimeError,
                     "expected numpy arrays to have tp_itemsize 0 (found %zd)",
                     basetype->tp_itemsize);
        return NULL;
    }
    if (!basetype->tp_new) {
        PyErr_SetString(PyExc_RuntimeError,
                        "numpy ndarray class does not have tp_new");
        return NULL;
    }

    Array_Type.tp_basicsize = basetype->tp_basicsize + sizeof(Fields);
    Array_Type.tp_base = (PyTypeObject *)baseclass;

    if (PyType_Ready(&Array_Type) < 0)
        return NULL;

    result = PyModule_AddObjectRef(mod, "Array", (PyObject*)&Array_Type);
    Py_DECREF(&Array_Type);
    if (result == -1) {
        Py_DECREF(mod);
        return NULL;
    }

    return mod;
}
