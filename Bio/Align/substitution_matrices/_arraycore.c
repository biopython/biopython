#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "_arraycore.h"

#define MISSING_LETTER -1

static PyTypeObject *basetype = NULL;
static PyObject *class_type = NULL;

typedef struct {
    PyObject* alphabet;
    Py_buffer mapping;
} Fields;

static Py_buffer* Array_get_mapping_buffer(PyObject* self) {
    if (!PyObject_IsInstance(self, class_type)) return NULL;
    Fields* fields = (Fields*)((intptr_t)self + basetype->tp_basicsize);
    Py_buffer* mapping = &fields->mapping;
    Py_XINCREF(mapping->obj);
    return mapping;
}

static PyObject *Array_get_alphabet(PyObject *self, void *closure) {
    Fields* fields = (Fields*)((intptr_t)self + basetype->tp_basicsize);
    PyObject* alphabet = fields->alphabet;
    if (!alphabet) Py_RETURN_NONE;
    Py_INCREF(alphabet);
    return alphabet;
}

static int Array_set_alphabet(PyObject *self, PyObject *arg, void *closure) {
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
    const Py_ssize_t length = PySequence_Size(arg);
    Py_buffer view;
    const int flag = PyBUF_STRIDES;
    if (PyObject_GetBuffer(self, &view, flag) != 0) {
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
        /* we cannot use mapping if alphabet is not a string */
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
                              self,
                              mapping,
                              mapping_size * sizeof(int),
                              0,
                              PyBUF_SIMPLE) == -1) {
            PyMem_Free(mapping);
            return -1;
        }
        fields->mapping.itemsize = sizeof(int);
    }
    else {
        Py_INCREF(self);
        fields->mapping.obj = self;
        fields->mapping.itemsize = sizeof(int);
        fields->mapping.len = length * sizeof(int);
    }
    Py_INCREF(arg);
    fields->alphabet = arg;
    return 0;
}

static void
Array_dealloc(PyObject *self)
{
    Fields* fields = (Fields*)((intptr_t)self + basetype->tp_basicsize);
    /* fields->alphabet may be NULL if this instance was created by numpy
     * and __array_finalize__ somehow failed.
     */
    Py_XDECREF(fields->alphabet);
    /* PyBuffer_Release won't do anything if fields->mapping is NULL. */
    PyBuffer_Release(&fields->mapping);
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

static PyGetSetDef Array_getset[] = {
    {"alphabet", (getter)Array_get_alphabet, (setter)Array_set_alphabet, "alphabet", NULL},
    {NULL}
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "_arraycore",
    "Base module defining the SubstitutionMatrix base class",
    -1,
    NULL,
    NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC
PyInit__arraycore(void)
{
    PyObject *mod = PyModule_Create(&module);
    if (!mod)
        return NULL;

    // Import the module containing the base class
    PyObject *basemodule = PyImport_ImportModule("numpy");
    if (!basemodule)
        return NULL;

    // Get the base class
    PyObject *baseclass = PyObject_GetAttrString(basemodule, "ndarray");
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
        basetype = NULL;
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

    // Create a tuple of bases
    PyObject *bases = PyTuple_Pack(1, baseclass);
    Py_DECREF(baseclass);
    if (!bases)
        return NULL;

    PyType_Slot Array_slots[] = {
        {Py_tp_dealloc, Array_dealloc},
        {Py_tp_methods, Array_methods},
        {Py_tp_getset, Array_getset},
        {0, 0}
    };

    // Define the type spec
    PyType_Spec class_spec = {
        .name = "_arraycore.SubstitutionMatrix",
        .basicsize = basetype->tp_basicsize + sizeof(Fields),
        .itemsize = 0,
        .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
        .slots = Array_slots,
    };

    // Create the type
    class_type = PyType_FromSpecWithBases(&class_spec, bases);
    Py_DECREF(bases);
    if (!class_type)
        return NULL;

    // Add it to the module
    if (PyModule_AddObject(mod, "SubstitutionMatrix", class_type) < 0) {
        Py_DECREF(class_type);
        Py_DECREF(mod);
        return NULL;
    }

    // Ensure that Array_get_mapping_buffer has the correct signature
    // by comparing to Array_get_mapping_buffer_signature:
    // (the compiler will complain if the signature of Array_get_mapping_buffer
    // is inconsistent with Array_get_mapping_buffer_signature as defined in
    // the header file):
    Array_get_mapping_buffer_signature function = Array_get_mapping_buffer;
    PyObject* capsule = PyCapsule_New((void*)function,
        "Bio.Align.substitution_matrices._arraycore.mapping_buffer_capsule",
        NULL);
    if (!capsule) {
        Py_DECREF(class_type);
        Py_DECREF(mod);
        return NULL;
    }
    if (PyModule_AddObject(mod, "mapping_buffer_capsule", capsule) < 0) {
        Py_DECREF(class_type);
        Py_DECREF(mod);
        Py_DECREF(capsule);
        return NULL;
    }

    return mod;
}
