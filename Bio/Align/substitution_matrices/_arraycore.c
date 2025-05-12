#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyTypeObject *basetype = NULL;

typedef struct {
    int value;  // your custom field
    PyObject* alphabet;
} Fields;

static int get_value(PyObject* self) {
    Fields* fields = (Fields*)((intptr_t)self + basetype->tp_basicsize);
    return fields->value;
}

static PyObject *Array_get_value(PyObject *self, void *closure) {
    Fields* fields = (Fields*)((intptr_t)self + basetype->tp_basicsize);
    return PyLong_FromLong(fields->value);
}

static int Array_set_value(PyObject *self, PyObject *arg, void *closure) {
    Fields* fields = (Fields*)((intptr_t)self + basetype->tp_basicsize);
    long val = PyLong_AsLong(arg);
    if (PyErr_Occurred())
        return -1;
    fields->value = (int)val;
    return 0;
}

static PyObject *Array_get_alphabet(PyObject *self, void *closure) {
    Fields* fields = (Fields*)((intptr_t)self + basetype->tp_basicsize);
    PyObject* alphabet = fields->alphabet;
    Py_INCREF(alphabet);
    return alphabet;
}

static int Array_set_alphabet(PyObject *self, PyObject *arg, void *closure) {
    Fields* fields = (Fields*)((intptr_t)self + basetype->tp_basicsize);
    PyObject* alphabet = fields->alphabet;
    if (!PySequence_Check(arg)) {
        PyErr_SetString(PyExc_TypeError,
                        "argument does not provide the sequence protocol");
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
    Py_XDECREF(alphabet);
    Py_INCREF(arg);
    fields->alphabet = arg;
    return 0;
}

static PyObject*
Array_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    if (!basetype) {
        PyErr_SetString(PyExc_RuntimeError, "base type was not intialized");
        return NULL;
    }
    if (!basetype->tp_new) {
        PyErr_SetString(PyExc_RuntimeError, "base type does not have tp_new");
        return NULL;
    }

    PyObject* alphabet = NULL;
    if (PyTuple_GET_SIZE(args) != 3) {
        // shape, dtype, alphabet
        PyErr_Format(PyExc_TypeError,
                     "Array() takes 3 positional arguments but %d were given",
                     PyTuple_GET_SIZE(args));
        return NULL;
    }

    alphabet = PyTuple_GET_ITEM(args, 2);
    if (!alphabet) {
        PyErr_SetString(PyExc_RuntimeError, "failed to get alphabet argument");
        return NULL;
    }

    args = PyTuple_GetSlice(args, 0, PyTuple_GET_SIZE(args) - 1);
    if (!args)
        return NULL;

    PyObject *obj = basetype->tp_new(type, args, kwds);
    Py_DECREF(args);

    if (!obj)
        return NULL;

    Fields* fields = (Fields*)((intptr_t)obj + basetype->tp_basicsize);
    Py_INCREF(alphabet);
    fields->alphabet = alphabet;

    return obj;
}

static void
Array_dealloc(PyObject *self)
{
    Fields* fields = (Fields*)((intptr_t)self + basetype->tp_basicsize);
    /* fields->alphabet may be NULL if this instance was created by numpy
     * and __array_finalize__ somehow failed.
     */
    Py_XDECREF(fields->alphabet);
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
    {"value", (getter)Array_get_value, (setter)Array_set_value, "int value", NULL},
    {"alphabet2", (getter)Array_get_alphabet, (setter)Array_set_alphabet, "alphabet", NULL},
    {NULL}
};

static PyObject *init_capsule(PyObject *self, PyObject *args) {
    return PyCapsule_New((void *)get_value, "_arraycore.get_value", NULL);
}

// Define the methods of your class (if any)
static PyMethodDef methods[] = {
    {"get_value_function", (PyCFunction)init_capsule, METH_NOARGS, "Get value function capsule"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "_arraycore",
    "Base module defining the SubstitutionMatrix base class",
    -1,
    methods,
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
        {Py_tp_new, Array_new},
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
    PyObject *class_type = PyType_FromSpecWithBases(&class_spec, bases);
    Py_DECREF(bases);
    if (!class_type)
        return NULL;

    // Add it to the module
    if (PyModule_AddObject(mod, "SubstitutionMatrix", class_type) < 0) {
        Py_DECREF(class_type);
        Py_DECREF(mod);
        return NULL;
    }

    return mod;
}
