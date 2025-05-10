#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyTypeObject *basetype = NULL;

typedef struct {
    int value;  // your custom field
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
    PyObject *obj = basetype->tp_new(type, args, kwds);
    if (!obj)
        return NULL;
    return obj;
}

static PyGetSetDef Array_getset[] = {
    {"value", (getter)Array_get_value, (setter)Array_set_value, "int value", NULL},
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
