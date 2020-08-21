#define PY_SSIZE_T_CLEAN

#include "Python.h"

static PyTypeObject SeqType;

typedef struct {
    PyObject_HEAD
    PyObject* data;
} SeqObject;

static SeqObject *
Seq_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    SeqObject *object;

    object = (SeqObject*)type->tp_alloc(type, 0);
    if (!object) return NULL;

    Py_INCREF(Py_None);
    object->data = Py_None;

    return object;
}

static void
Seq_dealloc(SeqObject *op)
{
    Py_DECREF(op->data);
    Py_TYPE(op)->tp_free(op);
}

static PyObject *
Seq_getattro(SeqObject *self, PyObject *name)
{
    PyObject* data = self->data;
    if (PyUnicode_CompareWithASCIIString(name, "data") == 0) {
        Py_INCREF(data);
        return data;
    }
    return PyObject_GenericGetAttr((PyObject*)self, name);
}

static int
Seq_setattro(SeqObject *self, PyObject *name, PyObject *value)
{
    if (PyUnicode_CompareWithASCIIString(name, "data") == 0) {
        Py_DECREF(self->data);
        Py_INCREF(value);
        self->data = value;
        return 0;
    }
    return PyObject_GenericSetAttr((PyObject*)self, name, value);
}

static int
Seq_buffer_getbuffer(SeqObject *self, Py_buffer *view, int flags)
{
    PyObject* data = self->data;
    return PyObject_GetBuffer(data, view, flags);
}

static PyBufferProcs Seq_as_buffer = {
    (getbufferproc)Seq_buffer_getbuffer,
    NULL,
};

static PyObject *
Seq_reduce(SeqObject *self)
{
    return Py_BuildValue("(O(O))", Py_TYPE(self), self->data);
}

static PyMethodDef
Seq_methods[] = {
    {"__reduce__", (PyCFunction)Seq_reduce,  METH_NOARGS},
    {NULL,     NULL}                         /* sentinel */
};

PyDoc_STRVAR(Seq_doc, "Seq() -> Seq");

static PyTypeObject SeqType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "Seq",                                      /* tp_name */
    sizeof(SeqObject),                          /* tp_basicsize */
    0,                                          /* tp_itemsize */
    (destructor)Seq_dealloc,                    /* tp_dealloc */
    0,                                          /* tp_print */
    0,                                          /* tp_getattr */
    0,                                          /* tp_setattr */
    0,                                          /* tp_reserved */
    0,                                          /* tp_repr */
    0,                                          /* tp_as_number */
    0,                                          /* tp_as_sequence */
    0,                                          /* tp_as_mapping */
    0,                                          /* tp_hash */
    0,                                          /* tp_call */
    0,                                          /* tp_str */
    (getattrofunc)Seq_getattro,                 /* tp_getattro */
    (setattrofunc)Seq_setattro,                 /* tp_setattro */
    &Seq_as_buffer,                             /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE |
        Py_TPFLAGS_BYTES_SUBCLASS,              /* tp_flags */
    Seq_doc,                                    /* tp_doc */
    0,                                          /* tp_traverse */
    0,                                          /* tp_clear */
    0,                                          /* tp_richcompare */
    0,                                          /* tp_weaklistoffset */
    0,                                          /* tp_iter */
    0,                                          /* tp_iternext */
    Seq_methods,                                /* tp_methods */
    0,                                          /* tp_members */
    0,                                          /* tp_getset */
    &PyBaseObject_Type,                         /* tp_base */
    0,                                          /* tp_dict */
    0,                                          /* tp_descr_get */
    0,                                          /* tp_descr_set */
    0,                                          /* tp_dictoffset */
    0,                                          /* tp_init */
    0,                                          /* tp_alloc */
    (newfunc)Seq_new,                           /* tp_new */
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_seqobject",
        NULL,
        -1,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL
};

PyObject *
PyInit__seqobject(void)
{
    PyObject *module;

    if (PyType_Ready(&SeqType) < 0)
        return NULL;

    module = PyModule_Create(&moduledef);
    if (module == NULL) return NULL;

    Py_INCREF(&SeqType);
    if (PyModule_AddObject(module, "Seq", (PyObject*) &SeqType) == 0)
        return module;

    Py_DECREF(module);
    Py_DECREF(&SeqType);
    return NULL;
}
