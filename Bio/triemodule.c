#include <Python.h>
#include <marshal.h>
#include "trie.h"

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

/* Must define Py_TYPE for Python 2.5 or older */
#ifndef Py_TYPE
# define Py_TYPE(o) ((o)->ob_type)
#endif

/* Must define PyVarObject_HEAD_INIT for Python 2.5 or older */
#ifndef PyVarObject_HEAD_INIT
#define PyVarObject_HEAD_INIT(type, size) \
PyObject_HEAD_INIT(type) size,
#endif

static PyTypeObject Trie_Type;

typedef struct {
    PyObject_HEAD
    Trie* trie;
} trieobject;

static PyObject*
trie_trie(PyObject* self, PyObject* args)
{
    trieobject* trieobj;
    Trie* trie;

    if (!PyArg_ParseTuple(args,":trie"))
        return NULL;
    if(!(trie = Trie_new()))
	return PyErr_NoMemory();
    if(!(trieobj = PyObject_New(trieobject, &Trie_Type)))
	return NULL;
    trieobj->trie = trie;
    return (PyObject*)trieobj;
}

static void
_decref_objects(const char *key, const void *value, void *data)
{
    Py_DECREF((PyObject *)value);
}

static void
trie_dealloc(PyObject* self)
{
    trieobject *mp = (trieobject *)self;
    Trie_iterate(mp->trie, _decref_objects, NULL);
    Trie_del(mp->trie);
    PyObject_Del(self);
}

static Py_ssize_t
trie_length(trieobject *mp)
{
    return Trie_len(mp->trie);
}

static PyObject *
trie_subscript(trieobject *mp, PyObject *py_key)
{
    const char *key;
    PyObject *py_value;
#ifdef IS_PY3K
    PyObject* bytes;
#endif

    /* Make sure key is a string. */
#ifdef IS_PY3K
    if(!PyUnicode_Check(py_key)) {
#else
    if(!PyString_Check(py_key)) {
#endif
	PyErr_SetString(PyExc_TypeError, "key must be a string");
	return NULL;
    }
#ifdef IS_PY3K
    bytes = PyUnicode_AsASCIIString(py_key);
    if(!bytes) {
	PyErr_SetString(PyExc_TypeError, "key must be an ASCII string");
	return NULL;
    }
    key = PyBytes_AsString(bytes);
#else
    key = PyString_AS_STRING(py_key);
#endif
    py_value = Trie_get(mp->trie, key);
#ifdef IS_PY3K
    Py_DECREF(bytes);
#endif
    if(py_value == NULL)
	PyErr_SetString(PyExc_KeyError, key);
    else
	Py_INCREF(py_value);
    return py_value;
}

static int
trie_ass_sub(trieobject *mp, PyObject *py_key, PyObject *py_value)
{
    int result = -1;
    const char *key;
    PyObject *py_prev;
#ifdef IS_PY3K
    PyObject* bytes;
#endif

    /* Make sure key is a string. */
#ifdef IS_PY3K
    if(!PyUnicode_Check(py_key)) {
#else
    if(!PyString_Check(py_key)) {
#endif
	PyErr_SetString(PyExc_TypeError, "key must be a string");
	return -1;
    }
#ifdef IS_PY3K
    bytes = PyUnicode_AsASCIIString(py_key);
    if(!bytes) {
	PyErr_SetString(PyExc_TypeError, "key must be an ASCII string");
	return -1;
    }
    key = PyBytes_AsString(bytes);
#else
    key = PyString_AS_STRING(py_key);
#endif

    /* Check to see whether something already exists at that key.  If
       there's already an object there, then I will have to remove it.
    */
    py_prev = Trie_get(mp->trie, key);
    if(py_prev) {
	Py_DECREF(py_prev);
    }

    /* The client wants to delete a key from a dictionary.  The Trie
       API doesn't support this, so I will just overwrite it with
       NULL. */
    if(!py_value) {
	/* If the key doesn't exist, raise a KeyError. */
	if(!py_prev)
	    PyErr_SetString(PyExc_KeyError, key);
        else {
	    Trie_set(mp->trie, key, NULL);
            result = 0;
        }
    }
    /* The client wants to set a key in the dictionary. */
    else {
	Py_INCREF(py_value);
	if(Trie_set(mp->trie, key, py_value))
	    PyErr_SetString(PyExc_AssertionError, "error setting trie");
        else
            result = 0;
    }
#ifdef IS_PY3K
    Py_DECREF(bytes);
#endif
    return result;
}

static int trie_contains(trieobject *mp, PyObject* py_key)
{
    int result;
#ifdef IS_PY3K
    PyObject* bytes;
#endif
    const char *key;
    /* Make sure key is a string. */
#ifdef IS_PY3K
    if(!PyUnicode_Check(py_key)) {
#else
    if(!PyString_Check(py_key)) {
#endif
	PyErr_SetString(PyExc_TypeError, "key must be a string");
	return -1;
    }
#ifdef IS_PY3K
    bytes = PyUnicode_AsASCIIString(py_key);
    if(!bytes) {
	PyErr_SetString(PyExc_TypeError, "key must be an ASCII string");
	return -1;
    }
    key = PyBytes_AsString(bytes);
#else
    key = PyString_AS_STRING(py_key);
#endif
    result = Trie_has_key(mp->trie, key);
#ifdef IS_PY3K
    Py_DECREF(bytes);
#endif
    return result;
}

static char has_key__doc__[] =
"D.has_key(k) -> 1 if D has a key k, else 0";

static PyObject *
trie_has_key(trieobject *mp, PyObject *py_key)
{
    int has_key = trie_contains(mp, py_key);
    if (has_key==-1) return NULL;
#ifdef IS_PY3K
    return PyLong_FromLong((long)has_key);
#else
    return PyInt_FromLong((long)has_key);
#endif
}

static PyObject *
trie_has_key_onearg(trieobject *mp, PyObject *py_args)
{
    PyObject *py_arg;
    if(!PyArg_ParseTuple(py_args, "O", &py_arg))
	return NULL;
    return trie_has_key(mp, py_arg);
}


static char has_prefix__doc__[] =
"D.has_prefix(k) -> 1 if D has a prefix k, else 0";

static PyObject *
trie_has_prefix(trieobject *mp, PyObject *py_prefix)
{
    const char *prefix;
    int has_prefix;
#ifdef IS_PY3K
    PyObject* bytes;
#endif

    /* Make sure prefix is a string. */
#ifdef IS_PY3K
    if(!PyUnicode_Check(py_prefix)) {
#else
    if(!PyString_Check(py_prefix)) {
#endif
	PyErr_SetString(PyExc_TypeError, "prefix must be a string");
	return NULL;
    }
#ifdef IS_PY3K
    bytes = PyUnicode_AsASCIIString(py_prefix);
    if(!bytes) {
	PyErr_SetString(PyExc_TypeError, "prefix must be an ASCII string");
	return NULL;
    }
    prefix = PyBytes_AsString(bytes);
#else
    prefix = PyString_AS_STRING(py_prefix);
#endif
    has_prefix = Trie_has_prefix(mp->trie, prefix);
#ifdef IS_PY3K
    Py_DECREF(bytes);
    return PyLong_FromLong((long)has_prefix);
#else
    return PyInt_FromLong((long)has_prefix);
#endif
}

static PyObject *
trie_has_prefix_onearg(trieobject *mp, PyObject *py_args)
{
    PyObject *py_arg;
    if(!PyArg_ParseTuple(py_args, "O", &py_arg))
	return NULL;
    return trie_has_prefix(mp, py_arg);
}

static char with_prefix__doc__[] =
"D.with_prefix(prefix) -> list of D's keys that begins with prefix";

static void
_trie_with_prefix_helper(const char *key, const void *value, void *data)
{
    PyObject *py_list = (PyObject *)data;
    PyObject *py_key;

    if(PyErr_Occurred())
	return;

#ifdef IS_PY3K
    if(!(py_key = PyUnicode_FromFormat(key)))
#else
    if(!(py_key = PyString_FromString(key)))
#endif
	return;
    PyList_Append(py_list, py_key);
    Py_DECREF(py_key);
}

static PyObject *
trie_with_prefix(trieobject *mp, PyObject *py_prefix)
{
    const char *prefix;
    PyObject *py_list;
#ifdef IS_PY3K
    PyObject* bytes;
#endif

    /* Make sure prefix is a string. */
#ifdef IS_PY3K
    if(!PyUnicode_Check(py_prefix)) {
#else
    if(!PyString_Check(py_prefix)) {
#endif
	PyErr_SetString(PyExc_TypeError, "prefix must be a string");
	return NULL;
    }
#ifdef IS_PY3K
    bytes = PyUnicode_AsASCIIString(py_prefix);
    if(!bytes) {
	PyErr_SetString(PyExc_TypeError, "prefix must be an ASCII string");
	return NULL;
    }
    prefix = PyBytes_AsString(bytes);
#else
    prefix = PyString_AS_STRING(py_prefix);
#endif

    if(!(py_list = PyList_New(0)))
	return NULL;
    Trie_with_prefix(mp->trie, prefix,
		     _trie_with_prefix_helper, (void *)py_list);
#ifdef IS_PY3K
    Py_DECREF(bytes);
#endif
    if(PyErr_Occurred()) {
	Py_DECREF(py_list);
	return NULL;
    }
    return py_list;
}

static PyObject *
trie_with_prefix_onearg(trieobject *mp, PyObject *py_args)
{
    PyObject *py_arg;
    if(!PyArg_ParseTuple(py_args, "O", &py_arg))
	return NULL;
    return trie_with_prefix(mp, py_arg);
}


static char keys__doc__[] =
"D.keys() -> list of D's keys";

static void
_trie_keys_helper(const char *key, const void *value, void *data)
{
    PyObject *py_list = (PyObject *)data;
    PyObject *py_key;

    if(PyErr_Occurred())
	return;

#ifdef IS_PY3K
    if(!(py_key = PyUnicode_FromFormat(key)))
#else
    if(!(py_key = PyString_FromString(key)))
#endif
	return;
    PyList_Append(py_list, py_key);
    Py_DECREF(py_key);
}

static PyObject *
trie_keys(trieobject *mp)
{
    PyObject *py_list;

    if(!(py_list = PyList_New(0)))
	return NULL;
    Trie_iterate(mp->trie, _trie_keys_helper, (void *)py_list);
    if(PyErr_Occurred()) {
	Py_DECREF(py_list);
	return NULL;
    }
    return py_list;
}

static PyObject *
trie_keys_noargs(trieobject *mp, PyObject *py_args)
{
    if(PyTuple_Size(py_args) != 0) {
	PyErr_SetString(PyExc_ValueError, "no args expected");
	return NULL;
    }
    return trie_keys(mp);
}

static char values__doc__[] =
"D.values() -> list of D's values";

static void
_trie_values_helper(const char *key, const void *value, void *data)
{
    PyObject *py_list = (PyObject *)data;
    if(PyErr_Occurred())
	return;
    PyList_Append(py_list, (PyObject *)value);
}

static PyObject *
trie_values(trieobject *mp)
{
    PyObject *py_list;

    if(!(py_list = PyList_New(0)))
	return NULL;
    Trie_iterate(mp->trie, _trie_values_helper, (void *)py_list);
    if(PyErr_Occurred()) {
	Py_DECREF(py_list);
	return NULL;
    }
    return py_list;
}

static PyObject *
trie_values_noargs(trieobject *mp, PyObject *py_args)
{
    if(PyTuple_Size(py_args) != 0) {
	PyErr_SetString(PyExc_ValueError, "no args expected");
	return NULL;
    }
    return trie_values(mp);
}

static char get__doc__[] =
"D.get(k[,d]) -> D[k] if D.has_key(k), else d.  d defaults to None.";

static PyObject *
trie_get(trieobject *mp, PyObject *args)
{
    const char *key;
    PyObject *py_value;
    PyObject *py_failobj = Py_None;

    if (!PyArg_ParseTuple(args, "s|O:get", &key, &py_failobj))
	return NULL;
    py_value = Trie_get(mp->trie, key);
    if(!py_value)
	py_value = py_failobj;
    Py_INCREF(py_value);
    return py_value;
}

static char get_approximate__doc__[] =
"D.get_approximate(key, k) -> List of (key, value, mismatches) in D, allowing up to k mismatches in key.";

static void
_trie_get_approximate_helper(const char *key, const void *value,
			     const int mismatches, void *data)
{
    /* Append a tuple of (key, value) to data, which is a PyList. */
    PyObject *py_list = (PyObject *)data,
	*py_value = (PyObject *)value,
	*py_key,
	*py_tuple,
	*py_mismatches;

    if(PyErr_Occurred())
	return;

#ifdef IS_PY3K
    if(!(py_key = PyUnicode_FromFormat(key)))
#else
    if(!(py_key = PyString_FromString(key)))
#endif
	return;
#ifdef IS_PY3K
    if(!(py_mismatches = PyLong_FromLong(mismatches))) {
#else
    if(!(py_mismatches = PyInt_FromLong(mismatches))) {
#endif
	Py_DECREF(py_key);
	return;
    }
    Py_INCREF(py_value);

    if(!(py_tuple = PyTuple_New(3))) {
	Py_DECREF(py_key);
	Py_DECREF(py_mismatches);
	Py_DECREF(py_value);
	return;
    }
    PyTuple_SetItem(py_tuple, 0, py_key);
    PyTuple_SetItem(py_tuple, 1, py_value);
    PyTuple_SetItem(py_tuple, 2, py_mismatches);
    PyList_Append(py_list, py_tuple);
    Py_DECREF(py_tuple);
}

static PyObject *
trie_get_approximate(trieobject *mp, PyObject *args)
{
    const char *key;
    int k;
    PyObject *py_list;

    if (!PyArg_ParseTuple(args, "si:get_approximate", &key, &k))
	return NULL;

    if(!(py_list = PyList_New(0)))
	return NULL;
    Trie_get_approximate(mp->trie, key, k,
			 _trie_get_approximate_helper, (void *)py_list);
    if(PyErr_Occurred()) {
	Py_DECREF(py_list);
	return NULL;
    }
    return py_list;
}

static long
trie_nohash(PyObject *self)
{
    PyErr_SetString(PyExc_TypeError, "trie objects are unhashable");
    return -1;
}

static PyMappingMethods trie_as_mapping = {
    (lenfunc)trie_length,        /*mp_length*/
    (binaryfunc)trie_subscript,  /*mp_subscript*/
    (objobjargproc)trie_ass_sub  /*mp_ass_subscript*/
};

static PySequenceMethods trie_as_sequence = {
    (lenfunc)trie_length,        /* sq_length */
    NULL,                        /* sq_concat */
    NULL,                        /* sq_repeat */
    NULL,                        /* sq_item */
    NULL,                        /* sq_slice */
    NULL,                        /* sq_ass_item */
    NULL,                        /* sq_ass_slice */
    (objobjproc)trie_contains,   /* sq_contains */
    NULL,                        /* sq_inplace_concat */
    NULL                         /* sq_inplace_repeat */
};

static PyMethodDef trieobj_methods[] = {
    {"has_key", (PyCFunction)trie_has_key_onearg,  METH_VARARGS,
     has_key__doc__},
    {"has_prefix", (PyCFunction)trie_has_prefix_onearg,  METH_VARARGS,
     has_prefix__doc__},
    {"with_prefix", (PyCFunction)trie_with_prefix_onearg,  METH_VARARGS,
     with_prefix__doc__},
    {"keys",    (PyCFunction)trie_keys_noargs,     METH_VARARGS,
     keys__doc__},
    {"values",  (PyCFunction)trie_values_noargs,   METH_VARARGS,
     values__doc__},

    {"get",     (PyCFunction)trie_get,      METH_VARARGS,
     get__doc__},
    {"get_approximate",  (PyCFunction)trie_get_approximate,  METH_VARARGS,
     get_approximate__doc__},
    {NULL, NULL}   /* sentinel */
};

static PyTypeObject Trie_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "trie",
    sizeof(trieobject),
    0,
    trie_dealloc,       /*tp_dealloc*/
    0,                  /*tp_print*/
    0,                  /*tp_getattr*/
    0,                  /*tp_setattr*/
    0,                  /*tp_compare*/
    0,                  /*tp_repr*/
    0,                  /*tp_as_number*/
    &trie_as_sequence,  /*tp_as_sequence*/
    &trie_as_mapping,   /*tp_as_mapping*/
    trie_nohash,        /*tp_hash */
    0,                  /*tp_call */
    0,                  /* tp_str */
    0,                  /* tp_getattro */
    0,                  /* tp_setattro */
    0,                  /* tp_as_buffer */
#ifdef IS_PY3K
    0,                  /* tp_flags */
#else
    Py_TPFLAGS_HAVE_SEQUENCE_IN,    /* tp_flags */
#endif
    0,                  /* tp_doc */
    0,                  /* tp_traverse */
    0,                  /* tp_clear */
    0,                  /* tp_richcompare */
    0,                  /* tp_weaklistoffset */
    0,                  /* tp_iter */
    0,                  /* tp_iternext */
    trieobj_methods,    /* tp_methods */
};

static int
_write_to_handle(const void *towrite, const int length, void *handle)
{
    PyObject *py_handle = (PyObject *)handle,
	*py_retval = NULL;
    int success = 0;

    if(!length)
	return 1;

    if(!(py_retval = PyObject_CallMethod(py_handle, "write", "s#",
					 towrite, length)))
	goto _write_to_handle_cleanup;
    success = 1;

 _write_to_handle_cleanup:
    if(py_retval) {
	Py_DECREF(py_retval);
    }
    return success;
}

static int _write_value_to_handle(const void *value, void *handle)
{
    PyObject *py_value = (PyObject *)value,
	*py_marshalled = NULL,
        *bytes = NULL;
    char *marshalled;
    Py_ssize_t length;
    int success = 0;

#ifdef Py_MARSHAL_VERSION
    if(!(py_marshalled =
	 PyMarshal_WriteObjectToString(py_value, Py_MARSHAL_VERSION)))
        goto _write_value_to_handle_cleanup;
#else
    if(!(py_marshalled = PyMarshal_WriteObjectToString(py_value)))
        goto _write_value_to_handle_cleanup;
#endif

#ifdef IS_PY3K
    if (!PyBytes_Check(py_marshalled)) {
        PyErr_SetString(PyExc_TypeError, "marshalled data expected to be bytes");
	goto _write_value_to_handle_cleanup;
    }
    if(PyBytes_AsStringAndSize(py_marshalled, &marshalled, &length) == -1)
#else
    if(PyString_AsStringAndSize(py_marshalled, &marshalled, &length) == -1)
#endif
	goto _write_value_to_handle_cleanup;
    if(!_write_to_handle(&length, sizeof(length), handle))
	goto _write_value_to_handle_cleanup;
    if (length != (int)length)
	goto _write_value_to_handle_cleanup;
    if(!_write_to_handle(marshalled, (int)length, handle))
	goto _write_value_to_handle_cleanup;
    success = 1;

 _write_value_to_handle_cleanup:
    Py_XDECREF(py_marshalled);
    Py_XDECREF(bytes);

    return success;
}

static PyObject *
trie_save(PyObject *self, PyObject *args)
{
    PyObject *py_handle,
	*py_trie;
    trieobject *mp;

    if(!PyArg_ParseTuple(args, "OO:save", &py_handle, &py_trie))
        return NULL;
    mp = (trieobject *)py_trie;
    if(!Trie_serialize(mp->trie, _write_to_handle, _write_value_to_handle,
		       (void *)py_handle)) {
	if(!PyErr_Occurred())
	    PyErr_SetString(PyExc_RuntimeError,
			    "saving failed for some reason");
	return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static int
_read_from_handle(void *wasread, const int length, void *handle)
{
    PyObject *py_handle = (PyObject *)handle;
    PyObject *py_retval = NULL;
#ifdef IS_PY3K
    PyObject* bytes = NULL;
#endif
    int success = 0;
    char* buffer;

    if(!length) {
        PyErr_SetString(PyExc_RuntimeError, "data length is zero");
	return 0;
    }

    py_retval = PyObject_CallMethod(py_handle, "read", "i", length);
#ifdef IS_PY3K
    if(!PyUnicode_Check(py_retval)) {
#else
    if(!PyString_Check(py_retval)) {
#endif
        PyErr_SetString(PyExc_TypeError, "expected a string");
        goto error;
    }
#ifdef IS_PY3K
    bytes = PyUnicode_AsASCIIString(py_retval);
    if(!bytes) {
	PyErr_SetString(PyExc_TypeError, "expected an ASCII string");
	goto error;
    }
    buffer = PyBytes_AsString(bytes);
#else
    buffer = PyString_AS_STRING(py_retval);
#endif
    memcpy(wasread, buffer, length);
    success = 1;
error:
#ifdef IS_PY3K
    Py_XDECREF(bytes);
#endif
    Py_XDECREF(py_retval);
    return success;
}

static void *
_read_value_from_handle(void *handle)
{
    Py_ssize_t length;
    char* KEY;
    PyObject* VALUE;

    if(!_read_from_handle(&length, sizeof(length), handle))
	return NULL;
    if(length < 0)
    {
        return NULL;
    }
    KEY = malloc(length);
    if(length < 0)
    {
        PyErr_SetString(PyExc_MemoryError, "insufficient memory to read value");
        return NULL;
    }
    VALUE = NULL;
    if(_read_from_handle(KEY, length, handle))
        VALUE = PyMarshal_ReadObjectFromString(KEY, length);
    free(KEY);
    return VALUE;
}


static PyObject *
trie_load(PyObject *self, PyObject *args)
{
    PyObject *py_handle;
    Trie* trie;
    trieobject *trieobj;

    if(!PyArg_ParseTuple(args, "O:load", &py_handle))
	return NULL;

    if(!(trie = Trie_deserialize(_read_from_handle, _read_value_from_handle,
				 py_handle))) {
	if(!PyErr_Occurred())
	    PyErr_SetString(PyExc_RuntimeError,
			    "loading failed for some reason");
	return NULL;
    }

    if(!(trieobj = PyObject_New(trieobject, &Trie_Type))) {
	Trie_del(trie);
	return NULL;
    }
    trieobj->trie = trie;
    return (PyObject *)trieobj;
}

static PyMethodDef trie_methods[] = {
    {"trie", trie_trie, METH_VARARGS,
     "trie() -> new Trie object."},
    {"load", trie_load, METH_VARARGS,
     "load(handle) -> trie object"},
    {"save", trie_save, METH_VARARGS,
     "save(handle, trie), save a trie object to a handle"},
    {NULL, NULL, 0, NULL}
};

static char trie__doc__[] =
"\
This module implements a trie data structure.  This allows an O(M)\n\
lookup of a string in a dictionary, where M is the length of the\n\
string.  It also supports approximate matches.\n\
\n\
Functions:\n\
trie    Create a new trie object.\n\
save    Save a trie to a handle.\n\
load    Load a trie from a handle.\n\
\n\
";

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "trie",
        trie__doc__,
        -1,
        trie_methods,
        NULL,
        NULL,
        NULL,
        NULL
};

PyObject *
PyInit_trie(void)

#else

void inittrie(void)
#endif
{
    Py_TYPE(&Trie_Type) = &PyType_Type;

    if (PyType_Ready(&Trie_Type) < 0)
#if PY_MAJOR_VERSION >= 3
        return NULL;
#else
        return;
#endif

#if PY_MAJOR_VERSION >= 3
    return PyModule_Create(&moduledef);
#else
    Py_InitModule4("trie",
                   trie_methods,
                   trie__doc__,
                   NULL,
                   PYTHON_API_VERSION);
#endif
}
