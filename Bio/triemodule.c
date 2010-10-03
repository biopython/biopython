#include <Python.h>
#include <marshal.h>
#include "trie.h"

#if PY_VERSION_HEX < 0x02050000
#define Py_ssize_t int
#endif



staticforward PyTypeObject Trie_Type;

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

    /* Make sure key is a string. */
    if(!PyString_Check(py_key)) {
	PyErr_SetString(PyExc_TypeError, "key must be a string");
	return NULL;
    }
    key = PyString_AS_STRING(py_key);
    py_value = Trie_get(mp->trie, key);
    if(py_value == NULL)
	PyErr_SetString(PyExc_KeyError, key);
    else
	Py_INCREF(py_value);
    return py_value;
}

static int
trie_ass_sub(trieobject *mp, PyObject *py_key, PyObject *py_value)
{
    const char *key;
    PyObject *py_prev;

    /* Make sure key is a string. */
    if(!PyString_Check(py_key)) {
	PyErr_SetString(PyExc_TypeError, "key must be a string");
	return -1;
    }
    key = PyString_AS_STRING(py_key);
    
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
	if(!py_prev) {
	    PyErr_SetString(PyExc_KeyError, key);
	    return -1;
	}
	Trie_set(mp->trie, key, NULL);
    }
    /* The client wants to set a key in the dictionary. */
    else {
	Py_INCREF(py_value);
	if(Trie_set(mp->trie, key, py_value)) {
	    PyErr_SetString(PyExc_AssertionError, "error setting trie");
	    return -1;
	}
    }
    return 0;
}

static char has_key__doc__[] =
"D.has_key(k) -> 1 if D has a key k, else 0";

static PyObject *
trie_has_key(trieobject *mp, PyObject *py_key)
{
    const char *key;
    int has_key;

    /* Make sure key is a string. */
    if(!PyString_Check(py_key)) {
	PyErr_SetString(PyExc_TypeError, "key must be a string");
	return NULL;
    }
    key = PyString_AS_STRING(py_key);
    has_key = Trie_has_key(mp->trie, key);
    return PyInt_FromLong((long)has_key);
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

    /* Make sure prefix is a string. */
    if(!PyString_Check(py_prefix)) {
	PyErr_SetString(PyExc_TypeError, "k must be a string");
	return NULL;
    }
    prefix = PyString_AS_STRING(py_prefix);
    has_prefix = Trie_has_prefix(mp->trie, prefix);
    return PyInt_FromLong((long)has_prefix);
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

    if(!(py_key = PyString_FromString(key)))
	return;
    PyList_Append(py_list, py_key);
    Py_DECREF(py_key);
}

static PyObject *
trie_with_prefix(trieobject *mp, PyObject *py_prefix)
{
    const char *prefix;
    PyObject *py_list;

    /* Make sure prefix is a string. */
    if(!PyString_Check(py_prefix)) {
	PyErr_SetString(PyExc_TypeError, "k must be a string");
	return NULL;
    }
    prefix = PyString_AS_STRING(py_prefix);

    if(!(py_list = PyList_New(0)))
	return NULL;
    Trie_with_prefix(mp->trie, prefix, 
		     _trie_with_prefix_helper, (void *)py_list);
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

    if(!(py_key = PyString_FromString(key)))
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

    if(!(py_key = PyString_FromString(key)))
	return;
    if(!(py_mismatches = PyInt_FromLong(mismatches))) {
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
/* The first member of PyMappingMethods was redefined in Python 2.5. */
#if PY_VERSION_HEX < 0x02050000
    (inquiry)trie_length,        /*mp_length*/
#else
    (lenfunc)trie_length,        /*mp_length*/
#endif
    (binaryfunc)trie_subscript,  /*mp_subscript*/
    (objobjargproc)trie_ass_sub  /*mp_ass_subscript*/
};

static PyMethodDef trieobj_methods[] = {
    /*  METH_O and METH_NOARGS require Python 2.2.
    {"has_key", (PyCFunction)trie_has_key,  METH_O,
     has_key__doc__},
    {"has_prefix", (PyCFunction)trie_has_prefix,  METH_O,
     has_prefix__doc__},
    {"with_prefix", (PyCFunction)trie_with_prefix,  METH_O,
     with_prefix__doc__},
    {"keys",    (PyCFunction)trie_keys,     METH_NOARGS,
     keys__doc__},
    {"values",  (PyCFunction)trie_values,   METH_NOARGS,
     values__doc__},
    */

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

static PyObject *trie_getattr(PyObject *obj, const char *name)
{
    return Py_FindMethod(trieobj_methods, obj, name);

}

static PyTypeObject Trie_Type = {
    PyObject_HEAD_INIT(NULL)
    0,
    "trie",
    sizeof(trieobject),
    0,
    trie_dealloc,       /*tp_dealloc*/
    0,                  /*tp_print*/
    (getattrfunc)trie_getattr,                  /*tp_getattr*/
    0,                  /*tp_setattr*/
    0,                  /*tp_compare*/
    0,                  /*tp_repr*/
    0,                  /*tp_as_number*/
    0,                  /*tp_as_sequence*/
    &trie_as_mapping,   /*tp_as_mapping*/
    trie_nohash,        /*tp_hash */
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
	*py_marshalled = NULL;
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
    if(PyString_AsStringAndSize(py_marshalled, &marshalled, &length) == -1)
	goto _write_value_to_handle_cleanup;
    if(!_write_to_handle(&length, sizeof(length), handle))
	goto _write_value_to_handle_cleanup;
    if (length != (int)length)
	goto _write_value_to_handle_cleanup;
    if(!_write_to_handle(marshalled, (int)length, handle))
	goto _write_value_to_handle_cleanup;
    success = 1;

 _write_value_to_handle_cleanup:
    if(py_marshalled) {
	Py_DECREF(py_marshalled);
    }

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
    PyObject *py_handle = (PyObject *)handle,
	*py_retval = NULL;
    void *retval;
    int success = 0;
    PyBufferProcs *buffer;
    int segment;
    int bytes_read, bytes_left;
    
    if(!length)
	return 1;

    if(!(py_retval = PyObject_CallMethod(py_handle, "read", "i", length)))
	goto _read_from_handle_cleanup;
    if(!py_retval->ob_type->tp_as_buffer) {
	PyErr_SetString(PyExc_ValueError, "read method should return buffer");
	goto _read_from_handle_cleanup;
    }
    if(!(py_retval->ob_type->tp_flags & Py_TPFLAGS_DEFAULT)) {
	PyErr_SetString(PyExc_ValueError, "no bf_getcharbuffer slot");
	goto _read_from_handle_cleanup;
    }
    buffer = py_retval->ob_type->tp_as_buffer;
    if(!buffer->bf_getreadbuffer) {
	PyErr_SetString(PyExc_ValueError, "no bf_getreadbuffer");
	goto _read_from_handle_cleanup;
    }

    bytes_left = length;
    segment = 0;
    while(bytes_left > 0) {
	if((bytes_read = buffer->bf_getreadbuffer(py_retval, 
						  segment, &retval)) == -1)
	    goto _read_from_handle_cleanup; 
	memcpy(wasread, retval, bytes_read);
	wasread = (void *)((char *)wasread + bytes_read);
	bytes_left -= bytes_read;
	segment += 1;
    }

    success = 1;
    
 _read_from_handle_cleanup:
    if(py_retval) {
	Py_DECREF(py_retval);
    }
    return success;
}

#define MAX_KEY_LENGTH 2000
static void *
_read_value_from_handle(void *handle)
{
    Py_ssize_t length;
    char KEY[MAX_KEY_LENGTH];

    if(!_read_from_handle((void *)&length, sizeof(length), (void *)handle))
	return NULL;
    if(length < 0 || length >= MAX_KEY_LENGTH)
	return NULL;
    if(!_read_from_handle((void *)KEY, length, (void *)handle))
	return NULL;
    return PyMarshal_ReadObjectFromString(KEY, length);
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

DL_EXPORT(void)
inittrie(void) 
{
    Trie_Type.ob_type = &PyType_Type;

    (void) Py_InitModule3("trie", trie_methods, trie__doc__);
}
