#define PY_SSIZE_T_CLEAN

#include "Python.h"

typedef struct {
    PyObject_HEAD
    char character;
    Py_ssize_t length;
} UndefinedSeqDataObject;

static void
UndefinedSeqDataObject_dealloc(UndefinedSeqDataObject *op)
{
    Py_TYPE(op)->tp_free(op);
}

static PyObject*
UndefinedSeqDataObject_str(UndefinedSeqDataObject* self)
{
    PyObject* text = PyUnicode_New(self->length, 127);
    memset(PyUnicode_DATA(text), self->character, self->length);
    return text;
}

static int
UndefinedSeqDataObject_bf_getbuffer(UndefinedSeqDataObject *self, Py_buffer *view, int flags)
{
   static Py_ssize_t stride = 0;
   if ((flags & PyBUF_STRIDES) != PyBUF_STRIDES) {
       /* Need to create a contiguous buffer */
       int result;
       PyObject* bytes = PyBytes_FromObject((PyObject*)self);
       if (!bytes) return -1;
       result = PyObject_GetBuffer(bytes, view, flags);
       Py_DECREF(bytes);
       return result;
   }
   view->obj = (PyObject*)self;
   Py_INCREF(self);
   view->buf = &self->character;
   view->len = self->length;
   view->readonly = 1;
   view->itemsize = 1;
   if (flags & PyBUF_FORMAT) view->format = "B";
   else view->format = NULL;
   view->ndim = 1;
   view->shape = NULL;
   if ((flags & PyBUF_ND) == PyBUF_ND)
       view->shape = &(view->len);
   view->strides = &stride;
   view->suboffsets = NULL;
   view->internal = NULL;
   return 0;
}


static PyBufferProcs UndefinedSeqDataObject_as_buffer = {
    (getbufferproc)UndefinedSeqDataObject_bf_getbuffer,
    NULL,
};

static PyTypeObject UndefinedSeqDataType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "UndefinedSeqDataType",                     /* tp_name */
    sizeof(UndefinedSeqDataObject),             /* tp_basicsize */
    0,                                          /* tp_itemsize */
    (destructor)UndefinedSeqDataObject_dealloc, /* tp_dealloc */
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
    (reprfunc)UndefinedSeqDataObject_str,       /* tp_str */
    0,                                          /* tp_getattro */
    0,                                          /* tp_setattro */
    &UndefinedSeqDataObject_as_buffer,          /* tp_as_buffer */
};

static PyTypeObject SeqType;

typedef struct {
    PyObject_HEAD
    PyObject* data;
    PyObject* id;
    PyObject* name;
    PyObject* description;
    PyObject* annotations;
    PyObject* features;
    PyObject* dbxrefs;
    PyObject* letter_annotations;
} SeqObject;

static SeqObject *
Seq_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    SeqObject *object;
    PyObject *data = Py_None;
    char* character = NULL;

    static char *kwlist[] = {"data", "character", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|Os", kwlist,
                                     &data, &character))
        return NULL;

    if (PyLong_Check(data)) {
        Py_ssize_t length = PyLong_AsSsize_t(data);
        if (length == -1 && PyErr_Occurred()) return NULL;
        if (length < 0) {
            PyErr_Format(PyExc_ValueError,
                         "expected sequence data or a positive integer "
                         "(received %zd)", length);
            return NULL;
        }
        if (character && strlen(character) != 1) {
            PyErr_SetString(PyExc_ValueError,
                "character should be a single letter");
            return NULL;
        }
        data = type->tp_alloc(&UndefinedSeqDataType, 0);
        if (!data) return NULL;
        ((UndefinedSeqDataObject*)data)->length = length;
        if (character) {
            ((UndefinedSeqDataObject*)data)->character = *character;
        }
        else
            ((UndefinedSeqDataObject*)data)->character = '?';
    }
    else {
        if (character != NULL) {
            PyErr_SetString(PyExc_ValueError,
                "character should be None if data is given");
            return NULL;
        }
        if (PyUnicode_Check(data)) {
            data = PyUnicode_AsASCIIString(data);
            if (!data) return NULL;
        }
        else if (PyBytes_Check(data)) {
            Py_INCREF(data);
        }
        else if (PyByteArray_Check(data)) {
            Py_INCREF(data);
        }
        else if (PyObject_IsInstance(data, (PyObject*)&SeqType)) {
            data = ((SeqObject*)data)->data;
            if (PyBytes_Check(data)) {
                Py_INCREF(data);
            }
            else if (PyByteArray_Check(data)) {
                data = PyByteArray_FromObject(data);
                if (!data) return NULL;
            }
            else {
                PyErr_Format(PyExc_TypeError,
                    "found sequence data of type %s "
                    "(expected None, bytes, or bytearray)",
                    Py_TYPE(data)->tp_name);
                return NULL;
            }
        }
        else if (data == Py_None) {
            Py_INCREF(data);
        }
        else {
            PyObject* bytes = PyBytes_FromObject(data);
            if (!bytes) {
                PyErr_Format(PyExc_TypeError,
                    "data should be a string, bytes, bytearray, Seq object, "
                    "or any other object that provides characters via the "
                    "buffer protocol, not %s",
                    Py_TYPE(data)->tp_name);
                return NULL;
            }
            data = bytes;
        }
    }

    object = (SeqObject*)type->tp_alloc(type, 0);
    if (!object) {
        Py_DECREF(data);
        return NULL;
    }

    object->data = data;
    object->id = NULL;
    object->description = NULL;
    object->annotations = NULL;
    object->features = NULL;
    object->dbxrefs = NULL;
    object->letter_annotations = NULL;

    return object;
}

static void
Seq_dealloc(SeqObject *op)
{
    Py_DECREF(op->data);
    Py_XDECREF(op->id);
    Py_XDECREF(op->description);
    Py_XDECREF(op->annotations);
    Py_XDECREF(op->features);
    Py_XDECREF(op->dbxrefs);
    Py_XDECREF(op->letter_annotations);
    Py_TYPE(op)->tp_free(op);
}

static PyObject*
Seq_str(SeqObject* self)
{
    PyObject* data = self->data;
    if (PyObject_IsInstance(data, (PyObject*)&UndefinedSeqDataType))
        return PyObject_Str(data);
    return PyUnicode_FromEncodedObject(data, NULL, NULL);
}

static PyObject *
Seq_getattro(SeqObject *self, PyObject *name)
{
    if (PyUnicode_CompareWithASCIIString(name, "id") == 0) {
        if (!self->id) return PyUnicode_New(0, 127);
        Py_INCREF(self->id);
        return self->id;
    }
    if (PyUnicode_CompareWithASCIIString(name, "name") == 0) {
        if (!self->name) return PyUnicode_New(0, 127);
        Py_INCREF(self->name);
        return self->name;
    }
    if (PyUnicode_CompareWithASCIIString(name, "description") == 0) {
        if (!self->description) {
            self->description = PyUnicode_New(0, 127);
            if (!self->description) return NULL;
        }
        Py_INCREF(self->description);
        return self->description;
    }
    if (PyUnicode_CompareWithASCIIString(name, "annotations") == 0) {
        if (!self->annotations) {
            self->annotations = PyDict_New();
            if (!self->annotations) return NULL;
        }
        Py_INCREF(self->annotations);
        return self->annotations;
    }
    if (PyUnicode_CompareWithASCIIString(name, "dbxrefs") == 0) {
        if (!self->dbxrefs) {
            self->dbxrefs = PyList_New(0);
            if (!self->dbxrefs) return NULL;
        }
        Py_INCREF(self->dbxrefs);
        return self->dbxrefs;
    }
    if (PyUnicode_CompareWithASCIIString(name, "features") == 0) {
        if (!self->features) {
            self->features = PyList_New(0);
            if (!self->features) return NULL;
        }
        Py_INCREF(self->features);
        return self->features;
    }
    if (PyUnicode_CompareWithASCIIString(name, "letter_annotations") == 0) {
        if (self->letter_annotations) {
            Py_INCREF(self->letter_annotations);
            return self->letter_annotations;
        }
        /* Let Bio.Seq.Seq create the restricted dictionary if needed. */
    }
    return PyObject_GenericGetAttr((PyObject*)self, name);
}

static int
Seq_setattro(SeqObject *self, PyObject *name, PyObject *value)
{
    if (PyUnicode_CompareWithASCIIString(name, "id") == 0) {
        if (value) {
            if (value != Py_None && !PyUnicode_Check(value)) {
                PyErr_Format(PyExc_TypeError,
                    "attribute id requires a string (received type %s)",
                    Py_TYPE(value)->tp_name);
                return -1;
            }
            Py_INCREF(value);
        }
        Py_XDECREF(self->id);
        self->id = value;
        return 0;
    }
    if (PyUnicode_CompareWithASCIIString(name, "name") == 0) {
        if (value) {
            if (value != Py_None && !PyUnicode_Check(value)) {
                PyErr_Format(PyExc_TypeError,
                    "attribute name requires a string (received type %s)",
                    Py_TYPE(value)->tp_name);
                return -1;
            }
            Py_INCREF(value);
        }
        Py_XDECREF(self->name);
        self->name = value;
        return 0;
    }
    if (PyUnicode_CompareWithASCIIString(name, "description") == 0) {
        if (value) {
            if (value != Py_None && !PyUnicode_Check(value)) {
                PyErr_Format(PyExc_TypeError,
                    "attribute description requires a string (received type %s)",
                    Py_TYPE(value)->tp_name);
                return -1;
            }
            Py_INCREF(value);
        }
        Py_XDECREF(self->description);
        self->description = value;
        return 0;
    }
    if (PyUnicode_CompareWithASCIIString(name, "annotations") == 0) {
        if (value) {
            if (!PyDict_Check(value)) {
                PyErr_Format(PyExc_TypeError,
                    "attribute annotations requires a dictionary (received type %s)",
                    Py_TYPE(value)->tp_name);
                return -1;
            }
            Py_INCREF(value);
        }
        Py_XDECREF(self->annotations);
        self->annotations = value;
        return 0;
    }
    if (PyUnicode_CompareWithASCIIString(name, "features") == 0) {
        if (value) {
            if (!PyList_Check(value)) {
                PyErr_Format(PyExc_TypeError,
                    "attribute features requires a list (received type %s)",
                    Py_TYPE(value)->tp_name);
                return -1;
            }
            Py_INCREF(value);
        }
        Py_XDECREF(self->features);
        self->features = value;
        return 0;
    }
    if (PyUnicode_CompareWithASCIIString(name, "dbxrefs") == 0) {
        if (value) {
            if (!PyList_Check(value)) {
                PyErr_Format(PyExc_TypeError,
                    "attribute dbxrefs requires a list (received type %s)",
                    Py_TYPE(value)->tp_name);
                return -1;
            }
            Py_INCREF(value);
        }
        Py_XDECREF(self->dbxrefs);
        self->dbxrefs = value;
        return 0;
    }
    if (PyUnicode_CompareWithASCIIString(name, "letter_annotations") == 0) {
        if (value) {
            if (!PyDict_Check(value)) {
                PyErr_Format(PyExc_TypeError,
                    "attribute letter_annotations requires a dictionary (received type %s)",
                    Py_TYPE(value)->tp_name);
                return -1;
            }
            Py_INCREF(value);
        }
        Py_XDECREF(self->letter_annotations);
        self->letter_annotations = value;
        return 0;
    }
    return PyObject_GenericSetAttr((PyObject*)self, name, value);
}

static Py_ssize_t
Seq_mapping_length(SeqObject* self)
{
    PyObject* data = self->data;
    if (PyObject_IsInstance(data, (PyObject*)&UndefinedSeqDataType)) {
        return ((UndefinedSeqDataObject*)data)->length;
    }
    return PyObject_Size(data);
}

static PyObject*
Seq_subscript(SeqObject *self, PyObject *key)
{
    PyObject* data = self->data;
    PyObject* result = NULL;
    Py_buffer view;
    const char* s;

    if (PyObject_GetBuffer(data, &view, PyBUF_STRIDES | PyBUF_FORMAT) < 0) {
        return PyObject_GetItem(data, key);
    }

    s = view.buf;

    if (PyIndex_Check(key)) {
        Py_ssize_t i = PyNumber_AsSsize_t(key, PyExc_IndexError);
        if (!(i == -1 && PyErr_Occurred())) {
            if (i < 0) i += view.len;
            if (i < 0 || i >= view.len)
                PyErr_SetString(PyExc_IndexError, "index out of range");
            else
                result = PyBytes_FromStringAndSize(s + i, 1);
        }
    }
    else if (PySlice_Check(key)) {
        Py_ssize_t start, stop, step, slicelength, cur, i;

        if (PySlice_Unpack(key, &start, &stop, &step) < 0) return NULL;

        slicelength = PySlice_AdjustIndices(view.len, &start, &stop, step);
        if (slicelength <= 0)
            result = PyBytes_FromStringAndSize("", 0);
        else if (step == 1)
            result = PyBytes_FromStringAndSize(s + start, slicelength);
        else {
            result = PyBytes_FromStringAndSize(NULL, slicelength);
            if (result) {
                char* target = PyBytes_AS_STRING(result);
                for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                    target[i] = s[cur];
                }
            }
        }
    }
    else {
        PyErr_Format(PyExc_TypeError,
                     "Seq indices must be integers or slices, not %.200s",
                     Py_TYPE(key)->tp_name);
    }

    PyBuffer_Release(&view);

    return result;
}

static int
Seq_ass_subscript(SeqObject *self, PyObject *key, PyObject *value)
{
    int result;
    PyObject* data = self->data;

    if (!PyByteArray_Check(data)) {
        PyErr_SetString(PyExc_ValueError, "sequence is immutable");
        return -1;
    }

    if (!value) return PyObject_DelItem(data, key);

    if (PyUnicode_Check(value)) {
        value = PyUnicode_AsASCIIString(value);
        if (!value) return -1;
    }
    else Py_INCREF(value);

    if (PyIndex_Check(key)) {
        Py_buffer view;
        char* s;
        if (PyObject_GetBuffer(value, &view, PyBUF_STRIDES | PyBUF_FORMAT) < 0) {
            Py_DECREF(value);
            return -1;
        }
        Py_DECREF(value);
        if (view.ndim != 1 || view.len != 1 || strcmp(view.format, "B") != 0) {
            PyBuffer_Release(&view);
            PyErr_SetString(PyExc_RuntimeError, "expected a single byte");
            return -1;
        }
        s = view.buf;
        value = PyLong_FromLong((long)(*s));
        PyBuffer_Release(&view);
        if (!value) return -1;
    }

    result = PyObject_SetItem(data, key, value);
    Py_DECREF(value);
    return result;
}

static PyMappingMethods Seq_as_mapping = {
    (lenfunc)Seq_mapping_length,
    (binaryfunc)Seq_subscript,
    (objobjargproc)Seq_ass_subscript,
};

static int
Seq_bf_getbuffer(SeqObject *self, Py_buffer *view, int flags)
{
    if (self->data == Py_None) { /* FIXME */
        PyObject* text = PyObject_Str((PyObject*)self);
        if (!text) return -1;
        if (PyUnicode_READY(text) == -1) {
            Py_DECREF(text);
            return -1;
        }
        if (PyUnicode_KIND(text) != PyUnicode_1BYTE_KIND) {
            PyErr_SetString(PyExc_ValueError,
                "only ASCII strings can be used in comparisons");
            Py_DECREF(text);
            return -1;
        }
        view->obj = text;
        view->buf = PyUnicode_1BYTE_DATA(text);
        view->ndim = 1;
        view->format = "B";
        view->len = PyUnicode_GET_LENGTH(text);
        view->strides = NULL;
        view->shape = NULL;
        view->suboffsets = NULL;
        return 0;
    }
    return PyObject_GetBuffer(self->data, view, flags);
}

static PyBufferProcs Seq_as_buffer = {
    (getbufferproc)Seq_bf_getbuffer,
    NULL,
};

static PyObject* compare_buffers(Py_buffer *view1, Py_buffer *view2, int op)
{
    int comparison;
    const char* s1 = view1->buf;
    const char* s2 = view2->buf;
    const Py_ssize_t n1 = view1->len;
    const Py_ssize_t n2 = view2->len;
    Py_ssize_t n;
    if (view1->strides && view1->strides[0] == 0 && view2->strides && view2->strides[0] == 0) {
        switch (op) {
        case Py_EQ:
            if (n1 != n2) Py_RETURN_FALSE;
            if (*s1 != *s2) Py_RETURN_FALSE;
            Py_RETURN_TRUE;
        case Py_NE:
            if (n1 != n2) Py_RETURN_TRUE;
            if (*s1 != *s2) Py_RETURN_TRUE;
            Py_RETURN_FALSE;
        case Py_LT:
            if (n2 == 0) Py_RETURN_FALSE;
            if (n1 == 0) Py_RETURN_TRUE;
            if (*s1 < *s2) Py_RETURN_TRUE;
            if (*s1 > *s2) Py_RETURN_FALSE;
            if (n1 < n2) Py_RETURN_TRUE;
            Py_RETURN_FALSE;
        case Py_LE:
            if (n1 == 0) Py_RETURN_TRUE;
            if (n2 == 0) Py_RETURN_FALSE;
            if (*s1 < *s2) Py_RETURN_TRUE;
            if (*s1 > *s2) Py_RETURN_FALSE;
            if (n1 <= n2) Py_RETURN_TRUE;
            Py_RETURN_FALSE;
        case Py_GT:
            if (n1 == 0) Py_RETURN_FALSE;
            if (n2 == 0) Py_RETURN_TRUE;
            if (*s1 < *s2) Py_RETURN_FALSE;
            if (*s1 > *s2) Py_RETURN_TRUE;
            if (n1 > n2) Py_RETURN_TRUE;
            Py_RETURN_FALSE;
        case Py_GE:
            if (n2 == 0) Py_RETURN_TRUE;
            if (n1 == 0) Py_RETURN_FALSE;
            if (*s1 > *s2) Py_RETURN_TRUE;
            if (*s1 < *s2) Py_RETURN_FALSE;
            if (n1 >= n2) Py_RETURN_TRUE;
            Py_RETURN_FALSE;
        }
    }
    n = (n1 < n2) ? n1 : n2;
    if (view1->strides && view1->strides[0] == 0) {
        Py_ssize_t i;
        const char c1 = *s1;
        switch (op) {
        case Py_EQ:
            if (n1 != n2) Py_RETURN_FALSE;
            for (i = 0; i < n; i++) if (c1 != s2[i]) Py_RETURN_FALSE;
            Py_RETURN_TRUE;
        case Py_NE:
            if (n1 != n2) Py_RETURN_TRUE;
            for (i = 0; i < n; i++) if (c1 != s2[i]) Py_RETURN_TRUE;
            Py_RETURN_FALSE;
        case Py_LT:
            if (n2 == 0) Py_RETURN_FALSE;
            if (n1 == 0) Py_RETURN_TRUE;
            for (i = 0; i < n; i++) {
                if (c1 < s2[i]) Py_RETURN_TRUE;
                if (c1 > s2[i]) Py_RETURN_FALSE;
            }
            if (n1 < n2) Py_RETURN_TRUE;
            Py_RETURN_FALSE;
        case Py_LE:
            if (n1 == 0) Py_RETURN_TRUE;
            if (n2 == 0) Py_RETURN_FALSE;
            for (i = 0; i < n; i++) {
                if (c1 < s2[i]) Py_RETURN_TRUE;
                if (c1 > s2[i]) Py_RETURN_FALSE;
            }
            if (n1 <= n2) Py_RETURN_TRUE;
            Py_RETURN_FALSE;
        case Py_GT:
            if (n1 == 0) Py_RETURN_FALSE;
            if (n2 == 0) Py_RETURN_TRUE;
            for (i = 0; i < n; i++) {
                if (c1 < s2[i]) Py_RETURN_FALSE;
                if (c1 > s2[i]) Py_RETURN_TRUE;
            }
            if (n1 <= n2) Py_RETURN_FALSE;
            Py_RETURN_TRUE;
        case Py_GE:
            if (n2 == 0) Py_RETURN_TRUE;
            if (n1 == 0) Py_RETURN_FALSE;
            for (i = 0; i < n; i++) {
                if (c1 < s2[i]) Py_RETURN_FALSE;
                if (c1 > s2[i]) Py_RETURN_TRUE;
            }
            if (n1 < n2) Py_RETURN_FALSE;
            Py_RETURN_TRUE;
        }
    }
    if (view2->strides && view2->strides[0] == 0) {
        Py_ssize_t i;
        const char c2 = *s2;
        switch (op) {
        case Py_EQ:
            if (n1 != n2) Py_RETURN_FALSE;
            for (i = 0; i < n; i++) if (s1[i] != c2) Py_RETURN_FALSE;
            Py_RETURN_TRUE;
        case Py_NE:
            if (n1 != n2) Py_RETURN_TRUE;
            for (i = 0; i < n; i++) if (s1[i] != c2) Py_RETURN_TRUE;
            Py_RETURN_FALSE;
        case Py_LT:
            if (n2 == 0) Py_RETURN_FALSE;
            if (n1 == 0) Py_RETURN_TRUE;
            for (i = 0; i < n; i++) {
                if (s1[i] < c2) Py_RETURN_TRUE;
                if (s1[i] > c2) Py_RETURN_FALSE;
            }
            if (n1 < n2) Py_RETURN_TRUE;
            Py_RETURN_FALSE;
        case Py_LE:
            if (n1 == 0) Py_RETURN_TRUE;
            if (n2 == 0) Py_RETURN_FALSE;
            for (i = 0; i < n; i++) {
                if (s1[i] < c2) Py_RETURN_TRUE;
                if (s1[i] > c2) Py_RETURN_FALSE;
            }
            if (n1 <= n2) Py_RETURN_TRUE;
            Py_RETURN_FALSE;
        case Py_GT:
            if (n1 == 0) Py_RETURN_FALSE;
            if (n2 == 0) Py_RETURN_TRUE;
            for (i = 0; i < n; i++) {
                if (s1[i] < c2) Py_RETURN_FALSE;
                if (s1[i] > c2) Py_RETURN_TRUE;
            }
            if (n1 <= n2) Py_RETURN_FALSE;
            Py_RETURN_TRUE;
        case Py_GE:
            if (n2 == 0) Py_RETURN_TRUE;
            if (n1 == 0) Py_RETURN_FALSE;
            for (i = 0; i < n; i++) {
                if (s1[i] < c2) Py_RETURN_FALSE;
                if (s1[i] > c2) Py_RETURN_TRUE;
            }
            if (n1 < n2) Py_RETURN_FALSE;
            Py_RETURN_TRUE;
        }
    }
    if ( (view1->strides && view1->strides[0] != 1)
      || (view2->strides && view2->strides[0] != 1) ) {
        PyErr_SetString(PyExc_ValueError, "data are not contiguous");
        return NULL;
    }
    comparison = memcmp(s1, s2, n);
    if (comparison == 0) {
        if (n1 < n2) comparison = -1;
        else if (n1 > n2) comparison = +1;
    }
    switch (op) {
    case Py_EQ:
        if (comparison == 0) Py_RETURN_TRUE;
        Py_RETURN_FALSE;
    case Py_NE:
        if (comparison == 0) Py_RETURN_FALSE;
        Py_RETURN_TRUE;
    case Py_LE:
        if (comparison <= 0) Py_RETURN_TRUE;
        Py_RETURN_FALSE;
    case Py_GE:
        if (comparison >= 0) Py_RETURN_TRUE;
        Py_RETURN_FALSE;
    case Py_LT:
        if (comparison < 0) Py_RETURN_TRUE;
        Py_RETURN_FALSE;
    case Py_GT:
        if (comparison > 0) Py_RETURN_TRUE;
        Py_RETURN_FALSE;
    default:
        PyErr_BadArgument();
        return NULL;
    }
}

static PyObject*
Seq_richcompare(SeqObject *self, PyObject *other, int op)
{
    Py_buffer view1;
    Py_buffer view2;
    PyObject *result;
    if ((PyObject*)self == other) {
        switch (op) {
        case Py_EQ:
        case Py_LE:
        case Py_GE:
            Py_RETURN_TRUE;
        case Py_NE:
        case Py_LT:
        case Py_GT:
            Py_RETURN_FALSE;
        default:
            PyErr_BadArgument();
            return NULL;
        }
    }

    if (PyObject_GetBuffer((PyObject*)self, &view1, PyBUF_STRIDES | PyBUF_FORMAT) < 0)
        return NULL;
    if (PyUnicode_Check(other)) {
        if (PyUnicode_KIND(other) != PyUnicode_1BYTE_KIND) {
            PyErr_SetString(PyExc_ValueError,
                "only ASCII strings can be used in comparisons");
            PyBuffer_Release(&view1);
            return NULL;
        }
        view2.obj = NULL;
        view2.buf = PyUnicode_1BYTE_DATA(other);
        view2.ndim = 1;
        view2.format = "B";
        view2.len = PyUnicode_GET_LENGTH(other);
        view2.strides = NULL;
        view2.suboffsets = NULL;
        view2.internal = NULL;
    }
    else if (PyObject_GetBuffer(other, &view2, PyBUF_STRIDES | PyBUF_FORMAT) < 0) {
        PyBuffer_Release(&view1);
        return NULL;
    }
    if (view1.ndim != 1
     || view2.ndim != 1
     || strcmp(view1.format, "B") != 0
     || strcmp(view2.format, "B") != 0) {
        PyBuffer_Release(&view1);
        PyBuffer_Release(&view2);
        PyErr_Format(PyExc_TypeError,
                     "comparison not supported between instances of %s and %s",
                     Py_TYPE(self)->tp_name, Py_TYPE(other)->tp_name);
        return NULL;
    }

    result = compare_buffers(&view1, &view2, op);
    PyBuffer_Release(&view1);
    PyBuffer_Release(&view2);

    return result;
}

PyDoc_STRVAR(Seq_reduce_doc,
"__reduce__($self, /)\n"
"--\n"
"\n"
"Return state information for pickling.");

static PyObject *
Seq_reduce(SeqObject *self, PyObject *Py_UNUSED(ignored))
{
    PyObject* id = self->id ? self->id : Py_None;
    PyObject* name = self->name ? self->name : Py_None;
    PyObject* description = self->description ? self->description : Py_None;
    PyObject* annotations = self->annotations ? self->annotations : Py_None;
    PyObject* features = self->features ? self->features : Py_None;
    PyObject* dbxrefs = self->dbxrefs ? self->dbxrefs : Py_None;
    PyObject* letter_annotations = self->letter_annotations ? self->letter_annotations : Py_None;

    switch (PyObject_IsInstance(self->data, (PyObject*)&UndefinedSeqDataType)) {
        case -1: return NULL;
        case 1: {
            UndefinedSeqDataObject *data = (UndefinedSeqDataObject*)self->data;
            return Py_BuildValue("O(nsOOOOOOO)",
                                 Py_TYPE(self), data->length, &data->character,
                                 id, name, description, annotations, features,
                                 dbxrefs, letter_annotations);
        }
        case 0: break;
    }

    return Py_BuildValue("O(OOOOOOOO)",
                         Py_TYPE(self), self->data, id, name, description,
                         annotations, features, dbxrefs, letter_annotations);
}

PyDoc_STRVAR(Seq_reverse_doc,
"Modify the mutable sequence to reverse itself.\n"
"\n"
"No return value.\n");

static PyObject *
Seq_reverse(SeqObject *self, PyObject *Py_UNUSED(ignored))
{
    Py_ssize_t i, n, m;
    char c;
    char* s;
    PyObject* data = self->data;

    if (!PyByteArray_Check(data)) {
        PyErr_SetString(PyExc_ValueError, "sequence is immutable");
        return NULL;
    }

    s = PyByteArray_AS_STRING(data);
    n = PyByteArray_GET_SIZE(data);
    m = n / 2;

    for (i = 0; i < m; i++) {
        c = s[i];
        s[i] = s[n-i-1];
        s[n-i-1] = c;
    }

    Py_RETURN_NONE;
}

PyDoc_STRVAR(Seq_complement_doc,
"Modify the mutable sequence into its RNA complement.\n"
"\n"
"No return value.\n");


static PyObject *
Seq_complement(SeqObject *self, PyObject *Py_UNUSED(ignored))
{
    Py_ssize_t i, n;
    char* s;
    PyObject* data = self->data;

    if (PyObject_IsInstance(data, (PyObject*)&UndefinedSeqDataType)) {
        Py_INCREF(self);
        return (PyObject*)self;
    }

    if (PyBytes_Check(data)) {
        /* make sure we create a new bytes object;
         * note that the bytes object may have been interned */
        n = PyBytes_GET_SIZE(data);
        s = PyBytes_AS_STRING(data);
        data = PyBytes_FromStringAndSize(NULL, n);
        if (!data) return NULL;
        memcpy(PyBytes_AS_STRING(data), s, n);
        s = PyBytes_AS_STRING(data);
    }
    else if (PyByteArray_Check(data)) {
        n = PyByteArray_GET_SIZE(data);
        s = PyByteArray_AS_STRING(data);
    }
    else {
        data = PyBytes_FromObject(data);
        if (!data) return NULL;
        n = PyBytes_GET_SIZE(data);
        s = PyBytes_AS_STRING(data);
    }

    for (i = 0; i < n; i++) {
        switch (s[i]) {
            case 'A': s[i] = 'T'; break;
            case 'B': s[i] = 'V'; break;
            case 'C': s[i] = 'G'; break;
            case 'D': s[i] = 'H'; break;
            case 'G': s[i] = 'C'; break;
            case 'H': s[i] = 'D'; break;
            case 'K': s[i] = 'M'; break;
            case 'M': s[i] = 'K'; break;
            case 'N': s[i] = 'N'; break;
            case 'R': s[i] = 'Y'; break;
            case 'S': s[i] = 'S'; break;
            case 'T': s[i] = 'A'; break;
            case 'U': s[i] = 'A'; break;
            case 'V': s[i] = 'B'; break;
            case 'W': s[i] = 'W'; break;
            case 'X': s[i] = 'X'; break;
            case 'Y': s[i] = 'R'; break;
            case 'a': s[i] = 't'; break;
            case 'b': s[i] = 'v'; break;
            case 'c': s[i] = 'g'; break;
            case 'd': s[i] = 'h'; break;
            case 'g': s[i] = 'c'; break;
            case 'h': s[i] = 'd'; break;
            case 'k': s[i] = 'm'; break;
            case 'm': s[i] = 'k'; break;
            case 'n': s[i] = 'n'; break;
            case 'r': s[i] = 'y'; break;
            case 's': s[i] = 's'; break;
            case 't': s[i] = 'a'; break;
            case 'u': s[i] = 'a'; break;
            case 'v': s[i] = 'b'; break;
            case 'w': s[i] = 'w'; break;
            case 'x': s[i] = 'x'; break;
            case 'y': s[i] = 'r'; break;
            default: break;
        }
    }

    if (PyByteArray_Check(data)) Py_RETURN_NONE;

    return data;
}

PyDoc_STRVAR(Seq_rna_complement_doc,
"Modify the mutable sequence into its RNA complement.\n"
"\n"
"No return value.\n");

static PyObject *
Seq_rna_complement(SeqObject *self, PyObject *Py_UNUSED(ignored))
{
    Py_ssize_t i, n;
    char* s;
    PyObject* data = self->data;

    if (PyObject_IsInstance(data, (PyObject*)&UndefinedSeqDataType)) {
        Py_INCREF(self);
        return (PyObject*)self;
    }

    if (PyBytes_Check(data)) {
        /* make sure we create a new bytes object;
         * note that the bytes object may have been interned */
        n = PyBytes_GET_SIZE(data);
        s = PyBytes_AS_STRING(data);
        data = PyBytes_FromStringAndSize(NULL, n);
        if (!data) return NULL;
        memcpy(PyBytes_AS_STRING(data), s, n);
        s = PyBytes_AS_STRING(data);
    }
    else if (PyByteArray_Check(data)) {
        n = PyByteArray_GET_SIZE(data);
        s = PyByteArray_AS_STRING(data);
    }
    else {
        data = PyBytes_FromObject(data);
        if (!data) return NULL;
        n = PyBytes_GET_SIZE(data);
        s = PyBytes_AS_STRING(data);
    }

    for (i = 0; i < n; i++) {
        switch (s[i]) {
            case 'A': s[i] = 'U'; break;
            case 'B': s[i] = 'V'; break;
            case 'C': s[i] = 'G'; break;
            case 'D': s[i] = 'H'; break;
            case 'G': s[i] = 'C'; break;
            case 'H': s[i] = 'D'; break;
            case 'K': s[i] = 'M'; break;
            case 'M': s[i] = 'K'; break;
            case 'N': s[i] = 'N'; break;
            case 'R': s[i] = 'Y'; break;
            case 'S': s[i] = 'S'; break;
            case 'T': s[i] = 'A'; break;
            case 'U': s[i] = 'A'; break;
            case 'V': s[i] = 'B'; break;
            case 'W': s[i] = 'W'; break;
            case 'X': s[i] = 'X'; break;
            case 'Y': s[i] = 'R'; break;
            case 'a': s[i] = 'u'; break;
            case 'b': s[i] = 'v'; break;
            case 'c': s[i] = 'g'; break;
            case 'd': s[i] = 'h'; break;
            case 'g': s[i] = 'c'; break;
            case 'h': s[i] = 'd'; break;
            case 'k': s[i] = 'm'; break;
            case 'm': s[i] = 'k'; break;
            case 'n': s[i] = 'n'; break;
            case 'r': s[i] = 'y'; break;
            case 's': s[i] = 's'; break;
            case 't': s[i] = 'a'; break;
            case 'u': s[i] = 'a'; break;
            case 'v': s[i] = 'b'; break;
            case 'w': s[i] = 'w'; break;
            case 'x': s[i] = 'x'; break;
            case 'y': s[i] = 'r'; break;
            default: break;
        }
    }

    if (PyByteArray_Check(data)) Py_RETURN_NONE;

    return data;
}

PyDoc_STRVAR(Seq_append_doc,
"Add a letter to the sequence object.\n"
"\n"
">>> my_seq = MutableSeq('ACTCGACGTCG')\n"
">>> my_seq.append('A')\n"
">>> my_seq\n"
"MutableSeq('ACTCGACGTCGA')\n"
"\n"
"No return value.\n"
"\n"
"A ValueError will be raised if the sequence is immutable.\n");

static PyObject *
Seq_append(SeqObject *self, PyObject *arg)
{
    char letter;
    PyObject *data = self->data;
    Py_ssize_t n = Py_SIZE(data);

    if (!PyByteArray_Check(data)) {
        PyErr_SetString(PyExc_ValueError, "sequence is immutable");
        return NULL;
    }

    if (n == PY_SSIZE_T_MAX) {
        PyErr_SetString(PyExc_OverflowError,
                        "cannot add more letters to sequence");
        return NULL;
    }

    if (PyUnicode_Check(arg)) {
        if (PyUnicode_READY(arg) == -1) return NULL;
        if (PyUnicode_GET_LENGTH(arg) == 1) {
            arg = PyUnicode_AsASCIIString(arg);
            if (!arg) return NULL;
            letter = PyBytes_AS_STRING(arg)[0];
            Py_DECREF(arg);
            if (PyByteArray_Resize(data, n + 1) < 0) return NULL;
            PyByteArray_AS_STRING(data)[n] = letter;
            Py_RETURN_NONE;
        }
    }

    PyErr_SetString(PyExc_ValueError, "expected a single letter");
    return NULL;
}

PyDoc_STRVAR(Seq_extend_doc,
"Extend a sequence object by a string or sequence.\n"
"\n"
">>> my_seq = MutableSeq('ACTCGACGTCG')\n"
">>> my_seq.extend('A')\n"
">>> my_seq\n"
"MutableSeq('ACTCGACGTCGA')\n"
">>> my_seq.extend('TTT')\n"
">>> my_seq\n"
"MutableSeq('ACTCGACGTCGATTT')\n"
"\n"
"No return value.\n"
"\n"
"A ValueError will be raised if the sequence is immutable.\n");

static PyObject *
Seq_extend(SeqObject *self, PyObject *arg)
{
    char *letters;
    PyObject *data = self->data;
    Py_ssize_t n = Py_SIZE(data);
    Py_ssize_t m;

    if (!PyByteArray_Check(data)) {
        PyErr_SetString(PyExc_ValueError, "sequence is immutable");
        return NULL;
    }

    if (PyUnicode_Check(arg)) {
        if (PyUnicode_READY(arg) == -1) return NULL;
        m = PyUnicode_GET_LENGTH(arg);
        if (n >= PY_SSIZE_T_MAX - m) {
            PyErr_SetString(PyExc_OverflowError,
                            "cannot add letters to sequence");
            return NULL;
        }

        arg = PyUnicode_AsASCIIString(arg);
        if (!arg) return NULL;
        letters = PyBytes_AS_STRING(arg);
        if (PyByteArray_Resize(data, n + m) < 0) {
            Py_DECREF(arg);
            return NULL;
        }
        memcpy(PyByteArray_AS_STRING(data) + n, letters, m);
        Py_DECREF(arg);
        Py_RETURN_NONE;
    }

    if ((PyObject*)self == arg) {
        /* special case; cannot be handled by the buffer code below */
        if (n >= PY_SSIZE_T_MAX - n) {
            PyErr_SetString(PyExc_OverflowError,
                            "cannot add letters to sequence");
            return NULL;
        }
        if (PyByteArray_Resize(data, 2*n) < 0) return NULL;
        letters = PyByteArray_AS_STRING(data);
        memcpy(letters + n, letters, n);
        Py_RETURN_NONE;
    }

    if (PyObject_IsInstance(arg, (PyObject*)&SeqType)) {
        Py_buffer view;
        if (PyObject_GetBuffer(arg, &view, PyBUF_STRIDES | PyBUF_FORMAT) < 0)
            return NULL;
        if (view.ndim != 1 || strcmp(view.format, "B") != 0) {
            PyBuffer_Release(&view);
            PyErr_SetString(PyExc_RuntimeError,
                            "expected a 1-dimensional sequence of bytes");
            return NULL;
        }
        m = view.len;
        if (n >= PY_SSIZE_T_MAX - m) {
            PyErr_SetString(PyExc_OverflowError,
                            "cannot add letters to sequence");
            return NULL;
        }
        letters = view.buf;
        if (PyByteArray_Resize(data, n + m) < 0) {
            PyBuffer_Release(&view);
            return NULL;
        }
        if (view.strides[0] == 1)
            memcpy(PyByteArray_AS_STRING(data) + n, letters, m);
        else if (view.strides[0] == 0)
            memset(PyByteArray_AS_STRING(data) + n, letters[0], m);
        else {
            PyErr_Format(PyExc_RuntimeError,
                         "unexpected stride %zd in Seq object",
                         view.strides[0]);
            PyBuffer_Release(&view);
            return NULL;
        }
        PyBuffer_Release(&view);
        Py_RETURN_NONE;
    }

    PyErr_SetString(PyExc_ValueError, "expected a string or a Seq object");
    return NULL;
}

PyDoc_STRVAR(Seq_insert_doc,
"Insert a letter into the sequence object at the specified index.\n"
"\n"
">>> my_seq = MutableSeq('ACTCGACGTCG')\n"
">>> my_seq.insert(0,'A')\n"
">>> my_seq\n"
"MutableSeq('AACTCGACGTCG')\n"
">>> my_seq.insert(8,'G')\n"
">>> my_seq\n"
"MutableSeq('AACTCGACGGTCG')\n"
"\n"
"No return value.\n"
"\n"
"A ValueError will be raised if the sequence is immutable.\n");

static PyObject *
Seq_insert(SeqObject *self, PyObject *args)
{
    char letter;
    char *buf;
    PyObject *data = self->data;
    Py_ssize_t i;
    Py_ssize_t n = Py_SIZE(data);

    if (!PyByteArray_Check(data)) {
        PyErr_SetString(PyExc_ValueError, "sequence is immutable");
        return NULL;
    }

    if (n == PY_SSIZE_T_MAX) {
        PyErr_SetString(PyExc_OverflowError,
                        "cannot add more letters to sequence");
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "ns:insert", &i, &buf)) return NULL;

    letter = buf[0];
    if (letter == '\0' || (int)letter >= 128) {
        PyErr_SetString(PyExc_ValueError, "only ASCII letters are allowed");
        return NULL;
    }
    if (buf[1] != '\0') {
        PyErr_SetString(PyExc_ValueError, "expected a single letter");
        return NULL;
    }

    if (PyByteArray_Resize(data, n + 1) < 0) return NULL;
    buf = PyByteArray_AS_STRING(data);
    if (i < 0) {
        i += n;
        if (i < 0) i = 0;
    }
    else if (i > n) i = n;
    memmove(buf + i + 1, buf + i, n - i);
    buf[i] = letter;

    Py_RETURN_NONE;
}

PyDoc_STRVAR(Seq_pop_doc,
"Remove the letter at given index and return it.\n"
"\n"
">>> my_seq = MutableSeq('ACTCGACGTCG')\n"
">>> my_seq.pop()\n"
"'G'\n"
">>> my_seq\n"
"MutableSeq('ACTCGACGTC')\n"
">>> my_seq.pop()\n"
"'C'\n"
">>> my_seq\n"
"MutableSeq('ACTCGACGT')\n");

static PyObject *
Seq_pop(SeqObject *self, PyObject *args)
{
    char letter;
    char* buf;
    PyObject *data = self->data;
    const Py_ssize_t n = Py_SIZE(data);
    Py_ssize_t i = n - 1;

    if (!PyByteArray_Check(data)) {
        PyErr_SetString(PyExc_ValueError, "sequence is immutable");
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "|n:pop", &i)) return NULL;

    if (i == -1 && PyErr_Occurred()) NULL;
    if (i < 0) i += n;
    if (i < 0 || i >= n) {
        PyErr_SetString(PyExc_IndexError, "sequence index out of range");
        return NULL;
    }

    buf = PyByteArray_AS_STRING(data);
    letter = buf[i];

    memmove(buf + i, buf + i + 1, n - i - 1);

    if (PyByteArray_Resize(data, n - 1) < 0) return NULL;

    return PyUnicode_FromStringAndSize(&letter, 1);
}

PyDoc_STRVAR(Seq_count_doc,
"Return a non-overlapping count, like that of a python string.\n"
"\n"
"This behaves like the python string method of the same name,\n"
"which does a non-overlapping count!\n"
"\n"
"For an overlapping search, use the count_overlap() method.\n"
"\n"
"Returns an integer, the number of occurrences of substring\n"
"argument sub in the (sub)sequence given by [start:end].\n"
"Optional arguments start and end are interpreted as in slice\n"
"notation.\n"
"\n"
"Arguments:\n"
" - sub - a string or another Seq object to look for\n"
" - start - optional integer, slice start\n"
" - end - optional integer, slice end\n"
"\n"
"e.g.\n"
"\n"
">>> from Bio.Seq import Seq\n"
">>> my_seq = Seq(\"AAAATGA\")\n"
">>> print(my_seq.count(\"A\"))\n"
"5\n"
">>> print(my_seq.count(\"ATG\"))\n"
"1\n"
">>> print(my_seq.count(Seq(\"AT\")))\n"
"1\n"
">>> print(my_seq.count(\"AT\", 2, -1))\n"
"1\n"
"\n"
"HOWEVER, please note because python strings and Seq objects (and\n"
"MutableSeq objects) do a non-overlapping search, this may not give\n"
"the answer you expect:\n"
"\n"
">>> \"AAAA\".count(\"AA\")\n"
"2\n"
">>> print(Seq(\"AAAA\").count(\"AA\"))\n"
"2\n"
"\n"
"An overlapping search, as implemented in .count_overlap(),\n"
"would give the answer as three!\n"
"\n");

static int
index_converter(PyObject* argument, void* pointer)
{
    Py_ssize_t index;
    if (argument == Py_None) return 1;
    index = PyNumber_AsSsize_t(argument, NULL);
    if (index == -1 && PyErr_Occurred()) return 0;
    *((Py_ssize_t*)pointer) = index;
    return 1;
}

static PyObject *
Seq_count(SeqObject *self, PyObject *args)
{
    PyObject *data = self->data;
    Py_ssize_t start = 0;
    Py_ssize_t end = PY_SSIZE_T_MAX;
    PyObject *sub;
    PyObject* result = NULL;

    if (!PyArg_ParseTuple(args, "O|O&O&:count", &sub,
                                index_converter, &start,
                                index_converter, &end)) return NULL;

    if (PyUnicode_Check(sub)) {
        sub = PyUnicode_AsASCIIString(sub);
        if (!sub) return NULL;
    }
    else Py_INCREF(sub);

    if (PyObject_IsInstance(data, (PyObject*)&UndefinedSeqDataType)) {
       Py_ssize_t length = ((UndefinedSeqDataObject*)data)->length;
       if (end < 0) end += length;
       else if (end > length) end = length;
       if (start < 0) start += length;
       if (start < 0) start = 0;
       if (start < end) {
           char* s;
           Py_ssize_t i;
           Py_ssize_t stride;
           Py_buffer view;
           char letter = ((UndefinedSeqDataObject*)data)->character;
           if (PyObject_GetBuffer(sub, &view, PyBUF_STRIDES | PyBUF_FORMAT) < 0)
               goto exit;
           s = view.buf;
           stride = view.strides[0];
           for (i = 0; i < view.len; i++) {
               if (s[i*stride] != letter)  {
                   result = PyLong_FromLong(0);
                   break;
               }
           }
           if (!result)
               result = PyLong_FromLong((end-start)/view.len);
           PyBuffer_Release(&view);
       }
       else
           result = PyLong_FromLong(0);
    }
    else {
        if (!PyBytes_Check(data) && ! PyByteArray_Check(data)) {
            if (!PySequence_Check(data)) {
                PyErr_SetString(PyExc_RuntimeError,
                                "data should support the sequence protocol");
                goto exit;
            }
            self->data = PySequence_GetSlice(data, 0, PY_SSIZE_T_MAX);
            if (!self->data) {
                self->data = data;
                goto exit;
            }
            Py_DECREF(data);
        }
        result = PyObject_CallMethod(data, "count", "Onn", sub, start, end);
    }

exit:
    Py_DECREF(sub);
    return result;
}

static PyMethodDef Seq_methods[] = {
    {"__reduce__", (PyCFunction)Seq_reduce, METH_NOARGS, Seq_reduce_doc},
    {"reverse", (PyCFunction)Seq_reverse, METH_NOARGS, Seq_reverse_doc},
    {"complement", (PyCFunction)Seq_complement, METH_NOARGS, Seq_complement_doc},
    {"rna_complement", (PyCFunction)Seq_rna_complement, METH_NOARGS, Seq_rna_complement_doc},
    {"append", (PyCFunction)Seq_append, METH_O, Seq_append_doc},
    {"extend", (PyCFunction)Seq_extend, METH_O, Seq_extend_doc},
    {"insert", (PyCFunction)Seq_insert, METH_VARARGS, Seq_insert_doc},
    {"pop", (PyCFunction)Seq_pop, METH_VARARGS, Seq_pop_doc},
    {"count", (PyCFunction)Seq_count, METH_VARARGS, Seq_count_doc},
    {NULL}  /* Sentinel */
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
    &Seq_as_mapping,                            /* tp_as_mapping */
    0,                                          /* tp_hash */
    0,                                          /* tp_call */
    (reprfunc)Seq_str,                          /* tp_str */
    (getattrofunc)Seq_getattro,                 /* tp_getattro */
    (setattrofunc)Seq_setattro,                 /* tp_setattro */
    &Seq_as_buffer,                             /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE |
        Py_TPFLAGS_BYTES_SUBCLASS,              /* tp_flags */
    Seq_doc,                                    /* tp_doc */
    0,                                          /* tp_traverse */
    0,                                          /* tp_clear */
    (richcmpfunc)Seq_richcompare,               /* tp_richcompare */
    0,                                          /* tp_weaklistoffset */
    0,                                          /* tp_iter */
    0,                                          /* tp_iternext */
    Seq_methods,                                /* tp_methods */
    0,                                          /* tp_members */
    0,                                          /* tp_getset */
    0,                                          /* tp_base */
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
    if (PyType_Ready(&UndefinedSeqDataType) < 0)
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
