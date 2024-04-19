#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"
#include <math.h>

void
integer_unpack_u8(PyObject *in, PyObject *out)
{
    npy_intp in_size = PyArray_SIZE(in);
    npy_intp in_index = 0;
    npy_intp out_index = 0;

    const npy_uint32 max_value = (npy_uint32) NPY_MAX_UINT8;

    while (in_index < in_size) {
        npy_uint8 *in_p = PyArray_GETPTR1(in, in_index);
        npy_uint32 sum = (npy_uint32) *in_p;

        if (sum == max_value) {
            while (in_index + 1 < in_size) {
                in_index += 1;
                in_p = PyArray_GETPTR1(in, in_index);
                sum += (npy_uint32) *in_p;

                if ((npy_uint32) *in_p != max_value) {
                    break;
                }
            }
        }

        npy_uint32 *out_p = PyArray_GETPTR1(out, out_index);
        *out_p = sum;
        in_index += 1;
        out_index += 1;
    }
}

void
integer_unpack_u16(PyObject *in, PyObject *out)
{
    npy_intp in_size = PyArray_SIZE(in);
    npy_intp in_index = 0;
    npy_intp out_index = 0;

    const npy_uint32 max_value = (npy_uint32) NPY_MAX_UINT16;

    while (in_index < in_size) {
        npy_uint16 *in_p = PyArray_GETPTR1(in, in_index);
        npy_uint32 sum = (npy_uint32) *in_p;

        if (sum == max_value) {
            while (in_index + 1 < in_size) {
                in_index += 1;
                in_p = PyArray_GETPTR1(in, in_index);
                sum += (npy_uint32) *in_p;

                if ((npy_uint32) *in_p != max_value) {
                    break;
                }
            }
        }

        npy_uint32 *out_p = PyArray_GETPTR1(out, out_index);
        *out_p = sum;
        in_index += 1;
        out_index += 1;
    }
}

void
integer_unpack_i8(PyObject *in, PyObject *out)
{
    npy_intp in_size = PyArray_SIZE(in);
    npy_intp in_index = 0;
    npy_intp out_index = 0;

    const npy_int8 min_value = (npy_int8) NPY_MIN_INT8;
    const npy_int32 max_value = (npy_int32) NPY_MAX_INT8;

    while (in_index < in_size) {
        npy_int8 *in_p = PyArray_GETPTR1(in, in_index);
        npy_int32 sum = (npy_int32) *in_p;

        if (sum == max_value || sum == min_value) {
            while (in_index + 1 < in_size) {
                in_index += 1;
                in_p = PyArray_GETPTR1(in, in_index);
                sum += (npy_int32) *in_p;
                npy_int32 value = (npy_int32) *in_p;

                if (value != max_value && value != min_value) {
                    break;
                }
            }
        }

        npy_int32 *out_p = PyArray_GETPTR1(out, out_index);
        *out_p = sum;
        in_index += 1;
        out_index += 1;
    }
}

void
integer_unpack_i16(PyObject *in, PyObject *out)
{
    npy_intp in_size = PyArray_SIZE(in);
    npy_intp in_index = 0;
    npy_intp out_index = 0;

    const npy_int32 min_value = (npy_int32) NPY_MIN_INT16;
    const npy_int32 max_value = (npy_int32) NPY_MAX_INT16;

    while (in_index < in_size) {
        npy_int16 *in_p = PyArray_GETPTR1(in, in_index);
        npy_int32 sum = (npy_int32) *in_p;

        if (sum == max_value || sum == min_value) {
            while (in_index + 1 < in_size) {
                in_index += 1;
                in_p = PyArray_GETPTR1(in, in_index);
                sum += (npy_int32) *in_p;
                npy_int32 value = (npy_int32) *in_p;

                if (value != max_value && value != min_value) {
                    break;
                }
            }
        }

        npy_int32 *out_p = PyArray_GETPTR1(out, out_index);
        *out_p = sum;
        in_index += 1;
        out_index += 1;
    }
}

static PyObject *
integer_unpack(PyObject *self, PyObject *args)
{
    PyObject *in = NULL;
    PyObject *out = NULL;

    if (!PyArg_ParseTuple(args, "O!O!", &PyArray_Type,
                          &in, &PyArray_Type, &out)) {
        return NULL;
    }

    npy_intp itemsize = PyArray_ITEMSIZE(in);

    if (PyArray_ISUNSIGNED(in)) {
        if (itemsize == 1) {
            integer_unpack_u8(in, out);
        }
        else {
            integer_unpack_u16(in, out);
        }
    }
    else {
        if (itemsize == 1) {
            integer_unpack_i8(in, out);
        }
        else {
            integer_unpack_i16(in, out);
        }
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef IntegerUnpackMethods[] = {
    {"integer_unpack", integer_unpack, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_bcif_helper",
    NULL,
    -1,
    IntegerUnpackMethods
};

PyMODINIT_FUNC
PyInit__bcif_helper(void)
{
    PyObject *m;

    import_array();
    import_umath();

    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    return m;
}
