#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdint.h>

void
integer_unpack_u8(Py_buffer *in_view, Py_buffer *out_view)
{
    Py_ssize_t in_size = in_view->shape[0];
    Py_ssize_t in_index = 0;
    Py_ssize_t out_index = 0;

    uint8_t *in_data = in_view->buf;
    uint32_t *out_data = out_view->buf;

    while (in_index < in_size) {
        uint32_t sum = in_data[in_index];

        if (sum == UINT8_MAX) {
            while (in_index + 1 < in_size) {
                in_index += 1;
                sum += in_data[in_index];

                if (in_data[in_index] != UINT8_MAX) {
                    break;
                }
            }
        }

        out_data[out_index] = sum;
        in_index += 1;
        out_index += 1;
    }
}

void
integer_unpack_u16(Py_buffer *in_view, Py_buffer *out_view)
{
    Py_ssize_t in_size = in_view->shape[0];
    Py_ssize_t in_index = 0;
    Py_ssize_t out_index = 0;

    uint16_t *in_data = in_view->buf;
    uint32_t *out_data = out_view->buf;

    while (in_index < in_size) {
        uint32_t sum = in_data[in_index];

        if (sum == UINT16_MAX) {
            while (in_index + 1 < in_size) {
                in_index += 1;
                sum += in_data[in_index];

                if (in_data[in_index] != UINT16_MAX) {
                    break;
                }
            }
        }

        out_data[out_index] = sum;
        in_index += 1;
        out_index += 1;
    }
}

void
integer_unpack_i8(Py_buffer *in_view, Py_buffer *out_view)
{
    Py_ssize_t in_size = in_view->shape[0];
    Py_ssize_t in_index = 0;
    Py_ssize_t out_index = 0;

    int8_t *in_data = in_view->buf;
    int32_t *out_data = out_view->buf;

    while (in_index < in_size) {
        int32_t sum = in_data[in_index];

        if (sum == INT8_MAX || sum == INT8_MIN) {
            while (in_index + 1 < in_size) {
                in_index += 1;
                sum += in_data[in_index];

                if (in_data[in_index] != INT8_MAX && in_data[in_index] != INT8_MIN) {
                    break;
                }
            }
        }

        out_data[out_index] = sum;
        in_index += 1;
        out_index += 1;
    }
}

void
integer_unpack_i16(Py_buffer *in_view, Py_buffer *out_view)
{
    Py_ssize_t in_size = in_view->shape[0];
    Py_ssize_t in_index = 0;
    Py_ssize_t out_index = 0;

    int16_t *in_data = in_view->buf;
    int32_t *out_data = out_view->buf;

    while (in_index < in_size) {
        int32_t sum = in_data[in_index];

        if (sum == INT16_MAX || sum == INT16_MIN) {
            while (in_index + 1 < in_size) {
                in_index += 1;
                sum += in_data[in_index];

                if (in_data[in_index] != INT16_MAX && in_data[in_index] != INT16_MIN) {
                    break;
                }
            }
        }

        out_data[out_index] = sum;
        in_index += 1;
        out_index += 1;
    }
}

static PyObject *
integer_unpack(PyObject *self, PyObject *args)
{
    PyObject *in = NULL;
    PyObject *out = NULL;

    if (!PyArg_ParseTuple(args, "OO", &in, &out)) {
        return NULL;
    }

    Py_buffer in_view, out_view;
    const int flags = PyBUF_ND | PyBUF_FORMAT;

    if (PyObject_GetBuffer(in, &in_view, flags) != 0) {
        return NULL;
    }
    if (PyObject_GetBuffer(out, &out_view, flags | PyBUF_WRITABLE) != 0) {
        PyBuffer_Release(&in_view);
        return NULL;
    }

    if (in_view.ndim != 1) {
        PyErr_SetString(PyExc_ValueError, "First argument should be one-dimensional.");
        goto exit;
    }
    if (out_view.ndim != 1) {
        PyErr_SetString(PyExc_ValueError, "Second argument should be one-dimensional.");
        goto exit;
    }

    const char format = in_view.format[0];

    if (format == 'B') {
        integer_unpack_u8(&in_view, &out_view);
    }
    else if (format == 'H') {
        integer_unpack_u16(&in_view, &out_view);
    }
    else if (format == 'b') {
        integer_unpack_i8(&in_view, &out_view);
    }
    else if (format == 'h') {
        integer_unpack_i16(&in_view, &out_view);
    }
    else {
        PyErr_Format(PyExc_ValueError,
            "Unexpected buffer format: %s",
            in_view.format);
    }

exit:
    PyBuffer_Release(&in_view);
    PyBuffer_Release(&out_view);
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

    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    return m;
}
