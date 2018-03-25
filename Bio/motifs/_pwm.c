#include <Python.h>


static void
calculate(const char sequence[], int s, Py_ssize_t m, double* matrix,
          Py_ssize_t n, float* scores)
{
    Py_ssize_t i, j;
    char c;
    double score;
    int ok;
    float* p = scores;
    float nan = 0.0;
    nan /= nan;
    for (i = 0; i < n; i++)
    {
        score = 0.0;
        ok = 1;
        for (j = 0; j < m; j++)
        {
            c = sequence[i+j];
            switch (c)
            {
              /* Handling mixed case input here rather than converting it to
                 uppercase in Python code first, since doing so could use too
                 much memory if sequence is too long (e.g. chromosome or
                 plasmid). */
                case 'A':
                case 'a':
                    score += matrix[j*4+0]; break;
                case 'C':
                case 'c':
                    score += matrix[j*4+1]; break;
                case 'G':
                case 'g':
                    score += matrix[j*4+2]; break;
                case 'T':
                case 't':
                    score += matrix[j*4+3]; break;
                default:
                    ok = 0;
            }
        }
        if (ok) *p = (float)score;
        else *p = nan;
        p++;
    }
}

static int
matrix_converter(PyObject* object, void* address)
{
    const int flags = PyBUF_C_CONTIGUOUS | PyBUF_FORMAT;
    char datatype;
    Py_buffer* view = address;
    if (PyObject_GetBuffer(object, view, flags) == -1) {
        PyErr_SetString(PyExc_RuntimeError,
                        "position-weight matrix is not an array");
        return 0;
    }
    datatype = view->format[0];
    switch (datatype) {
        case '@':
        case '=':
        case '<':
        case '>':
        case '!': datatype = view->format[1]; break;
        default: break;
    }
    if (datatype != 'd') {
        PyErr_Format(PyExc_RuntimeError,
            "position-weight matrix data format incorrect ('%c', expected 'd')",
            datatype);
        return 0;
    }
    if (view->ndim != 2) {
        PyErr_Format(PyExc_RuntimeError,
            "position-weight matrix has incorrect rank (%d expected 2)",
            view->ndim);
        return 0;
    }
    if (view->shape[1] != 4) {
        PyErr_Format(PyExc_RuntimeError,
            "position-weight matrix should have four columns "
            "(%zd columns found)", view->shape[1]);
        return 0;
    }
    return 1;
}

static int
scores_converter(PyObject* object, void* address)
{
    const int flags = PyBUF_C_CONTIGUOUS | PyBUF_FORMAT;
    char datatype;
    Py_buffer* view = address;
    if (PyObject_GetBuffer(object, view, flags) == -1)
        return 0;
    datatype = view->format[0];
    switch (datatype) {
        case '@':
        case '=':
        case '<':
        case '>':
        case '!': datatype = view->format[1]; break;
        default: break;
    }
    if (datatype != 'f') {
        PyErr_Format(PyExc_RuntimeError,
            "scores array has incorrect data format ('%c', expected 'f')",
            datatype);
        return 0;
    }
    if (view->ndim != 1) {
        PyErr_Format(PyExc_ValueError,
            "scores array has incorrect rank (%d expected 1)",
            view->ndim);
        return 0;
    }
    return 1;
}

static char calculate__doc__[] =
"    calculate(sequence, pwm) -> array of score values\n"
"\n"
"This function calculates the position-weight matrix scores for all\n"
"positions along the sequence, and returns them as a Numerical Python\n"
"array.\n";

static PyObject*
py_calculate(PyObject* self, PyObject* args, PyObject* keywords)
{
    const char* sequence;
    static char* kwlist[] = {"sequence", "matrix", "scores", NULL};
    Py_ssize_t m;
    Py_ssize_t n;
    int s;
    PyObject* result = NULL;
    Py_buffer scores;
    Py_buffer matrix;
    matrix.obj = NULL;
    scores.obj = NULL;
    if(!PyArg_ParseTupleAndKeywords(args, keywords, "s#O&O&", kwlist,
                                    &sequence,
                                    &s,
                                    matrix_converter, &matrix,
                                    scores_converter, &scores)) goto exit;
    m = matrix.shape[0];
    n = scores.shape[0];
    if (n != s - m + 1) {
        PyErr_SetString(PyExc_RuntimeError,
                        "size of scores array is inconsistent");
        goto exit;
    }
    calculate(sequence, s, m, matrix.buf, n, scores.buf);
    Py_INCREF(Py_None);
    result = Py_None;
exit:
    if (matrix.obj) PyBuffer_Release(&matrix);
    if (scores.obj) PyBuffer_Release(&scores);
    return result;
}

static struct PyMethodDef methods[] = {
   {"calculate", (PyCFunction)py_calculate, METH_VARARGS | METH_KEYWORDS, calculate__doc__},
   {NULL,          NULL, 0, NULL} /* sentinel */
};


#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_pwm",
        "Fast calculations involving position-weight matrices",
        -1,
        methods,
        NULL,
        NULL,
        NULL,
        NULL
};

PyObject*
PyInit__pwm(void)

#else

void init_pwm(void)
#endif
{
  PyObject *m;
#if PY_MAJOR_VERSION >= 3
  m = PyModule_Create(&moduledef);
  if (m==NULL) return NULL;
#else
  m = Py_InitModule4("_pwm",
                     methods,
                     "Fast calculations involving position-weight matrices",
                     NULL,
                     PYTHON_API_VERSION);
  if (m==NULL) return;
#endif

  if (PyErr_Occurred()) Py_FatalError("can't initialize module _pwm");
#if PY_MAJOR_VERSION >= 3
    return m;
#endif
}
