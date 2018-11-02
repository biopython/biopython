/* Copyright 2018 by Michiel de Hoon.  All rights reserved.
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */


#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "float.h"


#define HORIZONTAL 0x1
#define VERTICAL 0x2
#define DIAGONAL 0x4
#define ENDPOINT 0x8
#define M_MATRIX 0x1
#define Ix_MATRIX 0x2
#define Iy_MATRIX 0x4
#define DONE 0x8
#define NONE 0x7

#define OVERFLOW_ERROR -1
#define MEMORY_ERROR -2

#define CHARINDEX(s) (c = s, c >= 'a' ? c - 'a' : c - 'A')


typedef enum {NeedlemanWunschSmithWaterman,
              Gotoh,
              WatermanSmithBeyer,
              Unknown} Algorithm;

typedef enum {Global, Local} Mode;

typedef struct {
    unsigned char trace : 4;
    unsigned char path : 4;
} Trace;

typedef struct {
    unsigned char M : 4;
    unsigned char Ix : 4;
    unsigned char Iy : 4;
    unsigned char path : 4;
} Trace3;

typedef struct {
    double score;
    unsigned int trace : 4;
    unsigned char path : 4;
    unsigned int step;
} CellM; /* Used for the Waterman-Smith-Beyer algorithm. */

typedef struct {
    double score;
    int* traceM;
    int* traceXY;
    unsigned char path : 4;
    unsigned int step;
} CellXY; /* Used for the Waterman-Smith-Beyer algorithm. */

static int _convert_single_letter(PyObject* item)
{
    int i;
    char letter = '\0';
    Py_buffer view;
#if PY_MAJOR_VERSION >= 3
    if (PyUnicode_Check(item)) {
        Py_UCS1* data;
        if (PyUnicode_READY(item) == -1) return -1;
        switch (PyUnicode_KIND(item)) {
            case PyUnicode_1BYTE_KIND: break;
            case PyUnicode_2BYTE_KIND:
            case PyUnicode_4BYTE_KIND:
            case PyUnicode_WCHAR_KIND:
                PyErr_SetString(PyExc_ValueError,
                                "expected an ASCII character");
                return -1;
            default:
                PyErr_SetString(PyExc_SystemError,
                                "unknown PyUnicode kind constant");
                return -1;
        }
        data = PyUnicode_1BYTE_DATA(item);
        letter = *((char*)(data));
    }
    else
#endif
    {
        if (!PyObject_CheckBuffer(item)
         || PyObject_GetBuffer(item, &view, PyBUF_FORMAT) == -1
         || strcmp(view.format, "B") != 0
         || view.len != 1) {
            PyErr_SetString(PyExc_ValueError, "expected a single letter");
            return -1;
        }
        letter = *((char*)(view.buf));
    }
    if (letter >= 'a' && letter <= 'z') i = letter - 'a';
    else if (letter >= 'A' && letter <= 'Z') i = letter - 'A';
    else {
        PyErr_SetString(PyExc_ValueError, "expected an ASCII character");
        return -1;
    }
    return i;
}

static PyObject*
_create_path_needleman_wunsch_smith_waterman(Trace** M, int i, int j) {
    PyObject* tuple;
    PyObject* row;
    PyObject* value;
    int path;
    const int ii = i;
    const int jj = j;
    int n = 1;
    int direction = 0;

    while (1) {
        path = M[i][j].path;
        if (!path) break;
        if (path != direction) {
            n++;
            direction = path;
        }
        switch (path) {
            case HORIZONTAL: j++; break;
            case VERTICAL: i++; break;
            case DIAGONAL: i++; j++; break;
        }
    }
    i = ii;
    j = jj;

    direction = 0;
    tuple = PyTuple_New(n);
    if (!tuple) return NULL;
    n = 0;
    while (1) {
        path = M[i][j].path;
        if (path != direction) {
            row = PyTuple_New(2);
            if (!row) break;
#if PY_MAJOR_VERSION >= 3
            value = PyLong_FromLong(i);
#else
            value = PyInt_FromLong(i);
#endif
            if (!value) {
                Py_DECREF(row); /* all references were stolen */
                break;
            }
            PyTuple_SET_ITEM(row, 0, value);
#if PY_MAJOR_VERSION >= 3
            value = PyLong_FromLong(j);
#else
            value = PyInt_FromLong(j);
#endif
            if (!value) {
                Py_DECREF(row); /* all references were stolen */
                break;
            }
            PyTuple_SET_ITEM(row, 1, value);
            PyTuple_SET_ITEM(tuple, n, row);
            n++;
            direction = path;
        }
        switch (path) {
            case HORIZONTAL: j++; break;
            case VERTICAL: i++; break;
            case DIAGONAL: i++; j++; break;
            default: return tuple;
        }
    }
    Py_DECREF(tuple); /* all references were stolen */
    return NULL;
}

static PyObject*
_create_path_gotoh(Trace3** trace3, int i, int j) {
    PyObject* tuple;
    PyObject* row;
    PyObject* value;
    int path;
    const int ii = i;
    const int jj = j;
    int n = 1;
    int direction = 0;

    path = trace3[i][j].path;
    while (path) {
        if (path != direction) {
            n++;
            direction = path;
        }
        switch (path) {
            case HORIZONTAL: j++; break;
            case VERTICAL: i++; break;
            case DIAGONAL: i++; j++; break;
        }
        path = trace3[i][j].path;
    }
    i = ii;
    j = jj;

    direction = 0;
    tuple = PyTuple_New(n);
    if (!tuple) return NULL;
    n = 0;
    path = trace3[i][j].path;
    while (1) {
        if (path != direction) {
            row = PyTuple_New(2);
            if (!row) break;
#if PY_MAJOR_VERSION >= 3
            value = PyLong_FromLong(i);
#else
            value = PyInt_FromLong(i);
#endif
            if (!value) {
                Py_DECREF(row); /* all references were stolen */
                break;
            }
            PyTuple_SET_ITEM(row, 0, value);
#if PY_MAJOR_VERSION >= 3
            value = PyLong_FromLong(j);
#else
            value = PyInt_FromLong(j);
#endif
            if (!value) {
                Py_DECREF(row); /* all references were stolen */
                break;
            }
            PyTuple_SET_ITEM(row, 1, value);
            PyTuple_SET_ITEM(tuple, n, row);
            n++;
            direction = path;
        }
        switch (path) {
            case HORIZONTAL: j++; break;
            case VERTICAL: i++; break;
            case DIAGONAL: i++; j++; break;
            default: return tuple;
        }
        path = trace3[i][j].path;
    }
    Py_DECREF(tuple); /* all references were stolen */
    return NULL;
}

static PyObject*
_create_path_waterman_smith_beyer(CellM** M, CellXY** Ix, CellXY** Iy,
                                  int i, int j) {
    PyObject* tuple;
    PyObject* row;
    PyObject* value;
    int i1 = i;
    int j1 = j;
    int i2;
    int j2;
    int n = 1;
    int direction = 0;
    int path = M[i1][j1].path;
    int step = M[i1][j1].step;

    while (path) {
        switch (path) {
            case HORIZONTAL:
                i2 = i1;
                j2 = j1 + step;
                if (direction != HORIZONTAL) {
                    n++;
                    direction = HORIZONTAL;
                }
                i1 = i2;
                j1 = j2;
                path = Iy[i1][j1].path;
                step = Iy[i1][j1].step;
                break;
            case VERTICAL:
                i2 = i1 + step;
                j2 = j1;
                if (direction != VERTICAL) {
                    n++;
                    direction = VERTICAL;
                }
                i1 = i2;
                j1 = j2;
                path = Ix[i1][j1].path;
                step = Ix[i1][j1].step;
                break;
            case DIAGONAL:
                i2 = i1 + 1;
                j2 = j1 + 1;
                if (direction != DIAGONAL) {
                    n++;
                    direction = DIAGONAL;
                }
                i1 = i2;
                j1 = j2;
                path = M[i1][j1].path;
                step = M[i1][j1].step;
                break;
            default:
                PyErr_SetString(PyExc_RuntimeError,
                    "Unexpected path in _create_path_waterman_smith_beyer");
                return NULL;
        }
    }

    i1 = i;
    j1 = j;
    path = M[i1][j1].path;
    switch (path) {
        case HORIZONTAL:
            i2 = i1;
            j2 = j1 + M[i1][j1].step;
            break;
        case VERTICAL:
            i2 = i1 + M[i1][j1].step;
            j2 = j1;
            break;
        case DIAGONAL:
            i2 = i1 + 1;
            j2 = j1 + 1;
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected path in _create_path_waterman_smith_beyer");
            return NULL;
    }
    tuple = PyTuple_New(n);
    if (!tuple) return NULL;
    n = 0;
    direction = 0;
    while (1) {
        if ((i1==i2 && direction != HORIZONTAL)
         || (j1==j2 && direction != VERTICAL)
         || (direction != DIAGONAL)
         || (i2 == -1)) {
            row = PyTuple_New(2);
            if (!row) break;
#if PY_MAJOR_VERSION >= 3
            value = PyLong_FromLong(i1);
#else
            value = PyInt_FromLong(i1);
#endif
            if (!value) {
                Py_DECREF(row); /* all references were stolen */
                break;
            }
            PyTuple_SET_ITEM(row, 0, value);
#if PY_MAJOR_VERSION >= 3
            value = PyLong_FromLong(j1);
#else
            value = PyInt_FromLong(j1);
#endif
            PyTuple_SET_ITEM(row, 1, value);
            PyTuple_SET_ITEM(tuple, n, row);
            n++;
        }
        if (i2 < 0) return tuple;
        else if (i1==i2) {
            i1 = i2;
            j1 = j2;
            if (Iy[i1][j1].path == HORIZONTAL) {
                i2 = i1;
                j2 = j1 + Iy[i1][j1].step;
            }
            else if (Iy[i1][j1].path == VERTICAL) {
                i2 = i1 + Iy[i1][j1].step;
                j2 = j1;
            }
            else if (Iy[i1][j1].path == DIAGONAL) {
                i2 = i1 + 1;
                j2 = j1 + 1;
            }
            else if (Iy[i1][j1].path == 0) {
                i2 = -1;
            }
            else printf("RUNTIME ERROR 42\n");
            direction = HORIZONTAL;
        }
        else if (j1==j2) {
            i1 = i2;
            j1 = j2;
            if (Ix[i1][j1].path == HORIZONTAL) {
                i2 = i1;
                j2 = j1 + Ix[i1][j1].step;
            }
            else if (Ix[i1][j1].path == VERTICAL) {
                i2 = i1 + Ix[i1][j1].step;
                j2 = j1;
            }
            else if (Ix[i1][j1].path == DIAGONAL) {
                i2 = i1 + 1;
                j2 = j1 + 1;
            }
            else if (Ix[i1][j1].path == 0) {
                i2 = -1;
            }
            else printf("RUNTIME ERROR 43\n");
            direction = VERTICAL;
        }
        else {
            i1 = i2;
            j1 = j2;
            if (M[i1][j1].path == HORIZONTAL) {
                i2 = i1;
                j2 = j1 + M[i1][j1].step;
            }
            else if (M[i1][j1].path == VERTICAL) {
                i2 = i1 + M[i1][j1].step;
                j2 = j1;
            }
            else if (M[i1][j1].path == DIAGONAL) {
                i2 = i1 + 1;
                j2 = j1 + 1;
            }
            else if (M[i1][j1].path == 0) {
                i2 = -1;
            }
            else printf("RUNTIME ERROR 50\n");
            direction = DIAGONAL;
        }
    }
    Py_DECREF(tuple); /* all references were stolen */
    return NULL;
}

typedef struct {
    PyObject_HEAD
    union { Trace** trace1; CellM** general; } M;
    Trace3** trace3;
    CellXY** Ix; 
    CellXY** Iy; 
    int nA;
    int nB;
    int iA;
    int iB;
    Mode mode;
    Algorithm algorithm;
    Py_ssize_t length;
} PathGenerator;

typedef struct {
    PyObject_HEAD
    Mode mode;
    Algorithm algorithm;
    double match;
    double mismatch;
    double epsilon;
    double target_open_gap_score;
    double target_extend_gap_score;
    double target_left_open_gap_score;
    double target_left_extend_gap_score;
    double target_right_open_gap_score;
    double target_right_extend_gap_score;
    double query_open_gap_score;
    double query_extend_gap_score;
    double query_left_open_gap_score;
    double query_left_extend_gap_score;
    double query_right_open_gap_score;
    double query_right_extend_gap_score;
    PyObject* target_gap_function;
    PyObject* query_gap_function;
    double substitution_matrix[26][26]; /* 26 letters in the alphabet */
    int* letters;
} Aligner;

static int
Aligner_init(Aligner *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"match", "mismatch", NULL};
    int i, j;
    const int n = 26;

    self->mode = Global;
    self->match = 1.0;
    self->mismatch = 0.0;
    self->epsilon = 1.e-6;
    self->target_open_gap_score = 0;
    self->target_extend_gap_score = 0;
    self->query_open_gap_score = 0;
    self->query_extend_gap_score = 0;
    self->target_left_open_gap_score = 0;
    self->target_left_extend_gap_score = 0;
    self->target_right_open_gap_score = 0;
    self->target_right_extend_gap_score = 0;
    self->query_left_open_gap_score = 0;
    self->query_left_extend_gap_score = 0;
    self->query_right_open_gap_score = 0;
    self->query_right_extend_gap_score = 0;
    self->target_gap_function = NULL;
    self->query_gap_function = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dd", kwlist,
                                     &self->match, &self->mismatch))
        return -1;
    for (i = 0; i < n; i++) {
        self->substitution_matrix[i][i] = self->match;
        for (j = 0; j < i; j++) {
            self->substitution_matrix[i][j] = self->mismatch;
            self->substitution_matrix[j][i] = self->mismatch;
        }
    }
    i = 'X' - 'A';
    self->substitution_matrix[i][i] = 0.0;
    self->letters = NULL;
    self->algorithm = Unknown;
    return 0;
}

static void
Aligner_dealloc(Aligner* self)
{   if (self->letters) PyMem_Free(self->letters);
    Py_XDECREF(self->target_gap_function);
    Py_XDECREF(self->query_gap_function);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject*
Aligner_repr(Aligner* self)
{
  const char text[] = "Pairwise aligner, implementing the Needleman-Wunsch, Smith-Waterman, Gotoh, and Waterman-Smith-Beyer global and local alignment algorithms";
#if PY_MAJOR_VERSION >= 3
  return PyUnicode_FromString(text);
#else
  return PyString_FromString(text);
#endif
}

static PyObject*
Aligner_str(Aligner* self)
{
    int n;
    char text[1024];
    char* p = text;
    n = sprintf(text, "Pairwise sequence aligner with parameters\n");
    p += n;
    if (self->letters) {
        n = sprintf(p, "  match/mismatch score: <substitution matrix>\n");
        p += n;
    } else {
        n = sprintf(p, "  match score: %f\n", self->match);
        p += n;
        n = sprintf(p, "  mismatch score: %f\n", self->mismatch);
        p += n;
    }
    if (self->target_gap_function) {
#if PY_MAJOR_VERSION >= 3
        n = sprintf(p, "  target gap function: %%R\n");
        p += n;
#else
        char* s;
        PyObject* representation = PyObject_Repr(self->target_gap_function);
        if (!representation) {
            PyErr_SetString(PyExc_MemoryError, "Out of memory");
            return NULL;
        }
        s = PyString_AsString(representation);
        n = sprintf(p, "  target gap function: %s\n", s);
        p += n;
        Py_DECREF(representation);
#endif
    }
    else {
        n = sprintf(p, "  target open gap score: %f\n",
                       self->target_open_gap_score);
        p += n;
        n = sprintf(p, "  target extend gap score: %f\n",
                       self->target_extend_gap_score);
        p += n;
        n = sprintf(p, "  target left open gap score: %f\n",
                       self->target_left_open_gap_score);
        p += n;
        n = sprintf(p, "  target left extend gap score: %f\n",
                       self->target_left_extend_gap_score);
        p += n;
        n = sprintf(p, "  target right open gap score: %f\n",
                       self->target_right_open_gap_score);
        p += n;
        n = sprintf(p, "  target right extend gap score: %f\n",
                       self->target_right_extend_gap_score);
        p += n;
    }
    if (self->query_gap_function) {
#if PY_MAJOR_VERSION >= 3
        n = sprintf(p, "  query gap function: %%R\n");
        p += n;
#else
        char* s;
        PyObject* representation = PyObject_Repr(self->query_gap_function);
        if (!representation) {
            PyErr_SetString(PyExc_MemoryError, "Out of memory");
            return NULL;
        }
        s = PyString_AsString(representation);
        n = sprintf(p, "  query gap function: %s\n", s);
        p += n;
        Py_DECREF(representation);
#endif
    }
    else {
        n = sprintf(p, "  query open gap score: %f\n",
                       self->query_open_gap_score);
        p += n;
        n = sprintf(p, "  query extend gap score: %f\n",
                       self->query_extend_gap_score);
        p += n;
        n = sprintf(p, "  query left open gap score: %f\n",
                       self->query_left_open_gap_score);
        p += n;
        n = sprintf(p, "  query left extend gap score: %f\n",
                       self->query_left_extend_gap_score);
        p += n;
        n = sprintf(p, "  query right open gap score: %f\n",
                       self->query_right_open_gap_score);
        p += n;
        n = sprintf(p, "  query right extend gap score: %f\n",
                       self->query_right_extend_gap_score);
        p += n;
    }
    switch (self->mode) {
        case Global: n = sprintf(p, "  mode: global\n"); break;
        case Local: n = sprintf(p, "  mode: local\n"); break;
    }
    p += n;
#if PY_MAJOR_VERSION >= 3
    if (self->target_gap_function || self->query_gap_function)
        return PyUnicode_FromFormat(text, self->target_gap_function, self->query_gap_function);
    else if (self->target_gap_function)
        return PyUnicode_FromFormat(text, self->target_gap_function);
    else if (self->query_gap_function)
        return PyUnicode_FromFormat(text, self->query_gap_function);
    else
        return PyUnicode_FromString(text);
#else
    return PyString_FromString(text);
#endif
}

static char Aligner_mode__doc__[] = "alignment mode (global or local)";

static PyObject*
Aligner_get_mode(Aligner* self, void* closure)
{   const char* message = NULL;
    switch (self->mode) {
        case Global: message = "global"; break;
        case Local: message = "local"; break;
    }
#if PY_MAJOR_VERSION >= 3
    return PyUnicode_FromString(message);
#else
    return PyString_FromString(message);
#endif
}

static int
Aligner_set_mode(Aligner* self, PyObject* value, void* closure)
{
#if PY_MAJOR_VERSION >= 3
    if (!PyUnicode_Check(value)) {
#else
    char* mode;
    if (!PyString_Check(value)) {
#endif
        PyErr_SetString(PyExc_TypeError, "invalid mode");
        return -1;
    }
#if PY_MAJOR_VERSION >= 3
    if (PyUnicode_CompareWithASCIIString(value, "global") == 0)
       self->mode = Global;
    else if (PyUnicode_CompareWithASCIIString(value, "local") == 0)
       self->mode = Local;
#else
    mode = PyString_AsString(value);
    if (strcmp(mode, "global")==0) self->mode = Global;
    else if (strcmp(mode, "local")==0) self->mode = Local;
#endif
    else {
        PyErr_SetString(PyExc_ValueError, "invalid mode");
        return -1;
    }
    return 0;
}

static char Aligner_match__doc__[] = "match score";

static PyObject*
Aligner_get_match(Aligner* self, void* closure)
{   if (self->letters) {
        PyErr_SetString(PyExc_ValueError, "using a substitution matrix");
        return NULL;
    }
    return PyFloat_FromDouble(self->match);
}

static int
Aligner_set_match(Aligner* self, PyObject* value, void* closure)
{   int i;
    const int n = 26;
    const double match = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_ValueError, "invalid match score");
        return -1;
    }
    self->match = match;
    for (i = 0; i < n; i++) self->substitution_matrix[i][i] = match;
    i = 'X' - 'A';
    self->substitution_matrix[i][i] = 0.0;
    if (self->letters) {
        PyMem_Free(self->letters);
        self->letters = NULL;
    }
    return 0;
}

static char Aligner_mismatch__doc__[] = "mismatch score";

static PyObject*
Aligner_get_mismatch(Aligner* self, void* closure)
{   if (self->letters) {
        PyErr_SetString(PyExc_ValueError, "using a substitution matrix");
        return NULL;
    }
    return PyFloat_FromDouble(self->mismatch);
}

static int
Aligner_set_mismatch(Aligner* self, PyObject* value, void* closure)
{   int i, j;
    const int n = 26;
    const double mismatch = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_ValueError, "invalid match score");
        return -1;
    }
    self->mismatch = mismatch;
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            self->substitution_matrix[i][j] = mismatch;
            self->substitution_matrix[j][i] = mismatch;
        }
    }
    i = 'X' - 'A';
    for (j = 0; j < n; j++) {
        self->substitution_matrix[i][j] = 0;
        self->substitution_matrix[j][i] = 0;
    }
    if (self->letters) {
        PyMem_Free(self->letters);
        self->letters = NULL;
    }
    return 0;
}

static char Aligner_substitution_matrix__doc__[] = "substitution_matrix";

static PyObject*
Aligner_get_substitution_matrix(Aligner* self, void* closure)
{   if (!self->letters) {
        PyErr_SetString(PyExc_ValueError, "using affine gap scores");
        return NULL;
    }
    else {
        int i, j;
        const int n = 26;
        const int* letters = self->letters;
        PyObject* key = NULL;
        PyObject* value = NULL;
        PyObject* matrix = PyDict_New();
        if (!matrix) goto exit;
        for (i = 0; i < n; i++) {
            if (!letters[i]) continue;
            for (j = 0; j < n; j++) {
                if (!letters[j]) continue;
#if PY_MAJOR_VERSION >= 3
                key = Py_BuildValue("(CC)", 'A' + i, 'A' + j);
#else
                key = Py_BuildValue("(cc)", 'A' + i, 'A' + j);
#endif
                if (!key) goto exit;
                value = PyFloat_FromDouble(self->substitution_matrix[i][j]);
                if (!value) goto exit;
                if (PyDict_SetItem(matrix, key, value) == -1) goto exit;
            }
        }
        return matrix;
exit:
        Py_XDECREF(matrix);
        Py_XDECREF(key);
        Py_XDECREF(value);
        return NULL;
    }
}

static int
Aligner_set_substitution_matrix(Aligner* self, PyObject* values, void* closure)
{   int i, j;
    const int n = 26;
    PyObject* key;
    PyObject* value;
    Py_ssize_t pos = 0;
    double score;
    double substitution_matrix[26][26];
    int substitution_matrix_flags[26][26];
    int letters[26];
    PyObject* item;
    if (!PyDict_Check(values)) {
        PyErr_SetString(PyExc_ValueError, "expected a dictionary");
        return -1;
    }
    for (i = 0; i < n; i++) {
        letters[i] = 0;
        for (j = 0; j < n; j++) substitution_matrix_flags[i][j] = 0;
    }
    while (PyDict_Next(values, &pos, &key, &value)) {
        score = PyFloat_AsDouble(value);
        if (PyErr_Occurred()) {
            PyErr_SetString(PyExc_ValueError, "invalid score found");
            return -1;
        }
        if (!PyTuple_Check(key) || PyTuple_GET_SIZE(key) != 2) {
            PyErr_SetString(PyExc_ValueError,
                            "each key should be a tuple of two letters");
            return -1;
        }
        item = PyTuple_GET_ITEM(key, 0);
        i = _convert_single_letter(item);
        if (i < 0) return -1;
        item = PyTuple_GET_ITEM(key, 1);
        j = _convert_single_letter(item);
        if (j < 0) return -1;
        if (substitution_matrix_flags[i][j]) {
            i = 'A' + i;
            j = 'A' + j;
            PyErr_Format(PyExc_ValueError,
                         "score for (%c,%c) specified more than once (substitution matrix is case-insensitive)",
                         i, j);
            return -1;
        }
        substitution_matrix_flags[i][j] = 1;
        substitution_matrix[i][j] = score;
        letters[i] = 1;
        letters[j] = 1;
    }
    if (!self->letters) self->letters = PyMem_Malloc(n*sizeof(int));
    if (!self->letters) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory");
        return -1;
    }
    /* No errors - store the new substitution matrix */
    for (i = 0; i < n; i++) {
        self->letters[i] = letters[i];
        for (j = 0; j < n; j++) {
            if (!letters[i] || !letters[j]) continue;
            if (substitution_matrix_flags[i][j])
                score = substitution_matrix[i][j];
            else if (substitution_matrix_flags[j][i])
                score = substitution_matrix[j][i];
            else
                score = 0;
            self->substitution_matrix[i][j] = score;
        }
    }
    return 0;
}

static char Aligner_gap_score__doc__[] = "gap score";

static PyObject*
Aligner_get_gap_score(Aligner* self, void* closure)
{   
    if (self->target_gap_function || self->query_gap_function) {
        if (self->target_gap_function != self->query_gap_function) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        Py_INCREF(self->target_gap_function);
        return self->target_gap_function;
    }
    else {
        const double score = self->target_open_gap_score;
        if (score != self->target_extend_gap_score
         || score != self->target_left_open_gap_score
         || score != self->target_left_extend_gap_score
         || score != self->target_right_open_gap_score
         || score != self->target_right_extend_gap_score
         || score != self->query_open_gap_score
         || score != self->query_extend_gap_score
         || score != self->query_left_open_gap_score
         || score != self->query_left_extend_gap_score
         || score != self->query_right_open_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_gap_score(Aligner* self, PyObject* value, void* closure)
{   if (PyCallable_Check(value)) {
        Py_XDECREF(self->target_gap_function);
        Py_XDECREF(self->query_gap_function);
        Py_INCREF(value);
        Py_INCREF(value);
        self->target_gap_function = value;
        self->query_gap_function = value;
    }
    else {
        const double score = PyFloat_AsDouble(value);
        if (PyErr_Occurred()) return -1;
        if (self->target_gap_function) {
            Py_DECREF(self->target_gap_function);
            self->target_gap_function = NULL;
        }
        if (self->query_gap_function) {
            Py_DECREF(self->query_gap_function);
            self->query_gap_function = NULL;
        }
        self->target_open_gap_score = score;
        self->target_extend_gap_score = score;
        self->target_left_open_gap_score = score;
        self->target_left_extend_gap_score = score;
        self->target_right_open_gap_score = score;
        self->target_right_extend_gap_score = score;
        self->query_open_gap_score = score;
        self->query_extend_gap_score = score;
        self->query_left_open_gap_score = score;
        self->query_left_extend_gap_score = score;
        self->query_right_open_gap_score = score;
        self->query_right_extend_gap_score = score;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_open_gap_score__doc__[] = "internal and end open gap score";

static PyObject*
Aligner_get_open_gap_score(Aligner* self, void* closure)
{   
    if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_open_gap_score;
        if (score != self->target_left_open_gap_score
         || score != self->target_right_open_gap_score
         || score != self->query_open_gap_score
         || score != self->query_left_open_gap_score
         || score != self->query_right_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_open_gap_score = score;
    self->target_left_open_gap_score = score;
    self->target_right_open_gap_score = score;
    self->query_open_gap_score = score;
    self->query_left_open_gap_score = score;
    self->query_right_open_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_extend_gap_score__doc__[] = "extend gap score";

static PyObject*
Aligner_get_extend_gap_score(Aligner* self, void* closure)
{   
    if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_extend_gap_score;
        if (score != self->target_left_extend_gap_score
         || score != self->target_right_extend_gap_score
         || score != self->query_extend_gap_score
         || score != self->query_left_extend_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_extend_gap_score = score;
    self->target_left_extend_gap_score = score;
    self->target_right_extend_gap_score = score;
    self->query_extend_gap_score = score;
    self->query_left_extend_gap_score = score;
    self->query_right_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_internal_gap_score__doc__[] = "internal gap score";

static PyObject*
Aligner_get_internal_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_open_gap_score;
        if (score != self->target_extend_gap_score
         || score != self->query_open_gap_score
         || score != self->query_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_internal_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_open_gap_score = score;
    self->target_extend_gap_score = score;
    self->query_open_gap_score = score;
    self->query_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_internal_open_gap_score__doc__[] = "internal open gap score";

static PyObject*
Aligner_get_internal_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_open_gap_score;
        if (score != self->query_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_internal_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_open_gap_score = score;
    self->query_open_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_internal_extend_gap_score__doc__[] = "internal extend gap score";

static PyObject*
Aligner_get_internal_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_extend_gap_score;
        if (score != self->query_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_internal_extend_gap_score(Aligner* self, PyObject* value,
                                      void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_extend_gap_score = score;
    self->query_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_end_gap_score__doc__[] = "end gap score";

static PyObject*
Aligner_get_end_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_open_gap_score;
        if (score != self->target_left_extend_gap_score
         || score != self->target_right_open_gap_score
         || score != self->target_right_extend_gap_score
         || score != self->query_left_open_gap_score
         || score != self->query_left_extend_gap_score
         || score != self->query_right_open_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_end_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_left_open_gap_score = score;
    self->target_left_extend_gap_score = score;
    self->target_right_open_gap_score = score;
    self->target_right_extend_gap_score = score;
    self->query_left_open_gap_score = score;
    self->query_left_extend_gap_score = score;
    self->query_right_open_gap_score = score;
    self->query_right_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_end_open_gap_score__doc__[] = "end open gap score";

static PyObject*
Aligner_get_end_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_open_gap_score;
        if (score != self->target_right_open_gap_score
         || score != self->query_left_open_gap_score
         || score != self->query_right_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_end_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_left_open_gap_score = score;
    self->target_right_open_gap_score = score;
    self->query_left_open_gap_score = score;
    self->query_right_open_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_end_extend_gap_score__doc__[] = "end extend gap score";

static PyObject*
Aligner_get_end_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_extend_gap_score;
        if (score != self->target_right_extend_gap_score
         || score != self->query_left_extend_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_end_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_left_extend_gap_score = score;
    self->target_right_extend_gap_score = score;
    self->query_left_extend_gap_score = score;
    self->query_right_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_left_gap_score__doc__[] = "left gap score";

static PyObject*
Aligner_get_left_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_open_gap_score;
        if (score != self->target_left_extend_gap_score
         || score != self->query_left_open_gap_score
         || score != self->query_left_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_left_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_left_open_gap_score = score;
    self->target_left_extend_gap_score = score;
    self->query_left_open_gap_score = score;
    self->query_left_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_right_gap_score__doc__[] = "right gap score";

static PyObject*
Aligner_get_right_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_right_open_gap_score;
        if (score != self->target_right_extend_gap_score
         || score != self->query_right_open_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_right_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_right_open_gap_score = score;
    self->target_right_extend_gap_score = score;
    self->query_right_open_gap_score = score;
    self->query_right_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_left_open_gap_score__doc__[] = "left open gap score";

static PyObject*
Aligner_get_left_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_open_gap_score;
        if (score != self->query_left_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_left_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_left_open_gap_score = score;
    self->query_left_open_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_left_extend_gap_score__doc__[] = "left extend gap score";

static PyObject*
Aligner_get_left_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_extend_gap_score;
        if (score != self->query_left_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_left_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_left_extend_gap_score = score;
    self->query_left_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_right_open_gap_score__doc__[] = "right open gap score";

static PyObject*
Aligner_get_right_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_right_open_gap_score;
        if (score != self->query_right_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_right_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_right_open_gap_score = score;
    self->query_right_open_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_right_extend_gap_score__doc__[] = "right extend gap score";

static PyObject*
Aligner_get_right_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function || self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_right_extend_gap_score;
        if (score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_right_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->target_right_extend_gap_score = score;
    self->query_right_extend_gap_score = score;
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_open_gap_score__doc__[] = "target open gap score";

static PyObject*
Aligner_get_target_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_open_gap_score;
        if (score != self->target_left_open_gap_score
         || score != self->target_right_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_open_gap_score = score;
    self->target_left_open_gap_score = score;
    self->target_right_open_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_extend_gap_score__doc__[] = "target extend gap score";

static PyObject*
Aligner_get_target_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_extend_gap_score;
        if (score != self->target_left_extend_gap_score
         || score != self->target_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_extend_gap_score = score;
    self->target_left_extend_gap_score = score;
    self->target_right_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_gap_score__doc__[] = "target gap score";

static PyObject*
Aligner_get_target_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        Py_INCREF(self->target_gap_function);
        return self->target_gap_function;
    }
    else {
        const double score = self->target_open_gap_score;
        if (score != self->target_extend_gap_score
         || score != self->target_left_open_gap_score
         || score != self->target_left_extend_gap_score
         || score != self->target_right_open_gap_score
         || score != self->target_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_gap_score(Aligner* self, PyObject* value, void* closure)
{
    if (PyCallable_Check(value)) {
        Py_XDECREF(self->target_gap_function);
        Py_INCREF(value);
        self->target_gap_function = value;
    }
    else {
        const double score = PyFloat_AsDouble(value);
        if (PyErr_Occurred()) {
            PyErr_SetString(PyExc_ValueError,
                            "gap score should be numerical or callable");
            return -1;
        }
        self->target_open_gap_score = score;
        self->target_extend_gap_score = score;
        self->target_left_open_gap_score = score;
        self->target_left_extend_gap_score = score;
        self->target_right_open_gap_score = score;
        self->target_right_extend_gap_score = score;
        if (self->target_gap_function) {
            Py_DECREF(self->target_gap_function);
            self->target_gap_function = NULL;
        }
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_open_gap_score__doc__[] = "query gap open score";

static PyObject*
Aligner_get_query_open_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_open_gap_score;
        if (score != self->query_left_open_gap_score
         || score != self->query_right_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_open_gap_score = score;
    self->query_left_open_gap_score = score;
    self->query_right_open_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_extend_gap_score__doc__[] = "query gap extend score";

static PyObject*
Aligner_get_query_extend_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_extend_gap_score;
        if (score != self->query_left_extend_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_extend_gap_score = score;
    self->query_left_extend_gap_score = score;
    self->query_right_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_gap_score__doc__[] = "query gap score";

static PyObject*
Aligner_get_query_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        Py_INCREF(self->query_gap_function);
        return self->query_gap_function;
    }
    else {
        const double score = self->query_open_gap_score;
        if (score != self->query_left_open_gap_score
         || score != self->query_right_open_gap_score
         || score != self->query_extend_gap_score
         || score != self->query_left_extend_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_gap_score(Aligner* self, PyObject* value, void* closure)
{   if (PyCallable_Check(value)) {
        Py_XDECREF(self->query_gap_function);
        Py_INCREF(value);
        self->query_gap_function = value;
    }
    else {
        const double score = PyFloat_AsDouble(value);
        if (PyErr_Occurred()) {
            PyErr_SetString(PyExc_ValueError,
                            "gap score should be numerical or callable");
            return -1;
        }
        self->query_open_gap_score = score;
        self->query_extend_gap_score = score;
        self->query_left_open_gap_score = score;
        self->query_left_extend_gap_score = score;
        self->query_right_open_gap_score = score;
        self->query_right_extend_gap_score = score;
        if (self->query_gap_function) {
            Py_DECREF(self->query_gap_function);
            self->query_gap_function = NULL;
        }
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_internal_open_gap_score__doc__[] = "target internal open gap score";

static PyObject*
Aligner_get_target_internal_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->target_open_gap_score);
}

static int
Aligner_set_target_internal_open_gap_score(Aligner* self,
                                           PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_open_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_internal_extend_gap_score__doc__[] = "target internal extend gap score";

static PyObject*
Aligner_get_target_internal_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->target_extend_gap_score);
}

static int
Aligner_set_target_internal_extend_gap_score(Aligner* self,
                                             PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_internal_gap_score__doc__[] = "target internal gap score";

static PyObject*
Aligner_get_target_internal_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_open_gap_score;
        if (score != self->target_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_internal_gap_score(Aligner* self, PyObject* value,
                                      void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_open_gap_score = score;
    self->target_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_end_gap_score__doc__[] = "target end gap score";

static PyObject*
Aligner_get_target_end_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_open_gap_score;
        if (score != self->target_left_extend_gap_score
         || score != self->target_right_open_gap_score
         || score != self->target_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_end_gap_score(Aligner* self, PyObject* value, void* closure) {
    const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_left_open_gap_score = score;
    self->target_left_extend_gap_score = score;
    self->target_right_open_gap_score = score;
    self->target_right_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_end_open_gap_score__doc__[] = "target end open gap score";

static PyObject*
Aligner_get_target_end_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_open_gap_score;
        if (score != self->target_right_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_end_open_gap_score(Aligner* self, PyObject* value,
                                      void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_left_open_gap_score = score;
    self->target_right_open_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_end_extend_gap_score__doc__[] = "target end extend gap score";

static PyObject*
Aligner_get_target_end_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_extend_gap_score;
        if (score != self->target_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_end_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_left_extend_gap_score = score;
    self->target_right_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_left_open_gap_score__doc__[] = "target left open score";

static PyObject*
Aligner_get_target_left_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->target_left_open_gap_score);
}

static int
Aligner_set_target_left_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_left_open_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_left_extend_gap_score__doc__[] = "target left extend score";

static PyObject*
Aligner_get_target_left_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->target_left_extend_gap_score);
}

static int
Aligner_set_target_left_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_left_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_left_gap_score__doc__[] = "target left score";

static PyObject*
Aligner_get_target_left_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_left_open_gap_score;
        if (score != self->target_left_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_left_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_left_open_gap_score = score;
    self->target_left_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_right_gap_score_open__doc__[] = "target right open score";

static PyObject*
Aligner_get_target_right_open_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->target_right_open_gap_score);
}

static int
Aligner_set_target_right_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_right_open_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_right_extend_gap_score__doc__[] = "target right extend score";

static PyObject*
Aligner_get_target_right_extend_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->target_right_extend_gap_score);
}

static int
Aligner_set_target_right_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_right_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_target_right_gap_score__doc__[] = "target right score";

static PyObject*
Aligner_get_target_right_gap_score(Aligner* self, void* closure)
{   if (self->target_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->target_right_open_gap_score;
        if (score != self->target_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_target_right_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->target_right_open_gap_score = score;
    self->target_right_extend_gap_score = score;
    if (self->target_gap_function) {
        Py_DECREF(self->target_gap_function);
        self->target_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_end_gap_score__doc__[] = "query end score";

static PyObject*
Aligner_get_query_end_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_left_open_gap_score;
        if (score != self->query_left_extend_gap_score
         || score != self->query_right_open_gap_score
         || score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_end_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_left_open_gap_score = score;
    self->query_left_extend_gap_score = score;
    self->query_right_open_gap_score = score;
    self->query_right_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_end_open_gap_score__doc__[] = "query end open score";

static PyObject*
Aligner_get_query_end_open_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_left_open_gap_score;
        if (score != self->query_right_open_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_end_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_left_open_gap_score = score;
    self->query_right_open_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_end_extend_gap_score__doc__[] = "query end extend score";

static PyObject*
Aligner_get_query_end_extend_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_left_extend_gap_score;
        if (score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_end_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_left_extend_gap_score = score;
    self->query_right_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_internal_open_gap_score__doc__[] = "query internal open gap score";

static PyObject*
Aligner_get_query_internal_open_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->query_open_gap_score);
}

static int
Aligner_set_query_internal_open_gap_score(Aligner* self, PyObject* value,
                                          void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_open_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_internal_extend_gap_score__doc__[] = "query internal extend gap score";

static PyObject*
Aligner_get_query_internal_extend_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->query_extend_gap_score);
}

static int
Aligner_set_query_internal_extend_gap_score(Aligner* self, PyObject* value,
                                            void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_internal_gap_score__doc__[] = "query internal gap score";

static PyObject*
Aligner_get_query_internal_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_open_gap_score;
        if (score != self->query_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_internal_gap_score(Aligner* self, PyObject* value,
                                     void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_open_gap_score = score;
    self->query_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_left_open_gap_score__doc__[] = "query left open score";

static PyObject*
Aligner_get_query_left_open_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->query_left_open_gap_score);
}

static int
Aligner_set_query_left_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_left_open_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_left_extend_gap_score__doc__[] = "query left extend score";

static PyObject*
Aligner_get_query_left_extend_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->query_left_extend_gap_score);
}

static int
Aligner_set_query_left_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_left_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_left_gap_score__doc__[] = "query left score";

static PyObject*
Aligner_get_query_left_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_left_open_gap_score;
        if (score != self->query_left_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_left_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_left_open_gap_score = score;
    self->query_left_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_right_open_gap_score__doc__[] = "query right open score";

static PyObject*
Aligner_get_query_right_open_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->query_right_open_gap_score);
}

static int
Aligner_set_query_right_open_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_right_open_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_right_extend_gap_score__doc__[] = "query right extend score";

static PyObject*
Aligner_get_query_right_extend_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    return PyFloat_FromDouble(self->query_right_extend_gap_score);
}

static int
Aligner_set_query_right_extend_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_right_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_query_right_gap_score__doc__[] = "query right score";

static PyObject*
Aligner_get_query_right_gap_score(Aligner* self, void* closure)
{   if (self->query_gap_function) {
        PyErr_SetString(PyExc_ValueError, "using a gap score function");
        return NULL;
    }
    else {
        const double score = self->query_right_open_gap_score;
        if (score != self->query_right_extend_gap_score) {
            PyErr_SetString(PyExc_ValueError, "gap scores are different");
            return NULL;
        }
        return PyFloat_FromDouble(score);
    }
}

static int
Aligner_set_query_right_gap_score(Aligner* self, PyObject* value, void* closure)
{   const double score = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->query_right_open_gap_score = score;
    self->query_right_extend_gap_score = score;
    if (self->query_gap_function) {
        Py_DECREF(self->query_gap_function);
        self->query_gap_function = NULL;
    }
    self->algorithm = Unknown;
    return 0;
}

static char Aligner_epsilon__doc__[] = "roundoff epsilon";

static PyObject*
Aligner_get_epsilon(Aligner* self, void* closure)
{   return PyFloat_FromDouble(self->epsilon);
}

static int
Aligner_set_epsilon(Aligner* self, PyObject* value, void* closure)
{   const double epsilon = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) return -1;
    self->epsilon = epsilon;
    self->algorithm = Unknown;
    return 0;
}

static Algorithm _get_algorithm(Aligner* self)
{
    Algorithm algorithm = self->algorithm;
    if (algorithm == Unknown) {
        const double target_gap_open = self->target_open_gap_score;
        const double query_gap_open = self->query_open_gap_score;
        const double target_gap_extend = self->target_extend_gap_score;
        const double query_gap_extend = self->query_extend_gap_score;
        const double target_left_open = self->target_left_open_gap_score;
        const double target_left_extend = self->target_left_extend_gap_score;
        const double query_left_open = self->query_left_open_gap_score;
        const double target_right_open = self->target_right_open_gap_score;
        const double query_right_open = self->query_right_open_gap_score;
        const double target_right_extend = self->target_right_extend_gap_score;
        const double query_left_extend = self->query_left_extend_gap_score;
        const double query_right_extend = self->query_right_extend_gap_score;
        if (self->target_gap_function || self->query_gap_function)
            algorithm = WatermanSmithBeyer;
        else if (target_gap_open == target_gap_extend
              && query_gap_open == query_gap_extend
              && target_left_open == target_left_extend
              && target_right_open == target_right_extend
              && query_left_open == query_left_extend
              && query_right_open == query_right_extend)
            algorithm = NeedlemanWunschSmithWaterman;
        else
            algorithm = Gotoh;
        self->algorithm = algorithm;
    }
    return algorithm;
}


static char Aligner_algorithm__doc__[] = "alignment algorithm";

static PyObject*
Aligner_get_algorithm(Aligner* self, void* closure)
{
    const char* s = NULL;
    const Mode mode = self->mode;
    const Algorithm algorithm = _get_algorithm(self);
    switch (algorithm) {
        case NeedlemanWunschSmithWaterman:
            switch (mode) {
                case Global:
                    s = "Needleman-Wunsch";
                    break;
                case Local:
                    s = "Smith-Waterman";
                    break;
            }
            break;
        case Gotoh:
            switch (mode) {
                case Global:
                    s = "Gotoh global alignment algorithm";
                    break;
                case Local:
                    s = "Gotoh local alignment algorithm";
                    break;
            }
            break;
        case WatermanSmithBeyer:
            switch (mode) {
                case Global:
                    s = "Waterman-Smith-Beyer global alignment algorithm";
                    break;
                case Local:
                    s = "Waterman-Smith-Beyer local alignment algorithm";
                    break;
            }
            break;
        case Unknown:
        default:
            break;
    }
#if PY_MAJOR_VERSION >= 3
    return PyUnicode_FromString(s);
#else
    return PyString_FromString(s);
#endif
}

static PyGetSetDef Aligner_getset[] = {
    {"mode",
        (getter)Aligner_get_mode,
        (setter)Aligner_set_mode,
        Aligner_mode__doc__, NULL},
    {"match",
        (getter)Aligner_get_match,
        (setter)Aligner_set_match,
        Aligner_match__doc__, NULL},
    {"mismatch",
        (getter)Aligner_get_mismatch,
        (setter)Aligner_set_mismatch,
        Aligner_mismatch__doc__, NULL},
    {"substitution_matrix",
        (getter)Aligner_get_substitution_matrix,
        (setter)Aligner_set_substitution_matrix,
        Aligner_substitution_matrix__doc__, NULL},
    {"gap_score",
        (getter)Aligner_get_gap_score,
        (setter)Aligner_set_gap_score,
        Aligner_gap_score__doc__, NULL},
    {"open_gap_score",
        (getter)Aligner_get_open_gap_score,
        (setter)Aligner_set_open_gap_score,
        Aligner_open_gap_score__doc__, NULL},
    {"extend_gap_score",
        (getter)Aligner_get_extend_gap_score,
        (setter)Aligner_set_extend_gap_score,
        Aligner_extend_gap_score__doc__, NULL},
    {"internal_gap_score",
        (getter)Aligner_get_internal_gap_score,
        (setter)Aligner_set_internal_gap_score,
        Aligner_internal_gap_score__doc__, NULL},
    {"internal_open_gap_score",
        (getter)Aligner_get_internal_open_gap_score,
        (setter)Aligner_set_internal_open_gap_score,
        Aligner_internal_open_gap_score__doc__, NULL},
    {"internal_extend_gap_score",
        (getter)Aligner_get_internal_extend_gap_score,
        (setter)Aligner_set_internal_extend_gap_score,
        Aligner_internal_extend_gap_score__doc__, NULL},
    {"end_gap_score",
        (getter)Aligner_get_end_gap_score,
        (setter)Aligner_set_end_gap_score,
        Aligner_end_gap_score__doc__, NULL},
    {"end_open_gap_score",
        (getter)Aligner_get_end_open_gap_score,
        (setter)Aligner_set_end_open_gap_score,
        Aligner_end_open_gap_score__doc__, NULL},
    {"end_extend_gap_score",
        (getter)Aligner_get_end_extend_gap_score,
        (setter)Aligner_set_end_extend_gap_score,
        Aligner_end_extend_gap_score__doc__, NULL},
    {"left_gap_score",
        (getter)Aligner_get_left_gap_score,
        (setter)Aligner_set_left_gap_score,
        Aligner_left_gap_score__doc__, NULL},
    {"left_open_gap_score",
        (getter)Aligner_get_left_open_gap_score,
        (setter)Aligner_set_left_open_gap_score,
        Aligner_left_open_gap_score__doc__, NULL},
    {"left_extend_gap_score",
        (getter)Aligner_get_left_extend_gap_score,
        (setter)Aligner_set_left_extend_gap_score,
        Aligner_left_extend_gap_score__doc__, NULL},
    {"right_gap_score",
        (getter)Aligner_get_right_gap_score,
        (setter)Aligner_set_right_gap_score,
        Aligner_right_gap_score__doc__, NULL},
    {"right_open_gap_score",
        (getter)Aligner_get_right_open_gap_score,
        (setter)Aligner_set_right_open_gap_score,
        Aligner_right_open_gap_score__doc__, NULL},
    {"right_extend_gap_score",
        (getter)Aligner_get_right_extend_gap_score,
        (setter)Aligner_set_right_extend_gap_score,
        Aligner_right_extend_gap_score__doc__, NULL},
    {"target_open_gap_score",
        (getter)Aligner_get_target_open_gap_score,
        (setter)Aligner_set_target_open_gap_score,
        Aligner_target_open_gap_score__doc__, NULL},
    {"target_extend_gap_score",
        (getter)Aligner_get_target_extend_gap_score,
        (setter)Aligner_set_target_extend_gap_score,
        Aligner_target_extend_gap_score__doc__, NULL},
    {"target_gap_score",
        (getter)Aligner_get_target_gap_score,
        (setter)Aligner_set_target_gap_score,
        Aligner_target_gap_score__doc__, NULL},
    {"query_open_gap_score",
        (getter)Aligner_get_query_open_gap_score,
        (setter)Aligner_set_query_open_gap_score,
        Aligner_query_open_gap_score__doc__, NULL},
    {"query_extend_gap_score",
        (getter)Aligner_get_query_extend_gap_score,
        (setter)Aligner_set_query_extend_gap_score,
        Aligner_query_extend_gap_score__doc__, NULL},
    {"query_gap_score",
        (getter)Aligner_get_query_gap_score,
        (setter)Aligner_set_query_gap_score,
        Aligner_query_gap_score__doc__, NULL},
    {"target_end_gap_score",
        (getter)Aligner_get_target_end_gap_score,
        (setter)Aligner_set_target_end_gap_score,
        Aligner_target_end_gap_score__doc__, NULL},
    {"target_end_open_gap_score",
        (getter)Aligner_get_target_end_open_gap_score,
        (setter)Aligner_set_target_end_open_gap_score,
        Aligner_target_end_open_gap_score__doc__, NULL},
    {"target_end_extend_gap_score",
        (getter)Aligner_get_target_end_extend_gap_score,
        (setter)Aligner_set_target_end_extend_gap_score,
        Aligner_target_end_extend_gap_score__doc__, NULL},
    {"target_internal_open_gap_score",
        (getter)Aligner_get_target_internal_open_gap_score,
        (setter)Aligner_set_target_internal_open_gap_score,
        Aligner_target_internal_open_gap_score__doc__, NULL},
    {"target_internal_extend_gap_score",
        (getter)Aligner_get_target_internal_extend_gap_score,
        (setter)Aligner_set_target_internal_extend_gap_score,
        Aligner_target_internal_extend_gap_score__doc__, NULL},
    {"target_internal_gap_score",
        (getter)Aligner_get_target_internal_gap_score,
        (setter)Aligner_set_target_internal_gap_score,
        Aligner_target_internal_gap_score__doc__, NULL},
    {"target_left_open_gap_score",
        (getter)Aligner_get_target_left_open_gap_score,
        (setter)Aligner_set_target_left_open_gap_score,
        Aligner_target_left_open_gap_score__doc__, NULL},
    {"target_left_extend_gap_score",
        (getter)Aligner_get_target_left_extend_gap_score,
        (setter)Aligner_set_target_left_extend_gap_score,
        Aligner_target_left_extend_gap_score__doc__, NULL},
    {"target_left_gap_score",
        (getter)Aligner_get_target_left_gap_score,
        (setter)Aligner_set_target_left_gap_score,
        Aligner_target_left_gap_score__doc__, NULL},
    {"target_right_open_gap_score",
        (getter)Aligner_get_target_right_open_gap_score,
        (setter)Aligner_set_target_right_open_gap_score,
        Aligner_target_right_gap_score_open__doc__, NULL},
    {"target_right_extend_gap_score",
        (getter)Aligner_get_target_right_extend_gap_score,
        (setter)Aligner_set_target_right_extend_gap_score,
        Aligner_target_right_extend_gap_score__doc__, NULL},
    {"target_right_gap_score",
        (getter)Aligner_get_target_right_gap_score,
        (setter)Aligner_set_target_right_gap_score,
        Aligner_target_right_gap_score__doc__, NULL},
    {"query_end_gap_score",
        (getter)Aligner_get_query_end_gap_score,
        (setter)Aligner_set_query_end_gap_score,
        Aligner_query_end_gap_score__doc__, NULL},
    {"query_end_open_gap_score",
        (getter)Aligner_get_query_end_open_gap_score,
        (setter)Aligner_set_query_end_open_gap_score,
        Aligner_query_end_open_gap_score__doc__, NULL},
    {"query_end_extend_gap_score",
        (getter)Aligner_get_query_end_extend_gap_score,
        (setter)Aligner_set_query_end_extend_gap_score,
        Aligner_query_end_extend_gap_score__doc__, NULL},
    {"query_internal_open_gap_score",
        (getter)Aligner_get_query_internal_open_gap_score,
        (setter)Aligner_set_query_internal_open_gap_score,
        Aligner_query_internal_open_gap_score__doc__, NULL},
    {"query_internal_extend_gap_score",
        (getter)Aligner_get_query_internal_extend_gap_score,
        (setter)Aligner_set_query_internal_extend_gap_score,
        Aligner_query_internal_extend_gap_score__doc__, NULL},
    {"query_internal_gap_score",
        (getter)Aligner_get_query_internal_gap_score,
        (setter)Aligner_set_query_internal_gap_score,
        Aligner_query_internal_gap_score__doc__, NULL},
    {"query_left_open_gap_score",
        (getter)Aligner_get_query_left_open_gap_score,
        (setter)Aligner_set_query_left_open_gap_score,
        Aligner_query_left_open_gap_score__doc__, NULL},
    {"query_left_extend_gap_score",
        (getter)Aligner_get_query_left_extend_gap_score,
        (setter)Aligner_set_query_left_extend_gap_score,
        Aligner_query_left_extend_gap_score__doc__, NULL},
    {"query_left_gap_score",
        (getter)Aligner_get_query_left_gap_score,
        (setter)Aligner_set_query_left_gap_score,
         Aligner_query_left_gap_score__doc__, NULL},
    {"query_right_open_gap_score",
        (getter)Aligner_get_query_right_open_gap_score,
        (setter)Aligner_set_query_right_open_gap_score,
        Aligner_query_right_open_gap_score__doc__, NULL},
    {"query_right_extend_gap_score",
        (getter)Aligner_get_query_right_extend_gap_score,
        (setter)Aligner_set_query_right_extend_gap_score,
        Aligner_query_right_extend_gap_score__doc__, NULL},
    {"query_right_gap_score",
        (getter)Aligner_get_query_right_gap_score,
        (setter)Aligner_set_query_right_gap_score,
        Aligner_query_right_gap_score__doc__, NULL},
    {"epsilon",
        (getter)Aligner_get_epsilon,
        (setter)Aligner_set_epsilon,
        Aligner_epsilon__doc__, NULL},
    {"algorithm",
        (getter)Aligner_get_algorithm,
        (setter)NULL,
        Aligner_algorithm__doc__, NULL},
    {NULL}  /* Sentinel */
};

#define SELECT_SCORE_GLOBAL(score1, score2, score3) \
    score = score1; \
    temp = score2; \
    if (temp > score) score = temp; \
    temp = score3; \
    if (temp > score) score = temp;

#define SELECT_SCORE_WATERMAN_SMITH_BEYER(score1, score2) \
    temp = score1 + gapscore; \
    if (temp > score) score = temp; \
    temp = score2 + gapscore; \
    if (temp > score) score = temp;

#define SELECT_SCORE_GOTOH_LOCAL_ALIGN(score1, score2, score3, score4) \
    score = score1; \
    temp = score2; \
    if (temp > score) score = temp; \
    temp = score3; \
    if (temp > score) score = temp; \
    score += score4; \
    if (score < 0) score = 0; \
    else if (score > maximum) maximum = score;

#define SELECT_SCORE_LOCAL3(score1, score2, score3) \
    score = score1; \
    temp = score2; \
    if (temp > score) score = temp; \
    temp = score3; \
    if (temp > score) score = temp; \
    if (score < 0) score = 0; \
    else if (score > maximum) maximum = score;

#define SELECT_SCORE_LOCAL1(score1) \
    score = score1; \
    if (score < 0) score = 0; \
    else if (score > maximum) maximum = score;

#define SELECT_TRACE_NEEDLEMAN_WUNSCH(hgap, vgap) \
    score = temp + self->substitution_matrix[kA][kB]; \
    trace = DIAGONAL; \
    temp = scores[j-1] + hgap; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = HORIZONTAL; \
    } \
    else if (temp > score - epsilon) trace |= HORIZONTAL; \
    temp = scores[j] + vgap; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = VERTICAL; \
    } \
    else if (temp > score - epsilon) trace |= VERTICAL; \
    temp = scores[j]; \
    scores[j] = score; \
    M[i][j].trace = trace;

#define SELECT_TRACE_SMITH_WATERMAN_HVD(hgap, vgap) \
    trace = DIAGONAL; \
    score = temp + self->substitution_matrix[kA][kB]; \
    temp = scores[j-1] + hgap; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = HORIZONTAL; \
    } \
    else if (temp > score - epsilon) trace |= HORIZONTAL; \
    temp = scores[j] + vgap; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = VERTICAL; \
    } \
    else if (temp > score - epsilon) trace |= VERTICAL; \
    if (score < epsilon) { \
        score = 0; \
        trace = 0; \
    } \
    else if (trace & DIAGONAL && score > maximum - epsilon) { \
        if (score > maximum + epsilon) { \
            for ( ; im < i; im++, jm = 0) \
                for ( ; jm <= nB; jm++) M[im][jm].trace &= ~ENDPOINT; \
            for ( ; jm < j; jm++) M[im][jm].trace &= ~ENDPOINT; \
            im = i; \
            jm = j; \
        } \
        trace |= ENDPOINT; \
    } \
    M[i][j].trace = trace; \
    if (score > maximum) maximum = score; \
    temp = scores[j]; \
    scores[j] = score;

#define SELECT_TRACE_SMITH_WATERMAN_D \
    score = temp + self->substitution_matrix[kA][kB]; \
    trace = DIAGONAL; \
    if (score < epsilon) { \
        score = 0; \
        trace = 0; \
    } \
    else if (trace & DIAGONAL && score > maximum - epsilon) { \
        if (score > maximum + epsilon) { \
            for ( ; im < i; im++, jm = 0) \
                for ( ; jm <= nB; jm++) M[im][jm].trace &= ~ENDPOINT; \
            for ( ; jm < j; jm++) M[im][jm].trace &= ~ENDPOINT; \
            im = i; \
            jm = j; \
        } \
        trace |= ENDPOINT; \
    } \
    M[i][j].trace = trace; \
    if (score > maximum) maximum = score; \
    temp = scores[j]; \
    scores[j] = score

#define SELECT_TRACE_GOTOH_GLOBAL_GAP(matrix, score1, score2, score3) \
    trace = M_MATRIX; \
    score = score1; \
    temp = score2; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Ix_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Ix_MATRIX; \
    temp = score3; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Iy_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Iy_MATRIX; \
    trace3[i][j].matrix = trace;

#define SELECT_TRACE_GOTOH_GLOBAL_ALIGN \
    trace = M_MATRIX; \
    score = M_temp; \
    temp = Ix_temp; \
    if (temp > score + epsilon) { \
        score = Ix_temp; \
        trace = Ix_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Ix_MATRIX; \
    temp = Iy_temp; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Iy_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Iy_MATRIX; \
    trace3[i][j].M = trace;

#define SELECT_TRACE_GOTOH_LOCAL_ALIGN \
    trace = M_MATRIX; \
    score = M_temp; \
    if (Ix_temp > score + epsilon) { \
        score = Ix_temp; \
        trace = Ix_MATRIX; \
    } \
    else if (Ix_temp > score - epsilon) trace |= Ix_MATRIX; \
    if (Iy_temp > score + epsilon) { \
        score = Iy_temp; \
        trace = Iy_MATRIX; \
    } \
    else if (Iy_temp > score - epsilon) trace |= Iy_MATRIX; \
    score += self->substitution_matrix[kA][kB]; \
    if (score < epsilon) { \
        score = 0; \
        trace = 0; \
    } \
    else if (score > maximum - epsilon) { \
        if (score > maximum + epsilon) { \
            maximum = score; \
            for ( ; im < i; im++, jm = 0) \
                for ( ; jm <= nB; jm++) trace3[im][jm].M &= ~ENDPOINT; \
            for ( ; jm < j; jm++) trace3[im][jm].M &= ~ENDPOINT; \
            im = i; \
            jm = j; \
        } \
        trace |= ENDPOINT; \
    } \
    trace3[i][j].M = trace;

#define SELECT_TRACE_GOTOH_LOCAL_GAP(matrix, score1, score2, score3) \
    trace = M_MATRIX; \
    score = score1; \
    temp = score2; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Ix_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Ix_MATRIX; \
    temp = score3; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Iy_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Iy_MATRIX; \
    if (score < epsilon) { \
        score = -DBL_MAX; \
        trace = 0; \
    } \
    trace3[i][j].matrix = trace;

#define SELECT_TRACE_WATERMAN_SMITH_BEYER_GLOBAL_ALIGN(score4) \
    trace = M_MATRIX; \
    score = M[i-1][j-1].score; \
    temp = Ix[i-1][j-1].score; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Ix_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Ix_MATRIX; \
    temp = Iy[i-1][j-1].score; \
    if (temp > score + epsilon) { \
        score = temp; \
        trace = Iy_MATRIX; \
    } \
    else if (temp > score - epsilon) trace |= Iy_MATRIX; \
    M[i][j].score = score + score4; \
    M[i][j].trace = trace;

#define SELECT_TRACE_WATERMAN_SMITH_BEYER_GAP(score1, score2, gap) \
    temp = score1 + gapscore; \
    if (temp > score - epsilon) { \
        if (temp > score + epsilon) { \
            score = temp; \
            nm = 0; \
            ng = 0; \
        } \
        traceM[nm] = gap; \
        nm++; \
    } \
    temp = score2 + gapscore; \
    if (temp > score - epsilon) { \
        if (temp > score + epsilon) { \
            score = temp; \
            nm = 0; \
            ng = 0; \
        } \
        traceXY[ng] = gap; \
        ng++; \
    }

#define SELECT_TRACE_WATERMAN_SMITH_BEYER_ALIGN(cell, score1, score2, score3, score4) \
    trace = M_MATRIX; \
    score = score1; \
    if (score2 > score + epsilon) { \
        score = score2; \
        trace = Ix_MATRIX; \
    } \
    else if (score2 > score - epsilon) trace |= Ix_MATRIX; \
    if (score3 > score + epsilon) { \
        score = score3; \
        trace = Iy_MATRIX; \
    } \
    else if (score3 > score - epsilon) trace |= Iy_MATRIX; \
    score += score4; \
    if (score < epsilon) { \
        score = 0; \
        trace = 0; \
    } \
    else if (score > maximum - epsilon) { \
        if (score > maximum + epsilon) { \
            maximum = score; \
            for ( ; im < i; im++, jm = 0) \
                for ( ; jm <= nB; jm++) M[im][jm].trace &= ~ENDPOINT; \
            for ( ; jm < j; jm++) M[im][jm].trace &= ~ENDPOINT; \
            im = i; \
            jm = j; \
        } \
        trace |= ENDPOINT; \
    } \
    cell.score = score; \
    cell.trace = trace;

/* -------------- allocation & deallocation ------------- */

static void
_deallocate_watermansmithbeyer_matrices(Py_ssize_t nA, Py_ssize_t nB,
                                        CellM** M, CellXY** Ix, CellXY** Iy,
                                        Mode mode)
{
    int i, j;
    int* trace;
    if (M) {
        if (Ix) {
            if (Iy) {
                for (i = 0; i <= nA; i++) {
                    if (!M[i]) break;
                    PyMem_Free(M[i]);
                    if (Ix[i]) {
                        if (Iy[i]) {
                            for (j = 0; j <= nB; j++) {
                                trace = Ix[i][j].traceM;
                                if (trace) PyMem_Free(trace);
                                trace = Ix[i][j].traceXY;
                                if (trace) PyMem_Free(trace);
                                trace = Iy[i][j].traceM;
                                if (trace) PyMem_Free(trace);
                                trace = Iy[i][j].traceXY;
                                if (trace) PyMem_Free(trace);
                            }
                            PyMem_Free(Iy[i]);
                        } else break;
                        PyMem_Free(Ix[i]);
                    } else break;
                }
                PyMem_Free(Iy);
            }
            PyMem_Free(Ix);
        }
        PyMem_Free(M);
    }
}

static int
_allocate_watermansmithbeyer_matrices(Py_ssize_t nA, Py_ssize_t nB,
                                      CellM*** pM, CellXY*** pIx, CellXY*** pIy,
                                      Mode mode)
{
    int i, j;
    int* trace;
    CellM** M = NULL;
    CellXY** Ix = NULL;
    CellXY** Iy = NULL;
    M = PyMem_Malloc((nA+1)*sizeof(CellM*));
    if (!M) goto exit;
    Ix = PyMem_Malloc((nA+1)*sizeof(CellXY*));
    if (!Ix) goto exit;
    Iy = PyMem_Malloc((nA+1)*sizeof(CellXY*));
    if (!Iy) goto exit;
    for (i = 0; i <= nA; i++) {
        M[i] = PyMem_Malloc((nB+1)*sizeof(CellM));
        if (!M[i]) goto exit;
        Ix[i] = PyMem_Malloc((nB+1)*sizeof(CellXY));
        if (!Ix[i]) goto exit;
        Iy[i] = PyMem_Malloc((nB+1)*sizeof(CellXY));
        if (!Iy[i]) goto exit;
        for (j = 0; j <= nB; j++) {
            Ix[i][j].traceM = NULL;
            Ix[i][j].traceXY = NULL;
            Iy[i][j].traceM = NULL;
            Iy[i][j].traceXY = NULL;
        }
        M[i][0].trace = 0;
        M[i][0].path = 0;
        Ix[i][0].path = 0;
        Iy[i][0].path = 0;
        Iy[i][0].score = -DBL_MAX;
        if (i==0) {
            M[0][0].score = 0;
            Ix[0][0].score = -DBL_MAX;
        }
        else {
            switch (mode) {
                case Global:
                    M[i][0].score = -DBL_MAX;
                    Ix[i][0].score = 0;
                    trace = PyMem_Malloc(2*sizeof(int));
                    if (!trace) goto exit;
                    Ix[i][0].traceM = trace;
                    trace[0] = 0;
                    trace[1] = -1;
                    trace = PyMem_Malloc(sizeof(int));
                    if (!trace) goto exit;
                    Ix[i][0].traceXY = trace;
                    trace[0] = -1;
                    break;
                case Local:
                    M[i][0].score = 0;
                    Ix[i][0].score = -DBL_MAX;
                    Ix[i][0].traceM = NULL;
                    Ix[i][0].traceXY = NULL;
                    break;
            }
        }
    }
    for (i = 1; i <= nB; i++) {
        M[0][i].trace = 0;
        Ix[0][i].traceM = NULL;
        Ix[0][i].traceXY = NULL;
        Ix[0][i].score = -DBL_MAX;
        switch (mode) {
            case Global:
                M[0][i].score = -DBL_MAX;
                Iy[0][i].score = 0;
                trace = PyMem_Malloc(2*sizeof(int));
                if (!trace) goto exit;
                Iy[0][i].traceM = trace;
                trace[0] = 0;
                trace[1] = -1;
                trace = PyMem_Malloc(sizeof(int));
                if (!trace) goto exit;
                Iy[0][i].traceXY = trace;
                trace[0] = -1;
                break;
            case Local:
                M[0][i].score = 0;
                Iy[0][i].score = -DBL_MAX;
                Iy[0][i].traceM = NULL;
                Iy[0][i].traceXY = NULL;
                break;
        }
    }
    M[0][0].path = 0;
    *pM = M;
    *pIx = Ix;
    *pIy = Iy;
    return 1;
exit:
    _deallocate_watermansmithbeyer_matrices(nA, nB, M, Ix, Iy, mode);
    PyErr_SetString(PyExc_MemoryError, "Out of memory");
    return 0;
}

/* ----------------- alignment algorithms ----------------- */

static PyObject*
Aligner_needlemanwunsch_score(Aligner* self, const char* sA, Py_ssize_t nA,
                                             const char* sB, Py_ssize_t nB)
{
    char c;
    int i;
    int j;
    int kA;
    int kB;
    const double gap_extend_A = self->target_extend_gap_score;
    const double gap_extend_B = self->query_extend_gap_score;
    const double left_gap_extend_A = self->target_left_extend_gap_score;
    const double right_gap_extend_A = self->target_right_extend_gap_score;
    const double left_gap_extend_B = self->query_left_extend_gap_score;
    const double right_gap_extend_B = self->query_right_extend_gap_score;
    double score;
    double temp;
    PyObject* result = NULL;

    double* scores;

    /* Needleman-Wunsch algorithm */
    scores = PyMem_Malloc((nB+1)*sizeof(double));
    if (!scores) goto exit;

    /* The top row of the score matrix is a special case,
     * as there are no previously aligned characters.
     */
    scores[0] = 0.0;
    for (j = 1; j <= nB; j++) scores[j] = j * left_gap_extend_A;
    for (i = 1; i < nA; i++) {
        kA = CHARINDEX(sA[i-1]);
        temp = scores[0];
        scores[0] = i * left_gap_extend_B;
        for (j = 1; j < nB; j++) {
            kB = CHARINDEX(sB[j-1]);
            SELECT_SCORE_GLOBAL(temp + self->substitution_matrix[kA][kB],
                                scores[j] + gap_extend_B,
                                scores[j-1] + gap_extend_A);
            temp = scores[j];
            scores[j] = score;
        }
        kB = CHARINDEX(sB[nB-1]);
        SELECT_SCORE_GLOBAL(temp + self->substitution_matrix[kA][kB],
                            scores[nB] + right_gap_extend_B,
                            scores[nB-1] + gap_extend_A);
        temp = scores[nB];
        scores[nB] = score;
    }
    kA = CHARINDEX(sA[nA-1]);
    temp = scores[0];
    scores[0] = nA * right_gap_extend_B;
    for (j = 1; j < nB; j++) {
        kB = CHARINDEX(sB[j-1]);
        SELECT_SCORE_GLOBAL(temp + self->substitution_matrix[kA][kB],
                            scores[j] + gap_extend_B,
                            scores[j-1] + right_gap_extend_A);
        temp = scores[j];
        scores[j] = score;
    }
    kB = CHARINDEX(sB[nB-1]);
    SELECT_SCORE_GLOBAL(temp + self->substitution_matrix[kA][kB],
                        scores[nB] + right_gap_extend_B,
                        scores[nB-1] + right_gap_extend_A);
    result = PyFloat_FromDouble(score);
exit:
    if (scores) PyMem_Free(scores);
    if (!result) PyErr_SetString(PyExc_MemoryError, "Out of memory");
    return result;
}

static PyObject*
Aligner_smithwaterman_score(Aligner* self, const char* sA, Py_ssize_t nA,
                                           const char* sB, Py_ssize_t nB)
{
    char c;
    int i;
    int j;
    int kA;
    int kB;
    const double gap_extend_A = self->target_extend_gap_score;
    const double gap_extend_B = self->query_extend_gap_score;
    double score;
    double* scores;
    double temp;
    double maximum = 0;
    PyObject* result = NULL;

    /* Smith-Waterman algorithm */
    scores = PyMem_Malloc((nB+1)*sizeof(double));
    if (!scores) goto exit;

    /* The top row of the score matrix is a special case,
     * as there are no previously aligned characters.
     */
    for (j = 0; j <= nB; j++)
        scores[j] = 0;
    for (i = 1; i < nA; i++) {
        kA = CHARINDEX(sA[i-1]);
        temp = 0;
        for (j = 1; j < nB; j++) {
            kB = CHARINDEX(sB[j-1]);
            SELECT_SCORE_LOCAL3(temp + self->substitution_matrix[kA][kB],
                                scores[j] + gap_extend_B,
                                scores[j-1] + gap_extend_A);
            temp = scores[j];
            scores[j] = score;
        }
        kB = CHARINDEX(sB[nB-1]);
        SELECT_SCORE_LOCAL1(temp + self->substitution_matrix[kA][kB]);
        temp = scores[nB];
        scores[nB] = score;
    }
    kA = CHARINDEX(sA[nA-1]);
    temp = 0;
    for (j = 1; j < nB; j++) {
        kB = CHARINDEX(sB[j-1]);
        SELECT_SCORE_LOCAL1(temp + self->substitution_matrix[kA][kB]);
        temp = scores[j];
        scores[j] = score;
    }
    kB = CHARINDEX(sB[nB-1]);
    SELECT_SCORE_LOCAL1(temp + self->substitution_matrix[kA][kB]);
    result = PyFloat_FromDouble(maximum);
exit:
    if (scores) PyMem_Free(scores);
    if (!result) PyErr_SetString(PyExc_MemoryError, "Out of memory");
    return result;
}

static PyObject* PathGenerator_next_needlemanwunsch(PathGenerator* self)
{
    int i = 0;
    int j = 0;
    int path;
    int trace = 0;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M.trace1;

    path = M[i][j].path;
    if (path == DONE) return NULL;
    if (path == 0) {
        /* Generate the first path. */
        i = nA;
        j = nB;
    }
    else {
        /* We already have a path. Prune the path to see if there are
         * any alternative paths. */
        while (1) {
            if (path == HORIZONTAL) {
                trace = M[i][++j].trace;
                if (trace & VERTICAL) {
                    M[--i][j].path = VERTICAL;
                    break;
                }
                if (trace & DIAGONAL) {
                    M[--i][--j].path = DIAGONAL;
                    break;
                }
            }
            else if (path == VERTICAL) {
                trace = M[++i][j].trace;
                if (trace & DIAGONAL) {
                    M[--i][--j].path = DIAGONAL;
                    break;
                }
            }
            else /* DIAGONAL */ {
                i++;
                j++;
            }
            path = M[i][j].path;
            if (!path) {
                /* we reached the end of the alignment without finding
                 * an alternative path */
                M[0][0].path = DONE;
                return NULL;
            }
        }
    }
    /* Follow the traceback until we reach the origin. */
    while (1) {
        trace = M[i][j].trace;
        if (trace & HORIZONTAL) M[i][--j].path = HORIZONTAL;
        else if (trace & VERTICAL) M[--i][j].path = VERTICAL;
        else if (trace & DIAGONAL) M[--i][--j].path = DIAGONAL;
        else break;
    }
    return _create_path_needleman_wunsch_smith_waterman(M, 0, 0);
}

static PyObject* PathGenerator_next_smithwaterman(PathGenerator* self)
{
    int trace = 0;
    int i = self->iA;
    int j = self->iB;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M.trace1;
    int path = M[0][0].path;

    if (path == DONE || path == NONE) return NULL;

    path = M[i][j].path;
    if (path) {
        /* We already have a path. Prune the path to see if there are
         * any alternative paths. */
        while (1) {
            if (path == HORIZONTAL) {
                trace = M[i][++j].trace;
                if (trace & VERTICAL) {
                    M[--i][j].path = VERTICAL;
                    break;
                }
                else if (trace & DIAGONAL) {
                    M[--i][--j].path = DIAGONAL;
                    break;
                }
            }
            else if (path == VERTICAL) {
                trace = M[++i][j].trace;
                if (trace & DIAGONAL) {
                    M[--i][--j].path = DIAGONAL;
                    break;
                }
            }
            else /* DIAGONAL */ {
                i++;
                j++;
            }
            path = M[i][j].path;
            if (!path) break;
        }
    }

    if (path) {
        trace = M[i][j].trace;
    } else {
        /* Find a suitable end point for a path.
         * Only allow start points ending at the M matrix. */
        while (1) {
            if (j < nB) j++;
            else if (i < nA) {
                i++;
                j = 0;
            }
            else {
                /* we reached the end of the sequences without finding
                 * an alternative path */
                M[0][0].path = DONE;
                return NULL;
            }
            trace = M[i][j].trace;
            if (trace & ENDPOINT) {
                trace &= DIAGONAL; /* exclude paths ending in a gap */
                break;
            }
        }
        M[i][j].path = 0;
    }

    /* Follow the traceback until we reach the origin. */
    while (1) {
        if (trace & HORIZONTAL) M[i][--j].path = HORIZONTAL;
        else if (trace & VERTICAL) M[--i][j].path = VERTICAL;
        else if (trace & DIAGONAL) M[--i][--j].path = DIAGONAL;
        else break;
        trace = M[i][j].trace;
    }
    self->iA = i;
    self->iB = j;
    return _create_path_needleman_wunsch_smith_waterman(M, i, j);
}

static PyObject* PathGenerator_next_gotoh_global(PathGenerator* self)
{
    int i = 0;
    int j = 0;
    int m;
    int path;
    int trace = 0;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace3** trace3 = self->trace3;

    m = M_MATRIX;
    path = trace3[i][j].path;
    if (path == DONE) return NULL;
    if (path == 0) {
        i = nA;
        j = nB;
    }
    else {
        /* We already have a path. Prune the path to see if there are
         * any alternative paths. */
        while (1) {
            path = trace3[i][j].path;
            if (path == 0) {
                switch (m) {
                    case M_MATRIX: m = Ix_MATRIX; break;
                    case Ix_MATRIX: m = Iy_MATRIX; break;
                    case Iy_MATRIX: m = 0; break;
                }
                break;
            }
            switch (path) {
                case HORIZONTAL: trace = trace3[i][++j].Iy; break;
                case VERTICAL: trace = trace3[++i][j].Ix; break;
                case DIAGONAL: trace = trace3[++i][++j].M; break;
            }
            switch (m) {
                case M_MATRIX:
                    if (trace & Ix_MATRIX) {
                        m = Ix_MATRIX;
                        break;
                    }
                case Ix_MATRIX:
                    if (trace & Iy_MATRIX) {
                        m = Iy_MATRIX;
                        break;
                    }
                case Iy_MATRIX:
                default:
                    switch (path) {
                        case HORIZONTAL: m = Iy_MATRIX; break;
                        case VERTICAL: m = Ix_MATRIX; break;
                        case DIAGONAL: m = M_MATRIX; break;
                    }
                    continue;
            }
            switch (path) {
                case HORIZONTAL: j--; break;
                case VERTICAL: i--; break;
                case DIAGONAL: i--; j--; break;
            }
            trace3[i][j].path = path;
            break;
        }
    }

    if (path == 0) {
        /* Generate a new path. */
        switch (m) {
            case M_MATRIX:
                if (trace3[nA][nB].M) {
                   /* m = M_MATRIX; */
                   break;
                }
            case Ix_MATRIX:
                if (trace3[nA][nB].Ix) {
                   m = Ix_MATRIX;
                   break;
                }
            case Iy_MATRIX:
                if (trace3[nA][nB].Iy) {
                   m = Iy_MATRIX;
                   break;
                }
            default:
                /* exhausted this generator */
                trace3[0][0].path = DONE;
                return NULL;
        }
    }

    switch (m) {
        case M_MATRIX:
            trace = trace3[i][j].M;
            path = DIAGONAL;
            i--; j--;
            break;
        case Ix_MATRIX:
            trace = trace3[i][j].Ix;
            path = VERTICAL;
            i--;
            break;
        case Iy_MATRIX:
            trace = trace3[i][j].Iy;
            path = HORIZONTAL;
            j--;
            break;
    }

    while (1) {
        if (trace & M_MATRIX) {
            trace3[i][j].path = path;
            trace = trace3[i][j].M;
            path = DIAGONAL;
            i--; j--;
        }
        else if (trace & Ix_MATRIX) {
            trace3[i][j].path = path;
            trace = trace3[i][j].Ix;
            path = VERTICAL;
            i--;
        }
        else if (trace & Iy_MATRIX) {
            trace3[i][j].path = path;
            trace = trace3[i][j].Iy;
            path = HORIZONTAL;
            j--;
        }
        else break;
    }
    return _create_path_gotoh(trace3, 0, 0);
}

static PyObject* PathGenerator_next_gotoh_local(PathGenerator* self)
{
    int trace = 0;
    int i;
    int j;
    int m = M_MATRIX;
    int iA = self->iA;
    int iB = self->iB;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace3** trace3 = self->trace3;
    int path = trace3[0][0].path;

    if (path == DONE) return NULL;

    path = trace3[iA][iB].path;

    if (path) {
        i = iA;
        j = iB;
        while (1) {
            /* We already have a path. Prune the path to see if there are
             * any alternative paths. */
            path = trace3[i][j].path;
            if (path == 0) {
                m = M_MATRIX;
                iA = i;
                iB = j;
                break;
            }
            switch (path) {
                case HORIZONTAL: trace = trace3[i][++j].Iy; break;
                case VERTICAL: trace = trace3[++i][j].Ix; break;
                case DIAGONAL: trace = trace3[++i][++j].M; break;
            }
            switch (m) {
                case M_MATRIX:
                    if (trace & Ix_MATRIX) {
                        m = Ix_MATRIX;
                        break;
                    }
                case Ix_MATRIX:
                    if (trace & Iy_MATRIX) {
                        m = Iy_MATRIX;
                        break;
                    }
                case Iy_MATRIX:
                default:
                    switch (path) {
                        case HORIZONTAL: m = Iy_MATRIX; break;
                        case VERTICAL: m = Ix_MATRIX; break;
                        case DIAGONAL: m = M_MATRIX; break;
                    }
                    continue;
            }
            switch (path) {
                case HORIZONTAL: j--; break;
                case VERTICAL: i--; break;
                case DIAGONAL: i--; j--; break;
            }
            trace3[i][j].path = path;
            break;
        }
    }

    if (path == 0) {
        /* Find the end point for a new path. */
        while (1) {
            if (iB < nB) iB++;
            else if (iA < nA) {
                iA++;
                iB = 0;
            }
            else {
                /* we reached the end of the alignment without finding
                 * an alternative path */
                trace3[0][0].path = DONE;
                return NULL;
            }
            if (trace3[iA][iB].M & ENDPOINT) {
                trace3[iA][iB].path = 0;
                break;
            }
        }
        m = M_MATRIX;
        i = iA;
        j = iB;
    }

    while (1) {
        switch (m) {
            case M_MATRIX: trace = trace3[i][j].M; break;
            case Ix_MATRIX: trace = trace3[i][j].Ix; break;
            case Iy_MATRIX: trace = trace3[i][j].Iy; break;
        }
        if (trace==0) {
            self->iA = i;
            self->iB = j;
            return _create_path_gotoh(trace3, i, j);
        }
        switch (m) {
            case M_MATRIX:
                path = DIAGONAL;
                i--;
                j--;
                break;
            case Ix_MATRIX:
                path = VERTICAL;
                i--;
                break;
            case Iy_MATRIX:
                path = HORIZONTAL;
                j--;
                break;
        }
        if (trace & M_MATRIX) m = M_MATRIX;
        else if (trace & Ix_MATRIX) m = Ix_MATRIX;
        else if (trace & Iy_MATRIX) m = Iy_MATRIX;
/*
        switch (m) {
            case M_MATRIX: M[i][j].path = path; break;
            case Ix_MATRIX: Ix[i][j].path = path; break;
            case Iy_MATRIX: Iy[i][j].path = path; break;
        }
*/
        trace3[i][j].path = path;
    }
    return NULL;
}

static PyObject*
PathGenerator_next_waterman_smith_beyer_global(PathGenerator* self)
{
    int i, j;
    int iA, iB;
    int path;
    int step;
    int trace;
    int* traceXY;
    int* traceM;

    int m = M_MATRIX;
    const int nA = self->nA;
    const int nB = self->nB;
    CellM** M = self->M.general;
    CellXY** Ix = self->Ix;
    CellXY** Iy = self->Iy;

    if (M[0][0].path == DONE) return NULL;

    if (M[0][0].path) {
        /* We already have a path. Prune the path to see if there are
         * any alternative paths. */
        i = 0;
        j = 0;
        while (1) {
            switch (m) {
                case M_MATRIX:
                    path = M[i][j].path;
                    step = M[i][j].step;
                    break;
                case Ix_MATRIX:
                    path = Ix[i][j].path;
                    step = Ix[i][j].step;
                    break;
                case Iy_MATRIX:
                    path = Iy[i][j].path;
                    step = Iy[i][j].step;
                    break;
            }
            if (path == HORIZONTAL) {
                iA = i;
                iB = j + step;
            }
            else if (path == VERTICAL) {
                iA = i + step;
                iB = j;
            }
            else if (path == DIAGONAL) {
                iA = i + 1;
                iB = j + 1;
            }
            else if (path == 0) {
                m <<= 1;
                iA = -1;
                break;
            }
            else printf("RUNTIME ERROR\n");
            if (i == iA) { /* HORIZONTAL */
                traceXY = Iy[iA][iB].traceXY;
                traceM = Iy[iA][iB].traceM;
                if (m==M_MATRIX) {
                    while (*traceM != j) traceM++;
                    traceM++;
                    j = *traceM;
                    if (j >= 0) {
                        M[i][j].path = HORIZONTAL;
                        M[i][j].step = iB - j;
                        break;
                    }
                } else if (m==Ix_MATRIX) {
                    while (*traceXY != j) traceXY++;
                    traceXY++;
                }
                j = *traceXY;
                if (j >= 0) {
                    m = Ix_MATRIX;
                    Ix[i][j].path = HORIZONTAL;
                    Ix[i][j].step = iB - j;
                    break;
                }
                /* no alternative found; continue pruning */
                m = Iy_MATRIX;
                j = iB;
            }
            else if (j==iB) { /* VERTICAL */
                traceXY = Ix[iA][iB].traceXY;
                traceM = Ix[iA][iB].traceM;
                if (m==M_MATRIX) {
                    while (*traceM != i) traceM++;
                    traceM++;
                    i = *traceM;
                    if (i >= 0) {
                        M[i][j].path = VERTICAL;
                        M[i][j].step = iA - i;
                        break;
                    }
                } else if (m==Iy_MATRIX) {
                    while (*traceXY != i) traceXY++;
                    traceXY++;
                }
                i = *traceXY;
                if (i >= 0) {
                    m = Iy_MATRIX;
                    Iy[i][j].path = VERTICAL;
                    Iy[i][j].step = iA - i;
                    break;
                }
                /* no alternative found; continue pruning */
                m = Ix_MATRIX;
                i = iA;
            }
            else { /* DIAGONAL */
                i = iA - 1;
                j = iB - 1;
                trace = M[iA][iB].trace;
                switch (m) {
                    case M_MATRIX:
                        if (trace & Ix_MATRIX) {
                            m = Ix_MATRIX;
                            Ix[i][j].path = DIAGONAL;
                            break;
                        }
                    case Ix_MATRIX:
                        if (trace & Iy_MATRIX) {
                            m = Iy_MATRIX;
                            Iy[i][j].path = DIAGONAL;
                            break;
                        }
                    case Iy_MATRIX:
                    default:
                        /* no alternative found; continue pruning */
                        m = M_MATRIX;
                        i = iA;
                        j = iB;
                        continue;
                }
                /* alternative found; build path until starting point */
                break;
            }
        }
    } else iA = -1;

    if (iA < 0) {
        /* Find a suitable end point for a path. */
        switch (m) {
            case M_MATRIX:
                if (M[nA][nB].trace) {
                    /* m = M_MATRIX; */
                    break;
                }
            case Ix_MATRIX:
                if (Ix[nA][nB].traceM[0] != -1 || Ix[nA][nB].traceXY[0] != -1) {
                    m = Ix_MATRIX;
                    break;
                }
            case Iy_MATRIX:
                if (Iy[nA][nB].traceM[0] != -1 || Iy[nA][nB].traceXY[0] != -1) {
                    m = Iy_MATRIX;
                    break;
                }
            default:
                M[0][0].path = DONE;
                return NULL;
        }
        i = nA;
        j = nB;
    }

    /* Follow the traceback until we reach the origin. */
    while (1) {
        switch (m) {
            case M_MATRIX:
                iA = i-1;
                iB = j-1;
                trace = M[i][j].trace;
                if (trace & M_MATRIX) {
                    m = M_MATRIX;
                    M[iA][iB].path = DIAGONAL;
                }
                else if (trace & Ix_MATRIX) {
                    m = Ix_MATRIX;
                    Ix[iA][iB].path = DIAGONAL;
                }
                else if (trace & Iy_MATRIX) {
                    m = Iy_MATRIX;
                    Iy[iA][iB].path = DIAGONAL;
                } else {
                    return _create_path_waterman_smith_beyer(M, Ix, Iy, i, j);
                }
                i = iA;
                j = iB;
                continue;
            case Ix_MATRIX:
                traceXY = Ix[i][j].traceXY;
                traceM = Ix[i][j].traceM;
                iB = j;
                iA = *traceM;
                if (iA >= 0) m = M_MATRIX;
                else {
                    iA = *traceXY;
                    m = Iy_MATRIX;
                }
                break;
            case Iy_MATRIX:
                traceXY = Iy[i][j].traceXY;
                traceM = Iy[i][j].traceM;
                iA = i;
                iB = *traceM;
                if (iB >= 0) m = M_MATRIX;
                else {
                    iB = *traceXY;
                    m = Ix_MATRIX;
                }
                break;
        }
        switch (m) {
            case M_MATRIX:
                if (i==iA) {
                    M[iA][iB].path = HORIZONTAL;
                    M[iA][iB].step = j - iB;
                }
                else if (j==iB) {
                    M[iA][iB].path = VERTICAL;
                    M[iA][iB].step = i - iA;
                }
                else M[iA][iB].path = DIAGONAL;
                break;
            case Ix_MATRIX:
                Ix[iA][iB].path = HORIZONTAL;
                Ix[iA][iB].step = j - iB;
                break;
            case Iy_MATRIX:
                Iy[iA][iB].path = VERTICAL;
                Iy[iA][iB].step = i - iA;
                break;
        }
        i = iA;
        j = iB;
    }
}

static PyObject*
PathGenerator_next_waterman_smith_beyer_local(PathGenerator* self)
{
    int i, j, m;
    int trace = 0;
    int* traceXY;
    int* traceM;

    int iA = self->iA;
    int iB = self->iB;
    const int nA = self->nA;
    const int nB = self->nB;
    CellM** M = self->M.general;
    CellXY** Ix = self->Ix;
    CellXY** Iy = self->Iy;

    if (M[0][0].path == DONE) return NULL;
    m = 0;
    if (M[iA][iB].path) {
        /* We already have a path. Prune the path to see if there are
         * any alternative paths. */
        m = M_MATRIX;
        i = iA;
        j = iB;
        while (1) {
            switch (m) {
                case M_MATRIX:
                    if (M[i][j].path == HORIZONTAL) {
                        iA = i;
                        iB = j + M[i][j].step;
                    }
                    else if (M[i][j].path == VERTICAL) {
                        iA = i + M[i][j].step;
                        iB = j;
                    }
                    else if (M[i][j].path == DIAGONAL) {
                        iA = i + 1;
                        iB = j + 1;
                    }
                    else iA = -1;
                    break;
                case Ix_MATRIX:
                    if (Ix[i][j].path == HORIZONTAL) {
                        iA = i;
                        iB = j + Ix[i][j].step;
                    }
                    else if (Ix[i][j].path == VERTICAL) {
                        iA = i + Ix[i][j].step;
                        iB = j;
                    }
                    else if (Ix[i][j].path == DIAGONAL) {
                        iA = i + 1;
                        iB = j + 1;
                    }
                    else printf("RUNTIME ERROR 12\n");
                    break;
                case Iy_MATRIX:
                    if (Iy[i][j].path == HORIZONTAL) {
                        iA = i;
                        iB = j + Iy[i][j].step;
                    }
                    else if (Iy[i][j].path == VERTICAL) {
                        iA = i + Iy[i][j].step;
                        iB = j;
                    }
                    else if (Iy[i][j].path == DIAGONAL) {
                        iA = i + 1;
                        iB = j + 1;
                    }
                    else printf("RUNTIME ERROR 12\n");
                    break;
            }
            if (iA < 0) {
                m = 0;
                iA = i;
                iB = j;
                break;
            }
            if (i == iA) { /* HORIZONTAL */
                traceXY = Iy[iA][iB].traceXY;
                traceM = Iy[iA][iB].traceM;
                if (m==M_MATRIX) {
                    while (*traceM != j) traceM++;
                    traceM++;
                    j = *traceM;
                    if (j >= 0) {
                        M[i][j].path = HORIZONTAL;
                        M[i][j].step = iB - j;
                        break;
                    }
                } else if (m==Ix_MATRIX) {
                    while (*traceXY != j) traceXY++;
                    traceXY++;
                }
                j = *traceXY;
                if (j >= 0) {
                    m = Ix_MATRIX;
                    Ix[i][j].path = HORIZONTAL;
                    Ix[i][j].step = iB - j;
                    break;
                }
                /* no alternative found; continue pruning */
                m = Iy_MATRIX;
                j = iB;
            }
            else if (j == iB) { /* VERTICAL */
                traceXY = Ix[iA][iB].traceXY;
                traceM = Ix[iA][iB].traceM;
                if (m==M_MATRIX) {
                    while (*traceM != i) traceM++;
                    traceM++;
                    i = *traceM;
                    if (i >= 0) {
                        M[i][j].path = VERTICAL;
                        M[i][j].step = iA - i;
                        break;
                    }
                } else if (m==Iy_MATRIX) {
                    while (*traceXY != i) traceXY++;
                    traceXY++;
                }
                i = *traceXY;
                if (i >= 0) {
                    m = Iy_MATRIX;
                    Iy[i][j].path = VERTICAL;
                    Iy[i][j].step = iA - i;
                    break;
                }
                /* no alternative found; continue pruning */
                m = Ix_MATRIX;
                i = iA;
            }
            else { /* DIAGONAL */
                i = iA - 1;
                j = iB - 1;
                trace = M[iA][iB].trace;
                switch (m) {
                    case M_MATRIX:
                        if (trace & Ix_MATRIX) {
                            m = Ix_MATRIX;
                            Ix[i][j].path = DIAGONAL;
                            break;
                        }
                    case Ix_MATRIX:
                        if (trace & Iy_MATRIX) {
                            m = Iy_MATRIX;
                            Iy[i][j].path = DIAGONAL;
                            break;
                        }
                    case Iy_MATRIX:
                    default:
                        /* no alternative found; continue pruning */
                        m = M_MATRIX;
                        i = iA;
                        j = iB;
                        continue;
                }
                /* alternative found; build path until starting point */
                break;
            }
        }
    }
 
    if (m == 0) {
        /* We are at [nA][nB]. Find a suitable end point for a path. */
        while (1) {
            if (iB < nB) iB++;
            else if (iA < nA) {
                iA++;
                iB = 0;
            }
            else {
                /* exhausted this generator */
                M[0][0].path = DONE;
                return NULL;
            }
            if (M[iA][iB].trace & ENDPOINT) break;
        }
        M[iA][iB].path = 0;
        m = M_MATRIX;
        i = iA;
        j = iB;
    }

    /* Follow the traceback until we reach the origin. */
    while (1) {
        switch (m) {
            case Ix_MATRIX:
                traceXY = Ix[i][j].traceXY;
                traceM = Ix[i][j].traceM;
                iB = j;
                iA = *traceM;
                if (iA >= 0) m = M_MATRIX;
                else {
                    iA = *traceXY;
                    m = Iy_MATRIX;
                }
                break;
            case Iy_MATRIX:
                traceXY = Iy[i][j].traceXY;
                traceM = Iy[i][j].traceM;
                iA = i;
                iB = *traceM;
                if (iB >= 0) m = M_MATRIX;
                else {
                    iB = *traceXY;
                    m = Ix_MATRIX;
                }
                break;
            case M_MATRIX:
                iA = i-1;
                iB = j-1;
                trace = M[i][j].trace;
                if (trace & M_MATRIX) {
                    m = M_MATRIX;
                    M[iA][iB].path = DIAGONAL;
                }
                else if (trace & Ix_MATRIX) {
                    m = Ix_MATRIX;
                    Ix[iA][iB].path = DIAGONAL;
                }
                else if (trace & Iy_MATRIX) {
                    m = Iy_MATRIX;
                    Iy[iA][iB].path = DIAGONAL;
                } else {
                    self->iA = i;
                    self->iB = j;
                    return _create_path_waterman_smith_beyer(M, Ix, Iy, i, j);
                }
        }
        switch (m) {
            case M_MATRIX:
                if (iA == i) {
                    M[iA][iB].path = HORIZONTAL;
                    M[iA][iB].step = j - iB;
                }
                else if (iB == j) {
                    M[iA][iB].path = VERTICAL;
                    M[iA][iB].step = i - iA;
                }
                else M[iA][iB].path = DIAGONAL;
                break;
            case Ix_MATRIX:
                if (iA == i) {
                    Ix[iA][iB].path = HORIZONTAL;
                    Ix[iA][iB].step = j - iB;
                }
                else if (iB == j) {
                    Ix[iA][iB].path = VERTICAL;
                    Ix[iA][iB].step = i - iA;
                }
                else Ix[iA][iB].path = DIAGONAL;
                break;
            case Iy_MATRIX:
                if (iA == i) {
                    Iy[iA][iB].path = HORIZONTAL;
                    Iy[iA][iB].step = j - iB;
                }
                else if (iB == j) {
                    Iy[iA][iB].path = VERTICAL;
                    Iy[iA][iB].step = i - iA;
                }
                else Iy[iA][iB].path = DIAGONAL;
                break;
        }
        i = iA;
        j = iB;
    }
}

static PyObject *
PathGenerator_next(PathGenerator* self)
{
    const Mode mode = self->mode;
    const Algorithm algorithm = self->algorithm;
    switch (algorithm) {
        case NeedlemanWunschSmithWaterman:
            switch (mode) {
                case Global:
                    return PathGenerator_next_needlemanwunsch(self);
                case Local:
                    return PathGenerator_next_smithwaterman(self);
            }
        case Gotoh:
            switch (mode) {
                case Global:
                    return PathGenerator_next_gotoh_global(self);
                case Local:
                    return PathGenerator_next_gotoh_local(self);
            }
        case WatermanSmithBeyer:
            switch (mode) {
                case Global:
                    return PathGenerator_next_waterman_smith_beyer_global(self);
                case Local:
                    return PathGenerator_next_waterman_smith_beyer_local(self);
            }
        case Unknown:
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unknown algorithm");
            return NULL;
    }
}

static void
PathGenerator_dealloc(PathGenerator* self)
{
    int i;
    const int nA = self->nA;
    const Algorithm algorithm = self->algorithm;
    switch (algorithm) {
        case NeedlemanWunschSmithWaterman: {
            Trace** M = self->M.trace1;
            for (i = 0; i <= nA; i++) PyMem_Free(M[i]);
            PyMem_Free(M);
            break;
        }
        case Gotoh: {
            Trace3** trace3 = self->trace3;
            for (i = 0; i <= nA; i++) PyMem_Free(trace3[i]);
            PyMem_Free(trace3);
            break;
        }
        case WatermanSmithBeyer: {
            CellM** M = self->M.general;
            CellXY** Ix = self->Ix;
            CellXY** Iy = self->Iy;
            const Mode mode = self->mode;
            const int nB = self->nB;
            _deallocate_watermansmithbeyer_matrices(nA, nB, M, Ix, Iy, mode);
            break;
        }
        case Unknown:
        default:
            PyErr_WriteUnraisable((PyObject*)self);
            return;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static const char PathGenerator_reset__doc__[] = "reset the iterator";

static PyObject*
PathGenerator_reset(PathGenerator* self)
{
    switch (self->mode) {
        case Local:
            self->iA = 0;
            self->iB = 0;
        case Global:
            switch (self->algorithm) {
                case NeedlemanWunschSmithWaterman: {
                    Trace** M = self->M.trace1;
                    if (M[0][0].path != NONE) M[0][0].path = 0;
                    break;
                }
                case Gotoh: {
                    Trace3** trace3 = self->trace3;
                    if (trace3[0][0].path != NONE) trace3[0][0].path = 0;
                    break;
                }
                case WatermanSmithBeyer: {
                    CellM** M = self->M.general;
                    M[0][0].path = 0;
                    break;
                }
                case Unknown:
                default:
                    break;
        }
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef PathGenerator_methods[] = {
    {"reset",
     (PyCFunction)PathGenerator_reset,
     METH_NOARGS,
     PathGenerator_reset__doc__
    },
    {NULL}  /* Sentinel */
};

#define SAFE_ADD(t, s) \
{   term = t; \
    if (term > PY_SSIZE_T_MAX - s) { count = OVERFLOW_ERROR; goto exit; } \
    s += term; \
}

static Py_ssize_t
PathGenerator_needlemanwunsch_length(PathGenerator* self)
{
    int i;
    int j;
    int trace;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M.trace1;
    Py_ssize_t term;
    Py_ssize_t count;
    Py_ssize_t temp;
    Py_ssize_t* counts;
    counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!counts) return MEMORY_ERROR;
    counts[0] = 1;
    for (j = 1; j <= nB; j++) {
        trace = M[0][j].trace;
        count = 0;
        if (trace & HORIZONTAL) SAFE_ADD(counts[j-1], count);
        counts[j] = count;
    }
    for (i = 1; i <= nA; i++) {
        trace = M[i][0].trace;
        count = 0;
        if (trace & VERTICAL) SAFE_ADD(counts[0], count);
        temp = counts[0];
        counts[0] = count;
        for (j = 1; j <= nB; j++) {
            trace = M[i][j].trace;
            count = 0;
            if (trace & HORIZONTAL) SAFE_ADD(counts[j-1], count);
            if (trace & VERTICAL) SAFE_ADD(counts[j], count);
            if (trace & DIAGONAL) SAFE_ADD(temp, count);
            temp = counts[j];
            counts[j] = count;
        }
    }
exit:
    PyMem_Free(counts);
    return count;
}

static Py_ssize_t
PathGenerator_smithwaterman_length(PathGenerator* self)
{
    int i;
    int j;
    int trace;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace** M = self->M.trace1;
    Py_ssize_t term;
    Py_ssize_t count;
    Py_ssize_t total = 0;
    Py_ssize_t temp;
    Py_ssize_t* counts;
    counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!counts) return MEMORY_ERROR;
    counts[0] = 1;
    for (j = 1; j <= nB; j++) counts[j] = 1;
    for (i = 1; i <= nA; i++) {
        temp = counts[0];
        counts[0] = 1;
        for (j = 1; j <= nB; j++) {
            trace = M[i][j].trace;
            count = 0;
            if (trace & DIAGONAL) SAFE_ADD(temp, count);
            if (M[i][j].trace & ENDPOINT) SAFE_ADD(count, total);
            if (trace & HORIZONTAL) SAFE_ADD(counts[j-1], count);
            if (trace & VERTICAL) SAFE_ADD(counts[j], count);
            temp = counts[j];
            if (count==0) count = 1;
            counts[j] = count;
        }
    }
    count = total;
exit:
    PyMem_Free(counts);
    return count;
}

static Py_ssize_t
PathGenerator_gotoh_global_length(PathGenerator* self)
{
    int i;
    int j;
    int trace;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace3** trace3 = self->trace3;
    Py_ssize_t count = MEMORY_ERROR;
    Py_ssize_t term;
    Py_ssize_t M_temp;
    Py_ssize_t Ix_temp;
    Py_ssize_t Iy_temp;
    Py_ssize_t* M_counts = NULL;
    Py_ssize_t* Ix_counts = NULL;
    Py_ssize_t* Iy_counts = NULL;
    M_counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!M_counts) goto exit;
    Ix_counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!Ix_counts) goto exit;
    Iy_counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!Iy_counts) goto exit;
    M_counts[0] = 1;
    Ix_counts[0] = 0;
    Iy_counts[0] = 0;
    for (j = 1; j <= nB; j++) {
        M_counts[j] = 0;
        Ix_counts[j] = 0;
        Iy_counts[j] = 1;
    }
    for (i = 1; i <= nA; i++) {
        M_temp = M_counts[0];
        M_counts[0] = 0;
        Ix_temp = Ix_counts[0];
        Ix_counts[0] = 1;
        Iy_temp = Iy_counts[0];
        Iy_counts[0] = 0;
        for (j = 1; j <= nB; j++) {
            count = 0;
            trace = trace3[i][j].M;
            if (trace & M_MATRIX) SAFE_ADD(M_temp, count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_temp, count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_temp, count);
            M_temp = M_counts[j];
            M_counts[j] = count;
            count = 0;
            trace = trace3[i][j].Ix;
            if (trace & M_MATRIX) SAFE_ADD(M_temp, count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_counts[j], count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_counts[j], count);
            Ix_temp = Ix_counts[j];
            Ix_counts[j] = count;
            count = 0;
            trace = trace3[i][j].Iy;
            if (trace & M_MATRIX) SAFE_ADD(M_counts[j-1], count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_counts[j-1], count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_counts[j-1], count);
            Iy_temp = Iy_counts[j];
            Iy_counts[j] = count;
        }
    }
    count = 0;
    if (trace3[nA][nB].M) SAFE_ADD(M_counts[nB], count);
    if (trace3[nA][nB].Ix) SAFE_ADD(Ix_counts[nB], count);
    if (trace3[nA][nB].Iy) SAFE_ADD(Iy_counts[nB], count);
exit:
    if (M_counts) PyMem_Free(M_counts);
    if (Ix_counts) PyMem_Free(Ix_counts);
    if (Iy_counts) PyMem_Free(Iy_counts);
    return count;
}

static Py_ssize_t
PathGenerator_gotoh_local_length(PathGenerator* self)
{
    int i;
    int j;
    int trace;
    const int nA = self->nA;
    const int nB = self->nB;
    Trace3** trace3 = self->trace3;
    Py_ssize_t term;
    Py_ssize_t count = MEMORY_ERROR;
    Py_ssize_t total = 0;
    Py_ssize_t M_temp;
    Py_ssize_t Ix_temp;
    Py_ssize_t Iy_temp;
    Py_ssize_t* M_counts = NULL;
    Py_ssize_t* Ix_counts = NULL;
    Py_ssize_t* Iy_counts = NULL;
    M_counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!M_counts) goto exit;
    Ix_counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!Ix_counts) goto exit;
    Iy_counts = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
    if (!Iy_counts) goto exit;
    M_counts[0] = 1;
    Ix_counts[0] = 0;
    Iy_counts[0] = 0;
    for (j = 1; j <= nB; j++) {
        M_counts[j] = 1;
        Ix_counts[j] = 0;
        Iy_counts[j] = 0;
    }
    for (i = 1; i <= nA; i++) {
        M_temp = M_counts[0];
        M_counts[0] = 1;
        Ix_temp = Ix_counts[0];
        Ix_counts[0] = 0;
        Iy_temp = Iy_counts[0];
        Iy_counts[0] = 0;
        for (j = 1; j <= nB; j++) {
            count = 0;
            trace = trace3[i][j].M;
            if (trace & M_MATRIX) SAFE_ADD(M_temp, count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_temp, count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_temp, count);
            if (count==0) count = 1;
            M_temp = M_counts[j];
            M_counts[j] = count;
            if (trace3[i][j].M & ENDPOINT) SAFE_ADD(count, total);
            count = 0;
            trace = trace3[i][j].Ix;
            if (trace & M_MATRIX) SAFE_ADD(M_temp, count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_counts[j], count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_counts[j], count);
            Ix_temp = Ix_counts[j];
            Ix_counts[j] = count;
            count = 0;
            trace = trace3[i][j].Iy;
            if (trace & M_MATRIX) SAFE_ADD(M_counts[j-1], count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_counts[j-1], count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_counts[j-1], count);
            Iy_temp = Iy_counts[j];
            Iy_counts[j] = count;
        }
    }
    count = total;
exit:
    if (M_counts) PyMem_Free(M_counts);
    if (Ix_counts) PyMem_Free(Ix_counts);
    if (Iy_counts) PyMem_Free(Iy_counts);
    return count;
}

static Py_ssize_t
PathGenerator_waterman_smith_beyer_global_length(PathGenerator* self)
{
    int i;
    int j;
    int k;
    int trace;
    int* tracep;
    const int nA = self->nA;
    const int nB = self->nB;
    CellM** M = self->M.general;
    CellXY** Ix = self->Ix;
    CellXY** Iy = self->Iy;
    Py_ssize_t count = MEMORY_ERROR;
    Py_ssize_t term;
    Py_ssize_t** M_count = NULL;
    Py_ssize_t** Ix_count = NULL;
    Py_ssize_t** Iy_count = NULL;
    M_count = PyMem_Malloc((nA+1)*sizeof(Py_ssize_t*));
    if (!M_count) goto exit;
    Ix_count = PyMem_Malloc((nA+1)*sizeof(Py_ssize_t*));
    if (!Ix_count) goto exit;
    Iy_count = PyMem_Malloc((nA+1)*sizeof(Py_ssize_t*));
    if (!Iy_count) goto exit;
    for (i = 0; i <= nA; i++) {
        M_count[i] = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
        if (!M[i]) goto exit;
        Ix_count[i] = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
        if (!Ix[i]) goto exit;
        Iy_count[i] = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
        if (!Iy[i]) goto exit;
    }
    for (i = 0; i <= nA; i++) {
        for (j = 0; j <= nB; j++) {
            count = 0;
            trace = M[i][j].trace;
            if (trace & M_MATRIX) SAFE_ADD(M_count[i-1][j-1], count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_count[i-1][j-1], count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_count[i-1][j-1], count);
            if (count == 0) count = 1; /* happens at M[0][0] only */
            M_count[i][j] = count;
            count = 0;
            tracep = Ix[i][j].traceM;
            if (tracep) {
                while (1) {
                    k = *tracep;
                    if (k < 0) break;
                    SAFE_ADD(M_count[k][j], count);
                    tracep++;
                }
            }
            tracep = Ix[i][j].traceXY;
            if (tracep) {
                while (1) {
                    k = *tracep;
                    if (k < 0) break;
                    SAFE_ADD(Iy_count[k][j], count);
                    tracep++;
                }
            }
            Ix_count[i][j] = count;
            count = 0;
            tracep = Iy[i][j].traceM;
            if (tracep) {
                while (1) {
                    k = *tracep;
                    if (k < 0) break;
                    SAFE_ADD(M_count[i][k], count);
                    tracep++;
                }
            }
            tracep = Iy[i][j].traceXY;
            if (tracep) {
                while (1) {
                    k = *tracep;
                    if (k < 0) break;
                    SAFE_ADD(Ix_count[i][k], count);
                    tracep++;
                }
            }
            Iy_count[i][j] = count;
        }
    }
    count = 0;
    if (M[nA][nB].trace)
        SAFE_ADD(M_count[nA][nB], count);
    if (Ix[nA][nB].traceM[0] != -1 || Ix[nA][nB].traceXY[0] != -1)
        SAFE_ADD(Ix_count[nA][nB], count);
     if (Iy[nA][nB].traceM[0] != -1 || Iy[nA][nB].traceXY[0] != -1)
        SAFE_ADD(Iy_count[nA][nB], count);
exit:
    if (M_count) {
        if (Ix_count) {
            if (Iy_count) {
                for (i = 0; i <= nA; i++) {
                    if (!M_count[i]) break;
                    PyMem_Free(M_count[i]);
                    if (!Ix_count[i]) break;
                    PyMem_Free(Ix_count[i]);
                    if (!Iy_count[i]) break;
                    PyMem_Free(Iy_count[i]);
                }
                PyMem_Free(Iy_count);
            }
            PyMem_Free(Ix_count);
        }
        PyMem_Free(M_count);
    }
    return count;
}

static Py_ssize_t
PathGenerator_waterman_smith_beyer_local_length(PathGenerator* self)
{
    int i;
    int j;
    int k;
    int trace;
    int* tracep;
    const int nA = self->nA;
    const int nB = self->nB;
    CellM** M = self->M.general;
    CellXY** Ix = self->Ix;
    CellXY** Iy = self->Iy;
    Py_ssize_t term;
    Py_ssize_t count = MEMORY_ERROR;
    Py_ssize_t total = 0;
    Py_ssize_t** M_count = NULL;
    Py_ssize_t** Ix_count = NULL;
    Py_ssize_t** Iy_count = NULL;
    M_count = PyMem_Malloc((nA+1)*sizeof(Py_ssize_t*));
    if (!M_count) goto exit;
    Ix_count = PyMem_Malloc((nA+1)*sizeof(Py_ssize_t*));
    if (!Ix_count) goto exit;
    Iy_count = PyMem_Malloc((nA+1)*sizeof(Py_ssize_t*));
    if (!Iy_count) goto exit;
    for (i = 0; i <= nA; i++) {
        M_count[i] = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
        if (!M[i]) goto exit;
        Ix_count[i] = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
        if (!Ix[i]) goto exit;
        Iy_count[i] = PyMem_Malloc((nB+1)*sizeof(Py_ssize_t));
        if (!Iy[i]) goto exit;
    }
    for (i = 0; i <= nA; i++) {
        for (j = 0; j <= nB; j++) {
            count = 0;
            trace = M[i][j].trace;
            if (trace & M_MATRIX) SAFE_ADD(M_count[i-1][j-1], count);
            if (trace & Ix_MATRIX) SAFE_ADD(Ix_count[i-1][j-1], count);
            if (trace & Iy_MATRIX) SAFE_ADD(Iy_count[i-1][j-1], count);
            if (count==0) count = 1;
            M_count[i][j] = count;
            if (M[i][j].trace & ENDPOINT) SAFE_ADD(count, total);
            count = 0;
            tracep = Ix[i][j].traceM;
            if (tracep) {
                while (1) {
                    k = *tracep;
                    if (k < 0) break;
                    SAFE_ADD(M_count[k][j], count);
                    tracep++;
                }
            }
            tracep = Ix[i][j].traceXY;
            if (tracep) {
                while (1) {
                    k = *tracep;
                    if (k < 0) break;
                    SAFE_ADD(Iy_count[k][j], count);
                    tracep++;
                }
            }
            if (count == 0) count = 1;
            Ix_count[i][j] = count;
            count = 0;
            tracep = Iy[i][j].traceM;
            if (tracep) {
                while (1) {
                    k = *tracep;
                    if (k < 0) break;
                    SAFE_ADD(M_count[i][k], count);
                    tracep++;
                }
            }
            tracep = Iy[i][j].traceXY;
            if (tracep) {
                while (1) {
                    k = *tracep;
                    if (k < 0) break;
                    SAFE_ADD(Ix_count[i][k], count);
                    tracep++;
                }
            }
            if (count == 0) count = 1;
            Iy_count[i][j] = count;
        }
    }
    count = total;
exit:
    if (M_count) {
        if (Ix_count) {
            if (Iy_count) {
                for (i = 0; i <= nA; i++) {
                    if (!M_count[i]) break;
                    PyMem_Free(M_count[i]);
                    if (!Ix_count[i]) break;
                    PyMem_Free(Ix_count[i]);
                    if (!Iy_count[i]) break;
                    PyMem_Free(Iy_count[i]);
                }
                PyMem_Free(Iy_count);
            }
            PyMem_Free(Ix_count);
        }
        PyMem_Free(M_count);
    }
    return count;
}

static Py_ssize_t PathGenerator_length(PathGenerator* self) {
    Py_ssize_t length = self->length;
    if (length == 0) {
        switch (self->algorithm) {
            case NeedlemanWunschSmithWaterman:
                switch (self->mode) {
                    case Global:
                        length = PathGenerator_needlemanwunsch_length(self);
                        break;
                    case Local:
                        length = PathGenerator_smithwaterman_length(self);
                        break;
                    default:
                        /* should not happen, but some compilers complain that
                         * that length can be used uninitialized.
                         */
                        PyErr_SetString(PyExc_RuntimeError, "Unknown mode");
                        return -1;
                }
                break;
            case Gotoh:
                switch (self->mode) {
                    case Global:
                        length = PathGenerator_gotoh_global_length(self);
                        break;
                    case Local:
                        length = PathGenerator_gotoh_local_length(self);
                        break;
                    default:
                        /* should not happen, but some compilers complain that
                         * that length can be used uninitialized.
                         */
                        PyErr_SetString(PyExc_RuntimeError, "Unknown mode");
                        return -1;
                }
                break;
            case WatermanSmithBeyer:
                switch (self->mode) {
                    case Global:
                        length = PathGenerator_waterman_smith_beyer_global_length(self);
                        break;
                    case Local:
                        length = PathGenerator_waterman_smith_beyer_local_length(self);
                        break;
                    default:
                        /* should not happen, but some compilers complain that
                         * that length can be used uninitialized.
                         */
                        PyErr_SetString(PyExc_RuntimeError, "Unknown mode");
                        return -1;
                }
                break;
            case Unknown:
            default:
                PyErr_SetString(PyExc_RuntimeError, "Unknown algorithm");
                return -1;
        }
        self->length = length;
    }
    switch (length) {
        case OVERFLOW_ERROR:
            PyErr_Format(PyExc_OverflowError,
                         "number of optimal alignments is larger than %zd",
                         PY_SSIZE_T_MAX);
            break;
        case MEMORY_ERROR:
            PyErr_SetString(PyExc_MemoryError, "Out of memory");
            break;
        default:
            break;
    }
    return length;
}

static PySequenceMethods PathGenerator_as_sequence = {
    (lenfunc)PathGenerator_length,  /* sq_length */
    NULL,                           /* sq_concat */
    NULL,                           /* sq_repeat */
    NULL,                           /* sq_item */
    NULL,                           /* sq_ass_item */
    NULL,                           /* sq_contains */
    NULL,                           /* sq_inplace_concat */
    NULL,                           /* sq_inplace_repeat */
};

static PyTypeObject PathGenerator_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "Path generator",               /* tp_name */
    sizeof(PathGenerator),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)PathGenerator_dealloc,  /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    0,                              /* tp_repr */
    0,                              /* tp_as_number */
    &PathGenerator_as_sequence,     /* tp_as_sequence */
    0,                              /* tp_as_mapping */
    0,                              /* tp_hash */
    0,                              /* tp_call */
    0,                              /* tp_str */
    0,                              /* tp_getattro */
    0,                              /* tp_setattro */
    0,                              /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,             /* tp_flags */
    0,                              /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    PyObject_SelfIter,              /* tp_iter */
    (iternextfunc)PathGenerator_next,      /* tp_iternext */
    PathGenerator_methods,          /* tp_methods */
};

static PathGenerator*
_create_path_generator(const Aligner* aligner, int nA, int nB, double epsilon)
{
    Algorithm algorithm = aligner->algorithm;
    PathGenerator* generator;
    generator = (PathGenerator*)PyType_GenericAlloc(&PathGenerator_Type, 0);
    if (!generator) return NULL;

    generator->iA = 0;
    generator->iB = 0;
    generator->nA = nA;
    generator->nB = nB;
    switch (algorithm) {
        case NeedlemanWunschSmithWaterman:
        case Gotoh:
            generator->M.trace1 = NULL;
            generator->Ix = NULL;
            generator->Iy = NULL;
            generator->trace3 = NULL;
            break;
        case WatermanSmithBeyer:
            generator->M.general = NULL;
            generator->Ix = NULL;
            generator->Iy = NULL;
            break;
        case Unknown:
        default:
            Py_DECREF(generator);
            PyErr_SetString(PyExc_RuntimeError, "unknown algorithm");
            return NULL;
    }
    generator->algorithm = algorithm;
    generator->mode = aligner->mode;
    generator->length = 0;

    return generator;
}

static PyObject*
Aligner_needlemanwunsch_align(Aligner* self, const char* sA, Py_ssize_t nA,
                                             const char* sB, Py_ssize_t nB)
{
    char c;
    int i;
    int j;
    int kA;
    int kB;
    const double gap_extend_A = self->target_extend_gap_score;
    const double gap_extend_B = self->query_extend_gap_score;
    const double left_gap_extend_A = self->target_left_extend_gap_score;
    const double left_gap_extend_B = self->query_left_extend_gap_score;
    const double right_gap_extend_A = self->target_right_extend_gap_score;
    const double right_gap_extend_B = self->query_right_extend_gap_score;
    const double epsilon = self->epsilon;
    Trace** M = NULL;
    double score;
    int trace;
    double temp;
    double* scores = NULL;
    PathGenerator* paths = NULL;

    /* Needleman-Wunsch algorithm */
    M = PyMem_Malloc((nA+1)*sizeof(Trace*));
    if (!M) goto exit;
    for (i = 0; i <= nA; i++) {
        M[i] = PyMem_Malloc((nB+1)*sizeof(Trace));
        if (!M[i]) goto exit;
    }
    scores = malloc((nB+1)*sizeof(double));
    if (!scores) goto exit;
    M[0][0].trace = 0;
    scores[0] = 0;
    for (j = 1; j <= nB; j++) {
        M[0][j].trace = HORIZONTAL;
        scores[j] = j * left_gap_extend_A;
    }
    for (i = 1; i < nA; i++) {
        temp = scores[0];
        M[i][0].trace = VERTICAL;
        scores[0] = i * left_gap_extend_B;
        kA = CHARINDEX(sA[i-1]);
        for (j = 1; j < nB; j++) {
            kB = CHARINDEX(sB[j-1]);
            SELECT_TRACE_NEEDLEMAN_WUNSCH(gap_extend_A, gap_extend_B);
        }
        kB = CHARINDEX(sB[j-1]);
        SELECT_TRACE_NEEDLEMAN_WUNSCH(gap_extend_A, right_gap_extend_B);
    }
    temp = scores[0];
    M[nA][0].trace = VERTICAL;
    scores[0] = i * left_gap_extend_B;
    kA = CHARINDEX(sA[nA-1]);
    for (j = 1; j < nB; j++) {
        kB = CHARINDEX(sB[j-1]);
        SELECT_TRACE_NEEDLEMAN_WUNSCH(right_gap_extend_A, gap_extend_B);
    }
    kB = CHARINDEX(sB[j-1]);
    SELECT_TRACE_NEEDLEMAN_WUNSCH(right_gap_extend_A, right_gap_extend_B);
    M[0][0].path = 0;
    M[nA][nB].path = 0;
    PyMem_Free(scores);

    paths = _create_path_generator(self, nA, nB, epsilon);
    if (paths) {
        PyObject* result;
        paths->M.trace1 = M;
        result = Py_BuildValue("fO", score, paths);
        Py_DECREF(paths);
        return result;
    }

exit:
    if (M) {
        for (i = 0; i <= nA; i++) {
            if (!M[i]) break;
            PyMem_Free(M[i]);
        }
        PyMem_Free(M);
    }
    PyErr_SetString(PyExc_MemoryError, "Out of memory");
    return NULL;
}

static PyObject*
Aligner_smithwaterman_align(Aligner* self, const char* sA, Py_ssize_t nA,
                                           const char* sB, Py_ssize_t nB)
{
    char c;
    int i;
    int j;
    int im = nA;
    int jm = nB;
    int kA;
    int kB;
    const double gap_extend_A = self->target_extend_gap_score;
    const double gap_extend_B = self->query_extend_gap_score;
    const double epsilon = self->epsilon;
    Trace** M = NULL;
    double maximum = 0;
    double score = 0;
    double* scores = NULL;
    double temp;
    int trace;
    PathGenerator* paths = NULL;

    /* Smith-Waterman algorithm */
    M = PyMem_Malloc((nA+1)*sizeof(Trace*));
    if (!M) goto exit;
    for (i = 0; i <= nA; i++) {
        M[i] = PyMem_Malloc((nB+1)*sizeof(Trace));
        if (!M[i]) goto exit;
    }
    scores = malloc((nB+1)*sizeof(double));
    if (!scores) goto exit;
    M[0][0].path = 0;
    for (j = 0; j <= nB; j++) {
        M[0][j].trace = 0;
        scores[j] = 0;
    }
    for (i = 1; i < nA; i++) {
        temp = 0;
        M[i][0].trace = 0;
        kA = CHARINDEX(sA[i-1]);
        for (j = 1; j < nB; j++) {
            kB = CHARINDEX(sB[j-1]);
            SELECT_TRACE_SMITH_WATERMAN_HVD(gap_extend_A, gap_extend_B);
        }
        kB = CHARINDEX(sB[nB-1]);
        SELECT_TRACE_SMITH_WATERMAN_D;
    }
    M[nA][0].trace = 0;
    temp = 0;
    kA = CHARINDEX(sA[nA-1]);
    for (j = 1; j < nB; j++) {
        kB = CHARINDEX(sB[j-1]);
        SELECT_TRACE_SMITH_WATERMAN_D;
    }
    kB = CHARINDEX(sB[nB-1]);
    SELECT_TRACE_SMITH_WATERMAN_D;
    PyMem_Free(scores);

    paths = _create_path_generator(self, nA, nB, epsilon);
    if (paths) {
        PyObject* result;
        paths->M.trace1 = M;
        if (maximum==0) M[0][0].path = NONE;
        result = Py_BuildValue("fO", maximum, paths);
        Py_DECREF(paths);
        return result;
    }

exit:
    if (M) {
        for (i = 0; i <= nA; i++) {
            if (!M[i]) break;
            PyMem_Free(M[i]);
        }
        PyMem_Free(M);
    }
    PyErr_SetString(PyExc_MemoryError, "Out of memory");
    return NULL;
}

static PyObject*
Aligner_gotoh_global_score(Aligner* self, const char* sA, Py_ssize_t nA,
                                          const char* sB, Py_ssize_t nB)
{
    char c;
    int i;
    int j;
    int kA;
    int kB;
    const double gap_open_A = self->target_open_gap_score;
    const double gap_open_B = self->query_open_gap_score;
    const double gap_extend_A = self->target_extend_gap_score;
    const double gap_extend_B = self->query_extend_gap_score;
    const double left_gap_open_A = self->target_left_open_gap_score;
    const double left_gap_open_B = self->query_left_open_gap_score;
    const double left_gap_extend_A = self->target_left_extend_gap_score;
    const double left_gap_extend_B = self->query_left_extend_gap_score;
    const double right_gap_open_A = self->target_right_open_gap_score;
    const double right_gap_open_B = self->query_right_open_gap_score;
    const double right_gap_extend_A = self->target_right_extend_gap_score;
    const double right_gap_extend_B = self->query_right_extend_gap_score;
    double* M_scores = NULL;
    double* Ix_scores = NULL;
    double* Iy_scores = NULL;
    double score;
    double temp;
    double M_temp;
    double Ix_temp;
    double Iy_temp;
    PyObject* result = NULL;

    /* Gotoh algorithm with three states */
    M_scores = PyMem_Malloc((nB+1)*sizeof(double));
    if (!M_scores) goto exit;
    Ix_scores = PyMem_Malloc((nB+1)*sizeof(double));
    if (!Ix_scores) goto exit;
    Iy_scores = PyMem_Malloc((nB+1)*sizeof(double));
    if (!Iy_scores) goto exit;

    /* The top row of the score matrix is a special case,
     * as there are no previously aligned characters.
     */
    M_scores[0] = 0;
    Ix_scores[0] = -DBL_MAX;
    Iy_scores[0] = -DBL_MAX;
    for (j = 1; j <= nB; j++) {
        M_scores[j] = -DBL_MAX;
        Ix_scores[j] = -DBL_MAX;
        Iy_scores[j] = left_gap_open_A + left_gap_extend_A * (j-1);
    }

    for (i = 1; i < nA; i++) {
        M_temp = M_scores[0];
        Ix_temp = Ix_scores[0];
        Iy_temp = Iy_scores[0];
        M_scores[0] = -DBL_MAX;
        Ix_scores[0] = left_gap_open_B + left_gap_extend_B * (i-1);
        Iy_scores[0] = -DBL_MAX;
        kA = CHARINDEX(sA[i-1]);
        for (j = 1; j < nB; j++) {
            kB = CHARINDEX(sB[j-1]);
            SELECT_SCORE_GLOBAL(M_temp,
                                Ix_temp,
                                Iy_temp);
            M_temp = M_scores[j];
            M_scores[j] = score + self->substitution_matrix[kA][kB];
            SELECT_SCORE_GLOBAL(M_temp + gap_open_B,
                                Ix_scores[j] + gap_extend_B,
                                Iy_scores[j] + gap_open_B);
            Ix_temp = Ix_scores[j];
            Ix_scores[j] = score;
            SELECT_SCORE_GLOBAL(M_scores[j-1] + gap_open_A,
                                Ix_scores[j-1] + gap_open_A,
                                Iy_scores[j-1] + gap_extend_A);
            Iy_temp = Iy_scores[j];
            Iy_scores[j] = score;
        }
        kB = CHARINDEX(sB[nB-1]);
        SELECT_SCORE_GLOBAL(M_temp,
                            Ix_temp,
                            Iy_temp);
        M_temp = M_scores[nB];
        M_scores[nB] = score + self->substitution_matrix[kA][kB];
        SELECT_SCORE_GLOBAL(M_temp + right_gap_open_B,
                            Ix_scores[nB] + right_gap_extend_B,
                            Iy_scores[nB] + right_gap_open_B);
        Ix_scores[nB] = score;
        SELECT_SCORE_GLOBAL(M_scores[nB-1] + gap_open_A,
                            Iy_scores[nB-1] + gap_extend_A,
                            Ix_scores[nB-1] + gap_open_A);
        Iy_scores[nB] = score;
    }

    M_temp = M_scores[0];
    Ix_temp = Ix_scores[0];
    Iy_temp = Iy_scores[0];
    M_scores[0] = -DBL_MAX;
    Ix_scores[0] = left_gap_open_B + left_gap_extend_B * (i-1);
    Iy_scores[0] = -DBL_MAX;
    kA = CHARINDEX(sA[nA-1]);
    for (j = 1; j < nB; j++) {
        kB = CHARINDEX(sB[j-1]);
        SELECT_SCORE_GLOBAL(M_temp,
                            Ix_temp,
                            Iy_temp);
        M_temp = M_scores[j];
        M_scores[j] = score + self->substitution_matrix[kA][kB];
        SELECT_SCORE_GLOBAL(M_temp + gap_open_B,
                            Ix_scores[j] + gap_extend_B,
                            Iy_scores[j] + gap_open_B);
        Ix_temp = Ix_scores[j];
        Ix_scores[j] = score;
        SELECT_SCORE_GLOBAL(M_scores[j-1] + right_gap_open_A,
                            Iy_scores[j-1] + right_gap_extend_A,
                            Ix_scores[j-1] + right_gap_open_A);
        Iy_temp = Iy_scores[j];
        Iy_scores[j] = score;
    }

    kB = CHARINDEX(sB[nB-1]);
    SELECT_SCORE_GLOBAL(M_temp,
                        Ix_temp,
                        Iy_temp);
    M_temp = M_scores[nB];
    M_scores[nB] = score + self->substitution_matrix[kA][kB];
    SELECT_SCORE_GLOBAL(M_temp + right_gap_open_B,
                        Ix_scores[nB] + right_gap_extend_B,
                        Iy_scores[nB] + right_gap_open_B);
    Ix_temp = Ix_scores[nB];
    Ix_scores[nB] = score;
    SELECT_SCORE_GLOBAL(M_scores[nB-1] + right_gap_open_A,
                        Ix_scores[nB-1] + right_gap_open_A,
                        Iy_scores[nB-1] + right_gap_extend_A);
    Iy_temp = Iy_scores[nB];
    Iy_scores[nB] = score;

    SELECT_SCORE_GLOBAL(M_scores[nB], Ix_scores[nB], Iy_scores[nB]);
    result = PyFloat_FromDouble(score);

exit:
    if (M_scores) PyMem_Free(M_scores);
    if (Ix_scores) PyMem_Free(Ix_scores);
    if (Iy_scores) PyMem_Free(Iy_scores);
    if (!result) PyErr_SetString(PyExc_MemoryError, "Out of memory");
    return result;
}

static PyObject*
Aligner_gotoh_local_score(Aligner* self, const char* sA, Py_ssize_t nA,
                                         const char* sB, Py_ssize_t nB)
{
    char c;
    int i;
    int j;
    int kA;
    int kB;
    const double gap_open_A = self->target_open_gap_score;
    const double gap_open_B = self->query_open_gap_score;
    const double gap_extend_A = self->target_extend_gap_score;
    const double gap_extend_B = self->query_extend_gap_score;
    double* M_scores = NULL;
    double* Ix_scores = NULL;
    double* Iy_scores = NULL;
    double score;
    double temp;
    double M_temp;
    double Ix_temp;
    double Iy_temp;
    double maximum = 0.0;
    PyObject* result = NULL;

    /* Gotoh algorithm with three states */
    M_scores = PyMem_Malloc((nB+1)*sizeof(double));
    if (!M_scores) goto exit;
    Ix_scores = PyMem_Malloc((nB+1)*sizeof(double));
    if (!Ix_scores) goto exit;
    Iy_scores = PyMem_Malloc((nB+1)*sizeof(double));
    if (!Iy_scores) goto exit;

    /* The top row of the score matrix is a special case,
     * as there are no previously aligned characters.
     */
    M_scores[0] = 0;
    Ix_scores[0] = -DBL_MAX;
    Iy_scores[0] = -DBL_MAX;
    for (j = 1; j <= nB; j++) {
        M_scores[j] = -DBL_MAX;
        Ix_scores[j] = -DBL_MAX;
        Iy_scores[j] = 0;
    }

    for (i = 1; i < nA; i++) {
        M_temp = M_scores[0];
        Ix_temp = Ix_scores[0];
        Iy_temp = Iy_scores[0];
        M_scores[0] = -DBL_MAX;
        Ix_scores[0] = 0;
        Iy_scores[0] = -DBL_MAX;
        kA = CHARINDEX(sA[i-1]);
        for (j = 1; j < nB; j++) {
            kB = CHARINDEX(sB[j-1]);
            SELECT_SCORE_GOTOH_LOCAL_ALIGN(M_temp,
                                           Ix_temp,
                                           Iy_temp,
                                           self->substitution_matrix[kA][kB]);
            M_temp = M_scores[j];
            M_scores[j] = score;
            SELECT_SCORE_LOCAL3(M_temp + gap_open_B,
                                Ix_scores[j] + gap_extend_B,
                                Iy_scores[j] + gap_open_B);
            Ix_temp = Ix_scores[j];
            Ix_scores[j] = score;
            SELECT_SCORE_LOCAL3(M_scores[j-1] + gap_open_A,
                                Ix_scores[j-1] + gap_open_A,
                                Iy_scores[j-1] + gap_extend_A);
            Iy_temp = Iy_scores[j];
            Iy_scores[j] = score;
        }

        kB = CHARINDEX(sB[nB-1]);

        Ix_scores[nB] = 0;
        Iy_scores[nB] = 0;
        SELECT_SCORE_GOTOH_LOCAL_ALIGN(M_temp,
                                       Ix_temp,
                                       Iy_temp,
                                       self->substitution_matrix[kA][kB]);
        M_temp = M_scores[nB];
        M_scores[nB] = score;
    }

    M_temp = M_scores[0];
    Ix_temp = Ix_scores[0];
    Iy_temp = Iy_scores[0];
    M_scores[0] = -DBL_MAX;
    Ix_scores[0] = 0;
    Iy_scores[0] = -DBL_MAX;
    kA = CHARINDEX(sA[nA-1]);
    for (j = 1; j < nB; j++) {
        kB = CHARINDEX(sB[j-1]);
        SELECT_SCORE_GOTOH_LOCAL_ALIGN(M_temp,
                                       Ix_temp,
                                       Iy_temp,
                                       self->substitution_matrix[kA][kB]);
        M_temp = M_scores[j];
        M_scores[j] = score;
        Ix_temp = Ix_scores[j];
        Iy_temp = Iy_scores[j];
        Ix_scores[j] = 0;
        Iy_scores[j] = 0;
    }

    kB = CHARINDEX(sB[nB-1]);
    SELECT_SCORE_GOTOH_LOCAL_ALIGN(M_temp,
                                   Ix_temp,
                                   Iy_temp,
                                   self->substitution_matrix[kA][kB]);

    result = PyFloat_FromDouble(maximum);

exit:
    if (M_scores) PyMem_Free(M_scores);
    if (Ix_scores) PyMem_Free(Ix_scores);
    if (Iy_scores) PyMem_Free(Iy_scores);
    if (!result) PyErr_SetString(PyExc_MemoryError, "Out of memory");
    return result;
}

static PyObject*
Aligner_gotoh_global_align(Aligner* self, const char* sA, Py_ssize_t nA,
                                          const char* sB, Py_ssize_t nB)
{
    char c;
    int i;
    int j;
    int kA;
    int kB;
    const double gap_open_A = self->target_open_gap_score;
    const double gap_open_B = self->query_open_gap_score;
    const double gap_extend_A = self->target_extend_gap_score;
    const double gap_extend_B = self->query_extend_gap_score;
    const double left_gap_open_A = self->target_left_open_gap_score;
    const double left_gap_open_B = self->query_left_open_gap_score;
    const double left_gap_extend_A = self->target_left_extend_gap_score;
    const double left_gap_extend_B = self->query_left_extend_gap_score;
    const double right_gap_open_A = self->target_right_open_gap_score;
    const double right_gap_open_B = self->query_right_open_gap_score;
    const double right_gap_extend_A = self->target_right_extend_gap_score;
    const double right_gap_extend_B = self->query_right_extend_gap_score;
    const double epsilon = self->epsilon;
    Trace3** trace3 = NULL;
    double* M_scores = NULL;
    double* Ix_scores = NULL;
    double* Iy_scores = NULL;
    double score;
    int trace;
    double temp;
    double M_temp;
    double Ix_temp;
    double Iy_temp;

    PathGenerator* paths = NULL;

    M_scores = PyMem_Malloc((nB+1)*sizeof(double));
    if (!M_scores) goto exit;
    Ix_scores = PyMem_Malloc((nB+1)*sizeof(double));
    if (!Ix_scores) goto exit;
    Iy_scores = PyMem_Malloc((nB+1)*sizeof(double));
    if (!Iy_scores) goto exit;
    trace3 = PyMem_Malloc((nA+1)*sizeof(Trace3*));
    if (!trace3) goto exit;
    for (i = 0; i <= nA; i++) {
        trace3[i] = PyMem_Malloc((nB+1)*sizeof(Trace3));
        if (!trace3[i]) goto exit;
    }

    /* Gotoh algorithm with three states */
    M_scores[0] = 0;
    trace3[0][0].M = 0;
    Ix_scores[0] = -DBL_MAX;
    trace3[0][0].Ix = 0;
    Iy_scores[0] = -DBL_MAX;
    trace3[0][0].Iy = 0;
    for (i = 1; i <= nA; i++) {
        trace3[i][0].M = 0;
        trace3[i][0].Ix = Ix_MATRIX;
        trace3[i][0].Iy = 0;
    }
    trace3[1][0].Ix = M_MATRIX;

    for (j = 1; j <= nB; j++) {
        M_scores[j] = -DBL_MAX;
        trace3[0][j].M = 0;
        Ix_scores[j] = -DBL_MAX;
        trace3[0][j].Ix = 0;
        Iy_scores[j] = left_gap_open_A + left_gap_extend_A * (j-1);
        trace3[0][j].Iy = Iy_MATRIX;
    }
    trace3[0][1].Iy = M_MATRIX;

    for (i = 1; i < nA; i++) {
        kA = CHARINDEX(sA[i-1]);
        M_temp = M_scores[0];
        Ix_temp = Ix_scores[0];
        Iy_temp = Iy_scores[0];
        M_scores[0] = -DBL_MAX;
        Ix_scores[0] = left_gap_open_B + left_gap_extend_B * (i-1);
        Iy_scores[0] = -DBL_MAX;
        for (j = 1; j < nB; j++) {
            kB = CHARINDEX(sB[j-1]);
            SELECT_TRACE_GOTOH_GLOBAL_ALIGN;
            M_temp = M_scores[j];
            M_scores[j] = score + self->substitution_matrix[kA][kB];
            SELECT_TRACE_GOTOH_GLOBAL_GAP(Ix,
                                          M_temp + gap_open_B,
                                          Ix_scores[j] + gap_extend_B,
                                          Iy_scores[j] + gap_open_B);
            Ix_temp = Ix_scores[j];
            Ix_scores[j] = score;
            SELECT_TRACE_GOTOH_GLOBAL_GAP(Iy,
                                          M_scores[j-1] + gap_open_A,
                                          Ix_scores[j-1] + gap_open_A,
                                          Iy_scores[j-1] + gap_extend_A);
            Iy_temp = Iy_scores[j];
            Iy_scores[j] = score;
        }
        kB = CHARINDEX(sB[nB-1]);
        SELECT_TRACE_GOTOH_GLOBAL_ALIGN;
        M_temp = M_scores[nB];
        M_scores[nB] = score + self->substitution_matrix[kA][kB];
        SELECT_TRACE_GOTOH_GLOBAL_GAP(Ix,
                                      M_temp + right_gap_open_B,
                                      Ix_scores[nB] + right_gap_extend_B,
                                      Iy_scores[nB] + right_gap_open_B);
        Ix_temp = Ix_scores[nB];
        Ix_scores[nB] = score;
        SELECT_TRACE_GOTOH_GLOBAL_GAP(Iy,
                                      M_scores[nB-1] + gap_open_A,
                                      Ix_scores[nB-1] + gap_open_A,
                                      Iy_scores[nB-1] + gap_extend_A);
        Iy_temp = Iy_scores[nB];
        Iy_scores[nB] = score;
    }
    kA = CHARINDEX(sA[nA-1]);
    M_temp = M_scores[0];
    Ix_temp = Ix_scores[0];
    Iy_temp = Iy_scores[0];
    M_scores[0] = -DBL_MAX;
    Ix_scores[0] = left_gap_open_B + left_gap_extend_B * (nA-1);
    Iy_scores[0] = -DBL_MAX;
    for (j = 1; j < nB; j++) {
        kB = CHARINDEX(sB[j-1]);
        SELECT_TRACE_GOTOH_GLOBAL_ALIGN;
        M_temp = M_scores[j];
        M_scores[j] = score + self->substitution_matrix[kA][kB];
        SELECT_TRACE_GOTOH_GLOBAL_GAP(Ix,
                                      M_scores[j] + gap_open_B,
                                      Ix_scores[j] + gap_extend_B,
                                      Iy_scores[j] + gap_open_B);
        Ix_temp = Ix_scores[j];
        Ix_scores[j] = score;
        SELECT_TRACE_GOTOH_GLOBAL_GAP(Iy,
                                      M_scores[j-1] + right_gap_open_A,
                                      Ix_scores[j-1] + right_gap_open_A,
                                      Iy_scores[j-1] + right_gap_extend_A);
        Iy_temp = Iy_scores[j];
        Iy_scores[j] = score;
    }
    kB = CHARINDEX(sB[nB-1]);
    SELECT_TRACE_GOTOH_GLOBAL_ALIGN;
    M_temp = M_scores[j];
    M_scores[j] = score + self->substitution_matrix[kA][kB];
    SELECT_TRACE_GOTOH_GLOBAL_GAP(Ix,
                                  M_temp + right_gap_open_B,
                                  Ix_scores[j] + right_gap_extend_B,
                                  Iy_scores[j] + right_gap_open_B);
    Ix_scores[nB] = score;
    SELECT_TRACE_GOTOH_GLOBAL_GAP(Iy,
                                  M_scores[j-1] + right_gap_open_A,
                                  Ix_scores[j-1] + right_gap_open_A,
                                  Iy_scores[j-1] + right_gap_extend_A);
    Iy_scores[nB] = score;
    trace3[0][0].path = 0;
    trace3[nA][nB].path = 0;

    /* traceback */
    paths = _create_path_generator(self, nA, nB, epsilon);
    if (paths) {
        PyObject* result;
        SELECT_SCORE_GLOBAL(M_scores[nB],
                            Ix_scores[nB],
                            Iy_scores[nB]);
        paths->trace3 = trace3;
        result = Py_BuildValue("fO", score, paths);
        Py_DECREF(paths);
        score -= epsilon;
        if (M_scores[nB] < score) trace3[nA][nB].M = 0;
        if (Ix_scores[nB] < score) trace3[nA][nB].Ix = 0;
        if (Iy_scores[nB] < score) trace3[nA][nB].Iy = 0;
        return result;
    }
exit:
    if (trace3) {
        for (i = 0; i <= nA; i++) {
            if (!trace3[i]) break;
            PyMem_Free(trace3[i]);
        }
        PyMem_Free(trace3);
    }
    if (M_scores) PyMem_Free(M_scores);
    if (Ix_scores) PyMem_Free(Ix_scores);
    if (Iy_scores) PyMem_Free(Iy_scores);
    PyErr_SetString(PyExc_MemoryError, "Out of memory");
    return NULL;
}

static PyObject*
Aligner_gotoh_local_align(Aligner* self, const char* sA, Py_ssize_t nA,
                                         const char* sB, Py_ssize_t nB)
{
    char c;
    int i;
    int j;
    int im = nA;
    int jm = nB;
    int kA;
    int kB;
    const double gap_open_A = self->target_open_gap_score;
    const double gap_open_B = self->query_open_gap_score;
    const double gap_extend_A = self->target_extend_gap_score;
    const double gap_extend_B = self->query_extend_gap_score;
    const double epsilon = self->epsilon;
    Trace3** trace3 = NULL;
    double* M_scores = NULL;
    double* Ix_scores = NULL;
    double* Iy_scores = NULL;
    double score;
    int trace;
    double temp;
    double M_temp;
    double Ix_temp;
    double Iy_temp;
    double maximum = 0.0;

    PathGenerator* paths = NULL;

    /* Gotoh algorithm with three states */
    M_scores = PyMem_Malloc((nB+1)*sizeof(double));
    if (!M_scores) goto exit;
    Ix_scores = PyMem_Malloc((nB+1)*sizeof(double));
    if (!Ix_scores) goto exit;
    Iy_scores = PyMem_Malloc((nB+1)*sizeof(double));
    if (!Iy_scores) goto exit;
    trace3 = PyMem_Malloc((nA+1)*sizeof(Trace3*));
    if (!trace3) goto exit;
    for (i = 0; i <= nA; i++) {
        trace3[i] = PyMem_Malloc((nB+1)*sizeof(Trace3));
        if (!trace3[i]) goto exit;
    }

    M_scores[0] = 0;
    Ix_scores[0] = -DBL_MAX;
    Iy_scores[0] = -DBL_MAX;
    trace3[0][0].M = 0;
    trace3[0][0].Ix = 0;
    trace3[0][0].Iy = 0;
    trace3[0][0].path = 0;

    for (j = 1; j <= nB; j++) {
        M_scores[j] = 0;
        trace3[0][j].M = 0;
        Ix_scores[j] = -DBL_MAX;
        trace3[0][j].Ix = 0;
        Iy_scores[j] = -DBL_MAX;
        trace3[0][j].Iy = 0;
    }
    for (i = 1; i < nA; i++) {
        M_temp = M_scores[0];
        Ix_temp = Ix_scores[0];
        Iy_temp = Iy_scores[0];
        M_scores[0] = 0;
        Ix_scores[0] = -DBL_MAX;
        Iy_scores[0] = -DBL_MAX;
        trace3[i][0].M = 0;
        trace3[i][0].Ix = 0;
        trace3[i][0].Iy = 0;
        kA = CHARINDEX(sA[i-1]);
        for (j = 1; j < nB; j++) {
            kB = CHARINDEX(sB[j-1]);
            SELECT_TRACE_GOTOH_LOCAL_ALIGN
            M_temp = M_scores[j];
            M_scores[j] = score;
            SELECT_TRACE_GOTOH_LOCAL_GAP(Ix,
                                     M_temp + gap_open_B,
                                     Ix_scores[j] + gap_extend_B,
                                     Iy_scores[j] + gap_open_B);
            Ix_temp = Ix_scores[j];
            Ix_scores[j] = score;
            SELECT_TRACE_GOTOH_LOCAL_GAP(Iy,
                                     M_scores[j-1] + gap_open_A,
                                     Ix_scores[j-1] + gap_open_A,
                                     Iy_scores[j-1] + gap_extend_A);
            Iy_temp = Iy_scores[j];
            Iy_scores[j] = score;
        }
        kB = CHARINDEX(sB[nB-1]);
        SELECT_TRACE_GOTOH_LOCAL_ALIGN
        M_temp = M_scores[j];
        M_scores[j] = score;
        Ix_temp = Ix_scores[nB];
        Ix_scores[nB] = 0;
        trace3[i][nB].Ix = 0;
        Iy_temp = Iy_scores[nB];
        Iy_scores[nB] = 0;
        trace3[i][nB].Iy = 0;
    }
    M_temp = M_scores[0];
    M_scores[0] = 0;
    trace3[nA][0].M = 0;
    Ix_temp = Ix_scores[0];
    Ix_scores[0] = -DBL_MAX;
    trace3[nA][0].Ix = 0;
    trace3[nA][0].Iy = 0;
    Iy_temp = Iy_scores[0];
    Iy_scores[0] = -DBL_MAX;
    kA = CHARINDEX(sA[nA-1]);
    for (j = 1; j < nB; j++) {
        kB = CHARINDEX(sB[j-1]);
        SELECT_TRACE_GOTOH_LOCAL_ALIGN
        M_temp = M_scores[j];
        M_scores[j] = score;
        Ix_temp = Ix_scores[j];
        Ix_scores[j] = 0;
        trace3[nA][j].Ix = 0;
        Iy_temp = Iy_scores[j];
        Iy_scores[j] = 0;
        trace3[nA][j].Iy = 0;
    }
    kB = CHARINDEX(sB[nB-1]);
    SELECT_TRACE_GOTOH_LOCAL_ALIGN
    trace3[nA][nB].Ix = 0;
    Ix_scores[nB] = 0;
    trace3[nA][nB].Iy = 0;
    Iy_scores[nB] = 0;

    /* traceback */
    paths = _create_path_generator(self, nA, nB, epsilon);
    if (paths) {
        PyObject* result;
        paths->trace3 = trace3;
        if (maximum==0) trace3[0][0].path = DONE;
        result = Py_BuildValue("fO", maximum, paths);
        Py_DECREF(paths);
        return result;
    }
exit:
    if (trace3) {
        for (i = 0; i <= nA; i++) {
            if (!trace3[i]) break;
            PyMem_Free(trace3[i]);
        }
        PyMem_Free(trace3);
    }
    if (M_scores) PyMem_Free(M_scores);
    if (Ix_scores) PyMem_Free(Ix_scores);
    if (Iy_scores) PyMem_Free(Iy_scores);
    PyErr_SetString(PyExc_MemoryError, "Out of memory");
    return NULL;
}

static int
_call_query_gap_function(Aligner* aligner, int i, int j, double* score)
{
    double value;
    PyObject* result;
    PyObject* function = aligner->query_gap_function;
    if (!function)
        value = aligner->query_open_gap_score
              + (j-1) * aligner->query_extend_gap_score;
    else {
        result = PyObject_CallFunction(function, "ii", i, j);
        if (result == NULL) return 0;
        value = PyFloat_AsDouble(result);
        Py_DECREF(result);
        if (value == -1.0 && PyErr_Occurred()) return 0;
    }
    *score = value;
    return 1;
}

static int
_call_target_gap_function(Aligner* aligner, int i, int j, double* score)
{
    double value;
    PyObject* result;
    PyObject* function = aligner->target_gap_function;
    if (!function)
        value = aligner->target_open_gap_score
              + (j-1) * aligner->target_extend_gap_score;
    else {
        result = PyObject_CallFunction(function, "ii", i, j);
        if (result == NULL) return 0;
        value = PyFloat_AsDouble(result);
        Py_DECREF(result);
        if (value == -1.0 && PyErr_Occurred()) return 0;
    }
    *score = value;
    return 1;
}

static PyObject*
Aligner_waterman_smith_beyer_global_score(Aligner* self,
                                          const char* sA, Py_ssize_t nA,
                                          const char* sB, Py_ssize_t nB)
{
    char c;
    int i;
    int j;
    int k;
    int kA;
    int kB;
    double** M = NULL;
    double** Ix = NULL;
    double** Iy = NULL;
    double score = 0.0;
    double gapscore;
    double temp;
    int ok = 1;

    PyObject* result = NULL;

    /* Waterman-Smith-Beyer algorithm */
    M = PyMem_Malloc((nA+1)*sizeof(double*));
    if (!M) goto exit;
    Ix = PyMem_Malloc((nA+1)*sizeof(double*));
    if (!Ix) goto exit;
    Iy = PyMem_Malloc((nA+1)*sizeof(double*));
    if (!Iy) goto exit;
    for (i = 0; i <= nA; i++) {
        M[i] = PyMem_Malloc((nB+1)*sizeof(double));
        if (!M[i]) goto exit;
        Ix[i] = PyMem_Malloc((nB+1)*sizeof(double));
        if (!Ix[i]) goto exit;
        Iy[i] = PyMem_Malloc((nB+1)*sizeof(double));
        if (!Iy[i]) goto exit;
    }

    /* The top row of the score matrix is a special case,
     *  as there are no previously aligned characters.
     */
    M[0][0] = 0;
    Ix[0][0] = -DBL_MAX;
    Iy[0][0] = -DBL_MAX;
    for (i = 1; i <= nA; i++) {
        ok = _call_query_gap_function(self, 0, i, &score);
        if (!ok) goto exit;
        M[i][0] = -DBL_MAX;
        Ix[i][0] = score;
        Iy[i][0] = -DBL_MAX;
    }
    for (j = 1; j <= nB; j++) {
        ok = _call_target_gap_function(self, 0, j, &score);
        if (!ok) goto exit;
        M[0][j] = -DBL_MAX;
        Ix[0][j] = -DBL_MAX;
        Iy[0][j] = score;
    }
    for (i = 1; i <= nA; i++) {
        kA = CHARINDEX(sA[i-1]);
        for (j = 1; j <= nB; j++) {
            kB = CHARINDEX(sB[j-1]);
            SELECT_SCORE_GLOBAL(M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1]);
            M[i][j] = score + self->substitution_matrix[kA][kB];
            score = -DBL_MAX;
            for (k = 1; k <= i; k++) {
                ok = _call_query_gap_function(self, j, k, &gapscore);
                if (!ok) goto exit;
                SELECT_SCORE_WATERMAN_SMITH_BEYER(M[i-k][j], Iy[i-k][j]);
            }
            Ix[i][j] = score;
            score = -DBL_MAX;
            for (k = 1; k <= j; k++) {
                ok = _call_target_gap_function(self, i, k, &gapscore);
                if (!ok) goto exit;
                SELECT_SCORE_WATERMAN_SMITH_BEYER(M[i][j-k], Ix[i][j-k]);
            }
            Iy[i][j] = score;
        }
    }
    SELECT_SCORE_GLOBAL(M[nA][nB], Ix[nA][nB], Iy[nA][nB]);

    result = PyFloat_FromDouble(score);
exit:
    if (M) {
        /* If M is NULL, then Ix is also NULL. */
        if (Ix) {
            /* If Ix is NULL, then Iy is also NULL. */
            if (Iy) {
                /* If Iy is NULL, then M[i], Ix[i], and Iy[i] are also NULL. */ 
                for (i = 0; i <= nA; i++) {
                    if (!M[i]) break;
                    PyMem_Free(M[i]); 
                    if (!Ix[i]) break;
                    PyMem_Free(Ix[i]);
                    if (!Iy[i]) break;
                    PyMem_Free(Iy[i]);
                }
                PyMem_Free(Iy);
            }
            PyMem_Free(Ix);
        }
        PyMem_Free(M);
    }
    if (!ok) return NULL;
    if (!result) PyErr_SetString(PyExc_MemoryError, "Out of memory");
    return result;
}

static PyObject*
Aligner_waterman_smith_beyer_global_align(Aligner* self,
                                          const char* sA, Py_ssize_t nA,
                                          const char* sB, Py_ssize_t nB)
{
    char c;
    int i;
    int j;
    int k;
    int kA;
    int kB;
    const double epsilon = self->epsilon;
    CellM** M = NULL;
    CellXY** Ix = NULL;
    CellXY** Iy = NULL;
    int ng;
    int nm;
    double score;
    double gapscore;
    double temp;
    int trace;
    int* traceM;
    int* traceXY;
    int ok = 1;

    PathGenerator* paths = NULL;

    /* Waterman-Smith-Beyer algorithm */
    if (!_allocate_watermansmithbeyer_matrices(nA, nB, &M, &Ix, &Iy, Global))
        return NULL;
    for (i = 1; i <= nA; i++) {
        ok = _call_query_gap_function(self, 0, i, &score);
        if (!ok) goto exit;
        Ix[i][0].score = score;
    }
    for (j = 1; j <= nB; j++) {
        ok = _call_target_gap_function(self, 0, j, &score);
        if (!ok) goto exit;
        Iy[0][j].score = score;
    }
    for (i = 1; i <= nA; i++) {
        kA = CHARINDEX(sA[i-1]);
        for (j = 1; j <= nB; j++) {
            kB = CHARINDEX(sB[j-1]);
            SELECT_TRACE_WATERMAN_SMITH_BEYER_GLOBAL_ALIGN(self->substitution_matrix[kA][kB]);
            traceM = PyMem_Malloc((i+1)*sizeof(int));
            if (!traceM) goto exit;
            Ix[i][j].traceM = traceM;
            traceXY = PyMem_Malloc((i+1)*sizeof(int));
            if (!traceXY) goto exit;
            Ix[i][j].traceXY = traceXY;
            nm = 0;
            ng = 0;
            score = -DBL_MAX;
            for (k = 1; k <= i; k++) {
                ok = _call_query_gap_function(self, j, k, &gapscore);
                if (!ok) goto exit;
                SELECT_TRACE_WATERMAN_SMITH_BEYER_GAP(M[i-k][j].score,
                                                      Iy[i-k][j].score,
                                                      i-k);
            }
            traceM = PyMem_Realloc(traceM, (nm+1)*sizeof(int));
            if (!traceM) goto exit;
            Ix[i][j].traceM = traceM;
            traceM[nm] = -1;
            traceXY = PyMem_Realloc(traceXY, (ng+1)*sizeof(int));
            if (!traceXY) goto exit;
            Ix[i][j].traceXY = traceXY;
            traceXY[ng] = -1;
            Ix[i][j].score = score;
            traceM = PyMem_Malloc((j+1)*sizeof(int));
            if (!traceM) goto exit;
            Iy[i][j].traceM = traceM;
            traceXY = PyMem_Malloc((j+1)*sizeof(int));
            if (!traceXY) goto exit;
            Iy[i][j].traceXY = traceXY;
            nm = 0;
            ng = 0;
            score = -DBL_MAX;
            for (k = 1; k <= j; k++) {
                ok = _call_target_gap_function(self, i, k, &gapscore);
                if (!ok) goto exit;
                SELECT_TRACE_WATERMAN_SMITH_BEYER_GAP(M[i][j-k].score,
                                                      Ix[i][j-k].score,
                                                      j-k);
            }
            Iy[i][j].score = score;
            traceM = PyMem_Realloc(traceM, (nm+1)*sizeof(int));
            if (!traceM) goto exit;
            Iy[i][j].traceM = traceM;
            traceM[nm] = -1;
            traceXY = PyMem_Realloc(traceXY, (ng+1)*sizeof(int));
            if (!traceXY) goto exit;
            Iy[i][j].traceXY = traceXY;
            traceXY[ng] = -1;
        }
    }

    /* traceback */
    SELECT_SCORE_GLOBAL(M[nA][nB].score, Ix[nA][nB].score, Iy[nA][nB].score);
    M[nA][nB].path = 0;
    Ix[nA][nB].path = 0;
    Iy[nA][nB].path = 0;
    paths = _create_path_generator(self, nA, nB, epsilon);
    if (paths) {
        PyObject* result;
        paths->M.general = M;
        paths->Ix = Ix;
        paths->Iy = Iy;
        if (M[nA][nB].score < score - epsilon) M[nA][nB].trace = 0;
        if (Ix[nA][nB].score < score - epsilon) {
            traceM = PyMem_Realloc(Ix[nA][nB].traceM, sizeof(int));
            if (!traceM) goto exit;
            traceM[0] = -1;
            Ix[nA][nB].traceM = traceM;
            traceXY = PyMem_Realloc(Ix[nA][nB].traceXY, sizeof(int));
            if (!traceXY) goto exit;
            traceXY[0] = -1;
            Ix[nA][nB].traceXY = traceXY;
        }
        if (Iy[nA][nB].score < score - epsilon) {
            traceM = PyMem_Realloc(Iy[nA][nB].traceM, sizeof(int));
            if (!traceM) goto exit;
            traceM[0] = -1;
            Iy[nA][nB].traceM = traceM;
            traceXY = PyMem_Realloc(Iy[nA][nB].traceXY, sizeof(int));
            if (!traceXY) goto exit;
            traceXY[0] = -1;
            Iy[nA][nB].traceXY = traceXY;
        }
        result = Py_BuildValue("fO", score, paths);
        Py_DECREF(paths);
        score -= epsilon;
        return result;
    }

exit:
    if (ok) /* otherwise, an exception was already set */
        PyErr_SetString(PyExc_MemoryError, "Out of memory");
    _deallocate_watermansmithbeyer_matrices(nA, nB, M, Ix, Iy, Global);
    return NULL;
}

static PyObject*
Aligner_waterman_smith_beyer_local_score(Aligner* self,
                                         const char* sA, Py_ssize_t nA,
                                         const char* sB, Py_ssize_t nB)
{
    char c;
    int i;
    int j;
    int k;
    int kA;
    int kB;
    double** M = NULL;
    double** Ix = NULL;
    double** Iy = NULL;
    double score = 0.0;
    double gapscore = 0.0;
    double temp;
    int ok = 1;

    double maximum = 0.0;
    PyObject* result = NULL;

    /* Waterman-Smith-Beyer algorithm */
    M = PyMem_Malloc((nA+1)*sizeof(double*));
    if (!M) goto exit;
    Ix = PyMem_Malloc((nA+1)*sizeof(double*));
    if (!Ix) goto exit;
    Iy = PyMem_Malloc((nA+1)*sizeof(double*));
    if (!Iy) goto exit;
    for (i = 0; i <= nA; i++) {
        M[i] = PyMem_Malloc((nB+1)*sizeof(double));
        if (!M[i]) goto exit;
        Ix[i] = PyMem_Malloc((nB+1)*sizeof(double));
        if (!Ix[i]) goto exit;
        Iy[i] = PyMem_Malloc((nB+1)*sizeof(double));
        if (!Iy[i]) goto exit;
    }

    /* The top row of the score matrix is a special case,
     *  as there are no previously aligned characters.
     */
    M[0][0] = 0;
    Ix[0][0] = -DBL_MAX;
    Iy[0][0] = -DBL_MAX;
    for (i = 1; i <= nA; i++) {
        M[i][0] = -DBL_MAX;
        Ix[i][0] = 0;
        Iy[i][0] = -DBL_MAX;
    }
    for (j = 1; j <= nB; j++) {
        M[0][j] = -DBL_MAX;
        Ix[0][j] = -DBL_MAX;
        Iy[0][j] = 0;
    }
    for (i = 1; i <= nA; i++) {
        kA = CHARINDEX(sA[i-1]);
        for (j = 1; j <= nB; j++) {
            kB = CHARINDEX(sB[j-1]);
            SELECT_SCORE_GOTOH_LOCAL_ALIGN(M[i-1][j-1],
                                           Ix[i-1][j-1],
                                           Iy[i-1][j-1],
                                           self->substitution_matrix[kA][kB]);
            M[i][j] = score;
            if (i == nA || j == nB) {
                Ix[i][j] = 0;
                Iy[i][j] = 0;
                continue;
            }
            score = 0.0;
            for (k = 1; k <= i; k++) {
                ok = _call_query_gap_function(self, j, k, &gapscore);
                SELECT_SCORE_WATERMAN_SMITH_BEYER(M[i-k][j], Iy[i-k][j]);
                if (!ok) goto exit;
            }
            if (score > maximum) maximum = score;
            Ix[i][j] = score;
            score = 0.0;
            for (k = 1; k <= j; k++) {
                ok = _call_target_gap_function(self, i, k, &gapscore);
                if (!ok) goto exit;
                SELECT_SCORE_WATERMAN_SMITH_BEYER(M[i][j-k], Ix[i][j-k]);
            }
            if (score > maximum) maximum = score;
            Iy[i][j] = score;
        }
    }
    SELECT_SCORE_GLOBAL(M[nA][nB], Ix[nA][nB], Iy[nA][nB]);
    if (score > maximum) maximum = score;

    result = PyFloat_FromDouble(maximum);
exit:
    if (!result) {
        /* otherwise the alignment generator will PyMem_Free the matrices */
        if (M) {
            /* If M is NULL, then Ix is also NULL. */
            if (Ix) {
                /* If Ix is NULL, then Iy is also NULL. */
                if (Iy) {
                    /* If Iy is NULL, then M[i], Ix[i], and Iy[i] are
                     * also NULL. */
                    for (i = 0; i <= nA; i++) {
                        if (!M[i]) break;
                        PyMem_Free(M[i]);
                        if (!Ix[i]) break;
                        PyMem_Free(Ix[i]);
                        if (!Iy[i]) break;
                        PyMem_Free(Iy[i]);
                    }
                    PyMem_Free(Iy);
                }
                PyMem_Free(Ix);
            }
            PyMem_Free(M);
        }
    }
    if (!ok) return NULL;
    if (!result) PyErr_SetString(PyExc_MemoryError, "Out of memory");
    return result;
}

static PyObject*
Aligner_waterman_smith_beyer_local_align(Aligner* self,
                                         const char* sA, Py_ssize_t nA,
                                         const char* sB, Py_ssize_t nB)
{
    char c;
    int i;
    int j;
    int im = nA;
    int jm = nB;
    int k;
    int kA;
    int kB;
    const double epsilon = self->epsilon;
    CellM** M = NULL;
    CellXY** Ix = NULL;
    CellXY** Iy = NULL;
    double score;
    double gapscore;
    double temp;
    int trace;
    int* traceM;
    int* traceXY;
    int nm;
    int ng;
    int ok = 1;
    double maximum = 0;

    PathGenerator* paths = NULL;

    /* Waterman-Smith-Beyer algorithm */
    if (!_allocate_watermansmithbeyer_matrices(nA, nB, &M, &Ix, &Iy, Local))
        return NULL;

    for (i = 1; i <= nA; i++) {
        kA = CHARINDEX(sA[i-1]);
        for (j = 1; j <= nB; j++) {
            kB = CHARINDEX(sB[j-1]);
            nm = 0;
            ng = 0;
            SELECT_TRACE_WATERMAN_SMITH_BEYER_ALIGN(M[i][j],
                                           M[i-1][j-1].score,
                                           Ix[i-1][j-1].score,
                                           Iy[i-1][j-1].score,
                                           self->substitution_matrix[kA][kB]);
            M[i][j].path = 0;
            if (i == nA || j == nB) {
                Ix[i][j].score = score;
                Ix[i][j].traceM = NULL;
                Ix[i][j].traceXY = NULL;
                Iy[i][j].score = score;
                Iy[i][j].traceM = NULL;
                Iy[i][j].traceXY = NULL;
                continue;
            }
            traceM = PyMem_Malloc((i+1)*sizeof(int));
            if (!traceM) goto exit;
            Ix[i][j].traceM = traceM;
            traceXY = PyMem_Malloc((i+1)*sizeof(int));
            if (!traceXY) goto exit;
            Ix[i][j].traceXY = traceXY;
            score = -DBL_MAX;
            for (k = 1; k <= i; k++) {
                ok = _call_query_gap_function(self, j, k, &gapscore);
                if (!ok) goto exit;
                SELECT_TRACE_WATERMAN_SMITH_BEYER_GAP(M[i-k][j].score,
                                                      Iy[i-k][j].score,
                                                      i-k);
            }
            if (score < epsilon) {
                score = -DBL_MAX;
                nm = 0;
                ng = 0;
            }
            else if (score > maximum) maximum = score;
            traceM[nm] = -1;
            traceXY[ng] = -1;
            Ix[i][j].score = score;
            Ix[i][j].path = 0;
            traceM = PyMem_Realloc(traceM, (nm+1)*sizeof(int));
            if (!traceM) goto exit;
            Ix[i][j].traceM = traceM;
            traceM[nm] = -1;
            traceXY = PyMem_Realloc(traceXY, (ng+1)*sizeof(int));
            if (!traceXY) goto exit;
            Ix[i][j].traceXY = traceXY;
            traceXY[ng] = -1;
            traceM = PyMem_Malloc((j+1)*sizeof(int));
            if (!traceM) goto exit;
            Iy[i][j].traceM = traceM;
            traceXY = PyMem_Malloc((j+1)*sizeof(int));
            if (!traceXY) goto exit;
            Iy[i][j].traceXY = traceXY;
            nm = 0;
            ng = 0;
            score = -DBL_MAX;
            for (k = 1; k <= j; k++) {
                ok = _call_target_gap_function(self, i, k, &gapscore);
                if (!ok) goto exit;
                SELECT_TRACE_WATERMAN_SMITH_BEYER_GAP(M[i][j-k].score,
                                                      Ix[i][j-k].score,
                                                      j-k);
            }
            if (score < epsilon) {
                score = -DBL_MAX;
                nm = 0;
                ng = 0;
            }
            else if (score > maximum) maximum = score;
            traceM = PyMem_Realloc(traceM, (nm+1)*sizeof(int));
            if (!traceM) goto exit;
            Iy[i][j].traceM = traceM;
            traceXY = PyMem_Realloc(traceXY, (ng+1)*sizeof(int));
            if (!traceXY) goto exit;
            Iy[i][j].traceXY = traceXY;
            traceM[nm] = -1;
            traceXY[ng] = -1;
            Iy[i][j].score = score;
            Iy[i][j].path = 0;
        }
    }

    /* traceback */
    paths = _create_path_generator(self, nA, nB, epsilon);
    if (paths) {
        PyObject* result;
        paths->M.general = M;
        paths->Ix = Ix;
        paths->Iy = Iy;
        if (maximum==0) M[0][0].path = DONE;
        result = Py_BuildValue("fO", maximum, paths);
        Py_DECREF(paths);
        return result;
    }

exit:
    if (ok) /* otherwise, an exception was already set */
        PyErr_SetString(PyExc_MemoryError, "Out of memory");
    _deallocate_watermansmithbeyer_matrices(nA, nB, M, Ix, Iy, Local);
    return NULL;
}
 
static const char Aligner_score__doc__[] = "calculates the alignment score";

static PyObject*
Aligner_score(Aligner* self, PyObject* args, PyObject* keywords)
{
    const char* sA;
    const char* sB;
    Py_ssize_t nA;
    Py_ssize_t nB;
    const Mode mode = self->mode;
    const Algorithm algorithm = _get_algorithm(self);

    static char *kwlist[] = {"sequenceA", "sequenceB", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, keywords, "s#s#", kwlist,
                                    &sA, &nA, &sB, &nB))
        return NULL;

    switch (algorithm) {
        case NeedlemanWunschSmithWaterman:
            switch (mode) {
                case Global:
                    return Aligner_needlemanwunsch_score(self, sA, nA, sB, nB);
                case Local:
                    return Aligner_smithwaterman_score(self, sA, nA, sB, nB);
            }
        case Gotoh:
            switch (mode) {
                case Global:
                    return Aligner_gotoh_global_score(self, sA, nA, sB, nB);
                case Local:
                    return Aligner_gotoh_local_score(self, sA, nA, sB, nB);
            }
        case WatermanSmithBeyer:
            switch (mode) {
                case Global:
                    return Aligner_waterman_smith_beyer_global_score(self, sA, nA, sB, nB);
                case Local:
                    return Aligner_waterman_smith_beyer_local_score(self, sA, nA, sB, nB);
            }
        case Unknown:
        default:
            PyErr_SetString(PyExc_RuntimeError, "unknown algorithm");
            return NULL;
    }
}

static const char Aligner_align__doc__[] = "align two sequences";

static PyObject*
Aligner_align(Aligner* self, PyObject* args, PyObject* keywords)
{
    const char* sA;
    const char* sB;
    Py_ssize_t nA;
    Py_ssize_t nB;
    const Mode mode = self->mode;
    const Algorithm algorithm = _get_algorithm(self);
    static char *kwlist[] = {"sequenceA", "sequenceB", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, keywords, "s#s#", kwlist,
                                    &sA, &nA, &sB, &nB))
        return NULL;

    switch (algorithm) {
        case NeedlemanWunschSmithWaterman:
            switch (mode) {
                case Global:
                    return Aligner_needlemanwunsch_align(self, sA, nA, sB, nB);
                case Local:
                    return Aligner_smithwaterman_align(self, sA, nA, sB, nB);
            }
        case Gotoh:
            switch (mode) {
                case Global:
                    return Aligner_gotoh_global_align(self, sA, nA, sB, nB);
                case Local:
                    return Aligner_gotoh_local_align(self, sA, nA, sB, nB);
            }
        case WatermanSmithBeyer:
            switch (mode) {
                case Global:
                    return Aligner_waterman_smith_beyer_global_align(self, sA, nA, sB, nB);
                case Local:
                    return Aligner_waterman_smith_beyer_local_align(self, sA, nA, sB, nB);
            }
        case Unknown:
        default:
            PyErr_SetString(PyExc_RuntimeError, "unknown algorithm");
            return NULL;
    }
}

static char Aligner_doc[] =
"Aligner.\n";

static PyMethodDef Aligner_methods[] = {
    {"score",
     (PyCFunction)Aligner_score,
     METH_VARARGS | METH_KEYWORDS,
     Aligner_score__doc__
    },
    {"align",
     (PyCFunction)Aligner_align,
     METH_VARARGS | METH_KEYWORDS,
     Aligner_align__doc__
    },
    {NULL}  /* Sentinel */
};

static PyTypeObject AlignerType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_algorithms.PairwiseAligner", /* tp_name */
    sizeof(Aligner),               /* tp_basicsize */
    0,                             /* tp_itemsize */
    (destructor)Aligner_dealloc,   /* tp_dealloc */
    0,                             /* tp_print */
    0,                             /* tp_getattr */
    0,                             /* tp_setattr */
    0,                             /* tp_compare */
    (reprfunc)Aligner_repr,        /* tp_repr */
    0,                             /* tp_as_number */
    0,                             /* tp_as_sequence */
    0,                             /* tp_as_mapping */
    0,                             /* tp_hash */
    0,                             /* tp_call */
    (reprfunc)Aligner_str,         /* tp_str */
    0,                             /* tp_getattro */
    0,                             /* tp_setattro */
    0,                             /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    Aligner_doc,                   /* tp_doc */
    0,                             /* tp_traverse */
    0,                             /* tp_clear */
    0,                             /* tp_richcompare */
    0,                             /* tp_weaklistoffset */
    0,                             /* tp_iter */
    0,                             /* tp_iternext */
    Aligner_methods,               /* tp_methods */
    0,                             /* tp_members */
    Aligner_getset,                /* tp_getset */
    0,                             /* tp_base */
    0,                             /* tp_dict */
    0,                             /* tp_descr_get */
    0,                             /* tp_descr_set */
    0,                             /* tp_dictoffset */
    (initproc)Aligner_init,        /* tp_init */
};

/* Module definition */

static PyMethodDef _aligners_methods[] = {
    {NULL, NULL, 0, NULL}
};

static char _aligners__doc__[] =
"C extension module implementing pairwise alignment algorithms";

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_aligners",
        _aligners__doc__,
        -1,
        _aligners_methods,
        NULL,
        NULL,
        NULL,
        NULL
};

PyObject *
PyInit__aligners(void)

#else

void
init_aligners(void)
#endif

{
  PyObject* module;

  AlignerType.tp_new = PyType_GenericNew;

  if (PyType_Ready(&AlignerType) < 0
   || PyType_Ready(&PathGenerator_Type) < 0)
#if PY_MAJOR_VERSION >= 3
      return NULL;
#else
      return;
#endif

#if PY_MAJOR_VERSION >= 3
    module = PyModule_Create(&moduledef);
#else
    module = Py_InitModule3("_aligners", _aligners_methods, _aligners__doc__);
#endif

  Py_INCREF(&AlignerType);
  PyModule_AddObject(module, "PairwiseAligner", (PyObject*) &AlignerType);

#if PY_MAJOR_VERSION >= 3
  return module;
#endif
}
