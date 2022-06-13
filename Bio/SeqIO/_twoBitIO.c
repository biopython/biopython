#define PY_SSIZE_T_CLEAN
#include "Python.h"


static const char bases[][4] = {"TTTT",  /* 00 00 00 00 */
                                "TTTC",  /* 00 00 00 01 */
                                "TTTA",  /* 00 00 00 10 */
                                "TTTG",  /* 00 00 00 11 */
                                "TTCT",  /* 00 00 01 00 */
                                "TTCC",  /* 00 00 01 01 */
                                "TTCA",  /* 00 00 01 10 */
                                "TTCG",  /* 00 00 01 11 */
                                "TTAT",  /* 00 00 10 00 */
                                "TTAC",  /* 00 00 10 01 */
                                "TTAA",  /* 00 00 10 10 */
                                "TTAG",  /* 00 00 10 11 */
                                "TTGT",  /* 00 00 11 00 */
                                "TTGC",  /* 00 00 11 01 */
                                "TTGA",  /* 00 00 11 10 */
                                "TTGG",  /* 00 00 11 11 */
                                "TCTT",  /* 00 01 00 00 */
                                "TCTC",  /* 00 01 00 01 */
                                "TCTA",  /* 00 01 00 10 */
                                "TCTG",  /* 00 01 00 11 */
                                "TCCT",  /* 00 01 01 00 */
                                "TCCC",  /* 00 01 01 01 */
                                "TCCA",  /* 00 01 01 10 */
                                "TCCG",  /* 00 01 01 11 */
                                "TCAT",  /* 00 01 10 00 */
                                "TCAC",  /* 00 01 10 01 */
                                "TCAA",  /* 00 01 10 10 */
                                "TCAG",  /* 00 01 10 11 */
                                "TCGT",  /* 00 01 11 00 */
                                "TCGC",  /* 00 01 11 01 */
                                "TCGA",  /* 00 01 11 10 */
                                "TCGG",  /* 00 01 11 11 */
                                "TATT",  /* 00 10 00 00 */
                                "TATC",  /* 00 10 00 01 */
                                "TATA",  /* 00 10 00 10 */
                                "TATG",  /* 00 10 00 11 */
                                "TACT",  /* 00 10 01 00 */
                                "TACC",  /* 00 10 01 01 */
                                "TACA",  /* 00 10 01 10 */
                                "TACG",  /* 00 10 01 11 */
                                "TAAT",  /* 00 10 10 00 */
                                "TAAC",  /* 00 10 10 01 */
                                "TAAA",  /* 00 10 10 10 */
                                "TAAG",  /* 00 10 10 11 */
                                "TAGT",  /* 00 10 11 00 */
                                "TAGC",  /* 00 10 11 01 */
                                "TAGA",  /* 00 10 11 10 */
                                "TAGG",  /* 00 10 11 11 */
                                "TGTT",  /* 00 11 00 00 */
                                "TGTC",  /* 00 11 00 01 */
                                "TGTA",  /* 00 11 00 10 */
                                "TGTG",  /* 00 11 00 11 */
                                "TGCT",  /* 00 11 01 00 */
                                "TGCC",  /* 00 11 01 01 */
                                "TGCA",  /* 00 11 01 10 */
                                "TGCG",  /* 00 11 01 11 */
                                "TGAT",  /* 00 11 10 00 */
                                "TGAC",  /* 00 11 10 01 */
                                "TGAA",  /* 00 11 10 10 */
                                "TGAG",  /* 00 11 10 11 */
                                "TGGT",  /* 00 11 11 00 */
                                "TGGC",  /* 00 11 11 01 */
                                "TGGA",  /* 00 11 11 10 */
                                "TGGG",  /* 00 11 11 11 */
                                "CTTT",  /* 01 00 00 00 */
                                "CTTC",  /* 01 00 00 01 */
                                "CTTA",  /* 01 00 00 10 */
                                "CTTG",  /* 01 00 00 11 */
                                "CTCT",  /* 01 00 01 00 */
                                "CTCC",  /* 01 00 01 01 */
                                "CTCA",  /* 01 00 01 10 */
                                "CTCG",  /* 01 00 01 11 */
                                "CTAT",  /* 01 00 10 00 */
                                "CTAC",  /* 01 00 10 01 */
                                "CTAA",  /* 01 00 10 10 */
                                "CTAG",  /* 01 00 10 11 */
                                "CTGT",  /* 01 00 11 00 */
                                "CTGC",  /* 01 00 11 01 */
                                "CTGA",  /* 01 00 11 10 */
                                "CTGG",  /* 01 00 11 11 */
                                "CCTT",  /* 01 01 00 00 */
                                "CCTC",  /* 01 01 00 01 */
                                "CCTA",  /* 01 01 00 10 */
                                "CCTG",  /* 01 01 00 11 */
                                "CCCT",  /* 01 01 01 00 */
                                "CCCC",  /* 01 01 01 01 */
                                "CCCA",  /* 01 01 01 10 */
                                "CCCG",  /* 01 01 01 11 */
                                "CCAT",  /* 01 01 10 00 */
                                "CCAC",  /* 01 01 10 01 */
                                "CCAA",  /* 01 01 10 10 */
                                "CCAG",  /* 01 01 10 11 */
                                "CCGT",  /* 01 01 11 00 */
                                "CCGC",  /* 01 01 11 01 */
                                "CCGA",  /* 01 01 11 10 */
                                "CCGG",  /* 01 01 11 11 */
                                "CATT",  /* 01 10 00 00 */
                                "CATC",  /* 01 10 00 01 */
                                "CATA",  /* 01 10 00 10 */
                                "CATG",  /* 01 10 00 11 */
                                "CACT",  /* 01 10 01 00 */
                                "CACC",  /* 01 10 01 01 */
                                "CACA",  /* 01 10 01 10 */
                                "CACG",  /* 01 10 01 11 */
                                "CAAT",  /* 01 10 10 00 */
                                "CAAC",  /* 01 10 10 01 */
                                "CAAA",  /* 01 10 10 10 */
                                "CAAG",  /* 01 10 10 11 */
                                "CAGT",  /* 01 10 11 00 */
                                "CAGC",  /* 01 10 11 01 */
                                "CAGA",  /* 01 10 11 10 */
                                "CAGG",  /* 01 10 11 11 */
                                "CGTT",  /* 01 11 00 00 */
                                "CGTC",  /* 01 11 00 01 */
                                "CGTA",  /* 01 11 00 10 */
                                "CGTG",  /* 01 11 00 11 */
                                "CGCT",  /* 01 11 01 00 */
                                "CGCC",  /* 01 11 01 01 */
                                "CGCA",  /* 01 11 01 10 */
                                "CGCG",  /* 01 11 01 11 */
                                "CGAT",  /* 01 11 10 00 */
                                "CGAC",  /* 01 11 10 01 */
                                "CGAA",  /* 01 11 10 10 */
                                "CGAG",  /* 01 11 10 11 */
                                "CGGT",  /* 01 11 11 00 */
                                "CGGC",  /* 01 11 11 01 */
                                "CGGA",  /* 01 11 11 10 */
                                "CGGG",  /* 01 11 11 11 */
                                "ATTT",  /* 10 00 00 00 */
                                "ATTC",  /* 10 00 00 01 */
                                "ATTA",  /* 10 00 00 10 */
                                "ATTG",  /* 10 00 00 11 */
                                "ATCT",  /* 10 00 01 00 */
                                "ATCC",  /* 10 00 01 01 */
                                "ATCA",  /* 10 00 01 10 */
                                "ATCG",  /* 10 00 01 11 */
                                "ATAT",  /* 10 00 10 00 */
                                "ATAC",  /* 10 00 10 01 */
                                "ATAA",  /* 10 00 10 10 */
                                "ATAG",  /* 10 00 10 11 */
                                "ATGT",  /* 10 00 11 00 */
                                "ATGC",  /* 10 00 11 01 */
                                "ATGA",  /* 10 00 11 10 */
                                "ATGG",  /* 10 00 11 11 */
                                "ACTT",  /* 10 01 00 00 */
                                "ACTC",  /* 10 01 00 01 */
                                "ACTA",  /* 10 01 00 10 */
                                "ACTG",  /* 10 01 00 11 */
                                "ACCT",  /* 10 01 01 00 */
                                "ACCC",  /* 10 01 01 01 */
                                "ACCA",  /* 10 01 01 10 */
                                "ACCG",  /* 10 01 01 11 */
                                "ACAT",  /* 10 01 10 00 */
                                "ACAC",  /* 10 01 10 01 */
                                "ACAA",  /* 10 01 10 10 */
                                "ACAG",  /* 10 01 10 11 */
                                "ACGT",  /* 10 01 11 00 */
                                "ACGC",  /* 10 01 11 01 */
                                "ACGA",  /* 10 01 11 10 */
                                "ACGG",  /* 10 01 11 11 */
                                "AATT",  /* 10 10 00 00 */
                                "AATC",  /* 10 10 00 01 */
                                "AATA",  /* 10 10 00 10 */
                                "AATG",  /* 10 10 00 11 */
                                "AACT",  /* 10 10 01 00 */
                                "AACC",  /* 10 10 01 01 */
                                "AACA",  /* 10 10 01 10 */
                                "AACG",  /* 10 10 01 11 */
                                "AAAT",  /* 10 10 10 00 */
                                "AAAC",  /* 10 10 10 01 */
                                "AAAA",  /* 10 10 10 10 */
                                "AAAG",  /* 10 10 10 11 */
                                "AAGT",  /* 10 10 11 00 */
                                "AAGC",  /* 10 10 11 01 */
                                "AAGA",  /* 10 10 11 10 */
                                "AAGG",  /* 10 10 11 11 */
                                "AGTT",  /* 10 11 00 00 */
                                "AGTC",  /* 10 11 00 01 */
                                "AGTA",  /* 10 11 00 10 */
                                "AGTG",  /* 10 11 00 11 */
                                "AGCT",  /* 10 11 01 00 */
                                "AGCC",  /* 10 11 01 01 */
                                "AGCA",  /* 10 11 01 10 */
                                "AGCG",  /* 10 11 01 11 */
                                "AGAT",  /* 10 11 10 00 */
                                "AGAC",  /* 10 11 10 01 */
                                "AGAA",  /* 10 11 10 10 */
                                "AGAG",  /* 10 11 10 11 */
                                "AGGT",  /* 10 11 11 00 */
                                "AGGC",  /* 10 11 11 01 */
                                "AGGA",  /* 10 11 11 10 */
                                "AGGG",  /* 10 11 11 11 */
                                "GTTT",  /* 11 00 00 00 */
                                "GTTC",  /* 11 00 00 01 */
                                "GTTA",  /* 11 00 00 10 */
                                "GTTG",  /* 11 00 00 11 */
                                "GTCT",  /* 11 00 01 00 */
                                "GTCC",  /* 11 00 01 01 */
                                "GTCA",  /* 11 00 01 10 */
                                "GTCG",  /* 11 00 01 11 */
                                "GTAT",  /* 11 00 10 00 */
                                "GTAC",  /* 11 00 10 01 */
                                "GTAA",  /* 11 00 10 10 */
                                "GTAG",  /* 11 00 10 11 */
                                "GTGT",  /* 11 00 11 00 */
                                "GTGC",  /* 11 00 11 01 */
                                "GTGA",  /* 11 00 11 10 */
                                "GTGG",  /* 11 00 11 11 */
                                "GCTT",  /* 11 01 00 00 */
                                "GCTC",  /* 11 01 00 01 */
                                "GCTA",  /* 11 01 00 10 */
                                "GCTG",  /* 11 01 00 11 */
                                "GCCT",  /* 11 01 01 00 */
                                "GCCC",  /* 11 01 01 01 */
                                "GCCA",  /* 11 01 01 10 */
                                "GCCG",  /* 11 01 01 11 */
                                "GCAT",  /* 11 01 10 00 */
                                "GCAC",  /* 11 01 10 01 */
                                "GCAA",  /* 11 01 10 10 */
                                "GCAG",  /* 11 01 10 11 */
                                "GCGT",  /* 11 01 11 00 */
                                "GCGC",  /* 11 01 11 01 */
                                "GCGA",  /* 11 01 11 10 */
                                "GCGG",  /* 11 01 11 11 */
                                "GATT",  /* 11 10 00 00 */
                                "GATC",  /* 11 10 00 01 */
                                "GATA",  /* 11 10 00 10 */
                                "GATG",  /* 11 10 00 11 */
                                "GACT",  /* 11 10 01 00 */
                                "GACC",  /* 11 10 01 01 */
                                "GACA",  /* 11 10 01 10 */
                                "GACG",  /* 11 10 01 11 */
                                "GAAT",  /* 11 10 10 00 */
                                "GAAC",  /* 11 10 10 01 */
                                "GAAA",  /* 11 10 10 10 */
                                "GAAG",  /* 11 10 10 11 */
                                "GAGT",  /* 11 10 11 00 */
                                "GAGC",  /* 11 10 11 01 */
                                "GAGA",  /* 11 10 11 10 */
                                "GAGG",  /* 11 10 11 11 */
                                "GGTT",  /* 11 11 00 00 */
                                "GGTC",  /* 11 11 00 01 */
                                "GGTA",  /* 11 11 00 10 */
                                "GGTG",  /* 11 11 00 11 */
                                "GGCT",  /* 11 11 01 00 */
                                "GGCC",  /* 11 11 01 01 */
                                "GGCA",  /* 11 11 01 10 */
                                "GGCG",  /* 11 11 01 11 */
                                "GGAT",  /* 11 11 10 00 */
                                "GGAC",  /* 11 11 10 01 */
                                "GGAA",  /* 11 11 10 10 */
                                "GGAG",  /* 11 11 10 11 */
                                "GGGT",  /* 11 11 11 00 */
                                "GGGC",  /* 11 11 11 01 */
                                "GGGA",  /* 11 11 11 10 */
                                "GGGG",  /* 11 11 11 11 */
                               };

static int
extract(const unsigned char* bytes, uint32_t byteSize, uint32_t start, uint32_t end, char sequence[]) {
    uint32_t i;
    const uint32_t size = end - start;
    const uint32_t byteStart = start / 4;
    const uint32_t byteEnd = (end + 3) / 4;

    if (byteSize != byteEnd - byteStart) {
        PyErr_Format(PyExc_RuntimeError,
                     "unexpected number of bytes %u (expected %u)",
                     byteSize, byteEnd - byteStart);
        return -1;
    }

    start -= byteStart * 4;
    if (byteStart + 1 == byteEnd) {
        /* one byte only */
        memcpy(sequence, &(bases[*bytes][start]), size);
    }
    else {
        end -= byteEnd * 4;
        /* end is now a negative number equal to the distance to the byte end */
        memcpy(sequence, &(bases[*bytes][start]), 4 - start);
        bytes++;
        sequence += (4 - start);
        for (i = byteStart+1; i < byteEnd-1; i++, bytes++, sequence += 4)
            memcpy(sequence, bases[*bytes], 4);
        memcpy(sequence, bases[*bytes], end + 4);
        bytes++;
        bytes -= byteSize;
    }
    return 0;
}

static void
applyNs(char sequence[], uint32_t start, uint32_t end, Py_buffer *nBlocks)
{
    const Py_ssize_t nBlockCount = nBlocks->shape[0];
    const uint32_t* const nBlockPositions = nBlocks->buf;

    Py_ssize_t i;
    for (i = 0; i < nBlockCount; i++) {
        uint32_t nBlockStart = nBlockPositions[2*i];
        uint32_t nBlockEnd = nBlockPositions[2*i+1];
        if (nBlockEnd < start) continue;
        if (end < nBlockStart) break;
        if (nBlockStart < start) nBlockStart = start;
        if (end < nBlockEnd) nBlockEnd = end;
        memset(sequence + nBlockStart - start, 'N', nBlockEnd - nBlockStart);
    }
}

static void
applyMask(char sequence[], uint32_t start, uint32_t end, Py_buffer* maskBlocks)
{
    const Py_ssize_t maskBlockCount = maskBlocks->shape[0];
    const uint32_t* const maskBlockPositions = maskBlocks->buf;
    const char diff = 'a' - 'A';

    Py_ssize_t i;
    for (i = 0; i < maskBlockCount; i++) {
        uint32_t j;
        uint32_t maskBlockStart = maskBlockPositions[2*i];
        uint32_t maskBlockEnd = maskBlockPositions[2*i+1];
        if (maskBlockEnd < start) continue;
        if (end < maskBlockStart) break;
        if (maskBlockStart < start) maskBlockStart = start;
        if (end < maskBlockEnd) maskBlockEnd = end;
        for (j = maskBlockStart - start; j < maskBlockEnd - start; j++)
            sequence[j] += diff;
    }
}

static int
blocks_converter(PyObject* object, void* pointer)
{
    const int flag = PyBUF_ND | PyBUF_FORMAT;
    Py_buffer *view = pointer;

    if (object == NULL) goto exit;

    if (PyObject_GetBuffer(object, view, flag) == -1) {
        PyErr_SetString(PyExc_RuntimeError, "blocks have unexpected format.");
        return 0;
    }

    if (view->itemsize != sizeof(uint32_t)
     || (strcmp(view->format, "I") != 0 && strcmp(view->format, "L") != 0 )) {
        PyErr_Format(PyExc_RuntimeError,
                     "blocks have incorrect data type (itemsize %zd, format %s)",
                     view->itemsize, view->format);
        goto exit;
    }
    if (view->ndim != 2) {
        PyErr_Format(PyExc_RuntimeError,
                     "blocks have incorrect rank %d (expected 2)", view->ndim);
        goto exit;
    }
    if (view->shape[1] != 2) {
        PyErr_Format(PyExc_RuntimeError,
                     "blocks should have two columns (found %zd)",
                     view->shape[1]);
        goto exit;
    }
    return Py_CLEANUP_SUPPORTED;

exit:
    PyBuffer_Release(view);
    return 0;
}

static char TwoBit_convert__doc__[] = "convert twoBit data to the DNA sequence, apply blocks of N's (representing unknown sequences) and masked (lower case) blocks, and return the sequence as a bytes object";

static PyObject*
TwoBit_convert(PyObject* self, PyObject* args, PyObject* keywords)
{
    const unsigned char *data;
    Py_ssize_t start;
    Py_ssize_t end;
    Py_ssize_t step;
    Py_ssize_t size;
    Py_ssize_t length;
    Py_buffer nBlocks;
    Py_buffer maskBlocks;
    PyObject *object;
    char *sequence;

    static char* kwlist[] = {"data", "start", "end", "step",
                             "nBlocks", "maskBlocks", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywords, "y#nnnO&O&", kwlist,
                                     &data, &length, &start, &end, &step,
                                     &blocks_converter, &nBlocks,
                                     &blocks_converter, &maskBlocks))
        return NULL;

    size = (end - start) / step;
    object = PyBytes_FromStringAndSize(NULL, size);
    if (!object) goto exit;

    sequence = PyBytes_AS_STRING(object);

    if (step == 1) {
        if (extract(data, length, start, end, sequence) < 0) {
            Py_DECREF(object);
            object = NULL;
            goto exit;
        }
        applyNs(sequence, start, end, &nBlocks);
        applyMask(sequence, start, end, &maskBlocks);
    }
    else {
        Py_ssize_t current, i;
        Py_ssize_t full_start, full_end;
        char* full_sequence;
        if (start <= end) {
            full_start = start;
            full_end = end;
            current = 0; /* first position in sequence */
        }
        else {
            full_start = end + 1;
            full_end = start + 1;
            current = start - end - 1; /* last position in sequence */
        }
        full_sequence = PyMem_Malloc((full_end-full_start+1)*sizeof(char));
        full_sequence[full_end-full_start] = '\0';
        if (!full_sequence) {
            Py_DECREF(object);
            object = NULL;
            goto exit;
        }
        if (extract(data, length, full_start, full_end, full_sequence) < 0) {
            PyMem_Free(full_sequence);
            Py_DECREF(object);
            object = NULL;
            goto exit;
        }
        applyNs(full_sequence, full_start, full_end, &nBlocks);
        applyMask(full_sequence, full_start, full_end, &maskBlocks);
        for (i = 0; i < size; current += step, i++)
            sequence[i] = full_sequence[current];
        PyMem_Free(full_sequence);
    }

exit:
    blocks_converter(NULL, &nBlocks);
    blocks_converter(NULL, &maskBlocks);
    return object;
}

static struct PyMethodDef _twoBitIO_methods[] = {
    {"convert",
     (PyCFunction)TwoBit_convert,
     METH_VARARGS | METH_KEYWORDS,
     TwoBit_convert__doc__
    },
    {NULL,          NULL, 0, NULL} /* sentinel */
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_twoBitIO",
    "Parser for DNA sequence data in 2bit format",
    -1,
    _twoBitIO_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *
PyInit__twoBitIO(void)
{
    return PyModule_Create(&moduledef);
}
