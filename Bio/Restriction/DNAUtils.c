/* DNA IUPAC ambiguous                "ACBDGHKMNSRTWVYacbdghkmnsrtwvy"
// DNA IUPAC ambiguous complement     "TGVHCDMKNSYAWBRtgvhcdmknsyawbr"

DNA_complement from PyString_translate (stringobject.c Python2.3.4)
*/
#include "Python.h"

static PyObject *
DNA_checkbases(PyObject *bogus, PyStringObject *string)
{
        register char *input, *output;
        register int i;
        register char c;
        register const char  *table =
                 "         zz  z                  z             "
                 "  zzzzzzzzzz       ABCD  GH  K MN   RST VW Y z"
                 "     ABCD  GH  K MN   RST VW Y               "
                 "                                              "
                 "                                              "
                 "                         ";
        PyObject *input_obj = (PyObject *)string, *result;
        int len;
        i = PyString_Size(input_obj);
        input = PyString_AsString(input_obj);
        result = PyString_FromStringAndSize((char *)NULL, i+1);
        if (result == NULL)
                return NULL;

        output = PyString_AsString(result);
        output[0] = ' ';
        len = 1;
        for (;--i >= 0;) {
                c = table[Py_CHARMASK(*input++)];
                if (c == 'z') continue;
                else if (c == ' ') {
                        PyErr_SetString(PyExc_TypeError, "All bases must IUPAC letters");
                                 return NULL;
                         }
                output[len] = c;
                len++;
        }
        _PyString_Resize(&result, len);
        return result;
}

static PyObject *
DNA_complement(PyObject *bogus, PyStringObject *string)
{
        register char *input, *output;
        register int i;
        register const char  *table =
                 "                                              "
                 "                   TVGHEFCDIJMLKNOPQYSAUBWXRZ "
                 "     tvghefcdijmlknopqysaubwxrz              "
                 "                                              "
                 "                                              "
                 "                         ";
        PyObject *input_obj = (PyObject *)string, *result;

        i = PyString_Size(input_obj);
        result = PyString_FromStringAndSize((char *)NULL, i);
        if (result == NULL)
                return NULL;
        output = PyString_AsString(result);
        input = PyString_AsString(input_obj);

        for (;--i >= 0;) *output++ = table[Py_CHARMASK(*input++)];

        return result;
}

static PyObject *
DNA_rev_compl(PyObject *bogus, PyStringObject *string)
{
        register char *input, *output;
        register int i;
        register const char *table =
                "                                              "
                "                   TVGHEFCDIJMLKNOPQYSAUBWXRZ "
                "     tvghefcdijmlknopqysaubwxrz              "
                "                                              "
                "                                              "
                "                         ";
        PyObject *input_obj = (PyObject *)string, *result;

        i = PyString_Size(input_obj);
        result = PyString_FromStringAndSize((char *)NULL, i);
        if (result == NULL)
                return NULL;
        output = PyString_AsString(result);
        input = PyString_AsString(input_obj);

        for (;--i >= 0;) output[i] = table[Py_CHARMASK(*input++)];
        return result;
}


PyDoc_STRVAR(DNA_complement_doc,
"complement(string) -> string.\n\
\n\
the new string is the complementary strand of the \n\
input. string must be a IUPAC ambiguous alphabet.\n");

PyDoc_STRVAR(DNA_rev_compl_doc,
"antiparallel(string) -> string.\n\
\n\
the new string is the reverse complementary strand \n\
(antiparallel) of the input.\n\
 string must be a IUPAC ambiguous alphabet.\n");

 PyDoc_STRVAR(DNA_checkbases_doc,
"check_bases(string) -> string.\n\
\n\
check that the bases in string are in the IUPAC ambiguous alphabet.\n\
Remove digits and white space present in string.\n");

static PyMethodDef DNAUtils_functions[] = {
        {"check_bases",          (PyCFunction)DNA_checkbases, METH_O, DNA_checkbases_doc},
        {"complement",           (PyCFunction)DNA_complement, METH_O, DNA_complement_doc},
        {"antiparallel",         (PyCFunction)DNA_rev_compl, METH_O, DNA_rev_compl_doc},
        {NULL,                NULL}                /* sentinel */
};

PyMODINIT_FUNC
initDNAUtils(void)
{
        PyObject *m;

        m = Py_InitModule3("DNAUtils",
                           DNAUtils_functions,
                           "some functions for use with DNA object.");
        if (m == NULL)
                return;
}
