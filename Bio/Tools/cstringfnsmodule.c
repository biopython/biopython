/* Copyright 2000 by Jeffrey Chang.  All rights reserved.
 * This code is part of the Biopython distribution and governed by its
 * license.  Please see the LICENSE file that should have been included
 * as part of this package.
 *
 * cstringfnsmodule.c
 * Created 7 Jun 2000
 */

#include "Python.h"
#include <string.h>  /* memset */


/* Functions in this module. */

static char cstringfns_splitany__doc__[] = 
"splitany(str [,sep [,maxsplit [,negate]]]) -> list of strings\n\
\n\
Split a string.  Similar to string.split, except that this considers\n\
any one of the characters in sep to be a delimiter.  If negate is\n\
true, then everything but sep will be a separator.\n\
\n\
";

static PyObject *cstringfns_splitany(
     PyObject *self, PyObject *args, PyObject *keywds)
{
    int i, prev;
    int nsplit, maxsplit=0;
    /*int negate=0;*/
    PyObject *py_negate=NULL;
    PyObject *strlist, *newstr;
    unsigned char *str, 
	*sep=" \011\012\013\014\015";  /* whitespace */
    char tosplit[256];
    static char *kwlist[] = {"str", "sep", "maxsplit", "negate", NULL};

    if(!PyArg_ParseTupleAndKeywords(args, keywds, "s|siO", kwlist, 
				    &str, &sep, &maxsplit, &py_negate))
	return NULL;
    if(maxsplit < 0)
	maxsplit = 1;
    /* negate = (py_negate && PyObject_IsTrue(py_negate));*/
    /* XXX NO MORE NEGATE */

    /* Set the tosplit array to 1 for characters to split on. */
    memset(tosplit, 0, 256);
    while(*sep) {
	tosplit[(unsigned char)*sep++] = 1;
	}
    if(py_negate && PyObject_IsTrue(py_negate)) {
	for(i=0; i<256; i++)
	    tosplit[i] = !tosplit[i];
    }

    /* Create a new list to store the variables. */
    if(!(strlist = PyList_New(0))) {
	PyErr_SetString(PyExc_SystemError, "I could not create a new list");
	return NULL;
    }

    prev = 0;
    nsplit = 0;
    for(i=0; str[i] && (maxsplit == 0 || nsplit < maxsplit); i++) {
	/*if(!(tosplit[(int)str[i]] == !negate))
	  continue; */
	if(!tosplit[(int)str[i]])
	    continue;

	/* Split the string here. */
	if(!(newstr = PyString_FromStringAndSize(&str[prev], i-prev))) {
	    PyErr_SetString(PyExc_SystemError, 
			    "I could not create a new string");
	    break;
	}
	if(PyList_Append(strlist, newstr) == -1) {
	    Py_DECREF(newstr);
	    break;
	}
	Py_DECREF(newstr);
	prev = i+1;
	nsplit++;
    }
    if(!PyErr_Occurred()) {
	i = strlen(str);
	/* Add the last one. */
	if(!(newstr = PyString_FromStringAndSize(&str[prev], i-prev))) {
	    PyErr_SetString(PyExc_SystemError, 
			    "I could not create a new string");
	} else {
	    PyList_Append(strlist, newstr);
	    Py_DECREF(newstr);
	}
    } else {
	Py_DECREF(strlist);
	return NULL;
    }


    return strlist;
}



/* Module definition stuff */

static PyMethodDef cstringfnsMethods[] = {
  {"splitany", (PyCFunction)cstringfns_splitany, METH_VARARGS|METH_KEYWORDS, 
   cstringfns_splitany__doc__},
  {NULL, NULL}
};

static char cstringfns__doc__[] =
"This provides helper functions for the stringfns module.\n\
You should never import this module on its own.\n\
\n\
";

void initcstringfns(void)
{
  (void) Py_InitModule3("cstringfns", cstringfnsMethods, cstringfns__doc__);
}
