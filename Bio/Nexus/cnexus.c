/* Copyright 2005 by Frank Kauff.  All rights reserved.
 * This code is part of the Biopython distribution and governed by its
 * license. Please see the LICENSE file that should have been included
 * as part of this package.
 *
 * cnexus.c
 *
 * Parse input strings, cut out (nested) comments, deal with quoted text.
 * Input lines terminated with ; are separated by ASCII code 7 (something
 * that naturally doesn't occur in plain NEXUS files).
 *
 * Used by Nexus.py
 */

#include <Python.h>
#include <string.h>

static PyObject * cnexus_scanfile(PyObject *self, PyObject *args)
{
    PyObject *cleaninput;
    const char *input;
    char *scanned, *scanned_start;
    char t, quotelevel;
    int speciallevel, commlevel;

    quotelevel=0;
    speciallevel=0;
    commlevel=0;
    
    if (!PyArg_ParseTuple(args, "s", &input))
        return NULL;
    if (!(scanned=malloc(strlen(input)+1)))
        PyErr_NoMemory();
    scanned_start=scanned;
    for(t=*input;(t=*input);input++) 
    {
        /* end of standard quote */
        if (!(commlevel || speciallevel) && t==quotelevel)
            quotelevel=0;
        /* start of standard quote */
        else if (!quotelevel && !(commlevel || speciallevel) && (t=='\'' || t=='"'))
            quotelevel=t;
        /* start of comment outside quote */
        else if (!quotelevel  && t=='[')
        {
            /* check for special comments */
            /*if ((*(input+1)=='!' || *(input+1)=='&' || *(input+1)=='%' || 
                    *(input+1)=='/' || *(input+1)=='\\' || *(input+1)=='@')
                    && !(commlevel || speciallevel))
                speciallevel=1;
            */
            if ((*(input+1)=='&') && !(commlevel || speciallevel))
                speciallevel=1;
            else /* standard comment */
                commlevel++;
        }
        else if (!quotelevel && t==']')
        {
            /* does it end a special comment? */
            if (speciallevel)
                speciallevel=0;
            else
            {
                commlevel--;
                if (commlevel<0) /* error: unmatched ] */
                {
                    free(scanned_start);
                    return Py_BuildValue("s","]");
                }
                continue;
            }
        }
        if (!commlevel)
        {
            /* we replace the ; at the end of command lines with special
             * character to make subsequent parsing of blocks easier */
            if (t==';' && !(quotelevel || speciallevel))
                /* need an ascii code thats not part of a nexus file, used as
                 * separator */
                *(scanned++)=7;
            else
                *(scanned++)=t;
        }
        /* printf("t %c, commlevel %d, speciallevel %d, quotelevel '%c', scanned %d\n",
         * t,commlevel,speciallevel,quotelevel,scanned);
         */
    }               
    
    if (commlevel>0)
    {
        /* error: unmatched [ */
        free(scanned_start);
        return Py_BuildValue("s","[");
    }
    else
    {
        *scanned=0; /* end of string */
        cleaninput= Py_BuildValue("s",scanned_start);
        free(scanned_start);
        return cleaninput;
    }
}

static PyMethodDef cNexusMethods[]=
{
    {"scanfile",cnexus_scanfile,METH_VARARGS,"Scan file and deal with comments and quotes."},
    {NULL,NULL,0,NULL}
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "cnexus",
        NULL,
        -1,
        cNexusMethods,
        NULL,
        NULL,
        NULL,
        NULL
};

PyObject *
PyInit_cnexus(void)
{
    return PyModule_Create(&moduledef);
}


#else

PyMODINIT_FUNC initcnexus(void)
{
    (void) Py_InitModule("cnexus",cNexusMethods);
}
#endif
