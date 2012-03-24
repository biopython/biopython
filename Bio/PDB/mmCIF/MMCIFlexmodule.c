#include <Python.h>
#include "lex.yy.h"

FILE *fp;
char *classinst;

PyObject *MMCIFlex__init__(PyObject *self, PyObject *args) {
    char *filename;

    /* First member of args is a PyObject of the class instance */
    if (!PyArg_ParseTuple(args, "Os", &classinst, &filename))
        return NULL;

    /*printf("Opening file '%s'\n", filename);*/

    fp = fopen(filename, "r");  

    mmcif_set_file(fp);

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *MMCIFlex__del__(PyObject *self, PyObject *args) {
    /* verify only argument is PyObject of class instance */
    if (!PyArg_ParseTuple(args, "O", &classinst))
        return NULL;

    fclose(fp);

    Py_INCREF(Py_None);
    return Py_None;
}   

PyObject *MMCIFlex_get_token(PyObject *self, PyObject *args) {
    int flag = 0;
    char *value="";

    /* get token number */
    flag=mmcif_get_token();

    /* if flag==0 we are EOF */
    if(flag) {
        value=mmcif_get_string();
    }   

    /* return the (tokennumber, string) tuple */
    return Py_BuildValue("(is)", flag, value);
}

PyMethodDef MMCIFlexMethods[] =
{
    {"__init__", MMCIFlex__init__, METH_VARARGS, "Open file"},
    {"__del__", MMCIFlex__del__, METH_VARARGS, "Close file"},
    {"get_token", MMCIFlex_get_token, METH_VARARGS, "Emit tok (type,val)"},
    {NULL, NULL}  /* Sentinel */
};

PyMethodDef ModuleMethods[] = { {NULL} };

void initMMCIFlex(void) {
    PyMethodDef *def;

    /* create a new module and class */
    PyObject *module = Py_InitModule("MMCIFlex", ModuleMethods);
    PyObject *moduleDict = PyModule_GetDict(module);
    PyObject *classDict = PyDict_New();
    PyObject *className = PyString_FromString("MMCIFlex");
    PyObject *mmciflexClass = PyClass_New(NULL, classDict, className);
    PyDict_SetItemString(moduleDict, "MMCIFlex", mmciflexClass);
    Py_DECREF(classDict);
    Py_DECREF(className);
    Py_DECREF(mmciflexClass);
    
    /* add methods to class */
    for (def = MMCIFlexMethods; def->ml_name != NULL; def++) {
        PyObject *func = PyCFunction_New(def, NULL);
        PyObject *method = PyMethod_New(func, NULL, mmciflexClass);
        PyDict_SetItemString(classDict, def->ml_name, method);
        Py_DECREF(func);
        Py_DECREF(method);
    }
}
