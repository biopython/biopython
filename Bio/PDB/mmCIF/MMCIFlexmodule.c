#include <Python.h>
#include "lex.yy.h"

FILE *fp;

PyObject *MMCIFlex__init__(PyObject *self, PyObject *args) {
    char *filename;
    /*filename = "x";*/
    char *classname;

    if (!PyArg_ParseTuple(args, "s s", &classname &filename))
        return NULL;

    /*printf("%s", classname);*/
    /*printf("%s", filename);*/

    /*fp = fopen(filename, "r");	*/

    /*mmcif_set_file(fp);*/

	Py_INCREF(Py_None);
	return Py_None;
}

/*PyObject *MMCIFlex__del__(PyObject *self, PyObject *args) {*/
	/*[> verify no arguments <]*/
    /*[>if (!PyArg_ParseTuple(args, ""))<]*/
        /*[>return NULL;<]*/

	/*fclose(fp);*/

	/*Py_INCREF(Py_None);*/
	/*return Py_None;*/
/*}	*/

/*PyObject *MMCIFlex_get_token(PyObject *self, PyObject *args) {*/
	/*int flag;*/
	/*char *value="";*/

	/*[> get token number <]*/
	/*flag=mmcif_get_token();*/

	/*[> if flag==0 we are EOF <]*/
	/*if(flag) {*/
		/*value=mmcif_get_string();*/
	/*}	*/

	/*[> return the (tokennumber, string) tuple <]*/
	/*return Py_BuildValue("(is)", flag, value);*/
/*}*/

PyMethodDef MMCIFlexMethods[] =
{
	{"__init__",	MMCIFlex__init__, 	METH_VARARGS},
	/*{"__del__",	MMCIFlex__del__,	METH_VARARGS},*/
	/*{"get_token",  	MMCIFlex_get_token, 	METH_VARARGS},*/
	{NULL,      	NULL}        			/* Sentinel */
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
