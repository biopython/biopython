// Copyright 2004 by Harry Zuzan.  All rights reserved.
// This code is part of the Biopython distribution and governed by its
// license.  Please see the LICENSE file that should have been included
// as part of this package.

// This source code parses an Affymetrix cel file.
//
// Documentation to come.
//

#include "Python.h"
#include "Numeric/arrayobject.h"

using namespace std;

#include <string>
#include <deque>


static PyObject * cel_parse(PyObject * self, PyObject * args)
	{
	char * cel_file;
	int len_cel_file;

	if (!PyArg_ParseTuple(args, "s#", &cel_file, &len_cel_file)) return 0;

	deque<string>lines;
	string tmp_str = "";

	for (int i=0; i<len_cel_file; i++)
		{
		char c = cel_file[i];

		if ( c == '\n' )
			{
			lines.push_back(tmp_str);
			tmp_str = "";
			}

		else if ( c != '\r' ) tmp_str = tmp_str + c;
		}

	string HEADER = "[HEADER]";
	while (1)
		{
		string line = lines.front();
		lines.pop_front();
		if (line == HEADER) break;
		}

	size_t token_pos;

	token_pos = lines.front().find("=") + 1;
	int ncols = atoi(string(lines.front(),token_pos).c_str());
	lines.pop_front();

	token_pos = lines.front().find("=") + 1;
	int nrows = atoi(string(lines.front(),token_pos).c_str());
	lines.pop_front();


	string INTENSITY = "[INTENSITY]";
	while (lines.size())
		{
		string line = lines.front();
		lines.pop_front();
		if (line == INTENSITY) break;
		}

	lines.pop_front();
	lines.pop_front();

	int dims[2] = {nrows,ncols};


	PyArrayObject * mean_py = (PyArrayObject *)
		PyArray_FromDims(2, dims, PyArray_DOUBLE);

	PyArrayObject * stdev_py = (PyArrayObject *)
		PyArray_FromDims(2, dims, PyArray_DOUBLE);

	PyArrayObject * npix_py = (PyArrayObject *)
		PyArray_FromDims(2, dims, PyArray_LONG);

	for (int i=0; i<nrows*ncols; i++)
		{
		string col_str = "";
		string row_str = "";
		string mean_str = "";
		string stdev_str = "";
		string npix_str = "";

		string data_str = lines.front();
		lines.pop_front();

		unsigned int k = 0;
		while ( data_str[k] == ' ' || data_str[k] == '\t') k++;
		while ( data_str[k] != ' ' && data_str[k] != '\t')
			col_str += data_str[k++];

		while ( data_str[k] == ' ' || data_str[k] == '\t') k++;
		while ( data_str[k] != ' ' && data_str[k] != '\t')
			row_str += data_str[k++];

		while ( data_str[k] == ' ' || data_str[k] == '\t') k++;
		while ( data_str[k] != ' ' && data_str[k] != '\t')
			mean_str += data_str[k++]; k++;

		while ( data_str[k] == ' ' || data_str[k] == '\t') k++;
		while ( data_str[k] != ' ' && data_str[k] != '\t')
			stdev_str += data_str[k++];

		while ( data_str[k] == ' ' || data_str[k] == '\t') k++;
		while ( k < data_str.size()) npix_str += data_str[k++];

		int row = atoi(row_str.c_str());
		int col = atoi(col_str.c_str());

		double mean = atof(mean_str.c_str());
		double stdev = atof(stdev_str.c_str());

		int npix = atoi(npix_str.c_str());


		((double *)mean_py->data)[ncols*row + col] = mean;
		((double *)stdev_py->data)[ncols*row + col] = stdev;
		((long *)npix_py->data)[ncols*row + col] = npix;
		}

	PyObject * returnList = PyList_New(0);

	PyList_Append( returnList, (PyObject *)mean_py);
	PyList_Append( returnList, (PyObject *)stdev_py);
	PyList_Append( returnList, (PyObject *)npix_py);

	return returnList;
	}


extern "C" {

static PyMethodDef celMethods[] = {
	{"parse", cel_parse, METH_VARARGS },
	{NULL, NULL}
	};

void init_cel()
{    
Py_InitModule("_cel", celMethods);
import_array();  // need to when using Numeric Python
}

};
