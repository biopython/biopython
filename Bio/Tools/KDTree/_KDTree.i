%module _KDTree
%include constraints.i
%include typemaps.i 

%{
	#include "_KDTree.h"
	#include <Numeric/arrayobject.h> 
%}

// import_array() call initialises Numpy
%init 
%{
	import_array();
%}

// Return PyArrayObjects unchanged 
%typemap(python,out) PyArrayObject * 
{
  $target = (PyObject *) $source;
}      

// Return PyObject's unchanged 
%typemap(python,out) PyObject * 
{
  $target = $source;
}      

// Handle input of 0xD coordinates as a Numpy array
%typemap(python, in) float *onedim
{
	float *coord_data;
	PyArrayObject *array;
	long int n, i;

	array=(PyArrayObject *) $source;

	/* Check if it is an array */
	if (PyArray_Check($source))
	{
		if(array->nd!=1)
		{
			PyErr_SetString(PyExc_ValueError, "Array must be one dimensional.");
			return NULL;
		}	

		n=array->dimensions[0];

		// coord_data is deleted by the KDTree object
		coord_data=new float [n];

		for (i=0; i<n; i++)
		{
			coord_data[i]=*(float *) (array->data+i*array->strides[0]);
		}	
		$target=coord_data;
	}
	else
	{
		return NULL;
	}
}

// Handle input of NxD coordinates as a Numpy array
%typemap(python, in) float *twodim
{
	float *coord_data;
	PyArrayObject *array;
	long int n, m, i;

	array=(PyArrayObject *) $source;

	/* Check if it is an array */
	if (PyArray_Check($source))
	{
		if(array->nd!=2)
		{
			PyErr_SetString(PyExc_ValueError, "Array must be two dimensional.");
			return NULL;
		}	

		n=array->dimensions[0];
		m=array->dimensions[1];

		// coord_data is deleted by the KDTree object
		coord_data=new float [m*n];

		for (i=0; i<n; i++)
		{
			int j;

			for (j=0; j<m; j++)
			{
				coord_data[i*m+j]=*(float *) (array->data+i*array->strides[0]+j*array->strides[1]);
			}
		}	
		$target=coord_data;
	}
	else
	{
		return NULL;
	}
}	

// This is actually swigged

class KDTree
{
	public:
		KDTree(int POSITIVE, int POSITIVE);
		~KDTree();
		void set_data(float *twodim, unsigned long int POSITIVE);
		void search_center_radius(float *onedim, float POSITIVE);
		long int get_count(void);
		void neighbor_search(float POSITIVE);
		void neighbor_simple_search(float POSITIVE);
		long int neighbor_get_count(void);
};		


// Add a function to return indices of coordinates within radius
// as a Numpy array
%{
	PyObject *KDTree_get_indices(KDTree *kdtree)
	{
		int length[1];
		PyArrayObject *_array;
	 
		length[0]=kdtree->get_count();
		
		if (length[0]==0)
		{
			Py_INCREF(Py_None);
			return Py_None;
		}

		_array=(PyArrayObject *) PyArray_FromDims(1, length, PyArray_LONG);

		// copy the data into the Numpy data pointer
		kdtree->copy_indices((long int *) _array->data);
		return PyArray_Return(_array);
	}                                   
%}

// Add a function to return indices of coordinates within radius
// as a Numpy array
%{
	PyObject *KDTree_neighbor_get_indices(KDTree *kdtree)
	{
		int length[1];
		PyArrayObject *_array;
	 
		length[0]=2*kdtree->neighbor_get_count();
		
		if (length[0]==0)
		{
			Py_INCREF(Py_None);
			return Py_None;
		}

		_array=(PyArrayObject *) PyArray_FromDims(1, length, PyArray_LONG);

		// copy the data into the Numpy data pointer
		kdtree->neighbor_copy_indices((long int *) _array->data);
		return PyArray_Return(_array);
	}                                   
%}

// Add a function to return distances of coordinates within radius
// as a Numpy array
%{
	PyObject *KDTree_get_radii(KDTree *kdtree)
	{
		int length[1];
		PyArrayObject *_array;
	 
		length[0]=kdtree->get_count();

		if (length[0]==0)
		{
			Py_INCREF(Py_None);
			return Py_None;
		}
		
		_array=(PyArrayObject *) PyArray_FromDims(1, length, PyArray_FLOAT);

		// copy the data into the Numpy data pointer
		kdtree->copy_radii((float *) _array->data);
		return PyArray_Return(_array);
	}                                   
%}

// Add a function to return distances of coordinates within radius
// as a Numpy array
%{
	PyObject *KDTree_neighbor_get_radii(KDTree *kdtree)
	{
		int length[1];
		PyArrayObject *_array;
	 
		length[0]=kdtree->neighbor_get_count();

		if (length[0]==0)
		{
			Py_INCREF(Py_None);
			return Py_None;
		}
		
		_array=(PyArrayObject *) PyArray_FromDims(1, length, PyArray_FLOAT);

		// copy the data into the Numpy data pointer
		kdtree->neighbor_copy_radii((float *) _array->data);
		return PyArray_Return(_array);
	}                                   
%}

// Add above two methods to KDTree class
%addmethods KDTree 
{
	PyObject *get_indices();
	PyObject *get_radii();
	PyObject *neighbor_get_indices();
	PyObject *neighbor_get_radii();
}	

