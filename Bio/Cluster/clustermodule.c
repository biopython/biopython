#include "Python.h"
#include "Numeric/arrayobject.h"
#include <stdio.h>
#include <string.h>
#include "cluster.h"
 
static PyObject *ErrorObject;
static char buffer[512];
static char* message = NULL;

static const char known_distances[] = "ebhcauxsk";

/* ========================================================================== */
/* -- Helper routines ------------------------------------------------------- */
/* ========================================================================== */

/* -- data ------------------------------------------------------------------ */

static double**
parse_data(PyObject* object, PyArrayObject** array)
/* Takes the Python object from the argument list, and finds the microarray
 * data set. In case of an error, the array is DECREF'ed and set to NULL. */
{ int i, j;
  int nrows, ncols;
  double** data = NULL;
  if(!PyArray_Check (object)) /* Try to convert object to a 2D double array */
  { *array = (PyArrayObject*) PyArray_FromObject(object, PyArray_DOUBLE, 2, 2);
    if (*array==NULL)
    { strcpy (message, "data cannot be converted to needed array.");
      PyErr_SetString(ErrorObject, buffer);
      return NULL;
    }
  }
  else /* User passed an array */
  { *array = (PyArrayObject*) object;
    Py_INCREF(object);
    if ((*array)->descr->type_num != PyArray_DOUBLE) /* Cast to type double */
    { PyArrayObject* av = (PyArrayObject*) PyArray_Cast(*array, PyArray_DOUBLE);
      Py_DECREF((PyObject*) (*array));
      *array = av;
      if (!(*array))
      { strcpy (message, "data cannot be cast to needed type.");
        PyErr_SetString(ErrorObject, buffer);
        return NULL;
      }
    } 
    if ((*array)->nd != 2) /* Checking number of dimensions */
    { sprintf(message, "data has incorrect rank (%d expected 2)", (*array)->nd);
      PyErr_SetString (ErrorObject, buffer);
      Py_DECREF((PyObject*) (*array));
      *array = NULL;
      return NULL;
    }
  }
  nrows = (*array)->dimensions[0];
  ncols = (*array)->dimensions[1];
  if (nrows < 1 || ncols < 1)
  { strcpy (message, "data is an empty matrix");
    PyErr_SetString (ErrorObject, buffer);
    Py_DECREF((PyObject*) (*array));
    *array = NULL;
    return NULL;
  }
  data = malloc(nrows*sizeof(double*));
  if (((*array)->strides)[1]==sizeof(double)) /* Each row is contiguous */
  { const char* p = (char*) ((*array)->data);
    const int stride =  ((*array)->strides)[0];
    for (i=0; i < nrows; i++, p+=stride) data[i] = (double*)p;
  }
  else /* We need to create contiguous rows */
  { const char* p0 = (char*) (*array)->data;
    const int rowstride =  (*array)->strides[0];
    const int colstride =  (*array)->strides[1];
    for (i=0; i < nrows; i++)
    { const char* p = p0;
      data[i] = malloc(ncols*sizeof(double));
      for (j=0; j < ncols; j++, p+=colstride) data[i][j] = *((double*)p);
      p0 += rowstride;
    }
  }
  return data;
}

static void
free_data(PyArrayObject* array, double** data)
{ if(data[0]!=(double*)(array->data))
  { int i;
    const int nrows = array->dimensions[0];
    for (i=0; i<nrows; i++) free(data[i]);
  }
  free (data);
  Py_DECREF((PyObject*) array);
  return;
}

/* -- mask ------------------------------------------------------------------ */

static int**
parse_mask (PyObject* object, PyArrayObject** array, const int dimensions[2])
{ int i, j;
  const int nrows = dimensions[0];
  const int ncolumns = dimensions[1];
  int** mask;
  if (object==NULL) /* Return the default mask */
  { mask = malloc(nrows*sizeof(int*));
    for (i=0; i<nrows; i++)
    { mask[i] = malloc(ncolumns*sizeof(int));
      for (j=0; j<ncolumns; j++) mask[i][j] = 1;
    }
    *array = NULL;
    return mask;
  }
  if(!PyArray_Check (object)) /* Try to convert object to a 2D double array */
  { *array = (PyArrayObject*) PyArray_FromObject(object, PyArray_LONG, 2, 2);
    if (!(*array))
    { strcpy (message, "mask cannot be converted to needed array");
      PyErr_SetString(ErrorObject, buffer);
      return NULL;
    }
  }
  else /* User passed an array */
  { *array = (PyArrayObject*) object;
    Py_INCREF(object);
    if ((*array)->descr->type_num != PyArray_LONG)
    { PyArrayObject* av = (PyArrayObject*) PyArray_Cast (*array, PyArray_LONG);
      Py_DECREF((PyObject*) *array);
      *array = av;
      if (!(*array))
      { strcpy (message, "mask cannot be cast to needed type.");
        PyErr_SetString(ErrorObject, buffer);
        return NULL;
      }
    } 
  }
  if((*array)->nd != 2) /* Checking number of dimensions */
  { sprintf(message, "mask has incorrect rank (%d expected 2)", (*array)->nd);
    PyErr_SetString (ErrorObject, buffer);
    Py_DECREF((PyObject*)*array);
    *array = NULL;
    return NULL;
  }
  if((*array)->dimensions[0] != nrows) /* Checking number of rows */
  { sprintf(message,
      "mask has incorrect number of rows (%d expected %d)",
      (*array)->dimensions[0], nrows);
    PyErr_SetString (ErrorObject, buffer);
    Py_DECREF((PyObject*)*array);
    *array = NULL;
    return NULL;
  }
  /* no checking on last dimension of expected size 1 */
  if (ncolumns != 1 && (*array)->dimensions[1] != ncolumns)
  { sprintf(message,
      "mask incorrect number of columns (%d expected %d)",
      (*array)->dimensions[1], ncolumns);
    PyErr_SetString (ErrorObject, buffer);
    *array = NULL;
    return NULL;
  }
  /* All checks OK */
  mask = malloc(nrows*sizeof(int*));
  if ((*array)->strides[1]==sizeof(int)) /* Each row is contiguous */
  { const char* p = (char*) ((*array)->data);
    const int stride =  ((*array)->strides)[0]; /* to go to the next row */
    for (i=0; i < nrows; i++, p+=stride) mask[i] = (int*)p;
  }
  else /* We need to create contiguous rows */
  { const char* p0 = (char*) (*array)->data;
    const int rowstride =  (*array)->strides[0];
    const int colstride =  (*array)->strides[1];
    for (i=0; i < nrows; i++)
    { const char* p = p0;
      mask[i] = malloc(ncolumns*sizeof(int));
      for (j=0; j < ncolumns; j++, p+=colstride) mask[i][j] = *((int*)p);
      p0 += rowstride;
    }
  }
  return mask;
}

static void
free_mask(PyArrayObject* array, int** mask, int nrows)
{ int i;
  if (array)
  { if(mask[0]!=(int*)(array->data)) for (i=0; i<nrows; i++) free(mask[i]);
    Py_DECREF((PyObject*) array);
  } else for (i=0; i<nrows; i++) free(mask[i]);
  free(mask);
  return;
}

/* -- weight ---------------------------------------------------------------- */

static double*
parse_weight (PyObject* object, PyArrayObject** array, const int ndata)
{ int i;
  double* weight = NULL;
  if (object==NULL) /* Return the default weights */
  { weight = malloc(ndata*sizeof(double));
    for (i = 0; i < ndata; i++) weight[i] = 1.0;
    *array = NULL;
    return weight;
  }
  if(!PyArray_Check (object)) /* Try to convert object to a 1D double array */
  { *array = (PyArrayObject*) PyArray_FromObject(object, PyArray_DOUBLE, 1, 1);
    if (!(*array))
    { strcpy (message, "weight cannot be converted to needed array.");
      PyErr_SetString(ErrorObject, buffer);
      return NULL;
    }
  }
  else
  { *array = (PyArrayObject*) object;
    Py_INCREF(object);
  }
  if ((*array)->descr->type_num != PyArray_DOUBLE)
  { PyArrayObject* av = (PyArrayObject*)PyArray_Cast(*array, PyArray_DOUBLE);
    Py_DECREF((PyObject*) *array);
    *array = av;
    if (!(*array))
    { strcpy (message, "weight cannot be cast to needed type.");
      PyErr_SetString(ErrorObject, message);
      return NULL;
    }
  }
  if((*array)->nd == 1) /* Checking number of dimensions */
  { /* no checking on last dimension of expected size 1 */
    if (ndata!=1 && ndata!=(*array)->dimensions[0]) 
    { sprintf(message,
              "weight has incorrect extent (%d expected %d)",
              (*array)->dimensions[0], ndata);
      PyErr_SetString (ErrorObject, buffer);
      Py_DECREF(*array);
      *array = NULL;
      return NULL;
    }
  }
  else
  { if ((*array)->nd > 0 || ndata != 1)
    { sprintf(message,
             "weight has incorrect rank (%d expected 1)",
              (*array)->nd);
      PyErr_SetString (ErrorObject, buffer);
      Py_DECREF(*array);
      *array = NULL;
      return NULL;
    }
  }
  /* All checks OK */
  if ((*array)->flags & CONTIGUOUS) weight = (double*) ((*array)->data);
  else
  { const char* p = (char*) ((*array)->data);
    const int stride =  ((*array)->strides)[0];
    weight = malloc(ndata*sizeof(double));
    for (i = 0; i < ndata; i++, p += stride) weight[i] = *(double*)p;
  }
  return weight;
}

static void
free_weight(PyArrayObject* array, double* weight)
{ if (array)
  { if (weight!=(double*)(array->data)) free(weight);
    Py_DECREF((PyObject*) array);
  } else free(weight);
  return;
}

/* -- initialid ------------------------------------------------------------- */

static PyArrayObject*
parse_initialid(PyObject* object, int* nclusters, int nitems)
/* This function creates the clusterid variable for the kcluster and kmedoids
 * routines, and fills it with the initial clustering solution if specified
 * by the user in object. */
{ int i;
  int stride;
  const char* p;
  int* q;
  int* number;
  PyArrayObject* array;
  /* -- First we create the clusterid variable ------------------------ */
  PyArrayObject* clusterid =
    (PyArrayObject*) PyArray_FromDims(1, &nitems, PyArray_LONG);
  if (!clusterid)
  { strcpy(message, "Could not create clusterid array -- too big?");
    PyErr_SetString (ErrorObject, buffer);
    return NULL;
  }
  /* -- If the user didn't specify an initial clustering, we're done -- */
  if (object==NULL) return clusterid;
  /* -- Check if the specified object is an array --------------------- */
  if(!PyArray_Check (object))
  { array = (PyArrayObject*) PyArray_FromObject(object, PyArray_LONG,1,1);
    if (!array)
    { strcpy (message, "initialid cannot be converted to needed array.");
      PyErr_SetString(ErrorObject, buffer);
      Py_DECREF((PyObject*) clusterid);
      return NULL;
    }
  }
  else
  { array = (PyArrayObject*) object;
    Py_INCREF(object);
  }
  /* -- Check if the array contains integers -------------------------- */
  if (array->descr->type_num != PyArray_LONG)
  { PyArrayObject* av = (PyArrayObject*) PyArray_Cast(array, PyArray_LONG);
    Py_DECREF((PyObject*) array);
    array = av;
    if (!array)
    { strcpy (message, "initialid cannot be cast to needed type.");
      PyErr_SetString(ErrorObject, buffer);
      Py_DECREF((PyObject*) clusterid);
      return NULL;
    }
  } 
  /* -- Check the size of the array ----------------------------------- */
  if(array->nd == 1)
  { /* no checking on last dimension of expected size 1 */
    if (nitems!=1 && nitems!=array->dimensions[0]) 
    { sprintf(message, "initialid has incorrect extent (%d expected %d)",
        array->dimensions[0], nitems);
      PyErr_SetString (ErrorObject, buffer);
      Py_DECREF((PyObject*) array);
      Py_DECREF((PyObject*) clusterid);
      return NULL;
    }
  }
  else
  { if (array->nd > 0 || nitems != 1)
    { sprintf(message, "initialid has incorrect rank (%d expected 1)",
        array->nd);
      PyErr_SetString (ErrorObject, buffer);
      Py_DECREF((PyObject*) array);
      Py_DECREF((PyObject*) clusterid);
      return NULL;
    }
  }
  /* -- The array seems to be OK. Count the number of clusters -------- */
  *nclusters = -1;
  stride = array->strides[0];
  p = (const char*) (array->data);
  for (i = 0; i < nitems; i++, p+=stride)
  { const int j = *((int*)p);
    if (j > *nclusters) *nclusters = j;
    if (j < 0)
    { strcpy(message, "initialid contains a negative cluster number");
      PyErr_SetString (ErrorObject, buffer);
      Py_DECREF((PyObject*) array);
      Py_DECREF((PyObject*) clusterid);
      return NULL;
    }
  }
  (*nclusters)++; /* One more than the highest cluster index */
  /* Count the number of items in each cluster */
  number = calloc(*nclusters,sizeof(int));
  p = (const char*) (array->data);
  q = (int*) (clusterid->data);
  for (i = 0; i < nitems; i++, p+=stride, q++)
  { *q = *((int*)p);
    number[*q]++;
  }
  /* Check if any clusters are empty */
  for (i = 0; i < (*nclusters); i++) if(number[i]==0) break;
  free(number);
  Py_DECREF((PyObject*) array);
  if (i < (*nclusters)) /* Due to the break above */
  { sprintf (message, "argument initialid: Cluster %d is empty", i);
    PyErr_SetString (ErrorObject, buffer);
    Py_DECREF((PyObject*) clusterid);
    return NULL;
  }
  return clusterid;
}

/* -- clusterid ------------------------------------------------------------- */

static int*
parse_clusterid(PyObject* object, PyArrayObject** array, unsigned int nitems,
  int* nclusters)
/* This function reads the cluster assignments of all items from object */
{ int i;
  int stride;
  const char* p;
  int* number;
  int* clusterid;
  /* -- Default is to assign all items to the same cluster ------------ */
  if (object==NULL)
  { clusterid = calloc(nitems, sizeof(int));
    *array = NULL;
    *nclusters = 1;
    return clusterid;
  }
  /* -- The user specified something. Let's see if it is an array ----- */
  if(!PyArray_Check (object))
  { *array = (PyArrayObject*) PyArray_FromObject(object, PyArray_LONG,1,1);
    if (!(*array))
    { strcpy (message, "clusterid cannot be converted to needed array.");
      PyErr_SetString(ErrorObject, buffer);
      return NULL;
    }
  }
  else
  { *array = (PyArrayObject*) object;
    Py_INCREF(object);
  }
  /* -- Check if the array contains integers -------------------------- */
  if ((*array)->descr->type_num != PyArray_LONG)
  { PyArrayObject* av = (PyArrayObject*) PyArray_Cast(*array, PyArray_LONG);
    Py_DECREF((PyObject*) (*array));
    *array = av;
    if (!(*array))
    { strcpy (message, "clusterid cannot be cast to needed type.");
      PyErr_SetString(ErrorObject, buffer);
      return NULL;
    }
  } 
  /* -- Check the array size ------------------------------------------ */
  if((*array)->nd == 1)
  { /* no checking on last dimension of expected size 1 */
    if (nitems!=1 && nitems!=(*array)->dimensions[0]) 
    { sprintf(message,
              "clusterid has incorrect extent (%d expected %d)",
              (*array)->dimensions[0], nitems);
      PyErr_SetString (ErrorObject, buffer);
      Py_DECREF((PyObject*) (*array));
      return NULL;
    }
  }
  else
  { if ((*array)->nd > 0 || nitems != 1)
    { sprintf(message,
             "clusterid has incorrect rank (%d expected 1)",
              (*array)->nd);
      PyErr_SetString (ErrorObject, buffer);
      Py_DECREF((PyObject*) (*array));
      return NULL;
    }
  }
  /* -- The array seems to be OK. Count the number of clusters -------- */
  stride = (*array)->strides[0];
  p = (const char*) ((*array)->data);
  *nclusters = -1;
  for (i = 0; i < nitems; i++, p+=stride)
  { const int j = (*(int*)p);
    if (j > *nclusters) *nclusters = j;
    if (j < 0)
    { strcpy(message, "clusterid contains an invalid cluster number");
      PyErr_SetString (ErrorObject, buffer);
      Py_DECREF((PyObject*) (*array));
      return NULL;
    }
  }
  (*nclusters)++;
  /* -- Count the number of items in each cluster --------------------- */
  number = calloc(*nclusters, sizeof(int));
  p = (const char*) ((*array)->data);
  for (i = 0; i < nitems; i++, p+=stride)
  { int j = *((int*)p);
    number[j]++;
  }
  for (i = 0; i < (*nclusters); i++) if(number[i]==0) break;
  free(number);
  if (i < (*nclusters))
  { sprintf (message, "argument initialid: Cluster %d is empty", i);
    PyErr_SetString (ErrorObject, buffer);
    Py_DECREF((PyObject*) (*array));
    return NULL;
  }
  /* All checks OK */
  if ((*array)->flags & CONTIGUOUS) clusterid = (int*) ((*array)->data);
  else
  { const char* p = (char*) ((*array)->data);
    const int stride =  ((*array)->strides)[0];
    clusterid = malloc(nitems*sizeof(int));
    for (i = 0; i < nitems; i++, p += stride) clusterid[i] = *(int*)p;
  }
  return clusterid;
}

static void
free_clusterid(PyArrayObject* array, int* clusterid)
{ if (array)
  { if (clusterid!=(int*)(array->data)) free(clusterid);
    Py_DECREF((PyObject*) array);
  } else free(clusterid);
  return;
}

/* -- distance -------------------------------------------------------------- */

static double**
parse_distance(PyObject* object, PyArrayObject** array, int* n)
/* Takes the Python object from the argument list, and finds the distance
 * matrix. In case of an error, the array is DECREF'ed and set to NULL. */
{ int i, j;
  double** distance = NULL;
  if(!PyArray_Check (object))
  { /* Convert object to a 1D or 2D array of type double */
    *array = (PyArrayObject*) PyArray_FromObject(object, PyArray_DOUBLE, 1, 2);
    if (*array==NULL)
    { strcpy (message, "distance cannot be converted to needed array.");
      PyErr_SetString(ErrorObject, buffer);
      *array = NULL;
      *n = 0;
      return NULL;
    }
  }
  else
  { /* User passed an array */
    *array = (PyArrayObject*) object;
    Py_INCREF(object);
    if ((*array)->descr->type_num != PyArray_DOUBLE)
    { PyArrayObject* av = (PyArrayObject*) PyArray_Cast((*array), PyArray_DOUBLE);
      Py_DECREF((PyObject*) (*array));
      *array = av;
      if (!(*array))
      { strcpy (message, "distance cannot be cast to needed type.");
        PyErr_SetString(ErrorObject, buffer);
        *array = NULL;
        *n = 0;
        return NULL;
      }
    }
  }
  if ((*array)->nd == 1)
  { const int stride =  (*array)->strides[0];
    const char* p = (char*) ((*array)->data);
    const int m = (*array)->dimensions[0];
    *n = (int) ((1+sqrt(1+8*m))/2);
    if ((*n)*(*n)-(*n) != 2 * m)
    { strcpy(message,
        "Array size of distance is incompatible with a lower triangular matrix");
      PyErr_SetString(ErrorObject, buffer);
      Py_DECREF((PyObject*) (*array));
      *array = NULL;
      *n = 0;
      return NULL;
    }
    distance = malloc((*n)*sizeof(double*));
    distance[0] = NULL;
    if (stride==sizeof(double)) /* Data are contiguous */
      for (i=1; i < *n; p+=(i*stride), i++) distance[i] = (double*)p;
    else /* We need to create contiguous rows */
    { for (i=1; i < *n; i++)
      { distance[i] = malloc(i*sizeof(double));
        for (j=0; j < i; j++, p+=stride) distance[i][j] = *((double*)p);
      }
    }
  }
  else if ((*array)->nd == 2)
  { const char* p = (char*) ((*array)->data);
    *n = (*array)->dimensions[0];
    if ((*array)->dimensions[0]!=(*array)->dimensions[1])
    { strcpy(message,
        "The distance matrix should be square");
      PyErr_SetString(ErrorObject, buffer);
      Py_DECREF((PyObject*) (*array));
      *array = NULL;
      *n = 0;
      return NULL;
    }
    distance = malloc((*n)*sizeof(double*));
    distance[0] = NULL;
    if ((*array)->strides[1]==sizeof(double)) /* Each row is contiguous */
    { const int stride =  (*array)->strides[0];
      for (i=0; i < *n; i++, p+=stride) distance[i] = (double*)p;
    }
    else /* We need to create contiguous rows */
    { const int stride =  (*array)->strides[1];
      for (i=0; i < *n; i++)
      { distance[i] = malloc(i*sizeof(double));
        for (j=0; j < i; j++, p+=stride) distance[i][j] = *((double*)p);
      }
    }
  }
  else
  { sprintf(message,
            "distance has an incorrect rank (%d expected 1 or 2)",
            (*array)->nd);
    PyErr_SetString (ErrorObject, buffer);
    Py_DECREF((PyObject*) (*array));
    *array = NULL;
    *n = 0;
    return NULL;
  }
  return distance;
}

static void
free_distances(PyArrayObject* array, double** distance)
{ int i;
  if (array->nd == 1)
  { const int m = array->dimensions[0];
    const int n = (int) ((1+sqrt(1+8*m))/2);
    const int stride =  array->strides[0];
    if (stride!=sizeof(double))
      for (i=1; i < n; i++) free(distance[i]);
  }
  else
  { const int n = array->dimensions[0];
    const int stride =  array->strides[1];
    if (stride!=sizeof(double))
      for (i=1; i < n; i++) free(distance[i]);
  }
  Py_DECREF((PyObject*) array);
  free(distance);
  return;
}

/* -- celldata -------------------------------------------------------------- */

static double***
create_celldata(int nxgrid, int nygrid, int ndata, PyArrayObject** array)
{ int i;
  int shape[3];
  double* p;
  double** pp;
  double*** ppp;
  shape[0] = nxgrid;
  shape[1] = nygrid;
  shape[2] = ndata;
  *array = (PyArrayObject*) PyArray_FromDims(3, shape, PyArray_DOUBLE);
  pp = malloc(nxgrid*nygrid*sizeof(double*));
  ppp = malloc(nxgrid*sizeof(double**));
  if (!(*array) || !pp || !ppp)
  { Py_XDECREF((PyObject*)(*array));
    *array = NULL;
    if(pp) free(pp);
    if(ppp) free(ppp);
    strcpy(message, "Could not create celldata array -- too big?");
    PyErr_SetString (ErrorObject, buffer);
    return NULL;
  }
  p = (double*) ((*array)->data);
  for (i=0; i<nxgrid*nygrid; i++, p+=ndata) pp[i]=p;
  for (i=0; i<nxgrid; i++, pp+=nygrid) ppp[i]=pp;
  return ppp;
}

static void
free_celldata(double*** celldata)
{ double** pp = celldata[0];
  free(pp);
  free(celldata);
}

/* -- index ----------------------------------------------------------------- */

static int*
parse_index(PyObject* object, PyArrayObject** array, int* n)
{ int* index;
  /* Check if the user specified a single item as an integer */
  if(!object || PyInt_Check(object))
  { *array = NULL;
    index = malloc(sizeof(int));
    if (!object) index[0] = 0;
    else index[0] = PyInt_AS_LONG(object);
    *n = 1;
    return index;
  }
  /* Check if the user specified an array */
  if(!PyArray_Check (object)) /* Try to convert to an array of type long */
  { *array = (PyArrayObject*)
      PyArray_ContiguousFromObject(object, PyArray_LONG, 1, 1);
    if (!(*array))
    { strcpy(message, "index argument cannot be converted to needed type.");
      PyErr_SetString (ErrorObject, buffer);
      *n = 0;
      return NULL;
    }
  }
  /* If an array, make sure it contains integers */
  else if ((*array)->descr->type_num == PyArray_LONG)
  { *array = (PyArrayObject*) object;
    Py_INCREF(object);
  }
  else
  { strcpy(message, "index argument cannot be cast to needed type.");
    PyErr_SetString (ErrorObject, buffer);
    *array = NULL;
    *n = 0;
    return NULL;
  }
  /* We have an array */
  *n = (*array)->dimensions[0];
  if((*array)->nd != 1 && ((*array)->nd > 0 || (*array)->dimensions[0] != 1))
  { sprintf(message,
            "index argument has incorrect rank (%d expected 1)",
            (*array)->nd);
    PyErr_SetString (ErrorObject, buffer);
    Py_DECREF(object); /* can only happen if *array==(PyArrayObject*)object */
    *array = NULL;
    *n = 0;
    return NULL;
  }
  if (!(*array)->flags & CONTIGUOUS)
  { PyObject* av =
      PyArray_ContiguousFromObject((PyObject*) array, PyArray_LONG, 0, 0);
    Py_DECREF(object); /* can only happen if *array==(PyArrayObject*)object */
    if(!av)
    { strcpy(message, "Failed making argument index contiguous.");
      PyErr_SetString (ErrorObject, buffer);
      *array = NULL;
      *n = 0;
      return NULL;
    }
    *array = (PyArrayObject*) av;
  }
  index = (int*)((*array)->data);
  return index;
}

static void
free_index(PyArrayObject* array, int* index)
{ if (array) Py_DECREF((PyObject*) array);
  else free(index);
}

/* ========================================================================== */
/* -- Methods --------------------------------------------------------------- */
/* ========================================================================== */

/* kcluster */
static char kcluster__doc__[] =
"returns clusterid, centroids, error, nfound.\n"
"\n"
"This function implements k-means clustering.\n"
"The array data is a nrows x ncolumns array containing the expression data\n"
"The number of clusters is given by nclusters.\n"
"The array mask shows which data are missing. If mask[i][j]==0, then\n"
"data[i][j] is missing.\n"
"The array weight contains the weights to be used when calculating distances.\n"
"If transpose==0, then genes are clustered. If transpose==1, microarrays are\n"
"clustered.\n"
"The integer npass is the number of times the k-means clustering algorithm\n"
"is performed, each time with a different (random) initial condition.\n"
"The character method describes how the center of a cluster is found:\n"
"method=='a': arithmetic mean\n"
"method=='m': median\n"
"The character dist defines the distance function to be used:\n"
"dist=='e': Euclidean distance\n"
"dist=='b': City Block distance\n"
"dist=='h': Harmonically summed Euclidean distance\n"
"dist=='c': Pearson correlation\n"
"dist=='a': absolute value of the correlation\n"
"dist=='u': uncentered correlation\n"
"dist=='x': absolute uncentered correlation\n"
"dist=='s': Spearman's rank correlation\n"
"dist=='k': Kendall's tau\n"
"For other values of dist, the default (Euclidean distance) is used.\n"
"initialid specifies the initial clustering from which the algorithm\n"
"should start. If initialid is not given, the routine carries out npass\n"
"repetitions of the EM algorithm, each time starting from a different\n"
"random initial clustering. If initialid is given, the routine carries\n"
"out the EM algorithm only once, starting from the same initial\n"
"clustering and without randomizing the order in which items are assigned\n"
"to clusters (i.e., using the same order as in the data matrix). In that\n"
"case, the k-means algorithm is fully deterministic.\n"
"\n"
"Return values:\n"
"clusterid is an array containing the number of the cluster to which each\n"
"  gene/microarray was assigned;\n"
"centroids is an array containing the gene expression data for the cluster\n"
"  centroids;\n"
"error is the within-cluster sum of distances for the optimal k-means\n"
"  clustering solution;\n"
"nfound is the number of times the optimal solution was found.\n";

static PyObject*
py_kcluster (PyObject* self, PyObject* args, PyObject* keywords)
{ int NCLUSTERS = 2;
  int nrows, ncolumns;
  int nitems;
  int ndata;
  PyObject* DATA = NULL;
  PyArrayObject* aDATA = NULL;
  double** data = NULL;
  PyObject* MASK = NULL;
  PyArrayObject* aMASK = NULL;
  int** mask = NULL;
  PyObject* WEIGHT = NULL;
  PyArrayObject* aWEIGHT = NULL;
  double* weight = NULL;
  int TRANSPOSE = 0;
  int NPASS = 1;
  char METHOD = 'a';
  char DIST = 'e';
  PyObject* INITIALID = NULL;
  PyArrayObject* aCLUSTERID = NULL;
  PyArrayObject* aCDATA = NULL;
  double** cdata;
  int shape[2];
  double ERROR;
  int IFOUND;
  int i;

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = { "data",
                            "nclusters",
                            "mask",
                            "weight",
                            "transpose",
                            "npass",
                            "method",
                            "dist",
                            "initialid",
                             NULL };
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|lOOllccO", kwlist,
                                  &DATA,
                                  &NCLUSTERS,
                                  &MASK,
                                  &WEIGHT,
                                  &TRANSPOSE,
                                  &NPASS,
                                  &METHOD,
                                  &DIST,
                                  &INITIALID)) return NULL;
  /* Set the function name for error messages */
  strcpy (buffer, "kcluster: ");
  message = strchr(buffer, '\0');
  /* -- Check the method variable ---------------------------------------- */
  if (!strchr("am", METHOD))
  { sprintf(message, "method %c is unknown", METHOD);
    PyErr_SetString (ErrorObject, buffer);
    return NULL;
  }
  /* -- Check the dist variable ------------------------------------------ */
  if (!strchr(known_distances, DIST))
  { sprintf(message, "dist %c is an unknown distance function", DIST);
    PyErr_SetString (ErrorObject, buffer);
    return NULL;
  }
  /* -- Check the transpose variable ------------------------------------- */
  if (TRANSPOSE) TRANSPOSE = 1;
  /* -- Check the npass variable ----------------------------------------- */
  if (INITIALID) NPASS = 0;
  else if (NPASS <= 0)
  { strcpy(message, "npass should be a positive integer");
    PyErr_SetString (ErrorObject, buffer);
    return NULL;
  }
  /* -- Check the data input array --------------------------------------- */
  data = parse_data(DATA, &aDATA);
  if (!data) return NULL;
  nrows = aDATA->dimensions[0];
  ncolumns = aDATA->dimensions[1];
  /* -- Check the mask input --------------------------------------------- */
  mask = parse_mask(MASK, &aMASK, aDATA->dimensions);
  if (!mask)
  { free_data(aDATA, data);
    return NULL;
  }
  /* -- Create the clusterid output variable ----------------------------- */
  ndata = TRANSPOSE ? nrows : ncolumns;
  nitems = TRANSPOSE ? ncolumns : nrows;
  aCLUSTERID = parse_initialid(INITIALID, &NCLUSTERS, nitems);
  if (!aCLUSTERID)
  { free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    return NULL;
  }
  /* -- Check the number of clusters ------------------------------------- */
  if (NCLUSTERS < 1)
  { strcpy(message, "nclusters should be positive");
    PyErr_SetString (ErrorObject, buffer);
    free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    Py_DECREF((PyObject*) aCLUSTERID);
    return NULL;
  }
  if (nitems < NCLUSTERS)
  { strcpy(message, "More clusters than items to be clustered");
    PyErr_SetString (ErrorObject, buffer);
    free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    Py_DECREF((PyObject*) aCLUSTERID);
    return NULL;
  }
  /* -- Check the weight input ------------------------------------------- */
  weight = parse_weight(WEIGHT, &aWEIGHT, ndata);
  if (!weight)
  { free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    Py_DECREF((PyObject*) aCLUSTERID);
    return NULL;
  }
  /* -- Create the centroid data output variable ------------------------- */
  shape[0] = TRANSPOSE ? nrows : NCLUSTERS;
  shape[1] = TRANSPOSE ? NCLUSTERS : ncolumns;
  aCDATA = (PyArrayObject*) PyArray_FromDims(2, shape, PyArray_DOUBLE);
  if (!aCDATA)
  { strcpy(message, "Could not create centroids array -- too big?");
    PyErr_SetString (ErrorObject, buffer);
    free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    free_weight(aWEIGHT, weight);
    Py_DECREF((PyObject*) aCLUSTERID);
    return NULL;
  }
  cdata = malloc(shape[0]*sizeof(double*));
  for (i=0; i<shape[0]; i++)
    cdata[i] = ((double*) (aCDATA->data)) + i*shape[1];
  /* --------------------------------------------------------------------- */
  kcluster(NCLUSTERS, 
      nrows, 
      ncolumns, 
      data, 
      mask, 
      weight,
      TRANSPOSE, 
      NPASS, 
      METHOD, 
      DIST, 
      (int*) (aCLUSTERID->data), 
      cdata,
      &ERROR, 
      &IFOUND);
  /* --------------------------------------------------------------------- */
  free_data(aDATA, data);
  free_mask(aMASK, mask, nrows);
  free_weight(aWEIGHT, weight);
  free (cdata);
  /* --------------------------------------------------------------------- */

  return Py_BuildValue("NNdl",aCLUSTERID, aCDATA, ERROR, IFOUND);
} 
/* end of wrapper for kcluster */

/* kmedoids */
static char kmedoids__doc__[] =
"kmedoids(distance, nclusters=2, npass=1, initialid=None)\n"
"returns clusterid, error, nfound.\n"
"\n"
"This function implements k-medoids clustering.\n"
"The argument distance is a 2D array containing the distance matrix between\n"
"the elements. You can either pass a 2D Numerical Python array (in which\n"
"only the left-lower part of the array will be accessed), or you can pass\n"
"a 1D Numerical Python array containing the distances consecutively.\n" 
"Examples are:\n"
" distance = array([[0.0, 1.1, 2.3],\n"
"                   [1.1, 0.0, 4.5],\n"
"                   [2.3, 4.5, 0.0]])\n"
"and\n"
" distance = array([1.1, 2.3, 4.5])\n"
"These two correspond to the same distance matrix\n"
"The number of clusters is given by nclusters.\n"
"The integer npass is the number of times the k-medoids clustering algorithm\n"
"is performed, each time with a different (random) initial condition.\n"
"The argument initialid specifies the initial clustering from which the\n"
"algorithm should start. If initialid is not given, the routine carries\n"
"out npass repetitions of the EM algorithm, each time starting from a\n"
"different random initial clustering. If initialid is given, the routine\n"
"carries out the EM algorithm only once, starting from the initial\n"
"clustering specified by initialid and without randomizing the order in\n"
"which items are assigned to clusters (i.e., using the same order as in\n"
"the data matrix). In that case, the k-means algorithm is fully\n"
"deterministic.\n"
"\n"
"Return values:\n"
"clusterid is an array containing the number of the cluster to which each\n"
"  gene/microarray was assigned. The cluster number is equal to the number\n"
"  of the element which forms the cluster centroid.\n"
"error is the within-cluster sum of distances for the optimal k-means\n"
"  clustering solution;\n"
"nfound is the number of times the optimal solution was found.\n";

static PyObject*
py_kmedoids (PyObject* self, PyObject* args, PyObject* keywords)
{ int NCLUSTERS = 2;
  int nitems;
  PyObject* DISTANCES = NULL;
  PyArrayObject* aDISTANCES = NULL;
  double** distances = NULL;
  PyObject* INITIALID = NULL;
  PyArrayObject* aCLUSTERID = NULL;
  int NPASS = 1;
  double ERROR;
  int IFOUND;

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = { "distance",
                            "nclusters",
                            "npass",
                            "initialid",
                             NULL };
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|llO", kwlist,
                                  &DISTANCES,
                                  &NCLUSTERS,
                                  &NPASS,
                                  &INITIALID)) return NULL;
  /* Set the function name for error messages */
  strcpy (buffer, "kmedoids: ");
  message = strchr(buffer, '\0');
  /* -- Check the npass variable ----------------------------------------- */
  if (INITIALID) NPASS = 0;
  else if (NPASS < 0)
  { strcpy(message, "npass should be a positive integer");
    PyErr_SetString (ErrorObject, buffer);
    return NULL;
  }
  /* -- Check the distance matrix ---------------------------------------- */
  distances = parse_distance(DISTANCES, &aDISTANCES, &nitems);
  if (!distances) return NULL;
  /* -- Create the clusterid output variable ----------------------------- */
  aCLUSTERID = parse_initialid(INITIALID, &NCLUSTERS, nitems);
  if (!aCLUSTERID)
  { free_distances(aDISTANCES, distances);
    return NULL;
  }
  /* -- Check the nclusters variable ------------------------------------- */
  if (NCLUSTERS <= 0)
  { strcpy(buffer,"nclusters should be a positive integer");
    PyErr_SetString (ErrorObject, buffer);
    free_distances(aDISTANCES, distances);
    Py_DECREF((PyObject*) aCLUSTERID);
    return NULL;
  }
  if (nitems < NCLUSTERS)
  { strcpy(message, "More clusters than items to be clustered");
    PyErr_SetString (ErrorObject, buffer);
    free_distances(aDISTANCES, distances);
    Py_DECREF((PyObject*) aCLUSTERID);
    return NULL;
  }
  /* --------------------------------------------------------------------- */
  kmedoids(NCLUSTERS, 
      nitems, 
      distances, 
      NPASS, 
      (int*) (aCLUSTERID->data), 
      &ERROR, 
      &IFOUND);
  /* --------------------------------------------------------------------- */
  free_distances(aDISTANCES, distances);
  /* --------------------------------------------------------------------- */
  if(!IFOUND)
  { Py_DECREF((PyObject*) aCLUSTERID);
    strcpy(message, "Unknown error in kmedoids");
    return NULL;
  }
  return Py_BuildValue("Ndl",aCLUSTERID, ERROR, IFOUND);
} 
/* end of wrapper for kmedoids */

/* treecluster */
static char treecluster__doc__[] =
"returns tree, linkdist\n"
"\n"
"This function implements the pairwise single, complete, centroid, and\n"
"average linkage hierarchical clustering methods.\n"
"\n"
"The nrows x ncolumns array data contains the gene expression data.\n"
"The array mask declares missing data. If mask[i][j]==0, then data[i][j]\n"
"is missing.\n"
"The array weight contains the weights to be used for the distance\n"
"calculation.\n"
"If the integer applyscale is nonzero, then the distances in linkdist are\n"
"scaled such that all distances are between zero and two (as in case of the\n"
"Pearson distance).\n"
"The integer transpose defines if rows (genes) or columns (microarrays) are\n"
"clustered. If transpose==0, then genes are clustered. If transpose==1,\n"
"microarrays are clustered.\n"
"The character dist defines the distance function to be used:\n"
"dist=='e': Euclidean distance (default)\n"
"dist=='b': City Block distance\n"
"dist=='h': Harmonically summed Euclidean distance\n"
"dist=='c': Pearson correlation\n"
"dist=='a': absolute value of the Pearson correlation\n"
"dist=='u': uncentered correlation\n"
"dist=='x': absolute uncentered correlation\n"
"dist=='s': Spearman's rank correlation\n"
"dist=='k': Kendall's tau\n"
"For other values of dist, the default (Euclidean distance) is used.\n"
"The character method specifies which linkage method is used:\n"
"method=='s': Single pairwise linkage\n"
"method=='m': Complete (maximum) pairwise linkage (default)\n"
"method=='c': Centroid linkage\n"
"method=='a': Average pairwise linkage\n"
"The 2D array distancematrix, which is square and symmetric, is the distance\n"
"matrix. Either data or distancematrix should be None. If distancematrix==None,\n"
"the hierarchical clustering solution is calculated from the gene expression\n"
"data stored in the argument data. If data==None, the hierarchical clustering\n"
"solution is calculated from the distance matrix instead. Pairwise centroid-\n"
"linkage clustering can be calculated only from the gene expression data and\n"
"not from the distance matrix. Pairwise single-, maximum-, and average-linkage\n"
"clustering can be calculated from either the gene expression data or from\n"
"the distance matrix.\n"
"\n"
"Return values:\n"
"tree is an (nobject x 2) array describing the hierarchical clustering\n"
"  result. Each row in the array represents one node, with the two columns\n"
"  representing the two objects or nodes that are being joined. Objects are\n"
"  numbered 0 through (nobjects-1), while nodes are numbered -1 through\n"
"  -(nobjects-1).\n"
"linkdist is a vector with (nobjects-1) elements containing the distances\n"
"between the two subnodes that are joined at each node.\n"
"\n";

static PyObject*
py_treecluster (PyObject* self, PyObject* args, PyObject* keywords)
{ PyObject *DATA = NULL;
  PyObject *MASK = NULL;
  PyObject *WEIGHT = NULL;
  int APPLYSCALE = 0;
  int TRANSPOSE = 0;
  char DIST = 'e';
  char METHOD = 'm';
  PyObject *DISTANCEMATRIX = NULL;
  PyArrayObject* aRESULT = NULL;
  PyArrayObject* aLINKDIST = NULL;

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = { "data",
                            "mask",
                            "weight",
                            "applyscale",
                            "transpose",
                            "method",
                            "dist",
                            "distancematrix",
                             NULL };
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "|OOOllccO", kwlist,
                                  &DATA,
                                  &MASK,
                                  &WEIGHT,
                                  &APPLYSCALE,
                                  &TRANSPOSE,
                                  &METHOD,
                                  &DIST,
                                  &DISTANCEMATRIX)) return NULL;
  /* Set the function name for error messages */
  strcpy (buffer, "treecluster: ");
  message = strchr(buffer, '\0');

  /* -- Check if we are using the data matrix or the distance matrix ----- */
  if (DATA!=NULL && DISTANCEMATRIX!=NULL)
  { strcpy(message, "Use either data or distancematrix, do not use both");
    PyErr_SetString(ErrorObject, buffer);
    return NULL;
  }
  if (DATA==NULL && DISTANCEMATRIX==NULL)
  { strcpy(message, "Neither data nor distancematrix was given");
    PyErr_SetString(ErrorObject, buffer);
    return NULL;
  }

  if (DISTANCEMATRIX==NULL) /* DATA contains gene expression data */
  { int nrows;
    int ncolumns;
    int ndata;
    int nnodes;
    PyArrayObject* aDATA = NULL;
    PyArrayObject* aMASK = NULL;
    PyArrayObject* aWEIGHT = NULL;
    double** data = NULL;
    int** mask = NULL;
    double* weight = NULL;
    int shape[2];
    int (*result)[2];

    /* -- Check the method variable ---------------------------------------- */
    if (!strchr("csma", METHOD))
    { strcpy(message, "keyword method should be 'c', 's', 'm', or 'a'");
      PyErr_SetString(ErrorObject, buffer);
      return NULL;
    }
    /* -- Check the dist variable ------------------------------------------ */
    if (!strchr(known_distances, DIST))
    { sprintf(message, "unknown distance function specified (dist='%c')", DIST);
      PyErr_SetString(ErrorObject, buffer);
      return NULL;
    }
    /* -- Check the data input array --------------------------------------- */
    data = parse_data(DATA, &aDATA);
    if (!data) return NULL;
    nrows = aDATA->dimensions[0];
    ncolumns = aDATA->dimensions[1];
    ndata = TRANSPOSE ? nrows : ncolumns;
    /* -- Check the mask input --------------------------------------------- */
    mask = parse_mask(MASK, &aMASK, aDATA->dimensions);
    if (!mask)
    { free_data(aDATA, data);
      return NULL;
    }
    /* -- Check the weight input ------------------------------------------- */
    weight = parse_weight(WEIGHT, &aWEIGHT, ndata);
    if (!weight)
    { free_data(aDATA, data);
      free_mask(aMASK, mask, nrows);
      return NULL;
    }
    /* -- Create the output variable tree ---------------------------------- */
    nnodes = ((TRANSPOSE==0) ? nrows : ncolumns) - 1;
    shape[0] = nnodes;
    shape[1] = 2;
    aRESULT = (PyArrayObject*) PyArray_FromDims(2, shape, PyArray_LONG);
    if (!aRESULT)
    { strcpy(message, "Could not create array for return value -- too big?");
      PyErr_SetString(ErrorObject, buffer);
      free_data(aDATA, data);
      free_mask(aMASK, mask, nrows);
      free_weight(aWEIGHT, weight);
      return NULL;
    }
    result = (int(*)[2]) (aRESULT->data);
    /* -- Create the output variable linkdist ------------------------------ */
    aLINKDIST = (PyArrayObject*) PyArray_FromDims(1, &nnodes, PyArray_DOUBLE);
    if (!aLINKDIST)
    { strcpy(message, "Could not create array for return value -- too big?");
      PyErr_SetString(ErrorObject, buffer);
      free_data(aDATA, data);
      free_mask(aMASK, mask, nrows);
      free_weight(aWEIGHT, weight);
      Py_DECREF((PyObject*) aRESULT);
    }
    /* --------------------------------------------------------------------- */
    treecluster(nrows,
        ncolumns,
        data,
        mask,
        weight,
        APPLYSCALE,
        TRANSPOSE,
        DIST,
        METHOD,
        result,
        (double*) (aLINKDIST->data), 0);
    /* --------------------------------------------------------------------- */
    free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    free_weight(aWEIGHT, weight);
    /* -- Check if a memory allocation error occurred ---------------------- */
    if(result[0][0]==0 && result[0][1]==0)
    { strcpy(message, "insufficient memory to store the distance matrix");
      PyErr_SetString (ErrorObject, buffer);
      return NULL;
    }
  }
  else
  { double** distances = NULL;
    PyArrayObject* aDISTANCEMATRIX = NULL;
    int nitems;
    int nnodes;
    int shape[2];
    if (!strchr("sma", METHOD))
    { strcpy(message,
        "argument method should be 's', 'm', or 'a' when specifying the distance matrix");
      PyErr_SetString (ErrorObject, buffer);
      return NULL;
    }
    /* -- Check the distance matrix ---------------------------------------- */
    distances = parse_distance(DISTANCEMATRIX, &aDISTANCEMATRIX, &nitems);
    if (!distances) return NULL;
    /* -- Create the output variable tree ---------------------------------- */
    nnodes = nitems - 1;
    shape[0] = nnodes;
    shape[1] = 2;
    aRESULT = (PyArrayObject*) PyArray_FromDims(2, shape, PyArray_LONG);
    if (!aRESULT)
    { strcpy(message, "Could not create array for return value -- too big?");
      PyErr_SetString(ErrorObject, buffer);
      free_distances(aDISTANCEMATRIX, distances);
      return NULL;
    }
    /* -- Create the output variable linkdist ------------------------------ */
    aLINKDIST = (PyArrayObject*) PyArray_FromDims(1, &nnodes, PyArray_DOUBLE);
    if (!aLINKDIST)
    { strcpy(message, "Could not create array for return value -- too big?");
      PyErr_SetString(ErrorObject, buffer);
      free_distances(aDISTANCEMATRIX, distances);
      Py_DECREF((PyObject*) aRESULT);
    }
    /* --------------------------------------------------------------------- */
    treecluster(nitems,
        nitems,
        0, 
        0, 
        0, 
        APPLYSCALE, 
        TRANSPOSE, 
        DIST, 
        METHOD, 
        (int(*)[2]) (aRESULT->data),
        (double*) (aLINKDIST->data),
        distances);
    /* --------------------------------------------------------------------- */
    free_distances(aDISTANCEMATRIX, distances);
    /* --------------------------------------------------------------------- */
  }

  return Py_BuildValue("NN",PyArray_Return(aRESULT),PyArray_Return(aLINKDIST));
} 
/* end of wrapper for treecluster */

/* somcluster */
static char somcluster__doc__[] =
"This function implements a self-organizing map on a rectangular grid.\n"
"The nrows x ncolumns array data contains the measurement data\n"
"The array mask declares missing data. If mask[i][j]==0, then data[i][j]\n"
"is missing.\n"
"The array weights contains the weights to be used for the distance\n"
"calculation.\n"
"The integer transpose defines if rows (genes) or columns (microarrays) are\n"
"clustered. If transpose==0, then genes are clustered. If transpose==1,\n"
"microarrays are clustered.\n"
"The dimensions of the SOM map are nxgrid x nygrid.\n"
"The initial value of tau (the neighborbood function) is given by inittau.\n"
"The number of iterations is given by niter.\n"
"The character dist defines the distance function to be used:\n"
"dist=='e': Euclidean distance\n"
"dist=='b': City Block distance\n"
"dist=='h': Harmonically summed Euclidean distance\n"
"dist=='c': correlation\n"
"dist=='a': absolute value of the correlation\n"
"dist=='u': uncentered correlation\n"
"dist=='x': absolute uncentered correlation\n"
"dist=='s': Spearman's rank correlation\n"
"dist=='k': Kendall's tau\n"
"For other values of dist, the default (Euclidean distance) is used.\n"
"\n"
"Return values:\n"
"clusterid is an array with two columns, while the number of rows is equal to\n"
"  the number of genes or the number of microarrays depending on whether\n"
"  genes or microarrays are being clustered. Each row in the array contains\n"
"  the x and y coordinates of the cell in the rectangular SOM grid to which\n"
"  the gene or microarray was assigned.\n"
"celldata is an array with dimensions (nxgrid, nygrid, number of microarrays)\n"
"  if genes are being clustered, or (nxgrid, nygrid, number of genes) if\n"
"  microarrays are being clustered. Each element [ix][iy] of this array is\n"
"  a 1D vector containing the gene expression data for the centroid of the\n"
"  cluster in the SOM grid cell with coordinates (ix, iy).\n";

static PyObject*
py_somcluster (PyObject* self, PyObject* args, PyObject* keywords)
{ int nrows;
  int ncolumns;
  int nitems;
  int ndata;
  PyObject* DATA = NULL;
  PyArrayObject* aDATA = NULL;
  double** data = NULL;
  PyObject* MASK = NULL;
  PyArrayObject* aMASK = NULL;
  int** mask = NULL;
  PyObject* WEIGHT = NULL;
  PyArrayObject* aWEIGHT = NULL;
  double* weight = NULL;
  int TRANSPOSE = 0;
  int NXGRID = 2;
  int NYGRID = 1;
  double INITTAU = 0.02;
  int NITER = 1;
  char DIST = 'e';
  PyArrayObject* aCELLDATA = NULL;
  double*** celldata = NULL;
  PyArrayObject* aCLUSTERID = NULL;
  int shape[2];

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = { "data",
                            "mask",
                            "weight",
                            "transpose",
                            "nxgrid",
                            "nygrid",
                            "inittau",
                            "niter",
                            "dist",
                             NULL };
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|OOllldlc", kwlist,
                                  &DATA,
                                  &MASK,
                                  &WEIGHT,
                                  &TRANSPOSE,
                                  &NXGRID,
                                  &NYGRID,
                                  &INITTAU,
                                  &NITER,
                                  &DIST)) return NULL;
  /* Set the function name for error messages */
  strcpy (buffer, "somcluster: ");
  message = strchr(buffer, '\0');
  /* -- Check the nxgrid variable ---------------------------------------- */
  if (NXGRID < 1)
  { strcpy(message, "nxgrid should be a positive integer (default is 2)");
    PyErr_SetString (ErrorObject, buffer);
    return NULL;
  }
  /* -- Check the nygrid variable ---------------------------------------- */
  if (NYGRID < 1)
  { strcpy(message, "nygrid should be a positive integer (default is 1)");
    PyErr_SetString (ErrorObject, buffer);
    return NULL;
  }
  /* -- Check the niter variable ----------------------------------------- */
  if (NITER < 1)
  { strcpy(message, "number of iterations (niter) should be positive");
    PyErr_SetString (ErrorObject, buffer);
    return NULL;
  }
  /* -- Check the dist variable ------------------------------------------ */
  if (!strchr(known_distances, DIST))
  { sprintf(message, "dist %c is an unknown distance function", DIST);
    PyErr_SetString (ErrorObject, buffer);
    return NULL;
  }
  /* -- Check the transpose variable ------------------------------------- */
  if (TRANSPOSE) TRANSPOSE = 1;
  /* -- Check the data input array --------------------------------------- */
  data = parse_data(DATA, &aDATA);
  if (!data) return NULL;
  nrows = aDATA->dimensions[0];
  ncolumns = aDATA->dimensions[1];
  nitems = TRANSPOSE ? ncolumns : nrows;
  ndata = TRANSPOSE ? nrows : ncolumns;
  /* -- Check the mask input --------------------------------------------- */
  mask = parse_mask(MASK, &aMASK, aDATA->dimensions);
  if (!mask)
  { free_data(aDATA, data);
    return NULL;
  }
  /* -- Check the weight input ------------------------------------------- */
  weight = parse_weight(WEIGHT, &aWEIGHT, ndata);
  if (!weight)
  { free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    return NULL;
  }
  /* --------------------------------------------------------------------- */
  shape[0] = nitems;
  shape[1] = 2;
  aCLUSTERID = (PyArrayObject*) PyArray_FromDims(2, shape, PyArray_LONG);
  if (!aCLUSTERID)
  { strcpy(buffer, "somcluster: Could not create clusterid array -- too big?");
    PyErr_SetString (ErrorObject, buffer);
    free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    free_weight(aWEIGHT, weight);
    return NULL;
  }
  /* --------------------------------------------------------------------- */
  celldata = create_celldata(NXGRID, NYGRID, ndata, &aCELLDATA);
  if (!celldata)
  { free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    free_weight(aWEIGHT, weight);
    Py_DECREF((PyObject*) aCLUSTERID);
  }
  /* --------------------------------------------------------------------- */
  somcluster(nrows,
      ncolumns,
      data,
      mask,
      weight,
      TRANSPOSE,
      NXGRID,
      NYGRID,
      INITTAU,
      NITER,
      DIST,
      celldata,
      (int(*)[2]) (aCLUSTERID->data));
  /* --------------------------------------------------------------------- */
  free_data(aDATA, data);
  free_mask(aMASK, mask, nrows);
  free_weight(aWEIGHT, weight);
  free_celldata (celldata);
  /* --------------------------------------------------------------------- */
  return Py_BuildValue("NN",
                       PyArray_Return(aCLUSTERID),
                       PyArray_Return(aCELLDATA));
} 
/* end of wrapper for somcluster */

/* median */
static char median__doc__[] =
"median (data)\n"
"This function returns the median of the 1D array data.\n"
"Note: data will be partially ordered upon return.\n";

static PyObject*
py_median (PyObject* unused, PyObject* args)
{ double result;
  PyObject* DATA = NULL;
  PyArrayObject* aDATA = NULL;

  /* -- Read the input variables ----------------------------------------- */
  if(!PyArg_ParseTuple(args, "O", &DATA)) return NULL;

  /* -- Check the input variable ----------------------------------------- */
  if (PyFloat_Check(DATA) || PyInt_Check(DATA) || PyLong_Check(DATA))
  { Py_INCREF(DATA);
    return DATA;
  }
  if(!PyArray_Check (DATA))
  { aDATA = (PyArrayObject *) PyArray_ContiguousFromObject(DATA, PyArray_NOTYPE, 0, 0);
    if (!aDATA)
    { strcpy(buffer, "median: Argument cannot be converted to needed array.");
      PyErr_SetString (ErrorObject, buffer);
      return NULL;
    }
  }
  else
  { aDATA = (PyArrayObject*) DATA;
    Py_INCREF(DATA);
  }
  if (aDATA->descr->type_num != PyArray_DOUBLE)
  { PyObject* av = PyArray_Cast (aDATA, PyArray_DOUBLE);
    Py_DECREF((PyObject*) aDATA);
    aDATA = (PyArrayObject*) av;
    if (!aDATA)
    { strcpy(buffer, "median: Argument cannot be cast to needed type.");
      PyErr_SetString (ErrorObject, buffer);
      return NULL;
    }
  } 
  if (aDATA->nd != 1 && (aDATA->nd > 0 || aDATA->dimensions[0] != 1))
  { sprintf(buffer, "median: Argument has incorrect rank (%d expected 1).",
                    aDATA->nd);
    PyErr_SetString (ErrorObject, buffer);
    Py_DECREF((PyObject*) aDATA);
    return NULL;
  }
  if (!(aDATA->flags & CONTIGUOUS))
  { PyObject* av =
      PyArray_ContiguousFromObject((PyObject*) aDATA, aDATA->descr->type_num, 0, 0);
    Py_DECREF((PyObject*)aDATA);
    if(!av)
    { strcpy(buffer, "median: Failed making argument contiguous.");
      PyErr_SetString (ErrorObject, buffer);
    }
    aDATA = (PyArrayObject*) av;
  }
  /* --------------------------------------------------------------------- */
  result = median(aDATA->dimensions[0], (double*) (aDATA->data));
  /* --------------------------------------------------------------------- */
  Py_DECREF((PyObject*) aDATA);
  /* --------------------------------------------------------------------- */
  return PyFloat_FromDouble(result);
} 
/* end of wrapper for median */

/* mean */
static char mean__doc__[] =
"mean (data)\n"
"This function returns the mean of the 1D array data.\n";

static PyObject*
py_mean (PyObject* unused, PyObject* args)
{ double result;
  PyObject* DATA = NULL;
  PyArrayObject* aDATA = NULL;

  /* -- Read the input variables ----------------------------------------- */
  if(!PyArg_ParseTuple(args, "O", &DATA)) return NULL;

  /* -- Check the input variable ----------------------------------------- */
  if (PyFloat_Check(DATA) || PyInt_Check(DATA) || PyLong_Check(DATA))
  { Py_INCREF(DATA);
    return DATA;
  }
  if(!PyArray_Check (DATA))
  { aDATA = (PyArrayObject *) PyArray_ContiguousFromObject(DATA, PyArray_NOTYPE, 0, 0);
    if (!aDATA)
    { strcpy(buffer, "mean: Argument cannot be converted to needed array.");
      PyErr_SetString (ErrorObject, buffer);
      return NULL;
    }
  }
  else
  { aDATA = (PyArrayObject*) DATA;
    Py_INCREF(DATA);
  }
  if (aDATA->descr->type_num != PyArray_DOUBLE)
  { PyObject* av = PyArray_Cast (aDATA, PyArray_DOUBLE);
    Py_DECREF((PyObject*) aDATA);
    aDATA = (PyArrayObject*) av;
    if (!aDATA)
    { strcpy(buffer, "mean: Argument cannot be cast to needed type.");
      PyErr_SetString (ErrorObject, buffer);
      return NULL;
    }
  } 
  if (aDATA->nd != 1 && (aDATA->nd > 0 || aDATA->dimensions[0] != 1))
  { sprintf(buffer, "mean: Argument has incorrect rank (%d expected 1).",
                    aDATA->nd);
    PyErr_SetString (ErrorObject, buffer);
    Py_DECREF((PyObject*) aDATA);
    return NULL;
  }
  if (!(aDATA->flags & CONTIGUOUS))
  { PyObject* av =
      PyArray_ContiguousFromObject((PyObject*) aDATA, aDATA->descr->type_num, 0, 0);
    Py_DECREF((PyObject*)aDATA);
    if(!av)
    { strcpy(buffer, "mean: Failed making argument contiguous.");
      PyErr_SetString (ErrorObject, buffer);
    }
    aDATA = (PyArrayObject*) av;
  }
  /* --------------------------------------------------------------------- */
  result = mean(aDATA->dimensions[0], (double*) (aDATA->data));
  /* --------------------------------------------------------------------- */
  Py_DECREF((PyObject*) aDATA);
  /* --------------------------------------------------------------------- */
  return PyFloat_FromDouble(result);
} 
/* end of wrapper for mean */

/* clusterdistance */
static char clusterdistance__doc__[] =
"The distance between two clusters\n"
"\n"
"The array data is a nrows x ncolumns array containing the gene expression\n"
"data.\n"
"The array mask shows which data are missing. If mask[i][j]==0, then\n"
"data[i][j] is missing.\n"
"The array weight contains the weights to be used when calculating distances.\n"
"The list index1 identifies which genes/microarrays belong to the first\n"
"cluster. If the cluster contains only one gene, then index1 can also be\n"
"written as a single integer.\n"
"The list index2 identifies which genes/microarrays belong to the second\n"
"cluster. If the cluster contains only one gene, then index2 can also be\n"
"written as a single integer.\n"
"The character dist defines the distance function to be used:\n"
"dist=='e': Euclidean distance\n"
"dist=='b': City Block distance\n"
"dist=='h': Harmonically summed Euclidean distance\n"
"dist=='c': correlation\n"
"dist=='a': absolute value of the correlation\n"
"dist=='u': uncentered correlation\n"
"dist=='x': absolute uncentered correlation\n"
"dist=='s': Spearman's rank correlation\n"
"dist=='k': Kendall's tau\n"
"For other values of dist, the default (Euclidean distance) is used.\n"
"The character method specifies how the distance between two clusters is\n"
"defined:\n"
"method=='a': the distance between the arithmetic means of the two clusters\n"
"method=='m': the distance between the medians of the two clusters\n"
"method=='s': the smallest pairwise distance between members of the two\n"
"             clusters\n"
"method=='x': the largest pairwise distance between members of the two\n"
"             clusters\n"
"method=='v': average of the pairwise distances between members of the\n"
"             clusters\n"
"If transpose==0, then clusters of genes are considered. If transpose==1,\n"
"clusters of microarrays are considered.\n";

static PyObject*
py_clusterdistance (PyObject* self, PyObject* args, PyObject* keywords)
{ double result;
  int nrows;
  int ncolumns;
  int ndata;
  PyObject* DATA = NULL;
  PyArrayObject* aDATA = NULL;
  double** data;
  PyObject* MASK = NULL;
  PyArrayObject* aMASK = NULL;
  int** mask;
  PyObject* WEIGHT = NULL;
  PyArrayObject* aWEIGHT = NULL;
  double* weight;
  char DIST = 'e';
  char METHOD = 'a';
  int TRANSPOSE = 0;
  int N1;
  int N2;
  PyObject* INDEX1 = NULL;
  PyArrayObject* aINDEX1 = NULL;
  int* index1;
  PyObject* INDEX2 = NULL;
  PyArrayObject* aINDEX2 = NULL;
  int* index2;

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = { "data",
                            "mask",
                            "weight",
                            "index1",
                            "index2",
                            "method",
                            "dist",
                            "transpose",
                             NULL };
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|OOOOccl", kwlist,
                                  &DATA,
                                  &MASK,
                                  &WEIGHT,
                                  &INDEX1,
                                  &INDEX2,
                                  &METHOD,
                                  &DIST,
                                  &TRANSPOSE)) return NULL;
  /* Set the function name for error messages */
  strcpy (buffer, "clusterdistance: ");
  message = strchr(buffer, '\0');

  /* -- Check the data input array --------------------------------------- */
  data = parse_data(DATA, &aDATA);
  if (!data) return NULL;
  nrows = aDATA->dimensions[0];
  ncolumns = aDATA->dimensions[1];
  ndata = TRANSPOSE ? nrows : ncolumns;
  /* -- Check the mask input --------------------------------------------- */
  mask = parse_mask(MASK, &aMASK, aDATA->dimensions);
  if (!mask)
  { free_data(aDATA, data);
    return NULL;
  }
  /* -- Check the weight input ------------------------------------------- */
  weight = parse_weight(WEIGHT, &aWEIGHT, ndata);
  if (!weight)
  { free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    return NULL;
  }
  /* --------------------------------------------------------------------- */
  index1 = parse_index(INDEX1, &aINDEX1, &N1);
  if (index1==NULL)
  { free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    free_weight(aWEIGHT, weight);
    return NULL;
  }
  index2 = parse_index(INDEX2, &aINDEX2, &N2);
  if (index2==NULL)
  { free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    free_weight(aWEIGHT, weight);
    free_index(aINDEX1, index1);
    return NULL;
  }
  /* --------------------------------------------------------------------- */
  result = clusterdistance(nrows,
      ncolumns,
      data,
      mask,
      weight,
      N1,
      N2,
      index1,
      index2,
      DIST,
      METHOD,
      TRANSPOSE);
  /* --------------------------------------------------------------------- */
  free_data(aDATA, data);
  free_mask(aMASK, mask, nrows);
  free_weight(aWEIGHT, weight);
  free_index(aINDEX1, index1);
  free_index(aINDEX2, index2);
  /* --------------------------------------------------------------------- */

  return PyFloat_FromDouble(result);
} 
/* end of wrapper for clusterdistance */

/* clustercentroid */
static char clustercentroid__doc__[] =
"The clustercentroid routine calculates the cluster centroids, given to\n"
"which cluster each element belongs. The centroid is defined as either the\n"
"mean or the median over all elements for each dimension.\n"
"The ngenes x nmicroarrays array data contains the gene expression data.\n"
"The array mask declares missing data. If mask[i][j]==0, then data[i][j] is\n"
"missing.\n"
"The integer transpose defines if rows (genes) or columns (microarrays) are\n"
"clustered. If transpose==0, then genes are clustered. If transpose==1,\n"
"microarrays are clustered.\n"
"The array clusterid contains the cluster number for each gene or microarray.\n"
"The cluster number should be non-negative.\n"
"The parameter method specifies whether the centroid is calculated from the\n"
"arithmetic mean (method=='a', default) or the median (method=='m') over each\n"
"dimension.\n"
"This function returns an array cdata and an array cmask.\n"
"The array cdata contains the cluster centroids. If transpose==0, then the\n"
"dimensions of cdata are nclusters x nmicroarrays. If transpose==1, then the\n"
"dimensions of cdata are ngenes x nclusters.\n"
"The array cmask describes which elements in cdata, if any, are missing.\n";

static PyObject*
py_clustercentroid (PyObject* self, PyObject* args, PyObject* keywords)
{ int nrows;
  int ncolumns;
  unsigned int nitems;
  int nclusters;
  PyObject* DATA = NULL;
  PyArrayObject* aDATA = NULL;
  double** data;
  PyObject* MASK = NULL;
  PyArrayObject* aMASK = NULL;
  int** mask;
  PyObject* CLUSTERID = NULL;
  PyArrayObject* aCLUSTERID = NULL;
  int* clusterid;
  char METHOD = 'a';
  int shape[2];
  PyArrayObject* aCDATA = NULL;
  double** cdata;
  PyArrayObject* aCMASK = NULL;
  int** cmask;
  int TRANSPOSE = 0;
  int i;

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = { "data",
                            "mask",
                            "clusterid",
                            "method",
                            "transpose",
                             NULL };
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|OOcl", kwlist,
                                  &DATA,
                                  &MASK,
                                  &CLUSTERID,
                                  &METHOD,
                                  &TRANSPOSE)) return NULL;
  /* Set the function name for error messages */
  strcpy (buffer, "clustercentroid: ");
  message = strchr(buffer, '\0');
  /* -- Check the data input array --------------------------------------- */
  data = parse_data(DATA, &aDATA);
  if (!data) return NULL;
  nrows = aDATA->dimensions[0];
  ncolumns = aDATA->dimensions[1];
  nitems = TRANSPOSE ? ncolumns : nrows;
  /* -- Check the mask input --------------------------------------------- */
  mask = parse_mask(MASK, &aMASK, aDATA->dimensions);
  if (!mask)
  { free_data(aDATA, data);
    return NULL;
  }
  /* -- Check the cluster assignments ------------------------------------ */
  clusterid = parse_clusterid(CLUSTERID, &aCLUSTERID, nitems, &nclusters);
  if (!clusterid)
  { free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    return NULL;
  }
  /* -- Create the centroid data output variable ------------------------- */
  shape[0] = TRANSPOSE ? nrows : nclusters;
  shape[1] = TRANSPOSE ? nclusters : ncolumns;
  aCDATA = (PyArrayObject*) PyArray_FromDims(2, shape, PyArray_DOUBLE);
  if (!aCDATA)
  { strcpy(message, "Could not create centroids array -- too big?");
    PyErr_SetString (ErrorObject, buffer);
    free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    free_clusterid(aCLUSTERID, clusterid);
    return NULL;
  }
  cdata = malloc(shape[0]*sizeof(double*));
  for (i=0; i<shape[0]; i++)
    cdata[i] = ((double*) (aCDATA->data)) + i*shape[1];
  /* -- Create the centroid mask output variable ------------------------- */
  aCMASK = (PyArrayObject*) PyArray_FromDims(2, shape, PyArray_LONG);
  if (!aCMASK)
  { strcpy(message, "Could not create centroids array -- too big?");
    PyErr_SetString (ErrorObject, buffer);
    free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    free_clusterid(aCLUSTERID, clusterid);
    Py_DECREF((PyObject*) aCDATA);
    free (cdata);
    return NULL;
  }
  cmask = malloc(shape[0]*sizeof(int*));
  for (i=0; i<shape[0]; i++)
    cmask[i] = ((int*) (aCMASK->data)) + i*shape[1];
  /* --------------------------------------------------------------------- */
  if (METHOD=='m')
    getclustermedian(nclusters,
        nrows,
        ncolumns,
        data,
        mask,
        clusterid,
        cdata,
        cmask,
        TRANSPOSE);
  else
    getclustermean(nclusters,
        nrows,
        ncolumns,
        data,
        mask,
        clusterid,
        cdata,
        cmask,
        TRANSPOSE);
  /* --------------------------------------------------------------------- */
  free_data(aDATA, data);
  free_mask(aMASK, mask, nrows);
  free (cdata);
  free (cmask);
  free_clusterid(aCLUSTERID, clusterid);
  /* --------------------------------------------------------------------- */
  return Py_BuildValue("NN", PyArray_Return(aCDATA), PyArray_Return(aCMASK));
} 
/* end of wrapper for clustercentroid */

/* distancematrix */
static char distancematrix__doc__[] =
"returns matrix\n"
"\n"
"This function returns the distance matrix between gene expression data.\n"
"The array data is a nrows x ncolumns array containing the expression data\n"
"The array mask shows which data are missing. If mask[i][j]==0, then\n"
"data[i][j] is missing.\n"
"The array weight contains the weights to be used when calculating distances.\n"
"If transpose==0, then genes are clustered. If transpose==1, microarrays are\n"
"clustered.\n"
"The character dist defines the distance function to be used:\n"
"dist=='e': Euclidean distance\n"
"dist=='b': City Block distance\n"
"dist=='h': Harmonically summed Euclidean distance\n"
"dist=='c': Pearson correlation\n"
"dist=='a': absolute value of the correlation\n"
"dist=='u': uncentered correlation\n"
"dist=='x': absolute uncentered correlation\n"
"dist=='s': Spearman's rank correlation\n"
"dist=='k': Kendall's tau\n"
"For other values of dist, the default (Euclidean distance) is used.\n"
"\n"
"Return values:\n"
"matrix is a list of 1D arrays containing the distance matrix between the\n"
"gene expression data. The number of columns in each row is equal to the\n"
"row number. Hence, the first row has zero elements.\n"
"An example of the return value is\n"
"matrix = [[],\n"
"          array([1.]),\n"
"          array([7., 3.]),\n"
"          array([4., 2., 6.])]\n"
"This corresponds to the distance matrix\n"
" [0., 1., 7., 4.]\n"
" [1., 0., 3., 2.]\n"
" [7., 3., 0., 6.]\n"
" [4., 2., 6., 0.]\n";
 

static PyObject*
py_distancematrix (PyObject* self, PyObject* args, PyObject* keywords)
{ PyObject* result = NULL;
  PyObject* DATA = NULL;
  PyArrayObject* aDATA = NULL;
  double** data = NULL;
  PyObject* MASK = NULL;
  PyArrayObject* aMASK = NULL;
  int** mask = (int**) NULL;
  PyObject* WEIGHT = NULL;
  PyArrayObject* aWEIGHT = NULL;
  double* weight = NULL;
  int TRANSPOSE = 0;
  char DIST = 'e';
  double** distances = NULL;
  int nrows, ncolumns, nelements, ndata;
 
  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = { "data",
                            "mask",
                            "weight",
                            "transpose",
                            "dist",
                             NULL };
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|OOlc", kwlist,
                                  &DATA,
                                  &MASK,
                                  &WEIGHT,
                                  &TRANSPOSE,
                                  &DIST)) return NULL;
  /* Set the function name for error messages */
  strcpy (buffer, "distancematrix: ");
  message = strchr(buffer, '\0');
  /* -- Check the dist variable ------------------------------------------ */
  if (!strchr(known_distances, DIST))
  { sprintf(message, "dist %c is an unknown distance function", DIST);
    PyErr_SetString (ErrorObject, buffer);
    return NULL;
  }
  /* -- Check the transpose variable ------------------------------------- */
  if (TRANSPOSE) TRANSPOSE = 1;
  /* -- Check the data input array --------------------------------------- */
  data = parse_data(DATA, &aDATA);
  if (!data) return NULL;
  nrows = aDATA->dimensions[0];
  ncolumns = aDATA->dimensions[1];
  ndata = (TRANSPOSE==0) ? ncolumns : nrows;
  nelements = (TRANSPOSE==0) ? nrows : ncolumns;
  /* -- Check the mask input --------------------------------------------- */
  mask = parse_mask(MASK, &aMASK, aDATA->dimensions);
  if (!mask)
  { free_data(aDATA, data);
    return NULL;
  }
  /* -- Check the weight input ------------------------------------------- */
  weight = parse_weight(WEIGHT, &aWEIGHT, ndata);
  if (!weight)
  { free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    return NULL;
  }
  /* -- Create the matrix output variable -------------------------------- */
  result = PyList_New(nelements);
  if (result)
  { int i, j;
    /* ------------------------------------------------------------------- */
    distances = distancematrix (nrows,
                                ncolumns,
                                data,
                                mask,
                                weight,
                                DIST,
                                TRANSPOSE);
    /* ------------------------------------------------------------------- */
    if (distances)
    { for (i = 0; i < nelements; i++)
      { double* rowdata = NULL;
        PyObject* row = PyArray_FromDims(1, &i, PyArray_DOUBLE);
        if (!row)
        { strcpy(message, "Could not create distance matrix -- too big?");
          PyErr_SetString (ErrorObject, buffer);
          break;
        }
        rowdata = (double*) (((PyArrayObject*)row)->data);
        for (j = 0; j < i; j++) rowdata[j] = distances[i][j];
        free(distances[i]);
        PyList_SET_ITEM(result, i, row);
      }
      if (i < nelements)
      { for (j = 0; j < i; j++)
        { PyObject* row =  PyList_GET_ITEM(result, i);
          Py_DECREF(row);
        }
        for (j = i; j < nelements; j++) free(distances[j]);
        Py_DECREF(result);
        result = NULL;
      }
      free(distances);
    }
    else
    { Py_DECREF(result);
      result = NULL;
    }
  }
  if(result==NULL)
  { strcpy(message, "Could not create distance matrix -- too big?");
    PyErr_SetString (ErrorObject, buffer);
  }
  /* --------------------------------------------------------------------- */
  free_data(aDATA, data);
  free_mask(aMASK, mask, nrows);
  free_weight(aWEIGHT, weight);
  return result;
}

/* cuttree */
static char cuttree__doc__[] =
"clusterid = cuttree(tree, nclusters)\n"
"Given a hierarchical clustering result tree, the routine cuttree divides\n"
"the elements in the tree into clusters. The number of clusters is equal to\n"
"nclusters.\n";

static PyObject*
py_cuttree (PyObject* self, PyObject* args, PyObject* keywords)
{ int NELEMENTS;
  PyObject* TREE;
  PyArrayObject* aTREE = (PyArrayObject*) NULL;
  int NCLUSTERS = 1;
  PyArrayObject* aCLUSTERID = (PyArrayObject*) NULL;
  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = {"ctree", "nclusters", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|l", kwlist,
                                  &TREE, &NCLUSTERS)) return NULL;
  /* -- Check the tree variable (don't allow casting) -------------------- */
  if(!PyArray_Check (TREE))
  { aTREE = (PyArrayObject *) PyArray_ContiguousFromObject(TREE, PyArray_NOTYPE, 0, 0);
    if (!aTREE)
    { PyErr_SetString (ErrorObject,
        "cuttree: Failed converting input argument tree to needed array");
      return NULL;
    }
  }
  else
  { aTREE = (PyArrayObject*) TREE;
    Py_INCREF((PyObject*) aTREE);
  }
  if (aTREE->descr->type_num != PyArray_LONG)
  { PyErr_SetString (ErrorObject,
      "cuttree: Argument tree should contain integer values only");
    return NULL;
  }
  if(aTREE->nd != 2) {
     sprintf(buffer, "cuttree, argument tree: Incorrect rank (%d expected 2)",
                       aTREE->nd);
     PyErr_SetString (ErrorObject, buffer);
     Py_DECREF((PyObject*) aTREE);
     return NULL;
  }
  if (!(aTREE->flags & CONTIGUOUS)) {
    PyObject* av =
      PyArray_ContiguousFromObject((PyObject*) aTREE,
                                   aTREE->descr->type_num, 0, 0);
    Py_DECREF(aTREE);
    if(!av)
    { PyErr_SetString (ErrorObject,
        "cuttree: Failed making input argument tree contiguous");
      return NULL;
    }
    aTREE = (PyArrayObject*) av;
  }
  /* -- Check the nclusters variable ------------------------------------- */
  NELEMENTS = aTREE->dimensions[0] + 1;
  if (NCLUSTERS < 1)
  { PyErr_SetString (ErrorObject,
      "cuttree: Requested number of clusters should be positive");
    Py_DECREF((PyObject*) aTREE);
    return NULL;
  }
  if (NCLUSTERS > NELEMENTS)
  { PyErr_SetString (ErrorObject,
      "cuttree: More clusters requested than items available");
    Py_DECREF((PyObject*) aTREE);
    return NULL;
  }
  /* -- Create the clusterid output variable ----------------------------- */
  aCLUSTERID = (PyArrayObject*) PyArray_FromDims(1, &NELEMENTS, PyArray_LONG);
  if (!aCLUSTERID) {
    PyErr_SetString (ErrorObject,
      "cuttree: Could not create array for return value -- too big?");
    Py_DECREF((PyObject*) aTREE);
    return NULL;
  }
  /* --------------------------------------------------------------------- */
  cuttree(NELEMENTS,
    (int(*)[2]) (aTREE->data),
    NCLUSTERS,
    (int*) (aCLUSTERID->data));
  /* -- The aTREE variable is no longer needed --------------------------- */
  Py_DECREF((PyObject*) aTREE);
  /* -- Check for errors flagged by the C routine ------------------------ */
  if (((int*)(aCLUSTERID->data))[0]==-1)
  {  PyErr_SetString (ErrorObject,
                      "cuttree, argument tree: incompatible input");
     Py_DECREF((PyObject*) aCLUSTERID);
     return NULL;
  }
  /* --------------------------------------------------------------------- */
  return PyArray_Return(aCLUSTERID);
}
/* end of wrapper for cuttree */

/* ========================================================================== */
/* -- The methods table ----------------------------------------------------- */
/* ========================================================================== */


static struct PyMethodDef methods[] = {
   {"kcluster", (PyCFunction) py_kcluster, METH_KEYWORDS, kcluster__doc__},
   {"kmedoids", (PyCFunction) py_kmedoids, METH_KEYWORDS, kmedoids__doc__},
   {"treecluster", (PyCFunction) py_treecluster, METH_KEYWORDS, treecluster__doc__},
   {"somcluster", (PyCFunction) py_somcluster, METH_KEYWORDS, somcluster__doc__},
   {"median", (PyCFunction) py_median, METH_VARARGS, median__doc__},
   {"mean", (PyCFunction) py_mean, METH_VARARGS, mean__doc__},
   {"clusterdistance", (PyCFunction) py_clusterdistance, METH_KEYWORDS, clusterdistance__doc__},
   {"clustercentroid", (PyCFunction) py_clustercentroid, METH_KEYWORDS, clustercentroid__doc__},
   {"distancematrix", (PyCFunction) py_distancematrix, METH_KEYWORDS, distancematrix__doc__},
   {"cuttree", (PyCFunction) py_cuttree, METH_KEYWORDS, cuttree__doc__},
   {NULL,          NULL, 0, NULL}/* sentinel */
};

/* ========================================================================== */
/* -- Initialization -------------------------------------------------------- */
/* ========================================================================== */

void initcluster(void)
{
  PyObject *m, *d;

  import_array ();
  m = Py_InitModule4("cluster",
                     methods,
                     "C Clustering Library",
                     NULL,
                     PYTHON_API_VERSION);
  d = PyModule_GetDict(m);
  ErrorObject = PyString_FromString("cluster.error");
  PyDict_SetItemString(d, "error", ErrorObject);

  if (PyErr_Occurred()) Py_FatalError("can't initialize module cluster");
}
