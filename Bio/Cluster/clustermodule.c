#include "Python.h"
#include "Numeric/arrayobject.h"
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "cluster.h"

static char buffer[512];
static char* message = NULL;

static const char known_distances[] = "ebcauxsk";

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
      PyErr_SetString(PyExc_TypeError, buffer);
      return NULL;
    }
  }
  else /* User passed an array */
  { *array = (PyArrayObject*) object;
    /* Check number of dimensions */
    if ((*array)->nd == 2) Py_INCREF(object);
    else
    { sprintf(message, "data has incorrect rank (%d expected 2)", (*array)->nd);
      PyErr_SetString(PyExc_ValueError, buffer);
      *array = NULL;
      return NULL;
    }
    if ((*array)->descr->type_num != PyArray_DOUBLE) /* Cast to type double */
    { *array = (PyArrayObject*) PyArray_Cast(*array, PyArray_DOUBLE);
      Py_DECREF(object);
      if (!(*array))
      { strcpy (message, "data cannot be cast to needed type.");
        PyErr_SetString(PyExc_ValueError, buffer);
        return NULL;
      }
    } 
  }
  nrows = (*array)->dimensions[0];
  ncols = (*array)->dimensions[1];
  if (nrows < 1 || ncols < 1)
  { strcpy (message, "data is an empty matrix");
    PyErr_SetString(PyExc_ValueError, buffer);
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
  { *array = (PyArrayObject*) PyArray_FromObject(object, PyArray_INT, 2, 2);
    if (!(*array))
    { strcpy (message, "mask cannot be converted to needed array");
      PyErr_SetString(PyExc_TypeError, buffer);
      return NULL;
    }
  }
  else /* User passed an array */
  { *array = (PyArrayObject*) object;
    if((*array)->nd != 2) /* Checking number of dimensions */
    { sprintf(message, "mask has incorrect rank (%d expected 2)", (*array)->nd);
      PyErr_SetString(PyExc_ValueError, buffer);
      *array = NULL;
      return NULL;
    }
    if ((*array)->descr->type_num == PyArray_INT) Py_INCREF(object);
    else
    { *array = (PyArrayObject*) PyArray_Cast (*array, PyArray_INT);
      if (!(*array))
      { strcpy (message, "mask cannot be cast to needed type.");
        PyErr_SetString(PyExc_ValueError, buffer);
        return NULL;
      }
    } 
  }
  if((*array)->dimensions[0] != nrows) /* Checking number of rows */
  { sprintf(message,
      "mask has incorrect number of rows (%d expected %d)",
      (*array)->dimensions[0], nrows);
    PyErr_SetString (PyExc_ValueError, buffer);
    Py_DECREF((PyObject*)*array);
    *array = NULL;
    return NULL;
  }
  /* no checking on last dimension of expected size 1 */
  if (ncolumns != 1 && (*array)->dimensions[1] != ncolumns)
  { sprintf(message,
      "mask incorrect number of columns (%d expected %d)",
      (*array)->dimensions[1], ncolumns);
    PyErr_SetString (PyExc_ValueError, buffer);
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
      PyErr_SetString(PyExc_TypeError, buffer);
      return NULL;
    }
  }
  else
  { *array = (PyArrayObject*) object;
    if ((*array)->descr->type_num == PyArray_DOUBLE) Py_INCREF(object);
    else
    { *array = (PyArrayObject*)PyArray_Cast(*array, PyArray_DOUBLE);
      if (!(*array))
      { strcpy (message, "weight cannot be cast to needed type.");
        PyErr_SetString(PyExc_ValueError, message);
        return NULL;
      }
    }
  }
  if((*array)->nd == 1) /* Checking number of dimensions */
  { /* no checking on last dimension of expected size 1 */
    if (ndata!=1 && ndata!=(*array)->dimensions[0]) 
    { sprintf(message,
              "weight has incorrect extent (%d expected %d)",
              (*array)->dimensions[0], ndata);
      PyErr_SetString (PyExc_ValueError, buffer);
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
      PyErr_SetString (PyExc_ValueError, buffer);
      Py_DECREF(*array);
      *array = NULL;
      return NULL;
    }
  }
  /* All checks OK */
  if (PyArray_ISCONTIGUOUS(*array)) weight = (double*) ((*array)->data);
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
    (PyArrayObject*) PyArray_FromDims(1, &nitems, PyArray_INT);
  if (!clusterid)
  { strcpy(message, "could not create clusterid array");
    PyErr_SetString (PyExc_MemoryError, buffer);
    return NULL;
  }
  /* -- If the user didn't specify an initial clustering, we're done -- */
  if (object==NULL) return clusterid;
  /* -- Check if the specified object is an array --------------------- */
  if(!PyArray_Check (object))
  { array = (PyArrayObject*) PyArray_FromObject(object, PyArray_INT,1,1);
    if (!array)
    { strcpy (message, "initialid cannot be converted to needed array.");
      PyErr_SetString(PyExc_TypeError, buffer);
      Py_DECREF((PyObject*) clusterid);
      return NULL;
    }
  }
  else
  { array = (PyArrayObject*) object;
    /* -- Check if the array contains integers ------------------------ */
    if (array->descr->type_num == PyArray_INT) Py_INCREF(object);
    else
    { array = (PyArrayObject*) PyArray_Cast(array, PyArray_INT);
      if (!array)
      { strcpy (message, "initialid cannot be cast to needed type.");
        PyErr_SetString(PyExc_ValueError, buffer);
        Py_DECREF((PyObject*) clusterid);
        return NULL;
      }
    } 
  }
  /* -- Check the size of the array ----------------------------------- */
  if(array->nd == 1)
  { /* no checking on last dimension of expected size 1 */
    if (nitems!=1 && nitems!=array->dimensions[0]) 
    { sprintf(message, "initialid has incorrect extent (%d expected %d)",
        array->dimensions[0], nitems);
      PyErr_SetString (PyExc_ValueError, buffer);
      Py_DECREF((PyObject*) array);
      Py_DECREF((PyObject*) clusterid);
      return NULL;
    }
  }
  else
  { if (array->nd > 0 || nitems != 1)
    { sprintf(message, "initialid has incorrect rank (%d expected 1)",
        array->nd);
      PyErr_SetString (PyExc_ValueError, buffer);
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
      PyErr_SetString (PyExc_ValueError, buffer);
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
    PyErr_SetString (PyExc_ValueError, buffer);
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
  { *array = (PyArrayObject*) PyArray_FromObject(object, PyArray_INT,1,1);
    if (!(*array))
    { strcpy (message, "clusterid cannot be converted to needed array.");
      PyErr_SetString(PyExc_TypeError, buffer);
      return NULL;
    }
  }
  else
  { *array = (PyArrayObject*) object;
    /* -- Check if the array contains integers ------------------------ */
    if ((*array)->descr->type_num == PyArray_INT) Py_INCREF(object);
    else
    { *array = (PyArrayObject*) PyArray_Cast(*array, PyArray_INT);
      if (!(*array))
      { strcpy (message, "clusterid cannot be cast to needed type.");
        PyErr_SetString(PyExc_ValueError, buffer);
        return NULL;
      }
    } 
  }
  /* -- Check the array size ------------------------------------------ */
  if((*array)->nd == 1)
  { /* no checking on last dimension of expected size 1 */
    if (nitems!=1 && nitems!=(*array)->dimensions[0]) 
    { sprintf(message,
              "clusterid has incorrect extent (%d expected %d)",
              (*array)->dimensions[0], nitems);
      PyErr_SetString (PyExc_ValueError, buffer);
      Py_DECREF((PyObject*) (*array));
      return NULL;
    }
  }
  else if ((*array)->nd > 0 || nitems != 1)
  { sprintf(message,
           "clusterid has incorrect rank (%d expected 1)",
            (*array)->nd);
    PyErr_SetString (PyExc_ValueError, buffer);
    Py_DECREF((PyObject*) (*array));
    return NULL;
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
      PyErr_SetString (PyExc_ValueError, buffer);
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
    PyErr_SetString (PyExc_ValueError, buffer);
    Py_DECREF((PyObject*) (*array));
    return NULL;
  }
  /* All checks OK */
  if (PyArray_ISCONTIGUOUS(*array)) clusterid = (int*) ((*array)->data);
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

static void
free_distances(PyObject* object, PyArrayObject* array, double** distance, int n)
{ int i;
   return;
  if (array==NULL) /* User passed a lower-triangular matrix as a list of rows */
  { for (i = 1; i < n; i++)
    { PyObject* row = PyList_GET_ITEM(object, i);
      if (PyArray_Check(row))
      { PyArrayObject* a = (PyArrayObject*)row;
        if (distance[i] == (double*) (a->data))
        { Py_DECREF(row);
          continue;
        }
      }
      free(distance[i]);
    }
  }
  else
  { if (array->nd == 1)
    { const int stride =  array->strides[0];
      if (stride!=sizeof(double))
        for (i = 1; i < n; i++) free(distance[i]);
    }
    else
    { const int stride =  array->strides[1];
      if (stride!=sizeof(double))
        for (i = 1; i < n; i++) free(distance[i]);
    }
    Py_DECREF((PyObject*) array);
  }
  free(distance);
}

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
    { /* This is not necessarily an error; the user may have passed the
       * the lower-triangular matrix as a list of rows. Clear the error
       * indicator set by PyArrayFromObject first. */
      PyErr_Clear();
      if (!PyList_Check(object))
      { strcpy (message, "distance cannot be converted to needed array.");
        PyErr_SetString(PyExc_TypeError, buffer);
        *n = 0;
        return NULL;
      }
      *n = PyList_GET_SIZE(object);
      distance = malloc((*n)*sizeof(double*));
      if (!distance)
      { strcpy (message, "failed to store distance matrix.");
        PyErr_SetString(PyExc_MemoryError, buffer);
        *n = 0;
        return NULL;
      }
      for (i = 0; i < *n; i++)
      { PyObject* row = PyList_GET_ITEM(object, i);
        if (PyArray_Check(row))
        { PyArrayObject* a = (PyArrayObject*)row;
          if (a->nd != 1)
          { sprintf (message, "Row %d in the distance matrix is not one-dimensional.", i);
            PyErr_SetString(PyExc_ValueError, buffer);
            break;
          }
          if (a->dimensions[0]!=i)
          { sprintf (message, "Row %d in the distance matrix has incorrect size (%d, should be %d).", i, a->dimensions[0], i);
            PyErr_SetString(PyExc_ValueError, buffer);
            break;
          }
          if (i==0) continue;
          if (a->descr->type_num == PyArray_DOUBLE)
          { const int stride = a->strides[0];
            if (stride==sizeof(double)) /* Row is contiguous */
            { Py_INCREF(row);
              distance[i] = (double*) (a->data);
            }
            else
            { const char* p = (char*) (a->data);
              distance[i] = malloc(i*sizeof(double));
              if(!distance[i])
              { Py_DECREF((PyObject*)a);
                sprintf (message, "failed to store row %d in the distance matrix.", i);
                PyErr_SetString(PyExc_MemoryError, buffer);
                break;
              }
              for (j=0; j < i; j++, p+=stride) distance[i][j] = *((double*)p);
            }
          }
          else
          { row = PyArray_ContiguousFromObject(row, PyArray_DOUBLE, 1, 1);
            if (!row)
            { sprintf (message, "Failed to cast row %d in the distance matrix to double precision.", i);
              PyErr_SetString(PyExc_MemoryError, buffer);
              break;
            }
            else
            { const double* p;
              a = (PyArrayObject*)row;
              p = (double*) (a->data);
              distance[i] = malloc(i*sizeof(double));
              if(!distance[i])
              { Py_DECREF(row);
                sprintf (message, "failed to store row %d in the distance matrix.", i);
                PyErr_SetString(PyExc_MemoryError, buffer);
                break;
              }
              for (j=0; j < i; j++, p++) distance[i][j] = *p;
              Py_DECREF(row);
            }
          }
        }
        else
        { /* Convert row */
          const double* p;
          PyArrayObject* a = (PyArrayObject*)PyArray_ContiguousFromObject(row, PyArray_DOUBLE, 1, 1);
          if(!a)
          { sprintf (message, "Failed to convert row %d in the distance matrix.", i);
            PyErr_SetString(PyExc_TypeError, buffer);
            break;
          }
          if (a->dimensions[0]!=i)
          { sprintf (message, "Row %d in the distance matrix has incorrect size (%d, should be %d).", i, a->dimensions[0], i);
            PyErr_SetString(PyExc_ValueError, buffer);
            Py_DECREF((PyObject*)a);
            break;
          }
          if (i > 0)
          { distance[i] = malloc(i*sizeof(double));
            if(!distance[i])
            { Py_DECREF((PyObject*)a);
              sprintf (message, "failed to store row %d in the distance matrix.", i);
              PyErr_SetString(PyExc_MemoryError, buffer);
              break;
            }
            p = (double*) (a->data);
            for (j=0; j < i; j++) distance[i][j] = p[j];
          }
          Py_DECREF((PyObject*)a);
        }
      }
      if (i < *n) /* break encountered */
      { free_distances(object, NULL, distance, i);
        *n = 0;
        return NULL;
      }
      return distance;
    }
  }
  else
  { /* User passed an array */
    *array = (PyArrayObject*) object;
    if ((*array)->descr->type_num == PyArray_DOUBLE) Py_INCREF(object);
    else
    { *array = (PyArrayObject*) PyArray_Cast((*array), PyArray_DOUBLE);
      if (!(*array))
      { strcpy (message, "distance cannot be cast to needed type.");
        PyErr_SetString(PyExc_ValueError, buffer);
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
      PyErr_SetString(PyExc_ValueError, buffer);
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
    if ((*n) != (*array)->dimensions[1])
    { strcpy(message,
        "The distance matrix should be square");
      PyErr_SetString(PyExc_ValueError, buffer);
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
    PyErr_SetString (PyExc_ValueError, buffer);
    Py_DECREF((PyObject*) (*array));
    *array = NULL;
    *n = 0;
    return NULL;
  }
  return distance;
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
    PyErr_SetString(PyExc_MemoryError, buffer);
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
    else index[0] = (int) PyInt_AS_LONG(object);
    *n = 1;
    return index;
  }
  /* Check if the user specified an array */
  if(!PyArray_Check (object)) /* Try to convert to an array of type int */
  { *array = (PyArrayObject*)
      PyArray_ContiguousFromObject(object, PyArray_INT, 1, 1);
    if (!(*array))
    { strcpy(message, "index argument cannot be converted to needed type.");
      PyErr_SetString (PyExc_TypeError, buffer);
      *n = 0;
      return NULL;
    }
  }
  else
  { *array = (PyArrayObject*) object;
    /* -- Check if the array contains integers ------------------------ */
    if ((*array)->descr->type_num == PyArray_INT) Py_INCREF(object);
    else
    { object = PyArray_Cast(*array, PyArray_INT);
      if (!object)
      { strcpy (message, "index argument cannot be cast to needed type.");
        PyErr_SetString(PyExc_ValueError, buffer);
        *n = 0;
        return NULL;
      }
      *array = (PyArrayObject*) object;
    } 
  }
  /* We have an array */
  *n = (*array)->dimensions[0];
  if((*array)->nd != 1 && ((*array)->nd > 0 || (*array)->dimensions[0] != 1))
  { sprintf(message,
            "index argument has incorrect rank (%d expected 1)",
            (*array)->nd);
    PyErr_SetString (PyExc_ValueError, buffer);
    Py_DECREF(object); /* can only happen if *array==(PyArrayObject*)object */
    *array = NULL;
    *n = 0;
    return NULL;
  }
  if (!PyArray_ISCONTIGUOUS(*array))
  { *array = (PyArrayObject*) PyArray_ContiguousFromObject(object, PyArray_INT, 1, 1);
    Py_DECREF(object);
    if(!(*array))
    { strcpy(message, "Failed making argument index contiguous.");
      PyErr_SetString (PyExc_ValueError, buffer);
      *array = NULL;
      *n = 0;
      return NULL;
    }
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
/* -- Classes --------------------------------------------------------------- */
/* ========================================================================== */

typedef struct {
    PyObject_HEAD
    Node node;
} PyNode;

static int
PyNode_init(PyNode *self, PyObject *args, PyObject *kwds)
{
    int left, right;
    double distance = 0.0;
    static char *kwlist[] = {"left", "right", "distance", NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "ii|d", kwlist, 
                                      &left, &right, &distance))
        return -1; 
    self->node.left = left;
    self->node.right = right;
    self->node.distance = distance;

    return 0;
}

static PyObject*
PyNode_repr(PyNode* self)
{ char string[64];
  sprintf(string, "(%d, %d): %g",
                  self->node.left, self->node.right, self->node.distance);
  return PyString_FromString(string);
}

static char PyNode_left__doc__[] =
"integer representing the first member of this node";

static PyObject*
PyNode_getleft(PyNode* self, void* closure)
{ int left = self->node.left;
  return PyInt_FromLong((long)left);
}

static int
PyNode_setleft(PyNode* self, PyObject* value, void* closure)
{ long left = PyInt_AsLong(value);
  if (PyErr_Occurred()) return -1;
  self->node.left = (int) left;
  return 0;
}

static char PyNode_right__doc__[] =
"integer representing the second member of this node";

static PyObject*
PyNode_getright(PyNode* self, void* closure)
{ int right = self->node.right;
  return PyInt_FromLong((long)right);
}

static int
PyNode_setright(PyNode* self, PyObject* value, void* closure)
{ long right = PyInt_AsLong(value);
  if (PyErr_Occurred()) return -1;
  self->node.right = (int) right;
  return 0;
}

static PyObject*
PyNode_getdistance(PyNode* self, void* closure)
{ return PyFloat_FromDouble(self->node.distance);
}

static int
PyNode_setdistance(PyNode* self, PyObject* value, void* closure)
{ const double distance = PyFloat_AsDouble(value);
  if (PyErr_Occurred()) return -1;
  self->node.distance = distance;
  return 0;
}

static char PyNode_distance__doc__[] =
"the distance between the two members of this node\n";

static PyGetSetDef PyNode_getset[] = {
    {"left", (getter)PyNode_getleft, (setter)PyNode_setleft, PyNode_left__doc__, NULL},
    {"right", (getter)PyNode_getright, (setter)PyNode_setright, PyNode_right__doc__, NULL},
    {"distance", (getter)PyNode_getdistance, (setter)PyNode_setdistance, PyNode_distance__doc__, NULL},
    {NULL}  /* Sentinel */
};

static char PyNode_doc[] =
"A Node object describes a single node in a hierarchical clustering tree.\n"
"The integer attributes 'left' and 'right' represent the two members that\n"
"make up this node; the floating point attribute 'distance' contains the\n"
"distance between the two members of this node.\n";

static PyTypeObject PyNodeType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "cluster.Node",            /*tp_name*/
    sizeof(PyNode),            /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    0,                         /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    (reprfunc)PyNode_repr,     /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    PyNode_doc,                /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    0,                         /* tp_methods */
    0,                         /* tp_members */
    PyNode_getset,             /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)PyNode_init,     /* tp_init */
};

typedef struct {
    PyObject_HEAD
    Node* nodes;
    int n;
} PyTree;

static void
PyTree_dealloc(PyTree* self)
{ if (self->n) free(self->nodes);
  self->ob_type->tp_free((PyObject*)self);
}

static int
PyTree_init(PyTree* self, PyObject* args, PyObject* kwds)
{ int i;
  int n;
  Node* nodes;
  PyObject* arg;
  int* flag;

  if (! PyArg_ParseTuple(args, "O", &arg)) return -1;

  if (!PyList_Check(arg))
  { PyErr_SetString(PyExc_TypeError, "Argument should be a list of Node objects");
    return -1;
  }

  n = PyList_GET_SIZE(arg);
  if (n < 1)
  { PyErr_SetString(PyExc_ValueError, "List is empty");
    return -1;
  }
  nodes = malloc(n*sizeof(Node));
  for (i = 0; i < n; i++)
  { PyNode* p;
    PyObject* row = PyList_GET_ITEM(arg, i);
    if (row->ob_type != &PyNodeType)
    { free(nodes);
      sprintf(buffer, "Row %d in list is not a Node object", i);
      PyErr_SetString(PyExc_TypeError, buffer);
      return -1;
    }
    p = (PyNode*)row;
    nodes[i] = p->node;
  }
  /* --- Check if this is a bona fide tree ------------------------------- */
  flag = malloc((2*n+1)*sizeof(int));
  if(flag) /* Otherwise, we're in enough trouble already */
  { int j;
    for (i = 0; i < 2*n+1; i++) flag[i] = 0; 
    for (i = 0; i < n; i++)
    { j = nodes[i].left;
      if (j < 0)
      { j = -j-1;
        if (j>=i) break;
      }
      else j+=n;
      if (flag[j]) break;
      flag[j] = 1;
      j = nodes[i].right;
      if (j < 0)
      { j = -j-1;
        if (j>=i) break;
      }
      else j+=n;
      if (flag[j]) break;
      flag[j] = 1;
    }
    free(flag);
  }
  if (!flag || i < n) /* break encountered */
  { free(nodes);
    PyErr_SetString(PyExc_ValueError, "Inconsistent tree");
    return -1;
  }
  /* --------------------------------------------------------------------- */
  self->n = n;
  self->nodes = nodes;
  return 0;
}

static PyObject*
PyTree_str(PyTree* self)
{ int i;
  const int n = self->n;
  Node node;
  PyObject* line;
  PyObject* output = PyString_FromString("");
  char string[128];
  for (i = 0; i < n; i++)
  { node = self->nodes[i];
    sprintf(string, "(%d, %d): %g", node.left, node.right, node.distance);
    if (i < n-1) strcat(string, "\n");
    line = PyString_FromString(string);
    if(!line)
    { Py_DECREF(output);
      return NULL;
    }
    PyString_ConcatAndDel(&output, line);
    if(!output) return NULL;
  }
  return output;
}

static int
PyTree_length(PyTree *self)
{
  return self->n;
}

static PyObject*
PyTree_item(PyTree* self, int i)
{ PyNode* result;
  if (i < 0 || i >= self->n)
  { PyErr_SetString(PyExc_IndexError, "tree index out of range");
    return NULL;
  }
  result = (PyNode*) PyNodeType.tp_alloc(&PyNodeType, 0);
  if(!result)
  { PyErr_SetString(PyExc_MemoryError,
      "could not create node for return value");
    return NULL;
  }
  result->node = self->nodes[i];
  return (PyObject*) result;
}

static PyObject*
PyTree_slice(PyTree* self, int i, int j)
{ int row;
  const int n = self->n;
  PyObject* item;
  PyObject* result;
  if (i < 0) i = 0;
  if (j < 0) j = 0; /* Avoid signed/unsigned bug in next line */
  if (j > n) j = n;
  if (j < i) j = i;
  result = PyList_New(j-i);
  if(!result)
  { PyErr_SetString(PyExc_MemoryError,
      "could not create list for return value");
    return NULL;
  }
  for (row = 0; i < j; i++, row++)
  { item = PyTree_item(self, i);
    if(!item)
    { Py_DECREF(result);
      PyErr_SetString(PyExc_MemoryError,
        "could not create node for return value");
      return NULL;
    }
    PyList_SET_ITEM(result, row, item);
  }
  return result;
}

#if (PY_MAJOR_VERSION <= 2) & (PY_MINOR_VERSION <= 4)
#define lenfunc inquiry
#define ssizeargfunc intargfunc
#define ssizessizeargfunc intintargfunc
#endif

static PySequenceMethods PyTree_sequence = {
        (lenfunc)PyTree_length, /*sq_length*/
        NULL,                 /*sq_concat*/
        NULL,                 /*sq_repeat*/
        (ssizeargfunc)PyTree_item, /*sq_item*/
        (ssizessizeargfunc)PyTree_slice, /*sq_slice*/
        NULL,                 /*sq_ass_item*/
        NULL,                 /*sq_ass_slice*/
        NULL                  /*sq_contains*/
};

static char PyTree_scale__doc__[] =
"mytree.scale()\n"
"This method scales the node distances in the tree such that they are all\n"
"between one and zero.\n";

static PyObject*
PyTree_scale(PyTree* self)
{ int i;
  const int n = self->n;
  Node* nodes = self->nodes;
  double maximum = DBL_MIN;
  /* --------------------------------------------------------------------- */
  for (i = 0; i < n; i++)
  { double distance = nodes[i].distance;
    if (distance > maximum) maximum = distance;
  }
  if (maximum!=0.0)
    for (i = 0; i < n; i++) nodes[i].distance /= maximum;
  /* --------------------------------------------------------------------- */
  Py_INCREF(Py_None);
  return Py_None;
}
/* end of wrapper for cuttree */
static char PyTree_cut__doc__[] =
"clusterid = mytree.cut(nclusters=2)\n"
"Given a hierarchical clustering result mytree, cut() divides the elements\n"
"in the tree into clusters. The number of clusters is given by nclusters.\n";

static PyObject*
PyTree_cut(PyTree* self, PyObject* args)
{ int nclusters = 2;
  int n = self->n + 1;
  PyArrayObject* aCLUSTERID = (PyArrayObject*) NULL;
  int* clusterid = NULL;
  /* -- Read the input variables ----------------------------------------- */
  if(!PyArg_ParseTuple(args, "|i", &nclusters)) return NULL;
  /* -- Check the nclusters variable ------------------------------------- */
  if (nclusters < 1)
  { PyErr_SetString (PyExc_ValueError,
      "cut: Requested number of clusters should be positive");
    return NULL;
  }
  if (nclusters > n)
  { PyErr_SetString (PyExc_ValueError,
      "cut: More clusters requested than items available");
    return NULL;
  }
  /* -- Create the clusterid output variable ----------------------------- */
  aCLUSTERID = (PyArrayObject*) PyArray_FromDims(1, &n, PyArray_INT);
  if (!aCLUSTERID)
  { PyErr_SetString (PyExc_MemoryError,
      "cut: Could not create array for return value");
    return NULL;
  }
  clusterid = (int*) (aCLUSTERID->data);
  /* --------------------------------------------------------------------- */
  cuttree(n, self->nodes, nclusters, clusterid);
  /* -- Check for errors flagged by the C routine ------------------------ */
  if (clusterid[0]==-1)
  { PyErr_SetString (PyExc_MemoryError, "cut: error in the cuttree routine");
    Py_DECREF((PyObject*) aCLUSTERID);
    return NULL;
  }
  /* --------------------------------------------------------------------- */
  return PyArray_Return(aCLUSTERID);
}
/* end of wrapper for cuttree */

static PyMethodDef PyTree_methods[] = {
    {"scale", (PyCFunction)PyTree_scale, METH_NOARGS, PyTree_scale__doc__},
    {"cut", (PyCFunction)PyTree_cut, METH_VARARGS, PyTree_cut__doc__},
    {NULL}  /* Sentinel */
};

static char PyTree_doc[] =
"Tree objects store a hierarchical clustering solution.\n"
"Individual nodes in the tree can be accessed with tree[i], where i is\n"
"an integer. Whereas the tree itself is a read-only object, tree[:]\n"
"returns a list of all the nodes, which can then be modified. To create\n"
"a new Tree from this list, use Tree(list).\n"
"See the description of the Node class for more information.";

static PyTypeObject PyTreeType = {
    PyObject_HEAD_INIT(NULL)
    0,                           /*ob_size*/
    "cluster.Tree",              /*tp_name*/
    sizeof(PyTree),              /*tp_basicsize*/
    0,                           /*tp_itemsize*/
    (destructor)PyTree_dealloc,  /*tp_dealloc*/
    0,                           /*tp_print*/
    0,                           /*tp_getattr*/
    0,                           /*tp_setattr*/
    0,                           /*tp_compare*/
    0,                           /*tp_repr*/
    0,                           /*tp_as_number*/
    &PyTree_sequence,            /*tp_as_sequence*/
    0,                           /*tp_as_mapping*/
    0,                           /*tp_hash */
    0,                           /*tp_call*/
    (reprfunc)PyTree_str,        /*tp_str*/
    0,                           /*tp_getattro*/
    0,                           /*tp_setattro*/
    0,                           /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,          /*tp_flags*/
    PyTree_doc,                  /* tp_doc */
    0,		                 /* tp_traverse */
    0,		                 /* tp_clear */
    0,		                 /* tp_richcompare */
    0,		                 /* tp_weaklistoffset */
    0,		                 /* tp_iter */
    0,		                 /* tp_iternext */
    PyTree_methods,              /* tp_methods */
    NULL,                        /* tp_members */
    0,                           /* tp_getset */
    0,                           /* tp_base */
    0,                           /* tp_dict */
    0,                           /* tp_descr_get */
    0,                           /* tp_descr_set */
    0,                           /* tp_dictoffset */
    (initproc)PyTree_init,       /* tp_init */
};

/* ========================================================================== */
/* -- Methods --------------------------------------------------------------- */
/* ========================================================================== */

/* kcluster */
static char kcluster__doc__[] =
"kcluster(data, nclusters=2, mask=None, weight=None,\n"
"         transpose=0, npass=1, method='a', dist='e',\n"
"         initialid=None) -> clusterid, error, nfound\n"
"\n"
"This function implements k-means clustering.\n"
"data     : nrows x ncolumns array containing the expression data\n"
"nclusters: number of clusters (the 'k' in k-means)\n"
"mask     : nrows x ncolumns array of integers, showing which data are\n"
"           missing. If mask[i][j]==0, then data[i][j] is missing.\n"
"weight   : the weights to be used when calculating distances\n"
"transpose: if equal to 0, genes (rows) are clustered;\n"
"           if equal to 1, microarrays (columns) are clustered.\n"
"npass    : number of times the k-means clustering algorithm is performed,\n"
"           each time with a different (random) initial condition.\n"
"method   : specifies how the center of a cluster is found:\n"
"           method=='a': arithmetic mean\n"
"           method=='m': median\n"
"dist     : specifies the distance function to be used:\n"
"           dist=='e': Euclidean distance\n"
"           dist=='b': City Block distance\n"
"           dist=='c': Pearson correlation\n"
"           dist=='a': absolute value of the correlation\n"
"           dist=='u': uncentered correlation\n"
"           dist=='x': absolute uncentered correlation\n"
"           dist=='s': Spearman's rank correlation\n"
"           dist=='k': Kendall's tau\n"
"initialid: the initial clustering from which the algorithm should start.\n"
"           If initialid is None, the routine carries out npass repetitions\n"
"           of the EM algorithm, each time starting from a different random\n"
"           initial clustering. If initialid is given, the routine carries\n"
"           out the EM algorithm only once, starting from the given initial\n"
"           clustering and without randomizing the order in which items are\n"
"           assigned to clusters (i.e., using the same order as in the data\n"
"           matrix). In that case, the k-means algorithm is fully\n"
"           deterministic.\n"
"\n"
"Return values:\n"
"clusterid: array containing the number of the cluster to which each\n"
"           gene/microarray was assigned in the best k-means clustering\n"
"           solution that was found in the npass runs;\n"
"error:     the within-cluster sum of distances for the returned k-means\n"
"           clustering solution;\n"
"nfound:    the number of times this solution was found.\n";

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
  double ERROR;
  int IFOUND;

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
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|iOOiiccO", kwlist,
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
  /* -- Reset None variables to NULL ------------------------------------- */
  if(MASK==Py_None) MASK = NULL;
  if(WEIGHT==Py_None) WEIGHT = NULL;
  if(INITIALID==Py_None) INITIALID = NULL;
  /* -- Check the method variable ---------------------------------------- */
  if (!strchr("am", METHOD))
  { sprintf(message, "method %c is unknown", METHOD);
    PyErr_SetString (PyExc_ValueError, buffer);
    return NULL;
  }
  /* -- Check the dist variable ------------------------------------------ */
  if (!strchr(known_distances, DIST))
  { sprintf(message, "dist %c is an unknown distance function", DIST);
    PyErr_SetString (PyExc_ValueError, buffer);
    return NULL;
  }
  /* -- Check the transpose variable ------------------------------------- */
  if (TRANSPOSE) TRANSPOSE = 1;
  /* -- Check the npass variable ----------------------------------------- */
  if (INITIALID) NPASS = 0;
  else if (NPASS <= 0)
  { strcpy(message, "npass should be a positive integer");
    PyErr_SetString (PyExc_ValueError, buffer);
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
    PyErr_SetString (PyExc_ValueError, buffer);
    free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    Py_DECREF((PyObject*) aCLUSTERID);
    return NULL;
  }
  if (nitems < NCLUSTERS)
  { strcpy(message, "More clusters than items to be clustered");
    PyErr_SetString (PyExc_ValueError, buffer);
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
      &ERROR, 
      &IFOUND);
  /* --------------------------------------------------------------------- */
  free_data(aDATA, data);
  free_mask(aMASK, mask, nrows);
  free_weight(aWEIGHT, weight);
  /* --------------------------------------------------------------------- */

  return Py_BuildValue("Ndi", aCLUSTERID, ERROR, IFOUND);
} 
/* end of wrapper for kcluster */

/* kmedoids */
static char kmedoids__doc__[] =
"kmedoids(distance, nclusters=2, npass=1,\n"
"         initialid=None) -> clusterid, error, nfound.\n"
"\n"
"This function implements k-medoids clustering.\n"
"distance:  The distance matrix between the elements. There are three ways\n"
"           in which you can pass a distance matrix:\n"
"           #1: a 2D Numerical Python array (in which only the left-lower\n"
"               part of the array will be accessed);\n"
"           #2: a 1D Numerical Python array containing the distances\n"
"               consecutively;\n" 
"           #3: a list of rows containing the lower-triangular part of the\n"
"               distance matrix.\n"
"           Examples are:\n"
"           >>> distance = array([[0.0, 1.1, 2.3],\n"
"                                 [1.1, 0.0, 4.5],\n"
"                                 [2.3, 4.5, 0.0]])\n"
"           (option #1)\n"
"           >>> distance = array([1.1, 2.3, 4.5])\n"
"           (option #2)\n"
"           >>> distance = [array([]),\n"
"                           array([1.1]),\n"
"                           array([2.3, 4.5])\n"
"                          ]\n"
"           (option #3)\n"
"           These three correspond to the same distance matrix.\n"
"nclusters: number of clusters (the 'k' in k-medoids)\n"
"npass    : the number of times the k-medoids clustering algorithm is\n"
"           performed, each time with a different (random) initial\n"
"           condition.\n"
"initialid: the initial clustering from which the algorithm should start.\n"
"           If initialid is not given, the routine carries out npass\n"
"           repetitions of the EM algorithm, each time starting from a\n"
"           different random initial clustering. If initialid is given,\n"
"           the routine carries out the EM algorithm only once, starting\n"
"           from the initial clustering specified by initialid and without\n"
"           randomizing the order in which items are assigned to clusters\n"
"           (i.e., using the same order as in the data matrix). In that\n"
"           case, the k-means algorithm is fully deterministic.\n"
"\n"
"Return values:\n"
"clusterid: array containing the number of the cluster to which each\n"
"           gene/microarray was assigned in the best k-means clustering\n"
"           solution that was found in the npass runs;\n"
"error:     the within-cluster sum of distances for the returned k-means\n"
"           clustering solution;\n"
"nfound:    the number of times this solution was found.\n";

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
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|iiO", kwlist,
                                  &DISTANCES,
                                  &NCLUSTERS,
                                  &NPASS,
                                  &INITIALID)) return NULL;
  /* Set the function name for error messages */
  strcpy (buffer, "kmedoids: ");
  message = strchr(buffer, '\0');
  /* -- Reset None variables to NULL ------------------------------------- */
  if (INITIALID==Py_None) INITIALID = NULL;
  /* -- Check the npass variable ----------------------------------------- */
  if (INITIALID) NPASS = 0;
  else if (NPASS < 0)
  { strcpy(message, "npass should be a positive integer");
    PyErr_SetString (PyExc_ValueError, buffer);
    return NULL;
  }
  /* -- Check the distance matrix ---------------------------------------- */
  distances = parse_distance(DISTANCES, &aDISTANCES, &nitems);
  if (!distances) return NULL;
  /* -- Create the clusterid output variable ----------------------------- */
  aCLUSTERID = parse_initialid(INITIALID, &NCLUSTERS, nitems);
  if (!aCLUSTERID)
  { free_distances(DISTANCES, aDISTANCES, distances, nitems);
    return NULL;
  }
  /* -- Check the nclusters variable ------------------------------------- */
  if (NCLUSTERS <= 0)
  { strcpy(buffer,"nclusters should be a positive integer");
    PyErr_SetString (PyExc_ValueError, buffer);
    free_distances(DISTANCES, aDISTANCES, distances, nitems);
    Py_DECREF((PyObject*) aCLUSTERID);
    return NULL;
  }
  if (nitems < NCLUSTERS)
  { strcpy(message, "More clusters than items to be clustered");
    PyErr_SetString (PyExc_ValueError, buffer);
    free_distances(DISTANCES, aDISTANCES, distances, nitems);
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
  free_distances(DISTANCES, aDISTANCES, distances, nitems);
  /* --------------------------------------------------------------------- */
  if(IFOUND==0) /* should not occur */
  { Py_DECREF((PyObject*) aCLUSTERID);
    strcpy(message, "Error in kmedoids input arguments");
    PyErr_SetString (PyExc_RuntimeError, buffer);
    return NULL;
  }
  if(IFOUND==-1)
  { Py_DECREF((PyObject*) aCLUSTERID);
    strcpy(message, "Memory allocation error in kmedoids");
    PyErr_SetString (PyExc_MemoryError, buffer);
    return NULL;
  }
  return Py_BuildValue("Ndi",aCLUSTERID, ERROR, IFOUND);
} 
/* end of wrapper for kmedoids */

/* treecluster */
static char treecluster__doc__[] =
"treecluster(data=None, mask=None, weight=None, transpose=0, dist='e',\n"
"            method='m', distancematrix=None) -> Tree object\n"
"\n"
"This function implements the pairwise single, complete, centroid, and\n"
"average linkage hierarchical clustering methods.\n"
"data     : nrows x ncolumns array containing the gene expression data.\n"
"mask     : nrows x ncolumns array of integers, showing which data are\n"
"           missing. If mask[i][j]==0, then data[i][j] is missing.\n"
"weight   : the weights to be used when calculating distances.\n"
"transpose: if equal to 0, genes (rows) are clustered;\n"
"           if equal to 1, microarrays (columns) are clustered.\n"
"dist     : specifies the distance function to be used:\n"
"           dist=='e': Euclidean distance\n"
"           dist=='b': City Block distance\n"
"           dist=='c': Pearson correlation\n"
"           dist=='a': absolute value of the correlation\n"
"           dist=='u': uncentered correlation\n"
"           dist=='x': absolute uncentered correlation\n"
"           dist=='s': Spearman's rank correlation\n"
"           dist=='k': Kendall's tau\n"
"method   : specifies which linkage method is used:\n"
"           method=='s': Single pairwise linkage\n"
"           method=='m': Complete (maximum) pairwise linkage (default)\n"
"           method=='c': Centroid linkage\n"
"           method=='a': Average pairwise linkage\n"
"distancematrix:  The distance matrix between the elements. There are three\n"
"           ways in which you can pass a distance matrix:\n"
"           #1: a 2D Numerical Python array (in which only the left-lower\n"
"               part of the array will be accessed);\n"
"           #2: a 1D Numerical Python array containing the distances\n"
"               consecutively;\n" 
"           #3: a list of rows containing the lower-triangular part of the\n"
"               distance matrix.\n"
"           Examples are:\n"
"           >>> distance = array([[0.0, 1.1, 2.3],\n"
"                                 [1.1, 0.0, 4.5],\n"
"                                 [2.3, 4.5, 0.0]])\n"
"           (option #1)\n"
"           >>> distance = array([1.1, 2.3, 4.5])\n"
"           (option #2)\n"
"           >>> distance = [array([]),\n"
"                           array([1.1]),\n"
"                           array([2.3, 4.5])\n"
"                          ]\n"
"           (option #3)\n"
"           These three correspond to the same distance matrix.\n"
"           PLEASE NOTE:\n"
"           As the treecluster routine may shuffle the values in the distance\n"
"           matrix as part of the clustering algorithm, be sure to save this\n"
"           array in a different variable before calling treecluster if you\n"
"           need it later.\n"
"\n"
"Either data or distancematrix should be None. If distancematrix==None,\n"
"the hierarchical clustering solution is calculated from the gene expression\n"
"data stored in the argument data. If data==None, the hierarchical clustering\n"
"solution is calculated from the distance matrix instead. Pairwise centroid-\n"
"linkage clustering can be calculated only from the gene expression data and\n"
"not from the distance matrix. Pairwise single-, maximum-, and average-\n"
"linkage clustering can be calculated from either the gene expression data or\n"
"from the distance matrix.\n"
"\n"
"Return value:\n"
"treecluster returns a Tree object describing the hierarchical clustering\n"
"result. See the description of the Tree class for more information.\n";

static PyObject*
py_treecluster (PyObject* self, PyObject* args, PyObject* keywords)
{ PyObject *DATA = NULL;
  PyObject *MASK = NULL;
  PyObject *WEIGHT = NULL;
  int TRANSPOSE = 0;
  char DIST = 'e';
  char METHOD = 'm';
  PyObject *DISTANCEMATRIX = NULL;
  PyTree* tree;
  Node* nodes;
  int nitems;

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = { "data",
                            "mask",
                            "weight",
                            "transpose",
                            "method",
                            "dist",
                            "distancematrix",
                             NULL };
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "|OOOiccO", kwlist,
                                  &DATA,
                                  &MASK,
                                  &WEIGHT,
                                  &TRANSPOSE,
                                  &METHOD,
                                  &DIST,
                                  &DISTANCEMATRIX)) return NULL;
  /* Set the function name for error messages */
  strcpy (buffer, "treecluster: ");
  message = strchr(buffer, '\0');

  /* -- Reset None variables to NULL ------------------------------------- */
  if(DATA==Py_None) DATA = NULL;
  if(MASK==Py_None) MASK = NULL;
  if(WEIGHT==Py_None) WEIGHT = NULL;
  if(DISTANCEMATRIX==Py_None) DISTANCEMATRIX = NULL;

  /* -- Check if we are using the data matrix or the distance matrix ----- */
  if (DATA!=NULL && DISTANCEMATRIX!=NULL)
  { strcpy(message, "Use either data or distancematrix, do not use both");
    PyErr_SetString(PyExc_ValueError, buffer);
    return NULL;
  }
  if (DATA==NULL && DISTANCEMATRIX==NULL)
  { strcpy(message, "Neither data nor distancematrix was given");
    PyErr_SetString(PyExc_ValueError, buffer);
    return NULL;
  }

  if (DISTANCEMATRIX==NULL) /* DATA contains gene expression data */
  { int nrows;
    int ncolumns;
    int ndata;
    PyArrayObject* aDATA = NULL;
    PyArrayObject* aMASK = NULL;
    PyArrayObject* aWEIGHT = NULL;
    double** data = NULL;
    int** mask = NULL;
    double* weight = NULL;

    /* -- Check the method variable ---------------------------------------- */
    if (!strchr("csma", METHOD))
    { strcpy(message, "keyword method should be 'c', 's', 'm', or 'a'");
      PyErr_SetString(PyExc_ValueError, buffer);
      return NULL;
    }
    /* -- Check the dist variable ------------------------------------------ */
    if (!strchr(known_distances, DIST))
    { sprintf(message, "unknown distance function specified (dist='%c')", DIST);
      PyErr_SetString(PyExc_ValueError, buffer);
      return NULL;
    }
    /* -- Check the data input array --------------------------------------- */
    data = parse_data(DATA, &aDATA);
    if (!data) return NULL;
    nrows = aDATA->dimensions[0];
    ncolumns = aDATA->dimensions[1];
    ndata = TRANSPOSE ? nrows : ncolumns;
    nitems = TRANSPOSE ? ncolumns : nrows;
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
    /* -- Call treecluster to perform hierarchical clustering -------------- */
    nodes = treecluster(nrows,
                        ncolumns,
                        data,
                        mask,
                        weight,
                        TRANSPOSE,
                        DIST,
                        METHOD,
                        NULL);
    /* --------------------------------------------------------------------- */
    free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    free_weight(aWEIGHT, weight);
  }
  else
  { double** distances = NULL;
    PyArrayObject* aDISTANCEMATRIX = NULL;
    if (!strchr("sma", METHOD))
    { strcpy(message,
        "argument method should be 's', 'm', or 'a' when specifying the distance matrix");
      PyErr_SetString (PyExc_ValueError, buffer);
      return NULL;
    }
    /* -- Check the distance matrix ---------------------------------------- */
    distances = parse_distance(DISTANCEMATRIX, &aDISTANCEMATRIX, &nitems);
    if (!distances) return NULL;
    /* --------------------------------------------------------------------- */
    nodes = treecluster(nitems,
                        nitems,
                        0,
                        0,
                        0,
                        TRANSPOSE,
                        DIST,
                        METHOD,
                        distances);
    /* --------------------------------------------------------------------- */
    free_distances(DISTANCEMATRIX, aDISTANCEMATRIX, distances, nitems);
  }

  /* -- Check if a memory allocation error occurred ---------------------- */
  if(!nodes)
  { PyErr_SetString (PyExc_MemoryError, "error occurred in treecluster");
    return NULL;
  }
  tree = (PyTree*) PyTreeType.tp_alloc(&PyTreeType, 0);
  if(!tree)
  { PyErr_SetString (PyExc_MemoryError, "error occurred in treecluster");
    free(nodes);
    return NULL;
  }
  tree->nodes = nodes;
  tree->n = nitems-1;
  return (PyObject*) tree;
} 
/* end of wrapper for treecluster */

/* somcluster */
static char somcluster__doc__[] =
"somcluster(data, mask=None, weight=None, transpose=0, nxgrid=2, nygrid=1,\n"
"           inittau=0.02, niter=1, dist='e') -> clusterid, celldata\n"
"\n"
"This function implements a self-organizing map on a rectangular grid.\n"
"data     : nrows x ncolumns array containing the gene expression data\n"
"mask     : nrows x ncolumns array of integers, showing which data are\n"
"           missing. If mask[i][j]==0, then data[i][j] is missing.\n"
"weight   : the weights to be used when calculating distances\n"
"transpose: if equal to 0, genes (rows) are clustered;\n"
"           if equal to 1, microarrays (columns) are clustered.\n"
"nxgrid   : the horizontal dimension of the rectangular SOM map\n"
"nygrid   : the vertical dimension of the rectangular SOM map\n"
"inittau  : the initial value of tau (the neighborbood function)\n"
"niter    : the number of iterations\n"
"dist     : specifies the distance function to be used:\n"
"           dist=='e': Euclidean distance\n"
"           dist=='b': City Block distance\n"
"           dist=='c': Pearson correlation\n"
"           dist=='a': absolute value of the correlation\n"
"           dist=='u': uncentered correlation\n"
"           dist=='x': absolute uncentered correlation\n"
"           dist=='s': Spearman's rank correlation\n"
"           dist=='k': Kendall's tau\n"
"\n"
"Return values:\n"
"clusterid: array with two columns, while the number of rows is equal to the\n"
"           number of genes or the number of microarrays depending on whether\n"
"           genes or microarrays are being clustered. Each row in the array\n"
"           contains the x and y coordinates of the cell in the rectangular\n"
"           SOM grid to which the gene or microarray was assigned.\n"
"celldata:  an array with dimensions (nxgrid, nygrid, number of microarrays)\n"
"           if genes are being clustered, or (nxgrid, nygrid, number of\n"
"           genes) if microarrays are being clustered. Each element [ix][iy]\n"
"           of this array is a 1D vector containing the gene expression data\n"
"           for the centroid of the cluster in the SOM grid cell with\n"
"           coordinates (ix, iy).\n";

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
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|OOiiidic", kwlist,
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
  /* -- Reset None variables to NULL ------------------------------------- */
  if(WEIGHT==Py_None) WEIGHT = NULL;
  if(MASK==Py_None) MASK = NULL;
  /* -- Check the nxgrid variable ---------------------------------------- */
  if (NXGRID < 1)
  { strcpy(message, "nxgrid should be a positive integer (default is 2)");
    PyErr_SetString (PyExc_ValueError, buffer);
    return NULL;
  }
  /* -- Check the nygrid variable ---------------------------------------- */
  if (NYGRID < 1)
  { strcpy(message, "nygrid should be a positive integer (default is 1)");
    PyErr_SetString (PyExc_ValueError, buffer);
    return NULL;
  }
  /* -- Check the niter variable ----------------------------------------- */
  if (NITER < 1)
  { strcpy(message, "number of iterations (niter) should be positive");
    PyErr_SetString (PyExc_ValueError, buffer);
    return NULL;
  }
  /* -- Check the dist variable ------------------------------------------ */
  if (!strchr(known_distances, DIST))
  { sprintf(message, "dist %c is an unknown distance function", DIST);
    PyErr_SetString (PyExc_ValueError, buffer);
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
  aCLUSTERID = (PyArrayObject*) PyArray_FromDims(2, shape, PyArray_INT);
  if (!aCLUSTERID)
  { strcpy(buffer, "somcluster: Could not create clusterid array");
    PyErr_SetString (PyExc_MemoryError, buffer);
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
"median(data) -> median value of the 1D array data\n"
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
      PyErr_SetString (PyExc_TypeError, buffer);
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
      PyErr_SetString (PyExc_ValueError, buffer);
      return NULL;
    }
  } 
  if (aDATA->nd != 1 && (aDATA->nd > 0 || aDATA->dimensions[0] != 1))
  { sprintf(buffer, "median: Argument has incorrect rank (%d expected 1).",
                    aDATA->nd);
    PyErr_SetString (PyExc_ValueError, buffer);
    Py_DECREF((PyObject*) aDATA);
    return NULL;
  }
  if (!PyArray_ISCONTIGUOUS(aDATA))
  { PyObject* av =
      PyArray_ContiguousFromObject((PyObject*) aDATA, aDATA->descr->type_num, 0, 0);
    Py_DECREF((PyObject*)aDATA);
    if(!av)
    { strcpy(buffer, "median: Failed making argument contiguous.");
      PyErr_SetString (PyExc_ValueError, buffer);
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
"mean(data) -> arithmetic mean of the 1D array data.\n";

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
      PyErr_SetString (PyExc_TypeError, buffer);
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
      PyErr_SetString (PyExc_ValueError, buffer);
      return NULL;
    }
  } 
  if (aDATA->nd != 1 && (aDATA->nd > 0 || aDATA->dimensions[0] != 1))
  { sprintf(buffer, "mean: Argument has incorrect rank (%d expected 1).",
                    aDATA->nd);
    PyErr_SetString (PyExc_ValueError, buffer);
    Py_DECREF((PyObject*) aDATA);
    return NULL;
  }
  if (!PyArray_ISCONTIGUOUS(aDATA))
  { PyObject* av =
      PyArray_ContiguousFromObject((PyObject*) aDATA, aDATA->descr->type_num, 0, 0);
    Py_DECREF((PyObject*)aDATA);
    if(!av)
    { strcpy(buffer, "mean: Failed making argument contiguous.");
      PyErr_SetString (PyExc_ValueError, buffer);
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
"clusterdistance(data, mask=None, weight=None, index1, index2, dist='e',\n"
"                method='a', transpose=0) -> the distance between the\n"
"                                            two clusters\n"
"\n"
"data     : nrows x ncolumns array containing the expression data\n"
"mask     : nrows x ncolumns array of integers, showing which data are\n"
"           missing. If mask[i][j]==0, then data[i][j] is missing.\n"
"weight   : the weights to be used when calculating distances\n"
"index1   : 1D array identifying which genes/microarrays belong to the first\n"
"           cluster. If the cluster contains only one gene, then index1 can\n"
"           also be written as a single integer.\n"
"index2   : 1D array identifying which genes/microarrays belong to the second\n"
"           cluster. If the cluster contains only one gene, then index2 can\n"
"           also be written as a single integer.\n"
"transpose: if equal to 0, genes (rows) are clustered;\n"
"           if equal to 1, microarrays (columns) are clustered.\n"
"dist     : specifies the distance function to be used:\n"
"           dist=='e': Euclidean distance\n"
"           dist=='b': City Block distance\n"
"           dist=='c': Pearson correlation\n"
"           dist=='a': absolute value of the correlation\n"
"           dist=='u': uncentered correlation\n"
"           dist=='x': absolute uncentered correlation\n"
"           dist=='s': Spearman's rank correlation\n"
"           dist=='k': Kendall's tau\n"
"method   : specifies how the distance between two clusters is defined:\n"
"           method=='a': the distance between the arithmetic means of the\n"
"                        two clusters\n"
"           method=='m': the distance between the medians of the two clusters\n"
"           method=='s': the smallest pairwise distance between members of\n"
"                        the two clusters\n"
"           method=='x': the largest pairwise distance between members of\n"
"                        the two clusters\n"
"           method=='v': average of the pairwise distances between members\n"
"                        of the clusters\n"
"transpose: if equal to 0, clusters of genes (rows) are considered;\n"
"           if equal to 1, clusters of microarrays (columns) are considered.\n";

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
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|OOOOcci", kwlist,
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

  /* -- Reset None variables to NULL ------------------------------------- */
  if (MASK==Py_None) MASK = NULL;
  if (WEIGHT==Py_None) WEIGHT = NULL;
  if (INDEX1==Py_None) INDEX1 = NULL;
  if (INDEX2==Py_None) INDEX2 = NULL;
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
  if (result < -0.5) /* Actually -1.0; avoiding roundoff errors */
  { PyErr_SetString(PyExc_IndexError, "index out of range");
    return NULL;
  }
  return PyFloat_FromDouble(result);
} 
/* end of wrapper for clusterdistance */

/* clustercentroids */
static char clustercentroids__doc__[] =
"clustercentroids(data, mask=None, transport=0, clusterid,\n"
"                 method='a') -> cdata, cmask\n"
"\n"
"The clustercentroids routine calculates the cluster centroids, given to\n"
"which cluster each element belongs. The centroid is defined as either the\n"
"mean or the median over all elements for each dimension.\n"

"data     : nrows x ncolumns array containing the expression data\n"
"mask     : nrows x ncolumns array of integers, showing which data are\n"
"           missing. If mask[i][j]==0, then data[i][j] is missing.\n"
"transpose: if equal to 0, gene (row) clusters are considered;\n"
"           if equal to 1, microarray (column) clusters are considered.\n"
"clusterid: array containing the cluster number for each gene or microarray.\n"
"           The cluster number should be non-negative.\n"
"method   : specifies whether the centroid is calculated from the arithmetic\n"
"           mean (method=='a', default) or the median (method=='m') over each\n"
"           dimension.\n"
"\n"
"Return values:\n"
"cdata    : 2D array containing the cluster centroids. If transpose==0, then\n"
"           the dimensions of cdata are nclusters x ncolumns. If\n"
"           transpose==1, then the dimensions of cdata are nrows x nclusters.\n"
"cmask    : 2D array of integers describing which elements in cdata, if any,\n"
"           are missing.\n";

static PyObject*
py_clustercentroids(PyObject* self, PyObject* args, PyObject* keywords)
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
  int ok;

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = { "data",
                            "mask",
                            "clusterid",
                            "method",
                            "transpose",
                             NULL };
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|OOci", kwlist,
                                  &DATA,
                                  &MASK,
                                  &CLUSTERID,
                                  &METHOD,
                                  &TRANSPOSE)) return NULL;
  /* Set the function name for error messages */
  strcpy (buffer, "clustercentroids: ");
  message = strchr(buffer, '\0');
  /* -- Reset None variables to NULL ------------------------------------- */
  if (MASK==Py_None) MASK = NULL;
  if (CLUSTERID==Py_None) CLUSTERID = NULL;
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
  { strcpy(message, "could not create centroids array");
    PyErr_SetString (PyExc_MemoryError, buffer);
    free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    free_clusterid(aCLUSTERID, clusterid);
    return NULL;
  }
  cdata = malloc(shape[0]*sizeof(double*));
  for (i=0; i<shape[0]; i++)
    cdata[i] = ((double*) (aCDATA->data)) + i*shape[1];
  /* -- Create the centroid mask output variable ------------------------- */
  aCMASK = (PyArrayObject*) PyArray_FromDims(2, shape, PyArray_INT);
  if (!aCMASK)
  { strcpy(message, "could not create centroids array");
    PyErr_SetString (PyExc_MemoryError, buffer);
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
  ok = getclustercentroids(nclusters,
                           nrows,
                           ncolumns,
                           data,
                           mask,
                           clusterid,
                           cdata,
                           cmask,
                           TRANSPOSE,
                           METHOD);
  /* --------------------------------------------------------------------- */
  free_data(aDATA, data);
  free_mask(aMASK, mask, nrows);
  free (cdata);
  free (cmask);
  free_clusterid(aCLUSTERID, clusterid);
  /* --------------------------------------------------------------------- */
  if (!ok)
  { strcpy(message, "allocation error in clustercentroids");
    PyErr_SetString (PyExc_MemoryError, buffer);
  }
  return Py_BuildValue("NN", PyArray_Return(aCDATA), PyArray_Return(aCMASK));
} 
/* end of wrapper for clustercentroids */

/* distancematrix */
static char distancematrix__doc__[] =
"distancematrix(data, mask=None, weight=None, transpose=0,\n"
"               dist='e') -> distance matrix as a list of arrays\n"
"\n"
"This function returns the distance matrix between gene expression data.\n"
"data     : nrows x ncolumns array containing the expression data\n"
"mask     : nrows x ncolumns array of integers, showing which data are\n"
"           missing. If mask[i][j]==0, then data[i][j] is missing.\n"
"weight   : the weights to be used when calculating distances.\n"
"transpose: if equal to 0, the distances between genes (rows) are calculated;\n"
"           if equal to 1, the distances beteeen microarrays (columns) are\n"
"           calculated.\n"
"dist     : specifies the distance function to be used:\n"
"           dist=='e': Euclidean distance\n"
"           dist=='b': City Block distance\n"
"           dist=='c': Pearson correlation\n"
"           dist=='a': absolute value of the correlation\n"
"           dist=='u': uncentered correlation\n"
"           dist=='x': absolute uncentered correlation\n"
"           dist=='s': Spearman's rank correlation\n"
"           dist=='k': Kendall's tau\n"
"\n"
"Return value:\n"
"The distance matrix is returned as a list of 1D arrays containing the\n"
"distance matrix between the gene expression data. The number of columns in\n"
"each row is equal to the row number. Hence, the first row has zero elements.\n"
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
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|OOic", kwlist,
                                  &DATA,
                                  &MASK,
                                  &WEIGHT,
                                  &TRANSPOSE,
                                  &DIST)) return NULL;
  /* Set the function name for error messages */
  strcpy (buffer, "distancematrix: ");
  message = strchr(buffer, '\0');
  /* -- Reset None variables to NULL ------------------------------------- */
  if (MASK==Py_None) MASK = NULL;
  if (WEIGHT==Py_None) WEIGHT = NULL;
  /* -- Check the dist variable ------------------------------------------ */
  if (!strchr(known_distances, DIST))
  { sprintf(message, "dist %c is an unknown distance function", DIST);
    PyErr_SetString (PyExc_ValueError, buffer);
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
        { strcpy(message, "could not create distance matrix");
          PyErr_SetString (PyExc_MemoryError, buffer);
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
  { strcpy(message, "Could not create distance matrix");
    PyErr_SetString (PyExc_MemoryError, buffer);
  }
  /* --------------------------------------------------------------------- */
  free_data(aDATA, data);
  free_mask(aMASK, mask, nrows);
  free_weight(aWEIGHT, weight);
  return result;
}


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
   {"clustercentroids", (PyCFunction) py_clustercentroids, METH_KEYWORDS, clustercentroids__doc__},
   {"distancematrix", (PyCFunction) py_distancematrix, METH_KEYWORDS, distancematrix__doc__},
   {NULL,          NULL, 0, NULL}/* sentinel */
};

/* ========================================================================== */
/* -- Initialization -------------------------------------------------------- */
/* ========================================================================== */

void initcluster(void)
{
  PyObject *m;

  import_array();

  PyNodeType.tp_new = PyType_GenericNew;
  PyTreeType.tp_new = PyType_GenericNew;
  if (PyType_Ready(&PyNodeType) < 0) return;
  if (PyType_Ready(&PyTreeType) < 0) return;

  m = Py_InitModule4("cluster",
                     methods,
                     "C Clustering Library",
                     NULL,
                     PYTHON_API_VERSION);
  if (m==NULL) return;

  Py_INCREF(&PyTreeType);
  Py_INCREF(&PyNodeType);
  PyModule_AddObject(m, "Tree", (PyObject*) &PyTreeType);
  PyModule_AddObject(m, "Node", (PyObject*) &PyNodeType);

  if (PyErr_Occurred()) Py_FatalError("can't initialize module cluster");
}
