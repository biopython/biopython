#include "Python.h"
#include "Numeric/arrayobject.h"
#include <stdio.h>
#include <string.h>
#include "cluster.h"
 
static PyObject *ErrorObject;
static char buffer[512];
static char* message = NULL;

static const char known_distances[] = "ebhcauxsk";

static void
set_pyfort_error (char* routine, char* var, char* problem) {
    sprintf(buffer, "%s, argument %s: %s", routine, var, problem);
    PyErr_SetString (ErrorObject, buffer);
}

static PyArrayObject* make_contiguous(char* rname, char* vname, PyArrayObject* ap)
{
/* return an owned ref to a contiguous version of ap */
    PyArrayObject* result;
    if (ap->flags & CONTIGUOUS) {
        Py_INCREF (ap);
        return ap;
    } else {
        result = (PyArrayObject *) PyArray_ContiguousFromObject(
			(PyObject*) ap, ap->descr->type_num, 0, 0);
        if(!result) set_pyfort_error(rname, vname, "Failed making object contiguous.");
        return result;
    }
}

static int do_size_check (char* rname, char* vname, PyArrayObject *av, int rank,  int extents[])
{
    int size1;
    int i;

    size1 = av->nd;
    
    if( size1 == rank) {
        for(i=0; i < rank; ++i) {
            /* no checking on last dimension of expected size 1 */
            if (i == size1-1) {
               if (extents[i] == 1) break;
            }
            if(av->dimensions[i] != extents[i]) 
            {
               sprintf(buffer, "%s, argument %s: Incorrect extent in dimension %d (%d expected %d)", rname, vname, i+1, av->dimensions[i], extents[i]);
               PyErr_SetString (ErrorObject, buffer);
               return 0;
            }
        } 
    } else {
        if (rank != 1 || 
            size1 > 0 ||
            extents[0] != 1) 
        {    
           sprintf(buffer, "%s, argument %s: Incorrect rank (%d expected %d)", rname, vname, size1, rank);
           PyErr_SetString (ErrorObject, buffer);
           return 0;
        }
    }
    return 1; /* size ok */
}

static PyArrayObject*
do_array_in (char* rname, char* vname, PyObject *v, 
    enum PyArray_TYPES python_array_type)
{
    PyArrayObject* av;
    PyArrayObject* t;

    if(!PyArray_Check (v)) {
        t = (PyArrayObject *) PyArray_ContiguousFromObject(v, PyArray_NOTYPE, 0, 0);
        if (!t) {
            set_pyfort_error(rname, vname, "Argument cannot be converted to needed array.");
            return (PyArrayObject*) 0;
        }
    } else {
        t = (PyArrayObject*) v;
        Py_INCREF((PyObject*) t);
    }
    if (t->descr->type_num != python_array_type) {
        av = (PyArrayObject*) PyArray_Cast (t, python_array_type);
        Py_DECREF((PyObject*) t);
        t = av;
        if (!t) {
            set_pyfort_error(rname, vname, "Argument cannot be cast to needed type.");
            return (PyArrayObject*) 0;
        }
    } 
    return t;
}


static PyArrayObject*
do_array_create (char* rname, char* vname, enum PyArray_TYPES python_array_type, 
    int rank, int extents[])
{
    PyArrayObject* av =
        (PyArrayObject*) PyArray_FromDims(rank, extents, python_array_type);
    if (!av) {
        set_pyfort_error(rname, vname, "Could not create array -- too big?");
        return (PyArrayObject*) 0;
    }
    return av;
}

static double**
parse_data(PyObject* object, PyArrayObject** array)
/* Takes the Python object from the argument list, and finds the microarray
 * data set. In case of an error, the array is DECREF'ed and set to NULL. */
{ int i, j;
  int nrows, ncols;
  double** data = NULL;
  if(!PyArray_Check (object)) /* Convert object to a 2D array of type double */
  { *array = (PyArrayObject*) PyArray_FromObject(object, PyArray_DOUBLE, 2, 2);
    if (*array==NULL)
    { strcpy (message, "data cannot be converted to needed array.");
      PyErr_SetString(ErrorObject, buffer);
      return NULL;
    }
  }
  else
  { /* User passed an array */
    *array = (PyArrayObject*) object;
    Py_INCREF((PyObject*) *array);
    if ((*array)->descr->type_num != PyArray_DOUBLE)
    { PyArrayObject* av = (PyArrayObject*) PyArray_Cast(*array, PyArray_DOUBLE);
      Py_DECREF((PyObject*) (*array));
      *array = av;
      if (!(*array))
      { strcpy (message, "data cannot be cast to needed type.");
        PyErr_SetString(ErrorObject, buffer);
        return NULL;
      }
    } 
    if ((*array)->nd != 2)
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
  data = (double**)malloc((size_t)(nrows*sizeof(double*)));
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
      data[i] = (double*)malloc((size_t)(ncols*sizeof(double)));
      for (j=0; j < ncols; j++, p+=colstride) data[i][j] = *((double*)p);
      p0 += rowstride;
    }
  }
  return data;
}

static void
free_data(PyArrayObject* array, double** data)
{ int i;
  const int nrows = array->dimensions[0];
  if(data[0]!=(double*)(array->data)) for (i=0; i<nrows; i++) free(data[i]);
  free (data);
  Py_DECREF((PyObject*) array);
  return;
}

static int**
parse_mask (PyObject* object, PyArrayObject** array, const int dimensions[])
{ int i, j;
  const int nrows = dimensions[0];
  const int ncolumns = dimensions[1];
  int** mask = (int**)malloc((size_t)nrows*sizeof(int*));
  /* -- Return the default mask if the user didn't specify anything. -- */
  if (object==NULL)
  { for (i=0; i<nrows; i++)
    { mask[i] = (int*)malloc((size_t)ncolumns*sizeof(int));
      for (j=0; j<ncolumns; j++) mask[i][j] = 1;
    }
    *array = NULL;
    return mask;
  }
  /* -- The user specified something. Let's see if it is an array ----- */
  if(!PyArray_Check (object))
  { *array = (PyArrayObject*) PyArray_FromObject(object, PyArray_LONG, 2, 2);
    if (!(*array))
    { strcpy (message, "mask cannot be converted to needed array");
      PyErr_SetString(ErrorObject, buffer);
      return NULL;
    }
  }
  else
  { *array = (PyArrayObject*) object;
    Py_INCREF((PyObject*) *array);
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
  /* -- Now that we have an array, we need to check its size ---------- */
  if((*array)->nd != 2)
  { sprintf(message, "mask has incorrect rank (%d expected 2)", (*array)->nd);
    PyErr_SetString (ErrorObject, buffer);
    Py_DECREF((PyObject*)*array);
    *array = NULL;
    return NULL;
  }
  if((*array)->dimensions[0] != nrows) 
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
      mask[i] = (int*)malloc((size_t)(ncolumns*sizeof(int)));
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

static double*
parse_weight (PyObject* object, PyArrayObject** array, const int ndata)
{ int i;
  double* weight = NULL;
  if (object==NULL)
  { weight = (double*)malloc((size_t)ndata*sizeof(double));
    for (i = 0; i < ndata; i++) weight[i] = 1.0;
    *array = NULL;
    return weight;
  }
  if(!PyArray_Check (object))
  { *array = (PyArrayObject*) PyArray_FromObject(object, PyArray_DOUBLE, 1, 1);
    if (!(*array))
    { strcpy (message, "weight cannot be converted to needed array.");
      PyErr_SetString(ErrorObject, buffer);
      return NULL;
    }
  }
  else
  { *array = (PyArrayObject*) object;
    Py_INCREF((PyObject*) *array);
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
  if((*array)->nd == 1)
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
  if ((*array)->flags & CONTIGUOUS) weight = (double*) ((*array)->data);
  else
  { const char* p = (char*) ((*array)->data);
    const int stride =  ((*array)->strides)[0];
    weight = (double*)malloc((size_t)ndata*sizeof(double));
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

static int
parse_initial(PyObject* object, PyArrayObject** array, const PyArrayObject* clusterid, int nclusters)
{ int i;
  int stride;
  const char* p;
  int* q;
  int* number;
  const int nitems = clusterid->dimensions[0];
  if(!PyArray_Check (object))
  { *array = (PyArrayObject*) PyArray_FromObject(object, PyArray_LONG,1,1);
    if (!(*array))
    { strcpy (message, "initialid cannot be converted to needed array.");
      PyErr_SetString(ErrorObject, buffer);
      return 0;
    }
  }
  else
  { *array = (PyArrayObject*) object;
    Py_INCREF((PyObject*) (*array));
  }
  if ((*array)->descr->type_num != PyArray_LONG)
  { PyArrayObject* av = (PyArrayObject*) PyArray_Cast(*array, PyArray_LONG);
    Py_DECREF((PyObject*) (*array));
    *array = av;
    if (!(*array))
    { strcpy (message, "initialid cannot be cast to needed type.");
      PyErr_SetString(ErrorObject, buffer);
      return 0;
    }
  } 
  if((*array)->nd == 1)
  { /* no checking on last dimension of expected size 1 */
    if (nitems!=1 && nitems!=(*array)->dimensions[0]) 
    { sprintf(message,
              "initialid has incorrect extent (%d expected %d)",
              (*array)->dimensions[0], nitems);
      PyErr_SetString (ErrorObject, buffer);
      Py_DECREF((PyObject*) (*array));
      return 0;
    }
  }
  else
  { if ((*array)->nd > 0 || nitems != 1)
    { sprintf(message,
             "initialid has incorrect rank (%d expected 1)",
              (*array)->nd);
      PyErr_SetString (ErrorObject, buffer);
      Py_DECREF((PyObject*) (*array));
      return 0;
    }
  }
  stride = (*array)->strides[0];
  p = (const char*) ((*array)->data);
  q = (int*) (clusterid->data);
  number = (int*)malloc((size_t)nclusters*sizeof(int));
  for (i = 0; i < nclusters; i++) number[i] = 0;
  for (i = 0; i < nitems; i++, p+=stride, q++)
  { if ((*q) < 0 || (*q) >= nclusters)
    { strcpy(message, "initialid contains an invalid cluster number");
      PyErr_SetString (ErrorObject, buffer);
      return 0;
    }
    number[*q]++;
  }
  for (i = 0; i < nclusters; i++) if(number[i]==0) break;
  free(number);
   Py_DECREF((PyObject*) (*array));
  if (i < nclusters)
  { sprintf (message, "argument initialid: Cluster %d is empty", i);
    PyErr_SetString (ErrorObject, buffer);
    return 0;
  }
  return 1;
}

static double**
parse_distance(PyObject* object, PyArrayObject** array, int* n)
/* Takes the Python object from the argument list, and finds the distance
 * matrix. In case of an error, the array is DECREF'ed and set to NULL. */
{ int i, j;
  double** distance = NULL;
  if(!PyArray_Check (object)) /* Convert object to a 2D array of type double */
  { *array = (PyArrayObject*) PyArray_FromObject(object, PyArray_DOUBLE, 1, 2);
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
    Py_INCREF((PyObject*) (*array));
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
      return NULL;
    }
    distance = (double**)malloc((size_t)(*n)*sizeof(double*));
    distance[0] = NULL;
    if (stride==sizeof(double)) /* Data are contiguous */
      for (i=1; i < *n; p+=(i*stride), i++) distance[i] = (double*)p;
    else /* We need to create contiguous rows */
    { for (i=1; i < *n; i++)
      { distance[i] = (double*)malloc((size_t)i*sizeof(double));
        for (j=0; j < i; j++, p+=stride) distance[i][j] = *((double*)p);
      }
    }
  }
  else if ((*array)->nd == 2)
  { const char* p = (char*) ((*array)->data);
    *n = (*array)->dimensions[0];
    distance = (double**)malloc((size_t)(*n)*sizeof(double*));
    distance[0] = NULL;
    if ((*array)->strides[1]==sizeof(double)) /* Each row is contiguous */
    { const int stride =  (*array)->strides[0];
      for (i=0; i < *n; i++, p+=stride) distance[i] = (double*)p;
    }
    else /* We need to create contiguous rows */
    { const int stride =  (*array)->strides[1];
      for (i=0; i < *n; i++)
      { distance[i] = (double*)malloc((size_t)i*sizeof(double));
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
    const int n = (int) ((1+sqrt(1+4*m))/2);
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


/* Methods */

/* kcluster */
static char cluster_kcluster__doc__[] =
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
"method=='a': arithmic mean\n"
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
cluster_kcluster (PyObject* self, PyObject* args, PyObject* keywords) {

  int NCLUSTERS = 2;
  int nrows, ncolumns;
  int nitems;
  int ndata;
  PyObject* DATA = NULL;
  PyArrayObject* aDATA = NULL;
  double** data = NULL;
  PyObject* MASK = NULL;
  PyArrayObject* aMASK = NULL;
  int** mask = NULL;
  PyObject* WEIGHT = 0;
  PyArrayObject* aWEIGHT = NULL;
  double* weight = NULL;
  int TRANSPOSE = 0;
  int NPASS = 1;
  char METHOD = 'a';
  char DIST = 'e';
  PyObject* INITIAL = NULL;
  PyArrayObject* aINITIAL = NULL;
  PyArrayObject* aCLUSTERID = NULL;
  PyArrayObject* aCDATA = NULL;
  int shape[2];
  double ERROR;
  int IFOUND;
  double** cdata;
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
                                  &INITIAL)) return NULL;
  /* Set the function name for error messages */
  strcpy (buffer, "kcluster: ");
  message = strchr(buffer, '\0');
  /* -- Check the nclusters variable ------------------------------------- */
  if (NCLUSTERS < 1)
  { strcpy(message, "nclusters should be positive");
    PyErr_SetString (ErrorObject, buffer);
    return NULL;
  }
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
  if (NPASS < 0)
  { strcpy(message, "npass should be 0 or more");
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
  /* -- Check the number of clusters ------------------------------------- */
  ndata = TRANSPOSE ? nrows : ncolumns;
  nitems = TRANSPOSE ? ncolumns : nrows;
  if (nitems < NCLUSTERS)
  { strcpy(message, "More clusters than items to be clustered");
    PyErr_SetString (ErrorObject, buffer);
    free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    return NULL;
  }
  /* -- Check the weight input ------------------------------------------- */
  weight = parse_weight(WEIGHT, &aWEIGHT, ndata);
  if (!weight)
  { free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    return NULL;
  }
  /* -- Create the clusterid output variable ----------------------------- */
  aCLUSTERID = (PyArrayObject*) PyArray_FromDims(1, &nitems, PyArray_LONG);
  if (!aCLUSTERID)
  { strcpy(message, "Could not create clusterid array -- too big?");
    PyErr_SetString (ErrorObject, buffer);
    free_data(aDATA, data);
    free_mask(aMASK, mask, nrows);
    free_weight(aWEIGHT, weight);
    return NULL;
  }
  /* -- Check if the user specified an initial clustering ---------------- */
  if (INITIAL)
  { int result = parse_initial(INITIAL, &aINITIAL, aCLUSTERID, NCLUSTERS);
    if (!result)
    { free_data(aDATA, data);
      free_mask(aMASK, mask, nrows);
      free_weight(aWEIGHT, weight);
      Py_DECREF((PyObject*) aCLUSTERID);
      return NULL;
    }
    NPASS = 0;
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
    Py_DECREF((PyObject*) aINITIAL);
  }
  cdata = (double**)malloc((size_t)shape[0]*sizeof(double*));
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

  return Py_BuildValue("OOdl",aCLUSTERID, aCDATA, ERROR, IFOUND);
} 
/* end of wrapper for kcluster */

void CALL kmedoids (int nclusters, int nelements, double** distance,
  int npass, int clusterid[], double* error, int* ifound);

/* kmedoids */
static char cluster_kmedoids__doc__[] =
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
cluster_kmedoids (PyObject* self, PyObject* args, PyObject* keywords) {

  int NCLUSTERS = 2;
  int nitems;
  PyObject* DISTANCES = NULL;
  PyArrayObject* aDISTANCES = NULL;
  double** distances = NULL;
  PyObject* INITIAL = NULL;
  PyArrayObject* aINITIAL = NULL;
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
                                  &INITIAL)) return NULL;
  /* Set the function name for error messages */
  strcpy (buffer, "kmedoids: ");
  message = strchr(buffer, '\0');
  /* -- Check the npass variable ----------------------------------------- */
  if (NPASS < 0)
  { strcpy (message, "npass should be 0 or more");
    PyErr_SetString (ErrorObject, buffer);
    return NULL;
  }
  /* -- Check the nclusters variable ------------------------------------- */
  if (NCLUSTERS <= 0)
  { strcpy(buffer,"nclusters should be a positive integer");
    PyErr_SetString (ErrorObject, buffer);
    return NULL;
  }
  /* -- Check the distance matrix ---------------------------------------- */
  distances = parse_distance(DISTANCES, &aDISTANCES, &nitems);
  if (!distances) return NULL;
  /* -- Check the number of clusters ------------------------------------- */
  if (nitems < NCLUSTERS)
  { strcpy(message, "More clusters than items to be clustered");
    PyErr_SetString (ErrorObject, buffer);
    free_distances(aDISTANCES, distances);
    return NULL;
  }
  /* -- Create the clusterid output variable ----------------------------- */
  aCLUSTERID = (PyArrayObject*) PyArray_FromDims(1, &nitems, PyArray_LONG);
  if (!aCLUSTERID)
  { strcpy(message, "could not create clusterid array -- too big?");
    PyErr_SetString (ErrorObject, buffer);
    free_distances(aDISTANCES, distances);
    return NULL;
  }
  /* -- Check if the user specified an initial clustering ---------------- */
  if (INITIAL)
  { int result = parse_initial(INITIAL, &aINITIAL, aCLUSTERID, NCLUSTERS);
    if (!result)
    { free_distances(aDISTANCES, distances);
      Py_DECREF((PyObject*) aCLUSTERID);
      return NULL;
    }
    NPASS = 0;
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
  return Py_BuildValue("Odl",aCLUSTERID, ERROR, IFOUND);
} 
/* end of wrapper for kmedoids */

/* treecluster */
static char cluster_treecluster__doc__[] =
"hierarchical clustering\n"
"result, linkdist = treecluster(data,mask,weight,applyscale,transpose,dist,\n"
"                               method,distances)\n"
"This function implements the pairwise centroid-, single-, maximum-, and\n"
"average-linkage clustering algorithm.\n"
"The nrows x ncolumns array data contains the gene expression data.\n"
"The array mask declares missing data. If mask[i][j]==0, then data[i][j]\n"
"is missing.\n"
"The array weight contains the weights to be used for the distance calculation.\n"
"If the integer applyscale is nonzero, then the distances in linkdist are\n"
"scaled such that all distances are between zero and two (as in case of the\n"
"Pearson distance).\n"
"The integer transpose defines if rows (genes) or columns (microarrays) are\n"
"clustered. If transpose==0, then genes are clustered. If transpose==1,\n"
"microarrays are clustered.\n"
"The character dist defines the distance function to be used:\n"
"dist=='e': Euclidean distance\n"
"dist=='b': City-block distance\n"
"dist=='h': Harmonically summed Euclidean distance\n"
"dist=='c': correlation\n"
"dist=='a': absolute value of the correlation\n"
"dist=='u': uncentered correlation\n"
"dist=='x': absolute uncentered correlation\n"
"dist=='s': Spearman's rank correlation\n"
"dist=='k': Kendall's tau\n"
"For other values of dist, the default (Euclidean distance) is used.\n"
"The character method defines which hierarchical clustering method is used:\n"
"method=='s': Single-linkage\n"
"method=='m': Maximum- or complete-linkage\n"
"method=='a': Average-linkage\n"
"method=='c': Centroid-linkage\n"
"The integer distances denotes if the data array contains the original gene\n"
"expression data, or the distance matrix calculated from those data. If\n"
"distances==1, then data is interpreted as the distance matrix, and the\n"
"arguments mask, weight, transpose, and dist are ignored.\n"
"\n"
"Return values:\n"
"result is an (nobject x 2) array describing the hierarchical clustering\n"
"  result. Each row in the array represents one node, with the two columns\n"
"  representing the two objects or nodes that are being joined. Objects are\n"
"  numbered 0 through (nobjects-1), while nodes are numbered -1 through\n"
"  -(nobjects-1).\n"
"linkdist is a vector with (nobjects-1) elements containing the distances\n"
"between the two subnodes that are joined at each node.\n"
"\n";

static PyObject*
cluster_treecluster (PyObject* unused, PyObject* args) {

    PyObject *pyfort_result;
    int NROWS;
    int NCOLUMNS;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    int shape[2];
    PyObject* MASK;
    PyArrayObject* aMASK;
    PyObject* WEIGHT;
    PyArrayObject* aWEIGHT;
    int eWEIGHT;
    int APPLYSCALE;
    int TRANSPOSE;
    char DIST;
    char METHOD;
    int DISTANCES;
    PyArrayObject* aRESULT;
    PyObject* rRESULT;
    int eRESULT[2];
    PyArrayObject* aLINKDIST;
    PyObject* rLINKDIST;
    int eLINKDIST[1];
    int ii;
    aDATA = (PyArrayObject*) 0;
    aMASK = (PyArrayObject*) 0;
    aWEIGHT = (PyArrayObject*) 0;
    aRESULT = (PyArrayObject*) 0;
    aLINKDIST = (PyArrayObject*) 0;

    if(!PyArg_ParseTuple(args, "OOOllccl", &DATA, &MASK, &WEIGHT, &APPLYSCALE, &TRANSPOSE, &DIST, &METHOD, &DISTANCES)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("treecluster", "DATA", DATA, PyArray_DOUBLE))) goto err;
    NROWS = aDATA->dimensions[0];
    NCOLUMNS = aDATA->dimensions[1];
    shape[0] = NROWS;
    shape[1] = NCOLUMNS;
    eRESULT[0] = ((TRANSPOSE==1) ? NCOLUMNS : NROWS) - 1;
    eRESULT[1] = 2;
    eLINKDIST[0] = ((TRANSPOSE==1) ? NCOLUMNS : NROWS) - 1;
    if (!do_size_check ("treecluster", "DATA", aDATA, 2, shape)) goto err;
    pyarray_value = aDATA;
    aDATA = make_contiguous ("treecluster", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!(aRESULT = do_array_create ("treecluster", "RESULT", PyArray_LONG, 2, eRESULT))) goto err;
    if (!(aLINKDIST = do_array_create ("treecluster", "LINKDIST", PyArray_DOUBLE, 1, eLINKDIST))) goto err;
    if (DISTANCES==0) /* aDATA contains gene expression data */
    { double* paDATA = 0;
      double** ppaDATA = 0;
      int* paMASK = 0;
      int** ppaMASK = 0;
      if (!(aMASK = do_array_in ("treecluster", "MASK", MASK, PyArray_LONG))) goto err;
      if (!(aWEIGHT = do_array_in ("treecluster", "WEIGHT", WEIGHT, PyArray_DOUBLE))) goto err;
      eWEIGHT = TRANSPOSE ? NROWS : NCOLUMNS;
      if (!do_size_check ("treecluster", "MASK", aMASK, 2, shape)) goto err;
      pyarray_value = aMASK;
      aMASK = make_contiguous ("treecluster", "MASK", pyarray_value);
      Py_DECREF(pyarray_value);
      if(!aMASK) goto err;
      if (!do_size_check ("treecluster", "WEIGHT", aWEIGHT, 1, &eWEIGHT)) goto err;
      pyarray_value = aWEIGHT;
      aWEIGHT = make_contiguous ("treecluster", "WEIGHT", pyarray_value);
      Py_DECREF(pyarray_value);
      if(!aWEIGHT) goto err;
      ppaDATA = (double**)malloc((size_t)NROWS*sizeof(double*));
      paDATA = (double*) (aDATA->data);
      for (ii=0; ii<NROWS; ii++) ppaDATA[ii]=&(paDATA[ii*NCOLUMNS]);
      ppaMASK = (int**)malloc((size_t)NROWS*sizeof(int*));
      paMASK = (int*) (aMASK->data);
      for (ii=0; ii<NROWS; ii++) ppaMASK[ii]=&(paMASK[ii*NCOLUMNS]);
      treecluster(NROWS, 
          NCOLUMNS, 
          ppaDATA, 
          ppaMASK, 
          (double*) (aWEIGHT->data), 
          APPLYSCALE, 
          TRANSPOSE, 
          DIST, 
          METHOD, 
          (int(*)[2]) (aRESULT->data), 
          (double*) (aLINKDIST->data), 0);
      free (ppaDATA);
      free (ppaMASK);
    }
    else
    { int jj;
      double** distmatrix = 0;
      double* paDATA = 0;
      if(NROWS!=NCOLUMNS)
      { set_pyfort_error ("treecluster", "DATA", "matrix is not square");
        goto err;
      }
      distmatrix = (double**)malloc((size_t)NROWS*sizeof(double*));
      paDATA = (double*) (aDATA->data);
      for(ii=1; ii<NROWS; ii++)
      { distmatrix[ii] = (double*)malloc((size_t)ii*sizeof(double));
        for(jj=0; jj<ii; jj++)
          distmatrix[ii][jj] = paDATA[ii*NCOLUMNS+jj];
      }
      treecluster(NROWS, 
          NCOLUMNS, 
          0, 
          0, 
          0, 
          APPLYSCALE, 
          TRANSPOSE, 
          DIST, 
          METHOD, 
          (int(*)[2]) (aRESULT->data), 
          (double*) (aLINKDIST->data),
          distmatrix);
      for(ii=1; ii<NROWS; ii++) free(distmatrix[ii]);
      free(distmatrix);
    }
    rRESULT = PyArray_Return(aRESULT);
    rLINKDIST = PyArray_Return(aLINKDIST);
    Py_XDECREF((PyObject*) aDATA);
    Py_XDECREF((PyObject*) aMASK);
    Py_XDECREF((PyObject*) aWEIGHT);

    pyfort_result = Py_BuildValue("OO",rRESULT, rLINKDIST);

    Py_XDECREF(rRESULT);
    Py_XDECREF(rLINKDIST);
    return pyfort_result;
err:
    Py_XDECREF((PyObject*) aDATA);
    Py_XDECREF((PyObject*) aMASK);
    Py_XDECREF((PyObject*) aWEIGHT);
    Py_XDECREF((PyObject*) aRESULT);
    Py_XDECREF((PyObject*) aLINKDIST);
    return NULL;
} 
/* end of wrapper for treecluster */

/* somcluster */
static char cluster_somcluster__doc__[] =
"self-organizing map\n"
"somcluster (data,mask,weight,transpose,nxgrid,nygrid,inittau,niter,dist)\n"
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
"The initial value of tau (the neighborbood function)\n"
"The number of iterations is given by niter.\n"
"The character dist defines the distance function to be used:\n"
"dist=='e': Euclidean distance\n"
"dist=='b': City-block distance\n"
"dist=='h': Harmonically summed Euclidean distance\n"
"dist=='c': correlation\n"
"dist=='a': absolute value of the correlation\n"
"dist=='u': uncentered correlation\n"
"dist=='x': absolute uncentered correlation\n"
"dist=='s': Spearman's rank correlation\n"
"dist=='k': Kendall's tau\n"
"For other values of dist, the default (Euclidean distance) is used.\n"
"\n";

static PyObject*
cluster_somcluster (PyObject* unused, PyObject* args) {

    PyObject *pyfort_result;
    int NROWS;
    int NCOLUMNS;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    int shape[2];
    PyObject* MASK;
    PyArrayObject* aMASK;
    PyObject* WEIGHT;
    PyArrayObject* aWEIGHT;
    int eWEIGHT;
    int TRANSPOSE;
    int NXGRID;
    int NYGRID;
    double INITTAU;
    int NITER;
    char DIST;
    PyArrayObject* aCELLDATA;
    PyObject* rCELLDATA;
    int eCELLDATA[3];
    PyArrayObject* aCLUSTERID;
    PyObject* rCLUSTERID;
    int eCLUSTERID[2];
    int ii;
    double* paDATA;
    double** ppaDATA;
    int* paMASK;
    int** ppaMASK;
    double* paCELLDATA;
    double** ppaCELLDATA;
    double*** pppaCELLDATA;
    aDATA = (PyArrayObject*) 0;
    aMASK = (PyArrayObject*) 0;
    aWEIGHT = (PyArrayObject*) 0;
    aCELLDATA = (PyArrayObject*) 0;
    aCLUSTERID = (PyArrayObject*) 0;

    if(!PyArg_ParseTuple(args, "OOOllldlc", &DATA, &MASK, &WEIGHT, &TRANSPOSE, &NXGRID, &NYGRID, &INITTAU, &NITER, &DIST)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("somcluster", "DATA", DATA, PyArray_DOUBLE))) goto err;
    if (!(aMASK = do_array_in ("somcluster", "MASK", MASK, PyArray_LONG))) goto err;
    if (!(aWEIGHT = do_array_in ("somcluster", "WEIGHT", WEIGHT, PyArray_DOUBLE))) goto err;
    NROWS = aDATA->dimensions[0];
    NCOLUMNS = aDATA->dimensions[1];
    shape[0] = NROWS;
    shape[1] = NCOLUMNS;
    eWEIGHT = TRANSPOSE ? NROWS : NCOLUMNS;
    eCELLDATA[0] = NXGRID;
    eCELLDATA[1] = NYGRID;
    eCELLDATA[2] = TRANSPOSE ? NROWS : NCOLUMNS;
    eCLUSTERID[0] = TRANSPOSE ? NCOLUMNS : NROWS;
    eCLUSTERID[1] = 2;
    if (!do_size_check ("somcluster", "DATA", aDATA, 2, shape)) goto err;
    pyarray_value = aDATA;
    aDATA = make_contiguous ("somcluster", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!do_size_check ("somcluster", "MASK", aMASK, 2, shape)) goto err;
    pyarray_value = aMASK;
    aMASK = make_contiguous ("somcluster", "MASK", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aMASK) goto err;
    if (!do_size_check ("somcluster", "WEIGHT", aWEIGHT, 1, &eWEIGHT)) goto err;
    pyarray_value = aWEIGHT;
    aWEIGHT = make_contiguous ("somcluster", "WEIGHT", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aWEIGHT) goto err;
    if (!(aCELLDATA = do_array_create ("somcluster", "CELLDATA", PyArray_DOUBLE, 3, eCELLDATA))) goto err;
    if (!(aCLUSTERID = do_array_create ("somcluster", "CLUSTERID", PyArray_LONG, 2, eCLUSTERID))) goto err;
    ppaDATA = (double**)malloc((size_t)NROWS*sizeof(double*));
    ppaMASK = (int**)malloc((size_t)NROWS*sizeof(int*));
    ppaCELLDATA = (double**)malloc((size_t)NXGRID*NYGRID*sizeof(double*));
    pppaCELLDATA = (double***)malloc((size_t)NXGRID*sizeof(double**));
    paDATA = (double*) (aDATA->data);
    paMASK = (int*) (aMASK->data);
    paCELLDATA = (double*) (aCELLDATA->data);
    for (ii=0; ii<NROWS; ii++) ppaDATA[ii]=&(paDATA[ii*NCOLUMNS]);
    for (ii=0; ii<NROWS; ii++) ppaMASK[ii]=&(paMASK[ii*NCOLUMNS]);
    for (ii=0; ii<NXGRID*NYGRID; ii++) ppaCELLDATA[ii]=&(paCELLDATA[ii*((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NROWS)]);
    for (ii=0; ii<NXGRID; ii++) pppaCELLDATA[ii]=&(ppaCELLDATA[ii*NYGRID]);
    somcluster(NROWS, 
        NCOLUMNS, 
        ppaDATA, 
        ppaMASK, 
        (double*) (aWEIGHT->data), 
        TRANSPOSE, 
        NXGRID, 
        NYGRID, 
        INITTAU, 
        NITER, 
        DIST, 
        pppaCELLDATA, 
        (int(*)[2]) (aCLUSTERID->data));
    rCELLDATA = PyArray_Return(aCELLDATA);
    rCLUSTERID = PyArray_Return(aCLUSTERID);
    Py_XDECREF((PyObject*) aDATA);
    Py_XDECREF((PyObject*) aMASK);
    Py_XDECREF((PyObject*) aWEIGHT);
    free (ppaDATA);
    free (ppaMASK);
    free (ppaCELLDATA);
    free (pppaCELLDATA);

    pyfort_result = Py_BuildValue("OO", rCLUSTERID, rCELLDATA);

    Py_XDECREF(rCELLDATA);
    Py_XDECREF(rCLUSTERID);
    return pyfort_result;
err:
    Py_XDECREF((PyObject*) aDATA);
    Py_XDECREF((PyObject*) aMASK);
    Py_XDECREF((PyObject*) aWEIGHT);
    Py_XDECREF((PyObject*) aCELLDATA);
    Py_XDECREF((PyObject*) aCLUSTERID);
    return NULL;
} 
/* end of wrapper for somcluster */

/* median */
static char cluster_median__doc__[] =
"median (data)\n"
"This function returns the median of the 1D array data.\n"
"Note: data will be partially ordered upon return.\n";

static PyObject*
cluster_median (PyObject* unused, PyObject* args) {

    double result;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    aDATA = (PyArrayObject*) 0;

    if(!PyArg_ParseTuple(args, "O", &DATA)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("median", "DATA", DATA, PyArray_DOUBLE))) goto err;
    if (aDATA->nd != 1 && (aDATA->nd > 0 || aDATA->dimensions[0] != 1))
    { sprintf(buffer, "median, argument data: Incorrect rank (%d expected 1)",
                      aDATA->nd);
      PyErr_SetString (ErrorObject, buffer);
      goto err;
    }
    pyarray_value = aDATA;
    aDATA = make_contiguous ("median", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    result = median(aDATA->dimensions[0], (double*) (aDATA->data));
    Py_XDECREF((PyObject*) aDATA);

    return PyFloat_FromDouble(result);
err:
    Py_XDECREF((PyObject*) aDATA);
    return NULL;
} 
/* end of wrapper for median */

/* mean */
static char cluster_mean__doc__[] =
"mean (data)\n"
"This function returns the mean of the 1D array data.\n";

static PyObject*
cluster_mean (PyObject* unused, PyObject* args) {

    double result;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    aDATA = (PyArrayObject*) 0;

    if(!PyArg_ParseTuple(args, "O", &DATA)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("mean", "DATA", DATA, PyArray_DOUBLE))) goto err;
    
    if (aDATA->nd != 1 && (aDATA->nd > 0 || aDATA->dimensions[0] != 1))
    { sprintf(buffer, "mean, argument data: Incorrect rank (%d expected 1)",
                      aDATA->nd);
      PyErr_SetString (ErrorObject, buffer);
      goto err;
    }
    pyarray_value = aDATA;
    aDATA = make_contiguous ("mean", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    result = mean(aDATA->dimensions[0], (double*) (aDATA->data));
    Py_XDECREF((PyObject*) aDATA);

    return PyFloat_FromDouble(result);
err:
    Py_XDECREF((PyObject*) aDATA);
    return NULL;
} 
/* end of wrapper for mean */

/* clusterdistance */
static char cluster_clusterdistance__doc__[] =
"The distance between two clusters\n"
"The nrows x ncolumns array data contains the gene expression data.\n"
"The array mask shows which data are missing. If mask[i][j]==0, then\n"
"data[i][j] is missing.\n"
"The array weight contains the weights to be used when calculating distances.\n"
"If transpose==0, then genes are clustered. If transpose==1, microarrays are\n"
"clustered.\n"
"The vector index1 identifies which genes/microarrays belong to the first\n"
"cluster.\n"
"The vector index2 identifies which genes/microarrays belong to the second\n"
"cluster.\n"
"The character dist defines the distance function to be used:\n"
"dist=='e': Euclidean distance\n"
"dist=='b': City-block distance\n"
"dist=='h': Harmonically summed Euclidean distance\n"
"dist=='c': correlation\n"
"dist=='a': absolute value of the correlation\n"
"dist=='u': uncentered correlation\n"
"dist=='x': absolute uncentered correlation\n"
"dist=='s': Spearman's rank correlation\n"
"dist=='k': Kendall's tau\n"
"For other values of dist, the default (Euclidean distance) is used.\n"
"The character method defines how the distance between two clusters is defined:\n"
"method=='a': the distance between the arithmic means of the two clusters\n"
"method=='m': the distance between the medians of the two clusters\n"
"method=='s': the smallest pairwise distance between members of the two clusters\n"
"method=='x': the largest pairwise distance between members of the two clusters\n"
"method=='v': average of the pairwise distances between members of the clusters\n";

static PyObject*
cluster_clusterdistance (PyObject* unused, PyObject* args) {

    double result;
    int NROWS;
    int NCOLUMNS;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    int shape[2];
    PyObject* MASK;
    PyArrayObject* aMASK;
    PyObject* WEIGHT;
    PyArrayObject* aWEIGHT;
    int eWEIGHT;
    int N1;
    int N2;
    PyObject* INDEX1;
    PyArrayObject* aINDEX1;
    int eINDEX1[1];
    PyObject* INDEX2;
    PyArrayObject* aINDEX2;
    int eINDEX2[1];
    char DIST;
    char METHOD;
    int TRANSPOSE;
    int ii;
    double* paDATA;
    double** ppaDATA;
    int* paMASK;
    int** ppaMASK;
    aDATA = (PyArrayObject*) 0;
    aMASK = (PyArrayObject*) 0;
    aWEIGHT = (PyArrayObject*) 0;
    aINDEX1 = (PyArrayObject*) 0;
    aINDEX2 = (PyArrayObject*) 0;

    if(!PyArg_ParseTuple(args, "OOOOOccl", &DATA, &MASK, &WEIGHT, &INDEX1, &INDEX2, &DIST, &METHOD, &TRANSPOSE)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("clusterdistance", "DATA", DATA, PyArray_DOUBLE))) goto err;
    if (!(aMASK = do_array_in ("clusterdistance", "MASK", MASK, PyArray_LONG))) goto err;
    if (!(aWEIGHT = do_array_in ("clusterdistance", "WEIGHT", WEIGHT, PyArray_DOUBLE))) goto err;
    if (!(aINDEX1 = do_array_in ("clusterdistance", "INDEX1", INDEX1, PyArray_LONG))) goto err;
    if (!(aINDEX2 = do_array_in ("clusterdistance", "INDEX2", INDEX2, PyArray_LONG))) goto err;
    NROWS = aDATA->dimensions[0];
    NCOLUMNS = aDATA->dimensions[1];
    N1 = aINDEX1->dimensions[0];
    N2 = aINDEX2->dimensions[0];
    shape[0] = NROWS;
    shape[1] = NCOLUMNS;
    eWEIGHT = TRANSPOSE ? NROWS : NCOLUMNS;
    eINDEX1[0] = N1;
    eINDEX2[0] = N2;
    if (!do_size_check ("clusterdistance", "DATA", aDATA, 2, shape)) goto err;
    pyarray_value = aDATA;
    aDATA = make_contiguous ("clusterdistance", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!do_size_check ("clusterdistance", "MASK", aMASK, 2, shape)) goto err;
    pyarray_value = aMASK;
    aMASK = make_contiguous ("clusterdistance", "MASK", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aMASK) goto err;
    if (!do_size_check ("clusterdistance", "WEIGHT", aWEIGHT, 1, &eWEIGHT)) goto err;
    pyarray_value = aWEIGHT;
    aWEIGHT = make_contiguous ("clusterdistance", "WEIGHT", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aWEIGHT) goto err;
    if (!do_size_check ("clusterdistance", "INDEX1", aINDEX1, 1, eINDEX1)) goto err;
    pyarray_value = aINDEX1;
    aINDEX1 = make_contiguous ("clusterdistance", "INDEX1", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aINDEX1) goto err;
    if (!do_size_check ("clusterdistance", "INDEX2", aINDEX2, 1, eINDEX2)) goto err;
    pyarray_value = aINDEX2;
    aINDEX2 = make_contiguous ("clusterdistance", "INDEX2", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aINDEX2) goto err;
    ppaDATA = (double**)malloc((size_t)NROWS*sizeof(double*));
    ppaMASK = (int**)malloc((size_t)NROWS*sizeof(int*));
    paDATA = (double*) (aDATA->data);
    paMASK = (int*) (aMASK->data);
    for (ii=0; ii<NROWS; ii++) ppaDATA[ii]=&(paDATA[ii*NCOLUMNS]);
    for (ii=0; ii<NROWS; ii++) ppaMASK[ii]=&(paMASK[ii*NCOLUMNS]);
    result = clusterdistance(NROWS, 
        NCOLUMNS, 
        ppaDATA, 
        ppaMASK, 
        (double*) (aWEIGHT->data), 
        N1, 
        N2, 
        (int*) (aINDEX1->data), 
        (int*) (aINDEX2->data), 
        DIST, 
        METHOD, 
        TRANSPOSE);
    Py_XDECREF((PyObject*) aDATA);
    Py_XDECREF((PyObject*) aMASK);
    Py_XDECREF((PyObject*) aWEIGHT);
    Py_XDECREF((PyObject*) aINDEX1);
    Py_XDECREF((PyObject*) aINDEX2);
    free (ppaDATA);
    free (ppaMASK);

    return PyFloat_FromDouble(result);
err:
    Py_XDECREF((PyObject*) aDATA);
    Py_XDECREF((PyObject*) aMASK);
    Py_XDECREF((PyObject*) aWEIGHT);
    Py_XDECREF((PyObject*) aINDEX1);
    Py_XDECREF((PyObject*) aINDEX2);
    return NULL;
} 
/* end of wrapper for clusterdistance */

/* getclustermean */
static char cluster_getclustermean__doc__[] =
"cluster centroid, using mean\n"
"getclustermean(nclusters,data,mask,clusterid,transpose)\n"
"The getclustermean routine calculates the cluster centroids, given to which\n"
"cluster each element belongs. The centroid is defined as the mean over all\n"
"elements for each dimension.\n"
"The number of clusters is nclusters.\n"
"The nrows x ncolumns array data contains the gene expression data.\n"
"The array mask declares missing data. If mask[i][j]==0, then data[i][j]\n"
"is missing.\n"
"The integer transpose defines if rows (genes) or columns (microarrays) are\n"
"clustered. If transpose==0, then genes are clustered. If transpose==1,\n"
"microarrays are clustered.\n"
"The array clusterid contains the cluster number for each gene or microarray.\n"
"Upon return, the array cdata contains the cluster centroids. If\n"
"transpose==0, then the dimensions of cdata are nclusters x ncolumns. If\n"
"transpose==1, then the dimensions of cdata are nrows x nclusters.\n"
"The array cmask describes which elements in cdata, if any, are missing.\n"
"\n";

static PyObject*
cluster_getclustermean (PyObject* unused, PyObject* args) {

    PyObject *pyfort_result;
    int NCLUSTERS;
    int NROWS;
    int NCOLUMNS;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    int shape[2];
    PyObject* MASK;
    PyArrayObject* aMASK;
    PyObject* CLUSTERID;
    PyArrayObject* aCLUSTERID;
    int eCLUSTERID;
    PyArrayObject* aCDATA;
    PyObject* rCDATA;
    int eCshape[2];
    PyArrayObject* aCMASK;
    PyObject* rCMASK;
    int TRANSPOSE;
    int ii;
    double* paDATA;
    double** ppaDATA;
    int* paMASK;
    int** ppaMASK;
    double* paCDATA;
    double** ppaCDATA;
    int* paCMASK;
    int** ppaCMASK;
    aDATA = (PyArrayObject*) 0;
    aMASK = (PyArrayObject*) 0;
    aCLUSTERID = (PyArrayObject*) 0;
    aCDATA = (PyArrayObject*) 0;
    aCMASK = (PyArrayObject*) 0;

    if(!PyArg_ParseTuple(args, "lOOOl", &NCLUSTERS, &DATA, &MASK, &CLUSTERID, &TRANSPOSE)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("getclustermean", "DATA", DATA, PyArray_DOUBLE))) goto err;
    if (!(aMASK = do_array_in ("getclustermean", "MASK", MASK, PyArray_LONG))) goto err;
    if (!(aCLUSTERID = do_array_in ("getclustermean", "CLUSTERID", CLUSTERID, PyArray_LONG))) goto err;
    NROWS = aDATA->dimensions[0];
    NCOLUMNS = aDATA->dimensions[1];
    shape[0] = NROWS;
    shape[1] = NCOLUMNS;
    eCLUSTERID = TRANSPOSE ? NCOLUMNS : NROWS;
    eCshape[0] = TRANSPOSE ? NROWS : NCLUSTERS;
    eCshape[1] = TRANSPOSE ? NCLUSTERS : NCOLUMNS;
    if (!do_size_check ("getclustermean", "DATA", aDATA, 2, shape)) goto err;
    pyarray_value = aDATA;
    aDATA = make_contiguous ("getclustermean", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!do_size_check ("getclustermean", "MASK", aMASK, 2, shape)) goto err;
    pyarray_value = aMASK;
    aMASK = make_contiguous ("getclustermean", "MASK", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aMASK) goto err;
    if (!do_size_check ("getclustermean", "CLUSTERID", aCLUSTERID, 1, &eCLUSTERID)) goto err;
    pyarray_value = aCLUSTERID;
    aCLUSTERID = make_contiguous ("getclustermean", "CLUSTERID", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aCLUSTERID) goto err;
    if (!(aCDATA = do_array_create ("getclustermean", "CDATA", PyArray_DOUBLE, 2, eCshape))) goto err;
    if (!(aCMASK = do_array_create ("getclustermean", "CMASK", PyArray_LONG, 2, eCshape))) goto err;
    ppaDATA = (double**)malloc((size_t)NROWS*sizeof(double*));
    ppaMASK = (int**)malloc((size_t)NROWS*sizeof(int*));
    ppaCDATA = (double**)malloc((size_t)((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS)*sizeof(double*));
    ppaCMASK = (int**)malloc((size_t)((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS)*sizeof(int*));
    paDATA = (double*) (aDATA->data);
    paMASK = (int*) (aMASK->data);
    paCDATA = (double*) (aCDATA->data);
    paCMASK = (int*) (aCMASK->data);
    for (ii=0; ii<NROWS; ii++) ppaDATA[ii]=&(paDATA[ii*NCOLUMNS]);
    for (ii=0; ii<NROWS; ii++) ppaMASK[ii]=&(paMASK[ii*NCOLUMNS]);
    for (ii=0; ii<((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS); ii++) ppaCDATA[ii]=&(paCDATA[ii*((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS)]);
    for (ii=0; ii<((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS); ii++) ppaCMASK[ii]=&(paCMASK[ii*((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS)]);
    getclustermean(NCLUSTERS, 
        NROWS, 
        NCOLUMNS, 
        ppaDATA, 
        ppaMASK, 
        (int*) (aCLUSTERID->data), 
        ppaCDATA, 
        ppaCMASK, 
        TRANSPOSE);
    rCDATA = PyArray_Return(aCDATA);
    rCMASK = PyArray_Return(aCMASK);
    Py_XDECREF((PyObject*) aDATA);
    Py_XDECREF((PyObject*) aMASK);
    Py_XDECREF((PyObject*) aCLUSTERID);
    free (ppaDATA);
    free (ppaMASK);
    free (ppaCDATA);
    free (ppaCMASK);

    pyfort_result = Py_BuildValue("OO",rCDATA, rCMASK);

    Py_XDECREF(rCDATA);
    Py_XDECREF(rCMASK);
    return pyfort_result;
err:
    Py_XDECREF((PyObject*) aDATA);
    Py_XDECREF((PyObject*) aMASK);
    Py_XDECREF((PyObject*) aCLUSTERID);
    Py_XDECREF((PyObject*) aCDATA);
    Py_XDECREF((PyObject*) aCMASK);
    return NULL;
} 
/* end of wrapper for getclustermean */

/* getclustermedian */
static char cluster_getclustermedian__doc__[] =
"cluster centroid, using median\n"
"getclustermedian(nclusters,data,mask,clusterid,transpose)\n"
"The getclustermedian routine calculates the cluster centroids, given to which\n"
"cluster each element belongs. The centroid is defined as the median over all\n"
"elements for each dimension.\n"
"The number of clusters is nclusters.\n"
"The nrows x ncolumns array data contains the gene expression data.\n"
"The array mask declares missing data. If mask[i][j]==0, then data[i][j]\n"
"is missing.\n"
"The integer transpose defines if rows (genes) or columns (microarrays) are\n"
"clustered. If transpose==0, then genes are clustered. If transpose==1,\n"
"microarrays are clustered.\n"
"The array clusterid contains the cluster number for each gene or microarray.\n"
"Upon return, the array cdata contains the cluster centroids. If\n"
"transpose==0, then the dimensions of cdata are nclusters x ncolumns. If\n"
"transpose==1, then the dimensions of cdata are nrows x nclusters.\n"
"The array cmask describes which elements in cdata, if any, are missing.\n"
"\n";

static PyObject*
cluster_getclustermedian (PyObject* unused, PyObject* args) {

    PyObject *pyfort_result;
    int NCLUSTERS;
    int NROWS;
    int NCOLUMNS;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    int shape[2];
    PyObject* MASK;
    PyArrayObject* aMASK;
    PyObject* CLUSTERID;
    PyArrayObject* aCLUSTERID;
    int eCLUSTERID;
    PyArrayObject* aCDATA;
    PyObject* rCDATA;
    int eCshape[2];
    PyArrayObject* aCMASK;
    PyObject* rCMASK;
    int TRANSPOSE;
    int ii;
    double* paDATA;
    double** ppaDATA;
    int* paMASK;
    int** ppaMASK;
    double* paCDATA;
    double** ppaCDATA;
    int* paCMASK;
    int** ppaCMASK;
    aDATA = (PyArrayObject*) 0;
    aMASK = (PyArrayObject*) 0;
    aCLUSTERID = (PyArrayObject*) 0;
    aCDATA = (PyArrayObject*) 0;
    aCMASK = (PyArrayObject*) 0;

    if(!PyArg_ParseTuple(args, "lOOOl", &NCLUSTERS, &DATA, &MASK, &CLUSTERID, &TRANSPOSE)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("getclustermedian", "DATA", DATA, PyArray_DOUBLE))) goto err;
    if (!(aMASK = do_array_in ("getclustermedian", "MASK", MASK, PyArray_LONG))) goto err;
    if (!(aCLUSTERID = do_array_in ("getclustermedian", "CLUSTERID", CLUSTERID, PyArray_LONG))) goto err;
    NROWS = aDATA->dimensions[0];
    NCOLUMNS = aDATA->dimensions[1];
    shape[0] = NROWS;
    shape[1] = NCOLUMNS;
    eCLUSTERID = TRANSPOSE ? NCOLUMNS : NROWS;
    eCshape[0] = TRANSPOSE ? NROWS : NCLUSTERS;
    eCshape[1] = TRANSPOSE ? NCLUSTERS : NCOLUMNS;
    if (!do_size_check ("getclustermedian", "DATA", aDATA, 2, shape)) goto err;
    pyarray_value = aDATA;
    aDATA = make_contiguous ("getclustermedian", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!do_size_check ("getclustermedian", "MASK", aMASK, 2, shape)) goto err;
    pyarray_value = aMASK;
    aMASK = make_contiguous ("getclustermedian", "MASK", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aMASK) goto err;
    if (!do_size_check ("getclustermedian", "CLUSTERID", aCLUSTERID, 1, &eCLUSTERID)) goto err;
    pyarray_value = aCLUSTERID;
    aCLUSTERID = make_contiguous ("getclustermedian", "CLUSTERID", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aCLUSTERID) goto err;
    if (!(aCDATA = do_array_create ("getclustermedian", "CDATA", PyArray_DOUBLE, 2, eCshape))) goto err;
    if (!(aCMASK = do_array_create ("getclustermedian", "CMASK", PyArray_LONG, 2, eCshape))) goto err;
    ppaDATA = (double**)malloc((size_t)NROWS*sizeof(double*));
    ppaMASK = (int**)malloc((size_t)NROWS*sizeof(int*));
    ppaCDATA = (double**)malloc((size_t)((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS)*sizeof(double*));
    ppaCMASK = (int**)malloc((size_t)((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS)*sizeof(int*));
    paDATA = (double*) (aDATA->data);
    paMASK = (int*) (aMASK->data);
    paCDATA = (double*) (aCDATA->data);
    paCMASK = (int*) (aCMASK->data);
    for (ii=0; ii<NROWS; ii++) ppaDATA[ii]=&(paDATA[ii*NCOLUMNS]);
    for (ii=0; ii<NROWS; ii++) ppaMASK[ii]=&(paMASK[ii*NCOLUMNS]);
    for (ii=0; ii<((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS); ii++) ppaCDATA[ii]=&(paCDATA[ii*((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS)]);
    for (ii=0; ii<((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS); ii++) ppaCMASK[ii]=&(paCMASK[ii*((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS)]);
    getclustermedian(NCLUSTERS, 
        NROWS, 
        NCOLUMNS, 
        ppaDATA, 
        ppaMASK, 
        (int*) (aCLUSTERID->data), 
        ppaCDATA, 
        ppaCMASK, 
        TRANSPOSE);
    rCDATA = PyArray_Return(aCDATA);
    rCMASK = PyArray_Return(aCMASK);
    Py_XDECREF((PyObject*) aDATA);
    Py_XDECREF((PyObject*) aMASK);
    Py_XDECREF((PyObject*) aCLUSTERID);
    free (ppaDATA);
    free (ppaMASK);
    free (ppaCDATA);
    free (ppaCMASK);

    pyfort_result = Py_BuildValue("OO",rCDATA, rCMASK);

    Py_XDECREF(rCDATA);
    Py_XDECREF(rCMASK);
    return pyfort_result;
err:
    Py_XDECREF((PyObject*) aDATA);
    Py_XDECREF((PyObject*) aMASK);
    Py_XDECREF((PyObject*) aCLUSTERID);
    Py_XDECREF((PyObject*) aCDATA);
    Py_XDECREF((PyObject*) aCMASK);
    return NULL;
} 
/* end of wrapper for getclustermedian */

/* distancematrix */
static char cluster_distancematrix__doc__[] =
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
cluster_distancematrix (PyObject* self, PyObject* args, PyObject* keywords)
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
    for (i = 0; i < nelements; i++)
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
  { strcpy(message, "Could not create distance matrix -- too big?");
    PyErr_SetString (ErrorObject, buffer);
  }
  /* --------------------------------------------------------------------- */
  free_data(aDATA, data);
  free_mask(aMASK, mask, nrows);
  free_weight(aWEIGHT, weight);
  return Py_BuildValue("O",result);
}

/* cuttree */
static char cluster_cuttree__doc__[] =
"clusterid = cuttree(tree, nclusters)\n"
"Given a hierarchical clustering result tree, the routine cuttree divides\n"
"the elements in the tree into clusters. The number of clusters is equal to\n"
"nclusters.\n";

static PyObject*
cluster_cuttree (PyObject* unused, PyObject* args) {

    PyObject *pyfort_result;
    int NELEMENTS;
    PyArrayObject* pyarray_value;
    PyObject* TREE;
    PyArrayObject* aTREE = (PyArrayObject*) 0;
    int NCLUSTERS = 1;
    PyArrayObject* aCLUSTERID = (PyArrayObject*) 0;
    PyObject* rCLUSTERID;

    if(!PyArg_ParseTuple(args, "O|l", &TREE, &NCLUSTERS)) {
        return NULL;
    }
    aTREE = do_array_in ("cuttree", "TREE", TREE, PyArray_LONG);
    if (!aTREE) goto err;
    NELEMENTS = aTREE->dimensions[0] + 1;
    if(aTREE->nd != 2) {
       sprintf(buffer, "cuttree, argument tree: Incorrect rank (%d expected 2)",
                       aTREE->nd);
       PyErr_SetString (ErrorObject, buffer);
       goto err;
    }
    pyarray_value = aTREE;
    aTREE = make_contiguous ("cuttree", "TREE", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aTREE) goto err;
    if (!(aCLUSTERID = do_array_create ("cuttree", "CLUSTERID", PyArray_LONG, 1, &NELEMENTS))) goto err;
    cuttree(NELEMENTS,
        (int(*)[2]) (aTREE->data),
        NCLUSTERS,
        (int*) (aCLUSTERID->data));
    if (((int*)(aCLUSTERID->data))[0]==-1) /* indicates an error in tree */
    {  PyErr_SetString (ErrorObject,
                        "cuttree, argument tree: incompatible input");
       goto err;
    }
    
    rCLUSTERID = PyArray_Return(aCLUSTERID);
    Py_XDECREF((PyObject*) aTREE);

    pyfort_result = Py_BuildValue("O",rCLUSTERID);

    Py_XDECREF(rCLUSTERID);
    return pyfort_result;
err:
    Py_XDECREF((PyObject*) aTREE);
    Py_XDECREF((PyObject*) aCLUSTERID);
    return NULL;
}
/* end of wrapper for cuttree */


static struct PyMethodDef cluster_methods[] = {
   {"kcluster", (PyCFunction) cluster_kcluster, METH_KEYWORDS, cluster_kcluster__doc__},
   {"kmedoids", (PyCFunction) cluster_kmedoids, METH_KEYWORDS, cluster_kmedoids__doc__},
   {"treecluster", (PyCFunction) cluster_treecluster, METH_VARARGS, cluster_treecluster__doc__},
   {"somcluster", (PyCFunction) cluster_somcluster, METH_VARARGS, cluster_somcluster__doc__},
   {"median", (PyCFunction) cluster_median, METH_VARARGS, cluster_median__doc__},
   {"mean", (PyCFunction) cluster_mean, METH_VARARGS, cluster_mean__doc__},
   {"clusterdistance", (PyCFunction) cluster_clusterdistance, METH_VARARGS, cluster_clusterdistance__doc__},
   {"getclustermean", (PyCFunction) cluster_getclustermean, METH_VARARGS, cluster_getclustermean__doc__},
   {"getclustermedian", (PyCFunction) cluster_getclustermedian, METH_VARARGS, cluster_getclustermedian__doc__},
   {"distancematrix", (PyCFunction) cluster_distancematrix, METH_KEYWORDS, cluster_distancematrix__doc__},
   {"cuttree", (PyCFunction) cluster_cuttree, METH_VARARGS, cluster_cuttree__doc__},
   {NULL,          NULL, 0, NULL}/* sentinel */
};

static char cluster_module_documentation[] =
"C interface module cluster";

void initcluster(void)
{
        PyObject *m, *d;
 
        import_array ();
        m = Py_InitModule4("cluster", cluster_methods,
                cluster_module_documentation,
                (PyObject*)NULL,PYTHON_API_VERSION);
 
        d = PyModule_GetDict(m);
        ErrorObject = PyString_FromString("cluster.error");
        PyDict_SetItemString(d, "error", ErrorObject);

        if (PyErr_Occurred()) {
            Py_FatalError("can't initialize module cluster");
        }
}
