#include "Python.h"
#include "Numeric/arrayobject.h"
 
static PyObject *ErrorObject;

static void
set_pyfort_error (char* routine, char* var, char* problem) {
    char buf[512];
    sprintf(buf, "%s, argument %s: %s", routine, var, problem);
    PyErr_SetString (ErrorObject, buf);
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
    char buf[512];

    size1 = av->nd;
    
    if( size1 == rank) {
        for(i=0; i < rank; ++i) {
            /* no checking on last dimension of expected size 1 */
            if (i == size1-1) {
               if (extents[i] == 1) break;
            }
            if(av->dimensions[i] != extents[i]) 
            {
               sprintf(buf, "Incorrect extent in dimension %d (%d expected %d)",
                       i+1, av->dimensions[i], extents[i]);
               set_pyfort_error(rname, vname, buf);
               return 0;
            }
        } 
    } else {
        if (rank != 1 || 
            size1 > 0 ||
            extents[0] != 1) 
        {    
           sprintf(buf, "Incorrect rank (%d expected %d)",
                       size1, rank);
           set_pyfort_error(rname, vname, buf);
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
/* Methods */

/* kcluster */
static char cluster_kcluster__doc__[] =
"k-means clustering\n"
"kcluster(nclusters,data,mask,weight,transpose,npass,method,dist)\n"
"This function implements the k-means clustering algorithm.\n"
"The number of clusters is given by ncluster.\n"
"The array data is a nrows x ncolumns array containing the gene \n"
"expression data.\n"
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
"For other values of method, the arithmic mean is used.\n"
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

extern void kcluster(long NCLUSTERS, long NROWS, long NCOLUMNS, double** DATA, long** MASK, double* WEIGHT, long TRANSPOSE, long NPASS, char METHOD, char DIST, long* CLUSTERID, double** CDATA, double* ERROR, long* IFOUND);

static PyObject*
cluster_kcluster (PyObject* unused, PyObject* args) {

    PyObject *pyfort_result;
    long NCLUSTERS;
    long NROWS;
    long NCOLUMNS;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    int eDATA[2];
    PyObject* MASK;
    PyArrayObject* aMASK;
    int eMASK[2];
    PyObject* WEIGHT;
    PyArrayObject* aWEIGHT;
    int eWEIGHT[1];
    long TRANSPOSE;
    long NPASS;
    char METHOD;
    char DIST;
    PyArrayObject* aCLUSTERID;
    PyObject* rCLUSTERID;
    int eCLUSTERID[1];
    PyArrayObject* aCDATA;
    PyObject* rCDATA;
    int eCDATA[2];
    double ERROR;
    long IFOUND;
    int ii;
    double* paDATA;
    double** ppaDATA;
    long* paMASK;
    long** ppaMASK;
    double* paCDATA;
    double** ppaCDATA;
    aDATA = (PyArrayObject*) 0;
    aMASK = (PyArrayObject*) 0;
    aWEIGHT = (PyArrayObject*) 0;
    aCLUSTERID = (PyArrayObject*) 0;
    aCDATA = (PyArrayObject*) 0;

    if(!PyArg_ParseTuple(args, "lOOOllcc", &NCLUSTERS, &DATA, &MASK, &WEIGHT, &TRANSPOSE, &NPASS, &METHOD, &DIST)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("kcluster", "DATA", DATA, PyArray_DOUBLE))) goto err;
    if (!(aMASK = do_array_in ("kcluster", "MASK", MASK, PyArray_LONG))) goto err;
    if (!(aWEIGHT = do_array_in ("kcluster", "WEIGHT", WEIGHT, PyArray_DOUBLE))) goto err;
    NROWS = aDATA->dimensions[0];
    NCOLUMNS = aDATA->dimensions[1];
    eDATA[0] = NROWS;
    eDATA[1] = NCOLUMNS;
    eMASK[0] = NROWS;
    eMASK[1] = NCOLUMNS;
    eWEIGHT[0] = TRANSPOSE ? NROWS : NCOLUMNS;
    eCLUSTERID[0] = TRANSPOSE ? NCOLUMNS : NROWS;
    eCDATA[0] = TRANSPOSE ? NROWS : NCLUSTERS;
    eCDATA[1] = TRANSPOSE ? NCLUSTERS : NCOLUMNS;
    if (!do_size_check ("kcluster", "DATA", aDATA, 2, eDATA)) goto err;
    pyarray_value = aDATA;
    aDATA = make_contiguous ("kcluster", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!do_size_check ("kcluster", "MASK", aMASK, 2, eMASK)) goto err;
    pyarray_value = aMASK;
    aMASK = make_contiguous ("kcluster", "MASK", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aMASK) goto err;
    if (!do_size_check ("kcluster", "WEIGHT", aWEIGHT, 1, eWEIGHT)) goto err;
    pyarray_value = aWEIGHT;
    aWEIGHT = make_contiguous ("kcluster", "WEIGHT", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aWEIGHT) goto err;
    if (!(aCLUSTERID = do_array_create ("kcluster", "CLUSTERID", PyArray_LONG, 1, eCLUSTERID))) goto err;
    if (!(aCDATA = do_array_create ("kcluster", "CDATA", PyArray_DOUBLE, 2, eCDATA))) goto err;
    ppaDATA = (double**)malloc((size_t)NROWS*sizeof(double*));
    ppaMASK = (long**)malloc((size_t)NROWS*sizeof(long*));
    ppaCDATA = (double**)malloc((size_t)(TRANSPOSE ? NROWS : NCLUSTERS)*sizeof(double*));
    paDATA = (double*) (aDATA->data);
    paMASK = (long*) (aMASK->data);
    paCDATA = (double*) (aCDATA->data);
    for (ii=0; ii<NROWS; ii++) ppaDATA[ii]=&(paDATA[ii*NCOLUMNS]);
    for (ii=0; ii<NROWS; ii++) ppaMASK[ii]=&(paMASK[ii*NCOLUMNS]);
    for (ii=0; ii<(TRANSPOSE ? NROWS : NCLUSTERS); ii++) ppaCDATA[ii]=&(paCDATA[ii*(TRANSPOSE ? NCLUSTERS : NCOLUMNS)]);
    kcluster(NCLUSTERS, 
        NROWS, 
        NCOLUMNS, 
        ppaDATA, 
        ppaMASK, 
        (double*) (aWEIGHT->data), 
        TRANSPOSE, 
        NPASS, 
        METHOD, 
        DIST, 
        (long*) (aCLUSTERID->data), 
        ppaCDATA, 
        &ERROR, 
        &IFOUND);
    rCLUSTERID = PyArray_Return(aCLUSTERID);
    rCDATA = PyArray_Return(aCDATA);
    Py_XDECREF((PyObject*) aDATA);
    Py_XDECREF((PyObject*) aMASK);
    Py_XDECREF((PyObject*) aWEIGHT);
    free (ppaDATA);
    free (ppaMASK);
    free (ppaCDATA);

    pyfort_result = Py_BuildValue("OOdl",rCLUSTERID, rCDATA, ERROR, IFOUND);

    Py_XDECREF(rCLUSTERID);
    Py_XDECREF(rCDATA);
    return pyfort_result;
err:
    Py_XDECREF((PyObject*) aDATA);
    Py_XDECREF((PyObject*) aMASK);
    Py_XDECREF((PyObject*) aWEIGHT);
    Py_XDECREF((PyObject*) aCLUSTERID);
    Py_XDECREF((PyObject*) aCDATA);
    return NULL;
} 
/* end of wrapper for kcluster */

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

extern void treecluster(long NROWS, long NCOLUMNS, double** DATA, long** MASK, double* WEIGHT, long APPLYSCALE, long TRANSPOSE, char DIST, char METHOD, long** RESULT, double* LINKDIST, double** DISTMATRIX);

static PyObject*
cluster_treecluster (PyObject* unused, PyObject* args) {

    PyObject *pyfort_result;
    long NROWS;
    long NCOLUMNS;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    int eDATA[2];
    PyObject* MASK;
    PyArrayObject* aMASK;
    int eMASK[2];
    PyObject* WEIGHT;
    PyArrayObject* aWEIGHT;
    int eWEIGHT[1];
    long APPLYSCALE;
    long TRANSPOSE;
    char DIST;
    char METHOD;
    long DISTANCES;
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
    eDATA[0] = NROWS;
    eDATA[1] = NCOLUMNS;
    eRESULT[0] = ((TRANSPOSE==1) ? NCOLUMNS : NROWS) - 1;
    eRESULT[1] = 2;
    eLINKDIST[0] = ((TRANSPOSE==1) ? NCOLUMNS : NROWS) - 1;
    if (!do_size_check ("treecluster", "DATA", aDATA, 2, eDATA)) goto err;
    pyarray_value = aDATA;
    aDATA = make_contiguous ("treecluster", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!(aRESULT = do_array_create ("treecluster", "RESULT", PyArray_LONG, 2, eRESULT))) goto err;
    if (!(aLINKDIST = do_array_create ("treecluster", "LINKDIST", PyArray_DOUBLE, 1, eLINKDIST))) goto err;
    if (DISTANCES==0) /* aDATA contains gene expression data */
    { double* paDATA = 0;
      double** ppaDATA = 0;
      long* paMASK = 0;
      long** ppaMASK = 0;
      if (!(aMASK = do_array_in ("treecluster", "MASK", MASK, PyArray_LONG))) goto err;
      if (!(aWEIGHT = do_array_in ("treecluster", "WEIGHT", WEIGHT, PyArray_DOUBLE))) goto err;
      eMASK[0] = NROWS;
      eMASK[1] = NCOLUMNS;
      eWEIGHT[0] = ((TRANSPOSE==0) ? NCOLUMNS : NROWS);
      if (!do_size_check ("treecluster", "MASK", aMASK, 2, eMASK)) goto err;
      pyarray_value = aMASK;
      aMASK = make_contiguous ("treecluster", "MASK", pyarray_value);
      Py_DECREF(pyarray_value);
      if(!aMASK) goto err;
      if (!do_size_check ("treecluster", "WEIGHT", aWEIGHT, 1, eWEIGHT)) goto err;
      pyarray_value = aWEIGHT;
      aWEIGHT = make_contiguous ("treecluster", "WEIGHT", pyarray_value);
      Py_DECREF(pyarray_value);
      if(!aWEIGHT) goto err;
      ppaDATA = (double**)malloc((size_t)NROWS*sizeof(double*));
      paDATA = (double*) (aDATA->data);
      for (ii=0; ii<NROWS; ii++) ppaDATA[ii]=&(paDATA[ii*NCOLUMNS]);
      ppaMASK = (long**)malloc((size_t)NROWS*sizeof(long*));
      paMASK = (long*) (aMASK->data);
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
          (long**) (aRESULT->data), 
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
          (long**) (aRESULT->data), 
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

extern void somcluster(long NROWS, long NCOLUMNS, double** DATA, long** MASK, double* WEIGHT, long TRANSPOSE, long NXGRID, long NYGRID, double INITTAU, long NITER, char DIST, double*** CELLDATA, long** CLUSTERID);

static PyObject*
cluster_somcluster (PyObject* unused, PyObject* args) {

    PyObject *pyfort_result;
    long NROWS;
    long NCOLUMNS;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    int eDATA[2];
    PyObject* MASK;
    PyArrayObject* aMASK;
    int eMASK[2];
    PyObject* WEIGHT;
    PyArrayObject* aWEIGHT;
    int eWEIGHT[1];
    long TRANSPOSE;
    long NXGRID;
    long NYGRID;
    double INITTAU;
    long NITER;
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
    long* paMASK;
    long** ppaMASK;
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
    eDATA[0] = NROWS;
    eDATA[1] = NCOLUMNS;
    eMASK[0] = NROWS;
    eMASK[1] = NCOLUMNS;
    eWEIGHT[0] = ((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NROWS);
    eCELLDATA[0] = NXGRID;
    eCELLDATA[1] = NYGRID;
    eCELLDATA[2] = ((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NROWS);
    eCLUSTERID[0] = ((1-TRANSPOSE)*NROWS+TRANSPOSE*NCOLUMNS);
    eCLUSTERID[1] = 2;
    if (!do_size_check ("somcluster", "DATA", aDATA, 2, eDATA)) goto err;
    pyarray_value = aDATA;
    aDATA = make_contiguous ("somcluster", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!do_size_check ("somcluster", "MASK", aMASK, 2, eMASK)) goto err;
    pyarray_value = aMASK;
    aMASK = make_contiguous ("somcluster", "MASK", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aMASK) goto err;
    if (!do_size_check ("somcluster", "WEIGHT", aWEIGHT, 1, eWEIGHT)) goto err;
    pyarray_value = aWEIGHT;
    aWEIGHT = make_contiguous ("somcluster", "WEIGHT", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aWEIGHT) goto err;
    if (!(aCELLDATA = do_array_create ("somcluster", "CELLDATA", PyArray_DOUBLE, 3, eCELLDATA))) goto err;
    if (!(aCLUSTERID = do_array_create ("somcluster", "CLUSTERID", PyArray_LONG, 2, eCLUSTERID))) goto err;
    ppaDATA = (double**)malloc((size_t)NROWS*sizeof(double*));
    ppaMASK = (long**)malloc((size_t)NROWS*sizeof(long*));
    ppaCELLDATA = (double**)malloc((size_t)NXGRID*NYGRID*sizeof(double*));
    pppaCELLDATA = (double***)malloc((size_t)NXGRID*sizeof(double**));
    paDATA = (double*) (aDATA->data);
    paMASK = (long*) (aMASK->data);
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
        (long**) (aCLUSTERID->data));
    rCELLDATA = PyArray_Return(aCELLDATA);
    rCLUSTERID = PyArray_Return(aCLUSTERID);
    Py_XDECREF((PyObject*) aDATA);
    Py_XDECREF((PyObject*) aMASK);
    Py_XDECREF((PyObject*) aWEIGHT);
    free (ppaDATA);
    free (ppaMASK);
    free (ppaCELLDATA);
    free (pppaCELLDATA);

    pyfort_result = Py_BuildValue("OO",rCELLDATA, rCLUSTERID);

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

extern double median(long N, double* DATA);

static PyObject*
cluster_median (PyObject* unused, PyObject* args) {

    PyObject *pyfort_result;
    double fortran_result;
    long N;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    int eDATA[1];
    aDATA = (PyArrayObject*) 0;

    if(!PyArg_ParseTuple(args, "O", &DATA)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("median", "DATA", DATA, PyArray_DOUBLE))) goto err;
    N = aDATA->dimensions[0];
    eDATA[0] = N;
    if (!do_size_check ("median", "DATA", aDATA, 1, eDATA)) goto err;
    pyarray_value = aDATA;
    aDATA = make_contiguous ("median", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    fortran_result = median(N, 
        (double*) (aDATA->data));
    Py_XDECREF((PyObject*) aDATA);

    pyfort_result = Py_BuildValue("d",fortran_result);

    return pyfort_result;
err:
    Py_XDECREF((PyObject*) aDATA);
    return NULL;
} 
/* end of wrapper for median */

/* mean */
static char cluster_mean__doc__[] =
"mean (data)\n"
"This function returns the mean of the 1D array data.\n";

extern double mean(long N, double* DATA);

static PyObject*
cluster_mean (PyObject* unused, PyObject* args) {

    PyObject *pyfort_result;
    double fortran_result;
    long N;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    int eDATA[1];
    aDATA = (PyArrayObject*) 0;

    if(!PyArg_ParseTuple(args, "O", &DATA)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("mean", "DATA", DATA, PyArray_DOUBLE))) goto err;
    N = aDATA->dimensions[0];
    eDATA[0] = N;
    if (!do_size_check ("mean", "DATA", aDATA, 1, eDATA)) goto err;
    pyarray_value = aDATA;
    aDATA = make_contiguous ("mean", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    fortran_result = mean(N, 
        (double*) (aDATA->data));
    Py_XDECREF((PyObject*) aDATA);

    pyfort_result = Py_BuildValue("d",fortran_result);

    return pyfort_result;
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

extern double clusterdistance(long NROWS, long NCOLUMNS, double** DATA, long** MASK, double* WEIGHT, long N1, long N2, long* INDEX1, long* INDEX2, char DIST, char METHOD, long TRANSPOSE);

static PyObject*
cluster_clusterdistance (PyObject* unused, PyObject* args) {

    PyObject *pyfort_result;
    double fortran_result;
    long NROWS;
    long NCOLUMNS;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    int eDATA[2];
    PyObject* MASK;
    PyArrayObject* aMASK;
    int eMASK[2];
    PyObject* WEIGHT;
    PyArrayObject* aWEIGHT;
    int eWEIGHT[1];
    long N1;
    long N2;
    PyObject* INDEX1;
    PyArrayObject* aINDEX1;
    int eINDEX1[1];
    PyObject* INDEX2;
    PyArrayObject* aINDEX2;
    int eINDEX2[1];
    char DIST;
    char METHOD;
    long TRANSPOSE;
    int ii;
    double* paDATA;
    double** ppaDATA;
    long* paMASK;
    long** ppaMASK;
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
    eDATA[0] = NROWS;
    eDATA[1] = NCOLUMNS;
    eMASK[0] = NROWS;
    eMASK[1] = NCOLUMNS;
    eWEIGHT[0] = TRANSPOSE*NROWS+(1-TRANSPOSE)*NCOLUMNS;
    eINDEX1[0] = N1;
    eINDEX2[0] = N2;
    if (!do_size_check ("clusterdistance", "DATA", aDATA, 2, eDATA)) goto err;
    pyarray_value = aDATA;
    aDATA = make_contiguous ("clusterdistance", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!do_size_check ("clusterdistance", "MASK", aMASK, 2, eMASK)) goto err;
    pyarray_value = aMASK;
    aMASK = make_contiguous ("clusterdistance", "MASK", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aMASK) goto err;
    if (!do_size_check ("clusterdistance", "WEIGHT", aWEIGHT, 1, eWEIGHT)) goto err;
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
    ppaMASK = (long**)malloc((size_t)NROWS*sizeof(long*));
    paDATA = (double*) (aDATA->data);
    paMASK = (long*) (aMASK->data);
    for (ii=0; ii<NROWS; ii++) ppaDATA[ii]=&(paDATA[ii*NCOLUMNS]);
    for (ii=0; ii<NROWS; ii++) ppaMASK[ii]=&(paMASK[ii*NCOLUMNS]);
    fortran_result = clusterdistance(NROWS, 
        NCOLUMNS, 
        ppaDATA, 
        ppaMASK, 
        (double*) (aWEIGHT->data), 
        N1, 
        N2, 
        (long*) (aINDEX1->data), 
        (long*) (aINDEX2->data), 
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

    pyfort_result = Py_BuildValue("d",fortran_result);

    return pyfort_result;
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

extern void getclustermean(long NCLUSTERS, long NROWS, long NCOLUMNS, double** DATA, long** MASK, long* CLUSTERID, double** CDATA, long** CMASK, long TRANSPOSE);

static PyObject*
cluster_getclustermean (PyObject* unused, PyObject* args) {

    PyObject *pyfort_result;
    long NCLUSTERS;
    long NROWS;
    long NCOLUMNS;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    int eDATA[2];
    PyObject* MASK;
    PyArrayObject* aMASK;
    int eMASK[2];
    PyObject* CLUSTERID;
    PyArrayObject* aCLUSTERID;
    int eCLUSTERID[1];
    PyArrayObject* aCDATA;
    PyObject* rCDATA;
    int eCDATA[2];
    PyArrayObject* aCMASK;
    PyObject* rCMASK;
    int eCMASK[2];
    long TRANSPOSE;
    int ii;
    double* paDATA;
    double** ppaDATA;
    long* paMASK;
    long** ppaMASK;
    double* paCDATA;
    double** ppaCDATA;
    long* paCMASK;
    long** ppaCMASK;
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
    eDATA[0] = NROWS;
    eDATA[1] = NCOLUMNS;
    eMASK[0] = NROWS;
    eMASK[1] = NCOLUMNS;
    eCLUSTERID[0] = ((1-TRANSPOSE)*NROWS+TRANSPOSE*NCOLUMNS);
    eCDATA[0] = ((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS);
    eCDATA[1] = ((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS);
    eCMASK[0] = ((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS);
    eCMASK[1] = ((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS);
    if (!do_size_check ("getclustermean", "DATA", aDATA, 2, eDATA)) goto err;
    pyarray_value = aDATA;
    aDATA = make_contiguous ("getclustermean", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!do_size_check ("getclustermean", "MASK", aMASK, 2, eMASK)) goto err;
    pyarray_value = aMASK;
    aMASK = make_contiguous ("getclustermean", "MASK", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aMASK) goto err;
    if (!do_size_check ("getclustermean", "CLUSTERID", aCLUSTERID, 1, eCLUSTERID)) goto err;
    pyarray_value = aCLUSTERID;
    aCLUSTERID = make_contiguous ("getclustermean", "CLUSTERID", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aCLUSTERID) goto err;
    if (!(aCDATA = do_array_create ("getclustermean", "CDATA", PyArray_DOUBLE, 2, eCDATA))) goto err;
    if (!(aCMASK = do_array_create ("getclustermean", "CMASK", PyArray_LONG, 2, eCMASK))) goto err;
    ppaDATA = (double**)malloc((size_t)NROWS*sizeof(double*));
    ppaMASK = (long**)malloc((size_t)NROWS*sizeof(long*));
    ppaCDATA = (double**)malloc((size_t)((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS)*sizeof(double*));
    ppaCMASK = (long**)malloc((size_t)((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS)*sizeof(long*));
    paDATA = (double*) (aDATA->data);
    paMASK = (long*) (aMASK->data);
    paCDATA = (double*) (aCDATA->data);
    paCMASK = (long*) (aCMASK->data);
    for (ii=0; ii<NROWS; ii++) ppaDATA[ii]=&(paDATA[ii*NCOLUMNS]);
    for (ii=0; ii<NROWS; ii++) ppaMASK[ii]=&(paMASK[ii*NCOLUMNS]);
    for (ii=0; ii<((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS); ii++) ppaCDATA[ii]=&(paCDATA[ii*((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS)]);
    for (ii=0; ii<((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS); ii++) ppaCMASK[ii]=&(paCMASK[ii*((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS)]);
    getclustermean(NCLUSTERS, 
        NROWS, 
        NCOLUMNS, 
        ppaDATA, 
        ppaMASK, 
        (long*) (aCLUSTERID->data), 
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

extern void getclustermedian(long NCLUSTERS, long NROWS, long NCOLUMNS, double** DATA, long** MASK, long* CLUSTERID, double** CDATA, long** CMASK, long TRANSPOSE);

static PyObject*
cluster_getclustermedian (PyObject* unused, PyObject* args) {

    PyObject *pyfort_result;
    long NCLUSTERS;
    long NROWS;
    long NCOLUMNS;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    int eDATA[2];
    PyObject* MASK;
    PyArrayObject* aMASK;
    int eMASK[2];
    PyObject* CLUSTERID;
    PyArrayObject* aCLUSTERID;
    int eCLUSTERID[1];
    PyArrayObject* aCDATA;
    PyObject* rCDATA;
    int eCDATA[2];
    PyArrayObject* aCMASK;
    PyObject* rCMASK;
    int eCMASK[2];
    long TRANSPOSE;
    int ii;
    double* paDATA;
    double** ppaDATA;
    long* paMASK;
    long** ppaMASK;
    double* paCDATA;
    double** ppaCDATA;
    long* paCMASK;
    long** ppaCMASK;
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
    eDATA[0] = NROWS;
    eDATA[1] = NCOLUMNS;
    eMASK[0] = NROWS;
    eMASK[1] = NCOLUMNS;
    eCLUSTERID[0] = ((1-TRANSPOSE)*NROWS+TRANSPOSE*NCOLUMNS);
    eCDATA[0] = ((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS);
    eCDATA[1] = ((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS);
    eCMASK[0] = ((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS);
    eCMASK[1] = ((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS);
    if (!do_size_check ("getclustermedian", "DATA", aDATA, 2, eDATA)) goto err;
    pyarray_value = aDATA;
    aDATA = make_contiguous ("getclustermedian", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!do_size_check ("getclustermedian", "MASK", aMASK, 2, eMASK)) goto err;
    pyarray_value = aMASK;
    aMASK = make_contiguous ("getclustermedian", "MASK", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aMASK) goto err;
    if (!do_size_check ("getclustermedian", "CLUSTERID", aCLUSTERID, 1, eCLUSTERID)) goto err;
    pyarray_value = aCLUSTERID;
    aCLUSTERID = make_contiguous ("getclustermedian", "CLUSTERID", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aCLUSTERID) goto err;
    if (!(aCDATA = do_array_create ("getclustermedian", "CDATA", PyArray_DOUBLE, 2, eCDATA))) goto err;
    if (!(aCMASK = do_array_create ("getclustermedian", "CMASK", PyArray_LONG, 2, eCMASK))) goto err;
    ppaDATA = (double**)malloc((size_t)NROWS*sizeof(double*));
    ppaMASK = (long**)malloc((size_t)NROWS*sizeof(long*));
    ppaCDATA = (double**)malloc((size_t)((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS)*sizeof(double*));
    ppaCMASK = (long**)malloc((size_t)((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS)*sizeof(long*));
    paDATA = (double*) (aDATA->data);
    paMASK = (long*) (aMASK->data);
    paCDATA = (double*) (aCDATA->data);
    paCMASK = (long*) (aCMASK->data);
    for (ii=0; ii<NROWS; ii++) ppaDATA[ii]=&(paDATA[ii*NCOLUMNS]);
    for (ii=0; ii<NROWS; ii++) ppaMASK[ii]=&(paMASK[ii*NCOLUMNS]);
    for (ii=0; ii<((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS); ii++) ppaCDATA[ii]=&(paCDATA[ii*((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS)]);
    for (ii=0; ii<((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS); ii++) ppaCMASK[ii]=&(paCMASK[ii*((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS)]);
    getclustermedian(NCLUSTERS, 
        NROWS, 
        NCOLUMNS, 
        ppaDATA, 
        ppaMASK, 
        (long*) (aCLUSTERID->data), 
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

static struct PyMethodDef cluster_methods[] = {
   {"kcluster", (PyCFunction) cluster_kcluster, METH_VARARGS, cluster_kcluster__doc__},
   {"treecluster", (PyCFunction) cluster_treecluster, METH_VARARGS, cluster_treecluster__doc__},
   {"somcluster", (PyCFunction) cluster_somcluster, METH_VARARGS, cluster_somcluster__doc__},
   {"median", (PyCFunction) cluster_median, METH_VARARGS, cluster_median__doc__},
   {"mean", (PyCFunction) cluster_mean, METH_VARARGS, cluster_mean__doc__},
   {"clusterdistance", (PyCFunction) cluster_clusterdistance, METH_VARARGS, cluster_clusterdistance__doc__},
   {"getclustermean", (PyCFunction) cluster_getclustermean, METH_VARARGS, cluster_getclustermean__doc__},
   {"getclustermedian", (PyCFunction) cluster_getclustermedian, METH_VARARGS, cluster_getclustermedian__doc__},
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
