#ifdef __CPLUSPLUS__
extern "C" {
#endif
#include "Python.h"
#include "Numeric/arrayobject.h"
 
static PyObject *ErrorObject;

static int array_really_contiguous(PyArrayObject *ap) {
      int sd;
      int i;

      sd = ap->descr->elsize;
      for (i = ap->nd-1; i >= 0; --i) {
              if (ap->dimensions[i] == 0) return 1; /* contiguous by definition */
              if (ap->strides[i] != sd) return 0;
              sd *= ap->dimensions[i];
      }
      return 1;
}

struct fcomplex {
    float r;
    float i;
    };
typedef struct fcomplex Py_complex_float;
#define TRANSPOSE_OPTION 0
#define MIRROR_OPTION 0
#define get_fortran_dim(v,n) v->dimensions[(n)-1]

/* 
    Built by PyFort for C compiler.
*/

static int default_option = TRANSPOSE_OPTION;
static PyObject*
set_pyfort_option (PyObject* unused, PyObject* args) {
    if(!PyArg_ParseTuple(args, "i", &default_option)) return NULL;
    Py_INCREF(Py_None);
    return Py_None;
}

static void
set_pyfort_error (char* routine, char* var, char* problem) {
    char buf[512];
    sprintf(buf, "%s, argument %s: %s", routine, var, problem);
    PyErr_SetString (ErrorObject, buf);
}

static void set_transposed_strides (PyArrayObject* ar)
{
    int m1, n1, itmp; 
    n1 = ar->nd; 
    if (n1 < 2) return;
    m1 = ar->strides[n1-1];   /* stride for one element */ 
    for(itmp=0; itmp < n1 ; itmp++) { 
        ar->strides[itmp] = m1; 
        m1 *= ar->dimensions[itmp]; 
    } 
    ar->flags &= ~CONTIGUOUS; 
}


static PyArrayObject*
transpose_array (char* rname, char* vname, PyArrayObject* ap) {
/* return transpose of ap, possibly not contiguous so as to avoid copy if we
   are transposing a previously transposed contiguous array
   This means with the transpose option on the output of one call might
   not need any copying if used as input to another call. I.e., Fortran 
   arrays stay in row-major order.
*/
    int i, n;
    PyArrayObject *ret;
    n  = ap->nd;

    /* this allocates memory for dimensions and strides (but fills them
           incorrectly), sets up descr, and points data at ap->data. */
    ret = (PyArrayObject *)PyArray_FromDimsAndData(n, ap->dimensions,
                                                    ap->descr->type_num,
                                                    ap->data);
    if (!ret) {
        set_pyfort_error (rname, vname, "Could not create descriptors for transpose.");
        return NULL;
    }
       
    /* point at true owner of memory: */
    ret->base = (PyObject *)ap;
    Py_INCREF(ap);
    for(i=0; i<n; i++) {
        ret->dimensions[i] = ap->dimensions[n - i - 1];
        ret->strides[i] = ap->strides[n - i - 1];
    }
    if (array_really_contiguous(ret)) {
        ret->flags |= CONTIGUOUS;
    } else {
        ret->flags &= ~CONTIGUOUS;
    }
    return ret;
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

static int do_size_check (char* rname, char* vname, PyArrayObject *av, int rank,  int extents[], int mirror)
{
    int size1;
    int i, j;
    char buf[512];

    size1 = av->nd;
    
    if( size1 == rank) {
        for(i=0; i < rank; ++i) {
            /* no checking on last dimension of expected size 1 */
            if (i == size1-1) {
               if (extents[i] == 1) break;
            }
            j = mirror ? size1 - 1 - i : i;
            if(av->dimensions[j] != extents[i]) 
            {
               sprintf(buf, "Incorrect extent in dimension %d (%d expected %d)",
                       i+1, av->dimensions[j], extents[i]);
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
            goto err;
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
            goto err;
        }
    } 
    return t;

err:
   return (PyArrayObject*) 0;
}

static PyArrayObject*
do_array_inout (char* rname, char* vname, PyObject *v, 
    enum PyArray_TYPES python_array_type)
{
    PyArrayObject* av;

   if (!PyArray_Check (v)) {
        set_pyfort_error(rname, vname, "Argument intent(inout) must be an array.");
        goto err;
   }
   av = (PyArrayObject*) v;
   if (av->descr->type_num != python_array_type) {
        set_pyfort_error(rname, vname, "Argument intent(inout) must be of correct type.");
        goto err;
   }
   if (!(av->flags & CONTIGUOUS))  {
       set_pyfort_error(rname, vname, "Argument intent(inout) must be contiguous.");
       goto err;
   }
   Py_INCREF(av);
   return av;
err:
   return (PyArrayObject*) 0;
}

static PyArrayObject*
do_array_create (char* rname, char* vname, enum PyArray_TYPES python_array_type, 
    int rank, int extents[], int mirror)
{
    PyArrayObject* av;
    int i, dims[7];
    if (rank > 7) {
        set_pyfort_error(rname, vname, "Too many dimensions -- limit is 7.");
        goto err;
    }
    if(mirror) {
        for(i=0; i < rank; ++i) dims[i] = extents[rank-1-i];
    } else {
        for(i=0; i < rank; ++i) dims[i] = extents[i];
    }
    av = (PyArrayObject*) PyArray_FromDims(rank, dims, python_array_type);
    if (!av) {
        set_pyfort_error(rname, vname, "Could not create array -- too big?");
        goto err;
    }
    return av;
err:
    return (PyArrayObject*) 0;
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
    int keyoption;
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
    keyoption = default_option;

    if(!PyArg_ParseTuple(args, "lOOOllcc|i", &NCLUSTERS, &DATA, &MASK, &WEIGHT, &TRANSPOSE, &NPASS, &METHOD, &DIST, &keyoption)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("kcluster", "DATA", DATA, PyArray_DOUBLE))) goto err;
    if (!(aMASK = do_array_in ("kcluster", "MASK", MASK, PyArray_LONG))) goto err;
    if (!(aWEIGHT = do_array_in ("kcluster", "WEIGHT", WEIGHT, PyArray_DOUBLE))) goto err;
    NROWS = get_fortran_dim(aDATA, 1);
    NCOLUMNS = get_fortran_dim(aDATA, 2);
    eDATA[0] = NROWS;
    eDATA[1] = NCOLUMNS;
    eMASK[0] = NROWS;
    eMASK[1] = NCOLUMNS;
    eWEIGHT[0] = ((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NROWS);
    eCLUSTERID[0] = ((1-TRANSPOSE)*NROWS+TRANSPOSE*NCOLUMNS);
    eCDATA[0] = ((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS);
    eCDATA[1] = ((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS);
    if (!do_size_check ("kcluster", "DATA", aDATA, 2, eDATA, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aDATA->nd > 1)) {
        pyarray_value = aDATA;
        aDATA = transpose_array ("kcluster", "DATA", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aDATA) goto err;
    }
    pyarray_value = aDATA;
    aDATA = make_contiguous ("kcluster", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!do_size_check ("kcluster", "MASK", aMASK, 2, eMASK, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aMASK->nd > 1)) {
        pyarray_value = aMASK;
        aMASK = transpose_array ("kcluster", "MASK", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aMASK) goto err;
    }
    pyarray_value = aMASK;
    aMASK = make_contiguous ("kcluster", "MASK", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aMASK) goto err;
    if (!do_size_check ("kcluster", "WEIGHT", aWEIGHT, 1, eWEIGHT, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aWEIGHT->nd > 1)) {
        pyarray_value = aWEIGHT;
        aWEIGHT = transpose_array ("kcluster", "WEIGHT", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aWEIGHT) goto err;
    }
    pyarray_value = aWEIGHT;
    aWEIGHT = make_contiguous ("kcluster", "WEIGHT", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aWEIGHT) goto err;
    if (!(aCLUSTERID = do_array_create ("kcluster", "CLUSTERID", PyArray_LONG, 1, eCLUSTERID, keyoption&MIRROR_OPTION))) goto err;
    if (!(aCDATA = do_array_create ("kcluster", "CDATA", PyArray_DOUBLE, 2, eCDATA, keyoption&MIRROR_OPTION))) goto err;
    ppaDATA = (double**)malloc((size_t)NROWS*sizeof(double*));
    ppaMASK = (long**)malloc((size_t)NROWS*sizeof(long*));
    ppaCDATA = (double**)malloc((size_t)((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS)*sizeof(double*));
    paDATA = (double*) (aDATA->data);
    paMASK = (long*) (aMASK->data);
    paCDATA = (double*) (aCDATA->data);
    for (ii=0; ii<NROWS; ii++) ppaDATA[ii]=&(paDATA[ii*NCOLUMNS]);
    for (ii=0; ii<NROWS; ii++) ppaMASK[ii]=&(paMASK[ii*NCOLUMNS]);
    for (ii=0; ii<((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS); ii++) ppaCDATA[ii]=&(paCDATA[ii*((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS)]);
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
"result, linkdist = pclcluster(data,mask,weight,applyscale,transpose,dist,method)\n"
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
    int keyoption;
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
    PyArrayObject* aRESULT;
    PyObject* rRESULT;
    int eRESULT[2];
    PyArrayObject* aLINKDIST;
    PyObject* rLINKDIST;
    int eLINKDIST[1];
    int ii;
    double* paDATA;
    double** ppaDATA;
    long* paMASK;
    long** ppaMASK;
    aDATA = (PyArrayObject*) 0;
    aMASK = (PyArrayObject*) 0;
    aWEIGHT = (PyArrayObject*) 0;
    aRESULT = (PyArrayObject*) 0;
    aLINKDIST = (PyArrayObject*) 0;
    keyoption = default_option;

    if(!PyArg_ParseTuple(args, "OOOllcc|i", &DATA, &MASK, &WEIGHT, &APPLYSCALE, &TRANSPOSE, &DIST, &METHOD, &keyoption)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("treecluster", "DATA", DATA, PyArray_DOUBLE))) goto err;
    if (!(aMASK = do_array_in ("treecluster", "MASK", MASK, PyArray_LONG))) goto err;
    if (!(aWEIGHT = do_array_in ("treecluster", "WEIGHT", WEIGHT, PyArray_DOUBLE))) goto err;
    NROWS = get_fortran_dim(aDATA, 1);
    NCOLUMNS = get_fortran_dim(aDATA, 2);
    eDATA[0] = NROWS;
    eDATA[1] = NCOLUMNS;
    eMASK[0] = NROWS;
    eMASK[1] = NCOLUMNS;
    eWEIGHT[0] = ((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NROWS);
    eRESULT[0] = ((1-TRANSPOSE)*NROWS+TRANSPOSE*NCOLUMNS-1);
    eRESULT[1] = 2;
    eLINKDIST[0] = ((1-TRANSPOSE)*NROWS+TRANSPOSE*NCOLUMNS-1);
    if (!do_size_check ("treecluster", "DATA", aDATA, 2, eDATA, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aDATA->nd > 1)) {
        pyarray_value = aDATA;
        aDATA = transpose_array ("treecluster", "DATA", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aDATA) goto err;
    }
    pyarray_value = aDATA;
    aDATA = make_contiguous ("treecluster", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!do_size_check ("treecluster", "MASK", aMASK, 2, eMASK, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aMASK->nd > 1)) {
        pyarray_value = aMASK;
        aMASK = transpose_array ("treecluster", "MASK", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aMASK) goto err;
    }
    pyarray_value = aMASK;
    aMASK = make_contiguous ("treecluster", "MASK", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aMASK) goto err;
    if (!do_size_check ("treecluster", "WEIGHT", aWEIGHT, 1, eWEIGHT, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aWEIGHT->nd > 1)) {
        pyarray_value = aWEIGHT;
        aWEIGHT = transpose_array ("treecluster", "WEIGHT", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aWEIGHT) goto err;
    }
    pyarray_value = aWEIGHT;
    aWEIGHT = make_contiguous ("treecluster", "WEIGHT", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aWEIGHT) goto err;
    if (!(aRESULT = do_array_create ("treecluster", "RESULT", PyArray_LONG, 2, eRESULT, keyoption&MIRROR_OPTION))) goto err;
    if (!(aLINKDIST = do_array_create ("treecluster", "LINKDIST", PyArray_DOUBLE, 1, eLINKDIST, keyoption&MIRROR_OPTION))) goto err;
    ppaDATA = (double**)malloc((size_t)NROWS*sizeof(double*));
    ppaMASK = (long**)malloc((size_t)NROWS*sizeof(long*));
    paDATA = (double*) (aDATA->data);
    paMASK = (long*) (aMASK->data);
    for (ii=0; ii<NROWS; ii++) ppaDATA[ii]=&(paDATA[ii*NCOLUMNS]);
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
    rRESULT = PyArray_Return(aRESULT);
    rLINKDIST = PyArray_Return(aLINKDIST);
    Py_XDECREF((PyObject*) aDATA);
    Py_XDECREF((PyObject*) aMASK);
    Py_XDECREF((PyObject*) aWEIGHT);
    free (ppaDATA);
    free (ppaMASK);

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
    int keyoption;
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
    keyoption = default_option;

    if(!PyArg_ParseTuple(args, "OOOllldlc|i", &DATA, &MASK, &WEIGHT, &TRANSPOSE, &NXGRID, &NYGRID, &INITTAU, &NITER, &DIST, &keyoption)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("somcluster", "DATA", DATA, PyArray_DOUBLE))) goto err;
    if (!(aMASK = do_array_in ("somcluster", "MASK", MASK, PyArray_LONG))) goto err;
    if (!(aWEIGHT = do_array_in ("somcluster", "WEIGHT", WEIGHT, PyArray_DOUBLE))) goto err;
    NROWS = get_fortran_dim(aDATA, 1);
    NCOLUMNS = get_fortran_dim(aDATA, 2);
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
    if (!do_size_check ("somcluster", "DATA", aDATA, 2, eDATA, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aDATA->nd > 1)) {
        pyarray_value = aDATA;
        aDATA = transpose_array ("somcluster", "DATA", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aDATA) goto err;
    }
    pyarray_value = aDATA;
    aDATA = make_contiguous ("somcluster", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!do_size_check ("somcluster", "MASK", aMASK, 2, eMASK, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aMASK->nd > 1)) {
        pyarray_value = aMASK;
        aMASK = transpose_array ("somcluster", "MASK", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aMASK) goto err;
    }
    pyarray_value = aMASK;
    aMASK = make_contiguous ("somcluster", "MASK", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aMASK) goto err;
    if (!do_size_check ("somcluster", "WEIGHT", aWEIGHT, 1, eWEIGHT, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aWEIGHT->nd > 1)) {
        pyarray_value = aWEIGHT;
        aWEIGHT = transpose_array ("somcluster", "WEIGHT", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aWEIGHT) goto err;
    }
    pyarray_value = aWEIGHT;
    aWEIGHT = make_contiguous ("somcluster", "WEIGHT", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aWEIGHT) goto err;
    if (!(aCELLDATA = do_array_create ("somcluster", "CELLDATA", PyArray_DOUBLE, 3, eCELLDATA, keyoption&MIRROR_OPTION))) goto err;
    if (!(aCLUSTERID = do_array_create ("somcluster", "CLUSTERID", PyArray_LONG, 2, eCLUSTERID, keyoption&MIRROR_OPTION))) goto err;
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
    int keyoption;
    double fortran_result;
    long N;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    int eDATA[1];
    aDATA = (PyArrayObject*) 0;
    keyoption = default_option;

    if(!PyArg_ParseTuple(args, "O|i", &DATA, &keyoption)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("median", "DATA", DATA, PyArray_DOUBLE))) goto err;
    N = get_fortran_dim(aDATA, 1);
    eDATA[0] = N;
    if (!do_size_check ("median", "DATA", aDATA, 1, eDATA, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aDATA->nd > 1)) {
        pyarray_value = aDATA;
        aDATA = transpose_array ("median", "DATA", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aDATA) goto err;
    }
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
    int keyoption;
    double fortran_result;
    long N;
    PyArrayObject* pyarray_value;
    PyObject* DATA;
    PyArrayObject* aDATA;
    int eDATA[1];
    aDATA = (PyArrayObject*) 0;
    keyoption = default_option;

    if(!PyArg_ParseTuple(args, "O|i", &DATA, &keyoption)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("mean", "DATA", DATA, PyArray_DOUBLE))) goto err;
    N = get_fortran_dim(aDATA, 1);
    eDATA[0] = N;
    if (!do_size_check ("mean", "DATA", aDATA, 1, eDATA, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aDATA->nd > 1)) {
        pyarray_value = aDATA;
        aDATA = transpose_array ("mean", "DATA", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aDATA) goto err;
    }
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
    int keyoption;
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
    keyoption = default_option;

    if(!PyArg_ParseTuple(args, "OOOOOccl|i", &DATA, &MASK, &WEIGHT, &INDEX1, &INDEX2, &DIST, &METHOD, &TRANSPOSE, &keyoption)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("clusterdistance", "DATA", DATA, PyArray_DOUBLE))) goto err;
    if (!(aMASK = do_array_in ("clusterdistance", "MASK", MASK, PyArray_LONG))) goto err;
    if (!(aWEIGHT = do_array_in ("clusterdistance", "WEIGHT", WEIGHT, PyArray_DOUBLE))) goto err;
    if (!(aINDEX1 = do_array_in ("clusterdistance", "INDEX1", INDEX1, PyArray_LONG))) goto err;
    if (!(aINDEX2 = do_array_in ("clusterdistance", "INDEX2", INDEX2, PyArray_LONG))) goto err;
    NROWS = get_fortran_dim(aDATA, 1);
    NCOLUMNS = get_fortran_dim(aDATA, 2);
    N1 = get_fortran_dim(aINDEX1, 1);
    N2 = get_fortran_dim(aINDEX2, 1);
    eDATA[0] = NROWS;
    eDATA[1] = NCOLUMNS;
    eMASK[0] = NROWS;
    eMASK[1] = NCOLUMNS;
    eWEIGHT[0] = TRANSPOSE*NROWS+(1-TRANSPOSE)*NCOLUMNS;
    eINDEX1[0] = N1;
    eINDEX2[0] = N2;
    if (!do_size_check ("clusterdistance", "DATA", aDATA, 2, eDATA, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aDATA->nd > 1)) {
        pyarray_value = aDATA;
        aDATA = transpose_array ("clusterdistance", "DATA", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aDATA) goto err;
    }
    pyarray_value = aDATA;
    aDATA = make_contiguous ("clusterdistance", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!do_size_check ("clusterdistance", "MASK", aMASK, 2, eMASK, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aMASK->nd > 1)) {
        pyarray_value = aMASK;
        aMASK = transpose_array ("clusterdistance", "MASK", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aMASK) goto err;
    }
    pyarray_value = aMASK;
    aMASK = make_contiguous ("clusterdistance", "MASK", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aMASK) goto err;
    if (!do_size_check ("clusterdistance", "WEIGHT", aWEIGHT, 1, eWEIGHT, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aWEIGHT->nd > 1)) {
        pyarray_value = aWEIGHT;
        aWEIGHT = transpose_array ("clusterdistance", "WEIGHT", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aWEIGHT) goto err;
    }
    pyarray_value = aWEIGHT;
    aWEIGHT = make_contiguous ("clusterdistance", "WEIGHT", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aWEIGHT) goto err;
    if (!do_size_check ("clusterdistance", "INDEX1", aINDEX1, 1, eINDEX1, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aINDEX1->nd > 1)) {
        pyarray_value = aINDEX1;
        aINDEX1 = transpose_array ("clusterdistance", "INDEX1", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aINDEX1) goto err;
    }
    pyarray_value = aINDEX1;
    aINDEX1 = make_contiguous ("clusterdistance", "INDEX1", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aINDEX1) goto err;
    if (!do_size_check ("clusterdistance", "INDEX2", aINDEX2, 1, eINDEX2, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aINDEX2->nd > 1)) {
        pyarray_value = aINDEX2;
        aINDEX2 = transpose_array ("clusterdistance", "INDEX2", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aINDEX2) goto err;
    }
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
    int keyoption;
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
    keyoption = default_option;

    if(!PyArg_ParseTuple(args, "lOOOl|i", &NCLUSTERS, &DATA, &MASK, &CLUSTERID, &TRANSPOSE, &keyoption)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("getclustermean", "DATA", DATA, PyArray_DOUBLE))) goto err;
    if (!(aMASK = do_array_in ("getclustermean", "MASK", MASK, PyArray_LONG))) goto err;
    if (!(aCLUSTERID = do_array_in ("getclustermean", "CLUSTERID", CLUSTERID, PyArray_LONG))) goto err;
    NROWS = get_fortran_dim(aDATA, 1);
    NCOLUMNS = get_fortran_dim(aDATA, 2);
    eDATA[0] = NROWS;
    eDATA[1] = NCOLUMNS;
    eMASK[0] = NROWS;
    eMASK[1] = NCOLUMNS;
    eCLUSTERID[0] = ((1-TRANSPOSE)*NROWS+TRANSPOSE*NCOLUMNS);
    eCDATA[0] = ((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS);
    eCDATA[1] = ((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS);
    eCMASK[0] = ((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS);
    eCMASK[1] = ((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS);
    if (!do_size_check ("getclustermean", "DATA", aDATA, 2, eDATA, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aDATA->nd > 1)) {
        pyarray_value = aDATA;
        aDATA = transpose_array ("getclustermean", "DATA", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aDATA) goto err;
    }
    pyarray_value = aDATA;
    aDATA = make_contiguous ("getclustermean", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!do_size_check ("getclustermean", "MASK", aMASK, 2, eMASK, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aMASK->nd > 1)) {
        pyarray_value = aMASK;
        aMASK = transpose_array ("getclustermean", "MASK", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aMASK) goto err;
    }
    pyarray_value = aMASK;
    aMASK = make_contiguous ("getclustermean", "MASK", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aMASK) goto err;
    if (!do_size_check ("getclustermean", "CLUSTERID", aCLUSTERID, 1, eCLUSTERID, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aCLUSTERID->nd > 1)) {
        pyarray_value = aCLUSTERID;
        aCLUSTERID = transpose_array ("getclustermean", "CLUSTERID", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aCLUSTERID) goto err;
    }
    pyarray_value = aCLUSTERID;
    aCLUSTERID = make_contiguous ("getclustermean", "CLUSTERID", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aCLUSTERID) goto err;
    if (!(aCDATA = do_array_create ("getclustermean", "CDATA", PyArray_DOUBLE, 2, eCDATA, keyoption&MIRROR_OPTION))) goto err;
    if (!(aCMASK = do_array_create ("getclustermean", "CMASK", PyArray_LONG, 2, eCMASK, keyoption&MIRROR_OPTION))) goto err;
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
    int keyoption;
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
    keyoption = default_option;

    if(!PyArg_ParseTuple(args, "lOOOl|i", &NCLUSTERS, &DATA, &MASK, &CLUSTERID, &TRANSPOSE, &keyoption)) {
        return NULL;
    }
    if (!(aDATA = do_array_in ("getclustermedian", "DATA", DATA, PyArray_DOUBLE))) goto err;
    if (!(aMASK = do_array_in ("getclustermedian", "MASK", MASK, PyArray_LONG))) goto err;
    if (!(aCLUSTERID = do_array_in ("getclustermedian", "CLUSTERID", CLUSTERID, PyArray_LONG))) goto err;
    NROWS = get_fortran_dim(aDATA, 1);
    NCOLUMNS = get_fortran_dim(aDATA, 2);
    eDATA[0] = NROWS;
    eDATA[1] = NCOLUMNS;
    eMASK[0] = NROWS;
    eMASK[1] = NCOLUMNS;
    eCLUSTERID[0] = ((1-TRANSPOSE)*NROWS+TRANSPOSE*NCOLUMNS);
    eCDATA[0] = ((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS);
    eCDATA[1] = ((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS);
    eCMASK[0] = ((1-TRANSPOSE)*NCLUSTERS+TRANSPOSE*NROWS);
    eCMASK[1] = ((1-TRANSPOSE)*NCOLUMNS+TRANSPOSE*NCLUSTERS);
    if (!do_size_check ("getclustermedian", "DATA", aDATA, 2, eDATA, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aDATA->nd > 1)) {
        pyarray_value = aDATA;
        aDATA = transpose_array ("getclustermedian", "DATA", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aDATA) goto err;
    }
    pyarray_value = aDATA;
    aDATA = make_contiguous ("getclustermedian", "DATA", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aDATA) goto err;
    if (!do_size_check ("getclustermedian", "MASK", aMASK, 2, eMASK, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aMASK->nd > 1)) {
        pyarray_value = aMASK;
        aMASK = transpose_array ("getclustermedian", "MASK", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aMASK) goto err;
    }
    pyarray_value = aMASK;
    aMASK = make_contiguous ("getclustermedian", "MASK", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aMASK) goto err;
    if (!do_size_check ("getclustermedian", "CLUSTERID", aCLUSTERID, 1, eCLUSTERID, keyoption & MIRROR_OPTION)) goto err;
    if ((keyoption & TRANSPOSE_OPTION) && (aCLUSTERID->nd > 1)) {
        pyarray_value = aCLUSTERID;
        aCLUSTERID = transpose_array ("getclustermedian", "CLUSTERID", pyarray_value);
        Py_DECREF(pyarray_value);
        if(!aCLUSTERID) goto err;
    }
    pyarray_value = aCLUSTERID;
    aCLUSTERID = make_contiguous ("getclustermedian", "CLUSTERID", pyarray_value);
    Py_DECREF(pyarray_value);
    if(!aCLUSTERID) goto err;
    if (!(aCDATA = do_array_create ("getclustermedian", "CDATA", PyArray_DOUBLE, 2, eCDATA, keyoption&MIRROR_OPTION))) goto err;
    if (!(aCMASK = do_array_create ("getclustermedian", "CMASK", PyArray_LONG, 2, eCMASK, keyoption&MIRROR_OPTION))) goto err;
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
   {"set_pyfort_option", (PyCFunction) set_pyfort_option, METH_VARARGS, 
           "set_pyfort_option (value) sets default value of option keyword."},
   {NULL,          NULL, 0, NULL}/* sentinel */
};

static char cluster_module_documentation[] =
"Fortran interface module cluster";

void initcluster(void)
{
        PyObject *m, *d, *j;
 
        import_array ();
        m = Py_InitModule4("cluster", cluster_methods,
                cluster_module_documentation,
                (PyObject*)NULL,PYTHON_API_VERSION);
 
        d = PyModule_GetDict(m);
        ErrorObject = PyString_FromString("cluster.error");
        PyDict_SetItemString(d, "error", ErrorObject);
        j = PyInt_FromLong((long) TRANSPOSE_OPTION);
        PyDict_SetItemString(d, "TRANSPOSE", j);
        Py_XDECREF(j);
        j = PyInt_FromLong((long) MIRROR_OPTION);
        PyDict_SetItemString(d, "MIRROR", j);
        Py_XDECREF(j);
        j = PyInt_FromLong(0L);
        PyDict_SetItemString(d, "NONE", j);
        Py_XDECREF(j);

        if (PyErr_Occurred()) {
            Py_FatalError("can't initialize module cluster");
        }
}

/* C++ trailer */
#ifdef __CPLUSCPLUS__
}
#endif
