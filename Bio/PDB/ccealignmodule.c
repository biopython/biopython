//
//  Original Code
//      Copyright (C) Jason Vertrees
//  Modifications
//      Copyright (C) Joao Rodrigues.
//      Modifications include removal of RMSD calculation code and associated
//      dependencies. Output of the module is now the best paths.
//
//  All rights reserved.
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are
//  met:
//
//      * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
//      * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in
//      the documentation and/or other materials provided with the
//      distribution.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
//  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
//  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
//  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
//  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//  The following notice is provided since the code was adapted from
//  open-source Pymol.

//  Open-Source PyMOL Copyright Notice
//  ==================================

//  The Open-Source PyMOL source code is copyrighted, but you can freely
//  use and copy it as long as you don't change or remove any of the
//  Copyright notices. The Open-Source PyMOL product is made available
//  under the following open-source license terms:

//  ----------------------------------------------------------------------
//  Open-Source PyMOL is Copyright (C) Schrodinger, LLC.

//  All Rights Reserved

//  Permission to use, copy, modify, distribute, and distribute modified
//  versions of this software and its built-in documentation for any
//  purpose and without fee is hereby granted, provided that the above
//  copyright notice appears in all copies and that both the copyright
//  notice and this permission notice appear in supporting documentation,
//  and that the name of Schrodinger, LLC not be used in advertising or
//  publicity pertaining to distribution of the software without specific,
//  written prior permission.

//  SCHRODINGER, LLC DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
//  INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN
//  NO EVENT SHALL SCHRODINGER, LLC BE LIABLE FOR ANY SPECIAL, INDIRECT OR
//  CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
//  OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
//  OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE
//  USE OR PERFORMANCE OF THIS SOFTWARE.
//  ----------------------------------------------------------------------

//  PyMOL Trademark Notice
//  ======================

//  PyMOL(TM) is a trademark of Schrodinger, LLC. Derivative
//  software which contains PyMOL source code must be plainly
//  distinguished from any and all PyMOL products distributed by Schrodinger,
//  LLC in all publicity, advertising, and documentation.
//  The slogans, "Includes PyMOL(TM).", "Based on PyMOL(TM) technology.",
//  "Contains PyMOL(TM) source code.", and "Built using PyMOL(TM).", may
//  be used in advertising, publicity, and documentation of derivative
//  software provided that the notice, "PyMOL is a trademark of Schrodinger,
//  LLC.", is included in a footnote or at the end of the
//  document.

//  All other endorsements employing the PyMOL trademark require specific,
//  written prior permission.
//

#include "Python.h"

#define MAX_PATHS 20

// Typical XYZ point and array of points
typedef struct {
    double x;
    double y;
    double z;
} cePoint, *pcePoint;

// An AFP (aligned fragment pair), and list/pointer
typedef struct {
    int first;
    int second;
} afp, *path, **pathCache;

// Calculate distance matrix
static double **
calcDM(pcePoint coords, int len)
{
    double **dm = (double **)malloc(sizeof(double *) * len);

    for (int i = 0; i < len; i++) {
        dm[i] = (double *)malloc(sizeof(double) * len);
    }
    for (int row = 0; row < len; row++) {
        for (int col = row; col < len; col++) {
            double xd = coords[row].x - coords[col].x;
            double yd = coords[row].y - coords[col].y;
            double zd = coords[row].z - coords[col].z;
            double distance = sqrt(xd * xd + yd * yd + zd * zd);

            dm[row][col] = dm[col][row] = distance;
        }
    }

    return dm;
}

//
// This similarity corresponds to distance measure (i) or equation (6) in the paper.
//
static double
similarityI(
    double **dA,
    double **dB,
    const int iA,
    const int iB,
    const int jA,
    const int jB,
    const int fragmentSize)
{
    const int m = fragmentSize;
    double similarity = fabs(dA[iA][jA] - dB[iB][jB]) +
        fabs(dA[iA + m-1][jA + m-1] - dB[iB + m-1][jB + m-1]);

    for (int k = 1; k < m - 1; k++) {
        similarity += fabs(dA[iA + k][jA + m-1 - k] - fabs(dB[iB + k][jB + m-1 - k]));
    }

    return -similarity / m;
}

//
// This similarity corresponds to distance measure (ii) or equation (7) in the paper.
//
static double
similarityII(
    double **dA,
    double **dB,
    const int iA,
    const int iB,
    const int fragmentSize)
{
    double similarity = 0.0;
    // Term count is the closed form of the number of terms in the summation
    const int termCount = (fragmentSize - 1) * (fragmentSize - 2) / 2;

    for (int k = 0; k < fragmentSize - 2; k++) {
        for (int l = k + 2; l < fragmentSize; l++) {
            similarity +=
                fabs(dA[iA + k][iA + l] - dB[iB + k][iB + l]);
        }
    }

    return -similarity / termCount;
}

// Calculate similarity matrix
static double **
calcS(double **dA, double **dB, const int lenA, const int lenB, const int fragmentSize)
{
    // Initialize the 2D similarity matrix
    const int rowCount = lenA - fragmentSize + 1;
    const int colCount = lenB - fragmentSize + 1;
    double **S = (double **)malloc(sizeof(double *) * rowCount);

    for (int i = 0; i < rowCount; i++) {
        S[i] = (double *)malloc(sizeof(double) * colCount);
    }

    //
    // This is where the magic of CE comes out.  In the similarity matrix,
    // for each i and j, the value of S[i][j] is how well the residues
    // i - i+fragmentSize in protein A, match to residues j - j+fragmentSize in protein
    // B.  A value of 0 means absolute match; a value >> 1 means bad match.
    //
    for (int iA = 0; iA < rowCount; iA++) {
        for (int iB = 0; iB < colCount; iB++) {
            S[iA][iB] = similarityII(dA, dB, iA, iB, fragmentSize);
        }
    }

    return S;
}

static pcePoint
getCoords(PyObject *L, int length)
{
    // Make space for the current coords
    pcePoint coords = (pcePoint)malloc(sizeof(cePoint) * length);

    if (!coords)
        return NULL;

    // loop through the arguments, pulling out the
    // XYZ coordinates.
    for (int i = 0; i < length; i++) {
        PyObject *curCoord = PyList_GetItem(L, i);
        Py_INCREF(curCoord);

        PyObject *curVal = PyList_GetItem(curCoord, 0);
        Py_INCREF(curVal);
        coords[i].x = PyFloat_AsDouble(curVal);
        Py_DECREF(curVal);

        curVal = PyList_GetItem(curCoord, 1);
        Py_INCREF(curVal);
        coords[i].y = PyFloat_AsDouble(curVal);
        Py_DECREF(curVal);

        curVal = PyList_GetItem(curCoord, 2);
        Py_INCREF(curVal);
        coords[i].z = PyFloat_AsDouble(curVal);

        Py_DECREF(curVal);
        Py_DECREF(curCoord);
    }

    return coords;
}

// Find the best N alignment paths
static PyObject *
findPath(
    double **S,
    double **dA,
    double **dB,
    const int lenA,
    const int lenB,
    const int fragmentSize,
    const int gapMax)
{
    const double D0 = -3.0;
    const double D1 = -4.0;

    // Length of longest possible alignment
    const int smaller = (lenA < lenB) ? lenA : lenB;

    // For storing the best N paths
    int bufferIndex = MAX_PATHS - 1; // We want to insert the first item at 0
    int bufferSize = 0;
    int lenBuffer[MAX_PATHS];
    double similarityBuffer[MAX_PATHS];
    pathCache pathBuffer = (pathCache)malloc(sizeof(path *) * MAX_PATHS);

    for (int i = 0; i < MAX_PATHS; i++) {
        // Initialize the paths
        similarityBuffer[i] = -1e6;
        lenBuffer[i] = 0;
        pathBuffer[i] = 0;
    }

    //======================================================================
    // Start the search through the similarity matrix.
    //
    for (int iA = 0; iA <= lenA - fragmentSize; iA++) {
        if (iA > lenA - fragmentSize * (lenBuffer[bufferIndex] - 1))
            break;

        for (int iB = 0; iB <= lenB - fragmentSize; iB++) {
            if (S[iA][iB] <= D0)
                continue;
            if (iB > lenB - fragmentSize * (lenBuffer[bufferIndex] - 1))
                break;

            // Initialize current path
            path curPath = (path)malloc(sizeof(afp) * smaller);
            int curPathLength = 1;
            double curPathSimilarity = S[iA][iB];

            curPath[0].first = iA;
            curPath[0].second = iB;

            for (int i = 1; i < smaller; i++) {
                curPath[i].first = -1;
                curPath[i].second = -1;
            }

            //
            // Build the best path starting from iA, iB
            //
            while (1) {
                double gapBestSimilarity = -1e6;
                int gapBestIndex = -1;

                //
                // Check all possible gaps from here
                //
                for (int g = 0; g < (gapMax * 2) + 1; g++) {
                    int jA = curPath[curPathLength - 1].first + fragmentSize;
                    int jB = curPath[curPathLength - 1].second + fragmentSize;

                    if ((g + 1) % 2 == 0) {
                        jA += (g + 1) / 2;
                    } else { // ( g odd )
                        jB += (g + 1) / 2;
                    }

                    // Following are three heuristics to ensure high quality
                    // long paths and make sure we don't run over the end of
                    // the S, matrix.

                    // 1st: If jA or jB is at the end of the similarity matrix
                    if (jA > lenA - fragmentSize || jB > lenB - fragmentSize)
                        continue;
                    // 2nd: If this candidate AFP is bad, ignore it.
                    if (S[jA][jB] <= D0)
                        continue;

                    double curSimilarity = 0.0;

                    for (int s = 0; s < curPathLength; s++) {
                        const int iA = curPath[s].first;
                        const int iB = curPath[s].second;
                        curSimilarity += similarityI(dA, dB, iA, iB, jA, jB, fragmentSize);
                    }

                    curSimilarity /= curPathLength;

                    // store GAPPED best
                    if (curSimilarity > D1 && curSimilarity > gapBestSimilarity) {
                        curPath[curPathLength].first = jA;
                        curPath[curPathLength].second = jB;
                        gapBestSimilarity = curSimilarity;
                        gapBestIndex = g;
                    }
                } /// ROF -- END GAP SEARCHING

                if (gapBestIndex == -1) {
                    // if here, then there was no good candidate AFP,
                    // so quit and restart from starting point
                    break;
                }

                int gA, gB;
                int jGap = (gapBestIndex + 1) / 2;

                if ((gapBestIndex + 1) % 2 == 0) {
                    gA = curPath[curPathLength - 1].first + fragmentSize + jGap;
                    gB = curPath[curPathLength - 1].second + fragmentSize;
                }
                else {
                    gA = curPath[curPathLength - 1].first + fragmentSize;
                    gB =
                        curPath[curPathLength - 1].second + fragmentSize + jGap;
                }

                // TODO: This implementation calculates the average inter-residue distance difference.
                // Instead, the implementation should calculate the average similarity as described in the paper.
                const double n = (double) curPathLength;
                const int m = fragmentSize;
                const int w = (m - 1) * (m - 2) / 2; // w is the number of terms in the similarity summation

                const double curTermCount = n * w + m * n * (n - 1) / 2;
                const double addTermCount = w + m * n;
                const double addSimilarity = (gapBestSimilarity * m * n + S[gA][gB] * w) / addTermCount;
                const double newTermCount = curTermCount + addTermCount;
                const double newSimilarity = (curPathSimilarity * curTermCount + addSimilarity * addTermCount) / newTermCount;

                if (newSimilarity > D1) {
                    curPathSimilarity = newSimilarity;
                    curPathLength++;
                }
                else {
                    // Heuristic -- path is getting sloppy, stop looking
                    break;
                }
            } /// END WHILE

            //
            // At this point, we've found the best path starting at iA, iB.
            //
            if (curPathLength > lenBuffer[bufferIndex] ||
                (curPathLength == lenBuffer[bufferIndex] &&
                 curPathSimilarity > similarityBuffer[bufferIndex])) {
                bufferIndex += 1;
                bufferIndex %= MAX_PATHS;
                bufferSize =
                    (bufferSize < MAX_PATHS) ? bufferSize + 1 : MAX_PATHS;
                path pathCopy = (path)malloc(sizeof(afp) * smaller);

                for (int i = 0; i < smaller; i++) {
                    pathCopy[i].first = curPath[i].first;
                    pathCopy[i].second = curPath[i].second;
                }
                if (pathBuffer[bufferIndex]) {
                    free(pathBuffer[bufferIndex]);
                }
                pathBuffer[bufferIndex] = pathCopy;
                lenBuffer[bufferIndex] = curPathLength;
                similarityBuffer[bufferIndex] = curPathSimilarity;
            }

            free(curPath);
            curPath = 0;
        } // ROF -- end for iB
    }     // ROF -- end for iA

    // To make it simpler to use this code and more portable, we are decoupling
    // the path finding (the actual CEAlign innovation) from the RMSD
    // calculation.
    //
    // As such, we return the N best paths to Python-land. Since the paths are
    // encoded as structs, it's simpler to return the each path as a list of
    // lists with the corresponding atom indices. e.g. [path1, path2, path3,
    // ..., pathN], where pathN is defined as,
    // [[Ai, Aj, Ak, ...], [Bi, Bj, Bk, ...], where An and Bn are equivalent
    // coordinates for structures A and B.

    // List to store all paths
    PyObject *result = PyList_New(bufferSize);
    Py_INCREF(result);

    for (int o = 0; o < bufferSize; o++) {
        // Make a new list to store this path
        PyObject *pathAList = PyList_New(0);
        PyObject *pathBList = PyList_New(0);
        Py_INCREF(pathAList);
        Py_INCREF(pathBList);

        for (int j = 0; j < smaller; j++) {
            if (pathBuffer[o][j].first != -1) {
                const int idxA = pathBuffer[o][j].first;
                const int idxB = pathBuffer[o][j].second;

                for (int k = 0; k < fragmentSize; k++) {
                    PyObject *v = Py_BuildValue("i", idxA + k);
                    PyList_Append(pathAList, v);
                    Py_DECREF(v);
                    v = Py_BuildValue("i", idxB + k);
                    PyList_Append(pathBList, v);
                    Py_DECREF(v);
                }
            }
            else {
                break;
            }
        }

        PyObject *pairList = Py_BuildValue("[NN]", pathAList, pathBList);
        Py_INCREF(pairList);
        PyList_SET_ITEM(result, o, pairList);
    }

    // Free memory
    free(pathBuffer);

    return result;
}

// Main Function
PyObject *
PyCealign(PyObject *Py_UNUSED(self), PyObject *args)
{
    int fragmentSize = 8;
    int gapMax = 30;
    double **dmA, **dmB, **S;

    PyObject *listA, *listB, *result;

    /* Unpack the arguments from Python */
    PyArg_ParseTuple(args, "OO|ii", &listA, &listB, &fragmentSize, &gapMax);

    /* Get the list lengths */
    const int lenA = (int)PyList_Size(listA);
    const int lenB = (int)PyList_Size(listB);

    /* get the coodinates from the Python objects */
    pcePoint coordsA = (pcePoint)getCoords(listA, lenA);
    pcePoint coordsB = (pcePoint)getCoords(listB, lenB);

    /* calculate the distance matrix for each protein */
    dmA = (double **)calcDM(coordsA, lenA);
    dmB = (double **)calcDM(coordsB, lenB);

    /* calculate the CE Similarity matrix */
    S = (double **)calcS(dmA, dmB, lenA, lenB, fragmentSize);

    // Calculate Top N Paths
    result = (PyObject *)findPath(S, dmA, dmB, lenA, lenB, fragmentSize, gapMax);

    /* release memory */
    free(coordsA);
    free(coordsB);

    /* distance matrices	 */
    for (int i = 0; i < lenA; i++)
        free(dmA[i]);
    free(dmA);

    for (int i = 0; i < lenB; i++)
        free(dmB[i]);
    free(dmB);

    // Similarity matrix
    for (int i = 0; i <= lenA - fragmentSize; i++)
        free(S[i]);
    free(S);

    return result;
}

//
// Python Interface
//
PyDoc_STRVAR(method_doc,
"run_cealign(coordsA, coordsB, fragmentSize, gapMax) -> list\
\n\n\
Find the optimal alignments between two structures, using CEAlign.\
\n\n\
Arguments:\n\
- listA: List of lists with coordinates for structure A.\n\
- listB: List of lists with coordinates for structure B.\n\
- fragmentSize: Size of fragments to be used in alignment.\n\
- gapMax: Maximum gap allowed between two aligned fragment pairs.");

static PyMethodDef CEAlignMethods[] = {
    {"run_cealign", PyCealign, METH_VARARGS, method_doc},
    {NULL, NULL, 0, NULL}
};

PyDoc_STRVAR(module_doc,
"Pairwise structure alignment of 3D structures using combinatorial extension.\
\n\n\
This module implements a single function: run_cealign. \
Refer to its docstring for more documentation on usage and implementation.");

PyObject *PyInit_ccealign(void)
{
    static struct PyModuleDef moduledef = {PyModuleDef_HEAD_INIT,
                                           "ccealign",
                                           module_doc,
                                           -1,
                                           CEAlignMethods,
                                           NULL,
                                           NULL,
                                           NULL,
                                           NULL};
    return PyModule_Create(&moduledef);
}
