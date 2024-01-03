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

// Calculate Distance Matrix
double **
calcDM(pcePoint coords, int len) {

    double **dm = (double **)malloc(sizeof(double *) * len);
    for (int i = 0; i < len; i++)
        dm[i] = (double *)malloc(sizeof(double) * len);

    for (int row = 0; row < len; row++) {
        for (int col = row; col < len; col++) {
            double xd = coords[row].x - coords[col].x;
            double yd = coords[row].y - coords[col].y;
            double zd = coords[row].z - coords[col].z;
            double distsq = (xd * xd) + (yd * yd) + (zd * zd);
            dm[row][col] = dm[col][row] = sqrt(distsq);
        }
    }
    return dm;
}

// Calculate Score Matrix
double **
calcS(double **d1, double **d2, int lenA, int lenB, int wSize) {
    int i;
    double winSize = (double)wSize;

    // initialize the 2D similarity matrix
    double **S = (double **)malloc(sizeof(double *) * lenA);
    for (i = 0; i < lenA; i++) {
        S[i] = (double *)malloc(sizeof(double) * lenB);
    }

    double sumSize = (winSize - 1.0) * (winSize - 2.0) / 2.0;
    //
    // This is where the magic of CE comes out.  In the similarity matrix,
    // for each i and j, the value of ceSIM[i][j] is how well the residues
    // i - i+winSize in protein A, match to residues j - j+winSize in protein
    // B.  A value of 0 means absolute match; a value >> 1 means bad match.
    //
    int lenA_m_wSize = lenA - wSize;
    int lenB_m_wSize = lenB - wSize;
    int iA, iB, row, col;
    for (iA = 0; iA < lenA; iA++) {
        for (iB = 0; iB < lenB; iB++) {
            S[iA][iB] = -1.0;
            if (iA > lenA_m_wSize || iB > lenB_m_wSize)
                continue;

            double score = 0.0;

            //
            // We always skip the calculation of the distance from THIS
            // residue, to the next residue.  This is a time-saving heur-
            // istic decision.  Almost all alpha carbon bonds of neighboring
            // residues is 3.8 Angstroms.  Due to entropy, S = -k ln pi * pi,
            // this tell us nothing, so it doesn't help so ignore it.
            //
            for (row = 0; row < wSize - 2; row++) {
                for (col = row + 2; col < wSize; col++) {
                    score +=
                        fabs(d1[iA + row][iA + col] - d2[iB + row][iB + col]);
                }
            }

            S[iA][iB] = score / sumSize;
        }
    }
    return S;
}

pcePoint
getCoords(PyObject *L, int length) {
    // make space for the current coords
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
PyObject *
findPath(double **S, double **dA, double **dB, int lenA, int lenB,
                   int winSize, int gapMax) {

    const double D0 = 3.0;
    const double D1 = 4.0;
    int i, j;

    // Score of the best Path
    double bestPathScore = 1e6;
    int bestPathLength = 0;

    // Length of longest possible alignment
    int smaller = (lenA < lenB) ? lenA : lenB;
    int winSum = (winSize - 1) * (winSize - 2) / 2;

    path bestPath = (path)malloc(sizeof(afp) * smaller);
    for (i = 0; i < smaller; i++) {
        bestPath[i].first = -1;
        bestPath[i].second = -1;
    }

    //======================================================================
    // for storing the best N paths
    int bufferIndex = 0;
    int bufferSize = 0;
    int lenBuffer[MAX_PATHS];
    double scoreBuffer[MAX_PATHS];
    pathCache pathBuffer = (pathCache)malloc(sizeof(path *) * MAX_PATHS);

    for (i = 0; i < MAX_PATHS; i++) {
        // initialize the paths
        scoreBuffer[i] = 1e6;
        lenBuffer[i] = 0;
        pathBuffer[i] = 0;
    }

    // winCache
    // this array stores a list of residues seen.  We use it to calculate the
    // total score of a path from 1..M and then add it to M+1..N.
    int *winCache = (int *)malloc(sizeof(int) * smaller);
    for (i = 0; i < smaller; i++)
        winCache[i] = (i + 1) * i * winSize / 2 + (i + 1) * winSum;

    // allScoreBuffer
    // this 2D array keeps track of all partial gapped scores
    double **allScoreBuffer = (double **)malloc(sizeof(double *) * smaller);
    for (i = 0; i < smaller; i++) {
        allScoreBuffer[i] =
            (double *)malloc((gapMax * 2 + 1) * sizeof(double));
        // initialize the ASB
        for (j = 0; j < gapMax * 2 + 1; j++)
            allScoreBuffer[i][j] = 1e6;
    }

    int *tIndex = (int *)malloc(sizeof(int) * smaller);
    int gapBestIndex = -1;

    //======================================================================
    // Start the search through the CE matrix.
    //
    int iA, iB;
    for (iA = 0; iA < lenA; iA++) {
        if (iA > lenA - winSize * (bestPathLength - 1))
            break;

        for (iB = 0; iB < lenB; iB++) {
            if (S[iA][iB] >= D0)
                continue;

            if (S[iA][iB] == -1.0)
                continue;

            if (iB > lenB - winSize * (bestPathLength - 1))
                break;

            //
            // Restart curPath here.
            //
            path curPath = (path)malloc(sizeof(afp) * smaller);
            for (i = 0; i < smaller; i++) {
                curPath[i].first = -1;
                curPath[i].second = -1;
            }

            curPath[0].first = iA;
            curPath[0].second = iB;
            int curPathLength = 1;
            tIndex[curPathLength - 1] = 0;
            double curTotalScore = 0.0;

            //
            // Check all possible paths starting from iA, iB
            //
            int done = 0;
            while (!done) {
                double gapBestScore = 1e6;
                gapBestIndex = -1;

                //
                // Check all possible gaps [1..gapMax] from here
                //
                for (int g = 0; g < (gapMax * 2) + 1; g++) {
                    int jA = curPath[curPathLength - 1].first + winSize;
                    int jB = curPath[curPathLength - 1].second + winSize;

                    if ((g + 1) % 2 == 0) {
                        jA += (g + 1) / 2;
                    } else { // ( g odd )
                        jB += (g + 1) / 2;
                    }

                    //
                    // Following are three heuristics to ensure high quality
                    // long paths and make sure we don't run over the end of
                    // the S, matrix.

                    // 1st: If jA and jB are at the end of the matrix
                    if (jA > lenA - winSize || jB > lenB - winSize) {
                        // FIXME, was: jA > lenA-winSize-1 || jB >
                        // lenB-winSize-1
                        continue;
                    }
                    // 2nd: If this gapped octapeptide is bad, ignore it.
                    if (S[jA][jB] > D0)
                        continue;
                    // 3rd: if too close to end, ignore it.
                    if (S[jA][jB] == -1.0)
                        continue;

                    double curScore = 0.0;
                    for (int s = 0; s < curPathLength; s++) {
                        curScore += fabs(dA[curPath[s].first][jA] -
                                         dB[curPath[s].second][jB]);
                        curScore += fabs(dA[curPath[s].first + (winSize - 1)]
                                           [jA + (winSize - 1)] -
                                         dB[curPath[s].second + (winSize - 1)]
                                           [jB + (winSize - 1)]);
                        for (int k = 1; k < winSize - 1; k++)
                            curScore += fabs(dA[curPath[s].first + k]
                                               [jA + (winSize - 1) - k] -
                                             dB[curPath[s].second + k]
                                               [jB + (winSize - 1) - k]);
                    }

                    curScore /= (double)winSize * (double)curPathLength;

                    if (curScore >= D1) {
                        continue;
                    }

                    // store GAPPED best
                    if (curScore < gapBestScore) {
                        curPath[curPathLength].first = jA;
                        curPath[curPathLength].second = jB;
                        gapBestScore = curScore;
                        gapBestIndex = g;
                        allScoreBuffer[curPathLength - 1][g] = curScore;
                    }
                } /// ROF -- END GAP SEARCHING

                //
                // DONE GAPPING:
                //

                // calculate curTotalScore
                curTotalScore = 0.0;
                int jGap, gA, gB;
                double score1 = 0.0, score2 = 0.0;

                if (gapBestIndex != -1) {
                    jGap = (gapBestIndex + 1) / 2;
                    if ((gapBestIndex + 1) % 2 == 0) {
                        gA = curPath[curPathLength - 1].first + winSize + jGap;
                        gB = curPath[curPathLength - 1].second + winSize;
                    } else {
                        gA = curPath[curPathLength - 1].first + winSize;
                        gB =
                            curPath[curPathLength - 1].second + winSize + jGap;
                    }

                    // perfect
                    score1 = (allScoreBuffer[curPathLength - 1][gapBestIndex] *
                                  winSize * curPathLength +
                              S[gA][gB] * winSum) /
                             (winSize * curPathLength + winSum);

                    // perfect
                    score2 =
                        ((curPathLength > 1
                              ? (allScoreBuffer[curPathLength - 2]
                                               [tIndex[curPathLength - 1]])
                              : S[iA][iB]) *
                             winCache[curPathLength - 1] +
                         score1 * (winCache[curPathLength] -
                                   winCache[curPathLength - 1])) /
                        winCache[curPathLength];

                    curTotalScore = score2;

                    // heuristic -- path is getting sloppy, stop looking
                    if (curTotalScore > D1) {
                        done = 1;
                        gapBestIndex = -1;
                        break;
                    } else {
                        allScoreBuffer[curPathLength - 1][gapBestIndex] =
                            curTotalScore;
                        tIndex[curPathLength] = gapBestIndex;
                        curPathLength++;
                    }
                } else {
                    // if here, then there was no good gapped path
                    // so quit and restart from iA, iB+1
                    done = 1;
                    curPathLength--;
                    break;
                }

                //
                // test this gapped path against the best seen
                // starting from iA, iB
                //

                // if our currently best gapped path from iA and iB is LONGER
                // than the current best; or, it's equal length and the score's
                // better, keep the new path.
                if (curPathLength > bestPathLength ||
                    (curPathLength == bestPathLength &&
                     curTotalScore < bestPathScore)) {
                    bestPathLength = curPathLength;
                    bestPathScore = curTotalScore;
                    // deep copy curPath
                    path tempPath = (path)malloc(sizeof(afp) * smaller);

                    for (int i = 0; i < smaller; i++) {
                        tempPath[i].first = curPath[i].first;
                        tempPath[i].second = curPath[i].second;
                    }

                    free(bestPath);
                    bestPath = tempPath;
                }
            } /// END WHILE

            //
            // At this point, we've found the best path starting at iA, iB.
            //
            if (bestPathLength > lenBuffer[bufferIndex] ||
                (bestPathLength == lenBuffer[bufferIndex] &&
                 bestPathScore < scoreBuffer[bufferIndex])) {
                // we're going to add an entry to the ring-buffer.
                // Adjust maxSize values and curIndex accordingly.
                bufferIndex =
                    (bufferIndex == MAX_PATHS - 1) ? 0 : bufferIndex + 1;
                bufferSize =
                    (bufferSize < MAX_PATHS) ? (bufferSize) + 1 : MAX_PATHS;
                path pathCopy = (path)malloc(sizeof(afp) * smaller);

                int i;
                for (i = 0; i < smaller; i++) {
                    pathCopy[i].first = bestPath[i].first;
                    pathCopy[i].second = bestPath[i].second;
                }

                if (bufferIndex == 0 && (bufferSize) == MAX_PATHS) {
                    if (pathBuffer[MAX_PATHS - 1])
                        free(pathBuffer[MAX_PATHS - 1]);
                    pathBuffer[MAX_PATHS - 1] = pathCopy;
                    scoreBuffer[MAX_PATHS - 1] = bestPathScore;
                    lenBuffer[MAX_PATHS - 1] = bestPathLength;
                } else {
                    if (pathBuffer[bufferIndex - 1])
                        free(pathBuffer[bufferIndex - 1]);
                    pathBuffer[bufferIndex - 1] = pathCopy;
                    scoreBuffer[bufferIndex - 1] = bestPathScore;
                    lenBuffer[bufferIndex - 1] = bestPathLength;
                }
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

    PyObject *result = PyList_New(MAX_PATHS); // List to store all paths
    Py_INCREF(result);

    for (int o = 0; o < bufferSize; o++) {

        // Make a new list to store this path
        PyObject *pathAList = PyList_New(0);
        PyObject *pathBList = PyList_New(0);
        Py_INCREF(pathAList);
        Py_INCREF(pathBList);

        int j = 0;
        // Grab the current path
        while (j < smaller) {
            if (pathBuffer[o][j].first != -1) {
                int idxA = pathBuffer[o][j].first;
                int idxB = pathBuffer[o][j].second;

                for (int k = 0; k < winSize; k++) {
                    PyObject *v = Py_BuildValue("i", idxA + k);
                    PyList_Append(pathAList, v);
                    Py_DECREF(v);
                    v = Py_BuildValue("i", idxB + k);
                    PyList_Append(pathBList, v);
                    Py_DECREF(v);
                }
                j++;
            } else {
                break;
            }
        }
        PyObject *pairList = Py_BuildValue("[NN]", pathAList, pathBList);
        Py_INCREF(pairList);
        PyList_SET_ITEM(result, o, pairList);
    }

    // free memory
    for (i = 0; i < smaller; i++)
        free(allScoreBuffer[i]);
    free(allScoreBuffer);
    free(tIndex);
    free(winCache);
    free(bestPath);

    free(pathBuffer);

    return result;
}

// Main Function
PyObject *
PyCealign(PyObject *Py_UNUSED(self), PyObject *args) {

    int i = 0;
    int windowSize = 8;
    int gapMax = 30;
    double **dmA, **dmB, **S;

    PyObject *listA, *listB, *result;

    /* Unpack the arguments from Python */
    PyArg_ParseTuple(args, "OO|ii", &listA, &listB, &windowSize, &gapMax);

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
    S = (double **)calcS(dmA, dmB, lenA, lenB, windowSize);

    // Calculate Top N Paths
    result = (PyObject *)findPath(S, dmA, dmB, lenA, lenB, windowSize, gapMax);

    /* release memory */
    free(coordsA);
    free(coordsB);

    /* distance matrices	 */
    for (i = 0; i < lenA; i++)
        free(dmA[i]);
    free(dmA);

    for (i = 0; i < lenB; i++)
        free(dmB[i]);
    free(dmB);

    /* similarity matrix */
    for (i = 0; i < lenA; i++)
        free(S[i]);
    free(S);

    return result;
}

//
// Python Interface
//
PyDoc_STRVAR(method_doc,
"run_cealign(coordsA, coordsB, windowSize, gapMax) -> list\
\n\n\
Find the optimal alignments between two structures, using CEAlign.\
\n\n\
Arguments:\n\
- listA: List of lists with coordinates for structure A.\n\
- listB: List of lists with coordinates for structure B.\n\
- windowSize: Length of fragments to be used in alignment.\n\
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

PyObject *PyInit_ccealign(void) {
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
