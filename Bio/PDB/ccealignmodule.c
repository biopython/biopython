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
    int pA;
    int pB;
} afp, *path;

// Calculate distance matrix
static double **
calcDM(pcePoint coords, int len)
{
    double **dm = (double **)PyMem_RawMalloc(sizeof(double *) * len);

    for (int i = 0; i < len; i++) {
        dm[i] = (double *)PyMem_RawMalloc(sizeof(double) * len);
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
// This similarity corresponds to distance measure (i) or
// equation (6) in the paper.
//
static double
similarityI(
    double **dA,
    double **dB,
    const afp afpI,
    const afp afpJ,
    const int fragmentSize)
{
    const int iA = afpI.pA;
    const int iB = afpI.pB;
    const int jA = afpJ.pA;
    const int jB = afpJ.pB;
    const int m = fragmentSize;
    double similarity = fabs(dA[iA][jA] - dB[iB][jB]) +
        fabs(dA[iA + m-1][jA + m-1] - dB[iB + m-1][jB + m-1]);

    for (int k = 1; k < m - 1; k++) {
        similarity += fabs(dA[iA + k][jA + m-1 - k] -
            fabs(dB[iB + k][jB + m-1 - k]));
    }

    return -similarity / m;
}

//
// This similarity corresponds to distance measure (ii)
// or equation (7) in the paper.
//
static double
similarityII(
    double **dA,
    double **dB,
    const afp afpI,
    const int fragmentSize)
{
    const int iA = afpI.pA;
    const int iB = afpI.pB;
    double similarity = 0.0;
    // Term count is the number of terms in the summation
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
calcS(
    double **dA,
    double **dB,
    const int lenA,
    const int lenB,
    const int fragmentSize)
{
    // Initialize the 2D similarity matrix
    const int rowCount = lenA - fragmentSize + 1;
    const int colCount = lenB - fragmentSize + 1;
    double **S = (double **)PyMem_RawMalloc(sizeof(double *) * rowCount);

    for (int i = 0; i < rowCount; i++) {
        S[i] = (double *)PyMem_RawMalloc(sizeof(double) * colCount);
    }

    //
    // This is where the magic of CE comes out. In the similarity matrix,
    // for each i and j, the value of S[i][j] is how well the fragment starting
    // at i in protein A matches the fragment starting at j in protein
    // B. A value of 0 means absolute match; a value << -3 means bad match.
    //
    for (int iA = 0; iA < rowCount; iA++) {
        for (int iB = 0; iB < colCount; iB++) {
            S[iA][iB] = similarityII(dA, dB, (afp) {iA, iB}, fragmentSize);
        }
    }

    return S;
}

static const double tableZtoP[] = {
    1.0, 9.20e-01, 8.41e-01, 7.64e-01, 6.89e-01, 6.17e-01, 5.49e-01, 4.84e-01, 4.24e-01, 3.68e-01,
    3.17e-01, 2.71e-01, 2.30e-01, 1.94e-01, 1.62e-01, 1.34e-01, 1.10e-01, 8.91e-02, 7.19e-02, 5.74e-02,
    4.55e-02, 3.57e-02, 2.78e-02, 2.14e-02, 1.64e-02, 1.24e-02, 9.32e-03, 6.93e-03, 5.11e-03, 3.73e-03,
    2.70e-03, 1.94e-03, 1.37e-03, 9.67e-04, 6.74e-04, 4.65e-04, 3.18e-04, 2.16e-04, 1.45e-04, 9.62e-05,
    6.33e-05, 4.13e-05, 2.67e-05, 1.71e-05, 1.08e-05, 6.80e-06, 4.22e-06, 2.60e-06, 1.59e-06, 9.58e-07,
    5.73e-07, 3.40e-07, 1.99e-07, 1.16e-07, 6.66e-08, 3.80e-08, 2.14e-08, 1.20e-08, 6.63e-09, 3.64e-09,
    1.97e-09, 1.06e-09, 5.65e-10, 2.98e-10, 1.55e-10, 8.03e-11, 4.11e-11, 2.08e-11, 1.05e-11, 5.20e-12,
    2.56e-12, 1.25e-12, 6.02e-13, 2.88e-13, 1.36e-13, 6.38e-14, 2.96e-14, 1.36e-14, 6.19e-15, 2.79e-15,
    1.24e-15, 5.50e-16, 2.40e-16, 1.04e-16, 4.46e-17, 1.90e-17, 7.97e-18, 3.32e-18, 1.37e-18, 5.58e-19,
    2.26e-19, 9.03e-20, 3.58e-20, 1.40e-20, 5.46e-21, 2.10e-21, 7.99e-22, 3.02e-22, 1.13e-22, 4.16e-23,
    1.52e-23, 5.52e-24, 1.98e-24, 7.05e-25, 2.48e-25, 8.64e-26, 2.98e-26, 1.02e-26, 3.44e-27, 1.15e-27,
    3.82e-28, 1.25e-28, 4.08e-29, 1.31e-29, 4.18e-30, 1.32e-30, 4.12e-31, 1.27e-31, 3.90e-32, 1.18e-32,
    3.55e-33, 1.06e-33, 3.11e-34, 9.06e-35, 2.61e-35, 7.47e-36, 2.11e-36, 5.91e-37, 1.64e-37, 4.50e-38,
    1.22e-38, 3.29e-39, 8.77e-40, 2.31e-40, 6.05e-41, 1.56e-41, 4.00e-42, 1.02e-42, 2.55e-43, 6.33e-44,
    1.56e-44, 3.80e-45, 9.16e-46, 2.19e-46, 5.17e-47, 1.21e-47, 2.81e-48, 6.45e-49, 1.46e-49, 3.30e-50};

static const double tablePtoZ[] = {
    0.00, 0.73, 1.24, 1.64, 1.99, 2.30, 2.58, 2.83, 3.07, 3.29,
    3.50, 3.70, 3.89, 4.07, 4.25, 4.42, 4.58, 4.74, 4.89, 5.04,
    5.19, 5.33, 5.46, 5.60, 5.73, 5.86, 5.99, 6.11, 6.23, 6.35,
    6.47, 6.58, 6.70, 6.81, 6.92, 7.02, 7.13, 7.24, 7.34, 7.44,
    7.54, 7.64, 7.74, 7.84, 7.93, 8.03, 8.12, 8.21, 8.30, 8.40,
    8.49, 8.57, 8.66, 8.75, 8.84, 8.92, 9.01, 9.09, 9.17, 9.25,
    9.34, 9.42, 9.50, 9.58, 9.66, 9.73, 9.81, 9.89, 9.97, 10.04,
    10.12, 10.19, 10.27, 10.34, 10.41, 10.49, 10.56, 10.63, 10.70, 10.77,
    10.84, 10.91, 10.98, 11.05, 11.12, 11.19, 11.26, 11.32, 11.39, 11.46,
    11.52, 11.59, 11.66, 11.72, 11.79, 11.85, 11.91, 11.98, 12.04, 12.10,
    12.17, 12.23, 12.29, 12.35, 12.42, 12.48, 12.54, 12.60, 12.66, 12.72,
    12.78, 12.84, 12.90, 12.96, 13.02, 13.07, 13.13, 13.19, 13.25, 13.31,
    13.36, 13.42, 13.48, 13.53, 13.59, 13.65, 13.70, 13.76, 13.81, 13.87,
    13.92, 13.98, 14.03, 14.09, 14.14, 14.19, 14.25, 14.30, 14.35, 14.41,
    14.46, 14.51, 14.57, 14.62, 14.67, 14.72, 14.77, 14.83, 14.88, 14.93};

// Convert a z-score into a probability
static double zToP(const double z)
{
    int index = (int)(z / 0.1);

    if (index < 0) {
        index = 0;
    }
    if (index > 149) {
        index = 149;
    }

    return tableZtoP[index];
}

// Convert a probability into a z-score
static double pToZ(const double p)
{
    int index = (int)(-log10(p) * 3.0);

    if (index < 0) {
        index = 0;
    }
    if (index > 149) {
        index = 149;
    }

    return tablePtoZ[index];
}

// These empirical data are reproduced from the original CE source code.
static const double similarityAvgs[] =
    {2.54, 2.51, 2.72, 3.01, 3.31, 3.61, 3.90, 4.19, 4.47, 4.74,
        4.99, 5.22, 5.46, 5.70, 5.94, 6.13, 6.36, 6.52, 6.68, 6.91};
static const double similaritySDs[] =
    {1.33, 0.88, 0.73, 0.71, 0.74, 0.80, 0.86, 0.92, 0.98, 1.04,
        1.08, 1.10, 1.15, 1.19, 1.23, 1.25, 1.32, 1.34, 1.36, 1.45};

static double zScoreSimilarity(
    const int pathLength,
    const double similarity)
{
    // This method only works for the default fragment size
    if (pathLength < 1) {
        return 0.0;
    }

    double similarityAvg, similaritySD;

    // 20 is the number of stored statistics (averages and standard deviations)
    // in the arrays above.
    if (pathLength <= 20) {
        similarityAvg = similarityAvgs[pathLength - 1];
        similaritySD = similaritySDs[pathLength - 1];
    }
    else {
        similarityAvg = 0.209874 * pathLength + 2.944714;
        similaritySD = 0.039487 * pathLength + 0.675735;
    }
    if (similarity > similarityAvg) {
        return 0.0;
    }

    return (similarityAvg - similarity) / similaritySD;
}

static const double gapCountAvgs[] =
    {0.00, 11.50, 23.32, 35.95, 49.02, 62.44, 76.28, 90.26,
        104.86, 119.97, 134.86, 150.54, 164.86, 179.57, 194.39,
        209.38, 224.74, 238.96, 253.72, 270.79};
static const double gapCountSDs[] =
    {0.00, 9.88, 14.34, 17.99, 21.10, 23.89, 26.55, 29.00, 31.11,
        33.10, 35.02, 36.03, 37.19, 38.82, 41.04, 43.35, 45.45,
        48.41, 50.87, 52.27};

static double zScoreGapCount(
    const int pathLength,
    const int gapCount)
{
    if (pathLength < 1) {
        return 0.0;
    }

    double gapCountAvg, gapCountSD;

    // 20 is the number of stored statistics (averages and standard deviations)
    // in the arrays above.
    if (pathLength <= 20) {
        gapCountAvg = gapCountAvgs[pathLength - 1];
        gapCountSD = gapCountSDs[pathLength - 1];
    }
    else {
        gapCountAvg = 14.949173 * pathLength - 14.581193;
        gapCountSD = 2.045067 * pathLength + 13.191095;
    }
    if (gapCount > gapCountAvg) {
        return 0.0;
    }

    return (gapCountAvg - gapCount) / gapCountSD;
}

// The z-score calculation is adapted from the code in
// https://github.com/kad-ecoli/CE.
static double calcZScore(
    const int fragmentSize,
    const int pathLength,
    const double pathSimilarity,
    const int gapCount)
{
    if (fragmentSize != 8) {
        // Z-score calculation is only supported for the default fragment size.
        return 0.0;
    }

    const double z1 = zScoreSimilarity(pathLength, pathSimilarity);
    const double z2 = zScoreGapCount(pathLength, gapCount);

    return pToZ(zToP(z1) * zToP(z2));
}

static pcePoint
getCoords(PyObject *L, int length)
{
    // Make space for the current coords
    pcePoint coords = (pcePoint)PyMem_RawMalloc(sizeof(cePoint) * length);

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
    int bufferSize = 0;
    int lenBuffer[MAX_PATHS];
    double similarityBuffer[MAX_PATHS];
    path pathBuffer[MAX_PATHS];

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
        if (bufferSize > 0 &&
            iA > lenA - fragmentSize * (lenBuffer[bufferSize - 1] - 1))
            break;

        for (int iB = 0; iB <= lenB - fragmentSize; iB++) {
            if (S[iA][iB] <= D0)
                continue;
            if (bufferSize > 0 &&
                iB > lenB - fragmentSize * (lenBuffer[bufferSize - 1] - 1))
                break;

            // Initialize current path
            path curPath = (path)PyMem_RawMalloc(sizeof(afp) * smaller);
            int curPathLength = 1;
            double curPathSimilarity = S[iA][iB];

            curPath[0] = (afp) {iA, iB};

            for (int i = 1; i < smaller; i++) {
                curPath[i] = (afp) {-1, -1};
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
                    int jA = curPath[curPathLength - 1].pA + fragmentSize;
                    int jB = curPath[curPathLength - 1].pB + fragmentSize;

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

                    const afp afpJ = (afp) {jA, jB};
                    double curSimilarity = 0.0;

                    for (int s = 0; s < curPathLength; s++) {
                        curSimilarity +=
                            similarityI(
                                dA,
                                dB,
                                curPath[s],
                                afpJ,
                                fragmentSize);
                    }

                    curSimilarity /= curPathLength;

                    // store GAPPED best
                    if (curSimilarity > D1 &&
                        curSimilarity > gapBestSimilarity) {
                        curPath[curPathLength] = afpJ;
                        gapBestSimilarity = curSimilarity;
                        gapBestIndex = g;
                    }
                } /// ROF -- END GAP SEARCHING

                if (gapBestIndex == -1) {
                    // if here, then there was no good candidate AFP,
                    // so quit and restart from starting point
                    break;
                }

                // The current path has n AFPs, and we are considering adding
                // the (n+1)-th AFP.
                // Imagine a matrix where entry ij is D_ij of the i-th and j-th
                // AFPs in the path.
                // The path similarity is the average of the upper triangle of
                // this matrix.
                const afp afpJ = curPath[curPathLength];
                const double n = (double) curPathLength;
                const double curTermCount = n + n * (n - 1) / 2;
                const double newTermCount = n + 1 + n * (n + 1) / 2;
                // Notice that the new term count is
                // the current term count plus n + 1.
                const double newSimilarity =
                        (curTermCount * curPathSimilarity +
                        n * gapBestSimilarity +
                        S[afpJ.pA][afpJ.pB]) / newTermCount;

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
            for (int i = 0; i < bufferSize; i++) {
                if (curPathLength > lenBuffer[i] ||
                    (curPathLength == lenBuffer[i] &&
                     curPathSimilarity > similarityBuffer[i])) {
                     // Swap the current path with the path in the buffer
                     int tempLength = lenBuffer[i];
                     double tempSimilarity = similarityBuffer[i];
                     path tempPath = pathBuffer[i];

                     lenBuffer[i] = curPathLength;
                     similarityBuffer[i] = curPathSimilarity;
                     pathBuffer[i] = curPath;

                     curPathLength = tempLength;
                     curPathSimilarity = tempSimilarity;
                     curPath = tempPath;
                }
            }

            if (bufferSize < MAX_PATHS) {
                lenBuffer[bufferSize] = curPathLength;
                similarityBuffer[bufferSize] = curPathSimilarity;
                pathBuffer[bufferSize] = curPath;
                bufferSize += 1;
            }
            else {
                PyMem_RawFree(curPath);
            }
        } // ROF -- end for iB
    }     // ROF -- end for iA

    double zScoreBuffer[MAX_PATHS];

    for (int i = 0; i < bufferSize; i++) {
        const int pathLength = lenBuffer[i];
        const double pathSimilarity = similarityBuffer[i];
        int gapCount = 0;

        for (int j = 1; j < pathLength; j++) {
            gapCount += pathBuffer[i][j].pA - pathBuffer[i][j - 1].pA - 1;
            gapCount += pathBuffer[i][j].pB - pathBuffer[i][j - 1].pB - 1;
        }

        zScoreBuffer[i] = calcZScore(fragmentSize, pathLength, pathSimilarity, gapCount);
    }

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

        for (int j = 0; j < lenBuffer[o]; j++) {
            const int pA = pathBuffer[o][j].pA;
            const int pB = pathBuffer[o][j].pB;

            for (int k = 0; k < fragmentSize; k++) {
                PyObject *v = Py_BuildValue("i", pA + k);
                PyList_Append(pathAList, v);
                Py_DECREF(v);
                v = Py_BuildValue("i", pB + k);
                PyList_Append(pathBList, v);
                Py_DECREF(v);
            }
        }

        const double zScore = zScoreBuffer[o];
        const int length = lenBuffer[o];
        PyObject *pairList = Py_BuildValue("[NN]", pathAList, pathBList);
        Py_INCREF(pairList);

        PyStructSequence_Field namedtupleFields[] = {
            (PyStructSequence_Field) {
                "path",
                NULL,
            },
            (PyStructSequence_Field) {
                "z_score",
                NULL,
            },
            (PyStructSequence_Field) {
                "length",
                NULL,
            },
            {NULL},
        };
        PyStructSequence_Desc namedtupleDesc = (PyStructSequence_Desc) {
            "ccealign.CEAlignment",
            NULL,
            namedtupleFields,
            3,
        };
        PyTypeObject *namedtupleType =
            PyStructSequence_NewType(&namedtupleDesc);
        PyObject *namedtuple = PyStructSequence_New(namedtupleType);

        PyStructSequence_SetItem(namedtuple, 0, pairList);
        PyStructSequence_SetItem(namedtuple, 1, PyFloat_FromDouble(zScore));
        PyStructSequence_SetItem(namedtuple, 2, PyLong_FromLong(length * fragmentSize));

        PyList_SET_ITEM(result, o, namedtuple);
        Py_DECREF(namedtupleType);
    }

    return result;
}

// Main Function
PyObject *
PyCealign(PyObject *Py_UNUSED(self), PyObject *args)
{
    int fragmentSize = 8;
    int gapMax = 30;

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
    double **dA = (double **)calcDM(coordsA, lenA);
    double **dB = (double **)calcDM(coordsB, lenB);

    /* calculate the CE Similarity matrix */
    double **S = (double **)calcS(dA, dB, lenA, lenB, fragmentSize);

    // Calculate Top N Paths
    result = (PyObject *)findPath(S, dA, dB, lenA, lenB, fragmentSize, gapMax);

    /* release memory */
    PyMem_RawFree(coordsA);
    PyMem_RawFree(coordsB);

    /* distance matrices	 */
    for (int i = 0; i < lenA; i++)
        PyMem_RawFree(dA[i]);
    PyMem_RawFree(dA);

    for (int i = 0; i < lenB; i++)
        PyMem_RawFree(dB[i]);
    PyMem_RawFree(dB);

    // Similarity matrix
    for (int i = 0; i <= lenA - fragmentSize; i++)
        PyMem_RawFree(S[i]);
    PyMem_RawFree(S);

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
