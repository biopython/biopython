/**
 * The file contains implementation of the fast rmsd calculation method. The
 * sources is a marginally modified version of the implementation available
 * at theobald.brandeis.edu/qcp/ (qcprot.c). 
 *
 *********************************************************************************
 * NOTE: following changes have been made to the original (qcprot.c):
 *  - introduced python bindings: change in method signature and return types
 *  - removed methods that are not used by this module (only one method retained)
 *********************************************************************************
 *
 *  Copyright (c) 2009-2013 Pu Liu and Douglas L. Theobald
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without modification, are permitted
 *  provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this list of
 *    conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice, this list
 *    of conditions and the following disclaimer in the documentation and/or other materials
 *    provided with the distribution.
 *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to
 *    endorse or promote products derived from this software without specific prior written
 *    permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 *
 */
#include <math.h>
#include <Python.h>

static PyObject* py_FastCalcRMSDAndRotation(PyObject* self, PyObject* args) {
	double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
	double E0;
	double len;
	double minScore;
	double Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2, SyzSzymSyySzz2,
			Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2, SxzpSzx, SyzpSzy, SxypSyx,
			SyzmSzy, SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
	double C[4], rot[9];
	int i;
	double mxEigenV, rmsd;
	double oldg = 0.0;
	double b, a, delta, rms, qsqr;
	double q1, q2, q3, q4, normq;
	double a11, a12, a13, a14, a21, a22, a23, a24;
	double a31, a32, a33, a34, a41, a42, a43, a44;
	double a2, x2, y2, z2;
	double xy, az, zx, ay, yz, ax;
	double a3344_4334, a3244_4234, a3243_4233, a3143_4133, a3144_4134,
			a3142_4132;
	double evecprec = 1e-6;
	double evalprec = 1e-11;

	/* parse the arguments  */
	PyArg_ParseTuple(args, "dddddddddddd", &Sxx, &Sxy, &Sxz, &Syx, &Syy, &Syz, &Szx, &Szy, &Szz, &E0, &len, &minScore);

	Sxx2 = Sxx * Sxx;
	Syy2 = Syy * Syy;
	Szz2 = Szz * Szz;

	Sxy2 = Sxy * Sxy;
	Syz2 = Syz * Syz;
	Sxz2 = Sxz * Sxz;

	Syx2 = Syx * Syx;
	Szy2 = Szy * Szy;
	Szx2 = Szx * Szx;

	SyzSzymSyySzz2 = 2.0 * (Syz * Szy - Syy * Szz);
	Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

	C[2] = -2.0
			* (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
	C[1] = 8.0
			* (Sxx * Syz * Szy + Syy * Szx * Sxz + Szz * Sxy * Syx
					- Sxx * Syy * Szz - Syz * Szx * Sxy - Szy * Syx * Sxz);

	SxzpSzx = Sxz + Szx;
	SyzpSzy = Syz + Szy;
	SxypSyx = Sxy + Syx;
	SyzmSzy = Syz - Szy;
	SxzmSzx = Sxz - Szx;
	SxymSyx = Sxy - Syx;
	SxxpSyy = Sxx + Syy;
	SxxmSyy = Sxx - Syy;
	Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

	C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
			+ (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2)
					* (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
			+ (-(SxzpSzx) * (SyzmSzy) + (SxymSyx) * (SxxmSyy - Szz))
					* (-(SxzmSzx) * (SyzpSzy) + (SxymSyx) * (SxxmSyy + Szz))
			+ (-(SxzpSzx) * (SyzpSzy) - (SxypSyx) * (SxxpSyy - Szz))
					* (-(SxzmSzx) * (SyzmSzy) - (SxypSyx) * (SxxpSyy + Szz))
			+ (+(SxypSyx) * (SyzpSzy) + (SxzpSzx) * (SxxmSyy + Szz))
					* (-(SxymSyx) * (SyzmSzy) + (SxzpSzx) * (SxxpSyy + Szz))
			+ (+(SxypSyx) * (SyzmSzy) + (SxzmSzx) * (SxxmSyy - Szz))
					* (-(SxymSyx) * (SyzpSzy) + (SxzmSzx) * (SxxpSyy - Szz));

	mxEigenV = E0;
	for (i = 0; i < 50; ++i) {
		oldg = mxEigenV;
		x2 = mxEigenV * mxEigenV;
		b = (x2 + C[2]) * mxEigenV;
		a = b + C[1];
		delta = ((a * mxEigenV + C[0]) / (2.0 * x2 * mxEigenV + b + a));
		mxEigenV -= delta;
		/* printf("\n diff[%3d]: %16g %16g %16g", i, mxEigenV - oldg, evalprec*mxEigenV, mxEigenV); */
		if (fabs(mxEigenV - oldg) < fabs(evalprec * mxEigenV))
			break;
	}

	if (i == 50)
		PySys_WriteStderr("\nMore than %d iterations needed!\n", i);

	/* the fabs() is to guard against extremely small, but *negative* numbers due to floating point error */
	rms = sqrt(fabs(2.0 * (E0 - mxEigenV) / len));
	rmsd = rms;
	/* printf("\n\n %16g %16g %16g \n", rms, E0, 2.0 * (E0 - mxEigenV)/len); */

	if (minScore > 0) {
		if (rms < minScore)  {
			return Py_BuildValue("dddddddddddddd", -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
		}
	}

	a11 = SxxpSyy + Szz - mxEigenV;
	a12 = SyzmSzy;
	a13 = -SxzmSzx;
	a14 = SxymSyx;
	a21 = SyzmSzy;
	a22 = SxxmSyy - Szz - mxEigenV;
	a23 = SxypSyx;
	a24 = SxzpSzx;
	a31 = a13;
	a32 = a23;
	a33 = Syy - Sxx - Szz - mxEigenV;
	a34 = SyzpSzy;
	a41 = a14;
	a42 = a24;
	a43 = a34;
	a44 = Szz - SxxpSyy - mxEigenV;
	a3344_4334 = a33 * a44 - a43 * a34;
	a3244_4234 = a32 * a44 - a42 * a34;
	a3243_4233 = a32 * a43 - a42 * a33;
	a3143_4133 = a31 * a43 - a41 * a33;
	a3144_4134 = a31 * a44 - a41 * a34;
	a3142_4132 = a31 * a42 - a41 * a32;
	q1 = a22 * a3344_4334 - a23 * a3244_4234 + a24 * a3243_4233;
	q2 = -a21 * a3344_4334 + a23 * a3144_4134 - a24 * a3143_4133;
	q3 = a21 * a3244_4234 - a22 * a3144_4134 + a24 * a3142_4132;
	q4 = -a21 * a3243_4233 + a22 * a3143_4133 - a23 * a3142_4132;
	
	qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

	if (qsqr < evecprec) {
		q1 = a12 * a3344_4334 - a13 * a3244_4234 + a14 * a3243_4233;
		q2 = -a11 * a3344_4334 + a13 * a3144_4134 - a14 * a3143_4133;
		q3 = a11 * a3244_4234 - a12 * a3144_4134 + a14 * a3142_4132;
		q4 = -a11 * a3243_4233 + a12 * a3143_4133 - a13 * a3142_4132;
		qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

		if (qsqr < evecprec) {
			double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24
					- a14 * a22;
			double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24
					- a14 * a21;
			double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22
					- a12 * a21;

			q1 = a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
			q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
			q3 = a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
			q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
			qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

			if (qsqr < evecprec) {
				q1 = a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
				q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
				q3 = a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
				q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
				qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

				if (qsqr < evecprec) {
					rot[0] = rot[4] = rot[8] = 1.0;
					rot[1] = rot[2] = rot[3] = rot[5] = rot[6] = rot[7] = 0.0;

					return Py_BuildValue("dddddddddddddd", rmsd, rot[0], rot[1], rot[2], rot[3], rot[4], rot[5], rot[6], rot[7], rot[8], q1, q2, q3, q4);
				}
			}
		}
	}

	normq = sqrt(qsqr);
	q1 /= normq;
	q2 /= normq;
	q3 /= normq;
	q4 /= normq;

	a2 = q1 * q1;
	x2 = q2 * q2;
	y2 = q3 * q3;
	z2 = q4 * q4;

	xy = q2 * q3;
	az = q1 * q4;
	zx = q4 * q2;
	ay = q1 * q3;
	yz = q3 * q4;
	ax = q1 * q2;

	rot[0] = a2 + x2 - y2 - z2;
	rot[1] = 2 * (xy + az);
	rot[2] = 2 * (zx - ay);
	rot[3] = 2 * (xy - az);
	rot[4] = a2 - x2 + y2 - z2;
	rot[5] = 2 * (yz + ax);
	rot[6] = 2 * (zx + ay);
	rot[7] = 2 * (yz - ax);
	rot[8] = a2 - x2 - y2 + z2;

	return Py_BuildValue("dddddddddddddd", rmsd, rot[0], rot[1], rot[2], rot[3], rot[4], rot[5], rot[6], rot[7], rot[8], q1, q2, q3, q4);
}

static PyMethodDef qcprot_methods[] = {
        {"FastCalcRMSDAndRotation", (PyCFunction)py_FastCalcRMSDAndRotation, METH_VARARGS, "The method calculates the RMSD by solving for the most positive eigenvalue using the Newton-Raphson method. The rotation matrix is given by the corresponding eigenvector and is calculated by finding roots of the characteristic polynomial of the matrix. The method returns the rmsd, the rotation matrix and the 4 quaternions."},
        {NULL, NULL, 0, NULL} 
};


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT, "qcprotmodule", NULL,
        -1, qcprot_methods, NULL, NULL, NULL, NULL};


PyObject * PyInit_qcprotmodule(void) {
        return PyModule_Create(&moduledef);
}
