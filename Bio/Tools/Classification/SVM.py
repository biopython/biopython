# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This provides code for working with support vector machines.

Support vector machine (SVM) is a supervised machine learning
algorithm that has been shown to work well on a variety of
classification problems.

XXX describe kernel_fn, C, epsilon

For more information, see:
http://svm.first.gmd.de/
http://svm.research.bell-labs.com/

Classes:
SVM                  Holds data for a Support Vector Machine.
SMOTrainer           Trains an SVM using Sequential Minimal Optimization (SMO).
KeerthiTrainer       Trains an SVM using Keerthi's extensions to SMO.
TransductiveTrainer  Trains a transductive SVM.

LinearKernel
PolynomialKernel
RadialBasisFunctionKernel
HyperbolicTangentKernel


Functions:
train                Train a Support Vector Machine on some training data.
trans_transductive
classify             Use a Support Vector Machine to classify some data.

Usage:
The 'train' function is provided as a user-friendly interface module.
Use this function if all you want is a plain-vanilla SVM classifier.
However, if you're concerned about details such as the training
algorithm, you will need to use some of the classes provided.

"""

# To do:
# make update_fn better
#    take an event string, followed by data
#    or, make it an object that responds to various things

import math
import random

class LinearKernel:
    """k(x, y) = x*y"""
    def __init__(self):
        pass
    def __call__(self, x1, x2):
        return _dot(x1, x2)

# XXX document sparse kernels
# XXX try to generalize the kernels, so they work with either type of input
class SparseLinearKernel:
    """k(x, y) = x*y"""
    def __init__(self):
        pass
    def __call__(self, x1, x2):
        return _sparse_dot(x1, x2)

class PolynomialKernel:
    """k(x, y) = (x*y+1)**p"""
    _max_p = 20
    def __init__(self, p=3):
        if p < 1 or p > self._max_p:
            raise ValueError, \
                  "polynomial degree must be between 1 and %d" % self._max_p
        self.p = p
    def __call__(self, x1, x2):
        s = _dot(x1, x2) + 1.0
        return s ** self.p

class SparsePolynomialKernel:
    """k(x, y) = (x*y+1)**p"""
    _max_p = 20
    def __init__(self, p=3):
        if p < 1 or p > self._max_p:
            raise ValueError, \
                  "polynomial degree must be between 1 and %d" % self._max_p
        self.p = p
    def __call__(self, x1, x2):
        s = _sparse_dot(x1, x2) + 1.0
        return s ** self.p

class RadialBasisFunctionKernel:
    """      -(x-y)**2
         e** ----------
             2*sigma**2
    """
    def __init__(self, sigma=1.0):
        self.sigma = sigma
        self.sigma2 = sigma**2
    def __call__(self, x1, x2):
        x1_x2 = _subtract(x1, x2)
        num = -_dot(x1_x2, x1_x2)
        den = 2.0 * self.sigma2
        try:
            return math.exp(num/den)
        except OverflowError:
            return 0.0

class SparseRadialBasisFunctionKernel:
    """      -(x-y)**2
         e** ----------
             2*sigma**2
    """
    def __init__(self, sigma=1.0):
        self.sigma = sigma
        self.sigma2 = sigma**2
    def __call__(self, x1, x2):
        x1_x2 = _sparse_subtract(x1, x2)
        num = -_sparse_dot(x1_x2, x1_x2)
        den = 2.0 * self.sigma2
        try:
            return math.exp(num/den)
        except OverflowError:
            return 0.0

class HyperbolicTangentKernel:
    """tanh(x*y*kappa - delta)"""
    # not all kappa and delta obey Mercer's condition
    def __init__(self, kappa, delta):
        self.kappa = kappa
        self.delta = delta
    def __call__(self, x1, x2):
        x1x2 = _dot(x1, x2)
        return math.tanh(x1x2*self.kappa - self.delta)

class SVM:
    """Holds information for a non-linear Support Vector Machine.

    Members:
    xs          A list of the input vectors.
    ys          A list of the output values.  Must be either 1 or -1.
    alphas      A list of the Lagrange multipliers for each point.
    b           The threshold value for the machine.
    kernel_fn   Should take 2 vectors and return a distance.

    * xs, ys, and alphas should be a parallel list of vectors.
    ** There's a special read-only member variable 'w'.  For linear
    SVM's, this is the 'w' vector so that the decision hyperplane is at
    wx-b=0

    """
    def __init__(self, xs, ys, alphas=None, b=None, kernel_fn=LinearKernel()):
        if alphas is None:
            alphas = [0.0]*len(self.xs)
        if b is None:
            b = 0.0

        self.xs, self.ys = xs, ys
        self.alphas = alphas
        self.b = b
        self.kernel_fn = kernel_fn

        if len(xs) != len(ys):
            raise ValueError, "xs and ys must be parallel arrays"
        if len(xs) != len(alphas):
            raise ValueError, "xs and alphas must be parallel arrays"

        for i in range(len(alphas)):
            if self.ys[i] not in [1, -1]:
                raise ValueError, "Output must be either 1 or -1"

    def __getattr__(self, x):
        if x == 'w':
            w = [0.0] * len(self.xs[0])
            for i in range(len(w)):
                for j in range(len(self.xs)):
                    w[i] = w[i] + self.ys[j]*self.alphas[j]*self.xs[j][i]
            return w
        raise AttributeError, x
        
def classify(svm, x):
    """classify(svm, x) -> score

    Classify x based on an SVM object.

    """
    sum = 0.0
    for i in range(len(svm.xs)):
        if svm.alphas[i] == 0.0:
            continue
        sum = sum + svm.ys[i]*svm.alphas[i]*svm.kernel_fn(svm.xs[i], x)
    return sum - svm.b

class SMOTrainer:
    """Sequential Minimal Optimization Trainer

    This is an implementation of the SVM trainer called Sequential
    Minimal Optimization, described by John C. Platt:
    http://www.research.microsoft.com/~jplatt

    Methods:
    train     Train a new SVM.

    """
    def train(self, training_set, results, kernel_fn, C, epsilon,
              update_fn=None):
        """train(self, training_set, results, kernel_fn, C, epsilon,
        update_fn=None) -> SVM

        """
        # Initializing variables
        self.svm = SVM(training_set, results,
                       [0.0]*len(training_set), 0.0, kernel_fn)
        self.eps = epsilon
        self.C = C
        self._ecache = {}

        # Training algorithm starts here.
        num_changed = 0
        examine_all = 1

        while num_changed > 0 or examine_all:
            num_changed = 0
            if examine_all:
                for i2 in range(len(self.svm.alphas)):
                    num_changed = num_changed + self._examine_example(i2)
            else:
                for i2 in range(len(self.svm.alphas)):
                    # skip the bound alphas
                    if self._is_bound(self.svm.alphas[i2], self.C):
                        continue
                    num_changed = num_changed + self._examine_example(i2)
            if update_fn is not None:
                update_fn(num_changed, self.svm)
            if examine_all:
                examine_all = 0
            elif num_changed == 0:
                examine_all = 1

        return self.svm

    def _examine_example(self, i2):
        y2 = self.svm.ys[i2]
        alph2 = self.svm.alphas[i2]
        E2 = self._E(i2)
        r2 = E2*y2

        if (r2 < -self.eps and alph2 < self.C) or \
           (r2 > self.eps and alph2 > 0):
            # sort it by E1 - E2
            max_i1 = None
            max_error = None
            for i1 in range(len(self.svm.alphas)):
                if i1 == i2 or self._is_bound(self.svm.alphas[i1], self.C):
                    continue
                E1 = self._E(i1)
                error = math.fabs(E1 - E2)
                if max_i1 is None or error > max_error:
                    max_error = error
                    max_i1 = i1
            i1 = max_i1
            if i1 is not None:
                if self._take_step(i1, i2):
                    return 1
                
            i1 = random.randint(0, len(self.svm.alphas)-1)
            for i in range(len(self.svm.alphas)):
                if self._is_bound(self.svm.alphas[i1], self.C):
                    continue
                if self._take_step(i1, i2):
                    return 1
                i1 = (i1 + 1) % len(self.svm.alphas)

            i1 = random.randint(0, len(self.svm.alphas)-1)
            for i in range(len(self.svm.alphas)):
                if self._take_step(i1, i2):
                    return 1
                i1 = (i1 + 1) % len(self.svm.alphas)
        return 0

    def _take_step(self, i1, i2):
        if i1 == i2:
            return 0
        point = self.svm.xs
        alph1, alph2 = self.svm.alphas[i1], self.svm.alphas[i2]
        y1, y2 = self.svm.ys[i1], self.svm.ys[i2]
        E1 = self._E(i1)
        E2 = self._E(i2)
        s = y1*y2
        b = self.svm.b
        C = self.C
        eps = self.eps
        
        if y1 == y2:
            L = max(0.0, alph2+alph1-C)
            H = min(C, alph2+alph1)
        else:
            L = max(0.0, alph2-alph1)
            H = min(C, C+alph2-alph1)
        if L == H:
            return 0

        k11 = self.svm.kernel_fn(point[i1], point[i1])
        k12 = self.svm.kernel_fn(point[i1], point[i2])
        k22 = self.svm.kernel_fn(point[i2], point[i2])
        eta = k11+k22-2.0*k12
        if eta > 0.0:
            a2 = alph2 + y2*(E1-E2)/eta
            if a2 < L: a2 = L
            elif a2 > H: a2 = H
        else:
            f1 = y1*(E1+b) - alph1*k11 - s*alph2*k22
            f2 = y2*(E2+b) - s*alph1*k12 - alph2*k22
            L1 = alph1 + s*(alph2-L)
            H1 = alph1 + s*(alph2-H)
            Lobj = L1*f1 + L*f2 + 0.5*L1*L1*k11 + 0.5*L*L*k22 + s*L*L1*k12
            Hobj = H1*f1 + H*f2 + 0.5*H1*H1*k11 + 0.5*H*H*k22 + s*H*H1*k12
            
            if Lobj < Hobj-eps:
                a2 = L
            elif Lobj > Hobj+eps:
                a2 = H
            else:
                a2 = alph2

        if math.fabs(a2-alph2) < eps*(a2+alph2+eps):
            return 0

        a1 = alph1 + s*(alph2-a2)
        
        if not self._is_bound(a1, C):
            b1 = E1 + y1*(a1-alph1)*k11 + y2*(a2-alph2)*k12 + b
            b = b1
        elif not self._is_bound(a2, C):
            b2 = E2 + y1*(a1-alph1)*k12 + y2*(a2-alph2)*k22 + b
            b = b2
        else:
            b1 = E1 + y1*(a1-alph1)*k11 + y2*(a2-alph2)*k12 + b
            b2 = E2 + y1*(a1-alph1)*k12 + y2*(a2-alph2)*k22 + b
            b = (b1+b2)/2.0

        self.svm.alphas[i1], self.svm.alphas[i2] = a1, a2
        self.svm.b = b
        self._ecache = {}

        return 1

    def _E(self, i):
        if not self._ecache.has_key(i):
            self._ecache[i] = classify(self.svm, self.svm.xs[i])-self.svm.ys[i]
        return self._ecache[i]

    def _is_bound(self, alpha, C):
        return alpha <= 0 or alpha >= C

class KeerthiTrainer:
    """Sequential Minimal Optimization Trainer

    This is an implementation of Sequential Minimal Optimization plus
    the 2 modifications suggested by S.S. Keerthi, et al.:
    http://guppy.mpe.nus.edu.sg/~mpessk/

    Methods:
    train     Train a new SVM.

    """
    def train(self, training_set, results, kernel_fn, C, epsilon,
              update_fn=None):
        """train(self, training_set, results, kernel_fn, C, epsilon,
        update_fn=None) -> SVM

        """
        # Initializing variables
        self.svm = SVM(training_set, results,
                       [0.0]*len(training_set), 0.0, kernel_fn)
        self.eps = epsilon
        self.C = C
        self._fcache = {}

        # Initialize values.
        # b_low and b_up should be set to extreme values.
        # self.i_up and self.i_low should be set to a positive example,
        # and a negative example.
        self.b_low, self.b_up = 1, -1
        self.i_up = self.i_low = None
        i = 0
        while i < len(self.svm.alphas):
            if self.svm.ys[i] == 1 and self.i_up is None:
                self.i_up = i
            elif self.svm.ys[i] == -1 and self.i_low is None:
                self.i_low = i
            if self.i_up is not None and self.i_low is not None:
                break
            i = i + 1
        else:
            raise ValueError, \
                  "I could not find positive and negative training examples"

        # Training algorithm starts here.
        while 1:
            num_changed = 0
            for i2 in range(len(self.svm.alphas)):
                num_changed = num_changed + self._examine_example(i2)
            if not num_changed:
                break
            
            inner_loop_success = 1
            while (self.b_up <= self.b_low-2*self.eps) and \
                  inner_loop_success:
                i2 = self.i_low
                y2 = self.svm.ys[i2]
                alph2 = self.svm.alphas[i2]
                F2 = self._F(i2)
                i1 = self.i_up
                inner_loop_success = self._take_step(self.i_up, self.i_low)
                num_changed = num_changed + inner_loop_success
                
            if update_fn is not None:
                update_fn(num_changed, self.svm)

        # b_low and b_up define the range of allowable b's, within
        # a tolerance of epsilon.  I'll just choose the b in the middle.
        self.svm.b = (self.b_low+self.b_up)/2.0

        return self.svm

    def _examine_example(self, i2):
        y2 = self.svm.ys[i2]
        alph2 = self.svm.alphas[i2]

        F2 = self._F(i2)
        if self._in_I_0(i2):
            pass
        elif (self._in_I_1(i2) or self._in_I_2(i2)) and F2 < self.b_up:
            self.b_up = F2
            self.i_up = i2
        elif (self._in_I_3(i2) or self._in_I_4(i2)) and F2 > self.b_low:
            self.b_low = F2
            self.i_low = i2

        optimality = 1
        if self._in_I_0(i2) or self._in_I_1(i2) or self._in_I_2(i2):
            #if self.b_low-F2 > 2*self.eps:
            if self.b_low > F2:
                optimality = 0
                i1 = self.i_low
        if self._in_I_0(i2) or self._in_I_3(i2) or self._in_I_4(i2):
            #if F2-self.b_up > 2*self.eps:
            if F2 > self.b_up:
                optimality = 0
                i1 = self.i_up
        if optimality:
            return 0

        # for i2 in I_0 choose the better i1
        if self._in_I_0(i2):
            if self.b_low-F2 > F2-self.b_up:
                i1 = self.i_low
            else:
                i1 = self.i_up
        return self._take_step(i1, i2)

    def _take_step(self, i1, i2):
        if i1 == i2:
            return 0
        point = self.svm.xs
        alph1, alph2 = self.svm.alphas[i1], self.svm.alphas[i2]
        y1, y2 = self.svm.ys[i1], self.svm.ys[i2]
        F1, F2 = self._F(i1), self._F(i2)
        s = y1*y2
        b = self.svm.b
        C = self.C
        eps = self.eps
        
        if y1 == y2:
            L = max(0.0, alph2+alph1-C)
            H = min(C, alph2+alph1)
        else:
            L = max(0.0, alph2-alph1)
            H = min(C, C+alph2-alph1)
        if L == H:
            return 0

        k11 = self.svm.kernel_fn(point[i1], point[i1])
        k12 = self.svm.kernel_fn(point[i1], point[i2])
        k22 = self.svm.kernel_fn(point[i2], point[i2])
        eta = k11+k22-2.0*k12
        if eta > 0.0:
            a2 = alph2 + y2*(F1-F2)/eta
            if a2 < L: a2 = L
            elif a2 > H: a2 = H
        else:
            f1 = y1*F1 - alph1*k11 - s*alph2*k22
            f2 = y2*F2 - s*alph1*k12 - alph2*k22
            L1 = alph1 + s*(alph2-L)
            H1 = alph1 + s*(alph2-H)
            Lobj = L1*f1 + L*f2 + 0.5*L1*L1*k11 + 0.5*L*L*k22 + s*L*L1*k12
            Hobj = H1*f1 + H*f2 + 0.5*H1*H1*k11 + 0.5*H*H*k22 + s*H*H1*k12
            
            if Lobj < Hobj-eps:
                a2 = L
            elif Lobj > Hobj+eps:
                a2 = H
            else:
                a2 = alph2

        if math.fabs(a2-alph2) < eps*(a2+alph2+eps):
            return 0

        a1 = alph1 + s*(alph2-a2)
        
        self.svm.alphas[i1], self.svm.alphas[i2] = a1, a2
        self._fcache = {}

        i_low = i_up = b_low = b_up = None
        for i in range(len(self.svm.alphas)):
            if self._in_I_0(i) or self._in_I_3(i) or self._in_I_4(i):
                f = self._F(i)
                if b_low is None or f > b_low:
                    b_low = f
                    i_low = i
            if self._in_I_0(i) or self._in_I_1(i) or self._in_I_2(i):
                f = self._F(i)
                if b_up is None or f < b_up:
                    b_up = f
                    i_up = i
        self.i_low, self.i_up = i_low, i_up
        self.b_low, self.b_up = b_low, b_up

        return 1

    def _F(self, i):
        if not self._fcache.has_key(i):
            # assume b == 0
            self._fcache[i] = classify(self.svm, self.svm.xs[i])-self.svm.ys[i]
        return self._fcache[i]

    def _in_I_0(self, i):
        return self.svm.alphas[i] > 0 and self.svm.alphas[i] < self.C

    def _in_I_1(self, i):
        return self.svm.alphas[i] == 0 and self.svm.ys[i] == 1

    def _in_I_2(self, i):
        return self.svm.alphas[i] == self.C and self.svm.ys[i] == -1

    def _in_I_3(self, i):
        return self.svm.alphas[i] == self.C and self.svm.ys[i] == 1

    def _in_I_4(self, i):
        return self.svm.alphas[i] == 0 and self.svm.ys[i] == -1

class TransductiveTrainer:
    """Transductive Support Vector Machine Trainer

    This is used to train transductive support vector machines described
    in:
    Joachims, Thorsten.  Transductive Inference for Text Classification
    using Support Vector Machines.
    XXX url?

    The algorithm is modified from Platt's Sequential Minimal Optimization.

    Methods:
    train     Train a new SVM.

    """
    def train(self, training_set, results, num_test_cases,
              kernel_fn, C, C_negtest, C_postest, epsilon,
              update_fn=None):
        """train(self, training_set, results, num_test_cases,
        kernel_fn, C, C_negtest, C_postest, epsilon,
        update_fn=None) -> SVM

        """
        # Initializing variables
        self.svm = SVM(training_set, results,
                       [0.0]*len(training_set), 0.0, kernel_fn)
        self.eps = epsilon
        self.num_test_cases = num_test_cases
        self.C = C
        self.C_negtest, self.C_postest = C_negtest, C_postest
        self._ecache = {}

        # Training algorithm starts here.
        num_changed = 0
        examine_all = 1

        while num_changed > 0 or examine_all:
            num_changed = 0
            if examine_all:
                for i2 in range(len(self.svm.alphas)):
                    num_changed = num_changed + self._examine_example(i2)
            else:
                for i2 in range(len(self.svm.alphas)):
                    # skip the bound alphas
                    if self._is_bound(self.svm.alphas[i2], self._C(i2)):
                        continue
                    num_changed = num_changed + self._examine_example(i2)
            if update_fn is not None:
                update_fn(num_changed, self.svm)
            if examine_all:
                examine_all = 0
            elif num_changed == 0:
                examine_all = 1

        return self.svm

    def _examine_example(self, i2):
        y2 = self.svm.ys[i2]
        alph2 = self.svm.alphas[i2]
        E2 = self._E(i2)
        r2 = E2*y2

        if (r2 < -self.eps and alph2 < self._C(i2)) or \
           (r2 > self.eps and alph2 > 0):
            # sort it by E1 - E2
            max_i1 = None
            max_error = None
            for i1 in range(len(self.svm.alphas)):
                if i1 == i2 or self._is_bound(self.svm.alphas[i1],
                                              self._C(i1)):
                    continue
                E1 = self._E(i1)
                error = math.fabs(E1 - E2)
                if max_i1 is None or error > max_error:
                    max_error = error
                    max_i1 = i1
            i1 = max_i1
            if i1 is not None:
                if self._take_step(i1, i2):
                    return 1
                
            i1 = random.randint(0, len(self.svm.alphas)-1)
            for i in range(len(self.svm.alphas)):
                if self._is_bound(self.svm.alphas[i1], self._C(i1)):
                    continue
                if self._take_step(i1, i2):
                    return 1
                i1 = (i1 + 1) % len(self.svm.alphas)

            i1 = random.randint(0, len(self.svm.alphas)-1)
            for i in range(len(self.svm.alphas)):
                if self._take_step(i1, i2):
                    return 1
                i1 = (i1 + 1) % len(self.svm.alphas)
        return 0

    def _take_step(self, i1, i2):
        if i1 == i2:
            return 0
        point = self.svm.xs
        alph1, alph2 = self.svm.alphas[i1], self.svm.alphas[i2]
        y1, y2 = self.svm.ys[i1], self.svm.ys[i2]
        E1 = self._E(i1)
        E2 = self._E(i2)
        s = y1*y2
        b = self.svm.b
        eps = self.eps
        C1, C2 = self._C(i1), self._C(i2)

        if y1 == y2:
            L = max(0.0, alph2+alph1-C1)
            H = min(C2, alph2+alph1)
        else:
            L = max(0.0, alph2-alph1)
            H = min(C2, C1+alph2-alph1)
        if L == H:
            return 0

        k11 = self.svm.kernel_fn(point[i1], point[i1])
        k12 = self.svm.kernel_fn(point[i1], point[i2])
        k22 = self.svm.kernel_fn(point[i2], point[i2])
        eta = k11+k22-2.0*k12
        if eta > 0.0:
            a2 = alph2 + y2*(E1-E2)/eta
            if a2 < L: a2 = L
            elif a2 > H: a2 = H
        else:
            f1 = y1*(E1+b) - alph1*k11 - s*alph2*k22
            f2 = y2*(E2+b) - s*alph1*k12 - alph2*k22
            L1 = alph1 + s*(alph2-L)
            H1 = alph1 + s*(alph2-H)
            Lobj = L1*f1 + L*f2 + 0.5*L1*L1*k11 + 0.5*L*L*k22 + s*L*L1*k12
            Hobj = H1*f1 + H*f2 + 0.5*H1*H1*k11 + 0.5*H*H*k22 + s*H*H1*k12
            
            if Lobj < Hobj-eps:
                a2 = L
            elif Lobj > Hobj+eps:
                a2 = H
            else:
                a2 = alph2

        if math.fabs(a2-alph2) < eps*(a2+alph2+eps):
            return 0

        a1 = alph1 + s*(alph2-a2)

        if not self._is_bound(a1, C1):
            b1 = E1 + y1*(a1-alph1)*k11 + y2*(a2-alph2)*k12 + b
            b = b1
        elif not self._is_bound(a2, C2):
            b2 = E2 + y1*(a1-alph1)*k12 + y2*(a2-alph2)*k22 + b
            b = b2
        else:
            b1 = E1 + y1*(a1-alph1)*k11 + y2*(a2-alph2)*k12 + b
            b2 = E2 + y1*(a1-alph1)*k12 + y2*(a2-alph2)*k22 + b
            b = (b1+b2)/2.0

        self.svm.alphas[i1], self.svm.alphas[i2] = a1, a2
        self.svm.b = b
        self._ecache = {}

        return 1

    def _E(self, i):
        if not self._ecache.has_key(i):
            self._ecache[i] = classify(self.svm, self.svm.xs[i])-self.svm.ys[i]
        return self._ecache[i]

    def _C(self, i):
        if i < len(self.svm.ys)-self.num_test_cases:
            return self.C
        elif self.svm.ys[i] > 0:
            return self.C_postest
        return self.C_negtest

    def _is_bound(self, alpha, C):
        return alpha <= 0 or alpha >= C

def train(training_set, results,
          kernel_fn=LinearKernel(), C=1.0, epsilon=1E-3, update_fn=None):
    """train(training_set, results,
    kernel_fn=LinearKernel(), C=1.0, epsilon=1E-3, update_fn=None) -> SVM

    Train a new support vector machine.  training_set is a list
    containing a list of numbers, each which is a vector in the
    training set.  results is a parallel list that contains either
    1 or -1, which describes the class that the training vector belongs
    to.

    kernel_fn, C, and epsilon are optional parameters used to tune the
    support vector machine.

    update_fn is an optional callback function that called at the
    end of every round of training.  It should take 2 parameters:
        update_fn(num_changed, svm)

    """
    svm = KeerthiTrainer().train(
        training_set, results, kernel_fn, C, epsilon, update_fn=update_fn)
    return svm

def train_transductive(training_set, results, test_set, num_pos,
                       kernel_fn=LinearKernel(), C=1.0, C_test=0.5,
                       epsilon=1E-3, update_fn=None):
    """train_transductive(training_set, results, test_set, num_pos,
    kernel_fn=LinearKernel(), C=1.0, C_test=0.5,
    epsilon=1E-3, update_fn=None) -> SVM

    XXX document this

    """
    svm = train(training_set, results,
                kernel_fn=kernel_fn, C=C, epsilon=epsilon, update_fn=update_fn)

    # Assign the top num_pos results into the positive class.
    # Everything else is in the negative.
    test_scores = []   # list of tuples (score, index)
    for i in range(len(test_set)):
        score = classify(svm, test_set[i])
        test_scores.append((score, i))
    test_scores.sort()
    test_results = [None] * len(test_scores)
    for i in range(len(test_scores)):
        score, index = test_scores[i]
        if i < num_pos:
            test_results[index] = 1
        else:
            test_results[index] = -1
    
    trainer = TransductiveTrainer()
    C_neg = 1E-5
    C_pos = 1E-5 * num_pos/(len(test_set)-num_pos)
    while C_neg < C_test or C_pos < C_test:
        svm = trainer.train(training_set+test_set, results+test_results,
                            len(test_set), kernel_fn, C, C_neg, C_pos,
                            epsilon, update_fn=update_fn)
        changed = 1
        while changed:
            changed = 0
            # Look for an m and l that should be in different classes, and
            # that are both misclassified.  If I see such a case, then
            # switch them.
            #
            # The amount misclassified is:
            # gamma = 1 - y*score
            # If the sample is misclassified, y*score is negative.
            # If it's within the margin, y*score is positive and < 1.
            # Thus, I need to check to make sure gamma > 0, and the
            # sum of the gammas > 2.

            total = len(test_set)**2  # total number of combinations to try
            offset = random.randint(0, total-1)  # start at a random point
            for i in range(total):
                index = (i + offset) % total
                m = index / len(test_set)
                l = index % len(test_set)

                # Make sure m and l should be in different classes.
                if test_results[m] == test_results[l]:
                    continue
                # Make sure they're both misclassified.
                gamma_m = 1.0 - test_results[m]*classify(svm, test_set[m])
                gamma_l = 1.0 - test_results[l]*classify(svm, test_set[l])
                if gamma_m <= 0 or gamma_l <= 0 or (gamma_m+gamma_l) < 2:
                    continue
                
                # Switch them, and reclassify.
                test_results[l], test_results[m] = \
                                 test_results[m], test_results[l]
                
                svm = trainer.train(training_set+test_set,results+test_results,
                                    len(test_set), kernel_fn, C, C_neg, C_pos,
                                    epsilon, update_fn=update_fn)
                changed = 1
                break
                    
        C_neg = min(C_neg*2, C_test)
        C_pos = max(C_pos*2, C_test)
        
    return svm
    
def _dot(x, y):
    """_dot(x, y) -> x*y

    Return the dot product of x and y, where x and y are lists of numbers.

    """
    if len(x) != len(y):
        raise ValueError, "vectors must be same length"
    sum = 0.0
    for i in range(len(x)):
        sum = sum + x[i]*y[i]
    return sum

def _sparse_dot(x, y):
    sum = 0.0
    for k in x.keys():
        if y.has_key(k):
            sum = sum + x[k] * y[k]
    return sum

def _subtract(x, y):
    """_subtract(x, y) -> x-y

    Return x subtract y, where x and y are lists of numbers.

    """
    if len(x) != len(y):
        raise ValueError, "vectors must be same length"
    s = [None] * len(x)
    for i in range(len(x)):
        s[i] = x[i]-y[i]
    return s

def _sparse_subtract(x, y):
    done = {}  # have I handled this index yet?
    vec = {}
    for k in x.keys():
        if y.has_key(k):
            vec[k] = x[k] - y[k]
        else:
            vec[k] = x[k]
        done[k] = 1
    for k in y.keys():
        if done.has_key(k):
            continue
        vec[k] = -y[k]
    return vec

def _sign(x):
    """_sign(x) -> 1 or -1

    Return 1/-1 depending on the sign of x.

    """
    if x >= 0:
        return 1
    return -1

import cSVM
_sparse_dot = cSVM._sparse_dot
_dot = cSVM._dot
classify = cSVM.classify
