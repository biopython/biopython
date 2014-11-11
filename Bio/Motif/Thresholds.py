# Copyright 2008 by Norbert Dojer.  All rights reserved.
# Adapted by Bartek Wilczynski.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Approximate calculation of appropriate thresholds for motif finding
"""
import math
import random


class ScoreDistribution(object):
    """ Class representing approximate score distribution for a given motif.

    Utilizes a dynamic programming approch to calculate the distribution of
    scores with a predefined precision. Provides a number of methods for calculating
    thresholds for motif occurences.
    """
    def __init__(self, motif, precision=10**3):
        self.min_score=min(0.0, motif.min_score())
        self.interval=max(0.0, motif.max_score())-self.min_score
        self.n_points=precision*motif.length
        self.step=self.interval/(self.n_points-1)
        self.mo_density=[0.0]*self.n_points
        self.mo_density[-self._index_diff(self.min_score)]=1.0
        self.bg_density=[0.0]*self.n_points
        self.bg_density[-self._index_diff(self.min_score)]=1.0
        self.ic=motif.ic()
        for lo, mo in zip(motif.log_odds(), motif.pwm()):
            self.modify(lo, mo, motif.background)

    def _index_diff(self, x, y=0.0):
        return int((x-y+0.5*self.step)//self.step)

    def _add(self, i, j):
        return max(0, min(self.n_points-1, i+j))

    def modify(self, scores, mo_probs, bg_probs):
        mo_new=[0.0]*self.n_points
        bg_new=[0.0]*self.n_points
        for k, v in scores.items():
            d=self._index_diff(v)
            for i in range(self.n_points):
                mo_new[self._add(i, d)]+=self.mo_density[i]*mo_probs[k]
                bg_new[self._add(i, d)]+=self.bg_density[i]*bg_probs[k]
        self.mo_density=mo_new
        self.bg_density=bg_new

    def threshold_fpr(self, fpr):
        """
        Approximate the log-odds threshold which makes the type I error (false positive rate).
        """
        i=self.n_points
        prob=0.0
        while prob<fpr:
            i-=1
            prob+=self.bg_density[i]
        return self.min_score+i*self.step

    def threshold_fnr(self, fnr):
        """
        Approximate the log-odds threshold which makes the type II error (false negative rate).
        """
        i=-1
        prob=0.0
        while prob<fnr:
            i+=1
            prob+=self.mo_density[i]
        return self.min_score+i*self.step

    def threshold_balanced(self, rate_proportion=1.0, return_rate=False):
        """
        Approximate the log-odds threshold which makes FNR equal to FPR times rate_proportion
        """
        i=self.n_points
        fpr=0.0
        fnr=1.0
        while fpr*rate_proportion<fnr:
            i-=1
            fpr+=self.bg_density[i]
            fnr-=self.mo_density[i]
        if return_rate:
            return self.min_score+i*self.step, fpr
        else:
            return self.min_score+i*self.step

    def threshold_patser(self):
        """Threshold selection mimicking the behaviour of patser (Hertz, Stormo 1999) software.

        It selects such a threshold that the log(fpr)=-ic(M)
        note: the actual patser software uses natural logarithms instead of log_2, so the numbers
        are not directly comparable.
        """
        return self.threshold_fpr(fpr=2**-self.ic)
