# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import generators
from Bio.SubsMat import FreqTable

class Motif(object):
    """
    A class representing the sequence motifs.
    """
    def __init__(self):
        self.instances = []
        self.score = 0.0
        self.mask = []
        self._pwm_is_current = 0
        self._pwm = []
        self.alphabet=None
        self.length=None

    def _check_length(self, len):
        if self.length==None:
            self.length = len
        elif self.length != len:
            raise ValueError, "You can't change the length of the motif"

    def _check_alphabet(self,alphabet):
        if self.alphabet==None:
            self.alphabet=alphabet
        elif self.alphabet != alphabet:
                raise ValueError, "Wrong Alphabet"
        
    def add_instance(self,instance):
        """
        adds new instance to the motif
        """
        self._check_alphabet(instance.alphabet)
        self._check_length(len(instance))
        self.instances.append(instance)
        self._pwm_is_current = 0

    def set_mask(self,mask):
        """
        sets the mask for the motif

        The mask should be a string containing asterisks in the position of significant columns and spaces in other columns
        """
        self._check_length(len(mask))
        self.mask=[]
        for char in mask:
            if char=="*":
                self.mask.append(1)
            elif char==" ":
                self.mask.append(0)
            else:
                raise ValueError("Mask should contain only '*' or ' ' and not a '%s'"%char)

    def pwm(self):
        """
        returns the PWM computed for the set of instances
        """
        
        if self._pwm_is_current:
            return self._pwm
        #we need to compute new pwm
        self._pwm = []
        for i in xrange(len(self.mask)):
            dict = {}
            #filling the dict with 0's
            for letter in self.alphabet.letters:
                dict[letter]=0
            #counting the occurences of letters in instances
            for seq in self.instances:
                dict[seq[i]]=dict[seq[i]]+1
            self._pwm.append(FreqTable.FreqTable(dict,FreqTable.COUNT,self.alphabet)) 
        self._pwm_is_current=1
        return self._pwm

    def search_instances(self,sequence):
        """
        a generator function, returning found positions of instances of the motif in a given sequence
        """

        for pos in xrange(0,len(sequence)-self.length+1):
            for instance in self.instances:
                print "testing",instance.tostring(),sequence[pos:pos+self.length].tostring()
                if instance.tostring()==sequence[pos:pos+self.length].tostring():
                    yield(pos,instance)
                    break # no other instance will fit (we don't want to return multiple hits)

    def score_hit(self,sequence,position,normalized=1):
        """
        give the pwm score for a given position
        """
        score = 0.0
        for pos in xrange(self.length):
            score += self.pwm()[pos][sequence[position+pos]]
        if normalized:
            score/=self.length
        return score
    
    def search_pwm(self,sequence,threshold=0.0):
        """
        a generator function, returning found hits in a given sequence with the pwm score higher than the threshold
        """

        for pos in xrange(0,len(sequence)-self.length+1):
            score = self.score_hit(sequence,pos)
            if score > threshold:
                yield (pos,score)
    
            
            

        
        
