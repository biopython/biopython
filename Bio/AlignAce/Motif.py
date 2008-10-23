# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Implementation of sequence motifs.

Changes:
9.2007 (BW) : added the to_faste() and .weblogo() methods allowing to use the Berkeley weblogo server at http://weblogo.berkeley.edu/
"""

from Bio.SubsMat import FreqTable

class Motif(object):
    """
    A class representing sequence motifs.
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
            raise ValueError("You can't change the length of the motif")

    def _check_alphabet(self,alphabet):
        if self.alphabet==None:
            self.alphabet=alphabet
        elif self.alphabet != alphabet:
                raise ValueError("Wrong Alphabet")
        
    def add_instance(self,instance):
        """
        adds new instance to the motif
        """
        self._check_alphabet(instance.alphabet)
        self._check_length(len(instance))
        self.instances.append(instance)
        self._pwm_is_current = False

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
        self._pwm_is_current=True
        return self._pwm

    def search_instances(self,sequence):
        """
        a generator function, returning found positions of instances of the motif in a given sequence
        """
        for pos in xrange(0,len(sequence)-self.length+1):
            for instance in self.instances:
                if instance.tostring()==sequence[pos:pos+self.length].tostring():
                    yield(pos,instance)
                    break # no other instance will fit (we don't want to return multiple hits)

    def score_hit(self,sequence,position,normalized=1,masked=0):
        """
        give the pwm score for a given position
        """
        score = 0.0
        for pos in xrange(self.length):
            if not masked or self.mask[pos]:
                score += self.pwm()[pos][sequence[position+pos]]
        if normalized:
            if not masked:
                score/=self.length
            else:
                score/=len(filter(lambda x: x, self.mask))
        return score
    
    def search_pwm(self,sequence,threshold=0.0,normalized=1,masked=1):
        """
        a generator function, returning found hits in a given sequence with the pwm score higher than the threshold
        """

        for pos in xrange(0,len(sequence)-self.length+1):
            score = self.score_hit(sequence,pos,normalized,masked)
            if score > threshold:
                yield (pos,score)

    
    def sim(self, motif, masked = 0):
        """
        return the similarity score for the given motif against self.

        We use the Pearson's correlation of the respective probabilities.
        If the motifs have different length or mask raise the ValueError.
        """

        from math import sqrt
        
        if self.alphabet != motif.alphabet:
            raise ValueError("Wrong alphabet")
        
        if self.length != motif.length:
            raise ValueError("Wrong length")
        
        if masked and self.mask!=motif.mask:
            raise ValueError("Wrong mask")

        sxx = 0 # \sum x^2
        sxy = 0 # \sum x \cdot y
        sx = 0  # \sum x
        sy = 0  # \sum y
        syy = 0 # \sum x^2

        for pos in xrange(self.length):
            if not masked or self.mask:
                for l in self.alphabet.letters:
                    xi = self.pwm()[pos][l]
                    yi = motif.pwm()[pos][l]
                    sx = sx + xi
                    sy = sy + yi
                    sxx = sxx + xi * xi
                    syy = syy + yi * yi
                    sxy = sxy + xi * yi
                    
        if masked:
            norm = len(filter(lambda x: x,self.mask))
        else:
            norm = self.length
            
        norm *= len(self.alphabet.letters)
        s1 = (sxy - sx*sy*1.0/norm)
        s2 = (sxx - sx*sx*1.0/norm)*(syy- sy*sy*1.0/norm)
        
        return s1/sqrt(s2)
            
    def read(self,stream):
        """
        reads the motif from the stream

        the self.alphabet variable must be set before
        """
        from Bio.Seq import Seq
        while 1:
            ln = stream.readline()
            if "*" in ln:
                self.set_mask(ln.strip("\n\c"))
                break
            self.add_instance(Seq(ln.strip(),self.alphabet))
        
    def __str__(self):
        """
        string representation of motif
        """
        str = ""
        for inst in self.instances:
            str = str + inst.tostring() + "\n"

        for i in xrange(self.length):
            if self.mask[i]:
                str = str + "*"
            else:
                str = str + " "
        str = str + "\n"

        return str

    def write(self,stream):
        """
        writes the motif to the stream
        """

        stream.write(self.__str__())
            
            

    def to_fasta(self):
        """
        FASTA representation of motif
        """
        str = ""
        for i,inst in enumerate(self.instances):
            str = str + "> instance %d\n"%i + inst.tostring() + "\n"
            
        return str       
        
    def weblogo(self,fname,format="PNG",**kwds):
        """
        uses the Berkeley weblogo service to download and save a weblogo of itself
        
        requires an internet connection.
        The parameters from **kwds are passed directly to the weblogo server.
        """
        import urllib
        import urllib2
        #import Image
        al= self.to_fasta()

        url = 'http://weblogo.berkeley.edu/logo.cgi'
        values = {'sequence' : al,
                  'format' : format,
                  'logowidth' : '18',
                  'logoheight' : '5',
                  'logounits' : 'cm',
                  'kind' : 'AUTO',
                  'firstnum' : "1",
                  'command' : 'Create Logo',
                  'smallsamplecorrection' : "on",
                  'symbolsperline' : 32,
                  'res' : '96',
                  'res_units' : 'ppi',
                  'antialias' : 'on',
                  'title' : '',
                  'barbits' : '',
                  'xaxis': 'on',
                  'xaxis_label'  : '',
                  'yaxis': 'on',
                  'yaxis_label' : '',
                  'showends' : 'on',
                  'shrink' : '0.5',
                  'fineprint' : 'on',
                  'ticbits' : '1',
                  'colorscheme' : 'DEFAULT',
                  'color1' : 'green',
                  'color2' : 'blue',
                  'color3' : 'red',
                  'color4' : 'black',
                  'color5' : 'purple',
                  'color6' : 'orange',
                  'color1' : 'black',
                  }
        for k,v in kwds.items():
            values[k]=str(v)
            
        data = urllib.urlencode(values)
        req = urllib2.Request(url, data)
        response = urllib2.urlopen(req)
        f=open(fname,"w")
        im=response.read()
        
        f.write(im)
        f.close()
  
