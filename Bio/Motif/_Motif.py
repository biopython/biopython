# Copyright 2003-2009 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Implementation of sequence motifs (PRIVATE).
"""
from Bio.Seq import Seq
from Bio.SubsMat import FreqTable
from Bio.Alphabet import IUPAC
import math,random

class Motif(object):
    """
    A class representing sequence motifs.
    """
    def __init__(self,alphabet=IUPAC.unambiguous_dna):
        self.instances = []
        self.has_instances=False
        self.counts = {}
        self.has_counts=False
        self.mask = []
        self._pwm_is_current = False
        self._pwm = []
        self._log_odds_is_current = False
        self._log_odds = []
        self.alphabet=alphabet
        self.length=None
        self.background=dict((n, 1.0/len(self.alphabet.letters)) \
                             for n in self.alphabet.letters)
        self.beta=1.0
        self.info=None
        self.name=""

    def _check_length(self, len):
        if self.length==None:
            self.length = len
        elif self.length != len:
            print "len",self.length,self.instances, len
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
        if self.has_counts:
            for i in range(self.length):
                let=instance[i]
                self.counts[let][i]+=1

        if self.has_instances or not self.has_counts:
            self.instances.append(instance)
            self.has_instances=True
            
        self._pwm_is_current = False
        self._log_odds_is_current = False

 
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
    
    def pwm(self,laplace=True):
        """
        returns the PWM computed for the set of instances

        if laplace=True (default), pseudocounts equal to self.background multiplied by self.beta are added to all positions.
        """
        
        if self._pwm_is_current:
            return self._pwm
        #we need to compute new pwm
        self._pwm = []
        for i in xrange(self.length):
            dict = {}
            #filling the dict with 0's
            for letter in self.alphabet.letters:
                if laplace:
                    dict[letter]=self.beta*self.background[letter]
                else:
                    dict[letter]=0.0
            if self.has_counts:
                #taking the raw counts
                for letter in self.alphabet.letters:
                    dict[letter]+=self.counts[letter][i]
            elif self.has_instances:
                #counting the occurences of letters in instances
                for seq in self.instances:
                    #dict[seq[i]]=dict[seq[i]]+1
                    try:
                        dict[seq[i]]+=1
                    except KeyError: #we need to ignore non-alphabet letters
                        pass
            self._pwm.append(FreqTable.FreqTable(dict,FreqTable.COUNT,self.alphabet)) 
        self._pwm_is_current=1
        return self._pwm

    def log_odds(self,laplace=True):
        """
        returns the logg odds matrix computed for the set of instances
        """
        
        if self._log_odds_is_current:
            return self._log_odds
        #we need to compute new pwm
        self._log_odds = []
        pwm=self.pwm(laplace)
        for i in xrange(self.length):
            d = {}
            for a in self.alphabet.letters:
                    d[a]=math.log(pwm[i][a]/self.background[a],2)
            self._log_odds.append(d)
        self._log_odds_is_current=1
        return self._log_odds

    def ic(self):
        """Method returning the information content of a motif.
        """
        res=0
        pwm=self.pwm()
        for i in range(self.length):
            res+=2
            for a in self.alphabet.letters:
                if pwm[i][a]!=0:
                    res+=pwm[i][a]*math.log(pwm[i][a],2)
        return res

    def exp_score(self,st_dev=False):
        """
        Computes expected score of motif's instance and its standard deviation
        """
        exs=0.0
        var=0.0
        pwm=self.pwm()
        for i in range(self.length):
            ex1=0.0
            ex2=0.0
            for a in self.alphabet.letters:
                if pwm[i][a]!=0:
                    ex1+=pwm[i][a]*(math.log(pwm[i][a],2)-math.log(self.background[a],2))
                    ex2+=pwm[i][a]*(math.log(pwm[i][a],2)-math.log(self.background[a],2))**2
            exs+=ex1
            var+=ex2-ex1**2
        if st_dev:
            return exs,math.sqrt(var)
        else:
            return exs

    def search_instances(self,sequence):
        """
        a generator function, returning found positions of instances of the motif in a given sequence
        """
        if not self.has_instances:
            raise ValueError ("This motif has no instances")
        for pos in xrange(0,len(sequence)-self.length+1):
            for instance in self.instances:
                if instance.tostring()==sequence[pos:pos+self.length].tostring():
                    yield(pos,instance)
                    break # no other instance will fit (we don't want to return multiple hits)

    def score_hit(self,sequence,position,normalized=0,masked=0):
        """
        give the pwm score for a given position
        """
        lo=self.log_odds()
        score = 0.0
        for pos in xrange(self.length):
            a = sequence[position+pos]
            if not masked or self.mask[pos]:
                try:
                    score += lo[pos][a]
                except:
                    pass
        if normalized:
            if not masked:
                score/=self.length
            else:
                score/=len([x for x in self.mask if x])
        return score
    
    def search_pwm(self,sequence,normalized=0,masked=0,threshold=0.0,both=True):
        """
        a generator function, returning found hits in a given sequence with the pwm score higher than the threshold
        """
        if both:
            rc = self.reverse_complement()
            
        sequence=sequence.tostring().upper()
        for pos in xrange(0,len(sequence)-self.length+1):
            score = self.score_hit(sequence,pos,normalized,masked)
            if score > threshold:
                yield (pos,score)
            if both:
                rev_score = rc.score_hit(sequence,pos,normalized,masked)
                if rev_score > threshold:
                    yield (-pos,rev_score)

    def dist_pearson(self, motif, masked = 0):
        """
        return the similarity score based on pearson correlation for the given motif against self.

        We use the Pearson's correlation of the respective probabilities.
        """

        if self.alphabet != motif.alphabet:
            raise ValueError("Cannot compare motifs with different alphabets")

        max_p=-2
        for offset in range(-self.length+1,motif.length):
            if offset<0:
                p = self.dist_pearson_at(motif,-offset)
            else: #offset>=0
                p = motif.dist_pearson_at(self,offset)
            
            if max_p<p:
                max_p=p
                max_o=-offset
        return 1-max_p,max_o

    def dist_pearson_at(self,motif,offset):
        sxx = 0 # \sum x^2
        sxy = 0 # \sum x \cdot y
        sx = 0  # \sum x
        sy = 0  # \sum y
        syy = 0 # \sum x^2
        norm=max(self.length,offset+motif.length)
        
        for pos in range(max(self.length,offset+motif.length)):
            for l in self.alphabet.letters:
                xi = self[pos][l]
                yi = motif[pos-offset][l]
                sx = sx + xi
                sy = sy + yi
                sxx = sxx + xi * xi
                syy = syy + yi * yi
                sxy = sxy + xi * yi

        norm *= len(self.alphabet.letters)
        s1 = (sxy - sx*sy*1.0/norm)
        s2 = (norm*sxx - sx*sx*1.0)*(norm*syy- sy*sy*1.0)
        return s1/math.sqrt(s2)

    def dist_product(self,other):
        """
        A similarity measure taking into account a product probability of generating overlaping instances of two motifs
        """
        max_p=0.0
        for offset in range(-self.length+1,other.length):
            if offset<0:
                p = self.dist_product_at(other,-offset)
            else: #offset>=0
                p = other.dist_product_at(self,offset)
            if max_p<p:
                max_p=p
                max_o=-offset
        return 1-max_p/self.dist_product_at(self,0),max_o
            
    def dist_product_at(self,other,offset):
        s=0
        for i in range(max(self.length,offset+other.length)):
            f1=self[i]
            f2=other[i-offset]
            for n,b in self.background.iteritems():
                s+=b*f1[n]*f2[n]
        return s/i

    def dist_dpq(self,other):
        r"""Calculates the DPQ distance measure between motifs.

        It is calculated as a maximal value of DPQ formula (shown using LaTeX
        markup, familiar to mathematicians):
        
        \sqrt{\sum_{i=1}^{alignment.len()} \sum_{k=1}^alphabet.len() \
        \{ m1[i].freq(alphabet[k])*log_2(m1[i].freq(alphabet[k])/m2[i].freq(alphabet[k])) +
           m2[i].freq(alphabet[k])*log_2(m2[i].freq(alphabet[k])/m1[i].freq(alphabet[k]))
        }
        
        over possible non-spaced alignemts of two motifs.  See this reference:

        D. M Endres and J. E Schindelin, "A new metric for probability
        distributions", IEEE transactions on Information Theory 49, no. 7
        (July 2003): 1858-1860.
        """
        
        min_d=float("inf")
        min_o=-1
        d_s=[]
        for offset in range(-self.length+1,other.length):
            #print "%2.3d"%offset,
            if offset<0:
                d = self.dist_dpq_at(other,-offset)
                overlap = self.length+offset
            else: #offset>=0
                d = other.dist_dpq_at(self,offset)
                overlap = other.length-offset
            overlap = min(self.length,other.length,overlap)
            out = self.length+other.length-2*overlap
            #print d,1.0*(overlap+out)/overlap,d*(overlap+out)/overlap
            #d = d/(2*overlap)
            d = (d/(out+overlap))*(2*overlap+out)/(2*overlap)
            #print d
            d_s.append((offset,d))
            if min_d> d:
                min_d=d
                min_o=-offset
        return min_d,min_o#,d_s
            
    def dist_dpq_at(self,other,offset):
        """
        calculates the dist_dpq measure with a given offset.

        offset should satisfy 0<=offset<=len(self)
        """
        def dpq (f1,f2,alpha):
            s=0
            for n in alpha.letters:
                avg=(f1[n]+f2[n])/2
                s+=f1[n]*math.log(f1[n]/avg,2)+f2[n]*math.log(f2[n]/avg,2)
            return math.sqrt(s)
                
        s=0        
        for i in range(max(self.length,offset+other.length)):
            f1=self[i]
            f2=other[i-offset]
            s+=dpq(f1,f2,self.alphabet)
        return s
            
    def _read(self,stream):
        """Reads the motif from the stream (in AlignAce format).

        the self.alphabet variable must be set beforehand.
        If the last line contains asterisks it is used for setting mask
        """
        
        while 1:
            ln = stream.readline()
            if "*" in ln:
                self.set_mask(ln.strip("\n\c"))
                break
            self.add_instance(Seq(ln.strip(),self.alphabet))
        
    def __str__(self,masked=False):
        """ string representation of a motif.
        """
        str = ""
        for inst in self.instances:
            str = str + inst.tostring() + "\n"

        if masked:
            for i in xrange(self.length):
                if self.mask[i]:
                    str = str + "*"
                else:
                    str = str + " "
            str = str + "\n"
        return str

    def __len__(self):
        """return the length of a motif

        Please use this method (i.e. invoke len(m)) instead of refering to the m.length directly.
        """
        if self.length==None:
            return 0
        else:
            return self.length
        
    def _write(self,stream):
        """
        writes the motif to the stream
        """

        stream.write(self.__str__())
            
            

    def _to_fasta(self):
        """
        FASTA representation of motif
        """
        if not self.has_instances:
            self.make_instances_from_counts()
        str = ""
        for i,inst in enumerate(self.instances):
            str = str + ">instance%d\n"%i + inst.tostring() + "\n"
            
        return str       

    def reverse_complement(self):
        """
        Gives the reverse complement of the motif
        """
        res = Motif()
        if self.has_instances:
            for i in self.instances:
                res.add_instance(i.reverse_complement())
        else: # has counts
            res.has_counts=True
            res.counts["A"]=self.counts["T"][:]
            res.counts["T"]=self.counts["A"][:]
            res.counts["G"]=self.counts["C"][:]
            res.counts["C"]=self.counts["G"][:]
            res.counts["A"].reverse()
            res.counts["C"].reverse()
            res.counts["G"].reverse()
            res.counts["T"].reverse()
            res.length=self.length
        res.mask = self.mask
        return res

        
    def _from_jaspar_pfm(self,stream,make_instances=False):
        """
        reads the motif from Jaspar .pfm file

        The instances are fake, but the pwm is accurate.
        """
        return self._from_horiz_matrix(stream,letters="ACGT",make_instances=make_instances)

    def _from_vert_matrix(self,stream,letters=None,make_instances=False):
        """reads a vertical count matrix from stream and fill in the counts.
        """

        self.counts = {}
        self.has_counts=True
        if letters==None:
            letters=self.alphabet.letters
        self.length=0
        for i in letters:
            self.counts[i]=[]
        for ln in stream.readlines():
            rec=map(float,ln.strip().split())
            for k,v in zip(letters,rec):
                self.counts[k].append(v)
            self.length+=1
        self.set_mask("*"*self.length)
        if make_instances==True:
            self.make_instances_from_counts()
        return self
        
    def _from_horiz_matrix(self,stream,letters=None,make_instances=False):
        """reads a horizontal count matrix from stream and fill in the counts.
        """
        if letters==None:
            letters=self.alphabet.letters
        self.counts = {}
        self.has_counts=True

        for i in letters:
            ln = stream.readline().strip().split()
            #if there is a letter in the beginning, ignore it
            if ln[0]==i:
                ln=ln[1:]
            #print ln
            try:
                self.counts[i]=map(int,ln)
            except ValueError: #not integers
                self.counts[i]=map(float,ln) #map(lambda s: int(100*float(s)),ln)
            #print counts[i]
        
        s = sum(self.counts[nuc][0] for nuc in letters)
        l = len(self.counts[letters[0]])
        self.length=l
        self.set_mask("*"*l)
        if make_instances==True:
            self.make_instances_from_counts()
        return self
        

    def make_instances_from_counts(self):
        """Creates "fake" instances for a motif created from a count matrix.

        In case the sums of counts are different for different columnes, the
        shorter columns are padded with background.
        """
        alpha="".join(self.alphabet.letters)
        #col[i] is a column taken from aligned motif instances
        col=[]
        self.has_instances=True
        self.instances=[]
        s = sum(map(lambda nuc: self.counts[nuc][0],self.alphabet.letters))
        for i in range(self.length):
            col.append("")
            for n in self.alphabet.letters:
                col[i] = col[i]+ (n*(self.counts[n][i]))
            if len(col[i])<s:
                print "WARNING, column too short",len(col[i]),s
                col[i]+=(alpha*s)[:(s-len(col[i]))]
            #print i,col[i]
        #iterate over instances
        for i in range(s): 
            inst="" #start with empty seq
            for j in range(self.length): #iterate over positions
                inst+=col[j][i]
            #print i,inst
            inst=Seq(inst,self.alphabet)                
            self.add_instance(inst)
        return self.instances

    def make_counts_from_instances(self):
        """Creates the count matrix for a motif with instances.

        """
        #make strings for "columns" of motifs
        #col[i] is a column taken from aligned motif instances
        counts={}
        for a in self.alphabet.letters:
            counts[a]=[]
        self.has_counts=True
        s = len(self.instances)
        for i in range(self.length):
            ci = dict((a,0) for a in self.alphabet.letters)
            for inst in self.instances:
                ci[inst[i]]+=1
            for a in self.alphabet.letters:
                counts[a].append(ci[a])
        self.counts=counts
        return counts

    def _from_jaspar_sites(self,stream):
        """
        reads the motif from Jaspar .sites file

        The instances and pwm are OK.
        """
        
        while True:
            ln = stream.readline()# read the header "$>...."
            if ln=="" or ln[0]!=">":
                break
            
            ln=stream.readline().strip()#read the actual sequence
            i=0
            while ln[i]==ln[i].lower():
                i+=1
            inst=""
            while i<len(ln) and ln[i]==ln[i].upper():
                inst+=ln[i]
                i+=1
            inst=Seq(inst,self.alphabet)                
            self.add_instance(inst)

        self.set_mask("*"*len(inst))
        return self


    def __getitem__(self,index):
        """Returns the probability distribution over symbols at a given position, padding with background.

        If the requested index is out of bounds, the returned distribution comes from background.
        """
        if index in range(self.length):
            return self.pwm()[index]
        else:
            return self.background

    def consensus(self):
        """Returns the consensus sequence of a motif.
        """
        res=""
        for i in range(self.length):
            max_f=0
            max_n="X"
            for n in sorted(self[i]):
                if self[i][n]>max_f:
                    max_f=self[i][n]
                    max_n=n
            res+=max_n
        return Seq(res,self.alphabet)

    def anticonsensus(self):
        """returns the least probable pattern to be generated from this motif.
        """
        res=""
        for i in range(self.length):
            min_f=10.0
            min_n="X"
            for n in sorted(self[i]):
                if self[i][n]<min_f:
                    min_f=self[i][n]
                    min_n=n
            res+=min_n
        return Seq(res,self.alphabet)

    def max_score(self):
        """Maximal possible score for this motif.

        returns the score computed for the consensus sequence.
        """
        return self.score_hit(self.consensus(),0)
    
    def min_score(self):
        """Minimal possible score for this motif.

        returns the score computed for the anticonsensus sequence.
        """
        return self.score_hit(self.anticonsensus(),0)

    def weblogo(self,fname,format="PNG",**kwds):
        """
        uses the Berkeley weblogo service to download and save a weblogo of itself
        
        requires an internet connection.
        The parameters from **kwds are passed directly to the weblogo server.
        """
        import urllib
        import urllib2
        al= self._to_fasta()
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
        for k,v in kwds.iteritems():
            values[k]=str(v)
            
        data = urllib.urlencode(values)
        req = urllib2.Request(url, data)
        response = urllib2.urlopen(req)
        f=open(fname,"w")
        im=response.read()
        
        f.write(im)
        f.close()
  

    def _to_transfac(self):
        """Write the representation of a motif in TRANSFAC format
        """
        res="XX\nTY Motif\n" #header
        try:
            res+="ID %s\n"%self.name
        except:
            pass
        res+="BF undef\nP0"
        for a in self.alphabet.letters:
            res+=" %s"%a
        res+="\n"
        if not self.has_counts:
            self.make_counts_from_instances()
        for i in range(self.length):
            if i<9:
                res+="0%d"%(i+1)
            else:
                res+="%d"%(i+1)
            for a in self.alphabet.letters:
                res+=" %d"%self.counts[a][i]
            res+="\n"
        res+="XX\n"
        return res

    def _to_vertical_matrix(self,letters=None):
        """Return string representation of the motif as  a matrix.
        
        """
        if letters==None:
            letters=self.alphabet.letters
        self._pwm_is_current=False
        pwm=self.pwm(laplace=False)
        res=""
        for i in range(self.length):
            res+="\t".join([str(pwm[i][a]) for a in letters])
            res+="\n"
        return res
    
    def _to_horizontal_matrix(self,letters=None,normalized=True):
        """Return string representation of the motif as  a matrix.
        
        """
        if letters==None:
            letters=self.alphabet.letters
        res=""
        if normalized: #output PWM
            self._pwm_is_current=False
            mat=self.pwm(laplace=False)
            for a in letters:
                res+="\t".join([str(mat[i][a]) for i in range(self.length)])
                res+="\n"
        else: #output counts
            if not self.has_counts:
                self.make_counts_from_instances()
            mat=self.counts
            for a in letters:
                res+="\t".join([str(mat[a][i]) for i in range(self.length)])
                res+="\n"
        return res

    def _to_jaspar_pfm(self):
        """Returns the pfm representation of the motif
        """
        return self._to_horizontal_matrix(normalized=False,letters="ACGT")

    def format(self,format):
        """Returns a string representation of the Motif in a given format

        Currently supported fromats:
         - jaspar-pfm : JASPAR Position Frequency Matrix
         - transfac : TRANSFAC like files
         - fasta : FASTA file with instances
        """

        formatters={
            "jaspar-pfm":   self._to_jaspar_pfm,
            "transfac":     self._to_transfac,
            "fasta" :       self._to_fasta,
            }

        try:
            return formatters[format]()
        except KeyError:
            raise ValueError("Wrong format type")

    def scanPWM(self,seq):
        """Matrix of log-odds scores for a nucleotide sequence.
 
        scans a nucleotide sequence and returns the matrix of log-odds
        scores for all positions.

        - the result is a one-dimensional list or numpy array
        - the sequence can only be a DNA sequence
        - the search is performed only on one strand
        """
        if self.alphabet!=IUPAC.unambiguous_dna:
            raise ValueError("Wrong alphabet! Use only with DNA motifs")
        if seq.alphabet!=IUPAC.unambiguous_dna:
            raise ValueError("Wrong alphabet! Use only with DNA sequences")

        seq = seq.tostring()

        # check if the fast C code can be used
        try:
            import _pwm
        except ImportError:
            # use the slower Python code otherwise
            return self._pwm_calculate(seq)
        
        # get the log-odds matrix into a proper shape
        # (each row contains sorted (ACGT) log-odds values)
        logodds=[[y[1] for y in sorted(x.items())] for x in self.log_odds()]
        return _pwm.calculate(seq, logodds)

    def _pwm_calculate(self, sequence):
        logodds = self.log_odds()
        m = len(logodds)
        s = len(sequence)
        n = s - m + 1
        result = [None] * n
        for i in xrange(n):
            score = 0.0
            for j in xrange(m):
                c = sequence[i+j]
                temp = logodds[j].get(c)
                if temp==None:
                    break
                score += temp
            else:
                result[i] = score
        return result
