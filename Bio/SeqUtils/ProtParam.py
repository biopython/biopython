# Copyright Yair Benita Y.Benita@pharm.uu.nl
# Biopython (http://biopython.org) license applies

"""Simple protein analysis.

Example,

X = ProteinAnalysis("MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV")
print X.count_amino_acids()
print X.get_amino_acids_percent()
print X.molecular_weight()
print X.aromaticity()
print X.instability_index()
print X.flexibility()
print X.isoelectric_point()
print X.secondary_structure_fraction()
print X.protein_scale(ProtParamData.kd, 9, 0.4)
"""

import sys
import ProtParamData, IsoelectricPoint
from ProtParamData import kd  # Added by Iddo to enable the gravy method
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData
#from BioModule import 

class ProteinAnalysis(object):
    """Class containing methods for protein analysis.

    The class init method takes only one argument, the protein sequence as a
    string and builds a sequence object using the Bio.Seq module. This is done
    just to make sure the sequence is a protein sequence and not anything else.
    
    """
    def __init__(self, ProtSequence):
        if ProtSequence.islower():
            self.sequence = Seq(ProtSequence.upper(), IUPAC.protein)
        else:
            self.sequence = Seq(ProtSequence, IUPAC.protein)
        self.amino_acids_content = None
        self.amino_acids_percent = None
        self.length = len(self.sequence)
        
    def count_amino_acids(self):
        """Count standard amino acids, returns a dict.
            
        Simply counts the number times an amino acid is repeated in the protein
        sequence. Returns a dictionary {AminoAcid:Number} and also stores the
        dictionary in self.amino_acids_content.
        """
        ProtDic = dict([ (k, 0) for k in IUPACData.protein_letters])
        for i in ProtDic:
            ProtDic[i]=self.sequence.count(i)
        self.amino_acids_content = ProtDic
        return ProtDic
    
    def get_amino_acids_percent(self):
        """Calculate the amino acid content in percents.

        The same as count_amino_acids only returns the Number in percentage of
        entire sequence. Returns a dictionary and stores the dictionary in
        self.amino_acids_content_percent.
        
        input is the dictionary from CountAA.
        output is a dictionary with AA as keys.
        """
        if not self.amino_acids_content:
            self.count_amino_acids()
                
        PercentAA = {}
        for i in self.amino_acids_content:
            if self.amino_acids_content[i] > 0:
                PercentAA[i]=self.amino_acids_content[i]/float(self.length)
            else:
                PercentAA[i] = 0
        self.amino_acids_percent = PercentAA
        return PercentAA

    def molecular_weight (self):
        """Calculate MW from Protein sequence"""
        # make local dictionary for speed
        MwDict = {}
        # remove a molecule of water from the amino acid weight.
        for i in IUPACData.protein_weights:
            MwDict[i] = IUPACData.protein_weights[i] - 18.02
        MW = 18.02 # add just one water molecule for the whole sequence.
        for i in self.sequence:
            MW += MwDict[i]
        return MW

    def aromaticity(self):
        """Calculate the aromaticity according to Lobry, 1994.

        Calculates the aromaticity value of a protein according to Lobry, 1994.
        It is simply the relative frequency of Phe+Trp+Tyr.
        """
        if not self.amino_acids_percent:
            self.get_amino_acids_percent()
        
        Arom= self.amino_acids_percent['Y']+self.amino_acids_percent['W']+self.amino_acids_percent['F']
        return Arom

    def instability_index(self):
        """Calculate the instability index according to Guruprasad et al 1990.

        Implementation of the method of Guruprasad et al. 1990 to test a
        protein for stability. Any value above 40 means the protein is unstable
        (has a short half life). 
        
        See: Guruprasad K., Reddy B.V.B., Pandit M.W.
        Protein Engineering 4:155-161(1990).
        """
        #make the dictionary local for speed.
        DIWV=ProtParamData.DIWV.copy()
        score=0.0
        for i in range(self.length - 1):
            DiPeptide=DIWV[self.sequence[i]][self.sequence[i+1]]
            score += DiPeptide
        return (10.0/self.length) * score

    def flexibility(self):
        """Calculate the flexibility according to Vihinen, 1994.
        
        No argument to change window size because parameters are specific for a
        window=9. The parameters used are optimized for determining the flexibility.
        """
        Flex = ProtParamData.Flex.copy()
        Window=9
        Weights=[0.25,0.4375,0.625,0.8125,1]
        List=[]
        for i in range(self.length - Window):
            SubSeq=self.sequence[i:i+Window]
            score = 0.0
            for j in range(Window//2):
                score += (Flex[SubSeq[j]]+Flex[SubSeq[Window-j-1]]) * Weights[j]
            score += Flex[SubSeq[Window//2+1]]
            List.append(score/5.25)
        return List

    def gravy(self):
        """Calculate the gravy according to Kyte and Doolittle."""
        ProtGravy=0.0
        for i in self.sequence:
            ProtGravy += kd[i]
            
        return ProtGravy/self.length

    # this method is used to make a list of relative weight of the
    # window edges compared to the window center. The weights are linear.
    # it actually generates half a list. For a window of size 9 and edge 0.4
    # you get a list of [0.4, 0.55, 0.7, 0.85]. 
    def _weight_list(self, window, edge):
        unit = ((1.0-edge)/(window-1))*2
        list = [0.0]*(window//2)
        for i in range(window//2):
            list[i] = edge + unit * i
        return list
    
    # The weight list returns only one tail. If the list should be [0.4,0.7,1.0,0.7,0.4]
    # what you actually get from _weights_list is [0.4,0.7]. The correct calculation is done
    # in the loop.
    def protein_scale(self, ParamDict, Window, Edge=1.0):
        """Compute a profile by any amino acid scale.
        
        An amino acid scale is defined by a numerical value assigned to each type of
        amino acid. The most frequently used scales are the hydrophobicity or
        hydrophilicity scales and the secondary structure conformational parameters
        scales, but many other scales exist which are based on different chemical and
        physical properties of the amino acids.  You can set several  parameters that
        control the computation  of a scale profile, such as the window size and the
        window edge relative weight value.  WindowSize: The window size is the length
        of the interval to use for the profile computation. For a window size n, we
        use the i- ( n-1)/2 neighboring residues on each side of residue it compute
        the score for residue i. The score for residue is  the sum of the scale values
        for these amino acids,  optionally weighted according to their position in the
        window.  Edge: The central amino acid of the window always has a weight of 1.
        By default, the amino acids at the remaining window positions have the same
        weight, but  you can make the residue at the center of the window  have a
        larger weight than the others by setting the edge value for the  residues at
        the beginning and end of the interval to a value between 0 and 1. For
        instance, for Edge=0.4 and a window size of 5 the weights will be: 0.4, 0.7,
        1.0, 0.7, 0.4.  The method returns a list of values which can be plotted to
        view the change along a protein sequence.  Many scales exist. Just add your
        favorites to the ProtParamData modules.

        Similar to expasy's ProtScale: http://www.expasy.org/cgi-bin/protscale.pl
        """
        # generate the weights
        weight = self._weight_list(Window,Edge)
        list = []
        # the score in each Window is divided by the sum of weights
        sum_of_weights = 0.0
        for i in weight: sum_of_weights += i
        # since the weight list is one sided:
        sum_of_weights = sum_of_weights*2+1
        
        for i in range(self.length-Window+1):
            subsequence = self.sequence[i:i+Window]
            score = 0.0
            for j in range(Window//2):
                # walk from the outside of the Window towards the middle.
                # Iddo: try/except clauses added to avoid raising an exception on a non-standad amino acid
                    try:
                        score += weight[j] * ParamDict[subsequence[j]] + weight[j] * ParamDict[subsequence[Window-j-1]]
                    except KeyError:
                        sys.stderr.write('warning: %s or %s is not a standard amino acid.\n' %
                                 (subsequence[j],subsequence[Window-j-1]))

            # Now add the middle value, which always has a weight of 1.
            if subsequence[Window//2] in ParamDict:
                score += ParamDict[subsequence[Window//2]]
            else:
                sys.stderr.write('warning: %s  is not a standard amino acid.\n' % (subsequence[Window//2]))
        
            list.append(score/sum_of_weights)
        return list

    def isoelectric_point(self):
        """Calculate the isoelectric point.
        
        This method uses the module IsoelectricPoint to calculate the pI of a protein.
        """
        if not self.amino_acids_content:
            self.count_amino_acids()
        X = IsoelectricPoint.IsoelectricPoint(self.sequence, self.amino_acids_content)
        return X.pi()
        
    def secondary_structure_fraction (self):
        """Calculate fraction of helix, turn and sheet.
        
        This methods returns a list of the fraction of amino acids which tend
        to be in Helix, Turn or Sheet.
        
        Amino acids in helix: V, I, Y, F, W, L.
        Amino acids in Turn: N, P, G, S.
        Amino acids in sheet: E, M, A, L.
        
        Returns a tuple of three integers (Helix, Turn, Sheet).
        """
        if not self.amino_acids_percent:
            self.get_amino_acids_percent()
        Helix = self.amino_acids_percent['V'] + self.amino_acids_percent['I'] + self.amino_acids_percent['Y'] + self.amino_acids_percent['F'] + self.amino_acids_percent['W'] + self.amino_acids_percent['L']
        Turn = self.amino_acids_percent['N'] + self.amino_acids_percent['P'] + self.amino_acids_percent['G'] + self.amino_acids_percent['S']
        Sheet = self.amino_acids_percent['E'] + self.amino_acids_percent['M'] + self.amino_acids_percent['A'] + self.amino_acids_percent['L']
        return Helix, Turn, Sheet

