# Copyright Yair Benita Y.Benita@pharm.uu.nl
# Biopython (http://biopython.org) license applies

"""Simple protein analysis.

Example::

    X = ProteinAnalysis("MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV")
    print(X.count_amino_acids())
    print(X.get_amino_acids_percent())
    print(X.molecular_weight())
    print(X.aromaticity())
    print(X.instability_index())
    print(X.flexibility())
    print(X.isoelectric_point())
    print(X.secondary_structure_fraction())
    print(X.protein_scale(ProtParamData.kd, 9, 0.4))

"""

from __future__ import print_function

import sys
from Bio.SeqUtils import ProtParamData  # Local
from Bio.SeqUtils import IsoelectricPoint  # Local
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData
from Bio.SeqUtils import molecular_weight

__docformat__ = "restructuredtext en"


class ProteinAnalysis(object):
    """Class containing methods for protein analysis.

    The constructor takes two arguments.
    The first is the protein sequence as a string, which is then converted to a
    sequence object using the Bio.Seq module. This is done just to make sure
    the sequence is a protein sequence and not anything else.

    The second argument is optional. If set to True, the weight of the amino
    acids will be calculated using their monoisotopic mass (the weight of the
    most abundant isotopes for each element), instead of the average molecular
    mass (the averaged weight of all stable isotopes for each element).
    If set to false (the default value) or left out, the IUPAC average
    molecular mass will be used for the calculation.

    """
    def __init__(self, prot_sequence, monoisotopic=False):
        if prot_sequence.islower():
            self.sequence = Seq(prot_sequence.upper(), IUPAC.protein)
        else:
            self.sequence = Seq(prot_sequence, IUPAC.protein)
        self.amino_acids_content = None
        self.amino_acids_percent = None
        self.length = len(self.sequence)
        self.monoisotopic = monoisotopic

    def count_amino_acids(self):
        """Count standard amino acids, returns a dict.

        Counts the number times each amino acid is in the protein
        sequence. Returns a dictionary {AminoAcid:Number}.

        The return value is cached in self.amino_acids_content.
        It is not recalculated upon subsequent calls.
        """
        if self.amino_acids_content is None:
            prot_dic = dict((k, 0) for k in IUPACData.protein_letters)
            for aa in prot_dic:
                prot_dic[aa] = self.sequence.count(aa)

            self.amino_acids_content = prot_dic

        return self.amino_acids_content

    def get_amino_acids_percent(self):
        """Calculate the amino acid content in percentages.

        The same as count_amino_acids only returns the Number in percentage of
        entire sequence. Returns a dictionary of {AminoAcid:percentage}.

        The return value is cached in self.amino_acids_percent.

        input is the dictionary self.amino_acids_content.
        output is a dictionary with amino acids as keys.
        """
        if self.amino_acids_percent is None:
            aa_counts = self.count_amino_acids()

            percentages = {}
            for aa in aa_counts:
                percentages[aa] = aa_counts[aa] / float(self.length)

            self.amino_acids_percent = percentages

        return self.amino_acids_percent

    def molecular_weight(self):
        """Calculate MW from Protein sequence"""
        return molecular_weight(self.sequence, monoisotopic=self.monoisotopic)

    def aromaticity(self):
        """Calculate the aromaticity according to Lobry, 1994.

        Calculates the aromaticity value of a protein according to Lobry, 1994.
        It is simply the relative frequency of Phe+Trp+Tyr.
        """
        aromatic_aas = 'YWF'
        aa_percentages = self.get_amino_acids_percent()

        aromaticity = sum(aa_percentages[aa] for aa in aromatic_aas)

        return aromaticity

    def instability_index(self):
        """Calculate the instability index according to Guruprasad et al 1990.

        Implementation of the method of Guruprasad et al. 1990 to test a
        protein for stability. Any value above 40 means the protein is unstable
        (has a short half life).

        See: Guruprasad K., Reddy B.V.B., Pandit M.W.
        Protein Engineering 4:155-161(1990).
        """
        index = ProtParamData.DIWV
        score = 0.0

        for i in range(self.length - 1):
            this, next = self.sequence[i:i + 2]
            dipeptide_value = index[this][next]
            score += dipeptide_value

        return (10.0 / self.length) * score

    def flexibility(self):
        """Calculate the flexibility according to Vihinen, 1994.

        No argument to change window size because parameters are specific for a
        window=9. The parameters used are optimized for determining the flexibility.
        """
        flexibilities = ProtParamData.Flex
        window_size = 9
        weights = [0.25, 0.4375, 0.625, 0.8125, 1]
        scores = []

        for i in range(self.length - window_size):
            subsequence = self.sequence[i:i + window_size]
            score = 0.0

            for j in range(window_size // 2):
                front = subsequence[j]
                back = subsequence[window_size - j - 1]
                score += (flexibilities[front] + flexibilities[back]) * weights[j]

            middle = subsequence[window_size // 2 + 1]
            score += flexibilities[middle]

            scores.append(score / 5.25)

        return scores

    def gravy(self):
        """Calculate the gravy according to Kyte and Doolittle."""
        total_gravy = sum(ProtParamData.kd[aa] for aa in self.sequence)

        return total_gravy / self.length

    def _weight_list(self, window, edge):
        """Makes a list of relative weight of the
        window edges compared to the window center. The weights are linear.
        it actually generates half a list. For a window of size 9 and edge 0.4
        you get a list of [0.4, 0.55, 0.7, 0.85].
        """
        unit = 2 * (1.0 - edge) / (window - 1)
        weights = [0.0] * (window // 2)

        for i in range(window // 2):
            weights[i] = edge + unit * i

        return weights

    def protein_scale(self, param_dict, window, edge=1.0):
        """Compute a profile by any amino acid scale.

        An amino acid scale is defined by a numerical value assigned to each type of
        amino acid. The most frequently used scales are the hydrophobicity or
        hydrophilicity scales and the secondary structure conformational parameters
        scales, but many other scales exist which are based on different chemical and
        physical properties of the amino acids.  You can set several parameters that
        control the computation  of a scale profile, such as the window size and the
        window edge relative weight value.

        WindowSize: The window size is the length
        of the interval to use for the profile computation. For a window size n, we
        use the i-(n-1)/2 neighboring residues on each side to compute
        the score for residue i. The score for residue i is the sum of the scaled values
        for these amino acids, optionally weighted according to their position in the
        window.

        Edge: The central amino acid of the window always has a weight of 1.
        By default, the amino acids at the remaining window positions have the same
        weight, but you can make the residue at the center of the window  have a
        larger weight than the others by setting the edge value for the  residues at
        the beginning and end of the interval to a value between 0 and 1. For
        instance, for Edge=0.4 and a window size of 5 the weights will be: 0.4, 0.7,
        1.0, 0.7, 0.4.

        The method returns a list of values which can be plotted to
        view the change along a protein sequence.  Many scales exist. Just add your
        favorites to the ProtParamData modules.

        Similar to expasy's ProtScale: http://www.expasy.org/cgi-bin/protscale.pl
        """
        # generate the weights
        #   _weight_list returns only one tail. If the list should be [0.4,0.7,1.0,0.7,0.4]
        #   what you actually get from _weights_list is [0.4,0.7]. The correct calculation is done
        #   in the loop.
        weights = self._weight_list(window, edge)
        scores = []

        # the score in each Window is divided by the sum of weights
        # (* 2 + 1) since the weight list is one sided:
        sum_of_weights = sum(weights) * 2 + 1

        for i in range(self.length - window + 1):
            subsequence = self.sequence[i:i + window]
            score = 0.0

            for j in range(window // 2):
                # walk from the outside of the Window towards the middle.
                # Iddo: try/except clauses added to avoid raising an exception on a non-standard amino acid
                try:
                    front = param_dict[subsequence[j]]
                    back = param_dict[subsequence[window - j - 1]]
                    score += weights[j] * front + weights[j] * back
                except KeyError:
                    sys.stderr.write('warning: %s or %s is not a standard amino acid.\n' %
                             (subsequence[j], subsequence[window - j - 1]))

            # Now add the middle value, which always has a weight of 1.
            middle = subsequence[window // 2]
            if middle in param_dict:
                score += param_dict[middle]
            else:
                sys.stderr.write('warning: %s  is not a standard amino acid.\n' % (middle))

            scores.append(score / sum_of_weights)

        return scores

    def isoelectric_point(self):
        """Calculate the isoelectric point.

        Uses the module IsoelectricPoint to calculate the pI of a protein.
        """
        aa_content = self.count_amino_acids()

        ie_point = IsoelectricPoint.IsoelectricPoint(self.sequence, aa_content)
        return ie_point.pi()

    def secondary_structure_fraction(self):
        """Calculate fraction of helix, turn and sheet.

        Returns a list of the fraction of amino acids which tend
        to be in Helix, Turn or Sheet.

        Amino acids in helix: V, I, Y, F, W, L.
        Amino acids in Turn: N, P, G, S.
        Amino acids in sheet: E, M, A, L.

        Returns a tuple of three integers (Helix, Turn, Sheet).
        """
        aa_percentages = self.get_amino_acids_percent()

        helix = sum(aa_percentages[r] for r in 'VIYFWL')
        turn = sum(aa_percentages[r] for r in 'NPGS')
        sheet = sum(aa_percentages[r] for r in 'EMAL')

        return helix, turn, sheet
