# Copyright 2003 Yair Benita.  All rights reserved.
# Revisions copyright 2020 by Tianyi Shi.  All rights reserved.
# Isoelectric point extension 2021 by Lukasz P. Kozlowski. Public Domain.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package. If you find this LICENSE too restrictive you can use native code as
# last revision is based on IPC (https://www.isoelectric.org) and IPC2
# (http://ipc2-isoelectric-point.org) that are within Public Domain.
"""Calculate isoelectric points of polypeptides with 18 different methods.

pK values and the methods are taken from:

    * Kozlowski, L.P. (2021) IPC 2.0: prediction of isoelectric point and pKa
      dissociation constants. Nucleic Acids Res. 49 (W1): W285–W292
      DOI: 10.1093/nar/gkab295


Method       Cterm Asp   Glu   Cys   Tyr    Nterm  Lys    Arg    His
IPC2_protein 6.065 3.766 4.497 7.89  11.491  5.779  9.247 10.223 5.492
IPC2_peptide 2.977 3.969 4.504 9.454  9.153  7.947  8.165 11.493 6.439
IPC_protein  2.869 3.872 4.412 7.555 10.85   9.094  9.052 11.84  5.637
IPC_peptide  2.383 3.887 4.317 8.297 10.071  9.564 10.517 12.503 6.018
EMBOSS       3.6   3.9   4.1   8.5   10.1    8.6   10.8   12.5   6.5
DTASelect    3.1   4.4   4.4   8.5   10.0    8.0   10.0   12.0   6.5
Solomon      2.4   3.9   4.3   8.3   10.1    9.6   10.5   12.5   6.0
Sillero      3.2   4.0   4.5   9.0   10.0    8.2   10.4   12.0   6.4
Rodwell      3.1   3.68  4.25  8.33  10.07   8.0   11.5   11.5   6.0
Patrickios   4.2   4.2   4.2   0.0    0.0   11.2   11.2   11.2   0.0
Wikipedia    3.65  3.9   4.07  8.18  10.46   8.2   10.54  12.48  6.04
Grimsley     3.3   3.5   4.2   6.8   10.3    7.7   10.5   12.04  6.6
Lehninger    2.34  3.86  4.25  8.33  10.0    9.69  10.5   12.4   6.0
Bjellqvist   3.55  4.05  4.45  9.0   10.0    7.5   10.0   12.0   5.98
Toseland     3.19  3.6   4.29  6.87   9.61   8.71  10.45  12.0   6.33
Thurlkill    3.67  3.67  4.25  8.55   9.84   8.0   10.4   12.0   6.54
Nozaki       3.8   4.0   4.4   9.5    9.6    7.5   10.4   12.0   6.3
Dawson       3.2   3.9   4.3   8.3   10.1    8.2   10.5   12.0   6.0


    For references of individual methods listed in the table check:
    http://ipc2-isoelectric-point.org/theory.html

"""
pKa_scales = {
    "IPC2_peptide": [
        {"Cterm": 2.977, "D": 3.969, "E": 4.504, "C": 9.454, "Y": 9.153},
        {"H": 6.439, "Nterm": 7.947, "K": 8.165, "R": 11.493},
    ],
    "IPC2_protein": [
        {"Cterm": 6.065, "D": 3.766, "E": 4.497, "C": 7.89, "Y": 11.491},
        {"H": 5.492, "Nterm": 5.779, "K": 9.247, "R": 10.223},
    ],
    "IPC_peptide": [
        {"Cterm": 2.383, "D": 3.887, "E": 4.317, "C": 8.297, "Y": 10.071},
        {"H": 6.018, "Nterm": 9.564, "K": 10.517, "R": 12.503},
    ],
    "IPC_protein": [
        {"Cterm": 2.869, "D": 3.872, "E": 4.412, "C": 7.555, "Y": 10.85},
        {"H": 5.637, "Nterm": 9.094, "K": 9.052, "R": 11.84},
    ],
    "EMBOSS": [
        {"Cterm": 3.6, "D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1},
        {"H": 6.5, "Nterm": 8.6, "K": 10.8, "R": 12.5},
    ],
    "DTASelect": [
        {"Cterm": 3.1, "D": 4.4, "E": 4.4, "C": 8.5, "Y": 10.0},
        {"H": 6.5, "Nterm": 8.0, "K": 10.0, "R": 12.0},
    ],
    "Solomon": [
        {"Cterm": 2.4, "D": 3.9, "E": 4.3, "C": 8.3, "Y": 10.1},
        {"H": 6.0, "Nterm": 9.6, "K": 10.5, "R": 12.5},
    ],
    "Sillero": [
        {"Cterm": 3.2, "D": 4.0, "E": 4.5, "C": 9.0, "Y": 10.0},
        {"H": 6.4, "Nterm": 8.2, "K": 10.4, "R": 12.0},
    ],
    "Rodwell": [
        {"Cterm": 3.1, "D": 3.68, "E": 4.25, "C": 8.33, "Y": 10.07},
        {"H": 6.0, "Nterm": 8.0, "K": 11.5, "R": 11.5},
    ],
    "Patrickios": [
        {"Cterm": 4.2, "D": 4.2, "E": 4.2, "C": 0.0, "Y": 0.0},
        {"H": 0.0, "Nterm": 11.2, "K": 11.2, "R": 11.2},
    ],
    "Wikipedia": [
        {"Cterm": 3.65, "D": 3.9, "E": 4.07, "C": 8.18, "Y": 10.46},
        {"H": 6.04, "Nterm": 8.2, "K": 10.54, "R": 12.48},
    ],
    "Grimsley": [
        {"Cterm": 3.3, "D": 3.5, "E": 4.2, "C": 6.8, "Y": 10.3},
        {"H": 6.6, "Nterm": 7.7, "K": 10.5, "R": 12.04},
    ],
    "Lehninger": [
        {"Cterm": 2.34, "D": 3.86, "E": 4.25, "C": 8.33, "Y": 10.0},
        {"H": 6.0, "Nterm": 9.69, "K": 10.5, "R": 12.4},
    ],
    "Bjellqvist": [
        {"Cterm": 3.55, "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0},
        {"H": 5.98, "Nterm": 7.5, "K": 10.0, "R": 12.0},
    ],
    "Toseland": [
        {"Cterm": 3.19, "D": 3.6, "E": 4.29, "C": 6.87, "Y": 9.61},
        {"H": 6.33, "Nterm": 8.71, "K": 10.45, "R": 12},
    ],
    "Thurlkill": [
        {"Cterm": 3.67, "D": 3.67, "E": 4.25, "C": 8.55, "Y": 9.84},
        {"H": 6.54, "Nterm": 8.0, "K": 10.4, "R": 12.0},
    ],
    "Nozaki": [
        {"Cterm": 3.8, "D": 4.0, "E": 4.4, "C": 9.5, "Y": 9.6},
        {"H": 6.3, "Nterm": 7.5, "K": 10.4, "R": 12},
    ],
    "Dawson": [
        {"Cterm": 3.2, "D": 3.9, "E": 4.3, "C": 8.3, "Y": 10.1},
        {"H": 6.0, "Nterm": 8.2, "K": 10.5, "R": 12},
    ],
}

positive_pKs = {"Nterm": 7.5, "K": 10.0, "R": 12.0, "H": 5.98}
negative_pKs = {"Cterm": 3.55, "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
pKcterminal = {"D": 4.55, "E": 4.75}
pKnterminal = {
    "A": 7.59,
    "M": 7.0,
    "S": 6.93,
    "P": 8.36,
    "T": 6.82,
    "V": 7.44,
    "E": 7.7,
}
charged_aas = ("K", "R", "H", "D", "E", "C", "Y")


class IsoelectricPoint:
    """A class for calculating the IEP or charge at given pH of a protein.

    Parameters
    ----------
    :protein_sequence: A ``Bio.Seq`` or string object containing a protein
                       sequence.
    :aa_content: A dictionary with amino acid letters as keys and its
                 occurrences as integers, e.g. ``{"A": 3, "C": 0, ...}``.
                 Default: ``None``. If ``None``, the dic will be calculated
                 from the given sequence.
    :pKa_scale: Different pKa scales. Need to be one of: Bjellqvist, DTASelect,
            Dawson, EMBOSS, Grimsley, IPC2_peptide, IPC2_protein,
            IPC_peptide, IPC_protein, Lehninger, Nozaki, Patrickios, ProMoST,
            Rodwell, Sillero, Solomon, Thurlkill, Toseland, Wikipedia.
            Default scale is 'IPC2_protein' (we assume the protein as an input)

    Methods
    -------
    :charge_at_pH(pH):  Calculates the charge of the protein for a given pH
    :pi():              Calculates the isoelectric point


    Examples
    --------
    The methods of this class can either be accessed from the class itself
    or from a ``ProtParam.ProteinAnalysis`` object (with partially different
    names):

    >>> from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
    >>> protein = IP("ADDKNPLEECFRETDYEEFLEIARNGLKATSNPKRVV")
    >>> print(f"IEP of protein {protein.sequence} is {protein.pi():.2f}")
    IEP of protein ADDKNPLEECFRETDYEEFLEIARNGLKATSNPKRVV is 4.83
    >>> print(f"Its charge at pH 7.4 (cytoplasm) is {protein.charge_at_pH(7.4):.2f}")
    Its charge at pH 7.4 (cytoplasm) is -4.22

    >>> peptide = IP("EFYTVDGQK", "IPC2_peptide")
    >>> print(f"IEP of peptide {peptide.sequence} is {peptide.pi():.2f}")
    IEP of peptide EFYTVDGQK is 4.28
    >>> print(f"Its charge at pH 7.4 is {peptide.charge_at_pH(7.4):.2f}")
    Its charge at pH 7.4 is -1.38

    >>> from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
    >>> peptide = PA("EFYTVDGQK")
    >>> print(f"IEP of {peptide.sequence} is {peptide.isoelectric_point('Toseland'):.2f}")
    IEP of EFYTVDGQK is 4.06
    >>> print(f"Charge at pH 7.4 (cytoplasm) is {peptide.charge_at_pH(7.4, 'Toseland'):.2f}")
    Charge at pH 7.4 is -1.05

    """

    def __init__(self, aa_sequence, pKa_scale=None, aa_content=None):
        """Initialize the class."""
        self.sequence = str(aa_sequence).upper()
        if not aa_content:
            from Bio.SeqUtils.ProtParam import ProteinAnalysis as _PA

            aa_content = _PA(self.sequence).count_amino_acids()

        if not pKa_scale:
            self.pKa_scale = "IPC2_protein"
        else:
            self.pKa_scale = str(pKa_scale)
            if self.pKa_scale not in pKa_scales:
                self.pKa_scale = "IPC2_protein"

        self.charged_aas_content = self._select_charged(aa_content)

        self.pos_pKs, self.neg_pKs = self._update_pKs_tables()

    # This function creates a dictionary with the contents of each charged aa,
    # plus Cterm and Nterm.
    def _select_charged(self, aa_content):
        charged = {}
        for aa in charged_aas:
            charged[aa] = float(aa_content[aa])
        charged["Nterm"] = 1.0
        charged["Cterm"] = 1.0
        return charged

    def _update_pKs_tables(self):
        """Update pKs tables with seq specific values for N- and C-termini."""
        if self.pKa_scale == "Bjellqvist":
            neg_pKs = pKa_scales[self.pKa_scale][0].copy()
            pos_pKs = pKa_scales[self.pKa_scale][1].copy()
            nterm, cterm = self.sequence[0], self.sequence[-1]
            if nterm in pKnterminal:
                pos_pKs["Nterm"] = pKnterminal[nterm]
            if cterm in pKcterminal:
                neg_pKs["Cterm"] = pKcterminal[cterm]
            return pos_pKs, neg_pKs
        else:
            neg_pKs = pKa_scales[self.pKa_scale][0].copy()
            pos_pKs = pKa_scales[self.pKa_scale][1].copy()
            return pos_pKs, neg_pKs

    def charge_at_pH(self, pH):
        """Calculate the charge of a protein at given pH."""
        # derivation:
        #   Henderson Hasselbalch equation: pH = pKa + log([A-]/[HA])
        #   Rearranging: [HA]/[A-] = 10 ** (pKa - pH)
        #   partial_charge =
        #       [A-]/[A]total = [A-]/([A-] + [HA]) = 1 / { ([A-] + [HA])/[A-] } =
        #       1 / (1 + [HA]/[A-]) = 1 / (1 + 10 ** (pKa - pH)) for acidic residues;
        #                             1 / (1 + 10 ** (pH - pKa)) for basic residues
        positive_charge = 0.0
        for aa, pK in self.pos_pKs.items():
            partial_charge = 1.0 / (10 ** (pH - pK) + 1.0)
            positive_charge += self.charged_aas_content[aa] * partial_charge

        negative_charge = 0.0
        for aa, pK in self.neg_pKs.items():
            partial_charge = 1.0 / (10 ** (pK - pH) + 1.0)
            negative_charge += self.charged_aas_content[aa] * partial_charge

        return positive_charge - negative_charge

    # This is the action function, it tries different pH until the charge of
    # the protein is 0 (or close).
    def pi(self, pH=7.5, min_=1.0, max_=14):
        r"""Calculate and return the isoelectric point as float.

        This is a recursive function that uses bisection method.
        Wiki on bisection: https://en.wikipedia.org/wiki/Bisection_method

        Arguments:
         - pH: the pH at which the current charge of the protein is computed.
           This pH lies at the centre of the interval (mean of `min_` and `max_`).
         - min\_: the minimum of the interval. Initial value set to 1.0 which is
           below the theoretical maximum.
         - max\_: the maximum of the the interval. Initial value defaults to 14,
           which is above the theoretical maximum.
        """
        charge = self.charge_at_pH(pH)
        if max_ - min_ > 0.001:
            if charge > 0.0:
                min_ = pH
            else:
                max_ = pH
            next_pH = (min_ + max_) / 2
            return self.pi(next_pH, min_, max_)
        return pH
