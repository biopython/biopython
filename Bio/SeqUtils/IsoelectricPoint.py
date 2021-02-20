# Copyright 2003 Yair Benita.  All rights reserved.
# Revisions copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Calculate isoelectric points of polypeptides using methods of Bjellqvist.

pK values and the methos are taken from::

    * Bjellqvist, B.,Hughes, G.J., Pasquali, Ch., Paquet, N., Ravier, F.,
    Sanchez, J.-Ch., Frutiger, S. & Hochstrasser, D.F.
    The focusing positions of polypeptides in immobilized pH gradients can be
    predicted from their amino acid sequences. Electrophoresis 1993, 14,
    1023-1031.

    * Bjellqvist, B., Basse, B., Olsen, E. and Celis, J.E.
    Reference points for comparisons of two-dimensional maps of proteins from
    different human cell types defined in a pH scale where isoelectric points
    correlate with polypeptide compositions. Electrophoresis 1994, 15, 529-539.

I designed the algorithm according to a note by David L. Tabb, available at:
http://fields.scripps.edu/DTASelect/20010710-pI-Algorithm.pdf
"""

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
    >>> protein = IP("INGAR")
    >>> print(f"IEP of peptide {protein.sequence} is {protein.pi():.2f}")
    IEP of peptide INGAR is 9.75
    >>> print(f"Its charge at pH 7 is {protein.charge_at_pH(7.0):.2f}")
    Its charge at pH 7 is 0.76


    >>> from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
    >>> protein = PA("PETER")
    >>> print(f"IEP of {protein.sequence}: {protein.isoelectric_point():.2f}")
    IEP of PETER: 4.53
    >>> print(f"Charge at pH 4.53: {protein.charge_at_pH(4.53):.2f}")
    Charge at pH 4.53: 0.00

    """

    def __init__(self, protein_sequence, aa_content=None):
        """Initialize the class."""
        self.sequence = str(protein_sequence).upper()
        if not aa_content:
            from Bio.SeqUtils.ProtParam import ProteinAnalysis as _PA

            aa_content = _PA(self.sequence).count_amino_acids()
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
        pos_pKs = positive_pKs.copy()
        neg_pKs = negative_pKs.copy()
        nterm, cterm = self.sequence[0], self.sequence[-1]
        if nterm in pKnterminal:
            pos_pKs["Nterm"] = pKnterminal[nterm]
        if cterm in pKcterminal:
            neg_pKs["Cterm"] = pKcterminal[cterm]
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
    def pi(self, pH=7.775, min_=4.05, max_=12):
        r"""Calculate and return the isoelectric point as float.

        This is a recursive function that uses bisection method.
        Wiki on bisection: https://en.wikipedia.org/wiki/Bisection_method

        Arguments:
         - pH: the pH at which the current charge of the protein is computed.
           This pH lies at the centre of the interval (mean of `min_` and `max_`).
         - min\_: the minimum of the interval. Initial value defaults to 4.05,
           which is below the theoretical minimum, when the protein is composed
           exclusively of aspartate.
         - max\_: the maximum of the the interval. Initial value defaults to 12,
           which is above the theoretical maximum, when the protein is composed
           exclusively of arginine.
        """
        charge = self.charge_at_pH(pH)
        if max_ - min_ > 0.0001:
            if charge > 0.0:
                min_ = pH
            else:
                max_ = pH
            next_pH = (min_ + max_) / 2
            return self.pi(next_pH, min_, max_)
        return pH
