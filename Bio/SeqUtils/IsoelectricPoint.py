# Copyright 2003 Yair Benita.  All rights reserved.
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


class IsoelectricPoint(object):
    """A class for calculating the IEP or charge at given pH of a protein.

    Parameters
    ----------
    :protein_sequence: A ``Bio.Seq`` or string object containing a protein
                       sequence.
    :aa_content: A dictionary with amino acid letters as keys and it's
                 occurences as integers, e.g. ``{"A": 3, "C": 0, ...}``.
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
    >>> print("IEP of peptide {} is {:.2f}"
    ...       .format(protein.sequence, protein.pi()))
    IEP of peptide INGAR is 9.75
    >>> print("It's charge at pH 7 is {:.2f}"
    ...       .format(protein.charge_at_pH(7.0)))
    It's charge at pH 7 is 0.76


    >>> from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
    >>> protein = PA("PETER")
    >>> print("IEP of {}: {:.2f}".format(protein.sequence,
    ...                                  protein.isoelectric_point()))
    IEP of PETER: 4.53
    >>> print("Charge at pH 4.53: {:.2f}"
    ...       .format(protein.charge_at_pH(4.53)))
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

    # This function calculates the total charge of the protein at a given pH.
    def _chargeR(self, pH):
        positive_charge = 0.0
        for aa, pK in self.pos_pKs.items():
            CR = 10 ** (pK - pH)
            partial_charge = CR / (CR + 1.0)
            positive_charge += self.charged_aas_content[aa] * partial_charge

        negative_charge = 0.0
        for aa, pK in self.neg_pKs.items():
            CR = 10 ** (pH - pK)
            partial_charge = CR / (CR + 1.0)
            negative_charge += self.charged_aas_content[aa] * partial_charge

        return positive_charge - negative_charge

    def charge_at_pH(self, pH):
        """Calculate the charge of a protein at given pH."""
        return self._chargeR(pH)

    # This is the action function, it tries different pH until the charge of
    # the protein is 0 (or close).
    def pi(self):
        """Calculate and return the isoelectric point as float."""
        # Bracket between pH1 and pH2
        pH = 7
        charge = self._chargeR(pH)
        if charge > 0.0:
            pH1 = pH
            charge1 = charge
            while charge1 > 0.0:
                pH = pH1 + 1.0
                charge = self._chargeR(pH)
                if charge > 0.0:
                    pH1 = pH
                    charge1 = charge
                else:
                    pH2 = pH
                    charge2 = charge
                    break
        else:
            pH2 = pH
            charge2 = charge
            while charge2 < 0.0:
                pH = pH2 - 1.0
                charge = self._chargeR(pH)
                if charge < 0.0:
                    pH2 = pH
                    charge2 = charge
                else:
                    pH1 = pH
                    charge1 = charge
                    break

        # Bisection
        while pH2 - pH1 > 0.0001 and charge != 0.0:
            pH = (pH1 + pH2) / 2.0
            charge = self._chargeR(pH)
            if charge > 0.0:
                pH1 = pH
                charge1 = charge
            else:
                pH2 = pH
                charge2 = charge

        return pH
