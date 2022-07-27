# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Calculate Composition, Transition and Distribution descriptors of a protein.

Descriptors are calculated as described in:

    *  Inna Dubchak, Ilya Muchink, Stephen R.Holbrook and Sung-Hou Kim.
       Prediction of protein folding class using global description of amino
       acid sequence. Proc.Natl. Acad.Sci.USA, 1995, 92, 8700-8704.
       https://doi.org/10.1073/PNAS.92.19.8700

    *  Inna Dubchak, Ilya Muchink, Christopher Mayor, Igor Dralyuk and Sung-Hou
       Kim. Recognition of a Protein Fold in the Context of the SCOP
       classification. Proteins: Structure, Function and Genetics,1999,35,401-407.
       https://doi.org/10.1002/(SICI)1097-0134(19990601)35:4<401::AID-PROT3>3.0.CO;2-K

"""

from typing import Union, List, Tuple, Iterable
from Bio.Seq import Seq
from Bio.Data.IUPACData import protein_letters


class CTD_Property:
    """Auxiliary class for CTD calculations.

    Parameters
    ----------
    :groups: A List of strings, representing the classifications of each amino
             acid. For example, the list: ["DEKNQR", "AGHPSTY", "CLVIMFW"] means
             that the amino acids "DEKNQR", "AGHPSTY" and "CLVIMFW" will be mapped
             to the numbers 1, 2 and 3 respectively.

    Methods
    -------
    :transitions(): returns an iterator of the possible transitions between the
                    categories, without duplication with the reverse transition.
                    For a class with 4 groups it would return, in this order:
                    "12", "13", "14", "23", "24" and "34".

    Examples
    --------
    >>> from Bio.SeqUtils.ctd import CTD_Property
    >>> mock_property = CTD_Property(["K","R", "ANCQGHILMFPSTWYVDE"])

    It implements the __getitem__ and __len__ methods as well.

    >>> print("Amino Acid K is in group:", mock_property["K"])
    Amino Acid K is in group: 1
    >>> print("Amino Acid A is in group:", mock_property["A"])
    Amino Acid A is in group: 3
    >>> len(mock_property)
    3

    Groups must contain all standard amino acids.

    >>> mock_property = CTD_Property(["K", "R", "ANCQGHILMFPSTWYVD"])
    Traceback (most recent call last):
    ...
    ValueError: Given groups do not contain all Amino Acids.

    """

    def __init__(self, groups):
        """Initialize the class."""
        given_aas = "".join(groups)
        if any([(aa not in given_aas) for aa in protein_letters]):
            raise ValueError("Given groups do not contain all Amino Acids.")

        self.d = {}
        for i, group in enumerate(groups):
            for aa in group:
                self.d[aa] = str(i + 1)

        self.size = len(groups)

    def __getitem__(self, key):
        return self.d[key]

    def __len__(self):
        return self.size

    def transitions(self):
        """Return Iterable of the possible transitions between groups."""
        for i in range(self.size):
            for j in range(i + 1, self.size):
                yield f"{i+1}{j+1}"


# All values taken from the previous cited works of Dubchak, I. and Muchink, I.
# '1'stand for Polar; '2'stand for Neutral, '3' stand for Hydrophobicity
_Hydrophobicity = CTD_Property(groups=["DEKNQR", "AGHPSTY", "CLVIMFW"])

# '1'stand for (0-2.78); '2'stand for (2.95-4.0), '3' stand for (4.03-8.08)
_NormalizedVDWV = CTD_Property(groups=["GASTPDC", "NVEQIL", "MHKFRYW"])

# '1'stand for (4.9-6.2); '2'stand for (8.0-9.2), '3' stand for (10.4-13.0)
_Polarity = CTD_Property(groups=["LIFWCMVY", "PATGS", "HQRKNED"])

# '1'stand for Positive; '2'stand for Neutral, '3' stand for Negative
_Charge = CTD_Property(groups=["KR", "ANCQGHILMFPSTWYV", "DE"])

# '1'stand for Helix; '2'stand for Strand, '3' stand for coil
_SecondaryStr = CTD_Property(groups=["EALMQKRH", "VIYCWFT", "GNPSD"])

# '1'stand for Buried; '2'stand for Exposed, '3' stand for Intermediate
_SolventAccessibility = CTD_Property(groups=["ALFCGIVW", "RKQEND", "MPSTHY"])

# '1'stand for (0-0.108); '2'stand for (0.128-0.186), '3' stand for (0.219-0.409)
_Polarizability = CTD_Property(groups=["GASDT", "CPNVEQIL", "KMHFRYW"])

default_ctd_props = [
    _Hydrophobicity,
    _NormalizedVDWV,
    _Charge,
    _Polarity,
    _Polarizability,
    _SecondaryStr,
    _SolventAccessibility,
]


class CTD:
    """A Class for calculating the CTD descriptors of a protein.

    The default amino acid properties used for calculating the descriptors are:
    Hidrophobicity, Normalized Van der Waals Volume, Polarity, Charge, Secondary
    Structure, Solvent Acessibility and Polarizability. All of the properties
    are given as 3 different groups of Amino Acids, based on their similarity
    in the given property.

    Parameters
    ----------
    :protein_sequence: A ``Bio.Seq`` or string object containing a protein
                       sequence.

    Methods
    -------
    :all_descriptors(properties): returns a list with the values of the CTD
                                  descriptors. The size of this list is defined
                                  by the number of groups of each property, with
                                  the defaults it is: len(properties) * 21.
    :composition_descriptors(properties): returns the Compositon descriptors of
                                          CTD. The size of this list is defined
                                          by the number of groups of each property,
                                          with the defaults it is: len(properties) * 3.
    :transition_descriptors(properties): returns the Transition descriptors of
                                         CTD. The size of this list is defined
                                         by the number of groups of each property,
                                         with the defaults it is:
                                         len(properties) * 3.
    :distribution_descriptors(properties): returns the Distribution descriptors of
                                           CTD. The size of this list is defined
                                           by the number of groups of each property,
                                           with the defaults it is: len(properties) * 15.

    Examples
    --------
    >>> from Bio.SeqUtils.ctd import CTD
    >>> protein = CTD("ACKLAA")
    >>> desc = protein.all_descriptors()
    >>> print(len(desc))
    147

    Calculate only a part of the descriptors.

    >>> C = protein.composition_descriptors()
    >>> T = protein.transition_descriptors()
    >>> D = protein.distribution_descriptors()
    >>> print(len(C), len(T), len(D))
    21 21 105

    Using user defined property groups.

    >>> from Bio.SeqUtils.ctd import CTD_Property
    >>> mock_property = CTD_Property(["K","R", "ANCQGHILMFPSTWYVDE"])
    >>> desc = protein.all_descriptors([mock_property])
    >>> print(len(desc))
    21

    To add your custom defined property to the standard one:

    >>> from Bio.SeqUtils.ctd import default_ctd_props
    >>> desc = protein.all_descriptors(default_ctd_props + [mock_property])
    >>> print(len(desc))
    168

    This methods can either be accessed from the class itself or from a
    ``ProtParam.ProteinAnalysis`` object:

    >>> from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
    >>> protein = PA("ACKLAA")
    >>> desc = protein.get_CTD_descriptors()
    >>> print(len(desc))
    147
    """

    def __init__(self, protein_sequence):
        """Initialize the class."""
        self.seq = str(protein_sequence).upper()
        self.L = len(self.seq)

    def _get_prop_tuples(self, properties):
        """Auxiliary function to convert between aa and properties groups (PRIVATE).

        Returns a list of (numberSeq, p) tuples, where numberSeq is the protein sequence
        with each amino acid converted to a str representation of it's group, and p is
        a CTD_Property object that contains the groups information.
        """
        return [("".join([p[aa] for aa in self.seq]), p) for p in properties]

    def _calc_C(self, prop_tuples):
        """Calculate the Composition descriptors of CTD (PRIVATE)."""
        C_descriptors = []
        for string, prop in prop_tuples:
            C_descriptors += [
                string.count(str(i + 1)) / self.L for i in range(len(prop))
            ]

        return C_descriptors

    def _calc_T(self, prop_tuples):
        """Calculate the Transition descriptors of CTD (PRIVATE)."""
        T_descriptors = []
        for string, prop in prop_tuples:
            T_descriptors += [
                (string.count(t) + string.count(t[::-1])) / (self.L - 1)
                for t in prop.transitions()
            ]

        return T_descriptors

    def _calc_D(self, prop_tuples):
        """Calculate the Distribution descriptors of CTD (PRIVATE)."""
        D_descriptors = []
        for string, prop in prop_tuples:
            positions: dict = {str(k + 1): [] for k in range(len(prop))}

            # Get the positions of each group in sequence
            for i, c in enumerate(string):
                positions[c].append(i)

            # calculate the descriptor values
            for key in positions:
                pos_list = positions[key]

                # If group is not present on string, all 5 values are 0
                if pos_list == []:
                    D_descriptors += [0.0 for _ in range(5)]
                else:
                    l = len(pos_list) / 5
                    for i in range(5):
                        D_descriptors.append(pos_list[int(l * i)] / self.L)

        return D_descriptors

    def composition_descriptors(self, properties=default_ctd_props):
        """Return the Composition descriptors of CTD.

        Returns the relative composition of each group of given properties
        in the sequence.
        """
        prop_tuples = self._get_prop_tuples(properties)

        return self._calc_C(prop_tuples)

    def transition_descriptors(self, properties=default_ctd_props):
        """Return the Transition descriptors of CTD.

        For each property, returns the relative frequency of each transition
        between the defined groups of that property in the sequence.
        """
        prop_tuples = self._get_prop_tuples(properties)

        return self._calc_T(prop_tuples)

    def distribution_descriptors(self, properties=default_ctd_props):
        """Return the Distribution descriptors of CTD.

        For each property, calculates the % of the sequence that contains
        0, 25, 50, 75 and 100% of the amino acids of each defined group.
        """
        prop_tuples = self._get_prop_tuples(properties)

        return self._calc_D(prop_tuples)

    def all_descriptors(self, properties=default_ctd_props):
        """Return the CTD descriptors, defined in Dubchak et al, 1995.

        Each property in _properties_, will add G*7 values to the final result,
        where G is the number of groups defined in the property, with G
        Composition descriptors, G*(G-1)/2 Transition descriptors and 5*G
        Distribution descriptors.

        With the default properties returns 147 values.
        """
        prop_tuples = self._get_prop_tuples(properties)

        return (
            self._calc_C(prop_tuples)
            + self._calc_T(prop_tuples)
            + self._calc_D(prop_tuples)
        )


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
