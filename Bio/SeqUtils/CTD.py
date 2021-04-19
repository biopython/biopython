# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Calculate Composition, Transition and Distribution descriptors of a protein.

Descriptors are calculated as described in:

    *  Inna Dubchak, Ilya Muchink, Stephen R.Holbrook and Sung-Hou Kim.
       Prediction of protein folding class using global description of amino
       acid sequence. Proc.Natl. Acad.Sci.USA, 1995, 92, 8700-8704.

    *  Inna Dubchak, Ilya Muchink, Christopher Mayor, Igor Dralyuk and Sung-Hou
       Kim. Recognition of a Protein Fold in the Context of the SCOP
       classification. Proteins: Structure, Function and Genetics,1999,35,401-407.

"""

from typing import Union, List, Tuple
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
    >>> from Bio.SeqUtils.CTD import CTD_Property
    >>> mock_property = CTD_Property(["K","R", "ANCQGHILMFPSTWYVDE"])

    It implements the __getitem__ and __len__ methods as well.

    >>> mock_property["K"]
    "1"
    >>> mock_property["A"]
    "3"
    >>> len(mock_property)
    3

    """

    def __init__(self, groups: List[str]):
        """Initialize the class."""
        # TODO: assert all aa in groups
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
        """Return an iterator of the possible transitions."""
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

    Parameters
    ----------
    :protein_sequence: A ``Bio.Seq`` or string object containing a protein
                       sequence.

    Methods
    -------
    :CTD_descriptors(properties): returns a list with the values of the CTD
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
    >>> from Bio.SeqUtils.CTD import CTD
    >>> protein = CTD("ACKLAA")
    >>> ctd = protein.CTD_descriptors()
    >>> print(len(ctd))
    147

    This methods can either be accessed from the class itself or from a
    ``ProtParam.ProteinAnalysis`` object:

    >>> from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
    >>> protein = PA("ACKLAA")
    >>> ctd = protein.CTD_descriptors()
    147

    Calculate only a part of the descriptors.

    >>> C = protein.composition_descriptors()
    >>> T = protein.transition_descriptors()
    >>> D = protein.distribution_descriptors()
    >>> print(len(C), len(T), len(D))
    21 21 105

    Using user defined property groups.

    >>> from Bio.SeqUtils.CTD import CTD_Property
    >>> mock_property = CTD_Property(["K","R", "ANCQGHILMFPSTWYVDE"])
    >>> ctd = protein.CTD_descriptors([mock_property])
    >>> print(len(ctd))
    21

    To add your custom defined property to the standard one:

    >>> from Bio.SeqUtils.CTD import default_ctd_props
    >>> ctd = protein.CTD_descriptors(default_ctd_props + [mock_property])
    >>> print(len(ctd))
    168
    """

    def __init__(self, protein_sequence: Union[Seq, str]):
        """Initialize the class."""
        self.seq = str(protein_sequence).upper()
        self.L = len(self.seq)

    def _get_prop_tuples(
        self, properties: List[CTD_Property]
    ) -> List[Tuple[str, CTD_Property]]:
        """Auxiliary function to convert between aa and properties groups (PRIVATE).

        Returns a list of (numberSeq, p) tuples, where numberSeq is the protein sequence
        with each amino acid converted to a str representation of it's group, and p is
        a CTD_Property object that contains the groups information.
        """
        return [("".join([p[aa] for aa in self.seq]), p) for p in properties]

    def _calc_C(self, prop_tuples: list) -> List[float]:
        """Calculate the Composition descriptors of CTD (PRIVATE)."""
        C_descriptors = []
        for string, prop in prop_tuples:
            C_descriptors += [string.count(str(i)) / self.L for i in range(len(prop))]

        return C_descriptors

    def _calc_T(self, prop_tuples: list) -> List[float]:
        """Calculate the Transition descriptors of CTD (PRIVATE)."""
        T_descriptors = []
        for string, prop in prop_tuples:
            T_descriptors += [
                (string.count(t) + string.count(t[::-1])) / (self.L - 1)
                for t in prop.transitions()
            ]

        return T_descriptors

    def _calc_D(self, prop_tuples: list) -> List[float]:
        """Calculate the Distribution descriptors of CTD (PRIVATE)."""
        # TODO: probably a faster way of calculating this
        D_descriptors = []
        for s, p in prop_tuples:
            positions = {str(k + 1): [] for k in range(len(p))}

            # Get the positions of each aa
            for i, c in enumerate(s):
                positions[c].append(i)

            # calculate the descriptor values
            for key in positions:
                pos_list = positions[key]
                if pos_list == []:
                    D_descriptors += [0 for _ in range(5)]
                else:
                    l = len(pos_list) / 5
                    for i in range(5):
                        D_descriptors.append(pos_list[int(l * i)] / self.L)

        return D_descriptors

    def composition_descriptors(
        self, properties: List[CTD_Property] = default_ctd_props
    ) -> List[float]:
        """Return the Composition descriptors of CTD."""
        prop_tuples = self._get_prop_tuples(properties)

        return self._calc_T(prop_tuples)

    def transition_descriptors(
        self, properties: List[CTD_Property] = default_ctd_props
    ) -> List[float]:
        """Return the Transition descriptors of CTD."""
        prop_tuples = self._get_prop_tuples(properties)

        return self._calc_T(prop_tuples)

    def distribution_descriptors(
        self, properties: List[CTD_Property] = default_ctd_props
    ) -> List[float]:
        """Return the Distribution descriptors of CTD."""
        prop_tuples = self._get_prop_tuples(properties)

        return self._calc_D(prop_tuples)

    def CTD_descriptors(
        self, properties: List[CTD_Property] = default_ctd_props
    ) -> List[float]:
        """Return the CTD descriptors.

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
