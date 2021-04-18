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
       classification. Proteins: Structure, Function and
       Genetics,1999,35,401-407.

"""

# https://rdrr.io/bioc/Rcpi/src/R/409-extractProtCTDD.R
# https://rdrr.io/bioc/Rcpi/src/R/409-extractProtCTDT.R

from typing import Union, List, Tuple
from Bio.Seq import Seq
from Bio.Data.IUPACData import protein_letters

# TODO: references for groups
# TODO: add examples to CTD class
# TODO: improve all function docstrings
# TODO: add function to ProteinAnalysis


class CTD_Property:
    """Auxiliary class for CTD calculations.

    Parameters
    ----------
    :protein_sequence: A ``Bio.Seq`` or string object containing a protein
                       sequence.

    Methods
    -------
    :composition_descriptors(properties): returns the Compositon descriptors of CT
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

_default_ctd_props = [
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
    :composition_descriptors(properties): returns the Compositon descriptors of CTD.
    :transition_descriptors(properties): returns the Transition descriptors of CTD.
    :distribution_descriptors(properties): returns the Distribution descriptors of CTD.
    :CTD_descriptors(properties): returns the CTD descriptors.

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

    def _calc_C(self, prop_tuples) -> List[float]:
        """Calculate the Composition descriptors of CTD (PRIVATE)."""
        C_descriptors = []
        for string, prop in prop_tuples:
            C_descriptors += [string.count(str(i)) / self.L for i in range(len(prop))]

        return C_descriptors

    def _calc_T(self, prop_tuples) -> List[float]:
        """Calculate the Transition descriptors of CTD (PRIVATE)."""
        T_descriptors = []
        for string, prop in prop_tuples:
            T_descriptors += [
                (string.count(t) + string.count(t[::-1])) / (self.L - 1)
                for t in prop.transitions()
            ]

        return T_descriptors

    def _calc_D(self, prop_tuples) -> List[float]:
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
        self, properties: List[CTD_Property] = _default_ctd_props
    ) -> List[float]:
        """Return the Composition descriptors of CTD."""
        prop_tuples = self._seq2number_strings(properties)

        return self._calc_T(prop_tuples)

    def transition_descriptors(
        self, properties: List[CTD_Property] = _default_ctd_props
    ) -> List[float]:
        """Return the Transition descriptors of CTD."""
        prop_tuples = self._seq2number_strings(properties)

        return self._calc_T(prop_tuples)

    def distribution_descriptors(
        self, properties: List[CTD_Property] = _default_ctd_props
    ) -> List[float]:
        """Return the Distribution descriptors of CTD."""
        prop_tuples = self._seq2number_strings(properties)

        return self._calc_D(prop_tuples)

    def CTD_descriptors(
        self, properties: List[CTD_Property] = _default_ctd_props
    ) -> List[float]:
        """Return the CTD descriptors."""
        # Taking N as = len(properties)
        # First N*3 values are the composition of each group type
        # The next N*3 values are the Transition properties
        # The final N*5*3 values are the Distribution properties
        prop_tuples = self._get_prop_tuples(properties)

        return (
            self._calc_C(prop_tuples)
            + self._calc_T(prop_tuples)
            + self._calc_D(prop_tuples)
        )
