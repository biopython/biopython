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

from typing import Union
from Bio.Seq import Seq
from Bio.Data.IUPACData import protein_letters

# '1'stand for Polar; '2'stand for Neutral, '3' stand for Hydrophobicity
_Hydrophobicity = dict(
    [(aa, "1") for aa in "RKEDQN"]
    + [(aa, "2") for aa in "GASTPHY"]
    + [(aa, "3") for aa in "CLVIMFW"]
)

# '1'stand for (0-2.78); '2'stand for (2.95-4.0), '3' stand for (4.03-8.08)
_NormalizedVDWV = dict(
    [(aa, "1") for aa in "GASTPD"]
    + [(aa, "2") for aa in "NVEQIL"]
    + [(aa, "3") for aa in "MHKFRYW"]
)

# '1'stand for (4.9-6.2); '2'stand for (8.0-9.2), '3' stand for (10.4-13.0)
_Polarity = dict(
    [(aa, "1") for aa in "LIFWCMVY"]
    + [(aa, "2") for aa in "PATGS"]
    + [(aa, "3") for aa in "HQRKNED"]
)
# 'L', 'I', 'F', 'W', 'C', 'M', 'V', 'Y'
# 'P', 'A', 'T', 'G', 'S'
# 'H', 'Q', 'R', 'K', 'N', 'E', 'D'
# '1'stand for Positive; '2'stand for Neutral, '3' stand for Negative
_Charge = dict(
    [(aa, "1") for aa in "KR"]
    + [(aa, "2") for aa in "ANCQGHILMFPSTWYV"]
    + [(aa, "3") for aa in "DE"]
)

# '1'stand for Helix; '2'stand for Strand, '3' stand for coil
_SecondaryStr = dict(
    [(aa, "1") for aa in "EALMQKRH"]
    + [(aa, "2") for aa in "VIYCWFT"]
    + [(aa, "3") for aa in "GNPSD"]
)

# '1'stand for Buried; '2'stand for Exposed, '3' stand for Intermediate
_SolventAccessibility = dict(
    [(aa, "1") for aa in "ALFCGIVW"]
    + [(aa, "2") for aa in "RKQEND"]
    + [(aa, "3") for aa in "MPSTHY"]
)

# '1'stand for (0-0.108); '2'stand for (0.128-0.186), '3' stand for (0.219-0.409)
_Polarizability = dict(
    [(aa, "1") for aa in "GASDT"]
    + [(aa, "2") for aa in "CPNVEQIL"]
    + [(aa, "3") for aa in "KMHFRYW"]
)

_prop_groups = [
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

    def _seq2number_strings(self, seq, properties):
        """Auxiliary function to convert between aa and properties groups (PRIVATE)."""
        prop_strings = []
        for prop in properties:
            prop_strings.append("".join([prop[aa] for aa in seq]))
        return prop_strings

    def _calc_C(self, prop_strings):
        """Calculate the Composition descriptors of CTD (PRIVATE)."""
        C_descriptors = []
        for s in prop_strings:
            C_descriptors += [s.count(i) for i in "123"]

        return C_descriptors

    def _calc_T(self, prop_strings):
        """Calculate the Transition descriptors of CTD (PRIVATE)."""
        T_descriptors = []
        for s in prop_strings:
            for t in ("12", "13", "23"):
                T_descriptors.append((s.count(t) + s.count(t[::-1])) / (self.L - 1))

        return T_descriptors

    def _calc_D(self, prop_strings):
        """Calculate the Distribution descriptors of CTD (PRIVATE)."""
        # TODO: more descriptive var names
        # TODO: probably a faster way of calculating this
        D_descriptors = []
        for s in prop_strings:
            positions = {"1": [], "2": [], "3": []}
            # Get the positions of each aa
            for p, c in enumerate(s):
                positions[c].append(p)

            # calculate the descriptor values
            for key in positions:
                pos_list = positions[key]
                if pos_list == []:
                    D_descriptors += [0, 0, 0, 0, 0]
                    continue
                l = len(pos_list) / 5
                for i in range(5):
                    D_descriptors.append(pos_list[int(l * i)] / self.L)

        return D_descriptors

    def composition_descriptors(self, properties=_prop_groups):
        """Return the Composition descriptors of CTD."""
        prop_strings = self._seq2number_strings(self.seq, properties)

        return self._calc_T(prop_strings)

    def transition_descriptors(self, properties=_prop_groups):
        """Return the Transition descriptors of CTD."""
        prop_strings = self._seq2number_strings(self.seq, properties)

        return self._calc_T(prop_strings)

    def distribution_descriptors(self, properties=_prop_groups):
        """Return the Distribution descriptors of CTD."""
        prop_strings = self._seq2number_strings(self.seq, properties)

        return self._calc_D(prop_strings)

    def CTD_descriptors(self, properties=_prop_groups):
        """Return the CTD descriptors."""
        # Taking N as = len(properties)
        # First N*3 values are the composition of each group type
        # The next N*3 values are the Transition properties
        # The final N*5*3 values are the Distribution properties
        prop_strings = self._seq2number_strings(self.seq, properties)

        return (
            self._calc_C(prop_strings)
            + self._calc_T(prop_strings)
            + self._calc_D(prop_strings)
        )
