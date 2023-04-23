# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Classes and methods for tree construction."""

import itertools
import copy
import numbers
from Bio.Phylo import BaseTree
from Bio.Align import Alignment, MultipleSeqAlignment
from Bio.Align import substitution_matrices


# flake8: noqa


class _Matrix:
    """Base class for distance matrix or scoring matrix.

    Accepts a list of names and a lower triangular matrix.::

        matrix = [[0],
                  [1, 0],
                  [2, 3, 0],
                  [4, 5, 6, 0]]
        represents the symmetric matrix of
        [0,1,2,4]
        [1,0,3,5]
        [2,3,0,6]
        [4,5,6,0]

    :Parameters:
        names : list
            names of elements, used for indexing
        matrix : list
            nested list of numerical lists in lower triangular format

    Examples
    --------
    >>> from Bio.Phylo.TreeConstruction import _Matrix
    >>> names = ['Alpha', 'Beta', 'Gamma', 'Delta']
    >>> matrix = [[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]]
    >>> m = _Matrix(names, matrix)
    >>> m
    _Matrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], matrix=[[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]])

    You can use two indices to get or assign an element in the matrix.

    >>> m[1,2]
    3
    >>> m['Beta','Gamma']
    3
    >>> m['Beta','Gamma'] = 4
    >>> m['Beta','Gamma']
    4

    Further more, you can use one index to get or assign a list of elements related to that index.

    >>> m[0]
    [0, 1, 2, 4]
    >>> m['Alpha']
    [0, 1, 2, 4]
    >>> m['Alpha'] = [0, 7, 8, 9]
    >>> m[0]
    [0, 7, 8, 9]
    >>> m[0,1]
    7

    Also you can delete or insert a column&row of elements by index.

    >>> m
    _Matrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], matrix=[[0], [7, 0], [8, 4, 0], [9, 5, 6, 0]])
    >>> del m['Alpha']
    >>> m
    _Matrix(names=['Beta', 'Gamma', 'Delta'], matrix=[[0], [4, 0], [5, 6, 0]])
    >>> m.insert('Alpha', [0, 7, 8, 9] , 0)
    >>> m
    _Matrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], matrix=[[0], [7, 0], [8, 4, 0], [9, 5, 6, 0]])

    """

    def __init__(self, names, matrix=None):
        """Initialize matrix.

        Arguments are a list of names, and optionally a list of lower
        triangular matrix data (zero matrix used by default).
        """
        # check names
        if isinstance(names, list) and all(isinstance(s, str) for s in names):
            if len(set(names)) == len(names):
                self.names = names
            else:
                raise ValueError("Duplicate names found")
        else:
            raise TypeError("'names' should be a list of strings")

        # check matrix
        if matrix is None:
            # create a new one with 0 if matrix is not assigned
            matrix = [[0] * i for i in range(1, len(self) + 1)]
            self.matrix = matrix
        else:
            # check if all elements are numbers
            if (
                isinstance(matrix, list)
                and all(isinstance(row, list) for row in matrix)
                and all(
                    isinstance(item, numbers.Number) for row in matrix for item in row
                )
            ):
                # check if the same length with names
                if len(matrix) == len(names):
                    # check if is lower triangle format
                    if [len(row) for row in matrix] == list(range(1, len(self) + 1)):
                        self.matrix = matrix
                    else:
                        raise ValueError("'matrix' should be in lower triangle format")
                else:
                    raise ValueError("'names' and 'matrix' should be the same size")
            else:
                raise TypeError("'matrix' should be a list of numerical lists")

    def __getitem__(self, item):
        """Access value(s) by the index(s) or name(s).

        For a _Matrix object 'dm'::

            dm[i]                   get a value list from the given 'i' to others;
            dm[i, j]                get the value between 'i' and 'j';
            dm['name']              map name to index first
            dm['name1', 'name2']    map name to index first

        """
        # Handle single indexing
        if isinstance(item, (int, str)):
            index = None
            if isinstance(item, int):
                index = item
            elif isinstance(item, str):
                if item in self.names:
                    index = self.names.index(item)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if index > len(self) - 1:
                raise IndexError("Index out of range.")
            return [self.matrix[index][i] for i in range(0, index)] + [
                self.matrix[i][index] for i in range(index, len(self))
            ]
        # Handle double indexing
        elif len(item) == 2:
            row_index = None
            col_index = None
            if all(isinstance(i, int) for i in item):
                row_index, col_index = item
            elif all(isinstance(i, str) for i in item):
                row_name, col_name = item
                if row_name in self.names and col_name in self.names:
                    row_index = self.names.index(row_name)
                    col_index = self.names.index(col_name)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if row_index > len(self) - 1 or col_index > len(self) - 1:
                raise IndexError("Index out of range.")
            if row_index > col_index:
                return self.matrix[row_index][col_index]
            else:
                return self.matrix[col_index][row_index]
        else:
            raise TypeError("Invalid index type.")

    def __setitem__(self, item, value):
        """Set value by the index(s) or name(s).

        Similar to __getitem__::

            dm[1] = [1, 0, 3, 4]    set values from '1' to others;
            dm[i, j] = 2            set the value from 'i' to 'j'

        """
        # Handle single indexing
        if isinstance(item, (int, str)):
            index = None
            if isinstance(item, int):
                index = item
            elif isinstance(item, str):
                if item in self.names:
                    index = self.names.index(item)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if index > len(self) - 1:
                raise IndexError("Index out of range.")
            # check and assign value
            if isinstance(value, list) and all(
                isinstance(n, numbers.Number) for n in value
            ):
                if len(value) == len(self):
                    for i in range(0, index):
                        self.matrix[index][i] = value[i]
                    for i in range(index, len(self)):
                        self.matrix[i][index] = value[i]
                else:
                    raise ValueError("Value not the same size.")
            else:
                raise TypeError("Invalid value type.")
        # Handle double indexing
        elif len(item) == 2:
            row_index = None
            col_index = None
            if all(isinstance(i, int) for i in item):
                row_index, col_index = item
            elif all(isinstance(i, str) for i in item):
                row_name, col_name = item
                if row_name in self.names and col_name in self.names:
                    row_index = self.names.index(row_name)
                    col_index = self.names.index(col_name)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if row_index > len(self) - 1 or col_index > len(self) - 1:
                raise IndexError("Index out of range.")
            # check and assign value
            if isinstance(value, numbers.Number):
                if row_index > col_index:
                    self.matrix[row_index][col_index] = value
                else:
                    self.matrix[col_index][row_index] = value
            else:
                raise TypeError("Invalid value type.")
        else:
            raise TypeError("Invalid index type.")

    def __delitem__(self, item):
        """Delete related distances by the index or name."""
        index = None
        if isinstance(item, int):
            index = item
        elif isinstance(item, str):
            index = self.names.index(item)
        else:
            raise TypeError("Invalid index type.")
        # remove distances related to index
        for i in range(index + 1, len(self)):
            del self.matrix[i][index]
        del self.matrix[index]
        # remove name
        del self.names[index]

    def insert(self, name, value, index=None):
        """Insert distances given the name and value.

        :Parameters:
            name : str
                name of a row/col to be inserted
            value : list
                a row/col of values to be inserted

        """
        if isinstance(name, str):
            # insert at the given index or at the end
            if index is None:
                index = len(self)
            if not isinstance(index, int):
                raise TypeError("Invalid index type.")
            # insert name
            self.names.insert(index, name)
            # insert elements of 0, to be assigned
            self.matrix.insert(index, [0] * index)
            for i in range(index, len(self)):
                self.matrix[i].insert(index, 0)
            # assign value
            self[index] = value
        else:
            raise TypeError("Invalid name type.")

    def __len__(self):
        """Matrix length."""
        return len(self.names)

    def __repr__(self):
        """Return Matrix as a string."""
        return self.__class__.__name__ + "(names=%s, matrix=%s)" % tuple(
            map(repr, (self.names, self.matrix))
        )

    def __str__(self):
        """Get a lower triangular matrix string."""
        matrix_string = "\n".join(
            [
                self.names[i]
                + "\t"
                + "\t".join([format(n, "f") for n in self.matrix[i]])
                for i in range(0, len(self))
            ]
        )
        matrix_string = matrix_string + "\n\t" + "\t".join(self.names)
        return matrix_string.expandtabs(tabsize=4)


class DistanceMatrix(_Matrix):
    """Distance matrix class that can be used for distance based tree algorithms.

    All diagonal elements will be zero no matter what the users provide.
    """

    def __init__(self, names, matrix=None):
        """Initialize the class."""
        _Matrix.__init__(self, names, matrix)
        self._set_zero_diagonal()

    def __setitem__(self, item, value):
        """Set Matrix's items to values."""
        _Matrix.__setitem__(self, item, value)
        self._set_zero_diagonal()

    def _set_zero_diagonal(self):
        """Set all diagonal elements to zero (PRIVATE)."""
        for i in range(0, len(self)):
            self.matrix[i][i] = 0

    def format_phylip(self, handle):
        """Write data in Phylip format to a given file-like object or handle.

        The output stream is the input distance matrix format used with Phylip
        programs (e.g. 'neighbor'). See:
        http://evolution.genetics.washington.edu/phylip/doc/neighbor.html

        :Parameters:
            handle : file or file-like object
                A writeable text mode file handle or other object supporting
                the 'write' method, such as StringIO or sys.stdout.

        """
        handle.write(f"    {len(self.names)}\n")
        # Phylip needs space-separated, vertically aligned columns
        name_width = max(12, max(map(len, self.names)) + 1)
        value_fmts = ("{" + str(x) + ":.4f}" for x in range(1, len(self.matrix) + 1))
        row_fmt = "{0:" + str(name_width) + "s}" + "  ".join(value_fmts) + "\n"
        for i, (name, values) in enumerate(zip(self.names, self.matrix)):
            # Mirror the matrix values across the diagonal
            mirror_values = (self.matrix[j][i] for j in range(i + 1, len(self.matrix)))
            fields = itertools.chain([name], values, mirror_values)
            handle.write(row_fmt.format(*fields))


# Shim for compatibility with Biopython<1.70 (#1304)
_DistanceMatrix = DistanceMatrix


class DistanceCalculator:
    """Calculates the distance matrix from a DNA or protein sequence alignment.

    This class calculates the distance matrix from a multiple sequence alignment
    of DNA or protein sequences, and the given name of the substitution model.

    Currently only scoring matrices are used.

    :Parameters:
        model : str
            Name of the model matrix to be used to calculate distance.
            The attribute ``dna_models`` contains the available model
            names for DNA sequences and ``protein_models`` for protein
            sequences.

    Examples
    --------
    Loading a small PHYLIP alignment from which to compute distances::

      >>> from Bio.Phylo.TreeConstruction import DistanceCalculator
      >>> from Bio import AlignIO
      >>> aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')
      >>> print(aln)  # doctest:+NORMALIZE_WHITESPACE
      Alignment with 5 rows and 13 columns
      AACGTGGCCACAT Alpha
      AAGGTCGCCACAC Beta
      CAGTTCGCCACAA Gamma
      GAGATTTCCGCCT Delta
      GAGATCTCCGCCC Epsilon

    DNA calculator with 'identity' model::

      >>> calculator = DistanceCalculator('identity')
      >>> dm = calculator.get_distance(aln)
      >>> print(dm)  # doctest:+NORMALIZE_WHITESPACE
        Alpha   0.000000
        Beta    0.230769    0.000000
        Gamma   0.384615    0.230769    0.000000
        Delta   0.538462    0.538462    0.538462    0.000000
        Epsilon 0.615385    0.384615    0.461538    0.153846    0.000000
            Alpha   Beta    Gamma   Delta   Epsilon

    Protein calculator with 'blosum62' model::

      >>> calculator = DistanceCalculator('blosum62')
      >>> dm = calculator.get_distance(aln)
      >>> print(dm)  # doctest:+NORMALIZE_WHITESPACE
      Alpha   0.000000
      Beta    0.369048    0.000000
      Gamma   0.493976    0.250000    0.000000
      Delta   0.585366    0.547619    0.566265    0.000000
      Epsilon 0.700000    0.355556    0.488889    0.222222    0.000000
          Alpha   Beta    Gamma   Delta   Epsilon

    Same calculation, using the new Alignment object::

      >>> from Bio.Phylo.TreeConstruction import DistanceCalculator
      >>> from Bio import Align
      >>> aln = Align.read('TreeConstruction/msa.phy', 'phylip')
      >>> print(aln)  # doctest:+NORMALIZE_WHITESPACE
      Alpha             0 AACGTGGCCACAT 13
      Beta              0 AAGGTCGCCACAC 13
      Gamma             0 CAGTTCGCCACAA 13
      Delta             0 GAGATTTCCGCCT 13
      Epsilon           0 GAGATCTCCGCCC 13
      <BLANKLINE>

    DNA calculator with 'identity' model::

      >>> calculator = DistanceCalculator('identity')
      >>> dm = calculator.get_distance(aln)
      >>> print(dm)  # doctest:+NORMALIZE_WHITESPACE
      Alpha   0.000000
      Beta    0.230769    0.000000
      Gamma   0.384615    0.230769    0.000000
      Delta   0.538462    0.538462    0.538462    0.000000
      Epsilon 0.615385    0.384615    0.461538    0.153846    0.000000
          Alpha   Beta    Gamma   Delta   Epsilon

    Protein calculator with 'blosum62' model::

      >>> calculator = DistanceCalculator('blosum62')
      >>> dm = calculator.get_distance(aln)
      >>> print(dm)  # doctest:+NORMALIZE_WHITESPACE
      Alpha   0.000000
      Beta    0.369048    0.000000
      Gamma   0.493976    0.250000    0.000000
      Delta   0.585366    0.547619    0.566265    0.000000
      Epsilon 0.700000    0.355556    0.488889    0.222222    0.000000
          Alpha   Beta    Gamma   Delta   Epsilon

    """

    protein_alphabet = set("ABCDEFGHIKLMNPQRSTVWXYZ")

    dna_models = []
    protein_models = []

    # matrices available
    names = substitution_matrices.load()
    for name in names:
        matrix = substitution_matrices.load(name)
        if name == "NUC.4.4":
            # BLAST nucleic acid scoring matrix
            name = "blastn"
        else:
            name = name.lower()
        if protein_alphabet.issubset(set(matrix.alphabet)):
            protein_models.append(name)
        else:
            dna_models.append(name)

    del protein_alphabet
    del name
    del names
    del matrix

    models = ["identity"] + dna_models + protein_models

    def __init__(self, model="identity", skip_letters=None):
        """Initialize with a distance model."""
        # Shim for backward compatibility (#491)
        if skip_letters:
            self.skip_letters = skip_letters
        elif model == "identity":
            self.skip_letters = ()
        else:
            self.skip_letters = ("-", "*")

        if model == "identity":
            self.scoring_matrix = None
        elif model in self.models:
            if model == "blastn":
                name = "NUC.4.4"
            else:
                name = model.upper()
            self.scoring_matrix = substitution_matrices.load(name)
        else:
            raise ValueError(
                "Model not supported. Available models: " + ", ".join(self.models)
            )

    def _pairwise(self, seq1, seq2):
        """Calculate pairwise distance from two sequences (PRIVATE).

        Returns a value between 0 (identical sequences) and 1 (completely
        different, or seq1 is an empty string.)
        """
        score = 0
        max_score = 0
        if self.scoring_matrix is None:
            # Score by character identity, not skipping any special letters
            score = sum(
                l1 == l2
                for l1, l2 in zip(seq1, seq2)
                if l1 not in self.skip_letters and l2 not in self.skip_letters
            )
            max_score = len(seq1)
        else:
            max_score1 = 0
            max_score2 = 0
            for i in range(0, len(seq1)):
                l1 = seq1[i]
                l2 = seq2[i]
                if l1 in self.skip_letters or l2 in self.skip_letters:
                    continue
                try:
                    max_score1 += self.scoring_matrix[l1, l1]
                except IndexError:
                    raise ValueError(
                        f"Bad letter '{l1}' in sequence '{seq1.id}' at position '{i}'"
                    ) from None
                try:
                    max_score2 += self.scoring_matrix[l2, l2]
                except IndexError:
                    raise ValueError(
                        f"Bad letter '{l2}' in sequence '{seq2.id}' at position '{i}'"
                    ) from None
                score += self.scoring_matrix[l1, l2]
            # Take the higher score if the matrix is asymmetrical
            max_score = max(max_score1, max_score2)
        if max_score == 0:
            return 1  # max possible scaled distance
        return 1 - (score / max_score)

    def get_distance(self, msa):
        """Return a DistanceMatrix for an Alignment or MultipleSeqAlignment object.

        :Parameters:
            msa : Alignment or MultipleSeqAlignment object representing a
                DNA or protein multiple sequence alignment.

        """
        if isinstance(msa, Alignment):
            names = [s.id for s in msa.sequences]
            dm = DistanceMatrix(names)
            n = len(names)
            for i1 in range(n):
                for i2 in range(i1):
                    dm[names[i1], names[i2]] = self._pairwise(msa[i1], msa[i2])
        elif isinstance(msa, MultipleSeqAlignment):
            names = [s.id for s in msa]
            dm = DistanceMatrix(names)
            for seq1, seq2 in itertools.combinations(msa, 2):
                dm[seq1.id, seq2.id] = self._pairwise(seq1, seq2)
        else:
            raise TypeError(
                "Must provide an Alignment object or a MultipleSeqAlignment object."
            )

        return dm


class TreeConstructor:
    """Base class for all tree constructor."""

    def build_tree(self, msa):
        """Caller to build the tree from an Alignment or MultipleSeqAlignment object.

        This should be implemented in subclass.
        """
        raise NotImplementedError("Method not implemented!")


class DistanceTreeConstructor(TreeConstructor):
    """Distance based tree constructor.

    :Parameters:
        method : str
            Distance tree construction method, 'nj'(default) or 'upgma'.
        distance_calculator : DistanceCalculator
            The distance matrix calculator for multiple sequence alignment.
            It must be provided if ``build_tree`` will be called.

    Examples
    --------
    Loading a small PHYLIP alignment from which to compute distances, and then
    build a upgma Tree::

      >>> from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
      >>> from Bio.Phylo.TreeConstruction import DistanceCalculator
      >>> from Bio import AlignIO
      >>> aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')
      >>> constructor = DistanceTreeConstructor()
      >>> calculator = DistanceCalculator('identity')
      >>> dm = calculator.get_distance(aln)
      >>> upgmatree = constructor.upgma(dm)
      >>> print(upgmatree)
      Tree(rooted=True)
          Clade(branch_length=0, name='Inner4')
              Clade(branch_length=0.18749999999999994, name='Inner1')
                  Clade(branch_length=0.07692307692307693, name='Epsilon')
                  Clade(branch_length=0.07692307692307693, name='Delta')
              Clade(branch_length=0.11057692307692304, name='Inner3')
                  Clade(branch_length=0.038461538461538464, name='Inner2')
                      Clade(branch_length=0.11538461538461536, name='Gamma')
                      Clade(branch_length=0.11538461538461536, name='Beta')
                  Clade(branch_length=0.15384615384615383, name='Alpha')

    Build a NJ Tree::

      >>> njtree = constructor.nj(dm)
      >>> print(njtree)
      Tree(rooted=False)
          Clade(branch_length=0, name='Inner3')
              Clade(branch_length=0.18269230769230765, name='Alpha')
              Clade(branch_length=0.04807692307692307, name='Beta')
              Clade(branch_length=0.04807692307692307, name='Inner2')
                  Clade(branch_length=0.27884615384615385, name='Inner1')
                      Clade(branch_length=0.051282051282051266, name='Epsilon')
                      Clade(branch_length=0.10256410256410259, name='Delta')
                  Clade(branch_length=0.14423076923076922, name='Gamma')

    Same example, using the new Alignment class::

      >>> from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
      >>> from Bio.Phylo.TreeConstruction import DistanceCalculator
      >>> from Bio import Align
      >>> aln = Align.read(open('TreeConstruction/msa.phy'), 'phylip')
      >>> constructor = DistanceTreeConstructor()
      >>> calculator = DistanceCalculator('identity')
      >>> dm = calculator.get_distance(aln)
      >>> upgmatree = constructor.upgma(dm)
      >>> print(upgmatree)
      Tree(rooted=True)
          Clade(branch_length=0, name='Inner4')
              Clade(branch_length=0.18749999999999994, name='Inner1')
                  Clade(branch_length=0.07692307692307693, name='Epsilon')
                  Clade(branch_length=0.07692307692307693, name='Delta')
              Clade(branch_length=0.11057692307692304, name='Inner3')
                  Clade(branch_length=0.038461538461538464, name='Inner2')
                      Clade(branch_length=0.11538461538461536, name='Gamma')
                      Clade(branch_length=0.11538461538461536, name='Beta')
                  Clade(branch_length=0.15384615384615383, name='Alpha')

    Build a NJ Tree::

      >>> njtree = constructor.nj(dm)
      >>> print(njtree)
      Tree(rooted=False)
          Clade(branch_length=0, name='Inner3')
              Clade(branch_length=0.18269230769230765, name='Alpha')
              Clade(branch_length=0.04807692307692307, name='Beta')
              Clade(branch_length=0.04807692307692307, name='Inner2')
                  Clade(branch_length=0.27884615384615385, name='Inner1')
                      Clade(branch_length=0.051282051282051266, name='Epsilon')
                      Clade(branch_length=0.10256410256410259, name='Delta')
                  Clade(branch_length=0.14423076923076922, name='Gamma')

    """

    methods = ["nj", "upgma"]

    def __init__(self, distance_calculator=None, method="nj"):
        """Initialize the class."""
        if distance_calculator is None or isinstance(
            distance_calculator, DistanceCalculator
        ):
            self.distance_calculator = distance_calculator
        else:
            raise TypeError("Must provide a DistanceCalculator object.")
        if method in self.methods:
            self.method = method
        else:
            raise TypeError(
                "Bad method: "
                + method
                + ". Available methods: "
                + ", ".join(self.methods)
            )

    def build_tree(self, msa):
        """Construct and return a Tree, Neighbor Joining or UPGMA."""
        if self.distance_calculator:
            dm = self.distance_calculator.get_distance(msa)
            tree = None
            if self.method == "upgma":
                tree = self.upgma(dm)
            else:
                tree = self.nj(dm)
            return tree
        else:
            raise TypeError("Must provide a DistanceCalculator object.")

    def upgma(self, distance_matrix):
        """Construct and return an UPGMA tree.

        Constructs and returns an Unweighted Pair Group Method
        with Arithmetic mean (UPGMA) tree.

        :Parameters:
            distance_matrix : DistanceMatrix
                The distance matrix for tree construction.

        """
        if not isinstance(distance_matrix, DistanceMatrix):
            raise TypeError("Must provide a DistanceMatrix object.")

        # make a copy of the distance matrix to be used
        dm = copy.deepcopy(distance_matrix)
        # init terminal clades
        clades = [BaseTree.Clade(None, name) for name in dm.names]
        # init minimum index
        min_i = 0
        min_j = 0
        inner_count = 0
        while len(dm) > 1:
            min_dist = dm[1, 0]
            # find minimum index
            for i in range(1, len(dm)):
                for j in range(0, i):
                    if min_dist >= dm[i, j]:
                        min_dist = dm[i, j]
                        min_i = i
                        min_j = j

            # create clade
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            inner_count += 1
            inner_clade = BaseTree.Clade(None, "Inner" + str(inner_count))
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            # assign branch length
            if clade1.is_terminal():
                clade1.branch_length = min_dist / 2
            else:
                clade1.branch_length = min_dist / 2 - self._height_of(clade1)

            if clade2.is_terminal():
                clade2.branch_length = min_dist / 2
            else:
                clade2.branch_length = min_dist / 2 - self._height_of(clade2)

            # update node list
            clades[min_j] = inner_clade
            del clades[min_i]

            # rebuild distance matrix,
            # set the distances of new node at the index of min_j
            for k in range(0, len(dm)):
                if k != min_i and k != min_j:
                    dm[min_j, k] = (dm[min_i, k] + dm[min_j, k]) / 2

            dm.names[min_j] = "Inner" + str(inner_count)

            del dm[min_i]
        inner_clade.branch_length = 0
        return BaseTree.Tree(inner_clade)

    def nj(self, distance_matrix):
        """Construct and return a Neighbor Joining tree.

        :Parameters:
            distance_matrix : DistanceMatrix
                The distance matrix for tree construction.

        """
        if not isinstance(distance_matrix, DistanceMatrix):
            raise TypeError("Must provide a DistanceMatrix object.")

        # make a copy of the distance matrix to be used
        dm = copy.deepcopy(distance_matrix)
        # init terminal clades
        clades = [BaseTree.Clade(None, name) for name in dm.names]
        # init node distance
        node_dist = [0] * len(dm)
        # init minimum index
        min_i = 0
        min_j = 0
        inner_count = 0
        # special cases for Minimum Alignment Matrices
        if len(dm) == 1:
            root = clades[0]

            return BaseTree.Tree(root, rooted=False)
        elif len(dm) == 2:
            # minimum distance will always be [1,0]
            min_i = 1
            min_j = 0
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            clade1.branch_length = dm[min_i, min_j] / 2.0
            clade2.branch_length = dm[min_i, min_j] - clade1.branch_length
            inner_clade = BaseTree.Clade(None, "Inner")
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            clades[0] = inner_clade
            root = clades[0]

            return BaseTree.Tree(root, rooted=False)
        while len(dm) > 2:
            # calculate nodeDist
            for i in range(0, len(dm)):
                node_dist[i] = 0
                for j in range(0, len(dm)):
                    node_dist[i] += dm[i, j]
                node_dist[i] = node_dist[i] / (len(dm) - 2)

            # find minimum distance pair
            min_dist = dm[1, 0] - node_dist[1] - node_dist[0]
            min_i = 0
            min_j = 1
            for i in range(1, len(dm)):
                for j in range(0, i):
                    temp = dm[i, j] - node_dist[i] - node_dist[j]
                    if min_dist > temp:
                        min_dist = temp
                        min_i = i
                        min_j = j
            # create clade
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            inner_count += 1
            inner_clade = BaseTree.Clade(None, "Inner" + str(inner_count))
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            # assign branch length
            clade1.branch_length = (
                dm[min_i, min_j] + node_dist[min_i] - node_dist[min_j]
            ) / 2.0
            clade2.branch_length = dm[min_i, min_j] - clade1.branch_length

            # update node list
            clades[min_j] = inner_clade
            del clades[min_i]

            # rebuild distance matrix,
            # set the distances of new node at the index of min_j
            for k in range(0, len(dm)):
                if k != min_i and k != min_j:
                    dm[min_j, k] = (
                        dm[min_i, k] + dm[min_j, k] - dm[min_i, min_j]
                    ) / 2.0

            dm.names[min_j] = "Inner" + str(inner_count)
            del dm[min_i]

        # set the last clade as one of the child of the inner_clade
        root = None
        if clades[0] == inner_clade:
            clades[0].branch_length = 0
            clades[1].branch_length = dm[1, 0]
            clades[0].clades.append(clades[1])
            root = clades[0]
        else:
            clades[0].branch_length = dm[1, 0]
            clades[1].branch_length = 0
            clades[1].clades.append(clades[0])
            root = clades[1]

        return BaseTree.Tree(root, rooted=False)

    def _height_of(self, clade):
        """Calculate clade height -- the longest path to any terminal (PRIVATE)."""
        height = 0
        if clade.is_terminal():
            height = clade.branch_length
        else:
            height = height + max(self._height_of(c) for c in clade.clades)
        return height


# #################### Tree Scoring and Searching Classes #####################


class Scorer:
    """Base class for all tree scoring methods."""

    def get_score(self, tree, alignment):
        """Caller to get the score of a tree for the given alignment.

        This should be implemented in subclass.
        """
        raise NotImplementedError("Method not implemented!")


class TreeSearcher:
    """Base class for all tree searching methods."""

    def search(self, starting_tree, alignment):
        """Caller to search the best tree with a starting tree.

        This should be implemented in subclass.
        """
        raise NotImplementedError("Method not implemented!")


class NNITreeSearcher(TreeSearcher):
    """Tree searching with Nearest Neighbor Interchanges (NNI) algorithm.

    :Parameters:
        scorer : ParsimonyScorer
            parsimony scorer to calculate the parsimony score of
            different trees during NNI algorithm.

    """

    def __init__(self, scorer):
        """Initialize the class."""
        if isinstance(scorer, Scorer):
            self.scorer = scorer
        else:
            raise TypeError("Must provide a Scorer object.")

    def search(self, starting_tree, alignment):
        """Implement the TreeSearcher.search method.

        :Parameters:
           starting_tree : Tree
               starting tree of NNI method.
           alignment : Alignment or MultipleSeqAlignment object
               multiple sequence alignment used to calculate parsimony
               score of different NNI trees.

        """
        return self._nni(starting_tree, alignment)

    def _nni(self, starting_tree, alignment):
        """Search for the best parsimony tree using the NNI algorithm (PRIVATE)."""
        best_tree = starting_tree
        while True:
            best_score = self.scorer.get_score(best_tree, alignment)
            temp = best_score
            for t in self._get_neighbors(best_tree):
                score = self.scorer.get_score(t, alignment)
                if score < best_score:
                    best_score = score
                    best_tree = t
            # stop if no smaller score exist
            if best_score >= temp:
                break
        return best_tree

    def _get_neighbors(self, tree):
        """Get all neighbor trees of the given tree (PRIVATE).

        Currently only for binary rooted trees.
        """
        # make child to parent dict
        parents = {}
        for clade in tree.find_clades():
            if clade != tree.root:
                node_path = tree.get_path(clade)
                # cannot get the parent if the parent is root. Bug?
                if len(node_path) == 1:
                    parents[clade] = tree.root
                else:
                    parents[clade] = node_path[-2]
        neighbors = []
        root_childs = []
        for clade in tree.get_nonterminals(order="level"):
            if clade == tree.root:
                left = clade.clades[0]
                right = clade.clades[1]
                root_childs.append(left)
                root_childs.append(right)
                if not left.is_terminal() and not right.is_terminal():
                    # make changes around the left_left clade
                    # left_left = left.clades[0]
                    left_right = left.clades[1]
                    right_left = right.clades[0]
                    right_right = right.clades[1]
                    # neighbor 1 (left_left + right_right)
                    del left.clades[1]
                    del right.clades[1]
                    left.clades.append(right_right)
                    right.clades.append(left_right)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # neighbor 2 (left_left + right_left)
                    del left.clades[1]
                    del right.clades[0]
                    left.clades.append(right_left)
                    right.clades.append(right_right)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # change back (left_left + left_right)
                    del left.clades[1]
                    del right.clades[0]
                    left.clades.append(left_right)
                    right.clades.insert(0, right_left)
            elif clade in root_childs:
                # skip root child
                continue
            else:
                # method for other clades
                # make changes around the parent clade
                left = clade.clades[0]
                right = clade.clades[1]
                parent = parents[clade]
                if clade == parent.clades[0]:
                    sister = parent.clades[1]
                    # neighbor 1 (parent + right)
                    del parent.clades[1]
                    del clade.clades[1]
                    parent.clades.append(right)
                    clade.clades.append(sister)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # neighbor 2 (parent + left)
                    del parent.clades[1]
                    del clade.clades[0]
                    parent.clades.append(left)
                    clade.clades.append(right)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # change back (parent + sister)
                    del parent.clades[1]
                    del clade.clades[0]
                    parent.clades.append(sister)
                    clade.clades.insert(0, left)
                else:
                    sister = parent.clades[0]
                    # neighbor 1 (parent + right)
                    del parent.clades[0]
                    del clade.clades[1]
                    parent.clades.insert(0, right)
                    clade.clades.append(sister)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # neighbor 2 (parent + left)
                    del parent.clades[0]
                    del clade.clades[0]
                    parent.clades.insert(0, left)
                    clade.clades.append(right)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # change back (parent + sister)
                    del parent.clades[0]
                    del clade.clades[0]
                    parent.clades.insert(0, sister)
                    clade.clades.insert(0, left)
        return neighbors


# ######################## Parsimony Classes ##########################


class ParsimonyScorer(Scorer):
    """Parsimony scorer with a scoring matrix.

    This is a combination of Fitch algorithm and Sankoff algorithm.
    See ParsimonyTreeConstructor for usage.

    :Parameters:
        matrix : _Matrix
            scoring matrix used in parsimony score calculation.

    """

    def __init__(self, matrix=None):
        """Initialize the class."""
        if not matrix or isinstance(matrix, _Matrix):
            self.matrix = matrix
        else:
            raise TypeError("Must provide a _Matrix object.")

    def get_score(self, tree, alignment):
        """Calculate parsimony score using the Fitch algorithm.

        Calculate and return the parsimony score given a tree and the
        MSA using either the Fitch algorithm (without a penalty matrix)
        or the Sankoff algorithm (with a matrix).
        """
        # make sure the tree is rooted and bifurcating
        if not tree.is_bifurcating():
            raise ValueError("The tree provided should be bifurcating.")
        if not tree.rooted:
            tree.root_at_midpoint()
        # sort tree terminals and alignment
        terms = tree.get_terminals()
        terms.sort(key=lambda term: term.name)
        alignment.sort()
        if isinstance(alignment, MultipleSeqAlignment):
            if not all(t.name == a.id for t, a in zip(terms, alignment)):
                raise ValueError(
                    "Taxon names of the input tree should be the same with the alignment."
                )
        else:  # Alignment object
            if not all(t.name == s.id for t, s in zip(terms, alignment.sequences)):
                raise ValueError(
                    "Taxon names of the input tree should be the same with the alignment."
                )
        # term_align = dict(zip(terms, alignment))
        score = 0
        for i in range(len(alignment[0])):
            # parsimony score for column_i
            score_i = 0
            # get column
            column_i = alignment[:, i]
            # skip non-informative column
            if column_i == len(column_i) * column_i[0]:
                continue

            # start calculating score_i using the tree and column_i

            # Fitch algorithm without the penalty matrix
            if not self.matrix:
                # init by mapping terminal clades and states in column_i
                clade_states = dict(zip(terms, [{c} for c in column_i]))
                for clade in tree.get_nonterminals(order="postorder"):
                    clade_childs = clade.clades
                    left_state = clade_states[clade_childs[0]]
                    right_state = clade_states[clade_childs[1]]
                    state = left_state & right_state
                    if not state:
                        state = left_state | right_state
                        score_i += 1
                    clade_states[clade] = state
            # Sankoff algorithm with the penalty matrix
            else:
                inf = float("inf")
                # init score arrays for terminal clades
                alphabet = self.matrix.names
                length = len(alphabet)
                clade_scores = {}
                for j in range(len(column_i)):
                    array = [inf] * length
                    index = alphabet.index(column_i[j])
                    array[index] = 0
                    clade_scores[terms[j]] = array
                # bottom up calculation
                for clade in tree.get_nonterminals(order="postorder"):
                    clade_childs = clade.clades
                    left_score = clade_scores[clade_childs[0]]
                    right_score = clade_scores[clade_childs[1]]
                    array = []
                    for m in range(length):
                        min_l = inf
                        min_r = inf
                        for n in range(length):
                            sl = self.matrix[alphabet[m], alphabet[n]] + left_score[n]
                            sr = self.matrix[alphabet[m], alphabet[n]] + right_score[n]
                            if min_l > sl:
                                min_l = sl
                            if min_r > sr:
                                min_r = sr
                        array.append(min_l + min_r)
                    clade_scores[clade] = array
                # minimum from root score
                score_i = min(array)
                # TODO: resolve internal states
            score += score_i
        return score


class ParsimonyTreeConstructor(TreeConstructor):
    """Parsimony tree constructor.

    :Parameters:
        searcher : TreeSearcher
            tree searcher to search the best parsimony tree.
        starting_tree : Tree
            starting tree provided to the searcher.

    Examples
    --------
    We will load an alignment, and then load various trees which have already been computed from it::

      >>> from Bio import AlignIO, Phylo
      >>> aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')
      >>> print(aln)
      Alignment with 5 rows and 13 columns
      AACGTGGCCACAT Alpha
      AAGGTCGCCACAC Beta
      CAGTTCGCCACAA Gamma
      GAGATTTCCGCCT Delta
      GAGATCTCCGCCC Epsilon

    Load a starting tree::

      >>> starting_tree = Phylo.read('TreeConstruction/nj.tre', 'newick')
      >>> print(starting_tree)
      Tree(rooted=False, weight=1.0)
          Clade(branch_length=0.0, name='Inner3')
              Clade(branch_length=0.01421, name='Inner2')
                  Clade(branch_length=0.23927, name='Inner1')
                      Clade(branch_length=0.08531, name='Epsilon')
                      Clade(branch_length=0.13691, name='Delta')
                  Clade(branch_length=0.2923, name='Alpha')
              Clade(branch_length=0.07477, name='Beta')
              Clade(branch_length=0.17523, name='Gamma')

    Build the Parsimony tree from the starting tree::

      >>> scorer = Phylo.TreeConstruction.ParsimonyScorer()
      >>> searcher = Phylo.TreeConstruction.NNITreeSearcher(scorer)
      >>> constructor = Phylo.TreeConstruction.ParsimonyTreeConstructor(searcher, starting_tree)
      >>> pars_tree = constructor.build_tree(aln)
      >>> print(pars_tree)
      Tree(rooted=True, weight=1.0)
          Clade(branch_length=0.0)
              Clade(branch_length=0.19732999999999998, name='Inner1')
                  Clade(branch_length=0.13691, name='Delta')
                  Clade(branch_length=0.08531, name='Epsilon')
              Clade(branch_length=0.04194000000000003, name='Inner2')
                  Clade(branch_length=0.01421, name='Inner3')
                      Clade(branch_length=0.17523, name='Gamma')
                      Clade(branch_length=0.07477, name='Beta')
                  Clade(branch_length=0.2923, name='Alpha')

    Same example, using the new Alignment class::

      >>> from Bio import Align, Phylo
      >>> alignment = Align.read(open('TreeConstruction/msa.phy'), 'phylip')
      >>> print(alignment)
      Alpha             0 AACGTGGCCACAT 13
      Beta              0 AAGGTCGCCACAC 13
      Gamma             0 CAGTTCGCCACAA 13
      Delta             0 GAGATTTCCGCCT 13
      Epsilon           0 GAGATCTCCGCCC 13
      <BLANKLINE>

    Load a starting tree::

      >>> starting_tree = Phylo.read('TreeConstruction/nj.tre', 'newick')
      >>> print(starting_tree)
      Tree(rooted=False, weight=1.0)
          Clade(branch_length=0.0, name='Inner3')
              Clade(branch_length=0.01421, name='Inner2')
                  Clade(branch_length=0.23927, name='Inner1')
                      Clade(branch_length=0.08531, name='Epsilon')
                      Clade(branch_length=0.13691, name='Delta')
                  Clade(branch_length=0.2923, name='Alpha')
              Clade(branch_length=0.07477, name='Beta')
              Clade(branch_length=0.17523, name='Gamma')

    Build the Parsimony tree from the starting tree::

      >>> scorer = Phylo.TreeConstruction.ParsimonyScorer()
      >>> searcher = Phylo.TreeConstruction.NNITreeSearcher(scorer)
      >>> constructor = Phylo.TreeConstruction.ParsimonyTreeConstructor(searcher, starting_tree)
      >>> pars_tree = constructor.build_tree(alignment)
      >>> print(pars_tree)
      Tree(rooted=True, weight=1.0)
          Clade(branch_length=0.0)
              Clade(branch_length=0.19732999999999998, name='Inner1')
                  Clade(branch_length=0.13691, name='Delta')
                  Clade(branch_length=0.08531, name='Epsilon')
              Clade(branch_length=0.04194000000000003, name='Inner2')
                  Clade(branch_length=0.01421, name='Inner3')
                      Clade(branch_length=0.17523, name='Gamma')
                      Clade(branch_length=0.07477, name='Beta')
                  Clade(branch_length=0.2923, name='Alpha')

    """

    def __init__(self, searcher, starting_tree=None):
        """Initialize the class."""
        self.searcher = searcher
        self.starting_tree = starting_tree

    def build_tree(self, alignment):
        """Build the tree.

        :Parameters:
            alignment : MultipleSeqAlignment
                multiple sequence alignment to calculate parsimony tree.

        """
        # if starting_tree is none,
        # create a upgma tree with 'identity' scoring matrix
        if self.starting_tree is None:
            dtc = DistanceTreeConstructor(DistanceCalculator("identity"), "upgma")
            self.starting_tree = dtc.build_tree(alignment)
        return self.searcher.search(self.starting_tree, alignment)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
