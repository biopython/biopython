# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Classes and methods for tree construction"""

import itertools
from Bio.Align import MultipleSeqAlignment

class Matrix(object):
    """A base class for distance matrix or scoring matrix that accepts
     a list of names and a lower triangular matrix.

    matrix = [[0],
              [1, 0],
              [2, 3, 0],
              [4, 5, 6, 0]]
    represents the symmetric matrix of
    [0,1,2,4]
    [1,0,3,5]
    [2,3,0,6]
    [4,5,6,0]
    """

    def __init__(self, names, matrix=None):
        """Initialize matrix by a list of names and a list of
        lower triangular matrix data"""
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
            matrix = [[0]*i for i in range(1, len(self) + 1)]
            self.matrix = matrix
        else:
            # check if all elements are numbers
            if isinstance(matrix, list) and all(isinstance(l, list) for l in matrix) and all(isinstance(n, (int, long, float, complex)) for n in [item for sublist in matrix for item in sublist]):
                # check if the same length with names
                if len(matrix) == len(names):
                    # check if is lower triangle format
                    if [len(m) for m in matrix] == range(1, len(self) + 1):
                        self.matrix = matrix
                    else:
                        raise ValueError("'matrix' should be in lower triangle format")
                else:
                    raise ValueError("'names' and 'matrix' should be the same size")
            else:
                raise TypeError("'matrix' should be a list of numerical lists")

    def __getitem__(self, item):
        """Access value(s) by the index(s) or name(s).
        For a Matrix object 'dm':
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
            return [self.matrix[index][i] for i in range(0, index)] + [self.matrix[i][index] for i in range(index, len(self))]
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
        Similar to __getitem__
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
            if isinstance(value, list) and all(isinstance(n, (int, long, float, complex)) for n in value):
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
            if isinstance(value, (int, long, float, complex)):
                if row_index > col_index:
                    self.matrix[row_index][col_index] = value
                else:
                    self.matrix[col_index][row_index] = value
            else:
                raise TypeError("Invalid value type.")
        else:
            raise TypeError("Invalid index type.")

    def __len__(self):
        """Matrix length"""
        return len(self.names)

    def __str__(self):
        """Get a lower triangular matrix string"""
        matrix_string = '\n'.join([self.names[i] + "\t" + "\t".join([str(n) for n in self.matrix[i]]) for i in range(0, len(self))])
        matrix_string = matrix_string + "\n\t" + "\t".join(self.names)
        return matrix_string


class DistanceMatrix(Matrix):
    """Distance matrix class that can be used for distance based tree
     algorithms.
    All diagonal elements will be zero no matter what the users provide.
    """

    def __init__(self, names, matrix=None):
        Matrix.__init__(self, names, matrix)
        self._set_zero_diagonal()

    def __setitem__(self, item, value):
        Matrix.__setitem__(self, item, value)
        self._set_zero_diagonal()

    def __delitem__(self, item):
        """Delete related distances by the index or name"""
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
        """Insert distances given the name and value"""
        if isinstance(name, str):
            # insert at the given index or at the end
            if index is None or isinstance(index, int):
                index = index and index or len(self)
            else:
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

    def _set_zero_diagonal(self):
        """set all diagonal elements to zero"""
        for i in range(0, len(self)):
            self.matrix[i][i] = 0


class DistanceCaluculator(object):
    """Class to calculate the distance matrix from a DNA or Protein
    Multiple Sequence Alignment(MSA) and the given name of the
    substitution model.
    Currently only scoring matrices are used.
    """

    dna_alphabet = ['A', 'T', 'C', 'G']

    # Identity scoring matrix
    identity = [[ 1],
                [ 0,  1],
                [ 0,  0,  1],
                [ 0,  0,  0,  1]]
    # BLAST nucleic acid scoring matrix
    blastn = [[ 5],
              [-4,  5],
              [-4, -4,  5],
              [-4, -4, -4,  5]]

    # transition/transversion scoring matrix
    trans = [[ 6],
             [-5,  6],
             [-5, -1,  6],
             [-1, -5, -5,  6]]

    protein_alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z']

    blosum40 = [[ 5],
                [-1,  5],
                [-2, -2, 16],
                [-1,  6, -2,  9],
                [-1,  1, -2,  2,  7],
                [-3, -3, -2, -4, -3,  9],
                [ 1, -1, -3, -2, -3, -3,  8],
                [-2,  0, -4,  0,  0, -2, -2, 13],
                [-1, -3, -4, -4, -4,  1, -4, -3,  6],
                [-1,  0, -3,  0,  1, -3, -2, -1, -3,  6],
                [-2, -3, -2, -3, -2,  2, -4, -2,  2, -2,  6],
                [-1, -3, -3, -3, -2,  0, -2,  1,  1, -1,  3,  7],
                [-1,  4, -2,  2, -1, -3,  0,  1, -2,  0, -3, -2,  8],
                [-2, -2, -5, -2,  0, -4, -1, -2, -2, -1, -4, -2, -2, 11],
                [ 0,  0, -4, -1,  2, -4, -2,  0, -3,  1, -2, -1,  1, -2,  8],
                [-2, -1, -3, -1, -1, -2, -3,  0, -3,  3, -2, -1,  0, -3,  2,  9],
                [ 1,  0, -1,  0,  0, -2,  0, -1, -2,  0, -3, -2,  1, -1,  1, -1,  5],
                [ 0,  0, -1, -1, -1, -1, -2, -2, -1,  0, -1, -1,  0,  0, -1, -2,  2,  6],
                [ 0, -3, -2, -3, -3,  0, -4, -4,  4, -2,  2,  1, -3, -3, -3, -2, -1,  1,  5],
                [-3, -4, -6, -5, -2,  1, -2, -5, -3, -2, -1, -2, -4, -4, -1, -2, -5, -4, -3, 19],
                [ 0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1,  0, -1, -2, -1, -1,  0,  0, -1, -2, -1],
                [-2, -3, -4, -3, -2,  4, -3,  2,  0, -1,  0,  1, -2, -3, -1, -1, -2, -1, -1,  3, -1,  9],
                [-1,  2, -3,  1,  5, -4, -2,  0, -4,  1, -2, -2,  0, -1,  4,  0,  0, -1, -3, -2, -1, -2,  5]]

    blosum62 = [[ 4],
                [-2,  4],
                [ 0, -3,  9],
                [-2,  4, -3,  6],
                [-1,  1, -4,  2,  5],
                [-2, -3, -2, -3, -3,  6],
                [ 0, -1, -3, -1, -2, -3,  6],
                [-2,  0, -3, -1,  0, -1, -2,  8],
                [-1, -3, -1, -3, -3,  0, -4, -3,  4],
                [-1,  0, -3, -1,  1, -3, -2, -1, -3,  5],
                [-1, -4, -1, -4, -3,  0, -4, -3,  2, -2,  4],
                [-1, -3, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5],
                [-2,  3, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6],
                [-1, -2, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7],
                [-1,  0, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5],
                [-1, -1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5],
                [ 1,  0, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4],
                [ 0, -1, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5],
                [ 0, -3, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4],
                [-3, -4, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11],
                [ 0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1,  0,  0, -1, -2, -1],
                [-2, -3, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2, -1,  7],
                [-1,  1, -3,  1,  4, -3, -2,  0, -3,  1, -3, -1,  0, -1,  3,  0,  0, -1, -2, -3, -1, -2,  4],]

    blosum90 = [[ 5],
                [-2,  4],
                [-1, -4,  9],
                [-3,  4, -5,  7],
                [-1,  0, -6,  1,  6],
                [-3, -4, -3, -5, -5,  7],
                [ 0, -2, -4, -2, -3, -5,  6],
                [-2, -1, -5, -2, -1, -2, -3,  8],
                [-2, -5, -2, -5, -4, -1, -5, -4,  5],
                [-1, -1, -4, -1,  0, -4, -2, -1, -4,  6],
                [-2, -5, -2, -5, -4,  0, -5, -4,  1, -3,  5],
                [-2, -4, -2, -4, -3, -1, -4, -3,  1, -2,  2,  7],
                [-2,  4, -4,  1, -1, -4, -1,  0, -4,  0, -4, -3,  7],
                [-1, -3, -4, -3, -2, -4, -3, -3, -4, -2, -4, -3, -3,  8],
                [-1, -1, -4, -1,  2, -4, -3,  1, -4,  1, -3,  0,  0, -2,  7],
                [-2, -2, -5, -3, -1, -4, -3,  0, -4,  2, -3, -2, -1, -3,  1,  6],
                [ 1,  0, -2, -1, -1, -3, -1, -2, -3, -1, -3, -2,  0, -2, -1, -1,  5],
                [ 0, -1, -2, -2, -1, -3, -3, -2, -1, -1, -2, -1,  0, -2, -1, -2,  1,  6],
                [-1, -4, -2, -5, -3, -2, -5, -4,  3, -3,  0,  0, -4, -3, -3, -3, -2, -1,  5],
                [-4, -6, -4, -6, -5,  0, -4, -3, -4, -5, -3, -2, -5, -5, -3, -4, -4, -4, -3, 11],
                [-1, -2, -3, -2, -2, -2, -2, -2, -2, -1, -2, -1, -2, -2, -1, -2, -1, -1, -2, -3, -2],
                [-3, -4, -4, -4, -4,  3, -5,  1, -2, -3, -2, -2, -3, -4, -3, -3, -3, -2, -3,  2, -2,  8],
                [-1,  0, -5,  0,  4, -4, -3,  0, -4,  1, -4, -2, -1, -2,  4,  0, -1, -1, -3, -4, -1, -3,  4]]

    pam90 = [[ 5],
             [-2,  5],
             [-5, -9,  9],
             [-2,  5 -10,  7],
             [-1,  2 -10,  3,  7],
             [-6, -8, -9 -11 -10,  8],
             [ 0, -2, -7, -2, -2, -7,  6],
             [-5,  0, -6, -2, -3, -4, -6,  8],
             [-3, -4, -4, -5, -4, -1, -7, -6,  7],
             [-5, -1 -10, -2, -3 -10, -5, -4, -4,  6],
             [-4, -7 -11, -9, -7, -1, -8, -4,  0, -6,  6],
             [-3, -6 -10, -7, -5, -2, -6, -7,  1,  0,  2, 10],
             [-2,  5, -7,  2,  0, -6, -1,  1, -4,  0, -5, -6,  6],
             [ 0, -4, -6, -5, -3, -7, -4, -2, -6, -4, -5, -6, -4,  7],
             [-3, -1 -10, -1,  2, -9, -5,  2, -5, -1, -3, -2, -2, -1,  7],
             [-5, -5, -6, -6, -6, -7, -7,  0, -4,  2, -6, -2, -3, -2,  0,  8],
             [ 1,  0, -1, -2, -2, -5,  0, -4, -4, -2, -6, -4,  1,  0, -3, -2,  5],
             [ 1, -2, -5, -3, -4, -6, -3, -5, -1, -2, -5, -2, -1, -2, -4, -4,  1,  6],
             [-1, -5, -4, -6, -4, -5, -4, -5,  3, -6, -1,  0, -5, -4, -5, -5, -4, -1,  6],
             [10, -8 -12 -11 -12, -3 -11, -5 -10, -8, -4, -9, -6 -10, -9,  0, -4, -9 -11, 13],
             [-2, -3, -6, -3, -3, -5, -3, -3, -3, -3, -4, -3, -2, -3, -3, -4, -2, -2, -3, -8, -3],
             [-6, -5, -2, -8, -7,  3 -10, -2, -4, -7, -5, -7, -3 -10, -8, -8, -5, -5, -5, -3, -5,  9],
             [-2,  1 -10,  2,  5 -10, -3,  0, -4, -2, -5, -4, -1, -2,  6, -2, -3, -4, -5 -11, -3, -7,  5]]

    pam120 = [[ 3],
              [ 0,  4],
              [-3, -6,  9],
              [ 0,  4, -7,  5],
              [ 0,  3, -7,  3,  5],
              [-4, -5, -6, -7, -7,  8],
              [ 1,  0, -4,  0, -1, -5,  5],
              [-3,  1, -4,  0, -1, -3, -4,  7],
              [-1, -3, -3, -3, -3,  0, -4, -4,  6],
              [-2,  0, -7, -1, -1, -7, -3, -2, -3,  5],
              [-3, -4, -7, -5, -4,  0, -5, -3,  1, -4,  5],
              [-2, -4, -6, -4, -3, -1, -4, -4,  1,  0,  3,  8],
              [-1,  3, -5,  2,  1, -4,  0,  2, -2,  1, -4, -3,  4],
              [ 1, -2, -4, -3, -2, -5, -2, -1, -3, -2, -3, -3, -2,  6],
              [-1,  0, -7,  1,  2, -6, -3,  3, -3,  0, -2, -1,  0,  0,  6],
              [-3, -2, -4, -3, -3, -5, -4,  1, -2,  2, -4, -1, -1, -1,  1,  6],
              [ 1,  0,  0,  0, -1, -3,  1, -2, -2, -1, -4, -2,  1,  1, -2, -1,  3],
              [ 1,  0, -3, -1, -2, -4, -1, -3,  0, -1, -3, -1,  0, -1, -2, -2,  2,  4],
              [ 0, -3, -3, -3, -3, -3, -2, -3,  3, -4,  1,  1, -3, -2, -3, -3, -2,  0,  5],
              [-7, -6, -8, -8, -8, -1, -8, -3, -6, -5, -3, -6, -4, -7, -6,  1, -2, -6, -8, 12],
              [-1, -1, -4, -2, -1, -3, -2, -2, -1, -2, -2, -2, -1, -2, -1, -2, -1, -1, -1, -5, -2],
              [-4, -3, -1, -5, -5,  4, -6, -1, -2, -5, -2, -4, -2, -6, -5, -5, -3, -3, -3, -2, -3,  8],
              [-1,  2, -7,  3,  4, -6, -2,  1, -3, -1, -3, -2,  0, -1,  4, -1, -1, -2, -3, -7, -1, -5,  4]]

    pam250 = [[ 2],
              [ 0,  3],
              [-2, -4, 12],
              [ 0,  3, -5,  4],
              [ 0,  3, -5,  3,  4],
              [-3, -4, -4, -6, -5,  9],
              [ 1,  0, -3,  1,  0, -5,  5],
              [-1,  1, -3,  1,  1, -2, -2,  6],
              [-1, -2, -2, -2, -2,  1, -3, -2,  5],
              [-1,  1, -5,  0,  0, -5, -2,  0, -2,  5],
              [-2, -3, -6, -4, -3,  2, -4, -2,  2, -3,  6],
              [-1, -2, -5, -3, -2,  0, -3, -2,  2,  0,  4,  6],
              [ 0,  2, -4,  2,  1, -3,  0,  2, -2,  1, -3, -2,  2],
              [ 1, -1, -3, -1, -1, -5,  0,  0, -2, -1, -3, -2,  0,  6],
              [ 0,  1, -5,  2,  2, -5, -1,  3, -2,  1, -2, -1,  1,  0,  4],
              [-2, -1, -4, -1, -1, -4, -3,  2, -2,  3, -3,  0,  0,  0,  1,  6],
              [ 1,  0,  0,  0,  0, -3,  1, -1, -1,  0, -3, -2,  1,  1, -1,  0,  2],
              [ 1,  0, -2,  0,  0, -3,  0, -1,  0,  0, -2, -1,  0,  0, -1, -1,  1,  3],
              [ 0, -2, -2, -2, -2, -1, -1, -2,  4, -2,  2,  2, -2, -1, -2, -2, -1,  0,  4],
              [-6, -5, -8, -7, -7,  0, -7, -3, -5, -3, -2, -4, -4, -6, -5,  2, -2, -5, -6, 17],
              [ 0, -1, -3, -1, -1, -2, -1, -1, -1, -1, -1, -1,  0, -1, -1, -1,  0,  0, -1, -4, -1],
              [-3, -3,  0, -4, -4,  7, -5,  0, -1, -4, -1, -2, -2, -5, -4, -4, -3, -3, -2,  0, -2, 10],
              [ 0,  2, -5,  3,  3, -5,  0,  2, -2,  0, -3, -2,  1,  0,  3,  0,  0, -1, -2, -6, -1, -4,  3]]

    # matrices available
    dna_matrices = {'identity': identity, 'blastn': blastn, 'trans': trans}
    protein_matrices = {'blosum40': blosum40, 'blosum62': blosum62, 'blosum90': blosum90,
                        'pam90': pam90, 'pam120': pam120, 'pam250': pam250}

    def __init__(self, msa, model):
        """Initialize with a MSA object and the distance model to be
        used()"""
        if isinstance(msa, MultipleSeqAlignment):
            self.msa = msa
        else:
            raise ValueError("Must provide a MultipleSeqAlignment object.")

        dna_keys = self.dna_matrices.keys()
        protein_keys = self.protein_matrices.keys()
        if model in dna_keys:
            self.scoring_matrix = Matrix(self.dna_alphabet, self.dna_matrices[model])
        elif model in protein_keys:
            self.scoring_matrix = Matrix(self.protein_alphabet, self.protein_matrices[model])
        else:
            raise ValueError("Model not supported. Available models: " + ", ".join(dna_keys + protein_keys))

    def _pairwise(self, seq1, seq2):
        """Calculate pairwise distance from two sequences"""
        max_score1 = 0
        max_score2 = 0
        score = 0
        skip_letters = ['-', '*']
        for i in range(0, len(seq1)):
            l1 = seq1[i]
            l2 = seq2[i]
            if l1 in skip_letters or l2 in skip_letters:
                continue
            if l1 not in self.scoring_matrix.names:
                raise ValueError("Bad alphabet '%s' in sequence '%s' at position '%s'" % (l1, seq1.id, i))
            if l2 not in self.scoring_matrix.names:
                raise ValueError("Bad alphabet '%s' in sequence '%s' at position '%s'" % (l2, seq2.id, i))
            max_score1 += self.scoring_matrix[l1, l1]
            max_score2 += self.scoring_matrix[l2, l2]
            score += self.scoring_matrix[l1, l2]

        max_score = max_score1 > max_score2 and max_score1 or max_score2

        return 1 - (score * 1.0 / max_score)

    def get_distance(self):
        """Return a DistanceMatrix for MSA object"""
        names = [s.id for s in self.msa]
        dm = DistanceMatrix(names)
        for seq1, seq2 in itertools.combinations(self.msa, 2):
            dm[seq1.id, seq2.id] = self._pairwise(seq1, seq2)
        return dm


class TreeContructor(object):
    """Base class for all tree constructor."""


class DistanceTreeConstructor(TreeContructor):
    """Distance based tree constructor"""
    def __init__(self, distance_matrix):
        self.distance_matrix = distance_matrix

    def upgma(self):
        """Construct and return an UPGMA(Unweighted Pair Group Method
        with Arithmetic mean) tree."""
        pass

    def nj(self):
        """Construct and return an Neighbor Joining tree."""
        pass


class ParsimonyTreeConstructor(TreeContructor):
    """Parsimony tree constructor"""
    def __init__(self, alignment):
        self.alignment = alignment

    def mp(self):
        """Construct and return an Maximum Parsimony tree."""
        pass

    def _parsimony_score(self, tree):
        """Calculate and return the parsimony score given a tree and
        the MSA."""
        pass

    def _nni(self, tree):
        """Search the best parsimony tree by using the NNI(Nearest
        Neighbor Interchanges) algorithm"""
        pass
