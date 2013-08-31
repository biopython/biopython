# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Classes and methods for tree construction"""
__docformat__ = "restructuredtext en"

import itertools
import copy
from Bio.Phylo import BaseTree
from Bio.Align import MultipleSeqAlignment
from Bio.SubsMat.MatrixInfo import *

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

    :Parameters:
        names : list
            names of elements, used for indexing 
        matrix : list
            nested list of numerical lists in lower triangular format

    Example
    -------

    >>> from Bio.Phylo.TreeConstruction import Matrix
    >>> names = ['Alpha', 'Beta', 'Gamma', 'Delta']
    >>> matrix = [[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]]
    >>> m = Matrix(names, matrix)
    >>> m
    Matrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], matrix=[[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]])
    
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

    Also you can delete or insert a column&row of elemets by index.

    >>> m
    Matrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], matrix=[[0], [7, 0], [8, 4, 0], [9, 5, 6, 0]])
    >>> del m['Alpha']
    >>> m
    Matrix(names=['Beta', 'Gamma', 'Delta'], matrix=[[0], [4, 0], [5, 6, 0]])
    >>> m.insert('Alpha', [0, 7, 8, 9] , 0)
    >>> m
    Matrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], matrix=[[0], [7, 0], [8, 4, 0], [9, 5, 6, 0]])

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
        """Matrix length"""
        return len(self.names)
    
    def __repr__(self):
        return self.__class__.__name__ \
        + "(names=%s, matrix=%s)" \
        % tuple(map(repr, (self.names, self.matrix)))

    def __str__(self):
        """Get a lower triangular matrix string"""
        matrix_string = '\n'.join([self.names[i] + "\t" + 
            "\t".join([str(n) for n in self.matrix[i]]) for i in range(0, len(self))])
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

    def _set_zero_diagonal(self):
        """set all diagonal elements to zero"""
        for i in range(0, len(self)):
            self.matrix[i][i] = 0


class DistanceCalculator(object):
    """Class to calculate the distance matrix from a DNA or Protein
    Multiple Sequence Alignment(MSA) and the given name of the
    substitution model.

    Currently only scoring matrices are used.

    :Parameters:
        msa : MultipleSeqAlignment
            DNA or Protein multiple sequence alignment.  
        model : str
            Name of the model matrix to be used to calculate distance.
            The attribute `dna_matrices` contains the available model
            names for DNA sequences and `protein_matrices` for protein
            sequences.

    Example
    -------

    >>> from Bio.Phylo.TreeConstruction import DistanceCalculator
    >>> from Bio import AlignIO
    >>> aln = AlignIO.read(open('Tests/TreeConstruction/msa.phy'), 'phylip')
    >>> print aln
    SingleLetterAlphabet() alignment with 5 rows and 13 columns
    AACGTGGCCACAT Alpha
    AAGGTCGCCACAC Beta
    GAGATTTCCGCCT Delta
    GAGATCTCCGCCC Epsilon
    CAGTTCGCCACAA Gamma

    DNA calculator with 'identity' model:

    >>> calculator = DistanceCalculator('identity')
    >>> dm = calculator.get_distance(aln)
    >>> print dm
    Alpha   0
    Beta    0.230769230769  0
    Gamma   0.384615384615  0.230769230769  0
    Delta   0.538461538462  0.538461538462  0.538461538462  0
    Epsilon 0.615384615385  0.384615384615  0.461538461538  0.153846153846  0
            Alpha           Beta            Gamma           Delta           Epsilon

    Protein calculator with 'blosum62' model:
    >>> calculator = DistanceCalculator('blosum62')
    >>> dm = calculator.get_distance(aln)
    >>> print dm
    Alpha   0
    Beta    0.369047619048  0
    Gamma   0.493975903614  0.25            0
    Delta   0.585365853659  0.547619047619  0.566265060241  0
    Epsilon 0.7             0.355555555556  0.488888888889  0.222222222222  0
            Alpha           Beta            Gamma           Delta           Epsilon
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

    # matrices available
    dna_matrices = {'identity': identity, 'blastn': blastn, 'trans': trans}
    protein_matrices = {'benner6': benner6, 'benner22': benner22, 'benner74': benner74, 'blosum100': blosum100,
                      'blosum30': blosum30, 'blosum35': blosum35, 'blosum40': blosum40, 'blosum45': blosum45,
                      'blosum50': blosum50, 'blosum55': blosum55, 'blosum60': blosum60, 'blosum62': blosum62,
                      'blosum65': blosum65, 'blosum70': blosum70, 'blosum75': blosum75, 'blosum80': blosum80,
                      'blosum85': blosum85, 'blosum90': blosum90, 'blosum95': blosum95, 'feng': feng,
                      'fitch': fitch, 'genetic': genetic, 'gonnet': gonnet, 'grant': grant,
                      'ident': ident, 'johnson': johnson, 'levin': levin, 'mclach': mclach,
                      'miyata': miyata, 'nwsgappep': nwsgappep, 'pam120': pam120, 'pam180': pam180,
                      'pam250': pam250, 'pam30': pam30, 'pam300': pam300, 'pam60': pam60,
                      'pam90': pam90, 'rao': rao, 'risler': risler, 'structure': structure
                     }

    def __init__(self, model):
        """Initialize with a distance model"""

        dna_keys = self.dna_matrices.keys()
        protein_keys = self.protein_matrices.keys()
        if model in dna_keys:
            self.scoring_matrix = Matrix(self.dna_alphabet, self.dna_matrices[model])
        elif model in protein_keys:
            self.scoring_matrix = self._build_protein_matrix(self.protein_matrices[model])
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

    def get_distance(self, msa):
        """Return a DistanceMatrix for MSA object"""

        if not isinstance(msa, MultipleSeqAlignment):
            raise TypeError("Must provide a MultipleSeqAlignment object.")

        names = [s.id for s in msa]
        dm = DistanceMatrix(names)
        for seq1, seq2 in itertools.combinations(msa, 2):
            dm[seq1.id, seq2.id] = self._pairwise(seq1, seq2)
        return dm

    def _build_protein_matrix(self, subsmat):
        """Convert matrix from SubsMat format to Matrix object"""
        protein_matrix = Matrix(self.protein_alphabet)
        for k, v in subsmat.items():
            aa1, aa2 = k
            protein_matrix[aa1, aa2] = v
        return protein_matrix


class TreeContructor(object):
    """Base class for all tree constructor."""

    def build_tree(self):
        """Caller to built the tree. This should be implemented
        in subclass"""
        raise NotImplementedError("Method not implemented!")


class DistanceTreeConstructor(TreeContructor):
    """Distance based tree constructor.

    :Parameters:
        distance_matrix : DistanceMatrix
            The distance matrix for tree construction.

    Example
    --------

    >>> from TreeConstruction import DistanceTreeConstructor
    >>> constructor = DistanceTreeConstructor(dm)

    UPGMA Tree: 

    >>> upgmatree = constructor.upgma()
    >>> print upgmatree
    Tree(rooted=True)
        Clade(name='Inner4')
            Clade(branch_length=0.171955155115, name='Inner1')
                Clade(branch_length=0.111111111111, name='Epsilon')
                Clade(branch_length=0.111111111111, name='Delta')
            Clade(branch_length=0.0673103855608, name='Inner3')
                Clade(branch_length=0.0907558806655, name='Inner2')
                    Clade(branch_length=0.125, name='Gamma')
                    Clade(branch_length=0.125, name='Beta')
                Clade(branch_length=0.215755880666, name='Alpha')

    NJ Tree:

    >>> njtree = constructor.nj()
    >>> print njtree
    Tree(rooted=False)
        Clade(name='Inner3')
            Clade(branch_length=0.0142054862889, name='Inner2')
                Clade(branch_length=0.239265540676, name='Inner1')
                    Clade(branch_length=0.0853101915988, name='Epsilon')
                    Clade(branch_length=0.136912030623, name='Delta')
                Clade(branch_length=0.292306275042, name='Alpha')
            Clade(branch_length=0.0747705106139, name='Beta')
            Clade(branch_length=0.175229489386, name='Gamma')


    """
    def __init__(self, distance_matrix):
        if isinstance(distance_matrix, DistanceMatrix):
            self.distance_matrix = distance_matrix
        else:
            raise TypeError("Must provide a DistanceMatrix object.")

    def upgma(self):
        """Construct and return an UPGMA(Unweighted Pair Group Method
        with Arithmetic mean) tree."""
        # make a copy of the distance matrix to be used
        dm = copy.deepcopy(self.distance_matrix)
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
            #assign branch length
            if clade1.is_terminal():
                clade1.branch_length = min_dist * 1.0 / 2
            else:
                clade1.branch_length = min_dist * 1.0 / 2 - self._height_of(clade1)

            if clade2.is_terminal():
                clade2.branch_length = min_dist * 1.0 / 2
            else:
                clade2.branch_length = min_dist * 1.0 / 2 - self._height_of(clade2)

            # update node list
            clades[min_j] = inner_clade
            del clades[min_i]

            # rebuild distance matrix,
            # set the distances of new node at the index of min_j
            for k in range(0, len(dm)):
                if k != min_i and k != min_j:
                    dm[min_j, k] = (dm[min_i, k] + dm[min_j, k]) * 1.0 / 2

            dm.names[min_j] = "Inner" + str(inner_count)

            del dm[min_i]
        return BaseTree.Tree(inner_clade)

    def nj(self):
        """Construct and return an Neighbor Joining tree."""

        # make a copy of the distance matrix to be used
        dm = copy.deepcopy(self.distance_matrix)
        # init terminal clades
        clades = [BaseTree.Clade(None, name) for name in dm.names]
        # init node distance
        node_dist = [0] * len(dm)
        # init minimum index
        min_i = 0
        min_j = 0
        inner_count = 0
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
            #assign branch length
            clade1.branch_length = (dm[min_i, min_j] + node_dist[min_i] - node_dist[min_j]) / 2
            clade2.branch_length = dm[min_i, min_j] - clade1.branch_length

            # update node list
            clades[min_j] = inner_clade
            del clades[min_i]

            # rebuild distance matrix,
            # set the distances of new node at the index of min_j
            for k in range(0, len(dm)):
                if k != min_i and k != min_j:
                    dm[min_j, k] = (dm[min_i, k] + dm[min_j, k] - dm[min_i, min_j]) * 1.0 / 2

            dm.names[min_j] = "Inner" + str(inner_count)

            del dm[min_i]

        # set the last clade as one of the child of the inner_clade
        if clades[0] == inner_clade:
            clades[1].branch_length = dm[1, 0]
            clades[0].clades.append(clades[1])
        else:
            clades[0].branch_length = dm[1, 0]
            clades[1].clades.append(clades[0])

        return BaseTree.Tree(clades[0], rooted=False)

    def _height_of(self, clade):
        """calculate height of the clade -- the longest path to one of
        the terminals"""
        height = 0
        if clade.is_terminal():
            height = clade.branch_length
        else:
            height = height + max([self._height_of(c) for c in clade.clades])
        return height

######################### Tree Scoring and Searching Classes ##########################

class Scorer(object):
    """Base class for all tree scoring methods"""

    def get_score(self, tree, alignment):
        """Caller to get the score of a tree for the given alignment.
        This should be implemented in subclass"""
        raise NotImplementedError("Method not implemented!")


class TreeSearcher(object):
    """Base class for all tree searching methods"""

    def search(self, starting_tree):
        """Caller to search the best tree with a starting tree.
        This should be implemented in subclass"""
        raise NotImplementedError("Method not implemented!")


class NNITreeSearcher(TreeSearcher):
    """Tree searching class of NNI(Nearest Neighbor Interchanges)
     algorithm.

    :Parameters:
        scorer: ParsimonyScorer
            parsimony scorer to calculate the parsimony score of
            different trees during NNI algorithm.
     """

    def __init__(self, scorer):
        if isinstance(scorer, Scorer):
            self.scorer = scorer
        else:
            raise TypeError("Must provide a Scorer object.")

    def search(self, starting_tree, alignment):
        """Implement the TreeSearcher.search method.

        :Parameters:
           starting_tree : Tree
               starting tree of NNI method.
           alignment: MultipleSeqAlignment
               multiple sequence alignment used to calculate parsimony
               score of different NNI trees.
        """

        return self._nni(starting_tree, alignment)

    def _nni(self, starting_tree, alignment):
        """Search the best parsimony tree by using the NNI(Nearest
        Neighbor Interchanges) algorithm"""
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
            if best_score == temp:
                break
        return best_tree

    def _get_neighbors(self, tree):
        """Get all neighbor trees of the given tree(currently only
         for binary rooted tree)"""
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
                    #left_left = left.clades[0]
                    left_right = left.clades[1]
                    right_left = right.clades[0]
                    right_right = right.clades[1]
                    # neightbor 1 (left_left + right_right)
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

######################### Parsimony Classes ##########################


class ParsimonyScorer(Scorer):
    """Parsimony scorer with a scoring matrix.

    This is a combination of Fitch algorithm and Sankoff algorithm. 
    See ParsimonyTreeConstructor for usage.

    :Parameters:
        matrix: Matrix
            scoring matrix used in parsimony score calculation.
    """
    def __init__(self, matrix=None):
        if not matrix or isinstance(matrix, Matrix):
            self.matrix = matrix
        else:
            raise TypeError("Must provide a Matrix object.")

    def get_score(self, tree, alignment):
        """Calculate and return the parsimony score given a tree and
        the MSA using the Fitch algorithm without the penalty matrix
        the Sankoff algorithm with the matrix"""
        # make sure the tree is rooted and bifurcating
        if not tree.is_bifurcating():
            raise ValueError("The tree provided should be bifurcating.")
        if not tree.rooted:
            tree.root_at_midpoint()
        # sort tree terminals and alignment
        terms = tree.get_terminals()
        terms.sort(key=lambda term: term.name)
        alignment.sort()
        if not all([t.name == a.id for t, a in zip(terms, alignment)]):
            raise ValueError("Taxon names of the input tree should be the same with the alignment.")
        #term_align = dict(zip(terms, alignment))
        score = 0
        for i in range(len(alignment[0])):
            # parsimony score for column_i
            score_i = 0
            # get column
            column_i = alignment[:, i]
            # skip non-informative column
            if column_i == len(column_i) * column_i[0]:
                continue

            #### start calculating score_i using the tree and column_i

            ## Fitch algorithm without the penalty matrix
            if not self.matrix:
                # init by mapping terminal clades and states in column_i
                clade_states = dict(zip(terms, [set([c]) for c in column_i]))
                for clade in tree.get_nonterminals(order="postorder"):
                    clade_childs = clade.clades
                    left_state = clade_states[clade_childs[0]]
                    right_state = clade_states[clade_childs[1]]
                    state = left_state & right_state
                    if not state:
                        state = left_state | right_state
                        score_i = score_i + 1
                    clade_states[clade] = state
            ## Sankoff algorithm with the penalty matrix
            else:
                inf = float('inf')
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
            score = score + score_i
        return score

class ParsimonyTreeConstructor(TreeContructor):
    """Parsimony tree constructor.

    :Parameters:
        alignment: MultipleSeqAlignment
            multiple sequence alignment to calculate parsimony tree.
        searcher: TreeSearcher
            tree searcher to search the best parsimony tree.
        starting_tree: Tree
            starting tree provided to the searcher.

    Example
    --------

    >>> from Bio import AlignIO
    >>> from TreeConstruction import *
    >>> aln = AlignIO.read(open('Tests/TreeConstruction/msa.phy'), 'phylip')
    >>> print aln
    SingleLetterAlphabet() alignment with 5 rows and 13 columns
    AACGTGGCCACAT Alpha
    AAGGTCGCCACAC Beta
    GAGATTTCCGCCT Delta
    GAGATCTCCGCCC Epsilon
    CAGTTCGCCACAA Gamma
    >>> starting_tree = Phylo.read('Tests/TreeConstruction/nj.tre', 'newick')
    >>> print tree
    Tree(weight=1.0, rooted=False)
        Clade(branch_length=0.0, name='Inner3')
            Clade(branch_length=0.01421, name='Inner2')
                Clade(branch_length=0.23927, name='Inner1')
                    Clade(branch_length=0.08531, name='Epsilon')
                    Clade(branch_length=0.13691, name='Delta')
                Clade(branch_length=0.29231, name='Alpha')
            Clade(branch_length=0.07477, name='Beta')
            Clade(branch_length=0.17523, name='Gamma')
    >>> from TreeConstruction import *
    >>> scorer = ParsimonyScorer()
    >>> searcher = NNITreeSearcher(scorer)
    >>> constructor = ParsimonyTreeConstructor(aln, searcher, starting_tree)
    >>> pars_tree = constructor.build_tree()
    >>> print pars_tree
    Tree(weight=1.0, rooted=True)
        Clade(branch_length=0.0)
            Clade(branch_length=0.197335, name='Inner1')
                Clade(branch_length=0.13691, name='Delta')
                Clade(branch_length=0.08531, name='Epsilon')
            Clade(branch_length=0.041935, name='Inner2')
                Clade(branch_length=0.01421, name='Inner3')
                    Clade(branch_length=0.17523, name='Gamma')
                    Clade(branch_length=0.07477, name='Beta')
                Clade(branch_length=0.29231, name='Alpha')
    """
    def __init__(self, alignment, searcher, starting_tree=None):
        self.alignment = alignment
        self.searcher = searcher
        self.starting_tree = starting_tree

    def build_tree(self):
        return self.searcher.search(self.starting_tree, self.alignment)