# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""A superclass for reading [f]ixed-column type [f]lat-[f]ile records."""


class fff_rec(object):
    def __init__(self, inrec=''):
        """Initialize the class."""
        self.data = inrec

    def __repr__(self):
        return str(self.data)
    __str__ = __repr__

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.data[index]
        elif (isinstance(index, tuple) or isinstance(index, list)) \
        and len(index) == 2:
            # Not sure if this is needed anymore:
            return self.data[index[0]:index[1]]
        else:
            return self.data[index]


class align(object):
    """Definition of the align section in a FSSP file."""

    abs_res_num = (0, 4)
    pdb_res_num = (4, 9)
    chain_id = 10
    res_name = 12
    ss1 = 15
    turn3 = 17
    turn4 = 18
    turn5 = (20, 22)
    acc = (34, 37)
    start_aa_list = 42
