# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO objects to model homology search program outputs (PRIVATE).

"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


class Result(object):

    """Class representing search results from a single query.

    """


class Hit(object):

    """Class representing the entire database entry of a sequence match.

    """


class HSP(object):

    """Class representing high-scoring alignment regions of the query and hit.

    """
