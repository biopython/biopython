# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SearchIO objects to model similarity search program outputs.

The SearchIO object model consists of a hierarchy of four nested objects:

    * QueryResult, to represent a search query.

      This is the top-level object returned by the main SearchIO ``parse`` and
      ``read`` functions. QueryResult objects may contain zero or more Hit
      objects, each accessible by its ID string (like in Python dictionaries)
      or integer index (like in Python lists).

    * Hit, to represent a database entry containing a full or partial sequence
      match with the query sequence.

      Hit objects contain one or more HSP objects, each accessible by its integer
      index. They behave very similar to a Python list.

    * HSP, to represent a region of significant alignment(s) between the query
      and hit sequences.

      HSP objects contain one or more HSPFragment objects, each accessible by
      its integer index. In most cases, the HSP objects are where the bulk of
      search result statistics (e.g. e-value, bitscore) are stored. Like Hit
      objects, HSPs also behave very similar to a Python list.

    * HSPFragment, to represent a single contiguous alignment between the query
      and hit sequences.

      HSPFragment objects may store hit and query sequences resulting from the
      sequence search. If present, these sequences are stored as SeqRecord
      objects (see SeqRecord). If both of them are present, HSPFragment will
      create a MultipleSeqAlignment object from both sequences.

      Most search programs only have HSPs with one HSPFragment in them, making
      these two objects inseparable. However, there are programs (e.g. BLAT and
      Exonerate) which may have more than one HSPFragment objects in any given
      HSP. If you are not using these programs, you can safely consider HSP and
      HSPFragment as a single union.

"""

from .query import QueryResult
from .hit import Hit
from .hsp import HSP, HSPFragment


__all__ = ("QueryResult", "Hit", "HSP", "HSPFragment")


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
