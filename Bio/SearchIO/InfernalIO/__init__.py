# Copyright 2024 by Samuel Prince. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


"""Bio.SearchIO support for Infernal output formats.

This module adds support for parsing Infernal outputs. Infernal is a
suite of programs for searching DNA sequence databases for RNA structure
and sequence similarities using covariance models (CMs).

Bio.SearchIO.InfernalIO was tested on the following Infernal versions and flavors:

    - Infernal (1.0.0+): cmscan and cmsearch

More information on HMMER are available through these links:
  - Web page: http://eddylab.org/infernal/
  - User guide: http://eddylab.org/infernal/Userguide.pdf


Supported formats
=================

Bio.SearchIO.InfernalIO supports the following Infernal output formats:
    - Plain text   - 'infernal-text'    - parsing, indexing
    - Tabular      - 'infernal-tab'     - parsing, indexing

For all output formats, Infernal uses 'mdl' for 'query' and 'seq' for 'hit'.
InfernalIO is aware of this different naming scheme, and will use 'query'
and 'hit' to fit SearchIO's object model.

Infernal sometime reports 'local ends' (i.e., a large insertion or deletion
in the optimal alignment), which are expresented by a number in brackets in
the alignment (ex. AUUAC*[88]*GUAGU). In InfernalIO, these local alignment
are split into fragments of the same HSP.

infernal-text
=============

The Infernal plain text parser supports output files with alignment
blocks (default) or without (with the '-noali' flag). If the alignment blocks
are present, it can parse files with variable alignment width (using the
'-notextw' or '-textw' flag). Both CM or HMM searches (with the '--hmmonly'
flag) output are supported. The parser only supports non-verbose output formats.

The following SearchIO objects attributes are provided.

+-----------------+-------------------------+----------------------------------+
| Object          | Attribute               | Value                            |
+=================+=========================+==================================+
| QueryResult     | accession               | query accession                  |
|                 +-------------------------+----------------------------------+
|                 | description             | query sequence description       |
|                 +-------------------------+----------------------------------+
|                 | id                      | query sequence ID                |
|                 +-------------------------+----------------------------------+
|                 | program                 | Infernal flavor                  |
|                 +-------------------------+----------------------------------+
|                 | seq_len                 | full length of query sequence    |
|                 +-------------------------+----------------------------------+
|                 | target                  | target search database           |
|                 +-------------------------+----------------------------------+
|                 | version                 | Infernal version                 |
+-----------------+-------------------------+----------------------------------+
| Hit             | description             | hit sequence description         |
|                 +-------------------------+----------------------------------+
|                 | id                      | hit sequence ID                  |
+-----------------+-------------------------+----------------------------------+
| HSP             | evalue                  | hsp evalue                       |
|                 +-------------------------+----------------------------------+
|                 | bias                    | hsp bias                         |
|                 +-------------------------+----------------------------------+
|                 | bitscore                | hsp score                        |
|                 +-------------------------+----------------------------------+
|                 | gc                      | gc fraction                      |
|                 +-------------------------+----------------------------------+
|                 | is_included             | boolean, whether the hit of the  |
|                 |                         | hsp is in the inclusion          |
|                 |                         | threshold                        |
|                 +-------------------------+----------------------------------+
|                 | query_start             | query start position             |
|                 +-------------------------+----------------------------------+
|                 | query_end               | query end position               |
|                 +-------------------------+----------------------------------+
|                 | query_endtype           | query sequence end types (e.g.,  |
|                 |                         | '[]', '..', '[.', '.]', etc.)    |
|                 +-------------------------+----------------------------------+
|                 | hit_start               | hit start position               |
|                 +-------------------------+----------------------------------+
|                 | hit_end                 | hit end position                 |
|                 +-------------------------+----------------------------------+
|                 | hit_endtype             | hit sequence end types           |
|                 +-------------------------+----------------------------------+
|                 | acc_avg                 | expected accuracy per alignment  |
|                 |                         | residue (acc column)             |
|                 +-------------------------+----------------------------------+
|                 | model                   | type of model used (cm or hmm)   |
|                 +-------------------------+----------------------------------+
|                 | truncated               | indicate if the hit is truncated |
|                 |                         | (5', 3' or both) or not          |
+-----------------+-------------------------+----------------------------------+
| HSPFragment     | aln_annotation          | alignment similarity string and  |
|                 |                         | other annotations (PP, CS,       |
|                 |                         | similarity and NC (except for    |
|                 |                         | --hmmonly))                      |
|                 +-------------------------+----------------------------------+
|                 | aln_span                | length of alignment fragment     |
|                 +-------------------------+----------------------------------+
|                 | hit                     | hit sequence                     |
|                 +-------------------------+----------------------------------+
|                 | hit_start               | local alignment sequence start   |
|                 |                         | coordinate (seq from)            |
|                 +-------------------------+----------------------------------+
|                 | hit_end                 | local alignment sequence end     |
|                 |                         | coordinate (seq to)              |
|                 +-------------------------+----------------------------------+
|                 | hit_strand              | hit sequence strand              |
|                 +-------------------------+----------------------------------+
|                 | query                   | query sequence                   |
|                 +-------------------------+----------------------------------+
|                 | query_start             | local model alignment start      |
|                 |                         | coordinate (mdl from)            |
|                 +-------------------------+----------------------------------+
|                 | query_end               | local model alignment end        |
|                 |                         | coordinate (mdl to)              |
+-----------------+-------------------------+----------------------------------+

infernal-tab
============

The Infernal plain text parser supports the standard cmsearch tabular output and
cmscan tabular output files formats 1, 2 and 3 (inferred automatically from the
header).

Rows marked with '*' denotes attributes not available in the default format.

+-----------------+-------------------------+----------------------------------+
| Object          | Attribute               | Value                            |
+=================+=========================+==================================+
| QueryResult     | accession               | query accession                  |
|                 +-------------------------+----------------------------------+
|                 | id                      | query sequence ID                |
|                 +-------------------------+----------------------------------+
|                 | clan*                   | Rfam clan                        |
|                 +-------------------------+----------------------------------+
|                 | seq_len*                | query sequence length            |
+-----------------+-------------------------+----------------------------------+
| Hit             | description             | hit sequence description         |
|                 +-------------------------+----------------------------------+
|                 | id                      | hit sequence ID                  |
|                 +-------------------------+----------------------------------+
|                 | accession               | hit accession                    |
|                 +-------------------------+----------------------------------+
|                 | seq_len*                | hit sequence length              |
+-----------------+-------------------------+----------------------------------+
| HSP             | evalue                  | hsp evalue                       |
|                 +-------------------------+----------------------------------+
|                 | bias                    | hsp bias                         |
|                 +-------------------------+----------------------------------+
|                 | bitscore                | hsp score                        |
|                 +-------------------------+----------------------------------+
|                 | gc                      | gc fraction                      |
|                 +-------------------------+----------------------------------+
|                 | is_included             | boolean, whether the hit of the  |
|                 |                         | hsp is in the inclusion          |
|                 |                         | threshold                        |
|                 +-------------------------+----------------------------------+
|                 | model                   | type of model used (cm or hmm)   |
|                 +-------------------------+----------------------------------+
|                 | truncated               | indicate if the hit is truncated |
|                 |                         | (5', 3' or both) or not          |
|                 +-------------------------+----------------------------------+
|                 | pipeline_pass           | pipeline pass at which the hit   |
|                 |                         | was identified                   |
|                 +-------------------------+----------------------------------+
|                 | olp*                    | overlap status of this hit ('*', |
|                 |                         | '^', '$' or '=')                 |
|                 +-------------------------+----------------------------------+
|                 | anyidx*                 | index of the best scoring        |
|                 |                         | overlapping hit (or none if      |
|                 |                         | there are no overlap)            |
|                 +-------------------------+----------------------------------+
|                 | afrct1*                 | fraction of this hit that        |
|                 |                         | overlap with anyidx hit (or none |
|                 |                         | if there are no overlap)         |
|                 +-------------------------+----------------------------------+
|                 | afrct2*                 | fraction of anyidx hit           |
|                 |                         | with this hit (or none if there  |
|                 |                         | are no overlap)                  |
|                 +-------------------------+----------------------------------+
|                 | winidx*                 | index of the best scoring hit    |
|                 |                         | that overlaps with this hit that |
|                 |                         | is marked as '^' (or none if     |
|                 |                         | there are no overlap)            |
|                 +-------------------------+----------------------------------+
|                 | wfrct1*                 | fraction of this hit that        |
|                 |                         | overlap with winidx hit (or none |
|                 |                         | if there are no overlap)         |
|                 +-------------------------+----------------------------------+
|                 | wfrct2*                 | fraction of winidx hit           |
|                 |                         | with this hit (or none if there  |
|                 |                         | are no overlap)                  |
+-----------------+-------------------------+----------------------------------+
| HSPFragment     | hit_start               | local alignment sequence start   |
| (also via HSP)  |                         | coordinate (seq from)            |
|                 +-------------------------+----------------------------------+
|                 | hit_end                 | local alignment sequence end     |
|                 |                         | coordinate (seq to)              |
|                 +-------------------------+----------------------------------+
|                 | hit_strand              | hit sequence strand              |
|                 +-------------------------+----------------------------------+
|                 | query_start             | local model alignment start      |
|                 |                         | coordinate (mdl from)            |
|                 +-------------------------+----------------------------------+
|                 | query_end               | local model alignment end        |
|                 |                         | coordinate (mdl to)              |
+-----------------+-------------------------+----------------------------------+
"""


from .infernal_tab import InfernalTabParser
from .infernal_tab import InfernalTabIndexer
from .infernal_text import InfernalTextParser
from .infernal_text import InfernalTextIndexer


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
