# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SearchIO support for HMMER output formats.

This module adds support for parsing HMMER outputs. HMMER is a
suite of programs implementing the profile hidden Markov models to find
similarity across protein sequences.

Bio.SearchIO.HmmerIO was tested on the following HMMER versions and flavors:

    - HMMER3 flavors: hmmscan, hmmsearch, phmmer
    - HMMER2 flavors: hmmpfam, hmmsearch

More information on HMMER are available through these links:
  - Web page: http://hmmer.janelia.org/
  - User guide: ftp://selab.janelia.org/pub/software/hmmer3/3.0/Userguide.pdf


Supported formats
=================

Bio.SearchIO.HmmerIO supports the following HMMER output formats:

    - Plain text, v3.0   - 'hmmer3-text'    - parsing, indexing
    - Table, v3.0        - 'hmmer3-tab'     - parsing, indexing, writing
    - Domain table, v3.0 - 'hmmer3-domtab'* - parsing, indexing, writing
    - Plain text, v2.x   - 'hmmer2-text'    - parsing, indexing

* For the domain table output, due to the way HMMER outputs the sequence
  coordinates, you have to specify what HMMER flavor produced the output as the
  file format. So instead of using 'hmmer3-domtab', you have to use either
  'hmmscan3-domtab', 'hmmsearch3-domtab', or 'phmmer3-domtab' as the file format
  name.

Note that for all output formats, HMMER uses its own convention of input and
output coordinates. It does not use the term 'hit' or 'query', instead it
uses 'hmm' or 'ali'. For example, 'hmmfrom' is the start coordinate of the HMM
sequence while 'alifrom' is the start coordinate of the protein sequence.

HmmerIO is aware of this different naming scheme and will adjust them
accordingly to fit SearchIO's object model. If HmmerIO sees that the output file
to parse was written by hmmsearch or phmmer, all 'hmm' coordinates will be the
hit coordinates and 'ali' coordinates will be the query coordinates. Conversely,
if the HMMER flavor is hmmscan, 'hmm' will be query and 'ali' will be hit.

This is why the 'hmmer3-domtab' format has to be specified with the source HMMER
flavor. The parsers need to know which is the hit and which is the query.
'hmmer3-text' has its source program information present in the file, while
'hmmer3-tab' does not output any coordinates. That's why both of these formats
do not need direct flavor specification like 'hmmer3-domtab'.

Also note that when using the domain table format writers, it will use HMMER's
naming convention ('hmm' and 'ali') so the files you write will be similar to
files written by a real HMMER program.


hmmer2-text and hmmer3-text
===========================

The parser for HMMER 3.0 plain text output can parse output files with alignment
blocks (default) or without (with the '--noali' flag). If the alignment blocks
are present, you can also parse files with variable alignment width (using the
'--notextw' or '--textw' flag).

The following SearchIO objects attributes are provided. Rows marked with '*'
denotes attributes not available in the hmmer2-text format:

+-----------------+-------------------------+----------------------------------+
| Object          | Attribute               | Value                            |
+=================+=========================+==================================+
| QueryResult     | accession               | accession (if present)           |
|                 +-------------------------+----------------------------------+
|                 | description             | query sequence description       |
|                 +-------------------------+----------------------------------+
|                 | id                      | query sequence ID                |
|                 +-------------------------+----------------------------------+
|                 | program                 | HMMER flavor                     |
|                 +-------------------------+----------------------------------+
|                 | seq_len*                | full length of query sequence    |
|                 +-------------------------+----------------------------------+
|                 | target                  | target search database           |
|                 +-------------------------+----------------------------------+
|                 | version                 | BLAST version                    |
+-----------------+-------------------------+----------------------------------+
| Hit             | bias*                   | hit-level bias                   |
|                 +-------------------------+----------------------------------+
|                 | bitscore                | hit-level score                  |
|                 +-------------------------+----------------------------------+
|                 | description             | hit sequence description         |
|                 +-------------------------+----------------------------------+
|                 | domain_exp_num*         | expected number of domains in    |
|                 |                         | the hit (exp column)             |
|                 +-------------------------+----------------------------------+
|                 | domain_obs_num          | observed number of domains in    |
|                 |                         | the hit (N column)               |
|                 +-------------------------+----------------------------------+
|                 | evalue                  | hit-level e-value                |
|                 +-------------------------+----------------------------------+
|                 | id                      | hit sequence ID                  |
|                 +-------------------------+----------------------------------+
|                 | is_included*            | boolean, whether the hit is in   |
|                 |                         | the inclusion threshold or not   |
+-----------------+-------------------------+----------------------------------+
| HSP             | acc_avg*                | expected accuracy per alignment  |
|                 |                         | residue (acc column)             |
|                 +-------------------------+----------------------------------+
|                 | bias*                   | hsp-level bias                   |
|                 +-------------------------+----------------------------------+
|                 | bitscore                | hsp-level score                  |
|                 +-------------------------+----------------------------------+
|                 | domain_index            | the domain index set by HMMER    |
|                 +-------------------------+----------------------------------+
|                 | env_end*                | end coordinate of the envelope   |
|                 +-------------------------+----------------------------------+
|                 | env_endtype*            | envelope end types (e.g. '[]',   |
|                 |                         | '..', '[.', etc.)                |
|                 +-------------------------+----------------------------------+
|                 | env_start*              | start coordinate of the envelope |
|                 +-------------------------+----------------------------------+
|                 | evalue                  | hsp-level independent e-value    |
|                 +-------------------------+----------------------------------+
|                 | evalue_cond*            | hsp-level conditional e-value    |
|                 +-------------------------+----------------------------------+
|                 | hit_endtype             | hit sequence end types           |
|                 +-------------------------+----------------------------------+
|                 | is_included*            | boolean, whether the hit of the  |
|                 |                         | hsp is in the inclusion          |
|                 |                         | threshold                        |
|                 +-------------------------+----------------------------------+
|                 | query_endtype           | query sequence end types         |
+-----------------+-------------------------+----------------------------------+
| HSPFragment     | aln_annotation          | alignment similarity string and  |
| (also via HSP)  |                         | other annotations (e.g. PP, CS)  |
|                 +-------------------------+----------------------------------+
|                 | aln_span                | length of alignment fragment     |
|                 +-------------------------+----------------------------------+
|                 | hit                     | hit sequence                     |
|                 +-------------------------+----------------------------------+
|                 | hit_end                 | hit sequence end coordinate, may |
|                 |                         | be 'hmmto' or 'alito' depending  |
|                 |                         | on the HMMER flavor              |
|                 +-------------------------+----------------------------------+
|                 | hit_start               | hit sequence start coordinate,   |
|                 |                         | may be 'hmmfrom' or 'alifrom'    |
|                 |                         | depending on the HMMER flavor    |
|                 +-------------------------+----------------------------------+
|                 | hit_strand              | hit sequence strand              |
|                 +-------------------------+----------------------------------+
|                 | query                   | query sequence                   |
|                 +-------------------------+----------------------------------+
|                 | query_end               | query sequence end coordinate,   |
|                 |                         | may be 'hmmto' or 'alito'        |
|                 |                         | depending on the HMMER flavor    |
|                 +-------------------------+----------------------------------+
|                 | query_start             | query sequence start coordinate, |
|                 |                         | may be 'hmmfrom' or 'alifrom'    |
|                 |                         | depending on the HMMER flavor    |
|                 +-------------------------+----------------------------------+
|                 | query_strand            | query sequence strand            |
+-----------------+-------------------------+----------------------------------+


hmmer3-tab
==========
The following SearchIO objects attributes are provided:

+-----------------+-------------------------+----------------------------------+
| Object          | Attribute               | Column / Value                   |
+=================+=========================+==================================+
| QueryResult     | accession               | query accession (if present)     |
|                 +-------------------------+----------------------------------+
|                 | description             | query sequence description       |
|                 +-------------------------+----------------------------------+
|                 | id                      | query name                       |
+-----------------+-------------------------+----------------------------------+
| Hit             | accession               | hit accession                    |
|                 +-------------------------+----------------------------------+
|                 | bias                    | hit-level bias                   |
|                 +-------------------------+----------------------------------+
|                 | bitscore                | hit-level score                  |
|                 +-------------------------+----------------------------------+
|                 | description             | hit sequence description         |
|                 +-------------------------+----------------------------------+
|                 | cluster_num             | clu column                       |
|                 +-------------------------+----------------------------------+
|                 | domain_exp_num          | exp column                       |
|                 +-------------------------+----------------------------------+
|                 | domain_included_num     | inc column                       |
|                 +-------------------------+----------------------------------+
|                 | domain_obs_num          | dom column                       |
|                 +-------------------------+----------------------------------+
|                 | domain_reported_num     | rep column                       |
|                 +-------------------------+----------------------------------+
|                 | env_num                 | env column                       |
|                 +-------------------------+----------------------------------+
|                 | evalue                  | hit-level evalue                 |
|                 +-------------------------+----------------------------------+
|                 | id                      | target name                      |
|                 +-------------------------+----------------------------------+
|                 | overlap_num             | ov column                        |
|                 +-------------------------+----------------------------------+
|                 | region_num              | reg column                       |
+-----------------+-------------------------+----------------------------------+
| HSP             | bias                    | bias of the best domain          |
|                 +-------------------------+----------------------------------+
|                 | bitscore                | bitscore of the best domain      |
|                 +-------------------------+----------------------------------+
|                 | evalue                  | evalue of the best domain        |
+-----------------+-------------------------+----------------------------------+


hmmer3-domtab
=============
To parse domain table files, you must use the HMMER flavor that produced the
file. So instead of using 'hmmer3-domtab', use either 'hmmsearch3-domtab',
'hmmscan3-domtab', or 'phmmer3-domtab'.

The following SearchIO objects attributes are provided:

+-----------------+-------------------------+----------------------------------+
| Object          | Attribute               | Value                            |
+=================+=========================+==================================+
| QueryResult     | accession               | accession                        |
|                 +-------------------------+----------------------------------+
|                 | description             | query sequence description       |
|                 +-------------------------+----------------------------------+
|                 | id                      | query sequence ID                |
|                 +-------------------------+----------------------------------+
|                 | seq_len                 | full length of query sequence    |
+-----------------+-------------------------+----------------------------------+
| Hit             | accession               | accession                        |
|                 +-------------------------+----------------------------------+
|                 | bias                    | hit-level bias                   |
|                 +-------------------------+----------------------------------+
|                 | bitscore                | hit-level score                  |
|                 +-------------------------+----------------------------------+
|                 | description             | hit sequence description         |
|                 +-------------------------+----------------------------------+
|                 | evalue                  | hit-level e-value                |
|                 +-------------------------+----------------------------------+
|                 | id                      | hit sequence ID                  |
|                 +-------------------------+----------------------------------+
|                 | seq_len                 | length of hit sequence or HMM    |
+-----------------+-------------------------+----------------------------------+
| HSP             | acc_avg                 | expected accuracy per alignment  |
|                 |                         | residue (acc column)             |
|                 +-------------------------+----------------------------------+
|                 | bias                    | hsp-level bias                   |
|                 +-------------------------+----------------------------------+
|                 | bitscore                | hsp-level score                  |
|                 +-------------------------+----------------------------------+
|                 | domain_index            | the domain index set by HMMER    |
|                 +-------------------------+----------------------------------+
|                 | env_end                 | end coordinate of the envelope   |
|                 +-------------------------+----------------------------------+
|                 | env_start               | start coordinate of the envelope |
|                 +-------------------------+----------------------------------+
|                 | evalue                  | hsp-level independent e-value    |
|                 +-------------------------+----------------------------------+
|                 | evalue_cond             | hsp-level conditional e-value    |
+-----------------+-------------------------+----------------------------------+
| HSPFragment     | hit_end                 | hit sequence end coordinate, may |
| (also via HSP)  |                         | be 'hmmto' or 'alito' depending  |
|                 |                         | on the HMMER flavor              |
|                 +-------------------------+----------------------------------+
|                 | hit_start               | hit sequence start coordinate,   |
|                 |                         | may be 'hmmfrom' or 'alifrom'    |
|                 |                         | depending on the HMMER flavor    |
|                 +-------------------------+----------------------------------+
|                 | hit_strand              | hit sequence strand              |
|                 +-------------------------+----------------------------------+
|                 | query_end               | query sequence end coordinate,   |
|                 |                         | may be 'hmmto' or 'alito'        |
|                 |                         | depending on the HMMER flavor    |
|                 +-------------------------+----------------------------------+
|                 | query_start             | query sequence start coordinate, |
|                 |                         | may be 'hmmfrom' or 'alifrom'    |
|                 |                         | depending on the HMMER flavor    |
|                 +-------------------------+----------------------------------+
|                 | query_strand            | query sequence strand            |
+-----------------+-------------------------+----------------------------------+

"""

from .hmmer2_text import Hmmer2TextParser, Hmmer2TextIndexer
from .hmmer3_domtab import (
    Hmmer3DomtabParser,
    Hmmer3DomtabHmmhitParser,
    Hmmer3DomtabHmmqueryParser,
)
from .hmmer3_domtab import Hmmer3DomtabHmmhitIndexer, Hmmer3DomtabHmmqueryIndexer
from .hmmer3_domtab import Hmmer3DomtabHmmhitWriter, Hmmer3DomtabHmmqueryWriter
from .hmmer3_text import Hmmer3TextParser, Hmmer3TextIndexer
from .hmmer3_tab import Hmmer3TabParser, Hmmer3TabIndexer, Hmmer3TabWriter


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
