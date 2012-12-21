# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO support for BLAST+ output formats.

This module adds support for parsing BLAST+ outputs. BLAST+ is a rewrite of
NCBI's legacy BLAST (Basic Local Alignment Search Tool), based on the NCBI
C++ toolkit. The BLAST+ suite is available as command line programs or on
NCBI's web page.

Bio.SearchIO.BlastIO was tested on the following BLAST+ flavors and versions:

    - flavors: blastn, blastp, blastx, tblastn, tblastx
    - versions: 2.2.22+, 2.2.26+

You should also be able to parse outputs from a local BLAST+ search or from
NCBI's web interface. Although the module was not tested against all BLAST+,
it should still be able to parse these other versions' outputs. Please submit
a bug report if you stumble upon an unparseable file.

Some output formats from the BLAST legacy suite (BLAST+'s predecessor) may
still be parsed by this module. However, results are not guaranteed. You may
try to use the Bio.Blast module to parse them instead.

More information about BLAST are available through these links:
  - Publication: http://www.biomedcentral.com/1471-2105/10/421
  - Web interface: http://blast.ncbi.nlm.nih.gov/
  - User guide: http://www.ncbi.nlm.nih.gov/books/NBK1762/


Supported Formats
=================

Bio.SearchIO.BlastIO supports the following BLAST+ output formats:

  - XML        - 'blast-xml'  - parsing, indexing, writing
  - Tabular    - 'blast-tab'  - parsing, indexing, writing
  - Plain text - 'blast-text' - parsing


blast-xml
=========

The blast-xml parser follows the BLAST XML DTD written here:
http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.mod.dtd

It provides the following attributes for each SearchIO object:

+----------------+-------------------------+-----------------------------+
| Object         | Attribute               | XML Element                 |
+================+=========================+=============================+
| QueryResult    | target                  | BlastOutput_db              |
|                +-------------------------+-----------------------------+
|                | program                 | BlastOutput_program         |
|                +-------------------------+-----------------------------+
|                | reference               | BlastOutput_reference       |
|                +-------------------------+-----------------------------+
|                | version                 | BlastOutput_version*        |
|                +-------------------------+-----------------------------+
|                | description             | Iteration_query-def         |
|                +-------------------------+-----------------------------+
|                | id                      | Iteration_query-ID          |
|                +-------------------------+-----------------------------+
|                | seq_len                 | Iteration_query-len         |
|                +-------------------------+-----------------------------+
|                | param_evalue_threshold  | Parameters_expect           |
|                +-------------------------+-----------------------------+
|                | param_entrez_query      | Parameters_entrez-query     |
|                +-------------------------+-----------------------------+
|                | param_filter            | Parameters_filter           |
|                +-------------------------+-----------------------------+
|                | param_gap_extend        | Parameters_gap-extend       |
|                +-------------------------+-----------------------------+
|                | param_gap_open          | Parameters_gap-open         |
|                +-------------------------+-----------------------------+
|                | param_include           | Parameters_include          |
|                +-------------------------+-----------------------------+
|                | param_matrix            | Parameters_matrix           |
|                +-------------------------+-----------------------------+
|                | param_pattern           | Parameters_pattern          |
|                +-------------------------+-----------------------------+
|                | param_score_match       | Parameters_sc-match         |
|                +-------------------------+-----------------------------+
|                | param_score_mismatch    | Parameters_sc-mismatch      |
|                +-------------------------+-----------------------------+
|                | stat_db_num             | Statistics_db-num           |
|                +-------------------------+-----------------------------+
|                | stat_db_len             | Statistics_db-len           |
|                +-------------------------+-----------------------------+
|                | stat_eff_space          | Statistics_eff-space        |
|                +-------------------------+-----------------------------+
|                | stat_entropy            | Statistics_entropy          |
|                +-------------------------+-----------------------------+
|                | stat_hsp_len            | Statistics_hsp-len          |
|                +-------------------------+-----------------------------+
|                | stat_kappa              | Statistics_kappa            |
|                +-------------------------+-----------------------------+
|                | stat_lambda             | Statistics_lambda           |
+----------------+-------------------------+-----------------------------+
| Hit            | accession               | Hit_accession               |
|                +-------------------------+-----------------------------+
|                | description             | Hit_def                     |
|                +-------------------------+-----------------------------+
|                | id                      | Hit_id                      |
|                +-------------------------+-----------------------------+
|                | seq_len                 | Hit_len                     |
+----------------+-------------------------+-----------------------------+
| HSP            | bitscore                | Hsp_bit-score               |
|                +-------------------------+-----------------------------+
|                | density                 | Hsp_density                 |
|                +-------------------------+-----------------------------+
|                | evalue                  | Hsp_evalue                  |
|                +-------------------------+-----------------------------+
|                | gap_num                 | Hsp_gaps                    |
|                +-------------------------+-----------------------------+
|                | ident_num               | Hsp_identity                |
|                +-------------------------+-----------------------------+
|                | pos_num                 | Hsp_positive                |
|                +-------------------------+-----------------------------+
|                | bitscore_raw            | Hsp_score                   |
+----------------+-------------------------+-----------------------------+
| HSPFragment    | aln_span                | Hsp_align-len               |
| (also via      +-------------------------+-----------------------------+
| HSP)           | hit_frame               | Hsp_hit-frame               |
|                +-------------------------+-----------------------------+
|                | hit_start               | Hsp_hit-from                |
|                +-------------------------+-----------------------------+
|                | hit_end                 | Hsp_hit-to                  |
|                +-------------------------+-----------------------------+
|                | hit                     | Hsp_hseq                    |
|                +-------------------------+-----------------------------+
|                | aln_annotation          | Hsp_midline                 |
|                +-------------------------+-----------------------------+
|                | pattern_start           | Hsp_pattern-from            |
|                +-------------------------+-----------------------------+
|                | pattern_end             | Hsp_pattern-to              |
|                +-------------------------+-----------------------------+
|                | query_frame             | Hsp_query-frame             |
|                +-------------------------+-----------------------------+
|                | query_start             | Hsp_query-from              |
|                +-------------------------+-----------------------------+
|                | query_end               | Hsp_query-to                |
|                +-------------------------+-----------------------------+
|                | query                   | Hsp_qseq                    |
+----------------+-------------------------+-----------------------------+
* may be modified

You may notice that in BLAST XML files, sometimes BLAST replaces your true
sequence ID with its own generated ID. For example, the query IDs become
'Query_1', 'Query_2', and so on. While the hit IDs sometimes become
'gnl|BL_ORD_ID|1', 'gnl|BL_ORD_ID|2', and so on. In these cases, BLAST lumps the
true sequence IDs together with their descriptions.

The blast-xml parser is aware of these modifications and will attempt to extract
the true sequence IDs out of the descriptions. So when accessing QueryResult or
Hit objects, you will use the non-BLAST-generated IDs.

Conversely, the blast-xml writer will try to concatenate the true sequence IDs
with their descriptions and use the BLAST-generated IDs. This enables you to
write BLAST XML files using SearchIO as if they were written by a real BLAST
program.


blast-tab
=========

The default format for blast-tab support is the variant without comments (-m 6
flag). Commented BLAST tabular files may be parsed, indexed, or written using
the keyword argument 'comments' set to True:

    # blast-tab defaults to parsing uncommented files
    >>> from Bio import SearchIO
    >>> uncommented = 'Blast/tab_2226_tblastn_004.txt'
    >>> qresult = SearchIO.read(uncommented, 'blast-tab')
    >>> qresult
    QueryResult(id='gi|11464971:4-101', 5 hits)

    # set the keyword argument to parse commented files
    >>> commented = 'Blast/tab_2226_tblastn_008.txt'
    >>> qresult = SearchIO.read(commented, 'blast-tab', comments=True)
    >>> qresult
    QueryResult(id='gi|11464971:4-101', 5 hits)

For uncommented files, the parser defaults to using BLAST's default column
ordering: 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send
evalue bitscore'.

If you want to parse an uncommented file with a customized column order, you can
use the 'fields' keyword argument to pass the custom column order. The names of
the column follow BLAST's naming. For example, 'qseqid' is the column for the
query sequence ID. These names may be passed either as a Python list or as a
space-separated strings.

    # pass the custom column names as a Python list
    >>> fname = 'Blast/tab_2226_tblastn_009.txt'
    >>> custom_fields = ['qseqid', 'sseqid']
    >>> qresult = SearchIO.parse(fname, 'blast-tab', fields=custom_fields).next()
    >>> qresult
    QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)

    # pass the custom column names as a space-separated string
    >>> fname = 'Blast/tab_2226_tblastn_009.txt'
    >>> custom_fields = 'qseqid sseqid'
    >>> qresult = SearchIO.parse(fname, 'blast-tab', fields=custom_fields).next()
    >>> qresult
    QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)

You may also use the 'std' field name as an alias to BLAST's default 12 columns,
just like when you run a command line BLAST search.

Note that the 'fields' keyword argument will be ignored if the parsed file is
commented. Commented files have their column ordering stated explicitly in the
file, so there is no need to specify it again in SearchIO.

'comments' and 'fields' keyword arguments are both applicable for parsing,
indexing, and writing.

blast-tab provides the following attributes for each SearchIO objects:

+-------------+-------------------+--------------+
| Object      | Attribute         | Column name  |
+=============+===================+==============+
| QueryResult | accession         | qacc         |
|             +-------------------+--------------+
|             | accession_version | qaccver      |
|             +-------------------+--------------+
|             | gi                | qgi          |
|             +-------------------+--------------+
|             | seq_len           | qlen         |
|             +-------------------+--------------+
|             | id                | qseqid       |
+-------------+-------------------+--------------+
| Hit         | accession         | sacc         |
|             +-------------------+--------------+
|             | accession_version | sacc_ver     |
|             +-------------------+--------------+
|             | gi                | sgi          |
|             +-------------------+--------------+
|             | gi_all            | sallgi       |
|             +-------------------+--------------+
|             | id_all            | sallseqid    |
|             +-------------------+--------------+
|             | seq_len           | slen         |
|             +-------------------+--------------+
|             | id                | sseqid       |
+-------------+-------------------+--------------+
| HSP         | bitscore          | bitscore     |
|             +-------------------+--------------+
|             | btop              | btop         |
|             +-------------------+--------------+
|             | evalue            | evalue       |
|             +-------------------+--------------+
|             | gapopen_num       | gapopen      |
|             +-------------------+--------------+
|             | gap_num           | gaps         |
|             +-------------------+--------------+
|             | ident_pct         | nident       |
|             +-------------------+--------------+
|             | ident_num         | pident       |
|             +-------------------+--------------+
|             | mismatch_num      | mismatch     |
|             +-------------------+--------------+
|             | pos_pct           | ppos         |
|             +-------------------+--------------+
|             | pos_num           | positive     |
|             +-------------------+--------------+
|             | bitscore_raw      | score        |
+-------------+-------------------+--------------+
| HSPFragment | frames            | frames*      |
| (also via   +-------------------+--------------+
| HSP)        | aln_span          | length       |
|             +-------------------+--------------+
|             | query_end         | qend         |
|             +-------------------+--------------+
|             | query_frame       | qframe       |
|             +-------------------+--------------+
|             | query             | qseq         |
|             +-------------------+--------------+
|             | query_start       | qstart       |
|             +-------------------+--------------+
|             | hit_end           | send         |
|             +-------------------+--------------+
|             | hit_frame         | sframe       |
|             +-------------------+--------------+
|             | hit               | sseq         |
|             +-------------------+--------------+
|             | hit_start         | sstart       |
+-------------+-------------------+--------------+
* When 'frames' is present, both `query_frame` and `hit_frame` will be present
  as well. It is recommended that you use these instead of 'frames' directly.

If the parsed file is commented, the following attributes may be available as
well:

+--------------+---------------+----------------------------+
| Object       | Attribute     | Value                      |
+==============+===============+============================+
| QueryResult  | description   | query description          |
|              +---------------+----------------------------+
|              | fields        | columns in the output file |
|              +---------------+----------------------------+
|              | program       | BLAST flavor               |
|              +---------------+----------------------------+
|              | rid           | remote search ID           |
|              +---------------+----------------------------+
|              | target        | target database            |
|              +---------------+----------------------------+
|              | version       | BLAST version              |
+--------------+---------------+----------------------------+


blast-text
==========
The BLAST plain text output format has been known to change considerably between
BLAST versions. NCBI itself has recommended that users not rely on the plain
text output for parsing-related work.

However, in some cases parsing the plain text output may still be useful.
SearchIO provides parsing support for the plain text output, but guarantees only
a minimum level of support. Writing a parser that fully supports plain text
output for all BLAST versions is not a priority at the moment.

If you do have a BLAST plain text file that can not be parsed and would like to
submit a patch, we are more than happy to accept it.

The blast-text parser provides the following object attributes:

+-----------------+-------------------------+----------------------------------+
| Object          | Attribute               | Value                            |
+=================+=========================+==================================+
| QueryResult     | description             | query sequence description       |
|                 +-------------------------+----------------------------------+
|                 | id                      | query sequence ID                |
|                 +-------------------------+----------------------------------+
|                 | program                 | BLAST flavor                     |
|                 +-------------------------+----------------------------------+
|                 | seq_len                 | full length of query sequence    |
|                 +-------------------------+----------------------------------+
|                 | target                  | target database of the search    |
|                 +-------------------------+----------------------------------+
|                 | version                 | BLAST version                    |
+-----------------+-------------------------+----------------------------------+
| Hit             | evalue                  | hit-level evalue, from the hit   |
|                 |                         | table                            |
|                 +-------------------------+----------------------------------+
|                 | id                      | hit sequence ID                  |
|                 +-------------------------+----------------------------------+
|                 | description             | hit sequence description         |
|                 +-------------------------+----------------------------------+
|                 | score                   | hit-level score, from the hit    |
|                 |                         | table                            |
|                 +-------------------------+----------------------------------+
|                 | seq_len                 | full length of hit sequence      |
+-----------------+-------------------------+----------------------------------+
| HSP             | evalue                  | hsp-level evalue                 |
|                 +-------------------------+----------------------------------+
|                 | bitscore                | hsp-level bit score              |
|                 +-------------------------+----------------------------------+
|                 | bitscore_raw            | hsp-level score                  |
|                 +-------------------------+----------------------------------+
|                 | gap_num                 | number of gaps in alignment      |
|                 +-------------------------+----------------------------------+
|                 | ident_num               | number of identical residues     |
|                 |                         | in alignment                     |
|                 +-------------------------+----------------------------------+
|                 | pos_num                 | number of positive matches in    |
|                 |                         | alignment                        |
+-----------------+-------------------------+----------------------------------+
| HSPFragment     | aln_annotation          | alignment homology string        |
| (also via       +-------------------------+----------------------------------+
| HSP)            | aln_span                | length of alignment fragment     |
|                 +-------------------------+----------------------------------+
|                 | hit                     | hit sequence                     |
|                 +-------------------------+----------------------------------+
|                 | hit_end                 | hit sequence end coordinate      |
|                 +-------------------------+----------------------------------+
|                 | hit_frame               | hit sequence reading frame       |
|                 +-------------------------+----------------------------------+
|                 | hit_start               | hit sequence start coordinate    |
|                 +-------------------------+----------------------------------+
|                 | hit_strand              | hit sequence strand              |
|                 +-------------------------+----------------------------------+
|                 | query                   | query sequence                   |
|                 +-------------------------+----------------------------------+
|                 | query_end               | query sequence end coordinate    |
|                 +-------------------------+----------------------------------+
|                 | query_frame             | query sequence reading frame     |
|                 +-------------------------+----------------------------------+
|                 | query_start             | query sequence start coordinate  |
|                 +-------------------------+----------------------------------+
|                 | query_strand            | query sequence strand            |
+-----------------+-------------------------+----------------------------------+

"""

from blast_tab import *
from blast_xml import *
from blast_text import *


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
