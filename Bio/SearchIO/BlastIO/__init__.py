# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

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
a bug report if you stumble upon an unparsable file.

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
|                | version                 | BlastOutput_version [*]_    |
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

You may notice that in BLAST XML files, sometimes BLAST replaces your true
sequence ID with its own generated ID. For example, the query IDs become
'Query_1', 'Query_2', and so on. While the hit IDs sometimes become
'gnl|BL_ORD_ID|1', 'gnl|BL_ORD_ID|2', and so on. In these cases, BLAST lumps the
true sequence IDs together with their descriptions.

The blast-xml parser is aware of these modifications and will attempt to extract
the true sequence IDs out of the descriptions. So when accessing QueryResult or
Hit objects, you will use the non-BLAST-generated IDs.

This behavior on the query IDs can be disabled using the 'use_raw_query_ids'
parameter while the behavior on the hit IDs can be disabled using the
'use_raw_hit_ids' parameter. Both are boolean values that can be supplied
to SearchIO.read or SearchIO.parse, with the default values set to 'False'.

In any case, the raw BLAST IDs can always be accessed using the query or hit
object's 'blast_id' attribute.

The blast-xml write function also accepts 'use_raw_query_ids' and
'use_raw_hit_ids' parameters. However, note that the default values for the
writer are set to 'True'. This is because the writer is meant to mimic native
BLAST result as much as possible.


blast-tab
=========

The default format for blast-tab support is the variant without comments (-m 6
flag). Commented BLAST tabular files may be parsed, indexed, or written using
the keyword argument 'comments' set to True:

    >>> # blast-tab defaults to parsing uncommented files
    >>> from Bio import SearchIO
    >>> uncommented = 'Blast/tab_2226_tblastn_004.txt'
    >>> qresult = SearchIO.read(uncommented, 'blast-tab')
    >>> qresult
    QueryResult(id='gi|11464971:4-101', 5 hits)

    >>> # set the keyword argument to parse commented files
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

    >>> # pass the custom column names as a Python list
    >>> fname = 'Blast/tab_2226_tblastn_009.txt'
    >>> custom_fields = ['qseqid', 'sseqid']
    >>> qresult = next(SearchIO.parse(fname, 'blast-tab', fields=custom_fields))
    >>> qresult
    QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)

    >>> # pass the custom column names as a space-separated string
    >>> fname = 'Blast/tab_2226_tblastn_009.txt'
    >>> custom_fields = 'qseqid sseqid'
    >>> qresult = next(SearchIO.parse(fname, 'blast-tab', fields=custom_fields))
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
|             | ident_num         | nident       |
|             +-------------------+--------------+
|             | ident_pct         | pident       |
|             +-------------------+--------------+
|             | mismatch_num      | mismatch     |
|             +-------------------+--------------+
|             | pos_pct           | ppos         |
|             +-------------------+--------------+
|             | pos_num           | positive     |
|             +-------------------+--------------+
|             | bitscore_raw      | score        |
+-------------+-------------------+--------------+
| HSPFragment | frames            | frames [*]_  |
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


.. [*] may be modified

.. [*] When 'frames' is present, both ``query_frame`` and ``hit_frame`` will be
   present as well. It is recommended that you use these instead of 'frames' directly.

"""

from .blast_tab import BlastTabParser, BlastTabIndexer, BlastTabWriter
from .blast_xml import BlastXmlParser, BlastXmlIndexer, BlastXmlWriter

# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
