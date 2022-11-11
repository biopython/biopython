# Copyright 2018 by Adhemar Zerlotini. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SearchIO support for InterProScan output formats.

This module adds support for parsing InterProScan XML output.
The InterProScan is available as a command line program or on
EMBL-EBI's web page.
Bio.SearchIO.InterproscanIO was tested on the following version:

- versions: 5.26-65.0 (interproscan-model-2.1.xsd)

More information about InterProScan are available through these links:
- Publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3998142/
- Web interface: https://www.ebi.ac.uk/interpro/search/sequence-search
- Documentation: https://github.com/ebi-pf-team/interproscan/wiki


Supported format
================

Bio.SearchIO.InterproscanIO supports the following format:

- XML   - 'interproscan-xml' - parsing


interproscan-xml
================

The interproscan-xml parser follows the InterProScan XML described here:
https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats

+--------------+--------------------+------------------------------------------+
| Object       | Attribute          | XML Element                              |
+==============+====================+==========================================+
| QueryResult  | target             | ``InterPro``                             |
|              +--------------------+------------------------------------------+
|              | program            | ``InterProScan``                         |
|              +--------------------+------------------------------------------+
|              | version            | ``protein-matches.interproscan-version`` |
+--------------+--------------------+------------------------------------------+
| Hit          | accession          | ``signature.name``                       |
|              +--------------------+------------------------------------------+
|              | id                 | ``signature.ac``                         |
|              +--------------------+------------------------------------------+
|              | description        | ``signature.desc``                       |
|              +--------------------+------------------------------------------+
|              | dbxrefs            | ``IPR:entry.ac``                         |
|              |                    | ``go-xref.id``                           |
|              |                    | ``pathway-xref.db:pathway-xref.id``      |
|              +--------------------+------------------------------------------+
|              | attributes         |                                          |
|              | ['Target']         | ``*-match`` / ``*-location``             |
|              | ['Target version'] | ``signature-library-release.library``    |
|              | ['Hit type']       | ``signature-library-release.version``    |
+--------------+--------------------+------------------------------------------+
| HSP          | bitscore           | ``*-location.score``                     |
|              +--------------------+------------------------------------------+
|              | evalue             | ``*-location.evalue``                    |
+--------------+--------------------+------------------------------------------+
| HSPFragment  | query_start        | ``*-location.start``                     |
| (also via    +--------------------+------------------------------------------+
| HSP)         | query_end          | ``*-location.end``                       |
|              +--------------------+------------------------------------------+
|              | hit_start          | ``*-location.hmm-start``                 |
|              +--------------------+------------------------------------------+
|              | hit_end            | ``*-location.hmm-end``                   |
|              +--------------------+------------------------------------------+
|              | query              | ``sequence``                             |
+--------------+--------------------+------------------------------------------+

InterProScan XML files may contain a match with multiple locations or multiple
matches to the same protein with a single location. In both cases, the match
is uniquely stored as a HIT object and the locations as HSP objects.

``HSP.*start == *start - 1`` (Since every start position is 0-based in Biopython)

``HSP.aln_span ==  query-end - query-start``

The types of matches or locations (eg. hmmer3-match, hmmer3-location,
coils-match, panther-location) are stored in hit.attributes['Hit type'].
For instance, for every 'phobious-match', there will be a 'phobious-location'.
Therefore, Hit.type will store the string excluding '-match' or '-location'
('phobious', in this example).
"""

from .interproscan_xml import InterproscanXmlParser

# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
