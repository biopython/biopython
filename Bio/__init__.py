# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

__all__ = [
    "Align",
    "Alphabet",
    "Blast",
    "Clustalw",
    "Data",
    "Encodings",
    "Enzyme",
    "FSSP",
    "Fasta",
    "File",
    "GenBank",
    "Gobase",
    "Index",
    "InterPro",
    "KEGG",
    "Kabat",
    "Medline",
    "PDB",
    "ParserSupport",
    "PropertyManager",
    "Prosite",
    "Rebase",
    "RecordFile",
    "SCOP",
    "Seq",
    "SeqFeature",
    "SeqIO",
    "SeqRecord",
    "SubsMat",
    "SwissProt",
    "Tools",
    "UniGene",
    "WWW",
    "utils"
    ]

import sys
if getattr(sys, "version_info", (1, 5))[:2] >= (2, 1):
  import FormatRegistry
  formats = FormatRegistry.FormatRegistry("Bio")
  register_format = formats.register_format
  link_format = formats.link

  import DBRegistry
  db = DBRegistry.DBRegistry("Bio")
  register_db = db.register_db
  group_db = db.group
  del DBRegistry

  import config
  seqdatabase = config.SeqDatabase()
  del config
  
del sys
