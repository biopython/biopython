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
    "UniGene",
    "WWW",
    "utils"
    ]

def _load_registries():
    import sys, os
    from Bio.config.Registry import Registry
    
    if getattr(sys, "version_info", (1, 5))[:2] < (2, 1):
        return

    self = sys.modules[__name__]        # self refers to this module.
    # Load the registries.  Look in all the '.py' files in Bio.config
    # for Registry objects.  Save them all into the local namespace.
    x = os.listdir(
        os.path.dirname(__import__("Bio.config", {}, {}, ["Bio"]).__file__))
    x = filter(lambda x: not x.startswith("_") and x.endswith(".py"), x)
    x = map(lambda x: x[:-3], x)            # chop off '.py'
    for module in x:
        module = __import__("Bio.config.%s" % module, {}, {}, ["Bio","config"])
        for name, obj in module.__dict__.items():
            if name.startswith("_") or not isinstance(obj, Registry):
                continue
            setattr(self, name, obj)

# Put the registry loading code in a function so we don't polute the
# module namespace with local variables.
_load_registries()
del _load_registries
