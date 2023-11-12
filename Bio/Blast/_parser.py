# Copyright 2023 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code to parse the BLAST XML DTD file.

The BLAST XML DTD file is available on the NCBI site at:
https://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd
"""

import os.path
from xml.parsers import expat

from Bio import Entrez


class DTDHandler:
    """Parser for the BLAST XML DTD file."""

    def __init__(self):
        """Initialize the parser and parse the BLAST XML DTD file."""
        parser = expat.ParserCreate()
        parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_ALWAYS)
        parser.ExternalEntityRefHandler = self._externalEntityRefHandler
        self.parser = parser
        self.names = []
        self._externalEntityRefHandler(None, None, "NCBI_BlastOutput.dtd", None)

    def _elementDeclHandler(self, name, model):
        self.names.append(name)

    def _externalEntityRefHandler(self, context, base, systemId, publicId):
        assert context is None
        assert base is None
        directory = Entrez.__path__[0]
        path = os.path.join(directory, "DTDs", systemId)
        parser = self.parser.ExternalEntityParserCreate(None)
        parser.ElementDeclHandler = self._elementDeclHandler
        with open(path, "rb") as stream:
            parser.ParseFile(stream)
        return 1
