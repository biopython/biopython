#!/usr/bin/env python
#
# Copyright 2013 by Leighton Pritchard.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

""" Tests for online functionality of the KGML online component
"""

# Builtins
import os
import unittest

def check_internet():
    try:
        check.available
    except AttributeError:
        # I'm going to check for internet availability
        RELIABLE_DOMAIN = "biopython.org"
        import socket
        try:
            socket.getaddrinfo(RELIABLE_DOMAIN,
                               80,
                               socket.AF_UNSPEC,
                               socket.SOCK_STREAM)
        except socket.gaierror, x:
            check.available = False
        else:
            check.available = True
    if not check.available:
        raise MissingExternalDependencyError("internet not available")

# Biopython Bio.KEGG.KGML
from Bio.KEGG.KGML.KGML_scrape import *

class KGMLPathwayTest(unittest.TestCase):
    """ Import the ko01100 metabolic map from a local .xml KGML file, and from
        the KEGG site, and write valid KGML output for each
    """
    def setUp(self):
        # Does our output directory exist?  If not, create it
        if not os.path.isdir('KEGG'):
            os.mkdir('KEGG')

    def test_KEGG_download_and_write(self):
        """ Download a KEGG pathway from the KEGG server and write KGML.
        """
        # Download the KEGG ko01120 pathway and write to file as KGML
        retrieve_kgml_to_file("ko01120", os.path.join("KEGG", "ko01120.xml"))
        
    def test_KEGG_download_to_pathway(self):
        """ Download a KEGG pathway from the KEGG server and write KGML.
        """
        # Download the KEGG ko03070 pathway as a Pathway object
        retrieve_KEGG_pathway("ko03070")
        
    def test_KEGG_download_handle(self):
        """ Download a KEGG pathway from the KEGG server and write KGML.
        """
        # Download the KEGG ko03070 pathway as a filehandle
        retrieve_kgml_stream("ko03070")
        


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner = runner)
