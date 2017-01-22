# -*- coding: utf-8 -*-

# Copyright 2017 by Joanna Zbijewska, Agata Gruszczyńska, Michał Karlicki.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Note that:
# 1. RMSD scripts and license is available separately.
#    We added it in file: calculate_rmsd and README2.
#
# 2. RNAstructure license is available separately.
#    Please consult rna.urmc.rochester.edu .

"""
Author = "Joanna Zbijewska"
"""

import sys
from Bio.RNA_Structure.API.API_RNA_STRAND import *
from Bio.RNA_Structure.API.API_NDB import *


class RNA_API():
    """RNA_API is a class to combine two API classes for downloading RNA structures and sequneces."""

    def __init__(self, database, what, inpt, input_type = None, p = None):
        """Four arguments:
        database: to choose the database you want to access
        what: to choose type of data you want to get - sequence, structure or metadata
        inpt: string to search through the database - a string for sequence or PDB id
        input_type: if the inpt string is a PDB id you have to type "pdb_id"
        p: path to save your desired data, if none give - all files are save in current directory
        """
        self.db = database
        self.what = what
        self.inpt = inpt
        self.type = input_type
        self.p = p
        if p == None:
            self.p = ""

    def choose_db(self):
        """Function to choose the database. It gives the class to inpt."""
        db = self.db
        db = db.upper()
        if self.type is None:
            if db == "NDB":
                c = via_sequence(sequence = self.inpt, path = self.p)
            elif db == "RNA-STRAND" or db == "RNA-STRAND":
                c = RNA_STRAND(sequence = self.inpt, path = self.p)
            else:
                sys.stderr.write("Error: Invalid database name.")
                sys.exit(1)
        elif self.type.lower() == "pdb_id":
            c = Nucleic_acid_database(pdb_id = self.inpt, path = self.p)
        else:
            sys.stderr.write("Error: Invalid input type.")
            sys.exit(1)
        return(c)

    def use_API(self):
        """Downloads sequence, structure or metadata depending on what value"""
        what = self.what
        poss_what = ["sequence", "structure", "metadata"]
        if (what != poss_what[i] for i in xrange(2)):
            sys.stderr.write("Error: Invalid recquired output.")
            sys.exit(1)
        if what == poss_what[0]:
            c.download_fasta_sequence()
        elif what == poss_what[1]:
            if type(c) == "RNA_STRAND":
                c.download_bpseq_structure()
            else:
                c.download_pdb_structure()
        else:
            c.metadata_to_file()

#Przyklad uzycia:
#x = RNA_API("ndb","sequence","5SWE")
#x.use_API()
