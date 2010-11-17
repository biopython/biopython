# Copyright 2009 by David Winter.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import sys
import unittest
import subprocess

from Bio import MissingExternalDependencyError
from Bio import AlignIO
from Bio.Nexus import Trees # One day we should use planned TreeIO module

from Bio.Emboss.Applications import FDNADistCommandline, FNeighborCommandline
from Bio.Emboss.Applications import FSeqBootCommandline, FProtDistCommandline
from Bio.Emboss.Applications import FProtParsCommandline, FConsenseCommandline
from Bio.Emboss.Applications import FTreeDistCommandline, FDNAParsCommandline

#Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

exes_wanted = ['fdnadist', 'fneighbor', 'fprotdist','fprotpars','fconsense',
               'fseqboot', 'ftreedist', 'fdnapars']
exes = dict() #Dictionary mapping from names to exe locations

if "EMBOSS_ROOT" in os.environ:
    #Windows default installation path is C:\mEMBOSS which contains the exes.
    #EMBOSS also sets an environment variable which we will check for.
    path = os.environ["EMBOSS_ROOT"]
    if os.path.isdir(path):
        for name in exes_wanted:
            if os.path.isfile(os.path.join(path, name+".exe")):
                exes[name] = os.path.join(path, name+".exe")
    del path, name
if sys.platform!="win32":
    import commands
    for name in exes_wanted:
        #This will "just work" if installed on the path as normal on Unix
        output = commands.getoutput("%s -help" % name)
        if "not found" not in output and "not recognized" not in output:
            exes[name] = name
        del output
    del name

if len(exes) < len(exes_wanted):
    raise MissingExternalDependencyError(\
          "Install the Emboss package 'PhylipNew' if you want to use the "+\
          "Bio.Emboss.Applications wrappers for phylogenetic tools.")

 ###########################################################################

# A few top level functions that are called repeatedly in the test cases
def write_AlignIO_dna():
    """Convert opuntia.aln to a phylip file"""
    assert 1 == AlignIO.convert("Clustalw/opuntia.aln", "clustal",
                                "Phylip/opuntia.phy", "phylip")

def write_AlignIO_protein():
    """Convert hedgehog.aln to a phylip file"""
    assert 1 == AlignIO.convert("Clustalw/hedgehog.aln", "clustal",
                                "Phylip/hedgehog.phy", "phylip")

def clean_up():
    """Delete tests files (to be used as tearDown() function in test fixtures)"""
    for filename in ["test_file", "Phylip/opuntia.phy","Phylip/hedgehog.phy"]:
        if os.path.isfile(filename):
            os.remove(filename)

def parse_trees(filename):
    """Helper function until we have Bio.Phylo on trunk."""
    data = open("test_file", "r").read()
    for tree_str in data.split(";\n"):
        if tree_str:
            yield Trees.Tree(tree_str+";")

class DistanceTests(unittest.TestCase):
    """Tests for calculating distance based phylogenetic trees with phylip"""

    def tearDown(self):
        clean_up()

    test_taxa = ['Archaeohip', 'Calippus', 'Hypohippus', 'M._secundu',
                 'Merychippu', 'Mesohippus', 'Nannipus', 'Neohippari',
                 'Parahippus', 'Pliohippus']
    
    def distances_from_alignment(self, filename, DNA = True):
        """check we can make distance matrix from a given alignment"""
        self.assertTrue(os.path.isfile(filename), "Missing %s" % filename)
        if DNA:
            cline =  FDNADistCommandline(exes["fdnadist"],
                                         method = 'j',
                                         sequence= filename,
                                         outfile = "test_file",
                                         auto = True)
        else:
            cline = FProtDistCommandline(exes["fprotdist"],
                                         method = 'j',
                                         sequence= filename,
                                         outfile = "test_file",
                                         auto = True)
        stdout, strerr = cline()
        #biopython can't grok distance matrices, so we'll just check it exists
        self.assertTrue(os.path.isfile("test_file"))
    
    def tree_from_distances(self, filename):
        """Check we can estimate a tree from a distance matrix"""
        self.assertTrue(os.path.isfile(filename), "Missing %s" % filename)
        cline = FNeighborCommandline(exes["fneighbor"],
                                     datafile = filename,
                                     outtreefile = "test_file",
                                     auto= True, filter = True)
        stdout, stderr = cline()
        for tree in parse_trees("test_file"):
            tree_taxa = [t.replace(" ", "_") for t in tree.get_taxa()]
            self.assertEqual(self.test_taxa, sorted(tree_taxa))

    def test_distances_from_phylip_DNA(self):
        """Calculate a distance matrix from an phylip alignment"""
        self.distances_from_alignment("Phylip/horses.phy")

    def test_distances_from_AlignIO_DNA(self):
        """Calculate a distance matrix from an alignment written by AlignIO"""
        write_AlignIO_dna()
        self.distances_from_alignment("Phylip/opuntia.phy")

    #def test_distances_from_bootstrapped_phylip_DNA(self):
    #    """Calculate a set of distance matrices from phylip alignments"""
    #    self.distances_from_alignment("Phylip/bs_horses.phy")

    # fprotdist tests
    def test_distances_from_protein_phylip(self):
        """Calculate a distance matrix from phylip protein alignment"""
        self.distances_from_alignment("Phylip/interlaced.phy", DNA=False)

    def test_distances_from_protein_AlignIO(self):
        """Calculate distance matrix from an AlignIO written protein alignment"""
        write_AlignIO_protein()
        self.distances_from_alignment("Phylip/hedgehog.phy", DNA=False)

    #def test_distances_from_bootstrapped_phylip_protein(self):
    #    """Calculate distance matrices from a bootstrapped protein alignment"""
    #    self.distances_from_alignment("Clustalw/bs_interlaced.phy", DNA=False)

    # fneighbor tests
    #def test_tree_from_distances(self):
    #    """Estimate tree from distance matrix and parse it."""
    #    self.tree_from_distances("Phylip/horses.fdnadist")

    # This one won't work because of a bug in EMBOSS 6.0.1
    #def test_tree_from_bootstrapped_distances(self):
    #    """Estimate tree from bootstrapped distance matrix and parse it"""
    #    self.tree_from_distances("Phylip/bs_horses.fdnadist")

class ParsimonyTests(unittest.TestCase):
    """Tests for estimating parsimony based phylogenetic trees with phylip"""

    def tearDown(self):
        clean_up()

    def parsimony_tree(self, filename, format, DNA=True):
        """Estimate a parsimony tree from an alignment"""
        self.assertTrue(os.path.isfile(filename), "Missing %s" % filename)
        if DNA:
            cline = FDNAParsCommandline(exes["fdnapars"],
                                        sequence = filename,
                                        outtreefile = "test_file",
                                        auto= True, stdout=True)
        else:
            cline = FProtParsCommandline(exes["fprotpars"],
                                         sequence = filename,
                                         outtreefile = "test_file",
                                         auto= True, stdout=True)
        stdout, stderr = cline()
        a_taxa = [s.name.replace(" ", "_") for s in
                  AlignIO.parse(open(filename, "r"), format).next()]
        for tree in parse_trees("test_file"):
            t_taxa = [t.replace(" ", "_") for t in tree.get_taxa()]
            self.assertEqual(sorted(a_taxa), sorted(t_taxa))
    
    # fdnapars tests
    #def test_parsimony_tree_from_phylip_DNA(self):
    #    """Make a parsimony tree from a phylip DNA alignment"""
    #    self.parsimony_tree("Phylip/horses.phy", "phylip")

    def test_parsimony_tree_from_AlignIO_DNA(self):
        """Make a parsimony tree from an alignment written with AlignIO"""
        write_AlignIO_dna()
        self.parsimony_tree("Phylip/opuntia.phy", "phylip")

    #def test_parsimony_bootstrapped_phylip_DNA(self):
    #    """Make a parsimony tree from a bootstrapped phylip DNA alignment"""
    #    self.parsimony_tree("Phylip/bs_horses.phy", "phylip")

    # fprotpars tests
    #def test_parsimony_tree_from_phylip_protein(self):
    #    """Make a parsimony tree from a phylip DNA alignment"""
    #    self.parsimony_tree("Phylip/interlaced.phy", "phylip", DNA=False)

    def test_parsimony_from_AlignIO_protein(self):
        """Make a parsimony tree from protein alignment written with AlignIO"""
        write_AlignIO_protein()
        self.parsimony_tree("Phylip/interlaced.phy", "phylip", DNA=False)

    #def test_parsimony_tree_bootstrapped_phylip_protein(self):
    #    """Make a parsimony tree from a phylip DNA alignment"""
    #    self.parsimony_tree("Phylip/bs_interlaced.phy", "phylip", DNA=False)

class BootstrapTests(unittest.TestCase):
    """Tests for pseudosampling alignments with fseqboot"""

    def tearDown(self):
        clean_up()
     
    def check_bootstrap(self, filename, format, align_type="d"):
        """Check we can use fseqboot to pseudosample an alignment
        
        The align_type type argument is passed to the commandline object to
        set the output format to use (from [D]na,[p]rotein and [r]na )
        """
        self.assertTrue(os.path.isfile(filename), "Missing %s" % filename)
        cline = FSeqBootCommandline(exes["fseqboot"],
                                    sequence = filename,
                                    outfile =  "test_file",
                                    seqtype = align_type,
                                    reps = 2,
                                    auto = True, filter = True)
        stdout, stderr = cline()
        # the resultant file should have 2 alignments...
        bs = list(AlignIO.parse(open("test_file", "r" ), format))
        self.assertEqual(len(bs), 2)
        # ..and each name in the original alignment...
        a_names = [s.name.replace(" ", "_") for s in
                   AlignIO.read(open(filename, "r"), format)]
        # ...should be in each alignment in the bootstrapped file
        for a in bs:
            self.assertEqual(a_names, [s.name.replace(" ", "_") for s in a])

    def test_bootstrap_phylip_DNA(self):
        """Pseudosample a phylip DNA alignment"""
        self.check_bootstrap("Phylip/horses.phy", "phylip")

    def test_bootstrap_AlignIO_DNA(self):
        """Pseudosample a phylip DNA alignment written with AlignIO"""
        write_AlignIO_dna()
        self.check_bootstrap("Phylip/opuntia.phy", "phylip")

    def test_bootstrap_phylip_protein(self):
        """Pseudosample a phylip protein alignment"""
        self.check_bootstrap("Phylip/interlaced.phy", "phylip", "p")

    def test_bootstrap_AlignIO_protein(self):
        """Pseudosample a phylip protein alignment written with AlignIO"""
        write_AlignIO_protein()
        self.check_bootstrap("Phylip/hedgehog.phy", "phylip", "p")

class TreeComparisonTests(unittest.TestCase):
    """Tests for comparing phylogenetic trees with phylip tools"""

    def tearDown(self):
        clean_up()

    def test_fconsense(self):
        """Calculate a consensus tree with fconsense"""
        cline = FConsenseCommandline(exes["fconsense"],
                                     intreefile = "Phylip/horses.tree",
                                     outtreefile = "test_file",
                                     auto = True, filter = True)
        stdout, stderr = cline()
        #Split the next and get_taxa into two steps to help 2to3 work
        tree1 = parse_trees("test_file").next()
        taxa1 = tree1.get_taxa()
        for tree in parse_trees("Phylip/horses.tree"):
            taxa2 = tree.get_taxa()
            self.assertEqual(sorted(taxa1),sorted(taxa2))

    def test_ftreedist(self):
        """Calculate the distance between trees with ftreedist"""
        cline = FTreeDistCommandline(exes["ftreedist"],
                                     intreefile = "Phylip/horses.tree",
                                     outfile = "test_file",
                                     auto = True, filter = True)
        stdout, stderr = cline()
        self.assertTrue(os.path.isfile("test_file"))

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
    clean_up()
