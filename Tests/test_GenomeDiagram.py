#!/usr/bin/env python
"""Tests for GenomeDiagram general functionality.
"""

##########
# IMPORTS

# Builtins
import os
import sys
import unittest
import string

# Do we have ReportLab?  Raise error if not present.
try:
    from reportlab.lib import colors
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
except:
    raise MissingExternalDependencyError(\
            "Install reportlab if you want to use Bio.Graphics.")

# Biopython core
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio import SeqUtils

# Bio.Graphics.GenomeDiagram
from Bio.Graphics.GenomeDiagram.FeatureSet import FeatureSet
from Bio.Graphics.GenomeDiagram.GraphSet import GraphSet
from Bio.Graphics.GenomeDiagram.Track import Track
#from Bio.Graphics.GenomeDiagram.Utilities import *
from Bio.Graphics.GenomeDiagram import Diagram
from Bio.Graphics.GenomeDiagram.Colors import ColorTranslator
from Bio.Graphics.GenomeDiagram.Graph import GraphData

###############################################################################
# Utility functions for graph plotting, originally in GenomeDiagram.Utilities #
# See Bug 2705 for discussion on where to put these functions in Biopython... #
###############################################################################
def apply_to_window(sequence, window_size, function, step=None):
    """ apply_to_window(sequence, window_size, function) -> [(int, float),(int, float),...]

        o sequence      Bio.Seq.Seq object

        o window_size   Int describing the length of sequence to consider

        o step          Int describing the step to take between windows
                        (default = window_size/2)

        o function      Method or function that accepts a Bio.Seq.Seq object
                        as its sole argument and returns a single value

        Returns a list of (position, value) tuples for fragments of the passed
        sequence of length window_size (stepped by step), calculated by the
        passed function.  Returned positions are the midpoint of each window.
    """
    seqlen = len(sequence)      # Total length of sequence to be used
    if step is None:    # No step specified, so use half window-width or 1 if larger
        step = max(window_size/2, 1)
    else:               # Use specified step, or 1 if greater
        step = max(step, 1)

    results = []    # Holds (position, value) results

    # Perform the passed function on as many windows as possible, short of
    # overrunning the sequence
    pos = 0
    while pos < seqlen-window_size+1:
        # Obtain sequence fragment
        start, middle, end = pos, (pos+window_size+pos)/2, pos+window_size
        fragment = sequence[start:end]
        # Apply function to the sequence fragment
        value = function(fragment)
        results.append((middle, value)) # Add results to list
        # Advance to next fragment
        pos += step

    # Use the last available window on the sequence, even if it means
    # re-covering old ground
    if pos != seqlen - window_size:
        # Obtain sequence fragment
        pos = seqlen - window_size
        start, middle, end = pos, (pos+window_size+pos)/2, pos+window_size
        fragment = sequence[start:end]
        # Apply function to sequence fragment
        value = function(fragment)
        results.append((middle, value)) # Add results to list
        
    # Check on last sequence
    #print fragment
    #print seq[-100:]
    return results      # Return the list of (position, value) results

def calc_gc_content(sequence):
    """ calc_gc_content(sequence)

        o sequence  A Bio.Seq.Seq object

        Returns the % G+C content in a passed sequence
    """
    d = {}
    for nt in ['A','T','G','C']:
        d[nt] = sequence.count(nt) + sequence.count(nt.lower())
    gc = d.get('G',0) + d.get('C',0)

    if gc == 0: return 0
    #print gc*100.0/(d['A'] +d['T'] + gc)
    return gc*1./(d['A'] +d['T'] + gc)


def calc_at_content(sequence):
    """ calc_at_content(sequence)

        o sequence  A Bio.Seq.Seq object

        Returns the % A+T content in a passed sequence
    """
    seq = sequence.data
    d = {}
    for nt in ['A','T','G','C']:
        d[nt] = sequence.count(nt) + sequence.count(nt.lower())
    at = d.get('A',0) + d.get('T',0)

    if at == 0: return 0
    return at*1./(d['G'] +d['G'] + at)


def calc_gc_skew(sequence):
    """ calc_gc_skew(sequence)

        o sequence   A Bio.Seq.Seq object

        Returns the (G-C)/(G+C) GC skew in a passed sequence
    """
    g = sequence.count('G') + sequence.count('g')
    c = sequence.count('C') + sequence.count('c')
    if g+c == 0:
        return 0
    else:
        return (g-c)/float(g+c)


def calc_at_skew(sequence):
    """ calc_at_skew(sequence)

        o sequence   A Bio.Seq.Seq object

        Returns the (A-T)/(A+T) AT skew in a passed sequence
    """
    a = sequence.count('A')
    t = sequence.count('T')
    if a+t == 0:
        return 0
    else:
        return (a-t)/float(a+t)

###############################################################################
# End of utility functions for graph plotting                                 #
###############################################################################

# Tests
def run_tests(argv):
    test_suite = testing_suite()

    runner = unittest.TextTestRunner(sys.stdout, verbosity=2)
    runner.run(test_suite)

def testing_suite():
    test_suite = unittest.TestSuite()

    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'
    tests = [GraphTest, ColorsTest, TrackTest, DiagramTest]   

    for test in tests:
        current_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(current_suite)

    return test_suite        

class TrackTest(unittest.TestCase):
    # TODO Bring code from Track.py, unsure about what test does
    pass

class ColorsTest(unittest.TestCase):
    def t_color_conversions(self):
        """Test color translations.
        """
        translator = ColorTranslator()
        
        # Does the translate method correctly convert the passed argument?
        assert translator.float1_color((0.5, 0.5, 0.5)) == translator.translate((0.5, 0.5, 0.5)), \
            "Did not correctly translate colour from floating point RGB tuple"
        assert translator.int255_color((1, 75, 240)) == translator.translate((1, 75, 240)), \
            "Did not correctly translate colour from integer RGB tuple"
        assert translator.artemis_color(7) == translator.translate(7), \
            "Did not correctly translate colour from Artemis colour scheme"                        
        assert translator.scheme_color(2) == translator.translate(2), \
            "Did not correctly translate colour from user-defined colour scheme"

            
class GraphTest(unittest.TestCase):
    def setUp(self):
        self.data = [(1, 10), (5, 15), (20, 40)]
        
    def t_slicing(self):
        gd = GraphData()
        gd.set_data(self.data)
        gd.add_point((10, 20))
        
        assert gd[4:16] == [(5, 15), (10, 20)], \
                "Unable to insert and retrieve points correctly"


class DiagramTest(unittest.TestCase):
    def setUp(self):
        # This has to be done on Windows otherwise it fails with 
        # the following stack trace:
        # TODO Put stack trace
        pdfmetrics.registerFont(TTFont('Vera', 'Vera.ttf'))
    
    def load_sequence(self, filename):
        """Load a GenBank file as a SeqRecord."""
        handle = open(filename, 'r')
        genbank_entry = SeqIO.read(handle, "genbank")
        handle.close()
        return genbank_entry

    def t_diagram_pdf(self):
        """Output circular and linear diagrams of a GenBank sequence to PDF.
        """
        genbank_entry = self.load_sequence(os.path.join("GenBank","arab1.gb"))
        
        gdd = Diagram('Test Diagram')

        gdfs1 = FeatureSet(name='CDS features')
        gdfs2 = FeatureSet(name='gene features')
        gdfs3 = FeatureSet(name='misc_features')
        gdfs4 = FeatureSet(name='repeat regions')

        for feature in genbank_entry.features:
            if feature.type == 'source':
                start = str(feature.location._start)  # Feature start
                end = str(feature.location._end)      # Feature end

                # Remove extraneous leading chars
                while start[0] not in string.digits:  
                    start = start[1:]

                while end[0] not in string.digits:
                    end = end[1:]

                start, end = int(start), int(end)
                if end < start:
                    start, end = end, start
    	    
            gdd.start = start
            gdd.end = end

            if feature.type == 'CDS':
                gdfs1.add_feature(feature) #, name="Some feature or other")

            if feature.type == 'gene':
                gdfs2.add_feature(feature)

            if feature.type == 'misc_feature':
                gdfs3.add_feature(feature, color=colors.orange)

            if feature.type == 'repeat_region':
                gdfs4.add_feature(feature, color=colors.purple)


        gdfs1.set_all_features('label', 1)
        gdfs2.set_all_features('label', 1)
        gdfs3.set_all_features('label', 1)
        gdfs4.set_all_features('label', 1)

        gdfs3.set_all_features('hide', 0)
        gdfs4.set_all_features('hide', 0)

        gdfs1.set_all_features('color', colors.red)
        gdfs2.set_all_features('color', colors.blue)

        gdt1 = Track('CDS features', greytrack=1,
            scale_largetick_interval=1e4,
            scale_smalltick_interval=1e3,
            scale_format = "SInt")
        gdt1.add_set(gdfs1)

        gdt2 = Track('gene features', greytrack=1,
                   scale_largetick_interval=1e4)
        gdt2.add_set(gdfs2)
                
        gdt3 = Track('misc features', greytrack=1,
                   scale_largetick_interval=1e4)
        gdt3.add_set(gdfs3)
        gdt3.add_set(gdfs4)        

        step = (end-start)/2000
        gdgs = GraphSet('GC Content')
        
        graphdata1 = apply_to_window(genbank_entry.seq, step, calc_gc_skew, step)
        gdgs.new_graph(graphdata1, 'GC Skew', style='bar',
                color=colors.violet,
                altcolor=colors.purple)
        
        gdt4 = Track(\
                'GC Skew (bar), GCContent(green line), ATContent(red line)',
                height=1.94, greytrack=1,
                scale_largetick_interval=1e4)
        gdt4.add_set(gdgs)

        graphdata2 = apply_to_window(genbank_entry.seq, step, calc_gc_content, step)
        gdgs.new_graph(graphdata2, 'GC content', style='line', 
                color=colors.lightgreen,
                altcolor=colors.darkseagreen)

        graphdata3 = apply_to_window(genbank_entry.seq, step, calc_at_content, step)
        gdgs.new_graph(graphdata3, 'AT content', style='line', 
                color=colors.orange,
                altcolor=colors.red)    
        
        gdd.add_track(gdt1, 2)
        gdd.add_track(gdt2, 4)
        gdd.add_track(gdt4, 5)
        gdd.add_track(gdt3, 6)

        gdd.draw(format='circular', orientation='landscape',
             tracklines=0, pagesize='A0', circular=True,
             start=start, end=end)
        output_filename = os.path.join('Graphics', 'DiagramTestCircular.pdf')
        gdd.write(output_filename, 'PDF')
        
        gdd.draw(format='linear', orientation='landscape',
             tracklines=0, pagesize='A0', fragments=20,
             start=start, end=end)
        output_filename = os.path.join('Graphics', 'DiagramTestLinear.pdf')
        gdd.write(output_filename, 'PDF')

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
