# Copyright 2009-2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Lazy parsing and feature indexing of sequence files (PRIVATE).

You are not expected to access this module, or any of its code, directly. This
is all handled internally by the following functions when accessed with the
lazy=True kwarg.
Bio.SeqIO.parse(.... lazy=True)
Bio.SeqIO.index(..., lazy=True)
Bio.SeqIO.read(..., lazy=True) 


The lazy loading parsers will read an entire record efficiently storing and 
parsing only the mimimum amount of information to provide fast lookup. Records
returned from parse(), index(), and read() will act as a proxy to a regular 
SeqRecord class.

The lazy loading strategy is, for large sequences, faster to initialize than
SeqIO.parse() and more memory efficient. It is slower to initialize than
SeqIO.index() but also more memory efficient.

The lazy loader will partially parse a sequence file noting the sequence
span and file position of sequence annoations. These annotations will
be stored for efficient retrieval from the file when requested. Sequence 
data will also be efficiently pulled from structured regions
and file

Note that this means all parsing is on demand, so improperly formatted
records may not trigger an exception unless the problem region is requested.
This is by design.
"""

if __name__ == "__main__":
    #this junk allows relative importing for testing in situ
    #will not be included once all unit tests have been exported to
    #the external unit test suite for builds
    import imp
    import os
    biopath = os.path.split(os.getcwd())[0]
    srpath = os.path.join(biopath,"SeqRecord.py")
    SeqRecord = imp.load_source("SeqRecord", srpath)
    _RestrictedDict = SeqRecord._RestrictedDict
    SeqRecord = SeqRecord.SeqRecord
else:
    from ..SeqRecord import SeqRecord, _RestrictedDict
    
from copy import copy

class SeqRecordProxyBase(SeqRecord):
    """A SeqRecord object holds a sequence and information about it.

    Main attributes:
     - id          - Identifier such as a locus tag (string)
     - seq         - The sequence itself (Seq object or similar)

    Additional attributes:
     - name        - Sequence name, e.g. gene name (string)
     - description - Additional text (string)
     - dbxrefs     - List of database cross references (list of strings)
     - features    - Any (sub)features defined (list of SeqFeature objects)
     - annotations - Further information about the whole sequence (dictionary).
                     Most entries are strings, or lists of strings.
     - letter_annotations - Per letter/symbol annotation (restricted
                     dictionary). This holds Python sequences (lists, strings
                     or tuples) whose length matches that of the sequence.
                     A typical use would be to hold a list of integers
                     representing sequencing quality scores, or a string
                     representing the secondary structure.

    You will typically use Bio.SeqIO to read in sequences from files as
    SeqRecord objects.  However, you may want to create your own SeqRecord
    objects directly (see the __init__ method for further details):



    #subclassing to-do list
    attributes that must be managed or utilized 
     -   self._seq         :::   sequence type, set by _read_seq
     -   self.name         :::   this should be identical
     -   self.dbxrefs      :::   [] db cross references
     -   self.id           :::   id string
     -   self._features    :::   list of SeqFeatures
     -   self._annotations :::   dict of annotations
     -   self._letter_an...:::   per letter information


     -   self._index_begin :::   defines begin index w/r/t global
     -   self._index_end   :::   defines end index w/r/t global

     

    methods need implementation by derived class:
     -      self.__init__    ::: the init provided in the base class is for testing
     -      self._read_seq    ::: gets the sequence from file and sets it



    """
    def __init__(self):
        raise NotImplementedError( \
            "__init__ must be implemented in the derived class")
    
    
    def _return_seq(self):
        """this simple function removes getter logic from _read_seq"""
        if not self._seq:
            self._read_seq()
        return self._seq
    
    def _read_seq(self):
        """this is implemented to handle file access for setting _seq"""
        raise NotImplementedError( \
            "_read_seq must be implemented in the derived class")

    seq = property(fget=_return_seq,
                   doc="The sequence itself, as a Seq or MutableSeq object.")

    #TODO - letter annotations?
    letter_annotations = None

    def __getitem__(self, index):
        """Returns a sub-sequence or an individual letter.

        This will need to cleverly either return a sequence letter
        or a copy of itself with new marker indices.
        """
        if isinstance(index, int):
            #This mimics the behavior of the current SeqRecord
            #The recursive call is required to prevent full parsing while
            #only calling a single resiude
            return self[index:index+1].seq[0]
        
        elif isinstance(index, slice):
            parent_length = len(self)
            if parent_length <= 0:
                raise ValueError("If the sequence length is zero, we cannot slice it.")
            
            #this is what will be returned
            seq_proxy_copy = copy(self)
            #do some index math
            start, stop, step = index.indices(parent_length)
            seq_proxy_copy._index_begin = self._index_begin + start
            seq_proxy_copy._index_end = self._index_begin + stop
            #fix _seq property when set
            if self._seq:
                seq_proxy_copy._seq = self._seq[start:stop]

            return seq_proxy_copy
        raise ValueError("Invalid index")
    
    
    def __len__(self):
        """Returns the length of the sequence.
        
        The derived class must set attribute _seq_len
        """
        return self._index_end - self._index_begin


    def upper(self):
        """Returns a copy of the record with an upper case sequence.

        
        """
        if not self._seq:
            self._read_seq()
        newSelf = copy(self)
        newSelf._seq = self._seq.upper()
        return newSelf


    def lower(self):
        """Returns a copy of the record with a lower case sequence.
        """
        if not self._seq:
            self._read_seq()
        newSelf = copy(self)
        newSelf._seq = self._seq.lower()
        return newSelf

    def __repr__(self):
        """this is a shortened repr for the lazy loading parsers

        modification of the repr() magic method is required to prevent
        simple command line options from unnecessarily invoking full
        parsing of a file. 

        This strategy of redefining __repr__ contrasts with the __str__
        method that is left to the base class since invoking the str() is
        less used and more likely in output focused programs. 
        """
        idrepr = "id=%s" % str(self.id)
        namerepr = "name=%s" % str(self.name)
        descriptionrepr = "description=%s" % str(self.description)
        dbxrefsrepr = "dbxrefs=%s" % repr(self.dbxrefs)
        if self._seq is None:
            seqrepr = "seq=NOT_READ"
        else:
            seqrepr = "seq=%s" % repr(self.seq)                 
        return self.__class__.__name__ \
         + "(%s, %s, %s, %s, %s)" \
         % (seqrepr, idrepr, namerepr, descriptionrepr, dbxrefsrepr)


    # All methods tagged below are implemented in the base class
    #
    #def __bool__(self):
    #__nonzero__= __bool__
    #def reverse_complement(self, id=False, name=False, description=False,
    #def __radd__(self, other):
    #def format(self, format):
    #def __format__(self, format_spec):
    #def __add__(self, other):
    #    returns non-proxy class where possible
    #def __radd__(self, other):
    #    returns non-proxy class where possible
    #def __str__(self):
    #def __repr__(self):
    #def __contains__(self, char):
    #def __iter__(self):
if __name__ == "__main__":
    import unittest
    #from Bio import SeqRecord
    #from Bio.SeqIO import _lazy
    #from Bio.Seq import Seq
    seqpath = os.path.join(biopath,"Seq.py")
    Seq = imp.load_source("Seq", seqpath)
    Seq = Seq.Seq
    

    class TestSeqRecordBaseClass(SeqRecordProxyBase):
        """ this class implements the minimum funcationality for a working proxy
        """
        def __init__(self, handle, id = "<unknown id>", name = "<unknown name>",
                     description = "<unknown description>", dbxrefs = None,
                     features = None, annotations = None,
                     letter_annotations = None, index_begin=0,
                     index_end=None):
            """Create a SeqRecord.

            Arguments:
             - handle      - Sequence string, required (string)
             - id          - Sequence identifier, recommended (string)
             - name        - Sequence name, optional (string)
             - description - Sequence description, optional (string)
             - dbxrefs     - Database cross references, optional (list of strings)
             - features    - Any (sub)features, optional (list of SeqFeature objects)
             - annotations - Dictionary of annotations for the whole sequence
             - letter_annotations - Dictionary of per-letter-annotations, values
                                    should be strings, list or tuples of the same
                                    length as the full sequence.

            """
            self._index_begin = index_begin
            self._index_end = len(handle)
            self._seq = None

            if id is not None and not isinstance(id, basestring):
                #Lots of existing code uses id=None... this may be a bad idea.
                raise TypeError("id argument should be a string")
            if not isinstance(name, basestring):
                raise TypeError("name argument should be a string")
            if not isinstance(description, basestring):
                raise TypeError("description argument should be a string")
            self._handle = handle
            self.id = id
            self.name = name
            self.description = description

            # database cross references (for the whole sequence)
            if dbxrefs is None:
                dbxrefs = []
            elif not isinstance(dbxrefs, list):
                raise TypeError("dbxrefs argument should be a list (of strings)")
            self.dbxrefs = dbxrefs

            # annotations about the whole sequence
            if annotations is None:
                annotations = {}
            elif not isinstance(annotations, dict):
                raise TypeError("annotations argument should be a dict")
            self.annotations = annotations

            #TODO handle letter annotations properly
            self.letter_annotations = letter_annotations

            #TODO handle features annotations about parts of the sequence
            if features is None:
                features = []
            elif not isinstance(features, list):
                raise TypeError("features argument should be a list (of SeqFeature objects)")
            self.features = features
        
        def _read_seq(self):
            if self._seq:
                pass
            else:
                i = self._index_begin
                j = self._index_end 
                self._seq = Seq(self._handle[i:j])
    

    class SeqRecordProxyBaseClassTests(unittest.TestCase):

        def setUp(self):
            pass

        def test_nothing(self):
            """An addition test"""
            a = TestSeqRecordBaseClass("sequencefake", "fakeid")
            self.assertEqual(5, 5)
            self.assertRaises(NotImplementedError, SeqRecordProxyBase)

        def test_simple_tester(self):
            """checks on basic definitions in tester class"""
            a = TestSeqRecordBaseClass("sequencefake", "fakeid")
            self.assertEqual(0, a._index_begin)
            self.assertEqual(12, a._index_end)
            self.assertEqual(None, a._seq)
            self.assertEqual("sequencefake", a._handle)

        def test_seq_reader(self):
            """checks on basic worksing of _read_seq"""
            a = TestSeqRecordBaseClass("sequencefake", "fakeid")
            self.assertEqual(None, a._seq)
            self.assertEqual("sequencefake", str(a.seq))
            self.assertEqual("sequencefake", str(a._seq))

        def test_simple_index_grab(self):
            a = TestSeqRecordBaseClass("sequencefake", "fakeid")
            b = a[1:9]
            self.assertEqual(1, b._index_begin)
            self.assertEqual(9, b._index_end)

        def test_simple_index_grab_handle_identity(self):
            a = TestSeqRecordBaseClass("sequencefake", "fakeid")
            b = a[1:9]
            self.assertTrue(a._handle is b._handle)
            self.assertTrue(a is not b)

        def test_two_level_getter_indexes(self):
            a = TestSeqRecordBaseClass("sequencefake", "fakeid")
            b = a[1:9]
            c = b[1:8]
            self.assertEqual(2, c._index_begin)
            self.assertEqual(9, c._index_end)

        def test_two_level_getter_handle_identity(self):
            a = TestSeqRecordBaseClass("sequencefake", "fakeid")
            b = a[1:9]
            c = b[1:8]
            self.assertTrue(a._handle is c._handle)

        def test_upper(self):
            a = TestSeqRecordBaseClass("seQUencefake", "fakeid")
            b = a.upper()
            self.assertEqual("SEQUENCEFAKE", str(b._seq))

        def test_lower(self):
            a = TestSeqRecordBaseClass("seQUEncefake", "fakeid")
            b = a.lower()
            self.assertEqual("sequencefake", str(b._seq))

        def test_repr(self):
            out = r"""TestSeqRecordBaseClass(seq=NOT_READ, id=fakeid, """ + \
                  r"""name=<unknown name>, description=<unknown description>, dbxrefs=[])"""
            out2 = r"""TestSeqRecordBaseClass(seq=Seq('sequencefake', Alphabet()), id=fakeid, """ + \
                  r"""name=<unknown name>, description=<unknown description>, dbxrefs=[])"""
            a = TestSeqRecordBaseClass("sequencefake", "fakeid")
            self.assertEqual(out, repr(a))
            s = a.seq 
            self.assertEqual(out2, repr(a))

        def test_len_method(self):
            a = TestSeqRecordBaseClass("sequencefake", "fakeid")
            b = a[3:6]
            self.assertEqual(3, len(b))
            self.assertEqual(None, b._seq)
            self.assertEqual(len("sequencefake"), len(a))
    
    a = TestSeqRecordBaseClass("seQUencefake", "fakeid")
    unittest.main( exit=False )
    """if __name__ == "__main__":
        runner = unittest.TextTestRunner(verbosity = 2)
        unittest.main(testRunner=runner)"""
