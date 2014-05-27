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
    py3kpath = os.path.join(biopath,"_py3k", "__init__.py")
    _py3k = imp.load_source("_py3k", py3kpath)
    _is_int_or_long = _py3k._is_int_or_long
else:
    from Bio._py3k import _is_int_or_long
    from ..SeqRecord import SeqRecord, _RestrictedDict
    
from copy import copy
from math import floor, ceil, log


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

class FeatureBinCollection(object):
    """this class manages the creation and maintenance of feature indices

       This class is used to organize feature data in a quickly retrievable
       data structure. The feature data must be added as a tuple containing
       at least two indices: first annotated residue and the last as a half
       open half closed interval [first, last). The indices are assumed to be
       the first two elements of the stored tuple, but they may be re-assigned
       on instantiation via the beginindex and endindex kwarks.

       EXAMPLE
       -------
       defined below is a 3-tuple format of (beginindex, endindex, fileidx)
       three features are added to a newly initialized featurebin 

       >>> ft0 = (5574, 5613, 2300) 
       >>> ft1 = (0, 18141, 1300 )
       >>> ft2 = (5298, 6416, 3540)
       >>> featurebin = FeatureBinCollection()
       >>> featurebin.insert( ft0 )
       >>> featurebin.insert( ft1 )
       >>> featurebin.insert( ft2 )
       >>> len(featurebin)
       3
       
       Now that the 'featurebin' instance has some features, they can be
       retrieved with a standard getter using single integer indices or
       slice notation.

       >>> featurebin[1]
       [(0, 18141, 1300)]
       >>> sliceresult = featurebin[5200:5300]
       >>> sliceresult.sort()
       >>> sliceresult
       [(0, 18141, 1300), (5298, 6416, 3540)]


       BACKGROUND:
       -----------
       The basic idea of using feature bins is to group features into 
       bins organized by their span and sequence location. These bins then allow
       only likely candidate features to be queried rather than all features. The 
       example below illustrated with Figure 1 shows a similar scheme where feature1 
       is stored in bin-0, feature2 in bin-4 and feature3 in bin-2. Each sequence is
       stored in the smallest bin that will fully contain the sequence. A query of 
       all features in the region denoted by query1 could be quickly performed by 
       only searching through bins 0, 2, 5, and 6. Were this data structure many 
       levels deep, the performance savings would be large

       ___Figure 1_________________________________________________
       |                                                           |
       |    feature1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~              |
       |    feature2  |       ~~~~                  |              |
       |    feature3  |       |  |               ~~~~~~~~~~~~~     |
       |              |       |  |               |  |        |     |
       | bins:        |       |  |               |  |        |     |
       |    0_________|_______|__|_______._______|__|____.___|_    |
       |    1_________________|__|__   2_:_______|_______:___|_    |
       |    3__________  4____|__|__   5_:________  6____:_____    |
       |                                 :               :         |
       |                                 :               :         |
       |    query1                       [ ? ? ? ? ? ? ? ]         |
       |...........................................................|

       Further reading on the math behind the idea can be found in: 
           Journal:  Bioinformatics Vol 27 no. 5 2011, pages 718-719
           Article:  "Tabix: fast retrieval of sequence features from generic
                      tab delimited files"
           Author:   Heng Li

       The implementation by Li has, as its largest bin ~500 million (2^29) and its smallest
       bin ~16000 (2^14). Each level of binning is separated by a factor of 8 (2^3).
       The implementation herein abandons a static binning scheme and instead
       starts with the smallest and largest bins as 256 and 8 million respectively. 
       These bins can then be dynamically expanded increasing by a factor of 8
       every time new data is found to be larger than the largest bin. As a practical
       matter of sanity checking, bin sizes are capped at 2.2 trillion residues (2^41).

       Under some circumstances the exact size of a sequence and all related annotations
       is known beforehand. If this is the case the length kwarg allows the binning object
       to be solidified on instantiation at the correct length.
       
       This structure knows nothing about the global sequence index and is indexed 
       at zero. Any index transformation must be done at a higher level. It is important
       that all sequences and features stored here are indexed to zero.
       """
    
    def __init__(self, length = None, beginindex=0, endindex=1):
        """ initialize the class and set standard attributes

        kwargs:

          length:
            when length == None, the bins are dynamically sized.
            when length is a positive integer, the appropriate bin 
            size is selected and locked. Exceeding this value will
            cause exceptions when the max bin size is locked

          beginindex:
            the index of the first residue within the tuple that will
            be stored with the FeatureBinCollection.

          endindex:
            the index of the last residue (as a open interval) inside 
            the tuple that will be stored with FeatureBinCollection
        """ 
        #these should not be changed
        self._bin_level_count = 6
        self._bins = [[] for i in range(37449)]

        # this defines the indices of the begin and end sequence info
        # in the tuple structures stored in the bins
        self._beginindex = beginindex
        self._endindex = endindex

        #default action: start small (8M) and allow expansion
        self._sorted = False
        self._dynamic_size = True
        if length is None:
            self._set_max_bin_power(23)
            
        #alternate action if a sequence length is provided
        # set to smallest power able to fully contain
        elif _is_int_or_long(length) and length > 0:
            default_powers = [23,26,29,32,35,38,41]
            for power in default_powers:
                if length <= 2**power:
                    self._set_max_bin_power(power)
                    self._dynamic_size = False
                    break
            if self._dynamic_size: #this should have been set to False
                error_string = "Sequence length is {}: must be less than 2^41".format(length)
                raise ValueError(error_string)
        
    def _increase_bin_sizes(self):
        """increase max bin size 8x (2**3) and re-organize existing binned data
        
        In order to increase the total maximum bin size, the lowest set
        of bins must be merged up one level (first step) then the entire set
        must be moved down one level without disturbing the organization scheme.
        
        An assertion in this routine blocks sequences larger than 2**41 from
        being created.
        """  
        oldsizepower = self._max_bin_power
        newsizepower = oldsizepower + 3
        assert newsizepower <= 41
        self._set_max_bin_power(newsizepower)
        
        # first, remove the lowest level
        # by merging it up to the previous level
        level = 5
        oL = int((2**(3*level) - 1)/7.0)
        new_level = 4
        new_oL = int((2**(3*new_level) - 1)/7.0)
        old_size = 2**(oldsizepower - 3*level)
        new_size = 2**(oldsizepower - 3*new_level)
        for k in range(4681, 37449):    
            bin_begin_old = (k-oL)*old_size
            k_new = int(floor(new_oL + (bin_begin_old/new_size)))
            #extend required to save existing data
            self._bins[k_new].extend(self._bins[k])
            self._bins[k] = []
        
        #then, move everything down.
        for k_inverse in range(4681):
            k = 4680 - k_inverse
            level = int( floor( log((7*k + 1),2)/3.0 ) )
            new_level = level + 1
            oL = int((2**(3*level) - 1)/7.0)
            new_oL = int((2**(3*new_level) - 1)/7.0)
            
            new_index = k - oL + new_oL 
            self._bins[new_index] = self._bins[k]
            self._bins[k] = []
               
    def _set_max_bin_power(self, power):
        """sets the maximum bin power and fixes other necessary attributes"""
        
        self._max_bin_power = power
        self._min_bin_power = self._max_bin_power - 3*(self._bin_level_count-1)
        self._size_list = [2**(self._min_bin_power+3*n) for n in range(self._bin_level_count)]
        self._max_sequence_length = self._size_list[-1]
    
    def insert(self, feature_tuple):
        """inserts a tuple with a sequence range into the feature bins
        
        data is assumed to be somewhat scrubbed, coming from a parser
        or a parser consumer."""
        
        beginindex = self._beginindex
        endindex = self._endindex

        #reset sorted quality
        self._sorted = False

        begin = feature_tuple[beginindex]
        end = feature_tuple[endindex]
        #use _py3k module later
        assert _is_int_or_long(begin)
        assert _is_int_or_long(end)
        assert begin <= end
        span = end-begin
        
        bin_index = self._calculate_bin_index(begin, span)
        self._bins[bin_index].append(feature_tuple)
        
    def __len__(self):
        return sum(len(bin) for bin in self._bins) 

    def sort(self):
        """this performs bin-centric sorting, necessary for faster retrieval"""
        #bins must be sorted by the begin index, this is fastest
        if self._beginindex == 0:
            for i in range(len(self._bins)):
                self._bins[i].sort()
        #this is a bit slower but accomodates diverse data structures
        else:
            for i in range(len(self._bins)):
                self.bins[i].sort(key = lambda tup: tup[self._beginindex])
        #reset sorted quality
        self._sorted = True
            
    def __getitem__(self, key):
        """This getter efficiently retrieves the required entries
        
        This getter primarily works as expected and in a pythonic
        fashion, one exception it it's treatment of slices indices
        where the start is greater than the stop. Rather than just 
        throwing calculated output, an IndexError is raised.
        """    
        #set some locals
        beginindex = self._beginindex
        endindex = self._endindex
        
        #check that it is a slice and it has no step property (or step==1)
        if not isinstance(key, slice):
            if _is_int_or_long(key):
                key = slice(key, key+1)
            else:
                raise TypeError("lookups in the feature bin must use slice or int keys")
        
        #any integers are just converted to a 'len() == 1' slice
        if key.step is not None and key.step != 1:
            raise KeyError("lookups in the feature bin may not use slice stepping ex. bins[0:50:2]")
        
        #fix begin or end index for slicing of forms: bins[50:] or bins[:50] or even bins[:]
        keystart, keystop, keystep = key.indices(self._max_sequence_length)
        if keystart > keystop:
            raise IndexError("key not valid, slice.start > slice.stop")
        
        #pre-sort if necessary
        if not self._sorted:
            self.sort()
            
        #code taken from self._calculate_bin_index(), comments removed
        return_entries = []
        possible_entries = []
        bin_level_count = self._bin_level_count
        max_bin_power = self._max_bin_power 
        for l_inverse in range(bin_level_count):
            L = bin_level_count - 1 - l_inverse
            oL = (2**(3*L) - 1)/7
            sL = float(2**(max_bin_power-3*L))
            k1 = int(floor(oL + (keystart/sL)))
            #k2 is incremented since range is used
            k2 = int(ceil(oL - 1 + (keystop)/sL)) + 1
            if k2-k1 > 2:
                for bin in range(k1+1, k2-1):
                    return_entries.extend( self._bins[bin]) 
            for binn in set([k1,k2-1]):
                #for binn in range(k1,k2):
                for feature in self._bins[binn]:
                    #this covers fully bound sequence and left overlap
                    if keystart <= feature[beginindex] < keystop:
                        return_entries.append(feature)
                    #this covers left sequence right sequence overlap 
                    elif keystart < feature[endindex] <= keystop:
                        return_entries.append(feature)
                    #this covers seqyebces fully bound by a feature      
                    elif keystart > feature[beginindex] and\
                         keystop < feature[endindex]:
                        return_entries.append(feature)
                    if keystop < feature[beginindex]:
                        break
        return return_entries

    def _calculate_bin_index(self, begin,span):
        """ This function returns a bin index given a (begin, span) interval
        
        The equations for determination of bin index derived from this publication: 
           Journal:  Bioinformatics Vol 27 no. 5 2011, pages 718-719
           Article:  "Tabix: fast retrieval of sequence features from generic
                      tab delimited files"
           Author:   Heng Li
        
        This function should only be used privately in the context of having no
        easier relationship to assign bin index. Placing this in a loop for any
        task other than arbitrary assignments is a bad idea since many of the
        common tasks can provide bin index through other relationships.
        _increase_bin_sizes is an example of a routine that would suffer
        performance penalties were it to use this routine yet can be run
        efficiently using alternate bin index relationships.
        """
        
        #take care base cases with the span parameter
        assert span >= 0
        #this is required for determination of bin location for zero length seq's
        span = max(1, span)
        
        # fix bin sizes if needed. Also, if the bin size is
        # larger than expected, do some self-consistency checks
        while begin+span > self._max_sequence_length:
            if self._dynamic_size:
                assert begin+span <= 2**41   # len(seq) > 2.19 trillion is not reasonable
            elif not self._dynamic_size and begin+span > 2**self._max_bin_power:
                error_string = "feature index at {}: must be less than 2^{}".format \
                                                (begin+span, self._max_bin_power)
                raise ValueError(error_string)
            self._increase_bin_sizes()
            
        #run the assignment loop
        bin_level_count = self._bin_level_count
        max_bin_power = self._max_bin_power 
        for l_inverse in range(bin_level_count):
            level = bin_level_count - 1 - l_inverse
            # calculate offset (oL) of the list at level L
            offset_at_L = (2**(3*level) - 1)/7
            group_length = (2**(3*(level+1)) - 1)/7
            #calculate size (sL) of the list: the number of residues in width
            size_at_L = float(2**(max_bin_power-3*level))
            # interval[0] >= (k - oL)*sL
            # rearrange to form
            # k =< oL + (interval[0])/sL
            k1 = int(floor(offset_at_L + (begin/size_at_L)))
            # interval[1] < (k - oL + 1)*sL
            # rearrange to form
            #k > 1 + oL + (interval[1])/sL
            #k2 = oL - 1 + (begin+span)/sL  
            k2 = int(ceil(offset_at_L - 1 + (begin+span)/size_at_L))
            if k1 == k2 and k1 < group_length:
                return k1
        
        assert False # the assignment loop failed
              
    
    
    
    
    
    
if __name__ == "__main__":
    import unittest
    #from Bio import SeqRecord
    #from Bio.SeqIO import _lazy
    #from Bio.Seq import Seq
    seqpath = os.path.join(biopath,"Seq.py")
    Seq = imp.load_source("Seq", seqpath)
    Seq = Seq.Seq
    

    class TestSeqRecordBaseClass(SeqRecordProxyBase):
        """ this class implements the minimum functionality for a working proxy
        
        This class must be defined in the test suite to implement the basic hook
        methods used by the newer proxy class. Additional parsers should add rudundant
        code, but proper testing here will catch errors before they are invoked
        by a derived class.
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
            
    class TestFeatureBinCollection(unittest.TestCase):
        def setUp(self):
            self.bins = FeatureBinCollection()
        
        def test_initial_state_max_bin_power(self):
            self.assertEqual(self.bins._max_bin_power, 23)

        def test_initial_state_max_count_of_bins(self):
            self.assertEqual(len(self.bins._bins), 37449)

        def test_bin_index_finder_smallest_bin_and_leftmost(self):
            k_0_256 = self.bins._calculate_bin_index(0,256)
            self.assertEqual(k_0_256, 4681) #this is the leftmost level 5 bin
        
        def test_bin_index_finder_negative_span(self):
            #k_0_256 = self.bins._calculate_bin_index(0,-56)
            #self.assertEqual(k_0_256, 4681) #this is the leftmost level 5 bin
            self.assertRaises(AssertionError, self.bins._calculate_bin_index, \
                              0, -1*56)

        def test_bin_index_finder_negative_start(self):
            #k_0_256 = self.bins._calculate_bin_index(-20,56)
            self.assertRaises(AssertionError, self.bins._calculate_bin_index, \
                              -20, 56)
            #self.assertEqual(k_0_256, 4681) #this is the leftmost level 5 bin

        def test_bin_index_finder_smallest_bin_and_2ndleftmost(self):
            k_256_256 = self.bins._calculate_bin_index(256,256)
            #this is the second leftmost level 5 bin
            self.assertEqual(k_256_256, 4682) 
            
        def test_bin_index_finder_smallest_bin_and_2ndleftmost_zero_length(self):
            k_256_0 = self.bins._calculate_bin_index(256,0)
            self.assertEqual(k_256_0, 4682) 
        
        def test_bin_index_finder_smallest_bin_last_index(self):
            k_1_256 = self.bins._calculate_bin_index(1,256)
            self.assertEqual(k_1_256, 585) 
         
        def test_bin_index_finder_smallest_bin_last_index(self):
            k_big_256 = self.bins._calculate_bin_index(8388352,256)
            self.assertEqual(k_big_256, 37448) # this is the rightmost lv5 bin
        
        def test_bin_index_finder_largest_bin_and_ensure_no_recalculation(self):
            k_small_big = self.bins._calculate_bin_index(0, 8388608)
            self.assertEqual(k_small_big, 0)
            self.assertEqual(self.bins._max_bin_power, 23)
            #value below should not have changed
            
        def test_bin_index_finder_smallest_bin_largest_index(self):
            """ this should return the index of the last bin"""
            k_big_256 = self.bins._calculate_bin_index(8388352,256)
            self.assertEqual(k_big_256, 37448)
            self.assertEqual(self.bins._max_bin_power, 23)
          
        def test_binning_size_changer_once(self):
            """ test that the bin size chnges when list is 1-over"""            
            k_overFull = self.bins._calculate_bin_index(1,8388608)
            k_full = self.bins._calculate_bin_index(0,8388608)
            self.assertEqual(self.bins._max_sequence_length, 8388608*8)
            self.assertNotIn(256, self.bins._size_list)
            self.assertIn(2048, self.bins._size_list)
            self.assertEqual(k_overFull, 0)
            self.assertEqual(k_full, 1)
            self.assertEqual(self.bins._max_bin_power, 26)
        
        def test_binning_size_changer_multiple_steps(self):
            """trigger a size rearrangement twice"""
            k_overFull = self.bins._calculate_bin_index(1,8388608*8)
            self.assertEqual(self.bins._max_sequence_length, 8388608*8*8)
            k_full = self.bins._calculate_bin_index(0,8388608*8)
            self.assertEqual(k_overFull, 0)
            self.assertEqual(k_full, 1)
            self.assertEqual(self.bins._max_bin_power, 29)
            self.assertNotIn(2048, self.bins._size_list)
        
        def test_insertion_bad_negative_val(self):
            testTuple1 = (-1, 56)
            self.assertRaises(AssertionError,self.bins.insert,testTuple1)
            
        def test_insertion_once_smallest(self):
            testTuple1 = (0, 256)
            self.bins.insert(testTuple1)
            self.assertIn(testTuple1, self.bins._bins[4681])
            
        def test_insertion_once_smallest_but_overlaps(self):
            testTuple2 = (1, 257)
            self.bins.insert(testTuple2)
            self.assertIn(testTuple2, self.bins._bins[585])
            
        def test_insertion_once_smallestbin_rightmost(self):
            testTuple3 = (8388608-256, 8388608)
            self.bins.insert(testTuple3)
            self.assertIn(testTuple3, self.bins._bins[37448])

        def test_insertion_zero_length_smallestbin_rightmost(self):
            testTuple3 = (8388607, 8388607)
            self.bins.insert(testTuple3)
            self.assertIn(testTuple3, self.bins._bins[37448])
            
        
        def test_insertion_once_smallestbin_rightmost_lv4(self):
            testTuple4 = (8388608-257, 8388608)
            self.bins.insert(testTuple4)
            self.assertIn(testTuple4, self.bins._bins[4680])    

        def test_insertion_with_rearrangement(self):
            testTuple3 = (256, 256+256)
            self.bins.insert(testTuple3)
            self.assertIn(testTuple3, self.bins._bins[4682])
            #trigger a size rearrangement with an insertion
            testTuple5 = (0, 9000000)
            self.bins.insert(testTuple5)
            self.assertIn(testTuple3, self.bins._bins[4681])
            self.assertIn(testTuple5, self.bins._bins[0])
            self.assertEqual(self.bins._max_sequence_length, 8388608*8)
            self.assertNotIn(256, self.bins._size_list)
            self.assertIn(2048, self.bins._size_list)
            self.assertEqual(self.bins._max_bin_power, 26)
        
        def test_insertion_level3_and_level5_then_resize_to_fit_in_same_bin(self):
            test_tuple_lv3 = (16384 , 16384+16384)
            test_tuple_lv5 = (16384+16384-256 , 16384+16384)
            self.bins.insert(test_tuple_lv3)
            self.bins.insert(test_tuple_lv5)
            self.assertIn(test_tuple_lv3, self.bins._bins[74])
            self.assertIn(test_tuple_lv5, self.bins._bins[4681+127])
            #trigger rearrangement by 2 levels
            k_overFull = self.bins._calculate_bin_index(1,8388608*8)
            self.assertIn(test_tuple_lv3, self.bins._bins[4682])
            self.assertIn(test_tuple_lv5, self.bins._bins[4682])
            self.assertEqual(self.bins._max_bin_power, 29)
            
            
        def test_overflows_of_static_defined_lists(self):
            staticbins = FeatureBinCollection(length=67108864)
            #insert a chunk of data
            testTuple1 = (1, 2049)
            staticbins.insert(testTuple1)
            self.assertIn(testTuple1, staticbins._bins[585])
            #test that exceptions are raised if a over-sized bin is used
            overSizedTuple1 = (0, 67108865)
            overSizedTuple2 = (67108864, 67108865)
            self.assertRaises(ValueError, staticbins.insert, overSizedTuple1)
            self.assertRaises(ValueError, staticbins.insert, overSizedTuple2)
            
        def test_oversized_definition_of_the_collection(self):    
            reallyReallyoversizedEnd = 1+2**41
            self.assertRaises(ValueError, FeatureBinCollection, reallyReallyoversizedEnd)

        def test_getter_get_values_from_empty_set(self):
            resultsEmpty = self.bins[2**20]
            self.assertEqual([], resultsEmpty)
            
        def test_getter_get_values_edge_cases(self):
            resultsEdgeRight = self.bins[-1+2**23]
            resultsEdgeLeft = self.bins[0]
            self.assertEqual([], resultsEdgeRight)
            self.assertEqual([], resultsEdgeLeft)
            
        def test_getter_get_values_from_out_of_bounds(self):
            self.assertRaises(IndexError, self.bins.__getitem__, -1)
            self.assertRaises(IndexError, self.bins.__getitem__, 1+2**23)
            
        def test_getter_typeError_string(self):
            self.assertRaises(TypeError, self.bins.__getitem__, "hello")
            
        def test_getter_typeError_float(self):
            self.assertRaises(TypeError, self.bins.__getitem__, 5.45)          
            
        def test_getter_typeError_stepped_slice(self):
            self.assertRaises(KeyError, self.bins.__getitem__, slice(0,25,2))
        
        def test_getter_reversed_index(self):
            resultR = (20000,30000)
            self.bins.insert(resultR)
            self.assertRaises(IndexError, self.bins.__getitem__, slice(23000,21000))
            
        def test_getter_half_slices(self):
            emptyL = self.bins[:50]
            emptyR = self.bins[50:]
            emptyM = self.bins[:]
            self.assertEqual([], emptyL)
            self.assertEqual([], emptyR)
            self.assertEqual([], emptyM)

        def test_getter_half_slices_with_result_right_border(self):
            resultR = (20000,30000)
            self.bins.insert(resultR)
            emptyL = self.bins[:20000]
            emptyR = self.bins[20000:]
            emptyM = self.bins[:]
            self.assertEqual([], emptyL)
            self.assertIn(resultR, emptyR)
            self.assertIn(resultR, emptyM)
        
        def test_getter_half_slices_with_result_left_border(self):
            resultL = (10000,20000)
            self.bins.insert(resultL)
            emptyL = self.bins[:20000]
            emptyR = self.bins[20000:]
            emptyM = self.bins[:]
            self.assertEqual([], emptyR)
            self.assertIn(resultL, emptyL)
            self.assertIn(resultL, emptyM)
            
        def test_insertion_where_medium_sized_bin_is_out_of_bounds(self):
            resultL = (8388608, 8388608+2047)
            self.bins.insert(resultL)
            self.assertIn(resultL, self.bins[838840:8388610])
            self.assertEqual(self.bins._max_bin_power, 26)

        def test_getter_overlap_left(self):
            feature = (100000,200000)
            self.bins.insert(feature)
            result = self.bins[99000:101000]
            self.assertIn(feature, result)
        
        def test_getter_overlap_right(self):
            feature = (100000,200000)
            self.bins.insert(feature)
            result = self.bins[199000:201000]
            self.assertIn(feature, result)

        def test_getter_feature_inside_region(self):
            feature = (100000,200000)
            self.bins.insert(feature)
            result = self.bins[99000:201000]
            self.assertIn(feature, result) 
            
        def test_getter_region_inside_feature(self):
            feature = (100000,200000)
            self.bins.insert(feature)
            result = self.bins[101000:199000]
            self.assertIn(feature, result) 
        
        def test_getter_zero_length_inside(self):
            testTuple3 = (8388605, 8388605)
            self.bins.insert(testTuple3)
            self.assertIn(testTuple3,self.bins[8388604:8388606])
            self.assertIn(testTuple3,self.bins[8388604:8388605])
            self.assertIn(testTuple3,self.bins[8388605:8388606])
            self.assertEqual([],self.bins[8388603:8388604])
    
    a = TestSeqRecordBaseClass("seQUencefake", "fakeid")
    unittest.main( exit=False )
    """if __name__ == "__main__":
        runner = unittest.TextTestRunner(verbosity = 2)
        unittest.main(testRunner=runner)"""
