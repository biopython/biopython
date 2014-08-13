import unittest

import sys
import imp
import os
import tempfile

from Bio import SeqIO
from Bio.Alphabet import single_letter_alphabet
from Bio._py3k import _bytes_to_string, _as_bytes, _is_int_or_long, StringIO
from Bio.Seq import Seq
from Bio.SeqIO._lazy import *
from Bio.SeqIO import InsdcIO
from Bio._py3k import _bytes_to_string, basestring
from io import BytesIO
from collections import namedtuple


#
### Basic fasta unit tests
#

#a 2 record unix formatted fasta
tsu1 = ">sp|O15205|UBD_HUMAN Ubiquitin D OS=Homo sapiens GN=UBD PE=1\n" + \
        "MAPNASCLCVHVRSEEWDLMTFDANPYDSVKKIKEHVRSKTKVPVQDQVLLLGSKILKPR\n" + \
        "RSLSSYGIDKEKTIHLTLKVVKPSDEELPLFLVESGDEAKRHLLQVRRSSSVAQVKAMIE\n" + \
        "RSLSSYGIDKEKTIHLTLKVVKPSDEELPLFLVESGDEAKRHLLQVRRSSSVAQVKAMIE\n" + \
        "TKTGIIPETQIVTCNGKRLEDGKMMADYGIRKGNLLFLACYCIGG\n" + \
        ">sp|P15516|HIS3_HUMAN Histatin-3 OS=Homo sapiens GN=HTN3 PE=1\n" + \
        "MKFFVFALILALMLSMTGADSHAKRHHGYKRKFHEKHHSHRGYRSNYLYDN\n"
#extra newline after ..YCIGG in tsu2
tsu2 = ">sp|O15205|UBD_HUMAN Ubiquitin D OS=Homo sapiens GN=UBD PE=1\n" + \
        "MAPNASCLCVHVRSEEWDLMTFDANPYDSVKKIKEHVRSKTKVPVQDQVLLLGSKILKPR\n" + \
        "RSLSSYGIDKEKTIHLTLKVVKPSDEELPLFLVESGDEAKRHLLQVRRSSSVAQVKAMIE\n" + \
        "RSLSSYGIDKEKTIHLTLKVVKPSDEELPLFLVESGDEAKRHLLQVRRSSSVAQVKAMIE\n" + \
        "TKTGIIPETQIVTCNGKRLEDGKMMADYGIRKGNLLFLACYCIGG\n\n" + \
        ">sp|P15516|HIS3_HUMAN Histatin-3 OS=Homo sapiens GN=HTN3 PE=1\n" + \
        "MKFFVFALILALMLSMTGADSHAKRHHGYKRKFHEKHHSHRGYRSNYLYDN\n"
#the exact sequence of O15205 [1:164]
tsuseq = "MAPNASCLCVHVRSEEWDLMTFDANPYDSVKKIKEHVRSKTKVPVQDQVLLLGSKILKPR" + \
         "RSLSSYGIDKEKTIHLTLKVVKPSDEELPLFLVESGDEAKRHLLQVRRSSSVAQVKAMIE" + \
         "RSLSSYGIDKEKTIHLTLKVVKPSDEELPLFLVESGDEAKRHLLQVRRSSSVAQVKAMIE" + \
         "TKTGIIPETQIVTCNGKRLEDGKMMADYGIRKGNLLFLACYCIGG"
#sigle fasta record
tsu3 = ">sp|P15516|HIS3_HUMAN Histatin-3 OS=Homo sapiens GN=HTN3 PE=1\n" + \
        "MKFFVFALILALMLSMTGADSHAKRHHGYKRKFHEKHHSHRGYRSNYLYDN\n"
#a windows formatted fasta
tsw1 = ">sp|O15205|UBD_HUMAN Ubiquitin D OS=Homo sapiens GN=UBD PE=1\r\n" + \
        "MAPNASCLCVHVRSEEWDLMTFDANPYDSVKKIKEHVRSKTKVPVQDQVLLLGSKILKPR\r\n" + \
        "RSLSSYGIDKEKTIHLTLKVVKPSDEELPLFLVESGDEAKRHLLQVRRSSSVAQVKAMIE\r\n" + \
        "RSLSSYGIDKEKTIHLTLKVVKPSDEELPLFLVESGDEAKRHLLQVRRSSSVAQVKAMIE\r\n" + \
        "TKTGIIPETQIVTCNGKRLEDGKMMADYGIRKGNLLFLACYCIGG\r\n" + \
        ">sp|P15516|HIS3_HUMAN Histatin-3 OS=Homo sapiens GN=HTN3 PE=1\r\n" + \
        "MKFFVFALILALMLSMTGADSHAKRHHGYKRKFHEKHHSHRGYRSNYLYDN\r\n"

double_padded = ">sp|O15205|UBD_HUMAN Ubiquitin D\r\n" + \
                "MAPNASCLCVHVRSEEW\r\n" + \
                "MANSANCLLCPPPPPPP\r\n" + \
                "MANSAC\r\n\r\n" + \
                ">spfake some label\r\n" + \
                "HHHHHHHH\r\n\r\n" + \
                ">spfake2 some other label\r\n" + \
                "IIIIIIII\r\n\r\n"

bad_header = ">sp|O15205|UBD_HUMAN Ubiquitin D\n" + \
             "MAPNASCLCVHVRSEEW\n" + \
             "MANSANCLLCPPPPPPP\n" + \
             "MANSAC\n" + \
             "badlyformatted some label\n" + \
             "HHHHHHHH\n"

bad_line = ">sp|O15205|UBD_HUMAN Ubiquitin D\n" + \
             "MAPNCVHVRSEEW\n" + \
             "MANSLCPPPPPDEM\n" + \
             "ANSLCPPPPPPPM\n" + \
             "ANSAC\n" + \
             ">wellformatted label\n" + \
             "HHHHHHHH\n"

test_seq_unix = BytesIO(_as_bytes(tsu1))
test_seq_gaped = BytesIO(_as_bytes(tsu2))
test_seq_win = BytesIO(_as_bytes(tsw1))
test_double_padded = BytesIO(_as_bytes(double_padded))
test_bad_line = BytesIO(_as_bytes(bad_line))
test_bad_header = BytesIO(_as_bytes(bad_header))


class LazyFastaIOSimpleTests(unittest.TestCase):

    def setUp(self):
        self.parser = SeqIO.FastaIO.FastaLazyIterator

    def test_iterator_finishes_win_type(self):
        """Test lazy iterator finishes with correct number of records"""
        lazyiterator = self.parser(test_seq_win)
        for count, rec in enumerate(lazyiterator):
            pass
        self.assertEqual(count, 1)
        lazyiterator = self.parser(test_seq_unix)
        for count, rec in enumerate(lazyiterator):
            pass
        self.assertEqual(count, 1)

    def test_iterator_returns_a_lazy_record(self):
        """Checks the lazy iterator returns the correct record"""
        first_seq = next(self.parser(test_seq_win))
        self.assertTrue(isinstance(first_seq,SeqIO.FastaIO.FastaSeqRecProxy))

    def test_indexing_simple_record_win(self):
        """Test indexing of a record with windows line endings"""
        first_seq = next(self.parser(test_seq_win))
        idx = first_seq._index.record
        self.assertEqual(idx["sequenceletterwidth"], 60)
        self.assertEqual(idx["sequencelinewidth"], 62)
        self.assertEqual(idx["recordoffsetstart"], 0)

    def test_indexing_2nd_record_unix(self):
        """Test indexing of a record with unix line endings"""
        lazyiterator = self.parser(test_seq_unix)
        first_seq = next(lazyiterator)
        second_seq = next(lazyiterator)
        idx = second_seq._index.record
        self.assertEqual(idx["sequenceletterwidth"], 51)
        self.assertEqual(idx["sequencelinewidth"], 52)
        self.assertEqual(idx["recordoffsetstart"], 290)

    def test_len_win(self):
        """Test that the 'seqlen' index value is correct"""
        lazyiterator = self.parser(test_seq_win)
        first_seq = next(lazyiterator)
        self.assertEqual(first_seq._index.record["seqlen"], 225)

    def test_indexing_2nd_record_win(self):
        """Test several index values on the second record"""
        lazyiterator = self.parser(test_seq_win)
        first_seq = next(lazyiterator)
        second_seq = next(lazyiterator)
        idx = second_seq._index.record
        self.assertEqual(idx["sequenceletterwidth"], 51)
        self.assertEqual(idx["sequencelinewidth"], 53)
        self.assertEqual(idx["recordoffsetstart"], 295)

    def test_sequence_getter_zero_unix(self):
        """Test sequence getting from unix endlines"""
        firstseq = next(self.parser(test_seq_unix))
        s = firstseq[0:10]
        self.assertEqual(str(s.seq), "MAPNASCLCV")

    def test_sequence_getter_zero_win(self):
        """Test sequence getting from unix endlines"""
        firstseq = next(self.parser(test_seq_win))
        s = firstseq[0:10]
        self.assertEqual(str(s.seq), "MAPNASCLCV")

    def test_double_padded_seq_len(self):
        """Test sequence getting when padded with extra newline"""
        firstseq = next(self.parser(test_double_padded))
        self.assertEqual(len(firstseq), 40)

    def test_double_padded_padding_difference(self):
        """Sequence getting 2nd record padded with extra newline"""
        lazyiterator = self.parser(test_double_padded)
        firstseq = next(lazyiterator)
        secondseq = next(lazyiterator)
        idx = firstseq._index.record
        firstseqend = idx["recordoffsetstart"] + idx["recordoffsetlength"]
        padding = secondseq._index.record["recordoffsetstart"] - firstseqend
        self.assertEqual(2, padding)

    def test_one_bad_sequence_line(self):
        """Test corrupt record raises"""
        firstseq = next(self.parser(test_bad_line))
        self.assertRaises(ValueError, firstseq._read_seq)

    def test_end_iteration(self):
        """Test StopIteration happens properly"""
        lazyiterator = self.parser(test_seq_unix)
        firstseq = next(lazyiterator)
        secondseq = next(lazyiterator)
        self.assertRaises(StopIteration, next, lazyiterator)

    def test_one_bad_sequence_header(self):
        """Test poorly formatted seq-header raises ValueError"""
        firstseq = next(self.parser(test_bad_header))
        self.assertRaises(ValueError, firstseq._read_seq)

    def test_sequence_getter_unix(self):
        """Test sequence getter for unix formatted endlines"""
        firstseq = next(self.parser(test_seq_unix))
        s = firstseq[6:10]
        self.assertEqual(str(s.seq), "CLCV")

    def test_sequence_getter_win(self):
        """Test sequence getter for windows formatted endlines"""
        firstseq = next(self.parser(test_seq_win))
        s = firstseq[6:10]
        self.assertEqual(str(s.seq), "CLCV")

    def test_sequence_getter_unix_2_line_span1(self):
        """Test seq getter slicing in for unix endlines"""
        firstseq = next(self.parser(test_seq_unix))
        s = firstseq[59:62]
        self.assertEqual(str(s.seq), "RRS")

    def test_sequence_getter_unix_2_line_span2(self):
        """Test seq getter in 2nd record with unix endlines"""
        lazyiterator = self.parser(test_seq_gaped)
        firstseq = next(lazyiterator)
        secondseq = next(lazyiterator)
        s = secondseq[40:43]
        self.assertEqual(str(s.seq), "RGY")

    def test_sequence_getter_unix_2_line_span3(self):
        """Test sequence getter for win endlines 2line span + 1 char"""
        firstseq = next(self.parser(test_seq_win))
        s = firstseq[58:62]
        self.assertEqual(str(s.seq), tsuseq[58:62])

    def test_sequence_getter_win_2_line_span(self):
        """tet sequence getter win endlines 2 line span"""
        firstseq = next(self.parser(test_seq_win))
        s = firstseq[59:61]
        self.assertEqual(str(s.seq), tsuseq[59:61])

    def test_sequence_getter_explicit_full_span(self):
        """test near-full sequence getter unix endlines"""
        firstseq = next(self.parser(test_seq_unix))
        s = firstseq[3:225]
        self.assertEqual(str(s.seq), tsuseq[3:225])

    def test_sequence_getter_inset_full_span(self):
        """Test sequence getter fullspan with unix endlines"""
        firstseq = next(self.parser(test_seq_unix))
        s = firstseq[0:161]
        self.assertEqual(len(tsuseq[0:161]), 161)
        self.assertEqual(len(s), 161)
        self.assertEqual(len(s.seq), 161)
        self.assertEqual(str(s.seq), tsuseq[0:161])

#
### tests for base class wrapper
#

class TestWrapperWithGenBank(unittest.TestCase):

    filename = os.path.join('GenBank', "brca_FJ940752.gb")
    filename2 = os.path.join('GenBank',"cor6_6.gb")

    def setUp(self):
        #db setup: make tempfile and close os handle to help with cleanup
        oshandle, dbfilename = tempfile.mkstemp()
        self.dbfilename = dbfilename
        os.close(oshandle)
        self.returncls = InsdcIO.GenbankSeqRecProxy

    def tearDown(self):
        os.remove(self.dbfilename)

    def test_make_one_record_with_db(self):
        """Test lazy class wrapper with database input"""
        #make the iter
        lazy_iter = LazyIterator(files = [self.filename],
                                 return_class = self.returncls,
                                 index = self.dbfilename)
        #test that it iterates and assigns temprec as a record
        temprec = None
        for r in lazy_iter:
            temprec = r
            self.assertTrue(isinstance(r, InsdcIO.GenbankSeqRecProxy))
        self.assertTrue(isinstance(temprec, InsdcIO.GenbankSeqRecProxy))

    def test_make_two_records_with_same_db(self):
        """test making 2 lazy records sequentially accessing same DB"""
        lazy_iter = LazyIterator(files = [self.filename],
                                 return_class = self.returncls,
                                 index = self.dbfilename)
        lazy_iter2 = LazyIterator(files = [self.filename2],
                                 return_class = self.returncls,
                                 index = self.dbfilename)
        self.assertEqual(lazy_iter.keys(), ["FJ940752.1"])
        self.assertEqual(lazy_iter2.keys(), ['X55053.1', 'X62281.1', \
                'M81224.1', 'AJ237582.1', 'L31939.1', 'AF297471.1'])
        for r in lazy_iter2:
            self.assertTrue(isinstance(r, InsdcIO.GenbankSeqRecProxy))
        for r in lazy_iter:
            self.assertTrue(isinstance(r, InsdcIO.GenbankSeqRecProxy))

    def test_make_one_record_without_db(self):
        """Test using the lazy iter without a db"""
        lazy_iter = LazyIterator(files = [self.filename],
                                 return_class = self.returncls,
                                 index = True)
        #test that it iterates and assigns temprec as a record
        temprec = None
        for r in lazy_iter:
            temprec = r
            self.assertTrue(isinstance(r, InsdcIO.GenbankSeqRecProxy))
        self.assertTrue(isinstance(temprec, InsdcIO.GenbankSeqRecProxy))

    def test_getter_using_brca_id(self):
        """Test get-by-ID using the lazy-iter"""
        lazydict = LazyIterator(files = [self.filename],
                         return_class = self.returncls,
                         index = self.dbfilename)
        del(lazydict)
        d = LazyIterator(files = [self.filename],
                         return_class = self.returncls,
                         index = self.dbfilename)
        self.assertTrue(isinstance(d["FJ940752.1"], \
                        InsdcIO.GenbankSeqRecProxy))


#
### tests for base class
#

class MinimalLazySeqRecord(SeqRecordProxyBase):
    """ this class implements the minimum functionality for a SeqRecordProxy

    This class must be defined in the test suite to implement the basic hook
    methods used by the new proxy class. Additional parsers should add
    redundant code, but proper testing here will catch some errors before
    they present in a derived class.
    """

    _format = "testclass"

    def _make_record_index(self, new_index):
        self._handle.seek(0)
        new_index["id"] = _bytes_to_string(self._handle.readline()).strip()
        seqstart = self._handle.tell()
        new_index["sequencestart"] = seqstart
        new_index["seqlen"] = len(self._handle.readline())
        return new_index

    def _read_seq(self):
        self._handle.seek(self._index.record["sequencestart"])
        seqstr = _bytes_to_string(self._handle.readline())
        self._seq = Seq(seqstr, self._alphabet)

    def _load_non_lazy_values(self):
        self.id = self._index.record["id"]
        self.name = "<unknown name>"
        self.description = "<unknown description>"
        self.dbxrefs = []

    def _make_feature_index(self, new_list):
        return new_list

    def _read_features(self):
        return []


class SeqRecordProxyBaseClassTests(unittest.TestCase):

    def setUp(self):
        self.ab = single_letter_alphabet

    def test_simple_tester(self):
        """checks on basic definitions in tester class"""
        handle = "fakeid\nsequencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        handle.name = "fake"
        a = MinimalLazySeqRecord(handle, 0,
                                 alphabet = self.ab)
        self.assertEqual(0, a._index_begin)
        self.assertEqual(12, a._index_end)
        self.assertEqual(None, a._seq)
        #self.assertEqual("sequencefake", a._handle)

    def test_seq_reader(self):
        """checks on basic worksing of _read_seq"""
        handle = "fakeid\nsequencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        handle.name = "fake"
        a = MinimalLazySeqRecord(handle, 0,
                                   alphabet = self.ab)
        self.assertEqual(None, a._seq)
        self.assertEqual("sequencefake", str(a.seq))
        self.assertEqual("sequencefake", str(a._seq))

    def test_simple_index_grab(self):
        """Test simple index getting"""
        handle = "fakeid\nsequencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        handle.name = "fake"
        a = MinimalLazySeqRecord(handle, 0,
                                   alphabet = self.ab)
        b = a[1:9]
        self.assertEqual(1, b._index_begin)
        self.assertEqual(9, b._index_end)

    def test_simple_index_grab_handle_identity(self):
        """Test simple index slicing, handle identity"""
        handle = "fakeid\nsequencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        handle.name = "fake"
        a = MinimalLazySeqRecord(handle, 0,
                                 alphabet = self.ab)
        b = a[1:9]
        self.assertTrue(a._handle is b._handle)
        self.assertTrue(a is not b)

    def test_two_level_getter_indexes(self):
        """Test private index position has changed"""
        handle = "fakeid\nsequencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        handle.name = "fake"
        a = MinimalLazySeqRecord(handle, 0,
                                   alphabet = self.ab)
        b = a[1:9]
        c = b[1:8]
        self.assertEqual(2, c._index_begin)
        self.assertEqual(9, c._index_end)

    def test_two_level_getter_handle_identity(self):
        """Test multi level getting and handle identity"""
        handle = "fakeid\nsequencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        handle.name = "fake"
        a = MinimalLazySeqRecord(handle, 0,
                                   alphabet = self.ab)
        b = a[1:9]
        c = b[1:8]
        self.assertTrue(a._handle is c._handle)

    def test_upper(self):
        """Test upper function"""
        handle = "fakeid\nseQUencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        handle.name = "fake"
        a = MinimalLazySeqRecord(handle, 0,
                                   alphabet = self.ab)
        b = a.upper()
        self.assertEqual("SEQUENCEFAKE", str(b._seq))

    def test_lower(self):
        """Test lower function"""
        handle = "fakeid\nseQUEncefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        handle.name = "fake"
        a = MinimalLazySeqRecord(handle, 0,
                                   alphabet = self.ab)
        b = a.lower()
        self.assertEqual("sequencefake", str(b._seq))

    def test_repr(self):
        """Test modified lazy repr function"""
        out = r"""MinimalLazySeqRecord(seq=NOT_READ, id=fakeid, """ +\
              r"""name=<unknown name>, description=<unknown description>,""" +\
              r""" dbxrefs=[])"""
        out2 = r"""MinimalLazySeqRecord(seq=Seq('sequencefake', """ +\
               r"""SingleLetterAlphabet()), id=fakeid, """ +\
               r"""name=<unknown name>, description=<unknown description>""" +\
               r""", dbxrefs=[])"""
        handle = "fakeid\nsequencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        handle.name = "fake"
        a = MinimalLazySeqRecord(handle, 0,
                                   alphabet = self.ab)
        self.assertEqual(out, repr(a))
        s = a.seq
        self.assertEqual(out2, repr(a))

    def test_len_method(self):
        """Test lazy record length method"""
        handle = "fakeid\nsequencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        handle.name = "fake"
        a = MinimalLazySeqRecord(handle, 0,
                                   alphabet = self.ab)
        b = a[3:6]
        self.assertEqual(3, len(b))
        self.assertEqual(None, b._seq)
        self.assertEqual(len("sequencefake"), len(a))

class TestFeatureBinCollection(unittest.TestCase):
    def setUp(self):
        self.bins = FeatureBinCollection(bounded_only_returns = False)

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
        self.assertTrue(256 not in self.bins._size_list)
        self.assertTrue(2048 in self.bins._size_list)
        self.assertEqual(k_overFull, 0)
        self.assertEqual(k_full, 1)
        self.assertEqual(self.bins._max_bin_power, 26)

    def test_binning_size_changer_multiple_steps(self):
        k_overFull = self.bins._calculate_bin_index(1,8388608*8)
        self.assertEqual(self.bins._max_sequence_length, 8388608*8*8)
        k_full = self.bins._calculate_bin_index(0,8388608*8)
        self.assertEqual(k_overFull, 0)
        self.assertEqual(k_full, 1)
        self.assertEqual(self.bins._max_bin_power, 29)
        self.assertTrue(2048 not in self.bins._size_list)

    def test_insertion_bad_negative_val(self):
        testTuple1 = (-1, 56)
        self.assertRaises(AssertionError,self.bins.insert,testTuple1)

    def test_insertion_once_smallest(self):
        testTuple1 = (0, 256)
        self.bins.insert(testTuple1)
        self.assertTrue(testTuple1 in self.bins._bins[4681])

    def test_insertion_once_smallest_but_overlaps(self):
        testTuple2 = (1, 257)
        self.bins.insert(testTuple2)
        self.assertTrue(testTuple2 in self.bins._bins[585])

    def test_insertion_once_smallestbin_rightmost(self):
        testTuple3 = (8388608-256, 8388608)
        self.bins.insert(testTuple3)
        self.assertTrue(testTuple3 in self.bins._bins[37448])

    def test_insertion_zero_length_smallestbin_rightmost(self):
        testTuple3 = (8388607, 8388607)
        self.bins.insert(testTuple3)
        self.assertTrue(testTuple3 in self.bins._bins[37448])

    def test_insertion_once_smallestbin_rightmost_lv4(self):
        testTuple4 = (8388608-257, 8388608)
        self.bins.insert(testTuple4)
        self.assertTrue(testTuple4 in self.bins._bins[4680])

    def test_insertion_with_rearrangement(self):
        testTuple3 = (256, 256+256)
        self.bins.insert(testTuple3)
        self.assertTrue(testTuple3 in self.bins._bins[4682])
        #trigger a size rearrangement with an insertion
        testTuple5 = (0, 9000000)
        self.bins.insert(testTuple5)
        self.assertTrue(testTuple3 in self.bins._bins[4681])
        self.assertTrue(testTuple5 in self.bins._bins[0])
        self.assertEqual(self.bins._max_sequence_length, 8388608*8)
        self.assertTrue(256 not in self.bins._size_list)
        self.assertTrue(2048 in self.bins._size_list)
        self.assertEqual(self.bins._max_bin_power, 26)

    def test_insertion_level3_and_level5_then_resize_to_fit_in_same_bin(self):
        test_tuple_lv3 = (16384 , 16384+16384)
        test_tuple_lv5 = (16384+16384-256 , 16384+16384)
        self.bins.insert(test_tuple_lv3)
        self.bins.insert(test_tuple_lv5)
        self.assertTrue(test_tuple_lv3 in self.bins._bins[74])
        self.assertTrue(test_tuple_lv5 in self.bins._bins[4681+127])
        #trigger rearrangement by 2 levels
        k_overFull = self.bins._calculate_bin_index(1,8388608*8)
        self.assertTrue(test_tuple_lv3 in self.bins._bins[4682])
        self.assertTrue(test_tuple_lv5 in self.bins._bins[4682])
        self.assertEqual(self.bins._max_bin_power, 29)

    def test_overflows_of_static_defined_lists(self):
        staticbins = FeatureBinCollection(length=67108864)
        #insert a chunk of data
        testTuple1 = (1, 2049)
        staticbins.insert(testTuple1)
        self.assertTrue(testTuple1 in staticbins._bins[585])
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
        self.assertTrue(resultR in emptyR)
        self.assertTrue(resultR in emptyM)

    def test_getter_half_slices_with_result_left_border(self):
        resultL = (10000,20000)
        self.bins.insert(resultL)
        emptyL = self.bins[:20000]
        emptyR = self.bins[20000:]
        emptyM = self.bins[:]
        self.assertEqual([], emptyR)
        self.assertTrue(resultL in emptyL)
        self.assertTrue(resultL in emptyM)

    def test_insertion_where_medium_sized_bin_is_out_of_bounds(self):
        resultL = (8388608, 8388608+2047)
        self.bins.insert(resultL)
        self.assertTrue(resultL in self.bins[838840:8388610])
        self.assertEqual(self.bins._max_bin_power, 26)

    def test_getter_overlap_left(self):
        feature = (100000,200000)
        self.bins.insert(feature)
        result = self.bins[99000:101000]
        self.assertTrue(feature in result)

    def test_getter_overlap_right(self):
        feature = (100000,200000)
        self.bins.insert(feature)
        result = self.bins[199000:201000]
        self.assertTrue(feature in result)

    def test_getter_feature_inside_region(self):
        feature = (100000,200000)
        self.bins.insert(feature)
        result = self.bins[99000:201000]
        self.assertTrue(feature in result)

    def test_getter_region_inside_feature(self):
        feature = (100000,200000)
        self.bins.insert(feature)
        result = self.bins[101000:199000]
        self.assertTrue(feature in result)

    def test_getter_zero_length_inside(self):
        testTuple3 = (8388605, 8388605)
        self.bins.insert(testTuple3)
        self.assertTrue(testTuple3 in self.bins[8388604:8388606])
        self.assertTrue(testTuple3 in self.bins[8388604:8388605])
        self.assertTrue(testTuple3 in self.bins[8388605:8388606])
        self.assertEqual([], self.bins[8388603:8388604])

    def test_getter_bounded_only_returns(self):
        testTuple3 = (200, 400)
        self.bins.insert(testTuple3)
        #test non-promiscuous feature returns
        self.bins.bounded_only_returns = True
        self.assertEqual([testTuple3], self.bins[200:400])
        self.assertEqual([testTuple3], self.bins[199:410])
        self.assertEqual([], self.bins[200:399])
        self.assertEqual([], self.bins[201:400])
        #sanity check: can promiscuous returns be enforced
        self.bins.bounded_only_returns = False
        self.assertEqual([testTuple3], self.bins[200:400])
        self.assertEqual([testTuple3], self.bins[199:410])
        self.assertEqual([testTuple3], self.bins[200:399])
        self.assertEqual([testTuple3], self.bins[201:400])

#
## Test the IndexManager to confirm valid DB IO operations
#

class TestIndexManager(unittest.TestCase):

    fakeindex = {"id":"new",
                 "recordoffsetstart":180,
                 "sequencestart":200,
                 "recordoffsetstart":200,
                 "recordoffsetlength":300,
                 "nextrecordoffset":301,
                 "sequencelinewidth":10,
                 "sequenceletterwidth":8,
                 "seqlen":80}

    # format: (seqbegin, seqend, offetbegin, offsetend, qualifier)
    fakefeatures = [(10, 80, 210, 220, "sequence"),
                    (9, 10, 225, 250, "site"),
                    (30, 50, 260, 270, "sometext"),
                    (20, 80, 280, 300, "sometext")]

    def test_index_manager_works_without_db(self):
        oshandlelazy, db_name = tempfile.mkstemp()
        os.close(oshandlelazy)
        newindex = SeqProxyIndexManager("fasta", indexdb=None,
                                        recordkey=None, handlename=None)
        newindex.set_record_index(self.fakeindex)
        self.assertEqual(self.fakeindex, newindex.record)
        del newindex
        os.remove(db_name)

    def test_index_manager_works_with_db(self):
        oshandlelazy, db_name = tempfile.mkstemp()
        os.close(oshandlelazy)
        newindex = SeqProxyIndexManager("fasta", indexdb=db_name,
                                        recordkey=None, handlename="lame")
        newindex.set_record_index(self.fakeindex)
        self.assertEqual(self.fakeindex, newindex.record)
        del newindex
        os.remove(db_name)

    def test_index_manager_gets_from_set_db(self):
        oshandlelazy, db_name = tempfile.mkstemp()
        os.close(oshandlelazy)
        newindex = SeqProxyIndexManager("fasta", indexdb=db_name,
                                        recordkey=None, handlename="lame")
        newindex.set_record_index(self.fakeindex)
        #cleanup old manager
        del newindex
        # make new manager
        newindex2 = SeqProxyIndexManager("fasta", indexdb=db_name,
                                         recordkey="new", handlename="lame")
        self.assertEqual(self.fakeindex, newindex2.record)
        #use new manager to double assign the same value
        fxwillraise = lambda idx: newindex2.set_record_index(idx)
        self.assertRaises(ValueError, fxwillraise, self.fakeindex)
        os.remove(db_name)

    def test_feature_indexing_with_db(self):
        oshandlelazy, db_name = tempfile.mkstemp()
        os.close(oshandlelazy)
        newindex = SeqProxyIndexManager("fasta", indexdb=db_name,
                                        recordkey=None, handlename="lame")
        featurebin = FeatureBinCollection()
        for ft in self.fakefeatures:
            featurebin.insert(ft)
        newindex.set_record_index(self.fakeindex)
        newindex.set_feature_index(featurebin)

        recoveredFeatures = newindex.get_features(0,80)[:]
        self.assertEqual(recoveredFeatures[:], featurebin[:])
        del newindex
        os.remove(db_name)

    def test_feature_indexing_with_fetch_from_db(self):
        oshandlelazy, db_name = tempfile.mkstemp()
        os.close(oshandlelazy)
        newindex = SeqProxyIndexManager("fasta", indexdb=db_name,
                                        recordkey=None, handlename="lame")
        featurebin = FeatureBinCollection()
        for ft in self.fakefeatures:
            featurebin.insert(ft)
        newindex.set_record_index(self.fakeindex)
        newindex.set_feature_index(featurebin)

        recoveredFeatures = newindex.features[:]
        self.assertEqual(recoveredFeatures[:], featurebin[:])
        del newindex
        os.remove(db_name)

    def test_feature_indexing_with_fetch_from_db(self):
        oshandlelazy, db_name = tempfile.mkstemp()
        os.close(oshandlelazy)
        #setting in DB
        newindex = SeqProxyIndexManager("fasta", indexdb=db_name,
                                        recordkey=None, handlename="lame")
        featurebin = FeatureBinCollection()
        for ft in self.fakefeatures:
            featurebin.insert(ft)
        newindex.set_record_index(self.fakeindex)
        newindex.set_feature_index(featurebin)
        del newindex
        #making new index manager
        newindex2 = SeqProxyIndexManager("fasta", indexdb=db_name,
                                        recordkey="new", handlename="lame")
        #get features again
        newfeatures = newindex2.get_features(1,80)
        self.assertEqual(newfeatures[:], featurebin[:])
        self.assertTrue(newindex2._features is None)
        #try full feature get operation (should set index._features)
        newfeatures = newindex2.get_features(0,80)
        self.assertTrue(isinstance(newindex2._features, FeatureBinCollection))
        os.remove(db_name)

#
## Test the SeqIO.parse and SeqIO.index_db bindings
#

class TestSeqIOBindings(unittest.TestCase):

    def setUp(self):
        oshandlelazy, db_name = tempfile.mkstemp()
        os.close(oshandlelazy)
        self.dbfile = db_name
        self.returnclass = SeqIO.FastaIO.FastaSeqRecProxy
        #set files for later tests
        self.onefile = os.path.join("Fasta", "f002")
        duplicates = ["f002", "f002_duplicate_for_test.fasta"]
        self.duplicatefiles = [os.path.join("Fasta", f) for f in duplicates]

    def tearDown(self):
        os.remove(self.dbfile)

    def test_bad_key_function(self):
        """test that an aggressive key function throws an error"""
        dbf = self.dbfile
        fasta = self.onefile
        key_fx = lambda x: "returnval"
        indexer = lambda f: SeqIO.index_db(dbf, fasta, format="fasta", \
                                           key_function=f, lazy=True)
        self.assertRaises(KeyError, indexer, key_fx)

    def test_bad_file_set(self):
        """test that an aggressive key function throws an error"""
        dbf = self.dbfile
        fastas = self.duplicatefiles
        indexer = lambda f: SeqIO.index_db(dbf, f, format="fasta", lazy=True)
        self.assertRaises(KeyError, indexer, fastas)

    def test_seqio_parse_nodb(self):
        """test that an aggressive key function throws an error"""
        fasta = self.onefile
        key_fx = lambda x: "returnval"
        firstten = ["CGGACCAGAC",
                    "CGGAGCCAGC",
                    "GATCAAATCT"]
        parseriter = SeqIO.parse(fasta, "fasta", lazy=True)
        for index, record in enumerate(parseriter):
            self.assertEqual(str(record[:10].seq), firstten[index])

    def test_seqio_parse_get_raw(self):
        """test that an aggressive key function throws an error"""
        fasta = self.onefile
        parseriter = SeqIO.parse(fasta, "fasta", lazy=True)
        firstrec = next(parseriter)
        raw = _bytes_to_string(firstrec.get_raw())
        fakefile = StringIO(raw)
        fakerec = SeqIO.read(fakefile, 'fasta')
        self.assertEqual(str(firstrec.seq), str(fakerec.seq))
        self.assertEqual(str(firstrec[0:10].seq), "CGGACCAGAC")
        self.assertEqual(firstrec.id, fakerec.id)

    def test_seqio_indexdb_get_raw(self):
        """test that an aggressive key function throws an error"""
        fasta = self.onefile
        dbf = self.dbfile
        def key_fx(key):
            newkey = key.split("|")
            return newkey[3].strip()
        lazydict = SeqIO.index_db(dbf, fasta, "fasta",
                                  key_function=key_fx, lazy=True)
        raw = _bytes_to_string(lazydict.get_raw("G26685"))
        fakefile = StringIO(raw)
        fakerec = SeqIO.read(fakefile, 'fasta')
        firstrec = lazydict["G26685"]
        self.assertEqual(str(firstrec.seq), str(fakerec.seq))
        self.assertEqual(str(firstrec[0:10].seq), "CGGAGCCAGC")
        self.assertEqual(firstrec.id, fakerec.id)

    def test_seqio_parse_withdb(self):
        """test that an aggressive key function throws an error"""
        dbf = self.dbfile
        fasta = self.onefile
        key_fx = lambda x: "returnval"
        firstten = ["CGGACCAGAC",
                    "CGGAGCCAGC",
                    "GATCAAATCT"]
        parseriter = SeqIO.parse(fasta, "fasta", lazy=dbf)
        for index, record in enumerate(parseriter):
            self.assertEqual(str(record[:10].seq), firstten[index])

    def test_good_key_function(self):
        """test that an nice key function works as expected"""
        dbf = self.dbfile
        fasta = self.onefile
        def key_fx(key):
            newkey = key.split("|")
            return newkey[3].strip()
        index = SeqIO.index_db(dbf, fasta, format="fasta", \
                               key_function=key_fx, lazy=True)
        keys = list(index.keys())
        self.assertTrue(len(keys) == 3)
        self.assertTrue("G26680" in keys)
        self.assertTrue("G26685" in keys)
        self.assertTrue("G29385" in keys)
        onerec = index["G26685"]
        self.assertTrue(isinstance(onerec, self.returnclass))
        self.assertEqual(onerec.id, "gi|1348917|gb|G26685|G26685")
        self.assertEqual(str(onerec[0:10].seq), "CGGAGCCAGC")
        self.assertEqual(str(onerec[0:10].seq), "CGGAGCCAGC")

#a = MinimalLazySeqRecord("seQUencefake", "fakeid")
#unittest.main( exit=False )
if __name__ == "__main__":
    #Run the test cases
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
