import unittest

import sys
import imp
import os
import tempfile
from Bio import SeqIO
from Bio.Alphabet import single_letter_alphabet
from Bio._py3k import _bytes_to_string, _as_bytes, _is_int_or_long
from Bio.Seq import Seq
from Bio.SeqIO._lazy import *
from Bio._py3k import _bytes_to_string, basestring
from io import StringIO, BytesIO
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
        """The iterator must finish"""
        lazyiterator = self.parser(test_seq_win)
        for count, rec in enumerate(lazyiterator):
            pass
        self.assertEqual(count, 1)
        lazyiterator = self.parser(test_seq_unix)
        for count, rec in enumerate(lazyiterator):
            pass
        self.assertEqual(count, 1)
    
    def test_iterator_returns_a_lazy_record(self):
        """checks the iterator returns the correct record"""
        lazyiterator = self.parser(test_seq_win)
        first_seq = next(lazyiterator)
        self.assertTrue(isinstance(first_seq,SeqIO.FastaIO.FastaSeqRecProxy))

    def test_indexing_simple_record_win(self):
        lazyiterator = self.parser(test_seq_win)
        first_seq = next(lazyiterator)
        idx = first_seq._index
        self.assertEqual(idx["sequenceletterwidth"], 60)
        self.assertEqual(idx["sequencelinewidth"], 62)
        self.assertEqual(idx["recordoffsetstart"], 0)

    def test_indexing_2nd_record_unix(self):
        lazyiterator = self.parser(test_seq_unix)
        second_seq = next(lazyiterator)
        second_seq = next(lazyiterator)
        idx = second_seq._index
        self.assertEqual(idx["sequenceletterwidth"], 51)
        self.assertEqual(idx["sequencelinewidth"], 52)
        self.assertEqual(idx["recordoffsetstart"], 290)
    
    def test_len_win(self):
        lazyiterator = self.parser(test_seq_win)
        first_seq = next(lazyiterator)
        self.assertEqual(first_seq._index["seqlen"], 225)

    def test_indexing_2nd_record_win(self):
        lazyiterator = self.parser(test_seq_win)
        second_seq = next(lazyiterator)
        second_seq = next(lazyiterator)
        idx = second_seq._index
        self.assertEqual(idx["sequenceletterwidth"], 51)
        self.assertEqual(idx["sequencelinewidth"], 53)
        self.assertEqual(idx["recordoffsetstart"], 295)

    def test_sequence_getter_zero_unix(self):
        lazyiterator = self.parser(test_seq_unix)
        firstseq = next(lazyiterator)
        s = firstseq[0:10]
        self.assertEqual(str(s.seq), "MAPNASCLCV")

    def test_sequence_getter_zero_win(self):
        lazyiterator = self.parser(test_seq_win)
        firstseq = next(lazyiterator)
        s = firstseq[0:10]
        self.assertEqual(str(s.seq), "MAPNASCLCV")
   
    def test_double_padded_seq_len(self):
        lazyiterator = self.parser(test_double_padded)
        firstseq = next(lazyiterator)
        self.assertEqual(len(firstseq), 40)

    def test_double_padded_padding_difference(self):
        lazyiterator = self.parser(test_double_padded)
        firstseq = next(lazyiterator)
        secondseq = next(lazyiterator)
        firstseqi =firstseq._index
        firstseqend = firstseqi["recordoffsetstart"] + firstseqi["recordoffsetlength"]
        padding = secondseq._index["recordoffsetstart"] - firstseqend
        self.assertEqual(2, padding)

    def test_one_bad_sequence_line(self):
        lazyiterator = self.parser(test_bad_line)
        firstseq = next(lazyiterator)
        #self.assertEqual(45, len(firstseq))
        #self.assertEqual("DEMANS", str(firstseq[24:30].seq))
        self.assertRaises(ValueError, firstseq._read_seq)

    def test_end_iteration(self):
        lazyiterator = self.parser(test_seq_unix)
        firstseq = next(lazyiterator)
        secondseq = next(lazyiterator)
        self.assertRaises(StopIteration, next, lazyiterator)

    def test_one_bad_sequence_header(self):
        lazyiterator = self.parser(test_bad_header)
        firstseq = next(lazyiterator)
        #self.assertEqual(45, len(firstseq))
        #self.assertEqual("DEMANS", str(firstseq[24:30].seq))
        self.assertRaises(ValueError, firstseq._read_seq)

    def test_sequence_getter_unix(self):
        lazyiterator = self.parser(test_seq_unix)
        firstseq = next(lazyiterator)
        s = firstseq[6:10]
        self.assertEqual(str(s.seq), "CLCV")   

    def test_sequence_getter_win(self):
        lazyiterator = self.parser(test_seq_win)
        firstseq = next(lazyiterator)
        s = firstseq[6:10]
        self.assertEqual(str(s.seq), "CLCV")  

    def test_sequence_getter_unix_2_line_span1(self):
        lazyiterator = self.parser(test_seq_unix)
        firstseq = next(lazyiterator)
        s = firstseq[59:62]
        self.assertEqual(str(s.seq), "RRS")

    def test_sequence_getter_unix_2_line_span2(self):
        lazyiterator = self.parser(test_seq_gaped)
        firstseq = next(lazyiterator)
        secondseq = next(lazyiterator)
        s = secondseq[40:43]
        self.assertEqual(str(s.seq), "RGY")

    def test_sequence_getter_unix_2_line_span3(self):
        lazyiterator = self.parser(test_seq_win)
        firstseq = next(lazyiterator)
        s = firstseq[58:62]
        self.assertEqual(str(s.seq), tsuseq[58:62])

    def test_sequence_getter_win_2_line_span(self):
        lazyiterator = self.parser(test_seq_win)
        firstseq = next(lazyiterator)
        s = firstseq[59:61]
        self.assertEqual(str(s.seq), tsuseq[59:61])

    def test_sequence_getter_full_span(self):
        lazyiterator = self.parser(test_seq_unix)
        firstseq = next(lazyiterator)
        s = firstseq
        self.assertEqual(str(s.seq), tsuseq)
    
    def test_sequence_getter_explicit_full_span(self):
        lazyiterator = self.parser(test_seq_unix)
        firstseq = next(lazyiterator)
        s = firstseq[3:225]
        self.assertEqual(str(s.seq), tsuseq[3:225])

    def test_sequence_getter_inset_full_span(self):
        lazyiterator = self.parser(test_seq_unix)
        firstseq = next(lazyiterator)
        s = firstseq[0:161]
        self.assertEqual(len(tsuseq[0:161]), 161)
        self.assertEqual(len(s), 161)
        self.assertEqual(len(s.seq), 161)
        self.assertEqual(str(s.seq), tsuseq[0:161]) 

class LazyFastaIOSimpleTestsGenericIterator(LazyFastaIOSimpleTests):
    """Tests the generic iterator using tests from format specific iterator"""

    def setUp(self):
        returncls = SeqIO.FastaIO.FastaSeqRecProxy
        self.parser = lambda handle: SeqIO._lazy.lazy_iterator(handle, \
                                                returncls, 'fasta')

#
### tests for Fasta comparison
#

class TestFastaSeqRecord(unittest.TestCase):
    fastafile = "Fasta/f002"

    def setUp(self):
        self.standard_iter = SeqIO.parse(self.fastafile, 'fasta')
        self.file = open(self.fastafile, 'rb')
        lazyiter = SeqIO.FastaIO.FastaLazyIterator    
        self.lazy_iter = lazyiter(self.file)

    def tearDown(self):
        self.file.close()

    def test_same_length_and_end_behavior(self):
        for std, lzy in zip(self.standard_iter, self.lazy_iter):
            pass
        self.assertRaises(StopIteration, next, self.standard_iter)
        self.assertRaises(StopIteration, next, self.lazy_iter)

    def test_same_id(self):
        for std, lzy in zip(self.standard_iter, self.lazy_iter):
            self.assertEqual(std.id, lzy.id)

    def test_same_seq(self):
        for std, lzy in zip(self.standard_iter, self.lazy_iter):
            self.assertEqual(str(std.seq), str(lzy.seq))

    def test_same_name_description(self):
        for std, lzy in zip(self.standard_iter, self.lazy_iter):
            self.assertEqual(std.description, lzy.description)
            self.assertEqual(std.name, lzy.name)

    def test_slicing_behavior(self):
        for std, lzy in zip(self.standard_iter, self.lazy_iter):
            S1, L1 = std[50:], lzy[50:]
            S2, L2 = std[:50], lzy[:50]
            S3, L3 = std[261:380], lzy[261:380]
            self.assertEqual(str(S1.seq), str(L1.seq))
            self.assertEqual(str(S2.seq), str(L2.seq))
            self.assertEqual(str(S3.seq), str(L3.seq))

class TestFastaFromIndexDb( TestFastaSeqRecord):
    """Test iterating using a premade index database compare to standard"""
    fastafile = "Fasta/f002"
    def setUp(self):
        #db setup: make tempfile and close os handle to help with cleanup
        oshandle, self.dbfilename = tempfile.mkstemp()
        os.close(oshandle)
        #iter setup lazy
        self.file = open(self.fastafile, 'rb')
        returncls = SeqIO.FastaIO.FastaSeqRecProxy
        lazyiter = lambda handle: SeqIO._lazy.lazy_iterator(handle, \
                                        returncls, 'fasta', \
                                        index=self.dbfilename)
        self.lazy_iter = lazyiter(self.file)
        #iter setup regular
        self.standard_iter = SeqIO.parse(self.fastafile, 'fasta') 
        
    def tearDown(self):
        #delete db temp file
        os.remove(self.dbfilename)
        #close sequence file
        self.file.close()
    
    def test_setup_cleanup(self):
        for std, lzy in zip(self.standard_iter, self.lazy_iter):
            pass
        self.assertRaises(StopIteration, next, self.standard_iter)
        self.assertRaises(StopIteration, next, self.lazy_iter)

class TestFastaAlignAgainst(TestFastaSeqRecord):
    fastafile = "Fasta/fa01"

#
### tests for GenBank IO
#

from Bio.SeqIO.InsdcIO import GenbankSeqRecProxy

class TestGenbankLazy(unittest.TestCase):
    recordfile = "brca_FJ940752.gb" 

    def setUp(self):
        returncls = GenbankSeqRecProxy
        self.parser = lambda handle: SeqIO._lazy.lazy_iterator(handle, \
                                                returncls, 'genbank')
        self.handle.seek(0)

    @classmethod
    def setUpClass(cls):
        cls.handle = open(os.path.join('GenBank', cls.recordfile), 'rb')
    
    @classmethod
    def tearDownClass(cls):
        cls.handle.close()

    def test_parser_init(self):
        recordgen = self.parser(self.handle)
        record = next(recordgen)

    def test_id_name(self):
        record = self.parser(self.handle)
        record = next(record)
        self.assertEqual(record.id, 'FJ940752.1')
        self.assertEqual(record.name, 'FJ940752')

    def test_record_description(self):
        record = self.parser(self.handle)
        record = next(record)
        descr = "Homo sapiens BRCA1 (BRCA1) gene, exon 18 and partial cds."
        self.assertEqual(record.description, descr)

    def test_id_seq(self):
        record = self.parser(self.handle)
        record = next(record)
        self.assertEqual(str(record[0:5].seq), 'GGCTC')
        self.assertEqual(str(record[-5:].seq), 'GTCTC')
        self.assertEqual(str(record[70:75].seq), "TTCTG")
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
        self._index = new_index

    def _read_seq(self):
        self._handle.seek(self._index["sequencestart"])
        seqstr = _bytes_to_string(self._handle.readline())
        self._seq = Seq(seqstr, self._alphabet)

    def _load_non_lazy_values(self):
        self.id = self._index["id"]
        self.name = "<unknown name>"
        self.description = "<unknown description>"
        self.dbxrefs = []

    def _make_feature_index(self, new_list):
        pass

    def _read_features(self):
        return []


class SeqRecordProxyBaseClassTests(unittest.TestCase):

    def setUp(self):
        self.ab = single_letter_alphabet
        pass

    def test_simple_tester(self):
        """checks on basic definitions in tester class"""
        handle = "fakeid\nsequencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        a = MinimalLazySeqRecord(handle, 0, lenhandle,
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
        a = MinimalLazySeqRecord(handle, 0, lenhandle,
                                   alphabet = self.ab)
        self.assertEqual(None, a._seq)
        self.assertEqual("sequencefake", str(a.seq))
        self.assertEqual("sequencefake", str(a._seq))

    def test_simple_index_grab(self):
        handle = "fakeid\nsequencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        a = MinimalLazySeqRecord(handle, 0, lenhandle,
                                   alphabet = self.ab)
        b = a[1:9]
        self.assertEqual(1, b._index_begin)
        self.assertEqual(9, b._index_end)

    def test_simple_index_grab_handle_identity(self):
        handle = "fakeid\nsequencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        a = MinimalLazySeqRecord(handle, 0, lenhandle,
                                   alphabet = self.ab)
        b = a[1:9]
        self.assertTrue(a._handle is b._handle)
        self.assertTrue(a is not b)

    def test_two_level_getter_indexes(self):
        handle = "fakeid\nsequencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        a = MinimalLazySeqRecord(handle, 0, lenhandle,
                                   alphabet = self.ab)
        b = a[1:9]
        c = b[1:8]
        self.assertEqual(2, c._index_begin)
        self.assertEqual(9, c._index_end)

    def test_two_level_getter_handle_identity(self):
        handle = "fakeid\nsequencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        a = MinimalLazySeqRecord(handle, 0, lenhandle,
                                   alphabet = self.ab)
        b = a[1:9]
        c = b[1:8]
        self.assertTrue(a._handle is c._handle)

    def test_upper(self):
        handle = "fakeid\nseQUencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        a = MinimalLazySeqRecord(handle, 0, lenhandle,
                                   alphabet = self.ab)
        b = a.upper()
        self.assertEqual("SEQUENCEFAKE", str(b._seq))

    def test_lower(self):
        handle = "fakeid\nseQUEncefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        a = MinimalLazySeqRecord(handle, 0, lenhandle,
                                   alphabet = self.ab)
        b = a.lower()
        self.assertEqual("sequencefake", str(b._seq))

    def test_repr(self):
        out = r"""MinimalLazySeqRecord(seq=NOT_READ, id=fakeid, """ + \
              r"""name=<unknown name>, description=<unknown description>,"""+\
              r""" dbxrefs=[])"""
        out2 = r"""MinimalLazySeqRecord(seq=Seq('sequencefake', """ +\
               r"""SingleLetterAlphabet()), id=fakeid, """ + \
               r"""name=<unknown name>, description=<unknown description>"""+\
               r""", dbxrefs=[])"""
        handle = "fakeid\nsequencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        a = MinimalLazySeqRecord(handle, 0, lenhandle,
                                   alphabet = self.ab)
        self.assertEqual(out, repr(a))
        s = a.seq 
        self.assertEqual(out2, repr(a))

    def test_len_method(self):
        handle = "fakeid\nsequencefake"
        lenhandle = len(handle)
        handle = BytesIO(_as_bytes(handle))
        a = MinimalLazySeqRecord(handle, 0, lenhandle,
                                   alphabet = self.ab)
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
        self.assertTrue(256 not in self.bins._size_list)
        self.assertTrue(2048 in self.bins._size_list)
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
        self.assertEqual([],self.bins[8388603:8388604])



File = namedtuple('File', 'name')
testfile = File("fake.fasta")

class IndexDbIO(SeqRecordProxyBase):
    _handle = testfile
    _indexkey = None
    _format = 'fake'
    def __init__(self):
        pass
        
class IndexDbIOCreateDb(unittest.TestCase):

    def setUp(self):
        #db setup: make tempfile and close os handle to help with cleanup
        oshandle, self.dbfilename = tempfile.mkstemp()
        os.close(oshandle)
        
    def tearDown(self):
        os.remove(self.dbfilename)
        

    def test_index_getter_empty_database(self):
        """empty databases follwed by a get call should raise KeyError"""
        testIO = IndexDbIO()
        testIO._indexkey = "test"
        testIO._indexdb = self.dbfilename
        gettermethod = testIO._SeqRecordProxyBase__load_and_return_index
        self.assertRaises(KeyError, gettermethod)
            
    def test_index_creates_table_without_raising(self):
        testindex = {"recordoffsetstart":0, 
                 "recordoffsetlength":200,
                 "seqstart":43, "id":"test1"}
        testIO = IndexDbIO()
        testIO._indexdb = self.dbfilename
        testIO._handle = testfile
        con = testIO._SeqRecordProxyBase__connect_db_return_connection()
        testIO._SeqRecordProxyBase__create_tables(con, testindex)
        
        
    def test_index_creates_table_and_check(self):
        """tests creation of table and class that checks tables"""
        #index setup
        testindex = {"recordoffsetstart":0, 
                 "recordoffsetlength":200,
                 "seqstart":43, "id":"test1"}
        testIO = IndexDbIO()
        testIO._indexdb = self.dbfilename
        testIO._handle = testfile
        con = testIO._SeqRecordProxyBase__connect_db_return_connection()
        cursor = con.cursor()
        #no
        a = testIO._SeqRecordProxyBase__is_valid_db_with_tables(con)
        self.assertFalse(a)
        #create tables
        testIO._SeqRecordProxyBase__create_tables(con, testindex)
        #make sure a valid sequence is found after creating tables
        a = testIO._SeqRecordProxyBase__is_valid_db_with_tables(con)
        self.assertTrue(a)
        one = cursor.execute("SELECT rowid, * FROM main_index")
        onedata = one.fetchone()
        colnames = [ d[0] for d in one.description ]
        colnames.sort()
        expected = ['fileid', 'id', 'recordoffsetlength', \
                    'recordoffsetstart', 'rowid', 'seqstart',]
        self.assertEqual(colnames, expected)
        
        
    
    def test_fileid_getter_method(self):
        testindex = {"recordoffsetstart":0, 
                 "recordoffsetlength":200,
                 "seqstart":43, "id":"test1"}
        testIO = IndexDbIO()
        testIO._indexdb = self.dbfilename
        con = testIO._SeqRecordProxyBase__connect_db_return_connection()
        #make the tables first
        testIO._SeqRecordProxyBase__create_tables(con, testindex)
        cursor = con.cursor()
        #invoke file id finder
        testIO._SeqRecordProxyBase__get_fileid(con)
        name = cursor.execute("SELECT filename " +\
                       "FROM indexed_files WHERE " +\
                       "filename=?;",('fake.fasta',))
        self.assertEqual(name.fetchone(), None)
        found_id = testIO._SeqRecordProxyBase__get_fileid(con, write=True)
        name = cursor.execute("SELECT filename, fileid " +\
                       "FROM indexed_files WHERE " +\
                       "filename=?;",('fake.fasta',))
        name, fileid = name.fetchone()
        self.assertEqual(name, 'fake.fasta')
        self.assertTrue(_is_int_or_long(fileid))
        self.assertEqual(found_id, fileid)
        
    def test_write_index_filename_is_written(self):
        testindex = {"recordoffsetstart":0, 
                 "recordoffsetlength":200,
                 "seqstart":43, "id":"test1"}
        testIO = IndexDbIO()
        testIO._indexdb = self.dbfilename
        #invoke writer
        testIO._index = testindex
        self.assertEqual(None, testIO._indexkey)
        #get cursor
        con = testIO._SeqRecordProxyBase__connect_db_return_connection()
        cursor = con.cursor()
        name = cursor.execute("SELECT filename " +\
                       "FROM indexed_files WHERE " +\
                       "filename=?;",('fake.fasta',))
        name = name.fetchone()[0]
        self.assertEqual(name, 'fake.fasta')
        
    def test_write_index_contents_arewritten(self):
        testindex = {"recordoffsetstart":0, 
                 "recordoffsetlength":200,
                 "seqstart":43, "id":"test1"}
        testIO = IndexDbIO()
        testIO._indexdb = self.dbfilename
        #invoke writer
        testIO._index = testindex
        #get cursor
        con = testIO._SeqRecordProxyBase__connect_db_return_connection()
        cursor = con.cursor()
        recordindex = cursor.execute("SELECT main_index.* " +\
                       "FROM main_index " +\
                       "INNER JOIN indexed_files " +\
                       "ON main_index.fileid = indexed_files.fileid " +\
                       "WHERE indexed_files.filename=?;",('fake.fasta',))
        record_keys = [key[0] for key in recordindex.description]
        record_index = recordindex.fetchone()
        record_indexdict = dict((record_keys[i], record_index[i]) for \
                                i in range(len(record_keys)))
        del(record_indexdict["fileid"])
        self.assertEqual(testindex, record_indexdict)
    
    def test_write_retrieve_with_shortcut_from_memory(self):
        testindex = {"recordoffsetstart":0, 
                 "recordoffsetlength":200,
                 "seqstart":43, "id":"test1"}
        testIO = IndexDbIO()
        testIO._indexdb = self.dbfilename
        #invoke writer
        testIO._index = testindex
        self.assertEqual(testIO._index, testindex)
        
    def test_write_retrieve_from_db(self):
        testindex = {"recordoffsetstart":0, 
                 "recordoffsetlength":200,
                 "seqstart":43, "id":"test1"}
        testIO = IndexDbIO()
        testIO._indexdb = self.dbfilename
        testIO._indexkey = None
        con = testIO._SeqRecordProxyBase__connect_db_return_connection()
        #invoke writer
        testIO._index = testindex
        testIO._indexkey = "test1"
        #reset memory index to None so that getter touches db
        testIO._SeqRecordProxyBase__index = None
        isvalid = testIO._SeqRecordProxyBase__is_valid_db_with_tables(con)
        self.assertTrue(isvalid)
        #activate getter and compare
        self.assertEqual(testIO._index, testindex)
        
    def test_write_two_records_overlapping(self):
        testindex = {"recordoffsetstart":0, 
                 "recordoffsetlength":200,
                 "seqstart":43, "id":"test1"}
        testIO = IndexDbIO()
        testIO._indexdb = self.dbfilename
        #invoke writer
        testIO._index = testindex
        testIO2 = IndexDbIO()
        testIO2._indexdb = self.dbfilename
        #invoke writer
        setter = testIO2._SeqRecordProxyBase__set_and_save_index
        self.assertRaises(ValueError, setter, testindex)
        
    def test_write_two_records_overlapping_id(self):
        testindex = {"recordoffsetstart":0, 
                 "recordoffsetlength":200,
                 "seqstart":43, "id":"test1"}
        testindex2 = {"recordoffsetstart":220, 
                 "recordoffsetlength":3000,
                 "seqstart":43, "id":"test1"}
        testIO = IndexDbIO()
        testIO._indexdb = self.dbfilename
        #invoke writer
        testIO._index = testindex
        testIO2 = IndexDbIO()
        testIO2._indexdb = self.dbfilename
        #invoke writer
        setter = testIO2._SeqRecordProxyBase__set_and_save_index
        self.assertRaises(ValueError, setter, testindex2)
    
    def test_write_two_records_distinct(self):
        testindex = {"recordoffsetstart":0, 
                 "recordoffsetlength":200,
                 "seqstart":43, "id":"test1"}
        testindex2 = {"recordoffsetstart":220, 
                 "recordoffsetlength":3000,
                 "seqstart":43, "id":"test2"}
        testIO = IndexDbIO()
        testIO._indexdb = self.dbfilename
        #invoke writer
        testIO._index = testindex
        testIO2 = IndexDbIO()
        testIO2._indexdb = self.dbfilename
        #invoke writer
        testIO._index = testindex2     
        
    def test_read_two_records_distinct_idkeys(self):
        testindex = {"recordoffsetstart":0, 
                 "recordoffsetlength":200,
                 "seqstart":43, "id":"test1"}
        testindex2 = {"recordoffsetstart":220, 
                 "recordoffsetlength":3000,
                 "seqstart":43, "id":"test2"}
        self.test_write_two_records_distinct()
        testIO = IndexDbIO()
        testIO._indexdb = self.dbfilename
        testIO._indexkey = "test1"
        self.assertEqual(testIO._index, testindex)
        testIO2 = IndexDbIO()
        testIO._indexdb = self.dbfilename
        testIO._indexkey = "test2"
        self.assertEqual(testIO._index, testindex)
        
    def test_read_two_records_distinct_intkeys(self):
        testindex = {"recordoffsetstart":0, 
                 "recordoffsetlength":200,
                 "seqstart":43, "id":"test1"}
        testindex2 = {"recordoffsetstart":220, 
                 "recordoffsetlength":3000,
                 "seqstart":43, "id":"test2"}
        self.test_write_two_records_distinct()
        testIO = IndexDbIO()
        testIO._indexdb = self.dbfilename
        testIO._indexkey = 0
        self.assertEqual(testIO._index, testindex)
        testIO2 = IndexDbIO()
        testIO._indexdb = self.dbfilename
        testIO._indexkey = 1
        self.assertEqual(testIO._index, testindex)

#a = MinimalLazySeqRecord("seQUencefake", "fakeid")
#unittest.main( exit=False )
if __name__ == "__main__":
    #Run the test cases
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
