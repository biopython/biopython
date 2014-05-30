import unittest
#from Bio import SeqRecord
#from Bio.SeqIO import _lazy
#from Bio.Seq import Seq
if __name__=="__main__":
    import imp
    import os
    from Bio.Seq import Seq
    from Bio.SeqIO._lazy import *


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
