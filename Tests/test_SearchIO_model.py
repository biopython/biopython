# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO objects.

Tests the methods and behaviors of QueryResult, Hit, and HSP objects. All tests
are format-independent and are meant to check the fundamental behavior common
to all formats.

"""

import pickle
import unittest
from io import BytesIO
from copy import deepcopy

from search_tests_common import compare_search_obj

from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import single_letter_alphabet, generic_dna
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# mock HSPFragments
frag111 = HSPFragment('hit1', 'query1', hit='ATGCGCAT', query='ATGCGCAT')
frag112 = HSPFragment('hit1', 'query1', hit='ATG', query='GAT')
frag113 = HSPFragment('hit1', 'query1', hit='ATTCG', query='AT-CG')
frag113b = HSPFragment('hit1', 'query1', hit='ATTCG', query='AT-CG')
frag114 = HSPFragment('hit1', 'query1', hit='AT', query='AT')
frag114b = HSPFragment('hit1', 'query1', hit='ATCG', query='ATGG')
frag211 = HSPFragment('hit2', 'query1', hit='GGGCCC', query='GGGCC-')
frag311 = HSPFragment('hit3', 'query1', hit='GATG', query='GTTG')
frag312 = HSPFragment('hit3', 'query1', hit='ATATAT', query='ATATAT')
frag411 = HSPFragment('hit4', 'query1', hit='CC-ATG', query='CCCATG')
frag121 = HSPFragment('hit1', 'query2', hit='GCGAG', query='GCGAC')
# mock HSPs
hsp111 = HSP([frag111])
hsp112 = HSP([frag112])
hsp113 = HSP([frag113, frag113b])
hsp114 = HSP([frag114, frag114b])
hsp211 = HSP([frag211])
hsp311 = HSP([frag311])
hsp312 = HSP([frag312])
hsp411 = HSP([frag411])
hsp121 = HSP([frag121])
# mock Hits
hit11 = Hit([hsp111, hsp112, hsp113, hsp114])
hit21 = Hit([hsp211])
hit31 = Hit([hsp311, hsp312])
hit41 = Hit([hsp411])
hit12 = Hit([hsp121])


class QueryResultCases(unittest.TestCase):

    def setUp(self):
        self.qresult = QueryResult([hit11, hit21, hit31], 'query1')
        # set mock attributes
        self.qresult.seq_len = 1102
        self.qresult.target = 'refseq_rna'

    def test_pickle(self):
        """Test pickling and unpickling of QueryResult"""
        buf = BytesIO()
        pickle.dump(self.qresult, buf)
        unp = pickle.loads(buf.getvalue())
        self.assertTrue(compare_search_obj(self.qresult, unp))

    def test_order(self):
        # added hits should be ordered
        self.assertEqual(self.qresult[0], hit11)
        self.assertEqual(self.qresult[2], hit31)
        # removal of second item should bump item #2 to #1
        del self.qresult['hit2']
        self.assertEqual(self.qresult[0], hit11)
        self.assertEqual(self.qresult[1], hit31)

    def test_init_none(self):
        """Test QueryResult.__init__, no arguments"""
        qresult = QueryResult()
        self.assertEqual(None, qresult.id)
        self.assertEqual(None, qresult.description)

    def test_init_id_only(self):
        """Test QueryResult.__init__, with ID only"""
        qresult = QueryResult(id='query1')
        self.assertEqual('query1', qresult.id)
        self.assertEqual(None, qresult.description)

    def test_init_hits_only(self):
        """Test QueryResult.__init__, with hits only"""
        qresult = QueryResult([hit11, hit21, hit31])
        self.assertEqual('query1', qresult.id)
        self.assertEqual('<unknown description>', qresult.description)

    def test_repr(self):
        """Test QueryResult.__repr__"""
        self.assertEqual("QueryResult(id='query1', 3 hits)",
                         repr(self.qresult))

    def test_iter(self):
        """Test QueryResult.__iter__"""
        # iteration should return hits contained
        for counter, hit in enumerate(self.qresult):
            self.assertIn(hit, (hit11, hit21, hit31))
        self.assertEqual(2, counter)

    def test_hits(self):
        """Test QueryResult.hits"""
        # hits should return hits contained in qresult
        hits = [x for x in self.qresult.hits]
        self.assertEqual([hit11, hit21, hit31], hits)

    def test_hit_keys(self):
        """Test QueryResult.hit_keys"""
        # hit_keys should return hit keys (which default to hit ids)
        hit_keys = [x for x in self.qresult.hit_keys]
        self.assertEqual(['hit1', 'hit2', 'hit3'], hit_keys)

    def test_items(self):
        """Test QueryResult.items"""
        # items should return tuples of hit key, hit object pair
        items = [x for x in self.qresult.items]
        self.assertEqual([('hit1', hit11), ('hit2', hit21),
                         ('hit3', hit31)], items)

    def test_hsps(self):
        """Test QueryResult.hsps"""
        # hsps should return all hsps contained in qresult
        hsps = self.qresult.hsps
        self.assertEqual([hsp111, hsp112, hsp113, hsp114, hsp211, hsp311,
                         hsp312, ], hsps)

    def test_fragments(self):
        """Test QueryResult.fragments"""
        # fragments should return all fragments contained in qresult
        frags = self.qresult.fragments
        self.assertEqual([frag111, frag112, frag113, frag113b, frag114,
                         frag114b, frag211, frag311, frag312], frags)

    def test_contains(self):
        """Test QueryResult.__contains__"""
        # contains should work with hit ids or hit objects
        self.assertIn('hit1', self.qresult)
        self.assertIn(hit21, self.qresult)
        self.assertFalse('hit5' in self.qresult)
        self.assertFalse(hit41 in self.qresult)

    def test_contains_alt(self):
        """Test QueryResult.__contains__, with alternative IDs"""
        # contains should work with alternative hit IDs
        hit11._id_alt = ['alt1']
        query = QueryResult([hit11])
        self.assertIn('alt1', query)
        hit11._id_alt = []

    def test_len(self):
        """Test QueryResult.__len__"""
        # len() should return the number of hits contained
        self.assertEqual(3, len(self.qresult))

    def test_nonzero(self):
        """Test QueryResult.__nonzero__"""
        # nonzero should return true only if the qresult has hits
        self.assertTrue(self.qresult)
        blank_qresult = QueryResult()
        self.assertFalse(blank_qresult)

    def test_setitem_ok(self):
        """Test QueryResult.__setitem__"""
        # hit objects assignment should work with arbitrary string keys
        self.qresult['hit4'] = hit41
        self.assertEqual([hit11, hit21, hit31, hit41], list(self.qresult.hits))
        # and if the key already exist, the object should be overwritten
        self.qresult['hit4'] = hit11
        self.assertEqual([hit11, hit21, hit31, hit11], list(self.qresult.hits))

    def test_setitem_ok_alt(self):
        """Test QueryResult.__setitem__, checking alt hit IDs"""
        # hit objects assignment should make alt IDs visible
        hit11._id_alt = ['alt1', 'alt11']
        query = QueryResult()
        query['hit1'] = hit11
        self.assertEqual(hit11, query['hit1'])
        self.assertEqual(hit11, query['alt1'])
        self.assertEqual(hit11, query['alt11'])
        self.assertTrue(hit11.id != 'alt1')
        self.assertTrue(hit11.id != 'alt11')
        hit11._id_alt = []

    def test_setitem_ok_alt_existing(self):
        """Test QueryResult.__setitem__, existing key"""
        # hit objects assignment on existing hits should also update alt IDs
        hit11._id_alt = ['alt1']
        hit21._id_alt = ['alt2']
        query = QueryResult()
        query['hit'] = hit11
        self.assertEqual(hit11, query['hit'])
        self.assertEqual(hit11, query['alt1'])
        query['hit'] = hit21
        self.assertEqual(hit21, query['hit'])
        self.assertEqual(hit21, query['alt2'])
        self.assertRaises(KeyError, query.__getitem__, 'alt1')
        hit11._id_alt = []
        hit21._id_alt = []

    def test_setitem_ok_alt_ok_promote(self):
        """Test QueryResult.__setitem__, previously alt ID"""
        # hit objects assignment with ID previously existing as alternative
        # should make the ID primary
        hit11._id_alt = ['alt1']
        hit41._id_alt = ['alt4']
        hit31._id_alt = ['alt3']
        query = QueryResult([hit11, hit41])
        self.assertEqual(hit11, query['alt1'])
        self.assertEqual(hit41, query['alt4'])
        self.assertNotIn('alt1', query._items)
        self.assertIn('alt1', query._QueryResult__alt_hit_ids)
        query['alt1'] = hit31
        self.assertEqual(hit31, query['alt1'])
        self.assertEqual(hit41, query['alt4'])
        self.assertIn('alt1', query._items)
        self.assertNotIn('alt1', query._QueryResult__alt_hit_ids)
        hit11._id_alt = []
        hit41._id_alt = []
        hit31._id_alt = []

    def test_setitem_wrong_key_type(self):
        """Test QueryResult.__setitem__, wrong key type"""
        # item assignment should fail if the key is not string
        self.assertRaises(TypeError, self.qresult.__setitem__, 0, hit41)
        self.assertRaises(TypeError, self.qresult.__setitem__, slice(0, 2),
                          [hit41, hit31])

    def test_setitem_wrong_type(self):
        """Test QueryResult.__setitem__, wrong type"""
        # item assignment should fail if the object assigned is not a hit object
        self.assertRaises(TypeError, self.qresult.__setitem__, 'hit4', hsp111)
        self.assertRaises(TypeError, self.qresult.__setitem__, 'hit5', 'hit5')

    def test_setitem_wrong_query_id(self):
        """Test QueryResult.__setitem__, wrong query ID"""
        # item assignment should fail if the hit object does not have the same
        # query id
        self.assertRaises(ValueError, self.qresult.__setitem__, 'hit4', hit12)

    def test_setitem_from_empty(self):
        """Test QueryResult.__setitem__, from empty container"""
        qresult = QueryResult()
        # initial desc and id is None
        self.assertEqual(None, qresult.id)
        self.assertEqual(None, qresult.description)
        # but changes to the first item's after append
        qresult.append(hit11)
        self.assertEqual('query1', qresult.id)
        self.assertEqual('<unknown description>', qresult.description)
        # and remains the same after popping the last item
        qresult.pop()
        self.assertEqual('query1', qresult.id)
        self.assertEqual('<unknown description>', qresult.description)

    def test_getitem_default_ok(self):
        """Test QueryResult.__getitem__"""
        # hits should be retrievable by their keys (default to id)
        self.assertEqual(hit21, self.qresult['hit2'])
        self.assertEqual(hit11, self.qresult['hit1'])

    def test_getitem_int_ok(self):
        """Test QueryResult.__getitem__, with integer"""
        # hits should be retrievable by their index
        self.assertEqual(hit21, self.qresult[1])
        self.assertEqual(hit31, self.qresult[-1])

    def test_getitem_slice_ok(self):
        """Test QueryResult.__getitem__, with slice"""
        # if the index is a slice object, a new qresult object with the same
        # instance attributes should be returned
        self.assertEqual(1102, self.qresult.seq_len)
        self.assertEqual('refseq_rna', self.qresult.target)
        new_qresult = self.qresult[1:]
        self.assertEqual([hit21, hit31], list(new_qresult.hits))
        self.assertEqual(1102, new_qresult.seq_len)
        self.assertEqual('refseq_rna', new_qresult.target)

    def test_getitm_slice_alt_ok(self):
        """Test QueryResult.__getitem__, with slice and alt IDs"""
        # slicing should be reflected in the alt IDs as well
        hit31._id_alt = ['alt3']
        hit11._id_alt = ['alt1']
        query = QueryResult([hit31, hit11])
        self.assertEqual(hit11, query['hit1'])
        self.assertEqual(hit11, query['alt1'])
        self.assertEqual(hit31, query['hit3'])
        self.assertEqual(hit31, query['alt3'])
        query = query[:1]
        self.assertEqual(hit31, query['hit3'])
        self.assertEqual(hit31, query['alt3'])
        self.assertRaises(KeyError, query.__getitem__, 'hit1')
        self.assertRaises(KeyError, query.__getitem__, 'alt1')
        hit31._id_alt = []
        hit11._id_alt = []

    def test_getitem_alt_ok(self):
        """Test QueryResult.__getitem__, single item with alternative ID"""
        hit11._id_alt = ['alt1']
        query = QueryResult([hit11])
        self.assertEqual(hit11, query['hit1'])
        self.assertEqual(hit11, query['alt1'])
        self.assertTrue(hit11.id != 'alt1')
        hit11._id_alt = []

    def test_delitem_string_ok(self):
        """Test QueryResult.__getitem__, with string"""
        # delitem should work with string index
        del self.qresult['hit1']
        self.assertEqual(2, len(self.qresult))
        self.assertTrue([hit21, hit31], list(self.qresult.hits))

    def test_delitem_int_ok(self):
        """Test QueryResult.__delitem__"""
        # delitem should work with int index
        del self.qresult[-1]
        self.assertEqual(2, len(self.qresult))
        self.assertEqual([hit11, hit21], list(self.qresult.hits))
        del self.qresult[0]
        self.assertEqual(1, len(self.qresult))
        self.assertTrue([hit21], list(self.qresult.hits))

    def test_delitem_slice_ok(self):
        """Test QueryResult.__delitem__, with slice"""
        # delitem should work with slice objects
        del self.qresult[:-1]
        self.assertEqual(1, len(self.qresult))
        self.assertTrue([hit31], self.qresult.hits)

    def test_delitem_alt_ok(self):
        """Test QueryResult.__delitem__, with alt ID"""
        # delitem should work with alt IDs
        hit31._id_alt = ['alt3']
        qresult = QueryResult([hit31, hit41])
        self.assertEqual(2, len(qresult))
        del qresult['alt3']
        self.assertEqual(1, len(qresult))
        self.assertEqual(hit41, qresult['hit4'])
        self.assertRaises(KeyError, qresult.__getitem__, 'alt3')
        hit31._id_alt = []

    def test_description_set(self):
        """Test QueryResult.description setter"""
        # setting the description should change the query seqrecord description
        # of the contained hsps, if they have an alignment
        # test for default value
        qresult = deepcopy(self.qresult)
        new_desc = 'unicorn hox homolog'
        # test initial condition
        for hit in qresult:
            self.assertNotEqual(new_desc, hit.query_description)
            for hsp in hit:
                self.assertNotEqual(new_desc, hsp.query_description)
                for fragment in hsp:
                    self.assertNotEqual(new_desc, fragment.query_description)
                    self.assertNotEqual(new_desc, fragment.query.description)
        qresult.description = new_desc
        # test after setting
        for hit in qresult:
            self.assertEqual(new_desc, hit.query_description)
            for hsp in hit:
                self.assertEqual(new_desc, hsp.query_description)
                for fragment in hsp:
                    self.assertEqual(new_desc, fragment.query_description)
                    self.assertEqual(new_desc, fragment.query.description)

    def test_description_set_no_seqrecord(self):
        """Test QueryResult.description setter, without HSP SeqRecords"""
        frag1 = HSPFragment('hit1', 'query')
        frag2 = HSPFragment('hit1', 'query')
        frag3 = HSPFragment('hit2', 'query')
        hit1 = Hit([HSP([x]) for x in [frag1, frag2]])
        hit2 = Hit([HSP([frag3])])
        qresult = QueryResult([hit1, hit2])
        # test initial condition
        for hit in qresult:
            for hsp in hit.hsps:
                self.assertTrue(getattr(hsp, 'query') is None)
        qresult.description = 'unicorn hox homolog'
        # test after setting
        for hit in qresult:
            for hsp in hit.hsps:
                self.assertTrue(getattr(hsp, 'query') is None)

    def test_id_set(self):
        """Test QueryResult.id setter"""
        # setting an ID should change the query IDs of all contained Hit and HSPs
        qresult = deepcopy(self.qresult)
        self.assertEqual('query1', qresult.id)
        for hit in qresult:
            self.assertEqual('query1', hit.query_id)
            for hsp in hit:
                self.assertEqual('query1', hsp.query_id)
                for fragment in hsp:
                    self.assertEqual('query1', fragment.query_id)
                    self.assertEqual('query1', fragment.query.id)
        qresult.id = 'new_id'
        self.assertEqual('new_id', qresult.id)
        for hit in qresult:
            self.assertEqual('new_id', hit.query_id)
            for hsp in hit:
                self.assertEqual('new_id', hsp.query_id)
                for fragment in hsp:
                    self.assertEqual('new_id', fragment.query_id)
                    self.assertEqual('new_id', fragment.query.id)

    def test_absorb_hit_does_not_exist(self):
        """Test QueryResult.absorb, hit does not exist"""
        # absorb should work like append when the hit does not exist
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))
        self.qresult.absorb(hit41)
        self.assertEqual([hit11, hit21, hit31, hit41], list(self.qresult.hits))
        self.assertEqual(['hit1', 'hit2', 'hit3', 'hit4'],
                         list(self.qresult.hit_keys))

    def test_absorb_hit_exists(self):
        """Test QueryResult.absorb, hit with the same ID exists"""
        # absorb should combine the hit's hsps if an existing one is present
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))
        self.assertEqual(2, len(self.qresult['hit3']))
        hit = Hit([HSP([HSPFragment('hit3', 'query1')])])
        self.qresult.absorb(hit)
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))
        self.assertEqual(['hit1', 'hit2', 'hit3'], list(self.qresult.hit_keys))
        self.assertEqual(3, len(self.qresult['hit3']))
        # remove the mock hsp
        del self.qresult['hit3'][-1]

    def test_append_ok(self):
        """Test QueryResult.append"""
        # append should work with Hit objects
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))
        self.qresult.append(hit41)
        self.assertEqual([hit11, hit21, hit31, hit41], list(self.qresult.hits))
        self.assertEqual(['hit1', 'hit2', 'hit3', 'hit4'],
                         list(self.qresult.hit_keys))

    def test_append_custom_hit_key_function_ok(self):
        """Test QueryResult.append, with custom hit key function"""
        self.qresult._hit_key_function = lambda hit: hit.id + '_custom'
        # append should assign hit keys according to _hit_key_function
        self.assertEqual(['hit1', 'hit2', 'hit3'], list(self.qresult.hit_keys))
        self.qresult.append(hit41)
        self.assertEqual(['hit1', 'hit2', 'hit3', 'hit4_custom'],
                         list(self.qresult.hit_keys))

    def test_append_id_exists(self):
        """Test QueryResult.append, when ID exists"""
        # append should raise an error if hit_key already exists
        self.assertRaises(ValueError, self.qresult.append, hit11)

    def test_append_alt_id_exists(self):
        """Test QueryResult.append, when alt ID exists"""
        # append should raise an error if hit_key already exists as alt ID
        hit11._id_alt = ['alt']
        hit21._id_alt = ['alt']
        qresult = QueryResult([hit11])
        self.assertRaises(ValueError, qresult.append, hit21)
        hit11._id_alt = []
        hit21._id_alt = []

    def test_append_alt_id_exists_alt(self):
        """Test QueryResult.append, when alt ID exists as primary"""
        # append should raise an error if alt ID already exists as primary ID
        hit21._id_alt = ['hit1']
        qresult = QueryResult([hit11])
        self.assertRaises(ValueError, qresult.append, hit21)
        hit21._id_alt = []

    def test_hit_filter(self):
        """Test QueryResult.hit_filter"""
        # hit_filter should return a new QueryResult object (shallow copy),
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))
        # filter func: min hit length == 2
        # this would filter out hit21, since it only has 1 HSP
        filter_func = lambda hit: len(hit) >= 2
        filtered = self.qresult.hit_filter(filter_func)
        self.assertEqual([hit11, hit31], list(filtered.hits))
        # make sure all remaining hits return True for the filter function
        self.assertTrue(all(filter_func(hit) for hit in filtered))
        self.assertEqual(1102, filtered.seq_len)
        self.assertEqual('refseq_rna', filtered.target)

    def test_hit_filter_no_func(self):
        """Test QueryResult.hit_filter, without arguments"""
        # when given no arguments, hit_filter should create a new object with
        # the same contents
        filtered = self.qresult.hit_filter()
        self.assertTrue(compare_search_obj(filtered, self.qresult))
        self.assertNotEqual(id(filtered), id(self.qresult))
        self.assertEqual(1102, filtered.seq_len)
        self.assertEqual('refseq_rna', filtered.target)

    def test_hit_filter_no_filtered(self):
        """Test QueryResult.hit_filter, all hits filtered out"""
        # when the filter filters out all hits, hit_filter should return an
        # empty QueryResult object
        filter_func = lambda hit: len(hit) > 50
        filtered = self.qresult.hit_filter(filter_func)
        self.assertEqual(0, len(filtered))
        self.assertTrue(isinstance(filtered, QueryResult))
        self.assertEqual(1102, filtered.seq_len)
        self.assertEqual('refseq_rna', filtered.target)

    def test_hit_map(self):
        """Test QueryResult.hit_map"""
        # hit_map should apply the given function to all contained Hits
        # deepcopy the qresult since we'll change the objects within
        qresult = deepcopy(self.qresult)
        # map func: capitalize hit IDs

        def map_func(hit):
            hit.id = hit.id.upper()
            return hit

        # test before mapping
        self.assertEqual('hit1', qresult[0].id)
        self.assertEqual('hit2', qresult[1].id)
        self.assertEqual('hit3', qresult[2].id)
        mapped = qresult.hit_map(map_func)
        self.assertEqual('HIT1', mapped[0].id)
        self.assertEqual('HIT2', mapped[1].id)
        self.assertEqual('HIT3', mapped[2].id)
        # and make sure the attributes are transferred
        self.assertEqual(1102, mapped.seq_len)
        self.assertEqual('refseq_rna', mapped.target)

    def test_hit_map_no_func(self):
        """Test QueryResult.hit_map, without arguments"""
        # when given no arguments, hit_map should create a new object with
        # the same contents
        mapped = self.qresult.hit_map()
        self.assertTrue(compare_search_obj(mapped, self.qresult))
        self.assertNotEqual(id(mapped), id(self.qresult))
        self.assertEqual(1102, mapped.seq_len)
        self.assertEqual('refseq_rna', mapped.target)

    def test_hsp_filter(self):
        """Test QueryResult.hsp_filter"""
        # hsp_filter should return a new QueryResult object (shallow copy)
        # and any empty hits should be discarded
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))
        # filter func: no '-' in hsp query sequence
        # this would filter out hsp113 and hsp211, effectively removing hit21
        filter_func = lambda hsp: '-' not in str(hsp.fragments[0].query)
        filtered = self.qresult.hsp_filter(filter_func)
        self.assertIn('hit1', filtered)
        self.assertNotIn('hit2', filtered)
        self.assertIn('hit3', filtered)
        # test hsps in hit11
        self.assertTrue(all(hsp in filtered['hit1'] for hsp in
                            [hsp111, hsp112, hsp114]))
        # test hsps in hit31
        self.assertTrue(all(hsp in filtered['hit3'] for hsp in
                            [hsp311, hsp312]))

    def test_hsp_filter_no_func(self):
        """Test QueryResult.hsp_filter, no arguments"""
        # when given no arguments, hsp_filter should create a new object with
        # the same contents
        filtered = self.qresult.hsp_filter()
        self.assertTrue(compare_search_obj(filtered, self.qresult))
        self.assertNotEqual(id(filtered), id(self.qresult))
        self.assertEqual(1102, filtered.seq_len)
        self.assertEqual('refseq_rna', filtered.target)

    def test_hsp_filter_no_filtered(self):
        """Test QueryResult.hsp_filter, all hits filtered out"""
        # when the filter filters out all hits, hsp_filter should return an
        # empty QueryResult object
        filter_func = lambda hsp: len(hsp) > 50
        filtered = self.qresult.hsp_filter(filter_func)
        self.assertEqual(0, len(filtered))
        self.assertTrue(isinstance(filtered, QueryResult))
        self.assertEqual(1102, filtered.seq_len)
        self.assertEqual('refseq_rna', filtered.target)

    def test_hsp_map(self):
        """Test QueryResult.hsp_map"""
        # hsp_map should apply the given function to all contained HSPs
        # deepcopy the qresult since we'll change the objects within
        qresult = deepcopy(self.qresult)
        # apply mock attributes to hsp, for testing mapped hsp attributes
        for hit in qresult:
            for hsp in hit:
                setattr(hsp, 'mock', 13)

        # map func: remove first letter of all HSP.aln
        def map_func(hsp):
            mapped_frags = [x[1:] for x in hsp]
            return HSP(mapped_frags)

        mapped = qresult.hsp_map(map_func)
        # make sure old hsp attributes is not transferred to mapped hsps
        for hit in mapped:
            for hsp in hit.hsps:
                self.assertFalse(hasattr(hsp, 'mock'))
        # check hsps in hit1
        self.assertEqual('TGCGCAT', str(mapped['hit1'][0][0].hit.seq))
        self.assertEqual('TGCGCAT', str(mapped['hit1'][0][0].query.seq))
        self.assertEqual('TG', str(mapped['hit1'][1][0].hit.seq))
        self.assertEqual('AT', str(mapped['hit1'][1][0].query.seq))
        self.assertEqual('TTCG', str(mapped['hit1'][2][0].hit.seq))
        self.assertEqual('T-CG', str(mapped['hit1'][2][0].query.seq))
        self.assertEqual('TTCG', str(mapped['hit1'][2][1].hit.seq))
        self.assertEqual('T-CG', str(mapped['hit1'][2][1].query.seq))
        self.assertEqual('T', str(mapped['hit1'][3][0].hit.seq))
        self.assertEqual('T', str(mapped['hit1'][3][0].query.seq))
        self.assertEqual('TCG', str(mapped['hit1'][3][1].hit.seq))
        self.assertEqual('TGG', str(mapped['hit1'][3][1].query.seq))
        # check hsps in hit2
        self.assertEqual('GGCCC', str(mapped['hit2'][0][0].hit.seq))
        self.assertEqual('GGCC-', str(mapped['hit2'][0][0].query.seq))
        # check hsps in hit3
        self.assertEqual('ATG', str(mapped['hit3'][0][0].hit.seq))
        self.assertEqual('TTG', str(mapped['hit3'][0][0].query.seq))
        self.assertEqual('TATAT', str(mapped['hit3'][1][0].hit.seq))
        self.assertEqual('TATAT', str(mapped['hit3'][1][0].query.seq))
        # and make sure the attributes are transferred
        self.assertEqual(1102, mapped.seq_len)
        self.assertEqual('refseq_rna', mapped.target)

    def test_hsp_map_no_func(self):
        """Test QueryResult.hsp_map, without arguments"""
        # when given no arguments, hit_map should create a new object with
        # the same contents
        mapped = self.qresult.hsp_map()
        self.assertTrue(compare_search_obj(mapped, self.qresult))
        self.assertNotEqual(id(mapped), id(self.qresult))
        self.assertEqual(1102, mapped.seq_len)
        self.assertEqual('refseq_rna', mapped.target)

    def test_pop_ok(self):
        """Test QueryResult.pop"""
        self.assertEqual(3, len(self.qresult))
        hit = self.qresult.pop()
        self.assertEqual(hit, hit31)
        self.assertEqual([hit11, hit21], list(self.qresult.hits))

    def test_pop_int_index_ok(self):
        """Test QueryResult.pop, with integer index"""
        # pop should work if given an int index
        self.assertEqual(3, len(self.qresult))
        hit = self.qresult.pop(1)
        self.assertEqual(hit, hit21)
        self.assertEqual([hit11, hit31], list(self.qresult.hits))

    def test_pop_string_index_ok(self):
        """Test QueryResult.pop, with string index"""
        # pop should work if given a string index
        self.assertEqual(3, len(self.qresult))
        hit = self.qresult.pop('hit2')
        self.assertEqual(hit, hit21)
        self.assertEqual([hit11, hit31], list(self.qresult.hits))

    def test_pop_string_alt_ok(self):
        """Test QueryResult.pop, with alternative ID"""
        # pop should work with alternative index
        hit11._id_alt = ['alt1']
        hit21._id_alt = ['alt2']
        qresult = QueryResult([hit11, hit21])
        hit = qresult.pop('alt1')
        self.assertEqual(hit, hit11)
        self.assertEqual([hit21], list(qresult))
        self.assertNotIn('hit1', qresult)
        hit11._id_alt = []
        hit21._id_alt = []

    def test_index(self):
        """Test QueryResult.index"""
        # index should accept hit objects or hit key strings
        self.assertEqual(2, self.qresult.index('hit3'))
        self.assertEqual(2, self.qresult.index(hit31))

    def test_index_alt(self):
        """Test QueryResult.index, with alt ID"""
        # index should work with alt IDs
        hit11._id_alt = ['alt1']
        qresult = QueryResult([hit21, hit11])
        self.assertEqual(1, qresult.index('alt1'))
        hit11._id_alt = []

    def test_index_not_present(self):
        """Test QueryResult.index, when index is not present"""
        self.assertRaises(ValueError, self.qresult.index, 'hit4')
        self.assertRaises(ValueError, self.qresult.index, hit41)

    def test_sort_ok(self):
        """Test QueryResult.sort"""
        # sort without any arguments should keep the Hits in the same order
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))
        self.qresult.sort()
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))

    def test_sort_not_in_place_ok(self):
        """Test QueryResult.sort, not in place"""
        # sort without any arguments should keep the Hits in the same order
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))
        sorted_qresult = self.qresult.sort(in_place=False)
        self.assertEqual([hit11, hit21, hit31], list(sorted_qresult.hits))
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))

    def test_sort_reverse_ok(self):
        """Test QueryResult.sort, reverse"""
        # sorting with reverse=True should return a QueryResult with Hits reversed
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))
        self.qresult.sort(reverse=True)
        self.assertEqual([hit31, hit21, hit11], list(self.qresult.hits))

    def test_sort_reverse_not_in_place_ok(self):
        """Test QueryResult.sort, reverse, not in place"""
        # sorting with reverse=True should return a QueryResult with Hits reversed
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))
        sorted_qresult = self.qresult.sort(reverse=True, in_place=False)
        self.assertEqual([hit31, hit21, hit11], list(sorted_qresult.hits))
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))

    def test_sort_key_ok(self):
        """Test QueryResult.sort, with custom key"""
        # if custom key is given, sort using it
        key = lambda hit: len(hit)
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))
        self.qresult.sort(key=key)
        self.assertEqual([hit21, hit31, hit11], list(self.qresult.hits))

    def test_sort_key_not_in_place_ok(self):
        """Test QueryResult.sort, with custom key, not in place"""
        # if custom key is given, sort using it
        key = lambda hit: len(hit)
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))
        sorted_qresult = self.qresult.sort(key=key, in_place=False)
        self.assertEqual([hit21, hit31, hit11], list(sorted_qresult.hits))
        self.assertEqual([hit11, hit21, hit31], list(self.qresult.hits))


class HitCases(unittest.TestCase):

    def setUp(self):
        self.hit = Hit([hsp111, hsp112, hsp113])
        self.hit.evalue = 5e-10
        self.hit.name = 'test'

    def test_pickle(self):
        """Test pickling and unpickling of Hit"""
        buf = BytesIO()
        pickle.dump(self.hit, buf)
        unp = pickle.loads(buf.getvalue())
        self.assertTrue(compare_search_obj(self.hit, unp))

    def test_init_none(self):
        """Test Hit.__init__, no arguments"""
        hit = Hit()
        self.assertEqual(None, hit.id)
        self.assertEqual(None, hit.description)
        self.assertEqual(None, hit.query_id)
        self.assertEqual(None, hit.query_description)

    def test_init_id_only(self):
        """Test Hit.__init__, with ID only"""
        hit = Hit(id='hit1')
        self.assertEqual('hit1', hit.id)
        self.assertEqual(None, hit.description)
        self.assertEqual(None, hit.query_id)
        self.assertEqual(None, hit.query_description)

    def test_init_hsps_only(self):
        """Test Hit.__init__, with hsps only"""
        hit = Hit([hsp111, hsp112, hsp113])
        self.assertEqual('hit1', hit.id)
        self.assertEqual('<unknown description>', hit.description)
        self.assertEqual('query1', hit.query_id)  # set from the HSPs
        self.assertEqual('<unknown description>', hit.query_description)

    def test_repr(self):
        """Test Hit.__repr__"""
        # test for cases with 1 or other alignment numbers
        self.assertEqual("Hit(id='hit1', query_id='query1', 3 hsps)",
                         repr(self.hit))

    def test_hsps(self):
        """Test Hit.hsps"""
        # hsps should return the list of hsps contained
        self.assertEqual([hsp111, hsp112, hsp113], self.hit.hsps)

    def test_fragments(self):
        """Test Hit.fragments"""
        # fragments should return the list of fragments in each hsps
        # as a flat list
        self.assertEqual([frag111, frag112, frag113, frag113b],
                         self.hit.fragments)

    def test_iter(self):
        """Test Hit.__iter__"""
        # iteration should return hsps contained
        for counter, hsp in enumerate(self.hit):
            self.assertIn(hsp, [hsp111, hsp112, hsp113])
        self.assertEqual(2, counter)

    def test_len(self):
        """Test Hit.__len__"""
        # len() on Hit objects should return how many hsps it has
        self.assertEqual(3, len(self.hit))

    def test_nonzero(self):
        """Test Hit.__nonzero__"""
        # bool() on Hit objects should return True only if hsps is filled
        # which is always true
        self.assertTrue(self.hit)

    def test_setitem_single(self):
        """Test Hit.__setitem__, single item"""
        # test regular setitem overwrite
        self.hit[1] = hsp114
        self.assertEqual(self.hit.hsps, [hsp111, hsp114, hsp113])

    def test_item_multiple(self):
        """Test Hit.__setitem__, multiple items"""
        # test iterable setitem
        self.hit[:] = [hsp113, hsp112, hsp111]
        self.assertEqual(self.hit.hsps, [hsp113, hsp112, hsp111])

    def test_getitem_single(self):
        """Test Hit.__getitem__, single item"""
        # getitem using integer index should return a hsp object
        hsp1 = self.hit[0]
        self.assertEqual(hsp111, hsp1)
        hsp3 = self.hit[-1]
        self.assertEqual(hsp113, hsp3)

    def test_getitem_multiple(self):
        """Test Hit.__getitem__, multiple items"""
        # getitem using slices should return another hit object
        # with the hsps sliced accordingly, but other attributes preserved
        new_hit = self.hit[:2]
        self.assertEqual(2, len(new_hit))
        self.assertEqual([hsp111, hsp112], new_hit.hsps)
        self.assertEqual(self.hit.id, new_hit.id)
        self.assertEqual(self.hit.query_id, new_hit.query_id)
        self.assertEqual(5e-10, new_hit.evalue)
        self.assertEqual('test', new_hit.name)

    def test_delitem(self):
        """Test Hit.__delitem__"""
        # test delitem
        del self.hit[0]
        self.assertEqual(2, len(self.hit))
        self.assertEqual([hsp112, hsp113], self.hit.hsps)

    def test_validate_hsp_ok(self):
        """Test Hit._validate_hsp"""
        # validation should pass if item is an hsp object with matching
        # query and hit ids
        # if validation passes, None is returned
        self.assertEqual(None, self.hit._validate_hsp(hsp114))

    def test_validate_hsp_wrong_type(self):
        """Test Hit._validate_hsp, wrong type"""
        # validation should fail if item is not an hsp object
        self.assertRaises(TypeError, self.hit._validate_hsp, 1)
        self.assertRaises(TypeError, self.hit._validate_hsp, Seq(''))

    def test_validate_hsp_wrong_query_id(self):
        """Test Hit._validate_hsp, wrong query ID"""
        # validation should fail if query id does not match
        self.assertRaises(ValueError, self.hit._validate_hsp, hsp211)

    def test_validate_hsp_wrong_hit_id(self):
        """Test Hit._validate_hsp, wrong hit ID"""
        # validation should vail if hit id does not match
        self.assertRaises(ValueError, self.hit._validate_hsp, hsp121)

    def test_desc_set(self):
        """Test Hit.description setter"""
        # setting the description should change the hit seqrecord description
        # of the contained hsps, if they have an alignment
        # test for default value
        hit = deepcopy(self.hit)
        new_desc = 'unicorn hox homolog'
        # test initial condition
        for hsp in hit:
            self.assertNotEqual(new_desc, hsp.hit_description)
            for fragment in hsp:
                self.assertNotEqual(new_desc, fragment.hit_description)
                self.assertNotEqual(new_desc, fragment.hit.description)
        hit.description = new_desc
        # test after setting
        for hsp in hit:
            self.assertEqual(new_desc, hsp.hit_description)
            for fragment in hsp:
                self.assertEqual(new_desc, fragment.hit_description)
                self.assertEqual(new_desc, fragment.hit.description)

    def test_desc_set_no_seqrecord(self):
        """Test Hit.description setter, without HSP SeqRecords"""
        frag1 = HSPFragment('hit1', 'query')
        frag2 = HSPFragment('hit1', 'query')
        hit = Hit([HSP([x]) for x in [frag1, frag2]])
        new_desc = 'unicorn hox homolog'
        # test initial condition
        self.assertEqual(hit.description, '<unknown description>')
        for hsp in hit:
            self.assertEqual(hsp.hit_description, '<unknown description>')
            for fragment in hsp:
                self.assertEqual(hsp.hit_description, '<unknown description>')
        hit.description = new_desc
        # test after setting
        self.assertEqual(hit.description, new_desc)
        for hsp in hit:
            self.assertTrue(hsp.hit_description, new_desc)
            for fragment in hsp:
                self.assertEqual(hsp.hit_description, new_desc)

    def test_id_set(self):
        """Test Hit.id setter"""
        # setting an ID should change the query IDs of all contained HSPs
        hit = deepcopy(self.hit)
        self.assertEqual('hit1', hit.id)
        for hsp in hit.hsps:
            self.assertEqual('hit1', hsp.hit_id)
            for fragment in hsp:
                self.assertEqual(fragment.hit_id, 'hit1')
                self.assertEqual(fragment.hit.id, 'hit1')
        hit.id = 'new_id'
        self.assertEqual('new_id', hit.id)
        for hsp in hit.hsps:
            self.assertEqual('new_id', hsp.hit_id)
            for fragment in hsp:
                self.assertEqual(fragment.hit_id, 'new_id')
                self.assertEqual(fragment.hit.id, 'new_id')

    def test_append(self):
        """Test Hit.append"""
        # append should add hits to the last position
        self.hit.append(hsp114)
        self.assertEqual(4, len(self.hit))
        self.assertEqual(hsp114, self.hit[-1])

    def test_filter(self):
        """Test Hit.filter"""
        # filter should return a new QueryResult object (shallow copy),
        self.assertEqual([hsp111, hsp112, hsp113], self.hit.hsps)
        # filter func: min hsp length == 4
        filter_func = lambda hsp: len(hsp[0]) >= 4
        filtered = self.hit.filter(filter_func)
        self.assertEqual([hsp111, hsp113], filtered.hsps)
        # make sure all remaining hits return True for the filter function
        self.assertTrue(all(filter_func(hit) for hit in filtered))
        self.assertEqual(5e-10, filtered.evalue)
        self.assertEqual('test', filtered.name)

    def test_filter_no_func(self):
        """Test Hit.filter, without arguments"""
        # when given no arguments, filter should create a new object with
        # the same contents
        filtered = self.hit.filter()
        self.assertTrue(compare_search_obj(filtered, self.hit))
        self.assertNotEqual(id(filtered), id(self.hit))
        self.assertEqual(5e-10, filtered.evalue)
        self.assertEqual('test', filtered.name)

    def test_filter_no_filtered(self):
        """Test Hit.hit_filter, all hits filtered out"""
        # when the filter filters out all hits, it should return None
        filter_func = lambda hsp: len(hsp[0]) > 50
        filtered = self.hit.filter(filter_func)
        self.assertTrue(filtered is None)

    def test_index(self):
        """Test Hit.index"""
        # index should accept hsp objects
        self.assertEqual(1, self.hit.index(hsp112))

    def test_index_not_present(self):
        """Test Hit.index, when index is not present"""
        self.assertRaises(ValueError, self.hit.index, hsp114)

    def test_map(self):
        """Test Hit.hsp_map"""
        # map should apply the given function to all contained HSPs
        # deepcopy hit since we'll change the objects within
        hit = deepcopy(self.hit)
        # apply mock attributes to hsp, for testing mapped hsp attributes
        for hsp in hit:
            setattr(hsp, 'mock', 13)

        # map func: remove first letter of all HSP.alignment
        def map_func(hsp):
            mapped_frags = [x[1:] for x in hsp]
            return HSP(mapped_frags)

        mapped = hit.map(map_func)
        # make sure old hsp attributes is not transferred to mapped hsps
        for hsp in mapped:
            self.assertFalse(hasattr(hsp, 'mock'))
        # check hsps in hit1
        self.assertEqual('TGCGCAT', str(mapped[0][0].hit.seq))
        self.assertEqual('TGCGCAT', str(mapped[0][0].query.seq))
        self.assertEqual('TG', str(mapped[1][0].hit.seq))
        self.assertEqual('AT', str(mapped[1][0].query.seq))
        self.assertEqual('TTCG', str(mapped[2][0].hit.seq))
        self.assertEqual('T-CG', str(mapped[2][0].query.seq))
        self.assertEqual('TTCG', str(mapped[2][1].hit.seq))
        self.assertEqual('T-CG', str(mapped[2][1].query.seq))
        # and make sure the attributes are transferred
        self.assertEqual(5e-10, mapped.evalue)
        self.assertEqual('test', mapped.name)

    def test_hsp_map_no_func(self):
        """Test Hit.map, without arguments"""
        # when given no arguments, map should create a new object with
        # the same contents
        mapped = self.hit.map()
        self.assertTrue(compare_search_obj(mapped, self.hit))
        self.assertNotEqual(id(mapped), id(self.hit))
        self.assertEqual(5e-10, mapped.evalue)
        self.assertEqual('test', mapped.name)

    def test_pop(self):
        """Test Hit.pop"""
        # pop should return the last item by default
        self.assertEqual(hsp113, self.hit.pop())
        self.assertEqual(hsp111, self.hit.pop(0))

    def test_sort(self):
        """Test Hit.sort"""
        self.assertEqual([hsp111, hsp112, hsp113], self.hit.hsps)
        # sort by hsp length
        key = lambda batch_hsp: len(batch_hsp[0])
        self.hit.sort(key=key)
        self.assertEqual([hsp112, hsp113, hsp111], self.hit.hsps)

    def test_sort_not_in_place(self):
        """Test Hit.sort, not in place"""
        self.assertEqual([hsp111, hsp112, hsp113], self.hit.hsps)
        # sort by hsp length
        key = lambda hsp: len(hsp[0])
        sorted_hit = self.hit.sort(key=key, in_place=False)
        self.assertEqual([hsp112, hsp113, hsp111], sorted_hit.hsps)
        self.assertEqual([hsp111, hsp112, hsp113], self.hit.hsps)
        self.assertEqual(5e-10, sorted_hit.evalue)
        self.assertEqual('test', sorted_hit.name)


class HSPSingleFragmentCases(unittest.TestCase):

    def setUp(self):
        self.frag = HSPFragment('hit_id', 'query_id', 'ATCAGT', 'AT-ACT')
        self.frag.query_start = 0
        self.frag.query_end = 6
        self.frag.hit_start = 15
        self.frag.hit_end = 20
        self.hsp = HSP([self.frag])

    def test_init_no_fragment(self):
        """Test HSP.__init__ without fragments"""
        self.assertRaises(ValueError, HSP, [])

    def test_len(self):
        """Test HSP.__len__"""
        self.assertEqual(1, len(self.hsp))

    def test_fragment(self):
        """Test HSP.fragment property"""
        self.assertTrue(self.frag is self.hsp.fragment)

    def test_is_fragmented(self):
        """Test HSP.is_fragmented property"""
        self.assertFalse(self.hsp.is_fragmented)

    def test_seq(self):
        """Test HSP sequence properties"""
        self.assertEqual('ATCAGT', str(self.hsp.hit.seq))
        self.assertEqual('AT-ACT', str(self.hsp.query.seq))

    def test_alignment(self):
        """Test HSP.alignment property"""
        aln = self.hsp.aln
        self.assertTrue(isinstance(aln, MultipleSeqAlignment))
        self.assertEqual(2, len(aln))
        self.assertTrue('ATCAGT', str(aln[0].seq))
        self.assertTrue('AT-ACT', str(aln[1].seq))

    def test_aln_span(self):
        """Test HSP.aln_span property"""
        self.assertEqual(6, self.hsp.aln_span)

    def test_span(self):
        """Test HSP span properties"""
        self.assertEqual(5, self.hsp.hit_span)
        self.assertEqual(6, self.hsp.query_span)

    def test_range(self):
        """Test HSP range properties"""
        self.assertEqual((15, 20), self.hsp.hit_range)
        self.assertEqual((0, 6), self.hsp.query_range)

    def test_setters_readonly(self):
        """Test HSP read-only properties"""
        read_onlies = ('range', 'span', 'strand', 'frame', 'start', 'end')
        for seq_type in ('query', 'hit'):
            self.assertRaises(AttributeError, setattr,
                              self.hsp, seq_type, 'A')
            for attr in read_onlies:
                self.assertRaises(AttributeError, setattr,
                                  self.hsp, '%s_%s' % (seq_type, attr), 5)
        self.assertRaises(AttributeError, setattr,
                          self.hsp, 'aln', None)


class HSPMultipleFragmentCases(unittest.TestCase):

    def setUp(self):
        self.frag1 = HSPFragment('hit_id', 'query_id', 'ATCAGT', 'AT-ACT')
        self.frag1.query_start = 0
        self.frag1.query_end = 6
        self.frag1.hit_start = 15
        self.frag1.hit_end = 20
        self.frag2 = HSPFragment('hit_id', 'query_id', 'GGG', 'CCC')
        self.frag2.query_start = 10
        self.frag2.query_end = 13
        self.frag2.hit_start = 158
        self.frag2.hit_end = 161
        self.hsp = HSP([self.frag1, self.frag2])

    def test_pickle(self):
        """Test pickling and unpickling of HSP"""
        buf = BytesIO()
        pickle.dump(self.hsp, buf)
        unp = pickle.loads(buf.getvalue())
        self.assertTrue(compare_search_obj(self.hsp, unp))

    def test_len(self):
        """Test HSP.__len__"""
        self.assertEqual(2, len(self.hsp))

    def test_getitem(self):
        """Test HSP.__getitem__"""
        self.assertTrue(self.frag1 is self.hsp[0])
        self.assertTrue(self.frag2 is self.hsp[1])

    def test_setitem_single(self):
        """Test HSP.__setitem___, single item"""
        frag3 = HSPFragment('hit_id', 'query_id', 'AAA', 'AAT')
        self.hsp[1] = frag3
        self.assertEqual(2, len(self.hsp))
        self.assertTrue(self.frag1 is self.hsp[0])
        self.assertTrue(frag3 is self.hsp[1])

    def test_setitem_multiple(self):
        """Test HSP.__setitem__, multiple items"""
        frag3 = HSPFragment('hit_id', 'query_id', 'AAA', 'AAT')
        frag4 = HSPFragment('hit_id', 'query_id', 'GGG', 'GAG')
        self.hsp[:2] = [frag3, frag4]
        self.assertEqual(2, len(self.hsp))
        self.assertTrue(frag3 is self.hsp[0])
        self.assertTrue(frag4 is self.hsp[1])

    def test_delitem(self):
        """Test HSP.__delitem__"""
        del self.hsp[0]
        self.assertEqual(1, len(self.hsp))
        self.assertTrue(self.frag2 is self.hsp[0])

    def test_contains(self):
        """Test HSP.__contains__"""
        frag3 = HSPFragment('hit_id', 'query_id', 'AAA', 'AAT')
        self.assertIn(self.frag1, self.hsp)
        self.assertNotIn(frag3, self.hsp)

    def test_fragments(self):
        """Test HSP.fragments property"""
        self.assertEqual([self.frag1, self.frag2], self.hsp.fragments)

    def test_is_fragmented(self):
        """Test HSP.is_fragmented property"""
        self.assertTrue(self.hsp.is_fragmented)

    def test_seqs(self):
        """Test HSP sequence properties"""
        self.assertEqual(['ATCAGT', 'GGG'], [str(x.seq) for x in
                         self.hsp.hit_all])
        self.assertEqual(['AT-ACT', 'CCC'], [str(x.seq) for x in
                         self.hsp.query_all])

    def test_id_desc_set(self):
        """Test HSP query and hit id and description setters"""
        for seq_type in ('query', 'hit'):
            for attr in ('id', 'description'):
                attr_name = '%s_%s' % (seq_type, attr)
                value = getattr(self.hsp, attr_name)
                if attr == 'id':
                    # because we happen to have the same value for
                    # IDs and the actual attribute name
                    self.assertEqual(value, attr_name)
                    for fragment in self.hsp:
                        self.assertEqual(getattr(fragment, attr_name), attr_name)
                else:
                    self.assertEqual(value, '<unknown description>')
                    for fragment in self.hsp:
                        self.assertEqual(getattr(fragment, attr_name), '<unknown description>')
                new_value = 'new_' + value
                setattr(self.hsp, attr_name, new_value)
                self.assertEqual(getattr(self.hsp, attr_name), new_value)
                self.assertNotEqual(getattr(self.hsp, attr_name), value)
                for fragment in self.hsp:
                    self.assertEqual(getattr(fragment, attr_name), new_value)
                    self.assertNotEqual(getattr(fragment, attr_name), value)

    def test_alphabet(self):
        """Test HSP.alphabet getter"""
        self.assertTrue(self.hsp.alphabet is single_letter_alphabet)

    def test_alphabet_set(self):
        """Test HSP.alphabet setter"""
        # test initial values
        self.assertTrue(self.hsp.alphabet is single_letter_alphabet)
        for frag in self.hsp.fragments:
            self.assertTrue(frag.alphabet is single_letter_alphabet)
        self.hsp.alphabet = generic_dna
        # test values after setting
        self.assertTrue(self.hsp.alphabet is generic_dna)
        for frag in self.hsp.fragments:
            self.assertTrue(frag.alphabet is generic_dna)

    def test_range(self):
        """Test HSP range properties"""
        # range on HSP with multiple fragment should give the
        # min start and max end coordinates
        self.assertEqual((15, 161), self.hsp.hit_range)
        self.assertEqual((0, 13), self.hsp.query_range)

    def test_ranges(self):
        """Test HSP ranges properties"""
        self.assertEqual([(15, 20), (158, 161)], self.hsp.hit_range_all)
        self.assertEqual([(0, 6), (10, 13)], self.hsp.query_range_all)

    def test_span(self):
        """Test HSP span properties"""
        # span is always end - start
        self.assertEqual(146, self.hsp.hit_span)
        self.assertEqual(13, self.hsp.query_span)

    def test_setters_readonly(self):
        """Test HSP read-only properties"""
        read_onlies = ('range_all', 'strand_all', 'frame_all')
        for seq_type in ('query', 'hit'):
            for attr in read_onlies:
                self.assertRaises(AttributeError, setattr,
                                  self.hsp, '%s_%s' % (seq_type, attr), 5)
        self.assertRaises(AttributeError, setattr,
                          self.hsp, 'aln_all', None)
        self.assertRaises(AttributeError, setattr,
                          self.hsp, 'hit_all', None)
        self.assertRaises(AttributeError, setattr,
                          self.hsp, 'query_all', None)


class HSPFragmentWithoutSeqCases(unittest.TestCase):

    def setUp(self):
        self.fragment = HSPFragment('hit_id', 'query_id')

    def test_init(self):
        """Test HSPFragment.__init__ attributes"""
        fragment = HSPFragment('hit_id', 'query_id')
        for seq_type in ('query', 'hit'):
            self.assertTrue(getattr(fragment, seq_type) is None)
            for attr in ('strand', 'frame', 'start', 'end'):
                attr_name = '%s_%s' % (seq_type, attr)
                self.assertTrue(getattr(fragment, attr_name) is None)
        self.assertTrue(fragment.aln is None)
        self.assertTrue(fragment.alphabet is single_letter_alphabet)
        self.assertEqual(fragment.aln_annotation, {})

    def test_seqmodel(self):
        """Test HSPFragment sequence attributes, no alignments"""
        # all query, hit, and alignment objects should be None
        self.assertTrue(self.fragment.query is None)
        self.assertTrue(self.fragment.hit is None)
        self.assertTrue(self.fragment.aln is None)

    def test_len(self):
        """Test HSPFragment.__len__, no alignments"""
        self.assertRaises(TypeError, len, self)
        # len is a shorthand for .aln_span, and it can be set manually
        self.fragment.aln_span = 5
        self.assertEqual(5, len(self.fragment))

    def test_repr(self):
        """Test HSPFragment.__repr__, no alignments"""
        # test for minimum repr
        self.assertEqual("HSPFragment(hit_id='hit_id', query_id='query_id')",
                         repr(self.fragment))
        self.fragment.aln_span = 5
        self.assertEqual("HSPFragment(hit_id='hit_id', query_id='query_id', "
                         "5 columns)", repr(self.fragment))

    def test_getitem(self):
        """Test HSPFragment.__getitem__, no alignments"""
        # getitem not supported without alignment
        self.assertRaises(TypeError, self.fragment.__getitem__, 0)
        self.assertRaises(TypeError, self.fragment.__getitem__, slice(0, 2))

    def test_getitem_only_query(self):
        """Test HSPFragment.__getitem__, only query"""
        # getitem should work if only query is present
        self.fragment.query = 'AATCG'
        self.assertEqual('ATCG', str(self.fragment[1:].query.seq))

    def test_getitem_only_hit(self):
        """Test HSPFragment.__getitem__, only hit"""
        # getitem should work if only query is present
        self.fragment.hit = 'CATGC'
        self.assertEqual('ATGC', str(self.fragment[1:].hit.seq))

    def test_iter(self):
        """Test HSP.__iter__, no alignments"""
        # iteration not supported
        self.assertRaises(TypeError, iter, self)


class HSPFragmentCases(unittest.TestCase):

    def setUp(self):
        self.fragment = HSPFragment('hit_id', 'query_id', 'ATGCTAGCTACA',
                                    'ATG--AGCTAGG')

    def test_pickle(self):
        """Test pickling and unpickling of HSPFragment"""
        buf = BytesIO()
        pickle.dump(self.fragment, buf)
        unp = pickle.loads(buf.getvalue())
        self.assertTrue(compare_search_obj(self.fragment, unp))

    def test_init_with_seqrecord(self):
        """Test HSPFragment.__init__, with SeqRecord"""
        # init should work with seqrecords
        hit_seq = SeqRecord(Seq('ATGCTAGCTACA'))
        query_seq = SeqRecord(Seq('ATG--AGCTAGG'))
        hsp = HSPFragment('hit_id', 'query_id', hit_seq, query_seq)
        self.assertTrue(isinstance(hsp.query, SeqRecord))
        self.assertTrue(isinstance(hsp.hit, SeqRecord))
        self.assertTrue(isinstance(hsp.aln, MultipleSeqAlignment))

    def test_init_wrong_seqtypes(self):
        """Test HSPFragment.__init__, wrong sequence argument types"""
        # init should only work with string or seqrecords
        wrong_query = Seq('ATGC')
        wrong_hit = Seq('ATGC')
        self.assertRaises(TypeError, HSPFragment, 'hit_id', 'query_id',
                          wrong_hit, wrong_query)

    def test_seqmodel(self):
        """Test HSPFragment sequence attribute types and default values"""
        # check hit
        self.assertTrue(isinstance(self.fragment.hit, SeqRecord))
        self.assertEqual('<unknown description>', self.fragment.hit.description)
        self.assertEqual('aligned hit sequence', self.fragment.hit.name)
        self.assertEqual(single_letter_alphabet, self.fragment.hit.seq.alphabet)
        # check query
        self.assertTrue(isinstance(self.fragment.query, SeqRecord))
        self.assertEqual('<unknown description>', self.fragment.query.description)
        self.assertEqual('aligned query sequence', self.fragment.query.name)
        self.assertEqual(single_letter_alphabet, self.fragment.query.seq.alphabet)
        # check alignment
        self.assertTrue(isinstance(self.fragment.aln, MultipleSeqAlignment))
        self.assertEqual(single_letter_alphabet, self.fragment.aln._alphabet)

    def test_alphabet_no_seq(self):
        """Test HSPFragment alphabet property, query and hit sequences not present"""
        self.assertTrue(self.fragment.alphabet is single_letter_alphabet)
        self.fragment.alphabet = generic_dna
        self.assertTrue(self.fragment.alphabet is generic_dna)

    def test_alphabet_with_seq(self):
        """Test HSPFragment alphabet property, query or hit sequences present"""
        self.assertTrue(self.fragment.alphabet is single_letter_alphabet)
        self.fragment._hit = SeqRecord(Seq('AAA'))
        self.fragment._query = SeqRecord(Seq('AAA'))
        self.fragment.alphabet = generic_dna
        self.assertTrue(self.fragment.alphabet is generic_dna)
        self.assertTrue(self.fragment.hit.seq.alphabet is generic_dna)
        self.assertTrue(self.fragment.query.seq.alphabet is generic_dna)

    def test_seq_unequal_hit_query_len(self):
        """Test HSPFragment sequence setter with unequal hit and query lengths"""
        for seq_type in ('hit', 'query'):
            opp_type = 'query' if seq_type == 'hit' else 'hit'
            # reset values first
            fragment = HSPFragment('hit_id', 'query_id')
            # and test it against the opposite
            setattr(fragment, seq_type, 'ATGCACAACAGGA')
            self.assertRaises(ValueError, setattr, fragment, opp_type, 'ATGCGA')

    def test_len(self):
        """Test HSPFragment.__len__"""
        # len should equal alignment column length
        self.assertEqual(12, len(self.fragment))

    def test_repr(self):
        """Test HSPFragment.__repr__"""
        # test for minimum repr
        self.assertEqual("HSPFragment(hit_id='hit_id', query_id='query_id', "
                         "12 columns)", repr(self.fragment))

    def test_getitem(self):
        """Test HSPFragment.__getitem__"""
        # getitem is supported when alignment is present
        sliced_fragment = self.fragment[:5]
        self.assertTrue(isinstance(sliced_fragment, HSPFragment))
        self.assertEqual(5, len(sliced_fragment))
        self.assertEqual('ATGCT', str(sliced_fragment.hit.seq))
        self.assertEqual('ATG--', str(sliced_fragment.query.seq))

    def test_getitem_attrs(self):
        """Test HSPFragment.__getitem__, with attributes"""
        # attributes from the original instance should not be present in the new
        # objects, except for query, hit, and alignment - related attributes
        setattr(self.fragment, 'attr_original', 1000)
        setattr(self.fragment, 'hit_description', 'yeah')
        setattr(self.fragment, 'hit_strand', 1)
        setattr(self.fragment, 'query_frame', None)
        # test values prior to slicing
        self.assertEqual(1000, getattr(self.fragment, 'attr_original'))
        self.assertEqual('yeah', getattr(self.fragment, 'hit_description'))
        self.assertEqual(1, getattr(self.fragment, 'hit_strand'))
        self.assertEqual(None, getattr(self.fragment, 'query_frame'))
        new_hsp = self.fragment[:5]
        # test values after slicing
        self.assertFalse(hasattr(new_hsp, 'attr_original'))
        self.assertEqual(1000, getattr(self.fragment, 'attr_original'))
        self.assertEqual('yeah', getattr(self.fragment, 'hit_description'))
        self.assertEqual(1, getattr(self.fragment, 'hit_strand'))
        self.assertEqual(None, getattr(self.fragment, 'query_frame'))

    def test_getitem_alignment_annot(self):
        """Test HSPFragment.__getitem__, with alignment annotation"""
        # the alignment is annotated, it should be sliced accordingly
        # and transferred to the new object
        setattr(self.fragment, 'aln_annotation', {'test': '182718738172'})
        new_hsp = self.fragment[:5]
        self.assertEqual('18271', new_hsp.aln_annotation['test'])

    def test_default_attrs(self):
        """Test HSPFragment attributes' default values"""
        fragment = HSPFragment()
        self.assertEqual('<unknown id>', fragment.hit_id)
        self.assertEqual('<unknown id>', fragment.query_id)
        self.assertEqual('<unknown description>', fragment.hit_description)
        self.assertEqual('<unknown description>', fragment.query_description)
        self.assertEqual(None, fragment.hit)
        self.assertEqual(None, fragment.query)
        self.assertEqual(None, fragment.aln)
        self.assertEqual([], fragment.hit_features)
        self.assertEqual([], fragment.query_features)
        self.assertEqual(None, fragment.hit_strand)
        self.assertEqual(None, fragment.query_strand)
        self.assertEqual(None, fragment.hit_frame)
        self.assertEqual(None, fragment.query_frame)

    def test_id_desc_set(self):
        """Test HSPFragment query and hit id and description setters"""
        for seq_type in ('query', 'hit'):
            for attr in ('id', 'description'):
                attr_name = '%s_%s' % (seq_type, attr)
                value = getattr(self.fragment, attr_name)
                if attr == 'id':
                    # because we happen to have the same value for
                    # IDs and the actual attribute name
                    self.assertEqual(value, attr_name)
                else:
                    self.assertEqual(value, '<unknown description>')
                new_value = 'new_' + value
                setattr(self.fragment, attr_name, new_value)
                self.assertEqual(getattr(self.fragment, attr_name), new_value)
                self.assertNotEqual(getattr(self.fragment, attr_name), value)

    def test_frame_set_ok(self):
        """Test HSPFragment query and hit frame setters"""
        attr = 'frame'
        for seq_type in ('query', 'hit'):
            attr_name = '%s_%s' % (seq_type, attr)
            for value in (-3, -2, -1, 0, 1, 2, 3, None):
                setattr(self.fragment, attr_name, value)
                self.assertEqual(value, getattr(self.fragment, attr_name))

    def test_frame_set_error(self):
        """Test HSPFragment query and hit frame setters, invalid values"""
        attr = 'frame'
        for seq_type in ('query', 'hit'):
            func_name = '_%s_%s_set' % (seq_type, attr)
            func = getattr(self.fragment, func_name)
            for value in ('3', '+3', '-2', 'plus'):
                self.assertRaises(ValueError, func, value)

    def test_strand_set_ok(self):
        """Test HSPFragment query and hit strand setters"""
        attr = 'strand'
        for seq_type in ('query', 'hit'):
            attr_name = '%s_%s' % (seq_type, attr)
            for value in (-1, 0, 1, None):
                setattr(self.fragment, attr_name, value)
                self.assertEqual(value, getattr(self.fragment, attr_name))

    def test_strand_set_error(self):
        """Test HSPFragment query and hit strand setters, invalid values"""
        attr = 'strand'
        for seq_type in ('query', 'hit'):
            func_name = '_%s_%s_set' % (seq_type, attr)
            func = getattr(self.fragment, func_name)
            for value in (3, 'plus', 'minus', '-', '+'):
                self.assertRaises(ValueError, func, value)

    def test_strand_set_from_plus_frame(self):
        """Test HSPFragment query and hit strand getters, from plus frame"""
        for seq_type in ('query', 'hit'):
            attr_name = '%s_strand' % seq_type
            self.assertTrue(getattr(self.fragment, attr_name) is None)
            setattr(self.fragment, '%s_frame' % seq_type, 3)
            self.assertEqual(1, getattr(self.fragment, attr_name))

    def test_strand_set_from_minus_frame(self):
        """Test HSPFragment query and hit strand getters, from minus frame"""
        for seq_type in ('query', 'hit'):
            attr_name = '%s_strand' % seq_type
            self.assertTrue(getattr(self.fragment, attr_name) is None)
            setattr(self.fragment, '%s_frame' % seq_type, -2)
            self.assertEqual(-1, getattr(self.fragment, attr_name))

    def test_strand_set_from_zero_frame(self):
        """Test HSPFragment query and hit strand getters, from zero frame"""
        for seq_type in ('query', 'hit'):
            attr_name = '%s_strand' % seq_type
            self.assertTrue(getattr(self.fragment, attr_name) is None)
            setattr(self.fragment, '%s_frame' % seq_type, 0)
            self.assertEqual(0, getattr(self.fragment, attr_name))

    def test_coords_setters_getters(self):
        """Test HSPFragment query and hit coordinate-related setters and getters"""
        for seq_type in ('query', 'hit'):
            attr_start = '%s_%s' % (seq_type, 'start')
            attr_end = '%s_%s' % (seq_type, 'end')
            setattr(self.fragment, attr_start, 9)
            setattr(self.fragment, attr_end, 99)
            # check for span value
            span = getattr(self.fragment, '%s_span' % seq_type)
            self.assertEqual(90, span)
            # and range as well
            range = getattr(self.fragment, '%s_range' % seq_type)
            self.assertEqual((9, 99), range)

    def test_coords_setters_readonly(self):
        """Test HSPFragment query and hit coordinate-related read-only getters"""
        read_onlies = ('range', 'span')
        for seq_type in ('query', 'hit'):
            for attr in read_onlies:
                self.assertRaises(AttributeError, setattr,
                                  self.fragment, '%s_%s' % (seq_type, attr), 5)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
