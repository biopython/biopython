# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO objects.

Tests the methods and behaviors of QueryResult, Hit, and HSP objects. All tests
are format-independent and are meant to check the fundamental behavior common
to all formats.

"""

import unittest
from copy import deepcopy

from search_tests_common import compare_qresult, compare_hit

from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import single_letter_alphabet
from Bio.SearchIO._objects import QueryResult, Hit, HSP
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# create mock objects
hsp111 = HSP('hit1', 'query1', 'ATGCGCAT', 'ATGCGCAT')
hsp112 = HSP('hit1', 'query1', 'ATG', 'GAT')
hsp113 = HSP('hit1', 'query1', 'ATTCG', 'AT-CG')
hsp114 = HSP('hit1', 'query1', 'AT', 'AT')
hsp211 = HSP('hit2', 'query1', 'GGGCCC', 'GGGCC-')
hsp311 = HSP('hit3', 'query1', 'GATG', 'GTTG')
hsp312 = HSP('hit3', 'query1', 'ATATAT', 'ATATAT')
hsp411 = HSP('hit4', 'query1', 'CC-ATG', 'CCCATG')
hsp121 = HSP('hit1', 'query2', 'GCGAG', 'GCGAC')

hit11 = Hit('hit1', 'query1', [hsp111, hsp112, hsp113, hsp114])
hit21 = Hit('hit2', 'query1', [hsp211])
hit31 = Hit('hit3', 'query1', [hsp311, hsp312])
hit41 = Hit('hit4', 'query1', [hsp411])
hit12 = Hit('hit1', 'query2', [hsp121])


class QueryResultCases(unittest.TestCase):

    def setUp(self):
        self.qresult = QueryResult('query1', [hit11, hit21, hit31])
        # set mock attributes
        self.qresult.seq_len = 1102
        self.qresult.target = 'refseq_rna'

    def test_repr(self):
        """Test QueryResult.__repr__"""
        self.assertEqual("QueryResult(id='query1', 3 hits)", \
                repr(self.qresult))

    def test_iter(self):
        """Test QueryResult.__iter__"""
        # iteration should return hits contained
        counter = 0
        for hit in self.qresult:
            self.assertTrue(hit in [hit11, hit21, hit31])
            counter += 1
        self.assertEqual(3, counter)

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
        self.assertEqual([('hit1', hit11), ('hit2', hit21), \
                ('hit3', hit31)], items)

    def test_contains(self):
        """Test QueryResult.__contains__"""
        # contains should work with hit ids or hit objects
        self.assertTrue('hit1' in self.qresult)
        self.assertTrue(hit21 in self.qresult)
        self.assertFalse('hit5' in self.qresult)
        self.assertFalse(hit41 in self.qresult)

    def test_len(self):
        """Test QueryResult.__len__"""
        # len() should return the number of hits contained
        self.assertEqual(3, len(self.qresult))

    def test_nonzero(self):
        """Test QueryResult.__nonzero__"""
        # nonzero should return true only if the qresult has hits
        self.assertTrue(self.qresult)
        blank_qresult = QueryResult('queryX')
        self.assertFalse(blank_qresult)

    def test_reversed(self):
        """Test QueryResult.__reversed__"""
        # reversed should a return a new qresult with the same attributes
        # except the hits is reversed
        setattr(self.qresult, 'name', 'test')
        rev_qresult = reversed(self.qresult)
        self.assertEqual(self.qresult.hits[::-1], rev_qresult.hits[:])
        self.assertEqual('test', rev_qresult.name)

    def test_setitem_ok(self):
        """Test QueryResult.__setitem__"""
        # hit objects assignment should work with arbitrary string keys
        self.qresult['hit4'] = hit41
        self.assertEqual([hit11, hit21, hit31, hit41], self.qresult.hits)
        # and if the key already exist, the object should be overwritten
        self.qresult['hit4'] = hit11
        self.assertEqual([hit11, hit21, hit31, hit11], self.qresult.hits)

    def test_setitem_wrong_key_type(self):
        """Test QueryResult.__setitem__, wrong key type"""
        # item assignment should fail if the key is not string
        self.assertRaises(TypeError, self.qresult.__setitem__, 0, hit41)
        self.assertRaises(TypeError, self.qresult.__setitem__, slice(0, 2), \
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
        self.assertEqual([hit21, hit31], new_qresult.hits)
        self.assertEqual(1102, new_qresult.seq_len)
        self.assertEqual('refseq_rna', new_qresult.target)

    def test_delitem_string_ok(self):
        """Test QueryResult.__getitem__, with string"""
        # delitem should work with string index
        del self.qresult['hit1']
        self.assertEqual(2, len(self.qresult))
        self.assertTrue([hit21, hit31], self.qresult.hits)

    def test_delitem_int_ok(self):
        """Test QueryResult.__delitem__"""
        # delitem should work with int index
        del self.qresult[-1]
        self.assertEqual(2, len(self.qresult))
        self.assertEqual([hit11, hit21], self.qresult.hits)
        del self.qresult[0]
        self.assertEqual(1, len(self.qresult))
        self.assertTrue([hit21], self.qresult.hits)

    def test_delitem_slice_ok(self):
        """Test QueryResult.__delitem__, with slice"""
        # delitem should work with slice objects
        del self.qresult[:-1]
        self.assertEqual(1, len(self.qresult))
        self.assertTrue([hit31], self.qresult.hits)

    def test_desc_set(self):
        """Test QueryResult.description setter"""
        # setting the description should change the query seqrecord description
        # of the contained hsps, if they have an alignment
        # test for default value
        qresult = deepcopy(self.qresult)
        new_desc = 'unicorn hox homolog'
        # test initial condition
        for hit in qresult:
            for hsp in hit:
                self.assertNotEqual(new_desc, hsp.query.description)
        qresult.description = new_desc
        # test after setting
        for hit in qresult:
            for hsp in hit:
                self.assertEqual(new_desc, hsp.query.description)

    def test_desc_set_no_seqrecord(self):
        """Test QueryResult.description setter, without HSP SeqRecords"""
        hsp1 = HSP('hit1', 'query')
        hsp2 = HSP('hit1', 'query')
        hsp3 = HSP('hit2', 'query')
        hit1 = Hit('hit1', 'query', [hsp1, hsp2])
        hit2 = Hit('hit2', 'query', [hsp3])
        qresult = QueryResult('query', [hit1, hit2])
        # test initial condition
        for hit in qresult:
            for hsp in hit:
                self.assertTrue(not hasattr(hsp, 'query'))
        qresult.description = 'unicorn hox homolog'
        # test after setting
        for hit in qresult:
            for hsp in hit:
                self.assertTrue(not hasattr(hsp, 'query'))

    def test_id_set(self):
        """Test QueryResult.id setter"""
        # setting an ID should change the query IDs of all contained Hit and HSPs
        qresult = deepcopy(self.qresult)
        self.assertEqual('query1', qresult.id)
        for hit in qresult:
            self.assertEqual('query1', hit.query_id)
            for hsp in hit:
                self.assertEqual('query1', hsp.query_id)
        qresult.id = 'new_id'
        self.assertEqual('new_id', qresult.id)
        for hit in qresult:
            self.assertEqual('new_id', hit.query_id)
            for hsp in hit:
                self.assertEqual('new_id', hsp.query_id)

    def test_append_ok(self):
        """Test QueryResult.append"""
        # append should work with Hit objects
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)
        self.qresult.append(hit41)
        self.assertEqual([hit11, hit21, hit31, hit41], self.qresult.hits)
        self.assertEqual(['hit1', 'hit2', 'hit3', 'hit4'], self.qresult.hit_keys)

    def test_append_custom_hit_key_function_ok(self):
        """Test QueryResult.append, with custom hit key function"""
        self.qresult._hit_key_function = lambda hit: hit.id + '_custom'
        # append should assign hit keys according to _hit_key_function
        self.assertEqual(['hit1', 'hit2', 'hit3'], self.qresult.hit_keys)
        self.qresult.append(hit41)
        self.assertEqual(['hit1', 'hit2', 'hit3', 'hit4_custom'], \
                self.qresult.hit_keys)

    def test_append_id_exists(self):
        """Test QueryResult.append, when ID exists"""
        # append should raise an error if hit_key already exist
        self.assertRaises(ValueError, self.qresult.append, hit11)

    def test_hit_filter(self):
        """Test QueryResult.hit_filter"""
        # hit_filter should return a new QueryResult object (shallow copy),
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)
        # filter func: min hit length == 2
        # this would filter out hit21, since it only has 1 HSP
        filter_func = lambda hit: len(hit) >= 2
        filtered = self.qresult.hit_filter(filter_func)
        self.assertEqual([hit11, hit31], filtered.hits)
        # make sure all remaining hits return True for the filter function
        self.assertTrue(all([filter_func(hit) for hit in filtered]))
        self.assertEqual(1102, filtered.seq_len)
        self.assertEqual('refseq_rna', filtered.target)

    def test_hit_filter_no_func(self):
        """Test QueryResult.hit_filter, without arguments"""
        # when given no arguments, hit_filter should create a new object with
        # the same contents
        filtered = self.qresult.hit_filter()
        self.assertTrue(compare_qresult(filtered, self.qresult, 'mock'))
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
        self.assertTrue(compare_qresult(mapped, self.qresult, 'mock'))
        self.assertNotEqual(id(mapped), id(self.qresult))
        self.assertEqual(1102, mapped.seq_len)
        self.assertEqual('refseq_rna', mapped.target)

    def test_hsp_filter(self):
        """Test QueryResult.hsp_filter"""
        # hsp_filter should return a new QueryResult object (shallow copy)
        # and any empty hits should be discarded
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)
        # filter func: no '-' in hsp query sequence
        # this would filter out hsp113 and hsp211, effectively removing hit21
        filter_func = lambda hsp: '-' not in str(hsp.query)
        filtered = self.qresult.hsp_filter(filter_func)
        self.assertTrue('hit1' in filtered)
        self.assertTrue('hit2' not in filtered)
        self.assertTrue('hit3' in filtered)
        # test hsps in hit11
        self.assertTrue(all([hsp in filtered['hit1'] for hsp in \
                [hsp111, hsp112, hsp114]]))
        # test hsps in in hit31
        self.assertTrue(all([hsp in filtered['hit3'] for hsp in \
                [hsp311, hsp312]]))

    def test_hsp_filter_no_func(self):
        """Test QueryResult.hsp_filter, no arguments"""
        # when given no arguments, hsp_filter should create a new object with
        # the same contents
        filtered = self.qresult.hsp_filter()
        self.assertTrue(compare_qresult(filtered, self.qresult, 'mock'))
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
        # map func: remove first letter of all HSP.alignment
        def map_func(hsp):
            hsp = hsp[1:]
            return hsp
        mapped = qresult.hsp_map(map_func)
        # make sure old hsp attributes is not transferred to mapped hsps
        for hit in mapped:
            for hsp in hit:
                self.assertFalse(hasattr(hsp, 'mock'))
        # check hsps in hit1
        self.assertEqual('TGCGCAT', str(mapped['hit1'][0].hit.seq))
        self.assertEqual('TGCGCAT', str(mapped['hit1'][0].query.seq))
        self.assertEqual('TG', str(mapped['hit1'][1].hit.seq))
        self.assertEqual('AT', str(mapped['hit1'][1].query.seq))
        self.assertEqual('TTCG', str(mapped['hit1'][2].hit.seq))
        self.assertEqual('T-CG', str(mapped['hit1'][2].query.seq))
        self.assertEqual('T', str(mapped['hit1'][3].hit.seq))
        self.assertEqual('T', str(mapped['hit1'][3].query.seq))
        # check hsps in hit2
        self.assertEqual('GGCCC', str(mapped['hit2'][0].hit.seq))
        self.assertEqual('GGCC-', str(mapped['hit2'][0].query.seq))
        # check hsps in hit3
        self.assertEqual('ATG', str(mapped['hit3'][0].hit.seq))
        self.assertEqual('TTG', str(mapped['hit3'][0].query.seq))
        self.assertEqual('TATAT', str(mapped['hit3'][1].hit.seq))
        self.assertEqual('TATAT', str(mapped['hit3'][1].query.seq))
        # and make sure the attributes are transferred
        self.assertEqual(1102, mapped.seq_len)
        self.assertEqual('refseq_rna', mapped.target)

    def test_hsp_map_no_func(self):
        """Test QueryResult.hsp_map, without arguments"""
        # when given no arguments, hit_map should create a new object with
        # the same contents
        mapped = self.qresult.hsp_map()
        self.assertTrue(compare_qresult(mapped, self.qresult, 'mock'))
        self.assertNotEqual(id(mapped), id(self.qresult))
        self.assertEqual(1102, mapped.seq_len)
        self.assertEqual('refseq_rna', mapped.target)

    def test_pop_ok(self):
        """Test QueryResult.pop"""
        self.assertEqual(3, len(self.qresult))
        hit = self.qresult.pop()
        self.assertEqual(hit, hit31)
        self.assertEqual([hit11, hit21], self.qresult.hits)

    def test_pop_int_index_ok(self):
        """Test QueryResult.pop, with integer index"""
        # pop should work if given an int index
        self.assertEqual(3, len(self.qresult))
        hit = self.qresult.pop(1)
        self.assertEqual(hit, hit21)
        self.assertEqual([hit11, hit31], self.qresult.hits)

    def test_pop_string_index_ok(self):
        """Test QueryResult.pop, with string index"""
        # pop should work if given a string index
        self.assertEqual(3, len(self.qresult))
        hit = self.qresult.pop('hit2')
        self.assertEqual(hit, hit21)
        self.assertEqual([hit11, hit31], self.qresult.hits)

    def test_index(self):
        """Test QueryResult.index"""
        # index should accept hit objects or hit key strings
        self.assertEqual(2, self.qresult.index('hit3'))
        self.assertEqual(2, self.qresult.index(hit31))

    def test_index_not_present(self):
        """Test QueryResult.index, when index is not present"""
        # rank should return -1 if the hit key or hit object is not present
        self.assertEqual(-1, self.qresult.index('hit4'))
        self.assertEqual(-1, self.qresult.index(hit41))

    def test_sort_ok(self):
        """Test QueryResult.sort"""
        # sort without any arguments should keep the Hits in the same order
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)
        self.qresult.sort()
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)

    def test_sort_not_in_place_ok(self):
        """Test QueryResult.sort, not in place"""
        # sort without any arguments should keep the Hits in the same order
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)
        sorted_qresult = self.qresult.sort(in_place=False)
        self.assertEqual([hit11, hit21, hit31], sorted_qresult.hits)
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)

    def test_sort_reverse_ok(self):
        """Test QueryResult.sort, reverse"""
        # sorting with reverse=True should return a QueryResult with Hits reversed
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)
        self.qresult.sort(reverse=True)
        self.assertEqual([hit31, hit21, hit11], self.qresult.hits)

    def test_sort_reverse_not_in_place_ok(self):
        """Test QueryResult.sort, reverse, not in place"""
        # sorting with reverse=True should return a QueryResult with Hits reversed
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)
        sorted_qresult = self.qresult.sort(reverse=True, in_place=False)
        self.assertEqual([hit31, hit21, hit11], sorted_qresult.hits)
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)

    def test_sort_key_ok(self):
        """Test QueryResult.sort, with custom key"""
        # if custom key is given, sort using it
        key = lambda hit: len(hit)
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)
        self.qresult.sort(key=key)
        self.assertEqual([hit21, hit31, hit11], self.qresult.hits)

    def test_sort_key_not_in_place_ok(self):
        """Test QueryResult.sort, with custom key, not in place"""
        # if custom key is given, sort using it
        key = lambda hit: len(hit)
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)
        sorted_qresult = self.qresult.sort(key=key, in_place=False)
        self.assertEqual([hit21, hit31, hit11], sorted_qresult.hits)
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)


class HitCases(unittest.TestCase):

    def setUp(self):
        self.hit = Hit('hit1', 'query1', [hsp111, hsp112, hsp113])
        self.hit.evalue = 5e-10
        self.hit.name = 'test'

    def test_repr(self):
        """Test Hit.__repr__"""
        # test for cases with 1 or other alignment numbers
        self.assertEqual("Hit(id='hit1', query_id='query1', 3 hsps)", \
                repr(self.hit))

    def test_hsps(self):
        """Test Hit.hsps"""
        # hsps should return the list of hsps contained
        self.assertEqual([hsp111, hsp112, hsp113], self.hit.hsps)

    def test_iter(self):
        """Test Hit.__iter__"""
        # iteration should return hsps contained
        counter = 0
        for hsp in self.hit:
            self.assertTrue(hsp in [hsp111, hsp112, hsp113])
            counter += 1
        self.assertEqual(3, counter)

    def test_len(self):
        """Test Hit.__len__"""
        # len() on Hit objects should return how many hsps it has
        self.assertEqual(3, len(self.hit))

    def test_nonzero(self):
        """Test Hit.__nonzero__"""
        # bool() on Hit objects should return True only if hsps is filled
        # which is always true
        self.assertTrue(self.hit)

    def test_reversed(self):
        """Test Hit.__reversed__"""
        rev_hit = reversed(self.hit)
        self.assertEqual(self.hit.hsps[::-1], rev_hit.hsps[:])
        rev_hit = reversed(self.hit)
        self.assertEqual(5e-10, rev_hit.evalue)
        self.assertEqual('test', rev_hit.name)

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
            self.assertNotEqual(new_desc, hsp.hit.description)
        hit.description = new_desc
        # test after setting
        for hsp in hit:
            self.assertEqual(new_desc, hsp.hit.description)

    def test_desc_set_no_seqrecord(self):
        """Test Hit.description setter, without HSP SeqRecords"""
        hsp1 = HSP('hit1', 'query')
        hsp2 = HSP('hit1', 'query')
        hit = Hit('hit1', 'query', [hsp1, hsp2])
        # test initial condition
        for hsp in hit:
            self.assertTrue(not hasattr(hsp, 'hit'))
        hit.description = 'unicorn hox homolog'
        # test after setting
        for hsp in hit:
            self.assertTrue(not hasattr(hsp, 'hit'))

    def test_id_set(self):
        """Test Hit.id setter"""
        # setting an ID should change the query IDs of all contained HSPs
        hit = deepcopy(self.hit)
        self.assertEqual('hit1', hit.id)
        for hsp in hit:
            self.assertEqual('hit1', hsp.hit_id)
        hit.id = 'new_id'
        self.assertEqual('new_id', hit.id)
        for hsp in hit:
            self.assertEqual('new_id', hsp.hit_id)

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
        filter_func = lambda hsp: len(hsp) >= 4
        filtered = self.hit.filter(filter_func)
        self.assertEqual([hsp111, hsp113], filtered.hsps)
        # make sure all remaining hits return True for the filter function
        self.assertTrue(all([filter_func(hit) for hit in filtered]))
        self.assertEqual(5e-10, filtered.evalue)
        self.assertEqual('test', filtered.name)

    def test_filter_no_func(self):
        """Test Hit.filter, without arguments"""
        # when given no arguments, filter should create a new object with
        # the same contents
        filtered = self.hit.filter()
        self.assertTrue(compare_hit(filtered, self.hit, 'mock'))
        self.assertNotEqual(id(filtered), id(self.hit))
        self.assertEqual(5e-10, filtered.evalue)
        self.assertEqual('test', filtered.name)

    def test_filter_no_filtered(self):
        """Test Hit.hit_filter, all hits filtered out"""
        # when the filter filters out all hits, it should return None
        filter_func = lambda hsp: len(hsp) > 50
        filtered = self.hit.filter(filter_func)
        self.assertTrue(filtered is None)

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
            hsp = hsp[1:]
            return hsp
        mapped = hit.map(map_func)
        # make sure old hsp attributes is not transferred to mapped hsps
        for hsp in mapped:
            self.assertFalse(hasattr(hsp, 'mock'))
        # check hsps in hit1
        self.assertEqual('TGCGCAT', str(mapped[0].hit.seq))
        self.assertEqual('TGCGCAT', str(mapped[0].query.seq))
        self.assertEqual('TG', str(mapped[1].hit.seq))
        self.assertEqual('AT', str(mapped[1].query.seq))
        self.assertEqual('TTCG', str(mapped[2].hit.seq))
        self.assertEqual('T-CG', str(mapped[2].query.seq))
        # and make sure the attributes are transferred
        self.assertEqual(5e-10, mapped.evalue)
        self.assertEqual('test', mapped.name)

    def test_hsp_map_no_func(self):
        """Test Hit.map, without arguments"""
        # when given no arguments, map should create a new object with
        # the same contents
        mapped = self.hit.map()
        self.assertTrue(compare_hit(mapped, self.hit, 'mock'))
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
        key = lambda hsp: len(hsp)
        self.hit.sort(key=key)
        self.assertEqual([hsp112, hsp113, hsp111], self.hit.hsps)

    def test_sort_not_in_place(self):
        """Test Hit.sort, not in place"""
        self.assertEqual([hsp111, hsp112, hsp113], self.hit.hsps)
        # sort by hsp length
        key = lambda hsp: len(hsp)
        sorted_hit = self.hit.sort(key=key, in_place=False)
        self.assertEqual([hsp112, hsp113, hsp111], sorted_hit.hsps)
        self.assertEqual([hsp111, hsp112, hsp113], self.hit.hsps)
        self.assertEqual(5e-10, sorted_hit.evalue)
        self.assertEqual('test', sorted_hit.name)


class HSPCases(unittest.TestCase):

    def setUp(self):
        self.hsp = HSP('hit_id', 'query_id')

    def test_seq_objects(self):
        """Test HSP sequence attributes, no alignments"""
        # hsp should have no query, hit, and alignment objects
        self.assertFalse(hasattr(self.hsp, 'query'))
        self.assertFalse(hasattr(self.hsp, 'hit'))
        self.assertFalse(hasattr(self.hsp, 'alignment'))

    def test_len(self):
        """Test HSP.__len__, no alignments"""
        self.assertRaises(TypeError, len, self.hsp)

    def test_repr(self):
        """Test HSP.__repr__, no alignments"""
        # test for minimum repr
        self.assertEqual("HSP(hit_id='hit_id', query_id='query_id')", repr(self.hsp))

    def test_getitem(self):
        """Test HSP.__getitem__, no alignments"""
        # getitem not supported without alignment
        self.assertRaises(TypeError, self.hsp.__getitem__, 0)
        self.assertRaises(TypeError, self.hsp.__getitem__, slice(0, 2))

    def test_setitem(self):
        """Test HSP.__setitem__, no alignments"""
        # setitem not supported
        self.assertRaises(TypeError, self.hsp.__setitem__, 0, 'a')
        self.assertRaises(TypeError, self.hsp.__setitem__, slice(0, 2), [1, 2])

    def test_delitem(self):
        """Test HSP.__delitem__, no alignments"""
        # delitem not supported
        self.assertRaises(TypeError, self.hsp.__delitem__, 0)
        self.assertRaises(TypeError, self.hsp.__delitem__, slice(0, 2))

    def test_iter(self):
        """Test HSP.__iter__, no alignments"""
        # iteration not supported
        self.assertRaises(TypeError, iter, self.hsp)


class HSPWithAlignmentCases(unittest.TestCase):

    def setUp(self):
        self.hsp = HSP('hit_id', 'query_id', 'ATGCTAGCTACA', 'ATG--AGCTAGG')

    def test_init_with_seqrecord(self):
        """Test HSP.__init__, with SeqRecord"""
        # init should work with seqrecords
        hit_seq = SeqRecord(Seq('ATGCTAGCTACA'))
        query_seq = SeqRecord(Seq('ATG--AGCTAGG'))
        hsp = HSP('hit_id', 'query_id', hit_seq, query_seq)
        self.assertTrue(isinstance(hsp.query, SeqRecord))
        self.assertTrue(isinstance(hsp.hit, SeqRecord))
        self.assertTrue(isinstance(hsp.alignment, MultipleSeqAlignment))

    def test_init_wrong_seqtypes(self):
        """Test HSP.__init__, wrong sequence argument types"""
        # init should only work with string or seqrecords
        wrong_query = Seq('ATGC')
        wrong_hit = Seq('ATGC')
        self.assertRaises(TypeError, HSP, 'hit_id', 'query_id', wrong_hit, wrong_query)

    def test_seq_objects(self):
        """Test HSP sequence attribute types and default values"""
        # check hit
        self.assertTrue(isinstance(self.hsp.hit, SeqRecord))
        self.assertEqual('', self.hsp.hit.description)
        self.assertEqual('aligned hit sequence', self.hsp.hit.name)
        self.assertEqual(single_letter_alphabet, self.hsp.hit.seq.alphabet)
        # check query
        self.assertTrue(isinstance(self.hsp.query, SeqRecord))
        self.assertEqual('', self.hsp.query.description)
        self.assertEqual('aligned query sequence', self.hsp.query.name)
        self.assertEqual(single_letter_alphabet, self.hsp.query.seq.alphabet)
        # check alignment
        self.assertTrue(isinstance(self.hsp.alignment, MultipleSeqAlignment))
        self.assertEqual(single_letter_alphabet, self.hsp.alignment._alphabet)

    def test_len(self):
        """Test HSP.__len__"""
        # len should equal alignment column length
        self.assertEqual(12, len(self.hsp))

    def test_repr(self):
        """Test HSP.__repr__"""
        # test for minimum repr
        self.assertEqual("HSP(hit_id='hit_id', query_id='query_id', 12-column "
                "alignment)", repr(self.hsp))

    def test_getitem(self):
        """Test HSP.__getitem__"""
        # getitem is supported when alignment is present
        sliced_hsp = self.hsp[:5]
        self.assertTrue(isinstance(sliced_hsp, HSP))
        self.assertEqual(5, len(sliced_hsp))
        self.assertEqual('ATGCT', str(sliced_hsp.hit.seq))
        self.assertEqual('ATG--', str(sliced_hsp.query.seq))

    def test_getitem_attrs(self):
        """Test HSP.__getitem__, with attributes"""
        # attributes from the original instance should not be present in the new
        # objects, except for query, hit, and alignment
        setattr(self.hsp, 'attr_original', 1000)
        self.assertTrue(hasattr(self.hsp, 'attr_original'))
        new_hsp = self.hsp[:5]
        self.assertFalse(hasattr(new_hsp, 'attr_original'))

    def test_getitem_alignment_annot(self):
        """Test HSP.__getitem__, with alignment annotation"""
        # the alignment is annotated, it should be sliced accordingly
        # and transferred to the new object
        setattr(self.hsp, 'alignment_annotation', {'test': '182718738172'})
        new_hsp = self.hsp[:5]
        self.assertEqual('18271', new_hsp.alignment_annotation['test'])

    def test_setitem(self):
        """Test HSP.__setitem__"""
        # setitem not supported
        self.assertRaises(TypeError, self.hsp.__setitem__, 0, 'a')
        self.assertRaises(TypeError, self.hsp.__setitem__, slice(0, 2), [1, 2])

    def test_delitem(self):
        """Test HSP.__delitem__"""
        # delitem not supported
        self.assertRaises(TypeError, self.hsp.__delitem__, 0)
        self.assertRaises(TypeError, self.hsp.__delitem__, slice(0, 2))

    def test_iter(self):
        """Test HSP.__iter__"""
        # iteration not supported
        self.assertRaises(TypeError, iter, self.hsp)

    def test_default_attrs(self):
        """Test HSP attributes default values"""
        self.assertEqual(None, self.hsp.hit_strand)
        self.assertEqual(None, self.hsp.query_strand)
        self.assertEqual(None, self.hsp.hit_frame)
        self.assertEqual(None, self.hsp.query_frame)

    def test_query_strand_set_ok(self):
        """Test HSP.query_strand setter"""
        # only 1, 0, -1, and None is allowed as strands
        for value in [-1, 0, 1]:
            self.hsp.query_strand = value
            self.assertEqual(value, self.hsp.query_strand)

    def test_query_strand_set_error(self):
        """Test HSP.query_strand setter, wrong values"""
        for value in [3, 'plus', 'minus', '-', '+']:
            self.assertRaises(ValueError, self.hsp._query_strand_set, value)

    def test_query_strand_from_frame_plus(self):
        """Test HSP.query_strand getter, plus from frame"""
        self.hsp.query_frame = 3
        self.assertEqual(1, self.hsp.query_strand)

    def test_query_strand_from_frame_minus(self):
        """Test HSP.query_strand getter, minus from frame"""
        self.hsp.query_frame = -2
        self.assertEqual(-1, self.hsp.query_strand)

    def test_query_range_ok(self):
        """Test HSP.query_range getter"""
        self.hsp.query_start = 9
        self.hsp.query_end = 99
        self.assertEqual((9, 99), self.hsp.query_range)

    def test_query_span_ok(self):
        """Test HSP.query_span getter, smaller from"""
        # span is to - from
        self.hsp.query_start = 9
        self.hsp.query_end = 99
        self.assertEqual(90, self.hsp.query_span)

    def test_hit_strand_set_ok(self):
        """Test HSP.hit_strand setter"""
        # only 1, 0, -1, and None is allowed as strands
        for value in [-1, 0, 1]:
            self.hsp.hit_strand = value
            self.assertEqual(value, self.hsp.hit_strand)

    def test_hit_strand_set_error(self):
        """Test HSP.hit_strand setter, wrong values"""
        for value in [3, 'plus', 'minus', '-', '+']:
            self.assertRaises(ValueError, self.hsp._hit_strand_set, value)

    def test_hit_strand_from_frame_plus(self):
        """Test HSP.hit_strand getter, plus from frame"""
        self.hsp.hit_frame = 3
        self.assertEqual(1, self.hsp.hit_strand)

    def test_hit_strand_from_frame_minus(self):
        """Test HSP.hit_strand getter, minus from frame"""
        self.hsp.hit_frame = -2
        self.assertEqual(-1, self.hsp.hit_strand)

    def test_hit_range_ok(self):
        """Test HSP.hit_range getter"""
        self.hsp.hit_start = 9
        self.hsp.hit_end = 99
        self.assertEqual((9, 99), self.hsp.hit_range)

    def test_hit_span_ok(self):
        """Test HSP.hit_span getter, smaller from"""
        # span is to - from
        self.hsp.hit_start = 9
        self.hsp.hit_end = 99
        self.assertEqual(90, self.hsp.hit_span)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
