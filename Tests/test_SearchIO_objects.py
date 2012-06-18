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

from Bio.Align import MultipleSeqAlignment
from Bio.SearchIO._objects import BaseSearchObject, QueryResult, Hit, HSP
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# create mock objects
hsp111 = HSP('hit1', 'query1', 'ATGCGCAT', 'ATGCGCAT')
hsp112 = HSP('hit1', 'query1', 'ATG', 'GAT')
hsp113 = HSP('hit1', 'query1', 'ATCG', 'CGAT')
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


class BaseSearchObjectCases(unittest.TestCase):

    def setUp(self):
        self.base = BaseSearchObject()

    def test_setattr_string(self):
        """Test BaseSearchObject __setattr__ with string"""
        # string should stay as string if the attribute name is not listed in
        # _objects._{INTS,FLOATS}
        setattr(self.base, 'test', 'attribute')
        self.assertEqual('attribute', self.base.test)

    def test_setattr_int(self):
        """Test BaseSearchObject __setattr__ with integer"""
        # string should be cast into int if the attribute name is listed in
        # _objects._INTS
        setattr(self.base, 'seq_len', '29312')
        self.assertEqual(29312, self.base.seq_len)

    def test_setattr_float(self):
        """Test BaseSearchObject __setattr__ with float"""
        # string should be cast into float if the attribute name is listed in
        # _objects._FLOATS
        setattr(self.base, 'evalue', '22.23')
        self.assertEqual(22.23, self.base.evalue)


class QueryResultCases(unittest.TestCase):

    def setUp(self):
        self.qresult = QueryResult('query1', [hit11, hit21, hit31])

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
        setattr(self.qresult, 'runtime', '3hrs')
        new_qresult = self.qresult[1:]
        self.assertEqual([hit21, hit31], new_qresult.hits)
        self.assertEqual('3hrs', new_qresult.runtime)

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
        # if the hit objects do not have any evalue attributes
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)
        sorted_qresult = self.qresult.sort()
        self.assertEqual([hit11, hit21, hit31], sorted_qresult.hits)
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)

    def test_sort_reverse_ok(self):
        """Test QueryResult.sort, reverse"""
        # sorting with reverse=True should return a QueryResult with Hits reversed
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)
        sorted_qresult = self.qresult.sort(reverse=True)
        self.assertEqual([hit31, hit21, hit11], sorted_qresult.hits)
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)

    def test_sort_key_ok(self):
        """Test QueryResult.sort, with custom key"""
        # if custom key is given, sort using it
        key = lambda hit: len(hit)
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)
        sorted_qresult = self.qresult.sort(key=key)
        self.assertEqual([hit21, hit31, hit11], sorted_qresult.hits)
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)

    def test_sort_in_place_ok(self):
        """Test QueryResult.sort, in place"""
        # sort without any arguments should keep the Hits in the same order
        # if the hit objects do not have any evalue attributes
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)
        self.qresult.sort(in_place=True)
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)

    def test_sort_reverse_in_place_ok(self):
        """Test QueryResult.sort, reverse, in place"""
        # sorting with reverse=True should return a QueryResult with Hits reversed
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)
        self.qresult.sort(reverse=True, in_place=True)
        self.assertEqual([hit31, hit21, hit11], self.qresult.hits)

    def test_sort_key_in_place_ok(self):
        """Test QueryResult.sort, with custom key, in place"""
        # if custom key is given, sort using it
        key = lambda hit: len(hit)
        self.assertEqual([hit11, hit21, hit31], self.qresult.hits)
        self.qresult.sort(key=key, in_place=True)
        self.assertEqual([hit21, hit31, hit11], self.qresult.hits)


class HitCases(unittest.TestCase):

    def setUp(self):
        self.hit = Hit('hit1', 'query1', [hsp111, hsp112, hsp113])

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

    def test_reversed_attrs(self):
        """Test Hit.___reversed__, with attributes"""
        setattr(self.hit, 'evalue', 5e-10)
        setattr(self.hit, 'name', 'test')
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

    def test_getitem_attrs_multiple(self):
        """Test Hit.__getitem__, multiple items, with attributes"""
        # check if attributes are carried over
        setattr(self.hit, 'evalue', 1e-5)
        setattr(self.hit, 'name', 'test')
        new_hit = self.hit[:2]
        self.assertEqual(1e-5, new_hit.evalue)
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

    def test_append(self):
        """Test Hit.append"""
        # append should add hits to the last position
        self.hit.append(hsp114)
        self.assertEqual(4, len(self.hit))
        self.assertEqual(hsp114, self.hit[-1])

    def test_pop(self):
        """Test Hit.pop"""
        # pop should return the last item by default
        self.assertEqual(hsp113, self.hit.pop())
        self.assertEqual(hsp111, self.hit.pop(0))

    def test_reverse(self):
        """Test Hit.reverse"""
        self.assertEqual([hsp111, hsp112, hsp113], self.hit.hsps)
        self.hit.reverse()
        self.assertEqual([hsp113, hsp112, hsp111], self.hit.hsps)

    def test_sort(self):
        """Test Hit.sort"""
        self.assertEqual([hsp111, hsp112, hsp113], self.hit.hsps)
        # sort by hsp length
        key = lambda hsp: len(hsp)
        sorted_hit = self.hit.sort(key=key)
        self.assertEqual([hsp112, hsp113, hsp111], sorted_hit.hsps)
        self.assertEqual([hsp111, hsp112, hsp113], self.hit.hsps)

    def test_sort_in_place(self):
        """Test Hit.sort, in place"""
        self.assertEqual([hsp111, hsp112, hsp113], self.hit.hsps)
        # sort by hsp length
        key = lambda hsp: len(hsp)
        self.hit.sort(key=key, in_place=True)
        self.assertEqual([hsp112, hsp113, hsp111], self.hit.hsps)


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
        """Test HSP sequence attributes"""
        # test for query, hit, and alignment types
        self.assertTrue(isinstance(self.hsp.query, SeqRecord))
        self.assertTrue(isinstance(self.hsp.hit, SeqRecord))
        self.assertTrue(isinstance(self.hsp.alignment, MultipleSeqAlignment))

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
        self.assertEqual('ATGCT', sliced_hsp.hit.seq.tostring())
        self.assertEqual('ATG--', sliced_hsp.query.seq.tostring())

    def test_getitem_attrs(self):
        """Test HSP.__getitem__, with attributes"""
        # attributes from the original instance should not be present in the new
        # objects, except for query, hit, and alignment
        setattr(self.hsp, 'attr_original', 1000)
        self.assertTrue(hasattr(self.hsp, 'attr_original'))
        new_hsp = self.hsp[:5]
        self.assertFalse(hasattr(new_hsp, 'attr_original'))

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

    def test_query_strand_error(self):
        """Test HSP.query_strand getter, error"""
        self.assertRaises(AttributeError, self.hsp._query_strand_get, )

    def test_query_from_smaller(self):
        """Test HSP.query_from getter, smaller from"""
        # from is always smaller
        self.hsp.query_from = 1
        self.hsp.query_to = 10
        self.assertEqual(1, self.hsp.query_from)
        self.assertEqual(10, self.hsp.query_to)

    def test_query_from_bigger(self):
        """Test HSP.query_from getter, bigger from"""
        # from is always smaller
        self.hsp.query_from = 10
        self.hsp.query_to = 1
        self.assertEqual(1, self.hsp.query_from)
        self.assertEqual(10, self.hsp.query_to)

    def test_query_span_ok(self):
        """Test HSP.query_span getter, smaller from"""
        # span is to - from + 1
        self.hsp.query_from = 10
        self.hsp.query_to = 99
        self.assertEqual(90, self.hsp.query_span)

    def test_query_span_ok_2(self):
        """Test HSP.query_span getter, bigger from"""
        # span is to - from + 1
        self.hsp.query_from = 99
        self.hsp.query_to = 10
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

    def test_hit_strand_error(self):
        """Test HSP.hit_strand getter, error"""
        self.assertRaises(AttributeError, self.hsp._hit_strand_get, )

    def test_hit_from_smaller(self):
        """Test HSP.hit_strand getter, smaller from"""
        # from is always smaller
        self.hsp.hit_from = 1
        self.hsp.hit_to = 10
        self.assertEqual(1, self.hsp.hit_from)
        self.assertEqual(10, self.hsp.hit_to)

    def test_hit_from_bigger(self):
        """Test HSP.hit_strand getter, bigger from"""
        # from is always smaller
        self.hsp.hit_from = 10
        self.hsp.hit_to = 1
        self.assertEqual(1, self.hsp.hit_from)
        self.assertEqual(10, self.hsp.hit_to)

    def test_hit_span_ok(self):
        """Test HSP.hit_span getter, smaller from"""
        # span is to - from + 1
        self.hsp.hit_from = 10
        self.hsp.hit_to = 99
        self.assertEqual(90, self.hsp.hit_span)

    def test_hit_span_ok_2(self):
        """Test HSP.hit_span getter, bigger from"""
        # span is to - from + 1
        self.hsp.hit_from = 99
        self.hsp.hit_to = 10
        self.assertEqual(90, self.hsp.hit_span)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
