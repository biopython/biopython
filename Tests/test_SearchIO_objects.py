# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO objects.

Tests the methods and behaviors of Result, Hit, and HSP objects. All tests
are format-independent and are meant to check the fundamental behavior common
to all formats.

"""

import unittest

from Bio.Align import MultipleSeqAlignment
from Bio.SearchIO._objects import Result, Hit, HSP
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


class ResultCases(unittest.TestCase):

    def setUp(self):
        self.result = Result('query1', [hit11, hit21, hit31])

    def test_init_wrong_meta(self):
        # if meta argument is not a dictionary, an exception should be raised
        self.assertRaises(TypeError, Result, 'query_id', meta='test')

    def test_repr(self):
        self.assertEqual("Result(program='<unknown>', target='<unknown>', "
                "id='query1', 3 hits)", repr(self.result))

    def test_iter(self):
        # iteration should return hits contained
        counter = 0
        for hit in self.result:
            self.assertTrue(hit in [hit11, hit21, hit31])
            counter += 1
        self.assertEqual(3, counter)

    def test_hits(self):
        # hits should return hits contained in result
        hits = [x for x in self.result.hits]
        self.assertEqual([hit11, hit21, hit31], hits)

    def test_hit_ids(self):
        # hit_ids should return hit keys (which default to hit ids)
        hit_ids = [x for x in self.result.hit_ids]
        self.assertEqual(['hit1', 'hit2', 'hit3'], hit_ids)

    def test_items(self):
        # items should return tuples of hit key, hit object pair
        items = [x for x in self.result.items]
        self.assertEqual([('hit1', hit11), ('hit2', hit21), \
                ('hit3', hit31)], items)

    def test_contains(self):
        # contains should work with hit ids or hit objects
        self.assertTrue('hit1' in self.result)
        self.assertTrue(hit21 in self.result)
        self.assertFalse('hit5' in self.result)
        self.assertFalse(hit41 in self.result)

    def test_len(self):
        # len() should return the number of hits contained
        self.assertEqual(3, len(self.result))

    def test_nonzero(self):
        # nonzero should return true only if the result has hits
        self.assertTrue(self.result)
        blank_result = Result('queryX')
        self.assertFalse(blank_result)

    def test_reversed(self):
        # reversed should a return a new result with the same attributes
        # except the hits is reversed
        setattr(self.result, 'name', 'test')
        rev_result = reversed(self.result)
        self.assertEqual(self.result.hits[::-1], rev_result.hits[:])
        self.assertEqual('test', rev_result.name)

    def test_setitem_ok(self):
        # hit objects assignment should work with arbitrary string keys
        self.result['hit4'] = hit41
        self.assertEqual([hit11, hit21, hit31, hit41], self.result.hits)
        # and if the key already exist, the object should be overwritten
        self.result['hit4'] = hit11
        self.assertEqual([hit11, hit21, hit31, hit11], self.result.hits)

    def test_setitem_wrong_key_type(self):
        # item assignment should fail if the key is not string
        self.assertRaises(TypeError, self.result.__setitem__, 0, hit41)
        self.assertRaises(TypeError, self.result.__setitem__, slice(0, 2), \
                [hit41, hit31])

    def test_setitem_wrong_type(self):
        # item assignment should fail if the object assigned is not a hit object
        self.assertRaises(TypeError, self.result.__setitem__, 'hit4', hsp111)
        self.assertRaises(TypeError, self.result.__setitem__, 'hit5', 'hit5')

    def test_setitem_wrong_query_id(self):
        # item assignment should fail if the hit object does not have the same
        # query id
        self.assertRaises(ValueError, self.result.__setitem__, 'hit4', hit12)

    def test_getitem_default_ok(self):
        # hits should be retrievable by their keys (default to id)
        self.assertEqual(hit21, self.result['hit2'])
        self.assertEqual(hit11, self.result['hit1'])

    def test_getitem_int_ok(self):
        # hits should be retrievable by their index
        self.assertEqual(hit21, self.result[1])
        self.assertEqual(hit31, self.result[-1])

    def test_getitem_slice_ok(self):
        # if the index is a slice object, a new result object with the same
        # instance attributes should be returned
        setattr(self.result, 'runtime', '3hrs')
        new_result = self.result[1:]
        self.assertEqual([hit21, hit31], new_result.hits)
        self.assertEqual('3hrs', new_result.runtime)

    def test_getitem_tuple_ok(self):
        # hsps should be retrievable if the index is a tuple / list
        self.assertEqual(hsp113, self.result['hit1', 2])
        self.assertEqual([hsp111, hsp112], self.result['hit1', :2])
        self.assertEqual([hsp111, hsp112], self.result['hit1', :-2])
        self.assertEqual([hsp111, hsp113], self.result['hit1', ::2])
        self.assertEqual(hsp113, self.result[0, 2])
        self.assertEqual([hsp111, hsp112], self.result[0, :2])
        self.assertEqual([hsp111, hsp112], self.result[0, :-2])
        self.assertEqual([hsp111, hsp113], self.result[0, ::2])

    def test_getitem_tuple_wrong_length(self):
        # if the tuple index has length < 2, an error should be raised
        self.assertRaises(ValueError, self.result.__getitem__, ('hit1', ))

    def test_getitem_tuple_wrong_type(self):
        # if the tuple's first item is not a string or int, an error
        # should be raised
        self.assertRaises(TypeError, self.result.__getitem__, (slice(2), 0))
        self.assertRaises(TypeError, self.result.__getitem__, (slice(1), \
                slice(2)))

    def test_delitem_string_ok(self):
        # delitem should work with string index
        del self.result['hit1']
        self.assertEqual(2, len(self.result))
        self.assertTrue([hit21, hit31], self.result.hits)

    def test_delitem_int_ok(self):
        # delitem should work with int index
        del self.result[-1]
        self.assertEqual(2, len(self.result))
        self.assertEqual([hit11, hit21], self.result.hits)
        del self.result[0]
        self.assertEqual(1, len(self.result))
        self.assertTrue([hit21], self.result.hits)

    def test_delitem_slice_ok(self):
        # delitem should work with slice objects
        del self.result[:-1]
        self.assertEqual(1, len(self.result))
        self.assertTrue([hit31], self.result.hits)

    def test_delitem_tuple_ok(self):
        # create local mock objects to avoid conflict with other tests
        hit11 = Hit('hit1', 'query1', [hsp111, hsp112, hsp113, hsp114])
        hit21 = Hit('hit2', 'query1', [hsp211])
        hit31 = Hit('hit3', 'query1', [hsp311])
        result = Result('query1', [hit11, hit21, hit31])
        # delitem should work on hsp objects if index is a tuple
        self.assertEqual(3, len(result))
        self.assertEqual(4, len(result[0]))
        del result[0, 1]
        self.assertEqual(3, len(result))
        self.assertEqual([hsp111, hsp113, hsp114], result[0].hsps)
        del result[0, :-1]
        self.assertEqual(3, len(result))
        self.assertEqual([hsp114], result[0].hsps)

    def test_delitem_tuple_wrong_length(self):
        # if the tuple index has length < 2, an error should be raised
        self.assertRaises(ValueError, self.result.__delitem__, ('hit1', ))

    def test_delitem_tuple_wrong_type(self):
        # if the tuple's first item is not a string or int, an error
        # should be raised
        self.assertRaises(TypeError, self.result.__delitem__, (slice(2), 0))
        self.assertRaises(TypeError, self.result.__delitem__, (slice(1), \
                slice(2)))

    def test_append_ok(self):
        # append should work with Hit objects
        self.assertEqual([hit11, hit21, hit31], self.result.hits)
        self.result.append(hit41)
        self.assertEqual([hit11, hit21, hit31, hit41], self.result.hits)
        self.assertEqual(['hit1', 'hit2', 'hit3', 'hit4'], self.result.hit_ids)

    def test_append_custom_hit_key_function_ok(self):
        self.result._hit_key_function = lambda hit: hit.id + '_custom'
        # append should assign hit keys according to _hit_key_function
        self.assertEqual(['hit1', 'hit2', 'hit3'], self.result.hit_ids)
        self.result.append(hit41)
        self.assertEqual(['hit1', 'hit2', 'hit3', 'hit4_custom'], \
                self.result.hit_ids)

    def test_append_id_exists(self):
        # append should raise an error if hit_key already exist
        self.assertRaises(ValueError, self.result.append, hit11)

    def test_pop_ok(self):
        self.assertEqual(3, len(self.result))
        hit = self.result.pop()
        self.assertEqual(hit, hit31)
        self.assertEqual([hit11, hit21], self.result.hits)

    def test_pop_int_index_ok(self):
        # pop should work if given an int index
        self.assertEqual(3, len(self.result))
        hit = self.result.pop(1)
        self.assertEqual(hit, hit21)
        self.assertEqual([hit11, hit31], self.result.hits)

    def test_pop_string_index_ok(self):
        # pop should work if given a string index
        self.assertEqual(3, len(self.result))
        hit = self.result.pop('hit2')
        self.assertEqual(hit, hit21)
        self.assertEqual([hit11, hit31], self.result.hits)

    def test_rank(self):
        # rank should accept hit objects or hit key strings
        self.assertEqual(2, self.result.rank('hit3'))
        self.assertEqual(2, self.result.rank(hit31))

    def test_rank_not_present(self):
        # rank should return -1 if the hit key or hit object is not present
        self.assertEqual(-1, self.result.rank('hit4'))
        self.assertEqual(-1, self.result.rank(hit41))

    def test_sort_default_no_evalue(self):
        # sort without any arguments should keep the Hits in the same order
        # if the hit objects do not have any evalue attributes
        self.assertEqual([hit11, hit21, hit31], self.result.hits)
        self.result.sort()
        self.assertEqual([hit11, hit21, hit31], self.result.hits)

    def test_sort_default_with_evalue(self):
        # create local mock objects to avoid conflict with other tests
        hit11 = Hit('hit1', 'query1', [hsp111, hsp112, hsp113, hsp114])
        hit21 = Hit('hit2', 'query1', [hsp211])
        hit31 = Hit('hit3', 'query1', [hsp311])
        result = Result('query1', [hit11, hit21, hit31])
        # if the hit objects has evalues, sort should be based on that
        result.hits[0].evalue = 1
        result.hits[1].evalue = 10
        result.hits[2].evalue = 1e-10
        self.assertEqual([hit11, hit21, hit31], result.hits)
        result.sort()
        self.assertEqual([hit31, hit11, hit21], result.hits)

    def test_sort_reverse_ok(self):
        # sorting with reverse=True should return a Result with Hits reversed
        self.assertEqual([hit11, hit21, hit31], self.result.hits)
        self.result.sort(reverse=True)
        self.assertEqual([hit31, hit21, hit11], self.result.hits)

    def test_sort_key_ok(self):
        # if custom key is given, sort using it
        key = lambda hit: len(hit)
        self.assertEqual([hit11, hit21, hit31], self.result.hits)
        self.result.sort(key=key)
        self.assertEqual([hit21, hit31, hit11], self.result.hits)



class HitCases(unittest.TestCase):

    def setUp(self):
        self.hit = Hit('hit1', 'query1', [hsp111, hsp112, hsp113])

    def test_init_no_hsps(self):
        # init without hsp objects should fail
        self.assertRaises(ValueError, Hit, 'hit1', 'query1', [])

    def test_repr(self):
        # test for cases with 1 or other alignment numbers
        self.assertEqual("Hit(id='hit1', query_id='query1', 3 alignments)", \
                repr(self.hit))
        self.assertEqual("Hit(id='hit1', query_id='query1', 1 alignment)", \
                repr(self.hit[:1]))

    def test_hsps(self):
        # hsps should return the list of hsps contained
        self.assertEqual([hsp111, hsp112, hsp113], self.hit.hsps)

    def test_iter(self):
        # iteration should return hsps contained
        counter = 0
        for hsp in self.hit:
            self.assertTrue(hsp in [hsp111, hsp112, hsp113])
            counter += 1
        self.assertEqual(3, counter)

    def test_len(self):
        # len() on Hit objects should return how many hsps it has
        self.assertEqual(3, len(self.hit))

    def test_nonzero(self):
        # bool() on Hit objects should return True only if hsps is filled
        # which is always true
        self.assertTrue(self.hit)

    def test_reversed(self):
        rev_hit = reversed(self.hit)
        self.assertEqual(self.hit.hsps[::-1], rev_hit.hsps[:])

    def test_reversed_attrs(self):
        setattr(self.hit, 'evalue', 5e-10)
        setattr(self.hit, 'name', 'test')
        rev_hit = reversed(self.hit)
        self.assertEqual(5e-10, rev_hit.evalue)
        self.assertEqual('test', rev_hit.name)

    def test_setitem_single(self):
        # test regular setitem overwrite
        self.hit[1] = hsp114
        self.assertEqual(self.hit.hsps, [hsp111, hsp114, hsp113])

    def test_item_multiple(self):
        # test iterabie setitem
        self.hit[:] = [hsp113, hsp112, hsp111]
        self.assertEqual(self.hit.hsps, [hsp113, hsp112, hsp111])

    def test_getitem_single(self):
        # getitem using integer index should return a hsp object
        hsp1 = self.hit[0]
        self.assertEqual(hsp111, hsp1)
        hsp3 = self.hit[-1]
        self.assertEqual(hsp113, hsp3)

    def test_getitem_multiple(self):
        # getitem using slices should return another hit object
        # with the hsps sliced accordingly, but other attributes preserved
        new_hit = self.hit[:2]
        self.assertEqual(2, len(new_hit))
        self.assertEqual([hsp111, hsp112], new_hit.hsps)
        self.assertEqual(self.hit.id, new_hit.id)
        self.assertEqual(self.hit.query_id, new_hit.query_id)

    def test_getitem_attrs_multiple(self):
        # check if attributes are carried over
        setattr(self.hit, 'evalue', 1e-5)
        setattr(self.hit, 'name', 'test')
        new_hit = self.hit[:2]
        self.assertEqual(1e-5, new_hit.evalue)
        self.assertEqual('test', new_hit.name)

    def test_delitem(self):
        # test delitem
        del self.hit[0]
        self.assertEqual(2, len(self.hit))
        self.assertEqual([hsp112, hsp113], self.hit.hsps)

    def test_validate_hsp_ok(self):
        # validation should pass if item is an hsp object with matching
        # query and hit ids
        # if validation passes, None is returned
        self.assertEqual(None, self.hit._validate_hsp(hsp114))

    def test_validate_hsp_wrong_type(self):
        # validation should fail if item is not an hsp object
        self.assertRaises(TypeError, self.hit._validate_hsp, 1)
        self.assertRaises(TypeError, self.hit._validate_hsp, Seq(''))

    def test_validate_hsp_wrong_query_id(self):
        # validation should fail if query id does not match
        self.assertRaises(ValueError, self.hit._validate_hsp, hsp211)

    def test_validate_hsp_wrong_hit_id(self):
        # validation should vail if hit id does not match
        self.assertRaises(ValueError, self.hit._validate_hsp, hsp121)

    def test_append(self):
        # append should add hits to the last position
        self.hit.append(hsp114)
        self.assertEqual(4, len(self.hit))
        self.assertEqual(hsp114, self.hit[-1])

    def test_pop(self):
        # pop should return the last item by default
        self.assertEqual(hsp113, self.hit.pop())
        self.assertEqual(hsp111, self.hit.pop(0))

    def test_reverse(self):
        self.assertEqual([hsp111, hsp112, hsp113], self.hit.hsps)
        self.hit.reverse()
        self.assertEqual([hsp113, hsp112, hsp111], self.hit.hsps)

    def test_sort(self):
        self.assertEqual([hsp111, hsp112, hsp113], self.hit.hsps)
        # sort by hsp length
        key = lambda hsp: len(hsp)
        self.hit.sort(key=key)
        self.assertEqual([hsp112, hsp113, hsp111], self.hit.hsps)


class HSPCases(unittest.TestCase):

    def setUp(self):
        self.hsp = HSP('hit_id', 'query_id')

    def test_seq_objects(self):
        # hsp should have no query, hit, and alignment objects
        self.assertTrue(self.hsp.query is None)
        self.assertTrue(self.hsp.hit is None)
        self.assertTrue(self.hsp.alignment is None)

    def test_len(self):
        self.assertRaises(TypeError, len, self.hsp)

    def test_repr(self):
        # test for minimum repr
        self.assertEqual("HSP(hit_id='hit_id', query_id='query_id')", repr(self.hsp))
        # with evalue
        self.hsp.evalue = 1e-25
        self.assertEqual("HSP(hit_id='hit_id', query_id='query_id', evalue=1e-25)", repr(self.hsp))

    def test_getitem(self):
        # getitem not supported without alignment
        self.assertRaises(TypeError, self.hsp.__getitem__, 0)
        self.assertRaises(TypeError, self.hsp.__getitem__, slice(0, 2))

    def test_setitem(self):
        # setitem not supported
        self.assertRaises(TypeError, self.hsp.__setitem__, 0, 'a')
        self.assertRaises(TypeError, self.hsp.__setitem__, slice(0, 2), [1, 2])

    def test_delitem(self):
        # delitem not supported
        self.assertRaises(TypeError, self.hsp.__delitem__, 0)
        self.assertRaises(TypeError, self.hsp.__delitem__, slice(0, 2))

    def test_iter(self):
        # iteration not supported
        self.assertRaises(TypeError, iter, self.hsp)


class HSPWithAlignmentCases(unittest.TestCase):

    def setUp(self):
        self.hsp = HSP('hit_id', 'query_id', 'ATGCTAGCTACA', 'ATG--AGCTAGG')

    def test_init_with_seqrecord(self):
        # init should work with seqrecords
        hit_seq = SeqRecord(Seq('ATGCTAGCTACA'))
        query_seq = SeqRecord(Seq('ATG--AGCTAGG'))
        hsp = HSP('hit_id', 'query_id', hit_seq, query_seq)
        self.assertTrue(isinstance(hsp.query, SeqRecord))
        self.assertTrue(isinstance(hsp.hit, SeqRecord))
        self.assertTrue(isinstance(hsp.alignment, MultipleSeqAlignment))

    def test_init_wrong_seqtypes(self):
        # init should only work with string or seqrecords
        wrong_query = Seq('ATGC')
        wrong_hit = Seq('ATGC')
        self.assertRaises(TypeError, HSP, 'hit_id', 'query_id', wrong_hit, wrong_query)

    def test_seq_objects(self):
        # test for query, hit, and alignment types
        self.assertTrue(isinstance(self.hsp.query, SeqRecord))
        self.assertTrue(isinstance(self.hsp.hit, SeqRecord))
        self.assertTrue(isinstance(self.hsp.alignment, MultipleSeqAlignment))

    def test_len(self):
        # len should equal alignment column length
        self.assertEqual(12, len(self.hsp))

    def test_repr(self):
        # test for minimum repr
        self.assertEqual("HSP(hit_id='hit_id', query_id='query_id', 12-column "
                "alignment)", repr(self.hsp))
        # with evalue
        self.hsp.evalue = 1e-25
        self.assertEqual("HSP(hit_id='hit_id', query_id='query_id', "
                "evalue=1e-25, 12-column alignment)", repr(self.hsp))

    def test_getitem(self):
        # getitem is supported when alignment is present
        sliced_hsp = self.hsp[:5]
        self.assertTrue(isinstance(sliced_hsp, HSP))
        self.assertEqual(5, len(sliced_hsp))
        self.assertEqual('ATGCT', sliced_hsp.hit.seq.tostring())
        self.assertEqual('ATG--', sliced_hsp.query.seq.tostring())

    def test_getitem_attrs(self):
        # attributes from the original instance should be present in the new
        # objects, except for query, hit, and alignment
        setattr(self.hsp, 'attr_original', 1000)
        new_hsp = self.hsp[:5]
        self.assertEqual(1000, new_hsp.attr_original)

    def test_setitem(self):
        # setitem not supported
        self.assertRaises(TypeError, self.hsp.__setitem__, 0, 'a')
        self.assertRaises(TypeError, self.hsp.__setitem__, slice(0, 2), [1, 2])

    def test_delitem(self):
        # delitem not supported
        self.assertRaises(TypeError, self.hsp.__delitem__, 0)
        self.assertRaises(TypeError, self.hsp.__delitem__, slice(0, 2))

    def test_iter(self):
        # iteration not supported
        self.assertRaises(TypeError, iter, self.hsp)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
