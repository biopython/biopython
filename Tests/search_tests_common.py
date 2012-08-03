# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.SearchIO import _ITERATOR_MAP
from Bio.SeqRecord import SeqRecord


SUPPORTED = _ITERATOR_MAP.keys() + ['mock'] # mock format for other tests


def compare_qresult(qres_a, qres_b, fmt=None):
    """Compares attribute values of two QueryResult objects."""
    # check that file format is supported
    assert fmt in SUPPORTED, "File format %r not supported" % fmt

    # check the number of hits contained
    assert len(qres_a) == len(qres_b), "length: %r vs %r" % (len(qres_a), \
            len(qres_b))

    # set common attributes
    attrs = ['id', 'desc', 'program', 'target', 'version', 'hit_keys', ]
    # and add format-specific attributes
    if fmt == 'blast-xml':
        from Bio.SearchIO.BlastIO.blastxml import _ELEM_QRESULT_OPT, _ELEM_META
        attrs += [x[0] for x in _ELEM_QRESULT_OPT.values()] + \
                [x[0] for x in _ELEM_META.values()] + ['seq_len']
    elif 'blast-tab' in fmt:
        from Bio.SearchIO.BlastIO.blasttab import _COLUMN_QRESULT
        attrs += [x[0] for x in _COLUMN_QRESULT.values()] + ['rid', 'fields']
    elif 'hmmer-' in fmt:
        attrs += ['acc', 'seq_len']
    elif fmt == 'exonerate-text':
        attrs += ['model']

    # compare qresult attributes
    compare_attrs(qres_a, qres_b, attrs)

    # compare hit objects within
    for hit_a, hit_b in zip(qres_a, qres_b):
        assert compare_hit(hit_a, hit_b, fmt)

    return True


def compare_hit(hit_a, hit_b, fmt=None):
    """Compares attribute values of two Hit objects."""
    # check the number of hsps contained
    assert len(hit_a) == len(hit_b), "length: %r vs %r" % (len(hit_a), \
            len(hit_b))

    # set common attributes
    attrs = ['id', 'desc', 'query_id', ]
    # and add format-specific attributes
    if fmt == 'blast-xml':
        from Bio.SearchIO.BlastIO.blastxml import _ELEM_HIT
        attrs += [x[0] for x in _ELEM_HIT.values()]
    elif 'blast-tab' in fmt:
        from Bio.SearchIO.BlastIO.blasttab import _COLUMN_HIT
        attrs += [x[0] for x in _COLUMN_HIT.values()]
    elif fmt == 'hmmer-text':
        attrs += ['evalue', 'bitscore', 'bias', 'domain_exp_num', \
            'domain_obs_num', 'is_in_inclusion', ]
    elif fmt == 'hmmer-tab':
        attrs += ['acc', 'evalue', 'bitscore', 'bias', 'domain_exp_num', \
                'region_num', 'cluster_num', 'overlap_num', 'env_num', \
                'domain_obs_num', 'domain_reported_num', \
                'domain_included_num',]
    elif '-domtab' in fmt:
        attrs += ['id', 'acc', 'seq_len', 'evalue', 'bitscore' ,'bias', ]
    elif 'blat-' in fmt:
        attrs += ['seq_len']

    # compare hit attributes
    compare_attrs(hit_a, hit_b, attrs)

    # compare hit objects within
    for hsp_a, hsp_b in zip(hit_a, hit_b):
        assert compare_hsp(hsp_a, hsp_b, fmt)

    return True


def compare_hsp(hsp_a, hsp_b, fmt=None):
    """Compares attribute values of two HSP objects."""

    attrs = []
    if fmt == 'blast-xml':
        from Bio.SearchIO.BlastIO.blastxml import _ELEM_HSP
        attrs += [x[0] for x in _ELEM_HSP.values()]
    elif fmt == 'blast-tab':
        from Bio.SearchIO.BlastIO.blasttab import _COLUMN_HSP
        attrs += [x[0] for x in _COLUMN_HSP.values()]
    elif fmt == 'hmmer-text':
        attrs += ['domain_index', 'is_in_inclusion', 'bitscore', 'bias', \
                'evalue', 'evalue_cond', 'hit_endtype', 'query_endtype', \
                'env_from', 'env_to', 'env_endtype', 'acc_avg']
    elif fmt == 'hmmer-tab':
        attrs += ['evalue', 'bitscore', 'bias']
    elif '-domtab' in fmt:
        attrs += ['evalue_cond', 'evalue', 'bitscore', 'bias', 'env_start', \
                'env_end', 'acc_avg']
    elif fmt == 'fasta-m10':
        attrs += ['initn_score', 'init1_score', 'opt_score', 'z_score', \
                'bitscore', 'evalue', 'sw_score', 'ident_pct', 'pos_pct', ]
    elif 'blat-' in fmt:
        attrs += ['match_num', 'mismatch_num', 'match_rep_num', 'n_num', \
                'query_gapopen_num', 'query_gap_num', 'hit_gapopen_num', \
                'hit_gap_num', 'query_is_protein', 'ident_pct', 'score', \
                'query_strand', 'hit_strand', 'query_block_spans', \
                'hit_block_spans', 'block_num', 'query_coords', 'hit_coords']
    elif fmt in 'exonerate-':
        attrs += ['query_coords', 'hit_coords', 'query_intron_coords', \
                'hit_intron_coords', 'query_scodon_coords', \
                'hit_scodon_coords', 'query_ner_coords', 'hit_ner_coords']

    # compare hsp attributes
    compare_attrs(hsp_a, hsp_b, attrs)

    # compare fragment objects within
    for frag_a, frag_b in zip(hsp_a, hsp_b):
        assert compare_frag(frag_a, frag_b, fmt)

    return True


def compare_frag(frag_a, frag_b, fmt=None):
    """Compares attribute values of two HSPFragment objects."""

    # set common attributes
    # not comparing alignment attribute, since it's basically comparing
    # hit and query seqs
    attrs = ['hit_id', 'query_id', 'hit', 'query', 'hit_strand', \
            'hit_start', 'hit_end', 'hit_span', 'hit_range', 'query_strand',
            'query_start', 'query_end', 'query_span', 'query_span', 
            'alignment_annotation', ]

    # compare fragment attributes
    assert compare_attrs(frag_a, frag_b, attrs)

    return True


def compare_attrs(obj_a, obj_b, attrs):
    """Compares attribute values of two objects."""
    for attr in attrs:
        try:
            # get attribute values from each objects
            val_a = getattr(obj_a, attr)
            val_b = getattr(obj_b, attr)
        except AttributeError:
            # make sure if the attribute doesn't exist, it doesn't exist
            # in both objects
            # Note that this could cause the test to pass if we pass arbitrary
            # attr names ~ This is OK since we're not checking for attr
            # presence here (it's the job of the individual parser test) ~ we
            # are just checking for attr equality
            assert not hasattr(obj_a, attr)
            assert not hasattr(obj_b, attr)
            # set proxy values for comparison below, which will
            # always return True
            val_a, val_b = None, None

        # special case for HSP and HSPFragment {hit,query}
        # since they are seqrecords, we compare the strings only
        # comparing using compare_record is too slow
        if attr in ('hit', 'query') and (val_a is not None and val_b is \
                not None):
            # compare seq directly if it's a contiguous hsp
            if isinstance(val_a, SeqRecord) and isinstance(val_b, SeqRecord):
                assert str(val_a.seq) == str(val_b.seq), \
                        "%s: %r vs %r" % (attr, val_a, val_b)
            elif isinstance(val_a, list) and isinstance(val_b, list):
                for seq_a, seq_b in zip(val_a, val_b):
                    assert str(seq_a.seq) == str(seq_b.seq), \
                            "%s: %r vs %r" % (attr, seq_a, seq_b)
        # if it's a dictionary, compare values and keys
        elif isinstance(val_a, dict):
            assert isinstance(val_b, dict)
            keys_a, values_a = val_a.keys(), val_a.values()
            keys_b, values_b = val_b.keys(), val_b.values()
            # sort all values and keys
            [x.sort() for x in (keys_a, values_a, keys_b, values_b)]
            assert keys_a == keys_b
            assert values_a == values_b
        else:
            assert val_a == val_b, "%s: %r vs %r" % (attr, val_a, val_b)

    return True
