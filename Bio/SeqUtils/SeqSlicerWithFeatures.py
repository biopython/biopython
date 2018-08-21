from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition, BeforePosition, AfterPosition, BetweenPosition,\
    WithinPosition, CompoundLocation
# from Bio import SeqIO


def move_position(pos, move: int, keep_all: bool, slice_len: int):
    tt = type(pos)
    if isinstance(pos, (ExactPosition, BeforePosition, AfterPosition)):
        if int(pos) + move < 0:
            return BeforePosition(0)
        elif int(pos) + move > slice_len:
            return AfterPosition(slice_len)
        else:
            return tt(pos + move)
    elif isinstance(pos, (BetweenPosition, WithinPosition)):
        # This would need update if BetweenPosition or WithinPosition gets obsolete or changed
        # attributes accessed because they are not exposed in any other way
        diff_left = pos._left + move
        diff_right = pos._right + move

        if diff_left < 0 and diff_right <= 0:
            return BeforePosition(0)
        elif diff_left < 0 < diff_right:
            if keep_all:
                tt(0, 0, pos._right)
            else:
                return None
        elif diff_left > 0 and diff_right > 0:
            if diff_left < slice_len and diff_right < slice_len:
                return AfterPosition(slice_len)
            elif diff_left < slice_len < diff_right:
                if keep_all:
                    return tt(slice_len, diff_left, slice_len)
                else:
                    return None
            else:
                # This should not be reached
                raise ValueError
        else:
            # This should not be reached
            raise ValueError
    else:
        # raise for not known Position classes
        raise NotImplementedError('slicing for position {} not implemented'.format(type(pos)))


def start_violation_pos(x, sl: slice):
    return x.start.position < sl.start < x.end.position


def end_violation_pos(x, sl: slice):
    return x.start.position < sl.stop < x.end.position


def slice_feature(feature: FeatureLocation, sl: slice, keep_all_features: bool, seq_len: int):
    if isinstance(feature, CompoundLocation):
        of = []
        for ff in feature.parts:
            # this must be without location boundaries checking
            #  - only location objects adjacent to slice boundaries can cross it

            # However we need to check whether the location object is fully inside the slice
            #  (and discard those which are not)

            mps = move_position(ff.start, - sl.start, keep_all_features, seq_len)
            mpe = move_position(ff.end, - sl.start, keep_all_features, seq_len)
            if mps is not None and mpe is not None:
                fl = FeatureLocation(mps, mpe, ff.strand)

            # Feature locations outside of slice gets discarded
            if fl is not None:
                of.append(fl)

        return CompoundLocation(of, operator=feature.operator)
    else:
        if start_violation_pos(feature, sl) or end_violation_pos(feature, sl):

            mps = move_position(feature.start, - sl.start, keep_all_features, seq_len)
            mpe = move_position(feature.end, - sl.start, keep_all_features, seq_len)

            if mps is not None and mpe is not None:

                fl = FeatureLocation(mps, mpe, feature.strand)
                return fl
            else:
                return None
        else:
            return None


def slice_sequence_with_features(seq: SeqRecord, sl: slice, keep_all_features: bool = False) -> SeqRecord:
    """
    This is an attempt to make slicing function for Bio.Seq which retains features which are not fully in slice
      it is intended for visualization purposes where part of a feature also provides valuable information

    # location variants:
    - CompoundLocation - Collection of FeatureLocation objects (for joins etc).             ? -> solved
    - ExactPosition - Specify the position as being exact.                                  OK
    - WithinPosition - Specify a position occurring within some range.                      ? -> solved
    - BetweenPosition - Specify a position occurring between a range (OBSOLETE?).           ? -> solved
    - BeforePosition - Specify the position as being found before some base.                OK
    - AfterPosition - Specify the position as being found after some base.                  OK
    - OneOfPosition - Specify a position where the location can be multiple positions.      ? Not is spec ? -> raise
    - UnknownPosition - Represents missing information like '?' in UniProt.                 raise

    """
    ns = seq[sl.start:sl.stop]

    seq_len = len(ns)
    for feat in seq.features:
        fl = slice_feature(feat.location, sl, keep_all_features, seq_len)
        if fl is not None:
            nf = SeqFeature(
                fl,
                type=feat.type,
                id=feat.id,
                qualifiers=feat.qualifiers,
                ref=feat.ref,
                ref_db=feat.ref_db
            )
            ns.features.append(nf)

    return ns


if __name__ == '__main__':

    # to my understanding
    # the BetweenPosition should be OK (not DEPRECATED)
    #    - as the between position is only for (n, n+1) - slicing should only retain full feature
    # the within position (usage restricted)
    #    - the within position must be inspected
    # One-of does not exist in specification
    #    - I do not expect it to see it
    # Unknown position not considered (maybe do some basic testing that presence would not break the function)

    # this fails for unknown reason
    # from Bio import SeqIO
    # import os
    # if os.path.exists('../../Tests/GenBank/NC_005816.gb'):
    #     print('ok')

    # using example record from Biopython Tests
    # a = SeqIO.read('../../Tests/GenBank/NC_005816.gb', format='gb')
    # b = slice_sequence_with_features(a, slice(5933, 5935))
    # b = slice_sequence_with_features(a, slice(5932, 5934))
    # b = slice_sequence_with_features(a, slice(5933, 5934))
    # b = slice_sequence_with_features(a, slice(1370, 1670))
    # b = slice_sequence_with_features(a, slice(1430, 1460))
    # b = slice_sequence_with_features(a, slice(1438, 1619))
    # for f in b.features:
    #     print(f)

    print('a')
