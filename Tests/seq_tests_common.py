# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
from Bio.Seq import UnknownSeq
from Bio.SeqUtils.CheckSum import seguid
from Bio.SeqFeature import ExactPosition, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

def checksum_summary(record):
    if isinstance(record.seq, UnknownSeq):
        return repr(record.seq)
    if len(record.seq) < 25:
        short = record.seq.tostring()
    else:
        short = record.seq.tostring()[:19] \
              + "..." + record.seq.tostring()[-3:]
    return "%s [%s] len %i" \
           % (short, seguid(record.seq), len(record.seq))

def compare_reference(old_r, new_r):
    """Compare two Reference objects

    Note new_r is assumed to be a BioSQL DBSeqRecord, due to limitations
    of the BioSQL table structure.
    """
    assert old_r.title == new_r.title, \
           "%s vs %s" % (old_r.title, new_r.title)
    assert old_r.authors == new_r.authors, \
           "%s vs %s" % (old_r.authors, new_r.authors)
    assert old_r.journal == new_r.journal, \
           "%s vs %s" % (old_r.journal, new_r.journal)
    assert old_r.medline_id == new_r.medline_id, \
           "%s vs %s" % (old_r.medline_id, new_r.medline_id)

    if old_r.pubmed_id and new_r.pubmed_id:
        assert old_r.pubmed_id == new_r.pubmed_id
        #Looking at BioSQL/BioSeq.py function _retrieve_reference
        #it seems that it will get either the MEDLINE or PUBMED,
        #but not both.  I *think* the current schema does not allow
        #us to store both... must confirm this.
    
    #TODO - assert old_r.comment == new_r.comment
    #Looking at the tables, I *think* the current schema does not
    #allow us to store a reference comment.  Must confirm this.
    assert old_r.comment == new_r.comment or new_r.comment == "", \
                          "%r vs %r" % (old_r.comment, new_r.comment)

    #TODO - assert old_r.consrtm == new_r.consrtm
    #Looking at the tables, I *think* the current schema does not
    #allow us to store a consortium.
    assert old_r.consrtm == new_r.consrtm or new_r.consrtm == ""
    
    if len(old_r.location) == 0:
        assert len(new_r.location) == 0
    else:
        #BioSQL can only store ONE location!
        #TODO - Check BioPerl with a GenBank file with multiple ref locations
        assert isinstance(old_r.location[0], FeatureLocation)
        assert isinstance(new_r.location[0], FeatureLocation)
        assert old_r.location[0].start == new_r.location[0].start and \
               old_r.location[0].end == new_r.location[0].end

    return True

def compare_feature(old_f, new_f):
    """Compare two SeqFeature objects"""
    assert isinstance(old_f, SeqFeature)
    assert isinstance(new_f, SeqFeature)

    assert old_f.type == new_f.type, \
        "%s -> %s" % (old_f.type, new_f.type) 
    
    assert old_f.strand == new_f.strand, \
        "%s -> %s" % (old_f.strand, new_f.strand)

    assert old_f.ref == new_f.ref, \
        "%s -> %s" % (old_f.ref, new_f.ref)

    assert old_f.ref_db == new_f.ref_db, \
        "%s -> %s" % (old_f.ref_db, new_f.ref_db)

    #TODO - BioSQL does not store/retrieve feature's id (Bug 2526)   
    assert old_f.id == new_f.id or new_f.id == "<unknown id>"

    #TODO - Work out how the location_qualifier_value table should
    #be used, given BioPerl seems to ignore it (Bug 2766)
    #assert old_f.location_operator == new_f.location_operator, \
    #        "%s -> %s" % (old_f.location_operator, new_f.location_operator)
    
    # We dont store fuzzy locations:
    try:
        assert str(old_f.location) == str(new_f.location), \
           "%s -> %s" % (str(old_f.location), str(new_f.location))
    except AssertionError, e:
        if isinstance(old_f.location.start, ExactPosition) and \
            isinstance(old_f.location.end, ExactPosition):
            # Its not a problem with fuzzy locations, re-raise 
            raise e
        else:
            assert old_f.location.nofuzzy_start == \
                    new_f.location.nofuzzy_start, \
                    "%s -> %s" % (old_f.location.nofuzzy_start, \
                                  new_f.location.nofuzzy_start)
            assert old_f.location.nofuzzy_end == \
                    new_f.location.nofuzzy_end, \
                    "%s -> %s" % (old_f.location.nofuzzy_end, \
                                  new_f.location.nofuzzy_end)

    assert len(old_f.sub_features) == len(new_f.sub_features), \
        "number of sub_features: %s -> %s" % \
        (len(old_f.sub_features), len(new_f.sub_features))
    
    for old_sub, new_sub in zip(old_f.sub_features, new_f.sub_features):
        
        assert old_sub.type == new_sub.type, \
            "%s -> %s" % (old_sub.type, new_sub.type)
        
        assert old_sub.strand == new_sub.strand, \
            "%s -> %s" % (old_sub.strand, new_sub.strand)

        assert old_sub.ref == new_sub.ref, \
            "%s -> %s" % (old_sub.ref, new_sub.ref)

        assert old_sub.ref_db == new_sub.ref_db, \
            "%s -> %s" % (old_sub.ref_db, new_sub.ref_db)

        #TODO - Work out how the location_qualifier_value table should
        #be used, given BioPerl seems to ignore it (Bug 2766)
        #assert old_sub.location_operator == new_sub.location_operator, \
        #    "%s -> %s" % (old_sub.location_operator, new_sub.location_operator)

        # Compare sub-feature Locations:
        # 
        # BioSQL currently does not store fuzzy locations, but instead stores
        # them as FeatureLocation.nofuzzy_start FeatureLocation.nofuzzy_end.
        # The vast majority of cases will be comparisons of ExactPosition
        # class locations, so we'll try that first and catch the exceptions.

        try:
            assert str(old_sub.location) == str(new_sub.location), \
               "%s -> %s" % (str(old_sub.location), str(new_sub.location))
        except AssertionError, e:
            if isinstance(old_sub.location.start, ExactPosition) and \
                isinstance(old_sub.location.end, ExactPosition):
                # Its not a problem with fuzzy locations, re-raise 
                raise e
            else:
                #At least one of the locations is fuzzy
                assert old_sub.location.nofuzzy_start == \
                       new_sub.location.nofuzzy_start, \
                       "%s -> %s" % (old_sub.location.nofuzzy_start, \
                                     new_sub.location.nofuzzy_start)
                assert old_sub.location.nofuzzy_end == \
                       new_sub.location.nofuzzy_end, \
                       "%s -> %s" % (old_sub.location.nofuzzy_end, \
                                     new_sub.location.nofuzzy_end)

    assert len(old_f.qualifiers) == len(new_f.qualifiers)    
    assert set(old_f.qualifiers) == set(new_f.qualifiers)
    for key in old_f.qualifiers:
        if isinstance(old_f.qualifiers[key], str):
            if isinstance(new_f.qualifiers[key], str):
                assert old_f.qualifiers[key] == new_f.qualifiers[key]
            elif isinstance(new_f.qualifiers[key], list):
                #Maybe a string turning into a list of strings?
                assert [old_f.qualifiers[key]] == new_f.qualifiers[key], \
                        "%s -> %s" \
                        % (repr(old_f.qualifiers[key]),
                           repr(new_f.qualifiers[key]))
            else:
                assert False, "Problem with feature's '%s' qualifier" & key
        else:
            #Should both be lists of strings...
            assert old_f.qualifiers[key] == new_f.qualifiers[key], \
                "%s -> %s" % (old_f.qualifiers[key], new_f.qualifiers[key])
    return True

def compare_sequence(old, new):
    """Compare two Seq or DBSeq objects"""
    assert len(old) == len(new), "%i vs %i" % (len(old), len(new))
    assert old.tostring() == new.tostring()

    if isinstance(old, UnknownSeq):
        assert isinstance(new, UnknownSeq)
    else:
        assert not isinstance(new, UnknownSeq)

    ln = len(old)
    s = old.tostring()
    assert isinstance(s, str)

    #Don't check every single element; for long sequences
    #this takes far far far too long to run!
    #Test both positive and negative indices
    if ln < 50:
        indices = range(-ln,ln)
    else:
        #A selection of end cases, and the mid point
        indices = [-ln,-ln+1,-(ln//2),-1,0,1,ln//2,ln-2,ln-1]

    #Test element access,    
    for i in indices:
        expected = s[i]
        assert expected == old[i]
        assert expected == new[i]

    #Test slices
    indices.append(ln) #check copes with overflows
    indices.append(ln+1000) #check copes with overflows
    for i in indices:
        for j in indices:
            expected = s[i:j]
            assert expected == old[i:j].tostring(), \
                   "Slice %s vs %s" % (repr(expected), repr(old[i:j]))
            assert expected == new[i:j].tostring(), \
                   "Slice %s vs %s" % (repr(expected), repr(new[i:j]))
            #Slicing with step of 1 should make no difference.
            #Slicing with step 3 might be useful for codons.
            for step in [1,3]:
                expected = s[i:j:step]
                assert expected == old[i:j:step].tostring()
                assert expected == new[i:j:step].tostring()

        #Check automatic end points
        expected = s[i:]
        assert expected == old[i:].tostring()
        assert expected == new[i:].tostring()
                
        expected = s[:i]
        assert expected == old[:i].tostring()
        assert expected == new[:i].tostring()

    #Check "copy" splice
    assert s == old[:].tostring()
    assert s == new[:].tostring()
    return True

def compare_features(old_list, new_list):
    assert isinstance(old_list, list)
    assert isinstance(new_list, list)
    assert len(old_list) == len(new_list)
    for old_f, new_f in zip(old_list, new_list):
        if not compare_feature(old_f, new_f):
            return False
    return True
        
def compare_record(old, new):
    """Compare two SeqRecord or DBSeqRecord objects"""
    assert isinstance(old, SeqRecord)
    assert isinstance(new, SeqRecord)
    #Sequence:
    compare_sequence(old.seq, new.seq)
    #Basics:
    assert old.id == new.id
    assert old.name == new.name
    assert old.description == new.description
    assert old.dbxrefs == new.dbxrefs, \
           "dbxrefs mismatch\nOld: %s\nNew: %s" \
           % (old.dbxrefs, new.dbxrefs)
    #Features:
    if not compare_features(old.features, new.features):
        return False

    #Annotation:
    #We are expecting to see some "extra" annotations appearing,
    #such as 'cross_references', 'dates', 'data_file_division',
    #'ncbi_taxon' and 'gi'.
    #TODO - address these, see Bug 2681?
    new_keys = set(new.annotations).difference(old.annotations)
    new_keys = new_keys.difference(['cross_references', 'date', 
                                    'data_file_division', 'ncbi_taxid', 'gi'])
    assert not new_keys, "Unexpected new annotation keys: %s" \
           % ", ".join(new_keys)
    missing_keys = set(old.annotations).difference(new.annotations)
    missing_keys = missing_keys.difference(['ncbi_taxid', # Can't store chimeras
                                            ])
    assert not missing_keys, "Unexpectedly missing annotation keys: %s" \
           % ", ".join(missing_keys)
    
    #In the short term, just compare any shared keys:
    for key in set(old.annotations).intersection(new.annotations):
        if key == "references":
            assert len(old.annotations[key]) == len(new.annotations[key])
            for old_r, new_r in zip(old.annotations[key], new.annotations[key]):
                compare_reference(old_r, new_r)
        elif key == "comment":
            #Turn them both into containing strings for comparison - due to
            #line wrapping in GenBank etc we don't really expect the white
            #space to be 100% the same.
            if isinstance(old.annotations[key], list):
                old_comment = " ".join(old.annotations[key])
            else:
                old_comment = old.annotations[key]
            if isinstance(new.annotations[key], list):
                new_comment = " ".join(new.annotations[key])
            else:
                new_comment = new.annotations[key]
            old_comment = old_comment.replace("\n"," ").replace("  ", " ")
            new_comment = new_comment.replace("\n"," ").replace("  ", " ")
            assert old_comment == new_comment, \
                "Comment annotation changed by load/retrieve\n" \
                "Was:%s\nNow:%s" \
                % (repr(old_comment), repr(new_comment))
        elif key in ["taxonomy", "organism", "source"]:
            #If there is a taxon id recorded, these fields get overwritten
            #by data from the taxon/taxon_name tables.  There is no
            #guarantee that they will be identical after a load/retrieve.
            assert isinstance(new.annotations[key], basestring) \
                or isinstance(new.annotations[key], list)
        elif type(old.annotations[key]) == type(new.annotations[key]):
            assert old.annotations[key] == new.annotations[key], \
                "Annotation '%s' changed by load/retrieve\nWas:%s\nNow:%s" \
                % (key, old.annotations[key], new.annotations[key])
        elif isinstance(old.annotations[key], str) \
        and isinstance(new.annotations[key], list):
            #Any annotation which is a single string gets turned into
            #a list containing one string by BioSQL at the moment.
            assert [old.annotations[key]] == new.annotations[key], \
                "Annotation '%s' changed by load/retrieve\nWas:%s\nNow:%s" \
                % (key, old.annotations[key], new.annotations[key])
        elif isinstance(old.annotations[key], list) \
        and isinstance(new.annotations[key], str):
            assert old.annotations[key] == [new.annotations[key]], \
                "Annotation '%s' changed by load/retrieve\nWas:%s\nNow:%s" \
                % (key, old.annotations[key], new.annotations[key])
    return True

def compare_records(old_list, new_list):
    assert isinstance(old_list, list)
    assert isinstance(new_list, list)
    assert len(old_list) == len(new_list)
    for old_r, new_r in zip(old_list, new_list):
        if not compare_record(old_r, new_r):
            return False
    return True
