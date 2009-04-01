"""Testing BioSQL with BioSQL

Uses Bio.SeqIO to parse files, and then loads them into a BioSQL database,
and checks we can retreive them again.

Goals:
    Make sure that BioSQL preserves SeqRecord objects.
"""
import os

#TODO - Remove this work around once we drop python 2.3 support
try:
    set = set
except NameError:
    from sets import Set as set

from Bio import MissingExternalDependencyError
from Bio import SeqIO
from Bio.Seq import UnknownSeq
from StringIO import StringIO
from Bio.SeqUtils.CheckSum import seguid
from Bio.SeqFeature import ExactPosition

from BioSQL import BioSeqDatabase
from BioSQL import BioSeq

# This testing suite should try to detect whether a valid database
# installation exists on this computer.  Only run the tests if it
# does.
try :
    from setup_BioSQL import DBDRIVER, DBTYPE
    from setup_BioSQL import DBHOST, DBUSER, DBPASSWD, TESTDB
    from setup_BioSQL import DBSCHEMA, SQL_FILE
except (NameError, ImportError) :
    message = "Check settings in Tests/setup_BioSQL.py "\
              "if you plan to use BioSQL."
    raise MissingExternalDependencyError(message)

db_name = "biosql-seqio-test"
#####################################################################

#This list was based on a selection from test_SeqIO.py
test_files = [ \
#Following nucleic examples are also used in test_Fasta2.py
    ("fasta",  False, 'Nucleic/lupine.nu', 1),
    ("fasta",  False, 'Nucleic/elderberry.nu', 1),
    ("fasta",  False, 'Nucleic/phlox.nu', 1),
    ("fasta",  False, 'Nucleic/centaurea.nu', 1),
    ("fasta",  False, 'Nucleic/wisteria.nu', 1),
    ("fasta",  False, 'Nucleic/sweetpea.nu', 1),
    ("fasta",  False, 'Nucleic/lavender.nu', 1),
#Following protein examples are also used in test_Fasta2.py
    ("fasta",  False, 'Amino/aster.pro', 1),
    ("fasta",  False, 'Amino/loveliesbleeding.pro', 1),
    ("fasta",  False, 'Amino/rose.pro', 1),
    ("fasta",  False, 'Amino/rosemary.pro', 1),
#Following examples are also used in test_Fasta.py
    ("fasta",  False, 'Fasta/f001', 1), #Protein
    ("fasta",  False, 'Fasta/f002', 3), #DNA
    #("fasta", False, 'Fasta/f003', 2), #Protein with comments
    ("fasta",  False, 'Fasta/fa01', 2), #Protein with gaps
#Following examples are also used in test_GFF.py
    ("fasta",  False, 'GFF/NC_001802.fna', 1), #upper case
    ("fasta",  True,  'GFF/multi.fna', 3), #Trivial nucleotide alignment
#Following example is also used in test_registry.py
    ("fasta",  False, 'Registry/seqs.fasta', 2), #contains blank line
#Following examples are also used in test_SwissProt.py
    ("swiss",  False, 'SwissProt/sp001', 1),
    ("swiss",  False, 'SwissProt/sp002', 1),
    ("swiss",  False, 'SwissProt/sp003', 1),
    ("swiss",  False, 'SwissProt/sp004', 1),
    ("swiss",  False, 'SwissProt/sp005', 1),
    ("swiss",  False, 'SwissProt/sp006', 1),
    ("swiss",  False, 'SwissProt/sp007', 1),
    ("swiss",  False, 'SwissProt/sp008', 1),
    ("swiss",  False, 'SwissProt/sp009', 1),
    ("swiss",  False, 'SwissProt/sp010', 1),
    ("swiss",  False, 'SwissProt/sp011', 1),
    ("swiss",  False, 'SwissProt/sp012', 1),
    ("swiss",  False, 'SwissProt/sp013', 1),
    ("swiss",  False, 'SwissProt/sp014', 1),
    ("swiss",  False, 'SwissProt/sp015', 1),
    ("swiss",  False, 'SwissProt/sp016', 1),
#Following example is also used in test_registry.py
    ("swiss",  False, 'Registry/EDD_RAT.dat', 1),
#Following examples are also used in test_GenBank.py
    ("genbank",False, 'GenBank/noref.gb', 1),
    ("genbank",False, 'GenBank/cor6_6.gb', 6),
    ("genbank",False, 'GenBank/iro.gb', 1),
    ("genbank",False, 'GenBank/pri1.gb', 1),
    ("genbank",False, 'GenBank/arab1.gb', 1),
    #protein_refseq.gb had malformed db_xref, fixed in protein_refseq2.gb
    #("genbank",False, 'GenBank/protein_refseq.gb', 1), 
    ("genbank",False, 'GenBank/protein_refseq2.gb', 1),
    ("genbank",False, 'GenBank/extra_keywords.gb', 1),
    ("genbank",False, 'GenBank/one_of.gb', 1),
    ("genbank",False, 'GenBank/NT_019265.gb', 1),
    ("genbank",False, 'GenBank/origin_line.gb', 1),
    ("genbank",False, 'GenBank/blank_seq.gb', 1),
    ("genbank",False, 'GenBank/dbsource_wrap.gb', 1),
    ("genbank",False, 'GenBank/NC_005816.gb', 1),
# The next example is a truncated copy of gbvrl1.seq from
# ftp://ftp.ncbi.nih.gov/genbank/gbvrl1.seq.gz
# This includes an NCBI header, and the first three records:
    ("genbank",False, 'GenBank/gbvrl1_start.seq', 3),
#Following files are also used in test_GFF.py
    ("genbank",False, 'GFF/NC_001422.gbk', 1),
#Following files are currently only used here and test_SeqIO:
    ("embl",      False, 'EMBL/TRBG361.embl', 1),
    ("embl",      False, 'EMBL/DD231055_edited.embl', 1),
    ("embl",      False, 'EMBL/SC10H5.embl', 1), # Pre 2006 style ID line
    ("embl",      False, 'EMBL/U87107.embl', 1), # Old ID line with SV line
    ]

#####################################################################

def checksum_summary(record) :
    if isinstance(record.seq, UnknownSeq) :
        return repr(record.seq)
    if len(record.seq) < 25 :
        short = record.seq.tostring()
    else :
        short = record.seq.tostring()[:19] \
              + "..." + record.seq.tostring()[-3:]
    return "%s [%s] len %i" \
           % (short, seguid(record.seq), len(record.seq))

def compare_references(old_r, new_r) :
    """Compare two Reference objects"""
    assert old_r.title == new_r.title, \
           "%s vs %s" % (old_r.title, new_r.title)
    assert old_r.authors == new_r.authors, \
           "%s vs %s" % (old_r.authors, new_r.authors)
    assert old_r.journal == new_r.journal, \
           "%s vs %s" % (old_r.journal, new_r.journal)
    assert old_r.medline_id == new_r.medline_id, \
           "%s vs %s" % (old_r.medline_id, new_r.medline_id)

    #TODO assert old_r.pubmed_id == new_r.pubmed_id
    #Looking at BioSQL/BioSeq.py function _retrieve_reference
    #it seems that it will get either the MEDLINE or PUBMED,
    #but not both.  I *think* the current schema does not allow
    #us to store both... must confirm this.
    
    #TODO - assert old_r.comment == new_r.comment
    #Looking at the tables, I *think* the current schema does not
    #allow us to store a reference comment.  Must confirm this.
    assert new_r.comment == ""

    #TODO - assert old_r.consrtm == new_r.consrtm
    #Looking at the tables, I *think* the current schema does not
    #allow us to store a consortium.
    assert new_r.consrtm == ""
    
    #TODO - reference location?
    #The parser seems to give a location object (i.e. which
    #nucleotides from the file is the reference for), while the
    #we seem to use the database to hold the journal details (!)

def compare_features(old_f, new_f) :
    """Compare two SeqFeature objects"""

    assert old_f.type == new_f.type, \
        "%s -> %s" % (old_f.type, new_f.type) 
    
    assert old_f.strand == new_f.strand, \
        "%s -> %s" % (old_f.strand, new_f.strand)

    assert old_f.ref == new_f.ref, \
        "%s -> %s" % (old_f.ref, new_f.ref)

    assert old_f.ref_db == new_f.ref_db, \
        "%s -> %s" % (old_f.ref_db, new_f.ref_db)

    #TODO - BioSQL does not store/retrieve feature's id (Bug 2526)   
    #assert old_f.id == new_f.id
    assert new_f.id == "<unknown id>"

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
    
    for old_sub, new_sub in zip(old_f.sub_features, new_f.sub_features) :
        
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
    assert set(old_f.qualifiers.keys()) == set(new_f.qualifiers.keys())
    for key in old_f.qualifiers.keys() :
        if isinstance(old_f.qualifiers[key], str) :
            if isinstance(new_f.qualifiers[key], str) :
                assert old_f.qualifiers[key] == new_f.qualifiers[key]
            elif isinstance(new_f.qualifiers[key], list) :
                #Maybe a string turning into a list of strings?
                assert [old_f.qualifiers[key]] == new_f.qualifiers[key], \
                        "%s -> %s" \
                        % (repr(old_f.qualifiers[key]),
                           repr(new_f.qualifiers[key]))
            else :
                assert False, "Problem with feature's '%s' qualifier" & key
        else :
            #Should both be lists of strings...
            assert old_f.qualifiers[key] == new_f.qualifiers[key], \
                "%s -> %s" % (old_f.qualifiers[key], new_f.qualifiers[key])

def compare_sequences(old, new) :
    """Compare two Seq or DBSeq objects"""
    assert len(old) == len(new)
    assert old.tostring() == new.tostring()

    if isinstance(old, UnknownSeq) :
        assert isinstance(new, UnknownSeq)
    else :
        assert not isinstance(new, UnknownSeq)

    l = len(old)
    s = old.tostring()
    assert isinstance(s, str)

    #Don't check every single element; for long sequences
    #this takes far far far too long to run!
    #Test both positive and negative indices
    if l < 50 :
        indices = range(-l,l)
    else :
        #A selection of end cases, and the mid point
        indices = [-l,-1+1,-int(l/2),-1,0,1,int(l/2),l-2,l-1]

    #Test element access,    
    for i in indices :
        expected = s[i]
        assert expected == old[i]
        assert expected == new[i]

    #Test slices
    indices.append(l) #check copes with overflows
    indices.append(l+1000) #check copes with overflows
    for i in indices :
        for j in indices :
            expected = s[i:j]
            assert expected == old[i:j].tostring(), \
                   "Slice %s vs %s" % (repr(expected), repr(old[i:j]))
            assert expected == new[i:j].tostring(), \
                   "Slice %s vs %s" % (repr(expected), repr(new[i:j]))
            #Slicing with step of 1 should make no difference.
            #Slicing with step 3 might be useful for codons.
            for step in [1,3] :
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
                
def compare_records(old, new) :
    """Compare two SeqRecord or DBSeqRecord objects"""
    #Sequence:
    compare_sequences(old.seq, new.seq)
    #Basics:
    assert old.id == new.id
    assert old.name == new.name
    assert old.description == new.description
    assert old.dbxrefs == new.dbxrefs, \
           "dbxrefs mismatch\nOld: %s\nNew: %s" \
           % (old.dbxrefs, new.dbxrefs)
    #Features:
    assert len(old.features) == len(new.features)
    for old_f, new_f in zip(old.features, new.features) :
        compare_features(old_f, new_f)

    #Annotation:
    #We are expecting to see some "extra" annotations appearing,
    #such as 'cross_references', 'dates', 'data_file_division',
    #'ncbi_taxon' and 'gi'.
    #TODO - address these, see Bug 2681?
    new_keys = set(new.annotations).difference(old.annotations)
    new_keys = new_keys.difference(['cross_references', 'dates', 
                                    'data_file_division', 'ncbi_taxid', 'gi'])
    assert not new_keys, "Unexpected new annotation keys: %s" \
           % ", ".join(new_keys)
    missing_keys = set(old.annotations).difference(new.annotations)
    missing_keys = missing_keys.difference(['ncbi_taxid', # Can't store chimeras
                                            'contig', # See Bug 2745 comments
                                            ])
    assert not missing_keys, "Unexpectedly missing annotation keys: %s" \
           % ", ".join(missing_keys)
    
    #In the short term, just compare any shared keys:
    for key in set(old.annotations.keys()).intersection(new.annotations.keys()) :
        if key == "references" :
            assert len(old.annotations[key]) == len(new.annotations[key])
            for old_r, new_r in zip(old.annotations[key], new.annotations[key]) :
                compare_references(old_r, new_r)
        elif key == "comment":
            if isinstance(old.annotations[key], list):
                old_comment = [comm.replace("\n", " ") for comm in \
                    old.annotations[key]]
            else:
                old_comment = [old.annotations[key].replace("\n", " ")]
            assert len(old_comment) == len(new.annotations[key]), \
                "Number of annotation 'comment's changed by load/retrieve\n" \
                "Was:%s\nNow:%s" \
                % (len(old_comment), len(new.annotations[key]))
            for old_com, new_com in zip(old_comment, new.annotations[key]):
                assert old_com == new_com, \
                    "Annotation 'comment' changed by load/retrieve\n" \
                    "Was:%s\nNow:%s" % (old_com, new_com)
        elif key in ["taxonomy", "organism", "source"]:
            #If there is a taxon id recorded, these fields get overwritten
            #by data from the taxon/taxon_name tables.  There is no
            #guarantee that they will be identical after a load/retrieve.
            assert isinstance(new.annotations[key], str) \
                or isinstance(new.annotations[key], list)
        elif type(old.annotations[key]) == type(new.annotations[key]) :
            assert old.annotations[key] == new.annotations[key], \
                "Annotation '%s' changed by load/retrieve\nWas:%s\nNow:%s" \
                % (key, old.annotations[key], new.annotations[key])
        elif isinstance(old.annotations[key], str) \
        and isinstance(new.annotations[key], list) :
            #Any annotation which is a single string gets turned into
            #a list containing one string by BioSQL at the moment.
            assert [old.annotations[key]] == new.annotations[key], \
                "Annotation '%s' changed by load/retrieve\nWas:%s\nNow:%s" \
                % (key, old.annotations[key], new.annotations[key])
        elif isinstance(old.annotations[key], list) \
        and isinstance(new.annotations[key], str) :
            assert old.annotations[key] == [new.annotations[key]], \
                "Annotation '%s' changed by load/retrieve\nWas:%s\nNow:%s" \
                % (key, old.annotations[key], new.annotations[key])

        
#####################################################################

#TODO - Should we re-use the create_database() function currently
#       defined in test_BioSQL.py here too?  This would allow us
#       to deal with the error of an unknown database...
#
#print "Creating database"
#from setup_BioSQL import create_database
#create_database()

print "Connecting to database"
try :
    server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                      user = DBUSER, passwd = DBPASSWD,
                                      host = DBHOST, db = TESTDB)
except Exception, e :
    message = "Connection failed, check settings in Tests/setup_BioSQL.py "\
              "if you plan to use BioSQL: %s" % str(e)
    raise MissingExternalDependencyError(message)

print "Removing existing sub-database '%s' (if exists)" % db_name
if db_name in server.keys() :
    #Might exist from a failed test run...
    #db = server[db_name]
    server.remove_database(db_name)
    server.commit()

print "(Re)creating empty sub-database '%s'" % db_name
db = server.new_database(db_name)
     
for (t_format, t_alignment, t_filename, t_count) in test_files :
    print "Testing loading from %s format file %s" % (t_format, t_filename)
    assert os.path.isfile(t_filename)

    iterator = SeqIO.parse(handle=open(t_filename,"r"), format=t_format)
    count = db.load(iterator)
    assert count == t_count
    
    #print " - Committing %i records" % count
    server.commit()
    
    iterator = SeqIO.parse(handle=open(t_filename,"r"), format=t_format)
    for record in iterator :
        print " - %s, %s" % (checksum_summary(record), record.id)

        key = record.name
        print " - Retrieving by name/display_id '%s'," % key,
        db_rec = db.lookup(name=key)
        compare_records(record, db_rec)
        db_rec = db.lookup(display_id=key)
        compare_records(record, db_rec)
        print "OK"

        key = record.id
        if key.count(".")==1 and key.split(".")[1].isdigit() :
            print " - Retrieving by version '%s'," % key,
            db_rec = db.lookup(version=key)
            compare_records(record, db_rec)
            print "OK"
        
        if "accessions" in record.annotations :
            accs = set(record.annotations["accessions"])
            for key in accs :
                assert key, "Blank accession in annotation %s" % repr(accs)
                try :
                    print " - Retrieving by accession '%s'," % key,
                    db_rec = db.lookup(accession=key)
                    compare_records(record, db_rec)
                    print "OK"
                except IndexError :
                    print "Failed"
                    pass

        if "gi" in record.annotations :
            key = record.annotations['gi']
            if key != record.id :
                print " - Retrieving by GI '%s'," % key,
                db_rec = db.lookup(primary_id=key)
                compare_records(record, db_rec)
                print "OK"

print "Removing (deleting) '%s'" % db_name
server.remove_database(db_name)

print "Committing remaining changes"
server.commit()

print "Closing connection"
server.close()
