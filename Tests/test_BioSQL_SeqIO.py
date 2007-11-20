"""Testing BioSQL with BioSQL

Uses Bio.SeqIO to parse files, and then loads them into a BioSQL database,
and checks we can retreive them again.

Goals:
    Make sure that all BioSQL preserves SeqRecord objects.
"""

from sets import Set
import os
from Bio import SeqIO
from StringIO import StringIO
from Bio.SeqUtils.CheckSum import seguid

from BioSQL import BioSeqDatabase
from BioSQL import BioSeq

# This testing suite should try to detect whether a valid database
# installation exists on this computer.  Only run the tests if it
# does.
try :
    from setup_BioSQL import DBDRIVER, DBTYPE
    from setup_BioSQL import DBHOST, DBUSER, DBPASSWD, TESTDB
    from setup_BioSQL import DBSCHEMA, SQL_FILE
except NameError :
    message = "Enable tests in Tests/setup_BioSQL.py (not important if you do not plan to use BioSQL)."
    raise Bio.MissingExternalDependencyError(message)

db_name = "biosql-seqio-test"

#####################################################################

#This list was based on a selection from test_SeqIO.py
test_files = [ \
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
#TODO - Pin down the "Duplicate entry" IntegrityError from this:
#    ("genbank",False, 'GenBank/cor6_6.gb', 6),
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
    if len(record.seq) < 25 :
        short = record.seq.tostring()
    else :
        short = record.seq.tostring()[:19] \
              + "..." + record.seq.tostring()[-3:]
    return "%s [%s] len %i" \
           % (short, seguid(record.seq), len(record.seq))

def compare_records(old, new) :
    #Sequence:
    assert len(old.seq) == len(new.seq)
    assert old.seq.tostring() == new.seq.tostring()
    #Basics:
    assert old.id == new.id
    assert old.name == new.name
    assert old.description == new.description
    #database cross references:
    assert len(old.dbxrefs) == len(new.dbxrefs)
    #Features:
    assert len(old.features) == len(new.features)
    for old_f, new_f in zip(old.features, new.features) :
        assert old_f.type == new_f.type
        #TODO - Sort out zero/None in strand
        assert old_f.strand == new_f.strand \
            or not old_f.strand or not new_f.strand
        assert old_f.id == new_f.id
        #TODO: assert old_f.location_operator == new_f.location_operator
        #TODO: assert str(old_f.location) == str(new_f.location)
        assert len(old_f.sub_features) == len(new_f.sub_features)
        assert len(old_f.qualifiers) == len(new_f.qualifiers)    
        assert Set(old_f.qualifiers.keys()) == Set(new_f.qualifiers.keys())
        for key in old_f.qualifiers.keys() :
            if key == "db_xref":
                #Should be a list
                assert Set(old_f.qualifiers[key]) == Set(new_f.qualifiers[key]), \
                       "db_xref:\n%s\n%s" % (old_f.qualifiers[key], new_f.qualifiers[key])
                pass #TODO - Fix this, e.g. file GenBank/protein_refseq.gb
            else :
                assert old_f.qualifiers[key] == new_f.qualifiers[key]
    #Annotation:
    #TODO assert len(old.annotations) == len(new.annotations) #comments!
    for key in Set(old.annotations.keys()).intersection(new.annotations.keys()) :
        if key == "references" :
            assert len(old.annotations[key]) == len(new.annotations[key])
            for old_r, new_r in zip(old.annotations[key], new.annotations[key]) :
                assert old_r.title == new_r.title
                assert old_r.authors == new_r.authors
                assert old_r.consrtm == new_r.consrtm
                assert old_r.journal == new_r.journal
                assert old_r.medline_id == new_r.medline_id
                #TODO - assert old_r.pubmed_id == new_r.pubmed_id
                #TODO - assert old_r.comment == new_r.comment
                #TODO - reference location?
        else :
            assert old.annotations[key] == new.annotations[key]

        
#####################################################################

print "Connecting to database"
server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                      user = DBUSER, passwd = DBPASSWD,
                                      host = DBHOST, db = TESTDB)

print "Removing existing sub-database '%s' (if exists)" % db_name
if db_name in server.keys() :
    #Might exist from a failed test run...
    #db = server[db_name]
	server.remove_database(db_name)
	server.adaptor.conn.commit()

print "(Re)creating empty sub-database '%s'" % db_name
db = server.new_database(db_name)
     
for (t_format, t_alignment, t_filename, t_count) in test_files :
    print "Testing loading from %s format file %s" % (t_format, t_filename)
    assert os.path.isfile(t_filename)

    iterator = SeqIO.parse(handle=open(t_filename,"r"), format=t_format)
    count = db.load(iterator)
    assert count == t_count
    
    #print " - Committing %i records" % count
    server.adaptor.conn.commit()
    
    iterator = SeqIO.parse(handle=open(t_filename,"r"), format=t_format)
    for record in iterator :
        print " - %s, %s" % (checksum_summary(record), record.id)

        assert len(record.dbxrefs) == 0, "Update this unit test!"
        
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
            accs = Set(record.annotations["accessions"])
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
            if key <> record.id :
                print " - Retrieving by GI '%s'," % key,
                db_rec = db.lookup(primary_id=key)
                compare_records(record, db_rec)
                print "OK"

print "Removing (deleting) '%s'" % db_name
server.remove_database(db_name)

print "Committing remaining changes"
server.adaptor.conn.commit()

print "Closing connection"
server.adaptor.conn.close()
