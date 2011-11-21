# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing BioSQL with BioSQL

Uses Bio.SeqIO to parse files, and then loads them into a BioSQL database,
and checks we can retreive them again.

Goals:
    Make sure that BioSQL preserves SeqRecord objects.
"""
import os

from Bio import MissingExternalDependencyError
from Bio import SeqIO
from Bio.Seq import UnknownSeq
from StringIO import StringIO

from BioSQL import BioSeqDatabase
from BioSQL import BioSeq

# This testing suite should try to detect whether a valid database
# installation exists on this computer.  Only run the tests if it
# does.
try:
    from setup_BioSQL import DBDRIVER, DBTYPE
    from setup_BioSQL import DBHOST, DBUSER, DBPASSWD, TESTDB
    from setup_BioSQL import DBSCHEMA, SQL_FILE
except (NameError, ImportError):
    message = "Check settings in Tests/setup_BioSQL.py "\
              "if you plan to use BioSQL."
    raise MissingExternalDependencyError(message)

from seq_tests_common import checksum_summary, compare_record

db_name = "biosql-seqio-test"
#####################################################################

#This list was based on a selection from test_SeqIO.py
test_files = [ \
#Following nucleic examples are also used in test_SeqIO_FastaIO.py
    ("fasta",  False, 'Fasta/lupine.nu', 1),
    ("fasta",  False, 'Fasta/elderberry.nu', 1),
    ("fasta",  False, 'Fasta/phlox.nu', 1),
    ("fasta",  False, 'Fasta/centaurea.nu', 1),
    ("fasta",  False, 'Fasta/wisteria.nu', 1),
    ("fasta",  False, 'Fasta/sweetpea.nu', 1),
    ("fasta",  False, 'Fasta/lavender.nu', 1),
#Following protein examples are also used in test_SeqIO_FastaIO.py
    ("fasta",  False, 'Fasta/aster.pro', 1),
    ("fasta",  False, 'Fasta/loveliesbleeding.pro', 1),
    ("fasta",  False, 'Fasta/rose.pro', 1),
    ("fasta",  False, 'Fasta/rosemary.pro', 1),
#Following examples are also used in test_SeqIO.py
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

#TODO - Should we re-use the create_database() function currently
#       defined in test_BioSQL.py here too?  This would allow us
#       to deal with the error of an unknown database...
#
#print "Creating database"
#from setup_BioSQL import create_database
#create_database()

print "Connecting to database"
try:
    server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                      user = DBUSER, passwd = DBPASSWD,
                                      host = DBHOST, db = TESTDB)
except Exception, e:
    message = "Connection failed, check settings in Tests/setup_BioSQL.py "\
              "if you plan to use BioSQL: %s" % str(e)
    raise MissingExternalDependencyError(message)

print "Removing existing sub-database '%s' (if exists)" % db_name
if db_name in server.keys():
    #Might exist from a failed test run...
    #db = server[db_name]
    server.remove_database(db_name)
    server.commit()

print "(Re)creating empty sub-database '%s'" % db_name
db = server.new_database(db_name)

db_count = 0
for (t_format, t_alignment, t_filename, t_count) in test_files:
    print "Testing loading from %s format file %s" % (t_format, t_filename)
    assert os.path.isfile(t_filename), t_filename

    iterator = SeqIO.parse(handle=open(t_filename,"r"), format=t_format)
    count = db.load(iterator)
    assert count == t_count
    db_count += count
    
    #print " - Committing %i records" % count
    server.commit()
    
    iterator = SeqIO.parse(handle=open(t_filename,"r"), format=t_format)
    for record in iterator:
        print " - %s, %s" % (checksum_summary(record), record.id)

        key = record.name
        print " - Retrieving by name/display_id '%s'," % key,
        db_rec = db.lookup(name=key)
        compare_record(record, db_rec)
        db_rec = db.lookup(display_id=key)
        compare_record(record, db_rec)
        print "OK"

        key = record.id
        if key.count(".")==1 and key.split(".")[1].isdigit():
            print " - Retrieving by version '%s'," % key,
            db_rec = db.lookup(version=key)
            compare_record(record, db_rec)
            print "OK"
        
        if "accessions" in record.annotations:
            accs = sorted(set(record.annotations["accessions"]))
            for key in accs:
                assert key, "Blank accession in annotation %s" % repr(accs)
                try:
                    print " - Retrieving by accession '%s'," % key,
                    db_rec = db.lookup(accession=key)
                    compare_record(record, db_rec)
                    print "OK"
                except IndexError:
                    print "Failed"
                    pass

        if "gi" in record.annotations:
            key = record.annotations['gi']
            if key != record.id:
                print " - Retrieving by GI '%s'," % key,
                db_rec = db.lookup(primary_id=key)
                compare_record(record, db_rec)
                print "OK"

    assert db_count == len(db), "%i vs %i" % (count, len(db))
    assert db_count == len(db.keys())
    assert db_count == len(db.values())
    assert db_count == len(db.items())

for db_rec in db.itervalues():
    assert isinstance(db_rec, BioSeq.DBSeqRecord)
for key, db_rec in db.iteritems():
    assert isinstance(db_rec, BioSeq.DBSeqRecord)
    compare_record(db_rec, db[key])

print "Removing (deleting) '%s'" % db_name
server.remove_database(db_name)

print "Committing remaining changes"
server.commit()

print "Closing connection"
server.close()
