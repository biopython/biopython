"""Connect with a BioSQL database and load Biopython like objects from it.

This provides interfaces for loading biological objects from a relational
database, and is compatible with the BioSQL standards.
"""
import BioSeq

def open_database(driver = "MySQLdb", *args, **kwargs):
    """Main interface for loading an existing BioSQL-style database.

    This function is the easiest way to retrieve a connection to a
    database, doing something like:
        
        >>> from BioSeq import BioSeqDatabase
        >>> server = BioSeqDatabase.open_database(user = "root", db="minidb")

    the various options are:
    driver -> The name of the database driver to use for connecting. The
    driver should implement the python DB API. By default, the MySQLdb
    driver is used.
    user -> the username to connect to the database with.
    passwd -> the password to connect with
    host -> the hostname of the database
    db -> the name of the database
    """
    module = __import__(driver)
    connect = getattr(module, "connect")
    conn = connect(*args, **kwargs)
    return DBServer(conn)

class DBServer:
    def __init__(self, conn):
        self.conn = conn
        self.adaptor = Adaptor(self.conn)
    def __repr__(self):
        return self.__class__.__name__ + "(%r)" % self.conn
    def __getitem__(self, name):
        return BioSeqDatabase(self.adaptor, name)
    def keys(self):
        return self.adaptor.list_biodatabase_names()
    def values(self):
        return [self[key] for key in self.keys()]
    def items(self):
        return [(key, self[key]) for key in self.keys()]

class Adaptor:
    def __init__(self, conn):
        self.conn = conn
        self.cursor = conn.cursor()

    def fetch_dbid_by_dbname(self, dbname):
        count = self.cursor.execute(
            r"select biodatabase_id from biodatabase where name = %s",
            (dbname,))
        if count == 0:
            raise KeyError("Cannot find biodatabase with name %r" % dbname)
        assert count == 1, "More than one biodatabase with name %r" % dbname
        return self.cursor.fetchone()[0]

    def fetch_seqid_by_display_id(self, dbid, name):
        count = self.cursor.execute(
            r"select bioentry_id from bioentry where "
            r"    biodatabase_id = %s and display_id = %s",
            (dbid, name))
        if count == 0:
            raise IndexError("Cannot find display id %r" % name)
        assert count == 1, "More than one entry with display id of %r" % name
        seqid, = self.cursor.fetchone()
        return seqid

    def fetch_seqid_by_accession(self, dbid, name):
        count = self.cursor.execute(
            r"select bioentry_id from bioentry where "
            r"    biodatabase_id = %s and accession = %s",
            (dbid, name))
        if count == 0:
            raise IndexError("Cannot find accession %r" % name)
        assert count == 1, "More than one entry with accession of %r" % name
        seqid, = self.cursor.fetchone()
        return seqid

    def fetch_seqid_by_seqid(self, dbid, seqid):
        # XXX can't implement this right since it doesn't seem like the 
        # right id is stored in the database
        raise NotImplementedError("No retrieval by this id")

    def list_biodatabase_names(self):
        self.cursor.execute("select name from biodatabase")
        return [field[0] for field in self.cursor.fetchall()]

    def list_bioentry_ids(self, dbid):
        self.cursor.execute(
            r"select bioentry_id from bioentry where biodatabase_id = %s",
            (dbid,))
        return [field[0] for field in self.cursor.fetchall()]
        
    def list_bioentry_display_ids(self, dbid):
        self.cursor.execute(
            r"select display_ids from bioentry where biodatabase_id = %s",
            (dbid,))
        return [field[0] for field in self.cursor.fetchall()]

    def execute_one(self, sql, args):
        count = self.cursor.execute(sql, args)
        assert count == 1, "Expected 1 response, got %s" % count
        return self.cursor.fetchone()

    def get_subseq_as_string(self, seqid, start, end):
        length = end - start
        return self.execute_one(
            """select SUBSTRING(biosequence_str, %s, %s)
                     from biosequence where bioentry_id = %s""",
            (start+1, length, seqid))[0]

    def execute_and_fetch_col0(self, sql, args):
        self.cursor.execute(sql, args)
        return [field[0] for field in self.cursor.fetchall()]


    def execute_and_fetchall(self, sql, args):
        self.cursor.execute(sql, args)
        return self.cursor.fetchall()
    


_allowed_lookups = {
    # Lookup name / function name to get id, function to list all ids
    'primary_id': "fetch_seqid_by_seqid",
    'display_id': "fetch_seqid_by_display_id",
    'accession': "fetch_seqid_by_accession",
    }


class BioSeqDatabase:
    def __init__(self, adaptor, name):
        self.adaptor = adaptor
        self.name = name
        self.dbid = self.adaptor.fetch_dbid_by_dbname(name)
    def __repr__(self):
        return "BioSeqDatabase(%r, %r)" % (self.adaptor, self.name)
        
    def get_Seq_by_id(self, name):
        """Gets a Bio::Seq object by its name

        Example: seq = db.get_Seq_by_id('ROA1_HUMAN')
        
        """
        seqid = self.adaptor.fetch_seqid_by_display_id(self.dbid, name)
        return BioSeq.DBBioSeq(self.adaptor, seqid)

    def get_Seq_by_acc(self, name):
        """Gets a Bio::Seq object by accession number

        Example: seq = db.get_Seq_by_acc('X77802')

        """
        seqid = self.adaptor.fetch_seqid_by_accession(self.dbid, name)
        return BioSeq.DBBioSeq(self.adaptor, seqid)

    def get_PrimarySeq_stream(self):
        # my @array = $self->get_all_primary_ids;
        # my $stream = Bio::DB::BioDatabasePSeqStream->new(
        #         -adaptor => $self->_adaptor->db->get_PrimarySeqAdaptor,
        #         -idlist => \@array);
        raise NotImplementedError("waiting for Python 2.2's iter")

    def get_all_primary_ids(self):
        """ array of all the primary_ids of the sequences in the database.

        These maybe ids (display style) or accession numbers or
        something else completely different - they *are not*
        meaningful outside of this database implementation.
        """
        return self.adaptor.list_bioentry_ids(self.dbid)

    def __getitem__(self, key):
        return BioSeq.DBBioSeq(self.adaptor, key)
    def keys(self):
        return self.get_all_primary_ids()
    def values(self):
        return [self[key] for key in self.keys()]
    def items(self):
        return [(key, self[key]) for key in self.keys()]

    def lookup(self, **kwargs):
        if len(kwargs) != 1:
            raise TypeError("single key/value parameter expected")
        k, v = kwargs.items()[0]
        if not _allowed_lookups.has_key(k):
            raise TypeError("lookup() expects one of %s, not %r" % \
                            (repr(_allowed_lookups.keys())[1:-1], repr(k)))
        lookup_name = _allowed_lookups[k]
        lookup_func = getattr(self.adaptor, lookup_name)
        seqid = lookup_func(self.dbid, v)
        return BioSeq.DBBioSeq(self.adaptor, seqid)
        
    def get_Seq_by_primary_id(self, seqid):
        """Gets a Bio::Seq object by the primary (internal) id.

        The primary id in these cases has to come from
        $db->get_all_primary_ids.  There is no other way to get (or
        guess) the primary_ids in a database.
        """
        return self[seqid]
