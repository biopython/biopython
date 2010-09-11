# Copyright 2002 by Andrew Dalke.  All rights reserved.
# Revisions 2007-2009 copyright by Peter Cock.  All rights reserved.
# Revisions 2009 copyright by Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Note that BioSQL (including the database schema and scripts) is
# available and licensed separately.  Please consult www.biosql.org
"""Connect with a BioSQL database and load Biopython like objects from it.

This provides interfaces for loading biological objects from a relational
database, and is compatible with the BioSQL standards.
"""
import BioSeq
import Loader
import DBUtils

_POSTGRES_RULES_PRESENT = False # Hack for BioSQL Bug 2839

def open_database(driver = "MySQLdb", **kwargs):
    """Main interface for loading a existing BioSQL-style database.

    This function is the easiest way to retrieve a connection to a
    database, doing something like:
        
        >>> from BioSeq import BioSeqDatabase
        >>> server = BioSeqDatabase.open_database(user="root", db="minidb")

    the various options are:
    driver -> The name of the database driver to use for connecting. The
    driver should implement the python DB API. By default, the MySQLdb
    driver is used.
    user -> the username to connect to the database with.
    password, passwd -> the password to connect with
    host -> the hostname of the database
    database or db -> the name of the database
    """
    if driver == "psycopg":
        raise ValueError("Using BioSQL with psycopg (version one) is no "
                         "longer supported. Use psycopg2 instead.")

    module = __import__(driver)
    connect = getattr(module, "connect")

    # Different drivers use different keywords...
    kw = kwargs.copy()
    if driver == "MySQLdb":
        if "database" in kw:
            kw["db"] = kw["database"]
            del kw["database"]
        if "password" in kw:
            kw["passwd"] = kw["password"]
            del kw["password"]
    else:
        # DB-API recommendations
        if "db" in kw:
            kw["database"] = kw["db"]
            del kw["db"]
        if "passwd" in kw:
            kw["password"] = kw["passwd"]
            del kw["passwd"]
    if driver in ["psycopg2", "pgdb"] and not kw.get("database"):
        kw["database"] = "template1"
    # SQLite connect takes the database name as input
    if driver in ["sqlite3"]:
        conn = connect(kw["database"])
    else:
        try:
            conn = connect(**kw)
        except module.InterfaceError:
            # Ok, so let's try building a DSN
            # (older releases of psycopg need this)
            if "database" in kw:
                kw["dbname"] = kw["database"]
                del kw["database"]
            elif "db" in kw:
                kw["dbname"] = kw["db"]
                del kw["db"]
            dsn = ' '.join(['='.join(i) for i in kw.items()])
            conn = connect(dsn)

    server = DBServer(conn, module)

    # TODO - Remove the following once BioSQL Bug 2839 is fixed.
    # Test for RULES in PostgreSQL schema, see also Bug 2833.
    if driver in ["psycopg2", "pgdb"]:
        sql = "SELECT ev_class FROM pg_rewrite WHERE " + \
              "rulename='rule_bioentry_i1' OR " + \
              "rulename='rule_bioentry_i2';"
        if server.adaptor.execute_and_fetchall(sql):
            import warnings
            warnings.warn("Your BioSQL PostgreSQL schema includes some "
                          "rules currently required for bioperl-db but "
                          "which may cause problems loading data using "
                          "Biopython (see BioSQL Bug 2839). If you do not "
                          "use BioPerl, please remove these rules. "
                          "Biopython should cope with the rules present, "
                          "but with a performance penalty when loading "
                          "new records.")
            global _POSTGRES_RULES_PRESENT
            _POSTGRES_RULES_PRESENT = True

    return server

class DBServer:
    """Represents a BioSQL database continaing namespaces (sub-databases).
    
    This acts like a Python dictionary, giving access to each namespace
    (defined by a row in the biodatabase table) as a BioSeqDatabase object.
    """
    def __init__(self, conn, module, module_name=None):
        self.module = module
        if module_name is None:
            module_name = module.__name__
        self.adaptor = Adaptor(conn, DBUtils.get_dbutils(module_name))
        self.module_name = module_name
        
    def __repr__(self):
        return self.__class__.__name__ + "(%r)" % self.adaptor.conn

    def __getitem__(self, name):
        return BioSeqDatabase(self.adaptor, name)

    def __len__(self):
        """Number of namespaces (sub-databases) in this database."""
        sql = "SELECT COUNT(name) FROM biodatabase;"
        return int(self.adaptor.execute_and_fetch_col0(sql)[0])

    def __contains__(self, value):
        """Check if a namespace (sub-database) in this database."""
        sql = "SELECT COUNT(name) FROM biodatabase WHERE name=%s;"
        return bool(self.adaptor.execute_and_fetch_col0(sql, (value,))[0])
    
    def __iter__(self):
        """Iterate over namespaces (sub-databases) in the database."""
        #TODO - Iterate over the cursor, much more efficient
        return iter(self.adaptor.list_biodatabase_names())        

    if hasattr(dict, "iteritems"):
        #Python 2, use iteritems etc    
        def keys(self):
            """List of namespaces (sub-databases) in the database."""
            return self.adaptor.list_biodatabase_names()

        def values(self):
            """List of BioSeqDatabase objects in the database."""
            return [self[key] for key in self.keys()]
    
        def items(self):
            """List of (namespace, BioSeqDatabase) for entries in the database."""
            return [(key, self[key]) for key in self.keys()]
        
        def iterkeys(self):
            """Iterate over namespaces (sub-databases) in the database."""
            return iter(self)
    
        def itervalues(self):
            """Iterate over BioSeqDatabase objects in the database."""
            for key in self:
                yield self[key]
            
        def iteritems(self):
            """Iterate over (namespace, BioSeqDatabase) in the database."""
            for key in self:
                yield key, self[key]
    else:
        #Python 3, items etc are all iterators
        def keys(self):
            """Iterate over namespaces (sub-databases) in the database."""
            return iter(self)
            
        def values(self):
            """Iterate over BioSeqDatabase objects in the database."""
            for key in self:
                yield self[key]
    
        def items(self):
            """Iterate over (namespace, BioSeqDatabase) in the database."""
            for key in self:
                yield key, self[key]

    def __delitem__(self, name):
        """Remove a namespace and all its entries."""
        if name not in self:
            raise KeyError(name)
        self.remove_database(name)

    def remove_database(self, db_name):
        """Remove a namespace and all its entries (OBSOLETE).
        
        Try to remove all references to items in a database.
        
        server.remove_database(name)
        
        In keeping with the dictionary interface, you can now do this:
        
        del server[name]
        """
        import warnings
        warnings.warn("This method is obsolete.  In keeping with the dictionary interface, you can now use 'del server[name]' instead", PendingDeprecationWarning)
        db_id = self.adaptor.fetch_dbid_by_dbname(db_name)
        remover = Loader.DatabaseRemover(self.adaptor, db_id)
        remover.remove()

    def new_database(self, db_name, authority=None, description=None):
        """Add a new database to the server and return it.
        """
        # make the database
        sql = r"INSERT INTO biodatabase (name, authority, description)" \
              r" VALUES (%s, %s, %s)" 
        self.adaptor.execute(sql, (db_name,authority, description))
        return BioSeqDatabase(self.adaptor, db_name)

    def load_database_sql(self, sql_file):
        """Load a database schema into the given database.

        This is used to create tables, etc when a database is first created.
        sql_file should specify the complete path to a file containing
        SQL entries for building the tables.
        """
        # Not sophisticated enough for PG schema. Is it needed by MySQL?
        # Looks like we need this more complicated way for both. Leaving it
        # the default and removing the simple-minded approach.

        # read the file with all comment lines removed
        sql_handle = open(sql_file, "rU")
        sql = r""
        for line in sql_handle:
            if line.find("--") == 0: # don't include comment lines
                pass
            elif line.find("#") == 0: # ditto for MySQL comments
                pass
            elif line.strip(): # only include non-blank lines
                sql += line.strip()
                sql += ' '
        
        # two ways to load the SQL
        # 1. PostgreSQL can load it all at once and actually needs to
        # due to FUNCTION defines at the end of the SQL which mess up
        # the splitting by semicolons
        if self.module_name in ["psycopg2", "pgdb"]:
            self.adaptor.cursor.execute(sql)
        # 2. MySQL needs the database loading split up into single lines of
        # SQL executed one at a time
        elif self.module_name in ["MySQLdb", "sqlite3"]:
            sql_parts = sql.split(";") # one line per sql command
            for sql_line in sql_parts[:-1]: # don't use the last item, it's blank
                self.adaptor.cursor.execute(sql_line)
        else:
            raise ValueError("Module %s not supported by the loader." %
                    (self.module_name))

    def commit(self):
        """Commits the current transaction to the database."""
        return self.adaptor.commit()

    def rollback(self):
        """Rolls backs the current transaction."""
        return self.adaptor.rollback()

    def close(self):
        """Close the connection. No further activity possible."""
        return self.adaptor.close()

class Adaptor:
    def __init__(self, conn, dbutils):
        self.conn = conn
        self.cursor = conn.cursor()
        self.dbutils = dbutils

    def last_id(self, table):
        return self.dbutils.last_id(self.cursor, table)

    def autocommit(self, y=True):
        """Set the autocommit mode. True values enable; False value disable."""
        return self.dbutils.autocommit(self.conn, y)

    def commit(self):
        """Commits the current transaction."""
        return self.conn.commit()

    def rollback(self):
        """Rolls backs the current transaction."""
        return self.conn.rollback()

    def close(self):
        """Close the connection. No further activity possible."""
        return self.conn.close()

    def fetch_dbid_by_dbname(self, dbname):
        self.execute(
            r"select biodatabase_id from biodatabase where name = %s",
            (dbname,))
        rv = self.cursor.fetchall()
        if not rv:
            raise KeyError("Cannot find biodatabase with name %r" % dbname)
        # Cannot happen (UK)
##        assert len(rv) == 1, "More than one biodatabase with name %r" % dbname
        return rv[0][0]

    def fetch_seqid_by_display_id(self, dbid, name):
        sql = r"select bioentry_id from bioentry where name = %s"
        fields = [name]
        if dbid:
            sql += " and biodatabase_id = %s"
            fields.append(dbid)
        self.execute(sql, fields)
        rv = self.cursor.fetchall()
        if not rv:
            raise IndexError("Cannot find display id %r" % name)
        if len(rv) > 1:
            raise IndexError("More than one entry with display id %r" % name)
        return rv[0][0]

    def fetch_seqid_by_accession(self, dbid, name):
        sql = r"select bioentry_id from bioentry where accession = %s"
        fields = [name]
        if dbid:
            sql += " and biodatabase_id = %s"
            fields.append(dbid)
        self.execute(sql, fields)
        rv = self.cursor.fetchall()
        if not rv:
            raise IndexError("Cannot find accession %r" % name)
        if len(rv) > 1:
            raise IndexError("More than one entry with accession %r" % name)
        return rv[0][0]

    def fetch_seqids_by_accession(self, dbid, name):
        sql = r"select bioentry_id from bioentry where accession = %s"
        fields = [name]
        if dbid:
            sql += " and biodatabase_id = %s"
            fields.append(dbid)
        return self.execute_and_fetch_col0(sql, fields)

    def fetch_seqid_by_version(self, dbid, name):
        acc_version = name.split(".")
        if len(acc_version) > 2:
            raise IndexError("Bad version %r" % name)
        acc = acc_version[0]
        if len(acc_version) == 2:
            version = acc_version[1]
        else:
            version = "0"
        sql = r"SELECT bioentry_id FROM bioentry WHERE accession = %s" \
              r" AND version = %s"
        fields = [acc, version]
        if dbid:
            sql += " and biodatabase_id = %s"
            fields.append(dbid)
        self.execute(sql, fields)
        rv = self.cursor.fetchall()
        if not rv:
            raise IndexError("Cannot find version %r" % name)
        if len(rv) > 1:
            raise IndexError("More than one entry with version %r" % name)
        return rv[0][0]

    def fetch_seqid_by_identifier(self, dbid, identifier):
        # YB: was fetch_seqid_by_seqid
        sql = "SELECT bioentry_id FROM bioentry WHERE identifier = %s"
        fields = [identifier]
        if dbid:
            sql += " and biodatabase_id = %s"
            fields.append(dbid)
        self.execute(sql, fields)
        rv = self.cursor.fetchall()
        if not rv:
            raise IndexError("Cannot find display id %r" % identifier)
        return rv[0][0]

    def list_biodatabase_names(self):
        return self.execute_and_fetch_col0(
            "SELECT name FROM biodatabase")

    def list_bioentry_ids(self, dbid):
        return self.execute_and_fetch_col0(
            "SELECT bioentry_id FROM bioentry WHERE biodatabase_id = %s",
            (dbid,))

    def list_bioentry_display_ids(self, dbid):
        return self.execute_and_fetch_col0(
            "SELECT name FROM bioentry WHERE biodatabase_id = %s",
            (dbid,))

    def list_any_ids(self, sql, args):
        """Return ids given a SQL statement to select for them.
        
        This assumes that the given SQL does a SELECT statement that
        returns a list of items. This parses them out of the 2D list
        they come as and just returns them in a list.
        """
        return self.execute_and_fetch_col0(sql, args)

    def execute_one(self, sql, args=None):
        self.execute(sql, args or ())
        rv = self.cursor.fetchall()
        assert len(rv) == 1, "Expected 1 response, got %d" % len(rv)
        return rv[0]

    def execute(self, sql, args=None):
        """Just execute an sql command.
        """
        self.dbutils.execute(self.cursor, sql, args)

    def get_subseq_as_string(self, seqid, start, end):
        length = end - start
        # XXX Check this on MySQL and PostgreSQL. substr should be general,
        # does it need dbutils?
        #return self.execute_one(
        #    """select SUBSTRING(seq FROM %s FOR %s)
        #             from biosequence where bioentry_id = %s""",
        #    (start+1, length, seqid))[0]
        # 
        # Convert to a string on returning for databases that give back
        # unicode. Shouldn't need unicode for sequences so this seems safe.
        return str(self.execute_one(
            """select SUBSTR(seq, %s, %s)
                     from biosequence where bioentry_id = %s""",
            (start+1, length, seqid))[0])

    def execute_and_fetch_col0(self, sql, args=None):
        self.execute(sql, args or ())
        return [field[0] for field in self.cursor.fetchall()]

    def execute_and_fetchall(self, sql, args=None):
        self.execute(sql, args or ())
        return self.cursor.fetchall()

_allowed_lookups = {
    # Lookup name / function name to get id, function to list all ids
    'primary_id': "fetch_seqid_by_identifier",
    'gi':         "fetch_seqid_by_identifier",
    'display_id': "fetch_seqid_by_display_id",
    'name':       "fetch_seqid_by_display_id",
    'accession':  "fetch_seqid_by_accession",
    'version':    "fetch_seqid_by_version",
    }

class BioSeqDatabase:
    """Represents a namespace (sub-database) within the BioSQL database.
    
    i.e. One row in the biodatabase table, and all all rows in the bioentry
    table associated with it.
    """
    def __init__(self, adaptor, name):
        self.adaptor = adaptor
        self.name = name
        self.dbid = self.adaptor.fetch_dbid_by_dbname(name)

    def __repr__(self):
        return "BioSeqDatabase(%r, %r)" % (self.adaptor, self.name)
        
    def get_Seq_by_id(self, name):
        """Gets a DBSeqRecord object by its name

        Example: seq_rec = db.get_Seq_by_id('ROA1_HUMAN')
        
        The name of this method is misleading since it returns a DBSeqRecord
        rather than a DBSeq ojbect, and presumably was to mirror BioPerl.
        """
        seqid = self.adaptor.fetch_seqid_by_display_id(self.dbid, name)
        return BioSeq.DBSeqRecord(self.adaptor, seqid)

    def get_Seq_by_acc(self, name):
        """Gets a DBSeqRecord object by accession number

        Example: seq_rec = db.get_Seq_by_acc('X77802')

        The name of this method is misleading since it returns a DBSeqRecord
        rather than a DBSeq ojbect, and presumably was to mirror BioPerl.
        """
        seqid = self.adaptor.fetch_seqid_by_accession(self.dbid, name)
        return BioSeq.DBSeqRecord(self.adaptor, seqid)

    def get_Seq_by_ver(self, name):
        """Gets a DBSeqRecord object by version number

        Example: seq_rec = db.get_Seq_by_ver('X77802.1')

        The name of this method is misleading since it returns a DBSeqRecord
        rather than a DBSeq ojbect, and presumably was to mirror BioPerl.
        """
        seqid = self.adaptor.fetch_seqid_by_version(self.dbid, name)
        return BioSeq.DBSeqRecord(self.adaptor, seqid)

    def get_Seqs_by_acc(self, name):
        """Gets a list of DBSeqRecord objects by accession number

        Example: seq_recs = db.get_Seq_by_acc('X77802')

        The name of this method is misleading since it returns a list of
        DBSeqRecord objects rather than a list of DBSeq ojbects, and presumably
        was to mirror BioPerl.
        """
        seqids = self.adaptor.fetch_seqids_by_accession(self.dbid, name)
        return [BioSeq.DBSeqRecord(self.adaptor, seqid) for seqid in seqids]

    def get_all_primary_ids(self):
        """All the primary_ids of the sequences in the database (OBSOLETE).

        These maybe ids (display style) or accession numbers or
        something else completely different - they *are not*
        meaningful outside of this database implementation.
        
        Please use .keys() instead of .get_all_primary_ids()
        """
        import warnings
        warnings.warn("Use bio_seq_database.keys() instead of "
                      "bio_seq_database.get_all_primary_ids()",
                      PendingDeprecationWarning)
        return self.keys()

    def __getitem__(self, key):
        return BioSeq.DBSeqRecord(self.adaptor, key)

    def __delitem__(self, key):
        """Remove an entry and all its annotation."""
        if key not in self:
            raise KeyError(key)
        #Assuming this will automatically cascade to the other tables...
        sql = "DELETE FROM bioentry " + \
              "WHERE biodatabase_id=%s AND bioentry_id=%s;"
        self.adaptor.execute(sql, (self.dbid,key))

    def __len__(self):
        """Number of records in this namespace (sub database)."""
        sql = "SELECT COUNT(bioentry_id) FROM bioentry " + \
              "WHERE biodatabase_id=%s;"
        return int(self.adaptor.execute_and_fetch_col0(sql, (self.dbid,))[0])

    def __contains__(self, value):
        """Check if a primary (internal) id is this namespace (sub database)."""
        sql = "SELECT COUNT(bioentry_id) FROM bioentry " + \
              "WHERE biodatabase_id=%s AND bioentry_id=%s;"
        #The bioentry_id field is an integer in the schema.
        #PostgreSQL will throw an error if we use a non integer in the query.
        try:
            bioentry_id = int(value)
        except ValueError:
            return False
        return bool(self.adaptor.execute_and_fetch_col0(sql,
                                                  (self.dbid, bioentry_id))[0])
    
    def __iter__(self):
        """Iterate over ids (which may not be meaningful outside this database)."""
        #TODO - Iterate over the cursor, much more efficient
        return iter(self.adaptor.list_bioentry_ids(self.dbid))        

    if hasattr(dict, "iteritems"):
        #Python 2, use iteritems etc    
        def keys(self):
            """List of ids which may not be meaningful outside this database."""
            return self.adaptor.list_bioentry_ids(self.dbid)

        def values(self):
            """List of DBSeqRecord objects in the namespace (sub database)."""
            return [self[key] for key in self.keys()]
    
        def items(self):
            """List of (id, DBSeqRecord) for the namespace (sub database)."""
            return [(key, self[key]) for key in self.keys()]
        
        def iterkeys(self):
            """Iterate over ids (which may not be meaningful outside this database)."""
            return iter(self)
    
        def itervalues(self):
            """Iterate over DBSeqRecord objects in the namespace (sub database)."""
            for key in self:
                yield self[key]
            
        def iteritems(self):
            """Iterate over (id, DBSeqRecord) for the namespace (sub database)."""
            for key in self:
                yield key, self[key]
    else:
        #Python 3, items etc are all iterators
        def keys(self):
            """Iterate over ids (which may not be meaningful outside this database)."""
            return iter(self)
            
        def values(self):
            """Iterate over DBSeqRecord objects in the namespace (sub database)."""
            for key in self:
                yield self[key]
    
        def items(self):
            """Iterate over (id, DBSeqRecord) for the namespace (sub database)."""
            for key in self:
                yield key, self[key]

    def lookup(self, **kwargs):
        if len(kwargs) != 1:
            raise TypeError("single key/value parameter expected")
        k, v = kwargs.items()[0]
        if k not in _allowed_lookups:
            raise TypeError("lookup() expects one of %s, not %r" % \
                            (repr(_allowed_lookups.keys())[1:-1], repr(k)))
        lookup_name = _allowed_lookups[k]
        lookup_func = getattr(self.adaptor, lookup_name)
        seqid = lookup_func(self.dbid, v)
        return BioSeq.DBSeqRecord(self.adaptor, seqid)
        
    def get_Seq_by_primary_id(self, seqid):
        """Get a DBSeqRecord by the primary (internal) id (OBSOLETE).
        
        Rather than db.get_Seq_by_primary_id(my_id) use db[my_id]
        
        The name of this method is misleading since it returns a DBSeqRecord
        rather than a DBSeq ojbect, and presumably was to mirror BioPerl.
        """
        import warnings
        warnings.warn("Use bio_seq_database[my_id] instead of "
                      "bio_seq_database.get_Seq_by_primary_id(my_id)",
                      PendingDeprecationWarning)
        return self[seqid]

    def load(self, record_iterator, fetch_NCBI_taxonomy=False):
        """Load a set of SeqRecords into the BioSQL database.

        record_iterator is either a list of SeqRecord objects, or an
        Iterator object that returns SeqRecord objects (such as the
        output from the Bio.SeqIO.parse() function), which will be
        used to populate the database.

        fetch_NCBI_taxonomy is boolean flag allowing or preventing
        connection to the taxonomic database on the NCBI server
        (via Bio.Entrez) to fetch a detailed taxonomy for each
        SeqRecord.

        Example:
        from Bio import SeqIO
        count = db.load(SeqIO.parse(open(filename), format))

        Returns the number of records loaded.
        """
        db_loader = Loader.DatabaseLoader(self.adaptor, self.dbid, \
                                          fetch_NCBI_taxonomy)
        num_records = 0
        global _POSTGRES_RULES_PRESENT
        for cur_record in record_iterator:
            num_records += 1
            #Hack to work arround BioSQL Bug 2839 - If using PostgreSQL and
            #the RULES are present check for a duplicate record before loading
            if _POSTGRES_RULES_PRESENT:
                #Recreate what the Loader's _load_bioentry_table will do:
                if cur_record.id.count(".") == 1:
                    accession, version = cur_record.id.split('.')
                    try:
                        version = int(version)
                    except ValueError:
                        accession = cur_record.id
                        version = 0
                else:
                    accession = cur_record.id
                    version = 0
                gi = cur_record.annotations.get("gi", None)
                sql = "SELECT bioentry_id FROM bioentry WHERE (identifier " + \
                      "= '%s' AND biodatabase_id = '%s') OR (accession = " + \
                      "'%s' AND version = '%s' AND biodatabase_id = '%s')"
                self.adaptor.execute(sql % (gi, self.dbid, accession, version, self.dbid))
                if self.adaptor.cursor.fetchone():
                    raise self.adaptor.conn.IntegrityError("Duplicate record " 
                                     "detected: record has not been inserted")
            #End of hack
            db_loader.load_seqrecord(cur_record)
        return num_records
