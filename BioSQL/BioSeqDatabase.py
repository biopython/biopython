# Copyright 2002 by Andrew Dalke.  All rights reserved.
# Revisions 2007-2016 copyright by Peter Cock.  All rights reserved.
# Revisions 2009 copyright by Cymon J. Cox.  All rights reserved.
# Revisions 2013-2014 copyright by Tiago Antao.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# Note that BioSQL (including the database schema and scripts) is
# available and licensed separately.  Please consult www.biosql.org
"""Connect with a BioSQL database and load Biopython like objects from it.

This provides interfaces for loading biological objects from a relational
database, and is compatible with the BioSQL standards.
"""
import os

from . import BioSeq
from . import Loader
from . import DBUtils


_POSTGRES_RULES_PRESENT = False  # Hack for BioSQL Bug 2839


def open_database(driver="MySQLdb", **kwargs):
    """Load an existing BioSQL-style database.

    This function is the easiest way to retrieve a connection to a
    database, doing something like::

        from BioSQL import BioSeqDatabase
        server = BioSeqDatabase.open_database(user="root", db="minidb")

    Arguments:
     - driver - The name of the database driver to use for connecting. The
       driver should implement the python DB API. By default, the MySQLdb
       driver is used.
     - user -the username to connect to the database with.
     - password, passwd - the password to connect with
     - host - the hostname of the database
     - database or db - the name of the database

    """
    if driver == "psycopg":
        raise ValueError(
            "Using BioSQL with psycopg (version one) is no "
            "longer supported. Use psycopg2 instead."
        )

    if os.name == "java":
        from com.ziclix.python.sql import zxJDBC

        module = zxJDBC
        if driver in ["MySQLdb"]:
            jdbc_driver = "com.mysql.jdbc.Driver"
            url_pref = "jdbc:mysql://" + kwargs["host"] + "/"
        elif driver in ["psycopg2"]:
            jdbc_driver = "org.postgresql.Driver"
            url_pref = "jdbc:postgresql://" + kwargs["host"] + "/"

    else:
        module = __import__(driver, fromlist=["connect"])
    connect = module.connect

    # Different drivers use different keywords...
    kw = kwargs.copy()
    if driver in ["MySQLdb", "mysql.connector"] and os.name != "java":
        if "database" in kw:
            kw["db"] = kw["database"]
            del kw["database"]
        if "password" in kw:
            kw["passwd"] = kw["password"]
            del kw["password"]
        # kw["charset"] = "utf8"
        # kw["use_unicode"] = True
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
    if os.name == "java":
        if driver in ["MySQLdb"]:
            conn = connect(
                url_pref + kw.get("database", "mysql"),
                kw["user"],
                kw["password"],
                jdbc_driver,
            )
        elif driver in ["psycopg2"]:
            conn = connect(
                url_pref + kw.get("database", "postgresql") + "?stringtype=unspecified",
                kw["user"],
                kw["password"],
                jdbc_driver,
            )
    elif driver in ["sqlite3"]:
        conn = connect(kw["database"])
    else:
        conn = connect(**kw)

    if os.name == "java":
        server = DBServer(conn, module, driver)
    else:
        server = DBServer(conn, module)

    # Sets MySQL to allow double quotes, rather than only backticks
    if driver in ["MySQLdb", "mysql.connector"]:
        server.adaptor.execute("SET sql_mode='ANSI_QUOTES';")

    # TODO - Remove the following once BioSQL Bug 2839 is fixed.
    # Test for RULES in PostgreSQL schema, see also Bug 2833.
    if driver in ["psycopg2", "pgdb"]:
        sql = (
            "SELECT ev_class FROM pg_rewrite WHERE "
            "rulename='rule_bioentry_i1' OR "
            "rulename='rule_bioentry_i2';"
        )
        if server.adaptor.execute_and_fetchall(sql):
            import warnings
            from Bio import BiopythonWarning

            warnings.warn(
                "Your BioSQL PostgreSQL schema includes some rules "
                "currently required for bioperl-db but which may"
                "cause problems loading data using Biopython (see "
                "BioSQL's RedMine Bug 2839 aka GitHub Issue 4 "
                "https://github.com/biosql/biosql/issues/4). "
                "If you do not use BioPerl, please remove these "
                "rules. Biopython should cope with the rules "
                "present, but with a performance penalty when "
                "loading new records.",
                BiopythonWarning,
            )
            global _POSTGRES_RULES_PRESENT
            _POSTGRES_RULES_PRESENT = True

    elif driver == "sqlite3":
        # Tell SQLite that we want to use foreign keys
        # https://www.sqlite.org/foreignkeys.html#fk_enable
        server.adaptor.execute("PRAGMA foreign_keys = ON")

    return server


class DBServer:
    """Represents a BioSQL database containing namespaces (sub-databases).

    This acts like a Python dictionary, giving access to each namespace
    (defined by a row in the biodatabase table) as a BioSeqDatabase object.
    """

    def __init__(self, conn, module, module_name=None):
        """Create a DBServer object.

        Arguments:
         - conn - A database connection object
         - module - The module used to create the database connection
         - module_name - Optionally, the name of the module. Default: module.__name__

        Normally you would not want to create a DBServer object yourself.
        Instead use the open_database function, which returns an instance of DBServer.
        """
        self.module = module
        if module_name is None:
            module_name = module.__name__
        if module_name == "mysql.connector":
            wrap_cursor = True
        else:
            wrap_cursor = False
        # Get module specific Adaptor or the base (general) Adaptor
        Adapt = _interface_specific_adaptors.get(module_name, Adaptor)
        self.adaptor = Adapt(
            conn, DBUtils.get_dbutils(module_name), wrap_cursor=wrap_cursor
        )
        self.module_name = module_name

    def __repr__(self):
        """Return a short description of the class name and database connection."""
        return f"{self.__class__.__name__}({self.adaptor.conn!r})"

    def __getitem__(self, name):
        """Return a BioSeqDatabase object.

        Arguments:
            - name - The name of the BioSeqDatabase

        """
        return BioSeqDatabase(self.adaptor, name)

    def __len__(self):
        """Return number of namespaces (sub-databases) in this database."""
        sql = "SELECT COUNT(name) FROM biodatabase;"
        return int(self.adaptor.execute_and_fetch_col0(sql)[0])

    def __contains__(self, value):
        """Check if a namespace (sub-database) in this database."""
        sql = "SELECT COUNT(name) FROM biodatabase WHERE name=%s;"
        return bool(self.adaptor.execute_and_fetch_col0(sql, (value,))[0])

    def __iter__(self):
        """Iterate over namespaces (sub-databases) in the database."""
        # TODO - Iterate over the cursor, much more efficient
        return iter(self.adaptor.list_biodatabase_names())

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
        db_id = self.adaptor.fetch_dbid_by_dbname(name)
        remover = Loader.DatabaseRemover(self.adaptor, db_id)
        remover.remove()

    def new_database(self, db_name, authority=None, description=None):
        """Add a new database to the server and return it."""
        # make the database
        sql = (
            "INSERT INTO biodatabase (name, authority, description)"
            " VALUES (%s, %s, %s)"
        )
        self.adaptor.execute(sql, (db_name, authority, description))
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
        sql = ""
        with open(sql_file) as sql_handle:
            for line in sql_handle:
                if line.startswith("--"):  # don't include comment lines
                    pass
                elif line.startswith("#"):  # ditto for MySQL comments
                    pass
                elif line.strip():  # only include non-blank lines
                    sql += line.strip() + " "

        # two ways to load the SQL
        # 1. PostgreSQL can load it all at once and actually needs to
        # due to FUNCTION defines at the end of the SQL which mess up
        # the splitting by semicolons
        if self.module_name in ["psycopg2", "pgdb"]:
            self.adaptor.cursor.execute(sql)
        # 2. MySQL needs the database loading split up into single lines of
        # SQL executed one at a time
        elif self.module_name in ["mysql.connector", "MySQLdb", "sqlite3"]:
            sql_parts = sql.split(";")  # one line per sql command
            # don't use the last item, it's blank
            for sql_line in sql_parts[:-1]:
                self.adaptor.cursor.execute(sql_line)
        else:
            raise ValueError(f"Module {self.module_name} not supported by the loader.")

    def commit(self):
        """Commit the current transaction to the database."""
        return self.adaptor.commit()

    def rollback(self):
        """Roll-back the current transaction."""
        return self.adaptor.rollback()

    def close(self):
        """Close the connection. No further activity possible."""
        return self.adaptor.close()


class _CursorWrapper:
    """A wrapper for mysql.connector resolving bytestring representations."""

    def __init__(self, real_cursor):
        self.real_cursor = real_cursor

    def execute(self, operation, params=None, multi=False):
        """Execute a sql statement."""
        self.real_cursor.execute(operation, params, multi)

    def executemany(self, operation, params):
        """Execute many sql statements."""
        self.real_cursor.executemany(operation, params)

    def _convert_tuple(self, tuple_):
        """Decode any bytestrings present in the row (PRIVATE)."""
        tuple_list = list(tuple_)
        for i, elem in enumerate(tuple_list):
            if type(elem) is bytes:
                tuple_list[i] = elem.decode("utf-8")
        return tuple(tuple_list)

    def _convert_list(self, lst):
        ret_lst = []
        for tuple_ in lst:
            new_tuple = self._convert_tuple(tuple_)
            ret_lst.append(new_tuple)
        return ret_lst

    def fetchall(self):
        rv = self.real_cursor.fetchall()
        return self._convert_list(rv)

    def fetchone(self):
        tuple_ = self.real_cursor.fetchone()
        return self._convert_tuple(tuple_)


class Adaptor:
    """High level wrapper for a database connection and cursor.

    Most database calls in BioSQL are done indirectly though this adaptor
    class. This provides helper methods for fetching data and executing
    sql.
    """

    def __init__(self, conn, dbutils, wrap_cursor=False):
        """Create an Adaptor object.

        Arguments:
         - conn - A database connection
         - dbutils - A BioSQL.DBUtils object
         - wrap_cursor - Optional, whether to wrap the cursor object

        """
        self.conn = conn
        if wrap_cursor:
            self.cursor = _CursorWrapper(conn.cursor())
        else:
            self.cursor = conn.cursor()
        self.dbutils = dbutils

    def last_id(self, table):
        """Return the last row id for the selected table."""
        return self.dbutils.last_id(self.cursor, table)

    def autocommit(self, y=True):
        """Set the autocommit mode. True values enable; False value disable."""
        return self.dbutils.autocommit(self.conn, y)

    def commit(self):
        """Commit the current transaction."""
        return self.conn.commit()

    def rollback(self):
        """Roll-back the current transaction."""
        return self.conn.rollback()

    def close(self):
        """Close the connection. No further activity possible."""
        return self.conn.close()

    def fetch_dbid_by_dbname(self, dbname):
        """Return the internal id for the sub-database using its name."""
        self.execute(
            "select biodatabase_id from biodatabase where name = %s", (dbname,)
        )
        rv = self.cursor.fetchall()
        if not rv:
            raise KeyError(f"Cannot find biodatabase with name {dbname!r}")
        return rv[0][0]

    def fetch_seqid_by_display_id(self, dbid, name):
        """Return the internal id for a sequence using its display id.

        Arguments:
         - dbid - the internal id for the sub-database
         - name - the name of the sequence. Corresponds to the
           name column of the bioentry table of the SQL schema

        """
        sql = "select bioentry_id from bioentry where name = %s"
        fields = [name]
        if dbid:
            sql += " and biodatabase_id = %s"
            fields.append(dbid)
        self.execute(sql, fields)
        rv = self.cursor.fetchall()
        if not rv:
            raise IndexError(f"Cannot find display id {name!r}")
        if len(rv) > 1:
            raise IndexError(f"More than one entry with display id {name!r}")
        return rv[0][0]

    def fetch_seqid_by_accession(self, dbid, name):
        """Return the internal id for a sequence using its accession.

        Arguments:
         - dbid - the internal id for the sub-database
         - name - the accession of the sequence. Corresponds to the
           accession column of the bioentry table of the SQL schema

        """
        sql = "select bioentry_id from bioentry where accession = %s"
        fields = [name]
        if dbid:
            sql += " and biodatabase_id = %s"
            fields.append(dbid)
        self.execute(sql, fields)
        rv = self.cursor.fetchall()
        if not rv:
            raise IndexError(f"Cannot find accession {name!r}")
        if len(rv) > 1:
            raise IndexError(f"More than one entry with accession {name!r}")
        return rv[0][0]

    def fetch_seqids_by_accession(self, dbid, name):
        """Return a list internal ids using an accession.

        Arguments:
         - dbid - the internal id for the sub-database
         - name - the accession of the sequence. Corresponds to the
           accession column of the bioentry table of the SQL schema

        """
        sql = "select bioentry_id from bioentry where accession = %s"
        fields = [name]
        if dbid:
            sql += " and biodatabase_id = %s"
            fields.append(dbid)
        return self.execute_and_fetch_col0(sql, fields)

    def fetch_seqid_by_version(self, dbid, name):
        """Return the internal id for a sequence using its accession and version.

        Arguments:
         - dbid - the internal id for the sub-database
         - name - the accession of the sequence containing a version number.
           Must correspond to <accession>.<version>

        """
        acc_version = name.split(".")
        if len(acc_version) > 2:
            raise IndexError(f"Bad version {name!r}")
        acc = acc_version[0]
        if len(acc_version) == 2:
            version = acc_version[1]
        else:
            version = "0"
        sql = "SELECT bioentry_id FROM bioentry WHERE accession = %s AND version = %s"
        fields = [acc, version]
        if dbid:
            sql += " and biodatabase_id = %s"
            fields.append(dbid)
        self.execute(sql, fields)
        rv = self.cursor.fetchall()
        if not rv:
            raise IndexError(f"Cannot find version {name!r}")
        if len(rv) > 1:
            raise IndexError(f"More than one entry with version {name!r}")
        return rv[0][0]

    def fetch_seqid_by_identifier(self, dbid, identifier):
        """Return the internal id for a sequence using its identifier.

        Arguments:
         - dbid - the internal id for the sub-database
         - identifier - the identifier of the sequence. Corresponds to
           the identifier column of the bioentry table in the SQL schema.

        """
        # YB: was fetch_seqid_by_seqid
        sql = "SELECT bioentry_id FROM bioentry WHERE identifier = %s"
        fields = [identifier]
        if dbid:
            sql += " and biodatabase_id = %s"
            fields.append(dbid)
        self.execute(sql, fields)
        rv = self.cursor.fetchall()
        if not rv:
            raise IndexError(f"Cannot find display id {identifier!r}")
        return rv[0][0]

    def list_biodatabase_names(self):
        """Return a list of all of the sub-databases."""
        return self.execute_and_fetch_col0("SELECT name FROM biodatabase")

    def list_bioentry_ids(self, dbid):
        """Return a list of internal ids for all of the sequences in a sub-databae.

        Arguments:
         - dbid - The internal id for a sub-database

        """
        return self.execute_and_fetch_col0(
            "SELECT bioentry_id FROM bioentry WHERE biodatabase_id = %s", (dbid,)
        )

    def list_bioentry_display_ids(self, dbid):
        """Return a list of all sequence names in a sub-databae.

        Arguments:
         - dbid - The internal id for a sub-database

        """
        return self.execute_and_fetch_col0(
            "SELECT name FROM bioentry WHERE biodatabase_id = %s", (dbid,)
        )

    def list_any_ids(self, sql, args):
        """Return ids given a SQL statement to select for them.

        This assumes that the given SQL does a SELECT statement that
        returns a list of items. This parses them out of the 2D list
        they come as and just returns them in a list.
        """
        return self.execute_and_fetch_col0(sql, args)

    def execute_one(self, sql, args=None):
        """Execute sql that returns 1 record, and return the record."""
        self.execute(sql, args or ())
        rv = self.cursor.fetchall()
        if len(rv) != 1:
            raise ValueError(f"Expected 1 response, got {len(rv)}.")
        return rv[0]

    def execute(self, sql, args=None):
        """Just execute an sql command."""
        if os.name == "java":
            sql = sql.replace("%s", "?")
        self.dbutils.execute(self.cursor, sql, args)

    def executemany(self, sql, args):
        """Execute many sql commands."""
        if os.name == "java":
            sql = sql.replace("%s", "?")
        self.dbutils.executemany(self.cursor, sql, args)

    def get_subseq_as_string(self, seqid, start, end):
        """Return a substring of a sequence.

        Arguments:
         - seqid - The internal id for the sequence
         - start - The start position of the sequence; 0-indexed
         - end - The end position of the sequence

        """
        length = end - start
        # XXX Check this on MySQL and PostgreSQL. substr should be general,
        # does it need dbutils?
        # return self.execute_one(
        #    """select SUBSTRING(seq FROM %s FOR %s)
        #             from biosequence where bioentry_id = %s""",
        #    (start+1, length, seqid))[0]
        return self.execute_one(
            "SELECT SUBSTR(seq, %s, %s) FROM biosequence WHERE bioentry_id = %s",
            (start + 1, length, seqid),
        )[0]

    def execute_and_fetch_col0(self, sql, args=None):
        """Return a list of values from the first column in the row."""
        self.execute(sql, args or ())
        return [field[0] for field in self.cursor.fetchall()]

    def execute_and_fetchall(self, sql, args=None):
        """Return a list of tuples of all rows."""
        self.execute(sql, args or ())
        return self.cursor.fetchall()


class MysqlConnectorAdaptor(Adaptor):
    """A BioSQL Adaptor class with fixes for the MySQL interface.

    BioSQL was failing due to returns of bytearray objects from
    the mysql-connector-python database connector. This adaptor
    class scrubs returns of bytearrays and of byte strings converting
    them to string objects instead. This adaptor class was made in
    response to backwards incompatible changes added to
    mysql-connector-python in release 2.0.0 of the package.
    """

    @staticmethod
    def _bytearray_to_str(s):
        """If s is bytes or bytearray, convert to a string (PRIVATE)."""
        if isinstance(s, (bytes, bytearray)):
            return s.decode()
        return s

    def execute_one(self, sql, args=None):
        """Execute sql that returns 1 record, and return the record."""
        out = super().execute_one(sql, args)
        return tuple(self._bytearray_to_str(v) for v in out)

    def execute_and_fetch_col0(self, sql, args=None):
        """Return a list of values from the first column in the row."""
        out = super().execute_and_fetch_col0(sql, args)
        return [self._bytearray_to_str(column) for column in out]

    def execute_and_fetchall(self, sql, args=None):
        """Return a list of tuples of all rows."""
        out = super().execute_and_fetchall(sql, args)
        return [tuple(self._bytearray_to_str(v) for v in o) for o in out]


_interface_specific_adaptors = {
    # If SQL interfaces require a specific adaptor, use this to map the adaptor
    "mysql.connector": MysqlConnectorAdaptor,
    "MySQLdb": MysqlConnectorAdaptor,
}

_allowed_lookups = {
    # Lookup name / function name to get id, function to list all ids
    "primary_id": "fetch_seqid_by_identifier",
    "gi": "fetch_seqid_by_identifier",
    "display_id": "fetch_seqid_by_display_id",
    "name": "fetch_seqid_by_display_id",
    "accession": "fetch_seqid_by_accession",
    "version": "fetch_seqid_by_version",
}


class BioSeqDatabase:
    """Represents a namespace (sub-database) within the BioSQL database.

    i.e. One row in the biodatabase table, and all all rows in the bioentry
    table associated with it.
    """

    def __init__(self, adaptor, name):
        """Create a BioDatabase object.

        Arguments:
         - adaptor - A BioSQL.Adaptor object
         - name - The name of the sub-database (namespace)

        """
        self.adaptor = adaptor
        self.name = name
        self.dbid = self.adaptor.fetch_dbid_by_dbname(name)

    def __repr__(self):
        """Return a short summary of the BioSeqDatabase."""
        return f"BioSeqDatabase({self.adaptor!r}, {self.name!r})"

    def get_Seq_by_id(self, name):
        """Get a DBSeqRecord object by its name.

        Example: seq_rec = db.get_Seq_by_id('ROA1_HUMAN')

        The name of this method is misleading since it returns a DBSeqRecord
        rather than a Seq object, and presumably was to mirror BioPerl.
        """
        seqid = self.adaptor.fetch_seqid_by_display_id(self.dbid, name)
        return BioSeq.DBSeqRecord(self.adaptor, seqid)

    def get_Seq_by_acc(self, name):
        """Get a DBSeqRecord object by accession number.

        Example: seq_rec = db.get_Seq_by_acc('X77802')

        The name of this method is misleading since it returns a DBSeqRecord
        rather than a Seq object, and presumably was to mirror BioPerl.
        """
        seqid = self.adaptor.fetch_seqid_by_accession(self.dbid, name)
        return BioSeq.DBSeqRecord(self.adaptor, seqid)

    def get_Seq_by_ver(self, name):
        """Get a DBSeqRecord object by version number.

        Example: seq_rec = db.get_Seq_by_ver('X77802.1')

        The name of this method is misleading since it returns a DBSeqRecord
        rather than a Seq object, and presumably was to mirror BioPerl.
        """
        seqid = self.adaptor.fetch_seqid_by_version(self.dbid, name)
        return BioSeq.DBSeqRecord(self.adaptor, seqid)

    def get_Seqs_by_acc(self, name):
        """Get a list of DBSeqRecord objects by accession number.

        Example: seq_recs = db.get_Seq_by_acc('X77802')

        The name of this method is misleading since it returns a list of
        DBSeqRecord objects rather than a list of Seq objects, and presumably
        was to mirror BioPerl.
        """
        seqids = self.adaptor.fetch_seqids_by_accession(self.dbid, name)
        return [BioSeq.DBSeqRecord(self.adaptor, seqid) for seqid in seqids]

    def __getitem__(self, key):
        """Return a DBSeqRecord for one of the sequences in the sub-database.

        Arguments:
         - key - The internal id for the sequence

        """
        record = BioSeq.DBSeqRecord(self.adaptor, key)
        if record._biodatabase_id != self.dbid:
            raise KeyError(f"Entry {key!r} does exist, but not in current name space")
        return record

    def __delitem__(self, key):
        """Remove an entry and all its annotation."""
        if key not in self:
            raise KeyError(
                f"Entry {key!r} cannot be deleted. It was not found or is invalid"
            )
        # Assuming this will automatically cascade to the other tables...
        sql = "DELETE FROM bioentry WHERE biodatabase_id=%s AND bioentry_id=%s;"
        self.adaptor.execute(sql, (self.dbid, key))

    def __len__(self):
        """Return number of records in this namespace (sub database)."""
        sql = "SELECT COUNT(bioentry_id) FROM bioentry WHERE biodatabase_id=%s;"
        return int(self.adaptor.execute_and_fetch_col0(sql, (self.dbid,))[0])

    def __contains__(self, value):
        """Check if a primary (internal) id is this namespace (sub database)."""
        sql = (
            "SELECT COUNT(bioentry_id) FROM bioentry "
            "WHERE biodatabase_id=%s AND bioentry_id=%s;"
        )
        # The bioentry_id field is an integer in the schema.
        # PostgreSQL will throw an error if we use a non integer in the query.
        try:
            bioentry_id = int(value)
        except ValueError:
            return False
        return bool(
            self.adaptor.execute_and_fetch_col0(sql, (self.dbid, bioentry_id))[0]
        )

    def __iter__(self):
        """Iterate over ids (which may not be meaningful outside this database)."""
        # TODO - Iterate over the cursor, much more efficient
        return iter(self.adaptor.list_bioentry_ids(self.dbid))

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
        """Return a DBSeqRecord using an acceptable identifier.

        Arguments:
         - kwargs - A single key-value pair where the key is one
           of primary_id, gi, display_id, name, accession, version

        """
        if len(kwargs) != 1:
            raise TypeError("single key/value parameter expected")
        k, v = list(kwargs.items())[0]
        if k not in _allowed_lookups:
            raise TypeError(
                f"lookup() expects one of {list(_allowed_lookups.keys())!r}, not {k!r}"
            )
        lookup_name = _allowed_lookups[k]
        lookup_func = getattr(self.adaptor, lookup_name)
        seqid = lookup_func(self.dbid, v)
        return BioSeq.DBSeqRecord(self.adaptor, seqid)

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

        Example::

            from Bio import SeqIO
            count = db.load(SeqIO.parse(open(filename), format))

        Returns the number of records loaded.
        """
        db_loader = Loader.DatabaseLoader(self.adaptor, self.dbid, fetch_NCBI_taxonomy)
        num_records = 0
        global _POSTGRES_RULES_PRESENT
        for cur_record in record_iterator:
            num_records += 1
            # Hack to work around BioSQL Bug 2839 - If using PostgreSQL and
            # the RULES are present check for a duplicate record before loading
            if _POSTGRES_RULES_PRESENT:
                # Recreate what the Loader's _load_bioentry_table will do:
                if cur_record.id.count(".") == 1:
                    accession, version = cur_record.id.split(".")
                    try:
                        version = int(version)
                    except ValueError:
                        accession = cur_record.id
                        version = 0
                else:
                    accession = cur_record.id
                    version = 0
                gi = cur_record.annotations.get("gi")
                sql = (
                    "SELECT bioentry_id FROM bioentry "
                    "WHERE (identifier = '%s' AND biodatabase_id = '%s') "
                    "OR (accession = '%s' AND version = '%s' AND biodatabase_id = '%s')"
                )
                self.adaptor.execute(
                    sql % (gi, self.dbid, accession, version, self.dbid)
                )
                if self.adaptor.cursor.fetchone():
                    raise self.adaptor.conn.IntegrityError(
                        "Duplicate record detected: record has not been inserted"
                    )
            # End of hack
            db_loader.load_seqrecord(cur_record)
        return num_records
