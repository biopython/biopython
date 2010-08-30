# Copyright 2002 by Andrew Dalke.  All rights reserved.
# Revisions 2007-2010 copyright by Peter Cock.  All rights reserved.
# Revisions 2009 copyright by Brad Chapman.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Note that BioSQL (including the database schema and scripts) is
# available and licensed separately.  Please consult www.biosql.org

_dbutils = {}

class Generic_dbutils:
    """Default database utilities."""
    def __init__(self):
        pass

    def tname(self, table):
        if table != 'biosequence': return table
        else: return 'bioentry'

    def last_id(self, cursor, table):
        # XXX: Unsafe without transactions isolation
        table = self.tname(table)
        sql = r"select max(%s_id) from %s" % (table, table)
        cursor.execute(sql)
        rv = cursor.fetchone()
        return rv[0]
    
    def execute(self, cursor, sql, args=None):
        """Just execute an sql command.
        """
        cursor.execute(sql, args or ())

    def autocommit(self, conn, y = 1):
        # Let's hope it was not really needed
        pass


class Sqlite_dbutils(Generic_dbutils):
    """Custom database utilities for SQLite."""
    def execute(self, cursor, sql, args=None):
        """Execute SQL command, replacing %s with ? for variable substitution in sqlite3.
        """
        cursor.execute(sql.replace("%s", "?"), args or ())

_dbutils["sqlite3"] = Sqlite_dbutils


class Mysql_dbutils(Generic_dbutils):
    """Custom database utilities for MySQL."""
    def last_id(self, cursor, table):
        try:
            #This worked on older versions of MySQL
            return cursor.insert_id()
        except AttributeError:
            #See bug 2390
            #Google suggests this is the new way,
            #same fix also suggested by Eric Gibert:
            return cursor.lastrowid
        
_dbutils["MySQLdb"] = Mysql_dbutils


class _PostgreSQL_dbutils(Generic_dbutils):
    """Base class for any PostgreSQL adaptor."""
    def next_id(self, cursor, table):
        table = self.tname(table)
        sql = r"select nextval('%s_pk_seq')" % table
        cursor.execute(sql)
        rv = cursor.fetchone()
        return rv[0]
        
    def last_id(self, cursor, table):
        table = self.tname(table)
        sql = r"select currval('%s_pk_seq')" % table
        cursor.execute(sql)
        rv = cursor.fetchone()
        return rv[0]

class Psycopg2_dbutils(_PostgreSQL_dbutils):
    """Custom database utilities for Psycopg2 (PostgreSQL)."""
    def autocommit(self, conn, y = True):
        if y:
            conn.set_isolation_level(0)
        else:
            conn.set_isolation_level(1)

_dbutils["psycopg2"] = Psycopg2_dbutils


class Pgdb_dbutils(_PostgreSQL_dbutils):
    """Custom database utilities for Pgdb (aka PyGreSQL, for PostgreSQL)."""
    def autocommit(self, conn, y = True):
        raise NotImplementedError("pgdb does not support this!")

_dbutils["pgdb"] = Pgdb_dbutils


def get_dbutils(module_name):
    try:
        return _dbutils[module_name]()
    except KeyError:
        return Generic_dbutils()
