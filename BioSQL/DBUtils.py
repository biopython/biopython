# Copyright 2002 by Andrew Dalke.  All rights reserved.
# Revisions 2007-2009 copyright by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Note that BioSQL (including the database schema and scripts) is
# available and licensed separately.  Please consult www.biosql.org

_dbutils = {}

class Generic_dbutils:
    def __init__(self):
        pass

    def tname(self, table):
        if table != 'biosequence': return table
        else: return 'bioentry'

# Disabled: better safe than sorry
##    def next_id(self, cursor, table):
##        # XXX brain-dead! Hopefully, the database will enforce PK unicity..
##        table = self.tname(table)
##        sql = r"select 1+max(%s_id) from %s" % (table, table)
##        cursor.execute(sql)
##        rv = cursor.fetchone()
##        return rv[0]
        
    def last_id(self, cursor, table):
        # XXX: Unsafe without transactions isolation
        table = self.tname(table)
        sql = r"select max(%s_id) from %s" % (table, table)
        cursor.execute(sql)
        rv = cursor.fetchone()
        return rv[0]

    def autocommit(self, conn, y = 1):
        # Let's hope it was not really needed
        pass

class Mysql_dbutils(Generic_dbutils):
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

class Psycopg_dbutils(Generic_dbutils):
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

    def autocommit(self, conn, y = True):
        conn.autocommit(y)

_dbutils["psycopg"] = Psycopg_dbutils
 
class Psycopg2_dbutils(Psycopg_dbutils):
    def autocommit(self, conn, y = True):
        if y:
            conn.set_isolation_level(0)
        else:
            conn.set_isolation_level(1)

_dbutils["psycopg2"] = Psycopg2_dbutils

class Pgdb_dbutils(Generic_dbutils):
    """Add support for pgdb in the PyGreSQL database connectivity package.
    """
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

    def autocommit(self, conn, y = True):
        raise NotImplementedError("pgdb does not support this!")

_dbutils["pgdb"] = Pgdb_dbutils

def get_dbutils(module_name):
    try:
        return _dbutils[module_name]()
    except:
        return Generic_dbutils()
