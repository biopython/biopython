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
        sql = r"select max(%s_id) from %s" % table
        cursor.execute(sql)
        rv = cursor.fetchone()
        return rv[0]

    def autocommit(self, conn, y = 1):
        # Let's hope it was not really needed
        pass

class Mysql_dbutils(Generic_dbutils):
    def last_id(self, cursor, table):
        return cursor.insert_id()
_dbutils["MySQLdb"] = Mysql_dbutils

class Pg_dbutils(Generic_dbutils):
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
_dbutils["psycopg"] = Pg_dbutils

def get_dbutils(module_name):
    try:
        return _dbutils[module_name]()
    except:
        return Generic_dbutils()
