import os
from bsddb3 import db
import Location
import BaseDB
import Bio

_open = open  # rename for internal use -- gets redefined below

class CreateBerkeleyDB(BaseDB.CreateDB):
    def __init__(self, dbname, unique, data_fields, format = "sequence"):
        self.dbname = dbname        
        self.unique = unique
        self.data_fields = data_fields
        self.dbenv = None
        self.format = Bio.formats.normalize(format)

        os.mkdir(dbname)
        outfile = _open(os.path.join(dbname, "config.dat"), "wb")
        outfile.write("index\tBerkeleyDB/1\n")
        outfile.close()
        
        dbenv = db.DBEnv(0)
        envflags = db.DB_THREAD | db.DB_INIT_MPOOL
        dbenv.open(dbname, envflags | db.DB_CREATE)

        self.config = db.DB(dbenv)
        self.config.open("bioindex", "config",
                         db.DB_BTREE, db.DB_CREATE, 0660)
        self.config["namespace"] = unique
                         

        self.fileids = {}
        self.filemap = db.DB(dbenv)
        self.filemap.open("bioindex", "identifiers",
                          db.DB_BTREE, db.DB_CREATE, 0660)

        self.primary_table = db.DB(dbenv)
        self.primary_table.open("bioindex", "key_%s" % (self.unique,),
                                db.DB_BTREE, db.DB_CREATE, 0660)

        self.lookup_tables = {}
        for field in data_fields:
            x = db.DB(dbenv)
            x.open("bioindex", "id_%s" % (field,),
                   db.DB_BTREE, db.DB_CREATE, 0)
            self.lookup_tables[field] = x

        self.dbenv = dbenv

    def add_record(self, filetag, startpos, length, table):
        key_list = table[self.unique]
        if len(key_list) != 1:
            raise TypeError(
                "Field %s has %d entries but must have only one "
                "(must be unique)" % (repr(unique), len(key_list)))
        key = key_list[0]
        if self.primary_table.has_key(key):
            raise TypeError("Field %r = %r already exists" %
                            (self.unique, key))
        self.primary_table[key] = "%s\t%s\t%s" % (filetag, startpos, length)

        for field in self.data_fields:
            lookup = self.lookup_tables[field]
            # Get the list of secondary identifiers for this identifier
            for val in table.get(field, ()):
                # Go from secondary identifier to list of primary identifiers
                if lookup.has_key(val):
                    lookup[val] = lookup[val] + "\t" + key
                else:
                    lookup[val] = key

    def close(self):
        self.primary_table.close()
        self.filemap.close()
        self.config.close()
        [x.close() for x in self.lookup_tables.values()]
        self.dbenv.close()

        self.dbenv = self.primary_table = self.fileids = self.filemap = \
                     self.lookup_tables = None
    def __del__(self):
        if self.dbenv is not None:
            self.close()


def open(dbname):
    return OpenBerkeleyDB(dbname)

class OpenBerkeleyDB(BaseDB.OpenDB):
    def __init__(self, dbname):
        BaseDB.OpenDB.__init__(self, dbname)
        
        self.dbenv = None
        dbenv = db.DBEnv()
        envflags = db.DB_THREAD | db.DB_INIT_MPOOL
        dbenv.open(dbname, envflags | db.DB_CREATE)
        self.config = db.DB(dbenv)
        self.config.open("bioindex", "config",
                         db.DB_BTREE, 0, 0660)
        self.unique = self.config["namespace"]

        filemap = db.DB(dbenv)
        filemap.open("bioindex", "identifiers",
                          db.DB_BTREE, db.DB_CREATE, 0660)
        fileids = {}
        for k, v in filemap.items():
            filename, size = v.split("\t")
            size = long(size)
            if os.path.getsize(filename) != size:
                raise TypeError(
                    "File %s has changed size from %d to %d bytes!" %
                    (size, os.path.getsize(filename)))
            fileids[k] = filename
        self.fileids = fileids
        filemap.close()

        self.primary_table = db.DB(dbenv)
        self.primary_table.open("bioindex", "key_%s" % (self.unique,),
                                db.DB_BTREE, db.DB_CREATE, 0660)

        self.lookup_tables = {}
        self.dbenv = dbenv


    def lookup_query(self, namespace, name):
        if namespace == self.unique:
            loc = self.primary_table[name]
            filetag, start, length = loc.split("\t")
            filename = self.fileids[filetag]
            data = [
                Location.Location(namespace,
                                  name,
                                  filename,
                                  long(start),
                                  long(length))
                ]
            return data

        # There can be one or more aliases
        dbname = "id_%s" % namespace
        # Does it already exist?
        if not self.lookup_tables.has_key(dbname):
            # Nope, so load it up
            self.lookup_tables[dbname] = db.DB(self.dbenv)
            self.lookup_tables[dbname].open("bioindex", dbname,
                                            db.DB_BTREE, db.DB_CREATE, 0)

        table = self.lookup_tables[dbname]
        text = table.get(name, None)
        if text is None:
            raise KeyError("Cannot find %r key %r" % (namespace, name))
        data = []
        for key in text.split("\t"):
            loc = self.primary_table[key]
            filetag, start, length = loc.split("\t")
            filename = self.fileids[filetag]
            data.append(Location.Location(namespace,
                                          name,
                                          filename,
                                          long(start),
                                          long(length)))
        return data
                        
            
    def close(self):
        self.primary_table.close()
        self.config.close()
        [x.close() for x in self.lookup_tables.values()]
        self.dbenv.close()
        self.dbenv = self.primary_table = self.fileids = \
                     self.lookup_tables = None

        
    def __del__(self):
        if self.dbenv is not None:
            self.close()
