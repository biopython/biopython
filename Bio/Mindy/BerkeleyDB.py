import os
from bsddb3 import db
import Location
import BaseDB
import Bio

_open = open  # rename for internal use -- gets redefined below

class CreateBerkeleyDB(BaseDB.CreateDB):
    def __init__(self, dbname, primary_namespace, secondary_namespaces,
                 format = "sequence"):
        self.dbname = dbname        
        self.primary_namespace = primary_namespace
        self.secondary_namespaces = secondary_namespaces
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
        self.config["primary_namespace"] = primary_namespace
        self.config["secondary_namespaces"] = "\t".join(secondary_namespaces)

        self.fileids = {}
        self.filemap = db.DB(dbenv)
        self.filemap.open("bioindex", "identifiers",
                          db.DB_BTREE, db.DB_CREATE, 0660)

        self.primary_table = db.DB(dbenv)
        self.primary_table.open("bioindex", "key_%s" % (self.primary_namespace,),
                                db.DB_BTREE, db.DB_CREATE, 0660)

        self.lookup_tables = {}
        for namespace in secondary_namespaces:
            x = db.DB(dbenv)
            x.open("bioindex", "id_%s" % (namespace,),
                   db.DB_BTREE, db.DB_CREATE, 0)
            self.lookup_tables[namespace] = x

        self.dbenv = dbenv

    def add_record(self, filetag, startpos, length, table):
        key_list = table[self.primary_namespace]
        if len(key_list) != 1:
            raise TypeError(
                "Field %s has %d entries but must have only one "
                "(must be unique)" % (repr(self.primary_namespace),
                                      len(key_list)))
        key = key_list[0]
        if self.primary_table.has_key(key):
            raise TypeError("Field %r = %r already exists" %
                            (self.primary_namespace, key))
        self.primary_table[key] = "%s\t%s\t%s" % (filetag, startpos, length)

        for namespace in self.secondary_namespaces:
            lookup = self.lookup_tables[namespace]
            # Get the list of secondary identifiers for this identifier
            for val in table.get(namespace, ()):
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

class PrimaryNamespace:
    def __init__(self, db, namespace):
        self.db = db
        self.namespace = namespace
        assert namespace == db.primary_namespace

    def __getitem__(self, name):
        loc = self.db.primary_table[name]
        filetag, start, length = loc.split("\t")
        filename = self.db.fileids[filetag]
        data = [
            Location.Location(self.namespace,
                              name,
                              filename,
                              long(start),
                              long(length))
            ]
        return data

    def keys(self):
        return self.db.primary_table.keys()
    def values(self):
        return [self[key] for key in self.keys()]
    def items(self):
        return [ (key, self[key]) for key in self.keys()]
    

class SecondaryNamespace:
    def __init__(self, db, namespace):
        self.db = db
        self.namespace = namespace
        assert namespace in db.secondary_namespaces
        
    def __getitem__(self, name):
        table = self.db._load_namespace(self.namespace)
        text = table.get(name, None)
        if text is None:
            raise KeyError("Cannot find %r key %r" % (self.namespace, name))
        data = []
        for key in text.split("\t"):
            loc = self.db.primary_table[key]
            filetag, start, length = loc.split("\t")
            filename = self.db.fileids[filetag]
            data.append(Location.Location(self.namespace,
                                          name,
                                          filename,
                                          long(start),
                                          long(length)))
            
        return data
    
    def keys(self):
        table = self.db._load_namespace(self.namespace)
        return table.keys()
    def values(self):
        return [self[key] for key in self.keys()]
    def items(self):
        return [ (key, self[key]) for key in self.keys()]

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
        self.primary_namespace = self.config["primary_namespace"]
        self.secondary_namespaces = self.config["secondary_namespaces"].split("\t")

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
        self.primary_table.open("bioindex", "key_%s" % (self.primary_namespace,),
                                db.DB_BTREE, db.DB_CREATE, 0660)

        self.lookup_tables = {}
        self.dbenv = dbenv

    def _load_namespace(self, namespace):
        dbname = "id_%s" % namespace
        # Get the appropriate lookup table
        if not self.lookup_tables.has_key(dbname):
            # Nope, so load it up
            self.lookup_tables[dbname] = db.DB(self.dbenv)
            self.lookup_tables[dbname].open("bioindex", dbname,
                                            db.DB_BTREE, db.DB_CREATE, 0)
        return self.lookup_tables[dbname]

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

    def __getitem__(self, key):
        if key not in self.keys():
            raise KeyError(key)
        if key == self.primary_namespace:
            return PrimaryNamespace(self, key)
        else:
            return SecondaryNamespace(self, key)
