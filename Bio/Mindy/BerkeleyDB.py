import os
from bsddb3 import db
import Location
import BaseDB
import Bio

_open = open  # rename for internal use -- gets redefined below

def create(dbname, primary_namespace, secondary_namespaces):
    os.mkdir(dbname)
    outfile = _open(os.path.join(dbname, "config.dat"), "wb")
    outfile.write("index\tBerkeleyDB/1\n")
    outfile.close()

    dbenv = db.DBEnv(0)
    envflags = db.DB_THREAD | db.DB_INIT_MPOOL
    dbenv.open(dbname, envflags | db.DB_CREATE)

    config = db.DB(dbenv)
    config.open("bioindex", "config",
                     db.DB_BTREE, db.DB_CREATE, 0660)
    config["primary_namespace"] = primary_namespace
    config["secondary_namespaces"] = "\t".join(secondary_namespaces)

    fileids = {}
    filemap = db.DB(dbenv)
    filemap.open("bioindex", "identifiers",
                      db.DB_BTREE, db.DB_CREATE, 0660)

    primary_table = db.DB(dbenv)
    primary_table.open("bioindex", "key_%s" % (primary_namespace,),
                            db.DB_BTREE, db.DB_CREATE, 0660)

    lookup_tables = {}
    for namespace in secondary_namespaces:
        x = db.DB(dbenv)
        x.open("bioindex", "id_%s" % (namespace,),
               db.DB_BTREE, db.DB_CREATE, 0)
        lookup_tables[namespace] = x

    for x in lookup_tables.values():
        x.close()
    primary_table.close()
    filemap.close()
    config.close()
    dbenv.close()

    return open(dbname, "rw")
    

def open(dbname, mode = "r"):
    # XXX should handle modes somehow
    if mode not in ("r", "rw"):
        raise TypeError("Unknown mode: %r" % (mode,))
    return BerkeleyDB(dbname)

class PrimaryNamespace(BaseDB.DictLookup):
    def __init__(self, db, namespace):
        self.db = db
        self.namespace = namespace
        assert namespace == db.primary_namespace

    def __getitem__(self, name):
        loc = self.db.primary_table[name]
        filetag, startpos, length = loc.split("\t")
        filename = self.db.fileids[filetag]
        return [
            Location.Location(self.namespace,
                              name,
                              filename,
                              long(startpos),
                              long(length))
            ]

    def keys(self):
        return self.db.primary_table.keys()

class SecondaryNamespace(BaseDB.DictLookup):
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

class BerkeleyDB(BaseDB.OpenDB, BaseDB.CreateDB):
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

        self.filemap = db.DB(dbenv)
        self.filemap.open("bioindex", "identifiers",
                          db.DB_BTREE, db.DB_CREATE, 0660)
        fileids = {}
        for k, v in self.filemap.items():
            filename, size = v.split("\t")
            size = long(size)
            if os.path.getsize(filename) != size:
                raise TypeError(
                    "File %s has changed size from %d to %d bytes!" %
                    (size, os.path.getsize(filename)))
            fileids[k] = filename
        self.fileids = fileids

        self.primary_table = db.DB(dbenv)
        self.primary_table.open("bioindex",
                                "key_%s" % (self.primary_namespace,),
                                db.DB_BTREE, db.DB_CREATE, 0660)

        self.lookup_tables = {}
        self.dbenv = dbenv

    def _load_namespace(self, namespace):
        dbname = "id_%s" % namespace
        # Get the appropriate lookup table
        if not self.lookup_tables.has_key(namespace):
            # Nope, so load it up
            self.lookup_tables[namespace] = db.DB(self.dbenv)
            self.lookup_tables[namespace].open("bioindex", dbname,
                                               db.DB_BTREE, db.DB_CREATE, 0)
        return self.lookup_tables[namespace]


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
        self.primary_table[key] = "%s\t%s\t%s" % (filetag,
                                                  BaseDB._int_str(startpos),
                                                  BaseDB._int_str(length))

        for namespace in self.secondary_namespaces:
            lookup = self._load_namespace(namespace)
            # Get the list of secondary identifiers for this identifier
            for val in table.get(namespace, ()):
                # Go from secondary identifier to list of primary identifiers
                if lookup.has_key(val):
                    lookup[val] = lookup[val] + "\t" + key
                else:
                    lookup[val] = key

    def flush(self):
        pass

    def close(self):
        self.primary_table.close()
        self.config.close()
        [x.close() for x in self.lookup_tables.values()]
        self.filemap.close()
        self.dbenv.close()
        self.dbenv = self.primary_table = self.fileids = \
                     self.lookup_tables = self.filemap = None
        
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
