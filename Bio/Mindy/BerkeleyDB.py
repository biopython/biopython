import os
from bsddb3 import db
import Location
import BaseDB
import Bio

_open = open  # rename for internal use -- gets redefined below

INDEX_TYPE = "BerkeleyDB/1"

def create(dbname, primary_namespace, secondary_namespaces,
           formatname = "unknown"):
    os.mkdir(dbname)
    config_filename = os.path.join(dbname, "config.dat")
    BaseDB.write_config(config_filename = config_filename,
                        index_type = INDEX_TYPE,
                        primary_namespace = primary_namespace,
                        secondary_namespaces = secondary_namespaces,
                        fileid_info = {},
                        formatname = formatname
                        )

    dbenv = db.DBEnv(0)
    envflags = db.DB_THREAD | db.DB_INIT_MPOOL
    dbenv.open(dbname, envflags | db.DB_CREATE)

    primary_table = db.DB(dbenv)
    primary_table.open("key_%s" % (primary_namespace,), None,
                       db.DB_BTREE, db.DB_CREATE, 0660)

    secondary_tables = {}
    for namespace in secondary_namespaces:
        x = db.DB(dbenv)
        x.open("id_%s" % (namespace,), None, db.DB_BTREE, db.DB_CREATE, 0)
        secondary_tables[namespace] = x

    for x in secondary_tables.values():
        x.close()
    primary_table.close()
    dbenv.close()

    return open(dbname, "rw")
    

class PrimaryNamespace(BaseDB.DictLookup):
    def __init__(self, db, namespace):
        self.db = db
        self.namespace = namespace
        assert namespace == db.primary_namespace

    def __getitem__(self, name):
        loc = self.db.primary_table[name]
        filetag, startpos, length = loc.split("\t")
        filename = self.db.fileid_info[filetag][0]
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
            filename = self.db.fileid_info[filetag][0]
            data.append(Location.Location(self.namespace,
                                          name,
                                          filename,
                                          long(start),
                                          long(length)))
            
        return data
    
    def keys(self):
        table = self.db._load_namespace(self.namespace)
        return table.keys()

class BerkeleyDB(BaseDB.OpenDB, BaseDB.WriteDB):
    def __init__(self, dbname, mode = "r"):
        if mode not in ("r", "rw"):
            raise TypeError("Unknown mode: %r" % (mode,))
        self.__need_flush = 0
        BaseDB.OpenDB.__init__(self, dbname, INDEX_TYPE)
        
        self.dbenv = None
        dbenv = db.DBEnv()
        envflags = db.DB_THREAD | db.DB_INIT_MPOOL
        dbenv.open(dbname, envflags)
        if mode == "r":
            self._dbopen_flags = db.DB_RDONLY
        else:
            self._dbopen_flags = 0

        self.primary_table = db.DB(dbenv)
        self.primary_table.open("key_%s" % (self.primary_namespace,),
                                None,
                                db.DB_BTREE, self._dbopen_flags, 0660)

        self.secondary_tables = {}
        self.dbenv = dbenv

    def _load_namespace(self, namespace):
        dbname = "id_%s" % namespace
        # Get the appropriate lookup table
        if not self.secondary_tables.has_key(namespace):
            # Nope, so load it up
            self.secondary_tables[namespace] = db.DB(self.dbenv)
            self.secondary_tables[namespace].open(dbname, None,
                                               db.DB_BTREE,
                                               self._dbopen_flags, 0)
        return self.secondary_tables[namespace]


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
        self.__need_flush = 1

    def flush(self):
        if not self.__need_flush:
            return
        config_filename = os.path.join(self.dbname, "config.dat")
        BaseDB.write_config(config_filename = config_filename,
                            index_type = INDEX_TYPE,
                            primary_namespace = self.primary_namespace,
                            secondary_namespaces =
                                 self.secondary_tables.keys(),
                            fileid_info = self.fileid_info,
                            formatname = self.formatname,
                      )
        self.__need_flush = 0

    def close(self):
        self.flush()
        self.primary_table.close()
        [x.close() for x in self.secondary_tables.values()]
        self.dbenv.close()
        self.dbenv = self.primary_table = self.fileid_info = \
                     self.secondary_tables = self.fileid_info = None
        
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

open = BerkeleyDB

