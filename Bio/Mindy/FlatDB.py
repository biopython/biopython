
import os, bisect
import BaseDB, Location
import Bio

_open = open
INDEX_TYPE = "flat/1"

def _parse_primary_table_entry(s):
    name, filetag, startpos, length = s.rstrip().split("\t")
    return name, filetag, long(startpos), long(length)

def _read_primary_table(filename):
    infile = _open(filename, "rb")
    size = int(infile.read(4))
    table = {}
    while 1:
        s = infile.read(size)
        if not s:
            break
        assert len(s) == size, (repr(s), size)
        name, filetag, startpos, length = _parse_primary_table_entry(s)
        table[name] = filetag, startpos, length
    return table

def _write_primary_table(filename, primary_table):
    info = primary_table.items()
    info.sort()
    n = 1
    for k, v in info:
        # Find the longest width
        s = "%s\t%s" % (k, v)
        if len(s) > n:
            n = len(s)
            if n > 9999:
                raise AssertionError(
                    "Primary index record too large for format spec! " +
                    " %s bytes in %r" % (n, s))
    outfile = _open(filename, "wb")
    outfile.write("%04d" % n)
    for k, v in info:
        s = "%s\t%s" % (k, v)
        outfile.write(s.ljust(n))
    outfile.close()

def _parse_secondary_table_entry(s):
    return s.rstrip().split("\t")

def _read_secondary_table(filename):
    infile = _open(filename, "rb")
    size = int(infile.read(4))
    table = {}
    while 1:
        s = infile.read(size)
        if not s:
            break
        assert len(s) == size, (repr(s), size)
        alias, name = _parse_secondary_table_entry(s)
        table.setdefault(alias, []).append(name)
    infile.close()
    return table

def _write_secondary_table(filename, table):
    items = table.items()
    items.sort()
    # Find the largest field
    n = 0
    for k, v in items:
        for x in v:
            s = "%s\t%s" % (k, x)
            if len(s) > n:
                n = len(s)
                if n > 9999:
                    raise AssertionError(
               "Secondary index record too large for format spec! " +
               " %s bytes in %r" % (n, s))
    # And write the output
    outfile = _open(filename, "wb")
    outfile.write("%04d" % n)
    for k, v in items:
        for x in v:
            s = "%s\t%s" % (k, x)
            outfile.write(s.ljust(n))
    outfile.close()

class BaseFlatDB(BaseDB.OpenDB):
    def __init__(self, dbname):
        BaseDB.OpenDB.__init__(self, dbname, INDEX_TYPE)
        self.key_filename = os.path.join(dbname,
                                         "key_%s.key" % self.primary_namespace)
        
class PrimaryTable(BaseDB.DictLookup):
    def __init__(self, db, namespace, table):
        self.db = db
        self.namespace = namespace
        self.table = table
    def __getitem__(self, name):
        fileid, startpos, length = self.table[name]
        return [Location.Location(self.namespace,
                                  name,
                                  self.db.fileid_info[fileid][0],
                                  startpos,
                                  length)
                ]
    def keys(self):
        return self.table.keys()

class SecondaryTable(BaseDB.DictLookup):
    def __init__(self, db, namespace, table):
        self.db = db
        self.namespace = namespace
        self.table = table
    def __getitem__(self, name):
        data = []
        for entry in self.table[name]:
            fileid, startpos, length = self.db.primary_table[entry]
            data.append( Location.Location(self.namespace,
                                           name,
                                           self.db.fileid_info[fileid][0],
                                           startpos,
                                           length) )
        return data
    def keys(self):
        return self.table.keys()

class MemoryFlatDB(BaseDB.WriteDB, BaseFlatDB):
    def __init__(self, dbname):
        self.__in_constructor = 1
        self._need_flush = 0
        BaseFlatDB.__init__(self, dbname)

        primary_filename = os.path.join(self.dbname,
                             "key_%s.key" % (self.primary_namespace,) )
        self.primary_table = _read_primary_table(primary_filename)

        self.secondary_tables = {}
        for namespace in self.secondary_namespaces:
            filename = os.path.join(self.dbname, "id_%s.index" % namespace)
            self.secondary_tables[namespace] = _read_secondary_table(filename)

        self.__in_constructor = 0

    def add_record(self, filetag, startpos, length, table):
        key_list = table[self.primary_namespace]
        if len(key_list) != 1:
            raise TypeError(
                "Field %s has %d entries but must have only one "
                "(must be unique)" % (repr(unique), len(key_list)))
        key = key_list[0]
        if self.primary_table.has_key(key):
            raise TypeError("Field %r = %r already exists; must be unique" %
                            (self.primary_namespace, key))
        self.primary_table[key] = "%s\t%s\t%s" % (filetag,
                                                  BaseDB._int_str(startpos),
                                                  BaseDB._int_str(length))

        for namespace in self.secondary_namespaces:
            lookup = self.secondary_tables[namespace]
            # Get the list of secondary identifiers for this identifier
            for val in table.get(namespace, ()):
                # Go from secondary identifier to list of primary identifiers
                lookup.setdefault(val, []).append(key)
        self._need_flush = 1

    def flush(self):
        if not self._need_flush:
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

        primary_filename = os.path.join(self.dbname,
                           "key_%s.key" % (self.primary_namespace,) )
        _write_primary_table(filename = primary_filename,
                             primary_table = self.primary_table)


        # Write the secondary identifier information
        for namespace, table in self.secondary_tables.items():
            filename = os.path.join(self.dbname, "id_%s.index" % namespace)
            _write_secondary_table(filename = filename,
                                   table = table)

        self._need_flush = 0

    def close(self):
        self.flush()
        self.primary_table = self.fileid_info = self.filename_map = \
                             self.secondary_tables = None

    def __del__(self):
        if not self.__in_constructor:
            self.close()


    def __getitem__(self, namespace):
        """return the database table lookup for the given namespace"""
        if namespace == self.primary_namespace:
            return PrimaryTable(self, namespace, self.primary_table)
        return SecondaryTable(self, namespace, self.secondary_tables[namespace])


class BisectFile:
    def __init__(self, infile, size):
        self.infile = infile
        self.size = size
        infile.seek(0)
        self.record_size = int(infile.read(4))
        assert (size - 4) % self.record_size == 0, "record size is wrong"
    def __len__(self):
        if self.record_size == 0:
            return 0
        return int((self.size - 4) / self.record_size)
    def __getitem__(self, i):
        self.infile.seek(i * self.record_size + 4)
        return self.infile.read(self.record_size).split("\t")[0]
    def get_entry(self, i):
        self.infile.seek(i * self.record_size + 4)
        return self.infile.read(self.record_size)

def _find_entry(filename, wantword):
    size = os.path.getsize(filename)
    infile = _open(filename, "rb")

    bf = BisectFile(infile, size)
    left = bisect.bisect_left(bf, wantword)
    line = bf.get_entry(left)
    if not line.startswith(wantword):
        return None
    return line

def _find_range(filename, wantword):
    size = os.path.getsize(filename)
    infile = _open(filename, "rb")

    bf = BisectFile(infile, size)
    left = bisect.bisect_left(bf, wantword)
    line = bf.get_entry(left)
    if not line.startswith(wantword):
        return None
    
    right = bisect.bisect_right(bf, wantword)
    data = []
    for i in range(left, right):
        x = bf.get_entry(i)
        data.append(x)
    return data


def _lookup_location(key_filename, word):
    line = _find_entry(key_filename, word)
    if line is None:
        return None
    return _parse_primary_table_entry(line)[1:]

def _lookup_alias(id_filename, word):
    lines = _find_range(id_filename, word)
    if not lines:
        return None
    primary_keys = []
    for line in lines:
        alias, primary_key = _parse_secondary_table_entry(line)
        assert alias == word, (alias, word)
        primary_keys.append(primary_key)
    return primary_keys

def create(dbname, primary_namespace, secondary_namespaces,
           formatname = "unknown"):
    os.mkdir(dbname)
    config_filename = os.path.join(dbname, "config.dat")
    BaseDB.write_config(config_filename = config_filename,
                        index_type = INDEX_TYPE,
                        primary_namespace = primary_namespace,
                        secondary_namespaces = secondary_namespaces,
                        fileid_info = {},
                        formatname = formatname,
                        )

    primary_filename = os.path.join(dbname,
                       "key_%s.key" % (primary_namespace,) )
    _write_primary_table(filename = primary_filename,
                         primary_table = {})


    # Write the secondary identifier information
    for namespace in secondary_namespaces:
        filename = os.path.join(dbname, "id_%s.index" % namespace)
        _write_secondary_table(filename = filename,
                               table = {})
    return open(dbname, "rw")
    

def open(dbname, mode = "r"):
    if mode == "r":
        return DiskFlatDB(dbname)
    elif mode == "rw":
        return MemoryFlatDB(dbname)
    elif mode == "a":
        raise TypeError("Must call FlatDB.create to create the database")
    else:
        raise TypeError("Unknown mode: %r" % (mode,))

def _get_first_words(filename):
    infile = _open(filename, "rb")
    size = int(infile.read(4))
    data = []
    while 1:
        s = infile.read(size)
        if not s:
            break
        assert len(s) == size, (repr(s), size)
        s = s.split("\t")[0]
        if not data or data[-1] != s:
            data.append(s)
    return data

class PrimaryNamespace(BaseDB.DictLookup):
    def __init__(self, db, namespace):
        self.db = db
        self.namespace = namespace
    def __getitem__(self, name):
        loc = _lookup_location(self.db.key_filename, name)
        if loc is None:
            raise KeyError("Cannot find primary key %r" % (name,))
        data = [
            Location.Location(self.namespace,
                              name,
                              self.db.fileid_info[loc[0]][0],
                              loc[1],
                              loc[2])
        ]
        return data
    def keys(self):
        return _get_first_words(self.db.key_filename)
    
        
class SecondaryNamespace(BaseDB.DictLookup):
    def __init__(self, db, namespace):
        self.db = db
        self.namespace = namespace
    def __getitem__(self, name):
        id_filename = os.path.join(self.db.dbname,
                                   "id_%s.index" % self.namespace)
        primary_keys = _lookup_alias(id_filename, name)
        if primary_keys is None:
            raise KeyError("Cannot find %r key %r" % (self.namespace, name))

        data = []
        for key in primary_keys:
            loc = _lookup_location(self.db.key_filename, key)
            if loc is None:
                raise AssertionError("Cannot find primary key %r -- "
                                     "lost database integrety" % (key,))
            data.append(Location.Location(self.namespace, name,
                                          self.db.fileid_info[loc[0]][0],
                                          loc[1],
                                          loc[2]))
        return data

    def keys(self):
        id_filename = os.path.join(self.db.dbname,
                                   "id_%s.index" % self.namespace)
        return _get_first_words(id_filename)


class DiskFlatDB(BaseFlatDB):
    def __init__(self, dbname):
        BaseFlatDB.__init__(self, dbname)
        

    def __getitem__(self, namespace):
        """return the database table lookup for the given namespace"""
        if namespace == self.primary_namespace:
            return PrimaryNamespace(self, namespace)
        if namespace in self.secondary_namespaces:
            return SecondaryNamespace(self, namespace)
        raise KeyError(namespace)
