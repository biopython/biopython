import os, bisect
import BaseDB, Location
import Bio

_open = open

class CreateFlatDB(BaseDB.CreateDB):
    def __init__(self, dbname, unique, data_fields, format):
        self.dbname = dbname        
        self.unique = unique
        self.data_fields = data_fields
        self.format = Bio.formats.normalize(format)

        # This is used in __del__ to tell if the file has been closed.
        # Set it to None so that test can work if the mkdir or open fail.
        # After the directory is made, we'll set it to {}
        self.fileids = None
        self.filemap = {}

        self.primary_table = {}
        self.lookup_tables = {}
        for field in data_fields:
            self.lookup_tables[field] = {}

        os.mkdir(dbname)
        outfile = _open(os.path.join(dbname, "BIOINDEX.dat"), "w")
        outfile.write("index\tflat/1\n")
        outfile.close()

        self.fileids = {}

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
                lookup.setdefault(val, []).append(key)

    def close(self):
        # Write the fileid table
        outfile = _open(os.path.join(self.dbname, "fileids.dat"), "w")
        for tag, s in self.filemap.items():
            outfile.write("%s\t%s\n" % (tag, s))
        outfile.close()

        # Write the configuration
        outfile = _open(os.path.join(self.dbname, "config.dat"), "w")
        outfile.write("namespace\t%s\n" % self.unique)
        outfile.close()

        # Write the primary identifier information, which is fixed width
        filename = os.path.join(self.dbname, "key_%s.key" % (self.unique,) )
        outfile = _open(filename, "w")
        info = self.primary_table.items()
        info.sort()
        n = 1
        for k, v in info:
            # Find the longest width
            s = "%s\t%s" % (k, v)
            if len(s) > n:
                n = len(s)
        for k, v in info:
            s = "%s\t%s" % (k, v)
            outfile.write(s.ljust(n) + "\n")
        outfile.close()

        # Write the secondary identifier information
        for field, table in self.lookup_tables.items():
            
            filename = os.path.join(self.dbname, "id_%s.index" % field)
            outfile = _open(filename, "w")
            items = table.items()
            items.sort()
            # Find the largest field
            n = 0
            for k, v in items:
                for x in v:
                    s = "%s\t%s" % (k, x)
                    if len(s) > n:
                        n = len(s)
            # And write the output
            for k, v in items:
                for x in v:
                    s = "%s\t%s" % (k, x)
                    outfile.write(s.ljust(n) + "\n")
            outfile.close()

        self.primary_table = self.fileids = self.lookup_tables = None

    def __del__(self):
        if self.fileids is not None:
            self.close()


def _read_tab_dict(filename):
    d = {}
    for line in _open(filename).readlines():
        words = line.rstrip().split("\t")
        d[words[0]] = words[1:]
    return d

class BisectFile:
    def __init__(self, infile, size):
        self.infile = infile
        self.size = size
        infile.seek(0)
        s = infile.readline()
        if not s:
            self.record_size = None
        else:
            self.record_size = len(s)
        n = int(self.size / self.record_size)
        assert n * self.record_size == self.size, "record size is wrong"
    def __len__(self):
        if self.record_size is None:
            return 0
        return int(self.size / self.record_size)
    def __getitem__(self, i):
        self.infile.seek(i * self.record_size)
        return self.infile.readline().split("\t")[0]
    def getline(self, i):
        self.infile.seek(i * self.record_size)
        return self.infile.readline()

def _find_line(filename, wantword):
    size = os.path.getsize(filename)
    infile = _open(filename)

    bf = BisectFile(infile, size)
    left = bisect.bisect_left(bf, wantword)
    line = bf.getline(left)
    if not line.startswith(wantword):
        return None
    return line

def _find_range(filename, wantword):
    size = os.path.getsize(filename)
    infile = _open(filename)

    bf = BisectFile(infile, size)
    left = bisect.bisect_left(bf, wantword)
    line = bf.getline(left)
    if not line.startswith(wantword):
        return None
    right = bisect.bisect_right(bf, wantword)
    data = []
    for i in range(left, right):
        x = bf.getline(i)
        data.append(x)
    return data


def _lookup_location(key_filename, word):
    line = _find_line(key_filename, word)
    if line is None:
        return None
    line = line.rstrip()  # get rid of the packing and newline
    words = line.split("\t")[1:4]  # Get fileid, start, length
    return words[0], long(words[1]), long(words[2])

def _lookup_alias(id_filename, word):
    lines = _find_range(id_filename, word)
    if not lines:
        return None
    primary_keys = []
    for line in lines:
        # get rid of newline (and spaces? spec says there are none)
        line = line.rstrip()
        # There are only two fields
        found_word, primary_key = line.split("\t")
        assert found_word == word, (found_word, word)
        primary_keys.append(primary_key)
    return primary_keys

def open(dbname):
    return OpenFlatDB(dbname)

class OpenFlatDB(BaseDB.OpenDB):
    def __init__(self, dbname):
        BaseDB.OpenDB.__init__(self, dbname)
        
        config = _read_tab_dict(os.path.join(dbname, "config.dat"))
        fileids = _read_tab_dict(os.path.join(dbname, "fileids.dat"))
        for k, v in fileids.items():
            filename, size = v
            size = long(size)
            fileids[k] = filename, size
            if os.path.getsize(filename) != size:
                raise TypeError(
                    "File %s has changed size from %d to %d bytes!" %
                    (size, os.path.getsize(filename)))

        self.unique = config["namespace"][0]
        self.fileids = fileids
        self.dbname = dbname

        self.key_filename = os.path.join(dbname,
                                         "key_%s.key" % self.unique)

    def _turn_into_query(self, *args, **kwargs):
        if args:
            if kwargs:
                raise TypeError("Cannot specify both args and kwargs")
            if len(args) != 1:
                raise TypeError("Only one identifier handled")
            return self.unique, args[0]
        
        if len(kwargs) != 1:
            raise TypeError("lookup takes a single key")
        return kwargs.items()[0]

    def lookup(self, *args, **kwargs):
        query = self._turn_into_query(*args, **kwargs)
        namespace, name = query
        if namespace == self.unique:
            loc = _lookup_location(self.key_filename, name)
            if loc is None:
                raise KeyError("Cannot find primary key %r" % (name,))
            data = [
                Location.Location(namespace,
                                  name,
                                  self.fileids[loc[0]][0],
                                  loc[1],
                                  loc[2])
            ]
            return data

        # There can be one or more aliases
        id_filename = os.path.join(self.dbname,
                                   "id_%s.index" % namespace)
        primary_keys = _lookup_alias(id_filename, name)
        if primary_keys is None:
            raise KeyError("Cannot find %r key %r" % (namespace, name))

        data = []
        for key in primary_keys:
            loc = _lookup_location(self.key_filename, key)
            if loc is None:
                raise AssertionError("Cannot find primary key %r -- "
                                     "lost database integrety" % (key,))
            data.append(Location.Location(namespace, name,
                                          self.fileids[loc[0]][0],
                                          loc[1],
                                          loc[2]))
        return data



