
import os, bisect
import BaseDB, Location
import Bio

_open = open

class CreateFlatDB(BaseDB.CreateDB):
    def __init__(self, dbname, primary_namespace,
                 secondary_namespaces, format):
        self.dbname = dbname        
        self.primary_namespace = primary_namespace
        self.secondary_namespaces = secondary_namespaces
        self.format = Bio.formats.normalize(format)

        # This is used in __del__ to tell if the file has been closed.
        # Set it to None so that test can work if the mkdir or open fail.
        # After the directory is made, we'll set it to {}
        self.fileids = None
        self.filemap = {}

        self.primary_table = {}
        self.lookup_tables = {}
        for namespace in secondary_namespaces:
            self.lookup_tables[namespace] = {}

        os.mkdir(dbname)
        outfile = _open(os.path.join(dbname, "config.dat"), "wb")
        outfile.write("index\tflat/1\n")
        outfile.close()

        self.fileids = {}

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
        self.primary_table[key] = "%s\t%s\t%s" % (filetag, startpos, length)

        for namespace in secondary_namespaces:
            lookup = self.lookup_tables[namespace]
            # Get the list of secondary identifiers for this identifier
            for val in table.get(field, ()):
                # Go from secondary identifier to list of primary identifiers
                lookup.setdefault(val, []).append(key)

    def close(self):
        # Write the configuration
        configfile = _open(os.path.join(self.dbname, "config.dat"), "ab")

        # Write the namespace information
        configfile.write("primary_namespace\t%s\n" % self.primary_namespace)
        keys = self.lookup_tables.keys()
        keys.sort()
        configfile.write("secondary_namespaces\t")
        configfile.write("\t".join(keys) + "\n")
        
        # Write the fileid table
        items = self.filemap.items()
        items.sort()
        for fileid, s in items:
            configfile.write("fileid_%s\t%s\n" % (fileid, s))

        configfile.close()

        # Write the primary identifier information, which is fixed width
        filename = os.path.join(self.dbname, "key_%s.key" % (self.primary_namespace,) )
        info = self.primary_table.items()
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

        # Write the secondary identifier information
        for field, table in self.lookup_tables.items():
            
            filename = os.path.join(self.dbname, "id_%s.index" % field)
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

        self.primary_table = self.fileids = self.lookup_tables = None

    def __del__(self):
        if self.fileids is not None:
            self.close()



def _read_tab_dict(filename):
    d = {}
    for line in _open(filename, "rb").read().split("\n"):
        words = line.rstrip().split("\t")
        d[words[0]] = words[1:]
    return d

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
    def get_line(self, i):
        self.infile.seek(i * self.record_size + 4)
        return self.infile.read(self.record_size)

def _find_line(filename, wantword):
    size = os.path.getsize(filename)
    infile = _open(filename, "rb")

    bf = BisectFile(infile, size)
    left = bisect.bisect_left(bf, wantword)
    line = bf.get_line(left)
    if not line.startswith(wantword):
        return None
    return line

def _find_range(filename, wantword):
    size = os.path.getsize(filename)
    infile = _open(filename, "rb")

    bf = BisectFile(infile, size)
    left = bisect.bisect_left(bf, wantword)
    line = bf.get_line(left)
    if not line.startswith(wantword):
        return None
    
    right = bisect.bisect_right(bf, wantword)
    data = []
    for i in range(left, right):
        x = bf.get_line(i)
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
        if config["index"][0] != "flat/1":
            raise TypeError("FlatDB does not support %r index" %
                            (config["index"][0],))
        self.primary_namespace = config["primary_namespace"][0]
        self.secondary_namespaces = config["secondary_namespaces"]

        fileids = {}
        for k, v in config.items():
            if not k.startswith("fileid_"):
                continue
            fileid = k[7:]
            filename, size = v
            size = long(size)
            fileids[fileid] = filename, size
            if os.path.getsize(filename) != size:
                raise TypeError(
                    "File %s has changed size from %d to %d bytes!" %
                    (size, os.path.getsize(filename)))

        self.fileids = fileids
        self.dbname = dbname

        self.key_filename = os.path.join(dbname,
                                         "key_%s.key" % self.primary_namespace)

    def _turn_into_query(self, *args, **kwargs):
        if args:
            if kwargs:
                raise TypeError("Cannot specify both args and kwargs")
            if len(args) != 1:
                raise TypeError("Only one identifier handled")
            return self.primary_namespace, args[0]
        
        if len(kwargs) != 1:
            raise TypeError("lookup takes a single key")
        return kwargs.items()[0]

    def lookup(self, *args, **kwargs):
        query = self._turn_into_query(*args, **kwargs)
        namespace, name = query
        if namespace == self.primary_namespace:
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



