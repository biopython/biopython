import os
import Bio
import compression

def _int_str(i):
    s = str(i)
    if s[-1:] == "l":
        return s[:-1]
    return s

class WriteDB:
    # Must define 'self.filename_map' mapping from filename -> fileid
    # Must define 'self.fileid_info' mapping from fileid -> (filename,size)

    def add_filename(self, filename, size, fileid_info):
        fileid = self.filename_map.get(filename, None)
        if fileid is not None:
            return fileid
        s = str(len(self.filename_map))
        self.filename_map[filename] = s  # map from filename -> id
        assert s not in fileid_info.keys(), "Duplicate entry! %s" % (s,)
        self.fileid_info[s] = (filename, size)
        return s

    def load(self, filename, builder, fileid_info, record_tag = "record"):
        formatname = self.formatname
        size = os.path.getsize(filename)
        filetag = self.add_filename(filename, size, fileid_info)

        source = compression.open_file(filename, "rb")
        if formatname == "unknown":
            formatname = "sequence"
        
        format = Bio.formats.normalize(formatname).identifyFile(source)
        if format is None:
            raise TypeError("Cannot identify file as a %s format" %
                            (self.formatname,))
        if self.formatname == "unknown":
            expected_names = ["fasta", "embl", "swissprot", "genbank"]
            for node in format._parents:
                if node.name in expected_names:
                    self.formatname = node.name
                    break
            else:
                self.formatname = format.name
        
        iterator = format.make_iterator(
            record_tag,
            select_names = tuple(builder.uses_tags()) + (record_tag,),
            debug_level = 0)

        for record in iterator.iterate(source, cont_handler = builder):
            self.add_record(filetag,
                            iterator.start_position,
                            iterator.end_position - iterator.start_position,
                            record.document)

class DictLookup:
    def __getitem__(self, key):
        raise NotImplementedError("Must be implemented in subclass")
    def keys(self):
        raise NotImplementedError("Must be implemented in subclass")

    def values(self):
        return [self[key] for key in self.keys()]
    def items(self):
        return [(key, self[key]) for key in self.keys()]

    def get(self, key, default = None):
        try:
            return self[key]
        except KeyError:
            return default
    
        
class OpenDB(DictLookup):
    def __init__(self, dbname, index_type):
        self.dbname = dbname

        config = read_config(os.path.join(dbname, "config.dat"))
        if config["index"] != index_type:
            raise TypeError("FlatDB does not support %r index" %
                            (config["index"],))
        self.primary_namespace = config["primary_namespace"]
        self.secondary_namespaces = config["secondary_namespaces"]
        self.formatname = config["format"]

        filename_map = {}
        fileid_info = {}
        for k, v in config.items():
            if not k.startswith("fileid_"):
                continue
            fileid = k[7:]
            filename, size = v
            fileid_info[fileid] = v
            filename_map[filename] = fileid
            if os.path.getsize(filename) != size:
                raise TypeError(
                    "File %s has changed size from %d to %d bytes!" %
                    (size, os.path.getsize(filename)))

        self.filename_map = filename_map
        self.fileid_info = fileid_info
        

    def lookup(self, *args, **kwargs):
        if args:
            if kwargs:
                raise TypeError("Cannot specify both args and kwargs")
            if len(args) != 1:
                raise TypeError("Only one identifier handled")
            namespace, name = self.primary_namespace, args[0]
        
        else:
            if len(kwargs) != 1:
                raise TypeError("lookup takes a single key")
            namespace, name = kwargs.items()[0]
        return self[namespace][name]

    def __getitem__(self, namespace):
        """return the database table lookup for the given namespace"""
        raise NotImplementedError("must be implemented in the derived class")

    def keys(self):
        return [self.primary_namespace] + self.secondary_namespaces

# Write the configuration
def write_config(config_filename,
                 index_type,
                 primary_namespace,
                 secondary_namespaces,
                 fileid_info,
                 formatname):
    configfile = open(config_filename, "wb")

    # Write the header
    configfile.write("index\t" + index_type + "\n")

    # Write the namespace information
    configfile.write("primary_namespace\t%s\n" % primary_namespace)
    keys = secondary_namespaces[:]
    keys.sort()
    configfile.write("secondary_namespaces\t")
    configfile.write("\t".join(keys) + "\n")

    # Format name
    configfile.write("format\t" + formatname + "\n")

    # Write the fileid table
    items = fileid_info.items()
    items.sort()
    for fileid, (filename, size) in items:
        configfile.write("fileid_%s\t%s\t%s\n" % \
                         (fileid, filename, _int_str(size)))

    configfile.close()


def read_config(config_filename):
    d = {}
    for line in open(config_filename, "rb").read().split("\n"):
        words = line.rstrip().split("\t")
        assert not d.has_key(words[0]), \
               "Duplicate key %r in config file: old = %r, new = %r" % \
               (words[0], d[words[0]], line)
        if words[0] in ("index", "primary_namespace", "format"):
            if len(words) != 2:
                raise AssertionError(
                    "%s should only have one value, not %r" % \
                    (words[0], words[1:]))
            d[words[0]] = words[1]
            
        elif words[0].startswith("fileid_"):
            if len(words) != 3:
                raise AssertionError(
                    "%s should only have two values, not %r" % \
                    (words[0], words[1:]))
            d[words[0]] = (words[1], long(words[2]))
        
        elif words[0] in ("secondary_namespaces",):
            # This can have 0 or more values
            d[words[0]] = words[1:]
        
        else:
            # Unknown word, save as-is
            d[words[0]] = words[1:]
    
    return d
