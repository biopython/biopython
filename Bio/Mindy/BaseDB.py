import os
import Bio
import compression

def _int_str(i):
    s = str(i)
    if s[-1:] == "l":
        return s[:-1]
    return s

class CreateDB:
    def __init__(self, dbname, unique, data_fields, format = "sequence"):
        # Must define 'self.format' being a Format object
        # Must define 'self.fileids' mapping from filename -> fileid
        # Must define 'self.filemap' mapping from fileid -> filename \t size
        raise NotImplementedError

    def add_filename(self, filename, size):
        fileid = self.fileids.get(filename, None)
        if fileid is not None:
            return fileid
        s = str(len(self.fileids))
        self.fileids[filename] = s
        self.filemap[s] = "%s\t%s" % (filename, _int_str(size))
        return s

    def load(self, filename, builder, record_tag = "record"):
        size = os.path.getsize(filename)
        filetag = self.add_filename(filename, size)

        source = compression.open_file(filename, "rb")
        format = self.format.identifyFile(source)
        if format is None:
            raise TypeError("Cannot identify file as a %s format" %
                            (origformat,))
        
        iterator = format.make_iterator(
            record_tag,
            select_names = tuple(builder.uses_tags()) + (record_tag,),
            debug_level = 0)

        for record in iterator.iterate(source, cont_handler = builder):
            self.add_record(filetag,
                            iterator.start_position,
                            iterator.end_position - iterator.start_position,
                            record.document)
        
class OpenDB:
    def __init__(self, dbname):
        self.dbname = dbname

    def lookup(self, *args, **kwargs):
        if args:
            if kwargs:
                raise TypeError("Cannot specify both args and kwargs")
            if len(args) != 1:
                raise TypeError("Only one identifier handled")
            namespace, name = self.unique, args[0]
        
        else:
            if len(kwargs) != 1:
                raise TypeError("lookup takes a single key")
            namespace, name = kwargs.items()[0]
        return self[namespace][name]

    def __getitem__(self, name):
        raise NotImplementedError("must be implemented in the derived class")

    # Introspection
    def keys(self):
        return [self.primary_namespace] + self.secondary_namespaces
    def values(self):
        return [self[key] for key in self.keys()]
    def items(self):
        return [(key, self[key]) for key in self.keys()]
        
        

