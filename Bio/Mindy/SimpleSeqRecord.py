"""Index a file based on information in a SeqRecord object.

This indexer tries to make it simple to index a file of records (ie. like
a GenBank file full of entries) so that individual records can be 
readily retrieved.

The indexing in this file takes place by converting the elements in the
file into SeqRecord objects, and then indexing by some item in these 
SeqRecords. This is a slower method, but is very flexible.

We have two default functions to index by the id and name elements of a
SeqRecord (ie. LOCUS and accession number from GenBank). There is also 
a base class which you can derive from to create your own indexer
which allows you to index by anything you feel like using python code.
"""
from Bio.builders.SeqRecord.sequence import BuildSeqRecord

# --- base class for indexing using SeqRecord information
class BaseSeqRecordIndexer:
    """Base class for indexing using SeqRecord information.

    This is the class you should derive from to index using some type of
    information in a SeqRecord. This is an abstract base class, so it needs
    to be subclassed to be useful.
    """
    def __init__(self):
        pass

    def get_builder(self):
        tricky_builder = FixDocumentBuilder(self.get_id_dictionary)
        return tricky_builder

    def primary_key_name(self):
        raise NotImplementedError("Please implement in derived classes")
    
    def secondary_key_names(self):
        raise NotImplementedError("Please implement in derived classes")

    def get_id_dictionary(self, seq_record):
        raise NotImplementedError("Please implement in derived classes")

class SimpleIndexer(BaseSeqRecordIndexer):
    """Index a file based on .id and .name attributes of a SeqRecord.

    A simple-minded indexing scheme which should work for simple cases. The
    easiest way to use this is trhough the create_*db functions of this
    module.
    """
    def __init__(self):
        BaseSeqRecordIndexer.__init__(self)

    def primary_key_name(self):
        return "id"

    def secondary_key_names(self):
        return ["name", "aliases"]

    def get_id_dictionary(self, seq_record):
        # XXX implement aliases once we have this attribute in SeqRecords
        id_info = {"id" : [seq_record.id],
                   "name" : [seq_record.name],
                   "aliases" : []}
        return id_info

class FunctionIndexer(BaseSeqRecordIndexer):
    """Indexer to index based on values returned by a function.
    
    This class is passed a function which will return id, name and alias
    information from a SeqRecord object. It needs to return either one item, 
    which is an id from the title, or three items which are (in order), the id, 
    a list of names, and a list of aliases.

    This indexer allows indexing to be completely flexible based on passed
    functions.
    """
    def __init__(self, index_function):
        BaseSeqRecordIndexer.__init__(self)
        self.index_function = index_function

    def primary_key_name(self):
        return "id"

    def secondary_key_names(self):
        return ["name", "aliases"]

    def get_id_dictionary(self, seq_record):
        items = self.index_function(seq_record)
        if type(items) is not type([]) and type(items) is not type(()):
            items = [items]
        if len(items) == 1:
            seq_id = items[0]
            name = []
            aliases = []
        elif len(items) == 3:
            seq_id, name, aliases = items
        else:
            raise ValueError("Unexpected items from index function: %s" %
                    (items))
        
        return {"id" : [seq_id],
                "name" : name,
                "aliases" : aliases}

class FixDocumentBuilder(BuildSeqRecord):
    """A SAX builder-style class to make a parsed SeqRecord available.

    This class does a lot of trickery to make things fit in the SAX
    framework and still have the flexibility to use a built SeqRecord
    object. 

    You shouldn't really need to use this class unless you are doing 
    something really fancy-pants; otherwise, just use the 
    BaseSeqRecordIndexer interfaces.
    """
    def __init__(self, get_ids_callback):
        """Intialize with a callback function to gets id info from a SeqRecord.

        get_ids_callback should be a callable function that will take a 
        SeqRecord object and return a dictionary mapping id names to 
        the valid ids for these names.
        """
        BuildSeqRecord.__init__(self)
        self._ids_callback = get_ids_callback

    def end_record(self, tag):
        """Overrride the builder function to muck with the document attribute.
        """
        # first build up the SeqRecord
        BuildSeqRecord.end_record(self, tag)
        # now convert the SeqRecord into the dictionary that the indexer needs
        self.document = self._ids_callback(self.document)

# --- convenience functions for indexing
# you should just use these unless you are doing something fancy
def create_berkeleydb(files, db_name, indexer = SimpleIndexer()):
    from Bio.Mindy import BerkeleyDB
    unique_name = indexer.primary_key_name()
    alias_names = indexer.secondary_key_names()
    creator = BerkeleyDB.create(db_name, unique_name, alias_names)
    builder = indexer.get_builder()
    for filename in files:
        creator.load(filename, builder = builder, fileid_info = {})
    creator.close()

def create_flatdb(files, db_name, indexer = SimpleIndexer()):
    from Bio.Mindy import FlatDB
    unique_name = indexer.primary_key_name()
    alias_names = indexer.secondary_key_names()
    creator = FlatDB.create(db_name, unique_name, alias_names)
    builder = indexer.get_builder()
    for filename in files:
        creator.load(filename, builder = builder, fileid_info = {})
    creator.close()
