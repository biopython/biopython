"""I/O function wrappers for the RDF/CDAO file format.

"""
__docformat__ = "restructuredtext en"

from cStringIO import StringIO

from Bio.Phylo import Newick
import os
import RDF
#RDF.debug(1)


# Definitions retrieved from Bio.Nexus.Trees
NODECOMMENT_START = '[&'
NODECOMMENT_END = ']'


class CDAOError(Exception):
    """Exception raised when CDAO object construction cannot continue."""
    pass


# ---------------------------------------------------------
# Public API

def parse(handle, **kwargs):
    """Iterate over the trees in a CDAO file handle.

    :returns: generator of Bio.Phylo.Newick.Tree objects.
    """
    return Parser(handle).parse(**kwargs)


def write(trees, handle, plain=False, **kwargs):
    """Write a trees in CDAO format to the given file handle.

    :returns: number of trees written.
    """
    return Writer(trees).write(handle, plain=plain, **kwargs)


# ---------------------------------------------------------
# Input

class Parser(object):
    """Parse a CDAO tree given a file handle.

    Based on the parser in `Bio.Nexus.Trees`.
    """

    def __init__(self, handle):
        self.handle = handle
        self.model = None

    @classmethod
    def from_string(cls, treetext):
        handle = StringIO(treetext)
        return cls(handle)

    def parse(self, values_are_confidence=False, rooted=False,
              storage=None, mime_type='text/turtle'):
        """Parse the text stream this object was initialized with."""
        if storage is None:
            # store RDF model in memory for now
            storage = RDF.Storage(storage_name="hashes",
                                  name="serializer",
                                  options_string="new='yes',hash-type='memory',dir='.'")
            if storage is None:
                raise CDAOError("new RDF.Storage failed")

        if self.model is None:
            self.model = RDF.Model(storage)
            if self.model is None:
                raise CDAOError("new RDF.model failed")
        model = self.model
        
        self.values_are_confidence = values_are_confidence
        self.rooted = rooted
        
        parser = RDF.Parser(mime_type=mime_type)
        if parser is None:
            raise Exception('Failed to create RDF.Parser for MIME type %s' % mime_type)
        
        uri = RDF.Uri(string="file:"+self.handle.name)
        for s in parser.parse_string_as_stream(self.handle.read(), uri):
            model.append(s)
        
        # TODO: create a Tree object from RDF model


# ---------------------------------------------------------
# Output

class Writer(object):
    """Based on the writer in Bio.Nexus.Trees (str, to_string)."""
    urls = {
            'owl': 'http://www.w3.org/2002/07/owl#',
            'cdao': 'http://purl.obolibrary.org/obo/cdao.owl#',
            'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
            }

    def __init__(self, trees):
        self.trees = trees
        self.model = None
        self.count = 0

    def write(self, handle, **kwargs):
        """Write this instance's trees to a file handle.
        
        Keywords:
            mime_type: used to determine the serialization format.
                default is 'text/turtle'
        """
        try: mime_type = kwargs['mime_type']
        except KeyError: mime_type = 'text/turtle'
        
        self.add_trees_to_model()
        self.serialize_model(handle, mime_type=mime_type)
        
    def add_trees_to_model(self, trees=None, storage=None):
        """Add triples describing a set of trees to an RDF model."""
        Uri = RDF.Uri
        urls = self.urls
        
        def qUri(s):
            if ':' in s:
                s = s.split(':')
                s = urls[s[0]] + ':'.join(s[1:])
            return Uri(s)
            
        if trees is None:
            trees = self.trees
        
        if storage is None:
            # store RDF model in memory for now
            storage = RDF.Storage(storage_name="hashes",
                                  name="serializer",
                                  options_string="new='yes',hash-type='memory',dir='.'")
            if storage is None:
                raise CDAOError("new RDF.Storage failed")

        if self.model is None:
            self.model = RDF.Model(storage)
            if self.model is None:
                raise CDAOError("new RDF.model failed")
        model = self.model
        
        def add_statements(statements):
            for stmt in statements:
                model.append(RDF.Statement(*stmt))
        
        def process_clade(clade, parent=None):
            # TODO: what uri to use for clades?
            # currently, non-terminal nodes are all just "clade"
            clade.uri = (str(clade))
            statements = [
                          (Uri(clade.uri), qUri('rdf:type'), qUri('cdao:Node')),
                          ]
            if clade.name:
                # TODO: create TU
                pass
                          
            if not parent is None:
                edge_uri = parent.uri + "_" + clade.uri
                statements += [
                               (Uri(edge_uri), qUri('rdf:type'), qUri('cdao:Directed_Edge')),
                               (Uri(edge_uri), qUri('cdao:has_Parent_Node'), Uri(parent.uri)),
                               (Uri(edge_uri), qUri('cdao:has_Child_Node'), Uri(clade.uri)),
                               (Uri(clade.uri), qUri('cdao:belongs_to_Edge_as_Child'), Uri(edge_uri)),
                               (Uri(clade.uri), qUri('cdao:has_Parent'), Uri(parent.uri)),
                               (Uri(parent.uri), qUri('cdao:belongs_to_Edge_as_Parent'), Uri(edge_uri)),
                               ]
                          
            add_statements(statements)
            
            if not clade.is_terminal():
                for new_clade in clade.clades:
                    process_clade(new_clade, parent=clade)
        
        add_statements([
                        (Uri(urls['cdao']), qUri('rdf:type'), qUri('owl:Ontology')),
                        ])

        for tree in trees:
            first_clade = tree.clade
            process_clade(first_clade)
            
            self.count += 1
            
        model.sync()
            
    def serialize_model(self, handle, mime_type='text/turtle'):
        """Serialize RDF model to file handle"""        
        # serialize RDF model to output file
        serializer = RDF.Serializer(mime_type=mime_type)
        for prefix, url in self.urls.items():
            serializer.set_namespace(prefix, url)

        # TODO: this is going to be too memory intensive for large trees;
        # come up with something better
        print "Writing..."
        serializer.serialize_model_to_file(handle.name, self.model)
        print "Done."
        
        return self.count
