"""I/O function wrappers for the RDF/CDAO file format.

"""
__docformat__ = "restructuredtext en"

from cStringIO import StringIO

from Bio.Phylo import Newick
import os


class CDAOError(Exception):
    """Exception raised when CDAO object construction cannot continue."""
    pass


def import_rdf():
    try: import RDF
    except ImportError: raise CDAOError('Redland Python bindings are required for CDAO support.')
    #RDF.debug(1)
    return RDF
        
    
def new_storage():
    RDF = import_rdf()
    
    storage = RDF.Storage(storage_name="hashes",
                          name="serializer",
                          options_string="new='yes',hash-type='memory',dir='.'")
    if storage is None:
        raise CDAOError("new RDF.Storage failed")
    return storage


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
        RDF = import_rdf()

        if storage is None:
            # store RDF model in memory for now
            storage = new_storage()

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
        # get all cdao:RootedTree instances, then start tree creation at the 
        # node designated by cdao:has_root
        


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
        
        self.node_counter = 0
        self.edge_counter = 0
        self.tree_counter = 0
        self.tu_counter = 0

    def write(self, handle, **kwargs):
        """Write this instance's trees to a file handle.
        
        Keywords:
            mime_type: used to determine the serialization format.
                default is 'text/turtle'
        """
        RDF = import_rdf()
        
        try: mime_type = kwargs['mime_type']
        except KeyError: mime_type = 'text/turtle'
        
        self.add_trees_to_model()
        self.serialize_model(handle, mime_type=mime_type)
        
    def add_trees_to_model(self, trees=None, storage=None):
        """Add triples describing a set of trees to an RDF model."""
        # TODO: base uri?
        RDF = import_rdf()
        
        Uri = RDF.Uri
        urls = self.urls
        
        def qUri(s):
            '''returns the full URI from a namespaced URI string (i.e. rdf:type)'''
            for url in urls: 
                s = s.replace(url+':', urls[url])
            return Uri(s)
            
        if trees is None:
            trees = self.trees
        
        if storage is None:
            # store RDF model in memory for now
            storage = new_storage()

        if self.model is None:
            self.model = RDF.Model(storage)
            if self.model is None:
                raise CDAOError("new RDF.model failed")
        model = self.model
        
        def add_statements(statements):
            '''add RDF statements, represented by triples, to an RDF model'''
            for stmt in statements:
                model.append(RDF.Statement(*stmt))
        
        def process_clade(clade, parent=None, root=False):
            '''recursively add statements describing a tree of clades to the
            RDF model'''
            import uuid
            self.node_counter += 1
            clade.uri = 'node%s' % self.node_counter
            
            statements = []
            if root:
                # create a cdao:RootedTree with reference to the tree root
                self.tree_counter += 1
                tree_uri = 'tree%s' % self.tree_counter
                statements += [
                               (Uri(tree_uri), qUri('rdf:type'), qUri('cdao:RootedTree')),
                               (Uri(tree_uri), qUri('cdao:has_root'), Uri(clade.uri)),
                               ]
            
            if clade.name:
                # create TU
                self.tu_counter += 1
                tu_uri = 'tu%s' % self.tu_counter
                statements += [
                               (Uri(tu_uri), qUri('rdf:type'), qUri('cdao:TU')),
                               (Uri(clade.uri), qUri('cdao:represents_TU'), Uri(tu_uri)),
                               (Uri(tu_uri), qUri('rdf:label'), clade.name),
                               ]
                               
                # TODO: should be able to pass in an optional function for 
                # running each TU through TNRS, etc.
                
            # create this node
            node_type = 'cdao:TerminalNode' if clade.is_terminal() else 'cdao:AncestralNode'
            statements += [
                           (Uri(clade.uri), qUri('rdf:type'), qUri(node_type)),
                           ]
                          
            if not parent is None:
                # create edge from the parent node to this node
                self.edge_counter += 1
                edge_uri = 'edge%s' % self.edge_counter
                statements += [
                               (Uri(edge_uri), qUri('rdf:type'), qUri('cdao:Directed_Edge')),
                               (Uri(edge_uri), qUri('cdao:has_Parent_Node'), Uri(parent.uri)),
                               (Uri(edge_uri), qUri('cdao:has_Child_Node'), Uri(clade.uri)),
                               (Uri(clade.uri), qUri('cdao:belongs_to_Edge_as_Child'), Uri(edge_uri)),
                               (Uri(clade.uri), qUri('cdao:has_Parent'), Uri(parent.uri)),
                               (Uri(parent.uri), qUri('cdao:belongs_to_Edge_as_Parent'), Uri(edge_uri)),
                               ]
                # TODO: add edge annotations (i.e. edge lengths)
                          
            add_statements(statements)
            
            if not clade.is_terminal():
                for new_clade in clade.clades:
                    process_clade(new_clade, parent=clade, root=False)
        
        add_statements([
                        (Uri(urls['cdao']), qUri('rdf:type'), qUri('owl:Ontology')),
                        ])

        for tree in trees:
            first_clade = tree.clade
            process_clade(first_clade, root=True)
            
        model.sync()
            
    def serialize_model(self, handle, mime_type='text/turtle'):
        """Serialize RDF model to file handle"""        
        RDF = import_rdf()
        
        # serialize RDF model to output file
        serializer = RDF.Serializer(mime_type=mime_type)
        for prefix, url in self.urls.items():
            serializer.set_namespace(prefix, url)

        handle.write(serializer.serialize_model_to_string(self.model))
        
        return self.tree_counter
