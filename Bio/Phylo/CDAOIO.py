"""I/O function wrappers for the RDF/CDAO file format.

"""
__docformat__ = "restructuredtext en"

from cStringIO import StringIO

from Bio.Phylo import Newick
import os


RDF_NAMESPACES = {
                  'owl': 'http://www.w3.org/2002/07/owl#',
                  'cdao': 'http://purl.obolibrary.org/obo/cdao.owl#',
                  'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
                  }

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
                          name="test",
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
    """
    urls = RDF_NAMESPACES

    def __init__(self, handle):
        self.handle = handle
        self.model = None
        self.node_info = None
        self.children = {}

    @classmethod
    def from_string(cls, treetext):
        handle = StringIO(treetext)
        return cls(handle)

    def parse(self, **kwargs):
        """Parse the text stream this object was initialized with."""
        self.parse_handle_to_model(**kwargs)
        return self.parse_model()
        
    def parse_handle_to_model(self, rooted=False, storage=None, 
                              mime_type='text/turtle', **kwargs):
        '''Parse self.handle into RDF model self.model.'''
        RDF = import_rdf()

        if storage is None:
            # store RDF model in memory for now
            storage = new_storage()

        if self.model is None:
            self.model = RDF.Model(storage)
            if self.model is None:
                raise CDAOError("new RDF.model failed")
        model = self.model
        
        self.rooted = rooted
        
        parser = RDF.Parser(mime_type=mime_type)
        if parser is None:
            raise Exception('Failed to create RDF.Parser for MIME type %s' % mime_type)
        
        if 'base_uri' in kwargs: base_uri = kwargs['base_uri']
        else: base_uri = RDF.Uri(string="file://"+os.path.abspath(self.handle.name))
        
        statements = parser.parse_string_as_stream(self.handle.read(), base_uri)
        for s in statements:
            model.append(s)
            
        return self.parse_model(model)
            
            
    def parse_model(self, model=None):
        '''Construct a Newick.Tree from an RDF model.'''
        RDF = import_rdf()
        
        if model is None:
            model = self.model
        
        # look up branch lengths/TUs for all nodes
        self.get_node_info(model)
        
        for root_node in self.tree_roots:
            clade = self.parse_children(root_node)
            
            yield Newick.Tree(root=clade, rooted=self.rooted)
            
            
    def new_clade(self, node):
        '''Returns a Newick.Clade object for a given named node.'''
        RDF = import_rdf()
        
        result = self.node_info[node]
        
        kwargs = {}
        if 'branch_length' in result: kwargs['branch_length'] = result['branch_length']
        if 'label' in result: kwargs['name'] = result['label']
        
        clade = Newick.Clade(**kwargs)
        
        return clade
        
            
    def get_node_info(self, model):
        '''Creates a dictionary containing information about all nodes in the tree.'''
        RDF = import_rdf()
        
        Uri = RDF.Uri
        
        self.node_info = {}
        self.obj_info = {}
        self.children = {}
        self.nodes = set()
        self.tree_roots = set()
        
        for statement in model:
            # process each RDF triple in the model sequentially
            s, v, o = str(statement.subject), Uri(str(statement.predicate)), str(statement.object)
            
            if not s in self.obj_info: self.obj_info[s] = {}
            this = self.obj_info[s]
            
            assignments = {
                           qUri('cdao:has_Parent'): 'parent',
                           qUri('cdao:belongs_to_Edge_as_Child'): 'edge',
                           qUri('cdao:has_annotation'): 'annotation',
                           qUri('cdao:has_value'): 'value',
                           qUri('cdao:represents_TU'): 'tu',
                           qUri('rdf:label'): 'label',
                           }
            
            try:
                # if the predicate is one we care about, store information for later
                this[assignments[v]] = o
            except KeyError: pass
            
            if v == qUri('rdf:type'):
                if Uri(o) in (qUri('cdao:AncestralNode'), qUri('cdao:TerminalNode')):
                    # this is a tree node; store it in set of all nodes
                    self.nodes.add(s)
            if v == qUri('cdao:has_root'):
                # this is a tree; store it in set of all trees
                self.tree_roots.add(o)
                    
        for node in self.nodes:
            # for each node, look up all information needed to create a Newick.Clade
            self.node_info[node] = {}
            node_info = self.node_info[node]
            
            obj = self.obj_info[node]
            if 'edge' in obj:
                # if this object points to an edge, we need a branch length from
                # the annotation on that edge
                edge = self.obj_info[obj['edge']]
                if 'annotation' in edge:
                    annotation = self.obj_info[edge['annotation']]
                    if 'value' in annotation:
                        node_info['branch_length'] = float(annotation['value'])
            
            if 'tu' in obj:
                # if this object points to a TU, we need the label of that TU
                tu = self.obj_info[obj['tu']]
                if 'label' in tu:
                    node_info['label'] = tu['label']
            
            if 'parent' in obj:
                # store this node as a child of its parent, if it has one,
                # so that the tree can be traversed from parent to children
                parent = obj['parent']
                if not parent in self.children:
                    self.children[parent] = []
                self.children[parent].append(node)
                

    def parse_children(self, node):
        '''Return a list of clades representing all children nodes of the specified
        parent node.
        
        This function calls itself recursively for each child, traversing the 
        entire tree and creating a nested structure of Newick.Clade objects.'''
        
        clade = self.new_clade(node)
        
        children = self.children[node] if node in self.children else []
        clade.clades = [self.parse_children(child_node) for child_node in children]
        
        return clade


# ---------------------------------------------------------
# Output

class Writer(object):
    """Based on the writer in Bio.Nexus.Trees (str, to_string)."""
    urls = RDF_NAMESPACES

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
        
        try: base_uri = kwargs['base_uri']
        except KeyError: base_uri = ''
        
        self.add_trees_to_model(base_uri=base_uri)
        self.serialize_model(handle, mime_type=mime_type)
        
        
    def add_trees_to_model(self, trees=None, storage=None, base_uri='http://localhost/'):
        """Add triples describing a set of trees to an RDF model."""
        RDF = import_rdf()
        import Redland
        
        nUri = lambda s: namedUri(s, self.base_uri)
        Uri = RDF.Uri
        urls = self.urls
        
        self.base_uri = base_uri
            
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
                    
        Redland.librdf_model_transaction_start(model._model)
        
        for stmt in [(Uri(urls['cdao']), qUri('rdf:type'), qUri('owl:Ontology'))]:
            model.append(RDF.Statement(*stmt))

        for tree in trees:
            first_clade = tree.clade
            statements = self.process_clade(first_clade, root=True)
            for stmt in statements:
                model.append(stmt)
                
        Redland.librdf_model_transaction_commit(model._model)
            
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
                
                
    def process_clade(self, clade, parent=None, root=False):
        '''recursively generate statements describing a tree of clades'''
        RDF = import_rdf()
        
        self.node_counter += 1
        clade.uri = 'node%s' % self.node_counter
        
        nUri = lambda s: namedUri(s, self.base_uri)
        Uri = RDF.Uri
        urls = self.urls
        
        statements = []
        if root:
            # create a cdao:RootedTree with reference to the tree root
            self.tree_counter += 1
            tree_uri = 'tree%s' % self.tree_counter
            statements += [
                           (nUri(tree_uri), qUri('rdf:type'), qUri('cdao:RootedTree')),
                           (nUri(tree_uri), qUri('cdao:has_root'), nUri(clade.uri)),
                           ]
        
        if clade.name:
            # create TU
            self.tu_counter += 1
            tu_uri = 'tu%s' % self.tu_counter
            statements += [
                           (nUri(tu_uri), qUri('rdf:type'), qUri('cdao:TU')),
                           (nUri(clade.uri), qUri('cdao:represents_TU'), nUri(tu_uri)),
                           (nUri(tu_uri), qUri('rdf:label'), clade.name),
                           ]
                           
            # TODO: should be able to pass in an optional function for 
            # running each TU through TNRS, etc.
            
        # create this node
        node_type = 'cdao:TerminalNode' if clade.is_terminal() else 'cdao:AncestralNode'
        statements += [
                       (nUri(clade.uri), qUri('rdf:type'), qUri(node_type)),
                       ]
                      
        if not parent is None:
            # create edge from the parent node to this node
            self.edge_counter += 1
            edge_uri = 'edge%s' % self.edge_counter
            statements += [
                           (nUri(edge_uri), qUri('rdf:type'), qUri('cdao:Directed_Edge')),
                           (nUri(edge_uri), qUri('cdao:has_Parent_Node'), nUri(parent.uri)),
                           (nUri(edge_uri), qUri('cdao:has_Child_Node'), nUri(clade.uri)),
                           (nUri(clade.uri), qUri('cdao:belongs_to_Edge_as_Child'), nUri(edge_uri)),
                           (nUri(clade.uri), qUri('cdao:has_Parent'), nUri(parent.uri)),
                           (nUri(parent.uri), qUri('cdao:belongs_to_Edge_as_Parent'), nUri(edge_uri)),
                           ]
            # add branch length
            edge_ann_uri = 'edge_annotation%s' % self.edge_counter
            statements += [
                           (nUri(edge_ann_uri), qUri('rdf:type'), qUri('cdao:EdgeLength')),
                           (nUri(edge_uri), qUri('cdao:has_annotation'), nUri(edge_ann_uri)),
                           # TODO: does this type of numeric literal actually work?
                           (nUri(edge_ann_uri), qUri('cdao:has_value'), str(clade.branch_length)),
                           ]
            # TODO: annotate with confidences?
                      
        for stmt in statements:
            yield RDF.Statement(*stmt)
        
        if not clade.is_terminal():
            for new_clade in clade.clades:
                for stmt in self.process_clade(new_clade, parent=clade, root=False):
                    yield stmt
                    
                    
def qUri(s):
    '''returns the full URI from a namespaced URI string (i.e. rdf:type)'''
    RDF = import_rdf()
    
    for url in RDF_NAMESPACES: 
        s = s.replace(url+':', RDF_NAMESPACES[url])
    return RDF.Uri(s)
    
def namedUri(s, base_uri):
    '''append a URI to the base URI'''
    RDF = import_rdf()
    
    return RDF.Uri(base_uri + s)
