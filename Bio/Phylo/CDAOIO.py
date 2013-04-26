# Copyright (C) 2013 by Ben Morris (ben@bendmorris.com)
# Based on Bio.Nexus, copyright 2005-2008 by Frank Kauff & Cymon J. Cox
# and Bio.Phylo.Newick, copyright 2009 by Eric Talevich.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrappers for the RDF/CDAO file format.

This is an RDF format that conforms to the Comparative Data Analysis Ontology (CDAO).
See: http://www.evolutionaryontology.org/cdao

This module requires the librdf Python bindings (http://www.librdf.org)

The CDAOIO.Parser, in addition to parsing text files, can also parse directly
from a triple store that implements the Redland storage interface; similarly,
the CDAOIO.Writer can store triples in a triple store instead of serializing
them to a file.
"""

__docformat__ = "restructuredtext en"

from cStringIO import StringIO

from Bio.Phylo import CDAO
from _cdao_owl import cdao_elements, cdao_namespaces, resolve_uri
import os
import urlparse

class CDAOError(Exception):
    """Exception raised when CDAO object construction cannot continue."""
    pass

try: 
    import RDF
    import Redland
except ImportError:
    raise CDAOError('Support for CDAO tree format requires the librdf Python bindings.')

RDF_NAMESPACES = {
                  'owl':  'http://www.w3.org/2002/07/owl#',
                  'rdf':  'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
                  'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
                  }
RDF_NAMESPACES.update(cdao_namespaces)

def node_uri(graph, uri):
    '''Returns the full URI of a node by appending the node URI to the graph URI.'''
    if graph.endswith('/'):
        return urlparse.urljoin(graph, uri)
    elif graph:
        return urlparse.urljoin(graph, '#%s' % uri)
    else: return uri


def new_storage():
    '''Create a new in-memory Redland store for storing the RDF model.'''

    storage = RDF.Storage(storage_name="hashes",
                          name="test",
                          options_string="new='yes',hash-type='memory',dir='.'")
    if storage is None:
        raise CDAOError("Creation of new RDF.Storage failed.")
    return storage


# ---------------------------------------------------------
# Public API

def parse(handle, **kwargs):
    """Iterate over the trees in a CDAO file handle.

    :returns: generator of Bio.Phylo.CDAO.Tree objects.
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
    def __init__(self, handle=None):
        self.handle = handle
        self.model = None
        self.node_info = None
        self.children = {}
        self.rooted = False

    @classmethod
    def from_string(cls, treetext):
        handle = StringIO(treetext)
        return cls(handle)

    def parse(self, **kwargs):
        """Parse the text stream this object was initialized with."""
        self.parse_handle_to_model(**kwargs)
        return self.parse_model()
        
    def parse_handle_to_model(self, rooted=False, storage=None, 
                              parse_format='turtle', context=None, **kwargs):
        '''Parse self.handle into RDF model self.model.'''

        if storage is None:
            # store RDF model in memory for now
            storage = new_storage()

        if self.model is None:
            self.model = RDF.Model(storage)
            if self.model is None:
                raise CDAOError("new RDF.model failed")
        model = self.model
        
        self.rooted = rooted
        
        parser = RDF.Parser(name=parse_format)
        if parser is None:
            raise Exception('Failed to create RDF.Parser for MIME type %s' % mime_type)
        
        if 'base_uri' in kwargs: base_uri = kwargs['base_uri']
        else: base_uri = RDF.Uri(string="file://"+os.path.abspath(self.handle.name))
        
        statements = parser.parse_string_as_stream(self.handle.read(), base_uri)
        for s in statements:
            model.append(s)
            
        return self.parse_model(model, context=context)
            
            
    def parse_model(self, model=None, context=None):
        '''Generator that yields CDAO.Tree instances from an RDF model.'''
        
        if model is None:
            model = self.model

        if not context is None: context = RDF.Node(RDF.Uri(context))
        
        # look up branch lengths/TUs for all nodes
        self.get_node_info(model, context=context)
        
        for root_node in self.tree_roots:
            clade = self.parse_children(root_node)
            
            yield CDAO.Tree(root=clade, rooted=self.rooted)
            
            
    def new_clade(self, node):
        '''Returns a CDAO.Clade object for a given named node.'''
        
        result = self.node_info[node]
        
        kwargs = {}
        if 'branch_length' in result: kwargs['branch_length'] = result['branch_length']
        if 'label' in result: kwargs['name'] = result['label'].replace('_', ' ')
        if 'confidence' in result: kwargs['confidence'] = result['confidence']
        
        clade = CDAO.Clade(**kwargs)
        
        return clade
        
            
    def get_node_info(self, model, context=None):
        '''Creates a dictionary containing information about all nodes in the tree.'''
        
        Uri = RDF.Uri
        
        self.node_info = {}
        self.obj_info = {}
        self.children = {}
        self.nodes = set()
        self.tree_roots = set()
        
        for statement, context in model.as_stream_context(context):
            # process each RDF triple in the model sequentially
            s, v, o = str(statement.subject), Uri(str(statement.predicate)), str(statement.object)
            
            if not s in self.obj_info: self.obj_info[s] = {}
            this = self.obj_info[s]
            
            assignments = {
                           qUri('cdao:has_Parent'): 'parent',
                           qUri('cdao:belongs_to_Edge_as_Child'): 'edge',
                           qUri('cdao:has_Annotation'): 'annotation',
                           qUri('cdao:has_Value'): 'value',
                           qUri('cdao:represents_TU'): 'tu',
                           qUri('rdfs:label'): 'label',
                           qUri('cdao:has_Support_Value'): 'confidence',
                           }
            
            try:
                # if the predicate is one we care about, store information for later
                this[assignments[v]] = o
            except KeyError: pass
            
            if v == qUri('rdf:type'):
                if Uri(o) in (qUri('cdao:AncestralNode'), qUri('cdao:TerminalNode')):
                    # this is a tree node; store it in set of all nodes
                    self.nodes.add(s)
            if v == qUri('cdao:has_Root'):
                # this is a tree; store its root in set of all tree roots
                self.tree_roots.add(o)
                    
        for node in self.nodes:
            # for each node, look up all information needed to create a CDAO.Clade
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
        '''Return a CDAO.Clade, and calls itself recursively for each child, 
        traversing the  entire tree and creating a nested structure of CDAO.Clade 
        objects.'''
        
        clade = self.new_clade(node)
        
        children = self.children[node] if node in self.children else []
        clade.clades = [self.parse_children(child_node) for child_node in children]
        
        return clade


# ---------------------------------------------------------
# Output

class Writer(object):
    """Based on the writer in Bio.Nexus.Trees (str, to_string)."""
    prefixes = RDF_NAMESPACES

    def __init__(self, trees):
        self.trees = trees
        
        self.node_counter = 0
        self.edge_counter = 0
        self.tu_counter = 0
        self.tree_counter = 0

    def write(self, handle, tree_uri='', context=None, record_complete_ancestry=False, 
              rooted=False, **kwargs):
        """Write this instance's trees to a file handle.

        If the handle is a librdf model, statements will be added directly to it.
        """

        self.rooted = rooted
        self.record_complete_ancestry = record_complete_ancestry
        
        self.add_trees_to_handle(handle, tree_uri=tree_uri, context=context)
        
        
    def add_trees_to_handle(self, handle, trees=None, tree_uri='', context=None):
        """Add triples describing a set of trees to handle, which can be either 
        a file or a librdf model."""

        if tree_uri and not tree_uri.endswith('/'): tree_uri = tree_uri + '/'

        is_librdf_model = isinstance(handle, RDF.Model)

        if is_librdf_model and context: context = RDF.Node(RDF.Uri(context))
        else: context = None

        Uri = RDF.Uri
        prefixes = self.prefixes
        
        if trees is None:
            trees = self.trees
        
        if is_librdf_model: 
            Redland.librdf_model_transaction_start(handle._model)
        else:
            for prefix, url in prefixes.items():
                handle.write('@prefix %s: <%s> .\n' % (prefix, url))
        
        for stmt in [(Uri(prefixes['cdao']), qUri('rdf:type'), qUri('owl:Ontology'))]:
            self.add_stmt_to_handle(handle, RDF.Statement(*stmt), context)
        
        for tree in trees:
            self.tree_counter += 1
            self.tree_uri = node_uri(tree_uri, 'tree%s' % str(self.tree_counter).zfill(7))

            first_clade = tree.clade
            statements = self.process_clade(first_clade, root=tree_uri)
            for stmt in statements:
                self.add_stmt_to_handle(handle, stmt, context)
                
        if is_librdf_model: 
            Redland.librdf_model_transaction_commit(handle._model)
            handle.sync()


    def add_stmt_to_handle(self, handle, stmt, context):
        if isinstance(handle, RDF.Model): handle.append(stmt, context)
        else: 
            # apply URI prefixes
            stmt_parts = [stmt.subject, stmt.predicate, stmt.object]
            stmt_strings = []
            for n, part in enumerate(stmt_parts):
                if part.is_resource():
                    changed = False
                    node_uri = str(part)
                    if n > 0:
                        for prefix, uri in self.prefixes.items():
                            if node_uri.startswith(uri):
                                node_uri = node_uri.replace(uri, '%s:'%prefix, 1)
                                if node_uri == 'rdf:type': node_uri = 'a'
                                changed = True
                        if changed: stmt_strings.append(node_uri)
                        else: stmt_strings.append('<%s>' % node_uri)
                    else: stmt_strings.append('<%s>' % node_uri)

                elif part.is_literal():
                    stmt_strings.append(('"%s"' % str(part.literal[0])) + 
                                        (('^^<%s>' % str(part.literal[2])) if part.literal[2] else ''))
                    
                else:
                    stmt_strings.append(str(part))
            
            handle.write('%s .\n' % ' '.join(stmt_strings))
        
            
    def process_clade(self, clade, parent=None, root=False):
        '''recursively generate statements describing a tree of clades'''
        
        self.node_counter += 1
        clade.uri = 'node%s' % str(self.node_counter).zfill(7)
        if parent: clade.ancestors = parent.ancestors + [parent.uri]
        else: clade.ancestors = []
        
        Uri = RDF.Uri
        nUri = lambda s: Uri(node_uri(self.tree_uri, s))
        tree_id = nUri('')
        
        statements = []
        
        if not root is False:
            # create a cdao:RootedTree with reference to the tree root
            tree_type = qUri('cdao:RootedTree') if self.rooted else qUri('cdao:UnrootedTree')
            
            statements += [
                           (tree_id, qUri('rdf:type'), tree_type),
                           (tree_id, qUri('cdao:has_Root'), nUri(clade.uri)),
                           ]
        
        if clade.name:
            # create TU
            self.tu_counter += 1
            tu_uri = 'tu%s' % str(self.tu_counter).zfill(7)

            statements += [
                           (nUri(tu_uri), qUri('rdf:type'), qUri('cdao:TU')),
                           (nUri(clade.uri), qUri('cdao:represents_TU'), nUri(tu_uri)),
                           (nUri(tu_uri), qUri('rdfs:label'), RDF.Node(literal=clade.name.replace('_', ' '))),
                           ]
                           
            # TODO: should be able to pass in an optional function for 
            # running each TU through TNRS, etc.
            
        # create this node
        node_type = 'cdao:TerminalNode' if clade.is_terminal() else 'cdao:AncestralNode'
        statements += [
                       (nUri(clade.uri), qUri('rdf:type'), qUri(node_type)),
                       (nUri(clade.uri), qUri('cdao:belongs_to_Tree'), tree_id),
                       ]
                      
        if not parent is None:
            # create edge from the parent node to this node
            self.edge_counter += 1
            edge_uri = 'edge%s' % str(self.edge_counter).zfill(7)

            statements += [
                           (nUri(edge_uri), qUri('rdf:type'), qUri('cdao:DirectedEdge')),
                           (nUri(edge_uri), qUri('cdao:belongs_to_Tree'), tree_id),
                           (nUri(edge_uri), qUri('cdao:has_Parent_Node'), nUri(parent.uri)),
                           (nUri(edge_uri), qUri('cdao:has_Child_Node'), nUri(clade.uri)),
                           (nUri(clade.uri), qUri('cdao:belongs_to_Edge_as_Child'), nUri(edge_uri)),
                           (nUri(clade.uri), qUri('cdao:has_Parent'), nUri(parent.uri)),
                           (nUri(parent.uri), qUri('cdao:belongs_to_Edge_as_Parent'), nUri(edge_uri)),
                           ]

            if hasattr(clade, 'confidence') and not clade.confidence is None:
                confidence = RDF.Node(literal=str(clade.confidence), 
                                      datatype=RDF.Uri('http://www.w3.org/2001/XMLSchema#decimal'))

                statements += [(nUri(clade.uri), qUri('cdao:has_Support_Value'), confidence)]

            
            if self.record_complete_ancestry and len(clade.ancestors) > 0:
                statements += [(nUri(clade.uri), qUri('cdao:has_Ancestor'), nUri(ancestor))
                               for ancestor in clade.ancestors]
            
            # add branch length
            edge_ann_uri = 'edge_annotation%s' % str(self.edge_counter).zfill(7)

            branch_length = RDF.Node(literal=str(clade.branch_length), 
                                     datatype=RDF.Uri('http://www.w3.org/2001/XMLSchema#decimal'))
            statements += [
                           (nUri(edge_ann_uri), qUri('rdf:type'), qUri('cdao:EdgeLength')),
                           (nUri(edge_uri), qUri('cdao:has_Annotation'), nUri(edge_ann_uri)),
                           (nUri(edge_ann_uri), qUri('cdao:has_Value'), branch_length),
                           ]
                      
        for stmt in statements:
            yield RDF.Statement(*stmt)
        
        if not clade.is_terminal():
            for new_clade in clade.clades:
                for stmt in self.process_clade(new_clade, parent=clade, root=False):
                    yield stmt
                    

def qUri(s):
    '''returns the full URI from a namespaced URI string (i.e. rdf:type)'''

    return RDF.Uri(resolve_uri(s, namespaces=RDF_NAMESPACES))