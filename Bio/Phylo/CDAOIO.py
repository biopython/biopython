"""I/O function wrappers for the RDF/CDAO file format.

"""
__docformat__ = "restructuredtext en"

from cStringIO import StringIO

from Bio.Phylo import Newick
import os

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

'''    @classmethod
    def from_string(cls, treetext):
        handle = StringIO(treetext)
        return cls(handle)

    def parse(self, values_are_confidence=False, rooted=False):
        """Parse the text stream this object was initialized with."""
        self.values_are_confidence = values_are_confidence
        self.rooted = rooted    # XXX this attribue is useless
        buf = ''
        for line in self.handle:
            buf += line.rstrip()
            if buf.endswith(';'):
                yield self._parse_tree(buf, rooted)
                buf = ''
        if buf:
            # Last tree is missing a terminal ';' character -- that's OK
            yield self._parse_tree(buf, rooted)

    def _parse_tree(self, text, rooted):
        """Parses the text representation into an Tree object."""
        # XXX Pass **kwargs along from Parser.parse?
        return Newick.Tree(root=self._parse_subtree(text), rooted=self.rooted)

    def _parse_subtree(self, text):
        """Parse ``(a,b,c...)[[[xx]:]yy]`` into subcomponents, recursively."""
        text = text.strip().rstrip(';')
        if text.count('(')!=text.count(')'):
            raise NewickError("Parentheses do not match in (sub)tree: " + text)
        # Text is now "(...)..." (balanced parens) or "..." (leaf node)
        if text.count('(') == 0:
            # Leaf/terminal node -- recursion stops here
            return self._parse_tag(text)
        # Handle one layer of the nested subtree
        # XXX what if there's a paren in a comment or other string?
        close_posn = text.rfind(')')
        subtrees = []
        # Locate subtrees by counting nesting levels of parens
        plevel = 0
        prev = 1
        for posn in range(1, close_posn):
            if text[posn] == '(':
                plevel += 1
            elif text[posn] == ')':
                plevel -= 1
            elif text[posn] == ',' and plevel == 0:
                subtrees.append(text[prev:posn])
                prev = posn + 1
        subtrees.append(text[prev:close_posn])
        # Construct a new clade from trailing text, then attach subclades
        clade = self._parse_tag(text[close_posn+1:])
        clade.clades = [self._parse_subtree(st) for st in subtrees]
        return clade

    def _parse_tag(self, text):
        """Extract the data for a node from text.

        :returns: Clade instance containing any available data
        """
        # Extract the comment
        comment_start = text.find(NODECOMMENT_START)
        if comment_start != -1:
            comment_end = text.find(NODECOMMENT_END)
            if comment_end == -1:
                raise NewickError('Error in tree description: '
                                  'Found %s without matching %s'
                                  % (NODECOMMENT_START, NODECOMMENT_END))
            comment = text[comment_start+len(NODECOMMENT_START):comment_end]
            text = text[:comment_start] + text[comment_end+len(NODECOMMENT_END):]
        else:
            comment = None
        clade = Newick.Clade(comment=comment)
        # Extract name (taxon), and optionally support, branch length
        # Float values are support and branch length, the string is name/taxon
        values = []
        for part in (t.strip() for t in text.split(':')):
            if part:
                try:
                    values.append(float(part))
                except ValueError:
                    assert clade.name is None, "Two string taxonomies?"
                    clade.name = part
        if len(values) == 1:
            # Real branch length, or support as branch length
            if self.values_are_confidence:
                clade.confidence = values[0]
            else:
                clade.branch_length = values[0]
        elif len(values) == 2:
            # Two non-taxon values: support comes first. (Is that always so?)
            clade.confidence, clade.branch_length = values
        elif len(values) > 2:
            raise NewickError("Too many colons in tag: " + text)
        return clade'''


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

    def write(self, handle, **kwargs):
        """Write this instance's trees to a file handle.
        
        Keywords:
            mime_type: used to determine the serialization format.
                default is 'text/turtle'
        """
        try: mime_type = kwargs['mime_type']
        except KeyError: mime_type = 'text/turtle'
        
        model = self.add_trees_to_model()
        return self.serialize_model(model, mime_type=mime_type)
        
    def add_trees_to_model(self, storage=None):
        import RDF
        Uri = RDF.Uri
        urls = self.urls
        
        def qUri(s):
            if ':' in s:
                s = s.split(':')
                s = urls[s[0]] + ':'.join(s[1:])
            return Uri(s)
        
        if storage is None:
            # store RDF model in memory for now
            storage = RDF.Storage(storage_name="hashes",
                                  name="serializer",
                                  options_string="new='yes',hash-type='memory',dir='.'")
            if storage is None:
                raise CDAOError("new RDF.Storage failed")

        model = RDF.Model(storage)
        if model is None:
            raise CDAOError("new RDF.model failed")
        
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

        self.count = 0
        for tree in self.trees:
            first_clade = tree.clade
            process_clade(first_clade)
            
            self.count += 1
            
    def serialize_model(self, model, mime_type='text_turtle'):
        import RDF
        
        # serialize RDF model to output file
        serializer = RDF.Serializer(mime_type=mime_type)
        for prefix, url in urls.items():
            serializer.set_namespace(prefix, url)

        # TODO: this is going to be too memory intensive for large trees;
        # come up with something better
        print "Writing..."
        handle.write(serializer.serialize_model_to_string(model))
        print "Done."
        
        return self.count
