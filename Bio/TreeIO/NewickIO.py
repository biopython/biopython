# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrappers for the Newick file format.

See: U{ http://evolution.genetics.washington.edu/phylip/newick_doc.html }
"""
__docformat__ = "epytext en"

from cStringIO import StringIO

from Bio.Nexus import Trees


def parse(file):
    buf = ''
    for line in file:
        buf += line.rstrip()
        if buf.endswith(';'):
            yield Trees.Tree(tree=buf)
            buf = ''
    if buf:
        # Last tree is missing a terminal ';' character
        yield Trees.Tree(tree=buf)


def write(trees, file, plain=False, **kwargs):
    count = 0
    for tree in trees:
        file.write(tree.to_string(plain_newick=True, plain=plain, **kwargs)
                    + ';\n')
        count += 1
    return count

# Copying Bio.Nexus.Trees parser

class Parser(object):
    """Parse a Newick tree given a file handle."""

    def __init__(self, handle,
            # from Nexus.Trees
            tree=None, weight=1.0, rooted=False, name='', #data=NodeData,
            values_are_support=False, max_support=1.0
            ):
        Nodes.Chain.__init__(self)
        self.dataclass=data
        self.__values_are_support=values_are_support
        self.max_support=max_support
        self.weight=weight
        self.rooted=rooted
        self.name=name
        root=Nodes.Node(data())
        self.root = self.add(root)
        if tree:    # use the tree we have
            # if Tree is called from outside Nexus parser, we need to get rid of linebreaks, etc
            tree=tree.strip().replace('\n','').replace('\r','')
            # there's discrepancy whether newick allows semicolons et the end
            tree=tree.rstrip(';')
            subtree_info, base_info = self._parse(tree)
            root.data = self._add_nodedata(root.data, [[], base_info])
            self._add_subtree(parent_id=root.id,tree=subtree_info)

    @classmethod
    def from_string(cls, treetext):
        handle = StringIO(treetext)
        return cls(handle)

    def parse(self):
        """Parse the text stream this object was initialized with."""

    # --- XXX copied -----------------------------

    def _parse(self,tree):
        """Parses (a,b,c...)[[[xx]:]yy] into subcomponents and travels down recursively."""
        #Remove any leading/trailing white space - want any string starting
        #with " (..." should be recognised as a leaf, "(..."
        tree = tree.strip()
        if tree.count('(')!=tree.count(')'):
            raise TreeError('Parentheses do not match in (sub)tree: '+tree)
        if tree.count('(')==0: # a leaf
            #check if there's a colon, or a special comment, or both  after the taxon name
            nodecomment=tree.find(NODECOMMENT_START)
            colon=tree.find(':')
            if colon==-1 and nodecomment==-1: # none
                return [tree,[None]]
            elif colon==-1 and nodecomment>-1: # only special comment
                return [tree[:nodecomment],self._get_values(tree[nodecomment:])]
            elif colon>-1 and nodecomment==-1: # only numerical values
                return [tree[:colon],self._get_values(tree[colon+1:])]
            elif colon < nodecomment: # taxon name ends at first colon or with special comment
                return [tree[:colon],self._get_values(tree[colon+1:])]
            else:
                return [tree[:nodecomment],self._get_values(tree[nodecomment:])]
        else:
            closing=tree.rfind(')')
            val=self._get_values(tree[closing+1:])
            if not val:
                val=[None]
            subtrees=[]
            plevel=0
            prev=1
            for p in range(1,closing):
                if tree[p]=='(':
                    plevel+=1
                elif tree[p]==')':
                    plevel-=1
                elif tree[p]==',' and plevel==0:
                    subtrees.append(tree[prev:p])
                    prev=p+1
            subtrees.append(tree[prev:closing])
            subclades=[self._parse(subtree) for subtree in subtrees]
            return [subclades,val]

    def _add_subtree(self,parent_id=None,tree=None):
        """Adds leaf or tree (in newick format) to a parent_id. (self,parent_id,tree)."""
        if parent_id is None:
            raise TreeError('Need node_id to connect to.')
        for st in tree:
            nd=self.dataclass()
            nd = self._add_nodedata(nd, st)
            if type(st[0])==list: # it's a subtree
                sn=Nodes.Node(nd)
                self.add(sn,parent_id)
                self._add_subtree(sn.id,st[0])
            else: # it's a leaf
                nd.taxon=st[0]
                leaf=Nodes.Node(nd)
                self.add(leaf,parent_id)

    def _add_nodedata(self, nd, st):
        """Add data to the node parsed from the comments, taxon and support.
        """
        if isinstance(st[1][-1],str) and st[1][-1].startswith(NODECOMMENT_START):
            nd.comment=st[1].pop(-1)
        # if the first element is a string, it's the subtree node taxon
        elif isinstance(st[1][0], str):
            nd.taxon = st[1][0]
            st[1] = st[1][1:]
        if len(st)>1:
            if len(st[1])>=2: # if there's two values, support comes first. Is that always so?
                nd.support=st[1][0]
                if st[1][1] is not None:
                    nd.branchlength=st[1][1]
            elif len(st[1])==1: # otherwise it could be real branchlengths or support as branchlengths
                if not self.__values_are_support: # default
                    if st[1][0] is not None:
                        nd.branchlength=st[1][0]
                else:
                    nd.support=st[1][0]
        return nd

    def _get_values(self, text):
        """Extracts values (support/branchlength) from xx[:yyy], xx."""
       
        if text=='':
            return None
        nodecomment = None
        if NODECOMMENT_START in text: # if there's a [&....] comment, cut it out
            nc_start=text.find(NODECOMMENT_START)
            nc_end=text.find(NODECOMMENT_END)
            if nc_end==-1:
                raise TreeError('Error in tree description: Found %s without matching %s' \
                                % (NODECOMMENT_START, NODECOMMENT_END))
            nodecomment=text[nc_start:nc_end+1]
            text=text[:nc_start]+text[nc_end+1:]

        # pase out supports and branchlengths, with internal node taxa info
        values = []
        taxonomy = None
        for part in [t.strip() for t in text.split(":")]:
            if part:
                try:
                    values.append(float(part))
                except ValueError:
                    assert taxonomy is None, "Two string taxonomies?"
                    taxonomy = part
        if taxonomy:
            values.insert(0, taxonomy)
        if nodecomment:
            values.append(nodecomment)
        return values


# Copying Bio.Nexus.Trees writer (str, to_string)

class Writer(object):
    """
    """

    def __init__(self, trees):
        pass

    def write(self, handle):
        """Parse the text stream this object was initialized with."""

    def to_string(self):
        stream = StringIO()
        self.write(stream)
        return stream.read()

    # --- XXX copied -----------------------------

    def to_string(self,support_as_branchlengths=False,branchlengths_only=False,plain=True,plain_newick=False,ladderize=None):
        """Return a paup compatible tree line.
       
        to_string(self,support_as_branchlengths=False,branchlengths_only=False,plain=True)
        """
        # if there's a conflict in the arguments, we override plain=True
        if support_as_branchlengths or branchlengths_only:
            plain=False
        self.support_as_branchlengths=support_as_branchlengths
        self.branchlengths_only=branchlengths_only
        self.plain=plain

        def make_info_string(data,terminal=False):
            """Creates nicely formatted support/branchlengths."""
            # CHECK FORMATTING
            if self.plain: # plain tree only. That's easy.
                return ''
            elif self.support_as_branchlengths: # support as branchlengths (eg. PAUP), ignore actual branchlengths
                if terminal:    # terminal branches have 100% support
                    return ':%1.2f' % self.max_support
                else:
                    return ':%1.2f' % (data.support)
            elif self.branchlengths_only: # write only branchlengths, ignore support
                return ':%1.5f' % (data.branchlength)
            else:   # write suport and branchlengths (e.g. .con tree of mrbayes)
                if terminal:
                    return ':%1.5f' % (data.branchlength)
                else:
                    if data.branchlength is not None and data.support is not None:  # we have blen and suppport
                        return '%1.2f:%1.5f' % (data.support,data.branchlength)
                    elif data.branchlength is not None:                             # we have only blen
                        return '0.00000:%1.5f' % (data.branchlength)
                    elif data.support is not None:                                  # we have only support
                        return '%1.2f:0.00000' % (data.support)
                    else:
                        return '0.00:0.00000'
        def ladderize_nodes(nodes,ladderize=None):
            """Sorts node numbers according to the number of terminal nodes."""
            if ladderize in ['left','LEFT','right','RIGHT']:
                succnode_terminals=[(self.count_terminals(node=n),n) for n in nodes]
                succnode_terminals.sort()
                if (ladderize=='right' or ladderize=='RIGHT'):
                    succnode_terminals.reverse()
                if succnode_terminals:
                    succnodes=zip(*succnode_terminals)[1]
                else:
                    succnodes=[]
            else:
                succnodes=nodes
            return succnodes

        def newickize(node,ladderize=None):
            """Convert a node tree to a newick tree recursively."""

            if not self.node(node).succ:    #terminal
                return self.node(node).data.taxon+make_info_string(self.node(node).data,terminal=True)
            else:
                succnodes=ladderize_nodes(self.node(node).succ,ladderize=ladderize)
                subtrees=[newickize(sn,ladderize=ladderize) for sn in succnodes]
                return '(%s)%s' % (','.join(subtrees),make_info_string(self.node(node).data))
                    
        treeline=['tree']
        if self.name:
            treeline.append(self.name)
        else:
            treeline.append('a_tree')
        treeline.append('=')
        if self.weight != 1:
            treeline.append('[&W%s]' % str(round(float(self.weight),3)))
        if self.rooted:
            treeline.append('[&R]')
        succnodes=ladderize_nodes(self.node(self.root).succ)
        subtrees=[newickize(sn,ladderize=ladderize) for sn in succnodes]
        treeline.append('(%s)' % ','.join(subtrees))
        if plain_newick:
            return treeline[-1]
        else:
            return ' '.join(treeline)+';'
        
    def __str__(self):
        """Short version of to_string(), gives plain tree"""
        return self.to_string(plain=True)

