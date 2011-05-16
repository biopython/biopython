# Copyright 2005-2008 by Frank Kauff & Cymon J. Cox. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.
#
# Bug reports welcome: fkauff@biologie.uni-kl.de or on Biopython's bugzilla.
"""Tree class to handle phylogenetic trees.

Provides a set of methods to read and write newick-format tree descriptions,
get information about trees (monphyly of taxon sets, congruence between trees,
common ancestors,...) and to manipulate trees (reroot trees, split terminal
nodes).
"""

import sys, random, copy
import Nodes

PRECISION_BRANCHLENGTH=6
PRECISION_SUPPORT=6
NODECOMMENT_START='[&'
NODECOMMENT_END=']'

class TreeError(Exception): pass

class NodeData(object):
    """Stores tree-relevant data associated with nodes (e.g. branches or otus)."""
    def __init__(self,taxon=None,branchlength=0.0,support=None,comment=None):
        self.taxon=taxon
        self.branchlength=branchlength
        self.support=support
        self.comment=comment

class Tree(Nodes.Chain):
    """Represents a tree using a chain of nodes with on predecessor (=ancestor)
    and multiple successors (=subclades).
    """ 
    # A newick tree is parsed into nested list and then converted to a node list in two stages
    # mostly due to historical reasons. This could be done in one swoop). Note: parentheses ( ) and
    # colon : are not allowed in taxon names. This is against NEXUS standard, but makes life much
    # easier when parsing trees.
    
    ## NOTE: Tree should store its data class in something like self.dataclass=data,
    ## so that nodes that are generated have easy access to the data class
    ## Some routines use automatically NodeData, this needs to be more concise

    def __init__(self,tree=None,weight=1.0,rooted=False,name='',data=NodeData,values_are_support=False,max_support=1.0):
        """Ntree(self,tree)."""
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
   
    def _walk(self,node=None):
        """Return all node_ids downwards from a node."""
        
        if node is None:
            node=self.root
        for n in self.node(node).succ:
            yield n
            for sn in self._walk(n):
                yield sn

    def node(self,node_id):
        """Return the instance of node_id.
        
        node = node(self,node_id)
        """
        if node_id not in self.chain:
            raise TreeError('Unknown node_id: %d' % node_id)
        return self.chain[node_id]

    def split(self,parent_id=None,n=2,branchlength=1.0):
        """Speciation: generates n (default two) descendants of a node.
        
        [new ids] = split(self,parent_id=None,n=2,branchlength=1.0):
        """ 
        if parent_id is None:
            raise TreeError('Missing node_id.')
        ids=[]
        parent_data=self.chain[parent_id].data
        for i in range(n):
            node=Nodes.Node()
            if parent_data:
                node.data=self.dataclass()
                # each node has taxon and branchlength attribute
                if parent_data.taxon:
                    node.data.taxon=parent_data.taxon+str(i)
                node.data.branchlength=branchlength
            ids.append(self.add(node,parent_id))
        return ids

    def search_taxon(self,taxon):
        """Returns the first matching taxon in self.data.taxon. Not restricted to terminal nodes.
        
        node_id = search_taxon(self,taxon)
        """
        for id,node in self.chain.iteritems():
            if node.data.taxon==taxon:
                return id
        return None
   
    def prune(self,taxon):
        """Prunes a terminal taxon from the tree.
        
        id_of_previous_node = prune(self,taxon)
        If taxon is from a bifurcation, the connectiong node will be collapsed
        and its branchlength added to remaining terminal node. This might be no
        longer a meaningful value'
        """
        
        id=self.search_taxon(taxon)
        if id is None:
            raise TreeError('Taxon not found: %s' % taxon)
        elif id not in self.get_terminals():
            raise TreeError('Not a terminal taxon: %s' % taxon)
        else:
            prev=self.unlink(id)
            self.kill(id)
            if len(self.node(prev).succ)==1:
                if prev==self.root: # we deleted one branch of a bifurcating root, then we have to move the root upwards
                    self.root=self.node(self.root).succ[0]    
                    self.node(self.root).branchlength=0.0
                    self.kill(prev)
                else: 
                    succ=self.node(prev).succ[0]
                    new_bl=self.node(prev).data.branchlength+self.node(succ).data.branchlength
                    self.collapse(prev)
                    self.node(succ).data.branchlength=new_bl
            return prev
        
    def get_taxa(self,node_id=None):
        """Return a list of all otus downwards from a node (self, node_id).

        nodes = get_taxa(self,node_id=None)
        """

        if node_id is None:
            node_id=self.root
        if node_id not in self.chain:
            raise TreeError('Unknown node_id: %d.' % node_id)    
        if self.chain[node_id].succ==[]:
            if self.chain[node_id].data:
                return [self.chain[node_id].data.taxon]
            else:
                return None
        else:
            list=[]
            for succ in self.chain[node_id].succ:
                list.extend(self.get_taxa(succ))
            return list

    def get_terminals(self):
        """Return a list of all terminal nodes."""
        return [i for i in self.all_ids() if self.node(i).succ==[]]
    
    def is_terminal(self,node):
        """Returns True if node is a terminal node."""
        return self.node(node).succ==[]

    def is_internal(self,node):
        """Returns True if node is an internal node."""
        return len(self.node(node).succ)>0

    def is_preterminal(self,node):
        """Returns True if all successors of a node are terminal ones."""
        if self.is_terminal(node):
            return False not in [self.is_terminal(n) for n in self.node(node).succ]
        else:
            return False
    def count_terminals(self,node=None):
        """Counts the number of terminal nodes that are attached to a node."""
        if node is None:
            node=self.root
        return len([n for n in self._walk(node) if self.is_terminal(n)])

    def collapse_genera(self,space_equals_underscore=True):
        """Collapses all subtrees which belong to the same genus (i.e share the same first word in their taxon name."""

        while True:
            for n in self._walk():
                if self.is_terminal(n):
                    continue
                taxa=self.get_taxa(n)
                genera=[]
                for t in taxa:
                    if space_equals_underscore:
                        t=t.replace(' ','_')
                    try:
                        genus=t.split('_',1)[0]
                    except:
                        genus='None'
                    if genus not in genera:
                        genera.append(genus)
                if len(genera)==1:
                    self.node(n).data.taxon=genera[0]+' <collapsed>'
                    #now we kill all nodes downstream
                    nodes2kill=[kn for kn in self._walk(node=n)]
                    for kn in nodes2kill:
                        self.kill(kn)
                    self.node(n).succ=[]
                    break # break out of for loop because node list from _walk will be inconsistent
            else: # for loop exhausted: no genera to collapse left
                break # while


    def sum_branchlength(self,root=None,node=None):
        """Adds up the branchlengths from root (default self.root) to node.
        
        sum = sum_branchlength(self,root=None,node=None)
        """

        if root is None:
            root=self.root
        if node is None:
            raise TreeError('Missing node id.')
        blen=0.0
        while node is not None and node is not root: 
            blen+=self.node(node).data.branchlength
            node=self.node(node).prev
        return blen

    def set_subtree(self,node):
        """Return subtree as a set of nested sets.

        sets = set_subtree(self,node)
        """
        
        if self.node(node).succ==[]:
            return self.node(node).data.taxon
        else:
            try:
                return frozenset([self.set_subtree(n) for n in self.node(node).succ])
            except:
                print node
                print self.node(node).succ
                for n in self.node(node).succ:
                    print n, self.set_subtree(n)
                print [self.set_subtree(n) for n in self.node(node).succ]
                raise
            
    def is_identical(self,tree2):
        """Compare tree and tree2 for identity.

        result = is_identical(self,tree2)
        """
        return self.set_subtree(self.root)==tree2.set_subtree(tree2.root)

    def is_compatible(self,tree2,threshold,strict=True):
        """Compares branches with support>threshold for compatibility.
        
        result = is_compatible(self,tree2,threshold)
        """

        # check if both trees have the same set of taxa. strict=True enforces this.
        missing2=set(self.get_taxa())-set(tree2.get_taxa())
        missing1=set(tree2.get_taxa())-set(self.get_taxa())
        if strict and (missing1 or missing2):
            if missing1: 
                print 'Taxon/taxa %s is/are missing in tree %s' % (','.join(missing1) , self.name)
            if missing2:
                print 'Taxon/taxa %s is/are missing in tree %s' % (','.join(missing2) , tree2.name)
            raise TreeError('Can\'t compare trees with different taxon compositions.')
        t1=[(set(self.get_taxa(n)),self.node(n).data.support) for n in self.all_ids() if \
            self.node(n).succ and\
            (self.node(n).data and self.node(n).data.support and self.node(n).data.support>=threshold)]
        t2=[(set(tree2.get_taxa(n)),tree2.node(n).data.support) for n in tree2.all_ids() if \
            tree2.node(n).succ and\
            (tree2.node(n).data and tree2.node(n).data.support and tree2.node(n).data.support>=threshold)]
        conflict=[]
        for (st1,sup1) in t1:
            for (st2,sup2) in t2:
                if not st1.issubset(st2) and not st2.issubset(st1):                     # don't hiccup on upstream nodes
                    intersect,notin1,notin2=st1 & st2, st2-st1, st1-st2                 # all three are non-empty sets
                    # if notin1==missing1 or notin2==missing2  <==> st1.issubset(st2) or st2.issubset(st1) ???
                    if intersect and not (notin1.issubset(missing1) or notin2.issubset(missing2)):         # omit conflicts due to missing taxa
                        conflict.append((st1,sup1,st2,sup2,intersect,notin1,notin2))
        return conflict
        
    def common_ancestor(self,node1,node2):
        """Return the common ancestor that connects two nodes.
        
        node_id = common_ancestor(self,node1,node2)
        """
        
        l1=[self.root]+self.trace(self.root,node1)
        l2=[self.root]+self.trace(self.root,node2)
        return [n for n in l1 if n in l2][-1]


    def distance(self,node1,node2):
        """Add and return the sum of the branchlengths between two nodes.
        dist = distance(self,node1,node2)
        """
        
        ca=self.common_ancestor(node1,node2)
        return self.sum_branchlength(ca,node1)+self.sum_branchlength(ca,node2)

    def is_monophyletic(self,taxon_list):
        """Return node_id of common ancestor if taxon_list is monophyletic, -1 otherwise.
        
        result = is_monophyletic(self,taxon_list)
        """
        if isinstance(taxon_list,str):
            taxon_set=set([taxon_list])
        else:
            taxon_set=set(taxon_list)
        node_id=self.root
        while 1:
            subclade_taxa=set(self.get_taxa(node_id))
            if subclade_taxa==taxon_set:                                        # are we there?
                return node_id
            else:                                                               # check subnodes
                for subnode in self.chain[node_id].succ:
                    if set(self.get_taxa(subnode)).issuperset(taxon_set):  # taxon_set is downstream
                        node_id=subnode
                        break   # out of for loop
                else:
                    return -1   # taxon set was not with successors, for loop exhausted

    def is_bifurcating(self,node=None):
        """Return True if tree downstream of node is strictly bifurcating."""
        if node is None:
            node=self.root
        if node==self.root and len(self.node(node).succ)==3: #root can be trifurcating, because it has no ancestor
            return self.is_bifurcating(self.node(node).succ[0]) and \
                    self.is_bifurcating(self.node(node).succ[1]) and \
                    self.is_bifurcating(self.node(node).succ[2])
        if len(self.node(node).succ)==2:
            return self.is_bifurcating(self.node(node).succ[0]) and self.is_bifurcating(self.node(node).succ[1])
        elif len(self.node(node).succ)==0:
            return True
        else:
            return False

    def branchlength2support(self):
        """Move values stored in data.branchlength to data.support, and set branchlength to 0.0

        This is necessary when support has been stored as branchlength (e.g. paup), and has thus
        been read in as branchlength. 
        """

        for n in self.chain:
            self.node(n).data.support=self.node(n).data.branchlength
            self.node(n).data.branchlength=0.0

    def convert_absolute_support(self,nrep):
        """Convert absolute support (clade-count) to rel. frequencies.
        
        Some software (e.g. PHYLIP consense) just calculate how often clades appear, instead of
        calculating relative frequencies."""

        for n in self._walk():
            if self.node(n).data.support:
                self.node(n).data.support/=float(nrep)

    def has_support(self,node=None):
        """Returns True if any of the nodes has data.support != None."""
        for n in self._walk(node):
            if self.node(n).data.support:
                return True
        else:
            return False

    def randomize(self,ntax=None,taxon_list=None,branchlength=1.0,branchlength_sd=None,bifurcate=True):
        """Generates a random tree with ntax taxa and/or taxa from taxlabels.
    
        new_tree = randomize(self,ntax=None,taxon_list=None,branchlength=1.0,branchlength_sd=None,bifurcate=True)
        Trees are bifurcating by default. (Polytomies not yet supported).
        """

        if not ntax and taxon_list:
            ntax=len(taxon_list)
        elif not taxon_list and ntax:
            taxon_list=['taxon'+str(i+1) for i in range(ntax)]
        elif not ntax and not taxon_list:
            raise TreeError('Either numer of taxa or list of taxa must be specified.')
        elif ntax != len(taxon_list):
            raise TreeError('Length of taxon list must correspond to ntax.')
        # initiate self with empty root
        self.__init__()
        terminals=self.get_terminals()
        # bifurcate randomly at terminal nodes until ntax is reached
        while len(terminals)<ntax:
            newsplit=random.choice(terminals)
            new_terminals=self.split(parent_id=newsplit,branchlength=branchlength)
            # if desired, give some variation to the branch length
            if branchlength_sd:
                for nt in new_terminals:
                    bl=random.gauss(branchlength,branchlength_sd)
                    if bl<0:
                        bl=0
                    self.node(nt).data.branchlength=bl
            terminals.extend(new_terminals)
            terminals.remove(newsplit)
        # distribute taxon labels randomly
        random.shuffle(taxon_list)
        for (node,name) in zip(terminals,taxon_list):
            self.node(node).data.taxon=name

    def display(self):
        """Quick and dirty lists of all nodes."""
        table=[('#','taxon','prev','succ','brlen','blen (sum)','support','comment')]
        #Sort this to be consistent accross CPython, Jython, etc
        for i in sorted(self.all_ids()):
            n=self.node(i)
            if not n.data:
                table.append((str(i),'-',str(n.prev),str(n.succ),'-','-','-','-'))
            else:
                tx=n.data.taxon
                if not tx:
                    tx='-'
                blength="%0.2f" % n.data.branchlength
                if blength is None:
                    blength='-'
                    sum_blength='-'
                else:
                    sum_blength="%0.2f" % self.sum_branchlength(node=i)
                support=n.data.support
                if support is None:
                    support='-'
                else:
                    support="%0.2f" % support
                comment=n.data.comment
                if comment is None:
                    comment='-'
                table.append((str(i),tx,str(n.prev),str(n.succ),
                             blength, sum_blength, support, comment))
        print '\n'.join(['%3s %32s %15s %15s %8s %10s %8s %20s' % l for l in table])
        print '\nRoot: ',self.root

    def to_string(self,support_as_branchlengths=False,branchlengths_only=False,plain=True,plain_newick=False,ladderize=None,ignore_comments=True):
        """Return a paup compatible tree line.
       
        to_string(self,support_as_branchlengths=False,branchlengths_only=False,plain=True)
        """
        # if there's a conflict in the arguments, we override plain=True
        if support_as_branchlengths or branchlengths_only:
            plain=False
        self.support_as_branchlengths=support_as_branchlengths
        self.branchlengths_only=branchlengths_only
        self.ignore_comments=ignore_comments
        self.plain=plain

        def make_info_string(data,terminal=False):
            """Creates nicely formatted support/branchlengths."""
            # CHECK FORMATTING
            if self.plain: # plain tree only. That's easy.
                info_string= ''
            elif self.support_as_branchlengths: # support as branchlengths (eg. PAUP), ignore actual branchlengths
                if terminal:    # terminal branches have 100% support
                    info_string= ':%1.2f' % self.max_support
                elif data.support:
                    info_string= ':%1.2f' % (data.support)
                else:
                    info_string=':0.00'
            elif self.branchlengths_only: # write only branchlengths, ignore support
                info_string= ':%1.5f' % (data.branchlength)
            else:   # write suport and branchlengths (e.g. .con tree of mrbayes)
                if terminal:
                    info_string= ':%1.5f' % (data.branchlength)
                else:
                    if data.branchlength is not None and data.support is not None:  # we have blen and suppport
                        info_string= '%1.2f:%1.5f' % (data.support,data.branchlength)
                    elif data.branchlength is not None:                             # we have only blen
                        info_string= '0.00000:%1.5f' % (data.branchlength)
                    elif data.support is not None:                                  # we have only support
                        info_string= '%1.2f:0.00000' % (data.support)
                    else:
                        info_string= '0.00:0.00000'
            if not ignore_comments and hasattr(data,'nodecomment'):
                info_string=str(data.nodecomment)+info_string
            return info_string
            
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
        
    def unroot(self):
        """Defines a unrooted Tree structure, using data of a rooted Tree."""

        # travel down the rooted tree structure and save all branches and the nodes they connect

        def _get_branches(node):
            branches=[]
            for b in self.node(node).succ:
                branches.append([node,b,self.node(b).data.branchlength,self.node(b).data.support])
                branches.extend(_get_branches(b))
            return branches
    
        self.unrooted=_get_branches(self.root)
        # if root is bifurcating, then it is eliminated
        if len(self.node(self.root).succ)==2:
            # find the two branches that connect to root
            rootbranches=[b for b in self.unrooted if self.root in b[:2]]
            b1=self.unrooted.pop(self.unrooted.index(rootbranches[0]))
            b2=self.unrooted.pop(self.unrooted.index(rootbranches[1]))
            # Connect them two each other. If both have support, it should be identical (or one set to None?).
            # If both have branchlengths, they will be added
            newbranch=[b1[1],b2[1],b1[2]+b2[2]]
            if b1[3] is None:
                newbranch.append(b2[3]) # either None (both rootbranches are unsupported) or some support
            elif b2[3] is None:
                newbranch.append(b1[3]) # dito
            elif b1[3]==b2[3]:          
                newbranch.append(b1[3]) # identical support
            elif b1[3]==0 or b2[3]==0:
                newbranch.append(b1[3]+b2[3]) # one is 0, take the other
            else:
                raise TreeError('Support mismatch in bifurcating root: %f, %f' \
                                % (float(b1[3]),float(b2[3])))
            self.unrooted.append(newbranch)

    def root_with_outgroup(self,outgroup=None):
        
        def _connect_subtree(parent,child):
            """Hook subtree starting with node child to parent."""
            for i,branch in enumerate(self.unrooted):
                if parent in branch[:2] and child in branch[:2]:
                    branch=self.unrooted.pop(i)
                    break 
            else:
                raise TreeError('Unable to connect nodes for rooting: nodes %d and %d are not connected' \
                                % (parent,child))
            self.link(parent,child)
            self.node(child).data.branchlength=branch[2]
            self.node(child).data.support=branch[3]
            #now check if there are more branches connected to the child, and if so, connect them
            child_branches=[b for b in self.unrooted if child in b[:2]]
            for b in child_branches:
                if child==b[0]:
                    succ=b[1]
                else:
                    succ=b[0]
                _connect_subtree(child,succ) 
            
        # check the outgroup we're supposed to root with
        if outgroup is None:
            return self.root
        outgroup_node=self.is_monophyletic(outgroup)
        if outgroup_node==-1:
            return -1
        # if tree is already rooted with outgroup on a bifurcating root,
        # or the outgroup includes all taxa on the tree, then we're fine
        if (len(self.node(self.root).succ)==2 and outgroup_node in self.node(self.root).succ) or outgroup_node==self.root:
            return self.root
        
        self.unroot()
        # now we find the branch that connects outgroup and ingroup
        #print self.node(outgroup_node).prev
        for i,b in enumerate(self.unrooted):
            if outgroup_node in b[:2] and self.node(outgroup_node).prev in b[:2]:
                root_branch=self.unrooted.pop(i)
                break
        else:
            raise TreeError('Unrooted and rooted Tree do not match')
        if outgroup_node==root_branch[1]:
            ingroup_node=root_branch[0]
        else:
            ingroup_node=root_branch[1]
        # now we destroy the old tree structure, but keep node data. Nodes will be reconnected according to new outgroup
        for n in self.all_ids():
            self.node(n).prev=None
            self.node(n).succ=[]
        # now we just add both subtrees (outgroup and ingroup) branch for branch
        root=Nodes.Node(data=NodeData())            # new root    
        self.add(root)                              # add to tree description
        self.root=root.id                           # set as root
        self.unrooted.append([root.id,ingroup_node,root_branch[2],root_branch[3]])  # add branch to ingroup to unrooted tree
        self.unrooted.append([root.id,outgroup_node,0.0,0.0])   # add branch to outgroup to unrooted tree
        _connect_subtree(root.id,ingroup_node)      # add ingroup
        _connect_subtree(root.id,outgroup_node)     # add outgroup
        # if theres still a lonely node in self.chain, then it's the old root, and we delete it
        oldroot=[i for i in self.all_ids() if self.node(i).prev is None and i!=self.root]
        if len(oldroot)>1:
            raise TreeError('Isolated nodes in tree description: %s' \
                            % ','.join(oldroot))
        elif len(oldroot)==1:
            self.kill(oldroot[0])
        return self.root
        
    def merge_with_support(self,bstrees=None,constree=None,threshold=0.5,outgroup=None):
        """Merges clade support (from consensus or list of bootstrap-trees) with phylogeny.

        tree=merge_bootstrap(phylo,bs_tree=<list_of_trees>)
        or
        tree=merge_bootstrap(phylo,consree=consensus_tree with clade support)
        """

        if bstrees and constree:
            raise TreeError('Specify either list of boostrap trees or consensus tree, not both')
        if not (bstrees or constree):
            raise TreeError('Specify either list of boostrap trees or consensus tree.')
        # no outgroup specified: use the smallest clade of the root
        if outgroup is None:
            try:
                succnodes=self.node(self.root).succ
                smallest=min([(len(self.get_taxa(n)),n) for n in succnodes]) 
                outgroup=self.get_taxa(smallest[1])
            except:
                raise TreeError("Error determining outgroup.")
        else: # root with user specified outgroup
            self.root_with_outgroup(outgroup)

        if bstrees: # calculate consensus 
            constree=consensus(bstrees,threshold=threshold,outgroup=outgroup)
        else:
            if not constree.has_support():
                constree.branchlength2support()
            constree.root_with_outgroup(outgroup)
        # now we travel all nodes, and add support from consensus, if the clade is present in both
        for pnode in self._walk():
            cnode=constree.is_monophyletic(self.get_taxa(pnode))
            if cnode>-1:
                self.node(pnode).data.support=constree.node(cnode).data.support

         
def consensus(trees, threshold=0.5,outgroup=None):
    """Compute a majority rule consensus tree of all clades with relative frequency>=threshold from a list of trees."""

    total=len(trees)
    if total==0:
        return None
    # shouldn't we make sure that it's NodeData or subclass??
    dataclass=trees[0].dataclass
    max_support=trees[0].max_support
    clades={}
    #countclades={}
    alltaxa=set(trees[0].get_taxa())
    # calculate calde frequencies
    c=0
    for t in trees:
        c+=1
        #if c%100==0:
        #    print c
        if alltaxa!=set(t.get_taxa()):
            raise TreeError('Trees for consensus must contain the same taxa')
        t.root_with_outgroup(outgroup=outgroup)
        for st_node in t._walk(t.root):
            subclade_taxa=t.get_taxa(st_node)
            subclade_taxa.sort()
            subclade_taxa=str(subclade_taxa) # lists are not hashable
            if subclade_taxa in clades:
                clades[subclade_taxa]+=float(t.weight)/total
            else:
                clades[subclade_taxa]=float(t.weight)/total
            #if subclade_taxa in countclades:
            #    countclades[subclade_taxa]+=t.weight
            #else:
            #    countclades[subclade_taxa]=t.weight
    # weed out clades below threshold
    delclades=[c for c,p in clades.iteritems() if round(p,3)<threshold] # round can be necessary 
    for c in delclades:
        del clades[c]
    # create a tree with a root node
    consensus=Tree(name='consensus_%2.1f' % float(threshold),data=dataclass)
    # each clade needs a node in the new tree, add them as isolated nodes
    for c, s in clades.iteritems():
        node=Nodes.Node(data=dataclass())
        node.data.support=s
        node.data.taxon=set(eval(c))
        consensus.add(node)
    # set root node data
    consensus.node(consensus.root).data.support=None
    consensus.node(consensus.root).data.taxon=alltaxa
    # we sort the nodes by no. of taxa in the clade, so root will be the last
    consensus_ids=consensus.all_ids()
    consensus_ids.sort(lambda x,y:len(consensus.node(x).data.taxon)-len(consensus.node(y).data.taxon))
    # now we just have to hook each node to the next smallest node that includes all taxa of the current 
    for i,current in enumerate(consensus_ids[:-1]): # skip the last one which is the root
        #print '----'
        #print 'current: ',consensus.node(current).data.taxon
        # search remaining nodes
        for parent in consensus_ids[i+1:]:
            #print 'parent: ',consensus.node(parent).data.taxon
            if consensus.node(parent).data.taxon.issuperset(consensus.node(current).data.taxon):
                break
        else:
            sys.exit('corrupt tree structure?')
        # internal nodes don't have taxa
        if len(consensus.node(current).data.taxon)==1:
            consensus.node(current).data.taxon=consensus.node(current).data.taxon.pop()
            # reset the support for terminal nodes to maximum
            #consensus.node(current).data.support=max_support
        else:
            consensus.node(current).data.taxon=None
        consensus.link(parent,current)
    # eliminate root taxon name
    consensus.node(consensus_ids[-1]).data.taxon=None 
    if alltaxa != set(consensus.get_taxa()):
        raise TreeError('FATAL ERROR: consensus tree is corrupt') 
    return consensus

