#
# Trees.py 
#
# version 1.0
#
# Copyright 2005 by Frank Kauff & Cymon J. Cox. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.
#
# Tree class handles phylogenetic trees. Provides a set of methods to read and write newick-format tree
# descriptions, get information about trees (monphyly of taxon sets, congruence between trees, common ancestors,...)
# and to manipulate trees (reroot trees, split terminal nodes).
#
# Bug reports welcome: fkauff@duke.edu
#

import sys, random, sets
import Nodes

PRECISION_BRANCHLENGTH=6
PRECISION_SUPPORT=6

class TreeError(Exception): pass

class NodeData:
    """Stores tree-relevant data associated with nodes (e.g. branches or otus)."""
    def __init__(self,taxon=None,branchlength=0.0,support=None):
        self.taxon=taxon
        self.branchlength=branchlength
        self.support=support

class Tree(Nodes.Chain):
    """Represents a tree using a chain of nodes with on predecessor (=ancestor)
    and multiple successors (=subclades).
    """ 
    # A newick tree is parsed into nested list and then converted to a node list in two stages
    # mostly due to hsitorical reasons. This could be done in one swoop). Note: parentheses ( ) and
    # colon : are not allowed in taxon names. This is against NEXUS standard, but makes life much
    # easier when parsing trees.
    

    def __init__(self,tree=None,weight=1,rooted=False,name='',data=NodeData,values_are_support=False):
        """Ntree(self,tree)."""
        Nodes.Chain.__init__(self)
        self.__values_are_support=values_are_support
        self.weight=weight
        self.rooted=rooted
        self.name=name
        root=Nodes.Node(data())
        self.add(root)
        self.root=root.id
        if tree:    # use the tree we have
            self._add_subtree(parent_id=root.id,tree=self._parse(tree)[0],data=data)
        
    def _parse(self,tree):
        """Parses (a,b,c...)[[[xx]:]yy] into subcomponents and travels down recursively."""
        
        if tree.count('(')!=tree.count(')'):
            raise TreeError, 'Parentheses do not match in (sub)tree: '+tree
        if tree.count('(')==0: # a leaf
            colon=tree.rfind(':')   
            if colon>-1:
                return [tree[:colon],self._get_values(tree[colon+1:])]
            else:
                return [tree,[None]]
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
    
    def _add_subtree(self,parent_id=None,tree=None,data=NodeData):
        """Adds leaf or tree (in newick format) to a parent_id. (self,parent_id,tree)."""
        
        if parent_id is None:
            raise TreeError('Need node_id to connect to.')
        for st in tree:
            if type(st[0])==list: # it's a subtree
                nd=data()
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
                sn=Nodes.Node(nd)
                self.add(sn,parent_id)
                self._add_subtree(sn.id,st[0])
            else: # it's a leaf
                nd=data()
                nd.taxon=st[0]
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
                leaf=Nodes.Node(nd)
                self.add(leaf,parent_id)
    
    def _get_values(self, text):
        """Extracts values (support/branchlength) from xx[:yyy], xx."""
       
        if text=='':
            return None
        return [float(t) for t in text.split(':') if t.strip()] 
   
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
        parent_data=self.chain[parent_id].get_data()
        for i in range(n):
            node=Nodes.Node()
            if parent_data:
                # the data class is probably some subclass from NodeData we don't know
                # we need a new instance for that
                node.data=parent_data.__class__()
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
        for id,node in self.chain.items():
            if node.get_data().taxon==taxon:
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
            if not prev==self.root and len(self.node(prev).succ)==1:
                succ=self.node(prev).get_succ(0)
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
        if self.chain[node_id].get_succ()==[]:
            if self.chain[node_id].get_data():
                return [self.chain[node_id].get_data().taxon]
            else:
                return None
        else:
            list=[]
            for succ in self.chain[node_id].get_succ():
                list.extend(self.get_taxa(succ))
            return list

    def get_terminals(self):
        """Return a list of all terminal nodes."""
        return [i for i,n in self.chain.items() if n.succ==[]]

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
            blen+=self.chain[node].get_data().branchlength
            node=self.chain[node].get_prev()
        return blen

    def set_subtree(self,node):
        """Return subtree as a set of nested sets.

        sets = set_subtree(self,node)
        """
        
        if self.node(node).succ==[]:
            return self.node(node).data.taxon
        else:
            return sets.Set([self.set_subtree(n) for n in self.node(node).succ])
            
    def is_identical(self,tree2):
        """Compare tree and tree2 for identity.

        result = is_identical(self,tree2)
        """
        return self.set_subtree(self.root)==tree2.set_subtree(tree2.root)

    def is_compatible(self,tree2,threshold):
        """Compares branches with support>threshold for compatibility.
        
        result = is_compatible(self,tree2,threshold)
        """

        if sets.Set(self.get_taxa())!=sets.Set(tree2.get_taxa()):
            raise TreeError, 'Can\'t compare trees with different taxon compositions.'
        t1=[(sets.Set(self.get_taxa(n)),self.node(n).data.support) for n in self.chain.keys() if \
            self.node(n).succ and\
            (self.node(n).data and self.node(n).data.support and self.node(n).data.support>=threshold)]
        t2=[(sets.Set(tree2.get_taxa(n)),tree2.node(n).data.support) for n in tree2.chain.keys() if \
            tree2.node(n).succ and\
            (tree2.node(n).data and tree2.node(n).data.support and tree2.node(n).data.support>=threshold)]
        conflict=[]
        for (st1,sup1) in t1:
            for (st2,sup2) in t2:
                if not st1.issubset(st2) and not st2.issubset(st1):
                    intersect,notin1,notin2=st1 & st2, st2-st1, st1-st2
                    if intersect:
                        conflict.append((st1,sup1,st2,sup2,intersect,notin1,notin2))
        return conflict
        
    def common_ancestor(self,node1,node2):
        """Return the common ancestor that connects to nodes.
        
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
            taxon_set=sets.Set([taxon_list])
        else:
            taxon_set=sets.Set(taxon_list)
        node_id=self.root
        while 1:
            subclade_taxa=sets.Set(self.get_taxa(node_id))
            if subclade_taxa==taxon_set:                                        # are we there?
                return node_id
            else:                                                               # check subnodes
                for subnode in self.chain[node_id].get_succ():
                    if sets.Set(self.get_taxa(subnode)).issuperset(taxon_set):  # taxon_set is downstream
                        node_id=subnode
                        break   # out of for loop
                else:
                    return -1   # taxon set was not with successors, for loop exhausted

    def branchlength2support(self):
        """Move values stored in data.branchlength to data.support, and set branchlength to 0.0

        This is necessary when support has been stored as branchlength (e.g. paup), and has thus
        been read in as branchlength. 
        """

        for n in self.chain.keys():
            self.node(n).data.support=self.node(n).data.branchlength
            self.node(n).data.branchlength=0.0

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
        elif ntax<>len(taxon_list):
            raise TreeError('Length of taxon list must correspond to ntax.')
        # initiate self with empty root
        self.__init__(plain=False)
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
                    self.chain[nt].data.branchlength=bl
            terminals.extend(new_terminals)
            terminals.remove(newsplit)
        # distribute taxon labels randomly
        random.shuffle(taxon_list)
        for (node,name) in zip(terminals,taxon_list):
            self.chain[node].data.taxon=name

    def root_with_outgroup(self,outgroup):
        """Roots the tree with (monophyletic list of) outgroup taxa.
        
        new_root = root_with_outgroup(self,outgroup)
        """
        outgroup_node=self.is_monophyletic(outgroup)
        if outgroup_node==-1:
            return -1           # outgroup is not monophyletic
        elif outgroup_node in self.node(self.root).get_succ() or outgroup_node==self.root:
            return self.root    # don't try to root a tree where it's already rooted...
        #isolate the old root and get the subtrees
        subtrees=self.chain[self.root].get_succ()[:]
        newroot_subtree=False
        for st in subtrees:
            self.unlink(st)                           # cut the subtrees off
            if self.is_parent_of(st,outgroup_node):
                newroot_subtree=st
        if not newroot_subtree:                                # that's impossible!
            raise TreeError('Corrupt tree structure')
        # now kill the old root
        self.kill(self.root)                          # kill the now isolated root 
        new_root=Nodes.Node()                                    # create new root node
        outgroup_anc=self.chain[outgroup_node].get_prev()       # what's the outgroup connected to? 
        self.unlink(outgroup_node)                        # disconnect outgroup
        new_root_id=self.add(new_root)                    # add new root node to tree
        self.link(new_root_id,outgroup_node)              # add outgroup to new root
        # add the subtrees from the old root to the node which contains the subtree with the new root
        for s in subtrees:
            if s!=newroot_subtree:
                self.link(newroot_subtree,s)
        # now the branches between the old and the new root must be reoriented. 
        connection=[newroot_subtree]+self.trace(newroot_subtree,outgroup_anc)   # get all intermediate nodes
        connection.append(new_root_id)                                 # add the new root
        for i in range(1,len(connection)):                          
            self.unlink(connection[i])                        
            self.link(connection[i],connection[i-1])
        self.root=new_root_id
        return self.root
       
    def display(self):
        """Quick and dirty lists of all nodes."""
        table=[('#','taxon','prev','succ','brlen','blen (sum)','support')]
        for (i,n) in self.chain.items():
            if not n.get_data():
                table.append((str(i),'-',str(n.get_prev()),str(n.get_succ()),'-','-','-'))
            else:
                tx=n.get_data().taxon
                if not tx:
                    tx='-'
                blength=n.get_data().branchlength
                if blength is None:
                    blength='-'
                    sum_blength='-'
                else:
                    sum_blength=self.sum_branchlength(node=i)
                support=n.get_data().support
                if support is None:
                    support='-'
                table.append((str(i),tx,str(n.get_prev()),str(n.get_succ()),blength,sum_blength,support))
        print '\n'.join(['%3s %10s %15s %15s %8s %10s %8s' % l for l in table])
        print '\nRoot: ',self.root

    def to_string(self,support_as_branchlengths=False,branchlengths_only=False,plain=True):
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
            if self.plain: # plain tree only. That's easy.
                return ''
            elif self.support_as_branchlengths: # support as branchlengths (eg. PAUP), ignore actual branchlengths
                if terminal:    # terminal branches have 100% support
                    return ':100.00000'
                else:
                    return ':%1.5f' % (data.support)
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
                        return '%1.5f:0.00000' % (data.support)
                    else:
                        return '0.00000:0.00000'

        def newickize(node):
            """Convert a node tree to a newick tree recursively."""

            if not self.node(node).succ:    #terminal
                return self.node(node).data.taxon+make_info_string(self.node(node).data,terminal=True)
            else:
                return '(%s)%s' % (','.join(map(newickize,self.node(node).succ)),make_info_string(self.node(node).data))
            return subtree
                    
        treeline='tree '
        if self.name:
            treeline+=self.name
        else:
            treeline+='a_tree'
        if self.weight<>1:
            treeline+=' [&W%s]' % str(round(float(self.weight),3))
        if self.rooted:
            treeline+=' [&R]'
        treeline+=' = (%s);' % ','.join(map(newickize,self.node(self.root).succ))
        return treeline 
        
    def __str__(self):
        """Short version of to_string(), gives plain tree"""
        return self.to_string(plain=True)
            
