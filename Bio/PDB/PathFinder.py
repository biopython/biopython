from copy import copy

from PDBExceptions import PDBException


class PathFinder:
	"""
	Finds all connected paths in a unidirectional graph where each node
	is involved in at most two edges, and no cycles are present.
	"""

	def __init__(self):
		self.edges={}
		self.inverted_edges={}

	def _search(self, key, edges, inverted_edges):
		node_list=[]
		while 1:
			if edges.has_key(key):
				new_key=edges[key]
				node_list.append(new_key)
				del edges[key]
				del inverted_edges[new_key]
				key=new_key
			else:
				break
		return node_list

	def add_edge(self, node1, node2):
		"""
		Add an edge.
		"""
		if self.edges.has_key(node1):
			raise PDBException, "node1 already present."
		if self.inverted_edges.has_key(node2):
			raise PDBException, "node2 already present."
		if node1==node2:
			raise PDBException, "node1 equals node2."
		self.edges[node1]=node2
		self.inverted_edges[node2]=node1

	def get_paths(self):
		# we need a copy because the search method
		# destroys the dictionaries
		inverted_edges=copy(self.inverted_edges)
		edges=copy(self.edges)
		path_list=[]
		while 1:
			keys=edges.keys()
			if len(keys)==0:
				break
			else:
				key=keys[0]
				# note that search alters edges, inverted_edges
				fw_list=self._search(key, edges, inverted_edges)
				bw_list=self._search(key, inverted_edges, edges)
				bw_list.reverse()
				if len(fw_list)+len(bw_list)>0:
					node_list=bw_list+[key]+fw_list
					path_list.append(node_list)
		return path_list


if __name__=="__main__":

	from random import random

	pf=PathFinder()

	pf.add_edge(5,6)

	pf.add_edge(1,2)
	pf.add_edge(2,3)
	pf.add_edge(3,4)

	pf.add_edge(6,7)
	pf.add_edge(7,8)

	pf.add_edge(9,10)

	print pf.get_paths()

	# generate random graph
	pf=PathFinder()
	for i in range(0, 1000):
		a=int(100*random())
		b=int(100*random())
		try:
			pf.add_edge(a, b)
		except PDBException:
			pass
	l=pf.get_paths()
	
