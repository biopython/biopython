#!/usr/bin/env python
# Copyright 2002 by Tiaan Wessels.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This file implements a SOAP interface into the Bibliographic database of
the European Bioinformatics Institute. It is a low-level interface and is
intended to be used by higher-level objects to build object representations
from data retrieved by this interface. This file needs a version of the
pywebsvcs package LATER than 0.97 available from www.sourceforge.net.
"""

import warnings
warnings.warn("Bio.biblio is deprecated because it seems to be out of date, "
              "and no one came forward saying they use this module. If you "
              "use Bio.biblio, please join the biopython mailing list and "
              "email us.",
              DeprecationWarning)

import sys
import SOAP
import copy

#SOAP.Config.debug = 1
SOAP.Config.BuildWithNoType = 1
SOAP.Config.BuildWithNoNamespacePrefix = 1

namespace = 'http://industry.ebi.ac.uk/openBQS'

class Biblio:
	"""
	this class provides low-level access to the EBI Bibliographics services exported
	through SOAP. there exist an almost 1-to-1 mapping between the methods and the
	RPC's available on the SOAP server.
	"""

	def __init__(self, serverurl):
		self.serverurl = serverurl
		self.namespace = namespace

		self.server = SOAP.SOAPProxy(
			self.serverurl,
			namespace = self.namespace
		)

	def get_count(self, cid):
		if cid == -1:
			return self.server.getBibRefCount()
		else:
			return self.server.getBibRefCount(cid)

	def find(self, cid, keywords, attrs, criteria):
		if cid == -1:
			return self.server.find([keywords, attrs, criteria])
		else:
			return self.server.find(cid, [keywords, attrs, criteria])

	def reset_retrieval(self, cid):
		if cid == -1:
			raise 'no collection id'
		self.server.resetRetrieval(cid)

	def has_next(self, cid):
		if cid == -1:
			raise 'no collection id'
		return self.server.hasNext(cid)

	def get_next(self, cid):
		if cid == -1:
			raise 'no collection id'
		return self.server.getNext(cid)

	def get_more(self, cid, cnt):
		if cid == -1:
			raise 'no collection id'
		if cnt <= 0:
			raise 'invalid count' + cnt
		return self.server.getMore(cid, cnt)

	def get_all_ids(self, cid):
		if cid == -1:
			raise 'no collection id (large result safeguard)'
		return self.server.getAllIDs(cid)

	def get_all(self, cid):
		if cid == -1:
			raise 'no collection id (large result safeguard)'
		return self.server.getAllBibRefs(cid)

	def get_by_id(self, id):
		return self.server.getById(id)

	def exists(self, cid):
		if cid == -1:
			raise 'no collection id'
		return self.server.exists(cid)

	def destroy(self, cid):
		if cid == -1:
			raise 'no collection id'
		self.server.destroy(cid)

	def get_vocabulary_names(self):
		return self.server.getAllVocabularyNames()

	def get_all_values(self, vocab):
		return self.server.getAllValues(vocab)

	def get_entry_description(self, vocab, entry):
		return self.server.getEntryDescription(vocab, entry)

	def contains(self, vocab, entry):
		return self.server.contains(vocab, entry)


class BiblioCollection:
	"""
	this class attempts to hide the concept of a collection id from users. each
	find action's results are grouped in the server under a unique collection id.
	this id could be used in subsequent calls to refine its content more by
	entering more specific search criteria. it can also return a new collection
	by using the subcollection method. each collection has its own current
	collection id as returned by the SOAP server by using the lower level Biblio
	class's services and it takes care of freeing this collection in the server
	upon destruction.
	"""

	def __init__(self, biblio, cid = -1):
		self.biblio = biblio
		self.cid = cid

	def __del__(self):
		self.destroy()
		
	def get_collection_id(self):
		if self.cid == -1:
			raise 'no collection id (use find)'
		return self.cid

	def get_count(self):
		return self.biblio.get_count(self.cid)

	def refine(self, keywords, attrs, criteria):
		remembercid = self.cid
		self.cid = self.biblio.find(self.cid, keywords, attrs, criteria)
		if remembercid != -1:	
			self.biblio.destroy(remembercid)

	def subcollection(self, keywords, attrs, criteria):
		return BiblioCollection(
			self.biblio, self.biblio.find(self.cid, keywords, attrs, criteria)
		)

	def get_all_ids(self):
		return self.biblio.get_all_ids(self.cid)

	def get_all(self):
		return self.biblio.get_all(self.cid)

	def exists(self):
		if self.cid == -1:
			return 0
		return self.biblio.exists(self.cid)

	def destroy(self):
		if self.cid == -1:
			return
		self.biblio.destroy(self.cid)

def checkargv(idx, msg):
	if idx-1 > len(sys.argv):
		raise 'argument expected at position %d for option %s' % (idx, msg)

def main():
	"""
	this function implements a command-line utility using the classes implemented
	in this file above and serves as base-line access software to query the BQS.
	"""

	# the default server location at EBI
	serverurl = 'http://industry.ebi.ac.uk/soap/openBQS'

	# is help requested
	try:
		sys.argv.index('-h')
	except:
		pass
	else:
		print """
usage: biblio.py [options] [- [finds]]
where options may be:
	-l <server URL> to change the server URL
	   (the default URL is %s)
	-g <citation ID> to get the XML version of a citation
	-c to obtain the size of a citation collection with each refinement
	-a to retrieve the citations in a collection instead of showing only
	   their citation id's
	-f <prefix> to specify the location whereto dump citations (implies -a)
	   found in a collection
	-o get citations one-by-one i.e. each will end up in its own file if used
	   in conjunction with -f
	-Vn to get the vocabulary names in the database
	-Vv <vocabulary> to get all antries for vocabulary
	-Vd <vocabulary> <entry> to get description for vocabulary entry
	-Ve <vocabulary> <entry> to determine whether vocabulary entry exist
and finds are any number of successive occurrences of the following:
	-find <keyword> [-attr <attribute>]
where each new find occurrence refines the result of the previous
examples of using this script is:
	biblio.py -l http://192.168.0.163:8123 -g 21322295
	biblio.py -g 21322295
	biblio.py -a - -find study -find gene
	biblio.py -f genestudies - -find study -find gene
	biblio.py -f brazma - -find brazma -attr author
	biblio.py -Vn
	biblio.py -Vv MEDLINE/Person/properties
	biblio.py -Vd MEDLINE/Person/properties LAST_NAME
	biblio.py -Ve MEDLINE/Person/properties LAST_NAME
""" % serverurl
		sys.exit

	# make server. (see if different server URL specified with -l)
	idx = 0
	try:
		idx = sys.argv.index('-l')
	except:
		pass
	else:
		checkargv(idx+1, '-l')
		serverurl = sys.argv[idx+1]
	server = Biblio(serverurl)

	# handle all the possible command-line options:

	# get a citation by its id
	try:
		idx = sys.argv.index('-g')
	except:
		pass
	else:
		checkargv(idx+1, '-g')
		print server.get_by_id(sys.argv[idx+1])

	# size of citation collection
	showsize = 0
	try:
		idx = sys.argv.index('-c')
	except:
		pass
	else:
		print 'total number of citations ->', server.get_count()
		showsize = 1

	# get all citations in collection
	fetch = 0
	try:
		idx = sys.argv.index('-a')
	except:
		pass
	else:
		fetch = 1

	# dump to a file with prefix
	prefix = None
	try:
		idx = sys.argv.index('-f')
	except:
		pass
	else:
		checkargv(idx+1, '-f')
		prefix = sys.argv[idx+1]
		fetch = 1

	# get individualy ?
	indiv = 0
	try:
		idx = sys.argv.index('-o')
	except:
		pass
	else:
		checkargv(idx+1, '-o')
		indiv = 1

	# get vocabulary names
	try:
		idx = sys.argv.index('-Vn')
	except:
		pass
	else:
		checkargv(idx+1, '-Vn')
		vocab = server.get_vocabulary_names()
		if len(vocab) > 0:
			print 'the vocabulary names are:'
		else:
			print 'there is no names in the vocabulary.'
		for v in vocab:
			print v

	# get entries for vocabulary name
	try:
		idx = sys.argv.index('-Vv')
	except:
		pass
	else:
		checkargv(idx+1, '-Vv')
		values = server.get_all_values(sys.argv[idx+1])
		if len(values) > 0:
			print 'the vocabulary entries for %s are:' % sys.argv[idx+1]
		else:
			print 'there is no entries in the vocabulary %s.' % sys.argv[idx+1]
		for v in values:
			print v

	# get entry for vocabulary entry for name
	try:
		idx = sys.argv.index('-Vd')
	except:
		pass
	else:
		checkargv(idx+1, '-Vd name')
		checkargv(idx+2, '-Vd entry')
		print server.get_entry_description(sys.argv[idx+1], sys.argv[idx+2])

	# vocabulary entry for name exist ?
	try:
		idx = sys.argv.index('-Ve')
	except:
		pass
	else:
		checkargv(idx+1, '-Ve name')
		checkargv(idx+2, '-Ve entry')
		if server.contains(sys.argv[idx+1], sys.argv[idx+2]):
			print 'entry %s::%s exists.' % (sys.argv[idx+1], sys.argv[idx+2])
		else:
			print 'entry %s::%s doesn\'t exists.' % (sys.argv[idx+1], sys.argv[idx+2])

	# - separates find from rest so this is a Rubicon
	base = 0
	try:
		idx = sys.argv.index('-')
	except:
		sys.exit
	else:
		base = idx

	# handle the find's (each successive find refines previous)
	collection = BiblioCollection(server)
	while 1:
		attrs = ''
		keys = ''
		try:
			idx = sys.argv[base:].index('-find')
		except:
			break
		else:
			checkargv(base+idx+1, '-find')
			keys = sys.argv[base+idx+1]
			if len(sys.argv[base+idx+1:]) > 1:
				if sys.argv[base+idx+2] == '-attr':
					checkargv(base+idx+3, '-attr')
					attrs = sys.argv[base+idx+3]
			if fetch:
				collection.refine(keys, attrs, '')
			else:
				print 'search with:', keys, attrs
				collection.refine(keys, attrs, '')
				print 'collection ->', collection.get_collection_id()
				if showsize:
					print 'collection size is ->', collection.get_count()
				ids = collection.get_all_ids()
				if len(ids) > 0:
					print 'citations in collection ->'
				else:
					print 'no citations in collection.'
				for id in ids:
					print id
			base = base+idx+1

	if fetch:
		if prefix != None:
			if indiv:
				ids = collection.get_all_ids()
				for id in ids:
					print 'saving %s ...' % id
					fn = prefix + '-' + id + '.xml'
					try:
						f = open(fn, 'w')
					except:
						print 'failed to open %s.' % fn
					else:
						f.write(server.get_by_id(id))
						f.close()
			else:
				fn = prefix + '.xml'
				try:
					f = open(fn, 'w')
				except:
					print 'failed to open %s.' % fn
				else:
					f.write(collection.get_all())
					f.close()
		else:
			print collection.get_all()

if __name__ == "__main__":
	main()
