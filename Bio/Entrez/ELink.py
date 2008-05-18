# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eLink,
# as specified by NCBI's DTD file eLink_020511.dtd (2005-02-18 17:13:40)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

error = "ERROR"		# (#PCDATA)	<!-- .+ -->

booleans = (
)

integers = (
    "Score",		# (#PCDATA)	<!-- \d+ -->
    "Priority",		# Not in the DTD, but found in the XML output
)

strings = (
    "Attribute",	# (#PCDATA)	<!-- .+ -->
    "DbFrom",		# (#PCDATA)	<!-- \S+ -->
    "DbTo",		# (#PCDATA)	<!-- \S+ -->
    "IconUrl",		# (#PCDATA)	<!-- \S+ -->
    "Id",		# (#PCDATA)	<!-- \d+ -->
			# ATTLIST HasLinkOut  (Y|N)	#IMPLIED	
			#         HasNeighbor (Y|N)	#IMPLIED
    "Info",		# (#PCDATA)	<!-- .+ -->
    "LinkName",		# (#PCDATA)	<!-- \S+ -->
    "Name",		# (#PCDATA)	<!-- .+ -->
    "NameAbbr",		# (#PCDATA)	<!-- \S+ -->
    "SubjectType",	# (#PCDATA)	<!-- .+ -->
    "Url",		# (#PCDATA)	<!-- \S+ -->
    "MenuTag",		# Not in the DTD, but found in the XML output
    "HtmlTag",		# Not in the DTD, but found in the XML output
)

lists = (
    "eLinkResult",	# (LinkSet*, ERROR?)
    "IdCheckList",	# (Id*,ERROR?)
    "IdList",		# (Id*)
    "IdUrlList",	# (IdUrlSet*,ERROR?)
)

dictionaries = (
    "Link",		# (Id, Score?)
    "Provider",		# (Name, NameAbbr, Id, Url, IconUrl?)
    "LinkInfo",		# Not in the DTD, but found in the XML output
)

structures = {
    "IdUrlSet": ["ObjUrl"],	# (Id,(ObjUrl+|Info))
    "LinkSetDb": ["Link"],	# (DbTo, LinkName, (Link*|Info), ERROR?)
    "ObjUrl": ["SubjectType", "Attribute"],	# ( Url,
						#  IconUrl?,
						#  LinkName?,
						#  SubjectType*,
						#  Attribute*,
						#  Provider)
    "LinkSet": ["LinkSetDb"],	# (DbFrom, 
				#  ((IdList?, LinkSetDb*)
				#   | IdUrlList | IdCheckList | ERROR)  
    "IdLinkSet": ["LinkInfo"],	# Not in the DTD, but found in the XML output
}

items = ()
