# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eFetch
# from the PubMed database, as specified by NCBI's DTD file
# pubmed_080101.dtd (2007-11-30 16:19:51)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.


"""
This code is used to parse XML results returned by Entrez's PubMed
Article database.

The code is not meant to be used by itself, but by calling
Bio.Entrez.read() instead.  When parsing a PubMed Article XML file
in this way, the record returned is a list of dictionaries.

For example,

from Bio import Entrez
handle = Entrez.efetch(db="pubmed", id="17238260", retmode="XML")
articles = Entrez.read(handle)
assert len(articles)==1
"""

error = None

booleans = (
)

integers = (
)

strings = (
    "ArticleId",	# (#PCDATA) ATTLIST IdType %art.id.type; "pubmed"
    "ArticleTitle",	# <!ENTITY % ArticleTitle.Ref "ArticleTitle">
    "ISSN",		# (#PCDATA)
			# ATTLIST IssnType  (Electronic | Print
			#                    | Undetermined) #REQUIRED
    "Param",		# (#PCDATA)
			# ATTLIST Name CDATA #REQUIRED
    "PublicationStatus",	# (#PCDATA)
    "URL",		# (#PCDATA)
			# ATTLIST
			#  lang %iso.language.codes; #IMPLIED
			#  Type ( FullText | Summary | fulltext
			#       | summary) #IMPLIED
)

lists = (
    "ArticleIdList",	# (ArticleId+)
    "History",		# (PubMedPubDate+)
    "Object",		# (Param)*
			#  ATTLIST Type CDATA #REQUIRED
    "ObjectList",	# (Object)+
    "PubmedArticleSet",	# (PubmedArticle)+
)

dictionaries = (
    "PubmedArticle",	# ((NCBIArticle | MedlineCitation), PubmedData?)
    "PubDate",		# <!ENTITY % Pub.Date.Ref "PubDate?">
    "PubMedPubDate",	# (%normal.date;)
			# ATTLIST PubStatus %pub.status; #REQUIRED
)

structures = {
    "PubmedData": ["History"],	# (History*, PublicationStatus, ArticleIdList,
				# ObjectList?)
}

items = ()



if __name__ == "__main__" :
    print "Quick example/test"
    from Bio import Entrez
    handle = Entrez.efetch(db="pubmed", id=["14630660","17238260"], retmode="XML")
    articles = Entrez.read(handle)
    assert len(articles)==2
    for article in articles :
        print
        for key,value in article.iteritems() :
            if isinstance(value,list) :
                print "%s - %s" % (key,value[0])
                for entry in value[1:] :
                    print " "*len(key) + " - %s" % entry
            else :
                print "%s = %s" % (key,value)
