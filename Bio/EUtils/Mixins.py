"""implements functionality shared between HistoryClient and DBIdsClient"""
# These are mixins and use methods implemented by the main classes.

import parse

class LinkMixin:
    def prlinks(self,
                db = "pubmed",
                term = None,
                field = None,
                daterange = None):
        """get the prlinks as a Datatypes.LinksLinkSet"""
        infile = self.elink(db = db,
                            cmd = "prlinks",
                            term = term,
                            field = field,
                            daterange = daterange)
        return parse.parse_prlinks(infile)

    def llinks(self,
                db = "pubmed",
                term = None,
                field = None,
                daterange = None):
        """get the llinks as a Datatypes.LinksLinkSet"""
        infile = self.elink(db = db,
                            cmd = "llinks",
                            term = term,
                            field = field,
                            daterange = daterange)
        return parse.parse_llinks(infile)

    def lcheck(self,
                db = "pubmed",
                term = None,
                field = None,
                daterange = None):
        infile = self.elink(db = db,
                            cmd = "lcheck",
                            term = term,
                            field = field,
                            daterange = daterange)
        return parse.parse_lcheck(infile)

    def ncheck(self,
                db = "pubmed",
                term = None,
                field = None,
                daterange = None):
        infile = self.elink(db = db,
                            cmd = "ncheck",
                            term = term,
                            field = field,
                            daterange = daterange)
        return parse.parse_ncheck(infile)

    def neighbor_links(self,
                db = "pubmed",
                term = None,
                field = None,
                daterange = None):
        infile = self.elink(db = db,
                            cmd = "neighbor",
                            term = term,
                            field = field,
                            daterange = daterange)
        return parse.parse_neighbor_links(infile)

class SequenceFetchMixin(LinkMixin):
    pass


## Should turn this into a Bio.SeqRecord object
## Could do that using the XML code, but GenBank/GenPept is better.
##    def fetch(self,
##              seq_start = None, seq_stop = None, strand = None,
##              complexity = None):
##        return parse.parse_fetch_xml(self.efetch(seq_start = seq_start,
##                                                 seq_stop = seq_stop,
##                                                 strand = strand,
##                                                 complexity = complexity))

class PublicationFetchMixin(LinkMixin):
    def fetch(self):
        return parse.parse_fetch_publication_xml(self.efetch("xml"))
