# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eFetch
# from the Entrez Gene database. Specifically, this module was written for
# NCBI's DTD file NCBI_Entrezgene.mod.dtd (2007-01-18 23:07:18)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

def startElement(self, name, attrs):
    if name=="Entrezgene-Set":
        object = []
        self.path = []
        self.record = object
    elif name=="Entrezgene":
        object = {}
        self.path[-1].append(object)
    elif name in ("Entrezgene_track-info",
                  "Gene-track",
                  "Entrezgene_gene",
                  "Gene-ref"):
        object = {}
        self.path[-1][name] = {}
    elif name in ("Gene-track_status",
                  "Entrezgene_type"):
        object = [attrs["value"],None]
        self.path[-1][name] = object
    else:
        object = ""
    self.path.append(object)

def endElement(self, name):
    object = self.path.pop()
    if name in ("Entrezgene_summary",
                "Gene-ref_locus",
                "Gene-track_geneid"):
        self.path[-1][name] = str(self.content)
    elif name in ("Gene-track_status",
                  "Entrezgene_type"):
        self.path[-1][name][1] = int(self.content)
