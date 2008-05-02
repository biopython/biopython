# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eFetch
# from the Entrez Gene database. Specifically, this module was written for
# NCBI's DTD file NCBI_Entrezgene.mod.dtd (2007-01-18 23:07:18)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

def startElement(self, name, attrs):
    if self.element==["Entrezgene-Set"]:
        self.record = []
    elif self.element==["Entrezgene-Set", "Entrezgene"]:
        entrezgene = {}
        self.record.append(entrezgene)
    elif self.element==["Entrezgene-Set", "Entrezgene", "Entrezgene_track-info"]:
        entrezgene = self.record[-1]
        entrezgene["track-info"] = {}
    elif self.element==["Entrezgene-Set", "Entrezgene", "Entrezgene_track-info", "Gene-track"]:
        entrezgene = self.record[-1]
        entrezgene["track-info"]["Gene-track"] = {}
    elif self.element==["Entrezgene-Set", "Entrezgene", "Entrezgene_track-info", "Gene-track", "Gene-track_status"]:
        entrezgene = self.record[-1]
        entrezgene["track-info"]["Gene-track"]["status"] = [attrs["value"],None]
    elif self.element==["Entrezgene-Set", "Entrezgene", "Entrezgene_type"]:
        entrezgene = self.record[-1]
        entrezgene["Entrezgene_type"] = [attrs["value"],None]
    elif self.element==["Entrezgene-Set", "Entrezgene", "Entrezgene_gene"]:
        entrezgene = self.record[-1]
        entrezgene["gene"] = {}
    elif self.element==["Entrezgene-Set", "Entrezgene", "Entrezgene_gene", "Gene-ref"]:
        entrezgene = self.record[-1]
        entrezgene["gene"]["ref"] = {}

def endElement(self, name):
    if self.element==["Entrezgene-Set", "Entrezgene", "Entrezgene_track-info", "Gene-track", "Gene-track_geneid"]:
        entrezgene = self.record[-1]
        entrezgene["track-info"]["Gene-track"]["geneid"] = str(self.content)
    elif self.element==["Entrezgene-Set", "Entrezgene", "Entrezgene_track-info", "Gene-track", "Gene-track_status"]:
        entrezgene = self.record[-1]
        entrezgene["track-info"]["Gene-track"]["status"][1] = int(self.content)
    elif self.element==["Entrezgene-Set", "Entrezgene", "Entrezgene_type"]:
        entrezgene = self.record[-1]
        entrezgene["Entrezgene_type"][1] = int(self.content)
    elif self.element==["Entrezgene-Set", "Entrezgene", "Entrezgene_summary"]:
        entrezgene = self.record[-1]
        entrezgene["summary"] = self.content
    elif self.element==["Entrezgene-Set", "Entrezgene", "Entrezgene_gene", "Gene-ref", "Gene-ref_locus"]:
        entrezgene = self.record[-1]
        entrezgene["gene"]["ref"]["locus"] = self.content
