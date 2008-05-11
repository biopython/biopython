# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eFetch from
# the OMIM database, as specified by NCBI's DTD file NCBI_Mim.dtd
# (06/06/2006 23:03:48)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

def startElement(self, name, attrs):
    if name=="Mim-entries":
        object = []
        self.path = []
        self.record = object
    elif name in ("Mim-allelic-variant_aliases",
                  "Mim-allelic-variant_snpLinks",
                  "Mim-entry_aliases",
                  "Mim-entry_allelicVariants",
                  "Mim-entry_attribution",
                  "Mim-entry_clinicalSynopsis",
                  "Mim-entry_editHistory",
                  "Mim-entry_genomeLinks",
                  "Mim-entry_included",
                  "Mim-entry_nucleotideLinks",
                  "Mim-entry_medlineLinks",
                  "Mim-entry_proteinLinks",
                  "Mim-entry_references",
                  "Mim-entry_seeAlso",
                  "Mim-entry_structureLinks",
                  "Mim-entry_summary",
                  "Mim-entry_summaryAttribution",
                  "Mim-entry_summaryCreationDate",
                  "Mim-entry_summaryEditHistory",
                  "Mim-entry_synonyms",
                  "Mim-entry_synopsis",
                  "Mim-entry_synopsisAttribution",
                  "Mim-entry_synopsisCreationDate",
                  "Mim-entry_synopsisEditHistory",
                  "Mim-entry_text",
                  "Mim-entry_textfields",
                  "Mim-index-term_terms",
                  "Mim-reference_authors",
                  "Mim-text_neighbors"):
        object = []
        self.path[-1][name] = object
    elif name in ("Mim-allelic-variant",
                  "Mim-author",
                  "Mim-cit",
                  "Mim-entry",
                  "Mim-index-term",
                  "Mim-link",
                  "Mim-reference",
                  "Mim-text"):
        object = {}
        self.path[-1].append(object)
    elif name in ("Mim-date",
                  "Mim-edit-item_modDate",
                  "Mim-entry_creationDate",
                  "Mim-page",
                  "Mim-reference_pages",
                  "Mim-reference_pubDate"):
        object = {}
        self.path[-1][name] = object
    elif name=="Mim-entry_mimType":
        object = [str(attrs["value"]), None]
        self.path[-1][name] = object
    elif name=="Mim-reference_type":
        object = str(attrs["value"])
        self.path[-1][name] = object
    elif name in ("Mim-entry_hasSummary",
                  "Mim-entry_hasSynopsis",
                  "Mim-reference_ambiguous",
                  "Mim-reference_noLink"):
        if str(attrs["value"])=="true":
            object = True
        elif str(attrs["value"])=="false":
            object = False
        self.path[-1][name] = object
    elif name=="Mim-edit-item":
        object = {}
        if type(self.path[-1])==list:
            self.path[-1].append(object)
        elif type(self.path[-1])==dict:
            self.path[-1][name] = object
    else:
        object = ""
    self.path.append(object)

def endElement(self, name):
    self.path.pop()
    if name in ("Mim-allelic-variant_name",
                "Mim-allelic-variant_mutation",
                "Mim-allelic-variant_description",
                "Mim-author_name",
                "Mim-cit_author",
                "Mim-cit_others",
                "Mim-cit_year",
                "Mim-date_day",
                "Mim-date_month",
                "Mim-date_year",
                "Mim-edit-item_author",
                "Mim-entry_copyright",
                "Mim-entry_locus",
                "Mim-entry_mimNumber",
                "Mim-entry_symbol",
                "Mim-entry_title",
                "Mim-index-term_key",
                "Mim-link_uids",
                "Mim-reference_bookTitle",
                "Mim-reference_citationTitle",
                "Mim-reference_citationType",
                "Mim-reference_commNote",
                "Mim-reference_edition",
                "Mim-reference_editors",
                "Mim-reference_journal",
                "Mim-reference_miscInfo",
                "Mim-reference_number",
                "Mim-reference_origNumber",
                "Mim-reference_otherAuthors",
                "Mim-reference_place",
                "Mim-reference_primaryAuthor",
                "Mim-reference_publisher",
                "Mim-reference_pubmedUID",
                "Mim-reference_series",
                "Mim-reference_volume",
                "Mim-text_label",
                "Mim-text_text"):
        self.path[-1][name] = self.content
    elif name in ("Mim-author_index",
                  "Mim-allelic-variant_number",
                  "Mim-cit_number",
                  "Mim-entry_numGeneMaps",
                  "Mim-link_num",
                  "Mim-link_numRelevant",
                  "Mim-page_from",
                  "Mim-page_to"):
        self.path[-1][name] = int(self.content)
    elif name in ("Mim-index-term_terms_E",
                  "Mim-allelic-variant_aliases_E",
                  "Mim-entry_synonyms_E",
                  "Mim-entry_aliases_E",
                  "Mim-entry_included_E"):
        self.path[-1].append(self.content)
    elif name=="Mim-entry_mimType":
        self.path[-1][name][1] = self.content
