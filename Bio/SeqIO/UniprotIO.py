# Copyright 2010 by Andrea Pierleoni
# Revisions copyright 2010 by Peter Cock
# All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the "uniprot-xml" file format.

See also:

http://www.uniprot.org

The UniProt XML format essentially replaces the old plain text file format
originally introduced by SwissProt ("swiss" format in Bio.SeqIO).
"""
import sys

from Bio import Seq
from Bio import SeqFeature
from Bio import Alphabet
from Bio.SeqRecord import SeqRecord
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
import warnings

#For speed try to use cElementTree rather than ElemenTree
try:
    if (3,0,0) <= sys.version_info[:3] <= (3,1,3):
        #workaround for bug in python 3 to 3.1.3  see http://bugs.python.org/issue9257
        from xml.etree import ElementTree as ElementTree
    else:
        from xml.etree import cElementTree as ElementTree
except ImportError:
    from xml.etree import ElementTree as ElementTree

NS = "{http://uniprot.org/uniprot}"
REFERENCE_JOURNAL = "%(name)s %(volume)s:%(first)s-%(last)s(%(pub_date)s)"

def UniprotIterator(handle, alphabet=Alphabet.ProteinAlphabet(), return_raw_comments=False):
    """Generator function to parse UniProt XML as SeqRecord objects.
    
    parses an XML entry at a time from any UniProt XML file 
    returns a SeqRecord for each iteration
    
    This generator can be used in Bio.SeqIO
    
    return_raw_comments = True --> comment fields are returned as complete xml to allow further processing
    skip_parsing_errors = True --> if parsing errors are found, skip to next entry
    """
    if isinstance(alphabet, Alphabet.NucleotideAlphabet):
        raise ValueError, "Wrong alphabet %r" % alphabet
    if isinstance(alphabet, Alphabet.Gapped):
        if isinstance(alphabet.alphabet, Alphabet.NucleotideAlphabet):
            raise ValueError, "Wrong alphabet %r" % alphabet

    if not hasattr(handle, "read"):
        if type(handle)==type(''):
            handle=StringIO(handle)
        else:
            raise Exception('An XML-containing handler or an XML string must be passed')

    if ElementTree is None:
        from Bio import MissingExternalDependencyError
        raise MissingExternalDependencyError(
                "No ElementTree module was found. "
                "Use Python 2.5+, lxml or elementtree if you "
                "want to use Bio.SeqIO.UniprotIO.")
        
    for event, elem in ElementTree.iterparse(handle, events=("start", "end")):
        if event=="end" and elem.tag == NS + "entry":
            yield Parser(elem, alphabet=alphabet, return_raw_comments=return_raw_comments).parse()
            elem.clear()

class Parser(object):
    """Parse a UniProt XML entry to a SeqRecord.
    
    return_raw_comments=True to get back the complete comment field in XML format
    alphabet=Alphabet.ProteinAlphabet()    can be modified if needed, default is protein alphabet.
    """
    def __init__(self, elem, alphabet=Alphabet.ProteinAlphabet(), return_raw_comments=False):
        self.entry=elem
        self.alphabet=alphabet
        self.return_raw_comments=return_raw_comments
    
    def parse(self):
        """Parse the input."""
        assert self.entry.tag == NS + 'entry'
        
        def append_to_annotations(key, value):
            if key not in self.ParsedSeqRecord.annotations:
                self.ParsedSeqRecord.annotations[key]=[]
            if value not in self.ParsedSeqRecord.annotations[key]:
                self.ParsedSeqRecord.annotations[key].append(value)
            
        def _parse_name(element):
            self.ParsedSeqRecord.name=element.text
            self.ParsedSeqRecord.dbxrefs.append(self.dbname+':'+element.text)
        
        def _parse_accession(element):
            append_to_annotations('accessions', element.text)# to cope with SwissProt plain text parser
            self.ParsedSeqRecord.dbxrefs.append(self.dbname+':'+element.text)
        
        def _parse_protein(element):
            """Parse protein names (PRIVATE)."""
            descr_set=False
            for protein_element in element.getchildren():
                if protein_element.tag in [NS + 'recommendedName', NS + 'alternativeName']:#recommendedName tag are parsed before
                    #use protein fields for name and description
                    for rec_name in protein_element.getchildren():
                        ann_key='%s_%s' % (protein_element.tag.replace(NS,''), rec_name.tag.replace(NS,''))
                        append_to_annotations(ann_key, rec_name.text)
                        if (rec_name.tag==NS + 'fullName') and not descr_set:
                            self.ParsedSeqRecord.description=rec_name.text
                            descr_set=True
                elif protein_element.tag==NS + 'component':
                    pass #not parsed 
                elif protein_element.tag==NS + 'domain':
                    pass #not parsed 
        
        def _parse_gene(element):
            for genename_element in element.getchildren():  
                if 'type' in genename_element.attrib:
                    ann_key='gene_%s_%s' % (genename_element.tag.replace(NS,''), genename_element.attrib['type'])
                    if genename_element.attrib['type']=='primary':
                        self.ParsedSeqRecord.annotations[ann_key]=genename_element.text
                    else:
                        append_to_annotations(ann_key,genename_element.text)
        
        def _parse_geneLocation(element):
            append_to_annotations('geneLocation', element.attrib['type'])
        
        def _parse_organism(element):
            organism_name = com_name = sci_name = ''
            for organism_element in element.getchildren():  
                if organism_element.tag==NS + 'name':
                    if organism_element.text:
                        if organism_element.attrib['type'] == 'scientific':
                            sci_name = organism_element.text
                        elif organism_element.attrib['type'] == 'common':
                            com_name = organism_element.text
                        else:
                            #e.g. synonym
                            append_to_annotations("organism_name", organism_element.text)
                elif organism_element.tag==NS + 'dbReference':
                    self.ParsedSeqRecord.dbxrefs.append(organism_element.attrib['type']+':'+organism_element.attrib['id'])
                elif organism_element.tag==NS + 'lineage':
                    for taxon_element in organism_element.getchildren():
                        if taxon_element.tag==NS + 'taxon':
                            append_to_annotations('taxonomy',taxon_element.text)
            if sci_name and com_name:
                organism_name = '%s (%s)' % (sci_name, com_name)
            elif sci_name:
                organism_name = sci_name
            elif com_name:
                organism_name = com_name
            self.ParsedSeqRecord.annotations['organism']=organism_name
            
        def _parse_organismHost(element):
            for organism_element in element.getchildren():  
                if organism_element.tag==NS + 'name': 
                    append_to_annotations("organism_host", organism_element.text)
                        
        def _parse_keyword(element):      
            append_to_annotations('keywords',element.text)
        
        def _parse_comment(element):
            """Parse comments (PRIVATE).
            
            Comment fields are very heterogeneus. each type has his own (frequently mutated) schema.
            To store all the contained data, more complex data structures are needed, such as 
            annidated dictionaries. This is left to end user, by optionally setting:
            
            return_raw_comments=True 
            
            the orginal XMLs is returned in the annotation fields.
            
            available comment types at december 2009:
                "allergen"
                "alternative products"
                "biotechnology"
                "biophysicochemical properties"
                "catalytic activity"
                "caution"
                "cofactor"
                "developmental stage"
                "disease"
                "domain"
                "disruption phenotype"
                "enzyme regulation"
                "function"
                "induction"
                "miscellaneous"
                "pathway"
                "pharmaceutical"
                "polymorphism"
                "PTM"
                "RNA editing"
                "similarity"
                "subcellular location"
                "sequence caution"
                "subunit"
                "tissue specificity"
                "toxic dose"
                "online information"
                "mass spectrometry"
                "interaction"
            """
            
            simple_comments=["allergen",
                            "biotechnology",
                            "biophysicochemical properties",
                            "catalytic activity",
                            "caution",
                            "cofactor",
                            "developmental stage",
                            "disease",
                            "domain",
                            "disruption phenotype",
                            "enzyme regulation",
                            "function",
                            "induction",
                            "miscellaneous",
                            "pathway",
                            "pharmaceutical",
                            "polymorphism",
                            "PTM",
                            "RNA editing",#positions not parsed
                            "similarity",
                            "subunit",
                            "tissue specificity",
                            "toxic dose",
                             ]

            if element.attrib['type'] in simple_comments:
                ann_key='comment_%s' % element.attrib['type'].replace(' ','')
                for text_element in element.getiterator(NS + 'text'):
                    if text_element.text:
                        append_to_annotations(ann_key,text_element.text)
            elif element.attrib['type']=='subcellular location':
                for subloc_element in element.getiterator(NS + 'subcellularLocation'):
                    for el in subloc_element.getchildren():
                        if el.text:
                            ann_key='comment_%s_%s' % (element.attrib['type'].replace(' ',''), el.tag.replace(NS,''))
                            append_to_annotations(ann_key,el.text)
            elif element.attrib['type']=='interaction':
                for interact_element in element.getiterator(NS +'interactant'):
                    ann_key='comment_%s_intactId' % element.attrib['type']
                    append_to_annotations(ann_key,interact_element.attrib['intactId'])
            elif element.attrib['type']=='alternative products':
                for alt_element in element.getiterator(NS +'isoform'):
                    ann_key='comment_%s_isoform' % element.attrib['type'].replace(' ','')
                    for id_element in alt_element.getiterator(NS +'id'):
                        append_to_annotations(ann_key,id_element.text)
            elif element.attrib['type']=='mass spectrometry':
                ann_key='comment_%s' % element.attrib['type'].replace(' ','')
                start=end=0
                for loc_element in element.getiterator(NS +'location'):
                    pos_els=loc_element.getiterator(NS +'position')
                    pos_els=list(pos_els)
                    # this try should be avoided, maybe it is safer to skip postion parsing for mass spectrometry
                    try:
                        if pos_els:
                            end=int(pos_els[0].attrib['position'])
                            start=end-1
                        else:
                            start=int(loc_element.getiterator(NS +'begin')[0].attrib['position'])-1
                            end=int(loc_element.getiterator(NS +'end')[0].attrib['position'])
                    except :#undefined positions or erroneusly mapped
                        pass    
                mass=element.attrib['mass']
                method=element.attrib['mass'] #TODO - Check this, looks wrong!
                if start==end==0:  
                    append_to_annotations(ann_key,'undefined:%s|%s'%(mass,method))
                else:
                    append_to_annotations(ann_key,'%s..%s:%s|%s'%(start,end,mass,method))
            elif element.attrib['type']=='sequence caution':
                pass#not parsed: few information, complex structure
            elif element.attrib['type']=='online information':
                for link_element in element.getiterator(NS +'link'):
                    ann_key='comment_%s' % element.attrib['type'].replace(' ','')
                    for id_element in link_element.getiterator(NS +'link'):
                        append_to_annotations(ann_key,'%s@%s'%(element.attrib['name'],link_element.attrib['uri']))            
            
            #return raw XML comments if needed
            if self.return_raw_comments:
                ann_key='comment_%s_xml' % element.attrib['type'].replace(' ','')
                append_to_annotations(ann_key,ElementTree.tostring(element))
                
        
        def _parse_dbReference(element):
            self.ParsedSeqRecord.dbxrefs.append(element.attrib['type']+':'+element.attrib['id'])
            #e.g.
            # <dbReference type="PDB" key="11" id="2GEZ">
            #   <property value="X-ray" type="method"/>
            #   <property value="2.60 A" type="resolution"/>
            #   <property value="A/C/E/G=1-192, B/D/F/H=193-325" type="chains"/>
            # </dbReference>
            if 'type' in element.attrib:
                if element.attrib['type'] == 'PDB':
                        method=""
                        resolution=""
                        for ref_element in element.getchildren():  
                            if ref_element.tag==NS + 'property':
                                dat_type=ref_element.attrib['type']
                                if dat_type=='method':
                                    method=ref_element.attrib['value']
                                if dat_type=='resolution':
                                    resolution=ref_element.attrib['value']
                                if dat_type=='chains':
                                    pairs=ref_element.attrib['value'].split(',')
                                    for elem in pairs:
                                        pair=elem.strip().split('=')
                                        if pair[1]!='-':
                                            #TODO - How best to store these, do SeqFeatures make sense?
                                            feature=SeqFeature.SeqFeature()
                                            feature.type=element.attrib['type']
                                            feature.qualifiers['name']=element.attrib['id']
                                            feature.qualifiers['method']=method
                                            feature.qualifiers['resolution']=resolution
                                            feature.qualifiers['chains']=pair[0].split('/')
                                            start=int(pair[1].split('-')[0])-1
                                            end=int(pair[1].split('-')[1])
                                            feature.location=SeqFeature.FeatureLocation(start,end)
                                            #self.ParsedSeqRecord.features.append(feature)

            for ref_element in  element.getchildren():  
                if ref_element.tag==NS + 'property':
                    pass# this data cannot be fitted in a seqrecord object with a simple list. however at least ensembl and EMBL parsing can be improved to add entries in dbxrefs
            
        def _parse_reference(element):
            reference=SeqFeature.Reference()
            authors=[]
            scopes=[]
            tissues=[]
            journal_name=''
            pub_type=''
            pub_date=''
            for ref_element in element.getchildren():
                if ref_element.tag==NS + 'citation':
                    pub_type=ref_element.attrib['type']
                    if pub_type=='submission':
                        pub_type+=' to the '+ref_element.attrib['db']
                    if 'name' in ref_element.attrib:
                        journal_name=ref_element.attrib['name']
                    pub_date=ref_element.attrib.get('date','')
                    j_volume=ref_element.attrib.get('volume','')
                    j_first=ref_element.attrib.get('first','')
                    j_last=ref_element.attrib.get('last','')
                    for cit_element in ref_element.getchildren():
                        if cit_element.tag==NS + 'title':
                            reference.title=cit_element.text
                        elif cit_element.tag==NS + 'authorList':
                            for person_element in cit_element.getchildren():
                                authors.append(person_element.attrib['name'])
                        elif cit_element.tag==NS + 'dbReference':
                            self.ParsedSeqRecord.dbxrefs.append(cit_element.attrib['type']+':'+cit_element.attrib['id'])
                            if cit_element.attrib['type']=='PubMed':
                                reference.pubmed_id=cit_element.attrib['id']
                            elif ref_element.attrib['type']=='MEDLINE':
                                reference.medline_id=cit_element.attrib['id']
                elif ref_element.tag==NS + 'scope':
                    scopes.append(ref_element.text)
                elif ref_element.tag==NS + 'source':
                    for source_element in ref_element.getchildren():
                        if source_element.tag==NS + 'tissue':
                            tissues.append(source_element.text)
            if scopes:
                scopes_str='Scope: '+', '.join(scopes)
            else:
                scopes_str=''
            if tissues:
                tissues_str='Tissue: '+', '.join(tissues)
            else:
                tissues_str=''
            
            reference.location = [] #locations cannot be parsed since they are actually written in free text inside scopes so all the references are put in the annotation.
            reference.authors = ', '.join(authors) 
            if journal_name:
                if pub_date and j_volume and j_first and j_last:
                    reference.journal = REFERENCE_JOURNAL % dict(name=journal_name,
                        volume=j_volume, first=j_first, last=j_last, pub_date=pub_date)
                else:
                    reference.journal = journal_name 
            reference.comment = ' | '.join((pub_type,pub_date,scopes_str,tissues_str))
            append_to_annotations('references', reference)
        
        def _parse_position(element, offset=0):
            try:
                position=int(element.attrib['position']) + offset
            except KeyError, err:
                position=None
            status = element.attrib.get('status', '')
            if status == 'unknown':
                assert position is None
                return SeqFeature.UnknownPosition()
            elif not status:
                return SeqFeature.ExactPosition(position)
            elif status == 'greater than':
                return SeqFeature.AfterPosition(position)
            elif status == 'less than':
                return SeqFeature.BeforePosition(position)
            elif status == 'uncertain':
                return SeqFeature.UncertainPosition(position)
            else:
                raise NotImplementedError("Position status %r" % status)

        def _parse_feature(element):
            feature=SeqFeature.SeqFeature()
            for k,v in element.attrib.items():
                feature.qualifiers[k]=v
            feature.type=element.attrib.get('type','')
            if 'id' in element.attrib:
                feature.id=element.attrib['id']
            for feature_element in element.getchildren():
                if feature_element.tag==NS + 'location':
                    position_elements=feature_element.findall(NS + 'position')
                    if position_elements:
                        element = position_elements[0]
                        start_position = _parse_position(element, -1)
                        end_position = _parse_position(element)
                    else:
                        element = feature_element.findall(NS + 'begin')[0]
                        start_position=_parse_position(element, -1)
                        element = feature_element.findall(NS + 'end')[0]
                        end_position=_parse_position(element)
                    feature.location=SeqFeature.FeatureLocation(start_position,end_position)
                else:
                    try:
                        feature.qualifiers[feature_element.tag.replace(NS,'')]=feature_element.text
                    except:
                        pass#skip unparsable tag
            self.ParsedSeqRecord.features.append(feature)
            
        def _parse_proteinExistence(element):
            append_to_annotations('proteinExistence', element.attrib['type'])   
            
        def _parse_evidence(element):
            for k, v in  element.attrib.items():
                ann_key = k
                append_to_annotations(ann_key, v)   
        
        def  _parse_sequence(element):
            for k, v in element.attrib.items():
                if k in ("length", "mass", "version"):
                    self.ParsedSeqRecord.annotations['sequence_%s' % k] = int(v)
                else:
                    self.ParsedSeqRecord.annotations['sequence_%s' % k] = v
            seq=''.join((element.text.split()))
            self.ParsedSeqRecord.seq=Seq.Seq(seq,self.alphabet)
            
        #============================================#
        #Initialize SeqRecord
        self.ParsedSeqRecord=SeqRecord('', id='') 
        
        #Entry attribs parsing
        #Unknown dataset should not happen!
        self.dbname=self.entry.attrib.get('dataset', 'UnknownDataset')
        #add attribs to annotations
        for k, v in self.entry.attrib.items():
            if k in ("version"):
                #original
                #self.ParsedSeqRecord.annotations["entry_%s" % k] = int(v)
                #To cope with swissProt plain text parser. this can cause errors 
                #if the attrib has the same name of an other annotation
                self.ParsedSeqRecord.annotations[k] = int(v)
            else:
                #self.ParsedSeqRecord.annotations["entry_%s" % k] = v
                self.ParsedSeqRecord.annotations[k] = v # to cope with swissProt plain text parser

        #Top-to-bottom entry children parsing
        for element in self.entry.getchildren():
            if element.tag==NS + 'name':
                _parse_name(element)
            elif element.tag==NS + 'accession':
                _parse_accession(element)
            elif element.tag==NS + 'protein':
                _parse_protein(element)  
            elif element.tag==NS + 'gene':
                _parse_gene(element)
            elif element.tag==NS + 'geneLocation':
                _parse_geneLocation(element)
            elif element.tag==NS + 'organism':
                _parse_organism(element)          
            elif element.tag==NS + 'organismHost':
                _parse_organismHost(element)
            elif element.tag==NS + 'keyword':
                _parse_keyword(element)
            elif element.tag==NS + 'comment':
                _parse_comment(element)
            elif element.tag==NS + 'dbReference':
                _parse_dbReference(element)
            elif element.tag==NS + 'reference':
                _parse_reference(element)
            elif element.tag==NS + 'feature':
                _parse_feature(element)
            elif element.tag==NS + 'proteinExistence':
                _parse_proteinExistence(element)
            elif element.tag==NS + 'evidence':
                _parse_evidence(element)
            elif element.tag==NS + 'sequence':
                _parse_sequence(element)
            else:
                pass   
            
        self.ParsedSeqRecord.dbxrefs=list(set(self.ParsedSeqRecord.dbxrefs))#remove duplicate dbxrefs
        self.ParsedSeqRecord.dbxrefs.sort()

        # use first accession as id
        if not self.ParsedSeqRecord.id:
            self.ParsedSeqRecord.id=self.ParsedSeqRecord.annotations['accessions'][0]
        
        return self.ParsedSeqRecord
        
