# Copyright 2010 by Thomas Schmitt.  
# All rights reserved.
#
# This module is for reading and writing SeqXML format files as
# SeqRecord objects, and is expected to be used via the Bio.SeqIO API.
"""Bio.SeqIO support for the "seqxml" file format, SeqXML.

You are expected to use this module via the Bio.SeqIO functions.

SeqXML is a lightweight XML format which is supposed be an alternative for
FASTA files. For more Information see http://www.seqXML.org and Schmitt et al
(2011), http://dx.doi.org/10.1093/bib/bbr025 
"""

from xml.sax.saxutils import XMLGenerator
from xml.sax.xmlreader import AttributesImpl
from xml.dom import pulldom
from xml.sax import SAXParseException

from Bio import Alphabet
from Bio.Seq import Seq
from Bio.Seq import UnknownSeq
from Bio.SeqRecord import SeqRecord
from Interfaces import SequentialSequenceWriter


class XMLRecordIterator:
    """Base class for building iterators for record style XML formats. 
    
    It is assumed that all information for one record can be found within a
    record element or above. Two types of methods are called when the start
    tag of an element is reached. To receive only the attributes of an
    element before its end tag is reached implement _attr_TAGNAME.
    To get an element and its children as a DOM tree implement _elem_TAGNAME. 
    Everything that is part of the DOM tree will not trigger any further 
    method calls.
    """

    def __init__(self,handle,recordTag,namespace=None):
        """Creating the object and initializing the XML parser."""
        
        self._recordTag=recordTag
        self._namespace=namespace
        self._events=pulldom.parse(handle)
        
        
    def __iter__(self):
        """Iterate over the records in the XML file. 
        Returns the last parsed record.""" 
        
        record = None
        try:
            for event,node in self._events:
                
                if event == "START_ELEMENT" and node.namespaceURI == self._namespace:
                    
                    if node.localName == self._recordTag:
                        #create an empty SeqRecord
                        record = SeqRecord('', id='')

                    #call matching methods with attributes only                    
                    if hasattr(self,"_attr_" + node.localName):
                        getattr(self,"_attr_" + node.localName)(self._attributes(node),record)
    
                    #call matching methods with DOM tree 
                    if hasattr(self,"_elem_" + node.localName):
                        #read the element and all nested elements into a DOM tree
                        self._events.expandNode(node)
                        node.normalize()
                        
                        getattr(self,"_elem_" + node.localName)(node,record)
                    
                elif event == "END_ELEMENT" and node.namespaceURI == self._namespace and node.localName == self._recordTag:
                    yield record
                    
        except SAXParseException, e:
            
            if e.getLineNumber() == 1 and e.getColumnNumber() == 0:
                #empty file
                pass
            else:
                import os
                if e.getLineNumber() == 1 and e.getColumnNumber() == 1 \
                and os.name == "java":
                    #empty file, see http://bugs.jython.org/issue1774
                    pass
                else:
                    raise

    
    def _attributes(self,node):
        """Return the attributes of a DOM node as dictionary."""
        
        return dict( (node.attributes.item(i).name,node.attributes.item(i).value) for i in xrange(node.attributes.length) )
    
    

class SeqXmlIterator(XMLRecordIterator):
    """Breaks seqXML file into SeqRecords.
    
    Assumes valid seqXML please validate beforehand."""
    
    def __init__(self,handle):
        """Create the object."""
        XMLRecordIterator.__init__(self, handle,"entry")
        
        self._source = None
        self._source_version = None
        self._version = None
        self._speciesName = None
        self._ncbiTaxId = None

    def _attr_seqXML(self,attr_dict,record):
        """Parse the document metadata."""
        
        if "source" in attr_dict:
            self._source = attr_dict["source"]
        if "sourceVersion" in attr_dict:
            self._source_version = attr_dict["sourceVersion"]
        if "version" in attr_dict:
            self._version = attr_dict["seqXMLversion"]
        if "ncbiTaxID" in attr_dict:
            self._ncbiTaxId = attr_dict["ncbiTaxID"]
        if "speciesName" in attr_dict:
            self._speciesName = attr_dict["speciesName"]
    
    def _attr_property(self,attr_dict,record):
        """Parse key value pair properties and store them as annotations."""
        
        if "name" not in attr_dict:
            raise ValueError("Malformed property element.")
        
        value = attr_dict.get("value",None)
        
        if attr_dict["name"] not in record.annotations:
            record.annotations[attr_dict["name"]] = value
        elif isinstance(record.annotations[attr_dict["name"]],list):
            record.annotations[attr_dict["name"]].append(value)
        else:
            record.annotations[attr_dict["name"]] = [record.annotations[attr_dict["name"]],value]
            
        
    def _attr_species(self,attr_dict,record):
        """Parse the species information."""
        
        if "name" not in attr_dict or "ncbiTaxID" not in attr_dict:
            raise ValueError("Malformed species element!")
        
        #the keywords for the species annotation are taken from SwissIO   
        record.annotations["organism"] = attr_dict["name"]
        record.annotations["ncbi_taxid"] = attr_dict["ncbiTaxID"]
    
    def _attr_entry(self,attr_dict,record):
        """New entry set id and the optional entry source."""
        
        if "id" not in attr_dict:
            raise ValueError("Malformed entry! Identifier is missing.") 
        
        record.id = attr_dict["id"]
        if "source" in attr_dict:
            record.annotations["source"] = attr_dict["source"]
        elif self._source != None:
            record.annotations["source"] = self._source
            
        #initialize entry with global species definition
        #the keywords for the species annotation are taken from SwissIO   
        if self._ncbiTaxId != None:
            record.annotations["ncbi_taxid"] = self._ncbiTaxId
        if self._speciesName != None:
            record.annotations["organism"] = self._speciesName    


    def _elem_DNAseq(self,node,record):
        """Parse DNA sequence."""
        
        if not (node.hasChildNodes() and len(node.firstChild.data) > 0):
            raise ValueError("Sequence length should be greater than 0.")
            
        record.seq = Seq(node.firstChild.data,Alphabet.generic_dna)
        
        
    def _elem_RNAseq(self,node,record):
        """Parse RNA sequence."""
        
        if not (node.hasChildNodes() and len(node.firstChild.data) > 0):
            raise ValueError("Sequence length should be greater than 0.")
        
        record.seq = Seq(node.firstChild.data,Alphabet.generic_rna)
    
    def _elem_AAseq(self,node,record):
        """Parse protein sequence."""
        
        if not (node.hasChildNodes() and len(node.firstChild.data) > 0):
            raise ValueError("Sequence length should be greater than 0.")
        
        record.seq = Seq(node.firstChild.data,Alphabet.generic_protein)
        
        
    def _elem_description(self,node,record):
        """Parse the description."""
        
        if node.hasChildNodes() and len(node.firstChild.data) > 0:
            record.description = node.firstChild.data
        
    def _attr_DBRef(self,attr_dict,record):
        """Parse a database cross reference"""
        
        if "source" not in attr_dict or "id" not in attr_dict:
            raise ValueError("Invalid DB cross reference.")
        
        if "%s:%s" % (attr_dict["source"],attr_dict["id"]) not in record.dbxrefs:
            record.dbxrefs.append("%s:%s" % (attr_dict["source"],attr_dict["id"]) )



class SeqXmlWriter(SequentialSequenceWriter):
    """Writes SeqRecords into seqXML file.
    
    SeqXML requires the sequence alphabet be explicitly RNA, DNA or protein,
    i.e. an instance or subclass of Bio.Alphapet.RNAAlphabet,
    Bio.Alphapet.DNAAlphabet or Bio.Alphapet.ProteinAlphabet.
    """
    
    def __init__(self, handle,source=None,source_version=None,species=None,ncbiTaxId=None):
        """Create Object and start the xml generator."""
        
        SequentialSequenceWriter.__init__(self, handle)

        self.xml_generator = XMLGenerator(handle, "utf-8")
        self.xml_generator.startDocument()
        self.source = source
        self.source_version = source_version
        self.species = species
        self.ncbiTaxId = ncbiTaxId
            
    def write_header(self):
        """Write root node with document metadata."""
        SequentialSequenceWriter.write_header(self)
        
        attrs = {"xmlns:xsi":"http://www.w3.org/2001/XMLSchema-instance",
                 "xsi:noNamespaceSchemaLocation":"http://www.seqxml.org/0.4/seqxml.xsd",
                 "seqXMLversion":"0.4"}
        
        if self.source != None:
            attrs["source"] = self.source
        if self.source_version != None:
            attrs["sourceVersion"] = self.source_ersion
        if self.species != None:
            if not isinstance(species,basestring):
                raise TypeError("species should be of type string")
            attrs["speciesName"] = self.species
        if self.ncbiTaxId != None:
            if not isinstance(self.ncbiTaxId,(basestring,int)):
                raise TypeError("ncbiTaxID should be of type string or int")
            attrs["ncbiTaxID"] = self.ncbiTaxId
        
        self.xml_generator.startElement("seqXML", AttributesImpl(attrs))
        
    
    def write_record(self, record):
        """Write one record."""
        
        if not record.id or record.id == "<unknown id>":
            raise ValueError("SeqXML requires identifier")
        
        if not isinstance(record.id,basestring):
            raise TypeError("Identifier should be of type string")
        
        attrb = {"id" : record.id}
        
        if "source" in record.annotations and self.source != record.annotations["source"]:
            if not isinstance(record.annotations["source"],basestring):
                raise TypeError("source should be of type string")
            attrb["source"] = record.annotations["source"]
        
        self.xml_generator.startElement("entry", AttributesImpl(attrb))
        self._write_species(record)
        self._write_description(record)
        self._write_seq(record)
        self._write_dbxrefs(record)
        self._write_properties(record)
        self.xml_generator.endElement("entry")
    
    def write_footer(self):
        """Close the root node and finish the XML document."""
        
        SequentialSequenceWriter.write_footer(self)
        
        self.xml_generator.endElement("seqXML")
        self.xml_generator.endDocument()
    
    def _write_species(self,record):
        """Write the species if given."""
        
        if "organism" in record.annotations and "ncbi_taxid" in record.annotations:
            
            if not isinstance(record.annotations["organism"],basestring):
                raise TypeError("organism should be of type string")
            
            if not isinstance(record.annotations["ncbi_taxid"],(basestring,int)):
                raise TypeError("ncbiTaxID should be of type string or int")
            
            #The local species definition is only written if it differs from the global species definition
            if record.annotations["organism"] != self.species or record.annotations["ncbi_taxid"] != self.ncbiTaxId:
            
                attr = { "name" : record.annotations["organism"], "ncbiTaxID" :record.annotations["ncbi_taxid"] }   
                self.xml_generator.startElement("species",AttributesImpl(attr))
                self.xml_generator.endElement("species")
            
            
    def _write_description(self,record):
        """Write the description if given."""
        
        if record.description:
            
            if not isinstance(record.description,basestring):
                raise TypeError("Description should be of type string")
            
            description = record.description
            if description == "<unknown description>":
                description = ""
            
            if len(record.description) > 0:
                self.xml_generator.startElement("description",AttributesImpl( {} ))
                self.xml_generator.characters(description)
                self.xml_generator.endElement("description")
    
    def _write_seq(self,record):
        """Write the sequence.
        
        Note that SeqXML requires a DNA, RNA or protein alphabet.
        """
        
        if isinstance(record.seq,UnknownSeq):
            raise TypeError("Sequence type is UnknownSeq but SeqXML requires sequence")
        
        seq = record.seq.tostring()
        
        if not len(seq) > 0:
            raise ValueError("The sequence length should be greater than 0")
        
        #Get the base alphabet (underneath any Gapped or StopCodon encoding)
        alphabet = Alphabet._get_base_alphabet(record.seq.alphabet)
        if isinstance(alphabet,Alphabet.RNAAlphabet):
            seqElem = "RNAseq"
        elif isinstance(alphabet,Alphabet.DNAAlphabet):
            seqElem = "DNAseq"
        elif isinstance(alphabet,Alphabet.ProteinAlphabet):
            seqElem = "AAseq"
        else:
            raise ValueError("Need a DNA, RNA or Protein alphabet")
        
        self.xml_generator.startElement(seqElem,AttributesImpl( {} ))
        self.xml_generator.characters(seq)
        self.xml_generator.endElement(seqElem)    
        
    
    def _write_dbxrefs(self,record):
        """Write all database cross references."""
        if record.dbxrefs != None:
            
            for dbxref in record.dbxrefs:
                
                if not isinstance(dbxref,basestring):
                    raise TypeError("dbxrefs should be of type list of string")
                if dbxref.find(':') < 1:
                    raise ValueError("dbxrefs should be in the form ['source:id', 'source:id' ]")
                
                dbsource,dbid = dbxref.split(':',1)
                
                attr = { "source" : dbsource, "id" : dbid }
                self.xml_generator.startElement("DBRef",AttributesImpl(attr))
                self.xml_generator.endElement("DBRef")
    
    
    def _write_properties(self,record):
        """Write all annotations that are key value pairs with values of a primitive type or list of primitive types."""
        
        for key,value in record.annotations.items():
    
            if key not in ("organism","ncbi_taxid","source"):
            
                if value == None:
                    
                    attr = { "name" : key }
                    self.xml_generator.startElement("property",AttributesImpl(attr))
                    self.xml_generator.endElement("property")
                
                elif isinstance(value,list):
                    
                    for v in value:
                        if isinstance(value,(int,float,basestring)):
                            attr = { "name" : key , "value" : v }
                            self.xml_generator.startElement("property",AttributesImpl(attr))
                            self.xml_generator.endElement("property")
                    
                elif isinstance(value,(int,float,basestring)):
                
                    attr = { "name" : key , "value" : str(value) }
                    self.xml_generator.startElement("property",AttributesImpl(attr))
                    self.xml_generator.endElement("property")
    
if __name__ == "__main__":
    print "Running quick self test"
    
    from Bio import SeqIO
    import sys
    
    fileHandle = open("Tests/SeqXML/protein_example.xml","r")
    records = list(SeqIO.parse(fileHandle, "seqxml"))
    
    from StringIO import StringIO
    stringHandle = StringIO()

    SeqIO.write(records,stringHandle,"seqxml")
    SeqIO.write(records,sys.stdout,"seqxml")
    print
    
    stringHandle.seek(0)
    records = list(SeqIO.parse(stringHandle,"seqxml"))
    
    SeqIO.write(records,sys.stdout,"seqxml")
    print
