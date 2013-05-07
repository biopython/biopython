#!/usr/bin/env python
import sys
###########################################
"""
Parsers for the GAF, GPA and GPI formats
from UniProt-GOA.

Uniprot-GOA README:
ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/README

gp_association (GPA format) README:
ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gp_association_readme

gp_information (GPI format) README:
ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gp_information_readme

(c) 2013 Iddo Friedberg idoerg@gmail.com
http://iddo-friedberg.net
Distributed under Biopython license.
"""
###########################################

###########################################
# GAF: GO Annotation Format
#
# GAF version 2.0
GAF20FIELDS = ['DB' , 
        'DB_Object_ID' , 
        'DB_Object_Symbol' , 
        'Qualifier' , 
        'GO_ID' , 
        'DB:Reference' , 
        'Evidence' , 
        'With' , 
        'Aspect',
        'DB_Object_Name' , 
        'Synonym' , 
        'DB_Object_Type' , 
        'Taxon_ID' , 
        'Date' , 
        'Assigned_By' , 
        'Annotation_Extension' , 
        'Gene_Product_Form_ID']
# GAF version 1.0
GAF10FIELDS = ['DB' , 
        'DB_Object_ID' , 
        'DB_Object_Symbol' , 
        'Qualifier' , 
        'GO_ID' , 
        'DB:Reference' , 
        'Evidence' , 
        'With' , 
        'Aspect',
        'DB_Object_Name' , 
        'Synonym' , 
        'DB_Object_Type' , 
        'Taxon_ID' , 
        'Date' , 
        'Assigned_By'] 


# GPA version 1.0
GPA10FIELDS = [
      'DB',
      'DB_Object_ID',
      'Qualifier',
      'GO_ID',
      'DB:Reference',
      'Evidence code',
      'With',
      'Interacting_taxon_ID',
      'Date',
      'Assigned_by',
      'Annotation_Extension',
      'Spliceform_ID']

# GPA version 1.1
GPA11FIELDS = [
      'DB',
      'DB_Object_ID',
      'Qualifier',
      'GO_ID',
      'DB:Reference',
      'ECO_Evidence_code',
      'With',
      'Interacting_taxon_ID',
      'Date',
      'Assigned_by',
      'Annotation Extension',
      'Annotation_Properties']

# GPI version 1.1
GPI10FIELDS = [
      'DB',
      'DB_subset',
      'DB_Object_ID',
      'DB_Object_Symbol',
      'DB_Object_Name',
      'DB_Object_Synonym',
      'DB_Object_Type',
      'Taxon',
      'Annotation_Target_Set',
      'Annotation_Completed',
      'Parent_Object_ID']

def _gpi10iterator(handle):
    """ This is a generator function
    This iterator is used to read a gp_information.goa_uniprot
    file which is in the GPA 1.0 format."""
    for inline in handle:
        inrec = inline.rstrip('\n').split('\t')
        if len(inrec) == 1:
            continue
        gpi_rec = {}
        gpi_rec['DB'] = inrec[0]
        gpi_rec['DB_subset'] = inrec[1]
        gpi_rec['DB_Object_ID'] = inrec[2]
        gpi_rec['DB_Object_Symbol'] = inrec[3]
        gpi_rec['DB_Object_Name'] = inrec[4]
        gpi_rec['DB_Object_Synonym'] = inrec[5].split('|')
        gpi_rec['DB_Object_Type'] = inrec[6]
        gpi_rec['Taxon'] = inrec[7]
        gpi_rec['Annotation_Target_Set'] = inrec[8].split('|')
        gpi_rec['Annotation_Completed'] = inrec[9]
        gpi_rec['Parent_Object_ID'] = inrec[10]

        yield gpi_rec


def gpi_iterator(handle):
    """ This function should be called to read a
    gp_information.goa_uniprot file. At the moment, there is
    only one format, but this may change, so 
    this function is a placeholder a future wrapper.
    """
    return _gpi10iterator(handle)

def _gpa10iterator(handle):
    """ This is a generator function
    This iterator is used to read a gp_association.*
    file which is in the GPA 1.0 format."""

    for inline in handle:
        inrec = inline.rstrip('\n').split('\t')
        if len(inrec) == 1:
            continue
        gpa_rec = {}
        gpa_rec['DB'] = inrec[0]
        gpa_rec['DB_Object_ID'] = inrec[1]
        gpa_rec['Qualifier'] = inrec[2].split('|')
        gpa_rec['GO_ID'] = inrec[3]
        gpa_rec['DB:Reference'] = inrec[4].split('|')
        gpa_rec['Evidence code'] = inrec[5]
        gpa_rec['With'] = inrec[6].split('|')
        gpa_rec['Interacting_taxon_ID'] = inrec[7],split('|')
        gpa_rec['Date'] = inrec[8]
        gpa_rec['Assigned_by'] = inrec[9]
        gpa_rec['Annotation_Extension'] = inrec[10]
        gpa_rec['Spliceform_ID'] = inrec[11]

        yield gpa_rec

def _gpa11iterator(handle):
    """ This is a generator function
    This iterator is used to read a gp_association.goa_uniprot
    file which is in the GPA 1.1 format."""

    for inline in handle:
        inrec = inline.rstrip('\n').split('\t')
        if len(inrec) == 1:
            continue
        gpa_rec = {}
        gpa_rec['DB'] = inrec[0]
        gpa_rec['DB_Object_ID'] = inrec[1]
        gpa_rec['Qualifier'] = inrec[2]
        gpa_rec['GO_ID'] = inrec[3]
        gpa_rec['DB:Reference'] = inrec[4]
        gpa_rec['ECO_Evidence_code'] = inrec[5]
        gpa_rec['With'] = inrec[6]
        gpa_rec['Interacting_taxon_ID'] = inrec[7]
        gpa_rec['Date'] = inrec[8]
        gpa_rec['Assigned_by'] = inrec[9]
        gpa_rec['Annotation_Extension'] = inrec[10]
        gpa_rec['Spliceform_ID'] = inrec[11]

        gpa_rec['Qualifier'] = gpa_rec['Qualifier'].split('|')
        gpa_rec['DB:Reference'] = gpa_rec['DB:Reference'].split('|')
        gpa_rec['With'] = gpa_rec['With'].split('|')
        gpa_rec['Interacting_taxon_ID'] = \
               gpa_rec['Interacting_taxon_ID'].split('|')
        yield gpa_rec

def gpa_iterator(handle):
    """ This function should be called to read a
    gene_association.goa_uniprot file. Reads the first record and
    returns a gpa 1.1 or a gpa 1.0 iterator as needed
    """

    inline = handle.readline()
    if inline.strip() == '!gpa-version: 1.1':
        sys.stderr.write("gpa 1.1\n")
        return _gpa20iterator(handle)
    else:
        sys.stderr.write("gaf 1.1\n")
        return _gpa11iterator(handle)

def _gaf10iterator(handle):
    for inline in handle:
        inrec = inline.rstrip('\n').split('\t')
        if len(inrec) == 1:
            continue
        inrec[3] = inrec[3].split('|') #Qualifier
        inrec[5] = inrec[5].split('|') # DB:reference(s)
        inrec[7] = inrec[7].split('|') # With || From
        inrec[9] = inrec[9].split('|') # With || From
        inrec[10] = inrec[10].split('|') # Synonym
        inrec[12] = inrec[12].split('|') # Taxon
        yield dict(zip(GAF10FIELDS, inrec))

def _gaf20iterator(handle):
    for inline in handle:
        inrec = inline.rstrip('\n').split('\t')
        if len(inrec) == 1:
            continue
        inrec[3] = inrec[3].split('|') #Qualifier
        inrec[5] = inrec[5].split('|') # DB:reference(s)
        inrec[7] = inrec[7].split('|') # With || From
        inrec[9] = inrec[9].split('|') # With || From
        inrec[10] = inrec[10].split('|') # Synonym
        inrec[12] = inrec[12].split('|') # Taxon
        yield dict(zip(GAF20FIELDS, inrec))

def _gaf10iterator(handle):
    for inline in handle:
        inrec = inline.rstrip('\n').split('\t')
        if len(inrec) == 1:
            continue
        inrec[3] = inrec[3].split('|') #Qualifier
        inrec[5] = inrec[5].split('|') # DB:reference(s)
        inrec[7] = inrec[7].split('|') # With || From
        inrec[9] = inrec[9].split('|') # With || From
        inrec[10] = inrec[10].split('|') # Synonym
        inrec[12] = inrec[12].split('|') # Taxon
        yield dict(zip(GAF10FIELDS, inrec))

#def _gaf10iterator(handle):
#    """ This is a generator function
#    This iterator is used to read a gene_association.goa_uniprot
#    file which is in the GAF 1.0 format. Should not be called directly,
#    but rather through the GAFIterator function
#    This a shorter source code, but the nested enumerate loop makes it slower
#    hence, it is commented out"""
#
#    for inline in handle:
#        inrec = inline.rstrip('\n').split('\t')
#        if len(inrec) == 1:
#            continue
#        gaf_rec = {}
#        for (i, field) in enumerate(GAF10FIELDS):
#            gaf_rec[field] = inrec[i]
#        yield gaf_rec

def _gaf10iterator(handle):
    """ This is a generator function
    This iterator is used to read a gene_association.goa_uniprot
    file which is in the GAF 1.0 format.
    Should not be called directly, but rather through the GAFIterator
    function"""
    for inline in handle:
        inrec = inline.rstrip('\n').split('\t')
        if len(inrec) == 1:
            continue
        gaf_rec = {}
        gaf_rec['DB'] = inrec[0]
        gaf_rec['DB_Object_ID'] = inrec[1]
        gaf_rec['DB_Object_Symbol'] = inrec[2]
        gaf_rec['Qualifier'] = inrec[3]
        gaf_rec['GO_ID'] = inrec[4]
        gaf_rec['DB:Reference'] = inrec[5]
        gaf_rec['Evidence'] = inrec[6]
        gaf_rec['With'] = inrec[7]
        gaf_rec['Aspect']=inrec[8]
        gaf_rec['DB_Object_Name'] = inrec[9]
        gaf_rec['Synonym'] = inrec[10]
        gaf_rec['DB_Object_Type'] = inrec[11]
        gaf_rec['Taxon_ID'] = inrec[12]
        gaf_rec['Date'] = inrec[13]
        gaf_rec['Assigned_By'] = inrec[14]
        yield gaf_rec

def gafiterator(handle):
    """ This function should be called to read a
    gene_association.goa_uniprot file. Reads the first record and
    returns a gaf 2.0 or a gaf 1.0 iterator as needed
    """

    inline = handle.readline()
    if inline.strip() == '!gaf-version: 2.0':
        sys.stderr.write("gaf 2.0\n")
        return _gaf20iterator(handle)
    else:
        sys.stderr.write("gaf 1.0\n")
        return _gaf10iterator(handle)
    
def writerec(outrec,handle,FIELDS=GAF20FIELDS, header=None):
    """Write a single  record to an output stream. 
    Caller should know the  format version. Default: gaf-2.0
    If header has a value, then it is assumed this is the first record,
    a header is written and returned
    """
    if header:
        handle.write("%s\n" % header)
    else:
        outstr = ''
        for i in FIELDS[:-1]:
            outstr += outrec[i] + '\t'
        outstr += outrec[FIELDS[-1]] + '\n'
        handle.write("%s" % outstr)


def record_has(inrec, fieldvals = {}):
    """
    Accepts a record, and a dictionary of field values. The
    format is {'field_name': set([val1, val2])}.
    If any field in the record has  a matching value, the function returns
    True. Otherwise, returns False.
    """
    retval = False
    for field in fieldvals:
        if inrec[field] in fieldvals[field]:
            retval = True
            break
    return retval


if __name__ == '__main__':
    # Example: read a GAF file. Write only S. cerevisiae records, but
    # remove all records with IEA evidence
    banned = {'Evidence': set(['IEA','EXP'])}
    allowed = {'Taxon_ID': set(['taxon:4932'])}
    for inrec in gafiterator(open(sys.argv[1])):
        if (not record_has(inrec, banned)) and \
               record_has(inrec, allowed):
            writerec(inrec, sys.stdout,GAF10FIELDS)
        
        
        
         
