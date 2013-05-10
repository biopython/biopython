
import UniProtGOA as upg
""" Common utilities in parsing UniProtGOA files
"""
def split_to_ontologies(handle):
    """Splits a GAF file into three ontology files
    """
    header = "!gaf-version: 2.0\n"
    out_mfo = open("%s.MFO" % handle.name,'w')
    out_bpo = open("%s.BPO" % handle.name,'w')
    out_cco = open("%s.CCO" % handle.name,'w')
    out_bpo.write(header)
    out_bpo.write(header)
    out_cco.write(header)
    for inrec in upg.gafiterator(handle):
        if inrec['Aspect'] == 'F':
            upg.writerec(inrec,out_mfo)
        elif inrec['Aspect'] == 'P':
            upg.writerec(inrec,out_bpo)
        elif inrec['Aspect'] == 'C':
            upg.writerec(inrec,out_cco)
        else:
            raise ValueError, 'unknown ontology aspect %s' % inrec['Aspect']
    out_mfo.close()
    out_bpo.close()
    out_cco.close()

def species_stats(handle):
    """Statistics for species distributions in a gaf file"""
    taxa_count = {}
    for prot_rec in gafbyproteiniterator(handle):
        taxa_count[prot_rec[0]['Taxon_ID']] += 1
    return taxa_count
        
        
def exclusive_IEA(goa_reclist):
    f_only = [] # Molecular Function
    p_only = [] # Biological Process
    c_only = [] # Cellular Component

    for rec in goa_reclist:
        if rec['Evidence'] != 'IEA':
            continue

