from Bio.GenBank import FeatureParser
    
#This is a generator function!
def GenBankIterator(handle) :
    """Breaks up a Genbank file into SeqRecord objects

    Every section from the LOCUS line to the terminating // becomes
    a single SeqRecord with associated annotation and features.
    
    Note that for genomes or chromosomes, there is typically only
    one record."""
    parser_object = FeatureParser(debug_level=0) # Use defaults
    while True :
        try :
            record = parser_object.parse(handle)
            if record is None :
                break
            else :
                yield record
        except StopIteration :
            break
    #We have reached the end of the file.
    return
