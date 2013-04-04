# Copyright 2013 by Leighton Pritchard.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

""" This module provides a function to scrape specific KGML files from KEGG
"""

import urllib2
from KGML_parser import read

# List of 2010 IDs for metabolic pathways
metabolic = ["ko00010", "ko00020", "ko00030", "ko00040", "ko00051", "ko00052", 
             "ko00053", "ko00061", "ko00062", "ko00071", "ko00072", "ko00100", 
             "ko00120", "ko00121", "ko00130", "ko00140", "ko00190", "ko00195", 
             "ko00196", "ko00230", "ko00231", "ko00232", "ko00240", "ko00250", 
             "ko00253", "ko00260", "ko00270", "ko00280", "ko00281", "ko00290", 
             "ko00300", "ko00310", "ko00311", "ko00312", "ko00330", "ko00331", 
             "ko00340", "ko00350", "ko00351", "ko00360", "ko00361", "ko00362", 
             "ko00363", "ko00364", "ko00380", "ko00400", "ko00401", "ko00402", 
             "ko00410", "ko00430", "ko00440", "ko00450", "ko00460", "ko00471", 
             "ko00472", "ko00473", "ko00480", "ko00500", "ko00510", "ko00511", 
             "ko00512", "ko00513", "ko00514", "ko00520", "ko00521", "ko00522", 
             "ko00523", "ko00524", "ko00531", "ko00532", "ko00533", "ko00534", 
             "ko00540", "ko00550", "ko00561", "ko00562", "ko00563", "ko00564", 
             "ko00565", "ko00590", "ko00591", "ko00592", "ko00600", "ko00601", 
             "ko00603", "ko00604", "ko00620", "ko00621", "ko00622", "ko00623", 
             "ko00624", "ko00625", "ko00626", "ko00627", "ko00630", "ko00633", 
             "ko00640", "ko00642", "ko00643", "ko00650", "ko00660", "ko00670", 
             "ko00680", "ko00710", "ko00720", "ko00730", "ko00740", "ko00750", 
             "ko00760", "ko00770", "ko00780", "ko00785", "ko00790", "ko00791", 
             "ko00830", "ko00860", "ko00900", "ko00901", "ko00902", "ko00903", 
             "ko00904", "ko00905", "ko00906", "ko00908", "ko00909", "ko00910", 
             "ko00920", "ko00930", "ko00940", "ko00941", "ko00942", "ko00943", 
             "ko00944", "ko00945", "ko00950", "ko00960", "ko00965", "ko00966", 
             "ko00970", "ko00980", "ko00981", "ko00982", "ko00983", "ko01040", 
             "ko01051", "ko01053", "ko01055", "ko01056", "ko01057", "ko01058", 
             "ko01100", "ko01110", "ko01120", "ko04070"]

# List of 2010 IDs for non-metabolic pathways
non_metabolic = ["ko02010", "ko02020", "ko02030", "ko02040", "ko02060", 
                 "ko03008", "ko03010", "ko03013", "ko03015", "ko03018", 
                 "ko03020", "ko03022", "ko03030", "ko03040", "ko03050", 
                 "ko03060", "ko03070", "ko03320", "ko03410", "ko03420", 
                 "ko03430", "ko03440", "ko03450", "ko04010", "ko04011", 
                 "ko04012", "ko04013", "ko04020", "ko04060", "ko04062", 
                 "ko04070", "ko04075", "ko04080", "ko04110", "ko04111", 
                 "ko04112", "ko04113", "ko04114", "ko04115", "ko04120", 
                 "ko04122", "ko04130", "ko04140", "ko04141", "ko04142", 
                 "ko04144", "ko04145", "ko04146", "ko04150", "ko04210", 
                 "ko04260", "ko04270", "ko04310", "ko04320", "ko04330", 
                 "ko04340", "ko04350", "ko04360", "ko04370", "ko04380", 
                 "ko04510", "ko04512", "ko04514", "ko04520", "ko04530", 
                 "ko04540", "ko04610", "ko04612", "ko04614", "ko04620", 
                 "ko04621", "ko04622", "ko04623", "ko04626", "ko04630", 
                 "ko04640", "ko04650", "ko04660", "ko04662", "ko04664", 
                 "ko04666", "ko04670", "ko04672", "ko04710", "ko04711", 
                 "ko04712", "ko04720", "ko04722", "ko04730", "ko04740", 
                 "ko04742", "ko04744", "ko04745", "ko04810", "ko04910", 
                 "ko04912", "ko04914", "ko04916", "ko04920", "ko04930", 
                 "ko04940", "ko04950", "ko04960", "ko04961", "ko04962", 
                 "ko04964", "ko04966", "ko04970", "ko04971", "ko04972", 
                 "ko04973", "ko04974", "ko04975", "ko04976", "ko04977", 
                 "ko04978", "ko05010", "ko05012", "ko05014", "ko05016", 
                 "ko05020", "ko05100", "ko05110", "ko05111", "ko05120", 
                 "ko05130", "ko05131", "ko05140", "ko05142", "ko05143", 
                 "ko05144", "ko05145", "ko05146", "ko05150", "ko05152", 
                 "ko05160", "ko05162", "ko05200", "ko05210", "ko05211", 
                 "ko05212", "ko05213", "ko05214", "ko05215", "ko05216", 
                 "ko05217", "ko05218", "ko05219", "ko05220", "ko05221", 
                 "ko05222", "ko05223", "ko05310", "ko05320", "ko05322", 
                 "ko05323", "ko05330", "ko05332", "ko05340", "ko05410", 
                 "ko05412", "ko05414", "ko05416"]

# This string is the template for a KGML request via the web API
url_template = 'http://www.genome.jp/kegg-bin/download?entry=%s&format=kgml'

def retrieve_kgml(map_id):
    """ Returns the raw KGML response from KEGG for the passed KEGG map ID
    """
    f = urllib2.urlopen(url_template % map_id)
    return f.read()

def retrieve_kgml_stream(map_id):
    """ Returns a stream containing the raw KGML response from KEGG for 
        the passed KEGG map ID
    """
    f = urllib2.urlopen(url_template % map_id)
    return f


def retrieve_kgml_to_file(map_id, filename):
    """ Retrieves raw KGML from KEGG for the passed KEGG map ID, and writes it
        to a new file
    """
    outfh = open(filename, 'w')
    outfh.write(retrieve_kgml(map_id))
    outfh.close()

def retrieve_KEGG_pathway(map_id):
    """ Returns a KEGGPathway object, downloaded from KEGG, for the passed
        KEGG map ID
    """
    return read(retrieve_kgml_stream(map_id))


if __name__ == '__main__':
    retrieve_kgml_to_file('ddc00190', 'test_retrieve_ddc00190.kgml')

