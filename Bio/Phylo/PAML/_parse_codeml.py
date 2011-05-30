# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import re

line_floats_re = re.compile("-*\d+\.\d+")

def parse_basics(lines, results):
    version_re = re.compile(".+ \(in paml version (\d+\.\d+[a-z]*).*")
    model_re = re.compile("Model:\s+(.+)")
    codon_freq_re = re.compile("Codon frequency model:\s+(.+)")
    siteclass_re = re.compile("Site-class models:\s*(.*)")
    for line in lines:
        # Find all floating point numbers in this line
        line_floats_res = line_floats_re.findall(line)
        line_floats = [float(val) for val in line_floats_res]
        version_res = version_re.match(line)
        if version_res is not None:
            results["version"] = version_re.group(1)
            continue
        # Find the model description at the beginning of the file.
        model_res = model_re.match(line)
        if model_res is not None:
            results["model"] = model_re.group(1)
            continue
        codon_freq_res = codon_freq_re.match(line)
        if codon_freq_res is not None:
            results["codon model"] = codon_freq_res.group(1)
            continue
        # Find the site-class model name at the beginning of the file.
        # This exists only if a NSsites class other than 0 is used.
        # Example match: "Site-class models:  PositiveSelection"
        siteclass_res = siteclass_re.match("Site-class models:\s*(.*)", line)
        if siteclass_res is not None:
            siteclass_model = siteclass_re.group(1)
            if siteclass_model != "":
                results["site-class model"] = siteclass_model
        if "ln Lmax" in line and len(line_floats) > 0:
            results["lnL max"] = line_floats[0]
 
def parse_nssites(lines, results):
    pass

def parse_pairwise(lines, results):
    pass

def parse_distances(lines, results):
    pass

