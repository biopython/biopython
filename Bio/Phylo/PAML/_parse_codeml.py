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
                multi_models = False
            else:
                multi_models = True
        if "ln Lmax" in line and len(line_floats) > 0:
            results["lnL max"] = line_floats[0]
    return multi_models
 
def parse_nssites(lines, results, multi_models):
    """Determine which NSsites models are present and parse them."""

    results["NSsites"] = {}
    model_re = re.compile("Model (\d+):\s+(.+)")
    siteclass_model = results.get("site-class model")
    if not multi_models:
    # If there's only one model in the results, find out
    # which one it is and then parse it. 
        if siteclass_model is None:
            siteclass_model = "one-ratio"
        current_model = {"one-ratio" : 0,
                        "NearlyNeutral" : 1,
                        "PositiveSelection" : 2,
                        "discrete (4 categories)" : 3,
                        "beta (4 categories)" : 7,
                        "beta&w>1 (5 categories)" : 8}[siteclass_model]
        model_results = {"description" : siteclass_model}
        parse_model(lines, model_results)
        results["NSsites"][current_model] = model_results
    else:
    # If there are multiple models in the results, scan through
    # the file and send each model's text to be parsed individually.
        current_model = None
        model_start = None
        for line_num in range(len(lines)):
            # Find model names. If this is found on a line,
            # all of the following lines until the next time this matches
            # contain results for Model X.
            # Example match: "Model 1: NearlyNeutral (2 categories)"
            model_res = model_re.match(lines[line_num])
            if model_res is not None:
                if current_model:
                # We've already been tracking a model, so it's time
                # to send those lines off for parsing before beginning
                # a new one.
                    parse_model(lines[model_start:line_num], model_results)
                    results["NSsites"][current_model] = model_results
                model_start = line_num
                current_model = int(model_re.group(1))
                model_results = {"description":model_re.group(2)}
        if results["NSsites"].get(current_model) is None:
        # When we reach the end of the file, we'll still have one more
        # model to parse.
            parse_model(lines[model_start:], model_results)
            results["NSsites"][current_model] = model_results

def parse_model(lines, results):
    results["parameters"] = {}
    SEs_flag = False
    dS_tree_flag = False
    dN_tree_flag = False
    w_tree_flag = False
    tree_re = re.compile("\(\(+")
    for line in lines:
        # Find all floating point numbers in this line
        line_floats_res = line_floats_re.findall(line)
        line_floats = [float(val) for val in line_floats_res]
        # Find lnL values.
        # Example match (lnL = -2021.348300):
        # "lnL(ntime: 19  np: 22):  -2021.348300      +0.000000"
        if "lnL(ntime:" in line and len(line_floats) > 0:
            results["lnL"] = line_floats[0]
            np_res = re.match("lnL\(ntime:\s+\d+\s+np:\s+(\d+)\)",line)
            if np_res is not None:
                num_params = int(np_res.group(1))
        # Get parameter list. This can be useful for specifying starting
        # parameters in another run by copying the list of parameters
        # to a file called in.codeml. Since the parameters must be in
        # a fixed order and format, copying and pasting to the file is
        # best. For this reason, they are grabbed here just as a long
        # string and not as individual numbers.
        elif len(line_floats) == num_params and not SEs_flag:
            results["parameters"]["parameter list"] = line.strip()
        # Find SEs. The same format as parameters above is maintained
        # since there is a correspondance between the SE format and
        # the parameter format.
        # Example match:
        # "SEs for parameters:
        # -1.00000 -1.00000 -1.00000 801727.63247 730462.67590 -1.00000 
        elif "SEs for parameters:" in line:
            SEs_flag = True
        elif SEs_flag and len(line_floats) == num_params:
            results["parameters"]["SEs"] = line.strip()
            SEs_flag = False
        # Find tree lengths.
        # Example match: "tree length =   1.71931"
        elif "tree length =" in line and len(line_floats) > 0:
            results["tree length"] = line_floats[0]
        # Find the estimated trees only taking the tree if it has
        # lengths or rate estimates on the branches
        elif tree_re.match(line) is not None:
            if ":" in line:
                if dS_tree_flag:
                    results["dS tree"] = line.strip()
                    dS_tree_flag = False
                elif dN_tree_flag:
                    results["dN tree"] = line.strip()
                    dN_tree_flag = False
                elif w_tree_flag:
                    line_edit = line.replace(" '#", ": ")
                    line_edit = line_edit.replace("'", "")
                    line_edit = line_edit.replace(" ,", ",")
                    line_edit = line_edit.replace(" )", ")")
                    results["omega tree"] = line_edit.strip()
                    w_tree_flag = False
                else:
                    results["tree"] = line.strip()
        # Find rates for multiple genes
        # Example match: "rates for 2 genes:     1  2.75551"
        elif "rates for" in line and len(line_floats) > 0:
            line_floats.insert(0, 1.0)
            results["parameters"]["rates"] = line_floats
        # Find kappa values.
        # Example match: "kappa (ts/tv) =  2.77541"
        elif "kappa (ts/tv)" in line and len(line_floats) > 0:
            results["parameters"]["kappa"] = line_floats[0]
        # Find omega values.
        # Example match: "omega (dN/dS) =  0.25122"
        elif "omega (dN/dS)" in line and len(line_floats) > 0:
            results["parameters"]["omega"] = line_floats[0]
        elif "w (dN/dS)" in line and len(line_floats) > 0:
            results["parameters"]["omega"] = line_floats
        # Find omega and kappa values for multi-gene files
        # Example match: "gene # 1: kappa =   1.72615 omega =   0.39333"
        elif "gene # " in line:
            gene_num = int(re.match("gene # (\d+)", line).group(1))
            if results["parameters"].get("genes") is None:
                results["parameters"]["genes"] = {}
            results["parameters"]["genes"][gene_num] = {"kappa":line_floats[0],
                                                         "omega":line_floats[1]}
        # Find dN values.
        # Example match: "tree length for dN:       0.2990"
        elif "tree length for dN" in line and len(line_floats) > 0:
            results["parameters"]["dN"] = line_floats[0]
        # Find dS values
        # Example match: "tree length for dS:       1.1901"
        elif "tree length for dS" in line and len(line_floats) > 0:
            results["parameters"]["dS"] = line_floats[0]

 
def parse_pairwise(lines, results):
    pass

def parse_distances(lines, results):
    pass

