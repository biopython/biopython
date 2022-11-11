# Copyright (C) 2011, 2016 by Brandon Invergo (b.invergo@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Methods for parsing codeml results files."""

import re

line_floats_re = re.compile(r"-*\d+\.\d+")


def parse_basics(lines, results):
    """Parse the basic information that should be present in most codeml output files."""
    # multi_models is used to designate there being results for more than one
    # model in the file
    multi_models = False
    multi_genes = False
    version_re = re.compile(r".+ \(in paml version (\d+\.\d+[a-z]*).*")
    model_re = re.compile(r"Model:\s+(.+)")
    num_genes_re = re.compile(r"\(([0-9]+) genes: separate data\)")
    # In codeml 4.1, the codon substitution model is headed by:
    # "Codon frequencies:"
    # In codeml 4.3+ it is headed by:
    # "Codon frequency model:"
    codon_freq_re = re.compile(r"Codon frequenc[a-z\s]{3,7}:\s+(.+)")
    siteclass_re = re.compile(r"Site-class models:\s*([^\s]*)")
    for line in lines:
        # Find all floating point numbers in this line
        line_floats_res = line_floats_re.findall(line)
        line_floats = [float(val) for val in line_floats_res]
        # Get the program version number
        version_res = version_re.match(line)
        if version_res is not None:
            results["version"] = version_res.group(1)
            continue
        # Find the model description at the beginning of the file.
        model_res = model_re.match(line)
        if model_res is not None:
            results["model"] = model_res.group(1)
        # Find out if more than one genes are analyzed
        num_genes_res = num_genes_re.search(line)
        if num_genes_res is not None:
            results["genes"] = []
            num_genes = int(num_genes_res.group(1))
            for n in range(num_genes):
                results["genes"].append({})
            multi_genes = True
            continue
        # Get the codon substitution model
        codon_freq_res = codon_freq_re.match(line)
        if codon_freq_res is not None:
            results["codon model"] = codon_freq_res.group(1)
            continue
        # Find the site-class model name at the beginning of the file.
        # This exists only if a NSsites class other than 0 is used.
        # If there's no site class model listed, then there are multiple
        # models in the results file
        # Example match: "Site-class models:  PositiveSelection"
        siteclass_res = siteclass_re.match(line)
        if siteclass_res is not None:
            siteclass_model = siteclass_res.group(1)
            if siteclass_model != "":
                results["site-class model"] = siteclass_model
                multi_models = False
            else:
                multi_models = True
        # Get the maximum log-likelihood
        if "ln Lmax" in line and line_floats:
            results["lnL max"] = line_floats[0]
    return (results, multi_models, multi_genes)


def parse_nssites(lines, results, multi_models, multi_genes):
    """Determine which NSsites models are present and parse them."""
    ns_sites = {}
    model_re = re.compile(r"Model (\d+):\s+(.+)")
    gene_re = re.compile(r"Gene\s+([0-9]+)\s+.+")
    siteclass_model = results.get("site-class model")
    if not multi_models:
        # If there's only one model in the results, find out
        # which one it is and then parse it.
        if siteclass_model is None:
            siteclass_model = "one-ratio"
        current_model = {
            "one-ratio": 0,
            "NearlyNeutral": 1,
            "PositiveSelection": 2,
            "discrete": 3,
            "beta": 7,
            "beta&w>1": 8,
            "M2a_rel": 22,
        }[siteclass_model]
        if multi_genes:
            genes = results["genes"]
            current_gene = None
            gene_start = None
            model_results = None
            for line_num, line in enumerate(lines):
                gene_res = gene_re.match(line)
                if gene_res:
                    if current_gene is not None:
                        assert model_results is not None
                        parse_model(lines[gene_start:line_num], model_results)
                        genes[current_gene - 1] = model_results
                    gene_start = line_num
                    current_gene = int(gene_res.group(1))
                    model_results = {"description": siteclass_model}
            if len(genes[current_gene - 1]) == 0:
                model_results = parse_model(lines[gene_start:], model_results)
                genes[current_gene - 1] = model_results
        else:
            model_results = {"description": siteclass_model}
            model_results = parse_model(lines, model_results)
            ns_sites[current_model] = model_results
    else:
        # If there are multiple models in the results, scan through
        # the file and send each model's text to be parsed individually.
        current_model = None
        model_start = None
        for line_num, line in enumerate(lines):
            # Find model names. If this is found on a line,
            # all of the following lines until the next time this matches
            # contain results for Model X.
            # Example match: "Model 1: NearlyNeutral (2 categories)"
            model_res = model_re.match(line)
            if model_res:
                if current_model is not None:
                    # We've already been tracking a model, so it's time
                    # to send those lines off for parsing before beginning
                    # a new one.
                    parse_model(lines[model_start:line_num], model_results)
                    ns_sites[current_model] = model_results
                model_start = line_num
                current_model = int(model_res.group(1))
                model_results = {"description": model_res.group(2)}
        if ns_sites.get(current_model) is None:
            # When we reach the end of the file, we'll still have one more
            # model to parse.
            model_results = parse_model(lines[model_start:], model_results)
            ns_sites[current_model] = model_results
    # Only add the ns_sites dict to the results if we really have results.
    # Model M0 is added by default in some cases, so if it exists, make sure
    # it's not empty
    if len(ns_sites) == 1:
        m0 = ns_sites.get(0)
        if not m0 or len(m0) > 1:
            results["NSsites"] = ns_sites
    elif len(ns_sites) > 1:
        results["NSsites"] = ns_sites
    return results


def parse_model(lines, results):
    """Parse an individual NSsites model's results."""
    parameters = {}
    SEs_flag = False
    dS_tree_flag = False
    dN_tree_flag = False
    w_tree_flag = False
    num_params = None
    tree_re = re.compile(r"^\([\w #:',.()]*\);\s*$")
    branch_re = re.compile(r"\s+(\d+\.\.\d+)[\s+\d+\.\d+]+")
    model_params_re = re.compile(r"(?<!\S)([a-z]\d?)\s*=\s+(\d+\.\d+)")
    for line in lines:
        # Find all floating point numbers in this line
        line_floats_res = line_floats_re.findall(line)
        line_floats = [float(val) for val in line_floats_res]
        # Check if branch-specific results are in the line
        branch_res = branch_re.match(line)
        # Check if additional model parameters are in the line
        model_params = model_params_re.findall(line)
        # Find lnL values.
        # Example match (lnL = -2021.348300):
        # "lnL(ntime: 19  np: 22):  -2021.348300      +0.000000"
        if "lnL(ntime:" in line and line_floats:
            results["lnL"] = line_floats[0]
            np_res = re.match(r"lnL\(ntime:\s+\d+\s+np:\s+(\d+)\)", line)
            if np_res is not None:
                num_params = int(np_res.group(1))
        # Get parameter list. This can be useful for specifying starting
        # parameters in another run by copying the list of parameters
        # to a file called in.codeml. Since the parameters must be in
        # a fixed order and format, copying and pasting to the file is
        # best. For this reason, they are grabbed here just as a long
        # string and not as individual numbers.
        elif len(line_floats) == num_params and not SEs_flag:
            parameters["parameter list"] = line.strip()
        # Find SEs. The same format as parameters above is maintained
        # since there is a correspondence between the SE format and
        # the parameter format.
        # Example match:
        # "SEs for parameters:
        # -1.00000 -1.00000 -1.00000 801727.63247 730462.67590 -1.00000
        elif "SEs for parameters:" in line:
            SEs_flag = True
        elif SEs_flag and len(line_floats) == num_params:
            parameters["SEs"] = line.strip()
            SEs_flag = False
        # Find tree lengths.
        # Example match: "tree length =   1.71931"
        elif "tree length =" in line and line_floats:
            results["tree length"] = line_floats[0]
        # Find the estimated trees only taking the tree if it has
        # lengths or rate estimates on the branches
        elif tree_re.match(line) is not None:
            if ":" in line or "#" in line:
                if dS_tree_flag:
                    results["dS tree"] = line.strip()
                    dS_tree_flag = False
                elif dN_tree_flag:
                    results["dN tree"] = line.strip()
                    dN_tree_flag = False
                elif w_tree_flag:
                    results["omega tree"] = line.strip()
                    w_tree_flag = False
                else:
                    results["tree"] = line.strip()
        elif "dS tree:" in line:
            dS_tree_flag = True
        elif "dN tree:" in line:
            dN_tree_flag = True
        elif "w ratios as labels for TreeView:" in line:
            w_tree_flag = True
        # Find rates for multiple genes
        # Example match: "rates for 2 genes:     1  2.75551"
        elif "rates for" in line and line_floats:
            line_floats.insert(0, 1.0)
            parameters["rates"] = line_floats
        # Find kappa values.
        # Example match: "kappa (ts/tv) =  2.77541"
        elif "kappa (ts/tv)" in line and line_floats:
            parameters["kappa"] = line_floats[0]
        # Find omega values.
        # Example match: "omega (dN/dS) =  0.25122"
        elif "omega (dN/dS)" in line and line_floats:
            parameters["omega"] = line_floats[0]
        elif "w (dN/dS)" in line and line_floats:
            parameters["omega"] = line_floats
        # Find omega and kappa values for multi-gene files
        # Example match: "gene # 1: kappa =   1.72615 omega =   0.39333"
        elif "gene # " in line:
            gene_num = int(re.match(r"gene # (\d+)", line).group(1))
            if parameters.get("genes") is None:
                parameters["genes"] = {}
            parameters["genes"][gene_num] = {
                "kappa": line_floats[0],
                "omega": line_floats[1],
            }
        # Find dN values.
        # Example match: "tree length for dN:       0.2990"
        elif "tree length for dN" in line and line_floats:
            parameters["dN"] = line_floats[0]
        # Find dS values
        # Example match: "tree length for dS:       1.1901"
        elif "tree length for dS" in line and line_floats:
            parameters["dS"] = line_floats[0]
        # Find site class distributions.
        # Example match 1 (normal model, 2 site classes):
        # "p:   0.77615  0.22385"
        # Example match 2 (branch site A model, 4 site classes):
        # "proportion       0.00000  0.00000  0.73921  0.26079"
        elif line[0:2] == "p:" or line[0:10] == "proportion":
            site_classes = parse_siteclass_proportions(line_floats)
            parameters["site classes"] = site_classes
        # Find the omega value corresponding to each site class
        # Example match (2 site classes): "w:   0.10224  1.00000"
        elif line[0:2] == "w:":
            site_classes = parameters.get("site classes")
            site_classes = parse_siteclass_omegas(line, site_classes)
            parameters["site classes"] = site_classes
        # Find the omega values corresponding to a branch type from
        # the clade model C for each site class
        # Example match:
        # "branch type 0:    0.31022   1.00000   0.00000"
        elif "branch type " in line:
            branch_type = re.match(r"branch type (\d)", line)
            if branch_type:
                site_classes = parameters.get("site classes")
                branch_type_no = int(branch_type.group(1))
                site_classes = parse_clademodelc(
                    branch_type_no, line_floats, site_classes
                )
                parameters["site classes"] = site_classes
        # Find the omega values of the foreground branch for each site
        # class in the branch site A model
        # Example match:
        # "foreground w     0.07992  1.00000 134.54218 134.54218"
        elif line[0:12] == "foreground w":
            site_classes = parameters.get("site classes")
            site_classes = parse_branch_site_a(True, line_floats, site_classes)
            parameters["site classes"] = site_classes
        # Find the omega values of the background for each site
        # class in the branch site A model
        # Example match:
        # "background w     0.07992  1.00000  0.07992  1.00000"
        elif line[0:12] == "background w":
            site_classes = parameters.get("site classes")
            site_classes = parse_branch_site_a(False, line_floats, site_classes)
            parameters["site classes"] = site_classes
        # Find dN & dS for each branch, which is organized in a table
        # The possibility of NaNs forces me to not use the line_floats
        # method.
        # Example row (some spaces removed to make it smaller...).
        # " 6..7   0.000  167.7  54.3  0.0000  0.0000  0.0000  0.0  0.0"
        elif branch_res is not None and line_floats:
            branch = branch_res.group(1)
            if parameters.get("branches") is None:
                parameters["branches"] = {}
            params = line.strip().split()[1:]
            parameters["branches"][branch] = {
                "t": float(params[0].strip()),
                "N": float(params[1].strip()),
                "S": float(params[2].strip()),
                "omega": float(params[3].strip()),
                "dN": float(params[4].strip()),
                "dS": float(params[5].strip()),
                "N*dN": float(params[6].strip()),
                "S*dS": float(params[7].strip()),
            }
        # Find model parameters, which can be spread across multiple
        # lines.
        # Example matches:
        # "  p0=  0.99043  p=  0.36657 q=  1.04445
        # "  (p1=  0.00957) w=  3.25530"
        elif model_params:
            float_model_params = []
            for param in model_params:
                float_model_params.append((param[0], float(param[1])))
            parameters.update(dict(float_model_params))
    if parameters:
        results["parameters"] = parameters
    return results


def parse_siteclass_proportions(line_floats):
    """Find proportion of alignment assigned to each class.

    For models which have multiple site classes, find the proportion of the
    alignment assigned to each class.
    """
    site_classes = {}
    if line_floats:
        for n in range(len(line_floats)):
            site_classes[n] = {"proportion": line_floats[n]}
    return site_classes


def parse_siteclass_omegas(line, site_classes):
    """Find omega estimate for each class.

    For models which have multiple site classes, find the omega estimated
    for each class.
    """
    # The omega results are tabular with strictly 9 characters per column
    # (1 to 3 digits before the  decimal point and 5 after). This causes
    # numbers to sometimes run into each other, so we must use a different
    # regular expression to account for this. i.e.:
    # w:   0.00012  1.00000109.87121
    line_floats = re.findall(r"\d{1,3}\.\d{5}", line)
    if not site_classes or len(line_floats) == 0:
        return
    for n in range(len(line_floats)):
        site_classes[n]["omega"] = line_floats[n]
    return site_classes


def parse_clademodelc(branch_type_no, line_floats, site_classes):
    """Parse results specific to the clade model C."""
    if not site_classes or len(line_floats) == 0:
        return
    for n in range(len(line_floats)):
        if site_classes[n].get("branch types") is None:
            site_classes[n]["branch types"] = {}
        site_classes[n]["branch types"][branch_type_no] = line_floats[n]
    return site_classes


def parse_branch_site_a(foreground, line_floats, site_classes):
    """Parse results specific to the branch site A model."""
    if not site_classes or len(line_floats) == 0:
        return
    for n in range(len(line_floats)):
        if site_classes[n].get("branch types") is None:
            site_classes[n]["branch types"] = {}
        if foreground:
            site_classes[n]["branch types"]["foreground"] = line_floats[n]
        else:
            site_classes[n]["branch types"]["background"] = line_floats[n]
    return site_classes


def parse_pairwise(lines, results):
    """Parse results from pairwise comparisons."""
    # Find pairwise comparisons
    # Example:
    # 2 (Pan_troglo) ... 1 (Homo_sapie)
    # lnL = -291.465693
    #  0.01262 999.00000  0.00100
    #
    # t= 0.0126  S=    81.4  N=   140.6  dN/dS= 0.0010  dN= 0.0000  dS= 0.0115
    pair_re = re.compile(r"\d+ \((.+)\) ... \d+ \((.+)\)")
    pairwise = {}
    seq1 = None
    seq2 = None
    for line in lines:
        # Find all floating point numbers in this line
        line_floats_res = line_floats_re.findall(line)
        line_floats = [float(val) for val in line_floats_res]
        pair_res = pair_re.match(line)
        if pair_res:
            seq1 = pair_res.group(1)
            seq2 = pair_res.group(2)
            if seq1 not in pairwise:
                pairwise[seq1] = {}
            if seq2 not in pairwise:
                pairwise[seq2] = {}
        if len(line_floats) == 1 and seq1 is not None and seq2 is not None:
            pairwise[seq1][seq2] = {"lnL": line_floats[0]}
            pairwise[seq2][seq1] = pairwise[seq1][seq2]
        elif len(line_floats) == 6 and seq1 is not None and seq2 is not None:
            pairwise[seq1][seq2].update(
                {
                    "t": line_floats[0],
                    "S": line_floats[1],
                    "N": line_floats[2],
                    "omega": line_floats[3],
                    "dN": line_floats[4],
                    "dS": line_floats[5],
                }
            )
            pairwise[seq2][seq1] = pairwise[seq1][seq2]
    if pairwise:
        results["pairwise"] = pairwise
    return results


def parse_distances(lines, results):
    """Parse amino acid sequence distance results."""
    distances = {}
    sequences = []
    raw_aa_distances_flag = False
    ml_aa_distances_flag = False
    matrix_row_re = re.compile(r"(.+)\s{5,15}")
    for line in lines:
        # Find all floating point numbers in this line
        line_floats_res = line_floats_re.findall(line)
        line_floats = [float(val) for val in line_floats_res]
        if "AA distances" in line:
            raw_aa_distances_flag = True
            # In current versions, the raw distances always come
            # first but I don't trust this to always be true
            ml_aa_distances_flag = False
        elif "ML distances of aa seqs." in line:
            ml_aa_distances_flag = True
            raw_aa_distances_flag = False
        # Parse AA distances (raw or ML), in a lower diagonal matrix
        matrix_row_res = matrix_row_re.match(line)
        if matrix_row_res and (raw_aa_distances_flag or ml_aa_distances_flag):
            seq_name = matrix_row_res.group(1).strip()
            if seq_name not in sequences:
                sequences.append(seq_name)
            if raw_aa_distances_flag:
                if distances.get("raw") is None:
                    distances["raw"] = {}
                distances["raw"][seq_name] = {}
                for i in range(0, len(line_floats)):
                    distances["raw"][seq_name][sequences[i]] = line_floats[i]
                    distances["raw"][sequences[i]][seq_name] = line_floats[i]
            else:
                if distances.get("ml") is None:
                    distances["ml"] = {}
                distances["ml"][seq_name] = {}
                for i in range(0, len(line_floats)):
                    distances["ml"][seq_name][sequences[i]] = line_floats[i]
                    distances["ml"][sequences[i]][seq_name] = line_floats[i]
    if distances:
        results["distances"] = distances
    return results
