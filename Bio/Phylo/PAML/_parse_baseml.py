# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Methods for parsing baseml results files."""

import re


line_floats_re = re.compile(r"-*\d+\.\d+")


def parse_basics(lines, results):
    """Parse the basics that should be present in most baseml results files."""
    version_re = re.compile(r"BASEML \(in paml version (\d+\.\d+[a-z]*).*")
    np_re = re.compile(r"lnL\(ntime:\s+\d+\s+np:\s+(\d+)\)")
    num_params = -1
    for line in lines:
        # Find all floating point numbers in this line
        line_floats_res = line_floats_re.findall(line)
        line_floats = [float(val) for val in line_floats_res]
        # Find the version number
        # Example match:
        # "BASEML (in paml version 4.3, August 2009)  alignment.phylip"
        version_res = version_re.match(line)
        if version_res is not None:
            results["version"] = version_res.group(1)
        # Find max lnL
        # Example match:
        # ln Lmax (unconstrained) = -316.049385
        if "ln Lmax" in line and len(line_floats) == 1:
            results["lnL max"] = line_floats[0]
        # Find lnL values.
        # Example match (lnL = -2021.348300):
        # "lnL(ntime: 19  np: 22):  -2021.348300      +0.000000"
        elif "lnL(ntime:" in line and line_floats:
            results["lnL"] = line_floats[0]
            np_res = np_re.match(line)
            if np_res is not None:
                num_params = int(np_res.group(1))
        # Find tree lengths.
        # Example match: "tree length =   1.71931"
        elif "tree length" in line and len(line_floats) == 1:
            results["tree length"] = line_floats[0]
        # Find the estimated tree, only taking the tree if it has
        # branch lengths
        elif re.match(r"\(+", line) is not None:
            if ":" in line:
                results["tree"] = line.strip()
    return (results, num_params)


def parse_parameters(lines, results, num_params):
    """Parse the various parameters from the file."""
    parameters = {}
    parameters = parse_parameter_list(lines, parameters, num_params)
    parameters = parse_kappas(lines, parameters)
    parameters = parse_rates(lines, parameters)
    parameters = parse_freqs(lines, parameters)
    results["parameters"] = parameters
    return results


def parse_parameter_list(lines, parameters, num_params):
    """Parse the parameters list, which is just an unlabeled list of numeric values."""
    for line_num in range(len(lines)):
        line = lines[line_num]
        # Find all floating point numbers in this line
        line_floats_res = line_floats_re.findall(line)
        line_floats = [float(val) for val in line_floats_res]
        # Get parameter list. This can be useful for specifying starting
        # parameters in another run by copying the list of parameters
        # to a file called in.baseml. Since the parameters must be in
        # a fixed order and format, copying and pasting to the file is
        # best. For this reason, they are grabbed here just as a long
        # string and not as individual numbers.
        if len(line_floats) == num_params:
            parameters["parameter list"] = line.strip()
            # Find SEs. The same format as parameters above is maintained
            # since there is a correspondence between the SE format and
            # the parameter format.
            # Example match:
            # "SEs for parameters:
            # -1.00000 -1.00000 -1.00000 801727.63247 730462.67590 -1.00000
            if "SEs for parameters:" in lines[line_num + 1]:
                SEs_line = lines[line_num + 2]
                parameters["SEs"] = SEs_line.strip()
            break
    return parameters


def parse_kappas(lines, parameters):
    """Parse out the kappa parameters."""
    kappa_found = False
    for line in lines:
        # Find all floating point numbers in this line
        line_floats_res = line_floats_re.findall(line)
        line_floats = [float(val) for val in line_floats_res]
        # Find kappa parameter (F84, HKY85, T92 model)
        # Example match:
        # "Parameters (kappa) in the rate matrix (F84) (Yang 1994 J Mol Evol 39:105-111):
        #    3.00749"
        if "Parameters (kappa)" in line:
            kappa_found = True
        elif kappa_found and line_floats:
            branch_res = re.match(r"\s(\d+\.\.\d+)", line)
            if branch_res is None:
                if len(line_floats) == 1:
                    parameters["kappa"] = line_floats[0]
                else:
                    parameters["kappa"] = line_floats
                kappa_found = False
            else:
                if parameters.get("branches") is None:
                    parameters["branches"] = {}
                branch = branch_res.group(1)
                if line_floats:
                    parameters["branches"][branch] = {
                        "t": line_floats[0],
                        "kappa": line_floats[1],
                        "TS": line_floats[2],
                        "TV": line_floats[3],
                    }
        # Find kappa under REV
        # Example match:
        # kappa under REV: 999.00000 145.76453  0.00001  0.00001  0.00001
        elif "kappa under" in line and line_floats:
            if len(line_floats) == 1:
                parameters["kappa"] = line_floats[0]
            else:
                parameters["kappa"] = line_floats
    return parameters


def parse_rates(lines, parameters):
    """Parse the rate parameters."""
    Q_mat_found = False
    trans_probs_found = False
    for line in lines:
        # Find all floating point numbers in this line
        line_floats_res = line_floats_re.findall(line)
        line_floats = [float(val) for val in line_floats_res]
        # Find rate parameters
        # Example match:
        # "Rate parameters:   999.00000 145.59775  0.00001  0.00001  0.00001"
        if "Rate parameters:" in line and line_floats:
            parameters["rate parameters"] = line_floats
        # Find rates
        # Example match:
        # "rate:   0.90121  0.96051  0.99831  1.03711  1.10287"
        elif "rate: " in line and line_floats:
            parameters["rates"] = line_floats
        # Find Rate matrix Q & average kappa (REV model)
        # Example match:
        # Rate matrix Q, Average Ts/Tv =   3.0308
        #  -2.483179    1.865730    0.617449    0.000000
        #   2.298662   -2.298662    0.000000    0.000000
        #   0.335015    0.000000   -0.338059    0.003044
        #   0.000000    0.000000    0.004241   -0.004241
        elif "matrix Q" in line:
            parameters["Q matrix"] = {"matrix": []}
            if line_floats:
                parameters["Q matrix"]["average Ts/Tv"] = line_floats[0]
            Q_mat_found = True
        elif Q_mat_found and line_floats:
            parameters["Q matrix"]["matrix"].append(line_floats)
            if len(parameters["Q matrix"]["matrix"]) == 4:
                Q_mat_found = False
        # Find alpha (gamma shape parameter for variable rates)
        # Example match: "alpha (gamma, K=5) = 192.47918"
        elif "alpha" in line and line_floats:
            parameters["alpha"] = line_floats[0]
        # Find rho for auto-discrete-gamma model
        elif "rho" in line and line_floats:
            parameters["rho"] = line_floats[0]
        elif "transition probabilities" in line:
            parameters["transition probs."] = []
            trans_probs_found = True
        elif trans_probs_found and line_floats:
            parameters["transition probs."].append(line_floats)
            if len(parameters["transition probs."]) == len(parameters["rates"]):
                trans_probs_found = False
    return parameters


def parse_freqs(lines, parameters):
    """Parse the basepair frequencies."""
    root_re = re.compile(r"Note: node (\d+) is root.")
    branch_freqs_found = False
    base_freqs_found = False
    for line in lines:
        # Find all floating point numbers in this line
        line_floats_res = line_floats_re.findall(line)
        line_floats = [float(val) for val in line_floats_res]
        # Find base frequencies from baseml 4.3
        # Example match:
        # "Base frequencies:   0.20090  0.16306  0.37027  0.26577"
        if "Base frequencies" in line and line_floats:
            base_frequencies = {}
            base_frequencies["T"] = line_floats[0]
            base_frequencies["C"] = line_floats[1]
            base_frequencies["A"] = line_floats[2]
            base_frequencies["G"] = line_floats[3]
            parameters["base frequencies"] = base_frequencies
        # Find base frequencies from baseml 4.1:
        # Example match:
        # "base frequency parameters
        # "  0.20317  0.16768  0.36813  0.26102"
        elif "base frequency parameters" in line:
            base_freqs_found = True
        # baseml 4.4 returns to having the base frequencies on the next line
        # but the heading changed
        elif "Base frequencies" in line and not line_floats:
            base_freqs_found = True
        elif base_freqs_found and line_floats:
            base_frequencies = {}
            base_frequencies["T"] = line_floats[0]
            base_frequencies["C"] = line_floats[1]
            base_frequencies["A"] = line_floats[2]
            base_frequencies["G"] = line_floats[3]
            parameters["base frequencies"] = base_frequencies
            base_freqs_found = False
        # Find frequencies
        # Example match:
        # "freq:   0.90121  0.96051  0.99831  1.03711  1.10287"
        elif "freq: " in line and line_floats:
            parameters["rate frequencies"] = line_floats
        # Find branch-specific frequency parameters
        # Example match (note: I think it's possible to have 4 more
        # values per line, enclosed in brackets, so I'll account for
        # this):
        # (frequency parameters for branches)  [frequencies at nodes] (see Yang & Roberts 1995 fig 1)
        #
        # Node #1  ( 0.25824  0.24176  0.25824  0.24176 )
        # Node #2  ( 0.00000  0.50000  0.00000  0.50000 )
        elif "(frequency parameters for branches)" in line:
            parameters["nodes"] = {}
            branch_freqs_found = True
        elif branch_freqs_found:
            if line_floats:
                node_res = re.match(r"Node \#(\d+)", line)
                node_num = int(node_res.group(1))
                node = {"root": False}
                node["frequency parameters"] = line_floats[:4]
                if len(line_floats) > 4:
                    node["base frequencies"] = {
                        "T": line_floats[4],
                        "C": line_floats[5],
                        "A": line_floats[6],
                        "G": line_floats[7],
                    }
                parameters["nodes"][node_num] = node
            else:
                root_res = root_re.match(line)
                if root_res is not None:
                    root_node = int(root_res.group(1))
                    parameters["nodes"][root_node]["root"] = True
                    branch_freqs_found = False
    return parameters
