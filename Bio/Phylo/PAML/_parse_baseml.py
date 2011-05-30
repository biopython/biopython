# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import re

line_floats_re = re.compile("-*\d+\.\d+")

def parse_basics(lines, results):
    version_re = re.compile("BASEML \(in paml version (\d+\.\d+[a-z]*).*")
    np_re = re.compile("lnL\(ntime:\s+\d+\s+np:\s+(\d+)\)")
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
        elif "lnL(ntime:" in line and len(line_floats) > 0:
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
        elif re.match("\(+", line) is not None:
            if ":" in line:
                results["tree"] = line.strip()
    return num_params

def parse_parameters(lines, results, num_params): 
    results["parameters"] = {}
    parse_parameter_list(lines, results, num_params)
    parse_kappas(lines, results)
    parse_rates(lines, results)
    parse_freqs(lines, results)

def parse_parameter_list(lines, results, num_params):
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
           results["parameters"]["parameter list"] = line.strip()
        # Find SEs. The same format as parameters above is maintained
        # since there is a correspondance between the SE format and
        # the parameter format.
        # Example match:
        # "SEs for parameters:
        # -1.00000 -1.00000 -1.00000 801727.63247 730462.67590 -1.00000 
           if "SEs for parameters:" in lines[line_num + 1]:
                SEs_line = lines[line_num + 2]
                results["parameters"]["SEs"] = SEs_line.strip()
           break

def parse_kappas(lines, results):
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
        elif kappa_found and len(line_floats) > 0:
            branch_res = re.match("\s(\d+\.\.\d+)", line)
            if branch_res is None:
                if len(line_floats) == 1:
                    results["parameters"]["kappa"] = line_floats[0]
                else:
                    results["parameters"]["kappa"] = line_floats
                kappa_found = False
            else:
                if results["parameters"].get("branches") is None:
                    results["parameters"]["branches"] = {}
                branch = branch_res.group(1)
                if len(line_floats) > 0:
                    results["parameters"]["branches"][branch] = \
                        {"t":line_floats[0], "kappa":line_floats[1],
                        "TS":line_floats[2], "TV":line_floats[3]}
        # Find kappa under REV
        # Example match:
        # kappa under REV: 999.00000 145.76453  0.00001  0.00001  0.00001
        elif "kappa under" in line and len(line_floats) > 0:
            if len(line_floats) == 1:
                results["parameters"]["kappa"] = line_floats[0]
            else:
                results["parameters"]["kappa"] = line_floats

def parse_rates(lines, results):
    Q_mat_found = False
    for line in lines:
        # Find all floating point numbers in this line
        line_floats_res = line_floats_re.findall(line)
        line_floats = [float(val) for val in line_floats_res]
        # Find rate parameters
        # Example match: 
        # "Rate parameters:   999.00000 145.59775  0.00001  0.00001  0.00001"
        if "Rate parameters:" in line and len(line_floats) > 0:
            results["parameters"]["rate parameters"] = line_floats
        # Find rates
        # Example match: 
        # "rate:   0.90121  0.96051  0.99831  1.03711  1.10287"
        elif "rate: " in line and len(line_floats) > 0:
            results["parameters"]["rates"] = line_floats
        # Find Rate matrix Q & average kappa (REV model)
        # Example match:
        # Rate matrix Q, Average Ts/Tv =   3.0308
        #  -2.483179    1.865730    0.617449    0.000000
        #   2.298662   -2.298662    0.000000    0.000000
        #   0.335015    0.000000   -0.338059    0.003044
        #   0.000000    0.000000    0.004241   -0.004241
        elif "matrix Q" in line:
            results["parameters"]["Q matrix"] = {"matrix":[]}
            if len(line_floats) > 0:
                results["parameters"]["Q matrix"]["average Ts/Tv"] = \
                    line_floats[0]
            Q_mat_found = True
        elif Q_mat_found and len(line_floats) > 0:
            results["parameters"]["Q matrix"]["matrix"].append(line_floats)
            if len(results["parameters"]["Q matrix"]["matrix"]) == 4:
                Q_mat_found = False
        # Find alpha (gamma shape parameter for variable rates)
        # Example match: "alpha (gamma, K=5) = 192.47918"
        elif "alpha" in line and len(line_floats) > 0:
            results["parameters"]["alpha"] = line_floats[0]

def parse_freqs(lines, results):
    root_re = re.compile("Note: node (\d+) is root.")
    branch_freqs_found = False
    for line in lines:
        # Find all floating point numbers in this line
        line_floats_res = line_floats_re.findall(line)
        line_floats = [float(val) for val in line_floats_res]
        # Find base frequencies
        # Example match:
        # "Base frequencies:   0.20090  0.16306  0.37027  0.26577"  
        if "Base frequencies" in line and len(line_floats) > 0:
            results["parameters"]["base frequencies"] = {}
            results["parameters"]["base frequencies"]["T"] = line_floats[0]
            results["parameters"]["base frequencies"]["C"] = line_floats[1]
            results["parameters"]["base frequencies"]["A"] = line_floats[2]
            results["parameters"]["base frequencies"]["G"] = line_floats[3]
        # Find frequencies
        # Example match: 
        # "freq:   0.90121  0.96051  0.99831  1.03711  1.10287"
        elif "freq: " in line and len(line_floats) > 0:
            results["parameters"]["rate frequencies"] = line_floats
        # Find branch-specific frequency parameters
        # Example match (note: I think it's possible to have 4 more
        # values per line, enclosed in brackets, so I'll account for 
        # this):
        # (frequency parameters for branches)  [frequencies at nodes] (see Yang & Roberts 1995 fig 1)
        #
        # Node #1  ( 0.25824  0.24176  0.25824  0.24176 )
        # Node #2  ( 0.00000  0.50000  0.00000  0.50000 )
        elif "(frequency parameters for branches)" in line:
            results["parameters"]["nodes"] = {}
            branch_freqs_found = True
        elif branch_freqs_found is True:
            if len(line_floats) > 0:
                node_res = re.match("Node \#(\d+)", line)
                node_num = int(node_res.group(1))
                results["parameters"]["nodes"][node_num] = {"root":False}
                results["parameters"]["nodes"][node_num]\
                    ["frequency parameters"] = line_floats[:4]
                if len(line_floats) > 4:
                    results["parameters"]["nodes"][node_num]\
                        ["base frequencies"] = {"T":line_floats[4],
                                                "C":line_floats[5],
                                                "A":line_floats[6],
                                                "G":line_floats[7]}
            else:
                root_res = root_re.match(line)
                if root_res is not None:
                    root_node = int(root_res.group(1))
                    results["parameters"]["nodes"][root_node]["root"] =\
                        True
                    branch_freqs_found = False
