# Copyright (C) 2011, 2019 by Brandon Invergo (b.invergo@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Methods for parsing yn00 results files."""

import re


def parse_ng86(lines, results):
    """Parse the Nei & Gojobori (1986) section of the results.

    Nei_Gojobori results are organized in a lower
    triangular matrix, with the sequence names labeling
    the rows and statistics in the format:
    w (dN dS) per column
    Example row (2 columns):
    0.0000 (0.0000 0.0207) 0.0000 (0.0000 0.0421)
    """
    sequences = []
    for line in lines:
        # The purpose of this complex regex is to parse the NG86 section for
        # valid lines of data that are mixed in with citations and comments.
        # The data lines begin with a taxon name, followed by zero or more
        # fields containing numeric values, sometimes enclosed in parens.
        # Taxon names are from 1-30 characters and are usually separated from
        # the numeric portion of the line by space(s).  Long taxon names to are
        # truncated to 30 characters, and may run into the data fields without
        # any separator., e.g. some_long_name-1.0000
        # This regex is an attempt to cover more pathological cases while also
        # parsing all existing versions of yn00 output with shorter names.

        matrix_row_res = re.match(
            r"^([^\s]+?)(\s+-?\d+\.\d+.*$|\s*$|-1.0000\s*\(.*$)", line
        )
        if matrix_row_res is not None:
            # Find all floating point numbers in this line, accounting
            # for the fact that the sequence IDs might have bits that
            # look like floating point values.
            line_floats_res = re.findall(r"-*\d+\.\d+", matrix_row_res.group(2))
            line_floats = [float(val) for val in line_floats_res]

            seq_name = matrix_row_res.group(1).strip()
            sequences.append(seq_name)
            results[seq_name] = {}
            for i in range(0, len(line_floats), 3):
                NG86 = {}
                NG86["omega"] = line_floats[i]
                NG86["dN"] = line_floats[i + 1]
                NG86["dS"] = line_floats[i + 2]
                results[seq_name][sequences[i // 3]] = {"NG86": NG86}
                results[sequences[i // 3]][seq_name] = {"NG86": NG86}
    return (results, sequences)


def parse_yn00(lines, results, sequences):
    """Parse the Yang & Nielsen (2000) part of the results.

    Yang & Nielsen results are organized in a table with
    each row comprising one pairwise species comparison.
    Rows are labeled by sequence number rather than by
    sequence name.
    """
    # Example (header row and first table row):
    # seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE
    # 2    1    67.3   154.7   0.0136  3.6564  0.0000 -0.0000 +- 0.0000  0.0150
    # +- 0.0151
    for line in lines:
        # Find all floating point numbers in this line
        line_floats_res = re.findall(r"-*\d+\.\d+", line)
        line_floats = [float(val) for val in line_floats_res]
        row_res = re.match(r"\s+(\d+)\s+(\d+)", line)
        if row_res is not None:
            seq1 = int(row_res.group(1))
            seq2 = int(row_res.group(2))
            seq_name1 = sequences[seq1 - 1]
            seq_name2 = sequences[seq2 - 1]
            YN00 = {}
            YN00["S"] = line_floats[0]
            YN00["N"] = line_floats[1]
            YN00["t"] = line_floats[2]
            YN00["kappa"] = line_floats[3]
            YN00["omega"] = line_floats[4]
            YN00["dN"] = line_floats[5]
            YN00["dN SE"] = line_floats[6]
            YN00["dS"] = line_floats[7]
            YN00["dS SE"] = line_floats[8]
            results[seq_name1][seq_name2]["YN00"] = YN00
            results[seq_name2][seq_name1]["YN00"] = YN00
            seq_name1 = None
            seq_name2 = None
    return results


def parse_others(lines, results, sequences):
    """Parse the results from the other methods.

    The remaining methods are grouped together. Statistics
    for all three are listed for each of the pairwise
    species comparisons, with each method's results on its
    own line.
    The stats in this section must be handled differently
    due to the possible presence of NaN values, which won't
    get caught by my typical "line_floats" method used above.
    """
    # Example:
    # 2 (Pan_troglo) vs. 1 (Homo_sapie)

    # L(i):      143.0      51.0      28.0  sum=    222.0
    # Ns(i):    0.0000    1.0000    0.0000  sum=   1.0000
    # Nv(i):    0.0000    0.0000    0.0000  sum=   0.0000
    # A(i):     0.0000    0.0200    0.0000
    # B(i):    -0.0000   -0.0000   -0.0000
    # LWL85:  dS =  0.0227 dN =  0.0000 w = 0.0000 S =   45.0 N =  177.0
    # LWL85m: dS =    -nan dN =    -nan w =   -nan S =   -nan N =   -nan (rho = -nan)
    # LPB93:  dS =  0.0129 dN =  0.0000 w = 0.0000
    seq_name1 = None
    seq_name2 = None
    for line in lines:
        comp_res = re.match(r"\d+ \((.+)\) vs. \d+ \((.+)\)", line)
        if comp_res is not None:
            seq_name1 = comp_res.group(1)
            seq_name2 = comp_res.group(2)
        elif seq_name1 is not None and seq_name2 is not None:
            if "dS =" in line:
                stats = {}
                line_stats = line.split(":")[1].strip()
                # Find all of the xx = ###### values in a row
                # ie dS =  0.0227
                # For dN and dS, the values have 8 characters from the equals
                # sign, while the rest have 7 characters. On Windows,
                # NaNs take on weird values like -1.#IND, which might fill the
                # entire fixed column width.
                res_matches = re.findall(r"[dSNwrho]{1,3} =.{7,8}?", line_stats)
                for stat_pair in res_matches:
                    stat = stat_pair.split("=")[0].strip()
                    value = stat_pair.split("=")[1].strip()
                    try:
                        stats[stat] = float(value)
                    except ValueError:
                        stats[stat] = None
                if "LWL85:" in line:
                    results[seq_name1][seq_name2]["LWL85"] = stats
                    results[seq_name2][seq_name1]["LWL85"] = stats
                elif "LWL85m" in line:
                    results[seq_name1][seq_name2]["LWL85m"] = stats
                    results[seq_name2][seq_name1]["LWL85m"] = stats
                elif "LPB93" in line:
                    results[seq_name1][seq_name2]["LPB93"] = stats
                    results[seq_name2][seq_name1]["LPB93"] = stats
    return results
