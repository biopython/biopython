# Copyright 2019-2020 by moozeq, KacDom, erpeg and s0lczi. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


import argparse
import os
import sys
from typing import List, Tuple

from . import database, analyzer, visualizer


def print_success(msg: str):
    print(f'\033[92m[+] {msg}\033[0m')


sa_header = '''
███████╗███████╗ ██████╗  █████╗ ███╗   ██╗ █████╗ ██╗  ██╗   ██╗███████╗██╗███████╗
██╔════╝██╔════╝██╔═══██╗██╔══██╗████╗  ██║██╔══██╗██║  ╚██╗ ██╔╝██╔════╝██║██╔════╝
███████╗█████╗  ██║   ██║███████║██╔██╗ ██║███████║██║   ╚████╔╝ ███████╗██║███████╗
╚════██║██╔══╝  ██║▄▄ ██║██╔══██║██║╚██╗██║██╔══██║██║    ╚██╔╝  ╚════██║██║╚════██║
███████║███████╗╚██████╔╝██║  ██║██║ ╚████║██║  ██║███████╗██║   ███████║██║███████║
╚══════╝╚══════╝ ╚══▀▀═╝ ╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝╚══════╝╚═╝   ╚══════╝╚═╝╚══════╝
                                                                                    
Sequences analysis and plotting
'''

parser = argparse.ArgumentParser(description=sa_header, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('sequences', nargs='?', default='-', type=str, help='sequences as filename or list')
parser.add_argument('-d', '--directory', help='directory for analysis output', default='seq_analysis')
args = parser.parse_args()

seqs_data = None

# read from stdin if no filename provided
if args.sequences == '-':
    seqs_data = ','.join(line for line in sys.stdin)
else:
    with open(args.sequences, 'r') as seq_file:
        seqs_data = seq_file.read()

seqs_data = seqs_data.replace(':', ',').replace(';', ',').replace('\n', ',').replace('\t', ',').replace(' ', '').split(',')
seqs_data = [seq.strip() for seq in seqs_data if seq]

analysis_dir = args.directory
if not os.path.exists(analysis_dir):
    os.mkdir(analysis_dir)
with open(f'{analysis_dir}/sequences.txt', 'w') as seq_file:
    seq_file.write('\n'.join(seqs_data))


def create_db(sequences: List[str], analysis_name: str) -> dict:
    return database(sequences, analysis_name, analysis_dir).get()


def analyze(sequences: dict) -> Tuple[dict, dict]:
    return analyzer(sequences, analysis_dir).results()


def visualize(all_seqs: dict, separated_seqs: dict, out: str, show: bool):
    v = visualizer(all_seqs, separated_seqs, analysis_dir)
    legend = True if len(v.separated_sequences_dict) < 30 else False
    v.bar(parameter='sequence_length', out=out, show=show, add_legend=legend)
    v.bar(parameter='molecular_weight', out=out, show=show, add_legend=legend)
    v.bar(parameter='theoretical_pi', out=out, show=show, add_legend=legend)
    v.bar(parameter='aromaticity', out=out, show=show, add_legend=legend)
    v.bar(parameter='instability_index', out=out, show=show, add_legend=legend)
    v.bar(parameter='extinction_coefficient', out=out, show=show, add_legend=legend)
    v.bar(parameter='secondary_structure_fraction', out=out, show=show, add_legend=legend)
    v.bar(parameter='amino_acid_percent', out=out, show=show, add_legend=legend)
    v.box(parameter='extinction_coefficient', out=out, show=show, add_legend=legend)
    v.box(parameter='secondary_structure_fraction', out=out, show=show, add_legend=legend)
    v.scatter(x_parameter='sequence_length', y_parameter='gravy', z_parameter='theoretical_pi', size=(12, 7), out=out, show=show, add_legend=legend)
    v.scatter(x_parameter='sequence_length', y_parameter='theoretical_pi', size=(12, 7), out=out, show=show, add_legend=legend)
    v.scatter(x_parameter='molecular_weight', y_parameter='flexibility', z_parameter='sequence_length', size=(12, 7), out=out, show=show, add_legend=legend)
    v.scatter(x_parameter='theoretical_pi', y_parameter='instability_index', size=(12, 7), out=out, show=show, add_legend=legend)
    v.bar_interactive(parameter='amino_acid_comp', out=out, show=show)
    v.bar_interactive(parameter='flexibility', out=out, show=show)
    v.bar_interactive(parameter='c_term_freq', out=out, show=show)
    v.bar_interactive(parameter='gravy', out=out, show=show)
    v.bar_interactive(parameter='amino_acid_percent', out=out, show=show)
    v.bar_interactive(parameter='instability_index', out=out, show=show)
    v.bar_interactive(parameter='molecular_weight', out=out, show=show)


seqs = create_db(seqs_data, 'analysis')
print_success(f'Database created successfully, entries = {len(seqs)}')
dict_all_sequences, dict_each_sequence = analyze(seqs)
print_success(f'Analysis performed successfully')
visualize(dict_all_sequences, dict_each_sequence, 'def', False)
print_success(f'Visualizing performed successfully')
