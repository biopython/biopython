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
parser.add_argument('-a', '--analysis', help='name of analysis', default='analysis')
parser.add_argument('-o', action='store_true', help='save plots')
parser.add_argument('-v', action='store_true', help='show plots')
args = parser.parse_args()

seqs_data = None

# read from stdin if no filename provided
if args.sequences == '-':
    seqs_data = ','.join(line for line in sys.stdin)
else:
    with open(args.sequences, 'r') as seq_file:
        seqs_data = seq_file.read()

seqs_data = seqs_data.replace(':', ',').replace(';', ',').replace('\n', ',').replace(' ', '').split(',')
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
    v.bar(parameter='sequence_length', out=out, show=show)
    v.bar(parameter='molecular_weight', out=out, show=show)
    v.bar(parameter='theoretical_pi', out=out, show=show)
    v.bar(parameter='aromaticity', out=out, show=show)
    v.bar(parameter='instability_index', out=out, show=show)
    v.bar(parameter='extinction_coefficient', out=out, show=show)
    v.bar(parameter='secondary_structure_fraction', out=out, show=show)
    v.bar(parameter='amino_acid_percent', out=out, show=show)
    v.box(parameter='extinction_coefficient', out=out, show=show)
    v.box(parameter='secondary_structure_fraction', out=out, show=show)
    v.scatter(x_parameter='sequence_length', y_parameter='theoretical_pi', z_parameter='molecular_weight', size=(12, 7), out=out, show=show)
    v.scatter(x_parameter='sequence_length', y_parameter='theoretical_pi', size=(12, 7), out=out, show=show)
    v.scatter(x_parameter='molecular_weight', y_parameter='flexibility', z_parameter='flexibility', size=(12, 7), out=out, show=show)
    v.scatter(x_parameter='theoretical_pi', y_parameter='instability_index', size=(12, 7), out=out, show=show)
    v.bar_interactive(parameter='amino_acid_comp', out=out, show=show)
    v.bar_interactive(parameter='flexibility', out=out, show=show)
    v.bar_interactive(parameter='c_term_freq', out=out, show=show)
    v.bar_interactive(parameter='gravy', out=out, show=show)
    v.bar_interactive(parameter='amino_acid_percent', out=out, show=show)
    v.bar_interactive(parameter='instability_index', out=out, show=show)
    v.bar_interactive(parameter='molecular_weight', out=out, show=show)


seqs = create_db(seqs_data, args.analysis)
print_success(f'Database created successfully, entries = {len(seqs)}')
dict_all_sequences, dict_each_sequence = analyze(seqs)
print_success(f'Analysis performed successfully')
visualize(dict_all_sequences, dict_each_sequence, 'def' if args.o else None, True if args.v else False)
print_success(f'Visualizing performed successfully')
