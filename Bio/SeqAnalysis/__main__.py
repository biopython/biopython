import argparse
import os

from . import database, analyzer, visualizer


def print_success(msg: str):
    print(f'\033[92m[+] {msg}\033[0m')


sheader = '''
███████╗███████╗ ██████╗  █████╗ ███╗   ██╗ █████╗ ██╗  ██╗   ██╗███████╗██╗███████╗
██╔════╝██╔════╝██╔═══██╗██╔══██╗████╗  ██║██╔══██╗██║  ╚██╗ ██╔╝██╔════╝██║██╔════╝
███████╗█████╗  ██║   ██║███████║██╔██╗ ██║███████║██║   ╚████╔╝ ███████╗██║███████╗
╚════██║██╔══╝  ██║▄▄ ██║██╔══██║██║╚██╗██║██╔══██║██║    ╚██╔╝  ╚════██║██║╚════██║
███████║███████╗╚██████╔╝██║  ██║██║ ╚████║██║  ██║███████╗██║   ███████║██║███████║
╚══════╝╚══════╝ ╚══▀▀═╝ ╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝╚══════╝╚═╝   ╚══════╝╚═╝╚══════╝
                                                                                    
Sequences analysis and plotting
'''

parser = argparse.ArgumentParser(description=sheader, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-s', '--sequences', nargs='+', help='list of sequences')
parser.add_argument('-d', '--directory', help='directory for database')
parser.add_argument('-o', action='store_true', help='save plots')
parser.add_argument('-v', action='store_true', help='show plots')
args = parser.parse_args()

if not os.path.exists('Plots'):
    os.mkdir('Plots')
if not os.path.exists('Plots/interactive'):
    os.mkdir('Plots/interactive')
if not os.path.exists('Plots/static'):
    os.mkdir('Plots/static')

db = database(args.sequences, args.directory)
seqs = db.get()
print_success(f'Database created successfully, entries = {len(seqs)}')
an = analyzer(seqs)
results = an.results()
print_success(f'Analysis performed successfully')
v = visualizer(*results)

out = 'def' if args.o else None
show = True if args.v else False
v.bar(parameter='sequence_length', out=out, show=show)
v.bar(parameter='molecular_weight', out=out, show=show)
v.bar(parameter='theoretical_pi', out=out, show=show)
v.bar(parameter='aromaticity', out=out, show=show)
v.bar(parameter='instability_index', out=out, show=show)
v.bar(parameter='extinction_coefficient', out=out, show=show)
v.bar(parameter='secondary_structure_fraction', out=out, show=show)
v.bar(parameter='amino_acid_percent', out=out, show=show)
v.bar(parameter='amino_acid_percent', single_seq=list(seqs.keys())[0], out=out, show=show)
v.bar(parameter='theoretical_pi', single_seq=list(seqs.keys())[0], out=out, show=show)
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

print_success(f'Visualizing performed successfully, plots saved in Plots directory')
