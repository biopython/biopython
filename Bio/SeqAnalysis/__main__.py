import argparse
import doctest
import os
from pprint import pprint

from Bio.SeqAnalysis import test_analyzer, test_database, SeqDatabase, SeqAnalyzer, SeqVisualizer


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
parser.add_argument('-t', action='store_true', help='tests')
parser.add_argument('-o', action='store_true', help='save output')
parser.add_argument('-v', action='store_true', help='show output')
args = parser.parse_args()

if args.t:
    with open('test.txt', 'w') as test_file:
        test_file.write('A0A1A2:B0B1B2@C0C1C2!D0D1D2/E0E1E2|F0F1F2#G0G1G2\nH0H1H2\nI0I1I2$J0J1J2')
    print('======= Testing SeqDatabase =======')
    test_database()
    os.remove('test.txt')
    print('======= Testing SeqAnalyzer =======')
    test_analyzer()
    print('======= Testing SeqVisualizer =======')
    doctest.testfile('SeqVisualizer.py')
    exit(0)

if not os.path.exists('Plots'):
    os.mkdir('Plots')
if not os.path.exists('Plots/interactive'):
    os.mkdir('Plots/interactive')
if not os.path.exists('Plots/static'):
    os.mkdir('Plots/static')

database = SeqDatabase(args.sequences, args.directory)
seqs = database.get()
print_success(f'Database created successfully, entries = {len(seqs)}')
analyzer = SeqAnalyzer(seqs)
results = analyzer.results()
print_success(f'Analysis performed successfully')
v = SeqVisualizer(*results)

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

print_success(f'Visualization performed successfully, plots saved to Plots directory')
