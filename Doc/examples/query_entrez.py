"""Example code for querying Entrez and saving the results.

This code is meant to go with the 'Connecting with biological databases'
section of the Tutorial demonstrates Entrez connectivity.

See http://www.ncbi.nlm.nih.gov/entrez/query/static/linking.html for
more help understanding the parameters passed.

This also requires a web browser to run -- either netscape or lynx
are supported in this example."""
# standard library
import os

# biopython
from Bio.WWW import NCBI

search_command = 'Search'
search_database = 'Nucleotide'
return_format = 'FASTA'
search_term = 'Cypripedioideae'
my_browser = 'lynx'

result_handle = NCBI.query(search_command, search_database, term = search_term,
                           doptcmdl = return_format)

result_file_name = os.path.join(os.getcwd(), 'results.html')
result_file = open(result_file_name, 'w')
result_file.write(result_handle.read())
result_file.close()

if my_browser == 'lynx':
    os.system('lynx -force_html ' + result_file_name)
elif my_browser == 'netscape':
    os.system('netscape file:' + result_file_name)
