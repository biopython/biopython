# Local test code for pull request related to issue 1517
# https://github.com/biopython/biopython/issues/1517
# I also tested the revised Entrez code in a test Flask app deployed to AWS Lambda.
# I'm happy to share that.
# - Chad Parmet, @cparmet on GitHubg


from biopython.Bio import Entrez
Entrez.cache = '/tmp/'
from biopython.Bio.Entrez import Parser


print(Entrez.cache)

# Pull up an article on PubMed
PMID = 21340628
Entrez.email = "<my email address>"
handle = Entrez.esummary(db="pubmed", id=PMID)
record = Entrez.read(handle)[0]
handle.close()

# Did esummary work with the custom directory?
print(record['PubDate'])
print(record['Source'])

# Print parameters
print(Parser.DataHandler.local_xsd_dir)
print(Parser.DataHandler.local_dtd_dir)
print(Entrez.cache)
