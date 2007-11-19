"""Code for storing and retrieving biological objects from relational dbs.
"""
# Dictionary of database types, keyed by GenBank db_xref abbreviation
_db_dict = {'GeneID': 'Entrez',
           'GI': 'GeneIndex',
           'COG': 'COG',
           'CDD': 'CDD',
           'DDBJ': 'DNA Databank of Japan',
           'Entrez': 'Entrez',
           'GeneIndex': 'GeneIndex',
           'PUBMED': 'PubMed',
           'taxon': 'Taxon',
           'ATCC': 'ATCC',
           'ISFinder': 'ISFinder',
           'GOA': 'Gene Ontology Annotation',
           'ASAP': 'ASAP',
           'PSEUDO': 'PSEUDO',
           'InterPro': 'InterPro',
           'GEO': 'Gene Expression Omnibus',
           'EMBL': 'EMBL',
           'UniProtKB/Swiss-Prot': 'UniProtKB/Swiss-Prot',
           'ECOCYC': 'EcoCyc',
           'UniProtKB/TrEMBL': 'UniProtKB/TrEMBL'
           }
