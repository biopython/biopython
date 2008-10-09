# This is a Python module.
"""This module is DEPRECATED.

DBXref was used in building SeqRecords from Martel descriptions (see
Bio.builders.SeqRecord.sequence for more details).

Andrew Dalke is no longer maintaining Martel or Bio.Mindy, and these modules
and associate ones like Bio.DBXRef are now deprecated.  They are no longer
used in any of the current Biopython parsers, and are likely to be removed
in a future release.
"""

import warnings
warnings.warn("Martel and those parts of Biopython closely linked to it" \
              +" (such as Bio.DBXRef which is not used elsewhere) are now" \
              +" deprecated, and will be removed in a future release of"\
              +" Biopython.  If you want to continue to use this code,"\
              +" please get in contact with the Biopython developers via"\
              +" the mailing lists to avoid its permanent removal from"\
              +" Biopython.", \
              DeprecationWarning)

class DBXRef:
    def __init__(self, dbname, dbid, reftype = None, negate = 0):
        self.dbname = dbname
        self.dbid = dbid
        self.reftype = reftype
        self.negate = negate

    def __str__(self):
        if self.reftype is None:
            reftype = ""
        else:
            reftype = self.reftype + "="
        s = "%s/%s%s" % (self.dbname, reftype, self.dbid)
        if self.negate:
            s = "not(%s)" % s
        return s
    __repr__ = __str__

class BioformatDBName:
    def __getitem__(self, name):
        return name
class UnknownDBName:
    def __getitem__(self, name):
        return "x-unknown-" + name.lower()

dbname_conversions = {
    "bioformat": BioformatDBName(),
    "unknown": UnknownDBName(),
    "sp": {"AARHUS/GHENT-2DPAGE": "x-aarhus-ghent-2dpage",
           "CARBBANK": "x-carbbank",
           "DICTYDB": "x-dictydb",
           "ECO2DBASE": "x-eco2dbase",
           "ECOGENE": "x-ecogene",
           "EMBL": "embl",     # EMBL (in GO)
           "FLYBASE": "fb",    # Flybase (in GO)
           "GCRDB": "x-gcrdb",
           "HIV": "x-hiv",
           "HSC-2DPAGE": "x-hsc",
           "HSSP": "x-hssp",
           "MAIZE-2DPAGE": "x-maize",
           "MAIZEDB": "x-maizedb",
           "MENDEL": "x-mendel",
           "MGD": "mgd",       # (in GO)
           "MIM": "x-mim",
           "PDB": "x-pdb",       # Protein Data Bank
           "PFAM": "x-pfam",
           "PIR": "pir", # GO
           "PROSITE": "x-prosite",
           "REBASE": "x-rebase",
           "SGD": "sgd", # GO
           "STYGENE": "x-stygene",
           "SUBTILIST": "x-subtilist",
           "SWISS-2DPAGE": "x-swiss",
           "TIGR": "tigr", # GO
           "TRANSFAC": "x-transfac",
           "WORMPEP": "x-wormpep",
           "YEPD": "x-yepd",
           "ZFIN": "x-zfin",
           },
    "go": {"CGEN": "cgen",    # Compugen, Inc.
           "DDB": "ddb",      # DictyBase (Dictyostelium discoideum)
           "DDBJ": "ddbj",    # DNA Database of Japan
           "EC": "ec",        # Enzyme Commission
           "EMBL": "embl",    # EMBL Nucleotide Sequence Data Library
           "ENSEMBL": "ensembl", # ENSEMBL
           "ENZYME": "enzyme",   # ENZYME
           "FB": "fb",        # FlyBase
           "GB": "gb",        # GenBank
           "GO": "go",        # Gene Ontology
           "GXD": "gxd",      # Gene Expression Database (mouse)
           "IPR": "ipr",      # InterPro
           "ISBN": "isbn",    # International Standard Book Number
           "IUBMB": "iubmb",  # International Union of Biochemistry
                              #      and Molecular Biology
           "IUPAC": "iupac",  # International Union of Pure and Applied
                              #      Chemistry
           "MEDLINE": "medline", # MEDLINE
           "MGD": "mgd",     # Mouse Genome Database
           "MGI": "mgi",     # Mouse Genome Informatics
           "NC-IUBMB": "NC-IUBMB",
                             #  Nomenclature Committee of the International
                             #     Union of Biochemistry and Molecular Biology
           "PIR": "pir",     # PIR
           "PMID": "pmid",   # PubMed
           "Pombase": "pombase", # Schizosaccharomyces pombe
           "Pompep": "pompep",   # Schizosaccharomyces pombe Protein
                                 #    Sequence Database
           "RESID": "resid", # RESID (protein post-translational modifications)
           "SGD": "sgd",     # Saccharomyces Genome Database
           "SP": "sp",       # SWISS-PROT
           "SWALL": "swall", # SWISS-PROT + TrEMBL + TrEMBLnew
           "TAIR": "tair",   # The Arabidopsis Information Resource
           "taxonID": "taxonid", # Taxonomy ID
           "TC": "tc",       # Transport Commission
           "TIGR": "tigr",   # The Institute of Genome Research
           "TR": "tr",       # TrEMBL
           "WB": "wb",       # WormBase (Caenorhabditis elegans)
    },
    # http://www.ncbi.nlm.nih.gov/collab/db_xref.html
    "genbank": {
      "ATCC": "x-atcc",  # American Type Culture Collection database
                         #    /db_xref="ATCC:123456"
      "ATCC(in host)": "x-atcc-host", # See above
      "ATCC(dna)": "x-atcc-dna",      # See above

      "BDGP_EST": "x-bdgp-est", #  Berkeley Drosophila Genome Project
                                #       EST database
                                #   /db_xref="BDGP_EST:123456"
      
      "BDGP_INS": "x-bdgp-ins", #  Berkeley Drosophila Genome Project
                                #      database -- Insertion
                                #   /db_xref="BDGP_INS:123456"


      "dbEST": "x-dbest",  #  EST database maintained at the NCBI.
                           #  /db_xref="dbEST:123456"

      "dbSNP": "x-dbsnp",  #  Variation database maintained at the NCBI.
                           #  /db_xref="dbSNP:4647"

      "dbSTS": "x-dbsts",  # STS database maintained at the NCBI.
                           # /db_xref="dbSTS:456789"

      "ENSEMBL": "ensembl", #  Database of automatically annotated genomic data
                            # /db_xref="ENSEMBL:HUMAN-Clone-AC005612"
                            # /db_xref="ENSEMBL:HUMAN-Gene-ENSG00000007102" 

      "ESTLIB": "x-estlib", # EBI's EST library identifier  #'
                            # /db_xref="ESTLIB:1200"

      "FANTOM_DB": "x-fantom-db", # Database of Functional Annotation of Mouse
                                  # /db_xref="FANTOM_DB:0610005A07"

      "FLYBASE": "fb", # Database of Genetic and molecular data of Drosophila.
                       # /db_xref="FLYBASE:FBgn0000024"

      "GDB": "x-gdb",  # Human Genome Database accession numbers.
                       # /db_xref="GDB:G00-128-600"

      "GI": "x-gi",    # GenInfo identifier, used as a unique sequence
                       # identifier for nucleotide and proteins.
                       # /db_xref="GI:1234567890"

      "GO": "go",      # Gene Ontology Database identifier
                       # /db_xref="GO:123"

      "IMGT/LIGM": "x-imgt-ligm", #  Immunogenetics database, immunoglobulins
                                  #  and T-cell receptors
                                  # /db_xref="IMGT/LIGM:U03895"

      "IMGT/HLA": "x-imgt-hla",   # Immunogenetics database, human MHC
                                  # /db_xref="IMGT/HLA:HLA00031"


      "LocusID": "x-locus-id", # NCBI LocusLink ID.
                               # /db_xref="LocusID:51199"

      "MaizeDB": "x-maizedb",  # Maize Genome Database unique identifiers.
                               # /db_xref="MaizeDB:Probe/79847"

      "MGD": "mgd",  # Mouse Genome Database accession numbers.
                     # /db_xref="MGD:123456"

      "MGI": "mgi",  # Medicago Genome Initiative
                     # /db_xref="MGI:S:20819"

      "MIM": "x-mim", # Mendelian Inheritance in Man numbers.
                      # /db_xref="MIM:123456"

      "niaEST": "x-niaEST", # NIA Mouse cDNA Project
                            # /db_xref="niaEST:L0304H12-3"

      "PIR": "pir", # Protein Information Resource accession numbers.
                    # /db_xref="PIR:S12345"

      "PSEUDO": "x-pseudo-embl", #  EMBL pseudo protein identifier
                                 # /db_xref="PSEUDO:CAC44644.1"

      "RATMAP": "x-ratmap", #  Rat Genome Database 
                            # /db_xref="RATMAP:5"

      "RiceGenes": "x-ricegenes", #  Rice database accession numbers.
                                  # /db_xref="RiceGenes:AA231856"

      "REMTREMBL": "x-remtrembl",
              # Computer-annotated protein sequence database containing
              # the translations of those codings sequences (CDS) present
              # in the EMBL Nucleotide Sequence Database that won't be  '
              # included in SWISS-PROT. These include: immunoglobulins and
              # T-cell receptors, synthetic sequences, patent application
              # sequences, small fragments, CDS not coding for real 
              # proteins and truncated proteins.
              # example:      /db_xref="REMTREMBL:CAC01666"

      "RZPD": "x-rzpd", # Resource Centre Primary Database Clone Identifiers
                        # /db_xref="RZPD:IMAGp998I142450Q6"

      "SGD": "sgd",  # Saccharomyces Genome Database accession numbers.
                     # /db_xref="SGD:L0000470"

      "SoyBase": "x-soybase", #  Glycine max Genome Database 
                              # /db_xref="SoyBase:Satt005"

      "SPTREMBL": "x-sptrembl",  # is this the same as "swall" ?
              # Computer-annotated protein sequence database 
              # supplementing SWISS-PROT and containing the 
              # translations of all coding sequences (CDS) 
              # present in the EMBL Nucleotide Sequence 
              # Database not yet integrated in SWISS-PROT. 
              #   /db_xref="SPTREMBL:Q00177"                    

      "SWISS-PROT": "sp", # Swiss-Prot protein database accession numbers.
                          # /db_xref="SWISS-PROT:P12345"

      "taxon": "taxonid", #  NCBI taxonomic identifier.
                          # /db_xref="taxon:4932"
      },
}

def from_parser(dbname_style, dbname, idtype, dbid, negate):
    try:
        dbname = dbname_conversions[dbname_style][dbname]
    except KeyError:
        dbname = "x-unknown2-%s--%s" % (dbname_style, dbname)
    return DBXRef(dbname, dbid, idtype, negate)
