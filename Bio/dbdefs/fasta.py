"""Database definitions for retrieval of FASTA formatted sequences.
"""
from Bio.config.DBRegistry import DBGroup, EUtilsDB
from _support import *

from Martel import *

proxy_error_expr = has_expr(Alt(Str("500"), Str("502")) + Str(" Proxy Error"))
diagnostic_error_expr = has_str("WWW Error 500 Diagnostic")
error_expr = Str("ERROR")

ncbi_failures=[
    (proxy_error_expr, "proxy error"),
    (diagnostic_error_expr, "diagnostic error"),
    (error_expr, "ERROR"),
    (html_expr, "I got HTML and shouldn't have"),
    (Str("Please try again later"), "Please try again later"),
    (Str("The sequence has been intentionally withdrawn"),
     "Sequence withdrawn"),
    (blank_expr, "No data returned")
    ]

fasta_sequence_eutils = EUtilsDB(
        name = "fasta-sequence-eutils",
        doc = "Retrieve Fasta formatted sequence data from NCBI using EUtils",
        delay = 5.0,
        db = "sequences",
        rettype = "fasta",
        failure_cases = ncbi_failures
        )

fasta = DBGroup(
        name = "fasta",
        behavior = "serial")
fasta.add(fasta_sequence_eutils)
