"""Database definitions for retrieval of Medline information.
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

medline_eutils = EUtilsDB(
        name = "medline-eutils",
        doc = "Retrieve Medline data from NCBI using EUtils",
        delay = 5.0,
        db = "pubmed",
        rettype = "medline",
        failure_cases = ncbi_failures
        )

medline = DBGroup(
        name = "medline",
        behavior = "serial")
medline.add(medline_eutils)
