# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to work with Medline from the NCBI.

Classes:
 - Record           A dictionary holding Medline data.

Functions:
 - read             Reads one Medline record
 - parse            Allows you to iterate over a bunch of Medline records
"""

__docformat__ = "restructuredtext en"


class Record(dict):
    """A dictionary holding information from a Medline record.

    All data are stored under the mnemonic appearing in the Medline
    file. These mnemonics have the following interpretations:

    ========= ==============================
    Mnemonic  Description
    --------- ------------------------------
    AB        Abstract
    CI        Copyright Information
    AD        Affiliation
    IRAD      Investigator Affiliation
    AID       Article Identifier
    AU        Author
    FAU       Full Author
    CN        Corporate Author
    DCOM      Date Completed
    DA        Date Created
    LR        Date Last Revised
    DEP       Date of Electronic Publication
    DP        Date of Publication
    EDAT      Entrez Date
    GS        Gene Symbol
    GN        General Note
    GR        Grant Number
    IR        Investigator Name
    FIR       Full Investigator Name
    IS        ISSN
    IP        Issue
    TA        Journal Title Abbreviation
    JT        Journal Title
    LA        Language
    LID       Location Identifier
    MID       Manuscript Identifier
    MHDA      MeSH Date
    MH        MeSH Terms
    JID       NLM Unique ID
    RF        Number of References
    OAB       Other Abstract
    OCI       Other Copyright Information
    OID       Other ID
    OT        Other Term
    OTO       Other Term Owner
    OWN       Owner
    PG        Pagination
    PS        Personal Name as Subject
    FPS       Full Personal Name as Subject
    PL        Place of Publication
    PHST      Publication History Status
    PST       Publication Status
    PT        Publication Type
    PUBM      Publishing Model
    PMC       PubMed Central Identifier
    PMID      PubMed Unique Identifier
    RN        Registry Number/EC Number
    NM        Substance Name
    SI        Secondary Source ID
    SO        Source
    SFM       Space Flight Mission
    STAT      Status
    SB        Subset
    TI        Title
    TT        Transliterated Title
    VI        Volume
    CON       Comment on
    CIN       Comment in
    EIN       Erratum in
    EFR       Erratum for
    CRI       Corrected and Republished in
    CRF       Corrected and Republished from
    PRIN      Partial retraction in
    PROF      Partial retraction of
    RPI       Republished in
    RPF       Republished from
    RIN       Retraction in
    ROF       Retraction of
    UIN       Update in
    UOF       Update of
    SPIN      Summary for patients in
    ORI       Original report in
    ========= ==============================
    """


def parse(handle):
    """Read Medline records one by one from the handle.

    The handle is either is a Medline file, a file-like object, or a list
    of lines describing one or more Medline records.

    Typical usage::

        from Bio import Medline
        with open("mymedlinefile") as handle:
            records = Medline.parse(handle)
            for record in record:
                print(record['TI'])

    """
    # TODO - Turn that into a working doctest
    # These keys point to string values
    textkeys = ("ID", "PMID", "SO", "RF", "NI", "JC", "TA", "IS", "CY", "TT",
                "CA", "IP", "VI", "DP", "YR", "PG", "LID", "DA", "LR", "OWN",
                "STAT", "DCOM", "PUBM", "DEP", "PL", "JID", "SB", "PMC",
                "EDAT", "MHDA", "PST", "AB", "AD", "EA", "TI", "JT")
    handle = iter(handle)

    key = ""
    record = Record()
    for line in handle:
        line = line.rstrip()
        if line[:6] == "      ":  # continuation line
            if key == "MH":
                # Multi-line MESH term, want to append to last entry in list
                record[key][-1] += line[5:]  # including space using line[5:]
            else:
                record[key].append(line[6:])
        elif line:
            key = line[:4].rstrip()
            if key not in record:
                record[key] = []
            record[key].append(line[6:])
        elif record:
            # Join each list of strings into one string.
            for key in record:
                if key in textkeys:
                    record[key] = " ".join(record[key])
            yield record
            record = Record()
    if record:  # catch last one
        for key in record:
            if key in textkeys:
                record[key] = " ".join(record[key])
        yield record


def read(handle):
    """Read a single Medline record from the handle.

    The handle is either is a Medline file, a file-like object, or a list
    of lines describing a Medline record.

    Typical usage:

        >>> from Bio import Medline
        >>> with open("mymedlinefile") as handle:
        ...     record = Medline.read(handle)
        ...     print(record['TI'])

    """
    # TODO - Turn that into a working doctest
    records = parse(handle)
    return next(records)
