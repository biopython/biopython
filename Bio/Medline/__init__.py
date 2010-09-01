# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work with Medline.

Classes:
Record           A dictionary holding Medline data.

Functions:
read             Reads one Medline record
parse            Allows you to iterate over a bunch of Medline records
"""

class Record(dict):
    """A dictionary holding information from a Medline record.
    All data are stored under the mnemonic appearing in the Medline
    file. These mnemonics have the following interpretations:

    Mnemonic  Description
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
    """


def parse(handle):
    """Read Medline records one by one from the handle.

    The handle is either is a Medline file, a file-like object, or a list
    of lines describing one or more Medline records.

    Typical usage:

        from Bio import Medline
        handle = open("mymedlinefile")
        records = Medline.parse(handle)
        for record in record:
            print record['TI']

    """
    # These keys point to string values
    textkeys = ("ID", "PMID", "SO", "RF", "NI", "JC", "TA", "IS", "CY", "TT",
                "CA", "IP", "VI", "DP", "YR", "PG", "LID", "DA", "LR", "OWN",
                "STAT", "DCOM", "PUBM", "DEP", "PL", "JID", "SB", "PMC",
                "EDAT", "MHDA", "PST", "AB", "AD", "EA", "TI", "JT")
    handle = iter(handle)
    # First skip blank lines
    for line in handle:
        line = line.rstrip()
        if line:
            break
    else:
        return
    record = Record()
    finished = False
    while not finished:
        if line[:6]=="      ": # continuation line
            record[key].append(line[6:])
        elif line:
            key = line[:4].rstrip()
            if not key in record:
                record[key] = []
            record[key].append(line[6:])
        try:
            line = handle.next()
        except StopIteration:
            finished = True
        else:
            line = line.rstrip()
            if line:
                continue
        # Join each list of strings into one string.
        for key in textkeys:
            if key in record:
                record[key] = " ".join(record[key])
        if record:
            yield record
        record = Record()

def read(handle):
    """Read a single Medline records from the handle.

    The handle is either is a Medline file, a file-like object, or a list
    of lines describing a Medline record.

    Typical usage:

        from Bio import Medline
        handle = open("mymedlinefile")
        record = Medline.read(handle)
        print record['TI']

    """
    records = parse(handle)
    return records.next()
