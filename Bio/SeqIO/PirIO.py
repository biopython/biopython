# Copyright 2008-2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# This module is for reading and writing PIR or NBRF format files as
# SeqRecord objects.  The code is based on Bio.SeqIO.FastaIO

"""Bio.SeqIO support for the "pir" (aka PIR or NBRF) file format.

You are expected to use this module via the Bio.SeqIO functions, or if
the file contains a sequence alignment, optionally via Bio.AlignIO instead.

This format was introduced for the Protein Information Resource (PIR), a
project of the National Biomedical Research Foundation (NBRF).  The PIR
database itself is now part of UniProt.

The file format is described online at:
http://www.ebi.ac.uk/help/pir_frame.html
http://www.cmbi.kun.nl/bioinf/tools/crab_pir.html (currently down)

An example file in this format would be::

  >P1;CRAB_ANAPL
  ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN).
    MDITIHNPLI RRPLFSWLAP SRIFDQIFGE HLQESELLPA SPSLSPFLMR
    SPIFRMPSWL ETGLSEMRLE KDKFSVNLDV KHFSPEELKV KVLGDMVEIH
    GKHEERQDEH GFIAREFNRK YRIPADVDPL TITSSLSLDG VLTVSAPRKQ
    SDVPERSIPI TREEKPAIAG AQRK*

  >P1;CRAB_BOVIN
  ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN).
    MDIAIHHPWI RRPFFPFHSP SRLFDQFFGE HLLESDLFPA STSLSPFYLR
    PPSFLRAPSW IDTGLSEMRL EKDRFSVNLD VKHFSPEELK VKVLGDVIEV
    HGKHEERQDE HGFISREFHR KYRIPADVDP LAITSSLSSD GVLTVNGPRK
    QASGPERTIP ITREEKPAVT AAPKK*

Or, an example of a multiple sequence alignment::

  >P1;S27231
  rhodopsin - northern leopard frog
  MNGTEGPNFY IPMSNKTGVV RSPFDYPQYY LAEPWKYSVL AAYMFLLILL GLPINFMTLY
  VTIQHKKLRT PLNYILLNLG VCNHFMVLCG FTITMYTSLH GYFVFGQTGC YFEGFFATLG
  GEIALWSLVV LAIERYIVVC KPMSNFRFGE NHAMMGVAFT WIMALACAVP PLFGWSRYIP
  EGMQCSCGVD YYTLKPEVNN ESFVIYMFVV HFLIPLIIIS FCYGRLVCTV KEAAAQQQES
  ATTQKAEKEV TRMVIIMVIF FLICWVPYAY VAFYIFTHQG SEFGPIFMTV PAFFAKSSAI
  YNPVIYIMLN KQFRNCMITT LCCGKNPFGD DDASSAATSK TEATSVSTSQ VSPA*

  >P1;I51200
  rhodopsin - African clawed frog
  MNGTEGPNFY VPMSNKTGVV RSPFDYPQYY LAEPWQYSAL AAYMFLLILL GLPINFMTLF
  VTIQHKKLRT PLNYILLNLV FANHFMVLCG FTVTMYTSMH GYFIFGPTGC YIEGFFATLG
  GEVALWSLVV LAVERYIVVC KPMANFRFGE NHAIMGVAFT WIMALSCAAP PLFGWSRYIP
  EGMQCSCGVD YYTLKPEVNN ESFVIYMFIV HFTIPLIVIF FCYGRLLCTV KEAAAQQQES
  LTTQKAEKEV TRMVVIMVVF FLICWVPYAY VAFYIFTHQG SNFGPVFMTV PAFFAKSSAI
  YNPVIYIVLN KQFRNCLITT LCCGKNPFGD EDGSSAATSK TEASSVSSSQ VSPA*

  >P1;JN0120
  rhodopsin - Japanese lamprey
  MNGTEGDNFY VPFSNKTGLA RSPYEYPQYY LAEPWKYSAL AAYMFFLILV GFPVNFLTLF
  VTVQHKKLRT PLNYILLNLA MANLFMVLFG FTVTMYTSMN GYFVFGPTMC SIEGFFATLG
  GEVALWSLVV LAIERYIVIC KPMGNFRFGN THAIMGVAFT WIMALACAAP PLVGWSRYIP
  EGMQCSCGPD YYTLNPNFNN ESYVVYMFVV HFLVPFVIIF FCYGRLLCTV KEAAAAQQES
  ASTQKAEKEV TRMVVLMVIG FLVCWVPYAS VAFYIFTHQG SDFGATFMTL PAFFAKSSAL
  YNPVIYILMN KQFRNCMITT LCCGKNPLGD DE-SGASTSKT EVSSVSTSPV SPA*


As with the FASTA format, each record starts with a line beginning with ">"
character.  There is then a two letter sequence type (P1, F1, DL, DC, RL,
RC, or XX), a semi colon, and the identification code.  The second like is
free text description.  The remaining lines contain the sequence itself,
terminating in an asterisk.  Space separated blocks of ten letters as shown
above are typical.

Sequence codes and their meanings:

  - P1 - Protein (complete)
  - F1 - Protein (fragment)
  - D1 - DNA (e.g. EMBOSS seqret output)
  - DL - DNA (linear)
  - DC - DNA (circular)
  - RL - RNA (linear)
  - RC - RNA (circular)
  - N3 - tRNA
  - N1 - Other functional RNA
  - XX - Unknown
"""

from __future__ import print_function

from Bio.Alphabet import single_letter_alphabet, generic_protein, \
    generic_dna, generic_rna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

__docformat__ = "restructuredtext en"

_pir_alphabets = {"P1": generic_protein,
                  "F1": generic_protein,
                  "D1": generic_dna,
                  "DL": generic_dna,
                  "DC": generic_dna,
                  "RL": generic_rna,
                  "RC": generic_rna,
                  "N3": generic_rna,
                  "XX": single_letter_alphabet,
                  }


def PirIterator(handle):
    """Generator function to iterate over Fasta records (as SeqRecord objects).

    handle - input file
    alphabet - optional alphabet
    title2ids - A function that, when given the title of the FASTA
    file (without the beginning >), will return the id, name and
    description (in that order) for the record as a tuple of strings.

    If this is not given, then the entire title line will be used
    as the description, and the first word as the id and name.

    Note that use of title2ids matches that of Bio.Fasta.SequenceParser
    but the defaults are slightly different.
    """
    # Skip any text before the first record (e.g. blank lines, comments)
    while True:
        line = handle.readline()
        if line == "":
            return  # Premature end of file, or just empty?
        if line[0] == ">":
            break

    while True:
        if line[0] != ">":
            raise ValueError(
                "Records in PIR files should start with '>' character")
        pir_type = line[1:3]
        if pir_type not in _pir_alphabets or line[3] != ";":
            raise ValueError(
                "Records should start with '>XX;' "
                "where XX is a valid sequence type")
        identifier = line[4:].strip()
        description = handle.readline().strip()

        lines = []
        line = handle.readline()
        while True:
            if not line:
                break
            if line[0] == ">":
                break
            # Remove trailing whitespace, and any internal spaces
            lines.append(line.rstrip().replace(" ", ""))
            line = handle.readline()
        seq = "".join(lines)
        if seq[-1] != "*":
            # Note the * terminator is present on nucleotide sequences too,
            # it is not a stop codon!
            raise ValueError(
                "Sequences in PIR files should include a * terminator!")

        # Return the record and then continue...
        record = SeqRecord(Seq(seq[:-1], _pir_alphabets[pir_type]),
                           id=identifier, name=identifier,
                           description=description)
        record.annotations["PIR-type"] = pir_type
        yield record

        if not line:
            return  # StopIteration
    assert False, "Should not reach this line"

if __name__ == "__main__":
    print("Running quick self test")

    import os

    for name in ["clustalw", "DMA_nuc", "DMB_prot", "B_nuc", "Cw_prot"]:
        print(name)
        filename = "../../Tests/NBRF/%s.pir" % name
        if not os.path.isfile(filename):
            print("Missing %s" % filename)
            continue

        records = list(PirIterator(open(filename)))
        count = 0
        for record in records:
            count += 1
            parts = record.description.split()
            if "bases," in parts:
                assert len(record) == int(parts[parts.index("bases,") - 1])
        print("Could read %s (%i records)" % (name, count))
