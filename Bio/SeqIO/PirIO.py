# Copyright 2008-2015 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for the "pir" (aka PIR or NBRF) file format.

This module is for reading and writing PIR or NBRF format files as
SeqRecord objects.

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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import _get_seq_string
from .Interfaces import SequenceIterator
from .Interfaces import SequenceWriter


_pir_mol_type = {
    "P1": "protein",
    "F1": "protein",
    "D1": "DNA",
    "DL": "DNA",
    "DC": "DNA",
    "RL": "RNA",
    "RC": "RNA",
    "N3": "RNA",
    "XX": None,
}


class PirIterator(SequenceIterator):
    """Parser for PIR files."""

    def __init__(self, source):
        """Iterate over a PIR file and yield SeqRecord objects.

        source - file-like object or a path to a file.

        Examples
        --------
        >>> with open("NBRF/DMB_prot.pir") as handle:
        ...    for record in PirIterator(handle):
        ...        print("%s length %i" % (record.id, len(record)))
        HLA:HLA00489 length 263
        HLA:HLA00490 length 94
        HLA:HLA00491 length 94
        HLA:HLA00492 length 80
        HLA:HLA00493 length 175
        HLA:HLA01083 length 188

        """
        super().__init__(source, mode="t", fmt="Pir")

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        records = self.iterate(handle)
        return records

    def iterate(self, handle):
        """Iterate over the records in the PIR file."""
        # Skip any text before the first record (e.g. blank lines, comments)
        for line in handle:
            if line[0] == ">":
                break
        else:
            return  # Premature end of file, or just empty?

        while True:
            pir_type = line[1:3]
            if pir_type not in _pir_mol_type or line[3] != ";":
                raise ValueError(
                    "Records should start with '>XX;' where XX is a valid sequence type"
                )
            identifier = line[4:].strip()
            description = handle.readline().strip()

            lines = []
            for line in handle:
                if line[0] == ">":
                    break
                # Remove trailing whitespace, and any internal spaces
                lines.append(line.rstrip().replace(" ", ""))
            else:
                line = None
            seq = "".join(lines)
            if seq[-1] != "*":
                # Note the * terminator is present on nucleotide sequences too,
                # it is not a stop codon!
                raise ValueError(
                    "Sequences in PIR files should include a * terminator!"
                )

            # Return the record and then continue...
            record = SeqRecord(
                Seq(seq[:-1]), id=identifier, name=identifier, description=description
            )
            record.annotations["PIR-type"] = pir_type
            if _pir_mol_type[pir_type]:
                record.annotations["molecule_type"] = _pir_mol_type[pir_type]
            yield record

            if line is None:
                return  # StopIteration
        raise ValueError("Unrecognised PIR record format.")


class PirWriter(SequenceWriter):
    """Class to write PIR format files."""

    def __init__(self, handle, wrap=60, record2title=None, code=None):
        """Create a PIR writer.

        Arguments:
         - handle - Handle to an output file, e.g. as returned
           by open(filename, "w")
         - wrap - Optional line length used to wrap sequence lines.
           Defaults to wrapping the sequence at 60 characters
           Use zero (or None) for no wrapping, giving a single
           long line for the sequence.
         - record2title - Optional function to return the text to be
           used for the title line of each record.  By default
           a combination of the record.id, record.name and
           record.description is used.
         - code - Optional sequence code must be one of P1, F1,
           D1, DL, DC, RL, RC, N3 and XX. By default None is used,
           which means auto detection based on the molecule type
           in the record annotation.

        You can either use::

            handle = open(filename, "w")
            writer = PirWriter(handle)
            writer.write_file(myRecords)
            handle.close()

        Or, follow the sequential file writer system, for example::

            handle = open(filename, "w")
            writer = PirWriter(handle)
            writer.write_header() # does nothing for PIR files
            ...
            Multiple writer.write_record() and/or writer.write_records() calls
            ...
            writer.write_footer() # does nothing for PIR files
            handle.close()

        """
        super().__init__(handle)
        self.wrap = None
        if wrap:
            if wrap < 1:
                raise ValueError("wrap should be None, 0, or a positive integer")
        self.wrap = wrap
        self.record2title = record2title
        self.code = code

    def write_record(self, record):
        """Write a single PIR record to the file."""
        if self.record2title:
            title = self.clean(self.record2title(record))
        else:
            title = self.clean(record.id)

        if record.name and record.description:
            description = self.clean(record.name + " - " + record.description)
        elif record.name and not record.description:
            description = self.clean(record.name)
        else:
            description = self.clean(record.description)

        if self.code:
            code = self.code
        else:
            molecule_type = record.annotations.get("molecule_type")
            if molecule_type is None:
                code = "XX"
            elif "DNA" in molecule_type:
                code = "D1"
            elif "RNA" in molecule_type:
                code = "RL"
            elif "protein" in molecule_type:
                code = "P1"
            else:
                code = "XX"

        if code not in _pir_mol_type:
            raise TypeError(
                "Sequence code must be one of " + _pir_mol_type.keys() + "."
            )
        assert "\n" not in title
        assert "\r" not in description

        self.handle.write(f">{code};{title}\n{description}\n")

        data = _get_seq_string(record)  # Catches sequence being None

        assert "\n" not in data
        assert "\r" not in data

        if self.wrap:
            line = ""
            for i in range(0, len(data), self.wrap):
                line += data[i : i + self.wrap] + "\n"
            line = line[:-1] + "*\n"
            self.handle.write(line)
        else:
            self.handle.write(data + "*\n")


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
