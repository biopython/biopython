# Copyright 2000 Andrew Dalke.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Codon tables based on those from the NCBI.

These tables are based on parsing the NCBI file
ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
using Scripts/update_ncbi_codon_table.py

Last updated at Version 4.4 (May 2019)
"""

from Bio.Data import IUPACData


unambiguous_dna_by_name = {}
unambiguous_dna_by_id = {}
unambiguous_rna_by_name = {}
unambiguous_rna_by_id = {}
generic_by_name = {}  # unambiguous DNA or RNA
generic_by_id = {}  # unambiguous DNA or RNA

ambiguous_dna_by_name = {}
ambiguous_dna_by_id = {}
ambiguous_rna_by_name = {}
ambiguous_rna_by_id = {}
ambiguous_generic_by_name = {}  # ambiguous DNA or RNA
ambiguous_generic_by_id = {}  # ambiguous DNA or RNA

# standard IUPAC unambiguous codons
standard_dna_table = None
standard_rna_table = None


# In the future, the back_table could return a statistically
# appropriate distribution of codons, so do not cache the results of
# back_table lookups!


class TranslationError(Exception):
    """Container for translation specific exceptions."""

    pass


class CodonTable:
    """A codon-table, or genetic code."""

    forward_table = {}  # only includes codons which actually code
    back_table = {}  # for back translations
    start_codons = []
    stop_codons = []

    # Not always called from derived classes!
    def __init__(
        self,
        nucleotide_alphabet=None,
        protein_alphabet=None,
        forward_table=forward_table,
        back_table=back_table,
        start_codons=start_codons,
        stop_codons=stop_codons,
    ):
        """Initialize the class."""
        self.nucleotide_alphabet = nucleotide_alphabet
        self.protein_alphabet = protein_alphabet
        self.forward_table = forward_table
        self.back_table = back_table
        self.start_codons = start_codons
        self.stop_codons = stop_codons

    def __str__(self):
        """Return a simple text representation of the codon table.

        e.g.::

            >>> import Bio.Data.CodonTable
            >>> print(Bio.Data.CodonTable.standard_dna_table)
            Table 1 Standard, SGC0
            <BLANKLINE>
              |  T      |  C      |  A      |  G      |
            --+---------+---------+---------+---------+--
            T | TTT F   | TCT S   | TAT Y   | TGT C   | T
            T | TTC F   | TCC S   | TAC Y   | TGC C   | C
            ...
            G | GTA V   | GCA A   | GAA E   | GGA G   | A
            G | GTG V   | GCG A   | GAG E   | GGG G   | G
            --+---------+---------+---------+---------+--
            >>> print(Bio.Data.CodonTable.generic_by_id[1])
            Table 1 Standard, SGC0
            <BLANKLINE>
              |  U      |  C      |  A      |  G      |
            --+---------+---------+---------+---------+--
            U | UUU F   | UCU S   | UAU Y   | UGU C   | U
            U | UUC F   | UCC S   | UAC Y   | UGC C   | C
            ...
            G | GUA V   | GCA A   | GAA E   | GGA G   | A
            G | GUG V   | GCG A   | GAG E   | GGG G   | G
            --+---------+---------+---------+---------+--
        """
        if self.id:
            answer = "Table %i" % self.id
        else:
            answer = "Table ID unknown"
        if self.names:
            answer += " " + ", ".join([x for x in self.names if x])

        # Use the main four letters (and the conventional ordering)
        # even for ambiguous tables
        letters = self.nucleotide_alphabet
        if letters is not None and "T" in letters:
            letters = "TCAG"
        else:
            # Should be either RNA or generic nucleotides,
            # e.g. Bio.Data.CodonTable.generic_by_id[1]
            letters = "UCAG"

        # Build the table...
        answer += "\n\n"
        answer += "  |" + "|".join("  %s      " % c2 for c2 in letters) + "|"
        answer += "\n--+" + "+".join("---------" for c2 in letters) + "+--"
        for c1 in letters:
            for c3 in letters:
                line = c1 + " |"
                for c2 in letters:
                    codon = c1 + c2 + c3
                    line += " %s" % codon
                    if codon in self.stop_codons:
                        line += " Stop|"
                    else:
                        try:
                            amino = self.forward_table[codon]
                        except KeyError:
                            amino = "?"
                        except TranslationError:
                            amino = "?"
                        if codon in self.start_codons:
                            line += " %s(s)|" % amino
                        else:
                            line += " %s   |" % amino
                line += " " + c3
                answer += "\n" + line
            answer += "\n--+" + "+".join("---------" for c2 in letters) + "+--"
        return answer


def make_back_table(table, default_stop_codon):
    """Back a back-table (naive single codon mapping).

    ONLY RETURNS A SINGLE CODON, chosen from the possible alternatives
    based on their sort order.
    """
    # Do the sort so changes in the hash implementation won't affect
    # the result when one amino acid is coded by more than one codon.
    back_table = {}
    for key in sorted(table):
        back_table[table[key]] = key
    back_table[None] = default_stop_codon
    return back_table


class NCBICodonTable(CodonTable):
    """Codon table for generic nucleotide sequences."""

    nucleotide_alphabet = None
    protein_alphabet = IUPACData.protein_letters

    def __init__(self, id, names, table, start_codons, stop_codons):
        """Initialize the class."""
        self.id = id
        self.names = names
        self.forward_table = table
        self.back_table = make_back_table(table, stop_codons[0])
        self.start_codons = start_codons
        self.stop_codons = stop_codons

    def __repr__(self):
        """Represent the NCBI codon table class as a string for debugging."""
        return "%s(id=%r, names=%r, ...)" % (
            self.__class__.__name__,
            self.id,
            self.names,
        )


class NCBICodonTableDNA(NCBICodonTable):
    """Codon table for unambiguous DNA sequences."""

    nucleotide_alphabet = IUPACData.unambiguous_dna_letters


class NCBICodonTableRNA(NCBICodonTable):
    """Codon table for unambiguous RNA sequences."""

    nucleotide_alphabet = IUPACData.unambiguous_rna_letters


# ########  Deal with ambiguous forward translations


class AmbiguousCodonTable(CodonTable):
    """Base codon table for ambiguous sequences."""

    def __init__(
        self,
        codon_table,
        ambiguous_nucleotide_alphabet,
        ambiguous_nucleotide_values,
        ambiguous_protein_alphabet,
        ambiguous_protein_values,
    ):
        """Initialize the class."""
        CodonTable.__init__(
            self,
            ambiguous_nucleotide_alphabet,
            ambiguous_protein_alphabet,
            AmbiguousForwardTable(
                codon_table.forward_table,
                ambiguous_nucleotide_values,
                ambiguous_protein_values,
            ),
            codon_table.back_table,
            # These two are WRONG!  I need to get the
            # list of ambiguous codons which code for
            # the stop codons  XXX
            list_ambiguous_codons(
                codon_table.start_codons, ambiguous_nucleotide_values
            ),
            list_ambiguous_codons(codon_table.stop_codons, ambiguous_nucleotide_values),
        )
        self._codon_table = codon_table

    # Be sneaky and forward attribute lookups to the original table.
    # This lets us get the names, if the original table is an NCBI
    # table.
    def __getattr__(self, name):
        """Forward attribute lookups to the original table."""
        return getattr(self._codon_table, name)


def list_possible_proteins(codon, forward_table, ambiguous_nucleotide_values):
    """Return all possible encoded amino acids for ambiguous codon."""
    c1, c2, c3 = codon
    x1 = ambiguous_nucleotide_values[c1]
    x2 = ambiguous_nucleotide_values[c2]
    x3 = ambiguous_nucleotide_values[c3]
    possible = {}
    stops = []
    for y1 in x1:
        for y2 in x2:
            for y3 in x3:
                try:
                    possible[forward_table[y1 + y2 + y3]] = 1
                except KeyError:
                    # If tripping over a stop codon
                    stops.append(y1 + y2 + y3)
    if stops:
        if possible:
            raise TranslationError(
                "ambiguous codon %r codes for both proteins and stop codons" % codon
            )
        # This is a true stop codon - tell the caller about it
        raise KeyError(codon)
    return list(possible)


def list_ambiguous_codons(codons, ambiguous_nucleotide_values):
    """Extend a codon list to include all possible ambigous codons.

    e.g.::

         ['TAG', 'TAA'] -> ['TAG', 'TAA', 'TAR']
         ['UAG', 'UGA'] -> ['UAG', 'UGA', 'URA']

    Note that ['TAG', 'TGA'] -> ['TAG', 'TGA'], this does not add 'TRR'
    (which could also mean 'TAA' or 'TGG').
    Thus only two more codons are added in the following:

    e.g.::

        ['TGA', 'TAA', 'TAG'] -> ['TGA', 'TAA', 'TAG', 'TRA', 'TAR']

    Returns a new (longer) list of codon strings.
    """
    # Note ambiguous_nucleotide_values['R'] = 'AG' (etc)
    # This will generate things like 'TRR' from ['TAG', 'TGA'], which
    # we don't want to include:
    c1_list = sorted(
        letter
        for letter, meanings in ambiguous_nucleotide_values.items()
        if {codon[0] for codon in codons}.issuperset(set(meanings))
    )
    c2_list = sorted(
        letter
        for letter, meanings in ambiguous_nucleotide_values.items()
        if {codon[1] for codon in codons}.issuperset(set(meanings))
    )
    c3_list = sorted(
        letter
        for letter, meanings in ambiguous_nucleotide_values.items()
        if {codon[2] for codon in codons}.issuperset(set(meanings))
    )
    # candidates is a list (not a set) to preserve the iteration order
    candidates = []
    for c1 in c1_list:
        for c2 in c2_list:
            for c3 in c3_list:
                codon = c1 + c2 + c3
                if codon not in candidates and codon not in codons:
                    candidates.append(codon)
    answer = codons[:]  # copy
    # print "Have %i new candidates" % len(candidates)
    for ambig_codon in candidates:
        wanted = True
        # e.g. 'TRR' -> 'TAA', 'TAG', 'TGA', 'TGG'
        for codon in [
            c1 + c2 + c3
            for c1 in ambiguous_nucleotide_values[ambig_codon[0]]
            for c2 in ambiguous_nucleotide_values[ambig_codon[1]]
            for c3 in ambiguous_nucleotide_values[ambig_codon[2]]
        ]:
            if codon not in codons:
                # This ambiguous codon can code for a non-stop, exclude it!
                wanted = False
                # print "Rejecting %s" % ambig_codon
                continue
        if wanted:
            answer.append(ambig_codon)
    return answer


assert list_ambiguous_codons(["TGA", "TAA"], IUPACData.ambiguous_dna_values) == [
    "TGA",
    "TAA",
    "TRA",
]
assert list_ambiguous_codons(["TAG", "TGA"], IUPACData.ambiguous_dna_values) == [
    "TAG",
    "TGA",
]
assert list_ambiguous_codons(["TAG", "TAA"], IUPACData.ambiguous_dna_values) == [
    "TAG",
    "TAA",
    "TAR",
]
assert list_ambiguous_codons(["UAG", "UAA"], IUPACData.ambiguous_rna_values) == [
    "UAG",
    "UAA",
    "UAR",
]
assert list_ambiguous_codons(["TGA", "TAA", "TAG"], IUPACData.ambiguous_dna_values) == [
    "TGA",
    "TAA",
    "TAG",
    "TAR",
    "TRA",
]

# Forward translation is "onto", that is, any given codon always maps
# to the same protein, or it doesn't map at all.  Thus, I can build
# off of an existing table to produce the ambiguous mappings.
#
# This handles the general case.  Perhaps it's overkill?
#  >>> t = CodonTable.ambiguous_dna_by_id[1]
#  >>> t.forward_table["AAT"]
#  'N'
#  >>> t.forward_table["GAT"]
#  'D'
#  >>> t.forward_table["RAT"]
#  'B'
#  >>> t.forward_table["YTA"]
#  'L'


class AmbiguousForwardTable:
    """Forward table for translation of ambiguous nucleotide sequences."""

    def __init__(self, forward_table, ambiguous_nucleotide, ambiguous_protein):
        """Initialize the class."""
        self.forward_table = forward_table

        self.ambiguous_nucleotide = ambiguous_nucleotide
        self.ambiguous_protein = ambiguous_protein

        inverted = {}
        for name, val in ambiguous_protein.items():
            for c in val:
                x = inverted.get(c, {})
                x[name] = 1
                inverted[c] = x
        for name, val in inverted.items():
            inverted[name] = list(val)
        self._inverted = inverted

        self._cache = {}

    def __contains__(self, codon):
        """Check if codon works as key for ambiguous forward_table.

        Only returns 'True' if forward_table[codon] returns a value.
        """
        try:
            self.__getitem__(codon)
            return True
        except (KeyError, TranslationError):
            return False

    def get(self, codon, failobj=None):
        """Implement get for dictionary-like behaviour."""
        try:
            return self.__getitem__(codon)
        except KeyError:
            return failobj

    def __getitem__(self, codon):
        """Implement dictionary-like behaviour for AmbiguousForwardTable.

        forward_table[codon] will either return an amino acid letter,
        or throws a KeyError (if codon does not encode an amino acid)
        or a TranslationError (if codon does encode for an amino acid,
        but either is also a stop codon or does encode several amino acids,
        for which no unique letter is available in the given alphabet.
        """
        try:
            x = self._cache[codon]
        except KeyError:
            pass
        else:
            if x is TranslationError:
                raise TranslationError(codon)  # no unique translation
            if x is KeyError:
                raise KeyError(codon)  # it's a stop codon
            return x
        try:
            x = self.forward_table[codon]
            self._cache[codon] = x
            return x
        except KeyError:
            pass

        # XXX Need to make part of this into a method which returns
        # a list of all possible encodings for a codon!
        try:
            possible = list_possible_proteins(
                codon, self.forward_table, self.ambiguous_nucleotide
            )
        except KeyError:
            self._cache[codon] = KeyError
            raise KeyError(codon) from None  # stop codon
        except TranslationError:
            self._cache[codon] = TranslationError
            raise TranslationError(codon)  # does not code
        assert len(possible) > 0, "unambiguous codons must code"

        # Hah!  Only one possible protein, so use it
        if len(possible) == 1:
            self._cache[codon] = possible[0]
            return possible[0]

        # See if there's an ambiguous protein encoding for the multiples.
        # Find residues which exist in every coding set.
        ambiguous_possible = {}
        for amino in possible:
            for term in self._inverted[amino]:
                ambiguous_possible[term] = ambiguous_possible.get(term, 0) + 1

        n = len(possible)
        possible = []
        for amino, val in ambiguous_possible.items():
            if val == n:
                possible.append(amino)

        # No amino acid encoding for the results
        if len(possible) == 0:
            self._cache[codon] = TranslationError
            raise TranslationError(codon)  # no valid translation

        # All of these are valid, so choose one
        # To be unique, sort by smallet ambiguity then alphabetically
        # Can get this if "X" encodes for everything.
        # def _sort(x, y, table = self.ambiguous_protein):
        #    a = cmp(len(table[x]), len(table[y]))
        #    if a == 0:
        #        return cmp(x, y)
        #    return a

        # Sort by key is 2.x and 3.x compatible
        possible.sort(key=lambda x: (len(self.ambiguous_protein[x]), x))

        x = possible[0]
        self._cache[codon] = x
        return x


def register_ncbi_table(name, alt_name, id, table, start_codons, stop_codons):
    """Turn codon table data into objects (PRIVATE).

    The data is stored in the dictionaries.
    """
    # In most cases names are divided by "; ", however there is also
    # Table 11 'Bacterial, Archaeal and Plant Plastid Code', previously
    # 'Bacterial and Plant Plastid' which used to be just 'Bacterial'
    names = [
        x.strip() for x in name.replace(" and ", "; ").replace(", ", "; ").split("; ")
    ]

    dna = NCBICodonTableDNA(id, names + [alt_name], table, start_codons, stop_codons)
    ambig_dna = AmbiguousCodonTable(
        dna,
        IUPACData.ambiguous_dna_letters,
        IUPACData.ambiguous_dna_values,
        IUPACData.extended_protein_letters,
        IUPACData.extended_protein_values,
    )

    # replace all T's with U's for the RNA tables
    rna_table = {}
    generic_table = {}
    for codon, val in table.items():
        generic_table[codon] = val
        codon = codon.replace("T", "U")
        generic_table[codon] = val
        rna_table[codon] = val
    rna_start_codons = []
    generic_start_codons = []
    for codon in start_codons:
        generic_start_codons.append(codon)
        # We need to check if 'T' is in the codon, otherwise
        # generic_start_codons may contain duplicates
        if "T" in codon:
            codon = codon.replace("T", "U")
            generic_start_codons.append(codon)
        rna_start_codons.append(codon)
    rna_stop_codons = []
    generic_stop_codons = []
    for codon in stop_codons:
        generic_stop_codons.append(codon)
        if "T" in codon:
            codon = codon.replace("T", "U")
            generic_stop_codons.append(codon)
        rna_stop_codons.append(codon)

    generic = NCBICodonTable(
        id, names + [alt_name], generic_table, generic_start_codons, generic_stop_codons
    )

    # The following isn't very elegant, but seems to work nicely.
    _merged_values = dict(IUPACData.ambiguous_rna_values.items())
    _merged_values["T"] = "U"
    ambig_generic = AmbiguousCodonTable(
        generic,
        None,
        _merged_values,
        IUPACData.extended_protein_letters,
        IUPACData.extended_protein_values,
    )

    rna = NCBICodonTableRNA(
        id, names + [alt_name], rna_table, rna_start_codons, rna_stop_codons
    )

    ambig_rna = AmbiguousCodonTable(
        rna,
        IUPACData.ambiguous_rna_letters,
        IUPACData.ambiguous_rna_values,
        IUPACData.extended_protein_letters,
        IUPACData.extended_protein_values,
    )

    if id == 1:
        global standard_dna_table, standard_rna_table
        standard_dna_table = dna
        standard_rna_table = rna

    unambiguous_dna_by_id[id] = dna
    unambiguous_rna_by_id[id] = rna
    generic_by_id[id] = generic
    ambiguous_dna_by_id[id] = ambig_dna
    ambiguous_rna_by_id[id] = ambig_rna
    ambiguous_generic_by_id[id] = ambig_generic

    if alt_name is not None:
        names.append(alt_name)

    for name in names:
        unambiguous_dna_by_name[name] = dna
        unambiguous_rna_by_name[name] = rna
        generic_by_name[name] = generic
        ambiguous_dna_by_name[name] = ambig_dna
        ambiguous_rna_by_name[name] = ambig_rna
        ambiguous_generic_by_name[name] = ambig_generic


# The rest of this file is automatically generated, here we turn off
# black formatting in order to keep the codon tables compact.
#
# fmt: off

##########################################################################
# Start of auto-generated output from Scripts/update_ncbi_codon_table.py #
##########################################################################

# Data from NCBI genetic code table version 4.5

register_ncbi_table(
    name="Standard",
    alt_name="SGC0",
    id=1,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C",             "TGG": "W",   # noqa: E241
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG", "TGA"],
    start_codons=["TTG", "CTG", "ATG"],
)

register_ncbi_table(
    name="Vertebrate Mitochondrial",
    alt_name="SGC1",
    id=2,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S",                           # noqa: E241
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG", "AGA", "AGG"],
    start_codons=["ATT", "ATC", "ATA", "ATG", "GTG"],
)

register_ncbi_table(
    name="Yeast Mitochondrial",
    alt_name="SGC2",
    id=3,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "T", "CTC": "T", "CTA": "T", "CTG": "T",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG"],
    start_codons=["ATA", "ATG", "GTG"],
)

register_ncbi_table(
    name="Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate "
    "Mitochondrial; Mycoplasma; Spiroplasma",
    alt_name="SGC3",
    id=4,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG"],
    start_codons=["TTA", "TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"],
)

register_ncbi_table(
    name="Invertebrate Mitochondrial",
    alt_name="SGC4",
    id=5,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "S", "AGG": "S",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG"],
    start_codons=["TTG", "ATT", "ATC", "ATA", "ATG", "GTG"],
)

register_ncbi_table(
    name="Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear",
    alt_name="SGC5",
    id=6,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "Q", "TAG": "Q",
        "TGT": "C", "TGC": "C",             "TGG": "W",   # noqa: E241
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TGA"],
    start_codons=["ATG"],
)

register_ncbi_table(
    name="Echinoderm Mitochondrial; Flatworm Mitochondrial",
    alt_name="SGC8",
    id=9,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "N", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "S", "AGG": "S",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG"],
    start_codons=["ATG", "GTG"],
)

register_ncbi_table(
    name="Euplotid Nuclear",
    alt_name="SGC9",
    id=10,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "C", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG"],
    start_codons=["ATG"],
)

register_ncbi_table(
    name="Bacterial, Archaeal and Plant Plastid",
    alt_name=None,
    id=11,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C",             "TGG": "W",   # noqa: E241
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG", "TGA"],
    start_codons=["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"],
)

register_ncbi_table(
    name="Alternative Yeast Nuclear",
    alt_name=None,
    id=12,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C",             "TGG": "W",   # noqa: E241
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "S",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG", "TGA"],
    start_codons=["CTG", "ATG"],
)

register_ncbi_table(
    name="Ascidian Mitochondrial",
    alt_name=None,
    id=13,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "G", "AGG": "G",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG"],
    start_codons=["TTG", "ATA", "ATG", "GTG"],
)

register_ncbi_table(
    name="Alternative Flatworm Mitochondrial",
    alt_name=None,
    id=14,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "Y",               # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "N", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "S", "AGG": "S",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAG"],
    start_codons=["ATG"],
)

register_ncbi_table(
    name="Blepharisma Macronuclear",
    alt_name=None,
    id=15,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",             "TAG": "Q",   # noqa: E241
        "TGT": "C", "TGC": "C",             "TGG": "W",   # noqa: E241
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TGA"],
    start_codons=["ATG"],
)

register_ncbi_table(
    name="Chlorophycean Mitochondrial",
    alt_name=None,
    id=16,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",             "TAG": "L",   # noqa: E241
        "TGT": "C", "TGC": "C",             "TGG": "W",   # noqa: E241
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TGA"],
    start_codons=["ATG"],
)

register_ncbi_table(
    name="Trematode Mitochondrial",
    alt_name=None,
    id=21,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "N", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "S", "AGG": "S",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG"],
    start_codons=["ATG", "GTG"],
)

register_ncbi_table(
    name="Scenedesmus obliquus Mitochondrial",
    alt_name=None,
    id=22,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S",             "TCG": "S",   # noqa: E241
        "TAT": "Y", "TAC": "Y",             "TAG": "L",   # noqa: E241
        "TGT": "C", "TGC": "C",             "TGG": "W",   # noqa: E241
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TCA", "TAA", "TGA"],
    start_codons=["ATG"],
)

register_ncbi_table(
    name="Thraustochytrium Mitochondrial",
    alt_name=None,
    id=23,
    table={
        "TTT": "F", "TTC": "F",             "TTG": "L",   # noqa: E241
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C",             "TGG": "W",   # noqa: E241
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TTA", "TAA", "TAG", "TGA"],
    start_codons=["ATT", "ATG", "GTG"],
)

register_ncbi_table(
    name="Pterobranchia Mitochondrial",
    alt_name=None,
    id=24,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "S", "AGG": "K",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG"],
    start_codons=["TTG", "CTG", "ATG", "GTG"],
)

register_ncbi_table(
    name="Candidate Division SR1 and Gracilibacteria",
    alt_name=None,
    id=25,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "G", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG"],
    start_codons=["TTG", "ATG", "GTG"],
)

register_ncbi_table(
    name="Pachysolen tannophilus Nuclear",
    alt_name=None,
    id=26,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C",             "TGG": "W",   # noqa: E241
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "A",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG", "TGA"],
    start_codons=["CTG", "ATG"],
)

register_ncbi_table(
    name="Karyorelict Nuclear",
    alt_name=None,
    id=27,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "Q", "TAG": "Q",
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TGA"],
    start_codons=["ATG"],
)

register_ncbi_table(
    name="Condylostoma Nuclear",
    alt_name=None,
    id=28,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "Q", "TAG": "Q",
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG", "TGA"],
    start_codons=["ATG"],
)

register_ncbi_table(
    name="Mesodinium Nuclear",
    alt_name=None,
    id=29,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "Y", "TAG": "Y",
        "TGT": "C", "TGC": "C",             "TGG": "W",   # noqa: E241
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TGA"],
    start_codons=["ATG"],
)

register_ncbi_table(
    name="Peritrich Nuclear",
    alt_name=None,
    id=30,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "E", "TAG": "E",
        "TGT": "C", "TGC": "C",             "TGG": "W",   # noqa: E241
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TGA"],
    start_codons=["ATG"],
)

register_ncbi_table(
    name="Blastocrithidia Nuclear",
    alt_name=None,
    id=31,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "E", "TAG": "E",
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TAG"],
    start_codons=["ATG"],
)

register_ncbi_table(
    name="Balanophoraceae Plastid",
    alt_name=None,
    id=32,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",             "TAG": "W",   # noqa: E241
        "TGT": "C", "TGC": "C",             "TGG": "W",   # noqa: E241
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAA", "TGA"],
    start_codons=["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"],
)

register_ncbi_table(
    name="Cephalodiscidae Mitochondrial",
    alt_name=None,
    id=33,
    table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "Y",               # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "S", "AGG": "K",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    },
    stop_codons=["TAG"],
    start_codons=["TTG", "CTG", "ATG", "GTG"],
)

########################################################################
# End of auto-generated output from Scripts/update_ncbi_codon_table.py #
########################################################################
