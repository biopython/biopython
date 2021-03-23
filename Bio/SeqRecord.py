# Copyright 2000-2002 Andrew Dalke.  All rights reserved.
# Copyright 2002-2004 Brad Chapman.  All rights reserved.
# Copyright 2006-2020 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Represent a Sequence Record, a sequence with annotation."""
# NEEDS TO BE SYNCH WITH THE REST OF BIOPYTHON AND BIOPERL
# In particular, the SeqRecord and BioSQL.BioSeq.DBSeqRecord classes
# need to be in sync (this is the BioSQL "Database SeqRecord").
from io import StringIO

from Bio import StreamModeError
from Bio.Seq import UndefinedSequenceError


_NO_SEQRECORD_COMPARISON = "SeqRecord comparison is deliberately not implemented. Explicitly compare the attributes of interest."


class _RestrictedDict(dict):
    """Dict which only allows sequences of given length as values (PRIVATE).

    This simple subclass of the Python dictionary is used in the SeqRecord
    object for holding per-letter-annotations.  This class is intended to
    prevent simple errors by only allowing python sequences (e.g. lists,
    strings and tuples) to be stored, and only if their length matches that
    expected (the length of the SeqRecord's seq object).  It cannot however
    prevent the entries being edited in situ (for example appending entries
    to a list).

    >>> x = _RestrictedDict(5)
    >>> x["test"] = "hello"
    >>> x
    {'test': 'hello'}

    Adding entries which don't have the expected length are blocked:

    >>> x["test"] = "hello world"
    Traceback (most recent call last):
    ...
    TypeError: We only allow python sequences (lists, tuples or strings) of length 5.

    The expected length is stored as a private attribute,

    >>> x._length
    5

    In order that the SeqRecord (and other objects using this class) can be
    pickled, for example for use in the multiprocessing library, we need to
    be able to pickle the restricted dictionary objects.

    Using the default protocol, which is 3 on Python 3,

    >>> import pickle
    >>> y = pickle.loads(pickle.dumps(x))
    >>> y
    {'test': 'hello'}
    >>> y._length
    5

    Using the highest protocol, which is 4 on Python 3,

    >>> import pickle
    >>> z = pickle.loads(pickle.dumps(x, pickle.HIGHEST_PROTOCOL))
    >>> z
    {'test': 'hello'}
    >>> z._length
    5
    """

    def __init__(self, length):
        """Create an EMPTY restricted dictionary."""
        dict.__init__(self)
        self._length = int(length)

    def __setitem__(self, key, value):
        # The check hasattr(self, "_length") is to cope with pickle protocol 2
        # I couldn't seem to avoid this with __getstate__ and __setstate__
        if (
            not hasattr(value, "__len__")
            or not hasattr(value, "__getitem__")
            or (hasattr(self, "_length") and len(value) != self._length)
        ):
            raise TypeError(
                "We only allow python sequences (lists, tuples or strings) "
                f"of length {self._length}."
            )
        dict.__setitem__(self, key, value)

    def update(self, new_dict):
        # Force this to go via our strict __setitem__ method
        for (key, value) in new_dict.items():
            self[key] = value


class SeqRecord:
    """A SeqRecord object holds a sequence and information about it.

    Main attributes:
     - id          - Identifier such as a locus tag (string)
     - seq         - The sequence itself (Seq object or similar)

    Additional attributes:
     - name        - Sequence name, e.g. gene name (string)
     - description - Additional text (string)
     - dbxrefs     - List of database cross references (list of strings)
     - features    - Any (sub)features defined (list of SeqFeature objects)
     - annotations - Further information about the whole sequence (dictionary).
       Most entries are strings, or lists of strings.
     - letter_annotations - Per letter/symbol annotation (restricted
       dictionary). This holds Python sequences (lists, strings
       or tuples) whose length matches that of the sequence.
       A typical use would be to hold a list of integers
       representing sequencing quality scores, or a string
       representing the secondary structure.

    You will typically use Bio.SeqIO to read in sequences from files as
    SeqRecord objects.  However, you may want to create your own SeqRecord
    objects directly (see the __init__ method for further details):

    >>> from Bio.Seq import Seq
    >>> from Bio.SeqRecord import SeqRecord
    >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
    ...                    id="YP_025292.1", name="HokC",
    ...                    description="toxic membrane protein")
    >>> print(record)
    ID: YP_025292.1
    Name: HokC
    Description: toxic membrane protein
    Number of features: 0
    Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')

    If you want to save SeqRecord objects to a sequence file, use Bio.SeqIO
    for this.  For the special case where you want the SeqRecord turned into
    a string in a particular file format there is a format method which uses
    Bio.SeqIO internally:

    >>> print(record.format("fasta"))
    >YP_025292.1 toxic membrane protein
    MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
    <BLANKLINE>

    You can also do things like slicing a SeqRecord, checking its length, etc

    >>> len(record)
    44
    >>> edited = record[:10] + record[11:]
    >>> print(edited.seq)
    MKQHKAMIVAIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
    >>> print(record.seq)
    MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF

    """

    def __init__(
        self,
        seq,
        id="<unknown id>",
        name="<unknown name>",
        description="<unknown description>",
        dbxrefs=None,
        features=None,
        annotations=None,
        letter_annotations=None,
    ):
        """Create a SeqRecord.

        Arguments:
         - seq         - Sequence, required (Seq or MutableSeq)
         - id          - Sequence identifier, recommended (string)
         - name        - Sequence name, optional (string)
         - description - Sequence description, optional (string)
         - dbxrefs     - Database cross references, optional (list of strings)
         - features    - Any (sub)features, optional (list of SeqFeature objects)
         - annotations - Dictionary of annotations for the whole sequence
         - letter_annotations - Dictionary of per-letter-annotations, values
           should be strings, list or tuples of the same length as the full
           sequence.

        You will typically use Bio.SeqIO to read in sequences from files as
        SeqRecord objects.  However, you may want to create your own SeqRecord
        objects directly.

        Note that while an id is optional, we strongly recommend you supply a
        unique id string for each record.  This is especially important
        if you wish to write your sequences to a file.

        You can create a 'blank' SeqRecord object, and then populate the
        attributes later.
        """
        if id is not None and not isinstance(id, str):
            # Lots of existing code uses id=None... this may be a bad idea.
            raise TypeError("id argument should be a string")
        if not isinstance(name, str):
            raise TypeError("name argument should be a string")
        if not isinstance(description, str):
            raise TypeError("description argument should be a string")
        self._seq = seq
        self.id = id
        self.name = name
        self.description = description

        # database cross references (for the whole sequence)
        if dbxrefs is None:
            dbxrefs = []
        elif not isinstance(dbxrefs, list):
            raise TypeError("dbxrefs argument should be a list (of strings)")
        self.dbxrefs = dbxrefs

        # annotations about the whole sequence
        if annotations is None:
            annotations = {}
        elif not isinstance(annotations, dict):
            raise TypeError("annotations argument should be a dict")
        self.annotations = annotations

        if letter_annotations is None:
            # annotations about each letter in the sequence
            if seq is None:
                # Should we allow this and use a normal unrestricted dict?
                self._per_letter_annotations = _RestrictedDict(length=0)
            else:
                try:
                    self._per_letter_annotations = _RestrictedDict(length=len(seq))
                except TypeError:
                    raise TypeError(
                        "seq argument should be a Seq object or similar"
                    ) from None
        else:
            # This will be handled via the property set function, which will
            # turn this into a _RestrictedDict and thus ensure all the values
            # in the dict are the right length
            self.letter_annotations = letter_annotations

        # annotations about parts of the sequence
        if features is None:
            features = []
        elif not isinstance(features, list):
            raise TypeError(
                "features argument should be a list (of SeqFeature objects)"
            )
        self.features = features

    # TODO - Just make this a read only property?
    def _set_per_letter_annotations(self, value):
        if not isinstance(value, dict):
            raise TypeError(
                "The per-letter-annotations should be a (restricted) dictionary."
            )
        # Turn this into a restricted-dictionary (and check the entries)
        try:
            self._per_letter_annotations = _RestrictedDict(length=len(self.seq))
        except AttributeError:
            # e.g. seq is None
            self._per_letter_annotations = _RestrictedDict(length=0)
        self._per_letter_annotations.update(value)

    letter_annotations = property(
        fget=lambda self: self._per_letter_annotations,
        fset=_set_per_letter_annotations,
        doc="""Dictionary of per-letter-annotation for the sequence.

        For example, this can hold quality scores used in FASTQ or QUAL files.
        Consider this example using Bio.SeqIO to read in an example Solexa
        variant FASTQ file as a SeqRecord:

        >>> from Bio import SeqIO
        >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
        >>> print("%s %s" % (record.id, record.seq))
        slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
        >>> print(list(record.letter_annotations))
        ['solexa_quality']
        >>> print(record.letter_annotations["solexa_quality"])
        [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5]

        The letter_annotations get sliced automatically if you slice the
        parent SeqRecord, for example taking the last ten bases:

        >>> sub_record = record[-10:]
        >>> print("%s %s" % (sub_record.id, sub_record.seq))
        slxa_0001_1_0001_01 ACGTNNNNNN
        >>> print(sub_record.letter_annotations["solexa_quality"])
        [4, 3, 2, 1, 0, -1, -2, -3, -4, -5]

        Any python sequence (i.e. list, tuple or string) can be recorded in
        the SeqRecord's letter_annotations dictionary as long as the length
        matches that of the SeqRecord's sequence.  e.g.

        >>> len(sub_record.letter_annotations)
        1
        >>> sub_record.letter_annotations["dummy"] = "abcdefghij"
        >>> len(sub_record.letter_annotations)
        2

        You can delete entries from the letter_annotations dictionary as usual:

        >>> del sub_record.letter_annotations["solexa_quality"]
        >>> sub_record.letter_annotations
        {'dummy': 'abcdefghij'}

        You can completely clear the dictionary easily as follows:

        >>> sub_record.letter_annotations = {}
        >>> sub_record.letter_annotations
        {}

        Note that if replacing the record's sequence with a sequence of a
        different length you must first clear the letter_annotations dict.
        """,
    )

    def _set_seq(self, value):
        # TODO - Add a deprecation warning that the seq should be write only?
        if self._per_letter_annotations:
            if len(self) != len(value):
                # TODO - Make this a warning? Silently empty the dictionary?
                raise ValueError("You must empty the letter annotations first!")
            else:
                # Leave the existing per letter annotations unchanged:
                self._seq = value
        else:
            self._seq = value
            # Reset the (empty) letter annotations dict with new length:
            try:
                self._per_letter_annotations = _RestrictedDict(length=len(self.seq))
            except AttributeError:
                # e.g. seq is None
                self._per_letter_annotations = _RestrictedDict(length=0)

    seq = property(
        fget=lambda self: self._seq,
        fset=_set_seq,
        doc="The sequence itself, as a Seq or MutableSeq object.",
    )

    def __getitem__(self, index):
        """Return a sub-sequence or an individual letter.

        Slicing, e.g. my_record[5:10], returns a new SeqRecord for
        that sub-sequence with some annotation preserved as follows:

        * The name, id and description are kept as-is.
        * Any per-letter-annotations are sliced to match the requested
          sub-sequence.
        * Unless a stride is used, all those features which fall fully
          within the subsequence are included (with their locations
          adjusted accordingly). If you want to preserve any truncated
          features (e.g. GenBank/EMBL source features), you must
          explicitly add them to the new SeqRecord yourself.
        * With the exception of any molecule type, the annotations
          dictionary and the dbxrefs list are not used for the new
          SeqRecord, as in general they may not apply to the
          subsequence. If you want to preserve them, you must explicitly
          copy them to the new SeqRecord yourself.

        Using an integer index, e.g. my_record[5] is shorthand for
        extracting that letter from the sequence, my_record.seq[5].

        For example, consider this short protein and its secondary
        structure as encoded by the PDB (e.g. H for alpha helices),
        plus a simple feature for its histidine self phosphorylation
        site:

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> rec = SeqRecord(Seq("MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLAT"
        ...                     "EMMSEQDGYLAESINKDIEECNAIIEQFIDYLR"),
        ...                 id="1JOY", name="EnvZ",
        ...                 description="Homodimeric domain of EnvZ from E. coli")
        >>> rec.letter_annotations["secondary_structure"] = "  S  SSSSSSHHHHHTTTHHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHHHHHHHHHHHHHHTT  "
        >>> rec.features.append(SeqFeature(FeatureLocation(20, 21),
        ...                     type = "Site"))

        Now let's have a quick look at the full record,

        >>> print(rec)
        ID: 1JOY
        Name: EnvZ
        Description: Homodimeric domain of EnvZ from E. coli
        Number of features: 1
        Per letter annotation for: secondary_structure
        Seq('MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEE...YLR')
        >>> rec.letter_annotations["secondary_structure"]
        '  S  SSSSSSHHHHHTTTHHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHHHHHHHHHHHHHHTT  '
        >>> print(rec.features[0].location)
        [20:21]

        Now let's take a sub sequence, here chosen as the first (fractured)
        alpha helix which includes the histidine phosphorylation site:

        >>> sub = rec[11:41]
        >>> print(sub)
        ID: 1JOY
        Name: EnvZ
        Description: Homodimeric domain of EnvZ from E. coli
        Number of features: 1
        Per letter annotation for: secondary_structure
        Seq('RTLLMAGVSHDLRTPLTRIRLATEMMSEQD')
        >>> sub.letter_annotations["secondary_structure"]
        'HHHHHTTTHHHHHHHHHHHHHHHHHHHHHH'
        >>> print(sub.features[0].location)
        [9:10]

        You can also of course omit the start or end values, for
        example to get the first ten letters only:

        >>> print(rec[:10])
        ID: 1JOY
        Name: EnvZ
        Description: Homodimeric domain of EnvZ from E. coli
        Number of features: 0
        Per letter annotation for: secondary_structure
        Seq('MAAGVKQLAD')

        Or for the last ten letters:

        >>> print(rec[-10:])
        ID: 1JOY
        Name: EnvZ
        Description: Homodimeric domain of EnvZ from E. coli
        Number of features: 0
        Per letter annotation for: secondary_structure
        Seq('IIEQFIDYLR')

        If you omit both, then you get a copy of the original record (although
        lacking the annotations and dbxrefs):

        >>> print(rec[:])
        ID: 1JOY
        Name: EnvZ
        Description: Homodimeric domain of EnvZ from E. coli
        Number of features: 1
        Per letter annotation for: secondary_structure
        Seq('MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEE...YLR')

        Finally, indexing with a simple integer is shorthand for pulling out
        that letter from the sequence directly:

        >>> rec[5]
        'K'
        >>> rec.seq[5]
        'K'
        """
        if isinstance(index, int):
            # NOTE - The sequence level annotation like the id, name, etc
            # do not really apply to a single character.  However, should
            # we try and expose any per-letter-annotation here?  If so how?
            return self.seq[index]
        elif isinstance(index, slice):
            if self.seq is None:
                raise ValueError("If the sequence is None, we cannot slice it.")
            parent_length = len(self)
            try:
                from BioSQL.BioSeq import DBSeqRecord

                biosql_available = True
            except ImportError:
                biosql_available = False

            if biosql_available and isinstance(self, DBSeqRecord):
                answer = SeqRecord(
                    self.seq[index],
                    id=self.id,
                    name=self.name,
                    description=self.description,
                )
            else:
                answer = self.__class__(
                    self.seq[index],
                    id=self.id,
                    name=self.name,
                    description=self.description,
                )
            # TODO - The description may no longer apply.
            # It would be safer to change it to something
            # generic like "edited" or the default value.

            # Don't copy the annotation dict and dbxefs list,
            # they may not apply to a subsequence.
            # answer.annotations = dict(self.annotations.items())
            # answer.dbxrefs = self.dbxrefs[:]
            # TODO - Review this in light of adding SeqRecord objects?

            if "molecule_type" in self.annotations:
                # This will still apply, and we need it for GenBank/EMBL etc output
                answer.annotations["molecule_type"] = self.annotations["molecule_type"]

            # TODO - Cope with strides by generating ambiguous locations?
            start, stop, step = index.indices(parent_length)
            if step == 1:
                # Select relevant features, add them with shifted locations
                # assert str(self.seq)[index] == str(self.seq)[start:stop]
                for f in self.features:
                    if f.ref or f.ref_db:
                        # TODO - Implement this (with lots of tests)?
                        import warnings

                        warnings.warn(
                            "When slicing SeqRecord objects, any "
                            "SeqFeature referencing other sequences (e.g. "
                            "from segmented GenBank records) are ignored."
                        )
                        continue
                    if (
                        start <= f.location.nofuzzy_start
                        and f.location.nofuzzy_end <= stop
                    ):
                        answer.features.append(f._shift(-start))

            # Slice all the values to match the sliced sequence
            # (this should also work with strides, even negative strides):
            for key, value in self.letter_annotations.items():
                answer._per_letter_annotations[key] = value[index]

            return answer
        raise ValueError("Invalid index")

    def __iter__(self):
        """Iterate over the letters in the sequence.

        For example, using Bio.SeqIO to read in a protein FASTA file:

        >>> from Bio import SeqIO
        >>> record = SeqIO.read("Fasta/loveliesbleeding.pro", "fasta")
        >>> for amino in record:
        ...     print(amino)
        ...     if amino == "L": break
        X
        A
        G
        L
        >>> print(record.seq[3])
        L

        This is just a shortcut for iterating over the sequence directly:

        >>> for amino in record.seq:
        ...     print(amino)
        ...     if amino == "L": break
        X
        A
        G
        L
        >>> print(record.seq[3])
        L

        Note that this does not facilitate iteration together with any
        per-letter-annotation.  However, you can achieve that using the
        python zip function on the record (or its sequence) and the relevant
        per-letter-annotation:

        >>> from Bio import SeqIO
        >>> rec = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
        >>> print("%s %s" % (rec.id, rec.seq))
        slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
        >>> print(list(rec.letter_annotations))
        ['solexa_quality']
        >>> for nuc, qual in zip(rec, rec.letter_annotations["solexa_quality"]):
        ...     if qual > 35:
        ...         print("%s %i" % (nuc, qual))
        A 40
        C 39
        G 38
        T 37
        A 36

        You may agree that using zip(rec.seq, ...) is more explicit than using
        zip(rec, ...) as shown above.
        """
        return iter(self.seq)

    def __contains__(self, char):
        """Implement the 'in' keyword, searches the sequence.

        e.g.

        >>> from Bio import SeqIO
        >>> record = SeqIO.read("Fasta/sweetpea.nu", "fasta")
        >>> "GAATTC" in record
        False
        >>> "AAA" in record
        True

        This essentially acts as a proxy for using "in" on the sequence:

        >>> "GAATTC" in record.seq
        False
        >>> "AAA" in record.seq
        True

        Note that you can also use Seq objects as the query,

        >>> from Bio.Seq import Seq
        >>> Seq("AAA") in record
        True

        See also the Seq object's __contains__ method.
        """
        return char in self.seq

    def __str__(self):
        """Return a human readable summary of the record and its annotation (string).

        The python built in function str works by calling the object's ___str__
        method.  e.g.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
        ...                    id="YP_025292.1", name="HokC",
        ...                    description="toxic membrane protein, small")
        >>> print(str(record))
        ID: YP_025292.1
        Name: HokC
        Description: toxic membrane protein, small
        Number of features: 0
        Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')

        In this example you don't actually need to call str explicity, as the
        print command does this automatically:

        >>> print(record)
        ID: YP_025292.1
        Name: HokC
        Description: toxic membrane protein, small
        Number of features: 0
        Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')

        Note that long sequences are shown truncated.
        """
        lines = []
        if self.id:
            lines.append(f"ID: {self.id}")
        if self.name:
            lines.append(f"Name: {self.name}")
        if self.description:
            lines.append(f"Description: {self.description}")
        if self.dbxrefs:
            lines.append("Database cross-references: " + ", ".join(self.dbxrefs))
        lines.append(f"Number of features: {len(self.features)}")
        for a in self.annotations:
            lines.append(f"/{a}={str(self.annotations[a])}")
        if self.letter_annotations:
            lines.append(
                "Per letter annotation for: " + ", ".join(self.letter_annotations)
            )
        try:
            bytes(self.seq)
        except UndefinedSequenceError:
            lines.append(f"Undefined sequence of length {len(self.seq)}")
        else:
            # Don't want to include the entire sequence
            seq = repr(self.seq)
            lines.append(seq)
        return "\n".join(lines)

    def __repr__(self):
        """Return a concise summary of the record for debugging (string).

        The python built in function repr works by calling the object's ___repr__
        method.  e.g.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> rec = SeqRecord(Seq("MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKAT"
        ...                     "GEMKEQTEWHRVVLFGKLAEVASEYLRKGSQVYIEGQLRTRKWTDQ"
        ...                     "SGQDRYTTEVVVNVGGTMQMLGGRQGGGAPAGGNIGGGQPQGGWGQ"
        ...                     "PQQPQGGNQFSGGAQSRPQQSAPAAPSNEPPMDFDDDIPF"),
        ...                 id="NP_418483.1", name="b4059",
        ...                 description="ssDNA-binding protein",
        ...                 dbxrefs=["ASAP:13298", "GI:16131885", "GeneID:948570"])
        >>> print(repr(rec))
        SeqRecord(seq=Seq('MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKATGEMKEQTE...IPF'), id='NP_418483.1', name='b4059', description='ssDNA-binding protein', dbxrefs=['ASAP:13298', 'GI:16131885', 'GeneID:948570'])

        At the python prompt you can also use this shorthand:

        >>> rec
        SeqRecord(seq=Seq('MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKATGEMKEQTE...IPF'), id='NP_418483.1', name='b4059', description='ssDNA-binding protein', dbxrefs=['ASAP:13298', 'GI:16131885', 'GeneID:948570'])

        Note that long sequences are shown truncated. Also note that any
        annotations, letter_annotations and features are not shown (as they
        would lead to a very long string).
        """
        return (
            f"{self.__class__.__name__}(seq={self.seq!r}, id={self.id!r},"
            f" name={self.name!r}, description={self.description!r},"
            f" dbxrefs={self.dbxrefs!r})"
        )

    def format(self, format):
        r"""Return the record as a string in the specified file format.

        The format should be a lower case string supported as an output
        format by Bio.SeqIO, which is used to turn the SeqRecord into a
        string.  e.g.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
        ...                    id="YP_025292.1", name="HokC",
        ...                    description="toxic membrane protein")
        >>> record.format("fasta")
        '>YP_025292.1 toxic membrane protein\nMKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF\n'
        >>> print(record.format("fasta"))
        >YP_025292.1 toxic membrane protein
        MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
        <BLANKLINE>

        The Python print function automatically appends a new line, meaning
        in this example a blank line is shown.  If you look at the string
        representation you can see there is a trailing new line (shown as
        slash n) which is important when writing to a file or if
        concatenating multiple sequence strings together.

        Note that this method will NOT work on every possible file format
        supported by Bio.SeqIO (e.g. some are for multiple sequences only,
        and binary formats are not supported).
        """
        # See also the __format__ method
        # See also the Bio.Align.Generic.Alignment class and its format()
        return self.__format__(format)

    def __format__(self, format_spec):
        r"""Return the record as a string in the specified file format.

        This method supports the Python format() function and f-strings.
        The format_spec should be a lower case string supported by
        Bio.SeqIO as a text output file format. Requesting a binary file
        format raises a ValueError. e.g.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
        ...                    id="YP_025292.1", name="HokC",
        ...                    description="toxic membrane protein")
        ...
        >>> format(record, "fasta")
        '>YP_025292.1 toxic membrane protein\nMKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF\n'
        >>> print(f"Here is {record.id} in FASTA format:\n{record:fasta}")
        Here is YP_025292.1 in FASTA format:
        >YP_025292.1 toxic membrane protein
        MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
        <BLANKLINE>

        See also the SeqRecord's format() method.
        """
        if not format_spec:
            # Follow python convention and default to using __str__
            return str(self)
        from Bio import SeqIO

        # Easy case, can call string-building function directly
        if format_spec in SeqIO._FormatToString:
            return SeqIO._FormatToString[format_spec](self)

        # Harder case, make a temp handle instead
        handle = StringIO()
        try:
            SeqIO.write(self, handle, format_spec)
        except StreamModeError:
            raise ValueError(
                "Binary format %s cannot be used with SeqRecord format method"
                % format_spec
            ) from None
        return handle.getvalue()

    def __len__(self):
        """Return the length of the sequence.

        For example, using Bio.SeqIO to read in a FASTA nucleotide file:

        >>> from Bio import SeqIO
        >>> record = SeqIO.read("Fasta/sweetpea.nu", "fasta")
        >>> len(record)
        309
        >>> len(record.seq)
        309
        """
        return len(self.seq)

    def __lt__(self, other):
        """Define the less-than operand (not implemented)."""
        raise NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __le___(self, other):
        """Define the less-than-or-equal-to operand (not implemented)."""
        raise NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __eq__(self, other):
        """Define the equal-to operand (not implemented)."""
        raise NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __ne__(self, other):
        """Define the not-equal-to operand (not implemented)."""
        raise NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __gt__(self, other):
        """Define the greater-than operand (not implemented)."""
        raise NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __ge__(self, other):
        """Define the greater-than-or-equal-to operand (not implemented)."""
        raise NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __bool__(self):
        """Boolean value of an instance of this class (True).

        This behaviour is for backwards compatibility, since until the
        __len__ method was added, a SeqRecord always evaluated as True.

        Note that in comparison, a Seq object will evaluate to False if it
        has a zero length sequence.

        WARNING: The SeqRecord may in future evaluate to False when its
        sequence is of zero length (in order to better match the Seq
        object behaviour)!
        """
        return True

    def __add__(self, other):
        """Add another sequence or string to this sequence.

        The other sequence can be a SeqRecord object, a Seq object (or
        similar, e.g. a MutableSeq) or a plain Python string. If you add
        a plain string or a Seq (like) object, the new SeqRecord will simply
        have this appended to the existing data. However, any per letter
        annotation will be lost:

        >>> from Bio import SeqIO
        >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
        >>> print("%s %s" % (record.id, record.seq))
        slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
        >>> print(list(record.letter_annotations))
        ['solexa_quality']

        >>> new = record + "ACT"
        >>> print("%s %s" % (new.id, new.seq))
        slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNNACT
        >>> print(list(new.letter_annotations))
        []

        The new record will attempt to combine the annotation, but for any
        ambiguities (e.g. different names) it defaults to omitting that
        annotation.

        >>> from Bio import SeqIO
        >>> with open("GenBank/pBAD30.gb") as handle:
        ...     plasmid = SeqIO.read(handle, "gb")
        >>> print("%s %i" % (plasmid.id, len(plasmid)))
        pBAD30 4923

        Now let's cut the plasmid into two pieces, and join them back up the
        other way round (i.e. shift the starting point on this plasmid, have
        a look at the annotated features in the original file to see why this
        particular split point might make sense):

        >>> left = plasmid[:3765]
        >>> right = plasmid[3765:]
        >>> new = right + left
        >>> print("%s %i" % (new.id, len(new)))
        pBAD30 4923
        >>> str(new.seq) == str(right.seq + left.seq)
        True
        >>> len(new.features) == len(left.features) + len(right.features)
        True

        When we add the left and right SeqRecord objects, their annotation
        is all consistent, so it is all conserved in the new SeqRecord:

        >>> new.id == left.id == right.id == plasmid.id
        True
        >>> new.name == left.name == right.name == plasmid.name
        True
        >>> new.description == plasmid.description
        True
        >>> new.annotations == left.annotations == right.annotations
        True
        >>> new.letter_annotations == plasmid.letter_annotations
        True
        >>> new.dbxrefs == left.dbxrefs == right.dbxrefs
        True

        However, we should point out that when we sliced the SeqRecord,
        any annotations dictionary or dbxrefs list entries were lost.
        You can explicitly copy them like this:

        >>> new.annotations = plasmid.annotations.copy()
        >>> new.dbxrefs = plasmid.dbxrefs[:]
        """
        if not isinstance(other, SeqRecord):
            # Assume it is a string or a Seq.
            # Note can't transfer any per-letter-annotations
            return SeqRecord(
                self.seq + other,
                id=self.id,
                name=self.name,
                description=self.description,
                features=self.features[:],
                annotations=self.annotations.copy(),
                dbxrefs=self.dbxrefs[:],
            )
        # Adding two SeqRecord objects... must merge annotation.
        answer = SeqRecord(
            self.seq + other.seq, features=self.features[:], dbxrefs=self.dbxrefs[:]
        )
        # Will take all the features and all the db cross refs,
        length = len(self)
        for f in other.features:
            answer.features.append(f._shift(length))
        del length
        for ref in other.dbxrefs:
            if ref not in answer.dbxrefs:
                answer.dbxrefs.append(ref)
        # Take common id/name/description/annotation
        if self.id == other.id:
            answer.id = self.id
        if self.name == other.name:
            answer.name = self.name
        if self.description == other.description:
            answer.description = self.description
        for k, v in self.annotations.items():
            if k in other.annotations and other.annotations[k] == v:
                answer.annotations[k] = v
        # Can append matching per-letter-annotation
        for k, v in self.letter_annotations.items():
            if k in other.letter_annotations:
                answer.letter_annotations[k] = v + other.letter_annotations[k]
        return answer

    def __radd__(self, other):
        """Add another sequence or string to this sequence (from the left).

        This method handles adding a Seq object (or similar, e.g. MutableSeq)
        or a plain Python string (on the left) to a SeqRecord (on the right).
        See the __add__ method for more details, but for example:

        >>> from Bio import SeqIO
        >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
        >>> print("%s %s" % (record.id, record.seq))
        slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
        >>> print(list(record.letter_annotations))
        ['solexa_quality']

        >>> new = "ACT" + record
        >>> print("%s %s" % (new.id, new.seq))
        slxa_0001_1_0001_01 ACTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
        >>> print(list(new.letter_annotations))
        []
        """
        if isinstance(other, SeqRecord):
            raise RuntimeError(
                "This should have happened via the __add__ of "
                "the other SeqRecord being added!"
            )
        # Assume it is a string or a Seq.
        # Note can't transfer any per-letter-annotations
        offset = len(other)
        return SeqRecord(
            other + self.seq,
            id=self.id,
            name=self.name,
            description=self.description,
            features=[f._shift(offset) for f in self.features],
            annotations=self.annotations.copy(),
            dbxrefs=self.dbxrefs[:],
        )

    def upper(self):
        """Return a copy of the record with an upper case sequence.

        All the annotation is preserved unchanged. e.g.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> record = SeqRecord(Seq("acgtACGT"), id="Test",
        ...                    description = "Made up for this example")
        >>> record.letter_annotations["phred_quality"] = [1, 2, 3, 4, 5, 6, 7, 8]
        >>> print(record.upper().format("fastq"))
        @Test Made up for this example
        ACGTACGT
        +
        "#$%&'()
        <BLANKLINE>

        Naturally, there is a matching lower method:

        >>> print(record.lower().format("fastq"))
        @Test Made up for this example
        acgtacgt
        +
        "#$%&'()
        <BLANKLINE>
        """
        return SeqRecord(
            self.seq.upper(),
            id=self.id,
            name=self.name,
            description=self.description,
            dbxrefs=self.dbxrefs[:],
            features=self.features[:],
            annotations=self.annotations.copy(),
            letter_annotations=self.letter_annotations.copy(),
        )

    def lower(self):
        """Return a copy of the record with a lower case sequence.

        All the annotation is preserved unchanged. e.g.

        >>> from Bio import SeqIO
        >>> record = SeqIO.read("Fasta/aster.pro", "fasta")
        >>> print(record.format("fasta"))
        >gi|3298468|dbj|BAA31520.1| SAMIPF
        GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVG
        VTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI
        <BLANKLINE>
        >>> print(record.lower().format("fasta"))
        >gi|3298468|dbj|BAA31520.1| SAMIPF
        gghvnpavtfgafvggnitllrgivyiiaqllgstvaclllkfvtndmavgvfslsagvg
        vtnalvfeivmtfglvytvyataidpkkgslgtiapiaigfivgani
        <BLANKLINE>

        To take a more annotation rich example,

        >>> from Bio import SeqIO
        >>> old = SeqIO.read("EMBL/TRBG361.embl", "embl")
        >>> len(old.features)
        3
        >>> new = old.lower()
        >>> len(old.features) == len(new.features)
        True
        >>> old.annotations["organism"] == new.annotations["organism"]
        True
        >>> old.dbxrefs == new.dbxrefs
        True
        """
        return SeqRecord(
            self.seq.lower(),
            id=self.id,
            name=self.name,
            description=self.description,
            dbxrefs=self.dbxrefs[:],
            features=self.features[:],
            annotations=self.annotations.copy(),
            letter_annotations=self.letter_annotations.copy(),
        )

    def reverse_complement(
        self,
        id=False,
        name=False,
        description=False,
        features=True,
        annotations=False,
        letter_annotations=True,
        dbxrefs=False,
    ):
        """Return new SeqRecord with reverse complement sequence.

        By default the new record does NOT preserve the sequence identifier,
        name, description, general annotation or database cross-references -
        these are unlikely to apply to the reversed sequence.

        You can specify the returned record's id, name and description as
        strings, or True to keep that of the parent, or False for a default.

        You can specify the returned record's features with a list of
        SeqFeature objects, or True to keep that of the parent, or False to
        omit them. The default is to keep the original features (with the
        strand and locations adjusted).

        You can also specify both the returned record's annotations and
        letter_annotations as dictionaries, True to keep that of the parent,
        or False to omit them. The default is to keep the original
        annotations (with the letter annotations reversed).

        To show what happens to the pre-letter annotations, consider an
        example Solexa variant FASTQ file with a single entry, which we'll
        read in as a SeqRecord:

        >>> from Bio import SeqIO
        >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
        >>> print("%s %s" % (record.id, record.seq))
        slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
        >>> print(list(record.letter_annotations))
        ['solexa_quality']
        >>> print(record.letter_annotations["solexa_quality"])
        [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5]

        Now take the reverse complement, here we explicitly give a new
        identifier (the old identifier with a suffix):

        >>> rc_record = record.reverse_complement(id=record.id + "_rc")
        >>> print("%s %s" % (rc_record.id, rc_record.seq))
        slxa_0001_1_0001_01_rc NNNNNNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT

        Notice that the per-letter-annotations have also been reversed,
        although this may not be appropriate for all cases.

        >>> print(rc_record.letter_annotations["solexa_quality"])
        [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]

        Now for the features, we need a different example. Parsing a GenBank
        file is probably the easiest way to get an nice example with features
        in it...

        >>> from Bio import SeqIO
        >>> with open("GenBank/pBAD30.gb") as handle:
        ...     plasmid = SeqIO.read(handle, "gb")
        >>> print("%s %i" % (plasmid.id, len(plasmid)))
        pBAD30 4923
        >>> plasmid.seq
        Seq('GCTAGCGGAGTGTATACTGGCTTACTATGTTGGCACTGATGAGGGTGTCAGTGA...ATG')
        >>> len(plasmid.features)
        13

        Now, let's take the reverse complement of this whole plasmid:

        >>> rc_plasmid = plasmid.reverse_complement(id=plasmid.id+"_rc")
        >>> print("%s %i" % (rc_plasmid.id, len(rc_plasmid)))
        pBAD30_rc 4923
        >>> rc_plasmid.seq
        Seq('CATGGGCAAATATTATACGCAAGGCGACAAGGTGCTGATGCCGCTGGCGATTCA...AGC')
        >>> len(rc_plasmid.features)
        13

        Let's compare the first CDS feature - it has gone from being the
        second feature (index 1) to the second last feature (index -2), its
        strand has changed, and the location switched round.

        >>> print(plasmid.features[1])
        type: CDS
        location: [1081:1960](-)
        qualifiers:
            Key: label, Value: ['araC']
            Key: note, Value: ['araC regulator of the arabinose BAD promoter']
            Key: vntifkey, Value: ['4']
        <BLANKLINE>
        >>> print(rc_plasmid.features[-2])
        type: CDS
        location: [2963:3842](+)
        qualifiers:
            Key: label, Value: ['araC']
            Key: note, Value: ['araC regulator of the arabinose BAD promoter']
            Key: vntifkey, Value: ['4']
        <BLANKLINE>

        You can check this new location, based on the length of the plasmid:

        >>> len(plasmid) - 1081
        3842
        >>> len(plasmid) - 1960
        2963

        Note that if the SeqFeature annotation includes any strand specific
        information (e.g. base changes for a SNP), this information is not
        amended, and would need correction after the reverse complement.

        Note trying to reverse complement a protein SeqRecord raises an
        exception:

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> protein_rec = SeqRecord(Seq("MAIVMGR"), id="Test",
        ...                         annotations={"molecule_type": "protein"})
        >>> protein_rec.reverse_complement()
        Traceback (most recent call last):
           ...
        ValueError: Proteins do not have complements!

        If you have RNA without any U bases, it must be annotated as RNA
        otherwise it will be treated as DNA by default with A mapped to T:

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> rna1 = SeqRecord(Seq("ACG"), id="Test")
        >>> rna2 = SeqRecord(Seq("ACG"), id="Test", annotations={"molecule_type": "RNA"})
        >>> print(rna1.reverse_complement(id="RC", description="unk").format("fasta"))
        >RC unk
        CGT
        <BLANKLINE>
        >>> print(rna2.reverse_complement(id="RC", description="RNA").format("fasta"))
        >RC RNA
        CGU
        <BLANKLINE>

        Also note you can reverse complement a SeqRecord using a MutableSeq:

        >>> from Bio.Seq import MutableSeq
        >>> from Bio.SeqRecord import SeqRecord
        >>> rec = SeqRecord(MutableSeq("ACGT"), id="Test")
        >>> rec.seq[0] = "T"
        >>> print("%s %s" % (rec.id, rec.seq))
        Test TCGT
        >>> rc = rec.reverse_complement(id=True)
        >>> print("%s %s" % (rc.id, rc.seq))
        Test ACGA
        """
        from Bio.Seq import Seq, MutableSeq  # Lazy to avoid circular imports

        if "protein" in self.annotations.get("molecule_type", ""):
            raise ValueError("Proteins do not have complements!")
        if "RNA" in self.annotations.get("molecule_type", ""):
            if isinstance(self.seq, MutableSeq):
                # Does not currently have reverse_complement_rna method:
                answer = SeqRecord(Seq(self.seq).reverse_complement_rna())
            else:
                answer = SeqRecord(self.seq.reverse_complement_rna())
        else:
            # Default to DNA
            if isinstance(self.seq, MutableSeq):
                # Currently the MutableSeq reverse complement is in situ
                answer = SeqRecord(Seq(self.seq).reverse_complement())
            else:
                answer = SeqRecord(self.seq.reverse_complement())
        if isinstance(id, str):
            answer.id = id
        elif id:
            answer.id = self.id
        if isinstance(name, str):
            answer.name = name
        elif name:
            answer.name = self.name
        if isinstance(description, str):
            answer.description = description
        elif description:
            answer.description = self.description
        if isinstance(dbxrefs, list):
            answer.dbxrefs = dbxrefs
        elif dbxrefs:
            # Copy the old dbxrefs
            answer.dbxrefs = self.dbxrefs[:]
        if isinstance(features, list):
            answer.features = features
        elif features:
            # Copy the old features, adjusting location and string
            length = len(answer)
            answer.features = [f._flip(length) for f in self.features]
            # The old list should have been sorted by start location,
            # reversing it will leave it sorted by what is now the end position,
            # so we need to resort in case of overlapping features.
            # NOTE - In the common case of gene before CDS (and similar) with
            # the exact same locations, this will still maintain gene before CDS
            answer.features.sort(key=lambda x: x.location.start.position)
        if isinstance(annotations, dict):
            answer.annotations = annotations
        elif annotations:
            # Copy the old annotations,
            answer.annotations = self.annotations.copy()
        if isinstance(letter_annotations, dict):
            answer.letter_annotations = letter_annotations
        elif letter_annotations:
            # Copy the old per letter annotations, reversing them
            for key, value in self.letter_annotations.items():
                answer._per_letter_annotations[key] = value[::-1]
        return answer

    def translate(
        self,
        # Seq translation arguments:
        table="Standard",
        stop_symbol="*",
        to_stop=False,
        cds=False,
        gap=None,
        # SeqRecord annotation arguments:
        id=False,
        name=False,
        description=False,
        features=False,
        annotations=False,
        letter_annotations=False,
        dbxrefs=False,
    ):
        """Return new SeqRecord with translated sequence.

        This calls the record's .seq.translate() method (which describes
        the translation related arguments, like table for the genetic code),

        By default the new record does NOT preserve the sequence identifier,
        name, description, general annotation or database cross-references -
        these are unlikely to apply to the translated sequence.

        You can specify the returned record's id, name and description as
        strings, or True to keep that of the parent, or False for a default.

        You can specify the returned record's features with a list of
        SeqFeature objects, or False (default) to omit them.

        You can also specify both the returned record's annotations and
        letter_annotations as dictionaries, True to keep that of the parent
        (annotations only), or False (default) to omit them.

        e.g. Loading a FASTA gene and translating it,

        >>> from Bio import SeqIO
        >>> gene_record = SeqIO.read("Fasta/sweetpea.nu", "fasta")
        >>> print(gene_record.format("fasta"))
        >gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds
        CAGGCTGCGCGGTTTCTATTTATGAAGAACAAGGTCCGTATGATAGTTGATTGTCATGCA
        AAACATGTGAAGGTTCTTCAAGACGAAAAACTCCCATTTGATTTGACTCTGTGCGGTTCG
        ACCTTAAGAGCTCCACATAGTTGCCATTTGCAGTACATGGCTAACATGGATTCAATTGCT
        TCATTGGTTATGGCAGTGGTCGTCAATGACAGCGATGAAGATGGAGATAGCCGTGACGCA
        GTTCTACCACAAAAGAAAAAGAGACTTTGGGGTTTGGTAGTTTGTCATAACACTACTCCG
        AGGTTTGTT
        <BLANKLINE>

        And now translating the record, specifying the new ID and description:

        >>> protein_record = gene_record.translate(table=11,
        ...                                        id="phya",
        ...                                        description="translation")
        >>> print(protein_record.format("fasta"))
        >phya translation
        QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
        SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV
        <BLANKLINE>

        """
        if "protein" == self.annotations.get("molecule_type", ""):
            raise ValueError("Proteins cannot be translated!")
        answer = SeqRecord(
            self.seq.translate(
                table=table, stop_symbol=stop_symbol, to_stop=to_stop, cds=cds, gap=gap
            )
        )
        if isinstance(id, str):
            answer.id = id
        elif id:
            answer.id = self.id
        if isinstance(name, str):
            answer.name = name
        elif name:
            answer.name = self.name
        if isinstance(description, str):
            answer.description = description
        elif description:
            answer.description = self.description
        if isinstance(dbxrefs, list):
            answer.dbxrefs = dbxrefs
        elif dbxrefs:
            # Copy the old dbxrefs
            answer.dbxrefs = self.dbxrefs[:]
        if isinstance(features, list):
            answer.features = features
        elif features:
            # Does not make sense to copy old features as locations wrong
            raise TypeError("Unexpected features argument %r" % features)
        if isinstance(annotations, dict):
            answer.annotations = annotations
        elif annotations:
            # Copy the old annotations
            answer.annotations = self.annotations.copy()
        # Set/update to protein:
        answer.annotations["molecule_type"] = "protein"
        if isinstance(letter_annotations, dict):
            answer.letter_annotations = letter_annotations
        elif letter_annotations:
            # Does not make sense to copy these as length now wrong
            raise TypeError(
                "Unexpected letter_annotations argument %r" % letter_annotations
            )
        return answer


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
