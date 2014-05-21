# Copyright 2009-2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Lazy parsing and feature indexing of sequence files (PRIVATE).

You are not expected to access this module, or any of its code, directly. This
is all handled internally by the following functions when accessed with the
lazy=True kwarg.
Bio.SeqIO.parse(.... lazy=True)
Bio.SeqIO.index(..., lazy=True)
Bio.SeqIO.read(..., lazy=True) 


The lazy loading parsers will read an entire record efficiently storing and 
parsing only the mimimum amount of information to provide fast lookup. Records
returned from parse(), index(), and read() will act as a proxy to a regular 
SeqRecord class.

The lazy loading strategy is, for large sequences, faster to initialize than
SeqIO.parse() and more memory efficient. It is slower to initialize than
SeqIO.index() but also more memory efficient.

The lazy loader will partially parse a sequence file noting the sequence
span and file position of sequence annoations. These annotations will
be stored for efficient retrieval from the file when requested. Sequence 
data will also be efficiently pulled from structured regions
and file

Note that this means all parsing is on demand, so improperly formatted
records may not trigger an exception unless the problem region is requested.
This is by design.
"""

if __name__ == "__main__":
    #this junk allows relative importing for testing in situ
    #will not be included once all unit tests have been exported to
    #the external unit test suite for builds
    import imp
    import os
    biopath = os.path.split(os.getcwd())[0]
    srpath = os.path.join(biopath,"SeqRecord.py")
    SeqRecord = imp.load_source("SeqRecord", srpath)
    _RestrictedDict = SeqRecord._RestrictedDict
    SeqRecord = SeqRecord.SeqRecord
else:
    from ..SeqRecord import SeqRecord, _RestrictedDict

class SeqRecordProxyBase(SeqRecord):
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



    #subclassing to-do list
    attributes that must be filled
     -      self._seq_len    :::   int
     

    methods need implementation by derived class:
     -      self.__init__    ::: the init provided in the base class is for testing



    """
    
    def _set_seq(self, value):
        #TODO - Add a deprecation warning that the seq should be write only?
        if self._per_letter_annotations:
            #TODO - Make this a warning? Silently empty the dictionary?
            raise ValueError("You must empty the letter annotations first!")
        self._seq = value
        try:
            self._per_letter_annotations = _RestrictedDict(length=len(self.seq))
        except AttributeError:
            #e.g. seq is None
            self._per_letter_annotations = _RestrictedDict(length=0)

    seq = property(fget=lambda self: self._seq,
                   fset=_set_seq,
                   doc="The sequence itself, as a Seq or MutableSeq object.")

    def __getitem__(self, index):
        """Returns a sub-sequence or an individual letter.

        This will need to cleverly either return a sequence letter
        or a copy of itself with new marker indices.
        """
        if isinstance(index, int):
            #NOTE - The sequence level annotation like the id, name, etc
            #do not really apply to a single character.  However, should
            #we try and expose any per-letter-annotation here?  If so how?
            return self.seq[index]
        elif isinstance(index, slice):
            if self.seq is None:
                raise ValueError("If the sequence is None, we cannot slice it.")
            parent_length = len(self)
            answer = self.__class__(self.seq[index],
                                    id=self.id,
                                    name=self.name,
                                    description=self.description)
            #TODO - The desription may no longer apply.
            #It would be safer to change it to something
            #generic like "edited" or the default value.

            #Don't copy the annotation dict and dbxefs list,
            #they may not apply to a subsequence.
            #answer.annotations = dict(self.annotations.items())
            #answer.dbxrefs = self.dbxrefs[:]
            #TODO - Review this in light of adding SeqRecord objects?

            #TODO - Cope with strides by generating ambiguous locations?
            start, stop, step = index.indices(parent_length)
            if step == 1:
                #Select relevant features, add them with shifted locations
                #assert str(self.seq)[index] == str(self.seq)[start:stop]
                for f in self.features:
                    if f.ref or f.ref_db:
                        #TODO - Implement this (with lots of tests)?
                        import warnings
                        warnings.warn("When slicing SeqRecord objects, any "
                              "SeqFeature referencing other sequences (e.g. "
                              "from segmented GenBank records) are ignored.")
                        continue
                    if start <= f.location.nofuzzy_start \
                    and f.location.nofuzzy_end <= stop:
                        answer.features.append(f._shift(-start))

            #Slice all the values to match the sliced sequence
            #(this should also work with strides, even negative strides):
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
        """Implements the 'in' keyword, searches the sequence.

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
        >>> from Bio.Alphabet import generic_dna
        >>> Seq("AAA") in record
        True
        >>> Seq("AAA", generic_dna) in record
        True

        See also the Seq object's __contains__ method.
        """
        return char in self.seq

    def __str__(self):
        """A human readable summary of the record and its annotation (string).

        The python built in function str works by calling the object's ___str__
        method.  e.g.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Alphabet import IUPAC
        >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
        ...                         IUPAC.protein),
        ...                    id="YP_025292.1", name="HokC",
        ...                    description="toxic membrane protein, small")
        >>> print(str(record))
        ID: YP_025292.1
        Name: HokC
        Description: toxic membrane protein, small
        Number of features: 0
        Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF', IUPACProtein())

        In this example you don't actually need to call str explicity, as the
        print command does this automatically:

        >>> print(record)
        ID: YP_025292.1
        Name: HokC
        Description: toxic membrane protein, small
        Number of features: 0
        Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF', IUPACProtein())

        Note that long sequences are shown truncated.
        """
        lines = []
        if self.id:
            lines.append("ID: %s" % self.id)
        if self.name:
            lines.append("Name: %s" % self.name)
        if self.description:
            lines.append("Description: %s" % self.description)
        if self.dbxrefs:
            lines.append("Database cross-references: "
                         + ", ".join(self.dbxrefs))
        lines.append("Number of features: %i" % len(self.features))
        for a in self.annotations:
            lines.append("/%s=%s" % (a, str(self.annotations[a])))
        if self.letter_annotations:
            lines.append("Per letter annotation for: "
                         + ", ".join(self.letter_annotations))
        #Don't want to include the entire sequence,
        #and showing the alphabet is useful:
        lines.append(repr(self.seq))
        return "\n".join(lines)

    def __repr__(self):
        """A concise summary of the record for debugging (string).

        The python built in function repr works by calling the object's ___repr__
        method.  e.g.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Alphabet import generic_protein
        >>> rec = SeqRecord(Seq("MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKAT"
        ...                    +"GEMKEQTEWHRVVLFGKLAEVASEYLRKGSQVYIEGQLRTRKWTDQ"
        ...                    +"SGQDRYTTEVVVNVGGTMQMLGGRQGGGAPAGGNIGGGQPQGGWGQ"
        ...                    +"PQQPQGGNQFSGGAQSRPQQSAPAAPSNEPPMDFDDDIPF",
        ...                    generic_protein),
        ...                 id="NP_418483.1", name="b4059",
        ...                 description="ssDNA-binding protein",
        ...                 dbxrefs=["ASAP:13298", "GI:16131885", "GeneID:948570"])
        >>> print(repr(rec))
        SeqRecord(seq=Seq('MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKATGEMKEQTE...IPF', ProteinAlphabet()), id='NP_418483.1', name='b4059', description='ssDNA-binding protein', dbxrefs=['ASAP:13298', 'GI:16131885', 'GeneID:948570'])

        At the python prompt you can also use this shorthand:

        >>> rec
        SeqRecord(seq=Seq('MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKATGEMKEQTE...IPF', ProteinAlphabet()), id='NP_418483.1', name='b4059', description='ssDNA-binding protein', dbxrefs=['ASAP:13298', 'GI:16131885', 'GeneID:948570'])

        Note that long sequences are shown truncated. Also note that any
        annotations, letter_annotations and features are not shown (as they
        would lead to a very long string).
        """
        return self.__class__.__name__ \
         + "(seq=%s, id=%s, name=%s, description=%s, dbxrefs=%s)" \
         % tuple(map(repr, (self.seq, self.id, self.name,
                            self.description, self.dbxrefs)))

    def __len__(self):
        """Returns the length of the sequence.
        
        The derived class must set attribute _seq_len
        """
        return len(self._seq_len)


    def upper(self):
        """Returns a copy of the record with an upper case sequence.

           proxy class needs a new implementation
        """
        raise NotImplementedError(" this needs to be implemented")


    def lower(self):
        """Returns a copy of the record with a lower case sequence.

        possibly do something clever or get implemnted by the derived class
        """
        raise NotImplementedError(" this needs to be implemented")

    # All methods tagged below are implemented in the base class
    #
    #def __bool__(self):
    #__nonzero__= __bool__
    #def reverse_complement(self, id=False, name=False, description=False,
    #def __radd__(self, other):
    #def format(self, format):
    #def __format__(self, format_spec):
    #def __add__(self, other):
    #    returns non-proxy class where possible
    #def __radd__(self, other):
    #    returns non-proxy class where possible

class TestSeqRecordBaseClass(SeqRecordProxyBase):
    """ this class implements simple forms of required methods and attributes
    """
    def __init__(self, seq, id = "<unknown id>", name = "<unknown name>",
                 description = "<unknown description>", dbxrefs = None,
                 features = None, annotations = None,
                 letter_annotations = None):
        """Create a SeqRecord.

        Arguments:
         - seq         - Sequence, required (Seq, MutableSeq or UnknownSeq)
         - id          - Sequence identifier, recommended (string)
         - name        - Sequence name, optional (string)
         - description - Sequence description, optional (string)
         - dbxrefs     - Database cross references, optional (list of strings)
         - features    - Any (sub)features, optional (list of SeqFeature objects)
         - annotations - Dictionary of annotations for the whole sequence
         - letter_annotations - Dictionary of per-letter-annotations, values
                                should be strings, list or tuples of the same
                                length as the full sequence.

        """
        if id is not None and not isinstance(id, basestring):
            #Lots of existing code uses id=None... this may be a bad idea.
            raise TypeError("id argument should be a string")
        if not isinstance(name, basestring):
            raise TypeError("name argument should be a string")
        if not isinstance(description, basestring):
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
                #Should we allow this and use a normal unrestricted dict?
                self._per_letter_annotations = _RestrictedDict(length=0)
            else:
                try:
                    self._per_letter_annotations = \
                                              _RestrictedDict(length=len(seq))
                except:
                    raise TypeError("seq argument should be a Seq object or similar")
        else:
            #This will be handled via the property set function, which will
            #turn this into a _RestrictedDict and thus ensure all the values
            #in the dict are the right length
            self.letter_annotations = letter_annotations

        # annotations about parts of the sequence
        if features is None:
            features = []
        elif not isinstance(features, list):
            raise TypeError("features argument should be a list (of SeqFeature objects)")
        self.features = features

    #TODO - Just make this a read only property?
    def _set_per_letter_annotations(self, value):
        if not isinstance(value, dict):
            raise TypeError("The per-letter-annotations should be a "
                            "(restricted) dictionary.")
        #Turn this into a restricted-dictionary (and check the entries)
        try:
            self._per_letter_annotations = _RestrictedDict(length=len(self.seq))
        except AttributeError:
            #e.g. seq is None
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
        """)
    
if __name__ == "__main__":
    import unittest
    #from Bio import SeqRecord
    #from Bio.SeqIO import _lazy
    
    class SeqRecordProxyBaseClassTests(unittest.TestCase):

        def setUp(self):
            pass

        def test_nothing(self):
            """An addition test"""
            a = SeqRecordProxyBase("sequencefake", "fakeid")
            self.assertEqual(5, 5)
            
    unittest.main( exit=False )
    """if __name__ == "__main__":
        runner = unittest.TextTestRunner(verbosity = 2)
        unittest.main(testRunner=runner)"""