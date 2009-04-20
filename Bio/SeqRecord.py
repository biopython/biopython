# Copyright 2000-2002 Andrew Dalke.
# Copyright 2002-2004 Brad Chapman.
# Copyright 2006-2009 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Represent a Sequence Record, a sequence with annotation."""
__docformat__ = "epytext en" #Simple markup to show doctests nicely

# NEEDS TO BE SYNCH WITH THE REST OF BIOPYTHON AND BIOPERL
# In particular, the SeqRecord and BioSQL.BioSeq.DBSeqRecord classes
# need to be in sync (this is the BioSQL "Database SeqRecord", see
# also BioSQL.BioSeq.DBSeq which is the "Database Seq" class)

class _RestrictedDict(dict):
    """Dict which only allows sequences of given length as values (PRIVATE).

    This simple subclass of the python dictionary is used in the SeqRecord
    object for holding per-letter-annotations.  This class is intended to
    prevent simple errors by only allowing python sequences (e.g. lists,
    strings and tuples) to be stored, and only if their length matches that
    expected (the length of the SeqRecord's seq object).  It cannot however
    prevent the entries being edited in situ (for example appending entries
    to a list).
    """
    def __init__(self, length) :
        """Create an EMPTY restricted dictionary."""
        dict.__init__(self)
        self._length = int(length)
    def __setitem__(self, key, value) :
        if not hasattr(value,"__len__") or not hasattr(value,"__getitem__") \
        or len(value) != self._length :
            raise TypeError("We only allow python sequences (lists, tuples or "
                            "strings) of length %i." % self._length)
        dict.__setitem__(self, key, value)

class SeqRecord(object):
    """A SeqRecord object holds a sequence and information about it.

    Main attributes:
     - id          - Identifier such as a locus tag (string)
     - seq         - The sequence itself (Seq object)

    Additional attributes:
     - name        - Sequence name, e.g. gene name (string)
     - description - Additional text (string)
     - dbxrefs     - List of database cross references (list of strings)
     - features    - Any (sub)features defined (list of SeqFeature objects)
     - annotations - Further information about the whole sequence (dictionary)
                     Most entries are lists of strings.
     - letter_annotations - Per letter/symbol annotation (restricted
                     dictionary). This holds python sequences (lists, strings
                     or tuples) whose length matches that of the sequence.
                     A typical use would be to hold a list of integers
                     representing sequencing quality scores, or a string
                     representing the secondary structure.

    You will typically use Bio.SeqIO to read in sequences from files as
    SeqRecord objects.  However, you may want to create your own SeqRecord
    objects directly (see the __init__ method for further details):

    >>> from Bio.Seq import Seq
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.Alphabet import IUPAC
    >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
    ...                         IUPAC.protein),
    ...                    id="YP_025292.1", name="HokC",
    ...                    description="toxic membrane protein")
    >>> print record
    ID: YP_025292.1
    Name: HokC
    Description: toxic membrane protein
    Number of features: 0
    Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF', IUPACProtein())

    If you want to save SeqRecord objects to a sequence file, use Bio.SeqIO
    for this.  For the special case where you want the SeqRecord turned into
    a string in a particular file format there is a format method which uses
    Bio.SeqIO internally:

    >>> print record.format("fasta")
    >YP_025292.1 toxic membrane protein
    MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
    <BLANKLINE>
    """
    def __init__(self, seq, id = "<unknown id>", name = "<unknown name>",
                 description = "<unknown description>", dbxrefs = None,
                 features = None):
        """Create a SeqRecord.

        Arguments:
         - seq         - Sequence, required (Seq or Mutable object)
         - id          - Sequence identifier, recommended (string)
         - name        - Sequence name, optional (string)
         - description - Sequence description, optional (string)
         - dbxrefs     - Database cross references, optional (list of strings)
         - features    - Any (sub)features, optional (list of SeqFeature objects)

        You will typically use Bio.SeqIO to read in sequences from files as
        SeqRecord objects.  However, you may want to create your own SeqRecord
        objects directly.

        Note that while an id is optional, we strongly recommend you supply a
        unique id string for each record.  This is especially important
        if you wish to write your sequences to a file.

        If you don't have the actual sequence, but you do know its length,
        then using the UnknownSeq object from Bio.Seq is appropriate.

        You can create a 'blank' SeqRecord object, and then populate the
        attributes later.  Note that currently the annotations and the
        letter_annotations dictionaries cannot be specified when creating
        the SeqRecord.
        """
        if id is not None and not isinstance(id, basestring) :
            #Lots of existing code uses id=None... this may be a bad idea.
            raise TypeError("id argument should be a string")
        if not isinstance(name, basestring) :
            raise TypeError("name argument should be a string")
        if not isinstance(description, basestring) :
            raise TypeError("description argument should be a string")
        if dbxrefs is not None and not isinstance(dbxrefs, list) :
            raise TypeError("dbxrefs argument should be a list (of strings)")
        if features is not None and not isinstance(features, list) :
            raise TypeError("features argument should be a list (of SeqFeature objects)")
        self._seq = seq
        self.id = id
        self.name = name
        self.description = description
        if dbxrefs is None:
            dbxrefs = []
        self.dbxrefs = dbxrefs
        # annotations about the whole sequence
        self.annotations = {}

        # annotations about each letter in the sequence
        if seq is None :
            #Should we allow this and use a normal unrestricted dict?
            self._per_letter_annotations = _RestrictedDict(length=0)
        else :
            try :
                self._per_letter_annotations = _RestrictedDict(length=len(seq))
            except :
                raise TypeError("seq argument should be a Seq or MutableSeq")
        
        # annotations about parts of the sequence
        if features is None:
            features = []
        self.features = features

    #TODO - Just make this a read only property?
    def _set_per_letter_annotations(self, value) :
        if not isinstance(value, dict) :
            raise TypeError("The per-letter-annotations should be a "
                            "(restricted) dictionary.")
        #Turn this into a restricted-dictionary (and check the entries)
        try :
            self._per_letter_annotations = _RestrictedDict(length=len(self.seq))
        except AttributeError :
            #e.g. seq is None
            self._per_letter_annotations = _RestrictedDict(length=0)
        self._per_letter_annotations.update(value)
    letter_annotations = property( \
        fget=lambda self : self._per_letter_annotations,
        fset=_set_per_letter_annotations,
        doc="""Dictionary of per-letter-annotation for the sequence.

        For example, this can hold quality scores used in FASTQ or QUAL files.
        Consider this example using Bio.SeqIO to read in an example Solexa
        variant FASTQ file as a SeqRecord:

        >>> from Bio import SeqIO
        >>> handle = open("Quality/solexa.fastq", "rU")
        >>> record = SeqIO.read(handle, "fastq-solexa")
        >>> handle.close()
        >>> print record.id, record.seq
        slxa_0013_1_0001_24 ACAAAAATCACAAGCATTCTTATACACC
        >>> print record.letter_annotations.keys()
        ['solexa_quality']
        >>> print record.letter_annotations["solexa_quality"]
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -6, -1, -1, -4, -1, -4, -19, -10, -27, -18]

        The per-letter-annotaions get sliced automatically if you slice the
        parent SeqRecord, for example taking the last ten bases:

        >>> sub_record = record[-10:]
        >>> print sub_record.id, sub_record.seq
        slxa_0013_1_0001_24 CTTATACACC
        >>> print sub_record.letter_annotations["solexa_quality"]
        [-6, -1, -1, -4, -1, -4, -19, -10, -27, -18]

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

    def _set_seq(self, value) :
        #TODO - Add a deprecation warning that the seq should be write only?
        if self._per_letter_annotations :
            #TODO - Make this a warning? Silently empty the dictionary?
            raise ValueError("You must empty the letter annotations first!")
        self._seq = value
        try :
            self._per_letter_annotations = _RestrictedDict(length=len(self.seq))
        except AttributeError :
            #e.g. seq is None
            self._per_letter_annotations = _RestrictedDict(length=0)

    seq = property(fget=lambda self : self._seq,
                   fset=_set_seq,
                   doc="The sequence itself, as a Seq or MutableSeq object.")

    def __getitem__(self, index) :
        """Returns a sub-sequence or an individual letter.

        Splicing, e.g. my_record[5:10], returns a new SeqRecord for
        that sub-sequence with approriate annotation preserved.  The
        name, id and description are kept.

        Any per-letter-annotations are sliced to match the requested
        sub-sequence.  Unless a stride is used, all those features
        which fall fully within the subsequence are included (with
        their locations adjusted accordingly).

        However, the annotations dictionary and the dbxrefs list are
        not used for the new SeqRecord, as in general they may not
        apply to the subsequence.  If you want to preserve them, you
        must explictly copy them to the new SeqRecord yourself.

        Using an integer index, e.g. my_record[5] is shorthand for
        extracting that letter from the sequence, my_record.seq[5].

        For example, consider this short protein and its secondary
        structure as encoded by the PDB (e.g. H for alpha helices),
        plus a simple feature for its histidine self phosphorylation
        site:

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> from Bio.Alphabet import IUPAC
        >>> rec = SeqRecord(Seq("MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLAT"
        ...                     "EMMSEQDGYLAESINKDIEECNAIIEQFIDYLR",
        ...                     IUPAC.protein),
        ...                 id="1JOY", name="EnvZ",
        ...                 description="Homodimeric domain of EnvZ from E. coli")
        >>> rec.letter_annotations["secondary_structure"] = \
            "  S  SSSSSSHHHHHTTTHHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHHHHHHHHHHHHHHTT  "
        >>> rec.features.append(SeqFeature(FeatureLocation(20,21),
        ...                     type = "Site"))

        Now let's have a quick look at the full record,

        >>> print rec
        ID: 1JOY
        Name: EnvZ
        Description: Homodimeric domain of EnvZ from E. coli
        Number of features: 1
        Per letter annotation for: secondary_structure
        Seq('MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEE...YLR', IUPACProtein())
        >>> print rec.letter_annotations["secondary_structure"]
          S  SSSSSSHHHHHTTTHHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHHHHHHHHHHHHHHTT  
        >>> print rec.features[0].location
        [20:21]

        Now let's take a sub sequence, here chosen as the first (fractured)
        alpha helix which includes the histidine phosphorylation site:

        >>> sub = rec[11:41]
        >>> print sub
        ID: 1JOY
        Name: EnvZ
        Description: Homodimeric domain of EnvZ from E. coli
        Number of features: 1
        Per letter annotation for: secondary_structure
        Seq('RTLLMAGVSHDLRTPLTRIRLATEMMSEQD', IUPACProtein())
        >>> print sub.letter_annotations["secondary_structure"]
        HHHHHTTTHHHHHHHHHHHHHHHHHHHHHH
        >>> print sub.features[0].location
        [9:10]

        You can also of course omit the start or end values, for
        example to get the first ten letters only:

        >>> print rec[:10]
        ID: 1JOY
        Name: EnvZ
        Description: Homodimeric domain of EnvZ from E. coli
        Number of features: 0
        Per letter annotation for: secondary_structure
        Seq('MAAGVKQLAD', IUPACProtein())

        Or for the last ten letters:

        >>> print rec[-10:]
        ID: 1JOY
        Name: EnvZ
        Description: Homodimeric domain of EnvZ from E. coli
        Number of features: 0
        Per letter annotation for: secondary_structure
        Seq('IIEQFIDYLR', IUPACProtein())

        If you omit both, then you get a copy of the original record (although
        lacking the annotations and dbxrefs):

        >>> print rec[:]
        ID: 1JOY
        Name: EnvZ
        Description: Homodimeric domain of EnvZ from E. coli
        Number of features: 1
        Per letter annotation for: secondary_structure
        Seq('MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEE...YLR', IUPACProtein())

        Finally, indexing with a simple integer is shorthand for pulling out
        that letter from the sequence directly:

        >>> rec[5]
        'K'
        >>> rec.seq[5]
        'K'
        """
        if isinstance(index, int) :
            #NOTE - The sequence level annotation like the id, name, etc
            #do not really apply to a single character.  However, should
            #we try and expose any per-letter-annotation here?  If so how?
            return self.seq[index]
        elif isinstance(index, slice) :
            if self.seq is None :
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
            #answer.annotations = dict(self.annotations.iteritems())
            #answer.dbxrefs = self.dbxrefs[:]
            
            #TODO - Cope with strides by generating ambiguous locations?
            if index.step is None or index.step == 1 :
                #Select relevant features, add them with shifted locations
                if index.start is None :
                    start = 0
                else :
                    start = index.start
                if index.stop is None :
                    stop = -1
                else :
                    stop = index.stop
                if (start < 0 or stop < 0) and parent_length == 0 :
                    raise ValueError, \
                          "Cannot support negative indices without the sequence length"
                if start < 0 :
                    start = parent_length - start
                if stop < 0  :
                    stop  = parent_length - stop + 1
                #assert str(self.seq)[index] == str(self.seq)[start:stop]
                for f in self.features :
                    if start <= f.location.start.position \
                    and f.location.end.position < stop :
                        answer.features.append(f._shift(-start))

            #Slice all the values to match the sliced sequence
            #(this should also work with strides, even negative strides):
            for key, value in self.letter_annotations.iteritems() :
                answer._per_letter_annotations[key] = value[index]

            return answer
        raise ValueError, "Invalid index"

    def __iter__(self) :
        """Iterate over the letters in the sequence.

        For example, using Bio.SeqIO to read in a protein FASTA file:

        >>> from Bio import SeqIO
        >>> record = SeqIO.read(open("Amino/loveliesbleeding.pro"),"fasta")
        >>> for amino in record :
        ...     print amino
        ...     if amino == "L" : break
        X
        A
        G
        L
        >>> print record.seq[3]
        L

        This is just a shortcut for iterating over the sequence directly:

        >>> for amino in record.seq :
        ...     print amino
        ...     if amino == "L" : break
        X
        A
        G
        L
        >>> print record.seq[3]
        L
        
        Note that this does not facilitate iteration together with any
        per-letter-annotation.  However, you can achieve that using the
        python zip function on the record (or its sequence) and the relevant
        per-letter-annotation:
        
        >>> from Bio import SeqIO
        >>> rec = SeqIO.read(open("Quality/solexa.fastq", "rU"),
        ...                  "fastq-solexa")
        >>> print rec.id, rec.seq
        slxa_0013_1_0001_24 ACAAAAATCACAAGCATTCTTATACACC
        >>> print rec.letter_annotations.keys()
        ['solexa_quality']
        >>> for nuc, qual in zip(rec,rec.letter_annotations["solexa_quality"]) :
        ...     if qual < -10 :
        ...         print nuc, qual
        C -19
        C -27
        C -18

        You may agree that using zip(rec.seq, ...) is more explicit than using
        zip(rec, ...) as shown above.
        """
        return iter(self.seq)

    def __str__(self) :
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
        >>> print str(record)
        ID: YP_025292.1
        Name: HokC
        Description: toxic membrane protein, small
        Number of features: 0
        Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF', IUPACProtein())

        In this example you don't actually need to call str explicity, as the
        print command does this automatically:

        >>> print record
        ID: YP_025292.1
        Name: HokC
        Description: toxic membrane protein, small
        Number of features: 0
        Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF', IUPACProtein())

        Note that long sequences are shown truncated.
        """
        lines = []
        if self.id : lines.append("ID: %s" % self.id)
        if self.name : lines.append("Name: %s" % self.name)
        if self.description : lines.append("Description: %s" % self.description)
        if self.dbxrefs : lines.append("Database cross-references: " \
                                       + ", ".join(self.dbxrefs))
        lines.append("Number of features: %i" % len(self.features))
        for a in self.annotations:
            lines.append("/%s=%s" % (a, str(self.annotations[a])))
        if self.letter_annotations :
            lines.append("Per letter annotation for: " \
                         + ", ".join(self.letter_annotations.keys()))
        #Don't want to include the entire sequence,
        #and showing the alphabet is useful:
        lines.append(repr(self.seq))
        return "\n".join(lines)

    def __repr__(self) :
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
        >>> print repr(rec)
        SeqRecord(seq=Seq('MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKATGEMKEQTE...IPF', ProteinAlphabet()), id='NP_418483.1', name='b4059', description='ssDNA-binding protein', dbxrefs=['ASAP:13298', 'GI:16131885', 'GeneID:948570'])

        At the python prompt you can also use this shorthand:

        >>> rec
        SeqRecord(seq=Seq('MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKATGEMKEQTE...IPF', ProteinAlphabet()), id='NP_418483.1', name='b4059', description='ssDNA-binding protein', dbxrefs=['ASAP:13298', 'GI:16131885', 'GeneID:948570'])

        Note that long sequences are shown truncated.
        """
        return self.__class__.__name__ \
         + "(seq=%s, id=%s, name=%s, description=%s, dbxrefs=%s)" \
         % tuple(map(repr, (self.seq, self.id, self.name,
                            self.description, self.dbxrefs)))

    def format(self, format) :
        r"""Returns the record as a string in the specified file format.

        The format should be a lower case string supported as an output
        format by Bio.SeqIO, which is used to turn the SeqRecord into a
        string.  e.g.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Alphabet import IUPAC
        >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
        ...                         IUPAC.protein),
        ...                    id="YP_025292.1", name="HokC",
        ...                    description="toxic membrane protein")
        >>> record.format("fasta")
        '>YP_025292.1 toxic membrane protein\nMKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF\n'
        >>> print record.format("fasta")
        >YP_025292.1 toxic membrane protein
        MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
        <BLANKLINE>

        The python print command automatically appends a new line, meaning
        in this example a blank line is shown.  If you look at the string
        representation you can see there is a trailing new line (shown as
        slash n) which is important when writing to a file or if
        concatenating mutliple sequence strings together.

        Note that this method will NOT work on every possible file format
        supported by Bio.SeqIO (e.g. some are for multiple sequences only).
        """
        #See also the __format__ added for Python 2.6 / 3.0, PEP 3101
        #See also the Bio.Align.Generic.Alignment class and its format()
        return self.__format__(format)

    def __format__(self, format_spec) :
        """Returns the record as a string in the specified file format.

        This method supports the python format() function added in
        Python 2.6/3.0.  The format_spec should be a lower case
        string supported by Bio.SeqIO as an output file format.
        See also the SeqRecord's format() method.
        """
        if format_spec:
            from StringIO import StringIO
            from Bio import SeqIO
            handle = StringIO()
            SeqIO.write([self], handle, format_spec)
            return handle.getvalue()
        else :
            #Follow python convention and default to using __str__
            return str(self)    

    def __len__(self) :
        """Returns the length of the sequence.

        For example, using Bio.SeqIO to read in a FASTA nucleotide file:

        >>> from Bio import SeqIO
        >>> record = SeqIO.read(open("Nucleic/sweetpea.nu"),"fasta")
        >>> len(record)
        309
        >>> len(record.seq)
        309
        """
        return len(self.seq)

    def __nonzero__(self) :
        """Returns True regardless of the length of the sequence.

        This behaviour is for backwards compatibility, since until the
        __len__ method was added, a SeqRecord always evaluated as True.

        Note that in comparison, a Seq object will evaluate to False if it
        has a zero length sequence.

        WARNING: The SeqRecord may in future evaluate to False when its
        sequence is of zero length (in order to better match the Seq
        object behaviour)!
        """
        return True

def _test():
    """Run the Bio.SeqRecord module's doctests (PRIVATE).

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..","Tests")) :
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..","Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"

if __name__ == "__main__":
    _test()
