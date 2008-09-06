# Stores data about the sequence
"""Represent a Sequence Record, a sequence with annotation."""

# NEEDS TO BE SYNCH WITH THE REST OF BIOPYTHON AND BIOPERL
# In particular, the SeqRecord and BioSQL.BioSeq.DBSeqRecord classes
# need to be in sync (this is the BioSQL "Database SeqRecord", see
# also BioSQL.BioSeq.DBSeq which is the "Database Seq" class)

class SeqRecord:
    """A SeqRecord object holds a sequence and information about it.

    Main attributes:
    id          - Identifier such as a locus tag (string)
    seq         - The sequence itself (Seq object)

    Additional attributes:
    name        - Sequence name, e.g. gene name (string)
    description - Additional text (string)
    dbxrefs     - List of database cross references (list of strings)
    features    - Any (sub)features defined (list of SeqFeature objects)
    annotations - Further information about the whole sequence (dictionary)
                  Most entries are lists of strings.
    """
    def __init__(self, seq, id = "<unknown id>", name = "<unknown name>",
                 description = "<unknown description>", dbxrefs = None,
                 features = None):
        """Create a SeqRecord.

        Arguments:
        seq         - Sequence, required (Seq object)
        id          - Sequence identifier, recommended (string)
        name        - Sequence name, optional (string)
        description - Sequence description, optional (string)
        dbxrefs     - Database cross references, optional (list of strings)
        features    - Any (sub)features, optional (list of SeqFeature objects)

        Note that while an id is optional, we strongly recommend you supply a
        unique id string for each record.  This is especially important
        if you wish to write your sequences to a file.

        You can create a 'blank' SeqRecord object, and then populated the
        attributes later.  Note that currently the annotations dictionary
        cannot be specified when creating the SeqRecord."""
        self.seq = seq
        self.id = id
        self.name = name
        self.description = description
        if dbxrefs is None:
            dbxrefs = []
        self.dbxrefs = dbxrefs
        # annotations about the whole sequence
        self.annotations = {}
        
        # annotations about parts of the sequence
        if features is None:
            features = []
        self.features = features

    def __str__(self) :
        """A human readable summary of the record and its annotation."""
        lines = []
        if self.id : lines.append("ID: %s" % self.id)
        if self.name : lines.append("Name: %s" % self.name)
        if self.description : lines.append("Description: %s" % self.description)
        if self.dbxrefs : lines.append("Database cross-references: " \
                                       + ", ".join(self.dbxrefs))
        lines.append("Number of features: %i" % len(self.features))
        for a in self.annotations:
            lines.append("/%s=%s" % (a, str(self.annotations[a])))
        #Don't want to include the entire sequence,
        #and showing the alphabet is useful:
        lines.append(repr(self.seq))
        return "\n".join(lines)

    def __repr__(self) :
        """A concise summary of the record for debugging."""
        return self.__class__.__name__ \
         + "(seq=%s, id=%s, name=%s, description=%s, dbxrefs=%s)" \
         % tuple(map(repr, (self.seq, self.id, self.name,
                            self.description, self.dbxrefs)))

    def format(self, format) :
        """Returns the record as a string in the specified file format.

        The format should be a lower case string supported as an output
        format by Bio.SeqIO, which is used to turn the SeqRecord into a
        string.

        e.g.
        print my_record.to_format("fasta")

        This will NOT work on every possible file format supported by
        Bio.SeqIO (e.g. some are for multiple sequences only).
        """
        #See also the __format__ added for Python 2.6 / 3.0, PEP 3101
        #See also the Bio.Align.Generic.Alignment class and its format()
        return self.__format__(format)

    def __format__(self, format_spec) :
        """Returns the record as a string in the specified file format.

        This method supports the python format() function added in
        Python 2.6/3.0.  The format_spec should be a lower case
        string supported by Bio.SeqIO as an output file format.
        See also the SeqRecord's format() method."""
        if format_spec:
            from StringIO import StringIO
            from Bio import SeqIO
            handle = StringIO()
            SeqIO.write([self], handle, format_spec)
            handle.seek(0)
            return handle.read()
        else :
            #Follow python convention and default to using __str__
            return str(self)    

    def __len__(self) :
        """Returns the length of the sequence."""
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

if __name__ == "__main__" :
    #The following is a very quick example of how to create a SeqRecord object
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_protein
    record = SeqRecord(Seq("MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKAT" \
                          +"GEMKEQTEWHRVVLFGKLAEVASEYLRKGSQVYIEGQLRTRKWTDQ" \
                          +"SGQDRYTTEVVVNVGGTMQMLGGRQGGGAPAGGNIGGGQPQGGWGQ" \
                          +"PQQPQGGNQFSGGAQSRPQQSAPAAPSNEPPMDFDDDIPF",
                           generic_protein),
                       id="NP_418483.1", name="b4059",
                       description="ssDNA-binding protein",
                       dbxrefs=["ASAP:13298", "GI:16131885", "GeneID:948570"])

    #Note that annotations must be added AFTER creating the record
    record.annotations["note"] = "This annotation was added later"

    print str(record)
    print
    print repr(record)
    assert 178 == len(record)
    print
    print "Using .format('fasta'),", repr(record.format("fasta"))
    print
    print "Using .format('tab'),", repr(record.format("tab"))

    #One way to create a minimal record.
    record2 = SeqRecord(Seq(""))
    assert record2 #True eeven though length is zero
    assert not len(record2)
