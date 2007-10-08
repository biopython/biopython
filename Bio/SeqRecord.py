# Stores data about the sequence

# NEEDS TO BE SYNCH WITH THE REST OF BIOPYTHON AND BIOPERL

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
    """
    def __init__(self, seq, id = "<unknown id>", name = "<unknown name>",
                 description = "<unknown description>", dbxrefs = None,
                 features = None):
        """Create a SeqRecord

        Arguments:
        seq         - Sequence, required (Seq object)
        id          - Sequence identifier, recommended (string)
        name        - Seqeuence name, optional (string)
        description - Seqeuence description, optional (string)
        dbxrefs     - Database cross references, optional (list of strings)
        features    - Any (sub)features, optional (list of SeqFeature objects)

        Note that while an id is optional, we strongly recommend you supply a
        unique id string for each record.  This is especially important
        if you wish to write your sequences to a file.

        You can create a 'blank' SeqRecord object can then populated the
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
        lines = []
        if self.id : lines.append("ID: %s" % self.id)
        if self.name : lines.append("Name: %s" % self.name)
        if self.description : lines.append("Desription: %s" % self.description)
        if self.dbxrefs : lines.append("Database cross-references: " \
                                       + ", ".join(self.dbxrefs))
        for a in self.annotations:
            lines.append("/%s=%s" % (a, str(self.annotations[a])))
        lines.append(str(self.seq))
        return "\n".join(lines)

    def __repr__(self) :
        return "SeqRecord(seq=%s, id=%s, name=%s, description=%s, dbxrefs=%s)" \
        % tuple(map(repr, (self.seq, self.id, self.name,
                           self.description, self.dbxrefs)))
        
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

    print record

    #One way to create a minimal record.
    record2 = SeqRecord(Seq(""))
