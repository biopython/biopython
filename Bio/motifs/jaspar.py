from Bio.motifs import Motif, Instances
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


def read(handle, format):
    alphabet = IUPAC.unambiguous_dna
    counts = {}
    if format=="pfm":
        # reads the motif from Jaspar .pfm file
        letters = "ACGT"
        for letter, line in zip(letters, handle):
            words = line.split()
            #if there is a letter in the beginning, ignore it
            if words[0]==letter:
                words = words[1:]
            counts[letter] = map(float, words)
        motif = Motif(alphabet, counts=counts)
    elif format=="sites":
        # reads the motif from Jaspar .sites file
        instances = []
        for line in handle:
            if not line.startswith(">"):
                break
            # line contains the header ">...."
            # now read the actual sequence
            line = handle.next()
            instance = ""
            for c in line.strip():
                if c==c.upper():
                    instance += c
            instance = Seq(instance, alphabet)
            instances.append(instance)
        instances = Instances(instances, alphabet)
        motif = Motif(alphabet, instances=instances)
    else:
        raise ValueError("Unknown format %s" % format)
    motif.mask = "*"*motif.length
    return motif

def write(motif):
    """Returns the pfm representation of the motif
    """
    letters = "ACGT"
    counts = motif.counts
    lines = []
    for letter in letters:
        terms = map(str, counts[letter])
        line = "\t".join(terms) + "\n"
        lines.append(line)
    # Finished; glue the lines together
    text = "".join(lines)
    return text
