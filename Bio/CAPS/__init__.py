# Copyright Jonathan Taylor 2005
class CAPS:
  """A differential cutsite.

  Members:
  start       Where the CAPS lives in the alignment indexed from 1.
  enzyme      The enzyme that causes this CAPS.
  cuts_in     A list of sequences as indexes it cuts in.
  blocked_in  A list of sequences as indexes it is blocked in.

  """

  def __init__(self, **kwds):
    """Initialize a CAPS.

    Each member should be included as a keyword
    """
    
    self.start = int(kwds["start"])
    self.enzyme = kwds["enzyme"]
    self.cuts_in = kwds["cuts_in"]
    self.blocked_in = kwds["blocked_in"]

class AlignmentHasDifferentLengthsError(Exception):
  pass

class CAPSMap:
  """A map of an alignment showing all possible caps.

  Members:
  alignment  The alignment that is mapped.
  caps       A list of possible CAPS markers.
  """

  def __init__(self, alignment, enzymes = []):
    """Initialize the CAPSMap

    Required:
    alignment    The alignment to be mapped.

    Optional:
    enzymes      The enzymes to be used to create the map.
    """

    self.sequences = [rec.seq for rec in alignment.get_all_seqs()]
    self.size = len(self.sequences)
    self.length = len(self.sequences[0])
    for seq in self.sequences:
      if len(seq) != self.length:
        raise AlignmentHasDifferentLengthsError

    self.alignment = alignment
    self.enzymes = enzymes

    # look for caps
    self._digest()
  
  def _digest_with(self, enzyme):
    cuts = {}
    all = []

    # go through each sequence
    for seq in self.sequences:

      # FIXME: Possibly Bio.Restriction should be changed to handle gaps
      from Bio.Seq import Seq
      ungapped_seq = Seq(seq.tostring().replace("-", "n"))
      
      # grab all the cuts in the sequence
      cuts[seq] = [cut - enzyme.fst5 for cut in enzyme.search(ungapped_seq)]

      # maintain a list of all cuts in all sequences
      all.extend(cuts[seq])

    # we sort the all list and remove duplicates
    all.sort()
    
    last = -999
    new = []
    for cut in all:
      if cut != last:
        new.append(cut)
      last = cut

    all = new
    # all now has indices for all sequences in the alignment

    for cut in all:
      # test the caps

      cuts_in = []
      blocked_in = []

      for i in range(0, self.size):
        seq = self.sequences[i]
        if cut in cuts[seq]:
          cuts_in.append(i)
        else:
          blocked_in.append(i)

      if cuts_in != [] and blocked_in != []:
        self.caps.append(CAPS(start = cut, enzyme = enzyme, cuts_in = cuts_in, blocked_in = blocked_in))

  def _digest(self):
    """Finds all caps in the alignment"""

    # start with no caps
    self.caps = []

    # check each enzyme for caps
    for enzyme in self.enzymes:
      self._digest_with(enzyme)

  def __str__(self):
    """Gives a summary of the CAPS Map."""
    s = "Enzyme - Location - Cuts - Blocks"

    for c in self.caps:
      s += "\n%s - %d - %d - %d" % (c.enzyme, c.start, len(c.cuts_in), len(c.blocked_in))
    return s
