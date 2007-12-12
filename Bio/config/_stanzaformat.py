"""This module reads and writes (actually, write not implemented yet)
files in the OBF stanza format.

"""
class StanzaFormat:
    """Contains information from a stanza-formatted file."""
    def __init__(self, version, stanzas):
        self.version = version
        self.stanzas = stanzas

class Stanza:
    """Contains information about one stanza in a file."""
    def __init__(self, name, tag_value_pairs):
        self.name = name
        self.tag_value_pairs = tag_value_pairs
        
        dict = {}
        for tag, value in tag_value_pairs:
            dict[tag] = value
        self.tag_value_dict = dict

def load(handle):
    """load(handle) -> StanzaFormat object"""
    import ConfigParser
    parser = ConfigParser.ConfigParser()

    # Read the VERSION string.
    line = handle.readline()
    while line and not line.startswith("VERSION"):
        line = handle.readline()
    assert line, "I could not find the VERSION line"
    x, version = line.split("=")
    version = version.strip()

    try:
        parser.readfp(handle)
    except ConfigParser.Error, x:
        raise ValueError, x
    stanzas = []
    for section in parser.sections():
        pairs = []
        for tag in parser.options(section):
            value = parser.get(section, tag)
            pairs.append((tag, value))
        stanzas.append(Stanza(section, pairs))
    return StanzaFormat(version, stanzas)
