import xml.sax, re

from Bio import Std


# To help parse XPath queries
_name = "[a-zA-Z_:][-a-zA-Z0-9._:]*"
_pat_tag_re = re.compile(r"""^//(%s)(\[@(%s)=("[^"]*"|'[^']*')\])?$""" %
                         (_name, _name) )
                                                   #')  # emacs cruft


def parse_simple_xpath(s):
    # Only supports two formats
    # //tag
    # //tag[@attr="value"]
    m = _pat_tag_re.match(s)
    if m is None:
        raise TypeError("Cannot yet understand the XPath expression: %r" %
                        (s,))
    tag =  m.group(1)
    if m.group(3) is not None:
        varname = m.group(3)
        varvalue = m.group(4)[1:-1]
        node_matcher = (tag, [(varname, varvalue)])
    else:
        node_matcher = (tag, None)
    return node_matcher



def xpath_index(dbname,
                filenames,
                primary_namespace,
                extract_info,  # pair of (data_value, xpath)
                format = "sequence",
                record_tag = Std.record.tag,
                creator_factory = None,
                ):
    if creator_factory is None:
        import BerkeleyDB
        creator_factory = BerkeleyDB.create
    
    data_names = [x[0] for x in extract_info]
    if primary_namespace not in data_names:
        raise TypeError(
            "No way to get the %r field needed for the primary (unique) id" %
            (primary_namespace,))
    data_names.remove(primary_namespace)

    for prop, xpath in extract_info:
        if prop == primary_namespace:
            break
    else:
        raise TypeError("Property %r has no xpath definition" %
                        (primary_namespace,))

    creator = creator_factory(dbname, primary_namespace, data_names)
    builder = GrabXPathNodes(extract_info)
    for filename in filenames:
        creator.load(filename, builder = builder, record_tag = record_tag,
                     formatname = format)
    creator.close()


class GrabXPathNodes(xml.sax.ContentHandler):
    def __init__(self, extractinfo):
        self._fast_tags = _fast_tags = {}
        for property, xpath in extractinfo:
            tag, attrs = parse_simple_xpath(xpath)
            _fast_tags.setdefault(tag, []).append( (attrs, property) )

        # for doing the endElement in the correct order,
        # which is opposite to the input order
        self._rev_tags = _rev_tags = {}
        for k, v in self._fast_tags.items():
            v = v[:]
            v.reverse()
            self._rev_tags[k] = v

    def uses_tags(self):
        return self._fast_tags.keys()

    def startDocument(self):
        self._text = ""
        self._capture = []
        self.document = {}
        
    def startElement(self, tag, attrs):
        if not self._fast_tags.has_key(tag):
            return
        for want_attrs, prop in self._fast_tags[tag]:
            needed = []
            if want_attrs is None:
                needed.append(prop)
            else:
                for k, v in want_attrs:
                    if not attrs.has_key(k) or attrs[k] != v:
                        break
                else:
                    needed.append(prop)

            self.save_info(needed)

    def characters(self, s):
        if self._capture:
            self._text += s

    def save_info(self, needed):
        if not self._capture:
            self._text = ""
        self._capture.append( (needed, len(self._text) ) )

    def get_info(self):
        needed, n = self._capture.pop()
        s = self._text[n:]
        return s, needed

    def endElement(self, tag):
        if not self._rev_tags.has_key(tag):
            return
        text, needed = self.get_info()
        for need in needed:
            self.document.setdefault(need, []).append(text)
