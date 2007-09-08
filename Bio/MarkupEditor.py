"""Simplify adding markup to a piece of text."""

import warnings
warnings.warn("""\
Bio.MarkupEditor is deprecated.
If you use the code in Bio.MarkupEditor, Please get in touch on
the Biopython mailing lists to prevent permanent removal of this
module.""",
              DeprecationWarning)




from xml.sax import saxutils

# Helper functions to turn start/end tags into text
def _start_element(tag, attrs):
    s = "<" + tag
    for name, value in attrs.items():
        s += " %s=%s" % (name, saxutils.quoteattr(value))
    return s + ">"

def _end_element(tag):
    return "</%s>" % tag


class MidPoint:
    def __init__(self):
        self.left = []
        self.right = []

class MarkupEditor:
    def __init__(self, text):
        self.text = text
        self._max = len(text)
        self.midpoints = {}

    def _get(self, pos):
        assert 0 <= pos <= self._max, \
               "position (%d) out of range; 0 <= x <= %d" % self._max
        if self.midpoints.has_key(pos):
            return self.midpoints[pos]
        x = MidPoint()
        self.midpoints[pos] = x
        return x
        
    def insert_text(self, pos, text):
        point = self._get(pos)
        point.right.insert(0, saxutils.escape(text))
        
    def insert_raw_text(self, pos, text):
        point = self._get(pos)
        point.right.insert(0, text)

    def insert_element(self, leftpos, rightpos, tag, attrs = {}):
        strleft = _start_element(tag, attrs)
        strright = _end_element(tag)
        if leftpos == rightpos:
            point = self._get(leftpos)
            point.right.insert(0, strleft + strright)
        else:
            assert leftpos < rightpos, (leftpos, rightpos)
            left = self._get(leftpos)
            right = self._get(rightpos)

            left.right.insert(0, strleft)
            right.left.append(strright)

    def insert_singleton(self, pos, tag, attrs = {}):
        point = self._get(pos)
        s = _start_element(tag, attrs)
        point.right.insert(0, s)

    def to_file(self, outfile):
        items = self.midpoints.items()
        items.sort()
        text = self.text
        prevpos = 0
        for pos, midpoint in items:
            outfile.write(text[prevpos:pos])
            for s in midpoint.left:
                outfile.write(s)
            for s in midpoint.right:
                outfile.write(s)
            prevpos = pos
        outfile.write(text[prevpos:])

def _compare(markup, expect):
    from cStringIO import StringIO
    file = StringIO()
    markup.to_file(file)
    s = file.getvalue()
    assert s == expect, (s, expect)

def test():
    markup = MarkupEditor("01234")
    _compare(markup, "01234")
    markup.insert_text(0, "A")
    _compare(markup, "A01234")
    markup.insert_text(5, "Z")
    _compare(markup, "A01234Z")
    markup.insert_text(5, "Y")
    _compare(markup, "A01234YZ")
    markup.insert_text(5, "<")
    _compare(markup, "A01234&lt;YZ")
    markup.insert_raw_text(1, "BCD")
    _compare(markup, "A0BCD1234&lt;YZ")
    markup.insert_raw_text(3, "<P>")
    _compare(markup, "A0BCD12<P>34&lt;YZ")
    markup.insert_element(1, 3, "a", {"tag": "value"})
    _compare(markup, 'A0<a tag="value">BCD12</a><P>34&lt;YZ')
    markup.insert_element(1, 1, "Q", {"R": "S"})
    _compare(markup, 'A0<Q R="S"></Q><a tag="value">BCD12</a><P>34&lt;YZ')
    markup.insert_singleton(3, "br")
    _compare(markup, 'A0<Q R="S"></Q><a tag="value">BCD12</a><br><P>34&lt;YZ')

if __name__ == "__main__":
    test()
