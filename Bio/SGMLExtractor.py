# Copyright 2002 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code for more fancy file handles.


Classes:
SGMLExtractorHandle     File object that strips tags and returns content from specified
tags blocks.

SGMLExtractor   Object that scans for specified SGML tag pairs, removes any inner tags
and returns the raw content.
For example the object SGMLExtractor( [ 'h1' ] )on the following html file would return
'House that Jack built'
SGMLExtractor( [ 'dt' ] ) would return 'ratcatdogcowmaiden'
SGMLExtractor( [ 'dt', 'dd' ] ) would return 'rat that ate the malttcat ate  the rat' etc

<h1>House that Jack Built</h1>
<dl>
  <dt><big>rat</big></dt>
    <dd><big>ate the malt</big></dd>
  <dt><big>cat</big></dt>
    <dd><big>that ate the rat</big></dd>
  <dt><big>dog</big></dt>
    <dd><big>that worried the dats</big></dd>
  <dt><big>cow</big></dt>
    <dd><big>with crumpled horn</big></dd>
  <dt><big>maiden</big></dt>
    <dd><big>all forlorns</big></dd>
</dl>
"""
import os
import string
import StringIO
import sgmllib


class SGMLExtractorHandle:
    """A Python handle that automatically strips SGML tags and returns data from
    specified tag start and end pairs.

    """
    def __init__(self, handle, tags_of_interest = [] ):
        """SGMLExtractor(handle, tags_of_interest )

        handle is a file handle to SGML-formatted data.
        tags_of_interest is a list of root names for pairs of start and end tags

        """
        self._handle = handle
        self._stripper = SGMLExtractor( tags_of_interest )

    def read(self, *args, **keywds):
        data = apply(self._handle.read, args, keywds)
        return self._stripper.strip(data)

    def readline(self, *args, **keywds):
        line = apply(self._handle.readline, args, keywds)
        return self._stripper.strip(line)

    def readlines(self, *args, **keywds):
        lines = apply(self._handle.readlines, args, keywds)
        for i in range(len(lines)):
            lines[i] = self._stripper.strip(str)
        return lines

    def __getattr__(self, attr):
        return getattr(self._handle, attr)


def is_empty( items ):
    if( len( items ) > 0 ):
        return 0
    else:
        return 1

class SGMLExtractor:
    class LocalParser(sgmllib.SGMLParser):
        def __init__(self, tags_of_interest = [] ):
            sgmllib.SGMLParser.__init__(self)
            self.data = ''
            self._instack = []
            self._tags_of_interest = []
            for tag in tags_of_interest:
                self._tags_of_interest.append( tag.lower() )

        def handle_data(self, data):
            if( not is_empty( self._instack ) ):
                self.data = self.data + data

        def unknown_starttag(self, tag, attrs):
            lower_tag = tag.lower()
            if( lower_tag in self._tags_of_interest ):
                self._instack.append( lower_tag )

        def unknown_endtag(self, tag ):
            if( not is_empty( self._instack ) ):
                open_tag = self._instack.pop()
                try:
                    if( open_tag != tag.lower() ):
                        self._instack.append( open_tag )
                except:
                    print tag


    def __init__(self, tags_of_interest = [] ):
        self._parser = SGMLExtractor.LocalParser( tags_of_interest )

    def strip(self, str):
        """S.strip(str) -> string

        Strip the SGML tags from str.

        """
        if not str:  # empty string, don't do anything.
            return ''
        # I need to make sure that I don't return an empty string if
        # the buffer is not empty.  This can happen if there's a newline
        # character embedded within a tag.  Thus, I'll first check to
        # see if the last character is a newline.  If it is, and it's stripped
        # away, I'll add it back.
        is_newline = str[-1] in ['\n', '\r']

        self._parser.data = ''    # clear the parser's data (don't reset)
        self._parser.feed(str)
        if self._parser.data:
            str = self._parser.data
        elif is_newline:
            str = '\n'
        return str

