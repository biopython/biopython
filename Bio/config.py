# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This module is really ugly and things will be to be redone as the
# open-bio registry spec changes.

# XXX
# Are repeated tags allowed?
# What about key="" ?

import os
import urllib

import DBRegistry
import sources

class SeqDatabase(DBRegistry.DBRegistry):
    def __init__(self):
        DBRegistry.DBRegistry.__init__(self)
        
    def _autoload(self):
        if self._autoloaded or self._autoloading:
            return
        self._autoloading = 1
        self._load()
        self._autoloading = 0
        self._autoloaded = 1
    
    def _load(self, more_files=[]):
        files = _find_ini_files("seqdatabase.ini") + more_files
        for file in files:
            self._load_file(_open(file))

    def _load_file(self, ini_file):
        try:
            inidata = _parse_ini(ini_file)
        except SyntaxError, x:
            return
        serial_groups = []   # list of groupname, name
        for name, tagvalue_pairs in inidata:
            tagvalue_dict = _tagvalue_asdict(tagvalue_pairs)
            if tagvalue_dict['protocol'] in ('xembl', 'biofetch'):
                self._add_cgi_db(name, tagvalue_dict)
            else:
                # Ignore unknown protocols silently (?)
                pass
            if tagvalue_dict.has_key("fallback_group"):
                groupname = tagvalue_dict['fallback_group']
                serial_groups.append((groupname, name))
        for groupname, name in serial_groups:
            if self.get(groupname) is None:
                self.register_db(name=groupname, behavior="serial")
            self.group(groupname, name)

    def _add_cgi_db(self, name, tagvalue_dict):
        # Make the Source object for the registry.
        source_params = {}
        source_params['name'] = name
        source_params['cgi'] = tagvalue_dict['location']
        source = sources.CGI(**source_params)

        # Add the DB thing to the database.
        db_params = {}
        db_params['name'] = name
        db_params['source'] = source
        db_params['key'] = tagvalue_dict.get('key', None)
        db_params['params'] = _urldecode(tagvalue_dict.get('params', ''))
        # XXX ignore all other tags silently (?)

        self.register_db(**db_params)



def _pathsep():
    # I would use os.pathsep, except that it breaks on OS X.  It
    # return ":", which the OS doesn't accept...
    seps = ["\\", "/", ":"]
    for sep in seps:
        if os.path.exists(sep):
            return sep
    raise AssertionError, "I couldn't find a working path separator."

def _default_search_paths():
    # Return a list of the places to search for files.  These can be
    # URL's or filenames.
    path = [
        "http://www.open-bio.org/registry",
        os.path.join(_pathsep(), "etc", "bioinformatics"), #/etc/bioinformatics
        ]
    # $HOME/.bioinformatics
    if os.environ.has_key("HOME"):
        p = os.path.join(os.environ["HOME"], ".bioinformatics")
        path.append(p)
    return path

def _open(file_or_url):
    # Guess whether this is a file or url and open it.
    if file_or_url[:4].lower() == 'http':
        open_fn, args, keywds = urllib.urlopen, (file_or_url,), {}
    else:
        open_fn, args, keywds = open, (file_or_url,), {}
    return open_fn(*args, **keywds)

def _find_ini_files(filename, also_search=[]):
    files = []
    searchpath = _default_search_paths() + also_search
    for path in searchpath:
        # works for files and urls
        fullname = os.path.join(path, filename)
        # Check to see if this name works.  If so, add it to the list.
        try:
            _open(fullname)
        except IOError:
            pass
        else:
            files.append(fullname)
    return files

def _parse_ini(handle):
    import ConfigParser
    # return list of (section name, list of tag-value pairs)
    # raises SyntaxError
    parser = ConfigParser.ConfigParser()
    try:
        parser.readfp(handle)
    except ConfigParser.Error, x:
        raise SyntaxError, x
    filedata = []
    for section in parser.sections():
        pairs = []
        for tag in parser.options(section):
            value = parser.get(section, tag)
            pairs.append((tag, value))
        filedata.append((section, pairs))
    return filedata

def _tagvalue_asdict(tagvalue_pairs):
    # tag -> value (tag must be hashable)
    # duplicate tags are overwritten
    dict = {}
    for tag, value in tagvalue_pairs:
        dict[tag] = value
    return dict

def _urldecode(string):
    params = []
    pairs = string.split("&")
    for tagvalue in pairs:
        i = tagvalue.find("=")
        if i >= 0:
            tag, value = tagvalue[:i], tagvalue[i+1:]
        else:
            tag, value = "", tagvalue
        tag, value = urllib.unquote(tag), urllib.unquote(value)
        params.append((tag, value))
    return params
