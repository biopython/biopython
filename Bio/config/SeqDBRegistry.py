# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This module handles seqdatabase.INI file.

Classes:
SeqDBRegistry   Holds databases from seqdatabase.INI.

"""
import os
import urllib

from Bio.config import DBRegistry

# Functions:
# _load_file
# _add_cgi_db
# _add_biocorba_db
# _add_biosql_db
# _add_indexed_db
#
# _open                   Open a file or URL.
#
# _find_ini_files         Look for INI files in a set of standard places.
# _default_search_paths   Return a list of standard places for INI files.
# _pathsep                Return the path separator for this machine.

# _parse_ini              Parse an INI file to (name, list tag/value pairs).
# _urlparamdecode         Decode the parameters of a URL to tag/value pairs.
# _tagvalue_asdict        Convert list of tag/value to a dictionary.


class SeqDBRegistry(DBRegistry.DBRegistry):
    def __init__(self, name):
        DBRegistry.DBRegistry.__init__(self, name)
        
    def _load(self, path):
        # path is always None here
        self.sources = _find_ini_files("seqdatabase.ini")
        for file in self.sources:
            for obj in _load_file(_open(file)):
                self.register(obj)

seqdb = SeqDBRegistry("seqdb")

def _load_file(ini_file):
    # Return a list of RegisterableObjects.
    try:
        inidata = _parse_ini(ini_file)
    except SyntaxError, x:
        return

    protocol2handler = {
        'xembl' : _make_cgi_db,
        'biofetch' : _make_cgi_db,
        'biocorba' : _make_biocorba_db,
        'bsane-corba' : _make_biocorba_db,
        'biosql' : _make_biosql_db,
        'index-flat' : _make_indexed_db,
        'index-berkeleydb' : _make_indexed_db
        }

    registry_objects = []  # list of RegisterableObjects
    serial_groups = []     # list of groupname, obj in group
    for name, tagvalue_pairs in inidata:
        tagvalue_dict = _tagvalue_asdict(tagvalue_pairs)
        protocol = tagvalue_dict['protocol'].lower()
        handler = protocol2handler.get(protocol)
        if handler is None:
            import warnings
            warnings.warn("Protocol %s not recognized" % protocol)
            continue
        obj = handler(name, tagvalue_dict)
        registry_objects.append(obj)
        if tagvalue_dict.has_key("fallback_group"):
            groupname = tagvalue_dict['fallback_group']
            serial_groups.append((groupname, obj))

    # Now make the group objects.
    groups = {}  # name -> DBGroup object
    for groupname, obj in serial_groups:
        if not groups.has_key(groupname):
            groups[groupname] = DBRegistry.DBGroup(
                groupname, behavior="serial")
        groups[groupname].add(obj)
    registry_objects.extend(groups.values())
    return registry_objects

def _make_cgi_db(name, tagvalue_dict):
    # Make the CGIDB object for the registry.
    params = {}
    params['name'] = name
    params['cgi'] = tagvalue_dict['location']
    params['key'] = tagvalue_dict.get('key', '')
    params['params'] = _urlparamdecode(tagvalue_dict.get('params', ''))
    # XXX ignore all other tags silently (?)
    return DBRegistry.CGIDB(**params)

def _make_biocorba_db(name, tagvalue_dict):
    """Register a BioCorba style database defined in the registry.
    """
    params = {}
    params['name'] = name
    params['ior_ref'] = tagvalue_dict['location']
    return DBRegistry.BioCorbaDB(**params)

def _make_biosql_db(name, tagvalue_dict):
    """Register a BioSQL database defined in the registry.
    """
    params = {}
    params['name'] = name
    params['db_host'] = tagvalue_dict.get('host', 'localhost')
    params['db_port'] = tagvalue_dict.get('port', '')
    params['db_user'] = tagvalue_dict.get('user', 'root')
    params['db_passwd'] = tagvalue_dict.get('pass', '')
    params['sql_db'] = tagvalue_dict['dbname']
    params['namespace_db'] = tagvalue_dict['biodbname']
    params['db_type'] = tagvalue_dict.get('driver', 'mysql').lower()
    return DBRegistry.BioSQLDB(**params)

def _make_indexed_db(name, tagvalue_dict):
    """Register a Berkeley or Flat indexed file defined in the registry.
    """
    # get the source object
    params = {}
    params['name'] = name
    params['dbname'] = tagvalue_dict["indexdir"]
    return DBRegistry.IndexedFileDB(**params)
        
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
        return urllib.urlopen(file_or_url)
    # doesn't look like a URL, guess it's a file.
    return open(file_or_url)

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
    # return list of (section name, list of tag-value pairs)
    # raises SyntaxError
    import ConfigParser
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

def _urlparamdecode(string):
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
