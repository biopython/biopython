# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This module handles seqdatabase.INI file.

Classes:
SeqDBRegistry   Holds databases from seqdatabase.INI.

"""
import os

from Bio.config import DBRegistry

# Functions:
# _list_ini_paths    Return a list of standard places for INI files.
# _list_ini_files    Return a list of INI files that exist.
# 
# _load_ini_file     Load an INI file into a set of RegisterableObjects.
# _make_flat_db
# _make_biofetch_db
# _make_biosql_db
#
# _openfu            Open a file or URL.
# _warn              Register a warning message.

class SeqDBRegistry(DBRegistry.DBRegistry):
    """This object implements a dictionary-like interface to sequence
    databases.  To get a list of the databases available, do:
    Bio.seqdb.keys()

    Then, you can access the database using:
    Bio.seqdb[DATABASE_NAME][SEQUENCE_ID]

    """
    def __init__(self, name):
        DBRegistry.DBRegistry.__init__(self, name)
        
    def _load(self, path):
        # path is always None here
        sources = _list_ini_files("seqdatabase.ini")
        for file in sources:
            objects = _load_registry_objects(file)
            if objects:
                for obj in objects:
                    self.register(obj)
                break   # Use the first one that exists.
        else:
            _warn("I could not load any seqdatabase files.")

    def __getitem__(self, name):
        try:
            return DBRegistry.DBRegistry.__getitem__(self, name)
        except KeyError, x:
            raise KeyError, "Unknown database: %s" % str(x)
        

seqdb = SeqDBRegistry("seqdb")

def _warn(message):
    import warnings
    warnings.warn(message)
    
def _load_registry_objects(ini_file):
    # Return a list of RegisterableObjects.
    import _stanzaformat

    try:
        stanzas = _stanzaformat.load(_openfu(ini_file))
    except SyntaxError, x:
        _warn("Can't load seqdb.  Syntax error in %s: %s" % (ini_file, str(x)))
        return None

    # Make sure the file is the right version.
    if stanzas.version > "1.00":
        _warn("I can't handle stanza files with version %s" % stanzas.version)
        return None

    protocol2handler = {
        'flat' : _make_flat_db,
        'biofetch' : _make_biofetch_db,
        'biosql' : _make_biosql_db,
        }

    inidata = []  # list of (section name, section key, dict of tag->value)
    for stanza in stanzas.stanzas:
        section_name, tagvalue_dict = stanza.name, stanza.tag_value_dict
        section_key = section_name.lower()   # case insensitive
        inidata.append((section_name, section_key, tagvalue_dict))

    # Do some checking on each of the stanzas.  If there are errors,
    # then ignore them.
    seen = {}  # Which sections we have already seen.
    i = 0
    while i < len(inidata):
        section_name, section_key, tagvalue_dict = inidata[i]
        # Make sure the stanza has a "protocol".
        if "protocol" not in tagvalue_dict:
            _warn("%s stanza missing 'protocol'.  Skipping" % section_name)
            del inidata[i]
        # Make sure the stanza has a "location".
        elif "location" not in tagvalue_dict:
            _warn("%s stanza missing 'location'.  Skipping" % section_name)
            del inidata[i]
        # Make sure we can handle the "protocol".
        elif tagvalue_dict['protocol'] not in protocol2handler:
            _warn("%s protocol not handled.  Skipping" %
                  tagvalue_dict['protocol'])
            del inidata[i]
        # Make sure this stanza has not already been defined.
        elif section_key in seen:
            _warn("%s stanza already exists.  Skipping" %
                  section_key)
            del inidata[i]
        else:
            seen[section_key] = 1
            i += 1

    # serial_groups is a list of fallback groups.  This is an
    # undocumented feature unsupported in the OBDA spec 1.00!
    registry_objects = []  # list of RegisterableObjects
    serial_groups = []     # list of group_name, obj in group
    for section_name, section_key, tagvalue_dict in inidata:
        handler = protocol2handler.get(tagvalue_dict['protocol'])
        obj = handler(section_name, tagvalue_dict)
        registry_objects.append(obj)

        if tagvalue_dict.has_key("fallback_group"):
            group_name = tagvalue_dict['fallback_group']
            serial_groups.append((group_name, obj))

    # Now make the group objects.
    groups = {}  # name -> DBGroup object
    for group_name, obj in serial_groups:
        if not groups.has_key(group_name):
            groups[group_name] = DBRegistry.DBGroup(
                group_name, behavior="serial")
        groups[group_name].add(obj)
    registry_objects.extend(groups.values())
    return registry_objects

def _make_biofetch_db(name, tagvalue_dict):
    from Martel import Str
    # Make the CGIDB object for the registry.
    params = {}
    params['name'] = name
    params['cgi'] = tagvalue_dict['location']
    dbname = tagvalue_dict.get("dbname", "embl")
    params['params'] = [('style', 'raw'),
                        ('db', dbname),
                        ]
    params['key'] = 'id'
    params['doc'] = "Retrieve sequences from the %s database." % dbname

    params['failure_cases'] = [
        (Str("ERROR 1"), "Unknown database."),
        (Str("ERROR 2"), "Unknown style."),
        (Str("ERROR 3"), "Format not known for database."),
        (Str("ERROR 4"), "ID not found in database."),
        (Str("ERROR 5"), "Too many IDs."),
        ]
    
    # All other params are ignored silently.
    return DBRegistry.CGIDB(**params)

##def _make_biocorba_db(name, tagvalue_dict):
##    """Register a BioCorba style database defined in the registry."""
##    params = {}
##    params['name'] = name
##    params['ior_ref'] = tagvalue_dict['location']
##    return DBRegistry.BioCorbaDB(**params)

def _make_biosql_db(name, tagvalue_dict):
    """Register a BioSQL database defined in the registry."""
    import re
    params = {}
    params['name'] = name

    # Make sure the location has the right format.
    if not re.match(r"[a-zA-Z0-9_]+:\d+$", tagvalue_dict['location']):
        _warn("Invalid location string: %s.  I want <host:port>.  Skipping" %
              tagvalue_dict['location'])
    host, port = tagvalue_dict['location'].split(":")
    params['db_host'] = host
    params['db_port'] = port
    
    params['sql_db'] = tagvalue_dict['biodbname']
    params['db_type'] = tagvalue_dict.get('driver', 'mysql').lower()
    params['db_user'] = tagvalue_dict.get('user', 'root')
    params['db_passwd'] = tagvalue_dict.get('passwd', '')
    params['namespace_db'] = tagvalue_dict['dbname']

    params["doc"] = "Retrieve %s sequences from BioSQL hosted at %s" % (
        tagvalue_dict['dbname'], host)

    return DBRegistry.BioSQLDB(**params)

def _make_flat_db(name, tagvalue_dict):
    """Register a Berkeley or Flat indexed file defined in the registry."""
    params = {}
    params['name'] = name
    params['dbname'] = tagvalue_dict["dbname"]
    params['doc'] = "Retrieve %s sequences from a local database." % \
                    tagvalue_dict["dbname"]
    return DBRegistry.IndexedFileDB(**params)
        
def _openfu(file_or_url):
    """Guess whether this is a file or url and open it."""
    if file_or_url[:4].lower() == 'http':
        import urllib
        return urllib.urlopen(file_or_url)
    # doesn't look like a URL, guess it's a file.
    return open(file_or_url)

def _list_ini_paths():
    """_list_ini_paths() -> list of URL's or paths to search for files.

    The default places to look for registry files are:
    - ${HOME}/.bioinformatics
    - /etc/bioinformatics
    - http://www.open-bio.org/registry

    The OBDA_SEARCH_PATH environment variable, if specified, overrides
    the default.  This should be a "+" separated list of paths or
    URL's.

    """
    if os.environ.has_key("OBDA_SEARCH_PATH"):
        paths = os.environ["OBDA_SEARCH_PATH"].split("+")
    else:
        paths = [
            os.path.join(os.sep, "etc", "bioinformatics"), #/etc/bioinformatics
            "http://www.open-bio.org/registry",
            ]
        # $HOME/.bioinformatics
        if os.environ.has_key("HOME"):
            p = os.path.join(os.environ["HOME"], ".bioinformatics")
            paths.insert(0, p)
    return paths

def _list_ini_files(filename, also_search=[]):
    """_list_ini_files(filename) -> list of files to search (in order)"""
    files = []
    searchpath = _list_ini_paths() + also_search
    for path in searchpath:
        # works for files and urls
        fullname = os.path.join(path, filename)
        # Check to see if this name works.  If so, add it to the list.
        try:
            _openfu(fullname)
        except IOError:
            pass
        else:
            files.append(fullname)
    return files

##def _urlparamdecode(string):
##    """Return a list of (tag, value) from a URL's GET string"""
##    params = []
##    pairs = string.split("&")
##    for tagvalue in pairs:
##        i = tagvalue.find("=")
##        if i >= 0:
##            tag, value = tagvalue[:i], tagvalue[i+1:]
##        else:
##            tag, value = "", tagvalue
##        tag, value = urllib.unquote(tag), urllib.unquote(value)
##        params.append((tag, value))
##    return params
