# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This module is really ugly and has to be redone...
# - Should allow some sort of lazy parsing.
#   (e.g. shouldn't look for seqdatabase files unless people want to)
# - better handling of broken INI files
# - clean handling of INI files
# - cross-platform searching of paths
# - handling of HTTP files

import os,sys
import urllib, urlparse
import ConfigParser

import Bio
import DBRegistry
import sources

def _get_search_path():
    path = []
    # $HOME/.bioinformatics
    if os.environ.has_key("HOME"):
        p = os.path.join(os.environ["HOME"], ".bioinformatics")
        path.append(p)
        
    # /etc/bioinformatics
    # XXX this isn't cross-platform...
    path.append(os.path.join("/", "etc", "bioinformatics"))
    
    # http://www.open-bio.org/registry
    path.append("http://www.open-bio.org/registry")
    return path

def config():
    raise NotImplementedError

def _find_ini_files(filename, also_search=[]):
    files = []
    searchpath = _get_search_path() + also_search
    for path in searchpath:
        if path[:4].lower() == "http":
            # How do I know if something doesn't exist?
            pass
        else:   # assume it's a file
            fullname = os.path.join(path, filename)
            if os.path.exists(fullname):
                files.append(fullname)
    return files

def _config_seqdatabase():
    files = _find_ini_files("seqdatabase.ini")
    registry = DBRegistry.DBRegistry()
    for file in files:
        _config_one_seqdatabase(open(file), registry)
    return registry

# Are repeated tags allowed?
# What about key="" ?

def _config_one_seqdatabase(ini_file, registry):
    inidata = _parse_ini(ini_file)
    serial_groups = []   # list of groupname, name
    for name, tagvalue_pairs in inidata:
        tagvalue_dict = _tagvalue_asdict(tagvalue_pairs)
        if tagvalue_dict['protocol'] in ('xembl', 'biofetch'):
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
            
            registry.register_db(**db_params)
        else:
            # Ignore unknown protocols silently (?)
            pass
        if tagvalue_dict.has_key("fallback_group"):
            groupname = tagvalue_dict['fallback_group']
            serial_groups.append((groupname, name))
    for groupname, name in serial_groups:
        if registry.get(groupname) is None:
            registry.register_db(name=groupname, behavior="serial")
        registry.group(groupname, name)

def _parse_ini(handle):
    # return list of (section name, list of tag-value pairs)
    parser = ConfigParser.ConfigParser()
    parser.readfp(handle)
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

