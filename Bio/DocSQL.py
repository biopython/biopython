#!/usr/bin/env python
#
# Copyright 2002-2003 by Michael Hoffman.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Bio.DocSQL: easy access to DB API databases.

>>> import os
>>> import MySQLdb
>>> from Bio import DocSQL
>>> db=MySQLdb.connect(passwd='', db='test')
>>> class CreatePeople(DocSQL.Create):
...     '''
...     CREATE TEMPORARY TABLE people
...     (id INT UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
...     last_name TINYTEXT,
...     first_name TINYTEXT)
...     '''
...
>>> CreatePeople(connection=db)
CreatePeople(message=Success)
"""

__version__ = "$Revision: 1.13 $"
# $Source: /home/bartek/cvs2bzr/biopython_fastimport/cvs_repo/biopython/Bio/DocSQL.py,v $

import sys

from Bio import MissingPythonDependencyError

try:
    import MySQLdb
except:
    raise MissingPythonDependencyError("Install MySQLdb if you want to use "
                                       "Bio.DocSQL.")

connection = None

class NoInsertionError(Exception):
    pass

def _check_is_public(name):
    if name[:6] == "_names":
        raise AttributeError
    
class QueryRow(list):
    def __init__(self, cursor):
        try:
            row = cursor.fetchone()
            super(QueryRow, self).__init__(row)
        except TypeError:
            raise StopIteration

        object.__setattr__(self, "_names", [x[0] for x in cursor.description]) # FIXME: legacy
        object.__setattr__(self, "_names_hash", {})
        
        for i, name in enumerate(self._names):
            self._names_hash[name] = i

    def __getattr__(self, name):
        _check_is_public(name)
        try:
            return self[self._names_hash[name]]
        except (KeyError, AttributeError):
            raise AttributeError("'%s' object has no attribute '%s'" \
                                 % (self.__class__.__name__, name))

    def __setattr__(self, name, value):
        try:
            self._names_hash
        except AttributeError:
            return object.__setattr__(self, name, value)
            
        _check_is_public(name)
        try:
            index = self._names_hash[name]
            self[index] = value
        except KeyError:
            return object.__setattr__(self, name, value)

class Query(object):
    """
    SHOW TABLES
    """
    MSG_FAILURE = "Failure"
    MSG_SUCCESS = "Success"
    message = "not executed"
    error_message = ""
    prefix = ""
    suffix = ""
    row_class = QueryRow

    def __init__(self, *args, **keywds):
        try:
            self.connection = keywds['connection']
        except KeyError:
            self.connection = connection
        try:
            self.diagnostics = keywds['diagnostics']
        except KeyError:
            self.diagnostics = 0

        self.statement = self.prefix + self.__doc__ + self.suffix
        self.params = args

    def __iter__(self):
        return IterationCursor(self, self.connection)

    def __repr__(self):
        return "%s(message=%s)" % (self.__class__.__name__, self.message)

    def cursor(self):
        return iter(self).cursor

    def dump(self):
        for item in self:
            print item

class QueryGeneric(Query):
    def __init__(self, statement, *args, **keywds):
        Query.__init__(self, *args, **keywds)
        self.statement = statement,

class IterationCursor(object):
    def __init__(self, query, connection=connection):
        if connection is None:
            raise TypeError("database connection is None")
        self.cursor = connection.cursor()
        self.row_class = query.row_class
        if query.diagnostics:
            print >>sys.stderr, query.statement
            print >>sys.stderr, query.params
        self.cursor.execute(query.statement, query.params)

    def next(self):
        return self.row_class(self.cursor)

class QuerySingle(Query, QueryRow):
    ignore_warnings = 0
    def __init__(self, *args, **keywds):
        message = self.MSG_FAILURE
        Query.__init__(self, *args, **keywds)
        try:
            self.single_cursor = Query.cursor(self)
        except MySQLdb.Warning:
            if not self.ignore_warnings:
                raise
        self.row_class.__init__(self, self.cursor())
        object.__setattr__(self, "message", self.MSG_SUCCESS)

    def cursor(self):
        return self.single_cursor

class QueryAll(list, Query):
    def __init__(self, *args, **keywds):
        Query.__init__(self, *args, **keywds)
        list.__init__(self, map(self.process_row, self.cursor().fetchall()))

    def process_row(self, row):
        return row

class QueryAllFirstItem(QueryAll):
    def process_row(self, row):
        return row[0]

class Create(QuerySingle):
    def __init__(self, *args, **keywds):
        try:
            QuerySingle.__init__(self, *args, **keywds)
        except StopIteration:
            self.message = self.MSG_SUCCESS

class Update(Create):
    pass

class Insert(Create):
    MSG_INTEGRITY_ERROR = "Couldn't insert: %s. "
    
    def __init__(self, *args, **keywds):
        try:
            Create.__init__(self, *args, **keywds)
        except MySQLdb.IntegrityError, error_data:
            self.error_message += self.MSG_INTEGRITY_ERROR % error_data[1]
            try:
                self.total_count
            except AttributeError:
                self.total_count = 0
            
            raise MySQLdb.IntegrityError(self.error_message)
            
        self.id = self.cursor().insert_id()
        try:
            self.total_count += self.cursor().rowcount
        except AttributeError:
            self.total_count = self.cursor().rowcount

        if self.cursor().rowcount == 0:
            raise NoInsertionError

def _test(*args, **keywds):
    import doctest, sys
    doctest.testmod(sys.modules[__name__], *args, **keywds)

if __name__ == "__main__":
    if __debug__:
        _test()
