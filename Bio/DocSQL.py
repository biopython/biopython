#!/usr/bin/env python
#
# Copyright 2002 by Michael Hoffman.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Bio.DocSQL: easy access to DB API databases

>>> import DocSQL, MySQLdb, os
>>> db=MySQLdb.connect(passwd='', db='test')
>>> class CreatePeople(DocSQL.Create):
...     \"""
...     CREATE TEMPORARY TABLE people
...     (id INT UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
...     last_name TINYTEXT,
...     first_name TINYTEXT)
...     \"""
...
>>> CreatePeople(connection=db)
CreatePeople(message=Success)
"""

__version__ = "$Revision: 1.4 $"
# $Source: /home/bartek/cvs2bzr/biopython_fastimport/cvs_repo/biopython/Bio/DocSQL.py,v $

import exceptions
import MySQLdb
import sys

connection = None

class NoInsertionError(exceptions.Exception):
    pass

class QueryRow(list):
    def __init__(self, cursor):
        try:
            list.__init__(self, cursor.fetchone())
        except TypeError:
            raise StopIteration

        self._names = map(lambda x: x[0], cursor.description)

    def __getattr__(self, name):
        if name == "_names":
            raise AttributeError
        try:
            return self[self._names.index(name)]
        except ValueError:
            raise AttributeError, "'%s' object has no attribute '%s'" % (self.__class__.__name__, name)
        except AttributeError:
            raise AttributeError, "'%s' object has no attribute '%s'" % (self.__class__.__name__, name)

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
        self.message = self.MSG_SUCCESS

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
            
            raise MySQLdb.IntegrityError, self.error_message
            
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
