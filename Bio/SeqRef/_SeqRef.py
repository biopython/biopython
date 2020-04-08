# Copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


class _SeqId(object):
    """A helper class for representing the ID itself.

    Some types of references have IDs with multiple attributes, for example,
    a GenBank code may have a version and a PDB code may have a chain.

    Child classes should always have the attribute "id". If "id" is the only
    attribute, just inherit _SeqId and `pass`.
    """

    def __init__(self, id):
        self.id = id

    def __str__(self):
        return self.id


class SeqRef(object):
    """A type of sequence reference.

    This class defines common features of all sequence reference types, and is
    inherited by specific classes, such as PDB, GB, and EMBL.

    Each class reference is able to:
        - store the type and ID of the reference
        - show databases that use this reference, and construct URLs accordingly
    """

    name = "Type of reference"
    databases = tuple()

    def __init__(self, id):
        self.id = _SeqId(id)
        self.urls = self.get_urls()

    def get_urls(self):
        """Construct URLs for databases available for this type of sequence reference.

        Returns a dict { <name of the database> : <url> }
        """
        return {db.name: db.make_entry_url(self.id) for db in self.databases}

    def __repr__(self):
        """Debugging representation of the reference.

        The returned string can be evaluated to return the same object.
        """
        s = self.__class__.__name__ + "("
        for attr, val in self.id.__dict__.items():
            s += f"{attr}='{val}',"
        return s[:-1] + ")"

    def __str__(self):
        """String representation of the reference.

        Format:
        <type of reference e.g. PDB, GenBank>: <accession code/id>
            <database 1>: <url 1>
            <databse 2>: <url 2>
            ...
        <blank line>
        """
        s = f"{self.name}: {self.id}\n"
        for db, url in self.urls.items():
            s += f"    {db}: {url}\n"
        s += "\n"
        return s
