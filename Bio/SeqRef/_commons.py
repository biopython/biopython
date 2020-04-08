# Copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


class SeqId(object):
    def __init__(self, id):
        self.id = id


class SeqRef(object):
    """A type of sequence reference.

    This class defines common features of all sequence reference types, and is
    inherited by specific classes, such as PDB, GB, and EMBL.

    Each class reference is able to:
        - store the type and ID of the reference
        - show databases that use this reference, and construct URLs accordingly
    """

    databases = tuple()

    def __init__(self, id):
        self.id = SeqId(id)
        self.urls = self.get_urls()

    def get_urls(self):
        """Construct URLs for databases available for this type of sequence reference.

        Returns a dict { <name of the database> : <url> }
        """
        return {db.name: db.make_entry_url(self.id) for db in self.databases}
