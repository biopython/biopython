# Copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Contains a private class _SeqDb which is the parent class for a number of
classes that represent APIs for working with major databases.
"""

from urllib import request


class _SeqDb(object):
    """The parent class for a number of classes that represent APIs for working
    with major databases.

    All APIs are able to:
        - construct a URL with given accession code and type

    Some APIs are able to:
        - fetch a sequence with given accession code, type and format
    """

    name = ""
    base_url = ""
    entry_url = ""  # to be formatted; e.g. http://www.rcsb.org/structure/{} ; defaults to base_url + "/" + id
    fetch_url = ""
    fetch_file_format_map = dict(fasta="fasta", genbank="genbank")
    fetch_file_format_default = "genbank"

    @classmethod
    def make_identifier(cls, id_obj):
        return id_obj.id

    @classmethod
    def make_entry_url(cls, id_obj):
        identifier = cls.make_identifier(id_obj)
        if cls.entry_url:
            return cls.entry_url.format(identifier)
        return cls.base_url + "/" + str(identifier)

    @classmethod
    def fetch(cls, accession_code, file_format=None):
        if not cls.fetch_url:
            raise Exception("This database does not support fetching!")
        if not file_format:
            file_format = cls.fetch_file_format_default
        fmt = cls.fetch_file_format_map.get(file_format.lower())
        if not fmt:
            raise Exception(
                f"""Cannot fetch the format {file_format} from this database!
available formats: "{'", "'.join(cls.fetch_file_format_map.keys())}" """
            )
        url = cls.fetch_url.format(
            id=accession_code, format=cls.fetch_file_format_map[file_format],
        )
        return request.urlopen(url).read().decode("utf8")
