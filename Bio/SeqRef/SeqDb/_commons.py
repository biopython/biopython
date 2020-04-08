class _SeqDb(object):
    name = ""
    base_url = ""
    entry_url = ""  # to be formatted; e.g. http://www.rcsb.org/structure/{} ; defaults to base_url + "/" + id

    @classmethod
    def make_identifier(cls, id_obj):
        return id_obj.id

    @classmethod
    def make_entry_url(cls, id_obj):
        identifier = cls.make_identifier(id_obj)
        if cls.entry_url:
            return cls.entry_url.format(identifier)
        return cls.base_url + "/" + str(identifier)
