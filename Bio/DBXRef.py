class DBXRef:
    def __init__(self, dbname, dbid, reftype = None, negate = 0):
        self.dbname = dbname
        self.dbid = dbid
        self.reftype = reftype
        self.negate = negate

    def __str__(self):
        if self.reftype is None:
            reftype = ""
        else:
            reftype = self.reftype + "="
        s = "%s/%s%s" % (self.dbname, reftype, self.dbid)
        if self.negate:
            s = "not(%s)" % s
        return s
