# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's ePost,
# as specified by NCBI's DTD file ePost_020511.dtd (2004-09-28 18:47:37)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.



def startElement(self, name, attrs):
    if name=="ePostResult":
        object = {}
	self.path = []
        self.record = object
    else:
        previous = self.path[-1]
        if name=="InvalidIdList":
            object = []
        else:
            object = ""
        if object=="":
            pass
        else:
            previous[name] = object
    self.path.append(object)

def endElement(self, name):
    self.path = self.path[:-1]
    if name=="ERROR":
        raise ValueError(self.content)
    if name in ("Id", "QueryKey", "WebEnv"):
        value = self.content
    else:
        return
    previous = self.path[-1]
    if type(previous)==list:
        previous.append(value)
    elif type(previous)==dict:
        previous[name] = value
