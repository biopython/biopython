#!/usr/bin/python

# This borrrows heavily from the MultipartPostHandler module by Will Holcomb
# <wholcomb@gmail.com>


# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#

import urllib
import urllib2
import mimetools
import mimetypes
import os
import stat

class Callable:
    def __init__(self, anycallable):
        self.__call__ = anycallable


class MultipartPoster(urllib2.BaseHandler):
    handler_order = urllib2.HTTPHandler.handler_order - 10 # needs to run first

    def http_request(self, request):
        data = request.get_data()
        if data is not None and type(data) != str:
            req_files = []
            req_vars = []
            try:
                 for(key, value) in data.items():
                     if type(value) == file:
                        req_files.append((key, value))
                     else:
                         req_vars.append((key, value))
            except TypeError:
                raise TypeError(
                    "not a valid non-string sequence or mapping object"
                    )
            #if len(req_files) == 0:
            #    # Below, True == do the encoding sequentially
            #    data = urllib.urlencode(req_vars, True) 
            #else:
            boundary, data = self.multipart_encode(req_vars, req_files)
            contenttype = 'multipart/form-data; boundary=%s' % boundary
            request.add_unredirected_header('Content-Type', contenttype)
            request.add_data(data)
        return request

    def multipart_encode(vars, files, boundary = None, buffer = None):
        if boundary is None:
            boundary = mimetools.choose_boundary()
        if buffer is None:
            buffer = ''
        for(key, value) in vars:
            buffer += '--%s\r\n' % boundary
            buffer += 'Content-Disposition: form-data; name="%s"' % key
            buffer += '\r\n\r\n' + value + '\r\n'
        for(key, fd) in files:
            filename = fd.name.split('/')[-1]
            contenttype = mimetypes.guess_type(filename)[0] 
            if not contenttype:
                contenttype = 'application/octet-stream'
            buffer += '--%s\r\n' % boundary
            buffer += 'Content-Disposition: form-data; name="%s"; filename="%s"\r\n' % (key, filename)
            buffer += 'Content-Type: %s\r\n' % contenttype
            fd.seek(0)
            buffer += '\r\n' + fd.read() + '\r\n'
        buffer += '--%s--\r\n\r\n' % boundary
        return boundary, buffer
    multipart_encode = Callable(multipart_encode)

    https_request = http_request







def main():
    import tempfile, sys

    validatorURL = "http://validator.w3.org/check"
    opener = urllib2.build_opener(MultipartPoster)

    def validateFile(url):
        temp = tempfile.mkstemp(suffix=".html")
        os.write(temp[0], opener.open(url).read())
        params = { "ss" : "0",            # show source
                   "doctype" : "Inline",
                   "uploaded_file" : open(temp[1], "rb") }
        print opener.open(validatorURL, params).read()
        os.remove(temp[1])

    if len(sys.argv[1:]) > 0:
        for arg in sys.argv[1:]:
            validateFile(arg)
    else:
        validateFile("http://www.google.com")

if __name__=="__main__":
    main()
