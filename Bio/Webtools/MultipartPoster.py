# Copyright 2012 by Kevin Murray.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This borrrows heavily from the MultipartPostHandler module by Will Holcomb
# <wholcomb@gmail.com>, available at:
# http://pypi.python.org/pypi/MultipartPostHandler/
# MultipartPostHandler was licensed under the LGPL v2.1

"""Module to allow urllib2 to POST to multipart/form-data forms
"""

import urllib2
import mimetools
import mimetypes


class MultipartPoster(urllib2.BaseHandler):
    """
    Handler class to allow urllib2 to POST to multipart/form-data forms.
    """
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
            boundary, data = self.multipart_encode(req_vars, req_files)
            contenttype = 'multipart/form-data; boundary=%s' % boundary
            request.add_unredirected_header('Content-Type', contenttype)
            request.add_data(data)
        return request

    def multipart_encode(self, vars, files, boundary = None, buffer = None):
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
    https_request = http_request

