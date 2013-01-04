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
    handler_order = urllib2.HTTPHandler.handler_order - 10  # needs to run 1st

    def __init__(self):
        pass

    def http_request(self, request):
        """Processes request parameters and returns request object.
        """
        data = request.get_data()
        if data is not None and type(data) != str:
            req_files = []
            req_params = []
            try:
                for(key, value) in data.items():
                    if type(value) == file:
                        req_files.append((key, value))
                    else:
                        req_params.append((key, value))
            except TypeError:
                raise TypeError(
                    "not a valid non-string sequence or mapping object"
                    )
            boundary, data = self.multipart_encode(req_params, req_files)
            contenttype = 'multipart/form-data; boundary=%s' % boundary
            request.add_unredirected_header('Content-Type', contenttype)
            request.add_data(data)
        return request

    def multipart_encode(self, params, files, boundary=None, data=None):
        """Forms the multipart post request text.
        """
        if boundary is None:
            boundary = mimetools.choose_boundary()
        if data is None:
            data = ''
        for(key, value) in params:
            data += '--%s\r\n' % boundary
            data += 'Content-Disposition: form-data; name="%s"' % key
            data += '\r\n\r\n' + value + '\r\n'
        for(key, file_handle) in files:
            filename = file_handle.name.split('/')[-1]
            contenttype = mimetypes.guess_type(filename)[0]
            if not contenttype:
                contenttype = 'application/octet-stream'
            data += '--%s\r\n' % boundary
            data += 'Content-Disposition: form-data;'  # continued next line
            data += ' name="%s"; filename="%s"\r\n' % (key, filename)
            data += 'Content-Type: %s\r\n' % contenttype
            file_handle.seek(0)
            data += '\r\n' + file_handle.read() + '\r\n'
        data += '--%s--\r\n\r\n' % boundary
        return boundary, data
    https_request = http_request
