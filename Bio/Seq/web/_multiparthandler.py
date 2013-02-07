# Copyright (c) 2o12 Kevin Murray k.d.murray.91@gmail.com
# multiparthandler is licensed under the LGPL v3

# This borrows heavily from the MultipartPostHandler module by Will Holcomb
# <wholcomb@gmail.com>, available at:
# http://pypi.python.org/pypi/MultipartPostHandler/

"""Module to allow urllib2 to POST to multipart/form-data forms
"""
# Python 3 import hackery
try:
    import urllib2 as a_urllib
except ImportError:
    import urllib.request as a_urllib
import uuid
import mimetypes


class multiparthandler(a_urllib.BaseHandler):
    """
    Handler class to allow urllib2 to POST to multipart/form-data forms.
    """

    # Needs to run first
    handler_order = a_urllib.HTTPHandler.handler_order - 10

    def http_request(self, request):
        """Processes request parameters and returns request object.
        """
        data = request.get_data()
        if data is not None and type(data) != str:
            req_params = []
            for(key, value) in data.items():
                req_params.append((key, value))
            boundary, data = self.multipart_encode(req_params)
            contenttype = 'multipart/form-data; charset=UTF-8; boundary=%s'\
             % boundary
            request.add_unredirected_header('Content-Type', contenttype)
            request.add_data(data)
        return request

    def multipart_encode(self, params, boundary=None, data=None):
        """Forms the multipart post request text.
        """
        if boundary is None:
            boundary = uuid.uuid1()
            # using a uuid as boundary, mimetools.choose_boundary() is
            # deprecated
        if data is None:
            data = ''
        for(key, value) in params:
            try:
                # If this works, `value` should be a readable file which needs
                # posting!
                filename = value.name.split('/')[-1]
                contenttype = mimetypes.guess_type(filename)[0]
                if not contenttype:
                    contenttype = 'application/octet-stream'
                data += '--%s\r\n' % boundary
                data += 'Content-Disposition: form-data;'
                data += ' name="%s"; filename="%s"\r\n' % (key, filename)
                data += 'Content-Type: %s\r\n' % contenttype
                value.seek(0)
                data += '\r\n%s\r\n' % value.read()
            except AttributeError:
                # If it's not, then it must be a keyword parameter
                data += '--%s\r\n' % boundary
                data += 'Content-Disposition: form-data; name="%s"' % key
                data += '\r\n\r\n%s\r\n' % value
        # Remember to add a final boundary
        data += '--%s--\r\n\r\n' % boundary
        return boundary, data.encode("UTF-8")

    https_request = http_request
