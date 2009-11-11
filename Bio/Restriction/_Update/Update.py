#!/usr/bin/env python
#
#      Restriction Analysis Libraries.
#      Copyright (C) 2004. Frederic Sohm.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Update the Rebase emboss files used by Restriction to build the
Restriction_Dictionary.py module."""

import os
import sys
import sre
import time
import gzip
from urllib import FancyURLopener

from Bio.Restriction.RanaConfig import *


class RebaseUpdate(FancyURLopener):
    
    def __init__(self, e_mail='', ftpproxy=''):
        """RebaseUpdate([e_mail[, ftpproxy]]) -> new RebaseUpdate instance.

        if e_mail and ftpproxy are not given RebaseUpdate uses the corresponding
        variable from RanaConfig.

        e_mail is the password for the anonymous ftp connection to Rebase.
        ftpproxy is the proxy to use if any."""
        proxy = {'ftp' : ftpproxy or ftp_proxy}
        global Rebase_password
        Rebase_password = e_mail or Rebase_password
        if not Rebase_password:
            raise FtpPasswordError('Rebase')
        if not Rebase_name:
            raise FtpNameError('Rebase')
        FancyURLopener.__init__(self, proxy)

    def prompt_user_passwd(self, host, realm):
        return (Rebase_name, Rebase_password)

    def openRebase(self, name = ftp_Rebase):
        print '\n Please wait, trying to connect to Rebase\n'
        try:
            self.open(name)
        except:
            raise ConnectionError('Rebase')
        return

    def getfiles(self, *files):
        print '\n',
        for file in self.update(*files):
            print 'copying', file
            fn = os.path.basename(file)
            #filename = os.path.join(Rebase, fn)
            filename = os.path.join(os.getcwd(), fn)
            print 'to', filename
            self.retrieve(file, filename)
        self.close()
        return

    def localtime(self):
        t = time.gmtime()
        year = str(t.tm_year)[-1]
        month = str(t.tm_mon)
        if len(month) == 1 : month = '0'+month
        return year+month

    def update(self, *files):
        if not files:
            files = [ftp_emb_e, ftp_emb_s, ftp_emb_r]
        return [x.replace('###', self.localtime()) for x in files]

    def __del__(self):
        if hasattr(self, 'tmpcache') : self.close()
        #
        #   self.tmpcache is created by URLopener.__init__ method.
        #
        return


class FtpNameError(ValueError):

    def __init__(self, which_server):
        print " In order to connect to %s ftp server, you must provide a name.\
        \n Please edit Bio.Restriction.RanaConfig\n" % which_server
        sys.exit()

class FtpPasswordError(ValueError):

    def __init__(self, which_server):
        print "\n\
        \n In order to connect to %s ftp server, you must provide a password.\
        \n Use the --e-mail switch to enter your e-mail address.\
        \n\n" % which_server
        sys.exit()


class ConnectionError(IOError):

    def __init__(self, which_server):
        print '\
        \n Unable to connect to the %s ftp server, make sure your computer\
        \n is connected to the internet and that you have correctly configured\
        \n the ftp proxy.\
        \n Use the --proxy switch to enter the address of your proxy\
        \n' % which_server
        sys.exit()
        


__all__ = ['RebaseUpdate', 'FtpNameError', 'FtpPasswordError']
