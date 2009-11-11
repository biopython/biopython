#
#      Restriction Analysis Libraries.
#      Copyright (C) 2004. Frederic Sohm.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

import os

###############################################################################
#                   Configuration of the console.
#
#   Mainly used  by PrintFormat.PrintFormat
#
#   ConsoleWidth : width of the console used default to 80.
#                   should never be less than 60.
#   NameWidth    : space attributed to the name in PrintList method.
#   Indent       : Indent of the second line.
#   MaxSize      : Maximal size of the sequence (default=6:
#                                                -> 99 999 bp + 1 trailing ','
#                  people are unlikely to ask for restriction map of sequences
#                  bigger than 100.000 bp. This is needed to determine the
#                  space to be reserved for sites location.
#                   
#                  MaxSize = 5  =>   9.999 bp
#                  MaxSize = 6  =>  99.999 bp
#                  MaxSize = 7  => 999.999 bp
#   example:
#               
#   <------------ ConsoleWidth --------------->
#   <- NameWidth ->
#   EcoRI         :   1, 45, 50, 300, 400, 650,
#                         700, 1200, 2500.
#                     <-->
#                       Indent
#
ConsoleWidth        =     80                   
NameWidth           =     10
Indent              =     4
MaxSize             =     6 
###############################################################################
#                   Proxies
#
#   Enter here the address of your proxy if any.
#   If you don't use proxy use an empty string
#   i.e.
#   ftp_proxy       =   ''
#                   -> no proxy
#
#   ftp_proxy       =   'http://www.somewhere.something:one_number'
#                   -> www.somewhere.something is the address of the proxy.
#                      one_number is the port number.
#   
ftp_proxy           =   ''
###############################################################################
#                   Rebase ftp location
#
#   Do not modify the addresses.
#
ftp_Rebase          =   'ftp://ftp.neb.com/'
ftp_emb_e           =   ftp_Rebase+'pub/rebase/emboss_e.###'
ftp_emb_s           =   ftp_Rebase+'pub/rebase/emboss_s.###'
ftp_emb_r           =   ftp_Rebase+'pub/rebase/emboss_r.###'
###############################################################################
#                   ftp rebase account.
#
#   In order to update the rebase files, Rana need to connect to the
#   ftp server corresponding.
#
#   the general procedure for accessing a ftp server is generally to
#   connect as anonymous user (rebase_name) and providing your e-mail address
#   as password.
#
#   Therefore, you need to enter your e-mail address in rebase_password.
#   The address will not be send to anyone but is necessary to login the
#   ftp server of rebase when connecting as anonymous user.
#
#   Do not forget to enclose the address between "'".
#
Rebase_name         =   'anonymous'
Rebase_password     =   ''
#Rebase_password     =   'your_address@somewhere.something'


