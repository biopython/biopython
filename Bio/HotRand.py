# Copyright 2002 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""handles true random numbers supplied from the the web server of fourmilab. Based on
atmospheric noise.  The motivation is to support biosimulations that rely on random numbers.
"""

import os
import sys
import string
from urllib import FancyURLopener
from urllib import urlencode

from Bio.SGMLExtractor import SGMLExtractorHandle


def hex_convert( text ):
    num_val = 0l
    text = text.lower()
    for letter in text:
        hex_digit = string.hexdigits.find( letter )
        if( hex_digit < 0 ): raise ValueError
        num_val = ( num_val * 16 ) + hex_digit
    return num_val

class HotCache:

    def __init__( self  ):
#        self.url = 'http://www.fourmilab.ch/cgi-bin/uncgi/Hotbits?num=5000&min=1&max=6&col=1'
        self.url = 'http://www.random.org/cgi-bin/randbyte?nbytes=5000&format=h'
        self.query = { 'nbytes' : 128, 'fmt' : 'hex' }
        self.fill_hot_cache()

    def fill_hot_cache( self ):
        bases   = [ 'a', 'g', 'c', 't' ]
        url = self.url + urlencode( self.query )
        url_opener = FancyURLopener( )
        fh = url_opener.open( url )
        hot_rand_handle = SGMLExtractorHandle( fh, [ 'pre', ] )
        hot_cache = hot_rand_handle.read()
        lines = hot_cache.splitlines()
        hot_cache = ''.join( lines )
        hot_cache = hot_cache.replace( ' ', '' )
        self.hot_cache = hot_cache.strip()
        fh.close()
        return self.hot_cache

    def next_num( self, num_digits = 4 ):
        cache = self.hot_cache
        if( len( cache ) % num_digits != 0 ):
            print 'len_cache is %d' % len( cache )
            raise ValueError
        if( cache == '' ):
            self.fill_hot_cache()
            cache = self.hot_cache
        hexdigits = cache[ :num_digits ]
        self.hot_cache = cache[ num_digits: ]
        return hex_convert( hexdigits )


class HotRandom:

    def __init__( self ):
        self.hot_cache = HotCache( )

    def hot_rand( self, high, low = 0 ):
        span = high - low
        val = self.hot_cache.next_num()
        val = ( span * val ) >> 16
        val = val + low
        return val


if( __name__ == '__main__' ):
    hot_random = HotRandom()
    for j in range ( 0, 130 ):
        print hot_random.hot_rand( 25 )
    nums = [ '0000', 'abcd', '1234', '5555', '4321', 'aaaa', 'ffff' ]
    for num in nums:
        print hex_convert( num )



