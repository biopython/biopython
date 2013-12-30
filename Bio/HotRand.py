# Copyright 2002 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""handles true random numbers supplied from the web server of
   fourmilab. Based on atmospheric noise.  The motivation is to
   support biosimulations that rely on random numbers.
"""

from __future__ import print_function

#Importing these functions with leading underscore as not intended for reuse
from Bio._py3k import urlopen as _urlopen
from Bio._py3k import urlencode as _urlencode

from Bio import BiopythonDeprecationWarning
import warnings
warnings.warn("The HotRand module is deprecated and likely to be removed in a future release of Biopython. Please use an alternative RNG.", BiopythonDeprecationWarning)


def byte_concat(text):
    val = 0
    numbytes = len(text)
    for i in range(0, numbytes):
        val = val * 256
        # Slice trick for Python 2 and 3 to get single char (byte) string:
        val += ord(text[i:i+1])
    return val


class HotCache(object):

    def __init__( self  ):
#        self.url = 'http://www.fourmilab.ch/cgi-bin/uncgi/Hotbits?num=5000&min=1&max=6&col=1'
        self.url = 'http://www.random.org/cgi-bin/randbyte?'
        self.query = {'nbytes': 128, 'fmt': 'h'}
        self.fill_hot_cache()

    def fill_hot_cache(self):
        url = self.url + _urlencode(self.query)
        fh = _urlopen(url)
        self.hot_cache = fh.read()
        fh.close()

    def next_num(self, num_digits=4):
        cache = self.hot_cache
        # Must explicitly use integer division on python 3
        numbytes = num_digits // 2
        if len(cache) % numbytes != 0:
            print('len_cache is %d' % len(cache))
            raise ValueError
        if cache == '':
            self.fill_hot_cache()
            cache = self.hot_cache
        hexdigits = cache[:numbytes]
        self.hot_cache = cache[numbytes:]
        return byte_concat(hexdigits)


class HotRandom(object):

    def __init__(self):
        self.hot_cache = HotCache( )

    def hot_rand(self, high, low=0):
        span = high - low
        val = self.hot_cache.next_num()
        val = (span * val) >> 16
        val += low
        return val


if __name__ == '__main__':
    hot_random = HotRandom()
    for j in range(0, 130):
        print(hot_random.hot_rand(25))
    nums = ['0000', 'abcd', '1234', '5555', '4321', 'aaaa', 'ffff']
    for num in nums:
        print(int(num, 16))
