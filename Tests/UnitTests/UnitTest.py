# UnitTest
#
#
# MODIFICATION NOTES: See bottom of file.

# Copyright (c) 1999 Katharine Lindner
# This module is free software; you can redistribute it and/or modify
# it under the same terms as Python itself.

# IN NO EVENT SHALL THE AUTHOR BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
# SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OF
# THIS CODE, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
#
# THE AUTHOR SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE.  THE CODE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS,
# AND THERE IS NO OBLIGATION WHATSOEVER TO PROVIDE MAINTENANCE,
# SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

import sys

class UnitTest:

    #Title     : run()
    #Usage     : n/a, abstract function
    #Function  : Runs the test cases
    #          :
    #Returns   : test result
    #Argument  : none

    def run( self ):
        pass

    #Title     : count_test_cases()
    #Usage     : n/a, abstract function
    #Function  : Returns the number of test cases
    #          :
    #Returns   : number of test cases
    #Argument  : none

    def count_test_cases( self):
        pass