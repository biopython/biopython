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
import traceback
import UnitTestExceptions
import Defect
import UnitTestResults

try:
    raise UnitTestExceptions.FailedUnitTest( 'Test!' )
except UnitTestExceptions.FailedUnitTest, failedUnitTest:
    traceBack = sys.exc_info()[ 2 ]
    trace = traceback.extract_tb( traceBack )
    defect = Defect.Defect( 'failure', traceBack )
    defect.print_message()

    test_results = UnitTestResults.UnitTestResults()
    test_results.add_failure( defect )
    defect = Defect.Defect( 'another failure', traceBack )
    test_results.add_failure( defect )
    defect = Defect.Defect( 'yet another failure', traceBack )
    test_results.add_failure( defect )
    print test_results.count_failures()
    print '\n'
    test_results.print_failures()

    test_results.add_error( 'an error' )
    test_results.add_error( 'another error' )
    test_results.add_error( 'yet another error' )
    test_results.add_error( 'last error' )
    print test_results.count_errors()
    print '\n'
    test_results.print_errors()

    print 'Warnings\n';
    test_results.add_warning( 'a warning' )
    test_results.add_warning( 'another warning' )
    test_results.add_warning( 'yet another warning' )
    test_results.add_warning( 'last warning' )
    print test_results.count_warnings()
    print '\n';
    test_results.print_warnings()

except:
    print 'failed!!!!!!!'
else:
    print 'READY!!'