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


try:
    raise UnitTestExceptions.FailedUnitTest( 'Test!' )
except UnitTestExceptions.FailedUnitTest, failedUnitTest:
    traceBack = sys.exc_info()[ 2 ]
    trace = traceback.extract_tb( traceBack )
    print type( failedUnitTest )
    message = failedUnitTest.message
    print message
    for x in trace:
        for y in x:
	    print y
    defect = Defect.Defect( message, traceBack )

except:
    print 'failed!!!!!!!'
else:
    print 'READY!!'
defect.print_message()
defect.print_trace()
