IN NO EVENT SHALL THE AUTHOR BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OF
THIS CODE, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

THE AUTHOR SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE.  THE CODE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS,
AND THERE IS NO OBLIGATION WHATSOEVER TO PROVIDE MAINTENANCE,
SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
  
  Changes in PyUnit 0.5:

PyUnit 0.5 removes a bug in the inheritance of setup and tear_down.

PyUnit 0.5 includes code to report "success" or "failure" for each test.

PyUnit 0.5 includes new tests, AlwaysPass.py and NoPass.py, to test the bug.

To use:
 PyUnit 0.5 requires the installation of Python, available at www.python.org.
To test the end-to-end functionality of PyUnit:

 open a DOS window and type:
    "python NeverPass.py".

 open a DOS window and type:
    "python AlwaysPass.py".

 open a DOS window and type:
    "python NoPass.py".

NeverPass.py is a script that invokes NeverPassTestCase.py.  The output
of NeverPass.py should be an error message for every test in NeverPassTestCase.py.
The output of AlwaysPass.py should be a success for each test case in
AlwaysPassTestCase.  NoPass.py is similar to NeverPass.py, except that it
tests setup and tear_down.

In addition, each core module in PyUnit has a corresponding test script, with
the name "Test" plus the name of the module. For example, the script,
TestUnitTestSuite.py tests the module, UnitTestSuite.py.

  To create unit tests for a your own module, subclass UnitTestCase.py.
Use NeverPassTestCase.py as an example.  Add your test routines to the
subclassed module.  Make use of assert_condition and assert_equals from
UnitTestCase.py. The name of each test routine must begin with "test_". The
only parameter passed to the routine should be the object reference.  Copy
the invoke_suite routine, from NeverPassTestCase.py to the subclassed module.
I did not include this routine in  UnitTestCase.py, because it would have
created a circular  dependency and  Python handles circular dependencies poorly.
invoke_suite calls UnitTestSuite.build_suite.  UnitTestSuite.build_suite
searches the test class for routines that begin with "test_". build_suite
adds them to a test suite, that invoke_suite runs.

 Last, use NeverPass.py as the model for a script file, that creates and
 invokes the test module.  Substitute the name of your subclassed module for
 NeverPassTestCase.




