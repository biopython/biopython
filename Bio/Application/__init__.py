"""General mechanisms to access applications in biopython.
"""
import os, sys
import StringIO

from Bio import File

def generic_run(commandline):
    """Run an application with the given commandline.

    This expects a pre-built commandline that derives from 
    AbstractCommandline, and returns a ApplicationResult object
    to get results from a program, along with handles of the
    standard output and standard error.

    WARNING - This will read in the full program output into memory!
    This may be in issue when the program write a large amount of
    data to standard output.
    """
    # print str(commandline)

    #Try and use subprocess (available in python 2.4+)
    try :
        import subprocess, sys
        #We don't need to supply any piped input, but we setup the
        #standard input pipe anyway as a work around for a python
        #bug if this is called from a Windows GUI program.  For
        #details, see http://bugs.python.org/issue1124861
        child = subprocess.Popen(str(commandline),
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=(sys.platform!="win32"))
        child.stdin.close()
        r = child.stdout
        e = child.stderr 

        r_out = r.read()
        e_out = e.read()
        r.close()
        e.close()

        # capture error code
        error_code = child.wait()

    except ImportError :
        #For python 2.3 can't use subprocess, using popen2 instead
        #(deprecated in python 2.6)
        import popen2
        if sys.platform[:3]=='win':
            # Windows does not have popen2.Popen3
            r, w, e = popen2.popen3(str(commandline))
        
            r_out = r.read()
            e_out = e.read()
            w.close()
            r.close()
            e.close()

            # No way to get the error code; setting it to a dummy variable
            error_code = 0

        else:
            child = popen2.Popen3(str(commandline), 1)
            # get information and close the files, so if we call this function
            # repeatedly we won't end up with too many open files

            # here are the file descriptors
            r = child.fromchild
            w = child.tochild
            e = child.childerr
        
            r_out = r.read()
            e_out = e.read()
            w.close()
            r.close()
            e.close()
        
            # capture error code
            error_code = os.WEXITSTATUS(child.wait())

    return ApplicationResult(commandline, error_code), \
           File.UndoHandle(StringIO.StringIO(r_out)), \
           File.UndoHandle(StringIO.StringIO(e_out))

class ApplicationResult:
    """Make results of a program available through a standard interface.
    
    This tries to pick up output information available from the program
    and make it available programmatically.
    """
    def __init__(self, application_cl, return_code):
        """Intialize with the commandline from the program.
        """
        self._cl = application_cl

        # provide the return code of the application
        self.return_code = return_code

        # get the application dependent results we can provide
        # right now the only results we handle are output files
        self._results = {}

        for parameter in self._cl.parameters:
            if "file" in parameter.param_types and \
               "output" in parameter.param_types:
                if parameter.is_set:
                    self._results[parameter.names[-1]] = parameter.value

    def get_result(self, output_name):
        """Retrieve result information for the given output.
        """
        return self._results[output_name]

    def available_results(self):
        """Retrieve a list of all available results.
        """
        result_names = self._results.keys()
        result_names.sort()
        return result_names

class AbstractCommandline:
    """Generic interface for running applications from biopython.

    This class shouldn't be called directly; it should be subclassed to
    provide an implementation for a specific application.
    """
    def __init__(self):
        self.program_name = ""
        self.parameters = []
    
    def __str__(self):
        """Make the commandline with the currently set options.
        """
        commandline = "%s " % self.program_name
        for parameter in self.parameters:
            if parameter.is_required and not(parameter.is_set):
                raise ValueError("Parameter %s is not set." % parameter.names)
            if parameter.is_set:
                commandline += str(parameter)

        return commandline

    def set_parameter(self, name, value = None):
        """Set a commandline option for a program.
        """
        set_option = 0
        for parameter in self.parameters:
            if name in parameter.names:
                if value is not None:
                    self._check_value(value, name, parameter.checker_function)
                    parameter.value = value
                parameter.is_set = 1
                set_option = 1

        if set_option == 0:
            raise ValueError("Option name %s was not found." % name)

    def _check_value(self, value, name, check_function):
        """Check whether the given value is valid.

        This uses the passed function 'check_function', which can either
        return a [0, 1] (bad, good) value or raise an error. Either way
        this function will raise an error if the value is not valid, or
        finish silently otherwise.
        """
        if check_function is not None:
            is_good = check_function(value)
            if is_good in [0, 1]: # if we are dealing with a good/bad check
                if not(is_good):
                    raise ValueError(
                            "Invalid parameter value %r for parameter %s" %
                            (value, name))
                    
class _AbstractParameter:
    """A class to hold information about a parameter for a commandline.

    Do not use this directly, instead use one of the subclasses.

    Attributes:

    o names -- a list of string names by which the parameter can be
    referenced (ie. ["-a", "--append", "append"]). The first name in
    the list is considered to be the one that goes on the commandline,
    for those parameters that print the option. The last name in the list
    is assumed to be a "human readable" name describing the option in one
    word.

    o param_type -- a list of string describing the type of parameter, 
    which can help let programs know how to use it. Example descriptions
    include 'input', 'output', 'file'

    o checker_function -- a reference to a function that will determine
    if a given value is valid for this parameter. This function can either
    raise an error when given a bad value, or return a [0, 1] decision on
    whether the value is correct.

    o description -- a description of the option.

    o is_required -- a flag to indicate if the parameter must be set for
    the program to be run.

    o is_set -- if the parameter has been set

    o value -- the value of a parameter
    """
    def __init__(self, names = [], types = [], checker_function = None, 
                 is_required = 0, description = ""):
        self.names = names
        self.param_types = types
        self.checker_function = checker_function
        self.description = description
        self.is_required = is_required

        self.is_set = 0
        self.value = None

class _Option(_AbstractParameter):
    """Represent an option that can be set for a program.

    This holds UNIXish options like --append=yes and -a yes
    """
    def __str__(self):
        """Return the value of this option for the commandline.
        """
        # first deal with long options
        if self.names[0].find("--") >= 0:
            output = "%s" % self.names[0]
            if self.value is not None:
                output += "=%s " % self.value
            else:
                output += " "
        # now short options
        elif self.names[0].find("-") >= 0:
            output = "%s " % self.names[0]
            if self.value is not None:
                output += "%s " % self.value
        else:
            raise ValueError("Unrecognized option type: %s" % self.names[0])

        return output

class _Argument(_AbstractParameter):
    """Represent an argument on a commandline.
    """
    def __str__(self):
        if self.value is not None:
            return "%s " % self.value
        else:
            return " "
