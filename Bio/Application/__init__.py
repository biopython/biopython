# Copyright 2001-2004 Brad Chapman.
# Revisions copyright 2009-2010 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""General mechanisms to access applications in Biopython.

This module is not intended for direct use. It provides the basic objects which
are subclassed by our command line wrappers, such as:

 - Bio.Align.Applications
 - Bio.Blast.Applications
 - Bio.Emboss.Applications
 - Bio.Sequencing.Applications

These modules provide wrapper classes for command line tools to help you
construct command line strings by setting the values of each parameter.
The finished command line strings are then normally invoked via the built-in
Python module subprocess.
"""
import os, sys
import StringIO
import subprocess
import re

#TODO - Remove this hack once we drop Python 2.4 support.
try:
    from subprocess import CalledProcessError as _ProcessCalledError
except:
    #For Python 2.4 use Exception as base class
    _ProcessCalledError = Exception

from Bio import File

#Use this regular expresion to test the property names are going to
#be valid as Python properties or arguments
_re_prop_name = re.compile(r"[a-zA-Z][a-zA-Z0-9_]*")
assert _re_prop_name.match("t")
assert _re_prop_name.match("test")
assert _re_prop_name.match("_test") is None # we don't want private names
assert _re_prop_name.match("-test") is None
assert _re_prop_name.match("test_name")
assert _re_prop_name.match("test2")
#These are reserved names in Python itself,
_reserved_names = ["and", "del", "from", "not", "while", "as", "elif",
                   "global", "or", "with", "assert", "else", "if", "pass",
                   "yield", "break", "except", "import", "print", "class",
                   "exec", "in", "raise", "continue", "finally", "is",
                   "return", "def", "for", "lambda", "try"]
#These are reserved names due to the way the wrappers work
_local_reserved_names = ["set_parameter"]


class ApplicationError(_ProcessCalledError):
    """Raised when an application returns a non-zero exit status.
    
    The exit status will be stored in the returncode attribute, similarly
    the command line string used in the cmd attribute, and (if captured)
    stdout and stderr as strings.
    
    This exception is a subclass of subprocess.CalledProcessError
    (unless run on Python 2.4 where that does not exist).
    
    >>> err = ApplicationError(-11, "helloworld", "", "Some error text")
    >>> err.returncode, err.cmd, err.stdout, err.stderr
    (-11, 'helloworld', '', 'Some error text')
    >>> print err
    Command 'helloworld' returned non-zero exit status -11, 'Some error text'
    
    """
    def __init__(self, returncode, cmd, stdout="", stderr=""):
        self.returncode = returncode
        self.cmd = cmd
        self.stdout = stdout
        self.stderr = stderr
    
    def __str__(self):
        #get first line of any stderr message
        try:
            msg = self.stderr.lstrip().split("\n",1)[0].rstrip()
        except:
            msg = ""
        if msg:
            return "Command '%s' returned non-zero exit status %d, %r" \
                   % (self.cmd, self.returncode, msg)
        else:
            return "Command '%s' returned non-zero exit status %d" \
                   % (self.cmd, self.returncode)
    
    def __repr__(self):
        return "ApplicationError(%i, %s, %s, %s)" \
               % (self.returncode, self.cmd, self.stdout, self.stderr)


class AbstractCommandline(object):
    """Generic interface for constructing command line strings.

    This class shouldn't be called directly; it should be subclassed to
    provide an implementation for a specific application.

    For a usage example we'll show one of the EMBOSS wrappers.  You can set
    options when creating the wrapper object using keyword arguments - or
    later using their corresponding properties:

    >>> from Bio.Emboss.Applications import WaterCommandline
    >>> cline = WaterCommandline(gapopen=10, gapextend=0.5)
    >>> cline
    WaterCommandline(cmd='water', gapopen=10, gapextend=0.5)

    You can instead manipulate the parameters via their properties, e.g.

    >>> cline.gapopen
    10
    >>> cline.gapopen = 20
    >>> cline
    WaterCommandline(cmd='water', gapopen=20, gapextend=0.5)

    You can clear a parameter you have already added by 'deleting' the
    corresponding property:

    >>> del cline.gapopen
    >>> cline.gapopen
    >>> cline
    WaterCommandline(cmd='water', gapextend=0.5)

    Once you have set the parameters you need, turn the object into a string:

    >>> str(cline)
    Traceback (most recent call last):
    ...
    ValueError: You must either set outfile (output filename), or enable filter or stdout (output to stdout).

    In this case the wrapper knows certain arguments are required to construct
    a valid command line for the tool.  For a complete example,

    >>> from Bio.Emboss.Applications import WaterCommandline
    >>> water_cmd = WaterCommandline(gapopen=10, gapextend=0.5)
    >>> water_cmd.asequence = "asis:ACCCGGGCGCGGT"
    >>> water_cmd.bsequence = "asis:ACCCGAGCGCGGT"
    >>> water_cmd.outfile = "temp_water.txt"
    >>> print water_cmd
    water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5
    >>> water_cmd
    WaterCommandline(cmd='water', outfile='temp_water.txt', asequence='asis:ACCCGGGCGCGGT', bsequence='asis:ACCCGAGCGCGGT', gapopen=10, gapextend=0.5)

    You would typically run the command line via a standard Python operating
    system call using the subprocess module for full control. For the simple
    case where you just want to run the command and get the output:

    stdout, stderr = water_cmd()
    """
    #Note the call example above is not a doctest as we can't handle EMBOSS
    #(or any other tool) being missing in the unit tests.
    def __init__(self, cmd, **kwargs):
        """Create a new instance of a command line wrapper object."""
        # Init method - should be subclassed!
        # 
        # The subclass methods should look like this:
        # 
        # def __init__(self, cmd="muscle", **kwargs):
        #     self.parameters = [...]
        #     AbstractCommandline.__init__(self, cmd, **kwargs)
        # 
        # i.e. There should have an optional argument "cmd" to set the location
        # of the executable (with a sensible default which should work if the
        # command is on the path on Unix), and keyword arguments.  It should
        # then define a list of parameters, all objects derived from the base
        # class _AbstractParameter.
        # 
        # The keyword arguments should be any valid parameter name, and will
        # be used to set the associated parameter.
        self.program_name = cmd
        try:
            parameters = self.parameters
        except AttributeError:
            raise AttributeError("Subclass should have defined self.parameters")
        #Create properties for each parameter at run time
        aliases = set()
        for p in parameters:
            for name in p.names:
                if name in aliases:
                    raise ValueError("Parameter alias %s multiply defined" \
                                     % name)
                aliases.add(name)
            name = p.names[-1]
            if _re_prop_name.match(name) is None:
                raise ValueError("Final parameter name %s cannot be used as "
                                 "an argument or property name in python"
                                 % repr(name))
            if name in _reserved_names:
                raise ValueError("Final parameter name %s cannot be used as "
                                 "an argument or property name because it is "
                                 "a reserved word in python" % repr(name))
            if name in _local_reserved_names:
                raise ValueError("Final parameter name %s cannot be used as "
                                 "an argument or property name due to the "
                                 "way the AbstractCommandline class works"
                                 % repr(name))
            #Beware of binding-versus-assignment confusion issues
            def getter(name):
                return lambda x : x._get_parameter(name)
            def setter(name):
                return lambda x, value : x.set_parameter(name, value)
            def deleter(name):
                return lambda x : x._clear_parameter(name)
            doc = p.description
            if isinstance(p, _Switch):
                doc += "\n\nThis property controls the addition of the %s " \
                       "switch, treat this property as a boolean." % p.names[0]
            else:
                doc += "\n\nThis controls the addition of the %s parameter " \
                       "and its associated value.  Set this property to the " \
                       "argument value required." % p.names[0]
            prop = property(getter(name), setter(name), deleter(name), doc)
            setattr(self.__class__, name, prop) #magic!
        for key, value in kwargs.iteritems():
            self.set_parameter(key, value)
    
    def _validate(self):
        """Make sure the required parameters have been set (PRIVATE).

        No return value - it either works or raises a ValueError.

        This is a separate method (called from __str__) so that subclasses may
        override it.
        """
        for p in self.parameters:
            #Check for missing required parameters:
            if p.is_required and not(p.is_set):
                raise ValueError("Parameter %s is not set." \
                                 % p.names[-1])
            #Also repeat the parameter validation here, just in case?

    def __str__(self):
        """Make the commandline string with the currently set options.

        e.g.
        >>> from Bio.Emboss.Applications import WaterCommandline
        >>> cline = WaterCommandline(gapopen=10, gapextend=0.5)
        >>> cline.asequence = "asis:ACCCGGGCGCGGT"
        >>> cline.bsequence = "asis:ACCCGAGCGCGGT"
        >>> cline.outfile = "temp_water.txt"
        >>> print cline
        water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5
        >>> str(cline)
        'water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5'
        """
        self._validate()
        commandline = "%s " % self.program_name
        for parameter in self.parameters:
            if parameter.is_set:
                #This will include a trailing space:
                commandline += str(parameter)
        return commandline.strip() # remove trailing space

    def __repr__(self):
        """Return a representation of the command line object for debugging.

        e.g.
        >>> from Bio.Emboss.Applications import WaterCommandline
        >>> cline = WaterCommandline(gapopen=10, gapextend=0.5)
        >>> cline.asequence = "asis:ACCCGGGCGCGGT"
        >>> cline.bsequence = "asis:ACCCGAGCGCGGT"
        >>> cline.outfile = "temp_water.txt"
        >>> print cline
        water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5
        >>> cline
        WaterCommandline(cmd='water', outfile='temp_water.txt', asequence='asis:ACCCGGGCGCGGT', bsequence='asis:ACCCGAGCGCGGT', gapopen=10, gapextend=0.5)
        """
        answer = "%s(cmd=%s" % (self.__class__.__name__, repr(self.program_name))
        for parameter in self.parameters:
            if parameter.is_set:
                if isinstance(parameter, _Switch):
                    answer += ", %s=True" % parameter.names[-1]
                else:
                    answer += ", %s=%s" \
                              % (parameter.names[-1], repr(parameter.value))
        answer += ")"
        return answer

    def _get_parameter(self, name):
        """Get a commandline option value."""
        for parameter in self.parameters:
            if name in parameter.names:
                if isinstance(parameter, _Switch):
                    return parameter.is_set
                else:
                    return parameter.value
        raise ValueError("Option name %s was not found." % name)

    def _clear_parameter(self, name):
        """Reset or clear a commandline option value."""
        cleared_option = False
        for parameter in self.parameters:
            if name in parameter.names:
                parameter.value = None
                parameter.is_set = False
                cleared_option = True
        if not cleared_option:
            raise ValueError("Option name %s was not found." % name)
        
    def set_parameter(self, name, value = None):
        """Set a commandline option for a program.
        """
        set_option = False
        for parameter in self.parameters:
            if name in parameter.names:
                if isinstance(parameter, _Switch):
                    if value is None:
                        import warnings
                        warnings.warn("For a switch type argument like %s, "
                                      "we expect a boolean.  None is treated "
                                      "as FALSE!" % parameter.names[-1])
                    parameter.is_set = bool(value)
                    set_option = True
                else:
                    if value is not None:
                        self._check_value(value, name, parameter.checker_function)
                        parameter.value = value
                    parameter.is_set = True
                    set_option = True
        if not set_option:
            raise ValueError("Option name %s was not found." % name)

    def _check_value(self, value, name, check_function):
        """Check whether the given value is valid.

        No return value - it either works or raises a ValueError.

        This uses the passed function 'check_function', which can either
        return a [0, 1] (bad, good) value or raise an error. Either way
        this function will raise an error if the value is not valid, or
        finish silently otherwise.
        """
        if check_function is not None:
            is_good = check_function(value) #May raise an exception
            assert is_good in [0,1,True,False]
            if not is_good:
                raise ValueError("Invalid parameter value %r for parameter %s" \
                                 % (value, name))

    def __setattr__(self, name, value):
        """Set attribute name to value (PRIVATE).

        This code implements a workaround for a user interface issue.
        Without this __setattr__ attribute-based assignment of parameters
        will silently accept invalid parameters, leading to known instances
        of the user assuming that parameters for the application are set,
        when they are not.
        
        >>> from Bio.Emboss.Applications import WaterCommandline
        >>> cline = WaterCommandline(gapopen=10, gapextend=0.5, stdout=True)
        >>> cline.asequence = "a.fasta"
        >>> cline.bsequence = "b.fasta"
        >>> cline.csequence = "c.fasta"
        Traceback (most recent call last):
        ...
        ValueError: Option name csequence was not found.
        >>> print cline
        water -stdout -asequence=a.fasta -bsequence=b.fasta -gapopen=10 -gapextend=0.5

        This workaround uses a whitelist of object attributes, and sets the
        object attribute list as normal, for these.  Other attributes are
        assumed to be parameters, and passed to the self.set_parameter method
        for validation and assignment.
        """
        if name in ['parameters', 'program_name']: # Allowed attributes
            self.__dict__[name] = value
        else:
            self.set_parameter(name, value)  # treat as a parameter
    
    def __call__(self, stdin=None, stdout=True, stderr=True,
                 cwd=None, env=None):
        """Executes the command, waits for it to finish, and returns output.
        
        Runs the command line tool and waits for it to finish. If it returns
        a non-zero error level, an exception is raised. Otherwise two strings
        are returned containing stdout and stderr.
        
        The optional stdin argument should be a string of data which will be
        passed to the tool as standard input.

        The optional stdout and stderr argument are treated as a booleans, and
        control if the output should be captured (True, default), or ignored
        by sending it to /dev/null to avoid wasting memory (False). In the
        later case empty string(s) are returned.

        The optional cwd argument is a string giving the working directory to
        to run the command from. See Python's subprocess module documentation
        for more details.

        The optional env argument is a dictionary setting the environment
        variables to be used in the new process. By default the current
        process' environment variables are used. See Python's subprocess
        module documentation for more details.

        Default example usage:

        from Bio.Emboss.Applications import WaterCommandline
        water_cmd = WaterCommandline(gapopen=10, gapextend=0.5,
                                     stdout=True, auto=True,
                                     asequence="a.fasta", bsequence="b.fasta")
        print "About to run:\n%s" % water_cmd
        std_output, err_output = water_cmd()

        This functionality is similar to subprocess.check_output() added in
        Python 2.7. In general if you require more control over running the
        command, use subprocess directly.
        
        As of Biopython 1.56, when the program called returns a non-zero error
        level, a custom ApplicationError exception is raised. This includes
        any stdout and stderr strings captured as attributes of the exception
        object, since they may be useful for diagnosing what went wrong.
        """
        if stdout:
            stdout_arg = subprocess.PIPE
        else:
            stdout_arg = open(os.devnull)
        if stderr:
            stderr_arg = subprocess.PIPE
        else:
            stderr_arg = open(os.devnull)
        #We may not need to supply any piped input, but we setup the
        #standard input pipe anyway as a work around for a python
        #bug if this is called from a Windows GUI program.  For
        #details, see http://bugs.python.org/issue1124861
        #
        #Using universal newlines is important on Python 3, this
        #gives unicode handles rather than bytes handles.
        child_process = subprocess.Popen(str(self), stdin=subprocess.PIPE,
                                         stdout=stdout_arg, stderr=stderr_arg,
                                         universal_newlines=True,
                                         cwd=cwd, env=env,
                                         shell=(sys.platform!="win32"))
        #Use .communicate as can get deadlocks with .wait(), see Bug 2804
        stdout_str, stderr_str = child_process.communicate(stdin)
        if not stdout: assert not stdout_str
        if not stderr: assert not stderr_str
        return_code = child_process.returncode
        if return_code:
            raise ApplicationError(return_code, str(self),
                                   stdout_str, stderr_str)
        return stdout_str, stderr_str


class _AbstractParameter:
    """A class to hold information about a parameter for a commandline.

    Do not use this directly, instead use one of the subclasses.
    """
    def __init__(self):
        raise NotImplementedError

    def __str__(self):
        raise NotImplementedError

class _Option(_AbstractParameter):
    """Represent an option that can be set for a program.

    This holds UNIXish options like --append=yes and -a yes,
    where a value (here "yes") is generally expected.

    For UNIXish options like -kimura in clustalw which don't
    take a value, use the _Switch object instead.

    Attributes:

    o names -- a list of string names by which the parameter can be
    referenced (ie. ["-a", "--append", "append"]). The first name in
    the list is considered to be the one that goes on the commandline,
    for those parameters that print the option. The last name in the list
    is assumed to be a "human readable" name describing the option in one
    word.

    o description -- a description of the option.

    o filename -- True if this argument is a filename and should be
    automatically quoted if it contains spaces.

    o checker_function -- a reference to a function that will determine
    if a given value is valid for this parameter. This function can either
    raise an error when given a bad value, or return a [0, 1] decision on
    whether the value is correct.

    o equate -- should an equals sign be inserted if a value is used?

    o is_required -- a flag to indicate if the parameter must be set for
    the program to be run.

    o is_set -- if the parameter has been set

    o value -- the value of a parameter
    """
    def __init__(self, names, description, filename=False, checker_function=None,
                 is_required=False, equate=True):
        self.names = names
        assert isinstance(description, basestring), \
               "%r for %s" % (description, names[-1])
        self.is_filename = filename
        self.checker_function = checker_function
        self.description = description
        self.equate = equate
        self.is_required = is_required

        self.is_set = False
        self.value = None

    def __str__(self):
        """Return the value of this option for the commandline.

        Includes a trailing space.
        """
        # Note: Before equate was handled explicitly, the old
        # code would do either "--name " or "--name=value ",
        # or " -name " or " -name value ".  This choice is now
        # now made explicitly when setting up the option.
        if self.value is None:
            return "%s " % self.names[0]
        if self.is_filename:
            v = _escape_filename(self.value)
        else:
            v = str(self.value)
        if self.equate:
            return "%s=%s " % (self.names[0], v)
        else:
            return "%s %s " % (self.names[0], v)

class _Switch(_AbstractParameter):
    """Represent an optional argument switch for a program.

    This holds UNIXish options like -kimura in clustalw which don't
    take a value, they are either included in the command string
    or omitted.

    o names -- a list of string names by which the parameter can be
    referenced (ie. ["-a", "--append", "append"]). The first name in
    the list is considered to be the one that goes on the commandline,
    for those parameters that print the option. The last name in the list
    is assumed to be a "human readable" name describing the option in one
    word.

    o description -- a description of the option.

    o is_set -- if the parameter has been set

    NOTE - There is no value attribute, see is_set instead,
    """
    def __init__(self, names, description):
        self.names = names
        self.description = description
        self.is_set = False
        self.is_required = False

    def __str__(self):
        """Return the value of this option for the commandline.

        Includes a trailing space.
        """
        assert not hasattr(self, "value")
        if self.is_set:
            return "%s " % self.names[0]
        else:
            return ""

class _Argument(_AbstractParameter):
    """Represent an argument on a commandline.
    """
    def __init__(self, names, description, filename=False,
                 checker_function=None, is_required=False):
        self.names = names
        assert isinstance(description, basestring), \
               "%r for %s" % (description, names[-1])
        self.is_filename = filename
        self.checker_function = checker_function
        self.description = description
        self.is_required = is_required
        self.is_set = False
        self.value = None

    def __str__(self):
        if self.value is None:
            return " "
        elif self.is_filename:
            return "%s " % _escape_filename(self.value)
        else:
            return "%s " % self.value

def _escape_filename(filename):
    """Escape filenames with spaces by adding quotes (PRIVATE).

    Note this will not add quotes if they are already included:
    
    >>> print _escape_filename('example with spaces')
    "example with spaces"
    >>> print _escape_filename('"example with spaces"')
    "example with spaces"
    """
    #Is adding the following helpful
    #if os.path.isfile(filename):
    #    #On Windows, if the file exists, we can ask for
    #    #its alternative short name (DOS style 8.3 format)
    #    #which has no spaces in it.  Note that this name
    #    #is not portable between machines, or even folder!
    #    try:
    #        import win32api
    #        short = win32api.GetShortPathName(filename)
    #        assert os.path.isfile(short)
    #        return short
    #    except ImportError:
    #        pass
    if " " not in filename:
        return filename
    #We'll just quote it - works on Windows, Mac OS X etc
    if filename.startswith('"') and filename.endswith('"'):
        #Its already quoted
        return filename
    else:
        return '"%s"' % filename

def _test():
    """Run the Bio.Application module's doctests."""
    import doctest
    doctest.testmod(verbose=1)

if __name__ == "__main__":
    #Run the doctests
    _test()
