# Copyright 2001-2004 Brad Chapman.
# Revisions copyright 2009-2013 by Peter Cock.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
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
import os
import platform
import sys
import subprocess
import re


# Use this regular expression to test the property names are going to
# be valid as Python properties or arguments
_re_prop_name = re.compile(r"^[a-zA-Z][a-zA-Z0-9_]*$")
assert _re_prop_name.match("t")
assert _re_prop_name.match("test")
assert _re_prop_name.match("_test") is None  # we don't want private names
assert _re_prop_name.match("-test") is None
assert _re_prop_name.match("any-hyphen") is None
assert _re_prop_name.match("underscore_ok")
assert _re_prop_name.match("test_name")
assert _re_prop_name.match("test2")
# These are reserved names in Python itself,
_reserved_names = [
    "and",
    "del",
    "from",
    "not",
    "while",
    "as",
    "elif",
    "global",
    "or",
    "with",
    "assert",
    "else",
    "if",
    "pass",
    "yield",
    "break",
    "except",
    "import",
    "print",
    "class",
    "exec",
    "in",
    "raise",
    "continue",
    "finally",
    "is",
    "return",
    "def",
    "for",
    "lambda",
    "try",
]
# These are reserved names due to the way the wrappers work
_local_reserved_names = ["set_parameter"]


class ApplicationError(subprocess.CalledProcessError):
    """Raised when an application returns a non-zero exit status.

    The exit status will be stored in the returncode attribute, similarly
    the command line string used in the cmd attribute, and (if captured)
    stdout and stderr as strings.

    This exception is a subclass of subprocess.CalledProcessError.

    >>> err = ApplicationError(-11, "helloworld", "", "Some error text")
    >>> err.returncode, err.cmd, err.stdout, err.stderr
    (-11, 'helloworld', '', 'Some error text')
    >>> print(err)
    Non-zero return code -11 from 'helloworld', message 'Some error text'

    """

    def __init__(self, returncode, cmd, stdout="", stderr=""):
        """Initialize."""
        self.returncode = returncode
        self.cmd = cmd
        self.stdout = stdout
        self.stderr = stderr

    def __str__(self):
        """Format the error as a string."""
        # get first line of any stderr message
        try:
            msg = self.stderr.lstrip().split("\n", 1)[0].rstrip()
        except Exception:  # TODO, ValueError? AttributeError?
            msg = ""
        if msg:
            return "Non-zero return code %d from %r, message %r" % (
                self.returncode,
                self.cmd,
                msg,
            )
        else:
            return "Non-zero return code %d from %r" % (self.returncode, self.cmd)

    def __repr__(self):
        """Represent the error as a string."""
        return "ApplicationError(%i, %s, %s, %s)" % (
            self.returncode,
            self.cmd,
            self.stdout,
            self.stderr,
        )


class AbstractCommandline:
    r"""Generic interface for constructing command line strings.

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

    Once you have set the parameters you need, you can turn the object into
    a string (e.g. to log the command):

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
    >>> print(water_cmd)
    water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5
    >>> water_cmd
    WaterCommandline(cmd='water', outfile='temp_water.txt', asequence='asis:ACCCGGGCGCGGT', bsequence='asis:ACCCGAGCGCGGT', gapopen=10, gapextend=0.5)

    You would typically run the command line via a standard Python operating
    system call using the subprocess module for full control. For the simple
    case where you just want to run the command and get the output:

    stdout, stderr = water_cmd()

    Note that by default we assume the underlying tool is installed on the
    system $PATH environment variable. This is normal under Linux/Unix, but
    may need to be done manually under Windows. Alternatively, you can specify
    the full path to the binary as the first argument (cmd):

    >>> from Bio.Emboss.Applications import WaterCommandline
    >>> water_cmd = WaterCommandline(r"C:\Program Files\EMBOSS\water.exe",
    ...                              gapopen=10, gapextend=0.5,
    ...                              asequence="asis:ACCCGGGCGCGGT",
    ...                              bsequence="asis:ACCCGAGCGCGGT",
    ...                              outfile="temp_water.txt")
    >>> print(water_cmd)
    "C:\Program Files\EMBOSS\water.exe" -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5

    Notice that since the path name includes a space it has automatically
    been quoted.

    """

    # TODO - Replace the above example since EMBOSS doesn't work properly
    # if installed into a folder with a space like "C:\Program Files\EMBOSS"
    #
    # Note the call example above is not a doctest as we can't handle EMBOSS
    # (or any other tool) being missing in the unit tests.

    parameters = None  # will be a list defined in subclasses

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
            raise AttributeError(
                "Subclass should have defined self.parameters"
            ) from None
        # Create properties for each parameter at run time
        aliases = set()
        for p in parameters:
            if not p.names:
                if not isinstance(p, _StaticArgument):
                    raise TypeError("Expected %r to be of type _StaticArgument" % p)
                continue
            for name in p.names:
                if name in aliases:
                    raise ValueError("Parameter alias %s multiply defined" % name)
                aliases.add(name)
            name = p.names[-1]
            if _re_prop_name.match(name) is None:
                raise ValueError(
                    "Final parameter name %r cannot be used as "
                    "an argument or property name in python" % name
                )
            if name in _reserved_names:
                raise ValueError(
                    "Final parameter name %r cannot be used as "
                    "an argument or property name because it is "
                    "a reserved word in python" % name
                )
            if name in _local_reserved_names:
                raise ValueError(
                    "Final parameter name %r cannot be used as "
                    "an argument or property name due to the "
                    "way the AbstractCommandline class works" % name
                )

            # Beware of binding-versus-assignment confusion issues
            def getter(name):
                return lambda x: x._get_parameter(name)

            def setter(name):
                return lambda x, value: x.set_parameter(name, value)

            def deleter(name):
                return lambda x: x._clear_parameter(name)

            doc = p.description
            if isinstance(p, _Switch):
                doc += (
                    "\n\nThis property controls the addition of the %s "
                    "switch, treat this property as a boolean." % p.names[0]
                )
            else:
                doc += (
                    "\n\nThis controls the addition of the %s parameter "
                    "and its associated value.  Set this property to the "
                    "argument value required." % p.names[0]
                )
            prop = property(getter(name), setter(name), deleter(name), doc)
            setattr(self.__class__, name, prop)  # magic!
        for key, value in kwargs.items():
            self.set_parameter(key, value)

    def _validate(self):
        """Make sure the required parameters have been set (PRIVATE).

        No return value - it either works or raises a ValueError.

        This is a separate method (called from __str__) so that subclasses may
        override it.
        """
        for p in self.parameters:
            # Check for missing required parameters:
            if p.is_required and not (p.is_set):
                raise ValueError("Parameter %s is not set." % p.names[-1])
            # Also repeat the parameter validation here, just in case?

    def __str__(self):
        """Make the commandline string with the currently set options.

        e.g.

        >>> from Bio.Emboss.Applications import WaterCommandline
        >>> cline = WaterCommandline(gapopen=10, gapextend=0.5)
        >>> cline.asequence = "asis:ACCCGGGCGCGGT"
        >>> cline.bsequence = "asis:ACCCGAGCGCGGT"
        >>> cline.outfile = "temp_water.txt"
        >>> print(cline)
        water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5
        >>> str(cline)
        'water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5'
        """
        self._validate()
        commandline = "%s " % _escape_filename(self.program_name)
        for parameter in self.parameters:
            if parameter.is_set:
                # This will include a trailing space:
                commandline += str(parameter)
        return commandline.strip()  # remove trailing space

    def __repr__(self):
        """Return a representation of the command line object for debugging.

        e.g.

        >>> from Bio.Emboss.Applications import WaterCommandline
        >>> cline = WaterCommandline(gapopen=10, gapextend=0.5)
        >>> cline.asequence = "asis:ACCCGGGCGCGGT"
        >>> cline.bsequence = "asis:ACCCGAGCGCGGT"
        >>> cline.outfile = "temp_water.txt"
        >>> print(cline)
        water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5
        >>> cline
        WaterCommandline(cmd='water', outfile='temp_water.txt', asequence='asis:ACCCGGGCGCGGT', bsequence='asis:ACCCGAGCGCGGT', gapopen=10, gapextend=0.5)
        """
        answer = "%s(cmd=%r" % (self.__class__.__name__, self.program_name)
        for parameter in self.parameters:
            if parameter.is_set:
                if isinstance(parameter, _Switch):
                    answer += ", %s=True" % parameter.names[-1]
                else:
                    answer += ", %s=%r" % (parameter.names[-1], parameter.value)
        answer += ")"
        return answer

    def _get_parameter(self, name):
        """Get a commandline option value (PRIVATE)."""
        for parameter in self.parameters:
            if name in parameter.names:
                if isinstance(parameter, _Switch):
                    return parameter.is_set
                else:
                    return parameter.value
        raise ValueError("Option name %s was not found." % name)

    def _clear_parameter(self, name):
        """Reset or clear a commandline option value (PRIVATE)."""
        cleared_option = False
        for parameter in self.parameters:
            if name in parameter.names:
                parameter.value = None
                parameter.is_set = False
                cleared_option = True
        if not cleared_option:
            raise ValueError("Option name %s was not found." % name)

    def set_parameter(self, name, value=None):
        """Set a commandline option for a program (OBSOLETE).

        Every parameter is available via a property and as a named
        keyword when creating the instance. Using either of these is
        preferred to this legacy set_parameter method which is now
        OBSOLETE, and likely to be DEPRECATED and later REMOVED in
        future releases.
        """
        set_option = False
        for parameter in self.parameters:
            if name in parameter.names:
                if isinstance(parameter, _Switch):
                    if value is None:
                        import warnings

                        warnings.warn(
                            "For a switch type argument like %s, "
                            "we expect a boolean.  None is treated "
                            "as FALSE!" % parameter.names[-1]
                        )
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
        """Check whether the given value is valid (PRIVATE).

        No return value - it either works or raises a ValueError.

        This uses the passed function 'check_function', which can either
        return a [0, 1] (bad, good) value or raise an error. Either way
        this function will raise an error if the value is not valid, or
        finish silently otherwise.
        """
        if check_function is not None:
            is_good = check_function(value)  # May raise an exception
            if is_good not in [0, 1, True, False]:
                raise ValueError(
                    "Result of check_function: %r is of an unexpected value" % is_good
                )
            if not is_good:
                raise ValueError(
                    "Invalid parameter value %r for parameter %s" % (value, name)
                )

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
        >>> print(cline)
        water -stdout -asequence=a.fasta -bsequence=b.fasta -gapopen=10 -gapextend=0.5

        This workaround uses a whitelist of object attributes, and sets the
        object attribute list as normal, for these.  Other attributes are
        assumed to be parameters, and passed to the self.set_parameter method
        for validation and assignment.
        """
        if name in ["parameters", "program_name"]:  # Allowed attributes
            self.__dict__[name] = value
        else:
            self.set_parameter(name, value)  # treat as a parameter

    def __call__(self, stdin=None, stdout=True, stderr=True, cwd=None, env=None):
        """Execute command, wait for it to finish, return (stdout, stderr).

        Runs the command line tool and waits for it to finish. If it returns
        a non-zero error level, an exception is raised. Otherwise two strings
        are returned containing stdout and stderr.

        The optional stdin argument should be a string of data which will be
        passed to the tool as standard input.

        The optional stdout and stderr argument may be filenames (string),
        but otherwise are treated as a booleans, and control if the output
        should be captured as strings (True, default), or ignored by sending
        it to /dev/null to avoid wasting memory (False). If sent to a file
        or ignored, then empty string(s) are returned.

        The optional cwd argument is a string giving the working directory
        to run the command from. See Python's subprocess module documentation
        for more details.

        The optional env argument is a dictionary setting the environment
        variables to be used in the new process. By default the current
        process' environment variables are used. See Python's subprocess
        module documentation for more details.

        Default example usage::

            from Bio.Emboss.Applications import WaterCommandline
            water_cmd = WaterCommandline(gapopen=10, gapextend=0.5,
                                         stdout=True, auto=True,
                                         asequence="a.fasta", bsequence="b.fasta")
            print("About to run: %s" % water_cmd)
            std_output, err_output = water_cmd()

        This functionality is similar to subprocess.check_output(). In general
        if you require more control over running the command, use subprocess
        directly.

        When the program called returns a non-zero error level, a custom
        ApplicationError exception is raised. This includes any stdout and
        stderr strings captured as attributes of the exception object, since
        they may be useful for diagnosing what went wrong.
        """
        if not stdout:
            stdout_arg = open(os.devnull, "w")
        elif isinstance(stdout, str):
            stdout_arg = open(stdout, "w")
        else:
            stdout_arg = subprocess.PIPE

        if not stderr:
            stderr_arg = open(os.devnull, "w")
        elif isinstance(stderr, str):
            if stdout == stderr:
                stderr_arg = stdout_arg  # Write both to the same file
            else:
                stderr_arg = open(stderr, "w")
        else:
            stderr_arg = subprocess.PIPE

        # We may not need to supply any piped input, but we setup the
        # standard input pipe anyway as a work around for a python
        # bug if this is called from a Windows GUI program.  For
        # details, see http://bugs.python.org/issue1124861
        #
        # Using universal newlines is important on Python 3, this
        # gives unicode handles rather than bytes handles.

        # Windows 7, 8, 8.1 and 10 want shell = True
        if sys.platform != "win32":
            use_shell = True
        else:
            win_ver = platform.win32_ver()[0]
            if win_ver in ["7", "8", "post2012Server", "10"]:
                use_shell = True
            else:
                use_shell = False
        child_process = subprocess.Popen(
            str(self),
            stdin=subprocess.PIPE,
            stdout=stdout_arg,
            stderr=stderr_arg,
            universal_newlines=True,
            cwd=cwd,
            env=env,
            shell=use_shell,
        )
        # Use .communicate as can get deadlocks with .wait(), see Bug 2804
        stdout_str, stderr_str = child_process.communicate(stdin)
        if not stdout:
            assert not stdout_str, stdout_str
        if not stderr:
            assert not stderr_str, stderr_str
        return_code = child_process.returncode

        # Particularly important to close handles on Jython and PyPy
        # (where garbage collection is less predictable) and on Windows
        # (where cannot delete files with an open handle):
        if not stdout or isinstance(stdout, str):
            # We opened /dev/null or a file
            stdout_arg.close()
        if not stderr or (isinstance(stderr, str) and stdout != stderr):
            # We opened /dev/null or a file
            stderr_arg.close()

        if return_code:
            raise ApplicationError(return_code, str(self), stdout_str, stderr_str)
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
     - names -- a list of string names (typically two entries) by which
       the parameter can be set via the legacy set_parameter method
       (eg ["-a", "--append", "append"]). The first name in list is used
       when building the command line. The last name in the list is a
       "human readable" name describing the option in one word. This
       must be a valid Python identifier as it is used as the property
       name and as a keyword argument, and should therefore follow PEP8
       naming.
     - description -- a description of the option. This is used as
       the property docstring.
     - filename -- True if this argument is a filename (or other argument
       that should be quoted) and should be automatically quoted if it
       contains spaces.
     - checker_function -- a reference to a function that will determine
       if a given value is valid for this parameter. This function can either
       raise an error when given a bad value, or return a [0, 1] decision on
       whether the value is correct.
     - equate -- should an equals sign be inserted if a value is used?
     - is_required -- a flag to indicate if the parameter must be set for
       the program to be run.
     - is_set -- if the parameter has been set
     - value -- the value of a parameter

    """

    def __init__(
        self,
        names,
        description,
        filename=False,
        checker_function=None,
        is_required=False,
        equate=True,
    ):
        self.names = names
        if not isinstance(description, str):
            raise TypeError("Should be a string: %r for %s" % (description, names[-1]))
        # Note 'filename' is for any string with spaces that needs quoting
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

    Attributes:
     - names -- a list of string names (typically two entries) by which
       the parameter can be set via the legacy set_parameter method
       (eg ["-a", "--append", "append"]). The first name in list is used
       when building the command line. The last name in the list is a
       "human readable" name describing the option in one word. This
       must be a valid Python identifier as it is used as the property
       name and as a keyword argument, and should therefore follow PEP8
       naming.
     - description -- a description of the option. This is used as
       the property docstring.
     - is_set -- if the parameter has been set

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

    The names argument should be a list containing one string.
    This must be a valid Python identifier as it is used as the
    property name and as a keyword argument, and should therefore
    follow PEP8 naming.
    """

    def __init__(
        self,
        names,
        description,
        filename=False,
        checker_function=None,
        is_required=False,
    ):
        # if len(names) != 1:
        #    raise ValueError("The names argument to _Argument should be a "
        #                     "single entry list with a PEP8 property name.")
        self.names = names
        if not isinstance(description, str):
            raise TypeError("Should be a string: %r for %s" % (description, names[-1]))
        # Note 'filename' is for any string with spaces that needs quoting
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


class _ArgumentList(_Argument):
    """Represent a variable list of arguments on a command line, e.g. multiple filenames."""

    # TODO - Option to require at least one value? e.g. min/max count?

    def __str__(self):
        if not isinstance(self.value, list):
            raise TypeError("Arguments should be a list")
        if not self.value:
            raise ValueError("Requires at least one filename")
        # A trailing space is required so that parameters following the last filename
        # do not appear merged.
        # e.g.:  samtools cat in1.bam in2.bam-o out.sam  [without trailing space][Incorrect]
        #        samtools cat in1.bam in2.bam -o out.sam  [with trailing space][Correct]
        if self.is_filename:
            return " ".join(_escape_filename(v) for v in self.value) + " "
        else:
            return " ".join(self.value) + " "


class _StaticArgument(_AbstractParameter):
    """Represent a static (read only) argument on a commandline.

    This is not intended to be exposed as a named argument or
    property of a command line wrapper object.
    """

    def __init__(self, value):
        self.names = []
        self.is_required = False
        self.is_set = True
        self.value = value

    def __str__(self):
        return "%s " % self.value


def _escape_filename(filename):
    """Escape filenames with spaces by adding quotes (PRIVATE).

    Note this will not add quotes if they are already included:

    >>> print((_escape_filename('example with spaces')))
    "example with spaces"
    >>> print((_escape_filename('"example with spaces"')))
    "example with spaces"
    >>> print((_escape_filename(1)))
    1

    Note the function is more generic than the name suggests, since it
    is used to add quotes around any string arguments containing spaces.
    """
    # Is adding the following helpful
    # if os.path.isfile(filename):
    #    # On Windows, if the file exists, we can ask for
    #    # its alternative short name (DOS style 8.3 format)
    #    # which has no spaces in it.  Note that this name
    #    # is not portable between machines, or even folder!
    #    try:
    #        import win32api
    #        short = win32api.GetShortPathName(filename)
    #        assert os.path.isfile(short)
    #        return short
    #    except ImportError:
    #        pass
    if not isinstance(filename, str):
        # for example the NCBI BLAST+ -outfmt argument can be an integer
        return filename
    if " " not in filename:
        return filename
    # We'll just quote it - works on Windows, Mac OS X etc
    if filename.startswith('"') and filename.endswith('"'):
        # Its already quoted
        return filename
    else:
        return '"%s"' % filename


def _test():
    """Run the Bio.Application module's doctests (PRIVATE)."""
    import doctest

    doctest.testmod(verbose=1)


if __name__ == "__main__":
    # Run the doctests
    _test()
