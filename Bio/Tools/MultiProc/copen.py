"""copen.py

This implements a set of command classes that behave like
file objects.  This allows you fork off many different commands
and let the run concurrently.

Functions:
copen_sys     Open a file-like pipe to a system command.
copen_fn      Open a file-like pipe to a python function.

"""
import os
import sys
import traceback
import time
import signal
import select
try:
    import cPickle as pickle
except ImportError:
    import pickle

def copen_sys(syscmd, *args):
    """copen_sys(syscmd, *args) -> file-like object

    Open a file-like object that returns the output from a system
    command.

    """
    # python requires first element to be the path
    if not args or args[0] != syscmd:
        args = [syscmd] + list(args)

    r, w = os.pipe()
    er, ew = os.pipe()

    pid = os.fork()
    if pid == 0: # child process
        os.dup2(w, 1)
        os.dup2(ew, 2)
        try:
            os.execvp(syscmd, args)  # execute it!
        except:
            sys.stderr.write("%s could not be executed\n" % path)
            os._exit(-1)
        os._exit(0)

    # parent
    os.close(w)
    os.close(ew)
    return _CommandHandle(pid, os.fdopen(r, 'r'), os.fdopen(er, 'r'))

def copen_fn(func, *args, **keywords):
    """copen_fn(func, *args, **keywords) -> file-like object

    Open a file-like object that returns the output from function
    call.  The object's 'read' method returns the return value from
    the function.  The function is executed as a separate process, so
    any variables modified by the function does not affect the ones in
    the parent process.  The return value of the function must be
    pickle-able.

    """
    r, w = os.pipe()
    er, ew = os.pipe()

    pid = os.fork()
    if pid == 0: # child process
        cwrite, errwrite = os.fdopen(w, 'w'), os.fdopen(ew, 'w')
        try:
            output = apply(func, args, keywords)
            s = pickle.dumps(output, 1)
            cwrite.write(s)
            cwrite.flush()
        except pickle.PicklingError, x:
            errwrite.write(x)
            errwrite.flush()
            os._exit(-1)
        except:
            type, value, tb = sys.exc_info()   # get the traceback
            errwrite.writelines(traceback.format_tb(tb))
            del tb       # delete it, so no circular reference
            if value:
                errwrite.write("%s: %s" % (type, value))
            else:
                errwrite.write(type)
            errwrite.flush()
            os._exit(-1)
        os._exit(0)

    # parent
    os.close(w)
    os.close(ew)
    return _PickleHandle(pid, os.fdopen(r, 'r'), os.fdopen(er, 'r'))


# Keep a list of all the active child processes.  If the process is
# forcibly killed, e.g. by a SIGTERM, make sure the child processes
# die too.
_active = []   # list of _CommandHandle objects

class _CommandHandle:
    """This file-like object is a wrapper around a command.

    Members:
    pid         what is the PID of the subprocess?
    killsig     what signal killed the child process?
    status      what was the status of the command?
    error       if an error occurred, this describes it.

    Methods:
    close       Close this process, killing it if necessary.
    fileno      Return the fileno used to read from the process.
    wait        Wait for the process to finish.
    poll        Is the process finished?
    elapsed     How much time has this process taken?
    read
    readline
    readlines
    
    """
    def __init__(self, pid, cread, errread=None):
        """_CommandHandle(pid, cread[, errread]) -> instance

        Create a wrapper around a command.  pid should be the process
        ID of the command that was created, probably by a fork/exec.
        cread should be a file object used to read from the child.  If
        errread is given, then I will look there for messages
        pertaining to error conditions.

        """
        _active.append(self)
        
        self.pid = pid
        self.status = None
        self.killsig = None
        self.error = ""
        
        self._start, self._end = time.time(), None
        self._cread, self._errread = cread, errread
        self._output = []
        self._done = 0
        self._closed = 0

    def __del__(self):
        self.close()  # kill the process

    def close(self):
        """S.close()

        Close the process, killing it if I must.

        """
        # If this gets called in the middle of object initialization,
        # the _closed attribute will not exist.
        if not hasattr(self, '_closed') or self._closed:
            return
        # on cleanup, _active may not be defined!
        if _active and self in _active:
            _active.remove(self)
        if not self._done:
            try:
                pid, ind = os.waitpid(self.pid, os.WNOHANG)
                if pid != self.pid:    # not finished
                    os.kill(self.pid, signal.SIGTERM)
            except OSError:
                pass
            self._end = time.time()
            self.status = None
            self.killsig = signal.SIGTERM
            self._done = 1
        self._output = []
        self._closed = 1

    def fileno(self):
        """S.fileno() -> file descriptor

        Return the file descriptor associated with the pipe.

        """
        return self._cread.fileno()
        
    def readline(self):
        """S.readline() -> string

        Return the next line read, or '' if finished.

        """
        self.wait()
        if not self._output:
            return ''
        line = self._output[0]
        del self._output[0]
        return line

    def readlines(self):
        """S.readlines() -> list of strings

        Return the output as a list of strings.

        """
        self.wait()
        output = self._output
        self._output = []
        return output

    def read(self):
        """S.read() -> string

        Return the output as a string.

        """
        self.wait()
        output = self._output
        self._output = []
        return "".join(output)

    def wait(self):
        """S.wait()

        Wait until the process is finished.

        """
        if self._done:
            return
        # wait until stuff's ready to be read
        select.select([self], [], [])
        self._cleanup_child()

    def poll(self):
        """S.poll() -> boolean

        Is the process finished running?

        """
        if self._done:
            return 1
        # If I'm done, then read the results.
        if select.select([self], [], [], 0)[0]:
            self._cleanup_child()
        return input

    def elapsed(self):
        """S.elapsed() -> num seconds

        How much time has elapsed since the process began?

        """
        if self._end:  # if I've finished, return the total time
            return self._end - self._start
        return time.time() - self._start

    def _cleanup_child(self):
        """S._cleanup_child()

        Do necessary cleanup functions after child is finished running.

        """
        if self._done:
            return

        # read the output
        self._output = self._cread.readlines()
        self._cread.close()
        if self._errread:
            self.error = self._errread.read()
            self._errread.close()
        # Remove myself from the active list.
        if _active and self in _active:
            _active.remove(self)

        pid, ind = os.waitpid(self.pid, 0)
        self.status, self.killsig = ind >> 8, ind & 0xff
        self._end = time.time()
        self._done = 1

class _PickleHandle:
    """This is a decorator around a _CommandHandle.
    Instead of returning the results as a string, it returns a
    python object.  The child process must pickle its output to the pipe!

    Members:
    pid         what is the PID of the subprocess?
    killsig     what signal killed the child process?
    status      what was the status of the command?
    error       if an error occurred, this describes it.

    Methods:
    close       Close this process, killing it if necessary.
    fileno      Return the fileno used to read from the process.
    wait        Wait for the process to finish.
    poll        Is the process finished?
    elapsed     How much time has this process taken?
    read        Return a Python object.
    
    """
    def __init__(self, pid, cread, errread=None):
        """_PickleHandle(pid, cread[, errread])

        Create a wrapper around a command.  pid should be the process ID
        of the command that was created, probably by a fork/exec.
        cread should be a file object used to read from the child.
        If errread is given, then I will look there for messages pertaining to
        error conditions.

        """
        self._cmd_handle = _CommandHandle(pid, cread, errread)

    def __getattr__(self, attr):
        # This object does not support 'readline' or 'readlines'
        if attr.startswith('readline'):
            raise AttributeError, attr
        return getattr(self._cmd_handle, attr)

    def read(self):
        """S.read() -> python object, or None

        Returns None on error.  Most likely, the function returned
        an object that could not be pickled.

        """
        r = self._cmd_handle.read()
        if r:
            obj = pickle.loads(r)
        else:
            raise self.error
        return obj


# Handle SIGTERM below

def _handle_sigterm(signum, stackframe):
    """Handles a SIGTERM.  Cleans up."""
    _cleanup()
    # call the previous handler
    if _PREV_SIGTERM is not None:
        signal.signal(signal.SIGTERM, _PREV_SIGTERM)

def _cleanup():
    """_cleanup()

    Close all active commands.

    """
    for obj in _active[:]:
        obj.close()

_PREV_SIGTERM = signal.getsignal(signal.SIGTERM)
signal.signal(signal.SIGTERM, _handle_sigterm)
