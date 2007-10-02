"""
This implements a set of classes that wraps a file object interface
around code that executes in another process.  This allows you fork
many different commands and let the run concurrently.

Functions:
copen_sys     Open a file-like pipe to a system command.
copen_fn      Open a file-like pipe to a python function.

"""

import warnings
warnings.warn("Bio.MultiProc is deprecated. If you want to use this code, please let the Biopython developers know by sending an email to biopython-dev@biopython.org to avoid permanent removal of Bio.MultiProc.",
              DeprecationWarning)

import os
import sys
import time
import signal

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
        os.close(r)
        os.close(er)
        os.dup2(w, sys.stdout.fileno())
        os.dup2(ew, sys.stderr.fileno())
        try:
            os.execvp(syscmd, args)  # execute it!
        except:
            sys.stderr.write("%s could not be executed\n" % syscmd)
            os._exit(-1)
        os._exit(0)

    # parent
    os.close(w)
    os.close(ew)
    return _ProcHandle(pid, os.fdopen(r, 'r'), os.fdopen(er, 'r'))

def copen_fn(func, *args, **keywords):
    """copen_fn(func, *args, **keywords) -> file-like object

    Open a file-like object that returns the output from function
    call.  The object's 'read' method returns the return value from
    the function.  The function is executed as a separate process so
    any variables modified by the function does not affect the ones in
    the parent process.  The return value of the function must be
    pickle-able.

    """
    try:
        import cPickle as pickle
    except ImportError:
        import pickle
    r, w = os.pipe()
    er, ew = os.pipe()

    pid = os.fork()
    if pid == 0: # child process
        cwrite, errwrite = os.fdopen(w, 'w'), os.fdopen(ew, 'w')
        try:
            output = func(*args, **keywords)
            # Pickle may fail here is the object is not pickleable.
            s = pickle.dumps(output, 1)
        except:
            import traceback
            traceback.print_exc(file=errwrite)
            errwrite.flush()
            os._exit(-1)
        try:
            cwrite.write(s)
            cwrite.flush()
        except IOError, x:
            # There can be an IOError if the parent is no longer
            # listening.  Ignore it.
            pass
        os._exit(0)

    # parent
    os.close(w)
    os.close(ew)
    return _PickleHandle(pid, os.fdopen(r, 'r'), os.fdopen(er, 'r'))


class _ProcHandle:
    """This object provides a file-like interface to a running
    process.

    Members:
    pid         what is the PID of the subprocess?
    killsig     what signal killed the child process?
    status      what was the status of the command?

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
        """Create a wrapper around a running process.  pid is the
        process ID.  cread is the file object used to read from the
        child.  errread is an optional file object used to read errors
        from the child.

        """
        _active.append(self)
        
        self.pid = pid
        self.status = None
        self.killsig = None
        
        self._start, self._end = time.time(), None
        self._cread, self._errread = cread, errread
        self._output = []
        self._done = 0
        self._closed = 0

    def __del__(self):
        self.close()  # kill the process

    def _kill(self):
        """Kill the process and return killsig"""
        try:
            pid, ind = os.waitpid(self.pid, os.WNOHANG)
            if pid == self.pid:   # died
                return 0
            # First, try to kill it with a SIGTERM.
            os.kill(self.pid, signal.SIGTERM)
            # Wait .5 seconds for it to die.
            end = time.time() + 0.5
            while time.time() < end:
                pid, ind = os.waitpid(self.pid, os.WNOHANG)
                if pid == self.pid:
                    return ind & 0xff
                time.sleep(0.1)
            # It didn't die, so kill with a SIGKILL
            os.kill(self.pid, signal.SIGKILL)
            return signal.SIGKILL
        except OSError:
            pass
        
    def close(self):
        """Close the process, killing it if it is still running."""
        # If this gets called in the middle of object initialization,
        # the _closed attribute will not exist.
        if not hasattr(self, '_closed') or self._closed:
            return
        # on cleanup, _active may not be defined!
        if _active and self in _active:
            _active.remove(self)
        if not self._done:
            self.killsig = self._kill()
            self._end = time.time()
            self.status = None
            self.killsig = signal.SIGTERM
            self._done = 1
        self._output = []
        self._closed = 1

    def fileno(self):
        """Return the file descriptor used to read from the process."""
        return self._cread.fileno()
        
    def readline(self):
        """Return the next line or '' if finished."""
        self.wait()
        if not self._output:
            return ''
        line = self._output[0]
        del self._output[0]
        return line

    def readlines(self):
        """Return the output of the process as a list of strings."""
        self.wait()
        output = self._output
        self._output = []
        return output

    def read(self):
        """Return the output as a string."""
        self.wait()
        output = self._output
        self._output = []
        return "".join(output)

    def wait(self):
        """Wait for the process to finish."""
        import select
        if self._done:
            return
        # wait until stuff's ready to be read
        select.select([self], [], [])
        self._cleanup_child()

    def poll(self):
        """Return a boolean.  Is the process finished running?"""
        import select
        if self._done:
            return 1
        # If I'm done, then read the results.
        if select.select([self], [], [], 0)[0]:
            self._cleanup_child()
        return self._done

    def elapsed(self):
        """Return the number of seconds elapsed since the process began."""
        if self._end:  # if I've finished, return the total time
            return self._end - self._start
        return time.time() - self._start

    def _cleanup_child(self):
        """Do necessary cleanup functions after child is finished running."""
        if self._done:
            return

        # read the output
        self._output = self._cread.readlines()
        self._cread.close()
        if self._errread:
            error = self._errread.read()
            self._errread.close()
            if error:
                raise AssertionError, "Error in child process:\n\n%s" % error
                # It would be nice to be able to save the exception
                # and traceback somehow, and raise it in the parent.
                #raise etype, value
        # Remove myself from the active list.
        if _active and self in _active:
            _active.remove(self)

        pid, ind = os.waitpid(self.pid, 0)
        self.status, self.killsig = ind >> 8, ind & 0xff
        self._end = time.time()
        self._done = 1

class _PickleHandle:
    """

    Members:
    pid         what is the PID of the subprocess?
    killsig     what signal killed the child process?
    status      what was the status of the command?

    Methods:
    close       Close this process, killing it if necessary.
    fileno      Return the fileno used to read from the process.
    wait        Wait for the process to finish.
    poll        Is the process finished?
    elapsed     How much time has this process taken?
    read        Return a Python object.
    
    """
    def __init__(self, pid, cread, errread=None):
        """Create a wrapper around a running process.  pid is the
        process ID.  cread is the file object used to read from the
        child.  errread is an optional file object used to read errors
        from the child.

        """
        self._phandle = _ProcHandle(pid, cread, errread)

    def __getattr__(self, attr):
        # This object does not support 'readline' or 'readlines'
        if attr.startswith('readline'):
            raise AttributeError, attr
        return getattr(self._phandle, attr)

    def read(self):
        """Return a Python object or ''."""
        try:
            import cPickle as pickle
        except ImportError:
            import pickle
        r = self._phandle.read()
        if not r:
            return r
        return pickle.loads(r)


# Handle SIGTERM below

# Keep a list of all the active child processes.  If the process is
# forcibly killed, e.g. by a SIGTERM, make sure the child processes
# die too.
_active = []   # list of _ProcHandle objects

_HANDLING = 0
def _handle_sigterm(signum, stackframe):
    """Handles a SIGTERM.  Cleans up."""
    global _HANDLING
    if _HANDLING:
        return
    _HANDLING = 1
    _cleanup()
    # call the previous handler
    if _PREV_SIGTERM is not None:
        signal.signal(signal.SIGTERM, _PREV_SIGTERM)
    os.kill(os.getpid(), signum)

def _cleanup():
    """Close all active commands."""
    for obj in _active[:]:
        obj.close()

_PREV_SIGTERM = signal.getsignal(signal.SIGTERM)
signal.signal(signal.SIGTERM, _handle_sigterm)
