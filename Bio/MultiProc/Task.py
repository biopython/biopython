"""Task.py

Functionality for running multiple tasks simultaneously.

Classes:
Task    Processing that can be forked off in a separate process.

"""
import copen

# Copied from threading.py
# This is not thread safe!
_counter = 0
def _newname(template="Task-%d"):
    """_newname(template="Task-%d") -> name"""
    global _counter
    _counter = _counter + 1
    return template % _counter

class Task:
    """Contains information for one process.

    Implements part of the Thread interface.

    Methods:
    start      Start this task.  Should be called once.
    run        Called by start to really run the task.
    getName    Get the name of the task.
    setName    Set the name of the task.
    isAlive    Whether this Task is still running.

    Members:
    retval     Return value of the function.

    """
    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}):
        """Task([group][, target][, name][, args][, kwargs])

        Create a task object.  group should be None and is reserved
        for future expansion.  target is the function to be called.
        name is the name of the thread.  args and kwargs are the
        arguments to be passed to target.

        """
        self._target = target
        if name is None:
            name = _newname()
        self._name = str(name)
        self._args, self._kwargs = args, kwargs
        self._start = self._finished = None
        self._handle = None
        self.retval = None

    def __del__(self):
        # Close my handle if it's not finished running.
        if self._handle:
            self._handle.close()

    def start(self):
        """S.start()

        Start this task.  Should only be called once.

        """
        if self._start is not None:
            raise ValueError, "task %s already started" % self._name
        self._start = 1
        self.run()

    def run(self):
        """S.run()

        Run this task.  Should only be called by S.start().

        """
        # If the client didn't specify a target function, then don't
        # do any processing.
        if not self._target:
            self._finished = 1
        else:
            self._handle = copen.copen_fn(
                self._target, *self._args, **self._kwargs)

    def getName(self):
        """S.getName() -> name"""
        return self._name

    def setName(self, name):
        """S.setName(name)"""
        self._name = name

    def isAlive(self):
        """S.isAlive() -> boolean"""
        if not self._finished:
            if self._handle and self._handle.poll():
                self.retval = self._handle.read()
                self._finished = 1
        return not self._finished
