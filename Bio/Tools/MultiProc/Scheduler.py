"""Scheduler.py

Classes:
Scheduler   Schedules threads to be run.

"""
class Scheduler:
    """Schedules threads to be run.  No prioritization.  Nothing fancy.

    Methods:
    add           Add a thread to be run.
    num_left      Return the number of threads left.
    num_running   Return the number of threads currently running.
    run           Main loop.  Returns whether there's still threads left.
    
    """
    def __init__(self, max_threads, start_fn=None, finish_fn=None):
        """Scheduler(max_threads[, start_fn][, finish_fn]) -> object

        max_threads is the maximum number of threads to run at a time.
        start_fn and finish_fn are optional callbacks that take a
        thread as an argument.  They are called before and after each
        thread.

        """
        if max_threads <= 0:
            raise ValueError, "You must specify a positive number of threads!"
        self._max_threads = max_threads  # maximum active threads at a time
        self._start_fn = start_fn        # called before a thread starts
        self._finish_fn = finish_fn      # called when a thread is finished
        self._waiting = []               # list of threads waiting to run
        self._in_progress = []           # list of threads running

    def run(self):
        """S.run() -> boolean

        Execute the main loop.  Return a boolean indicating whether
        threads are still running.

        """
        # See if the running threads are finished.  Remove any that
        # are.
        i=0
        while i < len(self._in_progress):
            if not self._in_progress[i].isAlive():
                t = self._in_progress.pop(i)
                if self._finish_fn:
                    self._finish_fn(t)
            else:
                i = i + 1

        # If I have any more slots, start new threads.
        while self._waiting and len(self._in_progress) < self._max_threads:
            t = self._waiting.pop(0)
            if self._start_fn:
                self._start_fn(t)
            t.start()
            self._in_progress.append(t)
            
        return self._waiting or self._in_progress

    def add(self, thread):
        """S.add(thread)"""
        self._waiting.append(thread)

    def num_left(self):
        """S.num_left() -> number of threads left to run"""
        return len(self._waiting)

    def num_running(self):
        """S.num_running() -> number of threads currently running"""
        return len(self._in_progress)
