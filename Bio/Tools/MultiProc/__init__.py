# Copyright 2001 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

__all__ = [
    'Scheduler',
    'Task',
    'copen'
    ]

import Scheduler
import Task
import time

def run(nprocs, fn, fn_args=(), fn_keywds={},
        sleep=0.1, always_use_scheduler=0):
    """run(nprocs, fn[, fn_args][, fn_keywds][, sleep])

    Run nprocs copies of fn concurrently.  The first two parameters of
    fn should be a process number (0 to nprocs-1) and the total number
    of processes.  Additional arguments and keywords can be passed in
    the fn_args and fn_keywds parameters.  Sleep is the amount of time
    the scheduler should wait between process polls.

    """
    if nprocs < 1 or nprocs > 100:
        raise ValueError, "nprocs %d out of range" % nprocs
    # Undocumented: return a list of return values.  The return values
    # of the functions must be pickle-able types.  This is
    # undocumented and may go away in future releases.  I'm not sure
    # it's a good idea right now...
    retvals = [None] * nprocs
    
    # If they only want to run 1 process, then just run it without
    # going through the scheduler.
    if nprocs == 1 and not always_use_scheduler:
        args = (0, 1) + fn_args
        r = fn(*args, **fn_keywds)
        retvals[0] = r
    else:
        def save_retval(task, retvals=retvals):
            i = int(task.getName())
            retvals[i] = task.retval
        scheduler = Scheduler.Scheduler(nprocs, finish_fn=save_retval)
        for i in range(nprocs):
            args = (i, nprocs) + fn_args
            task = Task.Task(name=i, target=fn, args=args, kwargs=fn_keywds)
            scheduler.add(task)
        while scheduler.run():
            time.sleep(sleep)
    return retvals
