#!/usr/bin/env python
# Created: Sat Dec  2 16:02:17 2000
# Last changed: Time-stamp: <00/12/03 13:12:39 thomas>
# Thomas.Sicheritz@molbio.uu.se, http://evolution.bmc.uu.se/~thomas
# File: xbb_blastbg.py

import posixpath, posix
import os, sys, commands
sys.path.insert(0, '.')
import Queue, threading
import tempfile
from Tkinter import *
from xbb_utils import NotePad

class BlastDisplayer:
    def __init__(self, command, text_id = None):
        self.command = command
        self.tid = text_id
        
    def RunCommand(self):
        self.outfile = tempfile.mktemp()

        # make sure outfile exists and is empty
        fid = open(self.outfile,'w+'); fid.close()
        
        com = '%s > %s' % (self.command, self.outfile)

        self.worker = BlastWorker(com)
        self.worker.start()
        self.UpdateResults()
        
    def UpdateResults(self):
        # open the oufile and displays new appended text
        fid = open(self.outfile)
        size = 0
        while 1:
            if self.worker.finished: break
            fid.seek(size)
            txt = fid.read()
            size = os.stat(self.outfile)[6]
            try:
                self.tid.insert(END, txt)
                self.tid.update()
            except:
                # text widget is detroyed, we assume the search
                # has been cancelled
                break

        fid.close()
        self.Exit()            
        

    def Exit(self):
        if os.path.exists(self.outfile):
            os.remove(self.outfile)

        # do I need to stop the queue ?
        self.worker.shutdown()
        del self.worker
        
        
class BlastWorker(threading.Thread):
  
    def __init__(self, command):
        self.com = command
        queue = Queue.Queue(0)
        self.queue = queue
        threading.Thread.__init__(self)
        self.finished = 0
        print dir(queue)
        print queue.queue

    def shutdown(self):
        # GRRRR How do I explicitely kill a thread ???????
        #self.queue.put(None)
        del self.queue
        
    def run(self):
        print 'running', self.com
        os.system(self.com)
        self.finished = 1


if __name__ == '__main__':
    np = NotePad()
    com = "blastall -d nr -i test.fas -p blastx"
    test = BlastDisplayer(com, np.tid)
    test.RunCommand()
