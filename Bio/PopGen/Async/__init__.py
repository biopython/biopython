# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


'''
Support for asynchronous execution.

'''

import os
import threading

__docformat__ = "restructuredtext en"

class Async(object):
    '''Abstract Asynchronous execution class.

       This is the top abstract class.
       Concrete classes must implement the _run_program method.
    '''

    def __init__(self):
        '''Async constructor.

       Initializes the queues, among other things.
       Of notice, is the access_ds lock for controlling exclusive
       access to this object.
        '''
        self.running = {}
        self.waiting = []
        self.done = {}
        self.id = 0
        self.hooks = {}
        self.access_ds = threading.Lock()

    def run_program(self, program, parameters, input_files):
        '''Runs a program.

           Real _run_program to be implemented by concrete classes.

           parameters:
           program String identifying program.
           parameters List of String parameters.
           input_files Hash of Input file descriptors.

           returns:
           Task Id.

           The input_files hash key is the path that is passed
           to the program. It should always be relative.
           Value is a stream.
        '''
        if program in self.hooks:
            self.access_ds.acquire()
            self.id += 1
            id = self.id
            self.access_ds.release()
            self._run_program(id, self.hooks[program], parameters, input_files)
            return id

    def _run_program(self, id, program, parameters, input_files):
        """Actually run the program, handled by a subclass (PRIVATE).

        This method should be replaced by any derived class to do
        something useful. It will be called by the run_program method.
        """
        raise NotImplementedError("This object should be subclassed")

    def get_result(self, id):
        ''' Returns the results for a certain Id, the info for that Id is
            forgotten.

            parameters:
            id Id of the task.

            returns:
            (return_code, output_files) return code and file access
            object.

            The output_files hash key is a relative file name, and the value a
            output stream.
        '''
        self.access_ds.acquire()
        if id in self.done:
            returnCode, fileObject = self.done[id]
            del self.done[id]
            self.access_ds.release()
        else:
            self.access_ds.release()
            return None


class FileRetriever(object):
    '''An Abstract Support class to retrieve files.
    '''

    def __init__(self):
        self.file_list = []

    def get_File_list(self):
        '''Returns the list of available files.
        '''
        return self.file_list

    def get_file(self, name):
        raise NotImplementedError('Abstract method')


class DirectoryRetriever(FileRetriever):
    '''Retrieves a directory content.
    '''

    def __init__(self, directory):
        FileRetriever.__init__(self)
        self.directory = directory
        walk_list = os.walk(directory)
        for dir, dir_list, file_list in walk_list:
            for file in file_list:
                self.file_list.append(file[len(directory) + 1:])

    def get_file(self, name):
        return open(self.directory + os.sep + name)
