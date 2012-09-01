# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys


class BDFReader(object):
    def __init__(self, debug, log):
        self.relpath = True
        if sys.version_info < (2, 6):
            #version = sys.version_info
            self.relpath = False
            #raise RuntimeError("must use python 2.6 or greater...version=%s"
            #                   %(str(version)))

        if log is None:
            from pyNastran.general.logger import dummyLogger
            if debug:
                word = 'debug'
            else:
                word = 'info'
            loggerObj = dummyLogger()
            log = loggerObj.startLog(word)  # or info
        self.log = log

    def print_filename(self, filename):
        """
        Takes a path such as C:/work/fem.bdf and locates the file using
        relative paths.  If it's on another drive, the path is not modified.
        @param self the object pointer
        @param filename a filename string
        @retval filenameString a shortened representation of the filename
        """
        driveLetter = os.path.splitdrive(os.path.abspath(filename))[0]
        if driveLetter == os.path.splitdrive(os.curdir)[0] and self.relpath:
            return os.path.relpath(filename)
        else:
            return filename

    def open_file(self, infileName):
        """
        Takes a filename and opens the file.
        This method is used in order to support INCLUDE files.
        """
        #print self.isOpened
        if self.isOpened[infileName] == False:
            self._active_filenames.append(infileName)
            #self.log.info("*openFile bdf=|%s|  pwd=|%s|" %(infileName,
            #                                               os.getcwd()))
            if not os.path.exists(infileName):
                msg = "infileName=|%s| does not exist..." % (infileName)
                raise IOError(msg)
            infile = open(infileName, 'r')
            self.infilesPack.append(infile)
            self.lineNumbers.append(0)
            self.isOpened[infileName] = True
            self.linesPack.append([])
        else:
            pass
            #print "is already open...skipping"

    def get_file_stats(self):
        """
        gets information about the active BDF file being read
        @param self the object pointer
        @retval lineNumber the active file's line number
        """
        filename = self._active_filenames[-1]
        return (filename, self.get_line_number())

    def get_line_number(self):
        """
        Gets the line number of the active BDF (used for debugging).
        @param self the object pointer
        @retval returns the line number of the active BDF filename
        """
        lineNumber = self.lineNumbers[-1]
        return lineNumber

    def get_next_line(self, debug=False):
        """
        Gets the next line in the BDF
        @param self
          the BDF object
        @param debug
          developer debug
        @retval line
          the next line in the BDF or None if it's the end of the current file
        """
        self.lineNumbers[-1] += 1
        linesPack = self._make_lines_pack(debug=False)
        #print "len(linesPack) = ", len(linesPack)
        #for line in linesPack:
            #print("$  |%r|" %(line))

        if len(linesPack) == 0:
            self.close_file()
            return None
            #linesPack = self._make_lines_pack(debug=debug)
            #return lastLine
        #print linesPack[0]
        return linesPack.pop(0)

    def close_file(self, debug=False):
        """
        Closes the active file object.
        If no files are open, the function is skipped.
        This method is used in order to support INCLUDE files.
        @param self the object pointer
        @param debug developer debug
        """
        if len(self.infilesPack) == 0:
            return
        if debug:
            self.log.debug("*closing")
        infile = self.infilesPack.pop()
        infile.close()

        #if debug:
        #    print [os.path.relpath(fname) for fname in self._active_filenames]
        lineNumbers = self.lineNumbers.pop()
        active_filename = self._active_filenames.pop()
        linesPack = self.linesPack.pop()
        self.isOpened[active_filename] = False

        if len(self.linesPack) == 0:
            raise IOError('\nThe bdf closed unexpectedly...\n  an Executive '
                          'and Case Control Decks are required...'
                          'put a CEND and BEGIN BULK in the BDF')
        nlines = len(self.linesPack[-1])

        ## determines if self.activefilename should be closed at the next
        ## opportunity
        self.doneReading = False
        if debug:
            fnameA = self.print_filename(active_filename)
            fnameB = self.print_filename(self.bdf_filename)

            self.log.debug("active_filename=|%s| infilename=%s len(pack)=%s\n"
                           % (fnameA, fnameB, nlines))
        #print "\n\n"

    def _set_infile(self, bdf_filename, includeDir=None):
        """
        Sets up the basic file/lines/cardCounting operations
        @param self
          the BDF object
        @param bdf_filename
          the input BDF filename
        @param includeDir
          the location of include files if an absolute/relative path is
          not used (not supported in Nastran)
        """
        ## automatically rejects every parsable card (default=False)
        self._auto_reject = False
        ## is the active file done reading
        self.doneReading = False
        ## was an ENDDATA card found
        self.foundEndData = False

        if includeDir is None:
            includeDir = os.path.dirname(bdf_filename)
        ## the active filename (string)
        self.bdf_filename = bdf_filename
        ## the directory of the 1st BDF (include BDFs are relative to this one)
        self.includeDir = includeDir
        ## list of infile objects (needed for INCLUDE files)
        self.infilesPack = []
        ## list of lines from self.activeFilename that are stored
        self.linesPack = []
        ## list of all open filenames
        self._active_filenames = []
        ## stores the line number of self.activefilename that the parser is on
        ## very helpful when debugging
        self.lineNumbers = []
        ## dictionary that says whether self.bdf_filename is open/close
        ## (boolean)
        self.isOpened = {self.bdf_filename: False}
        ## list of all read in cards - useful in determining if
        ## entire BDF was read & really useful in debugging
        self.cardCount = {}
        ## stores the cardCount of cards that have been rejected
        self.rejectCount = {}
