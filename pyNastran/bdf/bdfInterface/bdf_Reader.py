## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
import os
import sys

from pyNastran.bdf.errors import *

class bdfReader(object):
    def __init__(self,debug,log):
        self.relpath = True
        if sys.version_info < (2,6):
            version = sys.version_info
            self.relpath = False
            #raise Exception("must use python 2.6 or greater...version=%s" %(str(version)))

        if log is None:
            from pyNastran.general.logger import dummyLogger
            loggerObj = dummyLogger()
            if debug:
                word = 'debug'
            else:
                word = 'info'
            log = loggerObj.startLog(word) # or info
        self.log = log

    def printFileName(self,filename):
        """
        Takes a path such as C:/work/fem.bdf and locates the file using relative paths
        If it's on another drive, the path is not modified.
        @param self the object pointer
        @param filename a filename string
        @retval filenameString a shortened representation of the filename
        """
        driveLetter = os.path.splitdrive(filename)[0]
        if driveLetter==os.path.splitdrive(os.curdir)[0] and self.relpath:
            return os.path.relpath(filename)
        else:
            return filename
        ###

    def openFile(self,infileName):
        """
        Takes a filename and opens the file.  
        This method is used in order to support INCLUDE files.
        """
        #print self.isOpened
        if self.isOpened[infileName]==False:
            self.activeFileNames.append(infileName)
            #self.log.info("*openFile bdf=|%s|  pwd=|%s|" %(infileName,os.getcwd()))
            if not os.path.exists(infileName):
                raise MissingFileError("infileName=|%s| does not exist..." %(infileName))
            infile = open(infileName,'r')
            self.infilesPack.append(infile)
            self.lineNumbers.append(0)
            self.isOpened[infileName]=True
            self.linesPack.append([])
        ###
        else:
            pass
            #print "is already open...skipping"
        ###

    def getFileStats(self):
        """
        gets information about the active BDF file being read
        @param self the object pointer
        @retval lineNumber the active file's line number
        """
        filename   = self.activeFileNames[-1]
        return (filename,self.getLineNumber())

    def getLineNumber(self):
        """
        Gets the line number of the active BDF (used for debugging).
        @param self the object pointer
        @retval returns the line number of the active BDF filename
        """
        lineNumber = self.lineNumbers[-1]
        return lineNumber

    def getNextLine(self,debug=False):
        """
        Gets the next line in the BDF
        @param self the object pointer
        @param debug developer debug
        @retval line the next line in the BDF or None if it's the end of a the current file
        """
        self.lineNumbers[-1]+=1
        linesPack = self.makeLinesPack(debug=False)
        #print "len(linesPack) = ",len(linesPack)
        #for line in linesPack:
            #print("$  |%r|" %(line))

        if len(linesPack)==0:
            self.closeFile()
            return None
            #linesPack = self.makeLinesPack(debug=debug)
            #return lastLine
        #print linesPack[0]
        return linesPack.pop(0)

    #def getNextLine2(self):
    #    self.lineNumbers[-1]+=1
    #    return self.infilesPack[-1].readline()

    #    infile = self.infilesPack[-1]
    #    #print "infile = |%s|" %(infile),type(infile)
    #    line = infile.readline()
    #    print "line = |%r|" %(line)
    #    return line

    def closeFile(self,debug=False):
        """
        Closes the active file object.
        If no files are open, the function is skipped.
        This method is used in order to support INCLUDE files.
        @param self the object pointer
        @param debug developer debug
        """
        if len(self.infilesPack)==0:
            return
        if debug:
            self.log.debug("*closing")
        infile = self.infilesPack.pop()
        infile.close()

        #if debug:
        #    print [os.path.relpath(fname) for fname in self.activeFileNames]
        lineNumbers = self.lineNumbers.pop()
        activeFileName = self.activeFileNames.pop()
        linesPack = self.linesPack.pop()
        self.isOpened[activeFileName] = False
        
        if len(self.linesPack)==0:
            raise ClosedBDFError('\nThe bdf closed unexpectedly...\n  an Executive and Case Control Decks are required...put a CEND and BEGIN BULK in the BDF')
        nlines = len(self.linesPack[-1])
        ## determines if self.activefilename should be closed at the next opportunity
        self.doneReading = False
        if debug:
            fnameA = self.printFileName(activeFileName)
            fnameB = self.printFileName(self.bdfFileName)

            self.log.debug("activeFileName=|%s| infilename=%s len(pack)=%s\n" %(fnameA,fnameB,nlines))
        ###
        #print "\n\n"

    def _setInfile(self,bdfFileName,includeDir=None):
        """
        sets up the basic file/lines/cardCounting operations
        @param self the object pointer
        @param bdfFileName  the input BDF filename
        @param includeDir   the location of include files if a absolute/relative path is not used (not supported in Nastran)
        """
        ## automatically rejects every parsable card (default=False)
        self.autoReject   = False
        ## is the active file done reading
        self.doneReading  = False
        ## was an ENDDATA card found
        self.foundEndData = False

        if includeDir is None:
            includeDir = os.path.dirname(bdfFileName)
        ## the active filename (string)
        self.bdfFileName = bdfFileName
        ## the directory of the 1st BDF (include BDFs are relative to this one)
        self.includeDir = includeDir
        ## list of infile objects (needed for INCLUDE files)
        self.infilesPack     = []
        ## list of lines from self.activeFilename that are stored
        self.linesPack       = []
        ## list of all open filenames
        self.activeFileNames = []
        ## stores the line number of self.activefilename that the parser is on
        ## very helpful when debugging
        self.lineNumbers     = []
        ## dictionary that says whether self.bdfFileName is open/close (boolean0
        self.isOpened = {self.bdfFileName: False}
        ## list of all read in cards - useful in determining if
        ## entire BDF was read & really useful in debugging
        self.cardCount = {}
        ## stores the cardCount of cards that have been rejected
        self.rejectCount = {}

