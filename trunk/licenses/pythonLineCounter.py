#!/usr/bin/python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# Name:       pytonLineCounter.py (originally plc.py)
# Purpose:    A Python line counter.  Reports how many lines are in the given 
#             input files, broken down into code, comment and blank lines.
#
# Author:     Wayne Koorts
# Created:    31/03/2009
# Copyright:  Copyright 2009 Wayne Koorts
# Licence:    This software is released under the terms of the
#             GNU General Public License.  Please see license.txt for details.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Per discussion with Wayne, he has agreed to release this under the 
#    LGPL license for pyNastran.  Updates for additional breakdowns were
#    added.  Wayne's script was also renamed for clarity.  It runs on a single
#    file, but creating a single file with all the pyNastran code allows plc to
#    run on the entire codebase.
#

import string
import sys

usage = \
"""Usage:
    Enter a list of files to count, one per line, in a file named
    inputs.txt in this script's path, or pass in a list of
    space-separated filenames to count, e.g:
        
        plc.py <file1.py> <file2.py> <file3.py>"""

def countLines(fnames=[]):
    """
    def f(x):                                        -> code
        '''                                          -> header line
        comment                                      -> header line
        '''                                          -> header line
        return y                                     -> code
                                                     -> blank
    def g(x):                                        -> code
        '''comment'''                                -> header line
        return y                                     -> code
                                                     -> blank
    print f(5)                                       -> code
    if fnames==[]:             # 1st line            -> code + comment line   -> code
        fnames = sys.argv[1:]  ## @todo finish code  -> code + doxygen marker -> code
    ###                                              -> formatting line
                                                     -> blank
    ## standard doxygen                              -> doxygen line
    comment_lines = 0          ## standard doxygen   -> code + doxygen line -> code
    blank_lines = 0                                  -> code
                                                     -> blank
    code_lines = 0                                   -> code
    """
    #if fnames==[]:
    #    fnames = sys.argv[1:]
    import os
    from applyLicense import getFoldersFiles
    folders,fnames = getFoldersFiles(os.path.join('..','pyNastran'))

    total_lines = 0
    comment_lines = 0
    blank_lines = 0
    code_lines = 0
    header_lines = 0
    format_lines = 0
    doxygen_lines = 0
    bad_doxygen_lines = 0
    quote_lines = 0
    doxygen_starter_lines = 0
    
    input_file_list = []

    print
    
    if 0 and len(sys.argv) < 2 and len(fnames)==0:
        try:
            inputs_file = open("inputs.txt", "r")
            for line in inputs_file.readlines():
                if line[-1] == "\n":
                    input_file_list.append(line[:-1])
                else:
                    input_file_list.append(line)
            inputs_file.close()
            if len(input_file_list) == 0:
                print usage
                sys.exit(1)
        except IOError:
            print usage
            sys.exit(1)
    else:
        input_file_list = fnames
        
    print "Input files:"
    
    for input_file_name in input_file_list:
        if not input_file_name.endswith('.py'):
            continue
        print input_file_name

        try:
            current_file = open(input_file_name, "r")
        except IOError:
            print "Error opening \"" + input_file_name + "\", skipping...\n"
            break
        except TypeError:        
            print "Error opening \"" + str(input_file_name) + "\", skipping...\n"
            break
        i = 0
        toggleComment = False
        for line in current_file.readlines():
            #if i==100:
            #    break
            total_lines += 1
            i+=1
            hasPound = False
            if '"""' in line or "'''" in line:
                nQuotes = line.count('"""')+line.count("'''")
                if nQuotes == 2:
                    header_lines += 1  # inline """comment"""
                    continue
                
                if toggleComment: # if already on
                    header_lines += 1
                toggleComment = not(toggleComment)
                
            
            if toggleComment:
                #line2 = line.strip(' \n\t\r').strip('"').strip("'")
                #if not line2:
                header_lines += 1
                #line = line2
                continue
            
            for char in line:
                if char == "\n":
                    blank_lines += 1
                    break
                if char == "#":
                    hasPound = True
                    #comment_lines += 1
                    break
                elif (char not in string.whitespace) and (char != "#"):
                    code_lines += 1
                    break
            if hasPound:
                n  = line.count('#')
                n2 = line.count('##')
                n3 = line.count('###')
                
                if n3==1:
                    format_lines += 1
                    print input_file_name, "###"
                elif n==1 and '@' in line and '<' not in line:
                    bad_doxygen_lines += 1
                    print input_file_name, "****bad****",i
                elif n2==1 and '@' in line:
                    doxygen_starter_lines += 1
                elif n2==1:
                    doxygen_lines += 1
                else:
                    comment_lines += 1
                ###
            ###
            
        current_file.close()
        
    total_lines2 = (blank_lines + comment_lines + code_lines +
                    header_lines + doxygen_lines + format_lines)
    
    print "Total lines:   " + str(total_lines)
    print "Total lines2:  " + str(total_lines2)
    print "Code lines:    " + str(code_lines)
    print "Header lines:  " + str(header_lines)
    print "Comment lines: " + str(comment_lines)
    print "Doxygen start: " + str(doxygen_starter_lines)
    print "Doxygen lines: " + str(doxygen_lines)
    print "Bad Doxy lines " + str(bad_doxygen_lines)
    print "Blank lines:   " + str(blank_lines)
    print "Format lines:  " + str(format_lines)

if __name__ == "__main__":
    countLines()
