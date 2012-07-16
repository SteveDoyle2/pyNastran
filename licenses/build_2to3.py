from __future__ import absolute_import

import sys
import os
from applyLicense import getFoldersFiles

#assert sys.version_info[:2]==(2,7),'this must be done in python 2.7'

def fix_object(files):
    for fname in files:
        if '.py' in fname:
            print "adding fname=%s" %(fname)
            f = open(fname,'r')
            lines = f.readlines()
            f.close()
            #print lines
            lines2 = []
            #print "len(lines) = ",len(lines)
            for codeLine in lines:
                codeLine.rstrip()
                if '#' in codeLine:
                    isPound = True
                    indx = codeLine.index('#')
                    sline = [codeLine[:indx],codeLine[indx:]]
                    try:
                        line = sline[0]
                        comment = sline[1]
                    except ValueError:
                        line = sline[0]
                        comment = ''
                        raise

                    nspaces = len(codeLine.strip())-len(codeLine.rstrip())
                    #indx = codeLine.index('#')
                else:
                    isPound = False
                    line = codeLine

                #line = ''
                #if '(object):' in line:
                    #line = line.replace('(object):','():')
                if 'nested_scopes' in line or print_function in line:
                    line = ''
                elif 'unicode' in line:
                    line = line.replace('unicode','str')
                elif '.iteritems()' in line:
                    line = line.replace('.iteritems()','.items()')
                elif '.itervalues()' in line:
                    line = line.replace('.itervalues()','.values()')
                elif '.iterkeys()' in line:
                    line = line.replace('.iterkeys()','.keys()')
                elif '.izip()' in line:
                    line = line.replace('.izip()','.zip()')

                #elif ",'wb')" in line:
                    #assert 'open' in line,codeLine
                    #line = line.replace(",'wb')",",'wb',encoding='ascii')")
                #elif ",'r')" in line:
                    #assert 'open' in line,codeLine
                    #line = line.replace(",'r')",",'r',encoding='ascii')")
                #elif 'import' in line and ('bdf' in line or 'op2' in line or 'f06' in line or 'gui' in line or 'op4' in line):
                    #assert 'pyNastran.' in line,line # validates that the import is correct
                #elif 'print '
                #elif 
                if isPound:
                    line = ' '*nspaces + line + comment #[indx:]
                #if line != codeLine:
                    #print line.rstrip()
                    #print codeLine.rstrip()
                    #adsf
                lines2.append(line)
            
            #print "len(lines2) = ",len(lines2)
            f = open(fname,'wb')
            #f = open('bbb.py','wb')
            for line in lines2:
                f.write(line.rstrip('\n\r')+'\n')
                #print line.rstrip()
            ###
            #sys.exit(fname)
            f.close()
        ###
    ###

if __name__ == "__main__":
    (folders,files) = getFoldersFiles('../pyNastran')
    fix_object(files)
