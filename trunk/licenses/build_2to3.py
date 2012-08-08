from __future__ import absolute_import

import os
import sys
from applyLicense import getFoldersFiles

#assert sys.version_info[:2]==(2,7),'this must be done in python 2.7'

def fix_object(files):
    isPython3 = True
    isPython2 = False
    assert isPython3 != isPython2
    
    for fname in files:

        if '.py' in fname and '.pyx' not in fname and 'licenses' not in fname:
            #print fname,
            
            #print os.path.dirname('licenses')==os.path.relpath(os.path.dirname(os.path.abspath(fname)))
            #sys.exit()

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

                if isPython2:
                    #if "utf-8" in line or "elementName" in line:
                    #    print line
                    if 'nested_scopes' in line:
                        line = ''
                    elif 'unicode_literals' in line:
                        line = 'from __future__ import division, print_function'
                    elif 'unicode' in line:
                        line = line.replace('unicode','str')
                    elif "encode('utf-8')" in line:
                        line = line.replace(".encode('utf-8')","")

                elif isPython3:
                    #line = ''
                    #if '(object):' in line:
                        #line = line.replace('(object):','():')
                    if 'nested_scopes' in line:
                        line = ''
                    elif 'unicode_literals' in line:
                        line = ''
                    elif '.iteritems()' in line:
                        line = line.replace('.iteritems()','.items()')
                    elif '.itervalues()' in line:
                        line = line.replace('.itervalues()','.values()')
                    elif '.iterkeys()' in line:
                        line = line.replace('.iterkeys()','.keys()')
                    elif '.izip()' in line:
                        line = line.replace('.izip()','.zip()')

                    elif '=' in line and ('keys()' in line or 'values()' in line) and 'sorted' not in line:
                        i = line.index('=')
                        before,after = line[:i],line[i+1:]
                        #after = after.strip('\r\n')
                        #print "before=|%s| after=|%s|" %(before,after)
                        after = after.strip()
                        line = before + '= list(' + after + ')'

                    elif ('keys()' in line or 'values()' in line) and 'sorted' in line:
                        i = line.index('sorted(')
                        before,after = line[:i+7],line[i+7:] # 7 is the length of 'sorted('
                        #after = after.strip('\r\n')
                        #print "before=|%s| after=|%s|" %(before,after)
                        after = after.strip()
                        line = before + 'list(' + after + ')'
                        print "***",line

                    elif "'wb')" in line:
                        line = line.replace("'wb')","'w')")
                    elif 'format1 = bytes(format1)' in line:
                        line = ''
                    elif 'format2 = bytes(format2)' in line:
                        line = ''
                    elif 'iFormat = bytes(iFormat)' in line:
                        line = ''
                    elif 'formatAll = bytes(formatAll)' in line:
                        line = ''
                    elif 'dataFormat = bytes(dataFormat)' in line:
                        line = ''
                    elif "unpack(b'" in line:
                        line = line.replace("unpack(b'","unpack('")
                    elif " b'" in line and 'data' not in line:
                        line = line.replace(" b'"," '")
                    elif 'unicode' in line:
                        line = line.replace('unicode','str')

                    elif 'keys()' in line:
                        print '??? ',line

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
            f = open(fname, 'wb')
            #f = open('bbb.py', 'wb')
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
