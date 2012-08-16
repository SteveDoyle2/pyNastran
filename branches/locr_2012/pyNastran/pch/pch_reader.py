import os
import sys


class PCH(object):
    def __init__(self):
        pass

    def readPCH(self, pchName):
        self.pchName = pchName
        pch = open(self.pchName, 'r')

        for line in pch.readline()[:72]:
            headerLines = []
            while '$' in line:
                headerLines.append(line)
                line = pch.readline()[:72]
            print "***line = ", line
            for line in headerLines:
                print line

            # read the headers
            headers = {}
            for line in headerLines:
                print "-----"
                print line
                if '=' in line:
                    i = line.index('=')
                    key = line[1:i].strip()
                    value = value = line[i + 1:].strip()
                else:
                    key = line[1:72].strip()
                    value = None
                if key:
                    headers[key] = value
                    print "key=|%s| value=|%s|" % (key, value)

            if 'REAL OUTPUT' in headers:  # MAGNITUDE-PHASE OUTPUT
                print "real",
                if 'DISPLACEMENTS' in headers:
                    print "displacements"
                    print pch.readline().strip(), '***'
                    line = readRealTable(pch, line)
                elif 'VELOCITY' in headers:
                    print "velocity"
                    line = readRealTable(pch, line)
                elif 'ACCELERATION' in headers:
                    print "acceleration"
                    line = readRealTable(pch, line)
                #elif 'OLOADS' in headers:
                    #print "oloads"
                else:
                    raise NotImplementedError(headers.keys())
            elif 'REAL-IMAGINARY OUTPUT':
                print "real-imaginary"
            else:
                raise NotImplementedError(headers.keys())
            sys.exit('done')


            if line == '':  # end of file
                break


def readRealTable(pch, line):
    """
    reads displacemnt, velocity, acceleration, spc/mpc forces,
    temperature, load vector
    """
    #print "?",line
    #a = pch.readline()
    #print "xxxxx",a
    line = ''
    lines = []
    while '$' not in line:
        line = pch.readline()
        print "?|%s|" % (line)
        lines.append(line)

    print '\n'.join(lines)
    return line

if __name__ == '__main__':
    pch = PCH()
    pchname = 'nltnln02.pch'
    pch.readPCH(pchname)

#$TITLE   = QUADR -- THICK CYLINDER ( NU = 0.4999 )                             1
#$SUBTITLE= QAJOB-Q404K003                                                      2
#$LABEL   =                                                                     3
#$DISPLACEMENTS                                                                 4
#$REAL OUTPUT                                                                   5
#$SUBCASE ID =           1                                                      6
#         1       G      5.219747E-03      0.000000E+00      0.000000E+00       7
#-CONT-                  0.000000E+00      0.000000E+00      0.000000E+00       8
#         3       G      4.530559E-03      0.000000E+00      0.000000E+00       9
#-CONT-                  0.000000E+00      0.000000E+00      0.000000E+00      10
#         5       G      3.854469E-03      0.000000E+00      0.000000E+00      11
