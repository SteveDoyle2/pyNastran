import sys
from datetime import date

import pyNastran

def makeStamp(Title):
    t = date.today()
    months = ['January','February','March','April','May','June','July','August','September','October','November','December']
    today = '%-9s %s, %s' %(months[t.month],t.day,t.year)
    
    releaseDate = '02/08/12'#pyNastran.__releaseDate__
    releaseDate = ''
    build = 'pyNastran v%s %s' %(pyNastran.__version__,releaseDate)
    out = '1    %-67s %20s  %-22s PAGE ' %(Title,today,build)
    return out

def makePyNastranTitle():
    n = ''
    lines1 = [
    n+'/* -------------------------------------------------------------------  */\n',
    n+'/*                              PYNASTRAN                               */\n',
    n+'/*                      - NASTRAN FILE INTERFACE -                      */\n',
    n+'/*                                                                      */\n',
    n+'/*              A Python reader/editor/writer for the various           */\n',
    n+'/*                        NASTRAN file formats.                         */\n',
    n+'/*                  Copyright (C) 2011-2012 Steven Doyle                */\n',
    n+'/*                                                                      */\n',
    n+'/*    This program is free software; you can redistribute it and/or     */\n',
    n+'/*    modify it under the terms of the GNU Lesser General Public        */\n',
    n+'/*    License as published by the Free Software Foundation;             */\n',
    n+'/*    version 3 of the License.                                         */\n',
    n+'/*                                                                      */\n',
    n+'/*    This program is distributed in the hope that it will be useful,   */\n',
    n+'/*    but WITHOUT ANY WARRANTY; without even the implied warranty of    */\n',
    n+'/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */\n',
    n+'/*    GNU Lesser General Public License for more details.               */\n',
    n+'/*                                                                      */\n',
    n+'/*    You should have received a copy of the GNU Lesser General Public  */\n',
    n+'/*    License along with this program; if not, write to the             */\n',
    n+'/*    Free Software Foundation, Inc.,                                   */\n',
    n+'/*    675 Mass Ave, Cambridge, MA 02139, USA.                           */\n',
    n+'/* -------------------------------------------------------------------  */\n',
    '\n']

    n = 46*' '
    lines2 = [
        n+'* * * * * * * * * * * * * * * * * * * *\n',
        n+'* * * * * * * * * * * * * * * * * * * *\n',
        n+'* *                                 * *\n',
        n+'* *            pyNastran            * *\n',
        n+'* *                                 * *\n',
        n+'* *                                 * *\n',
        n+'* *                                 * *\n',
        n+'* *        Version %8s       * *\n' %(pyNastran.__version__),
        n+'* *                                 * *\n',
        n+'* *                                 * *\n',
        n+'* *          %15s         * *\n' %(pyNastran.__releaseDate2__),
        n+'* *                                 * *\n',
        n+'* *            Questions            * *\n',
        n+'* *        mesheb82@gmail.com       * *\n',
        n+'* *                                 * *\n',
        n+'* *                                 * *\n',
        n+'* * * * * * * * * * * * * * * * * * * *\n',
        n+'* * * * * * * * * * * * * * * * * * * *\n\n\n']
    return ''.join(lines1+lines2)

def makeEnd():
    lines = [' \n'
             '1                                        * * * END OF JOB * * *\n'
             ' \n'
             ' \n']
    return ''.join(lines)

class F06Writer(object):
    def __init__(self,model='tria3'):
        self.Title = ''
        self.setF06Name(model)

    def setF06Name(self,model):
        self.model = model
        self.f06OutName = '%s.f06.out' %(self.model)

    def loadOp2(self,isTesting=False):
        if isTesting==False:  ## @todo implement in way that doesnt require a variable (e.g. check parent class)
            raise Exception("Don't call this method unless you're testing the F06Writer.  It breaks the F06 and OP2 classes.")
        from pyNastran.op2.op2 import OP2
        self.op2Name = model+'.op2'
        op2 = OP2(self.op2Name)
        op2.readOP2()

        # oug
        self.eigenvectors  = op2.eigenvectors
        self.displacements = op2.displacements
        self.temperatures  = op2.temperatures
        
        # oes
        #CBEAM
        #CSHEAR
        #CELASi
        self.rodStress            = op2.rodStress
        self.rodStrain            = op2.rodStrain
        self.barStress            = op2.barStress
        self.barStrain            = op2.barStrain
        self.plateStress          = op2.plateStress
        self.plateStrain          = op2.plateStrain
        self.compositePlateStress = op2.compositePlateStress
        self.compositePlateStrain = op2.compositePlateStrain
    
    def writeF06(self):
        f = open(self.f06OutName,'wb')
        f.write(makePyNastranTitle())
        
        #pageStamp = '1    MSC.NASTRAN JOB CREATED ON 10-DEC-07 AT 09:21:23                      NOVEMBER  14, 2011  MSC.NASTRAN  6/17/05   PAGE '
        
        #Title = 'MSC.NASTRAN JOB CREATED ON 10-DEC-07 AT 09:21:23'
        pageStamp = makeStamp(self.Title)
        #print "pageStamp = |%r|" %(pageStamp)
        #print "stamp     = |%r|" %(stamp)

        pageNum = 1
        for case,result in sorted(self.eigenvectors.items()): # has a special header
            header = ['0                                                                                                            SUBCASE %i'%(case)]
            (msg,pageNum) = result.writeF06(header,pageStamp,pageNum=pageNum)
            f.write(msg)
            pageNum +=1

        # subcase name, subcase ID, transient word & value
        header = ['     DEFAULT                                                                                                                        ',
                  '','']

        resTypes = [self.displacements,self.temperatures,
                    self.rodStrain,self.barStrain,self.beamStrain,self.shearStrain,self.plateStrain,self.compositePlateStrain,#self.solidStrain,  # beam, shear...not done
                    self.rodStress,self.barStress,self.beamStress,self.shearStress,self.plateStress,self.compositePlateStress,self.solidStress,
                    self.ctriaxStress,
                    self.ctriaxStrain,
                    self.gridPointForces,
                    self.loadVectors,
                    ]

        iSubcases = self.iSubcaseNameMap.keys()
        for iSubcase in iSubcases:
            (subtitle,label) = self.iSubcaseNameMap[iSubcase]
            subtitle = subtitle.strip()
            label = label.strip()
            header[0] = '     %s\n' %(subtitle)
            header[1] = '0                                                                                                            SUBCASE %i\n' %(iSubcase)
            for resType in resTypes:
                if resType.has_key(iSubcase):
                    result = resType[iSubcase]
                    try:
                        (msg,pageNum) = result.writeF06(header,pageStamp,pageNum=pageNum)
                    except:
                        print "result name = %s" %(result.name())
                        raise
                    f.write(msg)
                    pageNum +=1
                ###
        ###

        f.write(makeEnd())
        f.close()
        

if __name__=='__main__':
    #Title = 'MSC.NASTRAN JOB CREATED ON 10-DEC-07 AT 09:21:23'
    #stamp = makeStamp(Title) # doesnt have pageNum
    #print "|%s|" %(stamp+'22')
    
    model = sys.argv[1]
    f06 = F06Writer(model)
    f06.loadOp2(isTesting=True)
    f06.writeF06()


