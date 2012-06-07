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
import sys
from datetime import date

import pyNastran

def makeStamp(Title):
    #pageStamp = '1    MSC.NASTRAN JOB CREATED ON 10-DEC-07 AT 09:21:23                      NOVEMBER  14, 2011  MSC.NASTRAN  6/17/05   PAGE '
    #Title = 'MSC.NASTRAN JOB CREATED ON 10-DEC-07 AT 09:21:23'
    t = date.today()
    months = ['January','February','March','April','May','June','July','August','September','October','November','December']
    today = '%-9s %s, %s' %(months[t.month-1],t.day,t.year)
    
    releaseDate = '02/08/12'#pyNastran.__releaseDate__
    releaseDate = ''
    build = 'pyNastran v%s %s' %(pyNastran.__version__,releaseDate)
    out = '1    %-67s %20s  %-22s PAGE ' %(Title,today,build)
    return out

def makeF06Header():
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
        n+'* *          %15s        * *\n' %(pyNastran.__releaseDate2__),
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
    
    def makeF06Header(self):
        """If this class is inherited, the F06 Header may be overwritten"""
        return makeF06Header()

    def makeStamp(self,Title):
        """If this class is inherited, the PAGE stamp may be overwritten"""
        return makeStamp(Title)

    def writeF06(self,f06OutName,makeFile=True):
        """
        Writes an F06 file based on the data we have stored in the object
        @param self the object pointer
        @param f06OutName the name of the F06 file to write
        @param makeFile True->makes a file, False->makes a StringIO object for testing (default=True)
        """
        if makeFile:
            f = open(f06OutName,'w',encoding='ascii')
        else:
            import StringIO
            f = StringIO.StringIO()

        f.write(self.makeF06Header())
        pageStamp = self.makeStamp(self.Title)
        #print "pageStamp = |%r|" %(pageStamp)
        #print "stamp     = |%r|" %(stamp)

        pageNum = 1
        header = ['     DEFAULT                                                                                                                        \n',
                  '\n']
        for iSubcase,result in sorted(self.eigenvalues.items()): # goes first
            (subtitle,label) = self.iSubcaseNameMap[iSubcase]
            subtitle = subtitle.strip()
            header[0] = '     %s\n' %(subtitle)
            header[1] = '0                                                                                                            SUBCASE %i\n \n' %(iSubcase)
            (msg,pageNum) = result.writeF06(header,pageStamp,pageNum=pageNum)
            f.write(msg)
            pageNum +=1
        
        for iSubcase,result in sorted(self.eigenvectors.items()): # has a special header
            (subtitle,label) = self.iSubcaseNameMap[iSubcase]
            subtitle = subtitle.strip()
            header[0] = '     %s\n' %(subtitle)
            header[1] = '0                                                                                                            SUBCASE %i\n' %(iSubcase)
            (msg,pageNum) = result.writeF06(header,pageStamp,pageNum=pageNum)
            f.write(msg)
            pageNum +=1

        # subcase name, subcase ID, transient word & value
        header = ['     DEFAULT                                                                                                                        \n',
                  '\n',' \n']

        resTypes = [
                    self.displacements,self.temperatures,
                    self.loadVectors,

                    self.spcForces,self.mpcForces,
                    self.rodStrain,self.nonlinearRodStress,self.barStrain,self.plateStrain,self.nonlinearPlateStrain,self.compositePlateStrain,self.solidStrain,
                    self.beamStrain,self.ctriaxStrain,self.hyperelasticPlateStress,

                    self.rodStress,self.nonlinearRodStrain,self.barStress,self.plateStress,self.nonlinearPlateStress,self.compositePlateStress,self.solidStress,
                    self.beamStress,self.ctriaxStress,self.hyperelasticPlateStrain,
                    
                    
                    # beam, shear...not done
                    #self.shearStrain,self.shearStress,
                    
                    
                    self.gridPointForces,
                    ]

        if 1:
            iSubcases = self.iSubcaseNameMap.keys()
            for iSubcase in sorted(iSubcases):
                (subtitle,label) = self.iSubcaseNameMap[iSubcase]
                subtitle = subtitle.strip()
                label = label.strip()
                header[0] = '     %s\n' %(subtitle)
                header[1] = '0                                                                                                            SUBCASE %i\n' %(iSubcase)
                for resType in resTypes:
                    if iSubcase in resType:
                        result = resType[iSubcase]
                        self.log.info(result.__class__.__name__)
                        try:
                            #print result.__class__.__name__
                            (msg,pageNum) = result.writeF06(header,pageStamp,pageNum=pageNum)
                        except:
                            #print "result name = %s" %(result.name())
                            raise
                        f.write(msg)
                        pageNum +=1
                    ###
                ###
            ###
        if 0:
            for res in resTypes:
                for iSubcase,result in sorted(res.items()):
                    self.log.info(result.__class__.__name__)
                    (msg,pageNum) = result.writeF06(header,pageStamp,pageNum=pageNum)
                    f.write(msg)
                    pageNum +=1
                ###
            ###
        f.write(makeEnd())
        if not makeFile:
            print(f.getvalue())
        ###
        f.close()

if __name__=='__main__':
    #Title = 'MSC.NASTRAN JOB CREATED ON 10-DEC-07 AT 09:21:23'
    #stamp = makeStamp(Title) # doesnt have pageNum
    #print "|%s|" %(stamp+'22')
    
    model = sys.argv[1]
    f06 = F06Writer(model)
    f06.loadOp2(isTesting=True)
    f06.writeF06()


