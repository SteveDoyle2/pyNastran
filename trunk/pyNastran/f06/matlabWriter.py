from __future__ import print_function
import sys
import copy
from datetime import date

import pyNastran


def makeStamp(Title):
    #pageStamp = '1    MSC.NASTRAN JOB CREATED ON 10-DEC-07 AT 09:21:23                      NOVEMBER  14, 2011  MSC.NASTRAN  6/17/05   PAGE '
    #Title = 'MSC.NASTRAN JOB CREATED ON 10-DEC-07 AT 09:21:23'
    t = date.today()
    months = ['January', 'February', 'March', 'April', 'May', 'June',
              'July', 'August', 'September', 'October', 'November', 'December']
    today = '%-9s %s, %s' % (months[t.month - 1], t.day, t.year)

    releaseDate = '02/08/12'  # pyNastran.__releaseDate__
    releaseDate = ''
    build = 'pyNastran v%s %s' % (pyNastran.__version__, releaseDate)
    out = '1    %-67s %20s  %-22s PAGE ' % (Title, today, build)
    return out


def makeF06Header():
    n = ''
    lines1 = [
        n + '/* -------------------------------------------------------------------  */\n',
        n + '/*                              PYNASTRAN                               */\n',
        n + '/*                      - NASTRAN FILE INTERFACE -                      */\n',
        n + '/*                                                                      */\n',
        n + '/*              A Python reader/editor/writer for the various           */\n',
        n + '/*                        NASTRAN file formats.                         */\n',
        n + '/*                  Copyright (C) 2011-2012 Steven Doyle                */\n',
        n + '/*                                                                      */\n',
        n + '/*    This program is free software; you can redistribute it and/or     */\n',
        n + '/*    modify it under the terms of the GNU Lesser General Public        */\n',
        n + '/*    License as published by the Free Software Foundation;             */\n',
        n + '/*    version 3 of the License.                                         */\n',
        n + '/*                                                                      */\n',
        n + '/*    This program is distributed in the hope that it will be useful,   */\n',
        n + '/*    but WITHOUT ANY WARRANTY; without even the implied warranty of    */\n',
        n + '/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */\n',
        n + '/*    GNU Lesser General Public License for more details.               */\n',
        n + '/*                                                                      */\n',
        n + '/*    You should have received a copy of the GNU Lesser General Public  */\n',
        n + '/*    License along with this program; if not, write to the             */\n',
        n + '/*    Free Software Foundation, Inc.,                                   */\n',
        n + '/*    675 Mass Ave, Cambridge, MA 02139, USA.                           */\n',
        n + '/* -------------------------------------------------------------------  */\n',
        '\n']

    n = 46 * ' '
    lines2 = [
        n + '* * * * * * * * * * * * * * * * * * * *\n',
        n + '* * * * * * * * * * * * * * * * * * * *\n',
        n + '* *                                 * *\n',
        n + '* *            pyNastran            * *\n',
        n + '* *                                 * *\n',
        n + '* *                                 * *\n',
        n + '* *                                 * *\n',
        n + '* *        Version %8s       * *\n' % (pyNastran.__version__),
        n + '* *                                 * *\n',
        n + '* *                                 * *\n',
        n + '* *          %15s        * *\n' % (pyNastran.__releaseDate2__),
        n + '* *                                 * *\n',
        n + '* *            Questions            * *\n',
        n + '* *        mesheb82@gmail.com       * *\n',
        n + '* *                                 * *\n',
        n + '* *                                 * *\n',
        n + '* * * * * * * * * * * * * * * * * * * *\n',
        n + '* * * * * * * * * * * * * * * * * * * *\n\n\n']
    return ''.join(lines1 + lines2)


def makeEnd():
    lines = [' \n'
             '1                                        * * * END OF JOB * * *\n'
             ' \n'
             ' \n']
    return ''.join(lines)


class MatlabWriter(object):
    def __init__(self, model='tria3'):
        self.Title = ''
        self.setF06Name(model)

    def setF06Name(self, model):
        self.model = model
        self.f06OutName = '%s.f06.out' % (self.model)

    def loadOp2(self, isTesting=False):
        print("self.class = ",self.__class__.__name__)
        if isTesting == False:  # @todo implement in way that doesnt require a variable (e.g. check parent class)
            msg = ("Don't call this method unless you're testing the "
                   "F06Writer.  It breaks the F06 and OP2 classes.")
            raise RuntimeError(msg)
        from pyNastran.op2.op2 import OP2
        self.op2Name = model + '.op2'
        op2 = OP2(self.op2Name)
        op2.readOP2()

        # oug
        self.eigenvectors = op2.eigenvectors
        self.displacements = op2.displacements
        self.temperatures = op2.temperatures

        # oes
        #CBEAM
        #CSHEAR
        #CELASi
        self.rodStress = op2.rodStress
        self.rodStrain = op2.rodStrain
        self.barStress = op2.barStress
        self.barStrain = op2.barStrain
        self.plateStress = op2.plateStress
        self.plateStrain = op2.plateStrain
        self.compositePlateStress = op2.compositePlateStress
        self.compositePlateStrain = op2.compositePlateStrain

    def makeF06Header(self):
        """If this class is inherited, the F06 Header may be overwritten"""
        return makeF06Header()

    def makeStamp(self, Title):
        """If this class is inherited, the PAGE stamp may be overwritten"""
        return makeStamp(Title)

    def writeMatlab(self, mFileOutName, isMagPhase=False, deleteObjects=True):
        """
        Writes an F06 file based on the data we have stored in the object
        @param self
               the object pointer
        @param mFileOutName
               the name of the M (Matlab) file to write
        @param isMagPhase
               should complex data be written using
               Magnitude/Phase instead of Real/Imaginary (default=False; Real/Imag)
               Real objects don't use this parameter.
        """
        f = open(mFileOutName, 'wb')

        strLines = self.makeF06Header()
        lines = strLines.split('\n')
        for line in lines:
            f.write('%% %s\n' % (line))

        #pageStamp = self.makeStamp(self.Title)
        #print "pageStamp = |%r|" %(pageStamp)
        #print "stamp     = |%r|" %(stamp)

        #isMagPhase = False
        pageNum = 1
        header = ['     DEFAULT                                                                                                                        \n',
                  '\n']

        f.write("fem.title = '%s';\n" % (self.Title))

        subtitleMsg = 'fem.subtitles = {'
        labelMsg = 'fem.labels = {'
        iSubcaseMsg = 'fem.iSubcases = ['

        nSubcases = {}
        for n, (iSubcase, (subtitle, label)) in enumerate(sorted(self.iSubcaseNameMap.items())):
            labelMsg += "'%s'," % (label)
            subtitleMsg += "'%s'," % (subtitle)
            iSubcaseMsg += '%s,' % (iSubcase)
            nSubcases[iSubcase] = n + 1
        labelMsg += '};\n'
        subtitleMsg += '};\n'
        iSubcaseMsg += '];\n'
        f.write(labelMsg)
        f.write(subtitleMsg)
        f.write(iSubcaseMsg)

        for iSubcase, result in sorted(self.eigenvalues.iteritems()):  # goes first
            (subtitle, label) = self.iSubcaseNameMap[iSubcase]
            subtitle = subtitle.strip()
            header[0] = '     %s\n' % (subtitle)
            header[1] = '0                                                                                                            SUBCASE %i\n \n' % (iSubcase)
            print(result.__class__.__name__)
            nSubcase = nSubcases[iSubcase]
            result.writeMatlab(nSubcase, f=f, isMagPhase=isMagPhase)
           #result.writeMatlab(f=f)

        for iSubcase, result in sorted(self.eigenvectors.iteritems()):  # has a special header
            (subtitle, label) = self.iSubcaseNameMap[iSubcase]
            subtitle = subtitle.strip()
            header[0] = '     %s\n' % (subtitle)
            header[1] = '0                                                                                                            SUBCASE %i\n' % (iSubcase)
            print(result.__class__.__name__)
            nSubcase = nSubcases[iSubcase]
            result.writeMatlab(nSubcase, f=f, isMagPhase=isMagPhase)
           #result.writeMatlab(f=f)

        # subcase name, subcase ID, transient word & value
        headerOld = ['     DEFAULT                                                                                                                        \n',
                     '\n', ' \n']
        header = copy.deepcopy(headerOld)
        resTypes = [
                    self.displacements,  # self.displacementsPSD,self.displacementsATO,self.displacementsRMS,
                    #self.scaledDisplacements, # ???
                    self.velocities, self.accelerations,  # self.eigenvectors,
                    #self.temperatures,
                    #self.loadVectors,self.thermalLoadVectors,
                    #self.forceVectors,

                    self.spcForces, self.mpcForces,

                    #self.barForces,self.beamForces,self.springForces,self.damperForces,
                    #self.solidPressureForces,

                    #self.rodStrain,self.nonlinearRodStress,self.barStrain,self.plateStrain,self.nonlinearPlateStrain,self.compositePlateStrain,self.solidStrain,
                    #self.beamStrain,self.ctriaxStrain,self.hyperelasticPlateStress,

                    #self.rodStress,self.nonlinearRodStrain,self.barStress,self.plateStress,self.nonlinearPlateStress,self.compositePlateStress,self.solidStress,
                    #self.beamStress,self.ctriaxStress,self.hyperelasticPlateStrain,


                    # beam, shear...not done
                    #self.shearStrain,self.shearStress,


                    #self.gridPointStresses,self.gridPointVolumeStresses,#self.gridPointForces,
        ]

        if 1:
            iSubcases = self.iSubcaseNameMap.keys()
            for iSubcase in sorted(iSubcases):
                (subtitle, label) = self.iSubcaseNameMap[iSubcase]
                subtitle = subtitle.strip()
                label = label.strip()
                #print "label = ",label
                header[0] = '     %-127s\n' % (subtitle)
                header[1] = '0    %-72s                                SUBCASE %-15i\n' % (label, iSubcase)
                #header[1] = '0    %-72s                                SUBCASE %-15i\n' %('',iSubcase)
                for resType in resTypes:
                    if iSubcase in resType:
                        result = resType[iSubcase]
                        try:
                            print(result.__class__.__name__)
                            nSubcase = nSubcases[iSubcase]
                            result.writeMatlab(nSubcase, f=f, isMagPhase=False)
                        except:
                            #print "result name = %s" %(result.name())
                            raise
                        ###
                        if deleteObjects:
                            del result
                        ###
                    ###
                ###
            ###
        if 0:
            for res in resTypes:
                for iSubcase, result in sorted(res.iteritems()):
                    (msg, pageNum) = result.writeMatlab(
                        iSubcase, f=f, isMagPhase=False)
                    f.write(msg)
                    pageNum += 1
                ###
            ###
        ###
        f.close()

if __name__ == '__main__':
    #Title = 'MSC.NASTRAN JOB CREATED ON 10-DEC-07 AT 09:21:23'
    #stamp = makeStamp(Title) # doesnt have pageNum
    #print "|%s|" %(stamp+'22')

    model = sys.argv[1]
    f06 = MatlabWriter(model)
    f06.loadOp2(isTesting=True)
    f06.writeMatlab()
