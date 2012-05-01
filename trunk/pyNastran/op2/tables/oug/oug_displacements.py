import sys
import copy

# pyNastran
from pyNastran.op2.resultObjects.tableObject import TableObject,complexTableObject

class displacementObject(TableObject): # approachCode=1, sortCode=0, thermal=0
    def __init__(self,dataCode,iSubcase,dt=None):
        TableObject.__init__(self,dataCode,iSubcase,dt)

    def writeF06(self,header,pageStamp,pageNum=1):
        if self.dt is not None:
            return self.writeF06Transient(header,pageStamp,pageNum)
        msg = header+['                                             D I S P L A C E M E N T   V E C T O R\n',
               ' \n',
               '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        for nodeID,translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx,dy,dz) = translation
            (rx,ry,rz) = rotation
            vals = [dx,dy,dz,rx,ry,rz]
            (vals2,isAllZeros) = self.writeF06Floats13E(vals)
            [dx,dy,dz,rx,ry,rz] = vals2
            msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,dx,dy,dz,rx,ry,rz.rstrip()))
        ###
        msg.append(pageStamp+str(pageNum)+'\n')
        return (''.join(msg),pageNum)

    def writeF06Transient(self,header,pageStamp,pageNum=1):
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        msg = []
        for dt,translations in sorted(self.translations.iteritems()):
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            msg += header+words
            for nodeID,translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                gridType = self.gridTypes[nodeID]

                (dx,dy,dz) = translation
                (rx,ry,rz) = rotation
                vals = [dx,dy,dz,rx,ry,rz]
                (vals2,isAllZeros) = self.writeF06Floats13E(vals)
                [dx,dy,dz,rx,ry,rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,dx,dy,dz,rx,ry,rz.rstrip()))
                pageNum+=1
            ###
            msg.append(pageStamp+str(pageNum)+'\n')
        return (''.join(msg),pageNum-1)

    def __repr__(self):
        if self.dt is not None:
            return self.__reprTransient__()

        msg = '---DISPLACEMENTS---\n'
        msg += self.writeHeader()

        for nodeID,translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx,dy,dz) = translation
            (rx,ry,rz) = rotation

            msg += '%-10i %-8s ' %(nodeID,gridType)
            vals = [dx,dy,dz,rx,ry,rz]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%10s ' %(0)
                else:
                    msg += '%10.3e ' %(val)
                ###
            msg += '\n'
        return msg

    def __reprTransient__(self):
        msg = '---TRANSIENT DISPLACEMENTS---\n'
        msg += self.writeHeader()
        
        for dt,translations in sorted(self.translations.iteritems()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for nodeID,translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                gridType = self.gridTypes[nodeID]
                (dx,dy,dz) = translation
                (rx,ry,rz) = rotation

                msg += '%-10i %8s ' %(nodeID,gridType)
                vals = [dx,dy,dz,rx,ry,rz]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %(0)
                    else:
                        msg += '%10.3e ' %(val)
                    ###
                msg += '\n'
            ###
        return msg

class complexDisplacementObject(complexTableObject): # approachCode=1, sortCode=0, thermal=0
    def __init__(self,dataCode,iSubcase,dt=None):
        complexTableObject.__init__(self,dataCode,iSubcase,dt)

    def writeF06(self,header,pageStamp,pageNum=1):
        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n',
                 '                                                          (REAL/IMAGINARY)\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        msg = []
        for dt,translations in sorted(self.translations.iteritems()):
            header[2] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            msg += header+words
            for nodeID,translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                gridType = self.gridTypes[nodeID]

                (dx,dy,dz) = translation
                dxr=dx.real; dyr=dy.real; dzr=dz.real; 
                dxi=dx.imag; dyi=dy.imag; dzi=dz.imag

                (rx,ry,rz) = rotation
                rxr=rx.real; ryr=ry.real; rzr=rz.real
                rxi=rx.imag; ryi=ry.imag; rzi=rz.imag
                
                vals = [dxr,dyr,dzr,rxr,ryr,rzr,dxi,dyi,dzi,rxi,ryi,rzi]
                (vals2,isAllZeros) = self.writeF06Floats13E(vals)
                [dxr,dyr,dzr,rxr,ryr,rzr,dxi,dyi,dzi,rxi,ryi,rzi] = vals2
                msg.append('0 %12i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,dxr,dyr,dzr,rxr,ryr,rzr.rstrip()))
                msg.append('  %12s %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %('','',          dxi,dyi,dzi,rxi,ryi,rzi.rstrip()))
            ###
            msg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1
        return (''.join(msg),pageNum-1)

    def __repr__(self):
        return self.writeF06(['','',''],'PAGE ',1)[0]

        msg = '---COMPLEX DISPLACEMENTS---\n'
        #if self.dt is not None:
        #    msg += '%s = %g\n' %(self.dataCode['name'],self.dt)
        headers = ['DxReal','DxImag','DyReal','DyImag','DzReal','DyImag','RxReal','RxImag','RyReal','RyImag','RzReal','RzImag']
        msg += '%-10s ' %('nodeID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for freq,translations in sorted(self.translations.iteritems()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)

            for nodeID,translation in sorted(translations.iteritems()):
                rotation = self.rotations[freq][nodeID]

                msg += '%-10i ' %(nodeID)
                vals = translation+rotation
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %(0)
                    else:
                        msg += '%10.3e ' %(val)
                    ###
                msg += '\n'
            ###
        return msg
