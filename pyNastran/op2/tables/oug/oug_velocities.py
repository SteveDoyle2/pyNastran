import sys
import copy

# pyNastran
from pyNastran.op2.resultObjects.tableObject import TableObject

class velocityObject(TableObject): # approachCode=10, sortCode=0, thermal=0
    def __init__(self,dataCode,iSubcase,dt=None):
        TableObject.__init__(self,dataCode,iSubcase,dt)

    def writeF06(self,header,pageStamp,pageNum=1):
        if self.dt is not None:
            return self.writeF06Transient(header,pageStamp,pageNum)
        msg = ['                                                   V E L O C I T Y   V E C T O R\n',
               ' \n',
               '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        for nodeID,translation in sorted(self.translations.items()):
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
        words = ['                                                   V E L O C I T Y   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        msg = []
        for dt,translations in sorted(self.translations.items()):
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            msg += header+words
            for nodeID,translation in sorted(translations.items()):
                rotation = self.rotations[dt][nodeID]
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

    def __repr__(self):
        if self.dt is not None:
            return self.__reprTransient__()

        msg = '---VELOCITIES---\n'
        msg += self.writeHeader()

        for nodeID,translation in sorted(self.translations.items()):
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
        msg = '---TRANSIENT VELOCITY---\n'
        msg += self.writeHeader()
        
        for dt,translations in sorted(self.translations.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for nodeID,translation in sorted(translations.items()):
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

