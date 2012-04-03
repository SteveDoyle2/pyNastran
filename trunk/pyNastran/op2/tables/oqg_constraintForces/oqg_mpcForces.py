from numpy import array
from pyNastran.op2.resultObjects.tableObject import TableObject,complexTableObject

class mpcForcesObject(TableObject):
    def __init__(self,dataCode,iSubcase,dt=None):
        TableObject.__init__(self,dataCode,iSubcase,dt)

    def writeF06(self,header,pageStamp,pageNum=1):
        if self.dt is not None:
            return self.writeF06Transient(header,pageStamp,pageNum)
        msg = header+['                               F O R C E S   O F   M U L T I - P O I N T   C O N S T R A I N T\n',
               ' \n',
               '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        for nodeID,translation in sorted(self.translations.items()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx,dy,dz) = translation
            (rx,ry,rz) = rotation
            vals = [dx,dy,dz,rx,ry,rz]
            (vals2,isAllZeros) = self.writeF06Floats13E(vals)
            if not isAllZeros:
                [dx,dy,dz,rx,ry,rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,dx,dy,dz,rx,ry,rz.rstrip()))
            ###
        ###
        msg.append(pageStamp+str(pageNum)+'\n')
        return (''.join(msg),pageNum)

    def writeF06Transient(self,header,pageStamp,pageNum=1):
        words = ['                               F O R C E S   O F   M U L T I - P O I N T   C O N S T R A I N T\n',
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
                if not isAllZeros:
                    [dx,dy,dz,rx,ry,rz] = vals2
                    msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,dx,dy,dz,rx,ry,rz.rstrip()))
                ###
            ###
            msg.append(pageStamp+str(pageNum)+'\n')
        return (''.join(msg),pageNum)

    def __reprTransient__(self):
        msg = '---MPC FORCES---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' %(self.dt)

        headers = ['T1','T2','T3','R1','R2','R3']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,translations in sorted(self.translations.items()):
            msg += 'dt = %s' %(dt)
            for nodeID,translation in sorted(translations.items()):
                rotation = self.rotations[dt][nodeID]
                (Fx,Fy,Fz) = translation
                (Mx,My,Mz) = rotation

                msg += '%-8i ' %(nodeID)
                vals = [Fx,Fy,Fz,Mx,My,Mx]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %(0)
                    else:
                        msg += '%10.2f ' %(val)
                    ###
                msg += '\n'
            ###
        return msg

    def __repr__(self):
        if self.dt is not None:
            return self.__reprTransient__()

        msg = '---MPC FORCES---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' %(self.dt)

        headers = ['T1','T2','T3','R1','R2','R3']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for nodeID,translation in sorted(self.translations.items()):
            rotation = self.rotations[nodeID]
            (Fx,Fy,Fz) = translation
            (Mx,My,Mz) = rotation

            msg += '%-8i ' %(nodeID)
            vals = [Fx,Fy,Fz,Mx,My,Mx]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%10s ' %(0)
                else:
                    msg += '%10.2f ' %(val)
                ###
            msg += '\n'
        return msg

class complexMpcForcesObject(complexTableObject):
    def __init__(self,dataCode,iSubcase,dt=None):
        complexTableObject.__init__(self,dataCode,iSubcase,dt)

    def writeF06(self,header,pageStamp,pageNum=1):
        if self.dt is not None:
            return self.writeF06Transient(header,pageStamp,pageNum)
        msg = header+['                               F O R C E S   O F   M U L T I - P O I N T   C O N S T R A I N T\n',
               ' \n',
               '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        for nodeID,translation in sorted(self.translations.items()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx,dy,dz) = translation
            (rx,ry,rz) = rotation
            vals = [dx,dy,dz,rx,ry,rz]
            (vals2,isAllZeros) = self.writeF06Floats13E(vals)
            if not isAllZeros:
                [dx,dy,dz,rx,ry,rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,dx,dy,dz,rx,ry,rz.rstrip()))
        ###
        msg.append(pageStamp+str(pageNum)+'\n')
        return (''.join(msg),pageNum)

    def writeF06Transient(self,header,pageStamp,pageNum=1):
        words = ['                         C O M P L E X   F O R C E S   O F   M U L T I   P O I N T   C O N S T R A I N T\n',
                 '                                                          (REAL/IMAGINARY)\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        msg = []
        for dt,translations in sorted(self.translations.items()):
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            msg+= header+words
            for nodeID,translation in sorted(translations.items()):
                rotation = self.rotations[dt][nodeID]
                gridType = self.gridTypes[nodeID]

                #(dx,dy,dz) = translation
                #(rx,ry,rz) = rotation
                vals = translation+rotation #[dx,dy,dz,rx,ry,rz]
                (vals2,isAllZeros) = self.writeF06Floats13E(vals)
                if not isAllZeros:
                    [v1r,v1i,v2r,v2i,v3r,v3i,v4r,v4i,v5r,v5i,v6r,v6i] = vals2
                    msg.append('0%13i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,v1r,v2r,v3r,v4r,v5r,v6r.rstrip()))
                    msg.append(' %13i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,v1i,v2i,v3i,v4i,v5i,v6i.rstrip()))
                ###
            ###
            msg.append(pageStamp+str(pageNum)+'\n')
        return (''.join(msg),pageNum-1)

    def __reprTransient__(self):
        msg = '---COMPLEX MPC FORCES---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' %(self.dt)

        headers = ['T1','T2','T3','R1','R2','R3']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,translations in sorted(self.translations.items()):
            msg += 'dt = %s' %(dt)
            for nodeID,translation in sorted(translations.items()):
                rotation = self.rotations[dt][nodeID]
                msg += '%-8i ' %(nodeID)
                vals = translation+rotation
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %(0)
                    else:
                        msg += '%10.2f ' %(val)
                    ###
                msg += '\n'
            ###
        return msg

    def __repr__(self):
        if self.dt is not None:
            return self.__reprTransient__()

        msg = '---COMPLEX MPC FORCES---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' %(self.dt)

        headers = ['T1','T2','T3','R1','R2','R3']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for nodeID,translation in sorted(self.translations.items()):
            rotation = self.rotations[nodeID]
            (Fx,Fy,Fz) = translation
            (Mx,My,Mz) = rotation

            msg += '%-8i ' %(nodeID)
            vals = [Fx,Fy,Fz,Mx,My,Mx]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%10s ' %(0)
                else:
                    msg += '%10.2f ' %(val)
                ###
            msg += '\n'
        return msg

