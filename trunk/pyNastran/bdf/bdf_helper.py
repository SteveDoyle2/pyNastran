import sys
import copy
from pyNastran.bdf.fieldWriter import printCard
from pyNastran.bdf.bdf_writeMesh import writeMesh
from pyNastran.bdf.bdf_cardMethods import cardMethods
from pyNastran.bdf.crossReference import XrefMesh

class getMethods(object):
    def __init__(self):
        pass

    #--------------------
    # NODE CARDS

    def nNodes(self):
        return len(self.nodes)
        
    def nodeIDs(self):
        return self.nodes.keys()
        
    def getNodes(self):
        nodes = []
        for nid,node in sorted(self.nodes.items()):
            nodes.append(node)
        return nodes

    def getNodeIDsWithElement(self,eid):
        return self.getNodeIDsWithElements([eid])

    def getNodeIDsWithElements(self,eids):
        nids2 = set([])
        for eid in eids:
            element = self.Element(eid)
            self.log.debug("element.pid = %s" %(element.pid))
            nids = set(element.nids)
            nids2 = nids2.union(nids)
        ###
        return nids2

    def Node(self,nid):
        return self.nodes[nid]

    def Nodes(self,nids):
        """
        Returns a series of node objects given a list of node IDs
        """
        #print "nids",nids
        nodes = []
        for nid in nids:
            nodes.append(self.nodes[nid])
        return nodes

    #--------------------
    # ELEMENT CARDS

    def nElements(self):
        return len(self.elements)

    def elementIDs(self):
        return self.elements.keys()

    def getElementIDsWithPID(self,pid):
        return self.getElementIDsWithPIDs([pid])

    def getElementIDsWithPIDs(self,pids):
        #self.log.info("pids = %s" %(pids))

        eids = self.elementIDs()
        eids2 = []
        #print "eids = ",eids
        for eid in eids:
            element = self.Element(eid)
            if element.Pid() in pids:
                eids2.append(eid)
            ###
            #print ""
        ###
        return (eids2)

    def Element(self,eid):
        return self.elements[eid]

    def Elements(self,eids):
        elements = []
        for eid in eids:
            elements.append(self.elements[eid])
        return elements

    #--------------------
    # PROPERTY CARDS

    def propertyIDs(self):
        return self.properties.keys()

    def Property(self,pid):
        return self.properties[pid]

    def Properties(self,pids):
        properties = []
        for pid in pids:
            properties.append(self.properties[pid])
        return properties

    def Phbdy(self,pid):
        return self.phbdys[pid]

    #--------------------
    # MATERIAL CARDS

    def structuralMaterialIDs(self):
        return self.materials.keys()

    def materialIDs(self):
        return self.materials.keys()+self.thermalMaterials.keys()

    def thermalMaterialIDs(self):
        return self.thermalMaterials.keys()

    def Material(self,mid):
        if mid in self.materials:
            return self.materials[mid]
        elif mid in self.thermalMaterials:
            return self.thermalMaterials[mid]
        else:
            raise KeyError('Invalid Material ID:  mid=%s' %(mid))

    def StructuralMaterial(self,mid):
        return self.materials[mid]

    def ThermalMaterial(self,mid):
        return self.thermalMaterials[mid]

    def Materials(self,mids):
        materials = []
        for mid in mids:
            materials.append(self.Material(mid))
        return materials

    #--------------------
    # LOADS/CONSTRAINTS/COORDINATES CARDS

    def Load(self,lid):
        return self.loads[lid]

    def Coord(self,cid):
        return self.coords[cid]

    #--------------------
    # AERO CARDS

    def Aero(self,acsid):
        return self.aero[acsid]

    def Aeros(self,acsid):
        return self.aeros[acsid]

    def AEStat(self,eid):
        return self.aestats[eid]

    def CAero(self,eid):
        return self.caeros[eid]

    def Flfact(self,sid):
        return self.flfacts[sid]

    def Flutter(self,fid):
        return self.flutters[fid]

    def Gust(self,sid):
        return self.gusts[sid]

    def Spline(self,eid):
        return self.splines[eid]

    #--------------------
    # OPTIMIZATION CARDS

    def DConstr(self,oid):
        return self.dconstrs[oid]

    def Desvar(self,oid):
        return self.desvars[oid]

    def DDVal(self,oid):
        return self.ddvals[oid]

    #--------------------
    # NONLINEAR CARDS

    def NLParm(self,nid):
        return self.nlparms[nid]

    #--------------------
    # METHODS

    def MassProperties(self):
                 #Ixx Iyy Izz, Ixy, Ixz Iyz
        I = array(0., 0., 0.,  0.,  0., 0.,)
        for element in self.elements:
            p = e.Centroid()  # not really coded across the board
            m = e.Mass()
            (x,y,z) = p
            I[0] = m*x*x  # Ixx
            I[1] = m*y*y  # Iyy
            I[2] = m*z*z  # Izz
            I[3] = m*x*y  # Ixy
            I[4] = m*x*z  # Ixz
            I[5] = m*y*z  # Iyz
        ###
        return I
            
    def sumForces(self):
        for key,loadCase in self.loads.items():
            F = array([0.,0.,0.])
            #print "loadCase = ",loadCase
            for load in loadCase:
                #print "load = ",load
                if isinstance(load,Force):
                    f = load.mag*load.xyz
                    print "f = ",f
                    F += f
                ###
            self.log.info("case=%s F=%s\n\n" %(key,F))
        ###

    def sumMoments(self):
        p = array([0.,0.5,0.])
        for key,loadCase in self.loads.items():
            M = array([0.,0.,0.])
            F = array([0.,0.,0.])
            #print "loadCase = ",loadCase
            for load in loadCase:
                #print "load = ",load
                if isinstance(load,Force):
                    f = load.mag*load.xyz
                    node = self.Node(load.node)
                    #print "node = ",node
                    r = node.Position() - p
                    m = cross(r,f)
                    #print "m    = ",m
                    M += m
                    F += f
                elif isinstance(load,Moment):
                    m = load.mag*load.xyz
                    M += m
                ###
            print "case=%s F=%s M=%s\n\n" %(key,F,M)
        ###

class addMethods(object):
    def __init__(self):
        pass

    def addParam(self,param):
        key = param.key
        if key in self.params:
            if not param.isSameCard(self.params[key]):
                assert param.key not in self.params,'key=%s param=%s oldPARAM=%s' %(key,param,self.params[key])
        else:
            self.params[key] = param
        ###

    def addNode(self,node):
        #print node
        #assert node.nid not in self.nodes,'nid=%s\noldNode=\n%snewNode=\n%s' %(node.nid,self.nodes[node.nid],node)  ## @todo enable before release...
        assert node.nid>0,'nid=%s node=\n%s' %(node.nid,node)
        self.nodes[node.nid] = node

    def addElement(self,elem,allowOverwrites=False):
        key = elem.eid
        if key in self.elements and not allowOverwrites:
            print 'eid=%s\noldElement=\n%snewElement=\n%s' %(key,self.elements[key],elem)
            #assert elem.eid not in self.elements,'eid=%s\noldElement=\n%snewElement=\n%s' %(elem.eid,self.elements[elem.eid],elem)
            pass
        assert key>0,'eid=%s elem=%s' %(key,elem)
        self.elements[key] = elem

    def addThermalElement(self,elem):
        """same as addElement at the moment..."""
        self.addElement(elem)
        #assert elem.eid not in self.elements
        #assert elem.eid>0
        #self.elements[elem.eid] = elem

    def addProperty(self,prop,allowOverwrites=False):
        if not allowOverwrites:
            assert prop.pid not in self.properties,'pid=%s oldProperty=\n%snewProperty=\n%s' %(prop.pid,self.properties[prop.pid],prop)
        assert prop.pid>0,'property=\n%s' %(prop.pid,prop)
        self.properties[prop.pid] = prop

    def addMaterial(self,material,allowOverwrites=False):
        """
        only for adding structural materials
        @deprecated this method will be renamed in v0.3 to addStructuralMaterial.
        """
        self.addStructuralMaterial(self,material,allowOverwrites)

    def addStructuralMaterial(self,material,allowOverwrites=False):
        key = material.mid
        if key in self.materials:
            if not material.isSameCard(self.materials[key]):
                assert key not in self.materials,'mid=%s\noldMaterial=\n%snewMaterial=\n%s' %(key,self.materials[key],material)
        else:
            assert key>0,'mid=%s material=\n%s' %(key,material)
            self.materials[key] = material

    def addThermalMaterial(self,material,allowOverwrites=False):
        key = material.mid
        if key in self.thermalMaterials:
            if not material.isSameCard(self.thermalMaterials[key]):
                assert key not in self.thermalMaterials,'mid=%s\noldMaterial=\n%snewMaterial=\n%s' %(key,self.thermalMaterials[key],material)
        else:
            assert key>0,'mid=%s material=\n%s' %(key,material)
            self.thermalMaterials[key] = material

    def addCreepMaterial(self,material,allowOverwrites=False):
        """
        @note
            May be removed in the future.  Are CREEP cards materials?
            They have an MID, but reference structural materials.
        """
        key = material.mid
        if key in self.thermalMaterials:
            if not material.isSameCard(self.creepMaterials[key]):
                assert key not in self.creepMaterials,'mid=%s\noldMaterial=\n%snewMaterial=\n%s' %(key,self.creepMaterials[key],material)
        else:
            assert key>0,'mid=%s material=\n%s' %(key,material)
            self.creepMaterials[key] = material

    def addCoord(self,coord,allowOverwrites=False):
        key = coord.cid
        if not allowOverwrites:
            assert key not in self.coords,'cid=%s\noldElement=\n%snewElement=\n%s' %(key,self.coords[key],coord)
        assert coord.cid>-1,'cid=%s coord=\n%s' %(key,coord)
        self.coords[key] = coord

    def addLoad(self,load):
        key = load.lid
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]

    def addPHBDY(self,prop):
        assert prop.pid>0
        assert prop.pid not in self.phbdys
        self.phbdys[prop.pid] = prop

    def addConvectionProperty(self,prop):
        assert prop.pconid>0
        assert prop.pconid not in self.convectionProperties
        self.convectionProperties[prop.pconid] = prop

    #def addThermalProperty(self,prop):
    #    assert prop.pconid not in self.thermalProperties
    #    self.thermalProperties[prop.pconid] = prop

    def addThermalBC(self,bc,key):
        assert key>0
        if key in self.bcs:
            self.bcs[key].append(bc)
        else:
            self.bcs[key] = [bc]

    def addThermalLoad(self,load):  # same function at the moment...
        key = load.sid
        assert key>0
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]

    def addConstraint_MPCADD(self,constraint):
        self.mpcObject.add(constraint)

    def addConstraint_MPC(self,constraint):
        self.mpcObject.append(constraint)

    def addConstraint_SPCADD(self,constraint):
        self.spcObject.add(constraint)

    def addConstraint_SPC(self,constraint):
        self.spcObject.append(constraint)
        #key = constraint.cid
        #if self.constraints.has_key(key):
        #    self.constraints[key].append(constraint)
        #else:
        #    self.constraints[key] = [constraint]

    def addConstraint(self,constraint):
        #self.spcObject.append(constraint)
        key = constraint.cid
        assert key>0
        if self.constraints.has_key(key):
            self.constraints[key].append(constraint)
        else:
            self.constraints[key] = [constraint]

    def addSUPORT(self,suport):
        self.suports.append(suport)

    def addDArea(self,darea):
        key = (darea.sid,darea.p)
        if key in self.dareas:
            if not darea.isSameCard(self.dareas[key]):
                assert key not in self.dareas,'\ndarea=\n%s oldDArea=\n%s' %(darea,self.dareas[key])
        else:
            assert darea.sid>0
            self.dareas[key] = darea
        ###

    def addAero(self,aero):
        assert aero.acsid not in self.aero,'\naero=\n%s oldAERO=\n%s' %(aero,self.aero[aero.acsid])
        assert aero.acsid>=0
        self.aero[aero.acsid] = aero

    def addAeros(self,aero):
        assert aero.acsid not in self.aeros,'\naeros=\n%s oldAEROS=\n%s' %(aero,self.aeros[aero.acsid])
        assert aero.acsid>=0
        self.aeros[aero.acsid] = aero

    def addAEStat(self,aestat):
        assert aestat.id not in self.aestats,'\naestat=%s oldAESTAT=\n%s' %(aestat,self.aestats[aestat.id])
        assert aestat.id>=0
        self.aestats[aestat.id] = aestat

    def addCAero(self,caero):
        assert caero.eid not in self.caeros,'\nself.caeros=|%s| caero.eid=|%s|' %(self.caeros,caero.eid)
        assert caero.eid>0
        self.caeros[caero.eid] = caero

    def addGrav(self,grav):
        assert grav.sid not in self.gravs
        assert grav.sid>0
        self.gravs[grav.sid] = grav

    def addSpline(self,spline):
        assert spline.eid not in self.splines
        assert spline.eid>0
        self.splines[spline.eid] = spline

    def addGust(self,gust):
        assert gust.sid not in self.gusts
        assert gust.sid>0
        self.gusts[gust.sid] = gust

    def addFlutter(self,flutter):
        assert flutter.sid not in self.flutters
        assert flutter.sid>0
        self.flutters[flutter.sid] = flutter

    def addFLFACT(self,flfact):
        #assert flfact.sid not in self.flfacts
        assert flfact.sid>0
        self.flfacts[flfact.sid] = flfact # set id...
        #print "added flfact...flflact =\n"+flfact

    def addDConstr(self,dconstr):
        key = (dconstr.oid,dconstr.rid)
        assert key not in self.dconstrs
        assert dconstr.oid>0
        assert dconstr.rid>0
        self.dconstrs[key] = dconstr

    def addDesvar(self,desvar):
        assert desvar.oid not in self.desvars
        assert desvar.oid>0
        self.desvars[desvar.oid] = desvar

    def addDDVal(self,ddval):
        assert ddval.oid not in self.ddvals
        assert ddval.oid>0
        self.ddvals[ddval.oid] = ddval

    def addDResp(self,dresp):
        assert dresp.oid not in self.dresps
        assert dresp.oid>0
        self.dresps[dresp.oid] = dresp

    def addDvprel(self,dvprel):
        assert dvprel.oid not in self.dvprels
        assert dvprel.oid>0
        self.dvprels[dvprel.oid] = dvprel

    def addNLParm(self,nlparm):
        assert nlparm.nid not in self.nlparms
        assert nlparm.nid>0
        self.nlparms[nlparm.nid] = nlparm
