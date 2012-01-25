import sys
import copy
from pyNastran.bdf.fieldWriter     import printCard
#from pyNastran.bdf.bdf_writeMesh   import writeMesh
#from pyNastran.bdf.bdf_cardMethods import cardMethods
#from pyNastran.bdf.crossReference  import XrefMesh

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

    def RigidElement(self,eid):
        return self.rigidElements[eid]

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
    # LOADS/CONSTRAINTS

    def Load(self,lid):
        if lid in self.loads:
            load = self.loads[lid]
        if lid in self.gravs:
            return self.Grav(lid)
        raise Exception('cannot find LoadID=%s' %(lid))

    def Grav(self,sid):
        return self.gravs[sid]

    #--------------------
    # COORDINATES CARDS
    def Coord(self,cid):
        return self.coords[cid]

    #--------------------
    # AERO CARDS

    def nCAeros(self):
        return len(self.caeros)

    def Aero(self,acsid):
        return self.aero[acsid]

    def Aeros(self,acsid):
        return self.aeros[acsid]

    def Spline(self,eid):
        return self.splines[eid]

    def CAero(self,eid):
        return self.caeros[eid]

    def PAero(self,pid):
        return self.paeros[pid]

    def Gust(self,sid):
        return self.gusts[sid]

    #--------------------
    # AERO CONTROL SURFACE CARDS
    def AEStat(self,aid):
        return self.aestats[aid]

    def AELink(self,linkID):
        return self.aelinks[linkID]

    def AEParam(self,aid):
        return self.aeparams[aid]

    #--------------------
    # FLUTTER CARDS

    def Flfact(self,sid):
        return self.flfacts[sid]

    def Flutter(self,fid):
        return self.flutters[fid]

    #--------------------
    # OPTIMIZATION CARDS

    def DConstr(self,oid):
        return self.dconstrs[oid]

    def Desvar(self,oid):
        return self.desvars[oid]

    def DDVal(self,oid):
        return self.ddvals[oid]

    #--------------------
    # SET CARDS

    def Set(self,sid):
        return self.sets[sid]

    def SetSuper(self,seid):
        return self.setsSuper[seid]

    #--------------------
    # NONLINEAR CARDS

    def NLParm(self,nid):
        return self.nlparms[nid]
    ###

    #--------------------

class addMethods(object):
    def __init__(self):
        pass

    def addParam(self,param,allowOverwrites=False):
        key = param.key
        if key in self.params and allowOverwrites==False:
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
        self.elements[key] = elem  ## temporary
        return                     ## temporary

        if key in self.elements and allowOverwrites==False:
            if not elem.isSameCard(self.elements[key]):
                #print 'eid=%s\noldElement=\n%snewElement=\n%s' %(key,self.elements[key],elem)
                assert elem.eid not in self.elements,'eid=%s\noldElement=\n%snewElement=\n%s' %(elem.eid,self.elements[elem.eid],elem)
        else:
            assert key>0,'eid=%s elem=%s' %(key,elem)
            self.elements[key] = elem
        ###

    def addDamperElement(self,elem,allowOverwrites=False):
        """@warning can dampers have the same ID as a standard element?"""
        return self.addElement(elem,allowOverwrites)
        key = elem.eid
        if key in self.damperElements and allowOverwrites==False:
            if not elem.isSameCard(self.damperElements[key]):
                #print 'eid=%s\noldElement=\n%snewElement=\n%s' %(key,self.elements[key],elem)
                assert elem.eid not in self.damperElements,'eid=%s\noldDamperElement=\n%snewDamperElement=\n%s' %(elem.eid,self.damperElements[elem.eid],elem)
        else:
            assert key>0,'eid=%s elem=%s' %(key,elem)
            self.damperElements[key] = elem
        ###

    def addRigidElement(self,elem,allowOverwrites=False):
        key = elem.eid
        if key in self.rigidElements and not allowOverwrites:
            print 'eid=%s\noldElement=\n%snewElement=\n%s' %(key,self.rigidElements[key],elem)
            #assert elem.eid not in self.rigidElements,'eid=%s\noldElement=\n%snewElement=\n%s' %(elem.eid,self.rigidElements[elem.eid],elem)
            pass
        assert key>0,'eid=%s elem=%s' %(key,elem)
        self.rigidElements[key] = elem

    def addThermalElement(self,elem):
        """same as addElement at the moment..."""
        self.addElement(elem)
        #assert elem.eid not in self.elements
        #assert elem.eid>0
        #self.elements[elem.eid] = elem

    def addDEQATN(self,deqatn,allowOverwrites=False):
        key = deqatn.eqID
        #if not allowOverwrites:
        #    assert prop.pid not in self.properties,'pid=%s oldProperty=\n%snewProperty=\n%s' %(prop.pid,self.properties[prop.pid],prop)
        assert key>0,'deqatn\n%s' %(key,deqatn)
        self.dequations[key] = deqatn

    def addProperty(self,prop,allowOverwrites=False):   
        key = prop.pid
        self.properties[key] = prop  ## temporary
        return                       ## temporary

        if not allowOverwrites:
            assert key not in self.properties,'pid=%s oldProperty=\n%snewProperty=\n%s' %(key,self.properties[key],prop)
        assert key>0,'property=\n%s' %(key,prop)
        self.properties[key] = prop

    def addMaterial(self,material,allowOverwrites=False):
        """
        only for adding structural materials
        @deprecated this method will be renamed in v0.3 to addStructuralMaterial.
        """
        self.addStructuralMaterial(material,allowOverwrites)

    def addStructuralMaterial(self,material,allowOverwrites=False):
        key = material.mid
        if key in self.materials and allowOverwrites==False:
            if not material.isSameCard(self.materials[key]):
                assert key not in self.materials,'mid=%s\noldMaterial=\n%snewMaterial=\n%s' %(key,self.materials[key],material)
        else:
            assert key>0,'mid=%s material=\n%s' %(key,material)
            self.materials[key] = material

    def addThermalMaterial(self,material,allowOverwrites=False):
        key = material.mid
        if key in self.thermalMaterials and allowOverwrites==False:
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
        if key in self.thermalMaterials and allowOverwrites==False:
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

        if constraint.conid in self.mpcadds:
            raise Exception('must have unique MPCADD IDs')
        self.mpcadds[constraint.conid] = constraint

    def addConstraint_MPC(self,constraint):
        self.mpcObject.append(constraint)
        if constraint.conid in self.mpcs:
            self.mpcs[constraint.conid].append(constraint)
        else:
            self.mpcs[constraint.conid] = [constraint]
        

    def addConstraint_SPCADD(self,constraint):
        self.spcObject.add(constraint)
        if constraint.conid in self.spcadds:
            raise Exception('must have unique SPCADD IDs')
        self.spcadds[constraint.conid] = constraint

    def addConstraint_SPC(self,constraint):
        self.spcObject.append(constraint)

        if constraint.conid in self.spcs:
            self.spcs[constraint.conid].append(constraint)
        else:
            self.spcs[constraint.conid] = [constraint]
        #key = constraint.conid
        #if self.constraints.has_key(key):
        #    self.constraints[key].append(constraint)
        #else:
        #    self.constraints[key] = [constraint]

    def addConstraint(self,constraint):
        #self.spcObject.append(constraint)
        key = constraint.conid

        if constraint.conid in self.spcs:
            self.spcs[constraint.conid].append(constraint)
        else:
            self.spcs[constraint.conid] = [constraint]
        #key = constraint.conid

        assert key>0
        if self.constraints.has_key(key):
            self.constraints[key].append(constraint)
        else:
            self.constraints[key] = [constraint]

    def addSUPORT(self,suport):
        self.suports.append(suport)

    def addDArea(self,darea,allowOverwrites=False):
        key = (darea.sid,darea.p)
        if key in self.dareas and allowOverwrites==False:
            if not darea.isSameCard(self.dareas[key]):
                assert key not in self.dareas,'\ndarea=\n%s oldDArea=\n%s' %(darea,self.dareas[key])
        else:
            assert darea.sid>0
            self.dareas[key] = darea
        ###

    def addAero(self,aero):
        key = aero.acsid
        assert key not in self.aero,'\naero=\n%s oldAERO=\n%s' %(aero,self.aero[key])
        assert key>=0
        self.aero[key] = aero

    def addAeros(self,aero):
        key = aero.acsid
        assert key not in self.aeros,'\naeros=\n%s oldAEROS=\n%s' %(aero,self.aeros[key])
        assert key>=0
        self.aeros[key] = aero

    def addAESurf(self,aesurf):
        key = aesurf.aesid
        assert key not in self.aesurfs,'\naesurf=\n%s oldAESURF=\n%s' %(aesurf,self.aesurfs[key])
        assert key>=0
        self.aesurfs[key] = aesurf

    def addAEList(self,aelist):
        key = aelist.sid
        assert key not in self.aelists,'\naelist=\n%s oldAELIST=\n%s' %(aelist,self.aelists[key])
        assert key>=0
        self.aelists[key] = aelist

    def addAELink(self,aelink):
        key = aelink.id
        assert key>=0
        if key not in self.aelinks:
            self.aelinks[key] = []
        self.aelinks[key].append(aelink)
        #assert key not in self.aestats,'\naestat=%s oldAESTAT=\n%s' %(aestat,self.aestats[key])

    def addAEParam(self,aeparam):
        key = aeparam.id
        assert key not in self.aeparams,'\naeparam=\n%s oldAESTAT=\n%s' %(aeparam,self.aeparams[key])
        assert key>=0
        self.aeparams[key] = aeparam

    def addAEStat(self,aestat):
        key = aestat.id
        assert key not in self.aestats,'\naestat=\n%s oldAESTAT=\n%s' %(aestat,self.aestats[key])
        assert key>=0
        self.aestats[key] = aestat

    def addCAero(self,caero):
        key = caero.eid
        assert key not in self.caeros,'\ncaero=\n|%s| oldCAERO=\n|%s|' %(caero,self.caeros[key])
        assert key>0
        self.caeros[key] = caero

    def addPAero(self,paero):
        key = paero.pid
        assert key not in self.paeros,'\npaero=\n|%s| oldPAERO=\n|%s|' %(paero,self.paeros[key])
        assert key>0,'paero.pid = |%s|' %(key)
        self.paeros[key] = paero

    def addGrav(self,grav):
        assert grav.sid not in self.gravs
        assert grav.sid>0
        self.gravs[grav.sid] = grav

    def addSpline(self,spline):
        assert spline.eid not in self.splines
        assert spline.eid>0
        self.splines[spline.eid] = spline

    def addGust(self,gust):
        key = gust.sid
        assert key not in self.gusts
        assert key>0
        self.gusts[key] = gust

    def addTrim(self,trim,allowOverwrites=False):
        key = trim.sid
        if not allowOverwrites:
            assert key not in self.trims,'trim=%s oldTrim=\n%snewProperty=\n%s' %(key,self.trims[key],trim)
        assert key>0,'trim=\n%s' %(key,trim)
        self.trims[key] = trim

    def addFlutter(self,flutter):
        key = flutter.sid
        assert key not in self.flutters
        assert key>0
        self.flutters[key] = flutter

    def addFLFACT(self,flfact):
        key = flfact.sid
        #assert key not in self.flfacts
        assert key>0
        self.flfacts[key] = flfact # set id...
        #print "added flfact...flflact =\n"+flfact

    def addDConstr(self,dconstr):
        key = (dconstr.oid,dconstr.rid)
        assert key not in self.dconstrs
        assert dconstr.oid>0
        assert dconstr.rid>0
        self.dconstrs[key] = dconstr

    def addDesvar(self,desvar):
        key = desvar.oid
        assert key not in self.desvars
        assert key>0
        self.desvars[key] = desvar

    def addDDVal(self,ddval):
        key = ddval.oid
        assert key not in self.ddvals
        assert key>0
        self.ddvals[key] = ddval

    def addDResp(self,dresp):
        key = dresp.oid
        assert key not in self.dresps
        assert key>0
        self.dresps[key] = dresp

    def addDvprel(self,dvprel):
        key = dvprel.oid
        assert key not in self.dvprels
        assert key>0
        self.dvprels[key] = dvprel

    def addNLParm(self,nlparm):
        key = nlparm.nid
        assert key not in self.nlparms
        assert key>0
        self.nlparms[key] = nlparm

    def addFREQ(self,freq):
        key = freq.sid
        assert key>0
        if key in self.frequencies:
            self.frequencies[key].addFrequencyObject(freq)
        else:
            self.frequencies[key] = freq
        #assert key not in self.frequencies,'\nfreq=\n%s oldFreq=\n%s' %(freq,self.frequencies[key])

    def addSet(self,setObj):
        key = setObj.sid
        assert key not in self.sets,'\nSET=\n%s oldSET=\n%s' %(setObj,self.sets[key])
        assert key>=0
        self.sets[key] = setObj

    def addSetSuper(self,setObj):
        key = setObj.seid
        assert key not in self.setsSuper,'\nSESET=\n%s oldSESET=\n%s' %(setObj,self.setsSuper[key])
        assert key>=0
        self.setsSuper[key] = setObj

    def addTable(self,table):
        key = table.tid
        assert key not in self.tables
        assert key>0
        self.tables[key] = table
