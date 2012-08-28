# pylint: disable=E1101,C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)


class AddMethods(object):
    def __init__(self):
        pass

    def addDMI(self, dmi, allowOverwrites=False):
        name = dmi.name
        #if key in self.params and not allowOverwrites:
            #if not param.isSameCard(self.params[key]):
                #assert param.key not in self.params,'key=%s param=%s oldPARAM=%s' %(key,param,self.params[key])
                #self.log.warning('key=%s param=%s oldPARAM=%s' %(key,param,self.params[key]))
                #self.params[key] = param
        #else:
        self.dmis[name] = dmi
        ###

    def addDMIG(self, dmig, allowOverwrites=False):
        name = dmig.name
        #if key in self.params and not allowOverwrites:
            #if not param.isSameCard(self.params[key]):
                #assert param.key not in self.params,'key=%s param=%s oldPARAM=%s' %(key,param,self.params[key])
                #self.log.warning('key=%s param=%s oldPARAM=%s' %(key,param,self.params[key]))
                #self.params[key] = param
        #else:
        self.dmigs[name] = dmig
        ###

    def addDMIJ(self, dmij, allowOverwrites=False):
        name = dmij.name
        #if key in self.params and not allowOverwrites:
            #if not param.isSameCard(self.params[key]):
                #assert param.key not in self.params,'key=%s param=%s oldPARAM=%s' %(key,param,self.params[key])
                #self.log.warning('key=%s param=%s oldPARAM=%s' %(key,param,self.params[key]))
                #self.params[key] = param
        #else:
        self.dmijs[name] = dmij
        ###

    def addDMIJI(self, dmiji, allowOverwrites=False):
        name = dmiji.name
        #if key in self.params and not allowOverwrites:
            #if not param.isSameCard(self.params[key]):
                #assert param.key not in self.params,'key=%s param=%s oldPARAM=%s' %(key,param,self.params[key])
                #self.log.warning('key=%s param=%s oldPARAM=%s' %(key,param,self.params[key]))
                #self.params[key] = param
        #else:
        self.dmijis[name] = dmiji
        ###

    def addDMIK(self, dmik, allowOverwrites=False):
        name = dmik.name
        #if key in self.params and not allowOverwrites:
            #if not param.isSameCard(self.params[key]):
                #assert param.key not in self.params,'key=%s param=%s oldPARAM=%s' %(key,param,self.params[key])
                #self.log.warning('key=%s param=%s oldPARAM=%s' %(key,param,self.params[key]))
                #self.params[key] = param
        #else:
        self.dmiks[name] = dmik
        ###

    def addParam(self, param, allowOverwrites=False):
        key = param.key
        if key in self.params and not allowOverwrites:
            if not param.isSameCard(self.params[key]):
                #assert param.key not in self.params,'key=%s param=%s oldPARAM=%s' %(key,param,self.params[key])
                self.log.warning('key=%s param=%s oldPARAM=%s' %
                                (key, param, self.params[key]))
                self.params[key] = param
        else:
            self.params[key] = param
        ###

    def addNode(self, node, allowOverwrites=False):
        #print node
        #assert node.nid not in self.nodes,'nid=%s\noldNode=\n%snewNode=\n%s' %(node.nid,self.nodes[node.nid],node)  ## @todo enable before release...
        #assert node.nid>0,'nid=%s node=\n%s' %(node.nid,node)
        #self.nodes[key] = node

        key = node.nid
        if key in self.nodes and not allowOverwrites:
            if not node.isSameCard(self.nodes[key]):
                print('nid=%s\noldNode=\n%snewNode=\n%s' %
                     (key, self.nodes[key], node))
                assert node.nid not in self.nodes, 'nid=%s\noldNode=\n%snewNode=\n%s' % (node.nid, self.nodes[key], node)
            else:
                #print 'Node was duplicated...nid=%s\nnode=\n%s' %(key,node)
                pass
        else:
            assert key > 0, 'nid=%s node=%s' % (key, node)
            self.nodes[key] = node
        ###

    def addSPoint(self, spoint):
        if self.spoints is None:
            self.spoints = spoint
        else:
            self.spoints.addSPoints(spoint.spoints)
        ###

    def addElement(self, elem, allowOverwrites=False):
        key = elem.eid
        #self.elements[key] = elem  ## @todo temporary
        #return                     ## @todo temporary

        if key in self.elements and not allowOverwrites:
            if not elem.isSameCard(self.elements[key]):
                elem.isSameCard(self.elements[key], debug=True)
                #print 'eid=%s\noldElement=\n%snewElement=\n%s' %(key,self.elements[key],elem)
                assert elem.eid not in self.elements, 'eid=%s\noldElement=\n%snewElement=\n%s' % (elem.eid, self.elements[elem.eid], elem)
        else:
            assert key > 0, 'eid=%s elem=%s' % (key, elem)
            self.elements[key] = elem
        ###

    def addMassElement(self, elem, allowOverwrites=False):
        key = elem.eid
        #self.massElements[key] = elem  ## @todo temporary
        #return                         ## @todo temporary

        if key in self.massElements and not allowOverwrites:
            if not elem.isSameCard(self.massElements[key]):
                elem.isSameCard(self.elements[key], debug=True)
                #print 'eid=%s\noldElement=\n%snewElement=\n%s' %(key,self.elements[key],elem)
                assert elem.eid not in self.massElements, 'eid=%s\noldElement=\n%snewElement=\n%s' % (elem.eid, self.massElements[elem.eid], elem)
        else:
            assert key > 0, 'eid=%s elem=%s' % (key, elem)
            self.massElements[key] = elem
        ###

    def addDamperElement(self, elem, allowOverwrites=False):
        """@warning can dampers have the same ID as a standard element?"""
        return self.addElement(elem, allowOverwrites)
        key = elem.eid
        if key in self.damperElements and not allowOverwrites:
            if not elem.isSameCard(self.damperElements[key]):
                #print 'eid=%s\noldElement=\n%snewElement=\n%s' %(key,self.elements[key],elem)
                assert elem.eid not in self.damperElements, 'eid=%s\noldDamperElement=\n%snewDamperElement=\n%s' % (elem.eid, self.damperElements[elem.eid], elem)
        else:
            assert key > 0, 'eid=%s elem=%s' % (key, elem)
            self.damperElements[key] = elem

    def addRigidElement(self, elem, allowOverwrites=False):
        key = elem.eid
        if key in self.rigidElements and not allowOverwrites:
            print('eid=%s\noldElement=\n%snewElement=\n%s' % (
                key, self.rigidElements[key], elem))
            #assert elem.eid not in self.rigidElements,'eid=%s\noldElement=\n%snewElement=\n%s' %(elem.eid,self.rigidElements[elem.eid],elem)
            pass
        assert key > 0, 'eid=%s elem=%s' % (key, elem)
        self.rigidElements[key] = elem

    def addThermalElement(self, elem):
        """same as addElement at the moment..."""
        self.addElement(elem)
        #assert elem.eid not in self.elements
        #assert elem.eid>0
        #self.elements[elem.eid] = elem

    def addDEQATN(self, deqatn, allowOverwrites=False):
        key = deqatn.eqID
        #if not allowOverwrites:
        #    assert prop.pid not in self.properties,'pid=%s oldProperty=\n%snewProperty=\n%s' %(prop.pid,self.properties[prop.pid],prop)
        assert key > 0, 'ID=%s deqatn\n%s' % (key, deqatn)
        self.dequations[key] = deqatn

    def addProperty(self, prop, allowOverwrites=False):
        key = prop.pid
        #self.properties[key] = prop  ## @todo temporary
        #return                       ## @todo temporary

        if key in self.properties and not allowOverwrites:
            if not prop.isSameCard(self.properties[key]):
                #print 'pid=%s\noldProperty=\n%snewProperty=\n%s' %(key,self.properties[key],prop)
                assert key not in self.properties, 'pid=%s oldProperty=\n%snewProperty=\n%s' % (key, self.properties[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            self.properties[key] = prop
        ###

    def addMaterial(self, material, allowOverwrites=False):
        """
        only for adding structural materials
        @deprecated this method will be renamed in v0.3 to addStructuralMaterial.
        """
        self.addStructuralMaterial(material, allowOverwrites)

    def addStructuralMaterial(self, material, allowOverwrites=False):
        key = material.mid
        if key in self.materials and not allowOverwrites:
            if not material.isSameCard(self.materials[key]):
                assert key not in self.materials, 'mid=%s\noldMaterial=\n%snewMaterial=\n%s' % (key, self.materials[key], material)
        else:
            assert key > 0, 'mid=%s material=\n%s' % (key, material)
            self.materials[key] = material

    def addThermalMaterial(self, material, allowOverwrites=False):
        key = material.mid
        if key in self.thermalMaterials and not allowOverwrites:
            if not material.isSameCard(self.thermalMaterials[key]):
                assert key not in self.thermalMaterials, 'mid=%s\noldMaterial=\n%snewMaterial=\n%s' % (key, self.thermalMaterials[key], material)
        else:
            assert key > 0, 'mid=%s material=\n%s' % (key, material)
            self.thermalMaterials[key] = material

    def addMaterialDependence(self, material, allowOverwrites=False):
        key = material.mid
        if key in self.materialDeps and not allowOverwrites:
            if not material.isSameCard(self.materialDeps[key]):
                assert key not in self.materialDeps, 'mid=%s\noldMaterialDep=\n%snewMaterialDep=\n%s' % (key, self.materialDeps[key], material)
        else:
            assert key > 0, 'mid=%s material=\n%s' % (key, material)
            self.materialDeps[key] = material

    def addCreepMaterial(self, material, allowOverwrites=False):
        """
        Method addCreepMaterial:
        @note
         May be removed in the future.  Are CREEP cards materials? 
         They have an MID, but reference structural materials.
        """
        key = material.mid
        if key in self.thermalMaterials and not allowOverwrites:
            if not material.isSameCard(self.creepMaterials[key]):
                assert key not in self.creepMaterials, 'mid=%s\noldMaterial=\n%snewMaterial=\n%s' % (key, self.creepMaterials[key], material)
        else:
            assert key > 0, 'mid=%s material=\n%s' % (key, material)
            self.creepMaterials[key] = material

    def addCoord(self, coord, allowOverwrites=False):
        key = coord.cid
        if not allowOverwrites:
            assert key not in self.coords, 'cid=%s\noldElement=\n%snewElement=\n%s' % (key, self.coords[key], coord)
        assert coord.cid > -1, 'cid=%s coord=\n%s' % (key, coord)
        self.coords[key] = coord

    def addLoad(self, load):
        key = load.sid
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]

    def addLSeq(self, load):
        key = load.sid
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]

    def addThermalLoad(self, load):  # same function at the moment...
        key = load.sid
        assert key > 0
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]

    def addPHBDY(self, prop):
        assert prop.pid > 0
        assert prop.pid not in self.phbdys
        self.phbdys[prop.pid] = prop

    def addConvectionProperty(self, prop):
        assert prop.pconid > 0
        assert prop.pconid not in self.convectionProperties
        self.convectionProperties[prop.pconid] = prop

    #def addThermalProperty(self,prop):
    #    assert prop.pconid not in self.thermalProperties
    #    self.thermalProperties[prop.pconid] = prop

    def addThermalBC(self, bc, key):
        assert key > 0
        if key in self.bcs:
            self.bcs[key].append(bc)
        else:
            self.bcs[key] = [bc]

    def addConstraint_MPCADD(self, constraint):
        #self.mpcObject.add(constraint)
        if constraint.conid in self.mpcadds:
            raise RuntimeError('must have unique MPCADD IDs')
        self.mpcadds[constraint.conid] = constraint

    def addConstraint_MPC(self, constraint):
        #self.mpcObject.append(constraint)
        if constraint.conid in self.mpcs:
            self.mpcs[constraint.conid].append(constraint)
        else:
            self.mpcs[constraint.conid] = [constraint]

    def addConstraint_SPCADD(self, constraint):
        #self.spcObject.add(constraint)
        if constraint.conid in self.spcadds:
            raise RuntimeError('must have unique SPCADD IDs')
        self.spcadds[constraint.conid] = constraint

    def addConstraint_SPC(self, constraint):
        #self.spcObject.append(constraint)
        if constraint.conid in self.spcs:
            self.spcs[constraint.conid].append(constraint)
        else:
            self.spcs[constraint.conid] = [constraint]
        #key = constraint.conid
        #if self.constraints.has_key(key):
        #    self.constraints[key].append(constraint)
        #else:
        #    self.constraints[key] = [constraint]

    def addConstraint(self, constraint):
        #self.spcObject.append(constraint)
        key = constraint.conid

        if constraint.conid in self.spcs:
            self.spcs[key].append(constraint)
        else:
            self.spcs[key] = [constraint]

        #assert key > 0
        #if self.constraints.has_key(key):
        #    self.constraints[key].append(constraint)
        #else:
        #    self.constraints[key] = [constraint]

    def addSuport(self, suport):
        self.suports.append(suport)

    def addDArea(self, darea, allowOverwrites=False):
        key = (darea.sid, darea.p)
        if key in self.dareas and not allowOverwrites:
            if not darea.isSameCard(self.dareas[key]):
                assert key not in self.dareas, '\ndarea=\n%s oldDArea=\n%s' % (
                    darea, self.dareas[key])
        else:
            assert darea.sid > 0
            self.dareas[key] = darea
        ###

    def addAero(self, aero):
        key = aero.acsid
        assert key not in self.aero, '\naero=\n%s oldAERO=\n%s' % (
            aero, self.aero[key])
        assert key >= 0
        self.aero[key] = aero

    def addAeros(self, aero):
        key = aero.acsid
        assert key not in self.aeros, '\naeros=\n%s oldAEROS=\n%s' % (
            aero, self.aeros[key])
        assert key >= 0
        self.aeros[key] = aero

    def addAEFact(self, aefact, allowOverwrites=False):
        key = aefact.sid
        if key in self.aefacts and not allowOverwrites:
            if not aefact.isSameCard(self.aefacts[key]):
                assert key not in self.aefacts, 'sid=%s\noldAEFACT=\n%snewAEFACT=\n%s' % (key, self.aefacts[key], aefact)
        else:
            assert key > 0, 'sid=%s method=\n%s' % (key, aefact)
            self.aefacts[key] = aefact

    def addAEList(self, aelist):
        key = aelist.sid
        assert key not in self.aelists, '\naelist=\n%s oldAELIST=\n%s' % (
            aelist, self.aelists[key])
        assert key >= 0
        self.aelists[key] = aelist

    def addAELink(self, aelink):
        key = aelink.id
        assert key >= 0
        if key not in self.aelinks:
            self.aelinks[key] = []
        self.aelinks[key].append(aelink)
        #assert key not in self.aestats,'\naestat=%s oldAESTAT=\n%s' %(aestat,self.aestats[key])

    def addAEParam(self, aeparam):
        key = aeparam.id
        assert key not in self.aeparams, '\naeparam=\n%s oldAESTAT=\n%s' % (
            aeparam, self.aeparams[key])
        assert key >= 0
        self.aeparams[key] = aeparam

    def addAEStat(self, aestat):
        key = aestat.id
        assert key not in self.aestats, '\naestat=\n%s oldAESTAT=\n%s' % (
            aestat, self.aestats[key])
        assert key >= 0
        self.aestats[key] = aestat

    def addAESurf(self, aesurf):
        key = aesurf.aesid
        assert key not in self.aesurfs, '\naesurf=\n%s oldAESURF=\n%s' % (
            aesurf, self.aesurfs[key])
        assert key >= 0
        self.aesurfs[key] = aesurf

    def addCAero(self, caero):
        key = caero.eid
        assert key not in self.caeros, '\ncaero=\n|%s| oldCAERO=\n|%s|' % (
            caero, self.caeros[key])
        assert key > 0
        self.caeros[key] = caero

    def addPAero(self, paero):
        key = paero.pid
        assert key not in self.paeros, '\npaero=\n|%s| oldPAERO=\n|%s|' % (
            paero, self.paeros[key])
        assert key > 0, 'paero.pid = |%s|' % (key)
        self.paeros[key] = paero

    def addSpline(self, spline):
        assert spline.eid not in self.splines
        assert spline.eid > 0
        self.splines[spline.eid] = spline

    def addGUST(self, gust):
        key = gust.sid
        assert key not in self.gusts
        assert key > 0
        self.gusts[key] = gust

    def addTrim(self, trim, allowOverwrites=False):
        key = trim.sid
        if not allowOverwrites:
            assert key not in self.trims, 'trim=%s oldTrim=\n%snewProperty=\n%s' % (key, self.trims[key], trim)
        assert key > 0, 'trim=\n%s' % (key, trim)
        self.trims[key] = trim

    def addFlutter(self, flutter):
        key = flutter.sid
        assert key not in self.flutters
        assert key > 0
        self.flutters[key] = flutter

    def addFLFACT(self, flfact):
        key = flfact.sid
        #assert key not in self.flfacts
        assert key > 0
        self.flfacts[key] = flfact  # set id...
        #print "added flfact...flflact =\n"+flfact

    def addDConstr(self, dconstr):
        key = (dconstr.oid, dconstr.rid)
        assert key not in self.dconstrs
        assert dconstr.oid > 0
        assert dconstr.rid > 0
        self.dconstrs[key] = dconstr

    def addDesvar(self, desvar):
        key = desvar.oid
        assert key not in self.desvars
        assert key > 0
        self.desvars[key] = desvar

    def addDDVal(self, ddval):
        key = ddval.oid
        assert key not in self.ddvals
        assert key > 0
        self.ddvals[key] = ddval

    def addDLink(self, dlink):
        key = dlink.oid
        assert key not in self.dlinks
        assert key > 0
        self.dlinks[key] = dlink

    def addDResp(self, dresp):
        key = dresp.oid
        assert key not in self.dresps
        assert key > 0
        self.dresps[key] = dresp

    def addDvmrel(self, dvmrel):
        key = dvmrel.oid
        assert key not in self.dvmrels
        assert key > 0
        self.dvmrels[key] = dvmrel

    def addDvprel(self, dvprel):
        key = dvprel.oid
        assert key not in self.dvprels
        assert key > 0
        self.dvprels[key] = dvprel

    def addNLParm(self, nlparm):
        key = nlparm.nid
        assert key not in self.nlparms
        assert key > 0
        self.nlparms[key] = nlparm

    def addTSTEP(self, tstep, allowOverwrites=False):
        key = tstep.sid
        if key in self.tsteps and not allowOverwrites:
            if not tstep.isSameCard(self.tsteps[key]):
                assert key not in self.tsteps, 'sid=%s\noldTSTEP=\n%snewTSTEP=\n%s' % (key, self.tsteps[key], tstep)
        else:
            assert key > 0, 'sid=%s tstep=\n%s' % (key, tstep)
            self.tsteps[key] = tstep

    def addTSTEPNL(self, tstepnl, allowOverwrites=False):
        key = tstepnl.sid
        if key in self.tsteps and not allowOverwrites:
            if not tstepnl.isSameCard(self.tsteps[key]):
                assert key not in self.tstepnls, 'sid=%s\noldTSTEPNL=\n%snewTSTEPNL=\n%s' % (key, self.tstepnls[key], tstepnl)
        else:
            assert key > 0, 'sid=%s tstepnl=\n%s' % (key, tstepnl)
            self.tstepnls[key] = tstepnl

    def addFREQ(self, freq):
        key = freq.sid
        assert key > 0
        if key in self.frequencies:
            self.frequencies[key].addFrequencyObject(freq)
        else:
            self.frequencies[key] = freq
        #assert key not in self.frequencies,'\nfreq=\n%s oldFreq=\n%s' %(freq,self.frequencies[key])

    def addSet(self, setObj):
        key = setObj.sid
        assert key not in self.sets, '\nSET=\n%s oldSET=\n%s' % (
            setObj, self.sets[key])
        assert key >= 0
        self.sets[key] = setObj

    def addASet(self, setObj):
        self.asets.append(setObj)

    def addBSet(self, setObj):
        self.bsets.append(setObj)

    def addCSet(self, setObj):
        self.csets.append(setObj)

    def addQSet(self, setObj):
        self.qsets.append(setObj)

    def addSetSuper(self, setObj):
        key = setObj.seid
        assert key not in self.setsSuper, '\nSESET=\n%s oldSESET=\n%s' % (
            setObj, self.setsSuper[key])
        assert key >= 0
        self.setsSuper[key] = setObj

    def addTable(self, table):
        key = table.tid
        assert key not in self.tables, '\nTable=\n%s oldTable=\n%s' % (
            table, self.tables[key])
        assert key > 0
        self.tables[key] = table

    def addRandomTable(self, table):
        key = table.tid
        assert key not in self.randomTables, '\nTable=\n%s oldTable=\n%s' % (
            table, self.randomTables[key])
        assert key > 0
        self.randomTables[key] = table

    def addMethod(self, method, allowOverwrites=False):
        key = method.sid
        if key in self.methods and not allowOverwrites:
            if not method.isSameCard(self.methods[key]):
                assert key not in self.methods, 'sid=%s\noldMethod=\n%snewMethod=\n%s' % (key, self.methods[key], method)
        else:
            assert key > 0, 'sid=%s method=\n%s' % (key, method)
            self.methods[key] = method

    def addCMethod(self, cMethod, allowOverwrites=False):
        key = cMethod.sid
        if key in self.cMethods and not allowOverwrites:
            if not cMethod.isSameCard(self.cMethods[key]):
                assert key not in self.cMethods, 'sid=%s\noldCMethod=\n%snewCMethod=\n%s' % (key, self.cMethods[key], cMethod)
        else:
            assert key > 0, 'sid=%s cMethod=\n%s' % (key, cMethod)
            self.cMethods[key] = cMethod

    def addMKAero(self, mkaero):
        self.mkaeros.append(mkaero)
