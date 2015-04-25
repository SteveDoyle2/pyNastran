# pylint: disable=E1101,C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)


class AddMethods(object):
    def __init__(self):
        pass

    def add_DMI(self, dmi, allowOverwrites=False):
        name = dmi.name
        self.dmis[name] = dmi

    def add_DMIG(self, dmig, allowOverwrites=False):
        name = dmig.name
        self.dmigs[name] = dmig

    def add_DMIJ(self, dmij, allowOverwrites=False):
        name = dmij.name
        self.dmijs[name] = dmij

    def add_DMIJI(self, dmiji, allowOverwrites=False):
        name = dmiji.name
        self.dmijis[name] = dmiji

    def add_DMIK(self, dmik, allowOverwrites=False):
        name = dmik.name
        self.dmiks[name] = dmik

    def add_PARAM(self, param, allowOverwrites=False):
        key = param.key
        if key in self.params and not allowOverwrites:
            if not param._is_same_card(self.params[key]):
                #assert param.key not in self.params,'key=%s param=%s oldPARAM=%s' %(key,param,self.params[key])
                self.log.warning('key=%s param=%s oldPARAM=%s'
                                 % (key, param, self.params[key]))
                self.params[key] = param
        else:
            self.params[key] = param

    def add_node(self, node, allowOverwrites=False):
        key = node.nid
        if key in self.nodes and not allowOverwrites:
            if not node._is_same_card(self.nodes[key]):
                print('nid=%s\noldNode=\n%snewNode=\n%s'
                      % (key, self.nodes[key], node))
                assert node.nid not in self.nodes, 'nid=%s\noldNode=\n%snewNode=\n%s' % (node.nid, self.nodes[key], node)
            else:
                #print('Node was duplicated...nid=%s; node=\n%s' % (key, node))
                pass
        else:
            assert key > 0, 'nid=%s node=%s' % (key, node)
            self.nodes[key] = node

    def add_SPOINT(self, spoint):
        if self.spoints is None:
            self.spoints = spoint
        else:
            self.spoints.addSPoints(spoint.spoints)

    def add_element(self, elem, allowOverwrites=False):
        key = elem.eid
        assert key > 0, 'eid=%s must be positive; elem=\n%s' % (key, elem)
        if key in self.elements and not allowOverwrites:
            if not elem._is_same_card(self.elements[key]):
                self._duplicate_elements.append(elem)
        else:
            self.elements[key] = elem

    def add_mass(self, mass, allowOverwrites=False):
        key = mass.eid
        assert key > 0, 'eid=%s must be positive; mass=\n%s' % (key, mass)
        if key in self.masses and not allowOverwrites:
            if not mass._is_same_card(self.masses[key]):
                self._duplicate_masses.append(mass)
        else:
            self.masses[key] = mass

    def add_damper(self, elem, allowOverwrites=False):
        """.. warning:: can dampers have the same ID as a standard element?"""
        return self.add_element(elem, allowOverwrites)

    def add_rigid_element(self, elem, allowOverwrites=False):
        key = elem.eid
        if key in self.rigidElements and not allowOverwrites:
            print('eid=%s\noldElement=\n%snewElement=\n%s' % (
                key, self.rigidElements[key], elem))
            assert elem.eid not in self.rigidElements, 'eid=%s\noldElement=\n%snewElement=\n%s' % (elem.eid, self.rigidElements[elem.eid], elem)
        assert key > 0, 'eid=%s elem=%s' % (key, elem)
        self.rigidElements[key] = elem

    def add_thermal_element(self, elem):
        """same as add_element at the moment..."""
        self.add_element(elem)

    def add_DEQATN(self, deqatn, allowOverwrites=False):
        key = deqatn.eqID
        #if not allowOverwrites:
        #    assert prop.pid not in self.properties, 'pid=%s oldProperty=\n%snewProperty=\n%s' % (prop.pid,self.properties[prop.pid], prop)
        assert key > 0, 'ID=%s deqatn\n%s' % (key, deqatn)
        self.dequations[key] = deqatn

    def add_property(self, prop, allowOverwrites=False):
        key = prop.pid
        assert key > 0, 'pid=%s prop=%s' % (key, prop)

        if key in self.properties and not allowOverwrites:
            if not prop._is_same_card(self.properties[key]):
                self._duplicate_properties.append(prop)
        else:
            self.properties[key] = prop

    def add_property_mass(self, prop, allowOverwrites=False):
        key = prop.pid
        if key in self.properties_mass and not allowOverwrites:
            if not prop._is_same_card(self.properties_mass[key]):
                #print('pid=%s\noldProperty=\n%snewProperty=\n%s' %(key,self.properties_mass[key],prop))
                assert key not in self.properties_mass, 'pid=%s oldProperty=\n%snewProperty=\n%s' % (key, self.properties_mass[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            self.properties_mass[key] = prop

    def add_PBUSHT(self, prop, allowOverwrites=False):
        key = prop.pid

        if key in self.pbusht and not allowOverwrites:
            if not prop._is_same_card(self.pbusht[key]):
                #print('pid=%s\noldProperty=\n%snewProperty=\n%s' %(key,self.pbusht[key],prop))
                assert key not in self.pbusht, 'pid=%s oldProperty=\n%snewProperty=\n%s' % (key, self.pbusht[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            self.pbusht[key] = prop

    def add_BCRPARA(self, card, allowOverwrites=False):
        self.bcrparas[card.crid] = card

    def add_BCTADD(self, card, allowOverwrites=False):
        self.bctadds[card.csid] = card

    def add_BCTPARA(self, card, allowOverwrites=False):
        self.bctparas[card.csid] = card

    def add_BCTSET(self, card, allowOverwrites=False):
        self.bctsets[card.csid] = card

    def add_BSURF(self, card, allowOverwrites=False):
        self.bsurf[card.sid] = card

    def add_BSURFS(self, card, allowOverwrites=False):
        self.bsurfs[card.id] = card

    def add_PDAMPT(self, prop, allowOverwrites=False):
        key = prop.pid
        if key in self.pdampt and not allowOverwrites:
            if not prop._is_same_card(self.pdampt[key]):
                #print('pid=%s\noldProperty=\n%snewProperty=\n%s' %(key,self.pdampt[key],prop))
                assert key not in self.pdampt, 'pid=%s oldProperty=\n%snewProperty=\n%s' % (key, self.pdampt[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            self.pdampt[key] = prop

    def add_PELAST(self, prop, allowOverwrites=False):
        key = prop.pid
        assert key > 0, 'pid=%s prop=%s' % (key, prop)
        if key in self.pelast and not allowOverwrites:
            if not prop._is_same_card(self.pelast[key]):
                #print('pid=%s\noldProperty=\n%snewProperty=\n%s' %(key,self.pdampt[key],prop))
                assert key not in self.pelast, 'pid=%s oldProperty=\n%snewProperty=\n%s' % (key, self.pelast[key], prop)
        else:
            self.pdampt[key] = prop

    def add_structural_material(self, material, allowOverwrites=False):
        key = material.mid
        assert key > 0, 'mid=%s material=\n%s' % (key, material)
        if key in self.materials and not allowOverwrites:
            if not material._is_same_card(self.materials[key]):
                self._duplicate_materials.append(material)
        else:
            self.materials[key] = material

    def add_thermal_material(self, material, allowOverwrites=False):
        key = material.mid
        assert key > 0, 'mid=%s material=\n%s' % (key, material)
        if key in self.thermalMaterials and not allowOverwrites:
            if not material._is_same_card(self.thermalMaterials[key]):
                self._duplicate_thermal_materials.append(material)
        else:
            self.thermalMaterials[key] = material

    def add_hyperelastic_material(self, material, allowOverwrites=False):
        key = material.mid
        assert key > 0, 'mid=%s material=\n%s' % (key, material)
        if key in self.hyperelasticMaterials and not allowOverwrites:
            if not material._is_same_card(self.hyperelasticMaterials[key]):
                assert key not in self.hyperelasticMaterials, 'mid=%s\noldMaterial=\n%snewMaterial=\n%s' % (key, self.hyperelasticMaterials[key], material)
        else:
            self.hyperelasticMaterials[key] = material

    def add_material_dependence(self, material, allowOverwrites=False):
        Type = material.type
        key = material.mid
        mapper = {
            'MATS1' : self.MATS1,
            'MATS3' : self.MATS3,
            'MATS8' : self.MATS8,

            'MATT1' : self.MATT1,
            'MATT2' : self.MATT2,
            'MATT3' : self.MATT3,
            'MATT4' : self.MATT4,
            'MATT5' : self.MATT5,
            'MATT8' : self.MATT8,
            'MATT9' : self.MATT9,
        }
        slot = mapper[Type]
        if key in slot and not allowOverwrites:
            if not material._is_same_card(slot[key]):
                assert key not in slot, 'mid=%s Type=%r\noldMaterialDep=\n%snewMaterialDep=\n%s' % (key, Type, slot[key], material)
        else:
            assert key > 0, 'mid=%s material=\n%s' % (key, material)
            slot[key] = material

    def add_creep_material(self, material, allowOverwrites=False):
        """
        .. note:: May be removed in the future.  Are CREEP cards materials?
                  They have an MID, but reference structural materials.
        """
        key = material.mid
        if key in self.thermalMaterials and not allowOverwrites:
            if not material._is_same_card(self.creepMaterials[key]):
                assert key not in self.creepMaterials, 'mid=%s\noldMaterial=\n%snewMaterial=\n%s' % (key, self.creepMaterials[key], material)
        else:
            assert key > 0, 'mid=%s material=\n%s' % (key, material)
            self.creepMaterials[key] = material

    def add_coord(self, coord, allowOverwrites=False):
        key = coord.cid
        assert coord.cid > -1, 'cid=%s coord=\n%s' % (key, coord)
        if key in self.coords:
            #if not allowOverwrites:
            if not coord._is_same_card(self.coords[key]):
                self._duplicate_coords.append(coord)
        else:
            self.coords[key] = coord

    def add_load(self, load):
        key = load.sid
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]

    def add_dload(self, dload):
        key = dload.sid
        if key in self.dloads:
            self.dloads[key].append(dload)
        else:
            self.dloads[key] = [dload]

    def add_dload_entry(self, dload):
        print(dload)
        key = dload.sid
        if key in self.dload_entries:
            self.dload_entries[key].append(dload)
        else:
            self.dload_entries[key] = [dload]

    def add_LSEQ(self, load):
        key = load.sid
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]

    def add_thermal_load(self, load):  # same function at the moment...
        key = load.sid
        assert key > 0
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]

    def add_PHBDY(self, prop):
        assert prop.pid > 0
        assert prop.pid not in self.phbdys
        self.phbdys[prop.pid] = prop

    def add_convection_property(self, prop):
        assert prop.pconid > 0
        assert prop.pconid not in self.convectionProperties
        self.convectionProperties[prop.pconid] = prop

    def add_thermal_BC(self, bc, key):
        assert key > 0
        if key in self.bcs:
            self.bcs[key].append(bc)
        else:
            self.bcs[key] = [bc]

    def add_constraint_MPCADD(self, constraint):
        if constraint.conid in self.mpcadds:
            raise RuntimeError('must have unique MPCADD IDs')
        self.mpcadds[constraint.conid] = constraint

    def add_constraint_MPC(self, constraint):
        if constraint.conid in self.mpcs:
            self.mpcs[constraint.conid].append(constraint)
        else:
            self.mpcs[constraint.conid] = [constraint]

    def add_constraint_SPCADD(self, constraint):
        if constraint.conid in self.spcadds:
            raise RuntimeError('must have unique SPCADD IDs')
        self.spcadds[constraint.conid] = constraint

    def add_constraint_SPC(self, constraint):
        if constraint.conid in self.spcs:
            self.spcs[constraint.conid].append(constraint)
        else:
            self.spcs[constraint.conid] = [constraint]

    def add_constraint(self, constraint):
        key = constraint.conid

        if constraint.conid in self.spcs:
            self.spcs[key].append(constraint)
        else:
            self.spcs[key] = [constraint]

    def add_suport(self, suport):
        self.suports.append(suport)

    def add_DAREA(self, darea, allowOverwrites=False):
        key = (darea.sid, darea.p)
        if key in self.dareas and not allowOverwrites:
            if not darea._is_same_card(self.dareas[key]):
                assert key not in self.dareas, '\ndarea=\n%s oldDArea=\n%s' % (darea, self.dareas[key])
        else:
            assert darea.sid > 0
            self.dareas[key] = darea

    def add_AERO(self, aero):
        key = aero.acsid
        assert key not in self.aero, '\naero=\n%s oldAERO=\n%s' % (aero, self.aero[key])
        assert key >= 0
        self.aero[key] = aero

    def add_AEROS(self, aero):
        key = aero.acsid
        assert key not in self.aeros, '\naeros=\n%s oldAEROS=\n%s' % (aero, self.aeros[key])
        assert key >= 0
        self.aeros[key] = aero

    def add_AEFACT(self, aefact, allowOverwrites=False):
        key = aefact.sid
        if key in self.aefacts and not allowOverwrites:
            if not aefact._is_same_card(self.aefacts[key]):
                assert key not in self.aefacts, 'sid=%s\noldAEFACT=\n%snewAEFACT=\n%s' % (key, self.aefacts[key], aefact)
        else:
            assert key > 0, 'sid=%s method=\n%s' % (key, aefact)
            self.aefacts[key] = aefact

    def add_AELIST(self, aelist):
        key = aelist.sid
        assert key not in self.aelists, '\naelist=\n%s oldAELIST=\n%s' % (aelist, self.aelists[key])
        assert key >= 0
        self.aelists[key] = aelist

    def add_AELINK(self, aelink):
        key = aelink.id
        assert key >= 0
        if key not in self.aelinks:
            self.aelinks[key] = []
        self.aelinks[key].append(aelink)
        #assert key not in self.aestats,'\naestat=%s oldAESTAT=\n%s' %(aestat,self.aestats[key])

    def add_AEPARM(self, aeparam):
        key = aeparam.id
        assert key not in self.aeparams, '\naeparam=\n%s oldAEPARM=\n%s' % (aeparam, self.aeparams[key])
        assert key >= 0
        self.aeparams[key] = aeparam

    def add_AESTAT(self, aestat):
        key = aestat.id
        assert key not in self.aestats, '\naestat=\n%s oldAESTAT=\n%s' % (
            aestat, self.aestats[key])
        assert key >= 0
        self.aestats[key] = aestat

    def add_AESURF(self, aesurf):
        key = aesurf.aesid
        assert key not in self.aesurfs, '\naesurf=\n%s oldAESURF=\n%s' % (
            aesurf, self.aesurfs[key])
        assert key >= 0
        self.aesurfs[key] = aesurf

    def add_CAERO(self, caero):
        key = caero.eid
        assert key not in self.caeros, '\ncaero=\n|%s| oldCAERO=\n|%s|' % (
            caero, self.caeros[key])
        assert key > 0
        self.caeros[key] = caero

    def add_PAERO(self, paero):
        key = paero.pid
        assert key not in self.paeros, '\npaero=\n|%s| oldPAERO=\n|%s|' % (
            paero, self.paeros[key])
        assert key > 0, 'paero.pid = %r' % (key)
        self.paeros[key] = paero

    def add_SPLINE(self, spline):
        assert spline.eid not in self.splines
        assert spline.eid > 0
        self.splines[spline.eid] = spline

    def add_GUST(self, gust):
        key = gust.sid
        assert key not in self.gusts
        assert key > 0
        self.gusts[key] = gust

    def add_TRIM(self, trim, allowOverwrites=False):
        key = trim.sid
        if not allowOverwrites:
            assert key not in self.trims, 'trim=%s oldTrim=\n%snewProperty=\n%s' % (key, self.trims[key], trim)
        assert key > 0, 'key=%r trim=\n%s' % (key, trim)
        self.trims[key] = trim

    def add_FLUTTER(self, flutter):
        key = flutter.sid
        assert key not in self.flutters
        assert key > 0
        self.flutters[key] = flutter

    def add_FLFACT(self, flfact):
        key = flfact.sid
        #assert key not in self.flfacts
        assert key > 0
        self.flfacts[key] = flfact  # set id...

    def add_DCONSTR(self, dconstr):
        key = (dconstr.oid, dconstr.rid)
        assert key not in self.dconstrs
        assert dconstr.oid > 0
        assert dconstr.rid > 0
        self.dconstrs[key] = dconstr

    def add_DESVAR(self, desvar):
        key = desvar.oid
        assert key not in self.desvars
        assert key > 0
        self.desvars[key] = desvar

    def add_DDVAL(self, ddval):
        key = ddval.oid
        assert key not in self.ddvals
        assert key > 0
        self.ddvals[key] = ddval

    def add_DLINK(self, dlink):
        key = dlink.oid
        assert key not in self.dlinks
        assert key > 0
        self.dlinks[key] = dlink

    def add_DRESP(self, dresp):
        key = dresp.oid
        assert key not in self.dresps
        assert key > 0
        self.dresps[key] = dresp

    def add_DVMREL(self, dvmrel):
        key = dvmrel.oid
        assert key not in self.dvmrels
        assert key > 0
        self.dvmrels[key] = dvmrel

    def add_DVPREL(self, dvprel):
        key = dvprel.oid
        assert key not in self.dvprels
        assert key > 0
        self.dvprels[key] = dvprel

    def add_NLPARM(self, nlparm):
        key = nlparm.nid
        assert key not in self.nlparms
        assert key > 0
        self.nlparms[key] = nlparm

    def add_NLPCI(self, nlpci):
        key = nlpci.nlparm_id
        assert key not in self.nlpcis
        assert key > 0
        self.nlpcis[key] = nlpci

    def add_TSTEP(self, tstep, allowOverwrites=False):
        key = tstep.sid
        if key in self.tsteps and not allowOverwrites:
            if not tstep._is_same_card(self.tsteps[key]):
                assert key not in self.tsteps, 'sid=%s\noldTSTEP=\n%snewTSTEP=\n%s' % (key, self.tsteps[key], tstep)
        else:
            assert key > 0, 'sid=%s tstep=\n%s' % (key, tstep)
            self.tsteps[key] = tstep

    def add_TSTEPNL(self, tstepnl, allowOverwrites=False):
        key = tstepnl.sid
        if key in self.tstepnls and not allowOverwrites:
            if not tstepnl._is_same_card(self.tstepnls[key]):
                assert key not in self.tstepnls, 'sid=%s\noldTSTEPNL=\n%snewTSTEPNL=\n%s' % (key, self.tstepnls[key], tstepnl)
        else:
            assert key > 0, 'sid=%s tstepnl=\n%s' % (key, tstepnl)
            self.tstepnls[key] = tstepnl

    def add_FREQ(self, freq):
        key = freq.sid
        assert key > 0
        if key in self.frequencies:
            self.frequencies[key].add_frequency_object(freq)
        else:
            self.frequencies[key] = freq

    def add_SET(self, set_obj):
        key = set_obj.sid
        assert key >= 0
        if key in self.sets:
            self.sets[key].add_set(set_obj)
        else:
            self.sets[key] = set_obj

    def add_ASET(self, set_obj):
        self.asets.append(set_obj)

    def add_BSET(self, set_obj):
        self.bsets.append(set_obj)

    def add_CSET(self, set_obj):
        self.csets.append(set_obj)

    def add_QSET(self, set_obj):
        self.qsets.append(set_obj)

    def add_SESET(self, set_obj):
        key = set_obj.seid
        assert key >= 0
        if key in self.setsSuper:
            old_set = self.setsSuper[key]
            set_obj.add_SESET_Object(old_set)
        self.setsSuper[key] = set_obj

    def add_table(self, table):
        key = table.tid
        assert key not in self.tables, '\nTable=\n%s oldTable=\n%s' % (
            table, self.tables[key])
        assert key > 0
        self.tables[key] = table

    def add_random_table(self, table):
        key = table.tid
        assert key not in self.randomTables, '\nTable=\n%s oldTable=\n%s' % (
            table, self.randomTables[key])
        assert key > 0
        self.randomTables[key] = table

    def add_method(self, method, allowOverwrites=False):
        key = method.sid
        if key in self.methods and not allowOverwrites:
            if not method._is_same_card(self.methods[key]):
                assert key not in self.methods, 'sid=%s\noldMethod=\n%snewMethod=\n%s' % (key, self.methods[key], method)
        else:
            assert key > 0, 'sid=%s method=\n%s' % (key, method)
            self.methods[key] = method

    def add_cmethod(self, cMethod, allowOverwrites=False):
        key = cMethod.sid
        if key in self.cMethods and not allowOverwrites:
            if not cMethod._is_same_card(self.cMethods[key]):
                assert key not in self.cMethods, 'sid=%s\noldCMethod=\n%snewCMethod=\n%s' % (key, self.cMethods[key], cMethod)
        else:
            assert key > 0, 'sid=%s cMethod=\n%s' % (key, cMethod)
            self.cMethods[key] = cMethod

    def add_MKAERO(self, mkaero):
        self.mkaeros.append(mkaero)