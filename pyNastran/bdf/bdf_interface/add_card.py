# pylint: disable=R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from pyNastran.bdf.bdf_interface.attributes import BDFAttributes


class AddMethods(BDFAttributes):
    def __init__(self):
        BDFAttributes.__init__(self)

    def add_DMI(self, dmi, allow_overwrites=False):
        name = dmi.name
        self.dmis[name] = dmi
        self._type_to_id_map[dmi.type].append(name)

    def add_DMIG(self, dmig, allow_overwrites=False):
        name = dmig.name
        self.dmigs[name] = dmig
        self._type_to_id_map[dmig.type].append(name)

    def add_DMIJ(self, dmij, allow_overwrites=False):
        name = dmij.name
        self.dmijs[name] = dmij
        self._type_to_id_map[dmij.type].append(name)

    def add_DMIJI(self, dmiji, allow_overwrites=False):
        name = dmiji.name
        self.dmijis[name] = dmiji
        self._type_to_id_map[dmiji.type].append(name)

    def add_DMIK(self, dmik, allow_overwrites=False):
        name = dmik.name
        self.dmiks[name] = dmik
        self._type_to_id_map[dmik.type].append(name)

    def add_PARAM(self, param, allow_overwrites=False):
        key = param.key
        if key in self.params and not allow_overwrites:
            if not param._is_same_card(self.params[key]):
                #assert param.key not in self.params,'key=%s param=%s oldPARAM=%s' %(key,param,self.params[key])
                self.log.warning('key=%s param=%s oldPARAM=%s'
                                 % (key, param, self.params[key]))
                self.params[key] = param
        else:
            self.params[key] = param
            self._type_to_id_map[param.type].append(key)

    def add_node(self, node, allow_overwrites=False):
        key = node.nid
        if key in self.nodes and not allow_overwrites:
            if not node._is_same_card(self.nodes[key]):
                assert node.nid not in self.nodes, 'nid=%s\noldNode=\n%snewNode=\n%s' % (node.nid, self.nodes[key], node)
            else:
                #print('Node was duplicated...nid=%s; node=\n%s' % (key, node))
                pass
        else:
            assert key > 0, 'nid=%s node=%s' % (key, node)
            self.nodes[key] = node
            self._type_to_id_map[node.type].append(key)

    def add_spoint(self, spoints):
        if self.spoints is None:
            self.spoints = spoints
        else:
            self.spoints.add_points(spoints.points)

    def add_epoint(self, epoints):
        if self.epoints is None:
            self.epoints = epoints
        else:
            self.epoints.add_points(epoints.points)

    def add_plotel(self, elem, allow_overwrites=False):
        key = elem.eid
        assert key > 0, 'eid=%s must be positive; elem=\n%s' % (key, elem)
        if not allow_overwrites:
            if key in self.elements:
                if elem._is_same_card(self.elements[key]):
                    self._duplicate_elements.append(elem)
                    if self._stop_on_duplicate_error:
                        self.pop_parse_errors()
            elif key in self.plotels:
                if not elem._is_same_card(self.plotels[key]):
                    assert elem.eid not in self.plotels, 'eid=%s\nold_element=\n%snew_element=\n%s' % (elem.eid, self.plotels[elem.eid], elem)
        self.plotels[key] = elem
        self._type_to_id_map[elem.type].append(key)

    def add_element(self, elem, allow_overwrites=False):
        key = elem.eid
        assert key > 0, 'eid=%s must be positive; elem=\n%s' % (key, elem)
        if key in self.elements and not allow_overwrites:
            if not elem._is_same_card(self.elements[key]):
                self._duplicate_elements.append(elem)
                if self._stop_on_duplicate_error:
                    self.pop_parse_errors()
        else:
            self.elements[key] = elem
            self._type_to_id_map[elem.type].append(key)

    def add_mass(self, mass, allow_overwrites=False):
        key = mass.eid
        assert key > 0, 'eid=%s must be positive; mass=\n%s' % (key, mass)
        if key in self.masses and not allow_overwrites:
            if not mass._is_same_card(self.masses[key]):
                self._duplicate_masses.append(mass)
        else:
            self.masses[key] = mass
            self._type_to_id_map[mass.type].append(key)

    def add_damper(self, elem, allow_overwrites=False):
        """.. warning:: can dampers have the same ID as a standard element?"""
        return self.add_element(elem, allow_overwrites)

    def add_rigid_element(self, elem, allow_overwrites=False):
        key = elem.eid
        assert key > 0, 'eid=%s elem=%s' % (key, elem)
        if key in self.rigid_elements and not allow_overwrites:
            assert elem.eid not in self.rigid_elements, 'eid=%s\noldElement=\n%snewElement=\n%s' % (elem.eid, self.rigid_elements[elem.eid], elem)
        self.rigid_elements[key] = elem
        self._type_to_id_map[elem.type].append(key)

    def add_thermal_element(self, elem):
        """same as add_element at the moment..."""
        self.add_element(elem)

    def add_DEQATN(self, deqatn, allow_overwrites=False):
        key = deqatn.equation_id
        assert key > 0, 'ID=%s deqatn\n%s' % (key, deqatn)
        if key in self.dequations and not allow_overwrites:
            if not deqatn.write_card() == self.dequations[key].write_card():
                assert key not in self.dequations, 'id=%s old_eq=\n%snew_eq=\n%s' % (
                    key, self.dequations[key], deqatn)
        self.dequations[key] = deqatn
        self._type_to_id_map[deqatn.type].append(key)

    def add_property(self, prop, allow_overwrites=False):
        key = prop.pid
        assert key > 0, 'pid=%s prop=%s' % (key, prop)
        if key in self.properties and not allow_overwrites:
            if not prop._is_same_card(self.properties[key]):
                self._duplicate_properties.append(prop)
                if self._stop_on_duplicate_error:
                    self.pop_parse_errors()
        else:
            self.properties[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def add_property_mass(self, prop, allow_overwrites=False):
        key = prop.pid
        if key in self.properties_mass and not allow_overwrites:
            if not prop._is_same_card(self.properties_mass[key]):
                #print('pid=%s\noldProperty=\n%snewProperty=\n%s' %(key,self.properties_mass[key],prop))
                assert key not in self.properties_mass, 'pid=%s oldProperty=\n%snewProperty=\n%s' % (key, self.properties_mass[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            self.properties_mass[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def add_DTABLE(self, dtable, allow_overwrites=False):
        if self.dtable is not None:
            if not dtable._is_same_card(self.dtable):
                raise RuntimeError('DTABLE cannot be overwritten\nold:\n%s\nnew:\n%s',
                                   self.dtable, dtable)
        else:
            self.dtable = dtable
            #self._type_to_id_map[dtable.type].append(1)

    def add_BCRPARA(self, card, allow_overwrites=False):
        key = card.crid
        self.bcrparas[key] = card
        self._type_to_id_map[card.type].append(key)

    def add_BCTADD(self, card, allow_overwrites=False):
        key = card.csid
        self.bctadds[key] = card
        self._type_to_id_map[card.type].append(key)

    def add_BCTPARA(self, card, allow_overwrites=False):
        key = card.csid
        self.bctparas[key] = card
        self._type_to_id_map[card.type].append(key)

    def add_BCTSET(self, card, allow_overwrites=False):
        key = card.csid
        self.bctsets[key] = card
        self._type_to_id_map[card.type].append(key)

    def add_BSURF(self, card, allow_overwrites=False):
        key = card.sid
        self.bsurf[key] = card
        self._type_to_id_map[card.type].append(key)

    def add_BSURFS(self, card, allow_overwrites=False):
        key = card.id
        self.bsurfs[key] = card
        self._type_to_id_map[card.type].append(key)

    def add_DELAY(self, delay, allow_overwrites=False):
        key = delay.sid
        assert key > 0, 'sid=%s delay=%s' % (key, delay)
        if key in self.delays:
            self.delays[key].add(delay)
        else:
            self.delays[key] = delay
            self._type_to_id_map[delay.type].append(key)

    def add_TEMPD(self, tempd, allow_overwrites=False):
        key = tempd.sid
        if key in self.tempds and not allow_overwrites:
            if not tempd._is_same_card(self.tempds[key]):
                assert key not in self.tempds, 'TEMPD.sid=%s old=\n%snew=\n%s' % (
                    key, self.tempds[key], tempd)
        else:
            assert key > 0, 'sid=%s tempd=%s' % (key, tempd)
            self.tempds[key] = tempd
            self._type_to_id_map[tempd.type].append(key)

    def add_PBUSHT(self, prop, allow_overwrites=False):
        key = prop.pid
        if key in self.pbusht and not allow_overwrites:
            if not prop._is_same_card(self.pbusht[key]):
                assert key not in self.pbusht, 'PBUSHT.pid=%s old=\n%snew=\n%s' % (
                    key, self.pbusht[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            self.pbusht[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def add_PDAMPT(self, prop, allow_overwrites=False):
        key = prop.pid
        if key in self.pdampt and not allow_overwrites:
            if not prop._is_same_card(self.pdampt[key]):
                assert key not in self.pdampt, 'PDAMPT.pid=%s old=\n%snew=\n%s' % (
                    key, self.pdampt[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            self.pdampt[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def add_PELAST(self, prop, allow_overwrites=False):
        key = prop.pid
        assert key > 0, 'pid=%s prop=%s' % (key, prop)
        if key in self.pelast and not allow_overwrites:
            if not prop._is_same_card(self.pelast[key]):
                #print('pid=%s\noldProperty=\n%snewProperty=\n%s' % (key, self.pelast[key],prop))
                assert key not in self.pelast, 'PELAST.pid=%s old=\n%snew=\n%s' % (
                    key, self.pelast[key], prop)
        else:
            self.pelast[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def add_TF(self, tf, allow_overwrites=False):
        key = tf.sid
        assert key > 0, 'sid=%s tf=%s' % (key, tf)
        if key in self.transfer_functions:
            self.transfer_functions[key].append(tf)
        else:
            self.transfer_functions[key] = [tf]
            self._type_to_id_map[tf.type].append(key)

    def add_structural_material(self, material, allow_overwrites=False):
        key = material.mid
        assert key > 0, 'mid=%s material=\n%s' % (key, material)
        if key in self.materials and not allow_overwrites:
            if not material._is_same_card(self.materials[key]):
                self._duplicate_materials.append(material)
        else:
            self.materials[key] = material
            self._type_to_id_map[material.type].append(key)

    def add_thermal_material(self, material, allow_overwrites=False):
        key = material.mid
        assert key > 0, 'mid=%s material=\n%s' % (key, material)
        if key in self.thermal_materials and not allow_overwrites:
            if not material._is_same_card(self.thermal_materials[key]):
                self._duplicate_thermal_materials.append(material)
        else:
            self.thermal_materials[key] = material
            self._type_to_id_map[material.type].append(key)

    def add_hyperelastic_material(self, material, allow_overwrites=False):
        key = material.mid
        assert key > 0, 'mid=%s material=\n%s' % (key, material)
        if key in self.hyperelastic_materials and not allow_overwrites:
            if not material._is_same_card(self.hyperelastic_materials[key]):
                assert key not in self.hyperelastic_materials, 'mid=%s\nold=\n%snew=\n%s' % (key, self.hyperelastic_materials[key], material)
        else:
            self.hyperelastic_materials[key] = material
            self._type_to_id_map[material.type].append(key)

    def add_material_dependence(self, material, allow_overwrites=False):
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
        if key in slot and not allow_overwrites:
            if not material._is_same_card(slot[key]):
                assert key not in slot, 'dMATx.mid=%s Type=%r\nold=\n%snew=\n%s' % (key, Type, slot[key], material)
        else:
            assert key > 0, 'mid=%s material=\n%s' % (key, material)
            slot[key] = material
            self._type_to_id_map[material.type].append(key)

    def add_creep_material(self, material, allow_overwrites=False):
        """
        .. note:: May be removed in the future.  Are CREEP cards materials?
                  They have an MID, but reference structural materials.
        """
        key = material.mid
        if key in self.thermal_materials and not allow_overwrites:
            if not material._is_same_card(self.creep_materials[key]):
                assert key not in self.creep_materials, 'Material.mid=%s\nold=\n%snew=\n%s' % (key, self.creep_materials[key], material)
        else:
            assert key > 0, 'mid=%s material=\n%s' % (key, material)
            self.creep_materials[key] = material
            self._type_to_id_map[material.type].append(key)

    def add_coord(self, coord, allow_overwrites=False):
        key = coord.cid
        assert coord.cid > -1, 'cid=%s coord=\n%s' % (key, coord)
        if key in self.coords:
            #if not allow_overwrites:
            if not coord._is_same_card(self.coords[key]):
                self._duplicate_coords.append(coord)
        else:
            self.coords[key] = coord
            self._type_to_id_map[coord.type].append(key)

    def add_load(self, load):
        key = load.sid
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]
            self._type_to_id_map[load.type].append(key)

    def add_dload(self, load):
        key = load.sid
        if key in self.dloads:
            self.dloads[key].append(load)
        else:
            self.dloads[key] = [load]
            self._type_to_id_map[load.type].append(key)

    def add_dload_entry(self, load):
        key = load.sid
        if key in self.dload_entries:
            self.dload_entries[key].append(load)
        else:
            self.dload_entries[key] = [load]
            self._type_to_id_map[load.type].append(key)

    def add_LSEQ(self, load):
        key = load.sid
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]
            self._type_to_id_map[load.type].append(key)

    def add_thermal_load(self, load):  # same function at the moment...
        key = load.sid
        assert key > 0, 'key=%s; load=%s\n' % (key, load)
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]
            self._type_to_id_map[load.type].append(key)

    def add_PHBDY(self, prop):
        key = prop.pid
        if key in self.phbdys:
            if not prop._is_same_card(self.phbdys[key]):
                assert key not in self.phbdys, 'PHBDY.pid=%s\nold=\n%snew=\n%s' % (
                    key, self.phbdys[key], prop)
        else:
            assert key > 0, 'pid=%s prop=\n%s' % (key, prop)
            self.phbdys[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def add_convection_property(self, prop):
        key = prop.pconid
        assert key > 0, key
        assert key not in self.convection_properties, key
        self.convection_properties[key] = prop
        self._type_to_id_map[prop.type].append(key)

    def add_thermal_BC(self, bc, key):
        assert key > 0
        if key in self.bcs:
            self.bcs[key].append(bc)
        else:
            self.bcs[key] = [bc]
            self._type_to_id_map[bc.type].append(key)

    #def add_constraint_MPCADD(self, constraint):
        #raise RuntimeError('is this used?')
        #key = constraint.conid
        #if key in self.mpcadds:
            #raise RuntimeError('must have unique MPCADD IDs')
        #self.mpcadds[key] = constraint
        #self._type_to_id_map[constraint.type].append(key)

    def add_constraint_MPC(self, constraint):
        key = constraint.conid
        if key in self.mpcs:
            self.mpcs[key].append(constraint)
        else:
            self.mpcs[key] = [constraint]
            self._type_to_id_map[constraint.type].append(key)

    #def add_constraint_SPCADD(self, constraint):
        #raise RuntimeError('is this used?')
        #key = constraint.conid
        #if key in self.spcadds:
            #raise RuntimeError('must have unique SPCADD IDs')
        #self.spcadds[key] = constraint
        #self._type_to_id_map[constraint.type].append(key)

    def add_constraint_SPC(self, constraint):
        key = constraint.conid
        if key in self.spcs:
            self.spcs[key].append(constraint)
        else:
            self.spcs[key] = [constraint]
            self._type_to_id_map[constraint.type].append(key)

    def add_constraint(self, constraint):
        key = constraint.conid
        if key in self.spcs:
            self.spcs[key].append(constraint)
        else:
            self.spcs[key] = [constraint]
            self._type_to_id_map[constraint.type].append(key)

    def add_sesuport(self, se_suport):
        self._type_to_id_map[se_suport.type].append(len(self.se_suport))
        self.se_suport.append(se_suport)

    def add_suport(self, suport):
        self._type_to_id_map[suport.type].append(len(self.suport))
        self.suport.append(suport)

    def add_suport1(self, suport1):
        key = suport1.conid
        if key in self.suport1:
            self.suport1[key].add_suport1_to_set(suport1)
        else:
            assert suport1.conid > 0
            self.suport1[key] = suport1
            self._type_to_id_map[suport1.type].append(key)

    def add_DAREA(self, darea, allow_overwrites=False):
        key = (darea.sid, darea.p, darea.c)
        if key in self.dareas and not allow_overwrites:
            if not darea._is_same_card(self.dareas[key]):
                assert key not in self.dareas, '\ndarea=\n%s oldDArea=\n%s' % (darea, self.dareas[key])
        else:
            assert darea.sid > 0
            self.dareas[key] = darea
            self._type_to_id_map[darea.type].append(key)

    def add_DPHASE(self, dphase, allow_overwrites=False):
        #key = (dphase.sid, dphase.phase_leads, dphase.components)
        key = dphase.sid
        if key in self.dphases and not allow_overwrites:
            if not dphase._is_same_card(self.dphases[key]):
                assert key not in self.dphases, '\ndphase=\n%s old_DPHASE=\n%s' % (
                    dphase, self.dphases[key])
        else:
            assert dphase.sid > 0
            self.dphases[key] = dphase
            self._type_to_id_map[dphase.type].append(key)

    def add_AERO(self, aero):
        # only one AERRO card allowed
        assert self.aero is None, '\naero=\n%s old=\n%s' % (aero, self.aero)
        self.aero = aero
        #self._type_to_id_map[aero.type].append(key)

    def add_AEROS(self, aeros):
        # only one AERROS card allowed
        assert self.aeros is None, '\naeros=\n%s old=\n%s' % (aeros, self.aeros)
        self.aeros = aeros
        #self._type_to_id_map[aeros.type].append(key)

    def add_AEFACT(self, aefact, allow_overwrites=False):
        key = aefact.sid
        if key in self.aefacts and not allow_overwrites:
            if not aefact._is_same_card(self.aefacts[key]):
                assert key not in self.aefacts, 'AEFACT.sid=%s\nold=\n%snew=\n%s' % (key, self.aefacts[key], aefact)
        else:
            assert key > 0, 'sid=%s method=\n%s' % (key, aefact)
            self.aefacts[key] = aefact
            self._type_to_id_map[aefact.type].append(key)

    def add_AELIST(self, aelist):
        key = aelist.sid
        assert key not in self.aelists, 'AELIST.sid=%s\nold=\n%snew=\n%s' % (key, self.aelists[key], aelist)
        assert key >= 0
        self.aelists[key] = aelist
        self._type_to_id_map[aelist.type].append(key)

    def add_AELINK(self, aelink):
        key = aelink.id
        assert key >= 0
        if key not in self.aelinks:
            self.aelinks[key] = []
        self.aelinks[key].append(aelink)
        self._type_to_id_map[aelink.type].append(key)
        #assert key not in self.aestats,'\naestat=%s oldAESTAT=\n%s' %(aestat,self.aestats[key])

    def add_AECOMP(self, aecomp):
        key = aecomp.name
        assert key not in self.aecomps, '\naecomp=\n%s oldAECOMP=\n%s' % (aecomp, self.aecomps[key])
        self.aecomps[key] = aecomp
        self._type_to_id_map[aecomp.type].append(key)

    def add_AEPARM(self, aeparam):
        key = aeparam.id
        assert key not in self.aeparams, '\naeparam=\n%s oldAEPARM=\n%s' % (aeparam, self.aeparams[key])
        assert key >= 0
        self.aeparams[key] = aeparam
        self._type_to_id_map[aeparam.type].append(key)

    def add_AESTAT(self, aestat):
        key = aestat.id
        assert key not in self.aestats, '\naestat=\n%s old=\n%s' % (
            aestat, self.aestats[key])
        assert key >= 0
        self.aestats[key] = aestat
        self._type_to_id_map[aestat.type].append(key)

    def add_AESURF(self, aesurf):
        key = aesurf.aesid
        assert key not in self.aesurf, '\naesurf=\n%s old=\n%s' % (
            aesurf, self.aesurf[key])
        assert key >= 0
        self.aesurf[key] = aesurf
        self._type_to_id_map[aesurf.type].append(key)

    def add_AESURFS(self, aesurfs):
        key = aesurfs.aesid
        assert key not in self.aesurf, '\naesurfs=\n%s old=\n%s' % (
            aesurfs, self.aesurfs[key])
        assert key >= 0
        self.aesurfs[key] = aesurfs
        self._type_to_id_map[aesurfs.type].append(key)

    def add_CSSCHD(self, csschd):
        key = csschd.sid
        assert key not in self.csschds, '\nCSSCHD=\n%s old=\n%s' % (csschd, self.csschds[key])
        assert key >= 0
        self.csschds[key] = csschd
        self._type_to_id_map[csschd.type].append(key)

    def add_CAERO(self, caero):
        key = caero.eid
        assert key not in self.caeros, '\ncaero=\n|%s| oldCAERO=\n|%s|' % (
            caero, self.caeros[key])
        assert key > 0
        self.caeros[key] = caero
        self._type_to_id_map[caero.type].append(key)

    def add_PAERO(self, paero):
        key = paero.pid
        assert key not in self.paeros, '\npaero=\n|%s| oldPAERO=\n|%s|' % (
            paero, self.paeros[key])
        assert key > 0, 'paero.pid = %r' % (key)
        self.paeros[key] = paero
        self._type_to_id_map[paero.type].append(key)

    def add_MONPNT(self, monitor_point):
        key = monitor_point.name
        assert key not in self.monitor_points, '\nmonitor_point=\n%soldMNTPNT=\n%s' % (
            monitor_point, self.monitor_points[key])
        self.monitor_points.append(monitor_point)
        self._type_to_id_map[monitor_point.type].append(len(self.monitor_points) - 1)

    def add_SPLINE(self, spline):
        assert spline.eid not in self.splines
        assert spline.eid > 0
        key = spline.eid
        self.splines[key] = spline
        self._type_to_id_map[spline.type].append(key)

    def add_gust(self, gust):
        key = gust.sid
        assert key not in self.gusts
        assert key > 0
        self.gusts[key] = gust
        self._type_to_id_map[gust.type].append(key)

    def add_TRIM(self, trim, allow_overwrites=False):
        key = trim.sid
        if not allow_overwrites:
            assert key not in self.trims, 'TRIM=%s  old=\n%snew=\n%s' % (key, self.trims[key], trim)
        assert key > 0, 'key=%r trim=\n%s' % (key, trim)
        self.trims[key] = trim
        self._type_to_id_map[trim.type].append(key)

    def add_DIVERG(self, diverg, allow_overwrites=False):
        key = diverg.sid
        if not allow_overwrites:
            assert key not in self.divergs, 'DIVERG=%s  old=\n%snew=\n%s' % (key, self.divergs[key], diverg)
        assert key > 0, 'key=%r diverg=\n%s' % (key, diverg)
        self.divergs[diverg] = diverg
        self._type_to_id_map[diverg.type].append(key)

    def add_FLUTTER(self, flutter):
        key = flutter.sid
        assert key not in self.flutters, 'FLUTTER=%s old=\n%snew=\n%s' % (key, self.flutters[key], flutter)
        assert key > 0
        self.flutters[key] = flutter
        self._type_to_id_map[flutter.type].append(key)

    def add_FLFACT(self, flfact):
        key = flfact.sid
        #assert key not in self.flfacts
        assert key > 0
        self.flfacts[key] = flfact  # set id...
        self._type_to_id_map[flfact.type].append(key)

    def add_DCONSTR(self, dconstr):
        #key = (dconstr.oid, dconstr.rid)
        key = dconstr.oid
        #assert key not in self.dconstrs, 'key=%r DCONSTR/DCONADD=\n%s' % (key, dconstr)
        assert dconstr.oid > 0
        #assert dconstr.dresp_id > 0
        if key in self.dconstrs:
            self.dconstrs[key].append(dconstr)
        else:
            self.dconstrs[key] = [dconstr]
        self._type_to_id_map[dconstr.type].append(key)

    #def add_DCONADD(self, dconadd, allow_overwrites=False):
        #key = dconadd.oid
        #if key in self.dconstrs and not allow_overwrites:
            #if not dconadd._is_same_card(self.dconstrs[key]):
                #assert key not in self.dconstrs, 'DCONADD=%s old=\n%snew=\n%s' % (
                    #key, self.dconstrs[key], dconadd)
        #else:
            #assert key > 0, 'dcid=%s dconadd=%s' % (key, dconadd)
            #self.dconstrs[key] = dconadd
            #self._type_to_id_map[dconadd.type].append(key)

    def add_DESVAR(self, desvar):
        key = desvar.desvar_id
        assert key not in self.desvars, 'DESVAR=%s old=\n%snew=\n%s' % (
            key, self.desvars[key], desvar)
        assert key > 0
        self.desvars[key] = desvar
        self._type_to_id_map[desvar.type].append(key)

    def add_DDVAL(self, ddval):
        key = ddval.oid
        assert key not in self.ddvals, 'DDVAL=%s old=\n%snew=\n%s' % (
            key, self.ddvals[key], ddval)
        assert key > 0
        self.ddvals[key] = ddval
        self._type_to_id_map[ddval.type].append(key)

    def add_DLINK(self, dlink):
        key = dlink.oid
        assert key not in self.dlinks, 'DLINK=%s old=\n%snew=\n%s' % (
            key, self.dlinks[key], dlink)
        assert key > 0
        self.dlinks[key] = dlink
        self._type_to_id_map[dlink.type].append(key)

    def add_DRESP(self, dresp):
        key = dresp.dresp_id
        assert key not in self.dresps, 'DRESPx=%s old=\n%snew=\n%s' % (
                    key, self.dresps[key], dresp)
        assert key > 0
        self.dresps[key] = dresp
        self._type_to_id_map[dresp.type].append(key)

    def add_DVCREL(self, dvcrel):
        key = dvcrel.oid
        assert key not in self.dvprels, 'DVCRELx=%s old\n%snew=\n%s' % (
                    key, self.dvcrels[key], dvcrel)
        assert key > 0
        self.dvcrels[key] = dvcrel
        self._type_to_id_map[dvcrel.type].append(key)

    def add_DVMREL(self, dvmrel):
        key = dvmrel.oid
        assert key not in self.dvmrels, 'DVMRELx=%s old=\n%snew=\n%s' % (
                    key, self.dvmrels[key], dvmrel)
        assert key not in self.dvmrels
        assert key > 0
        self.dvmrels[key] = dvmrel
        self._type_to_id_map[dvmrel.type].append(key)

    def add_DVPREL(self, dvprel):
        key = dvprel.oid
        assert key not in self.dvprels, 'DVPRELx=%s old\n%snew=\n%s' % (
                    key, self.dvprels[key], dvprel)
        assert key > 0
        self.dvprels[key] = dvprel
        self._type_to_id_map[dvprel.type].append(key)

    def add_NLPARM(self, nlparm):
        key = nlparm.nlparm_id
        assert key not in self.nlparms
        assert key > 0, 'key=%s; nlparm=%s\n' % (key, nlparm)
        self.nlparms[key] = nlparm
        self._type_to_id_map[nlparm.type].append(key)

    def add_NLPCI(self, nlpci):
        key = nlpci.nlpci_id
        assert key not in self.nlpcis
        assert key > 0
        self.nlpcis[key] = nlpci
        self._type_to_id_map[nlpci.type].append(key)

    def add_TSTEP(self, tstep, allow_overwrites=False):
        key = tstep.sid
        if key in self.tsteps and not allow_overwrites:
            if not tstep._is_same_card(self.tsteps[key]):
                assert key not in self.tsteps, 'TSTEP=%s\nold=\n%snew=\n%s' % (key, self.tsteps[key], tstep)
        else:
            assert key > 0, 'sid=%s tstep=\n%s' % (key, tstep)
            self.tsteps[key] = tstep
            self._type_to_id_map[tstep.type].append(key)

    def add_TSTEPNL(self, tstepnl, allow_overwrites=False):
        key = tstepnl.sid
        if key in self.tstepnls and not allow_overwrites:
            if not tstepnl._is_same_card(self.tstepnls[key]):
                assert key not in self.tstepnls, 'TSTEPNL=%s\nold=\n%snew=\n%s' % (key, self.tstepnls[key], tstepnl)
        else:
            assert key > 0, 'sid=%s tstepnl=\n%s' % (key, tstepnl)
            self.tstepnls[key] = tstepnl
            self._type_to_id_map[tstepnl.type].append(key)

    def add_FREQ(self, freq):
        key = freq.sid
        assert key > 0
        if key in self.frequencies:
            self.frequencies[key].add_frequency_object(freq)
        else:
            self.frequencies[key] = freq
            self._type_to_id_map[freq.type].append(key)

    def add_SET(self, set_obj):
        key = set_obj.sid
        assert key >= 0
        if key in self.sets:
            self.sets[key].add_set(set_obj)
        else:
            self.sets[key] = set_obj
            self._type_to_id_map[set_obj.type].append(key)

    def add_ASET(self, set_obj):
        self.asets.append(set_obj)

    def add_BSET(self, set_obj):
        self.bsets.append(set_obj)

    def add_CSET(self, set_obj):
        self.csets.append(set_obj)

    def add_QSET(self, set_obj):
        self.qsets.append(set_obj)

    def add_USET(self, set_obj):
        key = set_obj.name
        if key in self.usets:
            self.usets[key].append(set_obj)
        else:
            self.usets[key] = [set_obj]

    def add_SEBSET(self, set_obj):
        self.se_bsets.append(set_obj)

    def add_SECSET(self, set_obj):
        self.se_csets.append(set_obj)

    def add_SEQSET(self, set_obj):
        self.se_qsets.append(set_obj)

    def add_SEUSET(self, set_obj):
        key = set_obj.name
        if key in self.se_usets:
            self.se_usets[key].append(set_obj)
        else:
            self.se_usets[key] = [set_obj]

    def add_SESET(self, set_obj):
        key = set_obj.seid
        assert key >= 0
        if key in self.se_sets:
            old_set = self.se_sets[key]
            set_obj.add_seset(old_set)
        self.se_sets[key] = set_obj
        self._type_to_id_map[set_obj.type].append(key)


    def add_table(self, table):
        key = table.tid
        assert key not in self.tables, '\nTable=\n%s oldTable=\n%s' % (
            table, self.tables[key])
        assert key > 0
        self.tables[key] = table
        self._type_to_id_map[table.type].append(key)

    def add_table_sdamping(self, table):
        key = table.tid
        assert key not in self.tables_sdamping, '\nTable=\n%s oldTable=\n%s' % (
            table, self.tables_sdamping[key])
        assert key > 0
        self.tables_sdamping[key] = table
        self._type_to_id_map[table.type].append(key)

    def add_random_table(self, table):
        key = table.tid
        assert key not in self.random_tables, '\nTable=\n%s old=\n%s' % (
            table, self.random_tables[key])
        assert key > 0
        self.random_tables[key] = table
        self._type_to_id_map[table.type].append(key)

    def add_method(self, method, allow_overwrites=False):
        key = method.sid
        if key in self.methods and not allow_overwrites:
            if not method._is_same_card(self.methods[key]):
                assert key not in self.methods, 'sid=%s\noldMethod=\n%snewMethod=\n%s' % (key, self.methods[key], method)
        else:
            assert key > 0, 'sid=%s method=\n%s' % (key, method)
            self.methods[key] = method
            self._type_to_id_map[method.type].append(key)

    def add_cmethod(self, method, allow_overwrites=False):
        key = method.sid
        if key in self.cMethods and not allow_overwrites:
            if not method._is_same_card(self.cMethods[key]):
                assert key not in self.cMethods, 'sid=%s\noldCMethod=\n%snewCMethod=\n%s' % (key, self.cMethods[key], method)
        else:
            assert key > 0, 'sid=%s cMethod=\n%s' % (key, method)
            self.cMethods[key] = method
            self._type_to_id_map[method.type].append(key)

    def add_MKAERO(self, mkaero):
        self.mkaeros.append(mkaero)
