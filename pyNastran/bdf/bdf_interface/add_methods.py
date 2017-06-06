# pylint: disable=R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from pyNastran.bdf.bdf_interface.attributes import BDFAttributes


class AddMethods(BDFAttributes):
    """defines methods to add card objects to the BDF"""
    def __init__(self):
        BDFAttributes.__init__(self)

    def _add_dmi_object(self, dmi, allow_overwrites=False):
        """adds a DMI object"""
        name = dmi.name
        self.dmis[name] = dmi
        self._type_to_id_map[dmi.type].append(name)

    def _add_dmig_object(self, dmig, allow_overwrites=False):
        """adds a DMIG object"""
        name = dmig.name
        self.dmigs[name] = dmig
        self._type_to_id_map[dmig.type].append(name)

    def _add_dmij_object(self, dmij, allow_overwrites=False):
        """adds a DMIJ object"""
        name = dmij.name
        self.dmijs[name] = dmij
        self._type_to_id_map[dmij.type].append(name)

    def _add_dmiji_object(self, dmiji, allow_overwrites=False):
        """adds a DMIJI object"""
        name = dmiji.name
        self.dmijis[name] = dmiji
        self._type_to_id_map[dmiji.type].append(name)

    def _add_dmik_object(self, dmik, allow_overwrites=False):
        """adds a DMIK object"""
        name = dmik.name
        self.dmiks[name] = dmik
        self._type_to_id_map[dmik.type].append(name)

    def _add_dti_object(self, dti, allow_overwrites=False):
        """adds an DTI object"""
        name = dti.name
        self.dti[name] = dti
        self._type_to_id_map[dti.type].append(name)

    def _add_param_object(self, param, allow_overwrites=False):
        """adds a PARAM object"""
        key = param.key
        if key in self.params and not allow_overwrites:
            if not param._is_same_card(self.params[key]):
                #if param.key in self.params:
                    #msg = 'key=%s param=%s old_param=%s' % (key, param, self.params[key])
                    #raise KeyError(msg)
                self.log.warning('key=%s param=%s old_param=%s' %
                                 (key, param, self.params[key]))
                self.params[key] = param
        else:
            self.params[key] = param
            self._type_to_id_map[param.type].append(key)

    def _add_node_object(self, node, allow_overwrites=False):
        """adds a GRID card"""
        key = node.nid
        if key in self.nodes and not allow_overwrites:
            if not node._is_same_card(self.nodes[key]):
                assert node.nid not in self.nodes, 'nid=%s\nold_node=\n%snew_node=\n%s' % (node.nid, self.nodes[key], node)
            else:
                #print('GRID was duplicated...nid=%s; node=\n%s' % (key, node))
                pass
        else:
            assert key > 0, 'nid=%s node=%s' % (key, node)
            self.nodes[key] = node
            self._type_to_id_map[node.type].append(key)

    def _add_seqgp_object(self, seqgp):
        """adds an SPOINT card"""
        if self.seqgp is None:
            self.seqgp = seqgp
        else:
            self.seqgp.append(seqgp)

    def _add_point_object(self, point, allow_overwrites=False):
        """adds a POINT card"""
        key = point.nid
        if key in self.points and not allow_overwrites:
            if not point._is_same_card(self.points[key]):
                assert point.nid not in self.points, 'nid=%s\nold_point=\n%snew_point=\n%s' % (point.nid, self.points[key], point)
            else:
                #print('POINT was duplicated...nid=%s; point=\n%s' % (key, point))
                pass
        else:
            assert key > 0, 'nid=%s point=%s' % (key, point)
            self.points[key] = point
            self._type_to_id_map[point.type].append(key)

    def _add_spoint_object(self, spoints):
        """adds an SPOINT card"""
        if self.spoints is None:
            self.spoints = spoints
        else:
            self.spoints.add_points(spoints.points)

    def _add_epoint_object(self, epoints):
        """adds an EPOINT card"""
        if self.epoints is None:
            self.epoints = epoints
        else:
            self.epoints.add_points(epoints.points)

    def _add_plotel_object(self, elem, allow_overwrites=False):
        """adds an PLOTEL object"""
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

    #def add_element(self, elem, allow_overwrites=False):
        #"""deprecated"""
        #return self._add_element_object(elem, allow_overwrites)

    def _add_element_object(self, elem, allow_overwrites=False):
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

    def _add_ao_object(self, elem_flag, allow_overwrites=False):
        """adds a CBARAO"""
        key = elem_flag.eid
        assert key > 0, 'eid=%s must be positive; elem_flag=\n%s' % (key, elem_flag)
        if key in self.ao_element_flags and not allow_overwrites:
            if not elem_flag._is_same_card(self.ao_element_flags[key]):
                #self._duplicate_elements.append(elem_flag)
                #if self._stop_on_duplicate_error:
                    #self.pop_parse_errors()
                assert elem_flag.eid not in self.ao_element_flags, 'eid=%s\nold_ao_element=\n%snew_ao_element=\n%s' % (
                    elem_flag.eid, self.ao_element_flags[elem_flag.eid], elem_flag)
        else:
            self.ao_element_flags[key] = elem_flag
            self._type_to_id_map[elem_flag.type].append(key)

    def _add_doptprm_object(self, doptprm, comment=''):
        """adds a DOPTPRM"""
        self.doptprm = doptprm

    def _add_nsm_object(self, nsm, allow_overwrites=False):
        """adds a nsm object to a nsm set"""
        key = nsm.sid
        assert key > 0, 'sid=%s must be positive; nsm=\n%s' % (key, nsm)
        if key in self.nsms:
            self.nsms[key].append(nsm)
        else:
            self.nsms[key] = [nsm]
            self._type_to_id_map[nsm.type].append(key)

    def _add_mass_object(self, mass, allow_overwrites=False):
        key = mass.eid
        assert key > 0, 'eid=%s must be positive; mass=\n%s' % (key, mass)
        if key in self.masses and not allow_overwrites:
            if not mass._is_same_card(self.masses[key]):
                self._duplicate_masses.append(mass)
        else:
            self.masses[key] = mass
            self._type_to_id_map[mass.type].append(key)

    def _add_damper_object(self, elem, allow_overwrites=False):
        """.. warning:: can dampers have the same ID as a standard element?"""
        return self._add_element_object(elem, allow_overwrites)

    def _add_rigid_element_object(self, elem, allow_overwrites=False):
        key = elem.eid
        assert key > 0, 'eid=%s elem=%s' % (key, elem)
        if key in self.rigid_elements and not allow_overwrites:
            assert elem.eid not in self.rigid_elements, 'eid=%s\noldElement=\n%snewElement=\n%s' % (elem.eid, self.rigid_elements[elem.eid], elem)
        self.rigid_elements[key] = elem
        self._type_to_id_map[elem.type].append(key)

    def _add_thermal_element_object(self, elem):
        """same as add_element at the moment..."""
        self._add_element_object(elem)

    def _add_deqatn_object(self, deqatn, allow_overwrites=False):
        """adds an DEQATN object"""
        key = deqatn.equation_id
        assert key > 0, 'ID=%s deqatn\n%s' % (key, deqatn)
        if key in self.dequations and not allow_overwrites:
            if not deqatn.write_card() == self.dequations[key].write_card():
                assert key not in self.dequations, 'id=%s old_eq=\n%snew_eq=\n%s' % (
                    key, self.dequations[key], deqatn)
        self.dequations[key] = deqatn
        self._type_to_id_map[deqatn.type].append(key)

    #def add_property(self, prop, allow_overwrites=False):
        #"""deprecated"""
        #return self._add_property_object(prop, allow_overwrites)

    def _add_property_object(self, prop, allow_overwrites=False):
        """
        adds one of the following objects:
          PELAS, PBUSH, PBUSH1D, PBUSH2D, PDAMP,
          PROD, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
          PSHELL, PCOMP, PCOMPG,
          PSOLID, PLSOLID
        """
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

    def _add_property_mass_object(self, prop, allow_overwrites=False):
        """adds an PMASS object"""
        key = prop.pid
        if key in self.properties_mass and not allow_overwrites:
            if not prop._is_same_card(self.properties_mass[key]):
                #print('pid=%s\noldProperty=\n%snewProperty=\n%s' %(key,self.properties_mass[key],prop))
                assert key not in self.properties_mass, 'pid=%s oldProperty=\n%snewProperty=\n%s' % (key, self.properties_mass[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            self.properties_mass[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def _add_dtable_object(self, dtable, allow_overwrites=False):
        """adds an DTABLE object"""
        if self.dtable is not None:
            if not dtable._is_same_card(self.dtable):
                raise RuntimeError('DTABLE cannot be overwritten\nold:\n%s\nnew:\n%s',
                                   self.dtable, dtable)
        else:
            self.dtable = dtable
            #self._type_to_id_map[dtable.type].append(1)

    def _add_bcrpara_object(self, card, allow_overwrites=False):
        """adds an BCRPARA object"""
        key = card.crid
        self.bcrparas[key] = card
        self._type_to_id_map[card.type].append(key)

    def _add_bctadd_object(self, card, allow_overwrites=False):
        """adds an BCTADD object"""
        key = card.csid
        self.bctadds[key] = card
        self._type_to_id_map[card.type].append(key)

    def _add_bctpara_object(self, card, allow_overwrites=False):
        """adds an BCTPARA object"""
        key = card.csid
        self.bctparas[key] = card
        self._type_to_id_map[card.type].append(key)

    def _add_bctset_object(self, card, allow_overwrites=False):
        """adds an BCTSET object"""
        key = card.csid
        self.bctsets[key] = card
        self._type_to_id_map[card.type].append(key)

    def _add_bsurf_object(self, card, allow_overwrites=False):
        """adds an BSURF object"""
        key = card.sid
        self.bsurf[key] = card
        self._type_to_id_map[card.type].append(key)

    def _add_bsurfs_object(self, card, allow_overwrites=False):
        """adds an BSURFS object"""
        key = card.id
        self.bsurfs[key] = card
        self._type_to_id_map[card.type].append(key)

    def _add_tempd_object(self, tempd, allow_overwrites=False):
        """adds an TEMPD object"""
        key = tempd.sid
        if key in self.tempds and not allow_overwrites:
            if not tempd._is_same_card(self.tempds[key]):
                assert key not in self.tempds, 'TEMPD.sid=%s old=\n%snew=\n%s' % (
                    key, self.tempds[key], tempd)
        else:
            assert key > 0, 'sid=%s tempd=%s' % (key, tempd)
            self.tempds[key] = tempd
            self._type_to_id_map[tempd.type].append(key)

    def _add_pbusht_object(self, prop, allow_overwrites=False):
        """adds an PBUSHT object"""
        key = prop.pid
        if key in self.pbusht and not allow_overwrites:
            if not prop._is_same_card(self.pbusht[key]):
                assert key not in self.pbusht, 'PBUSHT.pid=%s old=\n%snew=\n%s' % (
                    key, self.pbusht[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            self.pbusht[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def _add_pdampt_object(self, prop, allow_overwrites=False):
        """adds an PDAMPT object"""
        key = prop.pid
        if key in self.pdampt and not allow_overwrites:
            if not prop._is_same_card(self.pdampt[key]):
                assert key not in self.pdampt, 'PDAMPT.pid=%s old=\n%snew=\n%s' % (
                    key, self.pdampt[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            self.pdampt[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def _add_pelast_object(self, prop, allow_overwrites=False):
        """adds an PELAST object"""
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

    def _add_tf_object(self, tf, allow_overwrites=False):
        """adds an TF (transfer function) object"""
        key = tf.sid
        assert key > 0, 'sid=%s tf=%s' % (key, tf)
        if key in self.transfer_functions:
            self.transfer_functions[key].append(tf)
        else:
            self.transfer_functions[key] = [tf]
            self._type_to_id_map[tf.type].append(key)

    def _add_structural_material_object(self, material, allow_overwrites=False):
        """adds an MAT1, MAT2, MAT8 object"""
        key = material.mid
        assert key > 0, 'mid=%s material=\n%s' % (key, material)
        if key in self.materials and not allow_overwrites:
            if not material._is_same_card(self.materials[key]):
                self._duplicate_materials.append(material)
        else:
            self.materials[key] = material
            self._type_to_id_map[material.type].append(key)

    def _add_thermal_material_object(self, material, allow_overwrites=False):
        """adds an MAT4, MAT5 object"""
        key = material.mid
        assert key > 0, 'mid=%s material=\n%s' % (key, material)
        if key in self.thermal_materials and not allow_overwrites:
            if not material._is_same_card(self.thermal_materials[key]):
                self._duplicate_thermal_materials.append(material)
        else:
            self.thermal_materials[key] = material
            self._type_to_id_map[material.type].append(key)

    def _add_hyperelastic_material_object(self, material, allow_overwrites=False):
        """adds an MATHP, MATHE object"""
        key = material.mid
        assert key > 0, 'mid=%s material=\n%s' % (key, material)
        if key in self.hyperelastic_materials and not allow_overwrites:
            if not material._is_same_card(self.hyperelastic_materials[key]):
                assert key not in self.hyperelastic_materials, 'mid=%s\nold=\n%snew=\n%s' % (key, self.hyperelastic_materials[key], material)
        else:
            self.hyperelastic_materials[key] = material
            self._type_to_id_map[material.type].append(key)

    def _add_material_dependence_object(self, material, allow_overwrites=False):
        """
        adds the following objects:
            MATS1, MATS3, MATS8,
            MATT1, MATT2, MATT3,
            MATT4, MATT5, MATT8, MATT9
        """
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

    def _add_creep_material_object(self, material, allow_overwrites=False):
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

    #def add_coord(self, coord, allow_overwrites=False):
        #"""deprecated"""
        #self._add_coord_object(coord, allow_overwrites)

    def _add_coord_object(self, coord, allow_overwrites=False):
        key = coord.cid
        assert coord.cid > -1, 'cid=%s coord=\n%s' % (key, coord)
        if key in self.coords:
            #if not allow_overwrites:
            if not coord._is_same_card(self.coords[key]):
                self._duplicate_coords.append(coord)
        else:
            self.coords[key] = coord
            self._type_to_id_map[coord.type].append(key)

    def _add_load_object(self, load):
        """adds a load object to a load case"""
        key = load.sid
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]
            self._type_to_id_map[load.type].append(key)

    def _add_dload_object(self, load):
        """adds a dload object to a load case"""
        key = load.sid
        if key in self.dloads:
            self.dloads[key].append(load)
        else:
            self.dloads[key] = [load]
            self._type_to_id_map[load.type].append(key)

    def _add_dload_entry(self, load):
        """adds a sub-dload object to a load case"""
        key = load.sid
        if key in self.dload_entries:
            self.dload_entries[key].append(load)
        else:
            self.dload_entries[key] = [load]
            self._type_to_id_map[load.type].append(key)

    def _add_lseq_object(self, load):
        """adds a LSEQ object to a load case"""
        key = load.sid
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]
            self._type_to_id_map[load.type].append(key)

    def _add_thermal_load_object(self, load):  # same function at the moment...
        key = load.sid
        assert key > 0, 'key=%s; load=%s\n' % (key, load)
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]
            self._type_to_id_map[load.type].append(key)

    def _add_phbdy_object(self, prop):
        key = prop.pid
        if key in self.phbdys:
            if not prop._is_same_card(self.phbdys[key]):
                assert key not in self.phbdys, 'PHBDY.pid=%s\nold=\n%snew=\n%s' % (
                    key, self.phbdys[key], prop)
        else:
            assert key > 0, 'pid=%s prop=\n%s' % (key, prop)
            self.phbdys[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def _add_convection_property_object(self, prop):
        key = prop.pconid
        assert key > 0, key
        assert key not in self.convection_properties, key
        self.convection_properties[key] = prop
        self._type_to_id_map[prop.type].append(key)

    def _add_thermal_bc_object(self, bc, key):
        assert key > 0
        if key in self.bcs:
            self.bcs[key].append(bc)
        else:
            self.bcs[key] = [bc]
            self._type_to_id_map[bc.type].append(key)

    def _add_constraint_mpc_object(self, constraint):
        key = constraint.conid
        if key in self.mpcs:
            self.mpcs[key].append(constraint)
        else:
            self.mpcs[key] = [constraint]
            self._type_to_id_map[constraint.type].append(key)

    def _add_constraint_spc_object(self, constraint):
        key = constraint.conid
        if key in self.spcs:
            self.spcs[key].append(constraint)
        else:
            self.spcs[key] = [constraint]
            self._type_to_id_map[constraint.type].append(key)

    def _add_constraint_spcoff_object(self, constraint):
        """dumb key, but good enough..."""
        key = constraint.type
        if key in self.spcoffs:
            self.spcoffs[key].append(constraint)
        else:
            self.spcoffs[key] = [constraint]
            self._type_to_id_map[constraint.type].append(key)

    def _add_sesuport_object(self, se_suport):
        self._type_to_id_map[se_suport.type].append(len(self.se_suport))
        self.se_suport.append(se_suport)

    def _add_suport_object(self, suport):
        self._type_to_id_map[suport.type].append(len(self.suport))
        self.suport.append(suport)

    def _add_suport1_object(self, suport1):
        key = suport1.conid
        if key in self.suport1:
            self.suport1[key].add_suport1_to_set(suport1)
        else:
            assert suport1.conid > 0
            self.suport1[key] = suport1
            self._type_to_id_map[suport1.type].append(key)

    def _add_tic_object(self, tic, allow_overwrites=False):
        key = tic.sid
        if key in self.tics:
            self.tics[key].add(tic)
        else:
            assert tic.sid > 0
            self.tics[key] = tic
            self._type_to_id_map[tic.type].append(key)

    def _add_darea_object(self, darea, allow_overwrites=False):
        #key = (darea.sid, darea.p, darea.c)
        key = darea.sid
        if key in self.dareas:
            self.dareas[key].add(darea)
        else:
            assert darea.sid > 0
            self.dareas[key] = darea
            self._type_to_id_map[darea.type].append(key)

    def _add_dphase_object(self, dphase, allow_overwrites=False):
        #key = (dphase.sid, dphase.nid, dphase.component) # dphase.phase_lead,
        key = dphase.sid
        if key in self.dphases:
            self.dphases[key].add(dphase)
        else:
            assert dphase.sid > 0, key
            self.dphases[key] = dphase
            self._type_to_id_map[dphase.type].append(key)

    def _add_delay_object(self, delay, allow_overwrites=False):
        """adds an DELAY object"""
        #key = (delay.sid, delay.nid, delay.component)
        key = delay.sid
        assert key > 0, 'sid=%s delay=%s' % (key, delay)
        if key in self.delays:
            self.delays[key].add(delay)
        else:
            self.delays[key] = delay
            self._type_to_id_map[delay.type].append(key)

    def _add_aero_object(self, aero):
        """adds an AERO object"""
        # only one AERO card allowed
        assert self.aero is None, '\naero=\n%s old=\n%s' % (aero, self.aero)
        self.aero = aero
        #self._type_to_id_map[aero.type].append(key)

    def _add_aeros_object(self, aeros):
        """adds an AEROS object"""
        # only one AEROS card allowed
        assert self.aeros is None, '\naeros=\n%s old=\n%s' % (aeros, self.aeros)
        self.aeros = aeros
        #self._type_to_id_map[aeros.type].append(key)

    def _add_aefact_object(self, aefact, allow_overwrites=False):
        """adds an AEFACT object"""
        key = aefact.sid
        if key in self.aefacts and not allow_overwrites:
            if not aefact._is_same_card(self.aefacts[key]):
                assert key not in self.aefacts, 'AEFACT.sid=%s\nold=\n%snew=\n%s' % (key, self.aefacts[key], aefact)
        else:
            assert key > 0, 'sid=%s method=\n%s' % (key, aefact)
            self.aefacts[key] = aefact
            self._type_to_id_map[aefact.type].append(key)

    def _add_aelist_object(self, aelist):
        """adds an AELIST object"""
        key = aelist.sid
        assert key not in self.aelists, 'AELIST.sid=%s\nold=\n%snew=\n%s' % (key, self.aelists[key], aelist)
        assert key >= 0
        self.aelists[key] = aelist
        self._type_to_id_map[aelist.type].append(key)

    def _add_aelink_object(self, aelink):
        """adds an AELINK object"""
        key = aelink.id
        assert key >= 0
        if key not in self.aelinks:
            self.aelinks[key] = []
        self.aelinks[key].append(aelink)
        self._type_to_id_map[aelink.type].append(key)
        #assert key not in self.aestats,'\naestat=%s oldAESTAT=\n%s' %(aestat,self.aestats[key])

    def _add_aecomp_object(self, aecomp):
        """adds an AECOMP object"""
        key = aecomp.name
        assert key not in self.aecomps, '\naecomp=\n%s oldAECOMP=\n%s' % (aecomp, self.aecomps[key])
        self.aecomps[key] = aecomp
        self._type_to_id_map[aecomp.type].append(key)

    def _add_aeparm_object(self, aeparam):
        """adds an AEPARM object"""
        key = aeparam.id
        assert key not in self.aeparams, '\naeparam=\n%s oldAEPARM=\n%s' % (aeparam, self.aeparams[key])
        assert key >= 0
        self.aeparams[key] = aeparam
        self._type_to_id_map[aeparam.type].append(key)

    def _add_aestat_object(self, aestat):
        """adds an AESTAT object"""
        key = aestat.id
        assert key not in self.aestats, '\naestat=\n%s old=\n%s' % (
            aestat, self.aestats[key])
        assert key >= 0
        self.aestats[key] = aestat
        self._type_to_id_map[aestat.type].append(key)

    def _add_aesurf_object(self, aesurf):
        """adds an AESURF object"""
        key = aesurf.aesid
        assert key not in self.aesurf, '\naesurf=\n%s old=\n%s' % (
            aesurf, self.aesurf[key])
        assert key >= 0
        self.aesurf[key] = aesurf
        self._type_to_id_map[aesurf.type].append(key)

    def _add_aesurfs_object(self, aesurfs):
        """adds an AESURFS object"""
        key = aesurfs.aesid
        assert key not in self.aesurf, '\naesurfs=\n%s old=\n%s' % (
            aesurfs, self.aesurfs[key])
        assert key >= 0
        self.aesurfs[key] = aesurfs
        self._type_to_id_map[aesurfs.type].append(key)

    def _add_csschd_object(self, csschd):
        """adds an CSSCHD object"""
        key = csschd.sid
        assert key not in self.csschds, '\nCSSCHD=\n%s old=\n%s' % (csschd, self.csschds[key])
        assert key >= 0
        self.csschds[key] = csschd
        self._type_to_id_map[csschd.type].append(key)

    def _add_caero_object(self, caero):
        """adds an CAERO1/CAERO2/CAERO3/CAERO4/CAERO5 object"""
        key = caero.eid
        assert key not in self.caeros, '\ncaero=\n|%s| old_caero=\n|%s|' % (
            caero, self.caeros[key])
        assert key > 0
        self.caeros[key] = caero
        self._type_to_id_map[caero.type].append(key)

    def _add_paero_object(self, paero):
        key = paero.pid
        assert key not in self.paeros, '\npaero=\n|%s| old_paero=\n|%s|' % (
            paero, self.paeros[key])
        assert key > 0, 'paero.pid = %r' % (key)
        self.paeros[key] = paero
        self._type_to_id_map[paero.type].append(key)

    def _add_monpnt_object(self, monitor_point):
        """adds an MONPNT object"""
        key = monitor_point.name
        assert key not in self.monitor_points, '\nmonitor_point=\n%soldMOTPNT=\n%s' % (
            monitor_point, self.monitor_points[key])
        self.monitor_points.append(monitor_point)
        self._type_to_id_map[monitor_point.type].append(len(self.monitor_points) - 1)

    def _add_spline_object(self, spline):
        """adds an SPLINE1/SPLINE2/SPLINE3/SPLINE4/SPLINE5 object"""
        assert spline.eid not in self.splines
        assert spline.eid > 0
        key = spline.eid
        self.splines[key] = spline
        self._type_to_id_map[spline.type].append(key)

    def _add_gust_object(self, gust):
        """adds an GUST object"""
        key = gust.sid
        assert key not in self.gusts
        assert key > 0
        self.gusts[key] = gust
        self._type_to_id_map[gust.type].append(key)

    def _add_trim_object(self, trim, allow_overwrites=False):
        """adds an TRIM object"""
        key = trim.sid
        if not allow_overwrites:
            assert key not in self.trims, 'TRIM=%s  old=\n%snew=\n%s' % (key, self.trims[key], trim)
        assert key > 0, 'key=%r trim=\n%s' % (key, trim)
        self.trims[key] = trim
        self._type_to_id_map[trim.type].append(key)

    def _add_diverg_object(self, diverg, allow_overwrites=False):
        """adds an DIVERG object"""
        key = diverg.sid
        if not allow_overwrites:
            assert key not in self.divergs, 'DIVERG=%s  old=\n%snew=\n%s' % (key, self.divergs[key], diverg)
        assert key > 0, 'key=%r diverg=\n%s' % (key, diverg)
        self.divergs[key] = diverg
        self._type_to_id_map[diverg.type].append(key)

    def _add_flutter_object(self, flutter):
        """adds an FLUTTER object"""
        key = flutter.sid
        assert key not in self.flutters, 'FLUTTER=%s old=\n%snew=\n%s' % (key, self.flutters[key], flutter)
        assert key > 0
        self.flutters[key] = flutter
        self._type_to_id_map[flutter.type].append(key)

    def _add_flfact_object(self, flfact):
        """adds an FLFACT object"""
        key = flfact.sid
        #assert key not in self.flfacts
        assert key > 0
        self.flfacts[key] = flfact  # set id...
        self._type_to_id_map[flfact.type].append(key)

    def _add_dconstr_object(self, dconstr):
        """adds a DCONSTR object"""
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

    def _add_desvar_object(self, desvar):
        """adds a DESVAR object"""
        key = desvar.desvar_id
        assert key not in self.desvars, 'DESVAR=%s old=\n%snew=\n%s' % (
            key, self.desvars[key], desvar)
        assert key > 0
        self.desvars[key] = desvar
        self._type_to_id_map[desvar.type].append(key)

    def _add_ddval_object(self, ddval):
        """adds a DDVAL object"""
        key = ddval.oid
        assert key not in self.ddvals, 'DDVAL=%s old=\n%snew=\n%s' % (
            key, self.ddvals[key], ddval)
        assert key > 0
        self.ddvals[key] = ddval
        self._type_to_id_map[ddval.type].append(key)

    def _add_dlink_object(self, dlink):
        """adds a DLINK object"""
        key = dlink.oid
        assert key not in self.dlinks, 'DLINK=%s old=\n%snew=\n%s' % (
            key, self.dlinks[key], dlink)
        assert key > 0
        self.dlinks[key] = dlink
        self._type_to_id_map[dlink.type].append(key)

    def _add_dresp_object(self, dresp):
        """adds a DRESP1/DRESP2/DRESP3 object"""
        key = dresp.dresp_id
        assert key not in self.dresps, 'DRESPx=%s old=\n%snew=\n%s' % (
                    key, self.dresps[key], dresp)
        assert key > 0
        self.dresps[key] = dresp
        self._type_to_id_map[dresp.type].append(key)

    def _add_dvcrel_object(self, dvcrel):
        """adds a DVCREL1/DVCREL2 object"""
        key = dvcrel.oid
        assert key not in self.dvcrels, 'DVCRELx=%s old\n%snew=\n%s' % (
                    key, self.dvcrels[key], dvcrel)
        assert key > 0
        self.dvcrels[key] = dvcrel
        self._type_to_id_map[dvcrel.type].append(key)

    def _add_dvmrel_object(self, dvmrel):
        """adds a DVMREL1/DVMREL2 object"""
        key = dvmrel.oid
        assert key not in self.dvmrels, 'DVMRELx=%s old=\n%snew=\n%s' % (
                    key, self.dvmrels[key], dvmrel)
        assert key not in self.dvmrels
        assert key > 0
        self.dvmrels[key] = dvmrel
        self._type_to_id_map[dvmrel.type].append(key)

    def _add_dvprel_object(self, dvprel):
        """adds a DVPREL1/DVPREL2 object"""
        key = dvprel.oid
        assert key not in self.dvprels, 'DVPRELx=%s old\n%snew=\n%s' % (
                    key, self.dvprels[key], dvprel)
        assert key > 0
        self.dvprels[key] = dvprel
        self._type_to_id_map[dvprel.type].append(key)

    def _add_dvgrid_object(self, dvgrid):
        """adds a DVGRID object"""
        key = dvgrid.dvid
        assert key > 0
        if key not in self.dvgrids:
            self.dvgrids[key] = []
            self._type_to_id_map[dvgrid.type].append(key)
        self.dvgrids[key].append(dvgrid)

    def _add_nlparm_object(self, nlparm):
        """adds a NLPARM object"""
        key = nlparm.nlparm_id
        assert key not in self.nlparms
        assert key > 0, 'key=%s; nlparm=%s\n' % (key, nlparm)
        self.nlparms[key] = nlparm
        self._type_to_id_map[nlparm.type].append(key)

    def _add_rotor_object(self, rotor):
        """adds a ROTORG object"""
        key = rotor.sid
        assert key > 0, 'key=%s; rotor=%s\n' % (key, rotor)
        if key in self.rotors:
            rotor_old = self.rotors[key]
            assert rotor.type == rotor_old.type
            self.rotors[key].nids += rotor.nids
        else:
            self.rotors[key] = rotor
        self._type_to_id_map[rotor.type].append(key)

    def _add_nlpci_object(self, nlpci):
        """adds a NLPCI object"""
        key = nlpci.nlpci_id
        assert key not in self.nlpcis
        assert key > 0
        self.nlpcis[key] = nlpci
        self._type_to_id_map[nlpci.type].append(key)

    def _add_tstep_object(self, tstep, allow_overwrites=False):
        """adds a TSTEP object"""
        key = tstep.sid
        if key in self.tsteps and not allow_overwrites:
            if not tstep._is_same_card(self.tsteps[key]):
                assert key not in self.tsteps, 'TSTEP=%s\nold=\n%snew=\n%s' % (key, self.tsteps[key], tstep)
        else:
            assert key > 0, 'sid=%s tstep=\n%s' % (key, tstep)
            self.tsteps[key] = tstep
            self._type_to_id_map[tstep.type].append(key)

    def _add_tstepnl_object(self, tstepnl, allow_overwrites=False):
        """adds a TSTEPNL object"""
        key = tstepnl.sid
        if key in self.tstepnls and not allow_overwrites:
            if not tstepnl._is_same_card(self.tstepnls[key]):
                assert key not in self.tstepnls, 'TSTEPNL=%s\nold=\n%snew=\n%s' % (key, self.tstepnls[key], tstepnl)
        else:
            assert key > 0, 'sid=%s tstepnl=\n%s' % (key, tstepnl)
            self.tstepnls[key] = tstepnl
            self._type_to_id_map[tstepnl.type].append(key)

    def _add_freq_object(self, freq):
        key = freq.sid
        assert key > 0
        if key in self.frequencies:
            freq0 = self.frequencies[key][0]
            if freq0.type == 'FREQ' and freq.type == 'FREQ':
                freq0.add_frequency_object(freq)
            else:
                self.frequencies[key].append(freq)
        else:
            self.frequencies[key] = [freq]
            self._type_to_id_map[freq.type].append(key)

    def _add_set_object(self, set_obj):
        """adds an SET1/SET3 object"""
        key = set_obj.sid
        assert key >= 0
        if key in self.sets:
            self.sets[key].add_set(set_obj)
        else:
            self.sets[key] = set_obj
            self._type_to_id_map[set_obj.type].append(key)

    def _add_aset_object(self, set_obj):
        """adds an ASET/ASET1 object"""
        self.asets.append(set_obj)

    def _add_bset_object(self, set_obj):
        """adds an BSET/BSET1 object"""
        self.bsets.append(set_obj)

    def _add_cset_object(self, set_obj):
        """adds an CSET/USET1 object"""
        self.csets.append(set_obj)

    def _add_qset_object(self, set_obj):
        """adds an QSET/QSET1 object"""
        self.qsets.append(set_obj)

    def _add_uset_object(self, set_obj):
        """adds an USET/USET1 object"""
        key = set_obj.name
        if key in self.usets:
            self.usets[key].append(set_obj)
        else:
            self.usets[key] = [set_obj]

    def _add_sebset_object(self, set_obj):
        """adds an SEBSET/SEBSET1 object"""
        self.se_bsets.append(set_obj)

    def _add_secset_object(self, set_obj):
        """adds an SECSET/SECSTE1 object"""
        self.se_csets.append(set_obj)

    def _add_seqset_object(self, set_obj):
        """adds an SEQSET/SEQSET1 object"""
        self.se_qsets.append(set_obj)

    def _add_seuset_object(self, set_obj):
        """adds an SEUSET/SEUSET1 object"""
        key = set_obj.name
        if key in self.se_usets:
            self.se_usets[key].append(set_obj)
        else:
            self.se_usets[key] = [set_obj]

    def _add_seset_object(self, set_obj):
        """adds an SESET object"""
        key = set_obj.seid
        assert key >= 0
        if key in self.se_sets:
            old_set = self.se_sets[key]
            set_obj.add_seset(old_set)
        self.se_sets[key] = set_obj
        self._type_to_id_map[set_obj.type].append(key)

    def _add_table_object(self, table):
        """adds a TABLES1, TABLEST object"""
        key = table.tid
        assert key not in self.tables, '\nTable=\n%s old_table=\n%s' % (
            table, self.tables[key])
        assert key > 0
        self.tables[key] = table
        self._type_to_id_map[table.type].append(key)

    def _add_tabled_object(self, table):
        """adds a TABLED1, TABLED2, TABLED3, TABLED4 object"""
        key = table.tid
        assert key not in self.tables_d, '\ntabled=\n%s old_tabled=\n%s' % (
            table, self.tables_d[key])
        assert key > 0
        self.tables_d[key] = table
        self._type_to_id_map[table.type].append(key)

    def _add_tablem_object(self, table):
        """adds a TABLEM1, TABLEM2, TABLEM3, TABLEM4 object"""
        key = table.tid
        assert key not in self.tables_m, '\ntablem=\n%s old_tablem=\n%s' % (
            table, self.tables_m[key])
        assert key > 0
        self.tables_m[key] = table
        self._type_to_id_map[table.type].append(key)

    def _add_table_sdamping_object(self, table):
        """adds a TABDMP1 object"""
        key = table.tid
        assert key not in self.tables_sdamping, '\nTable=\n%s oldTable=\n%s' % (
            table, self.tables_sdamping[key])
        assert key > 0
        self.tables_sdamping[key] = table
        self._type_to_id_map[table.type].append(key)

    def _add_random_table_object(self, table):
        """adds a TABRND1, TABRNDG object"""
        key = table.tid
        assert key not in self.random_tables, '\nTable=\n%s old=\n%s' % (
            table, self.random_tables[key])
        assert key > 0
        self.random_tables[key] = table
        self._type_to_id_map[table.type].append(key)

    def _add_method_object(self, method, allow_overwrites=False):
        """adds a EIGR/EIGRL object"""
        key = method.sid
        if key in self.methods and not allow_overwrites:
            if not method._is_same_card(self.methods[key]):
                assert key not in self.methods, 'sid=%s\nold_method=\n%snew_method=\n%s' % (key, self.methods[key], method)
        else:
            assert key > 0, 'sid=%s method=\n%s' % (key, method)
            self.methods[key] = method
            self._type_to_id_map[method.type].append(key)

    def _add_cmethod_object(self, method, allow_overwrites=False):
        """adds a EIGB/EIGC object"""
        key = method.sid
        if key in self.cMethods and not allow_overwrites:
            if not method._is_same_card(self.cMethods[key]):
                assert key not in self.cMethods, 'sid=%s\nold_cmethod=\n%snew_cmethod=\n%s' % (key, self.cMethods[key], method)
        else:
            assert key > 0, 'sid=%s cMethod=\n%s' % (key, method)
            self.cMethods[key] = method
            self._type_to_id_map[method.type].append(key)

    def _add_mkaero_object(self, mkaero):
        """adds an MKAERO1/MKAERO2 object"""
        self.mkaeros.append(mkaero)
