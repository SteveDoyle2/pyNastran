from pyNastran.dev.bdf_vectorized.bdf_interface2.attributes import BDFAttributes


class AddCard(BDFAttributes):
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
        self.dmig[name] = dmig
        self._type_to_id_map[dmig.type].append(name)

    def _add_dmij_object(self, dmij, allow_overwrites=False):
        """adds a DMIJ object"""
        name = dmij.name
        self.dmij[name] = dmij
        self._type_to_id_map[dmij.type].append(name)

    def _add_dmiji_object(self, dmiji, allow_overwrites=False):
        """adds a DMIJI object"""
        name = dmiji.name
        self.dmiji[name] = dmiji
        self._type_to_id_map[dmiji.type].append(name)

    def _add_dmik_object(self, dmik, allow_overwrites=False):
        """adds a DMIK object"""
        name = dmik.name
        self.dmiks[name] = dmik
        self._type_to_id_map[dmik.type].append(name)

    def _add_param_object(self, param, allow_overwrites=False):
        """adds a PARAM object"""
        key = param.key
        if key in self.params and not allow_overwrites:
            if not param == self.params[key]:
                #if param.key in self.params:
                    #msg = 'key=%s param=%s old_param=%s' % (key, param, self.params[key])
                    #raise KeyError(msg)
                self.log.warning('key=%s param=%s old_param=%s' %
                                 (key, param, self.params[key]))
                self.params[key] = param
        else:
            self.params[key] = param
            self._type_to_id_map[param.type].append(key)

    def _add_plotel_object(self, elem, allow_overwrites=False):
        """adds an PLOTEL object"""
        key = elem.eid
        assert key > 0, 'eid=%s must be positive; elem=\n%s' % (key, elem)
        if not allow_overwrites:
            #if key in self.elements:
                #if elem ==self.elements[key]:
                    #self._duplicate_elements.append(elem)
                    #if self._stop_on_duplicate_error:
                        #self.pop_parse_errors()
            if key in self.plotels:
                if not elem == self.plotels[key]:
                    assert elem.eid not in self.plotels, 'eid=%s\nold_element=\n%snew_element=\n%s' % (elem.eid, self.plotels[elem.eid], elem)
        self.plotels[key] = elem
        self._type_to_id_map[elem.type].append(key)

    def _add_doptprm_object(self, doptprm, comment=''):
        """adds a DOPTPRM"""
        self.doptprm = doptprm


    def _add_rigid_element_object(self, elem, allow_overwrites=False):
        key = elem.eid
        assert key > 0, 'eid=%s elem=%s' % (key, elem)
        if key in self.rigid_elements and not allow_overwrites:
            assert elem.eid not in self.rigid_elements, 'eid=%s\noldElement=\n%snewElement=\n%s' % (elem.eid, self.rigid_elements[elem.eid], elem)
        self.rigid_elements[key] = elem
        self._type_to_id_map[elem.type].append(key)

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
            if not aefact ==self.aefacts[key]:
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
        key = aeparam.aeparm_id
        assert key not in self.aeparams, '\naeparam=\n%s oldAEPARM=\n%s' % (aeparam, self.aeparams[key])
        assert key >= 0
        self.aeparams[key] = aeparam
        self._type_to_id_map[aeparam.type].append(key)

    def _add_aestat_object(self, aestat):
        """adds an AESTAT object"""
        key = aestat.aestat_id
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
        assert key not in self.monitor_points, '\nmonitor_point=\n%soldMNTPNT=\n%s' % (
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
            if not tstep == self.tsteps[key]:
                assert key not in self.tsteps, 'TSTEP=%s\nold=\n%snew=\n%s' % (key, self.tsteps[key], tstep)
        else:
            assert key > 0, 'sid=%s tstep=\n%s' % (key, tstep)
            self.tsteps[key] = tstep
            self._type_to_id_map[tstep.type].append(key)

    def _add_tstepnl_object(self, tstepnl, allow_overwrites=False):
        """adds a TSTEPNL object"""
        key = tstepnl.sid
        if key in self.tstepnls and not allow_overwrites:
            if not tstepnl == self.tstepnls[key]:
                assert key not in self.tstepnls, 'TSTEPNL=%s\nold=\n%snew=\n%s' % (key, self.tstepnls[key], tstepnl)
        else:
            assert key > 0, 'sid=%s tstepnl=\n%s' % (key, tstepnl)
            self.tstepnls[key] = tstepnl
            self._type_to_id_map[tstepnl.type].append(key)


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

    def _add_omit_object(self, set_obj):
        """adds an OMIT/OMIT1 object"""
        self.omits.append(set_obj)

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

    def _add_method_object(self, method, allow_overwrites=False):
        """adds a EIGR/EIGRL object"""
        key = method.sid
        if key in self.methods and not allow_overwrites:
            if not method == self.methods[key]:
                assert key not in self.methods, 'sid=%s\nold_method=\n%snew_method=\n%s' % (key, self.methods[key], method)
        else:
            assert key > 0, 'sid=%s method=\n%s' % (key, method)
            self.methods[key] = method
            self._type_to_id_map[method.type].append(key)

    def _add_cmethod_object(self, method, allow_overwrites=False):
        """adds a EIGB/EIGC object"""
        key = method.sid
        if key in self.cMethods and not allow_overwrites:
            if not method == self.cMethods[key]:
                assert key not in self.cMethods, 'sid=%s\nold_cmethod=\n%snew_cmethod=\n%s' % (key, self.cMethods[key], method)
        else:
            assert key > 0, 'sid=%s cMethod=\n%s' % (key, method)
            self.cMethods[key] = method
            self._type_to_id_map[method.type].append(key)

    def _add_mkaero_object(self, mkaero):
        """adds an MKAERO1/MKAERO2 object"""
        self.mkaeros.append(mkaero)
