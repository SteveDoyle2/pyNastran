"""
defines some methods for cleaning up a model
 - model = remove_unused(bdf_filename, remove_nids=True, remove_cids=True,
                         remove_pids=True, remove_mids=True)

"""
from pyNastran.bdf.bdf import BDF, read_bdf
#from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber

def remove_unused(bdf_filename, remove_nids=True, remove_cids=True,
                  remove_pids=True, remove_mids=True, remove_spcs=True, remove_mpcs=True):
    """
    Takes an uncross-referenced bdf and removes unused data

    removes unused:
     - nodes
     - properties
     - materials
     - coords
     - spcs
     - mpcs

    """
    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        model = read_bdf(bdf_filename, xref=False)

    #nids = model.nodes.keys()
    #cids =
    #nids = set(list(model.nodes.keys()))
    #cids = set(list(model.coords.keys()))
    #pids = set(list(model.properties.keys()))

    nids_used = set()
    cids_used = set()
    eids_used = set()
    pids_used = set()
    pids_mass_used = set()
    mids_used = set()
    mids_thermal_used = set()
    sets_used = set()
    desvars_used = set()
    mpcs_used = set()
    spcs_used = set()
    #nsms_used = set()

    #card_types = list(model.card_count.keys())
    #card_map = model.get_card_ids_by_card_types(
        #card_types=card_types,
        #reset_type_to_slot_map=False,
        #stop_on_missing_card=True)

    #for nid, node in model.nodes.items():
        #cids_used.update([node.Cp(), node.Cd()])

    # ureferenced types aren't referenced by anything
    unreferenced_types = {
        'ENDDATA', 'PARAM', 'EIGR', 'EIGRL', 'EIGB', 'EIGP', 'EIGC',
        'SPOINT', 'EPOINT', 'DESVAR',
        'SET1', 'FREQ', 'FREQ1', 'FREQ2',
        'TSTEP', 'TSTEPNL', 'NLPCI',
        'NLPARM', 'ROTORG', 'ROTORD',
        'DAREA', 'DEQATN',
        'DMIG', 'DMI', 'DMIJ', 'DMIK', 'DMIJI', 'DMIAX',
        'POINT', 'EPOINT',
        'DELAY', 'DPHASE',
        'CBARAO', 'AEPARM',

        # properties
        'PELAS', 'PDAMP', 'PBUSH',
        'PELAST', 'PDAMPT', 'PBUSHT',
        'PGAP', 'PBUSH1D', 'PFAST', 'PVISC', 'PMASS',

        'FLFACT', 'FLUTTER', 'DLINK', 'DDVAL', 'DIVERG', 'GUST',
        'AELINK', 'AELIST', 'TRIM', 'TRIM2', 'PAERO1', 'AEFACT', 'AESTAT',

        # contact
        'BCTPARA', 'BCRPARA', 'BSURF', 'BSURFS', 'BCTADD',
        'BCTSET', 'BFRIC',

        'TABRNDG', 'DTI', 'TABLEH1',

    }

    # this are things that haven't been referenced yet
    not_implemented_types = {
        # not checked------------------------------------------
        'PHBDY', 'CHBDYG', 'CHBDYP', 'CHBDYE', 'RADBC', # 'CONV',
        'QVOL', 'PCONVM', # 'PCONV',
        #'PBCOMP', 'PDAMP5', 'CFAST',
        'AECOMP', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5',
        'PAERO2', 'PAERO3', 'PAERO4', 'PAERO5',
        'DCONADD',
        'GMCORD',
        'MONPNT1', 'MONPNT2', 'MONPNT3',
        'DSCREEN', 'DTI', 'NSMADD',
        'AESURFS', 'CSSCHD',
        'CGEN', 'NXSTRAT',

        # axisymmetric
        'FORCEAX',

        # acoustic
        'PACABS',

        # superelements
        'SELOC',

        # parametric
        'FEEDGE', 'FEFACE'
    }
    set_types_simple = [
        'SET1', 'SET3',
    ]
    set_types = {
        'ASET', 'ASET1', 'BSET', 'BSET1', 'CSET', 'CSET1',
        'QSET', 'QSET1', 'USET', 'USET1', 'OMIT', 'OMIT1',
    }
    #seset_types = {
        #'SESET',
    #}
    load_types = {
        'GRAV', 'RANDPS', 'FORCE', 'FORCE1', 'FORCE2',
        'MOMENT', 'MOMENT1', 'MOMENT2',
        'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4', 'SPCD',
        'GMLOAD', 'RFORCE', 'RFORCE1',
        'TEMP', 'QBDY1', 'QBDY2', 'QBDY3', 'QHBDY',
        'ACCEL', 'PLOADX1', 'SLOAD', 'ACCEL1', 'LOADCYN', 'LOAD', 'CLOAD',
        'LSEQ', 'DLOAD', 'QVECT', 'RADM', 'TEMPAX', 'DEFORM',
    }
    masses = {'CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4'}
    elements = {
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
        'CGAP', 'CBUSH', 'CBUSH1D', 'CBUSH2D', 'CFAST',
        'CVISC',
        'CSHEAR', 'CTUBE',
        'GENEL',
        'CTRIA3', 'CQUAD4', 'CTRIA6', 'CTRIAR', 'CQUAD8', 'CQUADR',
        'CTRIAX', 'CQUADX', 'CQUAD',

        'CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8',
        'CPLSTS3', 'CPLSTS4', 'CPLSTS6', 'CPLSTS8',
        'CQUADX4', 'CQUADX8', 'CTRIAX6',
        'CTRAX3', 'CTRAX6', 'CTRIAX6',
        'CTETRA', 'CHEXA', 'CPENTA', 'CPYRAM',
        'CROD', 'CRAC2D', 'CRAC3D',
        'CONROD', 'CCONEAX',
        'CBAR', 'CBEAM', 'CBEND', 'CBEAM3',

        # acoustic
        'CHACAB',
    }
    tableht_used = set([])
    tableh1_used = set([])
    pconv_used = set([])

    friction_ids_used = set([])

    # could remove some if we look at the rid_trace
    #for cid, coord in model.coords.items():
        #if coord.type in ['CORD1R', 'CORD1C', 'CORD1S']:
            #nids_used.update(node_ids)
        #elif coord.type in ['CORD1R', 'CORD1C', 'CORD1S']:
            #cids_used.update(coord.Rid())
        #else:
            #raise NotImplementedError(coord)

    for card_type, ids in model._type_to_id_map.items():
    #for card_type, ids in card_map.items():
        if card_type in ['CORD1R', 'CORD1C', 'CORD1S']:
            for cid in ids:
                coord = model.coords[cid]
                nids_used.update(coord.node_ids)
        elif card_type in ['CORD2R', 'CORD2C', 'CORD2S']:
            for cid in ids:
                coord = model.coords[cid]
                cids_used.add(coord.Rid())

        elif card_type in ['MAT1', 'MAT2', 'MAT3', 'MAT4', 'MAT5',
                           'MAT8', 'MAT9', 'MAT10', 'MAT11']:
            # todo: MATS1, MATT1, etc.
            pass
        elif card_type in ['MATS1', 'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5',
                           'MATT8', 'MATT9', 'MATHE', 'MATHP', 'CREEP']:
            mids_used.update(ids)

        elif card_type in masses:
            _store_masses(card_type, model, ids, nids_used, pids_mass_used, cids_used)
        elif card_type in elements:
            _store_elements(card_type, model, ids, nids_used, pids_used, mids_used, cids_used)

        elif card_type == 'PLPLANE':
            for pid in ids:
                prop = model.properties[pid]
                cids_used.add(prop.cid)
                mids_used.add(prop.Mid())
        elif card_type == 'PPLANE':
            for pid in ids:
                prop = model.properties[pid]
                mids_used.add(prop.Mid())

        elif card_type == 'PLOTEL':
            for eid in ids:
                elem = model.plotels[eid]
                nids_used.update(elem.node_ids)
        elif card_type in ['PSOLID', 'PLSOLID', 'PIHEX']:
            for pid in ids:
                prop = model.properties[pid]
                mids_used.add(prop.Mid())
        elif card_type == 'PDAMP5':
            for pid in ids:
                prop = model.properties[pid]
                mids_thermal_used.add(prop.Mid())
        elif card_type in ['PBAR', 'PBARL', 'PROD', 'PTUBE', 'PBEAM', 'PBEAML', 'PBEAM3',
                           'PSHEAR', 'PRAC2D', 'PRAC3D', 'PBEND']:
            for pid in ids:
                prop = model.properties[pid]
                mids_used.add(prop.Mid())
        elif card_type == 'PSHELL':
            for pid in ids:
                prop = model.properties[pid]
                mids = [mid for mid in prop.material_ids if mid is not None]
                mids_used.update(mids)
        elif card_type in ['PCOMP', 'PCOMPG']:
            for pid in ids:
                prop = model.properties[pid]
                mids = prop.material_ids
                mids_used.update(mids)
        elif card_type in ['PBCOMP']:
            for pid in ids:
                prop = model.properties[pid]
                mids = prop.Mids()
                mids_used.add(prop.Mid())
                mids_used.update(mids)
        elif card_type == 'PCOMPS':
            for pid in ids:
                prop = model.properties[pid]
                mids = prop.Mids()
                mids_used.update(mids)
                cids_used.add(prop.cordm)

        elif card_type == 'PCONEAX':
            for pid in ids:
                # MID1 T1 MID2 I MID3 T2 NSM
                prop = model.properties[pid]
                #print(prop.object_methods())
                mids = [mid for mid in prop.Mids() if mid not in (0, None)]
                prop = model.properties[pid]
                mids_used.update(mids)

        elif card_type in ['RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3', 'RROD', 'RSPLINE', 'RSSCON']:
            for eid in ids:
                elem = model.rigid_elements[eid]
                #print(elem.object_attributes())
                #print(elem.object_methods())
                nids_used.update(elem.independent_nodes)
                nids_used.update(elem.dependent_nodes)

        elif card_type in ['TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2', 'ACSRCE']:
            pass
        elif card_type in load_types:
            _store_loads(model, card_type, ids, nids_used, eids_used, cids_used)

        elif card_type == 'TEMPD':
            pass
            #for temp_id in ids:
                #tempd = self.tempds[temp_id]

        elif card_type == 'MPCADD':
            for mpcadds in model.mpcadds.values():
                for mpcadd in mpcadds:
                    mpcs_used.update(mpcadd.mpc_ids)
        elif card_type == 'MPC':
            for mpcs in model.mpcs.values():
                for mpc in mpcs:
                    nids_used.update(mpc.node_ids)

        elif card_type == 'SPCADD':
            for spcadds in model.spcadds.values():
                for spcadd in spcadds:
                    spcs_used.update(spcadd.spc_ids)
        elif card_type in ['SPC1', 'SPC', 'GMSPC', 'SPCAX']:
            for spcs in model.spcs.values():
                for spc in spcs:
                    if spc.type in ['GMSPC', 'SPCAX']:
                        pass
                    elif spc.type in ['SPC1', 'SPC']:
                        nids_used.update(spc.node_ids)
                    else:
                        raise NotImplementedError(spc)

        elif card_type in ['TABLED1', 'TABLED2', 'TABLED3', 'TABLED4',
                           'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',
                           'TABDMP1', 'TABRND1', 'TABLES1',]:
            pass
        elif card_type == 'SUPORT':
            for suport in model.suport:
                nids_used.update(suport.node_ids)
        elif card_type == 'SUPORT1':
            for suport1 in model.suport1.values():
                nids_used.update(suport1.node_ids)
        elif card_type == 'TIC':
            for tic in model.tics.values():
                nids_used.update(tic.node_ids)
        elif card_type == 'GRID':
            for unused_nid, node in model.nodes.items():
                cids_used.update([node.Cp(), node.Cd()])

        elif card_type in ['PBUSH']:
            pass
            #for pid in ids:
                #prop = model.properties[pid]
                #raise RuntimeError(prop)
        elif card_type == 'PBUSHT':
            # tables
            pass
        elif card_type == 'AESURF':
            #CID1  | ALID1 | CID2   | ALID2
            for aesurf in model.aesurf.values():
                cids_used.add(aesurf.Cid1())
                cid2 = aesurf.Cid2()
                if cid2 is not None:
                    cids_used.add(cid2)

        elif card_type in ['SPLINE1', 'SPLINE2', 'SPLINE3', 'SPLINE4', 'SPLINE5']:
            pass
            #for spline_id in ids:
                #spline = model.splines[spline_id]
                #if card_type in ['SPLINE1', 'SPLINE2', 'SPLINE4', 'SPLINE5']:
                    #sets_used.add(spline.Set())

        elif card_type == 'CAERO1':
            for eid in ids:
                caero = model.caeros[eid]
                # PID, LSPAN, LCHORD
                cids_used.add(caero.Cp())

        elif card_type in unreferenced_types:
            pass
        elif card_type in set_types_simple:
            # handled based on context in other blocks
            pass
        elif card_type in ['USET', 'USET1']:
            for set_cards in model.usets.values():
                for set_card in set_cards:
                    nids_used.update(set_card.ids)

        elif card_type in set_types:
            obj = card_type[:4].lower() + 's'
            sets = getattr(model, obj) # list of SETs
            for set_card in sets:
                nids_used.update(set_card.ids)

        elif card_type == 'SESET':
            sets = model.se_sets # list of SETs
            for unused_id, set_card in sorted(sets.items()):
                nids_used.update(set_card.ids)

        elif card_type in ['DCONSTR']:
            pass
        elif card_type == 'DRESP1':
            _store_dresp1(model, ids, nids_used, pids_used)
        elif card_type == 'DRESP2':
            pass
            #for dresp_id in ids:
                #dresp = model.dresps[dresp_id]
                #dresp.deqatn
                #if dresp.property_type in ['PSHELL', 'PCOMP', 'PBAR', 'PBARL', 'PBEAM', 'PROD']:
                    #pids_used.update(dresp.atti_values())
                #elif dresp.property_type is None:
                    #if dresp.response_type in ['WEIGHT', 'EIGN', 'VOLUME']:
                        #pass
                    #elif dresp.response_type in ['DISP']:
                        #nids_used.update(dresp.atti)
                    #else:
                        #msg = str(dresp) + 'response_type=%r' % dresp.response_type
                        #raise NotImplementedError(msg)
                #else:
                #raise NotImplementedError(dresp)
                #msg = str(dresp) + 'response_type=%r' % dresp.response_type
                #raise NotImplementedError(msg)
        elif card_type == 'DRESP3':
            pass

        elif card_type in ['DVPREL1', 'DVPREL2']:
            for dvprel_id in ids:
                dvprel = model.dvprels[dvprel_id]
                desvars_used.update(dvprel.desvar_ids)
                if dvprel.prop_type in ['PSHELL', 'PCOMP', 'PBAR', 'PBARL', 'PBEAM',
                                        'PROD', 'PELAS', 'PBUSH', 'PDAMP', 'PTUBE',
                                        'PSHEAR', 'PDAMP', 'PMASS', 'PBEAML', 'PCOMPG',
                                        'PVISC', 'PBUSHT', 'PELAST', 'PBUSH1D', 'PGAP']:
                    pids_used.add(dvprel.Pid())
                elif dvprel.prop_type in ['DISP']:
                    msg = str(dvprel) + 'dvprel.prop_type=%r' % dvprel.prop_type
                    raise NotImplementedError(msg)
                else:
                    raise NotImplementedError(dvprel)

        elif card_type in ['DVCREL1', 'DVCREL2']:
            for dvcrel_id in ids:
                dvcrel = model.dvcrels[dvcrel_id]
                desvars_used.update(dvcrel.desvar_ids)
                if dvcrel.element_type in ['CMASS2', 'CMASS4', 'CONM1', 'CONM2',
                                           'CELAS2', 'CELAS4', 'CBUSH',
                                           'CDAMP2', 'CQUAD4', 'CGAP', 'CBAR']:
                    #eids_used.add(dvcrel.Eid())  # we don't remove elements...for now
                    pass
                else:
                    msg = str(dvcrel) + 'element_type=%r' % dvcrel.element_type
                    raise NotImplementedError(msg)

        elif card_type in ['DVMREL1', 'DVMREL2']:
            for dvmrel_id in ids:
                dvmrel = model.dvmrels[dvmrel_id]
                desvars_used.update(dvmrel.desvar_ids)
                if dvmrel.mat_type in ['MAT1', 'MAT2', 'MAT8', 'MAT9', 'MAT10', 'MAT11']:
                    mids_used.add(dvmrel.Mid())
                else:
                    msg = str(dvmrel) + 'mat_type=%r' % dvmrel.mat_type
                    raise NotImplementedError(msg)
        elif card_type == 'DVGRID':
            for dvgrid_id in ids:
                dvgrids = model.dvgrids[dvgrid_id]
                for dvgrid in dvgrids:
                    desvars_used.add(dvgrid.desvar_id)
                    nids_used.add(dvgrid.nid)
                    cids_used.add(dvgrid.cid)
        elif card_type == 'TF':
            for tf_id in ids:
                tfs = model.transfer_functions[tf_id]
                for transfer_function in tfs:
                    nids_used.update(transfer_function.nids)
        elif card_type in ['NSM', 'NSM1', 'NSML', 'NSML1']:
            _store_nsm(model, ids, pids_used)

        elif card_type in ['POINTAX', 'AXIC', 'RINGAX']:
            pass
            #for eid in ids:
                #elem = model.plotels[eid]
                #nids_used.update(elem.node_ids)
        elif card_type in ['PBRSECT', 'PBMSECT']:
            for pid in ids:
                prop = model.properties[pid]
                if prop.outp:
                    sets_used.add(prop.outp)
                if prop.brps:
                    for unused_key, value in prop.brps.items():
                        sets_used.add(value)
                #if prop.cores:
                    #for key, value in prop.cores.items():
                        #pids_used.add(value)
        elif card_type == 'CYJOIN':
            for idi in ids:
                cyjoin = model.cyjoin[idi]
                nids_used.update(cyjoin.nids)
        elif card_type == 'TABLEHT':
            for idi in ids:
                table = model.tables[idi]
                tableh1_ids = table.y.tolist()
                tableh1_used.update(tableh1_ids)
                del tableh1_ids
        elif card_type == 'PCONV':
            for idi in ids:
                pconv = model.convection_properties[idi]
                if pconv.tid is not None:
                    tableht_used.add(pconv.tid)
        elif card_type == 'CONV':
            for idi in ids:
                bcs = model.bcs[idi]
                for conv in bcs:
                    if conv.type != 'CONV':
                        continue
                    pconv_used.add(conv.pconid)
        elif card_type == 'BLSEG':
            for idi in ids:
                blseg = model.blseg[idi]
                # line_id
                nids_used.update(blseg.nodes)
                #print(blseg.get_stats())
        elif card_type == 'BCONP':
            for idi in ids:
                bconp = model.bconp[idi]
                # master
                # slave
                # ???
                #print(bconp.get_stats())
                cids_used.add(bconp.cid)
                friction_ids_used.add(bconp.friction_id)
        #elif card_type == 'FORCEAX':
            #pass
            #ring_id
            #sid
        elif card_type in not_implemented_types:
            model.log.warning(f'skipping {card_type}')
        else:
            raise NotImplementedError(card_type)


    #for pid, prop in model.properties.items():
        #prop = model.properties[pid]
        #if prop.type in no_materials:
            #continue
        #elif prop.type == 'PSHELL':
            #mids_used.extend([mid for mid in prop.material_ids if mid is not None])
        #elif prop.type == 'PCONEAX':
            #mids_used.extend([mid for mid in model.Mids() if mid is not None])

        #elif prop.type in prop_mid:
            #mids_used.append(prop.Mid())
        #elif prop.type in ['PCOMP', 'PCOMPG', 'PCOMPS']:
            #mids_used.extend(prop.Mids())

        #elif prop.type == 'PBCOMP':
            #mids_used.append(prop.Mid())
            #mids_used.extend(prop.Mids())
        #else:
            #raise NotImplementedError(prop)
    remove_desvars = False
    _remove(
        model,
        nids_used, cids_used,
        pids_used, pids_mass_used,
        mids_used,
        spcs_used, mpcs_used,
        pconv_used, tableht_used, tableh1_used,
        desvars_used,
        remove_nids=remove_nids,
        remove_cids=remove_cids,
        remove_pids=remove_pids,
        remove_mids=remove_mids,
        remove_spcs=remove_spcs, remove_mpcs=remove_mpcs,
        unused_remove_desvars=remove_desvars,
    )


def _store_elements(card_type, model, ids, nids_used, pids_used, mids_used, cids_used):
    if card_type in ['CTETRA', 'CPENTA', 'CPYRAM', 'CHEXA', 'CHACAB']:
        for eid in ids:
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)
            pids_used.add(elem.Pid())
    elif card_type in ['CELAS1', 'CDAMP1', 'CVISC', 'CDAMP5']:
        for eid in ids:
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)
            pids_used.add(elem.Pid())
    elif card_type in ['CELAS2', 'CDAMP2']:
        for eid in ids:
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)
    elif card_type in ['CELAS3', 'CDAMP3']:
        for eid in ids:
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)
            pids_used.add(elem.Pid())
    elif card_type in ['CELAS4', 'CDAMP4', 'GENEL']:
        for eid in ids:
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)

    elif card_type in ['CTRIA3', 'CQUAD4', 'CTRIA6', 'CTRIAR', 'CQUAD8', 'CQUADR',
                       'CTRIAX', 'CQUADX', 'CQUAD']:
        for eid in ids:
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)
            pids_used.add(elem.Pid())
            if isinstance(elem.theta_mcid, int):
                cids_used.add(elem.theta_mcid)
    elif card_type in ['CTRIAX6', ]:
        for eid in ids:
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)
            mids_used.add(elem.Mid())
    elif card_type in ['CSHEAR', 'CTUBE']:
        for eid in ids:
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)
            pids_used.add(elem.Pid())
    elif card_type in ['CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8',
                       'CPLSTS3', 'CPLSTS4', 'CPLSTS6', 'CPLSTS8',
                       'CQUADX4', 'CQUADX8', 'CTRIAX6',
                       'CTRAX3', 'CTRAX6']:
        for eid in ids:
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)
            pids_used.add(elem.Pid())
    elif card_type in ['CROD', 'CRAC2D', 'CRAC3D']:
        for eid in ids:
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)
            pids_used.add(elem.Pid())
    elif card_type in ['CONROD']:
        for eid in ids:
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)
            pids_used.add(elem.Mid())
    elif card_type == 'CCONEAX':
        for eid in ids:
            elem = model.elements[eid]
            pids_used.add(elem.Pid())
    elif card_type in ['CBAR', 'CBEAM', 'CBEND']:
        for eid in ids:
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)
            pids_used.add(elem.Pid())
            if elem.g0 is not None:
                assert isinstance(elem.g0, int), elem.g0
                nids_used.add(elem.g0)
    elif card_type == 'CBEAM3':
        for eid in ids:
            elem = model.elements[eid]
            nids_used.add(elem.Ga())
            nids_used.add(elem.Gb())
            if elem.gc is not None:
                nids_used.add(elem.gc)
            pids_used.add(elem.Pid())
            if elem.g0 is not None:
                assert isinstance(elem.g0, int), elem.g0

    elif card_type == 'CFAST':
        for eid in ids:
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)
            pids_used.add(elem.Pid())
    elif card_type == 'CGAP':
        for eid in ids:
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)
            pids_used.add(elem.Pid())
            if elem.g0 is not None:
                assert isinstance(elem.g0, int), elem.g0
                nids_used.add(elem.G0())
    elif card_type == 'CBUSH':
        for eid in ids:
            elem = model.elements[eid]
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)
            if elem.g0 is not None:
                assert isinstance(elem.g0, int), elem.g0
                nids_used.add(elem.G0())
            pids_used.add(elem.Pid())
            if elem.cid is not None:
                cids_used.add(elem.Cid())
    elif card_type in ['CBUSH1D', 'CBUSH2D']:
        for eid in ids:
            elem = model.elements[eid]
            nids_used.update(elem.node_ids)
            pids_used.add(elem.Pid())
            cids_used.add(elem.Cid())
    else:
        raise NotImplementedError(card_type)

def _store_nsm(model, ids, pids_used):
    """helper for ``remove_unused``"""
    for nsm_id in ids:
        nsms = model.nsms[nsm_id]
        for nsm in nsms:
            idsi = nsm.ids

            if nsm.nsm_type in ['PROD', 'PBARL', 'PBEAML',
                                'PSHELL', 'PCOMP', ]:
                if len(idsi) == 1 and idsi[0] == 'ALL':
                    idsi = list(model.properties.keys())
                    #raise NotImplementedError('found ALL...\n%s' % str(nsm))
                pids_used.update(idsi)
            elif nsm.nsm_type in ['CONROD', 'ELEMENT']:
                # we skip this because we assume all elements are used
                #if len(idsi) == 1 and idsi[0] == 'ALL':
                    #raise NotImplementedError('found ALL...\n%s' % str(nsm))
                #eids_used.update(idsi)
                pass
            else:
                msg = 'found nsm_type=%r...\n%s' % (nsm.nsm_type, str(nsm))
                raise NotImplementedError(msg)

def _store_loads(model, unused_card_type, unused_ids, nids_used, eids_used, cids_used):
    """helper for ``remove_unused``"""
    for loads in model.loads.values():
        for load in loads:
            if load.type in ['FORCE', 'MOMENT']:
                nids_used.add(load.node_id)
                cids_used.add(load.Cid())
            elif load.type in ['FORCE1', 'FORCE2', 'MOMENT1', 'MOMENT2']:
                nids_used.update(load.node_ids)
            elif load.type == 'GRAV':
                cids_used.add(load.Cid())
            elif load.type == 'RANDPS':
                pass
            elif load.type == 'PLOAD':
                nids_used.update(load.node_ids)
            elif load.type == 'PLOAD1':
                #eid = integer(card, 2, 'eid')
                pass
            elif load.type == 'PLOAD2':
                #eids_used.update(load.element_ids)
                pass
            elif load.type == 'PLOAD4':
                # eids, g1, g34
                cids_used.add(load.Cid())
            elif load.type == 'DEFORM':
                eids_used.add(load.Eid())
            elif load.type == 'SPCD':
                nids_used.update(load.node_ids)
            elif load.type == 'GMLOAD':
                cids_used.add(load.Cid())
            elif load.type in ['RFORCE', 'RFORCE1']:
                nids_used.add(load.node_id)
                cids_used.add(load.Cid())
            elif load.type == 'TEMP':
                nids_used.update(list(load.temperatures.keys()))
            elif load.type == 'ACCEL':
                # nids?
                cids_used.add(load.Cid())
            elif load.type == 'ACCEL1':
                # nids?
                cids_used.add(load.Cid())
            elif load.type in ['QBDY1', 'QBDY2', 'QBDY3', 'QHBDY']:
                pass
            #'QBDY1', 'QBDY2', 'QBDY3', 'QHBDY', 'PLOADX1
            elif load.type in ['PLOADX1']:
                nids_used.update(load.node_ids)
            elif load.type in ['SLOAD']:
                nids_used.update(load.node_ids)
            elif load.type in ['LOAD', 'LSEQ', 'LOADCYN']:
                pass
            elif load.type in ['QVOL', 'TEMPRB']:
                # eids
                pass
            elif load.type in ['TEMPAX']:
                pass # not done...
            else:
                raise NotImplementedError(load)

def _store_dresp1(model, ids, nids_used, pids_used):
    """helper for ``remove_unused``"""
    for dresp_id in ids:
        dresp = model.dresps[dresp_id]
        if dresp.property_type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PBAR', 'PBARL', 'PBEAM',
                                   'PROD', 'PDAMP', 'PVISC', 'PTUBE', 'PSHEAR', 'PELAS',
                                   'PSOLID', 'PBEAML']:
            pids_used.update(dresp.atti_values())
        elif dresp.property_type == 'ELEM':
            if dresp.response_type in ['STRESS', 'FRSTRE',
                                       'CFAILURE',
                                       'TFORC', 'FRFORC']:
                #eids_used.update(dresp.atti_values())
                pass
            else:
                msg = (
                    str(dresp) + 'region=%r property_type=%r response_type=%r, '
                    'atta=%r attb=%s atti=%s' % (
                        dresp.region, dresp.property_type, dresp.response_type,
                        dresp.atta, dresp.attb, dresp.atti))
                raise NotImplementedError(msg)

        #elif dresp.property_type == 'STRESS':

        elif dresp.property_type is None:
            if dresp.response_type in ['WEIGHT', 'EIGN', 'VOLUME', 'LAMA', 'CEIG',
                                       'FREQ', 'STABDER']:
                pass
            elif dresp.response_type in ['DISP', 'FRDISP', 'TDISP', 'RMSDISP', 'PSDDISP',
                                         'TVELO', 'FRVELO', 'RMSVELO',
                                         'TACCL', 'FRACCL', 'RMSACCL',
                                         'SPCFORCE', 'TSPCF', 'FRSPCF',
                                         'FORCE', 'TFORC', 'FRFORC']:
                nids_used.update(dresp.atti)
            elif dresp.response_type in ['FLUTTER', 'TRIM', 'DIVERG']:
                # flutter_id / trim_id
                pass
            else:
                msg = (
                    str(dresp) + 'region=%r property_type=%r response_type=%r '
                    'atta=%r attb=%s atti=%s' % (
                        dresp.region, dresp.property_type, dresp.response_type,
                        dresp.atta, dresp.attb, dresp.atti))
                raise NotImplementedError(msg)
        else:
            msg = (
                str(dresp) + 'region=%r property_type=%r response_type=%r '
                'atta=%r attb=%s atti=%s' % (
                    dresp.region, dresp.property_type, dresp.response_type,
                    dresp.atta, dresp.attb, dresp.atti))
            raise NotImplementedError(msg)

def _store_masses(card_type, model, ids, nids_used, pids_mass_used, cids_used) -> None:
    """handles masses"""
    if card_type in ['CONM1', 'CONM2']:
        for eid in ids:
            elem = model.masses[eid]
            nids_used.add(elem.Nid())
            cids_used.add(elem.Cid())
    elif card_type in ['CMASS1', 'CMASS3']:
        for eid in ids:
            elem = model.masses[eid]
            pids_mass_used.add(elem.Pid())
            nids_used.update(elem.node_ids)
    elif card_type in ['CMASS2', 'CMASS4']:
        for eid in ids:
            elem = model.masses[eid]
            nids_used.update(elem.node_ids)
    else:
        raise NotImplementedError(card_type)

def _remove(model: BDF,
            nids_used, cids_used,
            pids_used, pids_mass_used, mids_used, spcs_used, mpcs_used,
            pconv_used, tableht_used, tableh1_used,
            unused_desvars_used,
            remove_nids=True, remove_cids=True,
            remove_pids=True, remove_mids=True,
            remove_spcs=True, remove_mpcs=True,
            unused_remove_desvars=True):
    """actually removes the cards"""
    nids = set(model.nodes.keys())
    pids = set(model.properties.keys())
    pids_mass = set(model.properties_mass.keys())
    cids = set(model.coords.keys())
    mids = set(model.materials.keys())
    spcs = set(model.spcs.keys())  # spcadds?
    mpcs = set(model.mpcs.keys()) # mpcadds?

    nids_to_remove = list(nids - nids_used)
    pids_to_remove = list(pids - pids_used)
    pids_mass_to_remove = list(pids_mass - pids_mass_used)
    mids_to_remove = list(mids - mids_used)
    cids_to_remove = list(cids - cids_used)

    #for subcase in model.subcases:
    #    if 'SPC' in subcase:
    #        value = subcase['SPC']
    #        spcs_used.add(value)
    #    if 'MPC' in subcase:
    #        value = subcase['MPC']
    #        mpcs_used.add(value)
    #    #if 'LOAD' in subcase:
    #        #value = subcase['LOAD']
    #        #loads_used.add(value)

    spcs_to_remove = list(spcs - spcs_used)
    mpcs_to_remove = list(mpcs - mpcs_used)

    if 0 in cids_to_remove:
        cids_to_remove.remove(0)

    if remove_nids and nids_to_remove:
        for nid in nids_to_remove:
            del model.nodes[nid]
        nids_to_remove.sort()
        model.log.debug('removed nodes %s' % nids_to_remove)

    if remove_cids and cids_to_remove:
        for cid in cids_to_remove:
            del model.coords[cid]
        cids_to_remove.sort()
        model.log.debug('removing coords %s' % cids_to_remove)

    if remove_pids and pids_to_remove:
        for pid in pids_mass_to_remove:
            del model.properties_mass[pid]
        pids_mass_to_remove.sort()
        model.log.debug('removing properties_mass %s' % pids_mass_to_remove)

        for pid in pids_to_remove:
            del model.properties[pid]
        pids_to_remove.sort()
        model.log.debug('removing properties %s' % pids_to_remove)

    if remove_mids and mids_to_remove:
        for mid in mids_to_remove:
            del model.materials[mid]
        mids_to_remove.sort()
        model.log.debug('removing materials %s' % mids_to_remove)

    #if remove_spcs and spcs_to_remove:
    #    for spc_id in spcs_to_remove:
    #        del model.spcs[spc_id]
    #    spcs_to_remove.sort()
    #    model.log.debug('removed spcs %s' % spcs_to_remove)
    #if remove_mpcs and mpcs_to_remove:
    #    for mpc_id in mpcs_to_remove:
    #        del model.mpcs[mpc_id]
    #    mpcs_to_remove.sort()
    #    model.log.debug('removed mpcs %s' % mpcs_to_remove)
    _remove_thermal(model, pconv_used, tableht_used, tableh1_used)

def _remove_thermal(model: BDF, pconv_used, tableht_used, tableh1_used) -> None:
    """removes some thermal cards"""
    pconv = {pid for pid, prop in model.convection_properties.items() if prop.type == 'PCONV'}
    tableht = {tid for tid, table in model.tables.items() if table.type == 'TABLEHT'}
    tableh1 = {tid for tid, table in model.tables.items() if table.type == 'TABLEH1'}

    pconv_to_remove = list(pconv - pconv_used)
    tableht_to_remove = list(tableht - tableht_used)
    tableh1_to_remove = list(tableh1 - tableh1_used)

    remove_pconv = True
    if remove_pconv and pconv_to_remove:
        for pid in pconv_to_remove:
            del model.convection_properties[pid]
        pconv_to_remove.sort()
        model.log.debug('removing convection_properties %s' % pconv_to_remove)

    remove_tableh1 = True
    remove_tableht = True
    if remove_tableh1 and tableh1_to_remove:
        for pid in tableh1_to_remove:
            del model.tables[pid]
        tableh1_to_remove.sort()
        model.log.debug('removing TABLEH1 %s' % tableh1_to_remove)
    if remove_tableht and tableht_to_remove:
        for pid in tableht_to_remove:
            del model.tables[pid]
        tableht_to_remove.sort()
        model.log.debug('removing TABLEH1 %s' % tableht_to_remove)
    return
