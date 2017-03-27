"""
defines some methods for cleaning up a model
 - remove_unassociated_nodes(...)
 - remove_unassociated_properties(...)
 - remove_unused_materials(...)
"""
from __future__ import print_function
from six import iteritems, itervalues, integer_types

from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber

def remove_unused(bdf_filename, remove_nids=True, remove_cids=True,
                  remove_pids=True, remove_mids=True):
    """
    Takes an uncross-referenced bdf and removes unused data

    removes unused:
     - nodes
     - properties
     - materials
     - coords
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

    nids_used = set([])
    cids_used = set([])
    pids_used = set([])
    pids_mass_used = set([])
    mids_used = set([])
    mids_thermal_used = set([])

    #card_types = list(model.card_count.keys())
    #card_map = model.get_card_ids_by_card_types(
        #card_types=card_types,
        #reset_type_to_slot_map=False,
        #stop_on_missing_card=True)

    #for nid, node in iteritems(model.nodes):
        #cids_used.update([node.Cp(), node.Cd()])

    skip_cards = [
        'ENDDATA', 'PARAM', 'EIGR', 'EIGRL', 'EIGB', 'EIGP', 'EIGC',
        'SPOINT', 'EPOINT', 'DESVAR',
        'SET1', 'FREQ', 'FREQ1', 'FREQ2',
        'TSTEP', 'TSTEPNL', 'NLPCI',
        #'LOAD', 'LSEQ', 'DLOAD', 'LOADCYN',
        'NLPARM', 'ROTORG', 'ROTORD',
        'DAREA', 'DEQATN',
        'DMIG', 'DMI', 'DMIJ', 'DMIK', 'DMIJI',
        'POINT', 'EPOINT',
        'DELAY', 'DPHASE',
        'CBARAO',


        # properties
        'PELAS', 'PDAMP', 'PBUSH',
        'PELAST', 'PDAMPT', 'PBUSHT',
        'PGAP', 'PBUSH1D', 'PFAST', 'PVISC', 'PMASS',

        'FLFACT', 'FLUTTER', 'DLINK', 'DDVAL', 'DIVERG', 'GUST',
        'AELINK', 'AELIST', 'TRIM', 'PAERO1', 'AEFACT', 'AESTAT',

        'BCTPARA', 'BCRPARA', 'BSURF', 'BSURFS', 'BCTADD',
        'BCTSET',

        # not checked------------------------------------------
        'PHBDY', 'CHBDYG', 'CHBDYP', 'CHBDYE', 'RADBC', 'CONV',
        'QVOL', 'PCONV', 'PCONVM',
        #'PBCOMP', 'PDAMP5',
        'AECOMP', 'CAERO2', 'CAERO3', 'CAERO4', 'PAERO3', 'PAERO4', #'CFAST',
        'DCONADD',
        'GMCORD',
        'MONPNT1', 'MONPNT2', 'MONPNT3',
    ]
    set_types = [
        'SET1', 'SET3',
        'ASET', 'ASET1', 'BSET', 'BSET1', 'CSET', 'CSET1',
        'QSET', 'SSET1', 'USET', 'USET1',
        'SESET',
    ]
    load_types = [
        'GRAV', 'RANDPS', 'FORCE', 'FORCE1', 'FORCE2',
        'MOMENT', 'MOMENT1', 'MOMENT2',
        'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4', 'SPCD',
        'GMLOAD', 'RFORCE', 'RFORCE1',
        'TEMP', 'QBDY1', 'QBDY2', 'QBDY3', 'QHBDY',
        'ACCEL', 'PLOADX1', 'SLOAD', 'ACCEL1', 'LOADCYN', 'LOAD',
        'LSEQ', 'DLOAD',
    ]

    # could remove some if we look at the rid_trace
    #for cid, coord in iteritems(model.coords):
        #if coord.type in ['CORD1R', 'CORD1C', 'CORD1S']:
            #nids_used.update(node_ids)
        #elif coord.type in ['CORD1R', 'CORD1C', 'CORD1S']:
            #cids_used.update(coord.Rid())
        #else:
            #raise NotImplementedError(coord)

    for card_type, ids in iteritems(model._type_to_id_map):
    #for card_type, ids in iteritems(card_map):
        if card_type in ['CORD1R', 'CORD1C', 'CORD1S']:
            #print(ids)
            for cid in ids:
                coord = model.coords[cid]
                nids_used.update(coord.node_ids)
        elif card_type in ['CORD2R', 'CORD2C', 'CORD2S']:
            #print(ids)
            for cid in ids:
                coord = model.coords[cid]
                cids_used.add(coord.Rid())

        elif card_type in ['MAT1', 'MAT2', 'MAT3', 'MAT4', 'MAT5',
                           'MAT8', 'MAT9', 'MAT10', 'MAT11']:
            # todo: MATS1, MATT1, etc.
            pass
        elif card_type in ['MATS1', 'MATT1', 'MATT2', 'MATT4', 'MATT5',
                           'MATHE', 'MATHP', 'CREEP']:
            mids_used.update(ids)

        elif card_type in ['CTETRA', 'CPENTA', 'CPYRAM', 'CHEXA']:
            for eid in ids:
                elem = model.elements[eid]
                nids_used.update(elem.node_ids)
                pids_used.add(elem.Pid())

        elif card_type in ['CONM1', 'CONM2']:
            for eid in ids:
                elem = model.masses[eid]
                nids_used.add(elem.Nid())
                cids_used.add(elem.Cid())
                #print(elem.object_attributes())
                #print(elem.object_methods())
                #aaa
        elif card_type in ['CMASS1', 'CMASS3']:
            for eid in ids:
                elem = model.masses[eid]
                pids_mass_used.add(elem.Pid())
                nids_used.update(elem.node_ids)
        elif card_type in ['CMASS2', 'CMASS4']:
            for eid in ids:
                elem = model.masses[eid]
                nids_used.update(elem.node_ids)

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
        elif card_type in ['CELAS4', 'CDAMP4']:
            for eid in ids:
                elem = model.elements[eid]
                nids_used.update(elem.node_ids)

        elif card_type in ['CTRIA3', 'CQUAD4', 'CTRIA6', 'CTRIAR', 'CQUAD8', 'CQUADR',
                           'CTRIAX', 'CQUADX', 'CQUAD']:
            for eid in ids:
                elem = model.elements[eid]
                nids_used.update(elem.node_ids)
                pids_used.add(elem.Pid())
                if isinstance(elem.theta_mcid, integer_types):
                    cids_used.add(elem.theta_mcid)
        elif card_type in ['CTRIAX6']:
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
        elif card_type == 'PLPLANE':
            for pid in ids:
                prop = model.properties[pid]
                cids_used.add(prop.cid)
                mids_used.add(prop.Mid())
        elif card_type == 'PPLANE':
            for pid in ids:
                prop = model.properties[pid]
                mids_used.add(prop.Mid())

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
        elif card_type in ['PLOTEL']:
            for eid in ids:
                elem = model.plotels[eid]
                nids_used.update(elem.node_ids)
        elif card_type in ['PSOLID', 'PLSOLID']:
            for pid in ids:
                prop = model.properties[pid]
                mids_used.add(prop.Mid())
        elif card_type in ['PDAMP5']:
            for pid in ids:
                prop = model.properties[pid]
                mids_thermal_used.add(prop.Mid())
        elif card_type in ['PBAR', 'PBARL', 'PROD', 'PTUBE', 'PBEAM', 'PBEAML', 'PSHEAR',
                           'PRAC2D', 'PRAC3D', 'PBEND']:
            for pid in ids:
                prop = model.properties[pid]
                mids_used.add(prop.Mid())
        elif card_type in ['PSHELL']:
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
        elif card_type in ['PCOMPS']:
            for pid in ids:
                prop = model.properties[pid]
                mids = prop.Mids()
                mids_used.update(mids)
                cids_used.update(prop.cordm)

        elif card_type in ['RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3', 'RROD', 'RSPLINE']:
            for eid in ids:
                elem = model.rigid_elements[eid]
                #print(elem.object_attributes())
                #print(elem.object_methods())
                nids_used.update(elem.independent_nodes)
                nids_used.update(elem.dependent_nodes)

        elif card_type in ['TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2', 'ACSRCE']:
            pass
        elif card_type in load_types:
            for loads in itervalues(model.loads):
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
                    elif load.type in ['QVOL']:
                        # eids
                        pass
                    else:
                        raise NotImplementedError(load)

        elif card_type == 'TEMPD':
            pass
            #for temp_id in ids:
                #tempd = self.tempds[temp_id]

        elif card_type in ['MPCADD', 'MPC']:
            for mpcs in itervalues(model.mpcs):
                for mpc in mpcs:
                    if mpc.type in ['MPCADD']:
                        pass
                    elif mpc.type in ['MPC']:
                        nids_used.update(mpc.node_ids)
                    else:
                        raise NotImplementedError(mpc)

        elif card_type in ['SPCADD', 'SPC1', 'SPC', 'GMSPC', 'SPCAX']:
            for spcs in itervalues(model.spcs):
                for spc in spcs:
                    if spc.type in ['SPCADD', 'GMSPC', 'SPCAX']:
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
            for suport1 in itervalues(model.suport1):
                nids_used.update(suport1.node_ids)
        elif card_type == 'GRID':
            for nid, node in iteritems(model.nodes):
                cids_used.update([node.Cp(), node.Cd()])

        elif card_type in ['CBAR', 'CBEAM', 'CBEND']:
            for eid in ids:
                elem = model.elements[eid]
                nids_used.update(elem.node_ids)
                pids_used.add(elem.Pid())
                if elem.g0 is not None:
                    assert isinstance(elem.g0, integer_types), elem.g0
                    nids_used.add(elem.g0)
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
                    assert isinstance(elem.G0(), integer_types), elem.G0()
                    nids_used.add(elem.G0())
        elif card_type in ['CBUSH1D', 'CBUSH2D']:
            for eid in ids:
                elem = model.elements[eid]
                nids_used.update(elem.node_ids)
                pids_used.add(elem.Pid())
                cids_used.add(elem.Cid())


        elif card_type in ['PBUSH']:
            pass
            #for pid in ids:
                #prop = model.properties[pid]
                #raise RuntimeError(prop)
        elif card_type == 'PBUSHT':
            # tables
            pass
        elif card_type in ['CBUSH']:
            for eid in ids:
                elem = model.elements[eid]
                nids_used.update(elem.node_ids)
                pids_used.add(elem.Pid())
                if elem.g0 is not None:
                    assert isinstance(elem.g0, integer_types), elem.g0
                    nids_used.add(elem.g0)
                # TODO: cid

        elif card_type == 'AESURF':
            #CID1  | ALID1 | CID2   | ALID2
            for aesurf in itervalues(model.aesurf):
                cids_used.add(aesurf.Cid1())
                cid2 = aesurf.Cid2()
                if cid2 is not None:
                    cids_used.add(cid2)
        elif card_type in ['SPLINE1', 'SPLINE2', 'SPLINE3', 'SPLINE4', 'SPLINE5']:
            pass
            #for spline_id in ids:
                #spline = model.splines[spline_id]
        elif card_type in ['CAERO1']:
            for eid in ids:
                caero = model.caeros[eid]
                # PID, LSPAN, LCHORD
                cids_used.add(caero.Cp())

        elif card_type in skip_cards:
            pass
        elif card_type in set_types:
            pass
        elif card_type in ['DCONSTR']:
            pass
        elif card_type == 'DRESP1':
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
        elif card_type == 'DRESP2':
            pass
            #for dresp_id in ids:
                #dresp = model.dresps[dresp_id]
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
                if dvprel.Type in ['PSHELL', 'PCOMP', 'PBAR', 'PBARL', 'PBEAM',
                                   'PROD', 'PELAS', 'PBUSH', 'PDAMP', 'PTUBE',
                                   'PSHEAR', 'PDAMP', 'PMASS', 'PBEAML', 'PCOMPG',
                                   'PVISC', 'PBUSHT', 'PELAST', 'PBUSH1D', 'PGAP']:
                    pids_used.add(dvprel.Pid())
                elif dvprel.Type in ['DISP']:
                    raise NotImplementedError(str(dvprel) + 'dvprel.Type=DISP')
                else:
                    raise NotImplementedError(dvprel)

        elif card_type in ['DVCREL1', 'DVCREL2']:
            for dvcrel_id in ids:
                dvcrel = model.dvcrels[dvcrel_id]
                if dvcrel.Type in ['CMASS2', 'CMASS4', 'CONM1', 'CONM2',
                                   'CELAS2', 'CELAS4',
                                   'CDAMP2', 'CQUAD4', 'CGAP', 'CBAR']:
                    pass
                    #pids_used.add(dvcrel.Eid())
                else:
                    raise NotImplementedError(str(dvcrel) + 'Type=%r' % dvcrel.Type)

        elif card_type in ['DVMREL1', 'DVMREL2']:
            for dvmrel_id in ids:
                dvmrel = model.dvmrels[dvmrel_id]
                if dvmrel.Type in ['MAT1', 'MAT2', 'MAT8', 'MAT9', 'MAT11']:
                    mids_used.add(dvmrel.Mid())
                else:
                    raise NotImplementedError(str(dvmrel) + 'Type=%r' % dvmrel.Type)
        elif card_type == 'DVGRID':
            for dvgrid_id in ids:
                dvgrids = model.dvgrids[dvgrid_id]
                for dvgrid in dvgrids:
                    nids_used.add(dvgrid.nid)
                    cids_used.add(dvgrid.cid)
        elif card_type == 'TF':
            for tf_id in ids:
                tfs = model.transfer_functions[tf_id]
                for transfer_function in tfs:
                    nids_used.update(transfer_function.nids)
        else:
            raise NotImplementedError(card_type)


    #for pid, prop in iteritems(model.properties):
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

    nids = set(model.nodes.keys())
    pids = set(model.properties.keys())
    pids_mass = set(model.properties_mass.keys())
    cids = set(model.coords.keys())
    mids = set(model.materials.keys())
    nids_to_remove = list(nids - nids_used)
    pids_to_remove = list(pids - pids_used)
    pids_mass_to_remove = list(pids_mass - pids_mass_used)
    mids_to_remove = list(mids - mids_used)
    cids_to_remove = list(cids - cids_used)

    if remove_nids:
        for nid in nids_to_remove:
            del model.nodes[nid]
        model.log.debug('removed nodes %s' % nids_to_remove)

    if remove_cids:
        for cid in cids_to_remove:
            del model.coords[cid]
        model.log.debug('removing coords %s' % cids_to_remove)

    if remove_pids:
        for pid in pids_mass_to_remove:
            del model.properties_mass[pid]
        model.log.debug('removing properties_mass %s' % pids_mass_to_remove)

        for pid in pids_to_remove:
            del model.properties[pid]
        model.log.debug('removing properties %s' % pids_to_remove)

    if remove_mids:
        for mid in mids_to_remove:
            del model.materials[mid]
        model.log.debug('removing materials %s' % mids_to_remove)
    return model
