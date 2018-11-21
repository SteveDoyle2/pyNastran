"""Defines various helper functions for loading a HDF5 BDF file"""
from __future__ import print_function
import numpy as np
from pyNastran.utils.dict_to_h5py import _cast
from pyNastran.bdf.bdf_interface.add_card import CARD_MAP

dict_int_obj_attrs = [
    # are handled explictly
    #'elements',
    #'nodes',
    #'coords',
    #'materials',
    #'properties',
    #'masses',
    #'tables',
    #'methods',
    #'creep_materials', 'csschds',
    #'flutters',
    #'gusts',
    #'trims',
    'MATS1', 'MATS3', 'MATS8', 'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9',

    # TODO: don't work
    #'mpcadds', 'nsmadds',
    #'rigid_elements',
    #'reject_count',

    'aecomps', 'aefacts', 'aelinks', 'aelists', 'aeparams',
    'aestats', 'aesurf', 'aesurfs', 'ao_element_flags', 'bconp', 'bcrparas', 'bcs', 'bctadds',
    'bctparas', 'bctsets', 'blseg', 'bsurf', 'bsurfs', 'cMethods', 'caeros',
    'convection_properties',
    'csuper', 'csupext', 'dareas',
    'dconadds', 'dconstrs', 'ddvals', 'delays', 'dequations', 'desvars', 'divergs', 'dlinks',
    'dload_entries', 'dloads', 'dmigs', 'dmijis', 'dmijs', 'dmiks', 'dmis', 'dphases', 'dresps',
    'dscreen', 'dti', 'dvcrels', 'dvgrids', 'dvmrels', 'dvprels',
    'epoints', 'flfacts',
    'frequencies', 'gridb',
    'hyperelastic_materials',
    'load_combinations', 'loads',
    'nlparms', 'nlpcis',
    'normals', 'nsms', 'nxstrats', 'paeros',
    'pbusht', 'pdampt', 'pelast', 'phbdys', 'plotels', 'points',
    'properties_mass',
    'radcavs', 'radmtx', 'random_tables',
    'ringaxs', 'ringfl',
    'rotors',
    'se_sets', 'se_usets', 'sebndry', 'sebulk', 'seconct', 'seelt',
    'seexcld', 'selabel', 'seload', 'seloc', 'sempln', 'senqset', 'setree', 'sets', 'spcadds',
    'spcoffs',
    'splines', 'spoints',
    'suport1',
    'tables_d', 'tables_m', 'tables_sdamping', 'tempds', 'thermal_materials',
    'tics', 'transfer_functions',
    'tstepnls', 'tsteps', 'usets',
    'view3ds', 'views',
]

def hdf5_load_coords(model, coords_group):
    """loads the coords from an HDF5 file"""
    for card_type in coords_group.keys():
        coords = coords_group[card_type]
        if card_type in ['CORD2R', 'CORD2C', 'CORD2S']:
            if card_type == 'CORD2R':
                func = model.add_cord2r
            elif card_type == 'CORD2C':
                func = model.add_cord2c
            elif card_type == 'CORD2S':
                func = model.add_cord2s

            cids = _cast(coords['cid'])
            rids = _cast(coords['rid'])
            e1s = _cast(coords['e1'])
            e2s = _cast(coords['e2'])
            e3s = _cast(coords['e3'])
            for cid, rid, origin, zaxis, xzplane in zip(
                    cids, rids, e1s, e2s, e3s):
                func(cid, origin, zaxis, xzplane, rid=rid, comment='')
        elif card_type in ['CORD1R', 'CORD1C', 'CORD1S']:
            if card_type == 'CORD1R':
                func = model.add_cord1r
            elif card_type == 'CORD1C':
                func = model.add_cord1c
            elif card_type == 'CORD1S':
                func = model.add_cord1s

            cids = _cast(coords['cid'])
            nodes = _cast(coords['nodes'])
            for cid, (n1, n2, n3) in zip(cid, nodes):
                func(cid, n1, n2, n3, comment='')
        else:
            load_cards_from_keys_values('coords/%s' % card_type, coords)
            #model.add_cord2r(cid, origin, zaxis, xzplane, rid=0, comment='')
            #model.add_cord2c(cid, origin, zaxis, xzplane, rid=0, comment='')
            #model.add_cord2s(cid, origin, zaxis, xzplane, rid=0, comment='')

            #model.add_cord1c(cid, g1, g2, g3, comment='')
            #model.add_cord1s(cid, g1, g2, g3, comment='')

def hdf5_load_tables(unused_model, group):
    for card_type in group.keys():
        sub_group = group[card_type]
        #if card_type == 'TABLES1':
            #pass
        load_cards_from_keys_values('tables/%s' % card_type, sub_group)

def hdf5_load_methods(unused_model, group):
    for card_type in group.keys():
        sub_group = group[card_type]
        #if card_type == 'TABLES1':
            #pass
        load_cards_from_keys_values('methods/%s' % card_type, sub_group)

def hdf5_load_masses(model, group):
    for card_type in group.keys():
        sub_group = group[card_type]
        if card_type == 'CONM2':
            eid = _cast(sub_group['eid'])
            nid = _cast(sub_group['nid'])
            cid = _cast(sub_group['cid'])
            X = _cast(sub_group['X'])
            I = _cast(sub_group['I'])
            mass = _cast(sub_group['mass'])
            for eidi, nidi, cidi, Xi, Ii, massi in zip(eid, nid, cid, X, I, mass):
                model.add_conm2(eidi, nidi, massi, cid=cidi, X=Xi, I=Ii, comment='')

        else:
            #model.add_conm1(eid, nid, mass_matrix, cid=0, comment='')
            load_cards_from_keys_values('masses/%s' % card_type, sub_group)


def hdf5_load_materials(model, group):
    for card_type in group.keys():
        sub_group = group[card_type]
        if card_type == 'MAT1':
            mid = _cast(sub_group['mid'])
            E = _cast(sub_group['E'])
            G = _cast(sub_group['G'])
            nu = _cast(sub_group['nu'])
            rho = _cast(sub_group['rho'])
            a = _cast(sub_group['A'])
            tref = _cast(sub_group['tref'])
            ge = _cast(sub_group['ge'])
            St = _cast(sub_group['St'])
            Sc = _cast(sub_group['Sc'])
            Ss = _cast(sub_group['Ss'])
            mcsid = _cast(sub_group['mcsid'])
            for midi, Ei, Gi, nui, rhoi, ai, trefi, gei, Sti, Sci, Ssi, mcsidi in zip(
                    mid, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid):
                model.add_mat1(midi, Ei, Gi, nui, rho=rhoi, a=ai, tref=trefi,
                               ge=gei, St=Sti, Sc=Sci, Ss=Ssi, mcsid=mcsidi, comment='')

        elif card_type == 'MAT3':
            mid = _cast(sub_group['mid'])
            ex = _cast(sub_group['Ex'])
            eth = _cast(sub_group['Eth'])
            ez = _cast(sub_group['Ez'])

            nuxth = _cast(sub_group['Nuxth'])
            nuzx = _cast(sub_group['Nuzx'])
            nuthz = _cast(sub_group['Nuthz'])
            gxz = _cast(sub_group['Gzx'])

            ax = _cast(sub_group['Ax'])
            ath = _cast(sub_group['Ath'])
            az = _cast(sub_group['Az'])

            rho = _cast(sub_group['rho'])
            tref = _cast(sub_group['tref'])
            ge = _cast(sub_group['ge'])
            for (midi, exi, ethi, ezi, nuxthi, nuzxi, nuthzi,
                 rhoi, gzxi, axi, athi, azi, trefi, gei) in zip(
                     mid, ex, eth, ez, nuxth, nuzx, nuthz, rho, gxz, ax, ath, az, tref, ge):
                model.add_mat3(midi, exi, ethi, ezi, nuxthi, nuthzi, nuzxi, rho=rhoi,
                               gzx=gzxi, ax=axi, ath=athi, az=azi, tref=trefi, ge=gei, comment='')

        #elif card_type == 'MAT8':
            #model.add_mat8(mid, e11, e22, nu12, g12=0.0, g1z=1e8, g2z=1e8, rho=0.,
                           #a1=0., a2=0., tref=0., Xt=0., Xc=None, Yt=0., Yc=None,
                           #S=0., ge=0., F12=0., strn=0., comment='')
        #elif card_type == 'MAT9':
            #model.add_mat9(mid, G11=0., G12=0., G13=0., G14=0., G15=0., G16=0.,
                           #G22=0., G23=0., G24=0., G25=0., G26=0., G33=0., G34=0., G35=0., G36=0.,
                           #G44=0., G45=0., G46=0., G55=0., G56=0., G66=0.,
                           #rho=0., A=None, tref=0., ge=0., comment='')
        elif card_type in ['MAT8', 'MAT9']:
            model.log.warning('skipping materials/%s because its vectorized '
                              'and needs a loader' % card_type)
        else:
            #model.add_mat2(mid, G11, G12, G13, G22, G23, G33, rho=0.,
                           #a1=None, a2=None, a3=None, tref=0., ge=0.,
                           #St=None, Sc=None, Ss=None, mcsid=None, comment='')
            #model.add_mat4(mid, k, cp=0.0, rho=1.0, H=None, mu=None, hgen=1.0,
                           #ref_enthalpy=None, tch=None, tdelta=None, qlat=None, comment='')
            #model.add_mat5(mid, kxx=0., kxy=0., kxz=0., kyy=0., kyz=0., kzz=0.,
                           #cp=0., rho=1., hgen=1., comment='')
            #model.add_mat10(mid, bulk, rho, c, ge=0.0, gamma=None,
                            #table_bulk=None, table_rho=None, table_ge=None,
                            #table_gamma=None, comment='')
            #model.add_mat11(mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23,
                            #rho=0.0, a1=0.0, a2=0.0, a3=0.0, tref=0.0, ge=0.0, comment='')
            load_cards_from_keys_values('materials/%s' % card_type, sub_group)

def hdf5_load_generic(unused_model, group, name):
    for card_type in group.keys():
        sub_group = group[card_type]
        #if card_type == 'TABLES1':
            #pass
        load_cards_from_keys_values('%s/%s' % (name, card_type), sub_group)



def hdf5_load_properties(model, properties_group):
    """loads the properties from an HDF5 file"""
    for card_type in properties_group.keys():
        properties = properties_group[card_type]
        if card_type == 'PSHELL':
            pids = _cast(properties['pid'])
            mids = _cast(properties['mids'])
            z = _cast(properties['z'])
            t = _cast(properties['t'])
            twelveIt3 = _cast(properties['twelveIt3'])
            tst = _cast(properties['tst'])
            nsm = _cast(properties['nsm'])
            for pid, (mid1, mid2, mid3, mid4), (z1, z2), ti, twelveIt3i, tsti, nsmi in zip(
                    pids, mids, z, t, twelveIt3, tst, nsm):
                model.add_pshell(pid, mid1=mid1, t=ti, mid2=mid2, twelveIt3=twelveIt3i,
                                 mid3=mid3, tst=tsti, nsm=nsmi, z1=z1, z2=z2, mid4=mid4,
                                 comment='')
        elif card_type == 'PROD':
            pid = _cast(properties['pid'])
            mid = _cast(properties['mid'])
            A = _cast(properties['A'])
            j = _cast(properties['J'])
            c = _cast(properties['c'])
            nsm = _cast(properties['nsm'])
            for pidi, midi, Ai, ji, ci, nsmi in zip(
                    pid, mid, A, j, c, nsm):
                model.add_prod(pidi, midi, Ai, j=ji, c=ci, nsm=nsm, comment='')

        elif card_type == 'PTUBE':
            pid = _cast(properties['pid'])
            mid = _cast(properties['mid'])
            OD = _cast(properties['OD'])
            t = _cast(properties['t'])
            nsm = _cast(properties['nsm'])
            for pidi, midi, (OD1, OD2), ti, nsmi in zip(
                    pid, mid, OD, t, nsm):
                model.add_ptube(pid, mid, OD1, t=ti, nsm=nsmi, OD2=OD2, comment='')

        elif card_type == 'PBAR':
            pid = _cast(properties['pid'])
            mid = _cast(properties['mid'])
            A = _cast(properties['A'])
            J = _cast(properties['J'])
            I = _cast(properties['I'])

            c = _cast(properties['c'])
            d = _cast(properties['d'])
            e = _cast(properties['e'])
            f = _cast(properties['f'])
            k = _cast(properties['k'])

            nsm = _cast(properties['nsm'])
            for (pidi, midi, Ai, Ji, (i1, i2, i12),
                 (c1, c2), (d1, d2), (e1, e2), (f1, f2), (k1, k2), nsmi) in zip(
                     pid, mid, A, J, I,
                     c, d, e, f, k, nsm):
                model.add_pbar(pidi, midi, A=Ai, i1=i1, i2=i2, i12=i12, j=Ji, nsm=nsmi,
                               c1=c1, c2=c2, d1=d1, d2=d2, e1=e1, e2=e2,
                               f1=f1, f2=f2, k1=k1, k2=k2, comment='')

        else:
            load_cards_from_keys_values('properties/%s' % card_type, properties)
            #model.add_pshear(pid, mid, t, nsm=0., f1=0., f2=0., comment='')
            #model.add_psolid(pid, mid, cordm=0, integ=None, stress=None,
                             #isop=None, fctn='SMECH', comment='')
            #model.add_pvisc(pid, ce, cr, comment='')
            #model.add_pelas(pid, k, ge=0., s=0., comment='')
            #model.add_pdamp(pid, b, comment='')
            #model.add_pcomp(pid, mids, thicknesses, thetas=None, souts=None, nsm=0., sb=0.,
                            #ft=None, tref=0., ge=0., lam=None, z0=None, comment='')
            #model.add_pcompg(pid, global_ply_ids, mids, thicknesses, thetas=None, souts=None,
                             #nsm=0.0, sb=0.0, ft=None, tref=0.0, ge=0.0, lam=None, z0=None,
                             #comment='')


def load_cards_from_keys_values(name, properties):
    try:
        keys = _cast(properties['keys'])
    except KeyError:
        print('name =', name)
        print(properties)
        raise
    values = properties['values']
    for key, keyi in zip(keys, values.keys()):
        value = values[keyi]
        keys_to_read = list(value.keys())
        card_type = _cast(value['type'])
        class_obj = CARD_MAP[card_type]
        if hasattr(class_obj, '_init_from_empty'):
            class_instance = class_obj._init_from_empty()
        else:
            try:
                class_instance = class_obj()
            except TypeError:
                print('error loading %r' % card_type)
                raise

        _properties = []
        if hasattr(class_obj, '_properties'):
            _properties = class_obj._properties
        for key_to_cast in keys_to_read:
            if key_to_cast in _properties:
                continue

            try:
                valuei = _cast(value[key_to_cast])
            except AttributeError:
                valuei = None

            try:
                setattr(class_instance, key_to_cast, valuei)
            except AttributeError:
                print('error loading %r' % card_type)
                print(_properties)
                print(key, key_to_cast, valuei)
                raise



def hdf5_load_elements(model, elements_group):
    """loads the elements from an HDF5 file"""
    for card_type in elements_group.keys():
        elements = elements_group[card_type]
        if card_type == 'CTETRA':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_ctetra(eid, pid, nids, comment='')
        elif card_type == 'CPENTA':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_chexa(eid, pid, nids, comment='')
        elif card_type == 'CHEXA':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_chexa(eid, pid, nids, comment='')

        elif card_type == 'CROD':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_crod(eid, pid, nids, comment='')
        elif card_type == 'CTUBE':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_ctube(eid, pid, nids, comment='')
        elif card_type == 'CONROD':
            eids = _cast(elements['eid'])
            mids = _cast(elements['mid'])
            nodes = _cast(elements['nodes']).tolist()
            A = _cast(elements['A'])
            J = _cast(elements['J'])
            c = _cast(elements['c'])
            nsm = _cast(elements['nsm'])
            for eid, mid, nids, ai, ji, ci, nsmi in zip(eids, mids, nodes, A, J, c, nsm):
                model.add_conrod(eid, mid, nids, A=ai, j=ji, c=ci, nsm=nsmi, comment='')

        elif card_type == 'CBAR':
            # TODO: support OFFT
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            g0 = _cast(elements['g0'])
            x = _cast(elements['x'])
            wa = _cast(elements['wa'])
            wb = _cast(elements['wb'])
            pa = _cast(elements['pa'])
            pb = _cast(elements['pb'])
            for eid, pid, nids, xi, g0i, pai, pbi, wai, wbi in zip(
                    eids, pids, nodes, x, g0, pa, pb, wa, wb):
                if g0i == -1:
                    g0i = None
                if xi[0] == np.nan:
                    xi = [None, None, None]
                model.add_cbar(eid, pid, nids, xi, g0i, offt='GGG',
                               pa=pai, pb=pbi, wa=wai, wb=wbi, comment='')

        elif card_type == 'CBEAM':
            # TODO: support OFFT
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            g0 = _cast(elements['g0'])
            x = _cast(elements['x'])
            sa = _cast(elements['sa'])
            sb = _cast(elements['sb'])
            wa = _cast(elements['wa'])
            wb = _cast(elements['wb'])
            pa = _cast(elements['pa'])
            pb = _cast(elements['pb'])
            for eid, pid, nids, xi, g0i, pai, pbi, wai, wbi, sai, sbi in zip(
                    eids, pids, nodes, x, g0, pa, pb, wa, wb, sa, sb):
                if g0i == -1:
                    g0i = None
                if xi[0] == np.nan:
                    xi = [None, None, None]
                model.add_cbeam(eid, pid, nids, xi, g0i, offt='GGG', bit=None,
                                pa=pai, pb=pbi, wa=wai, wb=wbi, sa=sai, sb=sbi, comment='')

        elif card_type == 'CBUSH':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            g0 = _cast(elements['g0'])
            x = _cast(elements['x'])
            cid = _cast(elements['cid'])
            ocid = _cast(elements['ocid'])
            s = _cast(elements['s'])
            si = _cast(elements['si'])
            for eid, pid, nids, xi, g0i, cidi, s2, ocidi, si2 in zip(
                    eids, pids, nodes, x, g0, cid, s, ocid, si):
                if g0i == -1:
                    g0i = None
                if xi[0] == np.nan:
                    xi = [None, None, None]
                if cidi == -1:
                    cidi = None

                if si2[0] == np.nan:
                    si2 = [None, None, None]
                model.add_cbush(eid, pid, nids, x, g0, cid=cidi, s=s2, ocid=ocidi, si=si2,
                                comment='')

        elif card_type == 'CTRIA3':
            # TODO: doesn't support tflag
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            thetas = _cast(elements['theta'])
            mcids = _cast(elements['mcid'])
            zoffsets = _cast(elements['zoffset'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, mcid, theta, zoffset in zip(
                    eids, pids, nodes, mcids, thetas, zoffsets):
                if mcid == -1:
                    theta_mcid = theta
                else:
                    theta_mcid = mcid
                model.add_ctria3(eid, pid, nids, zoffset=zoffset, theta_mcid=theta_mcid,
                                 tflag=0, T1=None, T2=None, T3=None, comment='')
        elif card_type == 'CQUAD4':
            # TODO: doesn't support tflag
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            thetas = _cast(elements['theta'])
            mcids = _cast(elements['mcid'])
            zoffsets = _cast(elements['zoffset'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, mcid, theta, zoffset in zip(
                    eids, pids, nodes, mcids, thetas, zoffsets):
                if mcid == -1:
                    theta_mcid = theta
                else:
                    theta_mcid = mcid
                model.add_cquad4(eid, pid, nids, zoffset=zoffset, theta_mcid=theta_mcid,
                                 tflag=0, T1=None, T2=None, T3=None, T4=None, comment='')

        elif card_type == 'CSHEAR':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, mcid, theta, zoffset in zip(eids, pids, nodes):
                model.add_cshear(eid, pid, nids, comment='')
        else:
            load_cards_from_keys_values('elements/%s' % card_type, elements)
            #model.add_celas1(eid, pid, nids, c1=0, c2=0, comment='')
            #model.add_celas2(eid, k, nids, c1=0, c2=0, ge=0., s=0., comment='')
            #model.add_celas3(eid, pid, nids, comment='')
            #model.add_celas4(eid, k, nids, comment='')
            #model.add_cdamp1(eid, pid, nids, c1=0, c2=0, comment='')
            #model.add_cdamp2(eid, b, nids, c1=0, c2=0, comment='')
            #model.add_cdamp3(eid, pid, nids, comment='')
            #model.add_cdamp4(eid, b, nids, comment='')
            #model.add_cvisc(eid, pid, nids, comment='')
            #model.add_cbush1d(eid, pid, nids, cid=None, comment='')
            #model.add_cbush2d(eid, pid, nids, cid=0, plane='XY', sptid=None, comment='')
            #model.add_cdamp5(eid, pid, nids, comment='')
            #model.add_cfast(eid, pid, Type, ida, idb, gs=None, ga=None, gb=None,
                            #xs=None, ys=None, zs=None, comment='')
            #model.add_cgap(eid, pid, nids, x, g0, cid=None, comment='')
            #model.add_cmass1(eid, pid, nids, c1=0, c2=0, comment='')
            #model.add_cmass2(eid, mass, nids, c1, c2, comment='')
            #model.add_cmass3(eid, pid, nids, comment='')
            #model.add_cmass4(eid, mass, nids, comment='')
            #model.add_cquad8(eid, pid, nids, theta_mcid=0., zoffset=0.,
                             #tflag=0, T1=None, T2=None, T3=None, T4=None, comment='')
            #model.add_cquad(eid, pid, nids, theta_mcid=0., comment='')
            #model.add_cquadr(eid, pid, nids, theta_mcid=0.0, zoffset=0.,
                             #tflag=0, T1=None, T2=None, T3=None, T4=None, comment='')
            #model.add_cquadx(eid, pid, nids, theta_mcid=0., comment='')
            #model.add_cquadx4(eid, pid, nids, theta=0., comment='')
            #model.add_cquadx8(eid, pid, nids, theta=0., comment='')
            #model.add_ctria6(eid, pid, nids, theta_mcid=0., zoffset=0.,
                             #tflag=0, T1=None, T2=None, T3=None, comment='')
            #model.add_ctriar(eid, pid, nids, theta_mcid=0.0, zoffset=0.0,
                             #tflag=0, T1=None, T2=None, T3=None, comment='')
            #model.add_ctriax(eid, pid, nids, theta_mcid=0., comment='')
            #model.add_ctriax6(eid, mid, nids, theta=0., comment='')
            print(card_type)
