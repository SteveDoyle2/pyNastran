"""
defines:
 - convert(model, units_to, units=None)
"""
from __future__ import print_function
from six import string_types
import numpy as np
from pyNastran.bdf.cards.base_card import break_word_by_trailing_parentheses_integer_ab

def convert(model, units_to, units=None):
    """
    Converts a model from a set of defined units

    Parameters
    ----------
    xref : bool
       cross references the model (default=True)
    units_to : List[str]
        [length, mass, time]
    units : list
        overwrites model.units

    """
    # units_start = 'in'
    # units_end = 'mm'
    # xyz_scale = in_to_mm
    # xyz_scale = 25.4
    if units is None:
        units = model.units
    xyz_scale, mass_scale, time_scale, weight_scale, gravity_scale = get_scale_factors(
        units, units_to)

    model.log.debug('xyz_scale = %s' % xyz_scale)
    model.log.debug('mass_scale = %s' % mass_scale)
    model.log.debug('time_scale = %s' % time_scale)
    model.log.debug('weight_scale = %s' % weight_scale)
    model.log.debug('gravity_scale = %s' % gravity_scale)
    #_set_wtmass(model, gravity_scale)

    _convert_nodes(model, xyz_scale)
    _convert_coordinates(model, xyz_scale)

    _convert_elements(model, xyz_scale, mass_scale, weight_scale)
    _convert_properties(model, xyz_scale, mass_scale, weight_scale)
    #_convert_masses(model)
    _convert_materials(model, xyz_scale, mass_scale, weight_scale)

    _convert_aero(model, xyz_scale, time_scale, weight_scale)
    _convert_constraints(model, xyz_scale)
    _convert_loads(model, xyz_scale, weight_scale)
    #_convert_sets(model)
    _convert_optimization(model, xyz_scale, mass_scale, weight_scale)


def _set_wtmass(model, gravity_scale):
    """
    set the PARAM,WTMASS

    ft-lbm-s : 1. / 32.2
    in-lbm-s : 1. / (32.2*12)
    m-kg-s   : 1.
    mm-Mg-s  : 1.

    """
    if 'WTMASS' in model.params:
        param = model.params['WTMASS']
        value0 = param.values[0]
        param.values = [value0 / gravity_scale]
    else:
        card = ['PARAM', 'WTMASS', 1. / gravity_scale]
        model.add_card(card, 'PARAM')

def _convert_nodes(model, xyz_scale):
    """converts the nodes"""
    for node in model.nodes.values():
        if node.cp == 0:
            node.xyz *= xyz_scale
        elif node.cp_ref.type in ['CORD1R', 'CORD2R']:
            node.xyz *= xyz_scale
        else:
            # only scale R
            node.xyz[0] *= xyz_scale

def _convert_coordinates(model, xyz_scale):
    """converts the coordinate systems"""
    for cid, coord in model.coords.items():
        if cid == 0:
            continue
        if coord.rid == 0:
            if np.abs(coord.e1).sum() == 0.:
                assert np.abs(coord.origin).sum() == 0., coord.origin
            else:
                coord.origin *= xyz_scale
                coord.e1 *= xyz_scale
                coord.e2 *= xyz_scale
                coord.e3 *= xyz_scale
        elif coord.rid.type in ['CORD1R', 'CORD2R']:
            coord.origin *= xyz_scale
            coord.e1 *= xyz_scale
            coord.e2 *= xyz_scale
            coord.e3 *= xyz_scale
        elif coord.rid.type in ['CORD1C', 'CORD1S', 'CORD2C', 'CORD2S']:
            coord.origin[0] *= xyz_scale
            coord.e1 *= xyz_scale
            raise NotImplementedError(coord)
        else:
            raise NotImplementedError(coord)
            #coord.e1 *= xyz_scale
        #elif coord.type in ['CORD1R', 'CORD2R']:
            #coord.origin *= xyz_scale
            #coord.e1 *= xyz_scale
            #coord.e2 *= xyz_scale
            #coord.e3 *= xyz_scale
        #else:
            #raise NotImplementedError(coord)

def _convert_elements(model, xyz_scale, mass_scale, weight_scale):
    """converts the elements"""
    time_scale = 1.
    force_scale = weight_scale
    area_scale = xyz_scale ** 2
    velocity_scale = xyz_scale / time_scale
    moi_scale = xyz_scale ** 4
    mass_moi_scale = mass_scale * xyz_scale ** 2
    #nsm_scale = mass_scale
    nsm_bar_scale = mass_scale / xyz_scale
    stiffness_scale = force_scale / xyz_scale
    damping_scale = force_scale / velocity_scale

    # these don't have any properties
    skip_elements = [
        'CCONEAX',
        'CELAS1', 'CELAS3', 'CDAMP3', 'CDAMP5', 'CVISC', 'CBUSH1D',
        'CROD', 'CTUBE',
        'CBUSH',
        'CQUAD', 'CSHEAR', 'CTRIAX', 'CTRIAX6',
        'CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM',

        # TODO: NX-verify
        'CTRAX3', 'CTRAX6',
        'CPLSTN3', 'CPLSTN6', 'CPLSTN4', 'CPLSTN8',
        'CQUADX4', 'CQUADX8',
    ]
    #skip_elements = ['CELAS2', 'CELAS4']

    tri_shells = ['CTRIA3', 'CTRIA6', 'CTRIAR']
    quad_shells = [
        'CQUAD4', 'CQUAD8', 'CQUADR']
    spring_elements = ['CELAS2', 'CELAS4']
    damper_elements = ['CDAMP2', 'CDAMP4']

    model.log.debug('--Element Scales--')
    model.log.debug('nsm_bar_scale = %g' % nsm_bar_scale)
    model.log.debug('moi_scale = %g' % moi_scale)
    model.log.debug('area_scale = %g' % area_scale)
    model.log.debug('stiffness_scale = %g\n' % stiffness_scale)
    model.log.debug('damping_scale = %g\n' % damping_scale)
    if len(model.masses):
        model.log.debug('mass_moi_scale = %g' % mass_moi_scale)

    for elem in model.elements.values():
        elem_type = elem.type
        if elem_type in skip_elements:
            continue

        if elem_type in spring_elements:
            elem.k *= stiffness_scale
        elif elem_type in damper_elements:
            elem.b *= damping_scale
        elif elem_type in tri_shells:
            # thickness
            elem.zoffset *= xyz_scale
            if elem.tflag == 0:
                if elem.T1 is not None:
                    elem.T1 *= xyz_scale
                    elem.T2 *= xyz_scale
                    elem.T3 *= xyz_scale
            # nsm
            #elem.nsm *= nsm_scale

        elif elem_type in quad_shells:
            # thickness
            # tflag=blank/0 - Ti value: Ti is specified thickness
            #                 Ti blank: Ti is PSHELL value
            # tflag=0
            # tflag=1 - thicknesses are relative (Ti default=1.0)
            #
            elem.zoffset *= xyz_scale
            if elem.tflag == 0:
                if elem.T1 is not None:
                    elem.T1 *= xyz_scale
                    elem.T2 *= xyz_scale
                    elem.T3 *= xyz_scale
                    elem.T4 *= xyz_scale
            # nsm
            #elem.nsm *= nsm_scale
        elif elem_type in ['CQUADX']:
            pass

        elif elem_type == 'CONROD':
            elem.A *= area_scale # area
            elem.nsm *= nsm_bar_scale
        elif elem_type == 'CGAP':
            if elem.x is not None:
                # vector
                elem.x = [x*xyz_scale for x in elem.x]

        elif elem_type == 'CBAR':
            if elem.x is not None:
                # vector
                #elem.x = [x*xyz_scale for x in elem.x]
                elem.wa *= xyz_scale
                elem.wb *= xyz_scale
        elif elem_type == 'CBEAM':
            if elem.x is not None:
                # vector
                #elem.x = [x*xyz_scale for x in elem.x]
                elem.wa *= xyz_scale
                elem.wb *= xyz_scale
        else:
            raise NotImplementedError('type=%r; elem:\n%s' % (elem.type, elem))

    for elem in model.masses.values():
        elem_type = elem.type
        if elem_type == 'CONM2':
            elem.mass *= mass_scale
            elem.X *= xyz_scale
            # I = m * r^2
            elem.I = [moi * mass_moi_scale for moi in elem.I]
        elif elem.type == 'CMASS1':
            pass
        elif elem.type == 'CMASS4':
            elem.mass *= mass_scale
        else:
            raise NotImplementedError(elem)

def _convert_properties(model, xyz_scale, mass_scale, weight_scale):
    """converts the properties"""
    time_scale = 1.
    force_scale = weight_scale
    moment_scale = force_scale * xyz_scale
    area_scale = xyz_scale ** 2
    moi_scale = xyz_scale ** 4
    velocity_scale = xyz_scale / time_scale

    # there are multiple nsm scales (CONM2, bar, plate)
    nsm_bar_scale = mass_scale / xyz_scale
    nsm_plate_scale = mass_scale / xyz_scale ** 2
    stiffness_scale = force_scale / xyz_scale
    damping_scale = force_scale / velocity_scale
    stress_scale = force_scale / xyz_scale ** 2

    model.log.debug('--Property Scales--')
    model.log.debug('nsm_bar_scale = %g' % nsm_bar_scale)
    model.log.debug('nsm_plate_scale = %g' % nsm_plate_scale)
    model.log.debug('stiffness_scale = %g' % stiffness_scale)
    model.log.debug('damping_scale = %g' % damping_scale)
    model.log.debug('stress_scale = %g\n' % stress_scale)

    skip_properties = (
        'PSOLID', 'PLSOLID', 'PLPLANE', 'PIHEX',

        # TODO: NX-verify
        'PPLANE',
    )
    for prop in model.properties.values():
        prop_type = prop.type
        if prop_type in skip_properties:
            continue
        elif prop_type == 'PELAS':
            prop.k *= stiffness_scale
        elif prop_type in ['PDAMP', 'PDAMP5']:
            prop.b *= damping_scale # force_scale / velocity_scale
        elif prop_type == 'PVISC':
            prop.ce *= force_scale / velocity_scale
            prop.cr *= moment_scale / velocity_scale

        elif prop_type == 'PROD':
            prop.A *= area_scale
            prop.j *= moi_scale
            #prop.c ???

        elif prop_type == 'PBAR':
            prop.A *= area_scale
            prop.i1 *= moi_scale
            prop.i2 *= moi_scale
            prop.i12 *= moi_scale
            prop.j *= moi_scale
            prop.nsm *= nsm_bar_scale
            prop.c1 *= xyz_scale
            prop.c2 *= xyz_scale
            prop.d1 *= xyz_scale
            prop.d2 *= xyz_scale
            prop.e1 *= xyz_scale
            prop.e2 *= xyz_scale
            prop.f1 *= xyz_scale
            prop.f2 *= xyz_scale

        elif prop_type == 'PBARL':
            prop.dim = [d * xyz_scale for d in prop.dim]
            prop.nsm *= nsm_bar_scale

        elif prop_type == 'PBEAM':
            prop.A *= area_scale
            prop.i1 *= moi_scale
            prop.i2 *= moi_scale
            prop.i12 *= moi_scale
            prop.j *= moi_scale
            prop.nsm *= nsm_bar_scale
            prop.c1 *= xyz_scale
            prop.c2 *= xyz_scale
            prop.d1 *= xyz_scale
            prop.d2 *= xyz_scale
            prop.e1 *= xyz_scale
            prop.e2 *= xyz_scale
            prop.f1 *= xyz_scale
            prop.f2 *= xyz_scale

            prop.m1a *= xyz_scale
            prop.m2a *= xyz_scale
            prop.m1b *= xyz_scale
            prop.m2b *= xyz_scale
            prop.n1a *= xyz_scale
            prop.n2a *= xyz_scale
            prop.n1b *= xyz_scale
            prop.n2b *= xyz_scale

        elif prop_type == 'PBEAML':
            prop.dim *= xyz_scale
            prop.nsm *= nsm_bar_scale

        elif prop_type == 'PSHELL':
            prop.t *= xyz_scale
            prop.nsm *= nsm_plate_scale
            prop.z1 *= xyz_scale
            prop.z2 *= xyz_scale
            prop.twelveIt3 /= xyz_scale ** 3
        elif prop_type == 'PSHEAR':
            prop.t *= xyz_scale
            prop.nsm *= nsm_plate_scale

        elif prop_type in ['PCOMP', 'PCOMPG']:
            prop.thicknesses = [t * xyz_scale for t in prop.thicknesses]
            prop.nsm *= nsm_plate_scale
            prop.z0 *= xyz_scale
            prop.sb *= stress_scale

        elif prop_type == 'PELAS':
            prop.k *= stiffness_scale
        elif prop_type == 'PTUBE':
            prop.OD1 *= xyz_scale
            prop.OD2 *= xyz_scale
            prop.t *= xyz_scale
            prop.nsm *= nsm_bar_scale

        elif prop_type == 'PBUSH':
            # can be length=0
            #assert len(prop.Ki) == 6, prop.Ki
            #assert len(prop.Bi) == 6, prop.Bi
            for var in prop.vars:
                # TODO: I think this needs to consider rotation
                if var == 'K':
                    prop.Ki = [ki*stiffness_scale if ki is not None else None
                               for ki in prop.Ki]
                elif var == 'B':
                    prop.Bi = [bi*velocity_scale if bi is not None else None
                               for bi in prop.Bi]
                elif var == 'GE':
                    pass
                else:
                    raise NotImplementedError(prop)

            #prop.rcv
            if prop.mass is not None:
                prop.mass *= mass_scale
            #rcv : List[float]; default=None -> (None, None, None, None)
                #[sa, st, ea, et] = rcv
                #length(mass_fields) = 4
            #mass : float; default=None
                #lumped mass of the CBUSH
                #This is an MSC only parameter.
        elif prop.type in ['PBUSH1D', 'PBUSH2D']:
            model.log.warning('skipping:\n%s' % str(prop))

        #elif prop.type == 'PCOMPS':
            #pass
        elif prop.type == 'PCONEAX':
            #T1 Membrane thickness. (Real > 0.0 if MID1 = 0)
            #T2 Transverse shear thickness. (Real > 0.0 if MID3 = 0)
            #I Moment of inertia per unit width. (Real)
            #NSM Nonstructural mass per unit area. (Real)
            #Z1, Z2 Fiber distances from the middle surface for stress recovery. (Real)
            #prop.mid1 = mid1
            if prop.t1 is not None:
                prop.t1 *= xyz_scale
            if prop.i is not None:
                prop.i *= moi_scale / xyz_scale
            if prop.t2 is not None:
                prop.t2 *= xyz_scale
            prop.nsm *= nsm_plate_scale
            prop.z1 *= xyz_scale
            prop.z2 *= xyz_scale

        #elif prop.type == 'PBCOMP':
            #pass
        #elif prop.type == 'PPLANE':
            #pass
        #elif prop.type == 'PRAC2D':
            #pass
        #elif prop.type == 'PRAC3D':
            #pass
        #elif prop.type == 'PFAST':
            #pass
        #elif prop.type == 'PDAMP':
            #pass
        #elif prop.type == 'PTUBE':
            #pass
        #elif prop.type == 'PELAST':
            #pass
        #elif prop.type == 'PBUSHT':
            #pass
        #elif prop.type == 'PBUSH1D':
            #pass
        #elif prop.type == 'PDAMPT':
            #pass
        #elif prop.type == 'PDAMP5':
            #pass
        elif prop.type == 'PGAP':
            #: initial gap opening
            prop.u0 *= xyz_scale
            #: preload
            prop.f0 *= force_scale

            #: axial stiffness of closed/open gap
            prop.ka *= stiffness_scale
            prop.kb *= stiffness_scale

            #: transverse stiffness of closed gap
            prop.kt *= stiffness_scale
        else:
            raise NotImplementedError(prop_type)

def _convert_materials(model, xyz_scale, mass_scale, weight_scale):
    """converts the materials"""
    density_scale = mass_scale / xyz_scale ** 3
    stress_scale = weight_scale / xyz_scale ** 2
    temp_scale = 1.
    a_scale = 1. / temp_scale # thermal expansion

    model.log.debug('--Material Scales--')
    model.log.debug('density_scale = %g' % density_scale)
    model.log.debug('stress_scale = %g\n' % stress_scale)

    for mat in model.materials.values():
        mat_type = mat.type
        if mat_type == 'MAT1':
            mat.e *= stress_scale
            mat.g *= stress_scale
            mat.a *= a_scale
            mat.tref *= temp_scale
            mat.rho *= density_scale
            mat.St *= stress_scale
            mat.Sc *= stress_scale
            mat.Ss *= stress_scale

        elif mat_type == 'MAT2':
            mat.G11 *= stress_scale
            mat.G12 *= stress_scale
            mat.G13 *= stress_scale
            mat.G22 *= stress_scale
            mat.G23 *= stress_scale
            mat.G33 *= stress_scale
            mat.rho *= density_scale
            # 1/dTemp
            if mat.a1 is not None:
                mat.a1 *= a_scale
            if mat.a2 is not None:
                mat.a2 *= a_scale
            if mat.a3 is not None:
                mat.a3 *= a_scale
            mat.tref *= temp_scale
            if mat.St is not None:
                mat.St *= stress_scale
            if mat.Sc is not None:
                mat.Sc *= stress_scale
            if mat.Ss is not None:
                mat.Ss *= stress_scale

        elif mat_type == 'MAT3':
            mat.ex *= stress_scale
            mat.eth *= stress_scale
            mat.ez *= stress_scale
            mat.rho *= density_scale
            if mat.gzx is not None:
                mat.gzx *= stress_scale
            mat.ax *= a_scale
            mat.ath *= a_scale
            mat.az *= a_scale
            mat.tref *= temp_scale

        #elif mat_type == 'MAT4':
            #mat.k
            #mat.cp
            #mat.rho *= density_scale
            #mat.H
            #mat.mu
            #mat.hgen
            #mat.refEnthalpy
            #mat.tch
            #mat.tdelta
            #mat.qlat

        elif mat_type == 'MAT8':
            mat.e11 *= stress_scale
            mat.e22 *= stress_scale
            mat.g12 *= stress_scale
            mat.g1z *= stress_scale
            mat.g2z *= stress_scale
            mat.rho *= density_scale
            mat.a1 *= a_scale
            mat.a2 *= a_scale
            mat.Xt *= stress_scale
            mat.Xc *= stress_scale
            mat.Yt *= stress_scale
            mat.Yc *= stress_scale
            mat.S *= stress_scale

        elif mat_type == 'MAT9':
            mat.G11 *= stress_scale
            mat.G12 *= stress_scale
            mat.G13 *= stress_scale
            mat.G14 *= stress_scale
            mat.G15 *= stress_scale
            mat.G16 *= stress_scale
            mat.G22 *= stress_scale
            mat.G23 *= stress_scale
            mat.G24 *= stress_scale
            mat.G25 *= stress_scale
            mat.G26 *= stress_scale
            mat.G33 *= stress_scale
            mat.G34 *= stress_scale
            mat.G35 *= stress_scale
            mat.G36 *= stress_scale
            mat.G44 *= stress_scale
            mat.G45 *= stress_scale
            mat.G46 *= stress_scale
            mat.G55 *= stress_scale
            mat.G56 *= stress_scale
            mat.G66 *= stress_scale
            mat.rho *= density_scale
            mat.A *= a_scale
            mat.tref *= temp_scale
        elif mat.type == 'MAT3D':
            mat.e1 *= stress_scale
            mat.e2 *= stress_scale
            mat.e3 *= stress_scale
            mat.g12 *= stress_scale
            mat.g13 *= stress_scale
            mat.g23 *= stress_scale
            mat.rho = density_scale
        elif mat.type == 'MAT11':
            mat.e1 *= stress_scale
            mat.e2 *= stress_scale
            mat.e3 *= stress_scale
            mat.g12 *= stress_scale
            mat.g13 *= stress_scale
            mat.g23 *= stress_scale
            mat.rho = density_scale
            #mat.a1 = a1
            #mat.a2 = a2
            #mat.a3 = a3
            mat.tref *= temp_scale
            #mat.ge = ge
        else:
            raise NotImplementedError(mat)

def _convert_constraints(model, xyz_scale):
    """converts the spc/mpcs"""
    for unused_spc_id, spcs in model.spcs.items():
        for spc in spcs:
            if spc.type in ['SPCADD', 'SPC1']:
                continue
            elif spc.type == 'SPC':
                spc.enforced = [enforcedi*xyz_scale for enforcedi in spc.enforced]
            elif spc.type == 'SPCAX':
                spc.d = spc.d * xyz_scale
            else:
                raise NotImplementedError(spc)

def _convert_loads(model, xyz_scale, weight_scale):
    """converts the loads"""
    time_scale = 1.
    frequency_scale = 1. / time_scale
    force_scale = weight_scale
    moment_scale = xyz_scale * weight_scale
    pressure_scale = weight_scale / xyz_scale ** 2
    accel_scale = weight_scale / xyz_scale
    velocity_scale = xyz_scale / time_scale

    if not model.loads:
        return
    model.log.debug('--Load Scales--')
    model.log.debug('force_scale = %s' % force_scale)
    model.log.debug('moment_scale = %s' % moment_scale)
    model.log.debug('pressure_scale = %s' % pressure_scale)
    model.log.debug('accel_scale = %s\n' % accel_scale)

    for dloads in model.dloads.values():
        assert isinstance(dloads, str), dloads  # TEMP

    for dloads in model.dload_entries.values():
        for dload in dloads:
            if dload.type == 'RLOAD1':
                #self.excite_id = excite_id
                #self.delay = delay
                #self.dphase = dphase
                #self.tc = tc
                #self.td = td
                # { P(f) }  = {A} [ C(f) + iD(f)] * e^{  i {\theta - 2 pi f tau } }
                if dload.Type == 'LOAD':
                    scale = force_scale  # moment_scale???
                elif dload.Type == 'DISP':
                    scale = xyz_scale
                elif dload.Type == 'VELO':
                    scale = velocity_scale
                elif dload.Type == 'ACCE':
                    scale = accel_scale
                else:
                    raise RuntimeError(dload)
            else:
                raise NotImplementedError(dload)

    for loads in model.loads.values():
        assert isinstance(loads, list), loads
        for load in loads: # list
            load_type = load.type
            if load_type == 'LOAD':
                pass
            elif load_type in ['FORCE', 'FORCE1', 'FORCE2']:
                load.mag *= force_scale
            elif load_type in ['MOMENT', 'MOMENT1', 'MOMENT2']:
                load.mag *= moment_scale
            elif load_type == 'GRAV':
                load.scale *= accel_scale
            elif load_type == 'ACCEL1':
                load.scale *= accel_scale
            elif load_type == 'PLOAD':
                load.pressure *= pressure_scale
            elif load_type == 'PLOAD1':
                # the errors should never hit
                if load.scale in ['LE', 'LEPR']:
                    if load.Type in ['FX', 'FY', 'FZ', 'FXE', 'FYE', 'FZE']:
                        load.p1 *= force_scale
                        load.p2 *= force_scale
                    elif load.Type in ['MX', 'MY', 'MZ', 'MXE', 'MYE', 'MZE']:
                        load.p1 *= moment_scale
                        load.p2 *= moment_scale
                    else:
                        raise RuntimeError(load)
                elif load.scale in ['FR', 'RFPR']:  # fractional
                    pass
                else:
                    raise RuntimeError(load)
            elif load_type == 'PLOAD2':
                load.pressure *= pressure_scale
            elif load_type == 'PLOAD4':
                load.pressures = [pressure*pressure_scale for pressure in load.pressures]
            elif load_type == 'RANDPS':
                table = load.tid # defines G(f)
                if table.type == 'TABRND1':
                    table.x *= frequency_scale # freq
                    table.y *= force_scale # G
                #elif table.type == 'TABRNDG':
                    #: Scale of turbulence divided by velocity (units of time; Real)
                    #self.LU = LU
                    #: Root-mean-square gust velocity. (Real)
                    #table.WG *= velocity_scale
                else:
                    raise NotImplementedError(table)
            else:
                raise NotImplementedError(load)

def _convert_aero(model, xyz_scale, time_scale, weight_scale):
    """
    Converts the aero cards
      - CAEROx, PAEROx, SPLINEx, AECOMP, AELIST, AEPARAM, AESTAT, AESURF, AESURFS
    """
    if not(model.aecomps or model.aefacts or model.aeparams or model.aelinks or
           model.aelists or model.aestats or model.aesurf or model.aesurfs or
           model.caeros or model.paeros or model.monitor_points or model.splines or
           model.aeros or model.trims or model.divergs):
        return

    area_scale = xyz_scale ** 2
    velocity_scale = xyz_scale / time_scale
    pressure_scale = weight_scale / xyz_scale ** 2
    density_scale = weight_scale / xyz_scale ** 3

    model.log.debug('--Aero Scales--')
    model.log.debug('area_scale = %s' % area_scale)
    model.log.debug('velocity_scale = %s' % velocity_scale)
    model.log.debug('pressure_scale = %s' % pressure_scale)
    model.log.debug('density_scale = %s\n' % density_scale)

    if model.aero:
        model.aero.cref *= xyz_scale
        if model.aero.velocity is not None:
            model.aero.velocity *= velocity_scale

        # we handle density on FLFACTs
        #assert np.allclose(model.aero.rho_ref, 1.0), model.aero

    if model.aeros:
        model.aeros.cref *= xyz_scale
        model.aeros.bref *= xyz_scale
        model.aeros.sref *= area_scale

    for caero in model.caeros.values():
        if caero.type == 'CAERO1':
            caero.p1 *= xyz_scale
            caero.p4 *= xyz_scale
            caero.x12 *= xyz_scale
            caero.x43 *= xyz_scale
        elif caero.type == 'CAERO2':
            caero.p1 *= xyz_scale
            caero.x12 *= xyz_scale
        else:
            raise NotImplementedError('\n' + str(caero))
    #for paero in model.paeros.values():
        #paero.cross_reference(model)
    for trim in model.trims.values():
        trim.q *= pressure_scale
    #for spline in model.splines.values():
        #spline.convert(model)
    #for aecomp in model.aecomps.values():
        #aecomp.cross_reference(model)
    #for aelist in model.aelists.values():
        #aelist.cross_reference(model)
    #for aeparam in model.aeparams.values():
        #aeparam.cross_reference(model)
    #for aestat in model.aestats.values():
        #aestat.cross_reference(model)
    #for aesurf in model.aesurf.values():
        #aesurf.cross_reference(model)
    #for aesurfs in model.aesurfs.values():
        #aesurfs.cross_reference(model)
    for monitor in model.monitor_points:
        if hasattr(monitor, 'xyz'):
            monitor.xyz *= xyz_scale
    # update only the FLFACTs corresponding to density
    flfact_ids = set([])
    for flutter in model.flutters.values():
        flfact = flutter.density
        flfact_ids.add(flfact.sid)
    for flfact_id in flfact_ids: # density
        flfact = model.flfacts[flfact_id]
        flfact.factors *= density_scale

def _convert_optimization(model, xyz_scale, mass_scale, weight_scale):
    """converts the optimization objects"""
    #time_scale = 1.
    #area_scale = xyz_scale ** 2
    #inertia_scale = xyz_scale ** 4
    #force_scale = weight_scale
    #velocity_scale = xyz_scale / time_scale
    pressure_scale = weight_scale / xyz_scale ** 2
    #stiffness_scale = force_scale / xyz_scale
    #damping_scale = force_scale / velocity_scale
    #for key, deqatn in model.dequations.items():
        #deqatn.cross_reference(model)
    #for key, dresp in model.dresps.items():
        #dresp.cross_reference(model)
    #for key, desvar in model.desvars.items():
        #desvar.xinit *= scale
        #desvar.xlb *= scale
        #desvar.xub *= scale
        #desvar.delx *= scale
        #raise NotImplementedError(desvar)

    for unused_key, dconstrs in model.dconstrs.items():
        for dconstr in dconstrs:
            # scale is appled to lid/uid
            _convert_dconstr(model, dconstr, pressure_scale)

    for unused_key, dvcrel in model.dvcrels.items():
        if dvcrel.type == 'DVCREL1':
            scale = _convert_dvcrel1(dvcrel, xyz_scale)
            desvars = dvcrel.dvids_ref
            assert len(desvars) == 1, len(desvars)
            scale_desvars(desvars, scale)
        else:
            raise NotImplementedError(dvcrel)

        desvars = dvcrel.dvids_ref
        assert len(desvars) == 1, len(desvars)
        scale_desvars(desvars, scale)

    for unused_key, dvmrel in model.dvmrels.items():
        raise NotImplementedError(dvmrel)

    for key, dvprel in model.dvprels.items():
        if dvprel.type == 'DVPREL1':
            scale = _convert_dvprel1(dvprel, xyz_scale, mass_scale, weight_scale)
            desvars = dvprel.dvids_ref
            assert len(desvars) == 1, len(desvars)
            scale_desvars(desvars, scale)
        else:
            raise NotImplementedError(dvprel)
        if dvprel.p_max != 1e20:
            dvprel.p_max *= scale
        if dvprel.p_min is not None:
            dvprel.p_min *= scale
        #print('------------')
        #print(dvprel)

def _convert_dconstr(model, dconstr, pressure_scale):
    """helper for ``_convert_optimization``"""
    otype = dconstr.type
    if otype == 'DCONSTR':
        dresp = dconstr.dresp_id_ref
        if dresp.type == 'DRESP1':
            #property_type = dresp.ptype
            response_type = dresp.rtype
            assert len(dresp.atti) == 1, dresp.atti
            for atti in dresp.atti:
                #label = dresp.label
                #atti_type = atti.type
                if response_type == 'STRESS':
                    scale = pressure_scale
                else:
                    raise RuntimeError(atti)
                #if atti
                #if property_type == 'PSHELL':
                    #if rst
        elif dresp.type == 'DRESP2':
            msg = 'skipping:\n%s%s' % (str(dconstr), str(dresp))
            model.log.warning(msg)
            return
        else:
            raise NotImplementedError(dresp)

        # lower bound
        dconstr.lid *= scale
        # upper bound
        dconstr.uid *= scale

        # low end of frequency range (Hz)
        self.lowfq = scale
        # high end of frequency range (Hz)
        self.highfq = scale
    else:
        raise NotImplementedError(dconstr)

def _convert_dvcrel1(dvcrel, xyz_scale):
    """helper for ``_convert_optimization``"""
    element_type = dvcrel.element_type
    if element_type == 'CBUSH':
        if dvcrel.cp_name in ['X1', 'X2', 'X3', 'S', 'S1', 'S2', 'S3']:
            scale = xyz_scale
        else:
            raise NotImplementedError(dvcrel)
    else:
        raise NotImplementedError(dvcrel)
    return scale

def _convert_dvprel1(dvprel, xyz_scale, mass_scale, weight_scale):
    """helper for ``_convert_optimization``"""
    time_scale = 1.
    area_scale = xyz_scale ** 2
    inertia_scale = xyz_scale ** 4
    force_scale = weight_scale
    velocity_scale = xyz_scale / time_scale
    #pressure_scale = weight_scale / xyz_scale ** 2
    stiffness_scale = force_scale / xyz_scale
    damping_scale = force_scale / velocity_scale

    #print(dvprel)
    prop_type = dvprel.prop_type
    var_to_change = dvprel.pname_fid

    scale = 1.
    if prop_type == 'PSHELL':
        if var_to_change == 'T':
            scale = xyz_scale
        elif var_to_change in [6, 8]: # 12I/t^3, ts/t
            scale = 1.
        else:
            raise NotImplementedError(dvprel)
    elif prop_type == 'PCOMP':
        if var_to_change.startswith('THETA'):
            return scale
        if var_to_change.startswith('T'):
            scale = xyz_scale
        else:
            raise NotImplementedError(dvprel)
    elif prop_type == 'PBARL':
        if var_to_change.startswith('DIM'):
            scale = xyz_scale
        else:
            raise NotImplementedError(dvprel)
    elif prop_type == 'PBEAM':
        if isinstance(prop_type, string_types):
            word, unused_num = break_word_by_trailing_parentheses_integer_ab(
                var_to_change)
            if word == 'A':
                scale = area_scale
            elif word in ['I1', 'I2', 'I12', 'J']:
                scale = inertia_scale
            elif word in ['C1', 'C2', 'D1', 'D2', 'E1', 'E2', 'F1', 'F2']:
                scale = xyz_scale
            else:
                raise NotImplementedError(dvprel)
        else:
            raise NotImplementedError(dvprel)

    elif prop_type == 'PBEAML':
        if var_to_change.startswith('DIM'):
            scale = xyz_scale
        else:
            raise NotImplementedError(dvprel)
    elif prop_type == 'PSHEAR':
        if var_to_change == 'T':
            scale = xyz_scale
        else:
            raise NotImplementedError(dvprel)
    elif prop_type == 'PBUSH':
        if var_to_change in ['K1', 'K2', 'K3', 'K4', 'K5', 'K6']:
            scale = stiffness_scale
        elif var_to_change in ['B1', 'B2', 'B3', 'B4', 'B5', 'B6']:
            scale = damping_scale
        else:
            raise NotImplementedError(dvprel)
    elif prop_type == 'PGAP':
        if var_to_change in ['KA', 'KB', 'KT']:
            scale = stiffness_scale
        else:
            raise NotImplementedError(dvprel)
        ##: initial gap opening
        #prop.u0 *= xyz_scale
        ##: preload
        #prop.f0 *= force_scale

    elif prop_type == 'PBUSH1D':
        if var_to_change in ['K']:
            scale = stiffness_scale
        elif var_to_change in ['C']:
            scale = damping_scale
        elif var_to_change in ['M']:
            scale = mass_scale
        else:
            raise NotImplementedError(dvprel)

    elif prop_type == 'PVISC':
        if var_to_change in ['CE1']:
            scale = force_scale / velocity_scale
        else:
            raise NotImplementedError(dvprel)

        #prop.ce *= force_scale / velocity_scale
        #prop.cr *= moment_scale / velocity_scale
    elif prop_type == 'PDAMP':
        if var_to_change in ['B1']:
            scale = force_scale / velocity_scale
        else:
            raise NotImplementedError(dvprel)

    else:
        raise NotImplementedError(dvprel)
    return scale

def scale_desvars(desvars, scale):
    """scales the DVPREL/DVCREL/DVMREL DESVAR values"""
    for desvar in desvars:
        desvar.xinit *= scale
        if desvar.xlb != -1e20:
            desvar.xlb *= scale
        if desvar.xub != 1e20:
            desvar.xub *= scale
        #if desvar.delx != 1e20:
        if desvar.delx is not None and desvar.delx != 1e20:
            desvar.delx *= scale
        if desvar.ddval is not None:
            msg = 'DESVAR id=%s DDVAL is not None\n%s' % str(desvar)
            raise RuntimeError(msg)
        assert desvar.ddval is None, desvar

def get_scale_factors(units_from, units_to):
    """
    [length, mass, time]
    [in, lb, s]
    """
    # in-lb-s
    # m-kg-s
    try:
        length_from, mass_from, time_from = units_from
        length_to, mass_to, time_to = units_to
    except ValueError:
        msg = 'units_from=%s\n' % units_from
        msg += 'units_to=%s' % units_to
        raise ValueError(msg)
    #length_from = units_from[0]
    #length_to = units_to[0]

    #mass_from = units_from[1]
    #mass_to = units_to[1]

    #time_from = units_from[2]
    #time_to = units_to[2]
    assert time_from == time_to, 'units_from=%s units_to=%s time_from=%r time_to=%r' % (
        units_from, units_to, time_from, time_to)
    time_scale = 1.0
    gravity_scale = 1.0
    xyz_scale, gravity_scale_length = convert_length(length_from, length_to)

    mass_scale, weight_scale, gravity_scale_mass = convert_mass(mass_from, mass_to)
    weight_scale /= time_scale ** 2   # doesn't consider cm/mm

    #print('gravity_scale_length=%s gravity_scale_mass=%s gravity_scale=%s' % (
        #gravity_scale_length, gravity_scale_mass, gravity_scale))

    ## doesn't consider cm/mm
    gravity_scale = gravity_scale_length * gravity_scale_mass / time_scale ** 2
    # 4.448N = 1 lbm
    # 1 slug = 14.5939 kg
    # 1g = 32.174 ft/s^2 = 386.088 = 9.80665 m/s^2
    return xyz_scale, mass_scale, time_scale, weight_scale, gravity_scale


def convert_length(length_from, length_to):
    """determines the length scale factor"""
    xyz_scale = 1.0
    gravity_scale = 1.0
    if length_from != length_to:
        if length_from == 'in':
            xyz_scale *= 0.0254
            gravity_scale /= 12.
        elif length_from == 'ft':
            xyz_scale *= 0.3048
            # ft/s^2 are the base units for english

        elif length_from == 'm':
            #xyz_scale = 1.0
            #gravity_scale = 1.0
            # m/s^2 are the base units for SI
            pass
        #elif length_from == 'cm':
            #xyz_scale *= 100.0
        elif length_from == 'mm':
            xyz_scale /= 1000.
            gravity_scale /= 1000.
        else:
            raise NotImplementedError('length from unit=%r; expected=[in, ft, m]' % length_from)

        if length_to == 'in':
            xyz_scale /= 0.0254
            gravity_scale *= 12.
        elif length_to == 'ft':
            xyz_scale /= 0.3048
        elif length_to == 'm':
            #xyz_scale /= 1.0
            pass
        #elif length_to == 'cm':
            #xyz_scale /= 100.0
        elif length_to == 'mm':
            xyz_scale *= 1000.
            gravity_scale *= 1000.
        else:
            raise NotImplementedError('length to unit=%r; expected=[in, ft, m]' % length_to)
    return xyz_scale, gravity_scale


def convert_mass(mass_from, mass_to):
    """determines the mass, weight, gravity scale factor"""
    mass_scale = 1.0
    weight_scale = 1.0
    gravity_scale = 1.0
    #gravity_english = 9.80665 / .3048 #= 32.174; exact
    gravity_english = 32.2

    if mass_from != mass_to:

        # convert to kg
        if mass_from == 'kg':
            pass
        elif mass_from == 'Mg': # mega-gram
            mass_scale *= 1000.
            gravity_scale *= 1000.
        elif mass_from == 'g':
            mass_scale *= 1./1000.
            gravity_scale *= 1/1000.
        elif mass_from == 'lbm':
            mass_scale *= 0.45359237
            weight_scale *= 4.4482216152605
            gravity_scale /= gravity_english # ft/s^2 to m/s^2
        else:
            raise NotImplementedError('mass from unit=%r; expected=[g, kg, Mg, lbm]' % mass_from)

        # convert from kg
        if mass_to == 'kg':
            pass
        elif mass_to == 'Mg': # mega-gram # what about weight_scale???
            mass_scale /= 1000.
            gravity_scale /= 1000.
        elif mass_to == 'g': # what about weight_scale???
            mass_scale *= 1000.
            gravity_scale *= 1000.
            #weight_scale *= 1000.
        elif mass_to == 'lbm':
            mass_scale /= 0.45359237
            weight_scale /= 4.4482216152605
            gravity_scale *= gravity_english
        else:
            raise NotImplementedError('mass to unit=%r; expected=[g, kg, Mg, lbm]' % mass_to)
    return mass_scale, weight_scale, gravity_scale
