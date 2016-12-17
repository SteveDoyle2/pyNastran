"""
defines:
 -
"""
from __future__ import print_function
from six import iteritems, itervalues
import numpy as np
from pyNastran.bdf.bdf import PBEAM, PBEAML

def convert(model, units_to, units=None):
    """
    Converts a model from a set of defined units
    Parameters
    ----------
    xref : bool
       cross references the model (default=True)
    units_to : dict
        ???
    units : dict
        overwrites model.units
    TODO: not done...
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

    if 'WTMASS' in model.params:
        param = model.params['WTMASS']
        value0 = param.values[0]
        param.values = [value0 * gravity_scale]
    else:
        card = ['PARAM', 'WTMASS', gravity_scale]
        model.add_card(card, 'PARAM')

    _convert_nodes(model, xyz_scale)
    _convert_coordinates(model, xyz_scale)
    #_convert_coords(model, xyz_scale)

    _convert_elements(model, xyz_scale, mass_scale, weight_scale)
    _convert_properties(model, xyz_scale, mass_scale, weight_scale)
    #_convert_masses(model)
    _convert_materials(model, xyz_scale, mass_scale, weight_scale)

    _convert_aero(model, xyz_scale, time_scale, weight_scale)
    #_convert_constraints(model)
    _convert_loads(model, xyz_scale, weight_scale)
    #_convert_sets(model)
    _convert_optimization(model, xyz_scale, mass_scale, weight_scale)
#def _convert_coords(model, xyz_scale):
    #for coord in model.coords:
        #if coord.type in ['CORD2R', 'CORD2C', 'CORD2S']:
            #coord.origin *= xyz_scale
            #coord.e1 *= xyz_scale
            #coord.e2 *= xyz_scale
            #coord.e3 *= xyz_scale
        #elif coord.type in ['CORD1R', 'CORD1C', 'CORD1S']:
            #coord.origin *= xyz_scale
        #else:
            #raise NotImplementedError(coord)

def _convert_nodes(model, xyz_scale):
    """converts the nodes"""
    for node in itervalues(model.nodes):
        if node.cp_ref.type in ['CORD1R', 'CORD2R']:
            node.xyz *= xyz_scale
        else:
            # only scale R
            node.xyz[0] *= xyz_scale

def _convert_coordinates(model, xyz_scale):
    for cid, coord in iteritems(model.coords):
        if cid == 0:
            continue
        #print(coord.object_methods())
        #print(dir(coord))
        #if coord.type in ['CORD1C', 'CORD1S', 'CORD2C', 'CORD2S']:
        #print('coord.rid =', coord.rid)
        if coord.rid == 0:
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
    area_scale = xyz_scale ** 2
    moi_scale = xyz_scale ** 4
    mass_moi_scale = mass_scale * xyz_scale ** 2
    #nsm_scale = mass_scale
    nsm_bar_scale = mass_scale / xyz_scale
    stiffness_scale = weight_scale / xyz_scale

    # these don't have any properties
    skip_elements = ['CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM', 'CROD', 'CELAS1']
    #skip_elements = ['CELAS2', 'CELAS4']

    tri_shells = ['CTRIA3', 'CTRIAX', 'CTRIAX6']
    quad_shells = ['CQUAD4', 'CQUAD', 'CQUAD8', 'CQUADX', 'CQUADX8']
    spring_elements = [ 'CELAS2', 'CELAS3', 'CELAS4']

    model.log.debug('--Element Scales--')
    model.log.debug('nsm_bar_scale = %g' % nsm_bar_scale)
    model.log.debug('moi_scale = %g' % moi_scale)
    model.log.debug('area_scale = %g' % area_scale)
    model.log.debug('stiffness_scale = %g\n' % stiffness_scale)
    if len(model.masses):
        model.log.debug('mass_moi_scale = %g' % mass_moi_scale)

    for elem in itervalues(model.elements):
        elem_type = elem.type
        if elem_type in skip_elements:
            continue
        if elem_type in spring_elements:
            elem.k *= stiffness_scale
            continue
        elif elem_type in tri_shells:
            # thickness
            elem.zOffset *= xyz_scale
            if elem.TFlag == 0:
                if elem.T1 is not None:
                    elem.T1 *= xyz_scale
                    elem.T2 *= xyz_scale
                    elem.T3 *= xyz_scale

            # nsm
            #elem.nsm *= nsm_scale

        elif elem_type in quad_shells:
            # thickness
            # TFlag=blank/0 - Ti value: Ti is specified thickness
            #                 Ti blank: Ti is PSHELL value
            # TFlag=0
            # TFlag=1 - thicknesses are relative (Ti default=1.0)
            #
            elem.zOffset *= xyz_scale
            if elem.TFlag == 0:
                if elem.T1 is not None:
                    elem.T1 *= xyz_scale
                    elem.T2 *= xyz_scale
                    elem.T3 *= xyz_scale
                    elem.T4 *= xyz_scale
            # nsm
            #elem.nsm *= nsm_scale

        elif elem_type == 'CONROD':
            elem.area *= area_scale
            elem.nsm *= nsm_bar_scale
        elif elem_type == 'CBAR':
            if elem.x is not None:
                # vector
                elem.x = [x*xyz_scale for x in elem.x]
        elif elem_type == 'CBEAM':
            if elem.x is not None:
                # vector
                elem.x = [x*xyz_scale for x in elem.x]
        else:
            raise NotImplementedError('type=%r; elem:\n%s' % (elem.type, elem))

    for elem in itervalues(model.masses):
        elem_type = elem.type
        if elem_type == 'CONM2':
            elem.mass *= mass_scale
            elem.X *= xyz_scale
                # I = m * r^2
            elem.I = [moi * mass_moi_scale for moi in elem.I]
        else:
            raise NotImplementedError(elem)

def _convert_properties(model, xyz_scale, mass_scale, weight_scale):
    """converts the properties"""
    area_scale = xyz_scale ** 2
    moi_scale = xyz_scale ** 4

    # there are multiple nsm scales (CONM2, bar, plate)
    nsm_scale = mass_scale
    nsm_bar_scale = mass_scale / xyz_scale
    nsm_plate_scale = mass_scale / xyz_scale ** 2
    stiffness_scale = weight_scale / xyz_scale
    stress_scale = weight_scale / xyz_scale ** 2


    model.log.debug('--Property Scales--')
    model.log.debug('nsm_bar_scale = %g' % nsm_bar_scale)
    model.log.debug('nsm_plate_scale = %g' % nsm_plate_scale)
    model.log.debug('stiffness_scale = %g' % stiffness_scale)
    model.log.debug('stress_scale = %g\n' % stress_scale)

    skip_properties = ['PSOLID']
    for prop in itervalues(model.properties):
        prop_type = prop.type
        if prop_type in skip_properties:
            continue
        elif prop_type == 'PELAS':
            prop.k *= stiffness_scale

        elif prop_type == 'PROD':
            prop.area *= area_scale
            prop.moi *= moi_scale

        elif prop_type == 'PBAR':
            prop.A *= area_scale
            prop.i1 *= moi_scale
            prop.i2 *= moi_scale
            prop.i12 *= moi_scale
            prop.j *= moi_scale
            prop.nsm *= nsm_bar_scale

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
            dim2 = []
            for dimi in prop.dim:
                dimi2 = [dim*xyz_scale for dim in dimi]
                dim2.append(dimi2)
            prop.dim = dim2
            prop.nsm = [nsm*nsm_bar_scale for nsm in prop.nsm]

        elif prop_type == 'PSHELL':
            prop.t *= xyz_scale
            prop.nsm *= nsm_plate_scale
            prop.z1 *= xyz_scale
            prop.z2 *= xyz_scale
            prop.twelveIt3 /= xyz_scale ** 3
        elif prop_type in ['PCOMP', 'PCOMPG']:
            prop.thicknesses = [t * xyz_scale for t in prop.thicknesses]
            prop.nsm *= nsm_plate_scale
            prop.z0 *= xyz_scale
            prop.sb *= stress_scale
        elif prop_type == 'PELAS':
            prop.k *= stiffness_scale
        else:
            raise NotImplementedError(prop_type)

def _convert_materials(model, xyz_scale, mass_scale, weight_scale):
    """converts the materials"""
    density_scale = mass_scale / xyz_scale ** 3
    stress_scale = weight_scale / xyz_scale ** 2

    model.log.debug('--Material Scales--')
    model.log.debug('density_scale = %g' % density_scale)
    model.log.debug('stress_scale = %g\n' % stress_scale)

    for mat in itervalues(model.materials):
        mat_type = mat.type
        if mat_type == 'MAT1':
            mat.e *= stress_scale
            mat.g *= stress_scale
            mat.rho *= density_scale
            mat.St *= stress_scale
            mat.Sc *= stress_scale
            mat.Ss *= stress_scale

        elif mat_type == 'MAT8':
            mat.e11 *= stress_scale
            mat.e22 *= stress_scale
            mat.g12 *= stress_scale
            mat.g1z *= stress_scale
            mat.g2z *= stress_scale
            mat.rho *= density_scale
            mat.Xt *= stress_scale
            mat.Xc *= stress_scale
            mat.Yt *= stress_scale
            mat.Yc *= stress_scale
            mat.S *= stress_scale
        else:
            raise NotImplementedError(mat)

def _convert_loads(model, xyz_scale, weight_scale):
    """converts the loads"""
    force_scale = weight_scale
    moment_scale = xyz_scale * weight_scale
    pressure_scale = weight_scale / xyz_scale ** 2
    accel_scale = weight_scale / xyz_scale

    if not model.loads:
        return
    model.log.debug('--Load Scales--')
    model.log.debug('force_scale = %s' % force_scale)
    model.log.debug('moment_scale = %s' % moment_scale)
    model.log.debug('pressure_scale = %s' % pressure_scale)
    model.log.debug('accel_scale = %s\n' % accel_scale)

    for loads in itervalues(model.loads):
        assert isinstance(loads, list), loads
        for load in loads: # list
            load_type = load.type
            if load_type == 'LOAD':
                pass
            elif load_type == 'FORCE':
                load.mag *= force_scale
            elif load_type == 'MOMENT':
                load.mag *= moment_scale
            elif load_type == 'GRAV':
                load.scale *= accel_scale
            #elif load_type == 'PLOAD1':
                #load.mag *= moment_scale
            #elif load_type == 'PLOAD2':
                #load.mag *= moment_scale
            elif load_type == 'PLOAD4':
                load.mag *= pressure_scale
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
    #if hasattr(model, 'aero'):
        #aero = model.aero
        print(model.aero.object_attributes())
        model.aero.refc *= xyz_scale
        model.aero.refb *= xyz_scale
        model.aero.sref *= area_scale
        model.aero.velocity *= velocity_scale
        assert np.allclose(model.aero.density, 1.0), model.aero
    if model.aeros:
        #print(aeros)
        #print(aeros.object_attributes())
        model.aeros.cref *= xyz_scale
        model.aeros.bref *= xyz_scale
        model.aeros.sref *= area_scale

    for caero in itervalues(model.caeros):
        if caero.type in ['CAERO1']:
            caero.p1 *= xyz_scale
            caero.p4 *= xyz_scale
            caero.x12 *= xyz_scale
            caero.x43 *= xyz_scale
        else:
            raise NotImplementedError(caero)
    #for paero in itervalues(model.paeros):
        #paero.cross_reference(model)
    for trim in itervalues(model.trims):
        trim.q *= pressure_scale
    #for spline in itervalues(model.splines):
        #spline.convert(model)
    #for aecomp in itervalues(model.aecomps):
        #aecomp.cross_reference(model)
    #for aelist in itervalues(model.aelists):
        #aelist.cross_reference(model)
    #for aeparam in itervalues(model.aeparams):
        #aeparam.cross_reference(model)
    #for aestat in itervalues(model.aestats):
        #aestat.cross_reference(model)
    #for aesurf in itervalues(model.aesurf):
        #aesurf.cross_reference(model)
    #for aesurfs in itervalues(model.aesurfs):
        #aesurfs.cross_reference(model)

    # update only the FLFACTs corresponding to density
    flfact_ids = set([])
    for flutter in itervalues(model.flutters):
        flfact = flutter.density
        flfact_ids.add(flfact.sid)
    for flfact_id in flfact_ids: # density
        flfact = model.flfacts[flfact_id]
        flfact.factors *= density_scale

def _convert_optimization(model, xyz_scale, mass_scale, weight_scale):
    """converts the optimization objects"""
    pressure_scale = weight_scale / xyz_scale ** 2
    #for key, deqatn in iteritems(model.dequations):
        #deqatn.cross_reference(model)
    #for key, dresp in iteritems(model.dresps):
        #dresp.cross_reference(model)
    #for key, desvar in iteritems(model.desvars):
        #desvar.xinit *= scale
        #desvar.xlb *= scale
        #desvar.xub *= scale
        #desvar.delx *= scale
        #raise NotImplementedError(desvar)
    for key, dconstr in iteritems(model.dconstrs):
        otype = dconstr.type
        if otype == 'DCONSTR':
            dresp = dconstr.rid
            if dresp.type == 'DRESP1':
                property_type = dresp.ptype
                response_type = dresp.rtype
                assert len(dresp.atti) == 1, dresp.atti
                for atti in dresp.atti:
                    label = dresp.label
                    atti_type = atti.type
                    if response_type == 'STRESS':
                        scale = pressure_scale
                    else:
                        raise RuntimeError(atti)
                    #if atti
                    #if property_type == 'PSHELL':
                        #if rst
            else:
                raise NotImplementedError(dresp)
            dconstr.lid *= scale
            dconstr.uid *= scale
        else:
            raise NotImplementedError(dconstr)

    for key, dvcrel in iteritems(model.dvcrels):
        raise NotImplementedError(dvcrel)
    for key, dvmrel in iteritems(model.dvmrels):
        raise NotImplementedError(dvmrel)
    for key, dvprel in iteritems(model.dvprels):
        if dvprel.type == 'DVPREL1':
            #print(dvprel)
            prop_type = dvprel.Type
            desvars = dvprel.dvids
            var_to_change = dvprel.pNameFid
            assert len(desvars) == 1, len(desvars)

            if prop_type == 'PSHELL':
                if var_to_change == 'T':
                    scale = xyz_scale
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
                else:
                    raise NotImplementedError(dvprel)
            else:
                raise NotImplementedError(dvprel)
        else:
            raise NotImplementedError(dvprel)
        if dvprel.pMax != 1e20:
            dvprel.pMax *= scale
        if dvprel.pMin is not None:
            dvprel.pMin *= scale
        #print('------------')
        #print(dvprel)
        #pass

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
    xyz_scale = convert_length(length_from, length_to)

    mass_scale, weight_scale, gravity_scale = convert_mass(mass_from, mass_to)
    weight_scale /= time_scale ** 2   # doesn't consider cm/mm
    gravity_scale /= time_scale ** 2  # doesn't consider cm/mm
    # 4.448N = 1 lbm
    # 1 slug = 14.5939 kg
    # 1g = 32.174 ft/s^2 = 386.088 = 9.80665 m/s^2
    return xyz_scale, mass_scale, time_scale, weight_scale, gravity_scale


def convert_length(length_from, length_to):
    """determines the length scale factor"""
    xyz_scale = 1.0
    if length_from != length_to:
        if length_from == 'in':
            xyz_scale *= 0.0254
        elif length_from == 'ft':
            xyz_scale *= 0.3048
        elif length_from == 'm':
            #xyz_scale *= 1.0
            pass
        #elif length_from == 'cm':
            #xyz_scale *= 100.0
        elif length_from == 'mm':
            xyz_scale /= 1000.0
        else:
            raise NotImplementedError('length from unit=%r; expected=[in, ft, m]' % length_from)

        if length_to == 'in':
            xyz_scale /= 0.0254
        elif length_to == 'ft':
            xyz_scale /= 0.3048
        elif length_to == 'm':
            #xyz_scale /= 1.0
            pass
        #elif length_to == 'cm':
            #xyz_scale /= 100.0
        elif length_to == 'mm':
            xyz_scale *= 1000.0
        else:
            raise NotImplementedError('length to unit=%r; expected=[in, ft, m]' % length_to)
    return xyz_scale


def convert_mass(mass_from, mass_to):
    """determines the mass, weight, gravity scale factor"""
    mass_scale = 1.0
    weight_scale = 1.0
    gravity_scale = 1.0
    if mass_from != mass_to:

        # convert to kg
        if mass_from == 'kg':
            pass
        elif mass_from == 'Mg': # mega-gram
            mass_scale *= 1000.
        elif mass_from == 'g':
            mass_scale *= 1./1000.
        elif mass_from == 'lbm':
            mass_scale *= 0.45359237
            weight_scale *= 4.4482216152605
            gravity_scale *= 0.3048
        else:
            raise NotImplementedError('mass from unit=%r; expected=[g, kg, Mg, lbm]' % mass_from)

        # convert from kg
        if mass_to == 'kg':
            pass
        elif mass_to == 'Mg': # mega-gram # what about weight_scale???
            mass_scale /= 1000.
        elif mass_to == 'g': # what about weight_scale???
            mass_scale *= 1000.
            #weight_scale *= 1000.
        elif mass_to == 'lbm':
            mass_scale /= 0.45359237
            weight_scale /= 4.4482216152605
            gravity_scale /= 0.3048
        else:
            raise NotImplementedError('mass to unit=%r; expected=[g, kg, Mg, lbm]' % mass_to)

    return mass_scale, weight_scale, gravity_scale
