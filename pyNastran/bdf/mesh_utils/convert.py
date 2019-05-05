# encoding: utf-8
"""
defines:
 - convert(model, units_to, units=None)

"""
from __future__ import print_function
from six import string_types
import numpy as np
from pyNastran.bdf.cards.base_card import break_word_by_trailing_parentheses_integer_ab
from pyNastran.bdf.bdf import read_bdf


def convert(model, units_to, units=None):
    """
    Converts a model from a set of defined units

    Parameters
    ----------
    model : BDF
       cross references the model (default=True)
    units_to : List[str]
        [length, mass, time]
        length = {in, ft, m, cm, mm}
        mass = {g, kg, Mg, lbm, slug, slinch}
        time = {s}

    units : list
        overwrites model.units

    Note
    ----
    mass refers to the elemental mass, which could be a weight

    """
    # units_start = 'in'
    # units_end = 'mm'
    # xyz_scale = in_to_mm
    # xyz_scale = 25.4
    if units is None:
        units = model.units
    xyz_scale, mass_scale, time_scale, weight_scale, gravity_scale = get_scale_factors(
        units, units_to, model.log)

    scale_model(model, xyz_scale, mass_scale, time_scale, weight_scale, gravity_scale)


def scale_by_terms(bdf_filename, terms, scales, bdf_filename_out=None,
                   encoding=None, log=None, debug=True):
    """
    Scales a BDF based on factors for 3 of the 6 independent terms

    Parameters
    ----------
    bdf_filename : str / BDF()
        a BDF filename
    terms : List[str]; length=3
        the names {M, L, T, F, P, V}
        mass, length, time, force, pressure, velocity
    scales : List[float]; length=3
        the scaling factors
    bdf_filename_out : str; default=None
        a BDF filename to write

    Returns
    -------
    model : BDF()
       the scaled BDF

    """
    mass_scale, xyz_scale, time_scale = _setup_scale_by_terms(scales, terms)

    #-------------------------------------------------
    weight_scale = mass_scale * xyz_scale / time_scale ** 2
    #gravity_scale = xyz_scale / time_scale ** 2
    gravity_scale = 1.0

    #cards_to_skip = [
        #'AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2',
        #'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    model = read_bdf(bdf_filename, validate=True, xref=True,
                     punch=False, save_file_structure=False,
                     skip_cards=None, read_cards=None,
                     encoding=encoding, log=log, debug=debug, mode='msc')
    scale_model(model, xyz_scale, mass_scale, time_scale, weight_scale, gravity_scale)

    if bdf_filename_out is not None:
        model.write_bdf(bdf_filename_out)
    return model

def _setup_scale_by_terms(scales, terms, quiet=False):
    """determines the mass, length, time scaling factors"""
    term_to_mlt_map = {
        #      M   L   T
        'M' : [1., 0., 0.],
        'L' : [0., 1., 0.],
        'T' : [0., 0., 1.],

        'F' : [1., 1., -2.],
        'P' : [1., -1., -2.],
        'V' : [0., 1., -1.],
    }
    assert len(terms) == 3, terms
    A = np.zeros((3, 3), dtype='float64')
    for i, term in enumerate(terms):
        mlt = term_to_mlt_map[term]
        #print(term, mlt)
        A[:, i] = mlt

    MLT = ['mass', 'length', 'time']
    for mlt, Ai in zip(MLT, A):
        if np.allclose(np.abs(Ai).max(), 0.):
            raise RuntimeError('%s is not solvable from [%s]' %  (mlt, ', '.join(terms)))

    detA = np.linalg.det(A)
    if detA == 0.0:
        raise RuntimeError('the equations are not independent '
                           '(e.g., length, time, and velocity) '
                           'and cannot determine mass, legnth, and time')
    M = np.linalg.solve(A, [1., 0., 0.])
    L = np.linalg.solve(A, [0., 1., 0.])
    T = np.linalg.solve(A, [0., 0., 1.])
    mass_scale, mass_msg = _scale_term('M', M, terms, scales)
    xyz_scale, xyz_msg = _scale_term('L', L, terms, scales)
    time_scale, time_msg = _scale_term('T', T, terms, scales)
    if not quiet:
        msg = ('MLT:\n'
               '%s\n'
               '%s\n'
               '%s\n' % (mass_msg, xyz_msg, time_msg))
        print(msg)
    return mass_scale, xyz_scale, time_scale

def _scale_term(name, coeffs, terms, scales):
    msg = '%s = ' % name
    value = 1.0
    for coeff, term, scale in zip(coeffs, terms, scales):
        if abs(coeff) > 0:
            msg += '%s^%s * ' % (term, coeff)
            value *= scale ** coeff
    msg = msg.strip('* ')
    return value, msg

def scale_model(model, xyz_scale, mass_scale, time_scale, weight_scale, gravity_scale,
                convert_nodes=True, convert_elements=True,
                convert_properties=True, convert_materials=True,
                convert_aero=True, convert_constraints=True,
                convert_loads=True, convert_optimization=True):
    """Performs the model scaling"""
    model.log.debug('xyz_scale = %s' % xyz_scale)
    model.log.debug('mass_scale = %s' % mass_scale)
    model.log.debug('time_scale = %s' % time_scale)
    model.log.debug('weight_scale = %s' % weight_scale)
    model.log.debug('gravity_scale = %s' % gravity_scale)
    _set_wtmass(model, gravity_scale)

    if convert_nodes:
        _convert_nodes(model, xyz_scale)
        _convert_coordinates(model, xyz_scale)

    if convert_elements:
        _convert_elements(model, xyz_scale, mass_scale, weight_scale)
    if convert_properties:
        _convert_properties(model, xyz_scale, mass_scale, weight_scale)
    if convert_materials:
        _convert_materials(model, xyz_scale, mass_scale, weight_scale)

    if convert_aero:
        _convert_aero(model, xyz_scale, time_scale, weight_scale)
    if convert_constraints:
        _convert_constraints(model, xyz_scale)
    if convert_loads:
        _convert_loads(model, xyz_scale, weight_scale)
    #_convert_sets(model)
    if convert_optimization:
        _convert_optimization(model, xyz_scale, mass_scale, weight_scale)


def _set_wtmass(model, gravity_scale):
    """
    set the PARAM,WTMASS

    ft-lbm-s-lbf-psf     : 1. / 32.2
    in-lbm-s-lbf-psi     : 1. / (32.2*12)
    m-kg-s-N-Pa          : 1.
    mm-Mg-s-N-Mpa        : 1.
    in-slinch-s-lbf-psi  : 1.
    ft-slug-s-lbf-psf    : 1.
    in-slug-s-lbf-psi    : 1 / 12.

    1 slug   * 1 ft/s^2 = 1 lbf
    1 slinch * 1 in/s^2 = 1 lbf
    1 slinch = 12 slug

    F = m*a
    386 lbf = 1 slinch * 386 * in/s^2
    32  lbf = 1 slug * 32 * ft/s^2
    1   N   = 1 kg * 9.8 m/s^2

    1 lbf = g_scale * 1 slinch  * 1 in/s^2
    1 lbf = g_scale * 12 slug   * 1 in/s^2
    --> g_scale = 1/12.

    """
    if 'WTMASS' in model.params:
        param = model.params['WTMASS']
        value0 = param.values[0]
        weight_mass = value0 / gravity_scale
        param.values = [weight_mass]
    else:
        weight_mass = 1. / gravity_scale
        card = ['PARAM', 'WTMASS', weight_mass]
        model.add_card(card, 'PARAM')
    model.log.debug('wtmass = %s' % weight_mass)

def _convert_nodes(model, xyz_scale):
    """
    Converts the nodes

    Supports: GRID
    """
    for node in model.nodes.values():
        if node.cp == 0:
            node.xyz *= xyz_scale
        elif node.cp_ref.type in ['CORD1R', 'CORD2R']:
            node.xyz *= xyz_scale
        else:
            # only scale R
            node.xyz[0] *= xyz_scale

def _convert_coordinates(model, xyz_scale):
    """
    Converts the coordinate systems

    Supports: CORD1x, CORD2x

    """
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
        elif coord.rid_ref.type in ['CORD1R', 'CORD2R']:
            coord.origin *= xyz_scale
            coord.e1 *= xyz_scale
            coord.e2 *= xyz_scale
            coord.e3 *= xyz_scale
        elif coord.rid_ref.type in ['CORD1C', 'CORD1S', 'CORD2C', 'CORD2S']:
            coord.origin[0] *= xyz_scale
            coord.e1 *= xyz_scale
            #raise NotImplementedError(coord)
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
    """
    Converts the elements

    Supports:  CTRIA3, CTRIA6, CTRIAR,  CQUAD4, CQUAD8, CQUADR,
               CELAS2, CELAS4, CDAMP2, CDAMP4, CBUSH,
               CONROD, CBAR, CBEAM, GENEL, CONM2, CMASS4

    Skips : CELAS1, CELAS3, CDAMP3, CDAMP5, CCONEAX,
            CROD, CTUBE, CVISC, CBUSH1D,
            CQUAD, CSHEAR, CTRIAX, CTRIAX6,
            CTETRA, CPENTA, CHEXA, CPYRAM,
            CMASS1, CMASS3,
    NX Skips: CTRAX3, CTRAX6, CPLSTN3, CPLSTN6, CPLSTN4', CPLSTN8, CQUADX4, CQUADX8
    *intentionally

    """
    time_scale = 1.
    force_scale = weight_scale
    area_scale = xyz_scale ** 2
    velocity_scale = xyz_scale / time_scale
    area_moi_scale = xyz_scale ** 4
    mass_moi_scale = mass_scale * xyz_scale ** 2
    #nsm_scale = mass_scale
    nsm_bar_scale = mass_scale / xyz_scale
    stiffness_scale = force_scale / xyz_scale
    flexibility_scale = 1. / stiffness_scale
    damping_scale = force_scale / velocity_scale

    # these don't have any properties
    skip_elements = {
        # nothing to convert (verified)
        'CCONEAX',
        'CELAS1', 'CELAS3',
        'CDAMP1', 'CDAMP3', 'CDAMP5',
        'CVISC', 'CBUSH1D',
        'CROD', 'CTUBE',
        'CSHEAR', 'CQUAD', 'CQUADX', 'CTRIAX', 'CTRIAX6',
        'CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM',
        'CRAC2D', 'CRAC3D',

        # TODO: NX-verify
        'CTRAX3', 'CTRAX6',
        'CPLSTN3', 'CPLSTN6', 'CPLSTN4', 'CPLSTN8',
        'CQUADX4', 'CQUADX8',
    }
    skip_masses = {'CMASS1', 'CMASS3'}

    tri_shells = {'CTRIA3', 'CTRIA6', 'CTRIAR'}
    quad_shells = {'CQUAD4', 'CQUAD8', 'CQUADR'}
    spring_elements = {'CELAS2', 'CELAS4'}
    damper_elements = {'CDAMP2', 'CDAMP4'}

    model.log.debug('--Element Scales--')
    model.log.debug('nsm_bar_scale (L) = %g' % nsm_bar_scale)
    model.log.debug('area_moi_scale = %g' % area_moi_scale)
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

        elif elem_type == 'CONROD':
            elem.A *= area_scale # area
            elem.nsm *= nsm_bar_scale
        elif elem_type == 'CGAP':
            if elem.x is not None:  # vector
                elem.x = [x*xyz_scale for x in elem.x]

        elif elem_type == 'CBAR':
            if elem.x is not None:  # vector
                elem.x = [x*xyz_scale for x in elem.x]
            elem.wa *= xyz_scale
            elem.wb *= xyz_scale
        elif elem_type == 'CBEAM':
            if elem.x is not None:  # vector
                elem.x = [x*xyz_scale for x in elem.x]
            elem.wa *= xyz_scale
            elem.wb *= xyz_scale
        elif elem_type == 'CBUSH':
            if elem.x[0] is not None:  # vector
                elem.x = [x*xyz_scale for x in elem.x]
            if elem.si[0] is not None:  # vector
                elem.si = [sii*xyz_scale for sii in elem.si]

        elif elem_type == 'GENEL':
            # I'm pretty sure [S] this is unitless
            if elem.k is not None:
                elem.k *= stiffness_scale
            if elem.z is not None:
                elem.z *= flexibility_scale
        else:
            raise NotImplementedError('type=%r; elem:\n%s' % (elem.type, elem))

    for elem in model.masses.values():
        elem_type = elem.type
        if elem_type in skip_masses:
            continue
        if elem_type == 'CONM2':
            elem.mass *= mass_scale
            elem.X *= xyz_scale
            # I = m * r^2
            elem.I = [moi * mass_moi_scale for moi in elem.I]
        elif elem.type == 'CONM1':
            elem.mass_matrix *= mass_scale
        elif elem.type == 'CMASS2':
            elem.mass *= mass_scale
        elif elem.type == 'CMASS4':
            elem.mass *= mass_scale
        else:
            raise NotImplementedError(elem)

def _convert_properties(model, xyz_scale, mass_scale, weight_scale):
    """
    Converts the properties

    Supports:  PELAS, PDAMP, PDAMP5, PVISC, PROD, PBAR, PBARL, PBEAM, PBEAML,
               PSHELL, PSHEAR, PCOMP, PCOMPG, PELAS, PTUBE, PBUSH,
               PCONEAX, PGAP,
    Skips : PSOLID, PLSOLID, PLPLANE, PIHEX

    Skips are unscaled (intentionally)

    """
    time_scale = 1.
    force_scale = weight_scale
    moment_scale = force_scale * xyz_scale
    area_scale = xyz_scale ** 2
    area_moi_scale = xyz_scale ** 4
    velocity_scale = xyz_scale / time_scale

    # there are multiple nsm scales (CONM2, bar, plate)
    nsm_bar_scale = mass_scale / xyz_scale
    nsm_plate_scale = mass_scale / xyz_scale ** 2
    stiffness_scale = force_scale / xyz_scale
    damping_scale = force_scale / velocity_scale
    stress_scale = force_scale / xyz_scale ** 2

    model.log.debug('--Property Scales--')
    model.log.debug('nsm_bar_scale (1/L) = %g' % nsm_bar_scale)
    model.log.debug('nsm_plate_scale (1/L^2) = %g' % nsm_plate_scale)
    model.log.debug('stiffness_scale = %g' % stiffness_scale)
    model.log.debug('damping_scale = %g' % damping_scale)
    model.log.debug('stress_scale = %g\n' % stress_scale)

    skip_properties = {
        'PSOLID', 'PLSOLID', 'PLPLANE', 'PIHEX',

        # TODO: NX-verify
        'PPLANE',
    }

    # we don't need to convert PBUSHT, PELAST, PDAMPT
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
            prop.j *= area_moi_scale
            #prop.c ???

        elif prop_type == 'PBAR':
            prop.A *= area_scale
            prop.i1 *= area_moi_scale
            prop.i2 *= area_moi_scale
            prop.i12 *= area_moi_scale
            prop.j *= area_moi_scale
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
            prop.i1 *= area_moi_scale
            prop.i2 *= area_moi_scale
            prop.i12 *= area_moi_scale
            prop.j *= area_moi_scale
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
                prop.i *= area_moi_scale / xyz_scale
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
    """
    Converts the materials

    Supports: MAT1, MAT2, MAT3, MAT8, MAT9, MAT10, MAT11

    """
    force_scale = weight_scale
    stress_scale = force_scale / xyz_scale ** 2

    density_scale = mass_scale / xyz_scale ** 3
    temp_scale = 1.
    a_scale = 1. / temp_scale # thermal expansion

    model.log.debug('--Material Scales--')
    model.log.debug('density_scale (L^3)= %g' % density_scale)
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
    """
    Converts the spc/mpcs

    Supports: SPC1, SPC, SPCAX
    Implicitly supports: MPC, MPCADD, SPCADD

    """
    for unused_spc_id, spcs in model.spcs.items():
        for spc in spcs:
            if spc.type in ['SPCADD', 'SPC1', 'GMSPC']:
                continue
            elif spc.type == 'SPC':
                spc.enforced = [enforcedi*xyz_scale for enforcedi in spc.enforced]
            elif spc.type == 'SPCAX':
                spc.enforced = spc.enforced * xyz_scale
            else:
                raise NotImplementedError(spc)

def _get_dload_scale(dload, xyz_scale, velocity_scale, accel_scale, force_scale):
    """
    LOAD asssumes force
    """
    if dload.Type == 'LOAD':
        scale = force_scale
    elif dload.Type == 'DISP':
        scale = xyz_scale
    elif dload.Type == 'VELO':
        scale = velocity_scale
    elif dload.Type == 'ACCE':
        scale = accel_scale
    else:
        raise RuntimeError(dload)
    return scale

def _convert_loads(model, xyz_scale, weight_scale):
    """
    Converts the loads

    Supports:
     - dloads: RLOAD1*, TLOAD1*
     - loads:  FORCE, FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2
               GRAV, ACCEL1, PLOAD, PLOAD1, PLOAD2, PLOAD4, RANDPS
     - combinations: DLOAD, LOAD

    * probably not done

    """
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

    tabled_scales = set()
    for dloads in model.dload_entries.values():
        for dload in dloads:
            if dload.type == 'RLOAD1':
                #self.excite_id = excite_id
                #self.delay = delay
                #self.dphase = dphase
                #self.tc = tc
                #self.td = td
                # { P(f) }  = {A} [ C(f) + iD(f)] * e^{  i {\theta - 2 pi f tau } }
                scale = _get_dload_scale(dload, xyz_scale, velocity_scale,
                                         accel_scale, force_scale)

                #delay : int/float; default=None
                    #the delay; if it's 0/blank there is no delay
                    #float : delay in units of time
                    #int : delay id
                #dphase : int/float; default=None
                    #the dphase; if it's 0/blank there is no phase lag
                    #float : delay in units of time
                    #int : delay id
                #tc : int/float; default=0
                    #TABLEDi id that defines C(f) for all degrees of freedom in
                    #EXCITEID entry
                #td : int/float; default=0
                    #TABLEDi id that defines D(f) for all degrees of freedom in
                    #EXCITEID entry
                if isinstance(dload.delay, float):
                    dload.delay *= time_scale
                #{P(f)}  = {A} [ C(f)+iD(f)] e^{i {theta - 2*pi*f*tau} }
                if dload.tc > 0:
                    tabled_scales.add((dload.tc, scale))
                if dload.td > 0:
                    tabled_scales.add((dload.td, scale))

                #darea = model.dareas[dload.excite_id]
                #print(darea.get_stats())
            elif dload.type == 'TLOAD1':
                if isinstance(dload.delay, float):
                    dload.delay *= time_scale
                dload.us0 *= xyz_scale
                dload.vs0 *= velocity_scale
                # {P(t)} = {A} ⋅F(t – τ)

                #darea = model.dareas[dload.excite_id]
                #print(darea.get_stats())
                scale = _get_dload_scale(dload, xyz_scale, velocity_scale,
                                         accel_scale, force_scale)
                tabled_scales.add((dload.tid, scale))
            else:
                raise NotImplementedError(dload)
    for tid, scale in tabled_scales:
        tabled = model.TableD(tid)
        tabled.y *= scale

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
            elif load_type == 'ACCEL':
                load.vals *= accel_scale
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
            elif load_type == 'DEFORM':
                load.deformation *= xyz_scale
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
      - CAEROx, PAEROx, SPLINEx, AECOMP, AELIST, AEPARAM, AESURF

    Supports: AERO, AEROS, CAERO1, CAERO2, TRIM*, MONPNT1, FLUTTER FLFACT-rho/vel
              GUST, AESURF, PAERO2
    Skips: PAERO1, AESTAT, AESURFS, AECOMP, AELIST
    Doesn't support:CAERO3-5, PAERO3-5, SPLINEx, AEPARAM,  AELINK, AEPRESS, AEFORCE,
    *probably not done

    """
    if not(model.aecomps or model.aefacts or model.aeparams or model.aelinks or
           model.aelists or model.aestats or model.aesurf or model.aesurfs or
           model.caeros or model.paeros or model.monitor_points or model.splines or
           model.aeros or model.trims or model.divergs or model.gusts):
        return

    area_scale = xyz_scale ** 2
    velocity_scale = xyz_scale / time_scale
    acceleration_scale = xyz_scale / time_scale ** 2
    moment_scale = weight_scale * xyz_scale
    pressure_scale = weight_scale / xyz_scale ** 2
    density_scale = weight_scale / xyz_scale ** 3

    angular_acceleration_scale = 1 / time_scale ** 2  # rad/s^2
    angular_velocity_scale = 1 / time_scale  # rad/s

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

    xyz_aefacts = set()
    for caero in model.caeros.values():
        if caero.type == 'CAERO1':
            caero.p1 *= xyz_scale
            caero.p4 *= xyz_scale
            caero.x12 *= xyz_scale
            caero.x43 *= xyz_scale
        elif caero.type == 'CAERO2':
            caero.p1 *= xyz_scale
            caero.x12 *= xyz_scale
            #: ID of an AEFACT Bulk Data entry for slender body division
            #: points; used only if NSB is zero or blank. (Integer >= 0)
            if caero.lsb > 0:
                xyz_aefacts.add(caero.lsb)

            #: ID of an AEFACT data entry containing a list of division
            #: points for interference elements; used only if NINT is zero
            #: or blank. (Integer > 0)
            if caero.lint > 0:
                xyz_aefacts.add(caero.lint)
        else:
            raise NotImplementedError('\n' + str(caero))

    for paero in model.paeros.values():
        if paero.type in ['PAERO1']:
            continue
        if paero.type == 'PAERO2':
            #: Reference half-width of body and the width of the constant width
            #: interference tube. (Real > 0.0)
            paero.width *= xyz_scale

            if paero.lrsb is not None:  # half-widths of slender body
                xyz_aefacts.add(paero.lrsb)
            if paero.lrib is not None:  # hafl-widths of interference elements
                xyz_aefacts.add(paero.lrib)
            #self.lrsb_ref = None
            #self.lrib_ref = None
    for aefact_id in xyz_aefacts:
        aefact = model.aefacts[aefact_id]
        aefact.fractions *= xyz_scale

    for trim in model.trims.values():
        trim.q *= pressure_scale
        uxs2 = []
        for label, ux in zip(trim.labels, trim.uxs):
            if label in ['URDD1', 'URDD2', 'URDD3']:
                ux *= acceleration_scale
            elif label in ['URDD4', 'URDD5', 'URDD6']:
                ux *= angular_acceleration_scale
            elif label in ['PITCH', 'ROLL', 'YAW']:
                ux *= angular_velocity_scale
            uxs2.append(ux)
        trim.uxs = uxs2

    for gust in model.gusts.values():
        gust.x0 *= xyz_scale
        gust.V *= velocity_scale
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

    if 'Q' in model.params:
        model.params['Q'].value *= pressure_scale

    q_scale_tables = set()
    for aesurf in model.aesurf.values():
        if aesurf.hmllim is not None:
            aesurf.hmllim *= moment_scale
        if aesurf.hmulim is not None:
            aesurf.hmulim *= moment_scale

        if aesurf.tqllim is not None:
            q_scale_tables.add(aesurf.tqllim)
        if aesurf.tqulim is not None:
            q_scale_tables.add(aesurf.tqulim)

    #for aesurfs in model.aesurfs.values():
        #aesurfs.cross_reference(model)
    for monitor in model.monitor_points:
        if hasattr(monitor, 'xyz'):
            monitor.xyz *= xyz_scale

    # update only the FLFACTs corresponding to density/velocity (not kferq)
    flfact_rho_ids = set()
    flfact_velocity_ids = set()
    for flutter in model.flutters.values():
        flfact = flutter.density_ref
        flfact_rho_ids.add(flfact.sid)
        if 'velocity' in flutter.headers: # we skip kfreq
            flfact = flutter.reduced_freq_velocity_ref
            flfact_velocity_ids.add(flfact.sid)
    for flfact_id in flfact_rho_ids: # density
        flfact = model.flfacts[flfact_id]
        flfact.factors *= density_scale
    for flfact_id in flfact_velocity_ids: # velocity
        flfact = model.flfacts[flfact_id]
        flfact.factors *= velocity_scale

def _convert_optimization(model, xyz_scale, mass_scale, weight_scale):
    """
    Converts the optimization objects

    Limited Support: DESVAR, DCONSTR, DVCREL1, DVPREL1

    """
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
            _convert_desvars(desvars, scale)
        else:
            raise NotImplementedError(dvcrel)

        desvars = dvcrel.dvids_ref
        assert len(desvars) == 1, len(desvars)
        _convert_desvars(desvars, scale)

    for unused_key, dvmrel in model.dvmrels.items():
        raise NotImplementedError(dvmrel)

    for unused_key, dvprel in model.dvprels.items():
        if dvprel.type == 'DVPREL1':
            scale = _convert_dvprel1(dvprel, xyz_scale, mass_scale, weight_scale)
            desvars = dvprel.dvids_ref
            assert len(desvars) == 1, len(desvars)
            _convert_desvars(desvars, scale)
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
        dconstr.lowfq = scale
        # high end of frequency range (Hz)
        dconstr.highfq = scale
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
    stress_scale = weight_scale / xyz_scale ** 2
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
        else:  # pragma: no cover
            raise NotImplementedError('cannot convert %r\n%s' % (var_to_change, dvprel))
    elif prop_type == 'PCOMP':
        if var_to_change.startswith('THETA') or var_to_change == 'GE':
            return scale
        elif var_to_change.startswith('T') or var_to_change == 'Z0':
            scale = xyz_scale
        elif var_to_change == 'SB': # Allowable shear stress of the bonding material
            scale = stress_scale
        else:  # pragma: no cover
            raise NotImplementedError('cannot convert %r\n%s' % (var_to_change, dvprel))
    elif prop_type == 'PBARL':
        if var_to_change.startswith('DIM'):
            scale = xyz_scale
        else:  # pragma: no cover
            raise NotImplementedError('cannot convert %r\n%s' % (var_to_change, dvprel))
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
        else:  # pragma: no cover
            raise NotImplementedError('cannot convert %r\n%s' % (var_to_change, dvprel))
    elif prop_type == 'PSHEAR':
        if var_to_change == 'T':
            scale = xyz_scale
        else:  # pragma: no cover
            raise NotImplementedError('cannot convert %r\n%s' % (var_to_change, dvprel))
    elif prop_type == 'PBUSH':
        if var_to_change in ['K1', 'K2', 'K3', 'K4', 'K5', 'K6']:
            scale = stiffness_scale
        elif var_to_change in ['B1', 'B2', 'B3', 'B4', 'B5', 'B6']:
            scale = damping_scale
        else:  # pragma: no cover
            raise NotImplementedError('cannot convert %r\n%s' % (var_to_change, dvprel))
    elif prop_type == 'PGAP':
        if var_to_change in ['KA', 'KB', 'KT']:
            scale = stiffness_scale
        else:  # pragma: no cover
            raise NotImplementedError('cannot convert %r\n%s' % (var_to_change, dvprel))
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
        else:  # pragma: no cover
            raise NotImplementedError('cannot convert %r\n%s' % (var_to_change, dvprel))

    elif prop_type == 'PVISC':
        if var_to_change in ['CE1']:
            scale = force_scale / velocity_scale
        else:  # pragma: no cover
            raise NotImplementedError('cannot convert %r\n%s' % (var_to_change, dvprel))

        #prop.ce *= force_scale / velocity_scale
        #prop.cr *= moment_scale / velocity_scale
    elif prop_type == 'PDAMP':
        if var_to_change in ['B1']:
            scale = force_scale / velocity_scale
        else:  # pragma: no cover
            raise NotImplementedError('cannot convert %r\n%s' % (var_to_change, dvprel))
    elif prop_type == 'PBAR':
        if var_to_change == 'A':
            scale = area_scale
        elif var_to_change in ['I1', 'I2', 'I12', 'J']:
            scale = inertia_scale
        else:  # pragma: no cover
            raise NotImplementedError('cannot convert %r\n%s' % (var_to_change, dvprel))
    elif prop_type == 'PROD':
        if var_to_change == 'A':
            scale = area_scale
        elif var_to_change == 'J':
            scale = inertia_scale
        else:  # pragma: no cover
            raise NotImplementedError('cannot convert %r\n%s' % (var_to_change, dvprel))
    elif prop_type == 'PTUBE':
        if var_to_change in ['OD', 'T']:
            scale = xyz_scale
        #elif var_to_change == 'J':
            #scale = inertia_scale
        else:  # pragma: no cover
            raise NotImplementedError('cannot convert %r\n%s' % (var_to_change, dvprel))
    else:  # pragma: no cover
        raise NotImplementedError('cannot convert %r\n%s' % (prop_type, dvprel))
    return scale

def _convert_desvars(desvars, scale):
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

def get_scale_factors(units_from, units_to, log):
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

    mass_scale, weight_scale, gravity_scale_mass = convert_mass(mass_from, mass_to, log)
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
    """
    Determines the length scale factor

    We crate a gravity_scale_length for any non-standard unit (ft, m)

    """
    xyz_scale = 1.0
    gravity_scale_length = 1.0

    #if length_from != length_to:
    if length_from == 'in':
        xyz_scale *= 0.0254
        gravity_scale_length /= 12.
    elif length_from == 'ft':
        xyz_scale *= 0.3048
        # ft/s^2 are the base units for english

    elif length_from == 'm':
        #xyz_scale = 1.0
        #gravity_scale_length = 1.0
        # m/s^2 are the base units for SI
        pass
    elif length_from == 'cm':
        xyz_scale *= 100.
        gravity_scale_length /= 100.
    elif length_from == 'mm':
        xyz_scale /= 1000.
        gravity_scale_length /= 1000.
    else:
        raise NotImplementedError('length from unit=%r; expected=[in, ft, m, cm, mm]' % length_from)

    if length_to == 'in':
        xyz_scale /= 0.0254
        gravity_scale_length *= 12.
    elif length_to == 'ft':
        xyz_scale /= 0.3048
    elif length_to == 'm':
        #xyz_scale /= 1.0
        pass
    elif length_to == 'cm':
        xyz_scale /= 100.
        gravity_scale_length *= 100.
    elif length_to == 'mm':
        xyz_scale *= 1000.
        gravity_scale_length *= 1000.
    else:
        raise NotImplementedError('length to unit=%r; expected=[in, ft, m, cm, mm]' % length_to)
    return xyz_scale, gravity_scale_length


def convert_mass(mass_from, mass_to, log):
    """
    determines the mass, weight, gravity scale factor

    We apply a gravity_scale_mass for any unit not {kg, slug}.
    Then we convert to N.

    So for SI, if we have kg, we have a base unit, and the length is assumed
    to be m, so we have a consistent system and gravity_scale_mass is 1.0.

    For lbm:
       F = m*a
       1 lbf = 1 lbm * 1 ft/s^2
       32 lbf = 1 slug * 32 ft/s^2
       gscale = 1/g
       F = gscale * m * a
       1 lbf = gscale * 1 lbm * 32 ft/s^2
       --> gscale = 1/32

    For slug:
       F = gscale * m * a
       32 lbf = gscale * 1 slug * 32 ft/s^2
       --> gscale = 1

    For slinch:
       F = gscale * m * a
       386 lbf = gscale * 1 slinch * 12*32 in/s^2
       1 slinch = 12 slug
       12 in = 1 ft
       386 lbf = gscale * 12 slug * 32 ft/s^2
       --> gscale = 1

    TODO: slinch/slug not validated

    """
    mass_scale = 1.0
    weight_scale = 1.0

    # gravity scale is only required if you have a force-based system
    gravity_scale_mass = 1.0

    #gravity_english = 9.80665 / .3048 #= 32.174; exact
    gravity_english_ft = 32.174
    #gravity_english_in = gravity_english_ft * 12. #32.2 * 12.

    slug_to_kg = 14.5939
    slinch_to_kg = 12. * slug_to_kg  # 1 slinch = 12 slug

    #slinch_to_lbf = gravity_english_in
    #slug_to_lbf = gravity_english_ft
    lbf_to_newton = 4.4482216152605
    lbm_to_kg = 0.45359237

    #slug_to_newton = slug_to_lbf * lbf_to_newton
    #slinch_to_newton = slinch_to_lbf * lbf_to_newton

    #if mass_from == mass_to:
        #return mass_scale, weight_scale, gravity_scale

    # convert to kg
    if mass_from == 'kg':
        pass
    elif mass_from == 'Mg': # mega-gram / ton
        mass_scale *= 1000.
        gravity_scale_mass *= 1000.
    elif mass_from == 'g':
        mass_scale *= 1./1000.
        gravity_scale_mass *= 1/1000.
    elif mass_from == 'lbm':
        mass_scale *= lbm_to_kg
        weight_scale *= lbf_to_newton
        gravity_scale_mass /= gravity_english_ft # ft/s^2 to m/s^2
    elif mass_from == 'slug':
        mass_scale *= slug_to_kg
        # assume lbf
        weight_scale *= lbf_to_newton #* gravity_english_ft)
        #log.warning("not scaling force/weight")
        #weight_scale *= lbf_to_newton
        #gravity_scale_mass /= gravity_english_in # in/s^2 to m/s^2
    elif mass_from == 'slinch':
        mass_scale *= slinch_to_kg
        weight_scale *= lbf_to_newton
        #weight_scale *= (slinch_to_lbf * gravity_english_ft)
        #log.warning("not scaling force/weight")
        #gravity_scale_mass /= gravity_english_in # in/s^2 to m/s^2

    else:
        raise NotImplementedError('mass from unit=%r; '
                                  'expected=[g, kg, Mg, lbm, slug, slinch]' % mass_from)

    # convert from kg
    if mass_to == 'kg':
        pass
    elif mass_to == 'Mg': # mega-gram # what about weight_scale???
        mass_scale /= 1000.
        gravity_scale_mass /= 1000.
    elif mass_to == 'g': # what about weight_scale???
        mass_scale *= 1000.
        gravity_scale_mass *= 1000.
        #weight_scale *= 1000.
    elif mass_to == 'lbm':
        mass_scale /= lbm_to_kg
        weight_scale /= lbf_to_newton
        gravity_scale_mass *= gravity_english_ft
    elif mass_to == 'slug':
        mass_scale /= slug_to_kg
        #weight_scale /= (slug_to_lbf * gravity_english_ft)
        weight_scale /= lbf_to_newton
        #gravity_scale_mass *= gravity_english_in # in/s^2 to m/s^2
    elif mass_to == 'slinch':
        mass_scale /= slinch_to_kg
        #weight_scale /= (slinch_to_lbf * gravity_english_ft)
        weight_scale /= lbf_to_newton
    #elif mass_to == 'slinch':
        #mass_scale /= 175.126836
        #gravity_scale_mass *= gravity_english_in # in/s^2 to m/s^2
    else:
        raise NotImplementedError('mass to unit=%r; '
                                  'expected=[g, kg, Mg, lbm, slug, slinch]' % mass_to)

    #print("weight_scale = ", mass_from, mass_to, weight_scale)
    return mass_scale, weight_scale, gravity_scale_mass
