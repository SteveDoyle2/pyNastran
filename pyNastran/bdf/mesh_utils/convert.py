# encoding: utf-8
"""
defines:
 - convert(model, units_to, units=None)

"""
from __future__ import annotations
from typing import Tuple, Optional, TYPE_CHECKING

import numpy as np
from pyNastran.bdf.cards.base_card import break_word_by_trailing_parentheses_integer_ab
from pyNastran.bdf.bdf import read_bdf
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.bdf.bdf import (BDF, DVCREL1, DVCREL2, DCONSTR,
                                   PBAR, PBEAM, PBEAM3, PBUSH, PBUSH1D)


def convert(model: BDF, units_to: List[str], units: Optional[List[str]]=None) -> None:
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
    units : List[str]
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


def scale_by_terms(bdf_filename: Union[BDF, str], terms: List[float], scales: List[float],
                   bdf_filename_out: Optional[str]=None,
                   encoding: Optional[str]=None, log=None, debug: bool=True) -> BDF:
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
    assert len(terms) == 3, f'terms={terms} scales={scales}'
    assert len(scales) == 3, f'terms={terms} scales={scales}'
    quiet = not debug
    mass_scale, xyz_scale, time_scale = _setup_scale_by_terms(scales, terms, quiet=quiet)

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

def _setup_scale_by_terms(scales: List[float],
                          terms: List[str], quiet: bool=False) -> Tuple[float, float, float]:
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

def _scale_term(name: str,
                coeffs: List[float],
                terms: List[float],
                scales: List[float]) -> Tuple[float, str]:
    msg = '%s = ' % name
    value = 1.0
    for coeff, term, scale in zip(coeffs, terms, scales):
        if abs(coeff) > 0:
            msg += '%s^%s * ' % (term, coeff)
            value *= scale ** coeff
    msg = msg.strip('* ')
    return value, msg

def scale_model(model: BDF,
                xyz_scale: float,
                mass_scale: float,
                time_scale: float,
                weight_scale: float,
                gravity_scale: float,
                convert_nodes: bool=True, convert_elements: bool=True,
                convert_properties: bool=True, convert_materials: bool=True,
                convert_aero: bool=True, convert_constraints: bool=True,
                convert_loads: bool=True, convert_optimization: bool=True):
    """Performs the model scaling"""
    model.log.debug('L, xyz_scale = %g' % xyz_scale)
    model.log.debug('M, mass_scale = %g' % mass_scale)
    model.log.debug('T, time_scale = %g' % time_scale)
    model.log.debug('F, weight_scale = %g' % weight_scale)
    model.log.debug('G, gravity_scale = %g' % gravity_scale)
    temperature_scale = 1.
    _set_wtmass(model, gravity_scale)

    if convert_nodes:
        _convert_nodes(model, xyz_scale)
        _convert_coordinates(model, xyz_scale)

    if convert_elements:
        _convert_elements(model, xyz_scale, time_scale, mass_scale, weight_scale)
    if convert_properties:
        _convert_properties(model, xyz_scale, time_scale, mass_scale, weight_scale,
                            temperature_scale)
    if convert_materials:
        _convert_materials(model, xyz_scale, mass_scale, weight_scale, temperature_scale)

    if convert_aero:
        _convert_aero(model, xyz_scale, time_scale, weight_scale)
    if convert_constraints:
        _convert_constraints(model, xyz_scale)
    if convert_loads:
        _convert_loads(model, xyz_scale, time_scale, weight_scale, temperature_scale)
    #_convert_sets(model)
    if convert_optimization:
        _convert_optimization(model, xyz_scale, mass_scale, weight_scale, time_scale)


def _set_wtmass(model: BDF, gravity_scale: float) -> None:
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

def _convert_nodes(model: BDF, xyz_scale: bool) -> None:
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

def _convert_coordinates(model: BDF, xyz_scale: float) -> None:
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

def _convert_elements(model: BDF,
                      xyz_scale: float,
                      time_scale: float,
                      mass_scale: float,
                      weight_scale: float) -> None:
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
        'CHBDYG', 'CHBDYE', 'CHBDYP',

        # TODO: NX-verify
        'CTRAX3', 'CTRAX6',
        'CPLSTN3', 'CPLSTN6', 'CPLSTN4', 'CPLSTN8',
        'CQUADX4', 'CQUADX8',

        # acoustic
        'CHACAB',
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
                if elem.T2 is not None:
                    elem.T2 *= xyz_scale
                if elem.T3 is not None:
                    elem.T3 *= xyz_scale
                if elem.T4 is not None:
                    elem.T4 *= xyz_scale
            # nsm
            #elem.nsm *= nsm_scale

        elif elem_type == 'CONROD':
            elem.A *= area_scale # area
            elem.nsm *= nsm_bar_scale
        elif elem_type == 'CGAP':
            if elem.g0 is None and None in elem.x:
                pass
            elif elem.x is not None:  # vector
                elem.x = [x*xyz_scale for x in elem.x]

            #g0 = None
            #x = [None, None, None]
            #if elem.g0 is None and :  # vector
                #print(elem.get_stats())
                #elem.x = [x*xyz_scale for x in elem.x]

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
        elif elem_type == 'CBEAM3':
            #print(elem.get_stats())
            if elem.x is not None:  # vector
                elem.x = [x*xyz_scale for x in elem.x]
            elem.wa *= xyz_scale
            elem.wb *= xyz_scale
            elem.wc *= xyz_scale
            # s
            # tw
        elif elem_type == 'CBEND':
            if elem.x is not None:  # vector
                elem.x = [x*xyz_scale for x in elem.x]
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

def _convert_properties(model: BDF,
                        xyz_scale: float,
                        time_scale: float,
                        mass_scale: float,
                        weight_scale: float,
                        temperature_scale: float) -> None:
    """
    Converts the properties

    Supports:  PELAS, PDAMP, PDAMP5, PVISC, PROD, PBAR, PBARL, PBEAM, PBEAML,
               PSHELL, PSHEAR, PCOMP, PCOMPG, PELAS, PTUBE, PBUSH,
               PCONEAX, PGAP, PBUSH1D
    Skips : PSOLID, PLSOLID, PLPLANE, PIHEX

    Skips are unscaled (intentionally)

    """
    if len(model.properties) == 0:
        return

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
    area_moi_length = area_moi_scale / xyz_scale

    # I, to the bending moment of inertia of a homogeneous shell, T3/12.
    # t^3 has a factor of L^3, so for a factor of 1
    # This is always 1.0.
    #plate_inertia = 1.

    log = model.log
    log.debug('--Property Scales--')
    log.debug('nsm_bar_scale (M/L) = %g' % nsm_bar_scale)
    log.debug('nsm_plate_scale (M/L^2) = %g' % nsm_plate_scale)
    log.debug('stiffness_scale (F/L) = %g' % stiffness_scale)
    log.debug('damping_scale (F/V) = %g' % damping_scale)
    log.debug('stress_scale (F/L^2) = %g\n' % stress_scale)

    skip_properties = {
        'PSOLID', 'PLSOLID', 'PLPLANE', 'PIHEX',

        # TODO: NX-verify
        'PPLANE',

        # acoustic
        'PACABS',
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
            _convert_pbar(prop, xyz_scale, area_scale, area_moi_scale, nsm_bar_scale)

        elif prop_type == 'PBARL':
            prop.dim = [d * xyz_scale for d in prop.dim]
            prop.nsm *= nsm_bar_scale

        elif prop_type == 'PBEAM':
            _convert_pbeam(prop, xyz_scale, area_scale, area_moi_scale, nsm_bar_scale)

        elif prop_type == 'PBEAML':
            prop.dim *= xyz_scale
            prop.nsm *= nsm_bar_scale
        elif prop_type == 'PBEAM3':
            _convert_pbeam3(prop, xyz_scale, area_scale, area_moi_scale, nsm_bar_scale)
        elif prop_type == 'PSHELL':
            prop.t *= xyz_scale
            prop.nsm *= nsm_plate_scale
            prop.z1 *= xyz_scale
            prop.z2 *= xyz_scale
            # prop.twelveIt3  # this is unchanged
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
            _convert_pbush(prop, velocity_scale, mass_scale, stiffness_scale, log)
        elif prop.type == 'PBUSH1D':
            _convert_pbush1d(model, prop, xyz_scale, area_scale,
                             mass_scale, damping_scale, stiffness_scale)

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
                prop.i *= area_moi_length
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

        elif prop_type in ['PBEND', 'PBCOMP', 'PBUSH2D']:
            model.log.warning('skipping %s convert' % prop_type)
        elif prop_type == 'PCOMPS':
            prop.thicknesses *= xyz_scale
            prop.tref *= temperature_scale
        else:
            raise NotImplementedError(f'{prop.get_stats()}\n{prop}')

def _convert_pbar(prop: PBAR,
                  xyz_scale: float,
                  area_scale: float,
                  area_moi_scale: float,
                  nsm_bar_scale: float) -> None:
    """converts a PBAR"""
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

def _convert_pbeam(prop: PBEAM,
                   xyz_scale: float,
                   area_scale: float,
                   area_moi_scale: float,
                   nsm_bar_scale: float) -> None:
    """converts a PBEAM"""
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

def _convert_pbeam3(prop: PBEAM3,
                    xyz_scale: float,
                    area_scale: float,
                    area_moi_scale: float,
                    nsm_bar_scale: float) -> None:
    """converts a PBEAM3"""
    prop.A = [areai * area_scale for areai in prop.A]
    # cw

    prop.cy = [ci * xyz_scale for ci in prop.cy]
    prop.cz = [ci * xyz_scale for ci in prop.cz]
    prop.dy = [ci * xyz_scale for ci in prop.dy]
    prop.dz = [ci * xyz_scale for ci in prop.dz]
    prop.ey = [ci * xyz_scale for ci in prop.ey]
    prop.ez = [ci * xyz_scale for ci in prop.ez]
    prop.fy = [ci * xyz_scale for ci in prop.fy]
    prop.fz = [ci * xyz_scale for ci in prop.fz]

    if prop.i1 is not None:
        prop.i1 = [i1i * area_moi_scale for i1i in prop.i1]
    if prop.i2 is not None:
        prop.i2 = [i2i * area_moi_scale for i2i in prop.i2]
    if prop.i12 is not None:
        prop.i12 = [i12i * area_moi_scale for i12i in prop.i12]
    prop.j = [ji * area_moi_scale for ji in prop.j]
    prop.nsm = [nsmi * nsm_bar_scale for nsmi in prop.nsm]

    #my
    #mz
    #nsiy
    #nsiyz
    #nsiz
    #ny
    #nz
    # w

def _convert_pbush(prop: PBUSH,
                   velocity_scale: float,
                   mass_scale: float,
                   stiffness_scale: float,
                   log: SimpleLogger) -> None:
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
        elif var == 'RCV':
            log.warning('Skipping RCV for PBUSH %i' % prop.pid)
        elif var == 'GE':
            pass
        else:  # pragma: no cover
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

def _convert_pbush1d(model: BDF,
                     prop: PBUSH1D,
                     xyz_scale: float,
                     area_scale: float,
                     mass_scale: float,
                     damping_scale: float,
                     stiffness_scale: float) -> None:
    prop.c *= damping_scale # Viscous damping (force/velocity)
    prop.k *= stiffness_scale
    prop.m *= mass_scale
    prop.sa /= area_scale  # Stress recovery coefficient [1/area]
    prop.se /= xyz_scale   # Strain recovery coefficient [1/length]
    spring_tables = set([])
    damper_tables = set([])
    for var in prop.vars:
        if var == 'SHOCKA':
            print(prop.get_stats())
             # Viscous damping coefficient (force/velocity)
            prop.shock_cvc = prop.shock_cvc * damping_scale if prop.shock_cvc is not None else None
            prop.shock_cvt = prop.shock_cvt * damping_scale if prop.shock_cvt is not None else None
            #shock_exp_vc : 1.0
            #shock_exp_vt : 1.0
            #shock_idecs : None
            #shock_idecsd : None
            #shock_idets : None
            #shock_idetsd : None
            #shock_idts : None
            #shock_type : 'TABLE'
        elif var == 'DAMPER':
            if prop.damper_type == 'TABLE':
                for key in ('damper_idc', 'damper_idcdv', 'damper_idt', 'damper_idtdv'):
                    value = getattr(prop, key)
                    if value is None:
                        continue
                    elif isinstance(value, int):
                        damper_tables.add(value)
                    else:
                        raise TypeError('key=%r value=%r' % (key, value))
            else:
                print(prop.get_stats())
                raise NotImplementedError(prop.damper_type)
        elif var == 'SPRING':
            if prop.spring_type == 'TABLE':
                for key in ('spring_idc', 'spring_idcdu', 'spring_idt', 'spring_idtdu'):
                    value = getattr(prop, key)
                    if value is None:
                        continue
                    elif isinstance(value, int):
                        spring_tables.add(value)
                    else:
                        raise TypeError('key=%r value=%r' % (key, value))
            else:
                raise NotImplementedError(prop.damper_type)
        else:
            print(prop.get_stats())
            raise RuntimeError('var=%r\n%s' % (var, str(prop)))
        if damper_tables:
            model.log.warning(f'scale PBUSH1D damper_tables={damper_tables}')
        if spring_tables:
            model.log.warning(f'scale PBUSH1D spring_tables={spring_tables}')
    #damper_idc : None
    #damper_idcdv : None
    #damper_idt : None
    #damper_idtdv : None
    #damper_type : None
    #type   : 'PBUSH1D'

def _convert_materials(model: BDF,
                       xyz_scale: float,
                       mass_scale: float,
                       weight_scale: float,
                       temperature_scale: float) -> None:
    """
    Converts the materials

    Supports: MAT1, MAT2, MAT3, MAT8, MAT9, MAT10, MAT11

    """
    nmaterials = len(model.materials)
    if nmaterials == 0:
        return

    force_scale = weight_scale
    stress_scale = force_scale / xyz_scale ** 2
    density_scale = mass_scale / xyz_scale ** 3
    a_scale = 1. / temperature_scale # thermal expansion

    model.log.debug('--Material Scales--')
    model.log.debug('density_scale (M/L^3)= %g' % density_scale)
    model.log.debug('stress_scale (F/L^2) = %g\n' % stress_scale)

    for mat in model.materials.values():
        mat_type = mat.type
        if mat_type == 'MAT1':
            mat.e *= stress_scale
            mat.g *= stress_scale
            mat.a *= a_scale
            mat.tref *= temperature_scale
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
            mat.tref *= temperature_scale
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
            mat.tref *= temperature_scale

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
            mat.tref *= temperature_scale

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
            mat.A = [ai * a_scale if ai is not None else None
                     for ai in mat.A]
            mat.tref *= temperature_scale
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
            mat.tref *= temperature_scale
            #mat.ge = ge
        elif mat.type == 'MAT10':
            model.log.warning('skipping %s convert' % mat.type)
        else:
            raise NotImplementedError(mat)

def _convert_constraints(model: BDF, xyz_scale: float) -> None:
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

def _get_dload_scale(dload,
                     xyz_scale: float,
                     velocity_scale: float,
                     accel_scale: float,
                     force_scale: float) -> None:
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

def _convert_loads(model: BDF,
                   xyz_scale: float,
                   time_scale: float,
                   weight_scale: float,
                   temperature_scale: float) -> None:
    """
    Converts the loads

    Supports:
     - dloads: RLOAD1*, TLOAD1*
     - loads:  FORCE, FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2
               GRAV, ACCEL1, PLOAD, PLOAD1, PLOAD2, PLOAD4, RANDPS
     - combinations: DLOAD, LOAD

    * probably not done

    """
    #if not model.loads:
        #return

    frequency_scale = 1. / time_scale
    force_scale = weight_scale
    moment_scale = xyz_scale * force_scale
    pressure_scale = force_scale / xyz_scale ** 2
    accel_scale = force_scale / xyz_scale
    velocity_scale = xyz_scale / time_scale

    model.log.debug('--Load Scales--')
    model.log.debug('force_scale = %s' % force_scale)
    model.log.debug('moment_scale = %s' % moment_scale)
    model.log.debug('pressure_scale = %s' % pressure_scale)
    model.log.debug('accel_scale = %s\n' % accel_scale)

    for dloads in model.dloads.values():
        assert isinstance(dloads, str), dloads  # TEMP

    skip_dloads = {'TLOAD2', 'RANDPS', 'RANDT1', 'QVECT'}
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
            elif dload.type in skip_dloads:
                model.log.warning(f'skipping {dload.type}')
            else:
                raise NotImplementedError(dload)
    for tid, scale in tabled_scales:
        tabled = model.TableD(tid)
        tabled.y *= scale

    skip_cards = {'PLOADX1', 'QVOL', 'QHBDY', 'QBDY1', 'QBDY2', 'QBDY3', 'GMLOAD'}
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
            elif load_type in {'GRAV', 'ACCEL1'}:
                load.scale *= accel_scale
            elif load_type == 'ACCEL':
                load.vals *= accel_scale
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
            #-------------------------------------------------------------
            # thermal
            elif load_type == 'TEMP':
                #print(load.get_stats())
                #temperatures : {1901: 100.0}
                for nid in load.temperatures:
                    load.temperatures[nid] *= temperature_scale
            elif load_type == 'FORCEAX':
                load.f_rtz *= force_scale
            elif load_type in  skip_cards:
                model.log.warning('skipping %s' % load)
            elif load_type == 'TEMPAX':
                load.temperature *= temperature_scale
            elif load_type == 'PRESAX':
                load.pressure *= pressure_scale
            elif load_type == 'TEMPRB':
                load.ta *= temperature_scale
                load.tb *= temperature_scale
                load.tbi = [tempi * temperature_scale for tempi in load.tbi]
                load.tai = [tempi * temperature_scale for tempi in load.tai]
                load.tp1 = [tempi * temperature_scale if tempi is not None else tempi for tempi in load.tp1]
                load.tp2 = [tempi * temperature_scale if tempi is not None else tempi for tempi in load.tp2]
            else:
                raise NotImplementedError(f'{load.get_stats()}\n{load}')

def _convert_aero(model: BDF,
                  xyz_scale: float,
                  time_scale: float,
                  weight_scale: float) -> None:
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
    force_scale = weight_scale
    moment_scale = force_scale * xyz_scale
    pressure_scale = force_scale / xyz_scale ** 2
    density_scale = force_scale / xyz_scale ** 3

    angular_acceleration_scale = 1 / time_scale ** 2  # rad/s^2
    angular_velocity_scale = 1 / time_scale  # rad/s

    model.log.debug('--Aero Scales--')
    model.log.debug('area_scale (L^2) = %s' % area_scale)
    model.log.debug('velocity_scale (L/T) = %s' % velocity_scale)
    model.log.debug('pressure_scale (F/L^2) = %s' % pressure_scale)
    model.log.debug('density_scale (F/L^3) = %s\n' % density_scale)

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
        _scale_caero(caero, xyz_scale, xyz_aefacts)

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

def _scale_caero(caero, xyz_scale: float, xyz_aefacts) -> None:
    try:
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
        elif caero.type == 'CAERO3':
            caero.p1 *= xyz_scale
            caero.p4 *= xyz_scale
            caero.x12 *= xyz_scale
            caero.x43 *= xyz_scale
        elif caero.type == 'CAERO5':
            caero.p1 *= xyz_scale
            caero.p4 *= xyz_scale
            caero.x12 *= xyz_scale
            caero.x43 *= xyz_scale
        else:
            raise NotImplementedError('\n' + str(caero))
    except TypeError:  # pragma: no cover
        print(caero.get_stats())
        raise

def _convert_optimization(model: BDF,
                          xyz_scale: float,
                          mass_scale: float,
                          weight_scale: float,
                          time_scale: float) -> None:
    """
    Converts the optimization objects

    Limited Support: DESVAR, DCONSTR, DVCREL1, DVPREL1

    """
    #time_scale = 1.
    #area_scale = xyz_scale ** 2
    #inertia_scale = xyz_scale ** 4
    force_scale = weight_scale
    #velocity_scale = xyz_scale / time_scale
    pressure_scale = force_scale / xyz_scale ** 2
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
            scale = _convert_dvprel1(dvprel, xyz_scale, mass_scale, weight_scale, time_scale)
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

def _convert_dconstr(model: BDF, dconstr: DCONSTR, pressure_scale: float) -> None:
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

def _convert_dvcrel1(dvcrel: Union[DVCREL1, DVCREL2], xyz_scale: float) -> float:
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

def _convert_dvprel1(dvprel, xyz_scale: float,
                     mass_scale: float,
                     weight_scale: float,
                     time_scale: float) -> float:
    """helper for ``_convert_optimization``"""
    area_scale = xyz_scale ** 2
    inertia_scale = xyz_scale ** 4
    force_scale = weight_scale
    velocity_scale = xyz_scale / time_scale
    #pressure_scale = force_scale / xyz_scale ** 2
    stress_scale = force_scale / xyz_scale ** 2
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
        if var_to_change.startswith('T') or var_to_change == 'Z0':
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
        if isinstance(prop_type, str):
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


def convert_mass(mass_from: str, mass_to: str,
                 log: SimpleLogger) -> Tuple[float, float, float]:
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
        raise NotImplementedError(f'mass from unit={mass_from!r}; '
                                  'expected=[g, kg, Mg, lbm, slug, slinch]')

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
        raise NotImplementedError(f'mass to unit={mass_to!r}; '
                                  'expected=[g, kg, Mg, lbm, slug, slinch]')

    #print("weight_scale = ", mass_from, mass_to, weight_scale)
    return mass_scale, weight_scale, gravity_scale_mass
