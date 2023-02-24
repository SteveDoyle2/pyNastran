from __future__ import annotations
from typing import Union, Optional, Any, TYPE_CHECKING

import numpy as np
from pyNastran.bdf.cards.aero.utils import (
    points_elements_from_quad_points)

from pyNastran.converters.avl.avl_helper import integer_types, get_spacing, save_wing_elements
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.converters.avl.body import Body

class Surface:
    def __init__(self,
                 name: str,
                 sections: list[Any],
                 nchord: int, chord_spacing: float,
                 nspan: int, span_spacing: float,
                 component: Optional[int]=None,
                 yduplicate: Optional[float]=None,
                 angle: Optional[float]=None,
                 nowake: Optional[bool]=False,
                 noload: Optional[bool]=False):
        self.name = name
        self.component = component
        self.yduplicate = yduplicate
        self.angle = angle

        self.sections = sections
        self.nchord = nchord
        self.chord_spacing = chord_spacing
        self.nspan = nspan
        self.span_spacing =span_spacing
        self.nowake = nowake
        self.noload = noload
        self.scale = np.ones(3)
        self.translate = np.zeros(3)

    def __contains__(self, variable_name: str) -> bool:
        if variable_name == 'name':
            return True
        if variable_name == 'component':
            return self.component is not None
        elif variable_name == 'yduplicate':
            return self.yduplicate is not None
        elif variable_name == 'angle':
            return self.angle is not None
        elif variable_name == 'scale':
            return self.has_scale()
        elif variable_name == 'translate':
            return self.has_translate()
        elif variable_name == 'nowake':
            return self.nowake is True
        elif variable_name == 'noload':
            return self.noload is True
        raise NotImplementedError(variable_name)

    def has_scale(self) -> bool:
        """does the surface have a SCALE field"""
        return not np.array_equal(self.scale, np.ones(3))

    def has_translate(self) -> bool:
        """does the surface have a TRANSLATE field"""
        return not np.array_equal(self.translate, np.zeros(3))

    def write(self) -> str:
        """writes the surface"""
        surface_msg = (
            '#--------------------------------------------------\n'
            'SURFACE\n'
            f'{self.name}\n'
        )
        surface_msg += self._write_chord_span()

        if self.nowake:
            surface_msg += 'NOWAKE\n'
        if self.noload:
            surface_msg += 'NOLOAD\n'

        if self.yduplicate is not None:
            yduplicate = self.yduplicate
            surface_msg += (
                'YDUPLICATE\n'
                f'  {yduplicate}\n'
            )

        if self.component is not None:
            surface_msg += (
                '\n'
                'COMPONENT\n'
                f'  {self.component}\n'
            )

        if self.has_scale():
            xscale, yscale, zscale = self.scale
            surface_msg += (
                '\n'
                'SCALE\n'
                f'  {xscale}   {yscale}   {zscale}\n'
            )

        #if 'body_file' in surface:
            #body_filename = surface['body_file']
            #surface_msg += (
                #'\n'
                #'BFIL\n'
                #f'{body_filename}\n'
            #)
            #del surface['body_file']

        if self.angle is not None:
            surface_msg += (
                '\n'
                'ANGLE\n'
                f'  {self.angle}\n'
            )

        if self.has_translate():
            dx, dy, dz = self.translate
            surface_msg += (
                '\n'
                'TRANSLATE\n'
                f'  {dx}  {dy}  {dz}\n'
            )

        surface_msg += self._write_sections()
        return surface_msg

    def _write_chord_span(self) -> str:
        """writes the size of the surface"""
        nchordwise = self.nchord
        c_space = self.chord_spacing
        nspanwise = self.nspan
        s_space = self.span_spacing
        if isinstance(nchordwise, integer_types) and isinstance(nspanwise, integer_types):
            surface_msg = ('!Nchordwise  Cspace  Nspanwise  Sspace\n'
                           f'{nchordwise}           {c_space}     {nspanwise}         {s_space}\n')

        elif isinstance(nchordwise, integer_types):
            surface_msg = ('!Nchordwise  Cspace  Nspanwise  Sspace\n'
                           f'{nchordwise}    {c_space}\n')
        else:
            raise NotImplementedError(self)
        return surface_msg

    def _write_sections(self) -> str:
        """writes a section"""
        surface_msg = ''
        for isection, section in enumerate(self.sections):
            section_msg = 'SECTION\n'

            xle, yle, zle = section['xyz_LE']
            chord, ainc, nspan, span_spacing = section['section']
            if nspan is None or span_spacing is None:
                section_msg += (
                    '#Xle    Yle    Zle     Chord   Ainc\n'
                    f'{xle}    {yle}    {zle}     {chord}    {ainc}\n'
                )
            else:
                section_msg += (
                    '#Xle    Yle    Zle     Chord   Ainc  Nspanwise  Sspace\n'
                    f'{xle}    {yle}    {zle}     {chord}    {ainc}     {nspan}    {span_spacing}\n'
                )

            for control in section['control']:
                section_msg += control.write()

            if 'is_afile' in section and section['is_afile']:
                afile = section['afile']
                section_msg += (
                    'AFIL\n'
                    f'{afile}\n'
                )

            #!AFILE
            #!a1.dat
            #!CONTROL
            #!flap     1.0  0.81  0. 0. 0.  +1
            surface_msg += section_msg + '#---------------------------\n'

        return surface_msg

    def get_nodes_elements(self,
                           isurface: int,
                           surfaces: list[Union[Surface, Body]],
                           dirname: str,
                           nodes: list[np.ndarray],
                           ipoint: int,
                           line_elements: list[np.ndarray],
                           quad_elements: list[np.ndarray],
                           is_cs_list: list[np.ndarray],
                           log: SimpleLogger) -> int:
        """builds the surface mesh"""
        xyz_scale = self.scale
        dxyz = self.translate
        assert isinstance(xyz_scale, np.ndarray)
        assert isinstance(dxyz, np.ndarray)

        yduplicate = self.yduplicate

        name = self.name
        log.debug("name=%r ipoint=%s" % (name, ipoint))
        #if 'chord' not in surface:
            #log.debug('no chord for %s...' % name)
            #return ipoint

        ipoint = self._get_wing(
            isurface, xyz_scale, dxyz, ipoint, nodes,
            quad_elements, surfaces, is_cs_list, yduplicate, log)
        return ipoint

    def _get_wing(self, isurface: int,
                  xyz_scale: np.ndarray,
                  dxyz: np.ndarray,
                  ipoint: int,
                  nodes: list[np.ndarray],
                  quad_elements: list[np.ndarray],
                  surfaces: list[Union[Surface, Body]],
                  is_cs_list: list[np.ndarray],
                  yduplicate: float,
                  log: SimpleLogger) -> int:
        log.debug('get_wing')
        name = self.name
        nchord = self.nchord
        chord_spacing = self.chord_spacing
        nspan = self.nspan
        span_spacing = self.span_spacing
        sections = self.sections

        span_stations, airfoil_sections, spanwise_distances, nspans = get_airfoils_from_sections(sections, log)
        nspan_total = nspan
        if nspan_total is None:
            nspan_total = sum(nspans)
        log.debug('span_stations %s' % span_stations)


        #for iairfoil, is_afile in enumerate(surface['is_afile']):
            #pass

        #surface['naca']
        #print('naca =', naca)
        #loft_sections = []
        #for naca in airfoils:
        #get_lofted_sections(None)

        assert nchord > 0, nchord
        #assert nspan > 0, nspan
        nsections = len(sections)
        if len(spanwise_distances) == 1:
            nspanwise_panels = [nspan_total]
        else:
            dy = spanwise_distances.sum()
            nspanwise_panels = (spanwise_distances / dy * nspan_total).astype('int32')
            izero = np.where(nspanwise_panels == 0)[0]
            nspanwise_panels[izero] = 1

        assert isinstance(nchord, int), f'name={name!r} nchord={nchord}'
        assert isinstance(nspan_total, int), f'name={name!r} nspan_total={nspan_total}'
        x = get_spacing(nchord, chord_spacing)

        #del surface['sections']
        #print('wing surface:', str(str))

        #print(surface.keys())
        #print('name =', surface['name'])

        for i in range(nsections-1):
            end = (i == nsections - 1)
            section0 = sections[i]

            nspan_global = nspanwise_panels[i]
            section_data = section0['section']
            nspani, span_spacingi = section_data[2:]

            #if 'afile' in section0:
                #del section0['afile']
            #if 'control' in section0:
                #del section0['control']

            section1 = sections[i+1]
            #if 'afile' in section1:
                #del section1['afile']
            #if 'control' in section1:
                #del section1['control']

            if nspan is not None:
                assert nspan_global >= 1, nspani
                y = get_spacing(nspan_global, span_spacing)
            else:
                # nspan / Sspace for last section are ignored
                assert isinstance(nspani, int), nspani
                assert nspani >= 1, nspani
                y = get_spacing(nspani, span_spacingi)

            ipoint = _section_get_nodes_elements(
                isurface, i,
                section0, section1, airfoil_sections,
                x, y, nspan, end, yduplicate,
                xyz_scale, dxyz,
                ipoint, nodes, quad_elements, surfaces, is_cs_list)

        nodes_temp = np.vstack(nodes)
        assert nodes_temp.shape[0] == ipoint, 'nodes_temp.shape=%s ipoint=%s' % (nodes_temp.shape, ipoint)
        return ipoint

    def __repr__(self) -> str:
        msg = (
            f'Surface(name={self.name})'
        )
        return msg

def _section_get_nodes_elements(isurface: int, i: int,
                                section0, section1,
                                airfoil_sections: list[Any],
                                x: np.ndarray,
                                y: np.ndarray,
                                nspan: int,
                                end: bool,
                                yduplicate: float,
                                xyz_scale: np.ndarray,
                                dxyz: np.ndarray,
                                ipoint: int,
                                nodes: list[np.ndarray],
                                quad_elements: list[np.ndarray],
                                surfaces: list[Union[Surface, Body]],
                                is_cs_list: list[np.ndarray]):
    #print(section0)
    #print('*****************')
    #print(section1)
    p1 = np.array(section0['xyz_LE'])
    p4 = np.array(section1['xyz_LE'])
    #Xle,Yle,Zle =  airfoil's leading edge location
    #Chord       =  the airfoil's chord  (trailing edge is at Xle+Chord,Yle,Zle)
    #Ainc        =  incidence angle, taken as a rotation (+ by RH rule) about
                   #the surface's spanwise axis projected onto the Y-Z plane.
    #Nspan       =  number of spanwise vortices until the next section [ optional ]
    #Sspace      =  controls the spanwise spacing of the vortices      [ optional ]

    #section = [chord, ainc]
    #section = [chord, ainc, nspan, span_spacing]
    chord0 = section0['section'][0]
    chord1 = section1['section'][0]

    #print('chords =', chord0, chord1)
    #print('xyz_scale =', xyz_scale)
    #incidence = section[1]
    p2 = p1 + np.array([chord0, 0., 0.])
    p3 = p4 + np.array([chord1, 0., 0.])

    alpha0 = section0['section'][1]
    alpha1 = section1['section'][1]

    if airfoil_sections[i] is not None:
        if not airfoil_sections[i].shape == airfoil_sections[i+1].shape:
            raise RuntimeError('airfoil_sections[%i]=%s airfoil_sections[%i]=%s' % (
                i, airfoil_sections[i].shape,
                i + 1, airfoil_sections[i+1].shape))
        interpolated_stations = interp_stations(
            y, nspan,
            airfoil_sections[i], chord0, alpha0, p1,
            airfoil_sections[i+1], chord1, alpha1, p4, end=end)

    #loft_sections.append(chord0*airfoil_sections[i])
    #loft_sections.append(chord1*airfoil_sections[i+1])

    assert len(x) > 1, x
    point, element = points_elements_from_quad_points(p1, p2, p3, p4,
                                                      x, y, dtype='int32')
    nelements = element.shape[0]
    is_cs = _section_get_is_cs(section0, section1, x, y, nelements)
    #dxyz[1] = 0.

    ipoint = save_wing_elements(
        isurface, point, element,
        xyz_scale, dxyz,
        nodes, quad_elements, surfaces,
        is_cs_list, is_cs,
        ipoint)

    nodes_temp = np.vstack(nodes)
    assert nodes_temp.shape[0] == ipoint, 'nodes_temp.shape=%s ipoint=%s' % (nodes_temp.shape, ipoint)

    #point2 = None
    #element2 = None
    #print("p1[%i]  = %s" % (i, p1[:2]))
    if yduplicate is not None:
        assert np.allclose(yduplicate, 0.0), 'yduplicate=%s' % yduplicate
        p1[1] *= -1
        p2[1] *= -1
        p3[1] *= -1
        p4[1] *= -1

        # dxyz2 has to be calculated like this because dxyz is global to the surface
        # and we need a mirrored dxyz
        dxyz2 = np.array([dxyz[0], -dxyz[1], dxyz[2]])

        point2, element2 = points_elements_from_quad_points(p1, p2, p3, p4,
                                                            x, y, dtype='int32')
        ipoint = save_wing_elements(
            isurface, point2, element2,
            xyz_scale, dxyz2,
            nodes, quad_elements, surfaces,
            is_cs_list, is_cs,
            ipoint)
        nodes_temp = np.vstack(nodes)
        assert nodes_temp.shape[0] == ipoint, 'nodes_temp.shape=%s ipoint=%s' % (nodes_temp.shape, ipoint)

    #for e in elements:
        #print("  ", e)
    #print('npoint=%s nelement=%s' % (npoint, nelement2))
    #break
    #if not section['afile']:
    #del section['afile']
    #if not section['control']:
    #del section['control']
    #print('  ', section)
    return ipoint

def get_airfoils_from_sections(sections, log: SimpleLogger) -> tuple[np.ndarray,
                                                                     list[Optional[np.ndarray]],
                                                                     np.ndarray,
                                                                     list[int]]:
    nspans = []
    airfoil_sections = []
    is_airfoil_defined = False
    span_stations = np.arange(len(sections))
    leading_edges = []
    trailing_edges = []
    for isection, section in enumerate(sections):
        leading_edge = section['xyz_LE']
        section_data = section['section']
        # chord, ainc, nspan, sspan
        chord = section_data[0]
        nspan = section_data[2]
        trailing_edge = leading_edge + np.array([chord, 0., 0.])
        leading_edges.append(leading_edge)
        trailing_edges.append(trailing_edge)
        nspans.append(nspan)

        log.debug(str(section))
        if 'is_afile' in section:
            is_afile = section['is_afile']
            is_airfoil_defined = True
        else:
            assert is_airfoil_defined is False, is_airfoil_defined
            airfoil_sections.append(None)
            continue

        if is_afile:
            xy = None
        else:
            naca = section['naca']
            xy = get_naca_4_series(log, naca=naca)
        airfoil_sections.append(xy)

    leading_edges = np.array(leading_edges)
    trailing_edges = np.array(trailing_edges)
    quarter_chord = 0.75 * leading_edges + 0.25 * trailing_edges
    # array([[-3.41 ,  0.   ,  0.   ],
    #        [-3.25 , 18.   ,  0.   ],
    #        [-2.5  , 41.66 ,  4.25 ],
    #        [-1.788, 55.75 ,  9.38 ],
    #        [-0.95 , 57.64 , 10.064],
    #        [ 0.   , 58.3  , 10.4  ]])
    dedge = quarter_chord[1:, :] - quarter_chord[:-1, :]
    # array([[ 0.16 , 18.   ,  0.   ],
    #        [ 0.75 , 23.66 ,  4.25 ],
    #        [ 0.712, 14.09 ,  5.13 ],
    #        [ 0.838,  1.89 ,  0.684],
    #        [ 0.95 ,  0.66 ,  0.336]])
    distances = np.linalg.norm(dedge, axis=1)
    # array([18.0007111 , 24.0503763 , 15.01172688,  2.17765929,  1.20457295])

    return span_stations, airfoil_sections, distances, nspans

def get_naca_4_series(log: SimpleLogger, naca: str='2412') -> np.ndarray:
    """
    m=max_camber=2%
    p=located at 40%
    t=max thickness=12%
    """
    log.debug('naca airfoil=%s' % naca)
    t = int(naca[2:]) / 100.
    m = int(naca[0]) / 100.
    p = int(naca[1]) / 10.
    log.debug('t=%s m=%s p=percent_of_max_camber=%s' % (t, m, p))

    # setup the x locations
    if p > 0.:
        # xa = x/chord before the location of max camber
        # xb = x/chord after the location of max camber
        xa = np.linspace(0., p, num=4, endpoint=False, retstep=False, dtype=None)
        xb = np.linspace(p, 1., num=6, endpoint=True, retstep=False, dtype=None)
        x = np.hstack([xa, xb])
    else:
        x = np.linspace(0., 1., num=100, endpoint=True, retstep=False, dtype=None)
        xa = x
        xb = np.array([])
    log.debug('x = %s' % x)

    # https://en.wikipedia.org/wiki/NACA_airfoil
    # t - max thickness in percent chord (the last 2 digits)
    y_half_thickness = 5*t * (0.2969*x**0.5 - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)
    # p - location of max camber (second digit)
    # m - max camber (1st digit)
    #xc2 = xc**2

    if p > 0.0:
        y_camber_a = m/p**2 * (2*p*xa - xa**2)  # if 0 <= x <= pc
        y_camber_b = m/(1-p)**2 * ((1-2*p) + 2*p*xb - xb**2) # pc <= x <= c
        y_camber = np.hstack([y_camber_a, y_camber_b])

        # we're just going to be lazy for now and set theta to 0
        #dy_dx_a = 2*m/p**2 * (p-xa) # if 0 <= x <= pc
        #dy_dx_b = 2*m/(1-p)**2 * (p-xb) # pc <= x <= c
        #theta = np.arctan(dy_dx)
        theta = 0.  #  TODO: wrong
    else:
        y_camber = np.zeros(x.shape)
        theta = 0.

    # thickness is applied perpendicular to the camber line
    x_upper = x - y_half_thickness * np.sin(theta)
    x_lower = x + y_half_thickness * np.sin(theta)
    y_upper = y_camber + y_half_thickness * np.cos(theta)
    y_lower = y_camber - y_half_thickness * np.cos(theta)

    xtotal = np.hstack([x_upper[::-1], x_lower[1:]])
    ytotal = np.hstack([y_upper[::-1], y_lower[1:]])
    #print('x_upper =', x_upper[::-1])
    #print('x_lower =', x_lower[1:])
    #print('xtotal =', xtotal)
    xy = np.vstack([xtotal, ytotal]).T
    #import matplotlib.pyplot as plt
    #plt.figure(1)
    #plt.plot(xtotal, ytotal)
    #plt.grid(True)
    #print(xy)
    return xy

def _section_get_is_cs(section0, section1,
                       x: np.ndarray, y: np.ndarray,
                       nelements: int) -> np.ndarray:
    is_cs = np.zeros(nelements, dtype='int32')
    if len(section0['control']) and len(section1['control']):
        control0 = section0['control'][0]
        control1 = section1['control'][0]
        control0_xhinge = control0.xhinge
        control1_xhinge = control1.xhinge

        # wing in x/c coordinates
        p1 = np.array([0., 0., 0.])
        p2 = np.array([1., 0., 0.])
        p3 = np.array([1., 1., 0.])
        p4 = np.array([0., 1., 0.])
        point_quad, element_quad = points_elements_from_quad_points(
            p1, p2, p3, p4,
            x, y, dtype='int32')
        n1 = element_quad[:, 0]
        n2 = element_quad[:, 1]
        n3 = element_quad[:, 2]
        n4 = element_quad[:, 3]
        # del point_quad
        xyz1 = point_quad[n1, :]
        xyz2 = point_quad[n2, :]
        xyz3 = point_quad[n3, :]
        xyz4 = point_quad[n4, :]
        centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4.
        xcentroid = centroid[:, 0]
        ycentroid = centroid[:, 1]
        xhinge_required = ycentroid * (control1_xhinge - control0_xhinge) + control0_xhinge
        i_cs = xcentroid > xhinge_required
        is_cs[i_cs] = 1
    return is_cs

def interp_stations(y, unused_nspan,
                    airfoil_section0, chord0, alpha0, xyz_le0,
                    airfoil_section1, chord1, alpha1, xyz_le1, end=True):
    """
    x is t/c (so x)
    y is t/c (so z)
    """
    if not airfoil_section0.shape == airfoil_section1.shape:  # pragma: no cover
        raise RuntimeError('airfoil_section0=%s airfoil_section1=%s' % (
            airfoil_section0.shape, airfoil_section1.shape))
    #import matplotlib.pyplot as plt
    y = np.array([0., 0.5, 1.0])
    # first we scale and rotate the section
    xy0 = airfoil_section0 * chord0
    x0 = xy0[:, 0]
    y0 = xy0[:, 1]

    #plt.figure(2)
    #plt.grid(True)
    #plt.plot(x0, y0, 'ro')
    x0_rotated = xyz_le0[0] + x0 * np.cos(alpha0) - y0 * np.sin(alpha0)
    y0_rotated = xyz_le0[2] + x0 * np.sin(alpha0) + y0 * np.cos(alpha0)
    #xy0_rotated = np.vstack([x0_rotated, y0_rotated])

    xy1 = airfoil_section1 * chord1
    x1 = xy1[:, 0]
    y1 = xy1[:, 1]
    #plt.plot(x1, y1, 'bo-')
    #plt.show()

    x1_rotated = xyz_le1[0] + x1 * np.cos(alpha1) - y1 * np.sin(alpha1)
    y1_rotated = xyz_le1[2] + x1 * np.sin(alpha1) + y1 * np.cos(alpha1)

    #plt.figure(4)
    #plt.plot(x0_rotated, y0_rotated)
    #plt.plot(x1_rotated, y1_rotated)
    #plt.grid(True)

    x0_rotated = xyz_le0[0] + x0
    y0_rotated = xyz_le0[2] + y0

    x1_rotated = xyz_le1[0] + x1
    y1_rotated = xyz_le1[2] + y1
    #xy1_rotated = np.vstack([x1_rotated, y1_rotated])

    end = True
    if not end:
        y = y[:-1]

    #print(y.shape) # 3
    #print(x0_rotated.shape)  #
    #y2 = y[np.newaxis, :] + 1
    #print(y2)
    # use linear interpolation to calculate the interpolated stations

    #x_final = y[:, np.newaxis] * x0_rotated * (1.-y[:, np.newaxis]) * x1_rotated
    #y_final = y[:, np.newaxis] * y0_rotated * (1.-y[:, np.newaxis]) * y1_rotated
    #print(x_final.shape)
    x_final = []
    y_final = []
    #plt.figure(5)
    #plt.grid(True)
    for yi in y:
        x_finali = yi * x0_rotated + (1.-yi) * x1_rotated
        y_finali = yi * y0_rotated + (1.-yi) * y1_rotated
        #plt.plot(x_finali, y_finali)
        x_final.append(x_finali)
        y_final.append(y_finali)
    x_final = np.array(x_final)
    y_final = np.array(y_final)

    #plt.show()
    # (nspan, nchord, 2) -> (2, nsan, )
    # (3, 11, 2) -> (2, 3, 11)
    interpolated_stations = np.dstack([x_final, y_final])#.swapaxes(0, 1).swapaxes(0, 2)
    #print(x_final.shape)
    #print(xy_final.shape)
    return interpolated_stations

