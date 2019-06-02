import os
import sys
import copy

import numpy as np
from cpylog import get_logger2
from pyNastran.bdf.cards.aero.utils import (
    points_elements_from_quad_points, create_axisymmetric_body)

AVL_KEYWORDS_LONG = [
    'SURFACE', 'COMPONENT', 'YDUPLICATE', 'BODY', 'ANGLE', 'BFILE',
    'NOWAKE', 'NOLOAD',
    'NACA', 'SCALE', 'TRANSLATE', 'SECTION', 'AFILE', 'CONTROL',
    # unsupported
    #'NOLABE', 'DESIGN', 'AIRFOIL', 'CLAF', 'CDCL',
]

AVL_KEYWORDS = [_avl_keyword[:4] for _avl_keyword in AVL_KEYWORDS_LONG]

def read_avl(avl_filename, log=None, debug=False):
    """reads a *.avl file"""
    avl = AVL(log=log, debug=debug)
    avl.read_avl(avl_filename)
    return avl

class AVL:
    """Interface to the AVL (Athena Vortex Lattice) code"""
    def __init__(self, log=None, debug=False):
        self.name = 'model_name'
        self.log = get_logger2(log=log, debug=debug, encoding='utf-8')
        self.avl_filename = ''
        self.mach = 0.

        self.sref = 0.
        self.cref = 0.
        self.bcref = 0.

        self.xref = 0.
        self.yref = 0.

        self.zref = 0.
        self.cd0 = 0.
        self.sections = []
        self.surfaces = None
        self.sym_iy = None
        self.sym_iz = None
        self.symz = None

    def read_avl(self, avl_filename):
        """only the first 4 chancters are read...not in this reader"""
        self.avl_filename = os.path.abspath(avl_filename)
        with open(avl_filename, 'r') as avl_file:
            lines = avl_file.readlines()

        # remove comments and whitespace lines
        lines = [line.split('#')[0].split('!')[0].rstrip() for line in lines
                 if line.split('#')[0].split('!')[0].strip()]

        self.name = lines[0].strip()
        sline1 = lines[1].strip().split()
        self.mach = float(sline1[0])
        #iYsym =  1  case is symmetric about Y=0    , (X-Z plane is a solid wall)
        #      = -1  case is antisymmetric about Y=0, (X-Z plane is at const. Cp)
        #      =  0  no Y-symmetry is assumed

        #iZsym =  1  case is symmetric about Z=Zsym    , (X-Y plane is a solid wall)
        #      = -1  case is antisymmetric about Z=Zsym, (X-Y plane is at const. Cp)
        #      =  0  no Z-symmetry is assumed (Zsym ignored)

        sline2 = lines[2].split()
        sym_iy, sym_iz, symz = sline2[:3]

        sline3 = lines[3].split()
        sref, cref, bcref = sline3[:3]

        sline4 = lines[4].split()
        xref, yref, zref = sline4[:3]

        line5 = lines[5].strip()
        self.sym_iy = int(sym_iy)
        self.sym_iz = int(sym_iz)

        self.symz = float(symz)

        self.sref = float(sref)
        self.cref = float(cref)
        self.bcref = float(bcref)

        self.xref = float(xref)
        self.yref = float(yref)
        self.zref = float(zref)

        i = 5
        if not line5.isalpha(): # not in AVL_KEYWORDS:
            sline5 = line5.split()
            self.cd0 = float(sline5[0])
            i += 1

        surfaces = []
        surface = {}
        sections = []
        while i < len(lines):
            line = lines[i]
            header = line[:4]
            #print("line = ", line)
            assert header in AVL_KEYWORDS, line
            if is_header(line, 'SURFACE'):
                #print('line: %r' % line)
                if surface:
                    surface = {}
                    sections = []
                    #xyz_les = []
                #print('---------------')
                #print
                name = lines[i+1]
                surface['name'] = name

                #print('name = %r' % name)
                #print(lines[i+2])
                sline2 = lines[i+2].split()
                if len(sline2) == 4:
                    nchord, chord_spacing, nspan, span_spacing = sline2
                    self.log.debug('name=%s nchord=%s chord_spacing=%s nspan=%s span_spacing=%s' % (
                        name, nchord, chord_spacing, nspan, span_spacing))
                    nspan = int(nspan)
                    span_spacing = float(span_spacing)
                    surface['span'] = (nspan, span_spacing)
                elif len(sline2) == 2:
                    nchord, chord_spacing = sline2
                    self.log.debug('name=%s nchord=%s chord_spacing=%s' % (
                        name, nchord, chord_spacing))
                else:
                    raise NotImplementedError(sline2)
                nchord = int(nchord)
                chord_spacing = float(chord_spacing)
                surface['chord'] = (nchord, chord_spacing)
                surface['sections'] = sections
                surfaces.append(surface)
                #surface['xyz_LE'] = xyz_les

                i += 3
                #print()
                #print('---------------')
                continue
            elif is_header(line, 'COMPONENT'):
                ncomponents = int(lines[i+1])
                assert ncomponents == 1, ncomponents
                assert 'component' not in surface
                surface['component'] = 1
                i += 1

            elif is_header(line, 'YDUPLICATE'):
                yduplicate = float(lines[i+1])
                assert 'yduplicate' not in surface
                surface['yduplicate'] = yduplicate
                i += 1
            elif line == 'BODY':
                if surface:
                    surface = {}
                    sections = []

                name = lines[i+1].strip()

                sline2 = lines[i+2].strip().split()
                nchord, chord_spacing = sline2

                nchord = int(nchord)
                chord_spacing = float(chord_spacing)

                #print('name = %r' % name)
                surface['name'] = name
                surface['chord'] = (nchord, chord_spacing)
                surface['sections'] = sections
                surfaces.append(surface)

                i += 2

            elif is_header(line, 'ANGLE'):
                angle = float(lines[i+1])
                assert 'angle' not in surface
                surface['angle'] = angle
                i += 1
            elif is_header(line, 'BFILE'):
                body_file = lines[i+1]
                assert 'body_file' not in surface
                surface['body_file'] = body_file
                i += 1

            elif is_header(line, 'NOWAKE'):
                assert 'nowake' not in surface
                surface['nowake'] = True
            elif is_header(line, 'NOLOAD'):
                assert 'noload' not in surface
                surface['noload'] = True

            elif is_header(line, 'SCALE'):
                xscale, yscale, zscale = lines[i+1].split()
                xscale = float(xscale)
                yscale = float(yscale)
                zscale = float(zscale)
                assert 'scale' not in surface
                surface['scale'] = (xscale, yscale, zscale)
                i += 1

            elif is_header(line, 'TRANSLATE'):
                xtranslate, ytranslate, ztranslate = lines[i+1].split()
                xtranslate = float(xtranslate)
                ytranslate = float(ytranslate)
                ztranslate = float(ztranslate)
                assert 'translate' not in surface
                surface['translate'] = (xtranslate, ytranslate, ztranslate)
                i += 1
            elif is_header(line, 'SECTION'):
                section_data = {
                    #'naca' : [],
                    #'afile' : [],
                    #'is_afile' : [],
                    'control' : [],
                    'xyz_LE' : None,
                }
                sline = lines[i+1].split()
                #if len(sline) == 6:
                    ## #Xle Yle Zle Chord    Nspanwise Sspace
                    #xle, yle, zle, chord, nspan, span_spacing = sline
                    #xle = float(xle)
                    #yle = float(yle)
                    #zle = float(zle)
                    #chord = float(chord)
                    #nspan = float(nspan)
                    #span_spacing = float(span_spacing)
                    #assert 'station' not in surface
                    #section = [chord, nspan, span_spacing]
                if len(sline) == 5:
                    # #Xle Yle Zle Chord Ainc
                    xle, yle, zle, chord, ainc = sline
                    xle = float(xle)
                    yle = float(yle)
                    zle = float(zle)
                    chord = float(chord)
                    ainc = float(ainc)
                    section = [chord, ainc]
                elif len(sline) == 7:
                    #Xle Yle  Zle  Chord  Ainc  Nspanwise  Sspace
                    xle, yle, zle, chord, ainc, nspan, span_spacing = sline
                    xle = float(xle)
                    yle = float(yle)
                    zle = float(zle)
                    ainc = float(ainc)
                    chord = float(chord)
                    nspan = float(nspan)
                    span_spacing = float(span_spacing)
                    assert 'station' not in surface
                    section = [chord, ainc, nspan, span_spacing]
                else:  # pragma: no cover
                    #print(sline, len(sline))
                    raise NotImplementedError(sline)
                section_data['section'] = section
                section_data['xyz_LE'] = [xle, yle, zle]
                sections.append(section_data)
                i += 1
            elif line == 'NACA':
                section_data['is_afile'] = False
                section_data['naca'] = lines[i+1]
                i += 1
            elif is_header(line, 'AFILE'):
                section_data['is_afile'] = True
                section_data['afile'] = lines[i+1]
                i += 1
            elif is_header(line, 'CONTROL'):
                sline = lines[i+1].split()
                #print(sline)
                name, afloat, bfloat, cfloat, dfloat, efloat, fint = sline
                section_data['control'].append([name, afloat, bfloat, cfloat, dfloat, efloat, fint])
                i += 1
            else:
                print(line)
            i += 1

        for surface in surfaces:
            sections = surface['sections']
            #if 'sections' in surface:
                #del surface['sections']
            #print('*', surface)

            #if 0:
                #sections = surface['sections']
                #del surface['sections']
                #print(surface)
                #for section in sections:
                    #if not section['afile']:
                        #del section['afile']
                    #if not section['control']:
                        #del section['control']
                    #print('  ', section)
                #print('')
        self.surfaces = surfaces

    def get_nodes_elements(self):
        self.log.debug('get_nodes_elements')
        dirname = os.path.dirname(self.avl_filename)
        nodes = []
        quad_elements = []
        line_elements = []
        ipoint = 0
        surfaces = []


        #print('----')
        for isurface, surface in enumerate(self.surfaces):
            self.log.debug('isurface = %s' % isurface)
            xyz_scale = np.ones(3)
            if 'scale' in surface:
                xyz_scale = np.array(surface['scale'])

            dxyz = np.zeros(3)
            if 'translate' in surface:
                dxyz = np.array(surface['translate'])

            yduplicate = None
            if 'yduplicate' in surface:
                yduplicate = surface['yduplicate']

            if 'name' not in surface:
                self.log.debug('no name...%s'  % surface)

            name = surface['name']
            self.log.debug("name=%r ipoint=%s" % (name, ipoint))
            if 'chord' not in surface:
                self.log.debug('no chord for %s...' % name)
                continue

            if 'span' not in surface:
                self.log.debug('fuselage surface: %s\n' % simplify_surface(surface))
                if nodes:
                    nodes_temp = np.vstack(nodes)
                    assert nodes_temp.shape[0] == ipoint, 'nodes_temp.shape=%s ipoint=%s' % (nodes_temp.shape, ipoint)

                ipoint = get_fuselage(dirname, isurface, surface, xyz_scale, dxyz, yduplicate,
                                      nodes, line_elements, quad_elements, surfaces, ipoint)
                #if npoint == 0:
                    #self.log.info('skipping %s because there are no sections' % surface)
                #ipoint += npoint
                #print("npoint = ", npoint)
                #print('-----------')
                nodes_temp = np.vstack(nodes)
                assert nodes_temp.shape[0] == ipoint, 'nodes_temp.shape=%s ipoint=%s' % (nodes_temp.shape, ipoint)
                #break
                continue

        #def get_wing(isurface, surface, xyz_scale, dxyz, nodes,
                     #quad_elements, surfaces, ipoint):
            self.log.debug('get_wing')
            nchord, unused_chord_spacing = surface['chord']
            nspan, unused_span_spacing = surface['span']
            sections = surface['sections']


            unused_span_stations, airfoil_sections = get_airfoils_from_sections(sections, self.log)
            self.log.debug('unused_span_stations %s' % unused_span_stations)


            #for iairfoil, is_afile in enumerate(surface['is_afile']):
                #pass

            #surface['naca']
            #print('naca =', naca)
            #loft_sections = []
            #for naca in airfoils:
            #get_lofted_sections(None)

            assert nchord > 0, nchord
            assert nspan > 0, nspan
            #nchord = 1 #  breaks b737 independently of nspan if we use the real value
            #nspan = 1 #  breaks b737 independently of nchord if we use the real value
            #nchord = 2
            #nspan = 2
            x = np.linspace(0., 1., num=nchord+1, endpoint=True, retstep=False, dtype=None)
            y = np.linspace(0., 1., num=nspan+1, endpoint=True, retstep=False, dtype=None)
            #print('x =', x)
            #print('y =', y)

            del surface['sections']

            print('wing surface:', simplify_surface(surface))

            #print(surface.keys())
            #print('name =', surface['name'])

            nsections = len(sections)
            for i in range(nsections-1):
                end = i == nsections - 1

                section0 = sections[i]
                #if 'afile' in section0:
                    #del section0['afile']
                #if 'control' in section0:
                    #del section0['control']

                section1 = sections[i+1]
                #if 'afile' in section1:
                    #del section1['afile']
                #if 'control' in section1:
                    #del section1['control']

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

                #dxyz[1] = 0.

                ipoint = save_wing_elements(
                    isurface, point, element,
                    xyz_scale, dxyz,
                    nodes, quad_elements, surfaces,
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
            #print('-----------')
            nodes_temp = np.vstack(nodes)
            assert nodes_temp.shape[0] == ipoint, 'nodes_temp.shape=%s ipoint=%s' % (nodes_temp.shape, ipoint)
            #print('')
            #break
        #print("end ipoint=%s" % (ipoint))

        nodes = np.vstack(nodes)
        #print(nodes.shape)
        quad_elements = np.vstack(quad_elements)
        if line_elements:
            line_elements = np.vstack(line_elements)
        surfaces = np.hstack(surfaces)
        assert len(surfaces) == quad_elements.shape[0]
        return nodes, quad_elements, line_elements, surfaces

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

def save_wing_elements(isurface, point, element,
                       xyz_scale, dxyz,
                       nodes, quad_elements, surfaces,
                       ipoint):
    npoint = point.shape[0]
    nelement2 = element.shape[0]
    #print("  npoint = ", npoint)

    # scaling is applied before translation
    point[:, 0] *= xyz_scale[0]
    point[:, 1] *= xyz_scale[1]
    point[:, 2] *= xyz_scale[2]
    #print(point)
    point += dxyz

    surfaces.append(np.ones(nelement2, dtype='int32') * isurface)

    nodes.append(point)
    quad_elements.append(ipoint + element)

    ipoint += npoint
    return ipoint

def get_lofted_sections(sections):
    """
    1
    1----2   5----6
    |    |   |    |
    4----3   8----7

    sections = [
        [1, 2, 3, 4],
        [5, 6, 7, 8],
    ]
    elements = get_lofted_sections(sections)
    elements
    [1, 2, 6, 5],
    [2, 3, 7, 6],
    [3, 4, 8, 7],
    [4, 1, 5, 8],
    """
    sections = [
        [1, 2, 3, 4],
        [5, 6, 7, 8],
    ]
    elements = []
    nsections = len(sections)
    for isection in range(nsections-1):
        section0 = sections[isection]
        section1 = sections[isection+1]
        ipoint_max = len(section0) - 1
        for ipoint in range(ipoint_max):
            p1 = section0[ipoint]
            p2 = section0[ipoint+1]
            p3 = section1[ipoint+1]
            p4 = section1[ipoint]
            elements.append([p1, p2, p3, p4])
        p1 = section0[ipoint_max]
        p2 = section0[0]
        p3 = section1[0]
        p4 = section1[ipoint_max]
        elements.append([p1, p2, p3, p4])

    elements = np.array(elements, dtype='int32')
    #print(elements)

def get_naca_4_series(log, naca='2412'):
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

def get_fuselage_from_file(dirname, isurface, surface, xyz_scale, dxyz,
                           nodes, quad_elements, surfaces, ipoint, nchord):
    #print(surface)
    body_file = os.path.join(dirname, surface['body_file'])
    # top view/side view
    # x/c vs
    assert os.path.exists(body_file), body_file
    xc_vs_h = np.loadtxt(body_file, skiprows=1)
    xc = xc_vs_h[:, 0]
    h = xc_vs_h[:, 1]

    # reverse the body "airfoil" to start at the nose
    #xc_reversed = xc[::-1]
    #h_reversed = h[::-1]

    #x = h

    h0 = 0.
    # find x/c for h=0
    #xc_nose = np.interp(h0, h_reversed, xc_reversed, left=None, right=None, period=None)
    xc_nose = xc.min()
    xc_end = xc[0]
    #h_end = h[0]

    # this is a likely buggy way to get the full range of x values
    # on a single side of the airfoil (as the x locations are different)
    # and shaped like a parabola
    #
    i_underside_truncated = np.where(h0 < h)[0].max() + 1
    xc_truncate = xc[:-i_underside_truncated]
    h_truncate = h[:-i_underside_truncated]
    xc_truncate_reversed = xc_truncate[::-1]
    h_truncate_reversed = h_truncate[::-1]

    # find x/c from 0 to 1 starting from the nose
    xc2 = np.linspace(xc_nose, xc_end, num=nchord+1, endpoint=True, retstep=False, dtype=None)
    h2 = np.interp(xc2, xc_truncate_reversed, h_truncate_reversed,
                   left=None, right=None, period=None)

    p1 = np.zeros(3) #dxyz
    xstation = xc2
    ystation = 0.
    zstation = 0.
    radii = h2
    aspect_ratio = 1.0
    point, element = create_axisymmetric_body(xstation, ystation, zstation, radii, aspect_ratio, p1)
    #print(point)
    #print(element)
    npoint = point.shape[0]
    nelement2 = element.shape[0]

    # scaling is applied before translation
    point[:, 0] *= xyz_scale[0]
    point[:, 1] *= xyz_scale[1]
    point[:, 2] *= xyz_scale[2]
    #print(point)
    #print('xyz_scale ', xyz_scale)
    #print('dxyz ', dxyz)
    point += dxyz

    surfaces.append(np.ones(nelement2, dtype='int32') * isurface)
    nodes.append(point)
    quad_elements.append(ipoint + element)
    ipoint += npoint
    return ipoint, nelement2

def get_fuselage(dirname, isurface, surface, xyz_scale, dxyz, yduplicate,
                 nodes, unused_line_elements, quad_elements, surfaces, ipoint):
    #print('----------------------------------------')
    #print(surface)
    nchord, unused_chord_spacing = surface['chord']
    #nchord = 1
    assert nchord >= 1, nchord
    x = np.linspace(0., 1., num=nchord+1, endpoint=True, retstep=False, dtype=None)
    y = np.array([0., 1.])
    #assert len(x) == len(y), 'x=%s y=%s' % (x, y)
    sections = surface['sections']
    del surface['sections']

    #print(surface)
    #print(sections)
    nsections = len(sections)
    if nsections == 0:
        # from file
        ipoint, unused_nelement2 = get_fuselage_from_file(
            dirname, isurface, surface, xyz_scale, dxyz, nodes, quad_elements, surfaces,
            ipoint, nchord)
        return ipoint

    for i in range(nsections-1):
        section0 = sections[i]
        if 'afile' in section0:
            del section0['afile']
        if 'control' in section0:
            del section0['control']

        section1 = sections[i+1]
        if 'afile' in section1:
            del section1['afile']
        if 'control' in section1:
            del section1['control']
        #print(section0)
        #print(section1)
        p1 = np.array(section0['xyz_LE'])
        p4 = np.array(section1['xyz_LE'])
        chord0 = section0['section'][0]
        chord1 = section1['section'][0]
        #print('chords =', chord0, chord1)
        #print('xyz_scale =', xyz_scale)
        #incidence = section[1]
        p2 = p1 + np.array([chord0, 0., 0.])
        p3 = p4 + np.array([chord1, 0., 0.])

        point, element = points_elements_from_quad_points(p1, p2, p3, p4, x, y, dtype='int32')
        ipoint = save_wing_elements(
            isurface, point, element,
            xyz_scale, dxyz,
            nodes, quad_elements, surfaces,
            ipoint)

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
                ipoint)
            nodes_temp = np.vstack(nodes)
            assert nodes_temp.shape[0] == ipoint, 'nodes_temp.shape=%s ipoint=%s' % (nodes_temp.shape, ipoint)

        #print('npoint=%s nelement=%s' % (npoint, nelement2))

    return ipoint

def get_airfoils_from_sections(sections, log):
    airfoil_sections = []
    is_airfoil_defined = False
    span_stations = np.arange(len(sections))
    for section in sections:
        log.debug(section)
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
    return span_stations, airfoil_sections

def is_header(line, name):
    """only the first 4 chancters are read, but we're going to ensure all the letters are correct"""
    return line.startswith(name[:4]) and line == name[:len(line)]

def simplify_surface(surface):
    """gets rid of extraneous data from the surface that makes it hard to read"""
    surface2 = copy.deepcopy(surface)

    if 'translate' in surface and np.allclose(surface['translate'], [0., 0., 0.]):
        del surface2['translate']

    if 'scale' in surface and np.allclose(surface['scale'], [1., 1., 1.]):
        del surface2['scale']

    if 'component' in surface and surface['component'] == 1:
        del surface2['component']

    if 'angle' in surface and surface['angle'] == 0.:
        del surface2['angle']

    if 'sections' not in surface2:
        return surface2

    sections = surface2['sections']
    if len(sections) == 0:
        del surface2['sections']
        return surface2

    #sections2 = {}
    for section in sections:
        if not section['control']:
            del section['control']
    surface2['sections'] = sections
    return surface2

def main():  # pragma: no cover
    avl_filename = sys.argv[1]
    model = read_avl(avl_filename)
    model.get_nodes_elements()

if __name__ == '__main__':  # pragma: no cover
    main()
