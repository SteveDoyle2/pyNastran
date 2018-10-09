from __future__ import print_function
import os
import sys
from collections import OrderedDict
import numpy as np
from pyNastran.bdf.cards.aero.utils import points_elements_from_quad_points, create_axisymmetric_body
from pyNastran.utils.log import get_logger2

AVL_KEYWORDS = [
    'SURFACE', 'COMPONENT', 'YDUPLICATE', 'BODY', 'ANGLE', 'BFIL', 'BFILE', 'NOWAKE', 'NOLOAD',
    'NACA', 'SCALE', 'TRANSLATE', 'SECTION', 'AFIL', 'AFILE', 'CONTROL',
]


def read_avl(avl_filename, log=None, debug=False):
    avl = AVL()
    avl.read_avl(avl_filename)
    return avl

class AVL(object):
    def __init__(self, log=None, debug=False):
        self.name = 'model_name'
        self.log = get_logger2(log=log, debug=debug, encoding='utf-8')
        self.mach = 0.

        self.sref = 0.
        self.cref = 0.
        self.bcref = 0.

        self.xref = 0.
        self.yref = 0.

        self.zref = 0.
        self.cd0 = 0.
        self.sections = []

    def read_avl(self, avl_filename):
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
        if line5 not in AVL_KEYWORDS:
            sline5 = line5.split()
            self.cd0 = float(sline5[0])
            i += 1

        surfaces = []
        surface = {}
        sections = []
        while i < len(lines):
            line = lines[i]
            #print("line = ", line)
            assert line in AVL_KEYWORDS, line
            if line == 'SURFACE':
                print('line: %r' % line)
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
                    self.log.debug('name=%s nchord=%s chord_spacing=%s' % (name, nchord, chord_spacing))
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
            elif line == 'COMPONENT':
                ncomponents = int(lines[i+1])
                assert ncomponents == 1, ncomponents
                assert 'component' not in surface
                surface['component'] = 1
                i += 1

            elif line == 'YDUPLICATE':
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

                print('name = %r' % name)
                surface['name'] = name
                surface['chord'] = (nchord, chord_spacing)
                surface['sections'] = sections
                surfaces.append(surface)

                i += 2

            elif line == 'ANGLE':
                angle = float(lines[i+1])
                assert 'angle' not in surface
                surface['angle'] = angle
                i += 1
            elif line in ['BFILE', 'BFIL']:
                body_file = lines[i+1]
                assert 'body_file' not in surface
                surface['body_file'] = body_file
                i += 1

            elif line == 'NOWAKE':
                assert 'nowake' not in surface
                surface['nowake'] = True
            elif line == 'NOLOAD':
                assert 'noload' not in surface
                surface['noload'] = True

            elif line == 'NACA':
                assert 'NACA' not in surface
                surface['naca'] = lines[i+1]
                i += 1

            elif line == 'SCALE':
                xscale, yscale, zscale = lines[i+1].split()
                xscale = float(xscale)
                yscale = float(yscale)
                zscale = float(zscale)
                assert 'scale' not in surface
                surface['scale'] = (xscale, yscale, zscale)
                i += 1

            elif line == 'TRANSLATE':
                xtranslate, ytranslate, ztranslate = lines[i+1].split()
                xtranslate = float(xtranslate)
                ytranslate = float(ytranslate)
                ztranslate = float(ztranslate)
                assert 'translate' not in surface
                surface['translate'] = (xtranslate, ytranslate, ztranslate)
                i += 1
            elif line == 'SECTION':
                section_data = {
                    'afile' : [],
                    'control' : [],
                    'xyz_LE' : None,
                }
                sline = lines[i+1].split()
                if len(sline) == 6:
                    # #Xle Yle Zle Chord    Nspanwise Sspace
                    xle, yle, zle, chord, nspan, span_spacing = sline
                    xle = float(xle)
                    yle = float(yle)
                    zle = float(zle)
                    chord = float(chord)
                    nspan = float(nspan)
                    span_spacing = float(span_spacing)
                    assert 'station' not in surface
                    section = [chord, nspan, span_spacing]
                elif len(sline) == 5:
                    # #Xle Yle Zle Chord Ainc
                    xle, yle, zle, chord, ainc = sline
                    xle = float(xle)
                    yle = float(yle)
                    zle = float(zle)
                    chord = float(chord)
                    ainc = float(ainc)
                    section = [chord, ainc]
                elif len(sline) == 7:
                    #Xle   Yle    Zle      Chord   Ainc  Nspanwise  Sspace
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
                else:
                    print(sline, len(sline))
                    asdf
                section_data['section'] = section
                section_data['xyz_LE'] = [xle, yle, zle]
                sections.append(section_data)
                i += 1
            elif line in ['AFILE', 'AFIL']:
                section_data['afile'].append(lines[i+1])
                i += 1
            elif line == 'CONTROL':
                sline = lines[i+1].split()
                #print(sline)
                name, afloat, bfloat, cfloat, dfloat, efloat, fint = sline
                section_data['control'].append([name, afloat, bfloat, cfloat, dfloat, efloat, fint])
                i += 1

            #elif line == 'CONTROL':


            else:
                print(line)
            i += 1

        for surface in surfaces:
            sections = surface['sections']
            #if 'sections' in surface:
                #del surface['sections']
            print('*', surface)

            if 0:
                sections = surface['sections']
                del surface['sections']
                print(surface)
                for section in sections:
                    if not section['afile']:
                        del section['afile']
                    if not section['control']:
                        del section['control']
                    print('  ', section)
                print('')
        self.surfaces = surfaces

    def get_nodes_elements(self):
        nodes = []
        quad_elements = []
        line_elements = []
        ipoint = 0
        ielement = 0
        surfaces = []


        for isurface, surface in enumerate(self.surfaces):
            xyz_scale = np.ones(3)
            if 'scale' in surface:
                xyz_scale = np.array(surface['scale'])
            dxyz = np.zeros(3)
            if 'translate' in surface:
                dxyz = np.array(surface['translate'])

            if 'name' not in surface:
                print('no name...%s'  % surface)

            name = surface['name']
            if 'chord' not in surface:
                print('no chord for %s...' % name)
                continue

            if 'span' not in surface:
                npoint, nelement2 = get_fuselage(isurface, surface, xyz_scale, dxyz, nodes,
                                                 line_elements, quad_elements, surfaces, ipoint)
                if npoint == 0:
                    self.log.info('skipping %s because there are no sections' % surface)
                ipoint += npoint
                ielement += nelement2
                #break
                continue

            nchord, chord_spacing = surface['chord']
            nspan, span_spacing = surface['span']

            x = np.linspace(0., 1., num=nchord, endpoint=True, retstep=False, dtype=None)
            y = np.linspace(0., 1., num=nspan, endpoint=True, retstep=False, dtype=None)


            sections = surface['sections']
            del surface['sections']

            print('wing surface:', surface)

            #print(surface.keys())
            #print('name =', surface['name'])

            nsections = len(sections)
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
                #section = [chord, nspan, span_spacing]
                #section = [chord, ainc, nspan, span_spacing]
                chord0 = section0['section'][0]
                chord1 = section1['section'][0]
                #print('chords =', chord0, chord1)
                #print('xyz_scale =', xyz_scale)
                #incidence = section[1]
                p2 = p1 + np.array([chord0, 0., 0.])
                p3 = p4 + np.array([chord1, 0., 0.])

                point, element = points_elements_from_quad_points(p1, p2, p3, p4, x, y, dtype='int32')
                npoint = point.shape[0]
                nelement2 = element.shape[0]

                # scaling is applied before translation
                point[:, 0] *= xyz_scale[0]
                point[:, 1] *= xyz_scale[1]
                point[:, 2] *= xyz_scale[2]
                #print(point)
                point += dxyz
                surfaces.append(np.ones(nelement2, dtype='int32') * isurface)

                nodes.append(point)
                quad_elements.append(ipoint + element)
                #print(element)

                #for e in elements:
                    #print("  ", e)
                #print('npoint=%s nelement=%s' % (npoint, nelement2))
                ipoint += npoint
                ielement += nelement2
                #break
                #if not section['afile']:
                #del section['afile']
                #if not section['control']:
                #del section['control']
                #print('  ', section)
            #print('')
            #break
        nodes = np.vstack(nodes)
        quad_elements = np.vstack(quad_elements)
        if line_elements:
            line_elements = np.vstack(line_elements)
        surfaces = np.hstack(surfaces)
        return nodes, quad_elements, line_elements, surfaces

def get_fuselage(isurface, surface, xyz_scale, dxyz, nodes, line_elements, quad_elements, surfaces, ipoint):
    #print('----------------------------------------')
    #print(surface)
    nchord, chord_spacing = surface['chord']
    #nchord = 2
    x = np.linspace(0., 1., num=nchord, endpoint=True, retstep=False, dtype=None)
    y = np.array([0., 1.])
    sections = surface['sections']
    del surface['sections']

    #print(surface)
    #print(sections)
    nsections = len(sections)
    if nsections == 0:
        print(surface)
        body_file = surface['body_file']
        # top view/side view
        # x/c vs
        assert os.path.exists(body_file), body_file
        xc_vs_h = np.loadtxt(body_file, skiprows=1)
        xc = xc_vs_h[:, 0]
        h = xc_vs_h[:, 1]
        xc = xc[::-1]
        h = h[::-1]

        #x = h

        h0 = 0.
        # find x/c for h=0
        xc0 = np.interp(h0, h, xc, left=None, right=None, period=None)
        print('xc0 =', xc0)
        print('h =', h)
        print('xc =', xc)

        h2 = np.linspace(xc0, xc[0], num=nchord, endpoint=True, retstep=False, dtype=None)

        # find x/c for a range of h
        xc2 = np.interp(h2, h, xc)
        print('h2 =', h2)
        print('xc2 =', xc2)

        p1 = np.zeros(3) #dxyz
        xstation = xc2
        ystation = 0.
        zstation = 0.
        radii = h2
        aspect_ratio = 1.0
        point, element = create_axisymmetric_body(xstation, ystation, zstation, radii, aspect_ratio, p1)

        npoint = point.shape[0]
        nelement2 = element.shape[0]

        # scaling is applied before translation
        point[:, 0] *= xyz_scale[0]
        point[:, 1] *= xyz_scale[1]
        point[:, 2] *= xyz_scale[2]
        #print(point)
        print('xyz_scale ', xyz_scale)
        print('dxyz ', dxyz)
        point += dxyz

        surfaces.append(np.ones(nelement2, dtype='int32') * isurface)
        nodes.append(point)
        quad_elements.append(ipoint + element)
        #npoint = 0
        #nelement2 = 0
        return npoint, nelement2

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

        print('p1 = ', p1)
        print('p2 = ', p2)
        print('p3 = ', p3)
        print('p4 = ', p4)

        if 0:
            ipointi = np.arange(0, nchord, dtype='int32')
            dp = np.linspace(0., 1., num=nchord, endpoint=True, retstep=False, dtype=None)
            point = []
            for dpi in dp:
                #print(p1, type(p1), p2, dpi)
                pointi = p1*(1-dpi)# + p2*dpi
                point.append(pointi)
            point = np.array(point)
            #point = p1[np.newaxis, :] + p2[np.newaxis, :] * dp

            #point, element = points_elements_from_quad_points(p1, p2, p3, p4, x, y, dtype='int32')
            nelement2 = nchord - 1
            element = np.zeros((nelement2, 2), dtype='int32')
            element[:, 0] = ipointi[:-1]
            element[:, 1] = ipointi[1:]
            npoint = point.shape[0]
            #nelement2 = element.shape[0]

            # scaling is applied before translation
            point[:, 0] *= xyz_scale[0]
            point[:, 1] *= xyz_scale[1]
            point[:, 2] *= xyz_scale[2]
            point += dxyz
            print(point)
        else:
            npoint = point.shape[0]
            nelement2 = element.shape[0]

            # scaling is applied before translation
            point[:, 0] *= xyz_scale[0]
            point[:, 1] *= xyz_scale[1]
            point[:, 2] *= xyz_scale[2]
            #print(point)
            point += dxyz
            surfaces.append(np.ones(nelement2, dtype='int32') * isurface)

            nodes.append(point)
            quad_elements.append(ipoint + element)

        surfaces.append(np.ones(nelement2, dtype='int32') * isurface)

        nodes.append(point)
        quad_elements.append(ipoint + element)

        #print(element)

        #for e in elements:
            #print("  ", e)
        #print('npoint=%s nelement=%s' % (npoint, nelement2))
    return npoint, nelement2

def main():
    avl_filename = sys.argv[1]
    model = read_avl(avl_filename)
    model.get_nodes_elements()

if __name__ == '__main__':  # pragma: no cover
    main()
