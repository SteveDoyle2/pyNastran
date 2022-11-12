import os
import sys
from typing import Union

import numpy as np
from cpylog import get_logger2
from pyNastran.converters.avl.control import Control
from pyNastran.converters.avl.surface import Surface
from pyNastran.converters.avl.body import Body

AVL_KEYWORDS_LONG = [
    'SURFACE', 'COMPONENT', 'YDUPLICATE', 'BODY', 'ANGLE', 'BFILE',
    'NOWAKE', 'NOLOAD',
    'NACA', 'SCALE', 'TRANSLATE', 'SECTION', 'AFILE', 'CONTROL',
    # unsupported
    #'NOLABE', 'DESIGN', 'AIRFOIL', 'CLAF', 'CDCL',
]

AVL_KEYWORDS = [_avl_keyword[:4] for _avl_keyword in AVL_KEYWORDS_LONG]

def read_avl(avl_filename, log=None, debug: Union[str, bool, None]=False):
    """reads a *.avl file"""
    avl = AVL(log=log, debug=debug)
    avl.read_avl(avl_filename)
    return avl


class Section:
    def __init__(self, xle, yle, zle, chord, ainc, nspan, span_spacing):
        self.xyz_le = np.array([xle, yle, zle])
        self.chord = chord
        self.ainc = ainc
        self.nspan = nspan
        self.span_spacing = span_spacing
        self.is_afile = None
        self.afile = ''
        self.naca = ''
        self.controls = []

    def __repr__(self) -> str:
        ainc = ''
        if self.ainc != 0.0:
            ainc = f', ainc={self.ainc}'
        msg = (
            f'Section(xyz_le={self.zyz_le}, chord={self.chord}, nspan={self.nspan}, span_spacing={self.span_spacing}{ainc})'
        )
        return msg

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
        self.cd0 = None
        self.sections = []
        self.surfaces = None
        self.sym_iy = None
        self.sym_iz = None
        self.symz = None

    def read_avl(self, avl_filename: str) -> None:
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
        sref, cref, bref = sline3[:3]

        sline4 = lines[4].split()
        xref, yref, zref = sline4[:3]

        line5 = lines[5].strip()
        self.sym_iy = int(sym_iy)
        self.sym_iz = int(sym_iz)

        self.symz = float(symz)

        self.sref = float(sref)
        self.cref = float(cref)
        self.bref = float(bref)

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
                elif len(sline2) == 2:
                    nchord, chord_spacing = sline2
                    self.log.debug('name=%s nchord=%s chord_spacing=%s' % (
                        name, nchord, chord_spacing))
                    #surface['span'] = (None, None)
                else:
                    raise NotImplementedError(sline2)
                nchord = int(nchord)
                chord_spacing = float(chord_spacing)
                surface['chord'] = (nchord, chord_spacing)
                surface['span'] = (nspan, span_spacing)
                surface['sections'] = sections
                surf = Surface(
                    name,
                    sections,
                    nchord, chord_spacing,
                    nspan, span_spacing)
                surface = surf
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
                surf.component = ncomponents
                i += 1

            elif is_header(line, 'YDUPLICATE'):
                yduplicate = float(lines[i+1])
                assert 'yduplicate' not in surface
                surf.yduplicate = yduplicate
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
                surface = Body(
                    name,
                    sections,
                    nchord, chord_spacing)
                surf = surface
                surfaces.append(surface)

                i += 2

            elif is_header(line, 'ANGLE'):
                angle = float(lines[i+1])
                assert 'angle' not in surface
                surf.angle = angle
                i += 1
            elif is_header(line, 'BFILE'):
                body_file = lines[i+1]
                assert 'body_file' not in surface
                surf.body_file = body_file
                i += 1

            elif is_header(line, 'NOWAKE'):
                assert 'nowake' not in surface
                surf.nowake = True
            elif is_header(line, 'NOLOAD'):
                assert 'noload' not in surface
                surf.noload = True

            elif is_header(line, 'SCALE'):
                xscale, yscale, zscale = lines[i+1].split()
                xscale = float(xscale)
                yscale = float(yscale)
                zscale = float(zscale)
                xyz_scale = np.array([xscale, yscale, zscale])
                assert 'scale' not in surface
                surf.scale = xyz_scale
                i += 1

            elif is_header(line, 'TRANSLATE'):
                xtranslate, ytranslate, ztranslate = lines[i+1].split()
                xtranslate = float(xtranslate)
                ytranslate = float(ytranslate)
                ztranslate = float(ztranslate)
                xyz_translate = np.array([xtranslate, ytranslate, ztranslate])
                assert 'translate' not in surface
                if isinstance(surface, (Body, Surface)):
                    assert np.array_equal(surface.translate, np.zeros(3))
                surf.translate = xyz_translate
                i += 1
            elif is_header(line, 'SECTION'):
                section, section_data, sectioni = _read_section(surface, lines, i)
                section_data['section'] = section
                sections.append(section_data)
                i += 1
            elif line == 'NACA':
                section_data['is_afile'] = False
                section_data['naca'] = lines[i+1]
                sectioni.is_afile = False
                sectioni.naca = lines[i+1]
                i += 1
            elif is_header(line, 'AFILE'):
                section_data['is_afile'] = True
                section_data['afile'] = lines[i+1]
                sectioni.is_afile = True
                sectioni.afile = lines[i+1]
                i += 1
            elif is_header(line, 'CONTROL'):
                sline = lines[i+1].split()
                #print(sline)
                try:
                    name, gain_str, xhinge_str, xhinge_vector, yhinge_vector, zhinge_vector, sign_deflection_duplicated_str = sline
                except Exception:
                    self.log.error(sline)
                    raise
                gain = float(gain_str)
                xhinge = float(xhinge_str)
                sign_deflection_duplicated = float(sign_deflection_duplicated_str)
                hinge_vector = np.array([xhinge_vector, yhinge_vector, zhinge_vector], dtype='float64')
                control = Control(name, gain, xhinge, hinge_vector, sign_deflection_duplicated)
                section_data['control'].append(control)
                i += 1
            else:
                print(line)
            i += 1

        for surface in surfaces:
            assert isinstance(surface, (Body, Surface)), surface
            name = surface.name
            sections = surface.sections
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


    def write_avl(self, avl_filename: str) -> None:
        """very preliminary and full of del commands that corrupt the model"""
        msg = (
            f'{self.name}\n'
            '\n'
            '#Mach\n'
            f' {self.mach}\n'
            '\n'
            '#IYsym   IZsym   Zsym\n'
            f' {self.sym_iy}       {self.sym_iz}       {self.symz}\n'
            '\n'
            '#Sref    Cref    Bref\n'
            f'{self.sref}   {self.cref}    {self.bref}\n'
            '\n'
            '#Xref    Yref    Zref\n'
            f'{self.xref}     {self.yref}     {self.zref}\n'
        )
        if self.cd0 is not None:
            msg += (
                '! CDo\n'
                f'{self.cd0}\n'
            )


        for isurface, surface in enumerate(self.surfaces):
            msg += surface.write()

        with open(avl_filename, 'w') as avl_file:
            avl_file.write(msg)

    def get_nodes_elements(self):
        log = self.log
        log.debug('get_nodes_elements')
        dirname = os.path.dirname(self.avl_filename)
        nodes = []
        quad_elements = []
        line_elements = []
        ipoint = 0
        surfaces = []
        is_cs_list = []

        #print('----')
        for isurface, surface in enumerate(self.surfaces):

            log.debug('isurface = %s' % isurface)
            if isinstance(surface, Body):
                ipoint = surface.get_nodes_elements(
                    isurface, surfaces,
                    dirname,
                    nodes, ipoint,
                    line_elements, quad_elements, is_cs_list,
                    log)
                continue
            elif isinstance(surface, Surface):
                ipoint = surface.get_nodes_elements(
                    isurface, surfaces,
                    dirname,
                    nodes, ipoint,
                    line_elements, quad_elements, is_cs_list,
                    log)
                continue
            raise RuntimeError(surface)

        #print("end ipoint=%s" % (ipoint))
        nodes = np.vstack(nodes)
        quad_elements = np.vstack(quad_elements)
        if line_elements:
            line_elements = np.vstack(line_elements)
        surfaces = np.hstack(surfaces)
        is_cs = np.hstack(is_cs_list)
        assert surfaces.shape == is_cs.shape
        assert len(surfaces) == quad_elements.shape[0]
        return nodes, quad_elements, line_elements, surfaces, is_cs


def get_lofted_sections(sections) -> None:
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

def is_header(line: str, name: str) -> bool:
    """only the first 4 chancters are read, but we're going to ensure all the letters are correct"""
    return line.startswith(name[:4]) and line == name[:len(line)]


def _read_section(surface, lines: list[str], i: int):
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
        nspan = None
        span_spacing = None
        section = [chord, ainc, None, None]
    elif len(sline) == 7:
        #Xle Yle  Zle  Chord  Ainc  Nspanwise  Sspace
        xle, yle, zle, chord, ainc, nspan, span_spacing = sline
        xle = float(xle)
        yle = float(yle)
        zle = float(zle)
        ainc = float(ainc)
        chord = float(chord)
        nspan = int(nspan)
        span_spacing = float(span_spacing)
        if isinstance(surface, dict):
            assert 'station' not in surface
        section = [chord, ainc, nspan, span_spacing]
    else:  # pragma: no cover
        #print(sline, len(sline))
        raise NotImplementedError(sline)

    section_data['xyz_LE'] = [xle, yle, zle]
    sectioni = Section(xle, yle, zle,
                       chord, ainc, nspan, span_spacing)
    return section, section_data, sectioni

def main():  # pragma: no cover
    avl_filename = sys.argv[1]
    model = read_avl(avl_filename)
    model.get_nodes_elements()

    import os
    base, ext = os.path.splitext(avl_filename)
    avl_filename_out = base + '_out' + ext
    model.write_avl(avl_filename_out)

if __name__ == '__main__':  # pragma: no cover
    main()
