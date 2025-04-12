from __future__ import annotations
import os
from typing import Optional, Any, TYPE_CHECKING

import numpy as np
from pyNastran.bdf.cards.aero.utils import (
    points_elements_from_quad_points, create_axisymmetric_body)
from pyNastran.converters.avl.avl_helper import get_spacing, save_wing_elements
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.converters.avl.surface import Surface


class Body:
    """
    BODY
    Fuse pod
    ! nchord schord
    12       1.0
    TRANSLATE
    ! x    y    z
    -12.5  0.0  -1.4
    BFIL
    fuseBD.dat
    """
    def __init__(self,
                 name: str,
                 sections: list[Any],
                 nchord: int, chord_spacing: float):

        self.name = name
        #self.component = component
        self.yduplicate = False
        #self.angle = angle
        self.body_file = None
        assert len(sections) == 0
        self.sections = sections
        self.nchord = nchord
        self.chord_spacing = chord_spacing
        #self.nspan = nspan
        #self.span_spacing =span_spacing
        #self.nowake = nowake
        self.scale = np.ones(3)
        self.translate = np.zeros(3)

    def __contains__(self, variable_name: str):
        if variable_name == 'name':
            return True
        if variable_name == 'scale':
            return self.has_scale()
        elif variable_name == 'translate':
            return self.has_translate()
        elif variable_name == 'body_file':
            return self.body_file is not None
        elif variable_name == 'yduplicate':
            return self.yduplicate is False
        raise NotImplementedError(variable_name)

    def write(self) -> str:
        surface_msg = (
            '#--------------------------------------------------\n'
            'BODY\n'
            f'{self.name}\n'
        )
        nchordwise = self.nchord
        c_space = self.chord_spacing
        surface_msg += ('!Nchordwise  Cspace  Nspanwise  Sspace\n'
                       f'{nchordwise}    {c_space}\n')

        if self.has_translate():
            dx, dy, dz = self.translate
            surface_msg += (
                '\n'
                'TRANSLATE\n'
                f'{dx}  {dy}  {dz}\n'
            )
        if self.has_scale():
            xscale, yscale, zscale = self.scale
            surface_msg += (
                '\n'
                'SCALE\n'
                f'{xscale}   {yscale}   {zscale}\n'
            )
        return surface_msg

    def has_scale(self) -> bool:
        return not np.array_equal(self.scale, np.ones(3))

    def has_translate(self) -> bool:
        return not np.array_equal(self.translate, np.zeros(3))

    def get_nodes_elements(self,
                           isurface: int, surfaces: list[Surface | Body],
                           dirname: str,
                           nodes: list[np.ndarray], ipoint: int,
                           line_elements: list[np.ndarray],
                           quad_elements: list[np.ndarray],
                           is_cs_list: list[np.ndarray],
                           log: SimpleLogger) -> int:
        log.debug('fuselage surface: %s\n' % str(self))
        if nodes:
            nodes_temp = np.vstack(nodes)
            assert nodes_temp.shape[0] == ipoint, 'nodes_temp.shape=%s ipoint=%s' % (nodes_temp.shape, ipoint)

        yduplicate = self.yduplicate
        ipoint = get_fuselage(dirname, isurface, self, yduplicate,
                              nodes, line_elements, quad_elements, surfaces, is_cs_list, ipoint)
        #if npoint == 0:
            #self.log.info('skipping %s because there are no sections' % surface)
        #ipoint += npoint
        #print("npoint = ", npoint)
        #print('-----------')
        nodes_temp = np.vstack(nodes)
        assert nodes_temp.shape[0] == ipoint, 'nodes_temp.shape=%s ipoint=%s' % (nodes_temp.shape, ipoint)
        #break
        return ipoint

    def __repr__(self) -> str:
        scale_translate = f', scale={self.scale}, translate={self.translate}'
        msg = (
            f'Body(name={self.name!r}, nchord={self.nchord}, chord_spacing={self.chord_spacing}{scale_translate})'
        )
        return msg

def get_fuselage(dirname: str, isurface: int,
                 surface: Body,
                 yduplicate: bool,
                 nodes: list[np.ndarray],
                 unused_line_elements: list[np.ndarray],
                 quad_elements: list[np.ndarray],
                 surfaces: list[Surface | Body],
                 is_cs_list: list[np.ndarray],
                 ipoint: int):
    """
    gets the fuselage

    TODO: doesn't support sine/cosine spacing on the fuselage
    """
    #print('----------------------------------------')
    #print(surface)
    nchord = surface.nchord
    chord_spacing = surface.chord_spacing

    assert nchord >= 1, nchord
    x = get_spacing(nchord, chord_spacing)
    y = np.array([0., 1.])
    #assert len(x) == len(y), 'x=%s y=%s' % (x, y)
    sections = surface.sections
    assert len(sections) == 0, surface

    #print(surface)
    #print(sections)
    xyz_scale = surface.scale
    dxyz = surface.translate
    nsections = len(sections)
    if nsections == 0:
        # from file
        ipoint, unused_nelement2 = get_fuselage_from_file(
            dirname, isurface, surface, nodes, quad_elements,
            surfaces, is_cs_list,
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
        p1 = section0['xyz_LE']
        p4 = section1['xyz_LE']
        chord0 = section0['section'][0]
        chord1 = section1['section'][0]
        #print('chords =', chord0, chord1)
        #print('xyz_scale =', xyz_scale)
        #incidence = section[1]
        p2 = p1 + np.array([chord0, 0., 0.])
        p3 = p4 + np.array([chord1, 0., 0.])

        point, element = points_elements_from_quad_points(p1, p2, p3, p4,
                                                          x, y, dtype='int32')
        nelements = element.shape[0]
        is_cs = np.zeros(nelements, dtype='int32')
        ipoint = save_wing_elements(
            isurface, point, element,
            xyz_scale, dxyz,
            nodes, quad_elements, surfaces,
            is_cs_list, is_cs,
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
                is_cs_list, is_cs,
                ipoint)
            nodes_temp = np.vstack(nodes)
            assert nodes_temp.shape[0] == ipoint, 'nodes_temp.shape=%s ipoint=%s' % (nodes_temp.shape, ipoint)

        #print('npoint=%s nelement=%s' % (npoint, nelement2))
    return ipoint

def get_fuselage_from_file(dirname: str,
                           isurface: int,
                           surface: Body,
                           nodes: list[np.ndarray],
                           quad_elements: list[np.ndarray],
                           surfaces: list[Surface | Body],
                           is_cs_list: list[np.ndarray],
                           ipoint: int,
                           nchord: int) -> tuple[int, int]:
    xyz_scale = surface.scale
    dxyz = surface.translate

    #print(surface)
    body_file = os.path.join(dirname, surface.body_file)
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
    point, element = create_axisymmetric_body(
        xstation, ystation, zstation, radii,
        aspect_ratio, p1)
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
    is_cs = np.zeros(nelement2, dtype='int32')
    assert len(is_cs) == nelement2
    is_cs_list.append(is_cs)

    nodes.append(point)
    quad_elements.append(ipoint + element)
    ipoint += npoint
    return ipoint, nelement2
