"""
defines NastranGuiAttributes, which defines
GUI specific geometry functions that don't involve PyQt/VTK
"""
# pylint: disable=E1101
from __future__ import print_function
from collections import defaultdict
import numpy as np
from numpy.linalg import norm

from pyNastran.utils import integer_types, iteritems


piover2 = np.pi / 2.
piover3 = np.pi / 3.


class NastranGuiAttributes(object):
    """GUI specific geometry functions that don't involve PyQt/VTK"""
    def __init__(self):
        # new options, no way to access them through the gui
        # they control results generation
        self.make_xyz = False
        self.make_offset_normals_dim = True
        self.make_nnodes_result = False  # make_offset_normals_dim must be True for this to work
        self.make_released_dofs1 = False
        self.make_released_dofs2 = False
        self.plot_applied_loads = True
        self.plot_pressures = True  # make_offset_normals_dim must be True for this to work

        # these are local variables
        self.model = None
        self.model_results = None
        self.bar_lines = None
        self.bar_eids = None
        self.xyz_cid0 = None
        self.normals = None
        self.save_data = True # was False
        #--------------------------------------

        #: flips the nastran CAERO subpaneling
        #:   False -> borders of CAEROs can be seen
        #:   True  -> individual subpanels can be seen
        self.show_caero_sub_panels = False

        #: coordinate systems can be messy, so this is the
        #: list of coords to show
        self.show_cids = []

        self.show_caero_actor = True  # show the caero mesh
        self.show_control_surfaces = True
        self.show_conm = True

        self.element_ids = None
        self.node_ids = None
        self.nid_map = None
        self.eid_map = None
        self.nNodes = None
        self.nElements = None
        self.model_type = None
        self.iSubcaseNameMap = None
        self.has_caero = False
        self.dependents_nodes = set([])
        self.i_transform = {}


class NastranGeometryHelper(NastranGuiAttributes):
    """
    Defines VTK/PyQt-less methods used by NastranIO-Geometry
    """
    def __init__(self):
        super(NastranGeometryHelper, self).__init__()

    def _get_bar_yz_arrays(self, model, bar_beam_eids, scale, debug):
        lines_bar_y = []
        lines_bar_z = []

        bar_types = {
            # PBAR
            'bar' : [],

            # PBEAML/PBARL
            "ROD": [],
            "TUBE": [],
            "TUBE2" : [],
            "I": [],
            "CHAN": [],
            "T": [],
            "BOX": [],
            "BAR": [],
            "CROSS": [],
            "H": [],
            "T1": [],
            "I1": [],
            "CHAN1": [],
            "Z": [],
            "CHAN2": [],
            "T2": [],
            "BOX1": [],
            "HEXA": [],
            "HAT": [],
            "HAT1": [],
            "DBOX": [],  # was 12

            # PBEAM
            'beam' : [],

            # PBEAML specfic
            "L" : [],
        }  # for GROUP="MSCBML0"
        allowed_types = [
            'BAR', 'BOX', 'BOX1', 'CHAN', 'CHAN1', 'CHAN2', 'CROSS', 'DBOX',
            'H', 'HAT', 'HAT1', 'HEXA', 'I', 'I1', 'L', 'ROD',
            'T', 'T1', 'T2', 'TUBE', 'TUBE2', 'Z', 'bar', 'beam',
        ]

        # bar_types['bar'] = [ [...], [...], [...] ]
        #bar_types = defaultdict(lambda : defaultdict(list))

        found_bar_types = set([])
        #neids = len(self.element_ids)
        for bar_type, data in iteritems(bar_types):
            eids = []
            lines_bar_y = []
            lines_bar_z = []
            bar_types[bar_type] = (eids, lines_bar_y, lines_bar_z)
            #bar_types[bar_type] = [eids, lines_bar_y, lines_bar_z]

        nid_release_map = defaultdict(list)

        #debug = True
        bar_nids = set([])
        #print('bar_beam_eids = %s' % bar_beam_eids)
        for eid in bar_beam_eids:
            if eid not in self.eid_map:
                self.log.error('eid=%s is not a valid bar/beam element...' % eid)
                if debug:
                    print('eid=%s is not a valid bar/beam element...' % eid)
                continue
            ieid = self.eid_map[eid]
            elem = model.elements[eid]
            pid = elem.pid
            assert not isinstance(pid, integer_types), elem
            if pid.type in ['PBAR', 'PBEAM']:
                bar_type = 'bar'
            elif pid.type in ['PBEAM']:
                bar_type = 'beam'
            elif pid.type in ['PBARL', 'PBEAML']:
                bar_type = pid.Type
            else:
                if debug:
                    print('NotImplementedError(pid)')
                raise NotImplementedError(pid)
            #print('bar_type =', bar_type)

            if debug:
                print('%s' % elem)
                print('  bar_type =', bar_type)
            found_bar_types.add(bar_type)

            (nid1, nid2) = elem.node_ids
            bar_nids.update([nid1, nid2])
            node1 = model.nodes[nid1]
            node2 = model.nodes[nid2]
            n1 = node1.get_position()
            n2 = node2.get_position()
            centroid = (n1 + n2) / 2.
            i = n2 - n1
            Li = norm(i)
            ihat = i / Li

            if elem.pa != 0:
                #assert elem.pa in [], elem.pa
                nid_release_map[nid1].append((eid, elem.pa))
            if elem.pb != 0:
                nid_release_map[nid2].append((eid, elem.pb))


            # OFFT flag
            # ---------
            # ABC or A-B-C (an example is G-G-G or B-G-G)
            # while the slots are:
            #  - A -> orientation; values=[G, B]
            #  - B -> End A; values=[G, O]
            #  - C -> End B; values=[G, O]
            #
            # and the values for A,B,C mean:
            #  - B -> basic
            #  - G -> global
            #  - O -> orientation
            #
            # so for example G-G-G, that's global for all terms.
            # BOG means basic orientation, orientation end A, global end B
            #
            # so now we're left with what does basic/global/orientation mean?
            # - basic -> the glboal coordinate system defined by cid=0
            # - global -> the local coordinate system defined by the
            #             CD field on the GRID card, but referenced by
            #             the CBAR/CBEAM
            # - orientation -> ???
            #
            if elem.g0:
                #debug = False
                msg = 'which is required by %s eid=%s\n%s' % (elem.type, elem.g0, str(elem))
                g0_ref = model.Node(elem.g0, msg=msg)
                if debug:
                    print('  g0 = %s' % elem.g0)
                    print('  g0_ref = %s' % g0_ref)
                n0 = g0_ref.get_position()
                v = n0 - n1
            else:
                #debug = False
                ga = model.nodes[elem.Ga()]
                cda = ga.Cd()
                cda_ref = model.Coord(cda)
                v = cda_ref.transform_node_to_global(elem.x)
                if debug:
                    print('  ga = %s' % elem.ga)
                    if cda != 0:
                        print('  cd = %s' % cda_ref)
                    else:
                        print('  cd = 0')

                    print('  x = %s' % elem.x)
                    print('  v = %s' % v)
                #v = elem.x

            offt_vector, offt_end_a, offt_end_b = elem.offt
            if debug:
                print('  offt vector,A,B=%r' % (elem.offt))
            # if offt_end_a == 'G' or (offt_end_a == 'O' and offt_vector == 'G'):

            cd1 = node1.Cd()
            cd2 = node2.Cd()
            cd1_ref = model.Coord(cd1)
            cd2_ref = model.Coord(cd2)
            # node1.cd_ref, node2.cd_ref

            if offt_vector == 'G':
                # end A
                # global - cid != 0
                if cd1 != 0:
                    v = cd1_ref.transform_node_to_global_assuming_rectangular(v)
                    #if node1.cd_ref.type not in ['CORD2R', 'CORD1R']:
                        #msg = 'invalid Cd type (%r) on Node %i; expected CORDxR' % (
                            #node1.cd_ref.type, node1.nid)
                        #self.log.error(msg)
                        #continue
                        #raise NotImplementedError(node1.cd)
            elif offt_vector == 'B':
                # basic - cid = 0
                pass
            else:
                msg = 'offt_vector=%r is not supported; offt=%s' % (offt_vector, elem.offt)
                self.log.error(msg)
                continue
                #raise NotImplementedError(msg)
            #print('v = %s' % v)

            # rotate wa
            wa = elem.wa
            if offt_end_a == 'G':
                if cd1 != 0:
                    #if node1.cd.type not in ['CORD2R', 'CORD1R']:
                        #continue # TODO: support CD transform
                    # TODO: fixme
                    wa = cd1_ref.transform_node_to_global_assuming_rectangular(wa)
            elif offt_end_a == 'B':
                pass
            elif offt_end_a == 'O':
                # TODO: fixme
                wa = cd1_ref.transform_node_to_global_assuming_rectangular(n1 - wa)
            else:
                msg = 'offt_end_a=%r is not supported; offt=%s' % (offt_end_a, elem.offt)
                self.log.error(msg)
                continue
                #raise NotImplementedError(msg)

            #print('wa = %s' % wa)
            # rotate wb
            wb = elem.wb
            if offt_end_b == 'G':
                if cd2 != 0:
                    #if cd2_ref.type not in ['CORD2R', 'CORD1R']:
                        #continue # TODO: MasterModelTaxi
                    # TODO: fixme
                    wb = cd2_ref.transform_node_to_global_assuming_rectangular(wb)

            elif offt_end_b == 'B':
                pass
            elif offt_end_b == 'O':
                # TODO: fixme
                wb = cd1_ref.transform_node_to_global(n2 - wb)
            else:
                msg = 'offt_end_b=%r is not supported; offt=%s' % (offt_end_b, elem.offt)
                model.log.error(msg)
                continue
                #raise NotImplementedError(msg)

            #print('wb =', wb)
            ## concept has a GOO
            #if not elem.offt in ['GGG', 'BGG']:
                #msg = 'offt=%r for CBAR/CBEAM eid=%s is not supported...skipping' % (
                    #elem.offt, eid)
                #self.log.error(msg)
                #continue

            vhat = v / norm(v) # j
            try:
                z = np.cross(ihat, vhat) # k
            except ValueError:
                msg = 'Invalid vector length\n'
                msg += 'n1  =%s\n' % str(n1)
                msg += 'n2  =%s\n' % str(n2)
                msg += 'nid1=%s\n' % str(nid1)
                msg += 'nid2=%s\n' % str(nid2)
                msg += 'i   =%s\n' % str(i)
                msg += 'Li  =%s\n' % str(Li)
                msg += 'ihat=%s\n' % str(ihat)
                msg += 'v   =%s\n' % str(v)
                msg += 'vhat=%s\n' % str(vhat)
                msg += 'z=cross(ihat, vhat)'
                print(msg)
                raise ValueError(msg)

            zhat = z / norm(z)
            yhat = np.cross(zhat, ihat) # j
            if debug:
                print('  centroid = %s' % centroid)
                print('  ihat = %s' % ihat)
                print('  yhat = %s' % yhat)
                print('  zhat = %s' % zhat)
                print('  scale = %s' % scale)
            #if eid == 5570:
                #print('  check - eid=%s yhat=%s zhat=%s v=%s i=%s n%s=%s n%s=%s' % (
                      #eid, yhat, zhat, v, i, nid1, n1, nid2, n2))

            if norm(ihat) == 0.0 or norm(yhat) == 0.0 or norm(z) == 0.0:
                print('  invalid_orientation - eid=%s yhat=%s zhat=%s v=%s i=%s n%s=%s n%s=%s' % (
                    eid, yhat, zhat, v, i, nid1, n1, nid2, n2))
            elif not np.allclose(norm(yhat), 1.0) or not np.allclose(norm(zhat), 1.0) or Li == 0.0:
                print('  length_error        - eid=%s Li=%s Lyhat=%s Lzhat=%s'
                      ' v=%s i=%s n%s=%s n%s=%s' % (
                          eid, Li, norm(yhat), norm(zhat), v, i, nid1, n1, nid2, n2))

            #print('adding bar %s' % bar_type)
            #print('   centroid=%s' % centroid)
            #print('   yhat=%s len=%s' % (yhat, np.linalg.norm(yhat)))
            #print('   zhat=%s len=%s' % (zhat, np.linalg.norm(zhat)))
            #print('   Li=%s scale=%s' % (Li, scale))
            if bar_type not in allowed_types:
                msg = 'bar_type=%r allowed=[%s]' % (bar_type, ', '.join(allowed_types))
                raise RuntimeError(msg)

            bar_types[bar_type][0].append(eid)
            bar_types[bar_type][1].append((centroid, centroid + yhat * Li * scale))
            bar_types[bar_type][2].append((centroid, centroid + zhat * Li * scale))
            #if len(bar_types[bar_type][0]) > 5:
                #break

        #print('bar_types =', bar_types)
        for bar_type in list(bar_types):
            bars = bar_types[bar_type]
            if len(bars[0]) == 0:
                del bar_types[bar_type]
                continue
            #bar_types[bar_type][1] = np.array(bars[1], dtype='float32')  # lines_bar_y
            #bar_types[bar_type][2] = np.array(bars[2], dtype='float32')  # lines_bar_z

        debug = False
        if debug:
            #np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
            for bar_type, data in sorted(iteritems(bar_types)):
                eids, lines_bar_y, lines_bar_z = data
                if len(eids):
                    #print('barsi =', barsi)
                    #print('bar_type = %r' % bar_type)
                    for eid, line_y, line_z  in zip(eids, lines_bar_y, lines_bar_z):
                        print('eid=%s centroid=%s cy=%s cz=%s' % (
                            eid, line_y[0], line_y[1], line_z[1]))

        #print('found_bar_types =', found_bar_types)
        #no_axial_torsion = (no_axial, no_torsion)
        #no_shear_bending = (no_shear_y, no_shear_z, no_bending_y, no_bending_z)
        #no_dofs = (no_bending, no_bending_bad, no_6_16, no_0_456,
                   #no_0_56, no_56_456, no_0_6, no_0_16)
        return bar_nids, bar_types, nid_release_map

    def _get_suport_node_ids(self, model, suport_id):
        """gets the nodes where SUPORTs and SUPORT1s are defined"""
        node_ids = []
        # list
        #for suport in model.suport:
            #node_ids += suport.IDs

        # dict
        if suport_id in model.suport1:
            suport1 = model.suport1[suport_id]
            node_ids += suport1.IDs
        else:
            for suport in model.suport:  # TODO: shouldn't this be included?
                if suport_id in suport.IDs:
                    node_ids.append(suport_id)
        return np.unique(node_ids)

    def _get_material_arrays(self, model, mids):
        """gets e11, e22, e33"""
        e11 = np.zeros(mids.shape, dtype='float32')
        e22 = np.zeros(mids.shape, dtype='float32')
        e33 = np.zeros(mids.shape, dtype='float32')
        rho = np.zeros(mids.shape, dtype='float32')
        bulk = np.zeros(mids.shape, dtype='float32')
        speed_of_sound = np.zeros(mids.shape, dtype='float32')

        has_mat8 = False
        has_mat9 = False
        #has_mat10 = False
        for umid in np.unique(mids):
            if umid == 0:
                continue
            e11i = e22i = e33i = 0.
            rhoi = 0.
            bulki = 0.
            speed_of_soundi = 0.
            try:
                mat = model.materials[umid]
            except KeyError:
                print('mids = %s' % mids)
                print('mids = %s' % model.materials.keys())
                continue
                #raise
            if mat.type == 'MAT1':
                e11i = e22i = e33i = mat.e
                rhoi = mat.rho
            elif mat.type == 'MAT8':
                e11i = e33i = mat.e11
                e22i = mat.e22
                has_mat8 = True
                rhoi = mat.rho
            elif mat.type in ['MAT11', 'MAT3D']:
                e11i = mat.e1
                e22i = mat.e2
                e33i = mat.e3
                has_mat9 = True
                rhoi = mat.rho
            elif mat.type == 'MAT10':
                bulki = mat.bulk
                rhoi = mat.rho
                speed_of_soundi = mat.c
                #has_mat10 = True
                #self.log.info('skipping\n%s' % mat)
                #continue
            else:
                print('skipping\n%s' % mat)
                continue
                #raise NotImplementedError(mat)
            #print('mid=%s e11=%e e22=%e' % (umid, e11i, e22i))
            i = np.where(umid == mids)[0]
            e11[i] = e11i
            e22[i] = e22i
            e33[i] = e33i
            rho[i] = rhoi
            bulk[i] = bulki
            speed_of_sound[i] = speed_of_soundi
        return has_mat8, has_mat9, e11, e22, e33

def tri_quality(p1, p2, p3):
    """gets the quality metrics for a tri"""
    e1 = (p1 + p2) / 2.
    e2 = (p2 + p3) / 2.
    e3 = (p3 + p1) / 2.

    #    3
    #    / \
    # e3/   \ e2
    #  /    /\
    # /    /  \
    # 1---/----2
    #    e1
    e21 = e2 - e1
    e31 = e3 - e1
    e32 = e3 - e2

    e3_p2 = e3 - p2
    e2_p1 = e2 - p1
    e1_p3 = e1 - p3

    v21 = p2 - p1
    v32 = p3 - p2
    v13 = p1 - p3
    length21 = np.linalg.norm(v21)
    length32 = np.linalg.norm(v32)
    length13 = np.linalg.norm(v13)
    min_edge_length = min(length21, length32, length13)
    areai = 0.5 * np.linalg.norm(np.cross(v21, v13))

    cos_skew1 = np.dot(e2_p1, e31) / (np.linalg.norm(e2_p1) * np.linalg.norm(e31))
    cos_skew2 = np.dot(e2_p1, -e31) / (np.linalg.norm(e2_p1) * np.linalg.norm(e31))
    cos_skew3 = np.dot(e3_p2, e21) / (np.linalg.norm(e3_p2) * np.linalg.norm(e21))
    cos_skew4 = np.dot(e3_p2, -e21) / (np.linalg.norm(e3_p2) * np.linalg.norm(e21))
    cos_skew5 = np.dot(e1_p3, e32) / (np.linalg.norm(e1_p3) * np.linalg.norm(e32))
    cos_skew6 = np.dot(e1_p3, -e32) / (np.linalg.norm(e1_p3) * np.linalg.norm(e32))
    max_skew = np.pi / 2. - np.abs(np.arccos(np.clip([
        cos_skew1, cos_skew2, cos_skew3,
        cos_skew4, cos_skew5, cos_skew6], -1., 1.))).min()
    lengths = np.linalg.norm([v21, v32, v13], axis=1)
    #assert len(lengths) == 3, lengths
    aspect_ratio = lengths.max() / lengths.min()

    cos_theta1 = np.dot(v21, -v13) / (length21 * length13)
    cos_theta2 = np.dot(v32, -v21) / (length32 * length21)
    cos_theta3 = np.dot(v13, -v32) / (length13 * length32)
    thetas = np.arccos(np.clip([cos_theta1, cos_theta2, cos_theta3], -1., 1.))
    min_thetai = thetas.min()
    max_thetai = thetas.max()
    dideal_thetai = max(max_thetai - piover3, piover3 - min_thetai)

    #theta_deg = np.degrees(np.arccos(max_cos_theta))
    #if theta_deg < 60.:
        #print('p1=%s' % xyz_cid0[p1, :])
        #print('p2=%s' % xyz_cid0[p2, :])
        #print('p3=%s' % xyz_cid0[p3, :])
        #print('theta1=%s' % np.degrees(np.arccos(cos_theta1)))
        #print('theta2=%s' % np.degrees(np.arccos(cos_theta2)))
        #print('theta3=%s' % np.degrees(np.arccos(cos_theta3)))
        #print('max_theta=%s' % theta_deg)
        #asdf
    return areai, max_skew, aspect_ratio, min_thetai, max_thetai, dideal_thetai, min_edge_length


def quad_quality(p1, p2, p3, p4):
    """gets the quality metrics for a quad"""
    v21 = p2 - p1
    v32 = p3 - p2
    v43 = p4 - p3
    v14 = p1 - p4
    length21 = np.linalg.norm(v21)
    length32 = np.linalg.norm(v32)
    length43 = np.linalg.norm(v43)
    length14 = np.linalg.norm(v14)
    min_edge_length = min(length21, length32, length43, length14)

    v42 = p4 - p2
    v31 = p3 - p1
    p12 = (p1 + p2) / 2.
    p23 = (p2 + p3) / 2.
    p34 = (p3 + p4) / 2.
    p14 = (p4 + p1) / 2.
    v31 = p3 - p1
    v42 = p4 - p2
    normal = np.cross(v31, v42)
    areai = 0.5 * np.linalg.norm(normal)

    # still kind of in development
    #
    # the ratio of the ideal area to the actual area
    # this is an hourglass check
    areas = [
        np.linalg.norm(np.cross(-v14, v21)), # v41 x v21
        np.linalg.norm(np.cross(v32, -v21)), # v32 x v12
        np.linalg.norm(np.cross(v43, -v32)), # v43 x v23
        np.linalg.norm(np.cross(v14, v43)),  # v14 x v43
    ]
    #
    # for:
    #   area=1; area1=0.5 -> area_ratioi1=2.0; area_ratio=2.0
    #   area=1; area1=2.0 -> area_ratioi2=2.0; area_ratio=2.0
    area_ratioi1 = areai / min(areas)
    area_ratioi2 = max(areas) / areai
    area_ratioi = max(area_ratioi1, area_ratioi2)

    area1 = 0.5 * np.linalg.norm(np.cross(-v14, v21)) # v41 x v21
    area2 = 0.5 * np.linalg.norm(np.cross(-v21, v32)) # v12 x v32
    area3 = 0.5 * np.linalg.norm(np.cross(v43, v32)) # v43 x v32
    area4 = 0.5 * np.linalg.norm(np.cross(v14, -v43)) # v14 x v34
    aavg = (area1 + area2 + area3 + area4) / 4.
    taper_ratioi = (abs(area1 - aavg) + abs(area2 - aavg) +
                    abs(area3 - aavg) + abs(area4 - aavg)) / aavg

    #    e3
    # 4-------3
    # |       |
    # |e4     |  e2
    # 1-------2
    #     e1
    e13 = p34 - p12
    e42 = p23 - p14
    cos_skew1 = np.dot(e13, e42) / (np.linalg.norm(e13) * np.linalg.norm(e42))
    cos_skew2 = np.dot(e13, -e42) / (np.linalg.norm(e13) * np.linalg.norm(e42))
    max_skew = np.pi / 2. - np.abs(np.arccos(
        np.clip([cos_skew1, cos_skew2], -1., 1.))).min()
    #aspect_ratio = max(p12, p23, p34, p14) / max(p12, p23, p34, p14)
    lengths = np.linalg.norm([v21, v32, v43, v14], axis=1)
    #assert len(lengths) == 3, lengths
    aspect_ratio = lengths.max() / lengths.min()

    cos_theta1 = np.dot(v21, -v14) / (length21 * length14)
    cos_theta2 = np.dot(v32, -v21) / (length32 * length21)
    cos_theta3 = np.dot(v43, -v32) / (length43 * length32)
    cos_theta4 = np.dot(v14, -v43) / (length14 * length43)
    #max_thetai = np.arccos([cos_theta1, cos_theta2, cos_theta3, cos_theta4]).max()

    # dot the local normal with the normal vector
    # then take the norm of that to determine the angle relative to the normal
    # then take the sign of that to see if we're pointing roughly towards the normal

    # np.sign(np.linalg.norm(np.dot(
    # a x b = ab sin(theta)
    # a x b / ab = sin(theta)
    # sin(theta) < 0. -> normal is flipped
    normal2 = np.sign(np.dot(np.cross(v21, v32), normal))
    normal3 = np.sign(np.dot(np.cross(v32, v43), normal))
    normal4 = np.sign(np.dot(np.cross(v43, v14), normal))
    normal1 = np.sign(np.dot(np.cross(v14, v21), normal))
    n = np.array([normal1, normal2, normal3, normal4])
    theta_additional = np.where(n < 0, 2*np.pi, 0.)

    theta = n * np.arccos(np.clip(
        [cos_theta1, cos_theta2, cos_theta3, cos_theta4], -1., 1.)) + theta_additional
    min_thetai = theta.min()
    max_thetai = theta.max()
    dideal_thetai = max(max_thetai - piover2, piover2 - min_thetai)
    #print('theta_max = ', theta_max)

    #if 0:
        ## warp
        #v31 = xyz_cid0[p3, :] - xyz_cid0[p1, :]
        #n1a = np.cross(v21, v31) # v21 x v31
        #n1b = np.cross(v31, -v14) # v31 x v41
        #warp1 = np.dot(n1a, n1b) / (np.linalg.norm(n1a) * np.linalg.norm(n1b))

        #v42 = xyz_cid0[p4, :] - xyz_cid0[p2, :]
        #n2a = np.cross(v32, v42) # v32 x v42
        #n2b = np.cross(v42, -v21) # v42 x v12
        #warp2 = np.dot(n2a, n2b) / (np.linalg.norm(n2a) * np.linalg.norm(n2b))
        #max_warp = max(np.arccos(warp1), np.arccos(warp2))
    out = (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
           min_thetai, max_thetai, dideal_thetai, min_edge_length)
    return out

def get_min_max_theta(faces, all_node_ids, nid_map, xyz_cid0):
    """get the min/max thetas for CTETRA, CPENTA, CHEXA, CPYRAM"""
    cos_thetas = []
    ideal_theta = []
    #print('faces =', faces)
    #assert len(faces) > 0, 'faces=%s nids=%s' % (faces, all_node_ids)
    for face in faces:
        if len(face) == 3:
            node_ids = all_node_ids[face[0]], all_node_ids[face[1]], all_node_ids[face[2]]
            n1, n2, n3 = [nid_map[nid] for nid in node_ids[:3]]
            v21 = xyz_cid0[n2, :] - xyz_cid0[n1, :]
            v32 = xyz_cid0[n3, :] - xyz_cid0[n2, :]
            v13 = xyz_cid0[n1, :] - xyz_cid0[n3, :]
            length21 = np.linalg.norm(v21)
            length32 = np.linalg.norm(v32)
            length13 = np.linalg.norm(v13)
            min_edge_length = min(length21, length32, length13)

            cos_theta1 = np.dot(v21, -v13) / (length21 * length13)
            cos_theta2 = np.dot(v32, -v21) / (length32 * length21)
            cos_theta3 = np.dot(v13, -v32) / (length13 * length32)
            cos_thetas.extend([cos_theta1, cos_theta2, cos_theta3])
            ideal_theta.extend([piover3, piover3, piover3])
        elif len(face) == 4:
            try:
                node_ids = (all_node_ids[face[0]], all_node_ids[face[1]],
                            all_node_ids[face[2]], all_node_ids[face[3]])
            except:
                print(face)
                print(node_ids)
                raise

            n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
            v21 = xyz_cid0[n2, :] - xyz_cid0[n1, :]
            v32 = xyz_cid0[n3, :] - xyz_cid0[n2, :]
            v43 = xyz_cid0[n4, :] - xyz_cid0[n3, :]
            v14 = xyz_cid0[n1, :] - xyz_cid0[n4, :]
            length21 = np.linalg.norm(v21)
            length32 = np.linalg.norm(v32)
            length43 = np.linalg.norm(v43)
            length14 = np.linalg.norm(v14)
            min_edge_length = min(length21, length32, length43, length14)
            cos_theta1 = np.dot(v21, -v14) / (length21 * length14)
            cos_theta2 = np.dot(v32, -v21) / (length32 * length21)
            cos_theta3 = np.dot(v43, -v32) / (length43 * length32)
            cos_theta4 = np.dot(v14, -v43) / (length14 * length43)
            cos_thetas.extend([cos_theta1, cos_theta2, cos_theta3, cos_theta4])
            ideal_theta.extend([piover2, piover2, piover2, piover2])
        else:
            raise NotImplementedError(face)
    thetas = np.arccos(cos_thetas)
    ideal_theta = np.array(ideal_theta)
    ideal_thetai = max((thetas - ideal_theta).max(), (ideal_theta - thetas).min())

    min_thetai = thetas.min()
    max_thetai = thetas.max()
    return min_thetai, max_thetai, ideal_thetai, min_edge_length
