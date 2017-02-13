# pylint: disable=E1101
from __future__ import print_function
import numpy as np
from numpy.linalg import norm

from pyNastran.utils import integer_types, iteritems
from pyNastran.bdf.cards.loads.static_loads import LOAD


piover2 = np.pi / 2.
piover3 = np.pi / 3.


class NastranGuiAttributes(object):
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
        #method = 0 # TODO: this should be reworked...
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

        no_axial = np.zeros(self.element_ids.shape, dtype='int32')
        no_torsion = np.zeros(self.element_ids.shape, dtype='int32')

        if self.make_released_dofs1:
            no_shear_y = np.zeros(self.element_ids.shape, dtype='int32')
            no_shear_z = np.zeros(self.element_ids.shape, dtype='int32')
            no_bending_y = np.zeros(self.element_ids.shape, dtype='int32')
            no_bending_z = np.zeros(self.element_ids.shape, dtype='int32')

        if self.make_released_dofs2:
            no_bending = np.zeros(self.element_ids.shape, dtype='int32')
            no_bending_bad = np.zeros(self.element_ids.shape, dtype='int32')
            no_6_16 = np.zeros(self.element_ids.shape, dtype='int32')
            no_0_56 = np.zeros(self.element_ids.shape, dtype='int32')
            no_0_456 = np.zeros(self.element_ids.shape, dtype='int32')
            no_56_456 = np.zeros(self.element_ids.shape, dtype='int32')
            no_0_6 = np.zeros(self.element_ids.shape, dtype='int32')
            no_0_16 = np.zeros(self.element_ids.shape, dtype='int32')

        #debug = True
        bar_nids = set([])
        print('bar_beam_eids = %s' % bar_beam_eids)
        for eid in bar_beam_eids:
            if eid not in self.eid_map:
                self.log.error('eid=%s is not a valid element...' % eid)
                if debug:
                    print('eid=%s is not a valid element...' % eid)
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

            if not self.make_released_dofs2:
                if elem.pa == 0 and elem.pb == 0:
                    pass

                if elem.pa == 1 or elem.pb == 1:
                    no_axial[ieid] = 1
                if elem.pa == 2 or elem.pb == 2:
                    no_axial[ieid] = 1
                if elem.pa == 3 or elem.pb == 3:
                    no_axial[ieid] = 1
                if elem.pa == 4 or elem.pb == 4:
                    no_torsion[ieid] = 1
                if elem.pa == 5 or elem.pb == 5:
                    no_axial[ieid] = 1
                if elem.pa == 6 or elem.pb == 6:
                    no_axial[ieid] = 1

            else:
                if elem.pa == 0 and elem.pb == 0:
                    pass
                elif (elem.pa == 6 and elem.pb == 16) or (elem.pa == 16 and elem.pb == 6):
                    no_axial[ieid] = 1
                    no_6_16[ieid] = 1
                elif (elem.pa == 56 and elem.pb == 0) or (elem.pa == 0 and elem.pb == 56):
                    no_bending[ieid] = 1
                    no_0_56[ieid] = 1
                    #print(elem)
                elif (elem.pa == 0 and elem.pb == 456) or (elem.pa == 456 and elem.pb == 0):
                    no_bending[ieid] = 1
                    no_torsion[ieid] = 1
                    no_0_456[ieid] = 1
                    #print(elem)
                elif (elem.pa == 456 and elem.pb == 56) or (elem.pa == 56 and elem.pb == 456):
                    no_torsion[ieid] = 1
                    no_56_456[ieid] = 1
                elif elem.pa == 6 and elem.pb == 0:
                    no_bending_bad[ieid] = 1
                    no_0_6[ieid] = 1
                    #print(elem)
                elif elem.pa == 0 and elem.pb == 16 or elem.pb == 0 and elem.pa == 16:
                    no_axial[ieid] = 1
                    no_bending_bad[ieid] = 1
                    # print(elem)
                    no_0_16[ieid] = 1
                elif elem.pa == 56 and elem.pb == 45 or elem.pb == 56 and elem.pa == 45:
                    no_torsion[ieid] = 1
                    no_bending[ieid] = 1
                else:
                    msg = 'pa=%r pb=%r; elem=\n%s' % (elem.pa, elem.pb, elem)
                    print(msg)
                    continue
                    #raise NotImplementedError(msg)


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
                v = ga.cd_ref.transform_node_to_global(elem.x)
                if debug:
                    print('  ga = %s' % elem.ga)
                    if ga.Cd() != 0:
                        print('  cd = %s' % ga.cd_ref)
                    else:
                        print('  cd = 0')

                    print('  x = %s' % elem.x)
                    print('  v = %s' % v)
                #v = elem.x

            offt_vector, offt_end_a, offt_end_b = elem.offt
            if debug:
                print('  offt vector,A,B=%r' % (elem.offt))
            # if offt_end_a == 'G' or (offt_end_a == 'O' and offt_vector == 'G'):

            if offt_vector == 'G':
                # end A
                # global - cid != 0
                if node1.Cd() != 0:
                    v = node1.cd_ref.transform_node_to_global_assuming_rectangular(v)
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
            #print('v =', v)

            # rotate wa
            wa = elem.wa
            if offt_end_a == 'G':
                if node1.Cd() != 0:
                    #if node1.cd.type not in ['CORD2R', 'CORD1R']:
                        #continue # TODO: support CD transform
                    # TODO: fixme
                    wa = node1.cd_ref.transform_node_to_global_assuming_rectangular(wa)
            elif offt_end_a == 'B':
                pass
            elif offt_end_a == 'O':
                # TODO: fixme
                wa = node1.cd_ref.transform_node_to_global_assuming_rectangular(n1 - wa)
            else:
                msg = 'offt_end_a=%r is not supported; offt=%s' % (offt_end_a, elem.offt)
                self.log.error(msg)
                continue
                #raise NotImplementedError(msg)

            #print('wa =', wa)
            # rotate wb
            wb = elem.wb
            if offt_end_b == 'G':
                if node2.Cd() != 0:
                    #if node2.cd_ref.type not in ['CORD2R', 'CORD1R']:
                        #continue # TODO: MasterModelTaxi
                    wb = node2.cd.transform_node_to_global_assuming_rectangular(wb)  # TODO: fixme

            elif offt_end_b == 'B':
                pass
            elif offt_end_b == 'O':
                wb = node1.cd.transform_node_to_global(n2 - wb)  # TODO: fixme
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
            if len(bar_types[bar_type][0]) == 0:
                del bar_types[bar_type]

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
        #no_axial, no_torsion
        out = None
        return bar_nids, bar_types, out

    def _get_rigid(self, model):
        """
        dependent = (lines[:, 0])
        independent = np.unique(lines[:, 1])
        """
        lines_rigid = []
        for eid, elem in iteritems(model.rigid_elements):
            if elem.type == 'RBE3':
                if elem.Gmi != []:
                    msg = 'UM is not supported; RBE3 eid=%s Gmi=%s' % (elem.eid, elem.Gmi)
                    raise RuntimeError(msg)
                #list_fields = ['RBE3', elem.eid, None, elem.ref_grid_id, elem.refc]
                n1 = elem.ref_grid_id
                assert isinstance(n1, int), 'RBE3 eid=%s ref_grid_id=%s' % (elem.eid, n1)
                for (_weight, ci, Gij) in elem.WtCG_groups:
                    Giji = elem._nodeIDs(nodes=Gij, allow_empty_nodes=True)
                    # list_fields += [wt, ci] + Giji
                    for n2 in Giji:
                        assert isinstance(n2, int), 'RBE3 eid=%s Giji=%s' % (elem.eid, Giji)
                        lines_rigid.append([n1, n2])
            elif elem.type == 'RBE2':
                #list_fields = ['RBE2', elem.eid, elem.Gn(), elem.cm
                               #] + elem.Gmi_node_ids + [elem.alpha]
                n2 = elem.Gn() # independent
                nids1 = elem.Gmi_node_ids # dependent
                for n1 in nids1:
                    lines_rigid.append([n1, n2])
            elif elem.type in ['RBAR', 'RBAR1', 'RROD']: ## TODO: these aren't quite right
                dependent = elem.Ga()
                independent = elem.Gb()
                lines_rigid.append([dependent, independent])
            else:
                print(str(elem))
        return lines_rigid

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
        e11 = np.zeros(mids.shape, dtype='float32')
        e22 = np.zeros(mids.shape, dtype='float32')
        e33 = np.zeros(mids.shape, dtype='float32')

        has_mat8 = False
        has_mat9 = False
        for umid in np.unique(mids):
            if umid == 0:
                continue
            try:
                mat = model.materials[umid]
            except KeyError:
                print('mids = %s' % mids)
                print('mids = %s' % model.materials.keys())
                continue
                #raise
            if mat.type == 'MAT1':
                e11i = e22i = e33i = mat.e
            elif mat.type == 'MAT8':
                e11i = e33i = mat.e11
                e22i = mat.e22
                has_mat8 = True
            elif mat.type in ['MAT11', 'MAT3D']:
                e11i = mat.e1
                e22i = mat.e2
                e33i = mat.e3
                has_mat9 = True
            else:
                print('skipping %s' % mat)
                continue
                #raise NotImplementedError(mat)
            #print('mid=%s e11=%e e22=%e' % (umid, e11i, e22i))
            i = np.where(umid == mids)[0]
            e11[i] = e11i
            e22[i] = e22i
            e33[i] = e33i
        return has_mat8, has_mat9, e11, e22, e33

    def get_pressure_array(self, model, load_case):
        eids = self.element_ids

        # account for scale factors
        loads2 = []
        scale_factors2 = []
        for load in load_case:
            if isinstance(load, LOAD):
                scale_factors, loads = load.get_reduced_loads()
                scale_factors2 += scale_factors
                loads2 += loads
            else:
                scale_factors2.append(1.)
                loads2.append(load)

        pressures = np.zeros(len(model.elements), dtype='float32')

        iload = 0
        nloads = len(loads2)
        show_nloads = nloads > 5000
        # loop thru scaled loads and plot the pressure
        for load, scale in zip(loads2, scale_factors2):
            if show_nloads and iload % 5000 == 0:
                self.log_debug('  NastranIOv iload=%s/%s' % (iload, nloads))
            if load.type == 'PLOAD4':
                #print(load.object_attributes())
                for elem in load.eids:
                    #elem = self.model.elements[eid]
                    if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',
                                     'CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        pressure = load.pressures[0] * scale

                        # single element per PLOAD
                        #eid = elem.eid
                        #pressures[eids.index(eid)] = pressure

                        # multiple elements
                        #for elem in load.eids:
                        ie = np.searchsorted(eids, elem.eid)
                        #pressures[ie] += p  # correct; we can't assume model orientation
                        nz = self.normals[ie, 2]  # considers normal of shell
                        pressures[ie] += pressure * nz

                    #elif elem.type in ['CTETRA', 'CHEXA', 'CPENTA']:
                        #A, centroid, normal = elem.get_face_area_centroid_normal(
                            #load.g34.nid, load.g1.nid)
                        #r = centroid - p
            iload += 1
        return pressures

    def get_load_array(self, model, subcase):
        found_load = False
        found_temperature = False

        load_keys = (
            'LOAD', 'TEMPERATURE(MATERIAL)', 'TEMPERATURE(INITIAL)',
            'TEMPERATURE(LOAD)', 'TEMPERATURE(BOTH)')
        temperature_keys = (
            'TEMPERATURE(MATERIAL)', 'TEMPERATURE(INITIAL)',
            'TEMPERATURE(LOAD)', 'TEMPERATURE(BOTH)')


        centroidal_pressures = None
        forces = None
        spcd = None
        temperature_key = None
        temperatures = None
        for key in load_keys:
            try:
                load_case_id = subcase.get_parameter(key)[0]
            except KeyError:
                # print('no %s for isubcase=%s' % (key, subcase_id))
                continue
            try:
                load_case = model.loads[load_case_id]
            except KeyError:
                self.log.warning('LOAD=%s not found' % load_case_id)
                continue

            if key == 'LOAD':
                p0 = np.array([0., 0., 0.], dtype='float32')
                centroidal_pressures, forces, spcd = self._get_forces_moments_array(
                    model, p0, load_case_id, include_grav=False)
                found_load = True
            elif key in temperature_keys:
                temperatures = self._get_temperatures_array(model, load_case_id)
                found_temperature = True
                temperature_key = key
            else:
                raise NotImplementedError(key)
        temperature_data = (temperature_key, temperatures)
        load_data = (centroidal_pressures, forces, spcd)
        return found_load, found_temperature, temperature_data, load_data

    def _get_dvprel_ndarrays(self, model, nelements, pids):
        """creates arrays for dvprel results"""
        dvprel_t_init = np.zeros(nelements, dtype='float32')
        dvprel_t_min = np.zeros(nelements, dtype='float32')
        dvprel_t_max = np.zeros(nelements, dtype='float32')
        design_region = np.zeros(nelements, dtype='int32')

        for key, dvprel in iteritems(model.dvprels):
            if dvprel.type == 'DVPREL1':
                prop_type = dvprel.Type
                desvars = dvprel.dvids
                coeffs = dvprel.coeffs
                pid = dvprel.pid.pid
                var_to_change = dvprel.pname_fid
                assert len(desvars) == 1, len(desvars)

                if prop_type == 'PSHELL':
                    i = np.where(pids == pid)
                    design_region[i] = dvprel.oid
                    assert len(i) > 0, i
                    if var_to_change == 'T':
                        #value = 0.
                        lower_bound = 0.
                        upper_bound = 0.
                        for desvar, coeff in zip(desvars, coeffs):
                            xiniti = desvar.xinit
                            if desvar.xlb != -1e20:
                                xiniti = max(xiniti, desvar.xlb)
                                lower_bound = desvar.xlb
                            if desvar.xub != 1e20:
                                xiniti = min(xiniti, desvar.xub)
                                upper_bound = desvar.xub
                            if desvar.delx is not None and desvar.delx != 1e20:
                                desvar.delx
                            if desvar.ddval is not None:
                                msg = 'DESVAR id=%s DDVAL is not None\n%s' % str(desvar)
                            assert desvar.ddval is None, desvar
                            xinit = coeff * xiniti
                        dvprel_t_init[i] = xinit
                        dvprel_t_min[i] = lower_bound
                        dvprel_t_max[i] = upper_bound
                    else:
                        raise NotImplementedError(dvprel)
                else:
                    raise NotImplementedError(dvprel)
            else:
                raise NotImplementedError(dvprel)
            if dvprel.p_max != 1e20:
                dvprel.p_max
            if dvprel.p_min is not None:
                dvprel.p_min
        return dvprel_t_init, dvprel_t_min, dvprel_t_max, design_region

    def _get_forces_moments_array(self, model, p0, load_case_id, include_grav=False):
        nids = sorted(model.nodes.keys())
        nnodes = len(nids)
        nid_map = self.nid_map
        eid_map = self.eid_map

        load_case = model.loads[load_case_id]
        loads2, scale_factors2 = self._get_loads_and_scale_factors(load_case)

        eids = sorted(model.elements.keys())
        centroidal_pressures = np.zeros(len(model.elements), dtype='float32')
        nodal_pressures = np.zeros(len(self.node_ids), dtype='float32')

        forces = np.zeros((nnodes, 3), dtype='float32')
        spcd = np.zeros((nnodes, 3), dtype='float32')
        # loop thru scaled loads and plot the pressure
        cards_ignored = {}

        assert self.normals is not None
        for load, scale in zip(loads2, scale_factors2):
            if load.type == 'FORCE':
                scale2 = load.mag * scale  # does this need a magnitude?
                nid = load.node
                if nid in self.dependents_nodes:
                    print('    nid=%s is a dependent node and has an FORCE applied\n%s' % (
                        nid, str(load)))
                forces[nid_map[nid]] += load.xyz * scale2

            elif load.type == 'PLOAD2':
                pressure = load.pressures[0] * scale  # there are 4 pressures, but we assume p0
                for eid in load.eids:
                    elem = self.elements[eid]
                    if elem.type in ['CTRIA3',
                                     'CQUAD4', 'CSHEAR']:
                        node_ids = elem.node_ids
                        nnodes = len(node_ids)
                        normal = elem.Normal()
                        area = elem.Area()
                        forcei = pressure * normal * area / nnodes
                        # r = elem.Centroid() - p0
                        # m = cross(r, f)
                        for nid in node_ids:
                            if nid in self.dependents_nodes:
                                print('    nid=%s is a dependent node and has an PLOAD2 applied\n'
                                      '%s' % (nid, str(load)))
                            forces[nid_map[nid]] += forcei
                        forces += forcei
                        # F += f
                        # M += m
                    else:
                        self.log.debug('    case=%s etype=%r loadtype=%r not supported' % (
                            load_case_id, elem.type, load.type))

            elif load.type == 'PLOAD4':
                # continue  ## TODO: should be removed
                # elem = load.eid
                # area = elem.get_area()
                if 0:
                    if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',]:
                        eid = elem.eid
                        node_ids = elem.node_ids
                        k = load.pressures[0] * scale / 3.
                        # TODO: doesn't consider load.eids for distributed pressures???
                        for nid in node_ids[3:]:
                            centroidal_pressures[nid_map[nid]] += k
                    elif elem.type in ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        eid = elem.eid
                        node_ids = elem.node_ids
                        k = load.pressures[0] * scale / 4.
                        # TODO: doesn't consider load.eids for distributed pressures???
                        for nid in node_ids[4:]:
                            if nid in self.dependents_nodes:
                                print('    nid=%s is a dependent node and has an PLOAD4 applied\n'
                                      '%s' % (nid, str(load)))
                            centroidal_pressures[nid_map[nid]] += k
                    else:
                        print('    PLOAD4 is unhandled\n%s' % str(load))

                else:
                    # single element per PLOAD
                    #eid = elem.eid
                    #pressures[eids.index(eid)] = p

                    pressure = load.pressures[0] * scale

                    # multiple elements
                    for elem in load.eids:
                        ie = eid_map[elem.eid]
                        nz = self.normals[ie, :]
                        # pressures[eids.index(elem.eid)] += p
                        if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',
                                         # TODO: this was split in bdf_methods...
                                         'CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                            area = elem.get_area()
                            elem_node_ids = elem.node_ids
                            elem_nnodes = len(elem_node_ids)
                            forcei = pressure * area / elem_nnodes
                            for nid in elem_node_ids:
                                if nid in self.dependents_nodes:
                                    print('    nid=%s is a dependent node and has a'
                                          ' PLOAD4 applied\n%s' % (nid, str(load)))
                                #forces[nids.index(nid)] += F
                                i = nid_map[nid]
                                try:
                                    forces[i, :] += forcei * nz
                                except IndexError:
                                    print('i = %s' % i)
                                    print('normals.shape = %s' %  str(self.normals.shape))
                                    print('forces.shape = %s' % str(forces.shape))
                                    print('nz = ', self.normals[i, :])
                                    print('forces[i, :] = ', forces[i, :])
                                    raise
                        else:
                            elem_node_ids = elem.node_ids
                            if elem.type == 'CTETRA':
                                #face1 = elem.get_face(load.g1.nid, load.g34.nid)
                                face, area, centroid, normal = elem.get_face_area_centroid_normal(
                                    load.g1.nid, load.g34.nid)
                                #assert face == face1
                                nface = 3
                            elif elem.type == 'CHEXA':
                                #face1 = elem.get_face(load.g34.nid, load.g1.nid)
                                face, area, centroid, normal = elem.get_face_area_centroid_normal(
                                    load.g34.nid, load.g1.nid)
                                #assert face == face1
                                nface = 4
                            elif elem.type == 'CPENTA':
                                g1 = load.g1.nid
                                if load.g34 is None:
                                    #face1 = elem.get_face(g1)
                                    face, area, centroid, normal = elem.get_face_area_centroid_normal(g1)
                                    nface = 3
                                else:
                                    #face1 = elem.get_face(g1, load.g34.nid)
                                    face, area, centroid, normal = elem.get_face_area_centroid_normal(
                                        g1, load.g34.nid)
                                    nface = 4
                                #assert face == face1
                            else:
                                msg = ('case=%s eid=%s etype=%r loadtype=%r not supported'
                                       % (load_case_id, eid, elem.type, load.type))
                                self.log.debug(msg)
                                continue
                            pressures = load.pressures[:nface]
                            assert len(pressures) == nface
                            if min(pressures) != max(pressures):
                                pressure = np.mean(pressures)
                                #msg = ('%s%s\npressure.min=%s != pressure.max=%s using average'
                                       #' of %%s; load=%s eid=%%s'  % (
                                           #str(load), str(elem), min(pressures), max(pressures),
                                           #load.sid)
                                #print(msg % (pressure, eid))
                            else:
                                pressure = pressures[0]
                            #centroidal_pressures
                            f = pressure * area * normal * scale
                            for inid in face:
                                inidi = nid_map[elem_node_ids[inid]]
                                nodal_pressures[inid] += pressure * scale / nface
                                forces[inidi, :] += f / nface
                            centroidal_pressures[ie] += pressure

                            #r = centroid - p
                            #load.cid.transformToGlobal()
                            #m = cross(r, f)
                            #M += m


                #elif elem.type in ['CTETRA', 'CHEXA', 'CPENTA']:
            elif load.type == 'SPCD':
                #self.gids = [integer(card, 2, 'G1'),]
                #self.constraints = [components_or_blank(card, 3, 'C1', 0)]
                #self.enforced = [double_or_blank(card, 4, 'D1', 0.0)]
                for nid, c1, d1 in zip(load.node_ids, load.constraints, load.enforced):
                    if nid in self.dependents_nodes:
                        self.log_warning('    nid=%s is a dependent node and has an'
                                         ' SPCD applied\n%s' % (nid, str(load)))
                    c1 = int(c1)
                    assert c1 in [1, 2, 3, 4, 5, 6], c1
                    if c1 < 4:
                        spcd[nid_map[nid], c1 - 1] = d1
            else:
                if load.type not in cards_ignored:
                    cards_ignored[load.type] = True
                    self.log_info('  NastranIOv _get_forces_moments_array - unsupported '
                                  'load.type = %s' % load.type)
        return centroidal_pressures, forces, spcd

    def _get_loads_and_scale_factors(self, load_case):
        # account for scale factors
        loads2 = []
        scale_factors2 = []
        for load in load_case:
            if isinstance(load, LOAD):
                scale_factors, loads = load.get_reduced_loads()
                scale_factors2 += scale_factors
                loads2 += loads
            else:
                scale_factors2.append(1.)
                loads2.append(load)
        return loads2, scale_factors2

    def _get_temperatures_array(self, model, load_case_id):
        """builds the temperature array based on thermal cards"""
        nids = sorted(model.nodes.keys())

        load_case = model.loads[load_case_id]
        loads2, scale_factors2 = self._get_loads_and_scale_factors(load_case)
        tempd = model.tempds[load_case_id].temperature if load_case_id in model.tempds else 0.
        temperatures = np.ones(len(model.nodes), dtype='float32') * tempd
        for load, scale in zip(loads2, scale_factors2):
            if load.type == 'TEMP':
                temps_dict = load.temperatures
                for nid, val in iteritems(temps_dict):
                    nidi = nids.index(nid)
                    temperatures[nidi] = val
            else:
                print(load.type)
        return temperatures


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

    cos_theta1 = np.dot(v21, -v13) / (np.linalg.norm(v21) * np.linalg.norm(v13))
    cos_theta2 = np.dot(v32, -v21) / (np.linalg.norm(v32) * np.linalg.norm(v21))
    cos_theta3 = np.dot(v13, -v32) / (np.linalg.norm(v13) * np.linalg.norm(v32))
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
    return areai, max_skew, aspect_ratio, min_thetai, max_thetai, dideal_thetai


def quad_quality(p1, p2, p3, p4):
    """gets the quality metrics for a quad"""
    v21 = p2 - p1
    v32 = p3 - p2
    v43 = p4 - p3
    v14 = p1 - p4

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

    cos_theta1 = np.dot(v21, -v14) / (np.linalg.norm(v21) * np.linalg.norm(v14))
    cos_theta2 = np.dot(v32, -v21) / (np.linalg.norm(v32) * np.linalg.norm(v21))
    cos_theta3 = np.dot(v43, -v32) / (np.linalg.norm(v43) * np.linalg.norm(v32))
    cos_theta4 = np.dot(v14, -v43) / (np.linalg.norm(v14) * np.linalg.norm(v43))
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
    return areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio, min_thetai, max_thetai, dideal_thetai

def get_min_max_theta(faces, all_node_ids, nid_map, xyz_cid0):
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

            cos_theta1 = np.dot(v21, -v13) / (np.linalg.norm(v21) * np.linalg.norm(v13))
            cos_theta2 = np.dot(v32, -v21) / (np.linalg.norm(v32) * np.linalg.norm(v21))
            cos_theta3 = np.dot(v13, -v32) / (np.linalg.norm(v13) * np.linalg.norm(v32))
            cos_thetas.extend([cos_theta1, cos_theta2, cos_theta3])
            ideal_theta.extend([piover3, piover3, piover3])
        elif len(face) == 4:
            try:
                node_ids = all_node_ids[face[0]], all_node_ids[face[1]], all_node_ids[face[2]], all_node_ids[face[3]]
            except:
                print(face)
                print(node_ids)
                raise

            n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
            v21 = xyz_cid0[n2, :] - xyz_cid0[n1, :]
            v32 = xyz_cid0[n3, :] - xyz_cid0[n2, :]
            v43 = xyz_cid0[n4, :] - xyz_cid0[n3, :]
            v14 = xyz_cid0[n1, :] - xyz_cid0[n4, :]
            cos_theta1 = np.dot(v21, -v14) / (np.linalg.norm(v21) * np.linalg.norm(v14))
            cos_theta2 = np.dot(v32, -v21) / (np.linalg.norm(v32) * np.linalg.norm(v21))
            cos_theta3 = np.dot(v43, -v32) / (np.linalg.norm(v43) * np.linalg.norm(v32))
            cos_theta4 = np.dot(v14, -v43) / (np.linalg.norm(v14) * np.linalg.norm(v43))
            cos_thetas.extend([cos_theta1, cos_theta2, cos_theta3, cos_theta4])
            ideal_theta.extend([piover2, piover2, piover2, piover2])
        else:
            raise NotImplementedError(face)
    thetas = np.arccos(cos_thetas)
    ideal_theta = np.array(ideal_theta)
    ideal_thetai = max((thetas - ideal_theta).max(), (ideal_theta - thetas).min())

    min_thetai = thetas.min()
    max_thetai = thetas.max()
    return min_thetai, max_thetai, ideal_thetai
