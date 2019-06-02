import numpy as np
from pyNastran.femutils.utils import unique2d
from pyNastran.dev.bdf_vectorized.utils import slice_to_iter
from pyNastran.dev.bdf_vectorized.cards.elements.solid.ctetra4 import volume4
from pyNastran.dev.bdf_vectorized.cards.elements.solid.chexa8 import quad_area_centroid
from pyNastran.dev.bdf_vectorized.cards.elements.solid.cpenta6 import tri_area_centroid

from pyNastran.dev.bdf_vectorized.cards.elements.shell.cquad4 import _cquad4_normal_A
from pyNastran.dev.bdf_vectorized.cards.elements.shell.ctria3 import _ctria3_normal_A
from pyNastran.dev.bdf_vectorized.cards.elements.utils import build_groups, asarray

from pyNastran.dev.bdf_vectorized.cards.vectorized_card import BaseMethods

class Elements(BaseMethods):
    def __init__(self, model):
        """
        Defines the Elements object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        self.model = model
        self.n = 0
        self.nelements = 0
        self.nproperties = 0
        self.element_ids = None
        self.property_ids = None
        self.element_groups = None
        self.property_groups = None


        #: stores PSHELL, PCOMP, PCOMPG
        self.properties_shell = model.properties_shell
        #: stores CTRIA3, CTRIA6, CQUAD4, CQUAD8
        self.elements_shell = model.elements_shell

        # shear
        #: stores CSHEAR
        self.cshear = model.cshear
        #: stores PSHEAR
        self.pshear = model.pshear

        # rigid
        #self.rbe2 = model.rbe2
        #self.rbe3 = model.rbe3

        # spring
        self.elements_spring = model.elements_spring
        self.pelas = model.pelas

        # bushings
        self.cbush = model.cbush
        self.cbush1d = model.cbush1d
        self.cbush2d = model.cbush2d
        self.pbush = model.pbush

        # rods
        self.conrod = model.conrod
        self.prod = model.prod
        self.crod = model.crod
        self.ctube = model.ctube
        self.ptube = model.ptube

        # mass
        #: stores CONM1, CONM2, CMASS1, CMASS2, CMASS3, CMASS4, CMASS5, PMASS
        self.mass = model.mass
        #self.conm1 = model.conm1
        #self.conm2 = model.conm2
        #self.cmass1 = self.cmass1
        #self.cmass1 = self.cmass1
        #self.cmass2 = self.cmass2
        #self.cmass3 = self.cmass3
        #self.cmass4 = self.cmass4
        #self.cmass5 = self.cmass5

        # bars
        #: stores CBAR
        self.cbar = model.cbar
        #: stores PBAR, PBARL
        self.properties_bar = model.properties_bar

        # beams
        #: stores CBEAM
        self.cbeam = model.cbeam
        #: stores PBEAM, PBEAML
        self.properties_beam = model.properties_beam

        # solids
        #: stores CTETRA4, CPENTA6, CHEXA8, CTETRA10, CPENTA15, CHEXA20
        self.elements_solid = model.elements_solid
        #: stores PSOLID, PLSOLID
        self.properties_solid = model.properties_solid


    def validate_nodes(self, elements):
        validate_nodes = False
        if not hasattr(elements, 'node_ids'):
            # this element isn't finished
            return
        if not validate_nodes:
            # no checks
            return

        grids = self.model.grid.node_id
        nids = np.unique(np.ravel(elements.node_ids))
        #nids.sort()
        diff = np.setdiff1d(nids, grids)
        if len(diff):
            eids = []
            # find the bad elements
            for i, eid in enumerate(elements.element_id):
                j = np.intersect1d(diff, elements.node_ids[i, :])
                if len(j):
                    eids.append(eid)

            # prevents really long arrays
            eids = np.array(eids)
            msg = "Couldn't find Node ID: %s, which is requried by %s %s" % (
                diff, elements.type, eids)
            raise RuntimeError(msg)

    def build(self):
        #print('elements')
        self.n = 0
        self._build_elements()
        self._build_properties()


        old_build = False
        if old_build:
            etypes = self._get_element_types(nlimit=False)
            ptypes = self._get_property_types(nlimit=False)

            for elems in etypes:
                #if elems is None:
                    #continue
                if hasattr(elems, 'type'):
                    if elems.type in self.model.card_count:
                        self.model.log.debug('building %s' % elems.__class__.__name__)
                else:
                    #if elems.n:
                    self.model.log.debug('building %s' % elems.__class__.__name__)
                elems.build()
                self.nelements += elems.n
                self.validate_nodes(elems)
                #print(nids - grids[i])

            for props in ptypes:
                #if props is None:
                    #continue
                if hasattr(props, 'type'):
                    if props.type in self.model.card_count:
                        self.model.log.debug('building %s' % props.__class__.__name__)
                else:
                    #if props.n:
                    self.model.log.debug('building %s' % props.__class__.__name__)
                props.build()
                self.nproperties += props.n
        else:
            etypes = self._get_element_types(nlimit=True)
            ptypes = self._get_property_types(nlimit=True)
            self.model.log.debug('etypes = %s' % etypes)
            self.model.log.debug('ptypes = %s' % ptypes)
            for elems in etypes:
                #if elems.type in ['CONROD']:
                    #self.nproperties += elems.n
                self.nelements += elems.n
                self.validate_nodes(elems)
            for props in ptypes:
                self.nproperties += props.n

        self.model.log.debug('finished building %s' % self.__class__.__name__)

        if self.nelements:
            eids = check_duplicate('element_id', etypes, self.model.log)
            self.element_ids = asarray(eids, dtype='int32')
            self.element_ids.sort()
            self.element_groups = build_groups(etypes, 'element_id', is_element=True)
            #self.model.log.info('self.element_groups = %s' % self.element_groups)
        else:
            self.model.log.warning('no elements...')

        if self.nproperties:
            pids = check_duplicate('property_id', ptypes, self.model.log)
            self.property_ids = asarray(pids, dtype='int32')
            self.property_ids.sort()
            self.property_groups = build_groups(ptypes, 'property_id')
            self.model.log.info('self.property_groups = %s' % self.property_groups)
        #print('*****self.element_ids =', self.element_ids)
        #print('*****self.property_ids =', self.property_ids)

    def get_elements(self, element_id):
        type_map = {
            'CELAS1'  : self.elements_spring.celas1,
            'CELAS2'  : self.elements_spring.celas2,
            'CELAS3'  : self.elements_spring.celas3,
            'CELAS4'  : self.elements_spring.celas4,

            #'CBUSH'  : self.elements_bush.cbush,
            #'CBUSH1D'  : self.elements_bush.cbush1d,
            #'CBUSH2D'  : self.elements_bush.cbush2d,

            'CBUSH'  : self.cbush,
            'CBUSH1D'  : self.cbush1d,
            'CBUSH2D'  : self.cbush2d,

            'CROD'  : self.crod,
            'CTUBE'  : self.ctube,
            'CONROD'  : self.conrod,
            'CSHEAR'  : self.cshear,

            'CTRIA3'  : self.elements_shell.ctria3,
            'CQUAD4'  : self.elements_shell.cquad4,
            'CTRIA6'  : self.elements_shell.ctria6,
            'CQUAD8'  : self.elements_shell.cquad8,

            'CTETRA4' : self.elements_solid.ctetra4,
            'CPENTA6' : self.elements_solid.cpenta6,
            'CHEXA8'  : self.elements_solid.chexa8,

            'CTETRA10' : self.elements_solid.ctetra10,
            'CPENTA15' : self.elements_solid.cpenta15,
            'CHEXA20'  : self.elements_solid.chexa20,
        }
        out = []
        for eid in element_id:
            obj = None
            for etype, eids in self.element_groups.items():
                if eid in eids:
                    i = np.where(eid == eids)[0]
                    obj = type_map[etype][i]
                    out.append(obj)
        return out

    def get_element_properties(self, exclude_types=None):
        if exclude_types is None:
            exclude_types = []
        element_objs = [
            self.elements_spring.celas1, self.elements_spring.celas2,
            self.elements_spring.celas3, self.elements_spring.celas4,
            self.cshear,
            self.crod, self.conrod, self.ctube,
            self.cbar, self.cbeam,
            self.cbush, self.cbush1d, self.cbush2d,
            self.elements_shell.ctria3, self.elements_shell.cquad4,
            self.elements_shell.ctria6, self.elements_shell.cquad8,
            self.elements_solid.ctetra4, self.elements_solid.cpenta6, self.elements_solid.chexa8,
            self.elements_solid.ctetra10, self.elements_solid.cpenta15, self.elements_solid.chexa20,
        ]
        if exclude_types is None:
            exclude_types = []

        element_objs2 = []
        for element_obj in element_objs:
            if element_obj.type not in exclude_types:
                element_objs2.append(element_obj)
        element_objs = element_objs2
        del element_objs2
        # this isn't working...
        #element_objs = [element_obj if element_obj.type not in exclude_types
                        #for element_obj in element_objs]
        elements_without_properties = ['CELAS2', 'CELAS4', 'CONROD']
        eids = np.hstack([element_obj.element_id for element_obj in element_objs])
        pids = np.hstack([np.zeros(element_obj.n, dtype='int32')
                          if element_obj.type in elements_without_properties
                          else element_obj.property_id for element_obj in element_objs])
        return element_objs, eids, pids

    def get_element_ids_by_property_type(self, element_ids, exclude_types=None):
        #self.model.log.debug('element_ids = %s' % element_ids)
        Types, eids, pids = self.get_element_properties(exclude_types)

        # remove undefined properties
        #existing_pids = setdiff1d(unique(pids), self.property_ids, assume_unique=True)
        #self.model.log.debug('pids = %s' % pids)
        #self.model.log.debug('self.property_ids = %s' % self.property_ids)
        #self.model.log.debug('existing_pids = %s' % existing_pids)

        # make sure oids is unique
        oids = np.hstack([Type.op2_id for Type in Types])
        oids2 = np.unique(oids)
        assert len(oids) == len(oids2), oids

        oids = np.hstack([Type.op2_id * np.ones(Type.n, dtype='int32') for Type in Types])
        i = np.argsort(eids)
        #self.model.log.debug('i = %s' % i)
        #self.model.log.debug('eids = %s len=%s' % (eids, len(eids)))
        #self.model.log.debug('pids = %s len=%s' % (pids, len(pids)))
        #self.model.log.debug('oids = %s len=%s' % (oids, len(oids)))
        assert len(eids) == len(pids), 'len(eids)=%i len(pids)=%i' % (len(eids), len(pids))
        assert len(eids) == len(oids), 'len(eids)=%i len(oids)=%i' % (len(eids), len(oids))
        eids = eids[i]
        pids = pids[i]
        oids = oids[i]

        data = np.vstack([eids, pids, oids]).T

        #self.model.log.debug(data)

        # drop extra elements
        # for eids greater than the max allowable eid located at data[-1,0],
        # we drop them
        i_less = np.where(data[-1, 0] >= element_ids)[0]
        element_ids = element_ids[i_less]

        # drop more extra elements
        # we're handling cases of skipped elements (e.g. CELASx cards)
        # that have a sorted location in data, but no unique value
        #print('++++++ %s' % element_ids)
        #print('++++++ %s' % data[:, 0])
        ie = np.unique(np.searchsorted(data[:, 0], element_ids))
        #print('ie = %s' % ie)
        #print('dataA \n%s' % data)
        return data[ie, :]
        #return data

    def get_nodes(self, node_id, xyz_cid0, msg=''):
        i = self.model.grid.get_node_index_by_node_id(node_id, msg=msg)
        return xyz_cid0[i, :]

    def _get_element_ids(self, element_ids_orig):
        if element_ids_orig is None:
            element_ids = self.element_ids
            element_ids_orig = element_ids
            #print('A %s' % element_ids)
        else:
            # remove elements that don't exist in the BDF
            #print("self.element_ids = \n%s" % str(self.element_ids))
            #print("element_ids = \n%s" % str(element_ids))
            element_ids_orig = asarray(element_ids_orig)
            element_ids = np.intersect1d(element_ids_orig, self.element_ids)

            # check for invalid IDs
            #n = len(element_ids)
            #i = where(element_ids < 1)[0]  #  < 100000000
            #j = where(element_ids[i] > 99999999)[0]
            #if len(i) + len(j):
                #eids = setdiff1d(element_ids[i[j]], element_ids)
                #eids = '???'
                #raise RuntimeError('0 < element_ids < 100000000 is invalid; eids=%s' % eids)
            #element_ids = element_ids[i[j]]
            #print('B %s' % element_ids)

        return element_ids, element_ids_orig

    def get_mass_by_element_id(self, element_ids_orig=None, xyz_cid0=None, sort_output=True):
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_position_by_node_index()

        element_ids, element_ids_orig = self._get_element_ids(element_ids_orig)
        if len(element_ids) == 0:
            nelements = len(element_ids_orig)
            mass = np.full(nelements, np.nan, 'float64')
            if sort_output:
                i = np.argsort(element_ids_orig)
                #print("i =", i, i.shape)
                #print("element_ids_orig =", element_ids_orig, element_ids_orig.shape)
                return element_ids_orig[i], mass
            else:
                return element_ids_orig, mass
        nelements_orig = len(element_ids_orig)
        #print('eids orig = %s' % element_ids_orig)

        type_map = {
            #'CELAS1' : self.elements_spring.celas1,
            'CELAS2' : self.elements_spring.celas2,
            'CELAS3' : self.elements_spring.celas3,
            'CELAS4' : self.elements_spring.celas4,

            'CBAR'  : self.cbar,
            'CBEAM' : self.cbeam,

            'CROD' : self.crod,
            'CONROD' : self.conrod,
            'CTUBE' : self.ctube,
            'CSHEAR' : self.cshear,

            'CQUAD4' : self.elements_shell.cquad4,
            'CTRIA3' : self.elements_shell.ctria3,

            'CTETRA4' : self.elements_solid.ctetra4,
            'CPENTA6' : self.elements_solid.cpenta6,
            'CHEXA8' : self.elements_solid.chexa8,

            'CTETRA10' : self.elements_solid.ctetra10,
            'CPENTA15' : self.elements_solid.cpenta15,
            'CHEXA20' : self.elements_solid.chexa20,
        }

        exclude_types = [
            'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
            'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
            'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',
            'CBUSH', 'CBUSH1D', 'CBUSH2D',
            ]
        pid_data = self.get_element_ids_by_property_type(element_ids, exclude_types=exclude_types)
        element_ids_to_analyze = pid_data[:, 0]
        data2 = self.group_elements_by_property_type_and_element_type(self, pid_data)

        #self.model.log.debug('**data2 = %s' % data2)
        nelements = len(element_ids)
        #self.model.log.debug('nelement_ids =', nelements)
        eids2 = np.zeros(nelements_orig, dtype='int32')
        #mass = full(nelements, nan, dtype='float64')
        mass = np.full(nelements_orig, np.nan, dtype='float64')
        #self.model.log.debug('mass.shape =', mass.shape)

        ni = 0
        self.model.log.debug('data2 = %s' % data2)
        for (pid, etype), element_ids in data2.items():
            #self.model.log.debug('pid=%s eType=%s element_ids=%s' % (pid, eType, element_ids))
            elements = type_map[etype]
            i = np.searchsorted(elements.element_id, element_ids)
            n = len(i)
            eids2[ni:ni+n] = elements.element_id[i]
            if pid == 0:
                # CONROD
                pass
            else:
                self.model.log.debug('*cat pid = %s' % pid)
                props = self.get_properties([pid])
                if len(props) == 0:
                    # the property doesn't exist
                    self.model.log.debug('Property %i does not exist and is needed by %s eid=%i' % (
                        pid, etype, element_ids[0]))
                    ni += n
                    #print('props = %s' % props)
                    continue

                # we only get one property at a time
                prop = props[0]
                #print('  prop = %s' % str(prop).rstrip())

            elements = type_map[etype]
            i = np.searchsorted(elements.element_id, element_ids)
            n = len(i)
            #print('ielements = %s' % i)

            if etype in ['CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',]:
                pass
            elif etype in ['CROD', 'CONROD', 'CBAR', 'CBEAM']:
                msg = ', which is required for %ss' % etype
                n1 = self.get_nodes(elements.node_ids[i, 0], xyz_cid0, msg=msg)
                n2 = self.get_nodes(elements.node_ids[i, 1], xyz_cid0, msg=msg)
                length = np.linalg.norm(n2 - n1, axis=1)
                #print('prop = %s' % prop)
                #print('  calling get_mass_per_area for pid=%s' % (pid))
                if etype == 'CONROD':
                    rho = self.model.materials.get_density_by_material_id(elements.material_id)
                    mass[ni:ni+n] = length * elements.A[i] * rho[i]  + elements.nsm[i]
                else:
                    mpl = prop.get_mass_per_length_by_property_id()
                    mass[ni:ni+n] = mpl * length
                    del prop
            elif etype in ['CTRIA3', 'CQUAD4', 'CSHEAR']:
                if etype == 'CTRIA3':
                    n1, n2, n3 = (elements.node_ids[i, 0], elements.node_ids[i, 1],
                                  elements.node_ids[i, 2])
                    n1 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n1), :]
                    n2 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n2), :]
                    n3 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n3), :]
                    area = _ctria3_normal_A(n1, n2, n3,
                                            calculate_area=True, normalize=False)[1]
                elif etype in ['CQUAD4', 'CSHEAR']:
                    n1, n2, n3, n4 = (elements.node_ids[i, 0], elements.node_ids[i, 1],
                                      elements.node_ids[i, 2], elements.node_ids[i, 3])
                    n1 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n1), :]
                    n2 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n2), :]
                    n3 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n3), :]
                    n4 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n4), :]
                    area = _cquad4_normal_A(n1, n2, n3, n4,
                                            calculate_area=True, normalize=False)[1]
                else:
                    self.model.log.debug("Element.get_mass doesn't support %s; try %s.get_mass" % (
                        etype, etype))
                    ni += n
                    continue
                #print('prop = %s' % prop)
                #print('  calling get_mass_per_area for pid=%s' % (pid))
                mpa = prop.get_mass_per_area_by_property_id()
                mass[ni:ni+n] = mpa * area
                del prop
            elif etype in ['CTETRA4', 'CTETRA10', 'CPENTA6', 'CPENTA15', 'CHEXA8', 'CHEXA20']:
                rho = prop.get_density_by_property_id()
                del prop
                if etype in ['CTETRA4', 'CTETRA10']:
                    n1, n2, n3, n4 = (elements.node_ids[i, 0], elements.node_ids[i, 1],
                                      elements.node_ids[i, 2], elements.node_ids[i, 3])
                    n1 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n1), :]
                    n2 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n2), :]
                    n3 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n3), :]
                    n4 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n4), :]

                    volumei = np.zeros(n, self.model.float_fmt)
                    i = 0
                    for n1i, n2i, n3i, n4i in zip(n1, n2, n3, n4):
                        volumei[i] = volume4(n1i, n2i, n3i, n4i)
                        i += 1
                elif etype in ['CPENTA6', 'CPENTA15']:
                    n1, n2, n3, n4, n5, n6 = (elements.node_ids[i, 0], elements.node_ids[i, 1],
                                              elements.node_ids[i, 2], elements.node_ids[i, 3],
                                              elements.node_ids[i, 4], elements.node_ids[i, 5], )
                    n1 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n1), :]
                    n2 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n2), :]
                    n3 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n3), :]
                    n4 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n4), :]
                    n5 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n5), :]
                    n6 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n6), :]
                    (area1, centroid1) = tri_area_centroid(n1, n2, n3)
                    (area2, centroid2) = tri_area_centroid(n4, n5, n6)
                    volumei = (area1 + area2) / 2. * np.linalg.norm(centroid1 - centroid2, axis=1)
                elif etype in ['CHEXA8', 'CHEXA20']:
                    n1, n2, n3, n4, n5, n6, n7, n8 = (
                        elements.node_ids[i, 0], elements.node_ids[i, 1],
                        elements.node_ids[i, 2], elements.node_ids[i, 3],
                        elements.node_ids[i, 4], elements.node_ids[i, 5],
                        elements.node_ids[i, 6], elements.node_ids[i, 7], )
                    n1 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n1), :]
                    n2 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n2), :]
                    n3 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n3), :]
                    n4 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n4), :]
                    n5 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n5), :]
                    n6 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n6), :]
                    n7 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n7), :]
                    n8 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n8), :]

                    (area1, centroid1) = quad_area_centroid(n1, n2, n3, n4)
                    (area2, centroid2) = quad_area_centroid(n5, n6, n7, n8)
                    volumei = (area1 + area2) / 2. * np.linalg.norm(centroid1 - centroid2, axis=1)
                else:
                    msg = "Element.get_mass doesn't support %s; try %s.get_mass" % (etype, etype)
                    self.model.log.debug(msg)
                    ni += n
                    continue
                mass[ni:ni+n] = volumei * rho
            else:
                msg = "  Element.get_mass doesn't support %s; try %s.get_mass" % (etype, etype)
                self.model.log.debug(msg)
                msg = "Element.get_mass doesn't support %s; try %s.get_mass" % (etype, etype)
                raise NotImplementedError(msg)
            #self.model.log.debug("")
            ni += n
        #self.model.log.debug('data2 = %s' % data2)
        diff = np.setdiff1d(element_ids_orig, element_ids_to_analyze)
        #print('A = %s' % element_ids_orig)
        #print('B = %s' % element_ids_to_analyze)
        #print('A-B = %s' % (diff))
        #print('eids2 = %s' % eids2)
        #if len(diff):
            #print('****diff = %s' % list(diff))
            #print('eids out = %s' % eids2)
            #print('len(eids2) =', len(eids2))
            #print('len(diff) =', len(diff))
        eids2[ni:] = diff
        if sort_output:
            i = np.argsort(eids2)
            eids2 = eids2[i]
            mass = mass[i]
        return eids2, mass

    def get_properties(self, property_id=None):
        #self.model.log.debug('property_idA = %s' % property_id)
        property_id, int_flag = slice_to_iter(property_id)
        #self.model.log.debug('property_idB = %s int_flag=%s' % (property_id, int_flag))
        type_map = self.get_property_typemap()
        out = []
        #self.model.log.debug('property_id = %s' % property_id)
        for pid in property_id:
            for Type, pids in self.property_groups.items():
                if not isinstance(pid, int):
                    self.model.log.debug('pids = %s' % pids)
                    self.model.log.debug('pid  = %s' % pid)

                if pid in pids:
                    i = np.where(pid == pids)[0]
                    #self.model.log.debug('i = %s, shape=%s' % (i, i.shape))
                    #self.model.log.debug('pids_extract = %s' % pids_extract)
                    #obj = TypeMap[Type].slice_by_index(i)

                    pids_extract = pids[i]
                    #i = TypeMap[Type].get_property_index_by_property_id(pids_extract)
                    #obj = TypeMap[Type][i]
                    obj = type_map[Type].slice_by_property_id(pids_extract)
                    out.append(obj)
        #return out
        self.model.log.info('out = %s' % out)
        #return out
        return out[0] if int_flag else out

    def allocate(self, card_count):
        self.model.log.debug('elements.allocate()')
        etypes = self._get_element_types(nlimit=True)
        for etype in etypes:
            self.model.log.info('etype=%r' % etype)
            netype = card_count[etype]
            if hasattr(etype, 'type'):
                if etype.type in card_count:
                    self.model.log.debug('allocate %s->%s' % (etype.type, netype))
                    etype.allocate(card_count[etype.type])
                    del card_count[etype.type]
                else:
                    self.model.log.debug('skipping %s->%s' % (etype.type, netype))
            else:
                # ElementsSpring
                self.model.log.debug('allocate %s' % etype.__class__.__name__)
                etype.allocate(card_count)

        ptypes = self._get_property_types(nlimit=False)
        for ptype in ptypes:
            if hasattr(ptype, 'type'):
                if ptype.type in card_count:
                    nptype = card_count[ptype]
                    self.model.log.info('allocate %s->%s' % (ptype.type, nptype))
                    ptype.allocate(nptype)
                    del card_count[ptype.type]
                #else:
                    #asdf
            else:
                self.model.log.debug('allocate %s' % ptype.__class__.__name__)
                ptype.allocate(card_count)
                #del card_count[ptype.type]

    def _build_elements(self):
        """builds the macro objects"""
        etypes = [
            #self.crod, self.conrod, self.ctube,
            #self.cbar, self.cbeam,
            #self.cshear,
            self.mass,
            self.elements_spring,
            self.elements_shell,
            self.elements_solid,
        ]
        for etype in etypes:
            etype.build()

    def _build_properties(self):
        """builds the macro objects"""
        ptypes = [
            #self.prod, self.ptube,
            self.properties_bar, self.properties_beam,
            self.properties_shell,
            self.properties_solid,
        ]
        for ptype in ptypes:
            ptype.build()

    def _get_element_types(self, nlimit=True):
        """
        Parameters
        ----------
        nlimit : bool; default=True
            limit the outputs to objects with data
        """
        types = [
            self.crod, self.conrod, self.ctube,
            self.cbar, self.cbeam,
            self.cshear,

            self.mass,

            self.elements_spring,
            #self.elements_spring.celas1,
            #self.elements_spring.celas2,
            #self.elements_spring.celas3,
            #self.elements_spring.celas4,

            #self.cbush,

            # rigid
            #self.rbe2,
            #self.rbe3,

            #self.elements_shell.ctria3,
            #self.elements_shell.cquad4,
            #self.elements_shell.ctria6,
            #self.elements_shell.cquad8,
            self.elements_shell,

            self.elements_solid,
            #self.elements_solid.ctetra4,
            #self.elements_solid.cpenta6,
            #self.elements_solid.chexa8,

            #self.elements_solid.ctetra10,
            #self.elements_solid.cpenta15,
            #self.elements_solid.chexa20,
        ]
        if nlimit:
            types2 = []
            for etype in types:
                if etype is None:
                    continue
                elif not nlimit:
                    types2.append(etype)
                elif etype.n > 0:
                    types2.append(etype)
            types = types2
        else:
            types2 = []
            for etype in types:
                if etype is None:
                    continue
                types2.append(etype)
            types = types2
        return types

    def _get_property_types(self, nlimit=True):
        """
        Parameters
        ----------
        nlimit : bool; default=True
            limit the outputs to objects with data
        """
        types = [
            # 0D
            self.pelas,
            #self.pbush,

            # 1D
            self.prod, self.ptube,
            #self.properties_bar.pbar, self.properties_bar.pbarl,
            #self.properties_beam.pbeam, self.properties_beam.pbeaml,
            self.properties_bar, self.properties_beam,
            self.pshear,

            self.properties_shell,
            #self.properties_shell.pshell,
            #self.properties_shell.pcomp,
            #self.properties_shell.pcompg,

            self.properties_solid,
            #self.properties_solid.psolid,
            #self.properties_solid.plsolid,
        ]
        if nlimit:
            types2 = []
            for etype in types:
                if etype.n > 0:
                    types2.append(etype)
            types = types2
        return types

    def get_element(self, element_id=None):
        #if element_id is None:
            #out = []
            #for Type, eids in self.element_groups.items():
            #    for i in range(len(eids)):
            #        obj = type_map[Type].slice_by_index(i)
            #        out.append(obj)
            #return self
            #return out
        return self.__getitem__(element_id)

    def get_element_typemap(self):
        type_map = {
            # 0D
            'CELAS1' : self.elements_spring.celas1,
            'CELAS2' : self.elements_spring.celas2,
            'CELAS3' : self.elements_spring.celas3,
            'CELAS4' : self.elements_spring.celas4,

            'CBUSH' : self.cbush,
            #'CONM1': self.conm1,
            #'CONM2': self.conm2,
            #'CMASS1': self.cmass1,
            #'CMASS2': self.cmass2,
            #'CMASS3': self.cmass3,
            #'CMASS4': self.cmass4,
            #'CMASS5': self.cmass5,

            # 1D
            'CROD' : self.crod,
            'CONROD' : self.conrod,

            'CBAR' : self.cbar,
            'CBEAM' : self.cbeam,

            # 2D
            'CSHEAR' : self.cshear,
            'CTRIA3'  : self.elements_shell.ctria3,
            'CQUAD4'  : self.elements_shell.cquad4,
            'CTRIA6'  : self.elements_shell.ctria6,
            'CQUAD8'  : self.elements_shell.cquad8,

            # 3D
            'CTETRA4' : self.elements_solid.ctetra4,
            'CPENTA6' : self.elements_solid.cpenta6,
            'CHEXA8'  : self.elements_solid.chexa8,

            'CTETRA10' : self.elements_solid.ctetra10,
            'CPENTA15' : self.elements_solid.cpenta15,
            'CHEXA20'  : self.elements_solid.chexa20,
        }
        return type_map

    def get_property_typemap(self):
        type_map = {
            # 0D
            'PELAS'  : self.pelas,
            'PBUSH'  : self.pbush,
            #'PMASS' : self.pmass,

            # 1D
            'PROD'  : self.prod,
            'PTUBE'  : self.ptube,

            'PBAR'  : self.properties_bar.pbar,
            'PBARL'  : self.properties_bar.pbarl,

            'PBEAM'  : self.properties_beam.pbeam,
            'PBEAML'  : self.properties_beam.pbeaml,

            # 2D
            'PSHEAR'  : self.pshear,
            'PSHELL'  : self.properties_shell.pshell,
            'PCOMP'  : self.properties_shell.pcomp,
            'PCOMPG'  : self.properties_shell.pcompg,

            'PSOLID' : self.properties_solid.psolid,
        }
        return type_map

    def write_card(self, bdf_file, size=8, is_double=True,
                   include_properties=False, interspersed=True):
        interspersed = False
        alternating = False
        include_properties = True
        if interspersed:
            self._write_interspersed_elements_properties(bdf_file, size, is_double)
            #self.conrod.write_card(bdf_file, size)
        else:
            if alternating:
                self._write_alternating_elements_properties(bdf_file, size, is_double)
            else:
                self._write_elements_properties(bdf_file, size, is_double)

    def _write_elements_properties(self, bdf_file, size, is_double):
        """
        writes the elements/properties grouped
        (e.g., PSHELL is near PCOMP, but not CQUAD4)
        """
        bdf_file.write('$ Elements-------------------------------------------------\n')
        self.elements_spring.write_card(bdf_file, is_double)
        self.conrod.write_card(bdf_file)
        self.crod.write_card(bdf_file)
        self.ctube.write_card(bdf_file)
        self.cbush.write_card(bdf_file, size)

        #self.elements_bars.write_card(bdf_file)
        self.cbar.write_card(bdf_file, size)
        self.cbeam.write_card(bdf_file, size)
        self.cshear.write_card(bdf_file, size)
        self.elements_shell.write_card(bdf_file)
        self.elements_solid.write_card(bdf_file)
        self.mass.write_card(bdf_file, size, is_double)

        bdf_file.write('$ Properties-----------------------------------------------\n')
        self.properties_bar.write_card(bdf_file, size)
        self.pelas.write_card(bdf_file, size, is_double)
        self.prod.write_card(bdf_file)
        self.ptube.write_card(bdf_file)
        self.pbush.write_card(bdf_file, size)
        self.properties_beam.write_card(bdf_file, size)
        self.pshear.write_card(bdf_file, size)
        self.properties_shell.write_card(bdf_file, size)
        self.properties_solid.write_card(bdf_file, size)

    def _write_interspersed_elements_properties(self, bdf_file, size, is_double):
        #raise NotImplementedError('interspersed=False')
        #print(self.object_attributes())
        #print(self.element_groups)
        #print(self.element_ids)
        element_objs = self.get_element_typemap()
        emap = {}
        if self.element_ids is None:
            self.build()
        neids = self.element_ids.size
        element_ids = np.ones((neids, 2)) * -1

        i = 0
        for key, element_obj in sorted(element_objs.items()):
            if element_obj.n == 0:
                i += 1
                continue
            emap[i] = element_obj
            i += 1

    def _write_alternating_elements_properties(self, bdf_file, size, is_double):
        """
        writes the elements/properties with properties near the element type
         (e.g., PSHELL is near CQUAD4, but not PCOMP)
        """
        self._write_alternating_elements_properties_0d(bdf_file, size, is_double)
        self._write_alternating_elements_properties_1d(bdf_file, size, is_double)
        self._write_alternating_elements_properties_2d(bdf_file, size, is_double)
        self._write_alternating_elements_properties_3d(bdf_file, size, is_double)
        bdf_file.write('$----------------------------------------------------------\n')

    def _write_alternating_elements_properties_0d(self, bdf_file, size, is_double):
        #self.properties_springs.write_card(bdf_file)
        self.elements_spring.write_card(bdf_file, is_double)
        self.pelas.write_card(bdf_file, size, is_double)

        #self.elements_damper.write_card(bdf_file)
        #self.pdamp.write_card(bdf_file, size)

        if self.mass.n:
            bdf_file.write('$ Mass-----------------------------------------------------\n')
            self.mass.write_card(bdf_file, size, is_double)

    def _write_alternating_elements_properties_1d(self, bdf_file, size, is_double):
        #self.properties_rods.write_card(f)
        #self.elements_rods.write_card(f)
        if self.conrod.n or self.crod.n or self.prod.n or self.ctube.n or self.ptube.n:
            bdf_file.write('$ Rods-----------------------------------------------------\n')
            self.conrod.write_card(bdf_file)
            self.crod.write_card(bdf_file)
            self.prod.write_card(bdf_file)
            self.ctube.write_card(bdf_file)
            self.ptube.write_card(bdf_file)

        if self.cbush.n or self.pbush.n:
            bdf_file.write('$ Bush-----------------------------------------------------\n')
            self.cbush.write_card(bdf_file, size)
            self.pbush.write_card(bdf_file, size)

        if self.cbar.n or self.properties_bar.n:
            bdf_file.write('$ Bars-----------------------------------------------------\n')
            #self.elements_bars.write_card(f)
            self.properties_bar.write_card(bdf_file, size)
            self.cbar.write_card(bdf_file, size)

        if self.cbeam.n or self.properties_beam.n:
            bdf_file.write('$ Beams----------------------------------------------------\n')
            self.properties_beam.write_card(bdf_file, size)
            self.cbeam.write_card(bdf_file, size)

    def _write_alternating_elements_properties_2d(self, bdf_file, size, is_double):
        if self.cshear.n or self.pshear.n:
            bdf_file.write('$ Shear----------------------------------------------------\n')
            self.pshear.write_card(bdf_file, size)
            self.cshear.write_card(bdf_file, size)

        if self.elements_shell.n or self.properties_shell.n:
            bdf_file.write('$ Shell----------------------------------------------------\n')
            self.properties_shell.write_card(bdf_file, size)
            self.elements_shell.write_card(bdf_file)

    def _write_alternating_elements_properties_3d(self, bdf_file, size, is_double):
        if self.elements_solid.n or self.properties_solid.n:
            bdf_file.write('$ Solid----------------------------------------------------\n')
            self.properties_solid.write_card(bdf_file, size)
            self.elements_solid.write_card(bdf_file)

    def __len__(self):
        return self.nelements

    def __getitem__(self, element_ids):
        type_map = self.get_element_typemap()
        element_ids2, int_flag = slice_to_iter(element_ids)

        out = []
        self.model.log.debug('element_ids = %s' % element_ids)
        self.model.log.debug('element_ids_getitem = %s' % element_ids2)
        for eid in element_ids2:
            #obj = None
            for Type, eids in self.element_groups.items():
                if eid in eids:
                    #self.model.log.debug('  found Type=%s' % Type)
                    i = np.where(eid == eids)[0]
                    #self.model.log.debug("    i = %s" % i)
                    obj = type_map[Type].slice_by_index(i)
                    self.model.log.debug("    found eid=%s " % obj.element_id)
                    out.append(obj)
                    #break
                #else:
                    #out.append(None)
            #print('*obj = %s' % obj)
            #out.append(obj)
            #print('element_ids = %s\n--%s' % (element_ids2, out))
        return out[0] if int_flag else out

    def group_elements_by_property_type_and_element_type(self, elements, pid_data):
        """
        group elements of the same type by property type
           same element type & different property id (e.g. CTRIA3 PSHELL/PCOMP)  -> different group
           different element type & same property id (e.g. CTRIA3/CQUAD4 PSHELL) -> different group
           same element & same type -> same group

        we do this in order to think about one property at a time and not
        have to do a lot of special work to handle different methods for
        getting the mass
        """
        # find unique groups
        #print("pid_data = \n%s\n" % str(pid_data))
        pid_etype = unique2d(pid_data[:, 1:])

        data2 = {}
        #print('model %s' % type(model))
        etype_map = self.model._element_type_to_element_name_mapper

        #print("pid_etype = \n%s\n" % str(pid_etype))
        for (pid, etype) in pid_etype:
            if pid not in elements.property_ids:
                self.model.log.debug('Property pid=%s does not exist' % pid)
                #continue
            i = np.where(pid_data[:, 1] == pid)[0]
            #self.model.log.debug("pid=%i etype=%s Step #1=> \n%s\n" % (
                #pid, etype, pid_data[i, :]))
            j = np.where(pid_data[i, 2] == etype)[0]
            eids = pid_data[i[j], 0]
            #self.model.log.debug("pid=%i etype=%s eids=%s Step #2=> \n%s\n" % (
                #pid, etype, eids, pid_data[i[j], :]))
            etype = etype_map[etype]
            data2[(pid, etype)] = eids
        return data2

    def __repr__(self):
        return '<%s object; n=%s>' % (self.type, self.n)


def check_duplicate(name, objs, log):
    unique_vals = set()
    #print("nobjs = %s" % len(objs))
    for obj in objs:
        #print('working on %s' % obj.__class__.__name__)
        if hasattr(obj, name):
            vals = getattr(obj, name)
            if len(vals):
                #log.debug("  %s vals = %s for class %s" % (
                    #name, np.array(vals), obj.__class__.__name__))
                unique_vals.update(list(vals))
            #else:
                #print("  %s has no %s"  % (obj.__class__.__name__, name))

            #print(unique_vals)
        else:
            #print("  %s has no %s"  % (obj.__class__.__name__, name))
            if hasattr(obj, '_get_types'):
                Types = obj._get_types()
                for Type in Types:
                    vals = getattr(Type, name)
                    if len(vals):
                        #print("    %s vals = %s for class %s" % (
                            #name, vals, Type.__class__.__name__))
                        unique_vals.update(list(vals))
                    #else:
                        #print("    %s has no %s"  % (Type.__class__.__name__, name))
            #else:
                #print("    %s has no _get_types"  % (obj.__class__.__name__))

    #print("unique %s = %s\n" %(name, unique_vals))
    if len(unique_vals) == 0:
        raise RuntimeError("unique %s = %s" %(name, unique_vals))
    #print('unique %s = %s' % (name, unique_vals))
    return unique_vals

