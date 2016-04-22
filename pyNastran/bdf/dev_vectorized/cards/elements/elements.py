from __future__ import print_function
from numpy import (array, zeros, searchsorted, unique, argsort,
                   hstack, where, vstack, ones, intersect1d, setdiff1d,
                   nan, full, ravel, ndarray, asarray)
from six import iteritems, integer_types
from six.moves import zip
from numpy.linalg import norm

from pyNastran.bdf.dev_vectorized.utils import slice_to_iter, unique2d
from pyNastran.bdf.dev_vectorized.cards.elements.solid.ctetra4 import volume4
from pyNastran.bdf.dev_vectorized.cards.elements.solid.chexa8 import quad_area_centroid
from pyNastran.bdf.dev_vectorized.cards.elements.solid.cpenta6 import tri_area_centroid

from pyNastran.bdf.dev_vectorized.cards.elements.shell.cquad4 import _cquad4_normal_A
from pyNastran.bdf.dev_vectorized.cards.elements.shell.ctria3 import _ctria3_normal_A

class Elements(object):
    def __init__(self, model):
        """
        Defines the ElementsShell object.

        :param model: the BDF object
        """
        self.model = model
        self.ne = 0
        self.np = 0

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
        self.rbe2 = model.rbe2
        self.rbe3 = model.rbe3

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
        nids = unique(ravel(elements.node_ids))
        #nids.sort()
        diff = setdiff1d(nids, grids)
        if len(diff):
            eids = []
            # find the bad elements
            for i, eid in enumerate(elements.element_id):
                j = intersect1d(diff, elements.node_ids[i, :])
                if len(j):
                    eids.append(eid)

            # prevents really long arrays
            eids = array(eids)
            msg = "Couldn't find Node ID: %s, which is requried by %s %s" % (diff, elements.type, eids)
            raise RuntimeError(msg)

    def build(self):
        #print('elements')
        etypes = self._get_element_types(nlimit=False)
        ptypes = self._get_property_types(nlimit=False)
        self.n = 0
        for elems in etypes:
            if hasattr(elems, 'type'):
                if elems.type in self.model.card_count:
                    self.model.log.debug('building %s' % elems.__class__.__name__)
            else:
                #if elems.n:
                self.model.log.debug('building %s' % elems.__class__.__name__)
            elems.build()
            self.ne += elems.n
            self.validate_nodes(elems)
            #print nids - grids[i]

        for props in ptypes:
            if hasattr(props, 'type'):
                if props.type in self.model.card_count:
                    self.model.log.debug('building %s' % props.__class__.__name__)
            else:
                #if props.n:
                self.model.log.debug('building %s' % props.__class__.__name__)
            props.build()
            self.np += props.n
        self.model.log.debug('finished building %s' % self.__class__.__name__)
        if self.ne:
            eids = check_duplicate('element_id', etypes)
            self.element_ids = array(list(eids), dtype='int32')
            self.element_ids.sort()
            self.element_groups = self.build_groups(etypes, 'element_id', is_element=True)
        if self.np:
            pids = check_duplicate('property_id', ptypes)
            self.property_ids = array(list(pids), dtype='int32')
            self.property_ids.sort()
            self.property_groups = self.build_groups(ptypes, 'property_id')
            self.model.log.info('self.property_groups = %s' % self.property_groups)
        #print('*****self.element_ids =', self.element_ids)
        #print('*****self.property_ids =', self.property_ids)

    def build_groups(self, objs, name, is_element=False):
        groups = {}
        Types = []
        for obj in objs:
            if hasattr(obj, '_get_types'):
                #print(obj.__class__.__name__)
                Types2 = obj._get_types(nlimit=False)
                Types += Types2
                #Types += [Type.type for Type in Types2]

            else:
                #Types += [obj.type]
                Types += [obj]
        for Type in Types:
            group_data = getattr(Type, name)
            if not isinstance(group_data, ndarray):
                raise RuntimeError('Type %s does not return an ndarray when %s is requested' % (Type.type, name))

            assert Type.__class__.__name__ == Type.type, 'class %s has a type of %r' % (Type.__class__.__name__, Type.type)
            #if is_element:
                #assert hasattr(Type, 'op2_id'), 'class %s has no attribute op2_id' % Type.type

            if len(group_data):
                groups[Type.type] = group_data
        #print("groups = %s" % groups)
        return groups

    def get_elements(self, element_id):
        TypeMap = {
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
            for Type, eids in iteritems(self.element_groups):
                if eid in eids:
                    i = where(eid == eids)[0]
                    obj = TypeMap[Type][i]
                    out.append(obj)
        return out

    def get_element_properties(self, exclude_types=None):
        if exclude_types is None:
            exclude_types = []
        Types = [
            self.elements_spring.celas1, self.elements_spring.celas2,
            self.elements_spring.celas3, self.elements_spring.celas4,
            self.cshear,
            self.crod, self.conrod, self.ctube,
            self.cbar, self.cbeam,
            self.elements_shell.ctria3, self.elements_shell.cquad4,
            self.elements_shell.ctria6, self.elements_shell.cquad8,
            self.elements_solid.ctetra4, self.elements_solid.cpenta6, self.elements_solid.chexa8,
            self.elements_solid.ctetra10, self.elements_solid.cpenta15, self.elements_solid.chexa20,
        ]
        if exclude_types is None:
            exclude_types = []

        Types2 = []
        for Type in Types:
            if Type.type not in exclude_types:
                Types2.append(Type)
        Types = Types2
        del Types2
        # this isn't working...
        #Types = [Type if Type.type not in exclude_types for Type in Types]
        elements_without_properties = ['CELAS2', 'CELAS4', 'CONROD']
        eids = hstack([Type.element_id for Type in Types])
        pids = hstack([zeros(Type.n, dtype='int32') if Type.type in elements_without_properties
                       else Type.property_id for Type in Types])
        return Types, eids, pids

    def get_element_ids_by_property_type(self, element_ids, exclude_types=None):
        #print('element_ids = %s' % element_ids)
        Types, eids, pids = self.get_element_properties(exclude_types)

        # remove undefined properties
        #existing_pids = setdiff1d(unique(pids), self.property_ids, assume_unique=True)
        #print('pids = %s' % pids)
        #print('self.property_ids = %s' % self.property_ids)
        #print('existing_pids = %s' % existing_pids)

        # make sure oids is unique
        oids = hstack([Type.op2_id for Type in Types])
        oids2 = unique(oids)
        assert len(oids) == len(oids2), oids

        oids = hstack([Type.op2_id * ones(Type.n, dtype='int32') for Type in Types])
        i = argsort(eids)
        #print('i = %s' % i)
        #print('eids = %s len=%s' % (eids, len(eids)))
        #print('pids = %s len=%s' % (pids, len(pids)))
        #print('oids = %s len=%s' % (oids, len(oids)))
        assert len(eids) == len(pids), 'len(eids)=%i len(pids)=%i' % (len(eids), len(pids))
        assert len(eids) == len(oids), 'len(eids)=%i len(oids)=%i' % (len(eids), len(oids))
        eids = eids[i]
        pids = pids[i]
        oids = oids[i]

        data = vstack([eids, pids, oids]).T

        #print(data)

        # drop extra elements
        # for eids greater than the max allowable eid located at data[-1,0],
        # we drop them
        i_less = where(data[-1, 0] >= element_ids)[0]
        element_ids = element_ids[i_less]

        # drop more extra elements
        # we're handling cases of skipped elements (e.g. CELASx cards)
        # that have a sorted location in data, but no unique value
        #print('++++++ %s' % element_ids)
        #print('++++++ %s' % data[:, 0])
        ie = unique(searchsorted(data[:, 0], element_ids))
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
            element_ids = intersect1d(element_ids_orig, self.element_ids)

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
            mass = full(nelements, nan, 'float64')
            if sort_output:
                i = argsort(element_ids_orig)
                #print("i =", i, i.shape)
                #print("element_ids_orig =", element_ids_orig, element_ids_orig.shape)
                return element_ids_orig[i], mass
            else:
                return element_ids_orig, mass
        nelements_orig = len(element_ids_orig)
        #print('eids orig = %s' % element_ids_orig)

        TypeMap = {
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
        eids2 = zeros(nelements_orig, dtype='int32')
        #mass = full(nelements, nan, dtype='float64')
        mass = full(nelements_orig, nan, dtype='float64')
        #self.model.log.debug('mass.shape =', mass.shape)

        ni = 0
        self.model.log.debug('data2 = %s' % data2)
        for (pid, eType), element_ids in iteritems(data2):
            #self.model.log.debug('pid=%s eType=%s element_ids=%s' % (pid, eType, element_ids))
            elements = TypeMap[eType]
            i = searchsorted(elements.element_id, element_ids)
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
                    self.model.log.debug('Property %i does not exist and is needed by %s eid=%i' % (pid, eType, element_ids[0]))
                    ni += n
                    #print('props = %s' % props)
                    continue

                # we only get one property at a time
                prop = props[0]
                #print('  prop = %s' % str(prop).rstrip())

            elements = TypeMap[eType]
            i = searchsorted(elements.element_id, element_ids)
            n = len(i)
            #print('ielements = %s' % i)

            if eType in ['CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',]:
                pass
            elif eType in ['CROD', 'CONROD', 'CBAR', 'CBEAM']:
                msg = 'which is required for %ss' % eType
                n1 = self.get_nodes(elements.node_ids[i, 0], xyz_cid0, msg=msg)
                n2 = self.get_nodes(elements.node_ids[i, 1], xyz_cid0, msg=msg)
                L = norm(n2 - n1, axis=1)
                #print('prop = %s' % prop)
                #print('  calling get_mass_per_area for pid=%s' % (pid))
                if eType in ['CONROD']:
                    rho = self.model.materials.get_density_by_material_id(elements.material_id)
                    mass[ni:ni+n] = L * elements.A[i] * rho[i]  + elements.nsm[i]
                else:
                    mpl = prop.get_mass_per_length_by_property_id()
                    mass[ni:ni+n] = mpl * L
                    del prop
            elif eType in ['CTRIA3', 'CQUAD4', 'CSHEAR']:
                if eType == 'CTRIA3':
                    n1, n2, n3 = elements.node_ids[i, 0], elements.node_ids[i, 1], elements.node_ids[i, 2]
                    n1 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n1), :]
                    n2 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n2), :]
                    n3 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n3), :]
                    normal, A = _ctria3_normal_A(n1, n2, n3, calculate_area=True, normalize=True)
                elif eType in ['CQUAD4', 'CSHEAR']:
                    n1, n2, n3, n4 = elements.node_ids[i, 0], elements.node_ids[i, 1], elements.node_ids[i, 2], elements.node_ids[i, 3]
                    n1 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n1), :]
                    n2 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n2), :]
                    n3 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n3), :]
                    n4 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n4), :]
                    normal, A = _cquad4_normal_A(n1, n2, n3, n4, calculate_area=True, normalize=True)
                else:
                    self.model.log.debug("Element.get_mass doesn't support %s; try %s.get_mass" % (eType, eType))
                    ni += n
                    continue
                #print('prop = %s' % prop)
                #print('  calling get_mass_per_area for pid=%s' % (pid))
                mpa = prop.get_mass_per_area_by_property_id()
                mass[ni:ni+n] = mpa * A
                del prop
            elif eType in ['CTETRA4', 'CTETRA10', 'CPENTA6', 'CPENTA15', 'CHEXA8', 'CHEXA20']:
                rho = prop.get_density_by_property_id()
                del prop
                if eType in ['CTETRA4', 'CTETRA10']:
                    n1, n2, n3, n4 = elements.node_ids[i, 0], elements.node_ids[i, 1], elements.node_ids[i, 2], elements.node_ids[i, 3]
                    n1 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n1), :]
                    n2 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n2), :]
                    n3 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n3), :]
                    n4 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n4), :]

                    Vi = zeros(n, self.model.float)
                    i = 0
                    for n1i, n2i, n3i, n4i in zip(n1, n2, n3, n4):
                        Vi[i] = volume4(n1i, n2i, n3i, n4i)
                        i += 1
                elif eType in ['CPENTA6', 'CPENTA15']:
                    n1, n2, n3, n4, n5, n6 = (elements.node_ids[i, 0], elements.node_ids[i, 1],
                                              elements.node_ids[i, 2], elements.node_ids[i, 3],
                                              elements.node_ids[i, 4], elements.node_ids[i, 5], )
                    n1 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n1), :]
                    n2 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n2), :]
                    n3 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n3), :]
                    n4 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n4), :]
                    n5 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n5), :]
                    n6 = xyz_cid0[self.model.grid.get_node_index_by_node_id(n6), :]
                    (A1, c1) = tri_area_centroid(n1, n2, n3)
                    (A2, c2) = tri_area_centroid(n4, n5, n6)
                    Vi = (A1 + A2) / 2. * norm(c1 - c2, axis=1)
                elif eType in ['CHEXA8', 'CHEXA20']:
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

                    (A1, c1) = quad_area_centroid(n1, n2, n3, n4)
                    (A2, c2) = quad_area_centroid(n5, n6, n7, n8)
                    Vi = (A1 + A2) / 2. * norm(c1 - c2, axis=1)
                else:
                    self.model.log.debug("Element.get_mass doesn't support %s; try %s.get_mass" % (eType, eType))
                    ni += n
                    continue
                mass[ni:ni+n] = Vi * rho
            else:
                self.model.log.debug("  Element.get_mass doesn't support %s; try %s.get_mass" % (eType, eType))
                raise NotImplementedError("Element.get_mass doesn't support %s; try %s.get_mass" % (eType, eType))
            #self.model.log.debug("")
            ni += n
        #self.model.log.debug('data2 = %s' % data2)
        diff = setdiff1d(element_ids_orig, element_ids_to_analyze)
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
            i = argsort(eids2)
            eids2 = eids2[i]
            mass = mass[i]
        return eids2, mass

    def get_properties(self, property_id=None):
        #self.model.log.debug('property_idA = %s' % property_id)
        property_id, int_flag = slice_to_iter(property_id)
        #self.model.log.debug('property_idB = %s int_flag=%s' % (property_id, int_flag))
        TypeMap = self.get_property_typemap()
        out = []
        #print('property_id = %s' % property_id)
        for pid in property_id:
            for Type, pids in iteritems(self.property_groups):
                if not isinstance(pid, integer_types):
                    self.model.log.debug('pids = %s' % pids)
                    self.model.log.debug('pid  = %s' % pid)

                if pid in pids:
                    i = where(pid == pids)[0]
                    #print('i = %s, shape=%s' % (i, i.shape))
                    #print('pids_extract = %s' % pids_extract)
                    #obj = TypeMap[Type].slice_by_index(i)

                    pids_extract = pids[i]
                    #i = TypeMap[Type].get_property_index_by_property_id(pids_extract)
                    #obj = TypeMap[Type][i]
                    obj = TypeMap[Type].slice_by_property_id(pids_extract)
                    out.append(obj)
        #return out
        #self.model.log.debug('out[0] = %s' % out[0])
        return out[0] if int_flag else out

    def allocate(self, card_count):
        etypes = self._get_element_types(nlimit=False)
        for etype in etypes:
            if hasattr(etype, 'type'):
                if etype.type in card_count:
                    self.model.log.debug('allocate %s' % etype.type)
                    etype.allocate(card_count[etype.type])
                    del card_count[etype.type]
            else:
                # ElementsSpring
                self.model.log.debug('allocate %s' % etype.__class__.__name__)
                etype.allocate(card_count)

        ptypes = self._get_property_types(nlimit=False)
        for ptype in ptypes:
            if hasattr(ptype, 'type'):
                if ptype.type in card_count:
                    self.model.log.debug('allocate %s' % ptype.type)
                    ptype.allocate(card_count[ptype.type])
                    del card_count[ptype.type]
                #else:
                    #asdf
            else:
                self.model.log.debug('allocate %s' % ptype.__class__.__name__)
                ptype.allocate(card_count)
                #del card_count[ptype.type]

    def _get_element_types(self, nlimit=True):
        """
        :param nlimit: limit the outputs to objects with data
        """
        types = [self.crod, self.conrod, self.ctube,
                 self.cbar, self.cbeam,
                 self.cshear,

                 self.mass,

                 self.elements_spring,
                 #self.elements_spring.celas1,
                 #self.elements_spring.celas2,
                 #self.elements_spring.celas3,
                 #self.elements_spring.celas4,

                 self.cbush,

                 # rigid
                 self.rbe2,
                 self.rbe3,

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
                if etype.n > 0:
                    types2.append(etype)
            types = types2
        return types

    def _get_property_types(self, nlimit=True):
        """
        :param nlimit: limit the outputs to objects with data
        """
        types = [
            # 0D
            self.pelas,
            self.pbush,

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
            #for Type, eids in iteritems(self.element_groups):
            #    for i in range(len(eids)):
            #        obj = TypeMap[Type].slice_by_index(i)
            #        out.append(obj)
            #return self
            #return out
        return self.__getitem__(element_id)

    def get_element_typemap(self):
        TypeMap = {
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
        return TypeMap

    def get_property_typemap(self):
        TypeMap = {
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
        return TypeMap

    def write_card(self, f, size=8, is_double=True, include_properties=False, interspersed=True):
        interspersed = False
        include_properties = True
        if interspersed:
            #raise NotImplementedError('interspersed=False')
            self._write_interspersed_elements_properties(f, size)
            self.conrod.write_card(f, size)
        else:
            self._write_alternating_elements_properties(f, size, is_double)

    def _write_alternating_elements_properties(self, f, size, is_double):
        self._write_alternating_elements_properties_0d(f, size, is_double)
        self._write_alternating_elements_properties_1d(f, size, is_double)
        self._write_alternating_elements_properties_2d(f, size, is_double)
        self._write_alternating_elements_properties_3d(f, size, is_double)
        f.write('$----------------------------------------------------------\n')

    def _write_alternating_elements_properties_0d(self, f, size, is_double):
        #self.properties_springs.write_card(f)
        self.elements_spring.write_card(f)
        self.pelas.write_card(f, size)

        #self.elements_damper.write_card(f)
        #self.pdamp.write_card(f, size)

        if self.mass.n:
            f.write('$ Mass-----------------------------------------------------\n')
            self.mass.write_card(f, size)

    def _write_alternating_elements_properties_1d(self, f, size, is_double):
        #self.properties_rods.write_card(f)
        #self.elements_rods.write_card(f)
        if self.conrod.n or self.crod.n or self.prod.n or self.ctube.n or self.ptube.n:
            f.write('$ Rods-----------------------------------------------------\n')
            self.conrod.write_card(f)
            self.crod.write_card(f)
            self.prod.write_card(f)
            self.ctube.write_card(f)
            self.ptube.write_card(f)

        if self.cbush.n or self.pbush.n:
            f.write('$ Bush-----------------------------------------------------\n')
            self.cbush.write_card(f, size)
            self.pbush.write_card(f, size)

        if self.cbar.n or self.properties_bar.n:
            f.write('$ Bars-----------------------------------------------------\n')
            #self.elements_bars.write_card(f)
            self.properties_bar.write_card(f, size)
            self.cbar.write_card(f, size)

        if self.cbeam.n or self.properties_beam.n:
            f.write('$ Beams----------------------------------------------------\n')
            self.properties_beam.write_card(f, size)
            self.cbeam.write_card(f, size)

    def _write_alternating_elements_properties_2d(self, f, size, is_double):
        if self.cshear.n or self.pshear.n:
            f.write('$ Shear----------------------------------------------------\n')
            self.pshear.write_card(f, size)
            self.cshear.write_card(f, size)

        if self.elements_shell.n or self.properties_shell.n:
            f.write('$ Shell----------------------------------------------------\n')
            self.properties_shell.write_card(f, size)
            self.elements_shell.write_card(f)

    def _write_alternating_elements_properties_3d(self, f, size, is_double):
        if self.elements_solid.n or self.properties_solid.n:
            f.write('$ Solid----------------------------------------------------\n')
            self.properties_solid.write_card(f, size)
            self.elements_solid.write_card(f)

    def __len__(self):
        return self.ne

    def __getitem__(self, element_ids):
        TypeMap = self.get_element_typemap()
        element_ids2, int_flag = slice_to_iter(element_ids)

        out = []
        #print('element_ids_getitem = %s' % element_ids2)
        for eid in element_ids2:
            #obj = None
            for Type, eids in iteritems(self.element_groups):
                if eid in eids:
                    #print('  found Type=%s' % Type)
                    i = where(eid == eids)[0]
                    #print("    i = %s" % i)
                    obj = TypeMap[Type].slice_by_index(i)
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
        pid_eType = unique2d(pid_data[:, 1:])

        data2 = {}
        #print('model %s' % type(model))
        eTypeMap = self.model._element_type_to_element_name_mapper

        #print("pid_eType = \n%s\n" % str(pid_eType))
        for (pid, eType) in pid_eType:
            if pid not in elements.property_ids:
                self.model.log.debug('Property pid=%s does not exist' % pid)
                #continue
            i = where(pid_data[:, 1] == pid)[0]
            #self.model.log.debug("pid=%i eType=%s Step #1=> \n%s\n" % (pid, eType, pid_data[i, :]))
            j = where(pid_data[i, 2] == eType)[0]
            eids = pid_data[i[j], 0]
            #self.model.log.debug("pid=%i eType=%s eids=%s Step #2=> \n%s\n" % (pid, eType, eids, pid_data[i[j], :]))
            eType = eTypeMap[eType]
            data2[(pid, eType)] = eids
        return data2

def check_duplicate(name, objs):
    unique_vals = set([])
    #print("nobjs = %s" % len(objs))
    for obj in objs:
        #print('working on %s' % obj.__class__.__name__)
        if hasattr(obj, name):
            vals = getattr(obj, name)
            if len(vals):
                print("  %s vals = %s for class %s" % (name, vals, obj.__class__.__name__))
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
                        #print("    %s vals = %s for class %s" % (name, vals, Type.__class__.__name__))
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

