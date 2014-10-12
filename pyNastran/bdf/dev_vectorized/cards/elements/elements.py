from numpy import (array, zeros, searchsorted, unique, concatenate, argsort,
                   hstack, where, vstack, ones, cross, intersect1d)
from numpy.linalg import norm
from itertools import izip

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

        :param self: the ElementsShell object
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

        
        # spring
        self.elements_spring = model.elements_spring
        self.pelas = model.pelas

        # rods
        self.conrod = model.conrod
        self.prod = model.prod
        self.crod = model.crod

        # mass
        #: stores CONM1, CONM2, CMASS1, CMASS2, CMASS3, CMASS4, CMASS5, PMASS
        self.mass = model.mass

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

    def build(self):
        #print('elements')
        etypes = self._get_element_types(nlimit=False)
        ptypes = self._get_property_types(nlimit=False)
        self.n = 0
        for elems in etypes:
            elems.build()
            self.ne += elems.n
        for props in ptypes:
            props.build()
            self.np += props.n

        eids = check_duplicate('element_id',  etypes)
        check_duplicate('property_id', ptypes)

        self.element_ids = array(list(eids), dtype='int32')
        self.element_ids.sort()

        self.element_groups = self.build_groups(etypes, 'element_id')
        self.property_groups = self.build_groups(ptypes, 'property_id')

    def build_groups(self, objs, name):
        groups = {}
        Types = []
        for obj in objs:
            if hasattr(obj, '_get_types'):
                #print obj.__class__.__name__
                Types2 = obj._get_types(nlimit=False)
                Types += Types2
                #Types += [Type.type for Type in Types2]

            else:
                #Types += [obj.type]
                Types += [obj]
        for Type in Types:
            group_data = getattr(Type, name)
            if len(group_data):
                groups[Type.type] = group_data
        #print("groups = %s" % groups)
        return groups

    def get_elements(self, element_ids):
        TypeMap = {
            'CELAS1'  : self.elements_spring.celas1,
            'CELAS2'  : self.elements_spring.celas2,
            'CELAS3'  : self.elements_spring.celas3,
            'CELAS4'  : self.elements_spring.celas4,

            'CROD'  : self.crod,
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
        for eid in element_ids:
            obj = None
            for Type, eids in self.element_groups.iteritems():
                if eid in eids:
                    i = where(eid == eids)[0]
                    obj = TypeMap[Type][i]
                    out.append(obj)
        return out

    def get_element_ids_by_property_type(self, element_ids):
        Types = [
            self.elements_spring.celas1, self.elements_spring.celas2, self.elements_spring.celas3, self.elements_spring.celas4,
            self.crod, #self.conrod,
            self.cbar, self.cbeam,
            self.elements_shell.ctria3, self.elements_shell.cquad4,
            self.elements_solid.ctetra4, self.elements_solid.cpenta6, self.elements_solid.chexa8,
            self.elements_solid.ctetra10, self.elements_solid.cpenta15, self.elements_solid.chexa20,
        ]
        eids = hstack([Type.element_id for Type in Types])
        pids = hstack([Type.property_id for Type in Types])

        # make sure oids is unique
        oids = hstack([Type.op2_id for Type in Types])
        oids2 = unique(oids)
        assert len(oids) == len(oids2), oids

        oids = hstack([Type.op2_id * ones(Type.n, dtype='int32') for Type in Types])
        i = argsort(eids)
        eids = eids[i]
        pids = pids[i]
        oids = oids[i]

        data = vstack([eids, pids, oids]).T
        #print data

        # drop extra elements
        ie = searchsorted(data[:, 0], element_ids)
        return data[ie, :]

    def get_mass(self, element_ids=None, xyz_cid0=None):
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_positions()

        TypeMap = {
            'CELAS1'  : self.elements_spring.celas1,
            'CELAS2'  : self.elements_spring.celas2,
            'CELAS3'  : self.elements_spring.celas3,
            'CELAS4'  : self.elements_spring.celas4,

            'CROD'  : self.crod,
            'CONROD'  : self.conrod,
            'CSHEAR'  : self.cshear,

            'CQUAD4'  : self.elements_shell.cquad4,
            'CTRIA3'  : self.elements_shell.ctria3,

            'CTETRA4' : self.elements_solid.ctetra4,
            'CPENTA6' : self.elements_solid.cpenta6,
            'CHEXA8'  : self.elements_solid.chexa8,

            'CTETRA10' : self.elements_solid.ctetra10,
            'CPENTA15' : self.elements_solid.cpenta15,
            'CHEXA20'  : self.elements_solid.chexa20,
        }

        if element_ids is None:
            element_ids = self.element_ids
        else:
            # remove invalid elements
            element_ids = intersect1d(element_ids, self.element_ids)
        if len(element_ids) == 0:
            return

        pid_data = self.get_element_ids_by_property_type(element_ids)
        #print(pid_data)
        #print unique(pid_data[:, 1:])

        # find unique groups
        pid_eType = unique2d(pid_data[:, 1:])

        data2 = {}
        eTypeMap = {
            1 : 'CROD', 5: 'CONROD',
            2 : 'CBEAM', 3 : 'CBAR',
            4 : 'CSHEAR',
            10 : 'CELAS1', 11 : 'CELAS2', 12 : 'CELAS3', 13 : 'CELAS4',
            73 : 'CTRIA3', 144 : 'CQUAD4',

            60 : 'CTETRA4', 61 : 'CTETRA10',
            62 : 'CPENTA6', 63 : 'CPENTA15',
            64 : 'CHEXA8', 65 : 'CHEXA20',
        }

        # group elements of the same type by property type
        #   same element type & different property id (e.g. CTRIA3 PSHELL/PCOMP)  -> different group
        #   different element type & same property id (e.g. CTRIA3/CQUAD4 PSHELL) -> different group
        #   same element & same type -> same group
        #
        # we do this in order to think about one property at a time and not
        # have to do a lot of special work to handle different methods for
        # getting the mass
        #print("pid_data = \n%s\n" % str(pid_data))
        #print("pid_eType = \n%s\n" % str(pid_eType))
        for (pid, eType) in pid_eType:
            i = where(pid_data[:, 1] == pid)[0]
            #print("pid=%i eType=%s Step #1=> \n%s\n" % (pid, eType, pid_data[i, :]))
            j = where(pid_data[i, 2] == eType)[0]
            eids = pid_data[i[j], 0]
            #print("pid=%i eType=%s eids=%s Step #2=> \n%s\n" % (pid, eType, eids, pid_data[i[j], :]))
            eType = eTypeMap[eType]
            data2[(pid, eType)] = eids

        #print('data2 = %s' % data2)
        nelements = len(element_ids)
        #print('nelement_ids =', nelements)
        mass = zeros(nelements, dtype='float64')

        ni = 0
        for (pid, eType), element_ids in data2.iteritems():
            #print('pid=%s eType=%s element_ids=%s' % (pid, eType, element_ids))
            prop = self.get_properties([pid])[0]
            #print('  prop = %s' % str(prop).rstrip())

            elements = TypeMap[eType]
            i = searchsorted(elements.element_id, element_ids)
            n = len(i)
            #print('ielements = %s' % i)

            if eType in ['CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                         'CROD', 'CONROD', 'CBAR', 'CBEAM']:
                n1, n2 = elements.node_ids[i, 0], elements.node_ids[i, 1]
                n1 = xyz_cid0[self.model.grid.index_map(n1), :]
                n2 = xyz_cid0[self.model.grid.index_map(n2), :]
                L = norm(n2 - n1, axis=0)
                #print('prop = %s' % prop)
                #print('  calling get_mass_per_area for pid=%s' % (pid))
                mpl = prop.get_mass_per_length()
                mass[ni:ni+n] = mpl * L
            elif eType in ['CTRIA3', 'CQUAD4']:
                if eType == 'CTRIA3':
                    n1, n2, n3 = elements.node_ids[i, 0], elements.node_ids[i, 1], elements.node_ids[i, 2]
                    n1 = xyz_cid0[self.model.grid.index_map(n1), :]
                    n2 = xyz_cid0[self.model.grid.index_map(n2), :]
                    n3 = xyz_cid0[self.model.grid.index_map(n3), :]
                    normal, A = _ctria3_normal_A(n1, n2, n3, calculate_area=True, normalize=True)
                elif eType == 'CQUAD4':
                    n1, n2, n3, n4 = elements.node_ids[i, 0], elements.node_ids[i, 1], elements.node_ids[i, 2], elements.node_ids[i, 3]
                    n1 = xyz_cid0[self.model.grid.index_map(n1), :]
                    n2 = xyz_cid0[self.model.grid.index_map(n2), :]
                    n3 = xyz_cid0[self.model.grid.index_map(n3), :]
                    n4 = xyz_cid0[self.model.grid.index_map(n4), :]
                    normal, A = _cquad4_normal_A(n1, n2, n3, n4, calculate_area=True, normalize=True)
                else:
                    print("Element.get_mass doesn't support %s; try %s.get_mass" % (eType, eType))
                    ni += n
                    continue
                #print('prop = %s' % prop)
                #print('  calling get_mass_per_area for pid=%s' % (pid))
                mpa = prop.get_mass_per_area()
                mass[ni:ni+n] = mpa * A
            elif eType in ['CTETRA4', 'CTETRA10', 'CPENTA6', 'CPENTA15', 'CHEXA8', 'CHEXA20']:
                rho = prop.get_density()
                if eType in ['CTETRA4', 'CTETRA10']:
                    n1, n2, n3, n4 = elements.node_ids[i, 0], elements.node_ids[i, 1], elements.node_ids[i, 2], elements.node_ids[i, 3]
                    n1 = xyz_cid0[self.model.grid.index_map(n1), :]
                    n2 = xyz_cid0[self.model.grid.index_map(n2), :]
                    n3 = xyz_cid0[self.model.grid.index_map(n3), :]
                    n4 = xyz_cid0[self.model.grid.index_map(n4), :]

                    Vi = zeros(n, self.model.float)
                    i = 0
                    for n1i, n2i, n3i, n4i in izip(n1, n2, n3, n4):
                        Vi[i] = volume4(n1i, n2i, n3i, n4i)
                        i += 1
                elif eType in ['CPENTA6', 'CPENTA15']:
                    n1, n2, n3, n4, n5, n6 = (elements.node_ids[i, 0], elements.node_ids[i, 1],
                                              elements.node_ids[i, 2], elements.node_ids[i, 3],
                                              elements.node_ids[i, 4], elements.node_ids[i, 5], )
                    n1 = xyz_cid0[self.model.grid.index_map(n1), :]
                    n2 = xyz_cid0[self.model.grid.index_map(n2), :]
                    n3 = xyz_cid0[self.model.grid.index_map(n3), :]
                    n4 = xyz_cid0[self.model.grid.index_map(n4), :]
                    n5 = xyz_cid0[self.model.grid.index_map(n5), :]
                    n6 = xyz_cid0[self.model.grid.index_map(n6), :]
                    (A1, c1) = tri_area_centroid(n1, n2, n3)
                    (A2, c2) = tri_area_centroid(n4, n5, n6)
                    Vi = (A1 + A2) / 2. * norm(c1 - c2, axis=1)
                elif eType in ['CHEXA8', 'CHEXA20']:
                    n1, n2, n3, n4, n5, n6, n7, n8 = (
                        elements.node_ids[i, 0], elements.node_ids[i, 1],
                        elements.node_ids[i, 2], elements.node_ids[i, 3],
                        elements.node_ids[i, 4], elements.node_ids[i, 5],
                        elements.node_ids[i, 6], elements.node_ids[i, 7], )
                    n1 = xyz_cid0[self.model.grid.index_map(n1), :]
                    n2 = xyz_cid0[self.model.grid.index_map(n2), :]
                    n3 = xyz_cid0[self.model.grid.index_map(n3), :]
                    n4 = xyz_cid0[self.model.grid.index_map(n4), :]
                    n5 = xyz_cid0[self.model.grid.index_map(n5), :]
                    n6 = xyz_cid0[self.model.grid.index_map(n6), :]
                    n7 = xyz_cid0[self.model.grid.index_map(n7), :]
                    n8 = xyz_cid0[self.model.grid.index_map(n8), :]

                    (A1, c1) = quad_area_centroid(n1, n2, n3, n4)
                    (A2, c2) = quad_area_centroid(n5, n6, n7, n8)
                    Vi = (A1 + A2) / 2. * norm(c1 - c2, axis=1)
                else:
                    print("Element.get_mass doesn't support %s; try %s.get_mass" % (eType, eType))
                    ni += n
                    continue
                mass[ni:ni+n] = Vi * rho
            else:
                print("  Element.get_mass doesn't support %s; try %s.get_mass" % (eType, eType))
                #raise NotImplementedError("Element.get_mass doesn't support %s; try %s.get_mass" % (eType, eType))
            #print("")
            ni += n
        print('data2 = %s' % data2)
        return mass

    def get_properties(self, property_ids=None):
        property_ids, int_flag = slice_to_iter(property_ids)
        TypeMap = self.get_property_typemap()
        out = []
        #print('property_ids = %s' % property_ids)
        for pid in property_ids:
            for Type, pids in self.property_groups.iteritems():
                #print('pids = %s' % pids)
                if pid in pids:
                    i = where(pid == pids)[0]
                    pids_extract = pids[i]
                    #obj = TypeMap[Type].slice_by_index(i)
                    obj = TypeMap[Type][pids_extract]
                    out.append(obj)
        return out
        return out[0] if int_flag else out

    def _get_element_types(self, nlimit=True):
        """
        :param nlimit: limit the outputs to objects with data
        """
        types = [self.crod, self.elements_spring, self.cbar, self.cbeam, self.elements_shell, self.elements_solid]
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
        types = [self.prod, self.pelas, self.properties_bar, self.properties_beam, self.properties_shell, self.properties_solid]
        if nlimit:
            types2 = []
            for etype in types:
                if etype.n > 0:
                    types2.append(etype)
            types = types2
        return types

    def get_elements(self, element_ids=None):
        if element_ids is None:
            out = []
            for Type, eids in self.element_groups.iteritems():
                for i in range(len(eids)):
                    obj = TypeMap[Type].slice_by_index(i)
                    out.append(obj)
            #return self
        return self.__getitem__(element_ids)

    def get_element_typemap(self):
        TypeMap = {
            'CROD' : self.crod,
            'CONROD' : self.conrod,

            'CBAR' : self.cbar,
            'CBEAM' : self.cbeam,
            'CSHEAR' : self.cshear,

            'CELAS1' : self.elements_spring.celas1,
            'CELAS2' : self.elements_spring.celas2,
            'CELAS3' : self.elements_spring.celas3,
            'CELAS4' : self.elements_spring.celas4,

            'CTRIA3'  : self.elements_shell.ctria3,
            'CQUAD4'  : self.elements_shell.cquad4,
            #'CTRIA6'  : self.elements_shell.ctria6,
            #'CQUAD8'  : self.elements_shell.cquad8,

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
            'PELAS'  : self.pelas,
            'PROD'  : self.prod,
            'PSHEAR'  : self.pshear,

            #'PBUSH'  : self.pbush,

            'PSHELL'  : self.properties_shell.pshell,
            'PCOMP'  : self.properties_shell.pcomp,
            'PCOMPG'  : self.properties_shell.pcompg,

            'PSOLID' : self.properties_solid.psolid,
        }
        return TypeMap

    def __getitem__(self, element_ids):
        TypeMap = self.get_element_typemap()
        element_ids2, int_flag = slice_to_iter(element_ids)

        out = []
        #print('element_ids = %s' % element_ids2)
        for eid in element_ids2:
            #obj = None
            for Type, eids in self.element_groups.iteritems():
                if eid in eids:
                    #print('  found Type=%s' % Type)
                    i = where(eid == eids)[0]
                    #print("    i = %s" % i)
                    obj = TypeMap[Type].slice_by_index(i)
                    #print("    found eid=%s " % obj.element_id)
                    out.append(obj)
                    break
                #else:
                    #out.append(None)
            #print('*obj = %s' % obj)
            #out.append(obj)
            #print('element_ids = %s\n--%s' % (element_ids2, out))
        return out[0] if int_flag else out

def check_duplicate(name, objs):
    unique_vals = set([])
    for obj in objs:
        if hasattr(obj, name):
            vals = getattr(obj, name)
            #print vals
            unique_vals.update(list(vals))
            #print unique_vals
        else:
            #print "  %s has no %s"  % (obj.__class__.__name__, name)
            pass
    #print "unique %s = %s\n" %(name, unique_vals)
    if len(unique_vals) == 0:
        raise RuntimeError("unique %s = %s" %(name, unique_vals))
    return unique_vals

