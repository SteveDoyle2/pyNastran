from numpy import array, zeros, searchsorted, unique, concatenate, argsort, hstack, where, vstack, ones, cross
from numpy.linalg import norm

from pyNastran.bdf.dev_vectorized.utils import slice_to_iter, unique2d


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

        check_duplicate('element_id',  etypes)
        check_duplicate('property_id', ptypes)

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

            'CQUAD4'  : self.elements_shell.cquad4,
            'CTRIA3'  : self.elements_shell.ctria3,

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

        pid_data = self.get_element_ids_by_property_type(element_ids)

        #print(pid_data)
        #print unique(pid_data[:, 1:])

        # find unique groups
        pid_eType = unique2d(pid_data[:, 1:])

        data2 = {}
        eTypeMap = {
            1 : 'CROD',
            2 : 'CBEAM',
            3 : 'CBAR',
            10 : 'CELAS1',
            11 : 'CELAS2',
            12 : 'CELAS3',
            13 : 'CELAS4',

            73 : 'CTRIA3',
            144 : 'CQUAD4',

            60 : 'CTETRA4',
            61 : 'CTETRA10',
            62 : 'CPENTA6',
            63 : 'CPENTA15',
            64 : 'CHEXA8',
            65 : 'CHEXA20',
        }

        # group elements of the same type by property type
        #   same element type & different property id (e.g. CTRIA3 PSHELL/PCOMP)  -> different group
        #   different element type & same property id (e.g. CTRIA3/CQUAD4 PSHELL) -> different group
        #   same element & same type -> same group
        #
        # we do this in order to think about one property at a time and not
        # have to do a lot of special work to handle different methods for
        # getting the mass
        print("pid_data = \n%s\n" % str(pid_data))
        print("pid_eType = \n%s\n" % str(pid_eType))
        for (pid, eType) in pid_eType:
            i = where(pid_data[:, 1] == pid)[0]
            #print("pid=%i eType=%s Step #1=> \n%s\n" % (pid, eType, pid_data[i, :]))
            j = where(pid_data[i, 2] == eType)[0]
            eids = pid_data[i[j], 0]
            print("pid=%i eType=%s eids=%s Step #2=> \n%s\n" % (pid, eType, eids, pid_data[i[j], :]))
            eType = eTypeMap[eType]
            data2[(pid, eType)] = eids
            assert eids.max() < 20, eids

        print('data2 = %s' % data2)
        for (pid, eType), element_ids in data2.iteritems():
            print('pid=%s eType=%s element_ids=%s' % (pid, eType, element_ids))
            prop = self.get_properties([pid])[0]
            print('prop = ' % prop)
            elements = TypeMap[eType]
            i = searchsorted(elements.element_id, element_ids)
            print('ielemenets = %s' % i)
            if eType in ['CTRIA3', 'CQUAD4']:
                if eType == 'CTRIA3':
                    print('all_pid/nodes =\n%s' % vstack([elements.element_id, elements.node_ids]))
                    n1, n2, n3 = elements.node_ids[i, 0], elements.node_ids[i, 1], elements.node_ids[i, 2]

                    n1 = xyz_cid0[self.model.grid.index_map(n1), :]
                    n2 = xyz_cid0[self.model.grid.index_map(n2), :]
                    n3 = xyz_cid0[self.model.grid.index_map(n3), :]
                    normal = cross(n2 - n1, n3 - n1)
                elif eType == 'CQUAD4':
                    n1, n2, n3, n4 = elements.node_ids[i, :]
                    n1 = xyz_cid0[self.model.grid.index_map(n1), :]
                    n2 = xyz_cid0[self.model.grid.index_map(n2), :]
                    n3 = xyz_cid0[self.model.grid.index_map(n3), :]
                    n4 = xyz_cid0[self.model.grid.index_map(n4), :]
                    normal = cross(n3 - n1, n4 - n2)
                else:
                    print("Element.get_mass doesn't support %s; try %s.get_mass" % (eType, eType))
                    continue
                n = norm(normal, axis=0)
                #normal /= n
                A = 0.5 * n
                print('prop = %s' % prop)
                print('calling get_mass_per_area for pid=%s' % (pid))
                mpa = prop.get_mass_per_area()
                mass = mpa * A
            else:
                print("Element.get_mass doesn't support %s; try %s.get_mass" % (eType, eType))
                #raise NotImplementedError("Element.get_mass doesn't support %s; try %s.get_mass" % (eType, eType))

        print data2
        raise NotImplementedError()

    def get_properties(self, property_ids=None):
        property_ids, int_flag = slice_to_iter(property_ids)
        TypeMap = self.get_property_typemap()
        out = []
        print('property_ids = %s' % property_ids)
        for pid in property_ids:
            for Type, pids in self.property_groups.iteritems():
                print('pids = %s' % pids)
                if pid in pids:
                    i = where(pid == pids)[0]
                    obj = TypeMap[Type][i]
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
            'CONROD'  : self.conrod,
            'CROD' : self.crod,
            'CBAR' : self.cbar,
            'CBEAM' : self.cbeam,

            'CQUAD4'  : self.elements_shell.cquad4,
            'CTRIA3'  : self.elements_shell.ctria3,

            'CTETRA4' : self.elements_solid.ctetra4,
            'CPENTA6' : self.elements_solid.cpenta6,
            'CHEXA8'  : self.elements_solid.chexa8,
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
        print('element_ids = %s' % element_ids2)
        for eid in element_ids2:
            #obj = None
            for Type, eids in self.element_groups.iteritems():
                if eid in eids:
                    #print('  found Type=%s' % Type)
                    i = where(eid == eids)[0]
                    print("    i = %s" % i)
                    obj = TypeMap[Type].slice_by_index(i)
                    print("    found eid=%s " % obj.element_id)
                    out.append(obj)
                    break
                #else:
                    #out.append(None)
            #print('*obj = %s' % obj)
            #out.append(obj)
            print('element_ids = %s\n--%s' % (element_ids2, out))
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

