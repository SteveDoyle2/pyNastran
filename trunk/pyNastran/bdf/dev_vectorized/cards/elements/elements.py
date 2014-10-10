from numpy import array, zeros, searchsorted, unique, concatenate, argsort, hstack, where


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
        print "groups", groups
        return groups

    def get_elements(self, element_ids):
        TypeMap = {
            'CQUAD4'  : self.elements_shell.cquad4,
            'CTRIA3'  : self.elements_shell.cquad4,

            'CTETRA4' : self.elements_solid.ctetra4,
            'CPENTA6' : self.elements_solid.cpenta6,
            'CHEXA8'  : self.elements_solid.chexa8,
        }
        out = []
        for eid in element_ids:
            for Type, eids in self.element_groups.iteritems():
                if eid in eids:
                    i = where(eid == eids)[0]
                    obj = TypeMap[Type][i]
                    out.append(obj)
        return out

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

    def __getitem__(self, element_ids):
        TypeMap = {
            'CQUAD4'  : self.elements_shell.cquad4,
            'CTRIA3'  : self.elements_shell.cquad4,

            'CTETRA4' : self.elements_solid.ctetra4,
            'CPENTA6' : self.elements_solid.cpenta6,
            'CHEXA8'  : self.elements_solid.chexa8,
        }
        int_flag = False
        if isinstance(element_ids, int):
            element_ids2 = [element_ids]
            int_flag = True
        elif isinstance(element_ids, slice):
            if element_ids.step is None:
                element_ids2 = xrange(element_ids.start, element_ids.stop)
            else:
                element_ids2 = xrange(element_ids.start, element_ids.stop, element_ids.step)
        else: # list, ndarray
            element_ids2 = element_ids

        out = []
        for eid in element_ids2:
            for Type, eids in self.element_groups.iteritems():
                if eid in eids:
                    i = where(eid == eids)[0]
                    obj = TypeMap[Type][i]
                    out.append(obj)
                else:
                    out.append(None)
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