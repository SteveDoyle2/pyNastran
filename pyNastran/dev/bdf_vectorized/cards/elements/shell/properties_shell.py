from numpy import (concatenate, hstack, unique,
                   array, nan, full, where, isnan)


class PropertiesShell:
    def __init__(self, model):
        """
        Defines the ShellProperties object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        self.model = model
        self.pshell = model.pshell
        self.pcomp = model.pcomp
        self.pcompg = model.pcompg
        self.n = 0
        self.property_id = []

    def allocate(self, card_count):
        ptypes = self._get_types(nlimit=False)
        self.model.log.debug('card_count = %s' % card_count)
        for ptype in ptypes:
            if ptype.type in card_count:
                self.model.log.debug('    allocate %s' % ptype.type)
                ptype.allocate(card_count[ptype.type])
                del card_count[ptype.type]
            #else:
                #assert hasattr(etype, 'allocate'), '%s doesnt support allocate' % ptype.type

    def build(self):
        for prop in [self.pshell, self.pcomp, self.pcompg]:
            prop.build()
            if hasattr(prop, 'n'):
                self.model.log.debug('    building %s; n=%s' % (prop.__class__.__name__, prop.n))
            self.n += prop.n

        npshell = self.pshell.n
        npcomp = self.pcomp.n
        npcompg = self.pcompg.n

        self.n = npshell + npcomp + npcompg
        self.property_id = pid = hstack([self.pshell.property_id, self.pcomp.property_id, self.pcompg.property_id])
        unique_pids = unique(self.property_id)
        #print unique_pids
        if len(unique_pids) != len(self.property_id):
            raise RuntimeError('There are duplicate PSHELL/PCOMP IDs...')

    def get_property_by_property_id(self, property_id):
        ipid = where(self.property_id == property_id)[0]
        if ipid > self.pshell.n:
            if ipid > self.pshell.n + self.pcomp.n:
                aaa
            else:
                bbb
        else:
            return self.pshell.slice_by_property_id(property_id)

    def rebuild(self):
        raise NotImplementedError()

    #=========================================================================
    def add_pshell(self, card, comment):
        self.pshell.add(card, comment)

    def add_pcomp(self, card, comment):
        self.pcomp.add(card, comment)

    def add_pcompg(self, card, comment):
        self.pcompg.add(card, comment)

    #=========================================================================
    def get_material_id_by_property_id(self, property_id):
        types = self._get_types(nlimit=True)
        _property_id = concatenate([ptype.property_id for ptype in types])
        #print _property_id
        return _property_id

    def get_nonstructural_mass_by_property_id(self, property_id=None):
        return self._getmethod(property_id, 'get_nonstructural_mass_by_property_id')

    def get_mass_per_area_by_property_id(self, property_id=None):
        return self._getmethod(property_id, 'get_mass_per_area_by_property_id')

    def get_thickness_by_property_id(self, property_id=None):
        """
        Gets the thickness of the PSHELLs/PCOMPs.

        :param property_id: the property IDs to consider (default=None -> all)
        """
        return self._getmethod(property_id, 'get_thickness_by_property_id')

    def get_density_by_property_id(self, property_id=None):
        """
        Gets the density of the PSHELLs/PCOMPs.

        :param property_ids: the property IDs to consider (default=None -> all)
        """
        return self._getmethod(property_id, 'get_density_by_property_id')

    def _getattr(self, property_id, method):
        TypeMap = {
            'PSHELL': (self.pshell, self.pshell.property_id),
            'PCOMP' : (self.pcomp, self.pcomp.property_id),
            'PCOMPG': (self.pcompg, self.pcompg.property_id),
        }
        out = array([])
        for Type, (ptype, pids) in sorted(TypeMap.items()):
            for pid in pids:
                if pid in property_id:
                    value = getattr(ptype, method)
                    #assert len(value) == 1, value
                    #print('pid =', pid, value)
                    out = hstack([out, value])
        #print("out = %s" % out)
        return out

    def _getmethod(self, property_id, method):
        #types = self._get_types(nlimit=True)
        #data = hstack([getattr(ptype, method)() for ptype in types])
        #if property_ids is None:
            #return data

        TypeMap = {
            'PSHELL': (self.pshell, self.pshell.property_id),
            'PCOMP' : (self.pcomp, self.pcomp.property_id),
            'PCOMPG': (self.pcompg, self.pcompg.property_id),
        }
        #self.model.log.debug('property_id = %s' % property_id)
        n = len(property_id)
        out = full(n, nan, dtype='float64')
        for Type, (ptype, pids) in sorted(TypeMap.items()):
            #self.model.log.debug('Type=%s pids=%s' % (Type, pids))
            for pid in pids:
                if pid in property_id:
                    value = getattr(ptype, method)([pid])
                    j = where(pid == property_id)[0]
                    #assert len(value) == 1, value
                    #self.model.log.debug('pid = %s %s' % (pid, value))
                    out[j] = value
        #self.model.log.debug("out = %s" % out)
        assert out.shape == (n, ), out.shape
        assert isnan(out) == False, out
        return out
        #_property_ids = hstack([ptype.property_id for ptype in types])
        #assert isinstance(property_ids, ndarray), type(property_ids)
        #i = argsort(_property_ids)
        #j = searchsorted(property_ids, _property_ids[i])
        #data_short = data[j]
        #return data_short

    #=========================================================================
    def _get_types(self, nlimit=True):
        types = [self.pshell, self.pcomp, self.pcompg]
        if nlimit:
            types2 = []
            for ptype in types:
                if ptype.n:
                    types2.append(ptype)
            types = types2
        else:
            #self.model.log.debug('ptypes=%s' % types)
            pass
        return types

    def get_stats(self):
        msg = []
        types = self._get_types()
        for prop in types:
            nprop = prop.n
            if nprop:
                msg.append('  %-8s: %i' % (prop.type, nprop))
        return msg

    def write_card(self, bdf_file, size=8, is_double=False, property_id=None):
        bdf_file.write('$PROPERTIES_SHELL\n')
        types = self._get_types()
        for prop in types:
            #print('*SHELL', prop.type)
            prop.write_card(bdf_file, size=size, is_double=is_double, property_id=property_id)

    def __repr__(self):
        msg = '<%s object; n=%s>\n' % (self.__class__.__name__, self.n)
        types = self._get_types()
        for prop in types:
            msg += '  <%s object; n=%s>\n' % (prop.type, prop.n)
        return msg

