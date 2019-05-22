from numpy import searchsorted, hstack, full, nan, unique, where


class PropertiesSolid:
    def __init__(self, model):
        """
        Defines the PropertiesSolid object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        self.model = model

        self.psolid = model.psolid
        self.plsolid = model.plsolid
        self.n = 0

    #def allocate(self, card_count):
        #ptypes = self._get_types(nlimit=False)
        #for ptype in ptypes:
            #if ptype.type in card_count:
                #ptype.allocate(card_count[ptype.type])
                #del card_count[ptype.type]
            #else:
                #assert hasattr(etype, 'allocate'), '%s doesnt support allocate' % ptype.type

    def build(self):
        self.model.log.debug("building solid properties")
        types = self._get_types(nlimit=False)
        for prop in types:
            if prop.n:
                self.model.log.debug('    building %s' % prop.__class__.__name__)
            prop.build()
            self.n += prop.n

        #npsolid = self.psolid.n

        self.property_id = hstack([self.psolid.property_id, self.plsolid.property_id])
        #self.model.log.debug('dtype property_id=%s' % str(self.property_id.dtype))
        #print('npsolid =', npsolid)
        #assert npsolid > 0
        #nplsolid = self.plsolid.n

        #if npshell and npcomp and pshear:
            #asf
        #elif npshell and npcomp:
            #pid = concatenate(self.pshell.pid, self.pcomp.pid)
            #unique_pids = unique(pid)
            #print unique_pids
            #if len(unique_pids) != len(pid):
                #raise RuntimeError('There are duplicate PSHELL/PCOMP IDs...')
        #else:
            #pass

    def rebuild(self):
        raise NotImplementedError()

    def add_psolid(self, card, comment):
        self.psolid.add(card, comment)

    def add_plsolid(self, card, comment):
        self.plsolid.add(card, comment)

    #=========================================================================
    def get_material_id_by_property_id(self, property_id):
        #i = self.get_property_index_by_property_id(property_id)
        #material_id = self.get_material_id_by_property_index(i)
        #return material_id
        assert nan not in property_id, property_id
        n = len(property_id)
        ptypes = self._get_types(nlimit=True)
        material_id = full(n, nan, dtype='int32')
        upids = unique(property_id)
        for ptype in ptypes:
            for upid in upids:
                if upid in ptype.property_id:
                    i = searchsorted(ptype.property_id, upid)
                    j = where(upid == property_id)[0]
                    material_id[j] = ptype.material_id[i]
                    #break

        self.model.log.debug('property_id = %s' % property_id)
        self.model.log.debug('material_id = %s' % material_id)
        assert material_id.shape == (n,), 'material_id.shape=%s; n=%s' % (str(material_id.shape), n)
        return material_id

    #def get_material_id_by_property_id(self, property_id):
        #types = self._get_types(nlimit=True)
        #_material_ids = concatenate([ptype.material_id for ptype in types])
        #print _property_id
        #return _property_id

    def _get_types(self, nlimit=True):
        types = [self.psolid, self.plsolid]
        if nlimit:
            types2 = []
            for ptype in types:
                if ptype.n:
                    types2.append(ptype)
            types = types2
        return types

    #=========================================================================
    def get_stats(self):
        msg = []
        types = self._get_types()
        for prop in types:
            nprop = prop.n
            if nprop:
                msg.append('  %-8s: %i' % (prop.type, nprop))
        return msg

    def write_card(self, bdf_file, size=8, property_id=None):
        bdf_file.write('$PROPERTIES_SOLID\n')
        types = self._get_types()
        for prop in types:
            prop.write_card(bdf_file, size=size, property_id=property_id)

    def __repr__(self):
        msg = '<%s object; n=%s>\n' % (self.__class__.__name__, self.n)
        types = self._get_types()
        for prop in types:
            msg += '  <%s object; n=%s>\n' % (prop.type, prop.n)
        return msg
