
class ElementsDamper:
    def __init__(self, model):
        """
        Defines the ElementsDamper object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        self.model = model

        self.n = 0
        self.cdamp1 = model.cdamp1
        self.cdamp2 = model.cdamp2
        self.cdamp3 = model.cdamp3
        self.cdamp4 = model.cdamp4

    def allocate(self, card_count):
        etypes = self._get_types(nlimit=False)
        for etype in etypes:
            if etype.type in card_count:
                self.model.log.debug('etype.type = %r n=%s' % (etype.type, etype.n))
                etype.allocate(card_count[etype.type])
            #else:
                #assert hasattr(ptype, 'allocate'), '%s doesnt support allocate' % ptype.type

    def build(self):
        types = self._get_types(nlimit=False)
        for elems in types:
            elems.build()
            self.n += elems.n

        #eid = concatenate(pshell.pid, pcomp.pid)
        #unique_eids = unique(eid)
        #if unique_eids != len(eid):
            #raise RuntimeError('There are duplicate CTRIA3/CQUAD4 IDs...')

    def rebuild(self):
        raise NotImplementedError()

    def add_cdamp1(self, card, comment):
        self.cdamp1.add(card, comment)

    def add_cdamp2(self, card, comment):
        self.cdamp2.add(card, comment)

    def add_cdamp3(self, card, comment):
        self.cdamp3.add(card, comment)
        raise NotImplementedError()

    def add_cdamp4(self, card, comment):
        self.cdamp4.add(card, comment)
        raise NotImplementedError()

    def write_card(self, bdf_file, size=8, eids=None):
        bdf_file.write('$ELEMENTS\n')
        types = self._get_types()
        for element in types:
            element.write_card(bdf_file, size=size, eids=eids)

    def _get_types(self, nlimit=True):
        types = [self.cdamp1,
                 self.cdamp2,
                 self.cdamp3,
                 self.cdamp4,
                 ]
        if nlimit:
            types2 = []
            for etype in types:
                if etype.n > 0:
                    #print("etype.Type =", etype.Type)
                    types2.append(etype)
            types = types2
        return types

    def get_stats(self):
        msg = []
        types = self._get_types()
        for element in types:
            nele = element.n
            if nele:
                msg.append('  %-8s: %i' % (element.type, nele))
        return msg

    def _verify(self, xref=True):
        types = self._get_types()
        for elems in types:
            elems._verify(xref=xref)

    #def __repr__(self):
        #return '<%s object; n=%s>' % (self.__class__.__name__, self.n)

    def __repr__(self):
        msg = '<%s object; n=%s>\n' % (self.__class__.__name__, self.n)
        types = self._get_types()
        for elem in types:
            msg += '  <%s object; n=%s>\n' % (elem.type, elem.n)
        return msg
