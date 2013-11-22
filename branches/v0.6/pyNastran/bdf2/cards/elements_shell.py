from .ctria3 import CTRIA3


class ElementsShell(object):
    def __init__(self, model):
        """
        Defines the ShellProperties object.

        :param self: the ShellProperties object
        :param model: the BDF object
        """
        self.model = model

        self.ctria3 = CTRIA3(self.model)
        #self.cquad4 = CQUAD4(self.model)
        #self.cquad8 = CQUAD8(self.model)

        self._ctria3 = []
        self._cquad4 = []
        self._cquad8 = []

        self._ctria3_comment = []
        self._cquad4_comment = []
        self._cquad8_comment = []

    def build(self):
        if self._ctria3:
            self.ctria3.build(self._ctria3)
            self._ctria3 = []
            self._ctria3_comment = []

        if self._cquad4:
            self.cquad4.build(self._cquad4)
            self._cquad4 = []
            self._cquad4_comment = []

        if self._cquad8:
            self.cquad8.build(self._cquad8)
            self._cquad8 = []
            self._cquad8_comment = []

        #eid = concatenate(pshell.pid, pcomp.pid)
        #unique_eids = unique(eid)
        #if unique_eids != len(eid):
        #    raise RuntimeError('There are duplicate CTRIA3/CQUAD4 IDs...')

    def rebuild(self):
        raise NotImplementedError()

    def add_ctria3(self, card, comment):
        self._ctria3.append(card)
        self._ctria3_comment.append(comment)

    #def add_pcomp(self, card, comment):
        #self._pcomp.append(card)
        #self._pcomp_comment.append(comment)

    #def add_pcomp(self, card, comment):
        #self._pshear.append(card)
        #self._pshear_comment.append(comment)

    def write_bdf(self, f, size=8, eids=None):
        f.write('$ELEMENTS\n')
        self.ctria3.write_bdf(f, size=size, eids=eids)
        #self.cquad4.write_bdf(f, size=size, eids=eids)
        #self.cquad8.write_bdf(f, size=size, eids=eids)

    def get_stats(self):
        msg = []
        types = [self.ctria3]
        for element in types:
            try:
                nele = len(element.eid)
                if nele:
                    msg.append('  %-8s: %i' % (element.type, nele))
            except AttributeError:  # eid is undefined
                pass
        return msg

    def _verify(self):
        self.ctria3._verify()
        #self.ctria6._verify()
        #self.cquad4._verify()
        #self.cquad8._verify()