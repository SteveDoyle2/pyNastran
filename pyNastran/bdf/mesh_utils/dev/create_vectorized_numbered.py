"""
defines:
 - FakeBDFVectorized

"""
from pyNastran.bdf.bdf import BDF as BDF

class FakeBDFVectorized(BDF):
    """
    Renumbers the element ids to be in the same order as the BDFVectorized

    Intended for GUI testing, not anything serious.

    """
    #def __init__(self, debug=True, log=None, mode='msc'):
        #"""see ``BDF.read_bdf``"""
        #BDF.__init__(self, debug=debug, log=log, mode=mode)

    def read_bdf(self, bdf_filename=None, validate=True, xref=True,
                 punch=False, read_includes=True,
                 encoding=None):
        """see ``BDF.read_bdf``"""
        BDF.read_bdf(self, bdf_filename=bdf_filename, validate=validate, xref=False,
                     punch=punch, read_includes=read_includes,
                     encoding=encoding)

        # not done
        etypes = [
            'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
            'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5', 'CVISC',
            'CBUSH', 'CBUSH1D', 'CBUSH2D',
            'CONROD', 'CROD', 'CTUBE', 'PLOTEL',
            'CBAR', 'CBEAM', 'CBEND',
            'CSHEAR',
            'CTRIA3', 'CQUAD4', 'CTRIA6', 'CQUAD8', 'CTRIAR', 'CQUADR',
            'CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM',
            'CHBDYG', 'CHBDYE', 'CHBDYP',
        ]
        eid0 = 1
        #eid_map = {}
        elements2 = {}
        for etype in etypes:
            if etype in self._type_to_id_map:
                eids = self._type_to_id_map[etype]
                print(etype, eids)
                for eid in eids:
                    #eid_map[eid0] = eid
                    if etype == 'PLOTEL':
                        element = self.plotels[eid]
                    else:
                        element = self.elements[eid]
                    element.eid = eid0
                    print(element)
                    elements2[eid0] = element
                    #self.elements[eid] =
                    eid0 += 1

        failed_types = set()
        for elem in self.elements.values():
            if elem.type not in etypes:
                failed_types.add(elem.type)
        self.elements = elements2
        self.plotels = {}

        if failed_types:
            msg = 'The following types were not converted and may result in a bad deck\n'
            msg += '  %s' % failed_types
            self.log.warning(msg)

        # loads are not currently supported
        self.loads = {}  # type: Dict[int, List[Any]]
        self.load_combinations = {}  # type: Dict[int, List[Any]]

def create_vectorized_numbered(bdf_filename_in, bdf_filename_out, debug=True):
    model = FakeBDFVectorized()
    model.read_bdf(bdf_filename=bdf_filename_in, validate=True, xref=False, punch=False,
                   read_includes=True, encoding=None)
    model.write_bdf(bdf_filename_out)
