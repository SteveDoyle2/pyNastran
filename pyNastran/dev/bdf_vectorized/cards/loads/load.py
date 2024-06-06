from numpy import asarray

#from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import integer, double


class LOAD:
    type = 'LOAD'

    def allocate(self, ncards):
        pass

    def __init__(self, model):
        """
        Defines the LOAD object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        self.model = model
        self.load_id = None
        self.scale = None

        #: individual scale factors (corresponds to load_ids)
        self.scale_factors = []
        #: individual load_ids (corresponds to scale_factors)
        self.load_ids = []

    def add_from_bdf(self, card, comment):
        """
        Fills the LOAD object from the BDF reader

        :param card:  the BDFCard object
        :param comment: a comment
        """
        if comment:
            self.comment = comment

        #: load ID
        self.load_id = integer(card, 1, 'sid')

        #: overall scale factor
        self.scale = double(card, 2, 'scale')

        # alternating of scale factor & load set ID
        nLoads = len(card) - 3
        assert nLoads % 2 == 0
        for i in range(nLoads // 2):
            n = 2 * i + 3
            self.scale_factors.append(double(card, n, 'scale_factor'))
            self.load_ids.append(integer(card, n + 1, 'load_id'))

    def build(self):
        pass

    def resolve_loads(self):
        raise NotImplementedError()
        types = [self.model.load, self.model.dload]
        types2 = [self.model.force,  #self.model.force1,  self.model.force2,
                  self.model.moment, ]#self.model.moment1, self.model.moment2]
        #for loadtype in types:
            #for load_id, loads in loadtype.items():
                #self.resolve_loads()
        #for loadtype in types2:
            #loadtype.resolve_loads()

    def comment(self):
        if hasattr(self, '_comment'):
            return self.comment
        return []

    def get_stats(self) -> list[str]:
        msg = ['  %-8s: %i' % ('LOAD[%s]' % self.load_id)]
        return msg

    def write_card(self, bdf_file, size=8, is_double=False, load_id=None):
        if load_id is None:
            load_id = self.load_id
        load_id = asarray(load_id)
        list_fields = ['LOAD', self.load_id, self.scale]
        for (scale_factor, lid) in zip(self.scale_factors, self.load_ids):
            list_fields += [scale_factor, lid]
            if size == 8:
                bdf_file.write(print_card_8(list_fields))
            else:
                bdf_file.write(print_card_16(list_fields))
