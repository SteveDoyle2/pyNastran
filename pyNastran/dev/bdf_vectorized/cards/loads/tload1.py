import numpy as np
from numpy import zeros, unique

#from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double_or_blank, integer_string_or_blank,
    integer_double_or_blank)

from pyNastran.dev.bdf_vectorized.cards.loads.vectorized_load import VectorizedLoad

class TLOAD1(VectorizedLoad):
    r"""
    Transient Response Dynamic Excitation, Form 1

    Defines a time-dependent dynamic load or enforced motion of the form:

    .. math::
      \left\{ P(t) \right\} = \left\{ A \right\} \cdot F(t-\tau)

    for use in transient response analysis.
    """
    type = 'TLOAD1'
    def __init__(self, model):
        """
        Defines the TLOAD1 object.

        Parameters
        ----------
        model : BDF
           the BDF object

        .. todo:: collapse loads
        """
        VectorizedLoad.__init__(self, model)

    @property
    def load_id(self):
        return self.sid

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt

            #: SID must be unique for all TLOAD1, TLOAD2, RLOAD1, RLOAD2, and ACSRCE entries.
            self.sid = zeros(ncards, 'int32')

            #: Identification number of DAREA or SPCD entry set or a thermal load
            #: set (in heat transfer analysis) that defines {A}. (Integer > 0)
            self.excite_id = zeros(ncards, 'int32')

            #: If it is a non-zero integer, it represents the
            #: identification number of DELAY Bulk Data entry that defines .
            self.delay_id = zeros(ncards, 'int32')
            #: If it is real, then it directly defines the value of that will
            #: be used for all degrees-of-freedom that are excited by this
            #: dynamic load entry.  See also Remark 9. (Integer >= 0,
            #: real or blank)
            self.delay_tau = zeros(ncards, float_fmt)
            self.is_delay_table = zeros(ncards, 'bool')

            #: Defines the type of the dynamic excitation. (Integer; character
            #: or blank; Default = 0)
            self.Type = zeros(ncards, '|U4')

            #: Identification number of TABLEDi entry that gives F(t). (Integer > 0)
            self.table_id = zeros(ncards, 'int32')

            #: Factor for initial displacements of the enforced degrees-of-freedom.
            #: (Real; Default = 0.0)
            self.us0 = zeros(ncards, float_fmt)

            #: Factor for initial velocities of the enforced degrees-of-freedom
            #: (Real; Default = 0.0)
            self.vs0 = zeros(ncards, float_fmt)

    def get_load_ids(self):
        return unique(self.load_id)

    #def __getitem__(self, i):
        #unique_lid = unique(self.load_id)
        #if len(i):
            #f = FORCE1(self.model)
            #f.load_id = self.load_id[i]
            #f.node_id = self.node_id[i]
            #f.coord_id = self.coord_id[i]
            #f.mag = self.mag[i]
            #f.xyz = self.xyz[i]
            #f.n = len(i)
            #return f
        #raise RuntimeError('len(i) = 0')

    #def __mul__(self, value):
        #obj = TLOAD1(self.model)

        #obj.sid = self.sid
        #obj.excite_id = self.excite_id
        #obj.delay_id = self.delay_id
        #obj.Type = self.Type
        #obj.tau = self.tau
        ##obj.T1 = self.T1
        ##obj.T2 = self.T2
        #obj.frequency = self.frequency
        #obj.phase = self.phase
        #obj.c = self.c
        #obj.b = self.b
        #obj.us0 = self.us0
        #obj.vs0 = self.vs0
        #obj.n = self.n
        #return obj

    #def __rmul__(self, value):
        #return self.__mul__(value)

    def add_card(self, card: BDFCard, comment: str=''):
        i = self.i

        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        delay = integer_double_or_blank(card, 3, 'delay', 0.0)
        Type = integer_string_or_blank(card, 4, 'Type', 'LOAD')
        table_id = integer(card, 5, 'table_id')
        us0 = double_or_blank(card, 6, 'us0', 0.0)
        vs0 = double_or_blank(card, 7, 'vs0', 0.0)

        if Type in [0, 'L', 'LO', 'LOA', 'LOAD']:
            Type = 'LOAD'
        elif Type in [1, 'D', 'DI', 'DIS', 'DISP']:
            Type = 'DISP'
        elif Type in [2, 'V', 'VE', 'VEL', 'VELO']:
            Type = 'VELO'
        elif Type in [3, 'A', 'AC', 'ACC', 'ACCE']:
            Type = 'ACCE'
        elif Type in [5, 6, 7, 12, 13]:
            pass
        else:
            msg = 'invalid TLOAD2 type  Type=%r' % Type
            raise RuntimeError(msg)
        assert len(card) <= 13, 'len(TLOAD2 card) = %i\ncard=%s' % (len(card), card)

        self.sid[i] = sid
        self.excite_id[i] = excite_id
        if isinstance(delay, int):
            self.delay_id[i] = delay
            self.is_delay_table[i] = 1
        elif isinstance(delay, float):
            self.delay_tau[i] = delay
            self.is_delay_table[i] = 0
        else:
            raise NotImplementedError(delay)

        self.Type[i] = Type
        self.table_id[i] = table_id
        self.us0[i] = us0
        self.vs0[i] = vs0
        self.i += 1

    def build(self):
        if self.n:
            i = self.load_id.argsort()
            self.sid = self.sid[i]
            self.excite_id = self.excite_id[i]
            self.delay_id = self.delay_id[i]
            self.delay_tau = self.delay_tau[i]
            self.is_delay_table = self.is_delay_table[i]
            self.Type = self.Type[i]
            self.table_id = self.table_id[i]
            self.us0 = self.us0[i]
            self.vs0 = self.vs0[i]

    def write_card_by_index(self, bdf_file, size=8, is_double=False, i=None):
        for (sidi, excite_idi, delay_idi, delay_taui, is_delay_tablei, load_typei,
             table_idi, us0i, vs0i) in zip(
                 self.sid[i], self.excite_id[i], self.delay_id[i], self.delay_tau[i],
                 self.is_delay_table[i], self.Type[i], self.table_id[i],
                 self.us0[i], self.vs0[i]):

            if is_delay_tablei:
                list_fields = ['TLOAD1', sidi, excite_idi, delay_idi, load_typei,
                               table_idi, us0i, vs0i]
            else:
                list_fields = ['TLOAD1', sidi, excite_idi, delay_taui, load_typei,
                               table_idi, us0i, vs0i]

            if size == 8:
                bdf_file.write(print_card_8(list_fields))
            else:
                bdf_file.write(print_card_16(list_fields))
