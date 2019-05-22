from itertools import count

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16

from pyNastran.bdf.bdf_interface.assign_type import (integer,
    double, double_or_blank, string, string_or_blank)


class TRIM:
    def __init__(self, model):
        self.model = model

    def add_card(self, card, comment=''):
        #: Trim set identification number. (Integer > 0)
        self.trim_id = integer(card, 1, 'trim_id')

        #: Mach number. (Real > 0.0 and != 1.0)
        self.mach = double(card, 2, 'mach')
        assert self.mach >= 0.0, 'mach = %r' % self.mach
        assert self.mach != 1.0, 'mach = %r' % self.mach

        #: Dynamic pressure. (Real > 0.0)
        self.q = double(card, 3, 'q')
        assert self.q > 0.0, 'q=%s' % self.q

        #: The label identifying aerodynamic trim variables defined on an
        #: AESTAT or AESURF entry.
        self.labels = []

        #: The magnitude of the aerodynamic extra point degree-of-freedom.
        #: (Real)
        self.uxs = []

        label = string_or_blank(card, 4, 'label1')
        if label:
            ux = double(card, 5, 'ux1')
            self.uxs.append(ux)
            self.labels.append(label)

        label = string_or_blank(card, 6, 'label2')
        if label:
            ux = double(card, 7, 'ux1')
            self.uxs.append(ux)
            self.labels.append(label)

        #: Flag to request a rigid trim analysis (Real > 0.0 and < 1.0;
        #: Default = 1.0. A value of 0.0 provides a rigid trim analysis,
        #: not supported
        self.aeqr = double_or_blank(card, 8, 'aeqr', 1.0)

        i = 9
        n = 3
        while i < len(card):
            label = string(card, i, 'label%i' % n)
            ux = double(card, i + 1, 'ux%i' % n)
            self.labels.append(label)
            self.uxs.append(ux)
            i += 2

    def repr_fields(self):
        return self._repr_fields()

    def _repr_fields(self):
        card = ['TRIM', self.trim_id, self.mach, self.q]
        for (i, label, ux) in zip(count(), self.labels, self.uxs):
            card += [label, ux]
            if i == 1:
                card += [self.aeqr]
        return card

    def write_card(self, bdf_file, size=8):
        card = self._repr_fields()
        if size == 8:
            bdf_file.write(print_card_8(card))
        else:
            bdf_file.write(print_card_16(card))
