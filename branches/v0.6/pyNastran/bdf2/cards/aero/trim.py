from itertools import izip, count
from numpy import array, pi, linspace

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import (BaseCard, expand_thru,
                                          wipe_empty_fields)
from pyNastran.bdf.bdfInterface.assign_type import (fields,
    integer, integer_or_blank,
    double, double_or_blank, 
    string, string_or_blank,
    integer_or_string, double_string_or_blank,
    blank)


class TRIM(object):
    def __init__(self, model):
        self.model = model
    
    def add(self, card, comment):
        #: Trim set identification number. (Integer > 0)
        self.trim_id = integer(card, 1, 'trim_id')

        #: Mach number. (Real > 0.0 and != 1.0)
        self.mach = double(card, 2, 'mach')
        assert self.mach >= 0.0 and self.mach != 1.0, 'mach = %s' % self.mach

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

    def write_bdf(self, f, size=8):
        card = ['TRIM', self.trim_id, self.mach, self.q]
        for (i, label, ux) in izip(count(), self.labels, self.uxs):
            card += [label, ux]
            if i == 1:
                card += [self.aeqr]
        f.write(print_card(card))
