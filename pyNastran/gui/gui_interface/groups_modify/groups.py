"""
defines:
 - Group
"""
# -*- coding: utf-8 -*-
from __future__ import print_function, unicode_literals
from six import iteritems

from numpy import setdiff1d, unique, hstack, ndarray

from pyNastran.gui.qt_version import qt_version
if qt_version == 4:
    from PyQt4 import QtCore, QtGui
elif qt_version == 5:
    from PyQt5 import QtCore, QtGui
elif qt_version == 'pyside':
    from PySide import QtCore, QtGui
else:
    raise NotImplementedError('qt_version = %r' % qt_version)

from pyNastran.bdf.utils import parse_patran_syntax #, parse_patran_syntax_dict
from pyNastran.bdf.cards.collpase_card import collapse_colon_packs


class Group(object):
    def __init__(self, name, element_str, elements_pound, editable=True):
        if len(name):
            assert len(name) > 0, name
            assert name[-1] != ' ', name
            assert '\n' not in name, name
            assert '\r' not in name, name
            assert '\t' not in name, name
        self.name = name
        #self.cids = [0]
        if isinstance(element_str, list):
            element_str = ' '.join(str(s) for s in element_str)
        else:
            assert isinstance(element_str, (str, unicode)), 'element_str=%r type=%s' % (element_str, type(element_str))
        self.element_str = element_str
        self.elements_pound = elements_pound
        self.editable = editable

    @property
    def element_ids(self):
        return parse_patran_syntax(self.element_str, pound=self.elements_pound)

    @element_ids.setter
    def element_ids(self, eids):
        assert isinstance(eids, ndarray), eids
        self.element_str = _get_collapsed_text(eids).strip()

    def __repr__(self):
        msg = 'Group:\n'
        msg += '  name: %s\n' % self.name
        msg += '  editable: %s\n' % self.editable
        #msg += '  cids: [%s]\n' % _get_collapsed_text(self.cids).strip()
        msg += '  element_str: [%s]\n' % self.element_str
        #msg += '  elements: [%s]\n' % _get_collapsed_text(self.elements).strip()
        msg += '  elements_pound: %s\n' % self.elements_pound
        return msg


def _get_collapsed_text(values):
    singles, doubles = collapse_colon_packs(values)
    text = ' '.join([str(s) for s in singles]) + ' '
    text += ' '.join([''.join([str(doublei) for doublei in double]) for double in doubles])
    return text
