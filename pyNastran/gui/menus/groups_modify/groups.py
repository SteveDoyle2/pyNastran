"""
defines:
 - Group
"""
# -*- coding: utf-8 -*-
from __future__ import print_function, unicode_literals

from numpy import ndarray

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


class NodeGroup(object):
    def __init__(self, name, node_str, nodes_pound, editable=True):
        if len(name):
            assert len(name) > 0, name
            assert name[-1] != ' ', name
            assert '\n' not in name, name
            assert '\r' not in name, name
            assert '\t' not in name, name
        self.name = name
        #self.cids = [0]
        if isinstance(node_str, list):
            node_str = ' '.join(str(s) for s in node_str)
        else:
            assert isinstance(node_str, string_types), 'node_str=%r type=%s' % (node_str, type(node_str))
        self.node_str = node_str
        self.nodes_pound = nodes_pound
        self.editable = editable

    @property
    def node_ids(self):
        return parse_patran_syntax(self.node_str, pound=self.nodes_pound)

    @node_ids.setter
    def node_ids(self, nids):
        assert isinstance(nids, ndarray), nids
        self.node_str = _get_collapsed_text(nids).strip()

    def __repr__(self):
        msg = 'NodeGroup:\n'
        msg += '  name: %s\n' % self.name
        msg += '  editable: %s\n' % self.editable
        msg += '  node_str: [%s]\n' % self.node_str
        msg += '  nodes_pound: %s\n' % self.nodes_pound
        return msg

def _get_collapsed_text(values):
    singles, doubles = collapse_colon_packs(values)
    text = ' '.join([str(s) for s in singles]) + ' '
    text += ' '.join([''.join([str(doublei) for doublei in double]) for double in doubles])
    return text
