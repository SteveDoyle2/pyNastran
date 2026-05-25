"""
defines:
 - Group

"""
# -*- coding: utf-8 -*-
import numpy as np

from pyNastran.bdf.utils import parse_patran_syntax  # parse_patran_syntax_dict
from pyNastran.bdf.cards.collpase_card import collapse_colon_packs


class Group:
    def __init__(self, name: str,
                 element_str: str, elements_pound: int,
                 node_str: str, nodes_pound: int,
                 editable: bool=True) -> None:
        """
        Creates an element group

        Parameters
        ----------
        name : str
            the name of the group
        element_str / node_str: str
            a shortened form for storing ids
        elements_pound / nodes_pound: int
            the max id
        editable : bool; default=True
            not sure what this is used for

        """
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
            assert isinstance(element_str, str), 'element_str=%r type=%s' % (element_str, type(element_str))

        if isinstance(node_str, list):
            node_str = ' '.join(str(s) for s in node_str)
        else:
            assert isinstance(node_str, str), 'node_str=%r type=%s' % (node_str, type(node_str))
        self._element_str = element_str
        self.elements_pound = elements_pound
        self._node_str = node_str
        self.nodes_pound = nodes_pound
        self.editable = editable
        self._element_ids_cache: np.ndarray | None = None
        self._node_ids_cache: np.ndarray | None = None

    @property
    def element_str(self) -> str:
        return self._element_str

    @element_str.setter
    def element_str(self, value: str) -> None:
        self._element_str = value
        self._element_ids_cache = None

    @property
    def node_str(self) -> str:
        return self._node_str

    @node_str.setter
    def node_str(self, value: str) -> None:
        self._node_str = value
        self._node_ids_cache = None

    @property
    def element_ids(self) -> np.ndarray:
        if self._element_ids_cache is None:
            self._element_ids_cache = parse_patran_syntax(self._element_str, pound=self.elements_pound)
        return self._element_ids_cache

    @element_ids.setter
    def element_ids(self, eids: np.ndarray) -> None:
        assert isinstance(eids, np.ndarray), eids
        self._element_str = _get_collapsed_text(eids).strip()
        self._element_ids_cache = eids

    @property
    def node_ids(self) -> np.ndarray:
        if self._node_ids_cache is None:
            self._node_ids_cache = parse_patran_syntax(self._node_str, pound=self.nodes_pound)
        return self._node_ids_cache

    @node_ids.setter
    def node_ids(self, nids: np.ndarray) -> None:
        assert isinstance(nids, np.ndarray), nids
        self._node_str = _get_collapsed_text(nids).strip()
        self._node_ids_cache = nids

    def __repr__(self) -> str:
        msg = 'Group:\n'
        msg += '  name: %s\n' % self.name
        msg += '  editable: %s\n' % self.editable
        #msg += '  cids: [%s]\n' % _get_collapsed_text(self.cids).strip()
        msg += '  element_str: [%s]\n' % self.element_str
        #msg += '  elements: [%s]\n' % _get_collapsed_text(self.elements).strip()
        msg += '  elements_pound: %s\n' % self.elements_pound
        return msg


class NodeGroup:
    def __init__(self, name: str, node_str: str,
                 nodes_pound: int, editable: bool=True) -> None:
        """
        Creates a node group

        Parameters
        ----------
        name : str
            the name of the group
        node_str : str
            a shortened form for storing node ids
        nodes_pound : int
            the max id
        editable : bool; default=True
            not sure what this is used for
        """
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
            assert isinstance(node_str, str), 'node_str=%r type=%s' % (node_str, type(node_str))
        self._node_str = node_str
        self.nodes_pound = nodes_pound
        self.editable = editable
        self._node_ids_cache: np.ndarray | None = None

    @property
    def node_str(self) -> str:
        return self._node_str

    @node_str.setter
    def node_str(self, value: str) -> None:
        self._node_str = value
        self._node_ids_cache = None

    @property
    def node_ids(self) -> np.ndarray:
        if self._node_ids_cache is None:
            self._node_ids_cache = parse_patran_syntax(self._node_str, pound=self.nodes_pound)
        return self._node_ids_cache

    @node_ids.setter
    def node_ids(self, nids: np.ndarray):
        assert isinstance(nids, np.ndarray), nids
        self._node_str = _get_collapsed_text(nids).strip()
        self._node_ids_cache = nids

    def __repr__(self) -> str:
        msg = 'NodeGroup:\n'
        msg += '  name: %s\n' % self.name
        msg += '  editable: %s\n' % self.editable
        msg += '  node_str: [%s]\n' % self.node_str
        msg += '  nodes_pound: %s\n' % self.nodes_pound
        return msg


def _get_collapsed_text(values: list[int]) -> str:
    """writes the collapsed text for ``Group`` and ``NodeGroup``"""
    singles, doubles = collapse_colon_packs(values)
    text = ' '.join([str(s) for s in singles]) + ' '
    text += ' '.join([''.join([str(doublei) for doublei in double]) for double in doubles])
    return text
