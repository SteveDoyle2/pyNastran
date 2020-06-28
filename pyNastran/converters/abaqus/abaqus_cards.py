"""
defines:
 - SolidSection
 - Material
 - Part

"""
from __future__ import annotations
from typing import Dict, Optional, Any, TYPE_CHECKING
import numpy as np
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger

allowed_element_types = [
    'r2d2', 'conn2d2',
    'cpe3', 'cpe4', 'cpe4r',
    'cps3', 'cps4', 'cps4r',

    'coh2d4', 'c3d10h', 'cohax4',
    'cax3', 'cax4r', 'mass', 'rotaryi', 't2d2', 'c3d8r',
]

class SolidSection:
    """a SolidSection defines depth and a material"""
    def __init__(self, param_map, data_lines,
                 log: SimpleLogger):
        self.param_map = param_map
        self.data_lines = data_lines
        self.material = param_map['material']
        if len(data_lines) == 0:
            pass
        elif len(data_lines) == 1:
            assert len(data_lines) == 1, data_lines
            line0 = data_lines[0]
            assert len(line0) == 1, data_lines

            try:
                self.thickness = float(line0[0])
            except ValueError:
                self.thickness = 0.

        for line in data_lines:
            log.info('solid - %r' % line)

    def __repr__(self):
        """prints a summary for the solid section"""
        msg = 'SolidSection(\n'
        msg += '    param_map = %r,\n' % self.param_map
        msg += '    thickness = %s,\n' % self.thickness
        msg += ')\n'
        return msg


class Material:
    """a Material object is a series of nodes & elements (of various types)"""
    def __init__(self, name: str,
                 sections: Dict[str, float],
                 density: Optional[float]=None,
                 ndepvars: Optional[int]=None,
                 ndelete: Optional[int]=None):
        self.name = name
        self.density = density

        #self.depvar = None
        self.ndelete = ndelete
        self.ndepvars = ndepvars

        self.user_material = None
        #print(sections)
        #if 'density' in sections:
            #self.density = sections['density']
        #if 'depvar' in sections:
            #self.depvar = sections['depvar']
        #if 'user_material' in sections:
            #self.user_material = sections['user_material']
        self.sections = sections

    def __repr__(self) -> str:
        """prints a summary for the material"""
        msg = 'Material(\n'
        msg += '  name=%r,\n' % self.name
        for key, value in self.sections.items():
            msg += '  %r : %r,\n' % (key, value)
        msg += ')\n'
        return msg

    def write(self, abq_file) -> None:
        """
        *Material, name=Glassy828DEA
        *Density
        1180.,
        *Elastic
            2.14078e+09, 0.42718
        *Material, name=MAT1_828DEA_Dam
        *Density
        1180.,
        *Depvar, delete=4
            20,
        *User Material, constants=16
        ** K      CTELIN      C10          C01       DAM_FORM   FUNC_FORM     EVOLF         EVMF
            3.2e+09, 5.667e-05,  3.75e+08,        0.,        2.,        1.,    50000.,      0.05
            **EVM0ISO    EVM0VOL      EVM0VM DAM_METHOD      ALPHA         A          B          C
                0.,       0.5,       0.5,        1.,        0.,        0.,       0.5,       0.6
        *Material, name=Steel
        *Density
        7800.,
        *Elastic
            2e+11, 0.3
        """
        abq_file.write('*Material, name=%s\n' % write_name(self.name))
        if self.density is not None:
            abq_file.write(f'*Density\n  {self.density},\n')
        if self.ndepvars:
            ndelete = '' if self.ndelete is None else ', delete=%s' % self.ndelete
            abq_file.write(f'*Depvar{ndelete}\n  {self.ndepvars},\n')
        if self.user_material:
            nconstants = ''
            abq_file.write('*User Material%s\n  %s,\n' % (nconstants, self.user_material))
        #abq_file.write('** skipping Material %s\n' % self.name)

class Assembly:
    def __init__(self, element_types, node_sets, element_sets):
        self.element_types = element_types
        self.node_sets = node_sets
        self.element_sets = element_sets

    def write(self, abq_file):
        abq_file.write('** skipping Assembly\n')

    def __repr__(self):
        """summary for the Assembly"""
        etypes = list(self.element_types.keys())
        nsets = list(self.node_sets.keys())
        esets = list(self.element_sets.keys())
        msg = (
            'Assembly:\n'
            '  element_types = %s\n'
            '  node_sets = %s\n'
            '  element_sets = %s\n' % (etypes, nsets, esets)
        )
        return msg

class Part:
    """a Part object is a series of nodes & elements (of various types)"""
    def __init__(self, name: str,
                 nids: np.ndarray,
                 nodes: np.ndarray,
                 element_types: Dict[str, np.ndarray],
                 node_sets: Dict[str, np.ndarray],
                 element_sets: Dict[str, np.ndarray],
                 solid_sections: List[SolidSection],
                 log: SimpleLogger):
        """
        creates a Part object

        Parameters
        ----------
        name : str
            the name
        element_types : Dict[element_type] : node_ids
            element_type : str
                the element type
            bars:
                r2d2 : (nelements, 2) int ndarray
            shells:
                cpe3 : (nelements, 3) int ndarray
                cpe4 : (nelements, 4) int ndarray
                cpe4r : (nelements, 4) int ndarray
                cps3 : (nelements, 3) int ndarray
                cps4 : (nelements, 4) int ndarray
                cps4r : (nelements, 4) int ndarray
                coh2d4 : (nelements, 4) int ndarray
                cohax4 : (nelements, 4) int ndarray
                cax3 : (nelements, 3) int ndarray
                cax4r : (nelements, 4) int ndarray
            solids:
                c3d10h : (nelements, 10) int ndarray

        """
        self.name = name
        self.log = log
        self.node_sets = node_sets
        self.element_sets = element_sets
        self.solid_sections = solid_sections

        try:
            self.nids = np.array(nids, dtype='int32')
        except ValueError:
            msg = 'nids=%s is not integers' % nids
            raise ValueError(msg)
        nnodes = len(self.nids)

        node0 = nodes[0]
        node_shape = len(node0)

        if node_shape == 3:
            self.nodes = np.array(nodes, dtype='float32')
        elif node_shape == 2:
            # abaqus can have only x/y coordinates, so we fake the z coordinate
            self.nodes = np.zeros((nnodes, 3), dtype='float32')
            nodes2 = np.array(nodes, dtype='float32')
            #print(nodes2.shape, self.nodes.shape)
            self.nodes[:, :2] = nodes2
        else:
            raise NotImplementedError(node0)

        # bars
        self.r2d2 = None

        # ---shells---
        # plane strain
        self.cpe3 = None
        self.cpe4 = None
        self.cpe4r = None

        # plane stress
        self.cps3 = None
        self.cps4 = None
        self.cps4r = None

        # other
        self.coh2d4 = None
        self.cohax4 = None
        self.cax3 = None
        self.cax4r = None

        # solids
        self.c3d10h = None
        self.c3d8r = None
        #-----------------------------------
        # eids
        self.r2d2_eids = None

        self.cpe3_eids = None
        self.cpe4_eids = None
        self.cpe4r_eids = None

        self.cps3_eids = None
        self.cps4_eids = None
        self.cps4r_eids = None

        self.coh2d4_eids = None
        self.cohax4_eids = None
        self.cax3_eids = None
        self.cax4r_eids = None

        # rigid elements
        self.c3d10h_eids = None
        self.c3d8r_eids = None
        self._store_elements(element_types)

    def _etypes_nnodes(self):
        """internal helper method"""
        etypes_nnodes = [
            ('r2d2', 2),  #  similar to a CBAR

            #  shells
            ('cpe3', 3),
            ('cpe4', 4),
            ('cpe4r', 4),

            ('cps3', 3),
            ('cps4', 4),
            ('cps4r', 4), # quad, plane stress, reduced

            ('coh2d4', 4), #  cohesive zone
            ('cohax4', 4), #  cohesive zone
            ('cax3', 3),
            ('cax4r', 4),

            #  solids
            ('c3d10h', 10),  #  tet10
            ('c3d8r', 8),  #  hexa8
        ]
        return etypes_nnodes

    def _store_elements(self, element_types):
        """helper method for the init"""
        etypes_nnodes = self._etypes_nnodes()
        for etype, nnodes in etypes_nnodes:
            if etype in element_types:
                etype_eids = '%s_eids' % etype
                elements = element_types[etype]
                eids_elements = np.array(elements, dtype='int32')
                setattr(self, etype, eids_elements)  # r2d2
                setattr(self, etype_eids, eids_elements[:,  0]) #  r2d2_eids
                assert eids_elements.shape[1] == nnodes + 1, eids_elements.shape


    def element(self, eid):
        """gets a specific element of the part"""
        elem = None
        etypes_nnodes = self._etypes_nnodes()
        for etype, nnodes in etypes_nnodes:
            etype_eids = '%s_eids' % etype
            eids = getattr(self, etype_eids)  # r2d2_eids
            if eids is not None:
                ieid = np.where(eid == eids)[0]
                if len(ieid):
                    ieidi = ieid[0]
                    elems = getattr(self, etype)  # r2d2
                    elem = elems[ieid, :]
                    return etype, ieid, elem
        return None, None, None

    def check_materials(self, materials):
        """validates the materials"""
        for section in self.solid_sections:
            key = section.material
            if key in materials:
                self.log.debug('material=%r for part=%r exists' % (key, self.name))
            else:
                self.log.warning('key=%r is an invalid material' % key)

    @property
    def nelements(self):
        """Gets the total number of elements"""
        n_r2d2 = self.r2d2.shape[0] if self.r2d2 is not None else 0

        # plane strain
        n_cpe3 = self.cpe3.shape[0] if self.cpe3 is not None else 0
        n_cpe4 = self.cpe4.shape[0] if self.cpe4 is not None else 0
        n_cpe4r = self.cpe4r.shape[0] if self.cpe4r is not None else 0

        # plane stress
        n_cps3 = self.cps3.shape[0] if self.cps3 is not None else 0
        n_cps4 = self.cps4.shape[0] if self.cps4 is not None else 0
        n_cps4r = self.cps4r.shape[0] if self.cps4r is not None else 0

        n_coh2d4 = self.coh2d4.shape[0] if self.coh2d4 is not None else 0
        n_c3d10h = self.c3d10h.shape[0] if self.c3d10h is not None else 0

        n_cohax4 = self.cohax4.shape[0] if self.cohax4 is not None else 0
        n_cax3 = self.cax3.shape[0] if self.cax3 is not None else 0
        n_cax4r = self.cax4r.shape[0] if self.cax4r is not None else 0

        n_c3d8r = self.c3d8r.shape[0] if self.c3d8r is not None else 0

        neids = (n_r2d2 +
                 n_cpe3 + n_cpe4 + n_cpe4r +  # plane strain
                 n_cps3 + n_cps4 + n_cps4r +  # plane stress
                 n_coh2d4 +
                 n_c3d10h + n_cohax4 + n_cax3 + n_cax4r +
                 n_c3d8r)
        assert neids > 0, neids
        return neids

    def __repr__(self):
        """prints a summary for the part"""
        nnodes = self.nodes.shape[0]
        n_r2d2 = self.r2d2.shape[0] if self.r2d2 is not None else 0

        # plane strain
        n_cpe3 = self.cpe3.shape[0] if self.cpe3 is not None else 0
        n_cpe4 = self.cpe4.shape[0] if self.cpe4 is not None else 0
        n_cpe4r = self.cpe4r.shape[0] if self.cpe4r is not None else 0

        # plane stress
        n_cps3 = self.cps3.shape[0] if self.cps3 is not None else 0
        n_cps4 = self.cps4.shape[0] if self.cps4 is not None else 0
        n_cps4r = self.cps4r.shape[0] if self.cps4r is not None else 0

        n_coh2d4 = self.coh2d4.shape[0] if self.coh2d4 is not None else 0
        n_c3d10h = self.c3d10h.shape[0] if self.c3d10h is not None else 0

        n_cohax4 = self.cohax4.shape[0] if self.cohax4 is not None else 0
        n_cax3 = self.cax3.shape[0] if self.cax3 is not None else 0
        n_cax4r = self.cax4r.shape[0] if self.r2d2 is not None else 0

        n_c3d8r = self.c3d8r.shape[0] if self.c3d8r is not None else 0

        neids = (n_r2d2 +
                 n_cpe3 + n_cpe4 + n_cpe4r +  # plane strain
                 n_cps3 + n_cps4 + n_cps4r +  # plane stress
                 n_coh2d4 +
                 n_c3d10h + n_cohax4 + n_cax3 + n_cax4r +
                 n_c3d8r)
        assert neids == self.nelements, 'something is out of date...'
        msg = (
            f'Part(name={self.name}, nnodes={nnodes:d}, neids={neids:d},\n'
            f'     n_r2d2={n_r2d2}, n_cps3={n_cps3}, n_cpe3={n_cpe3}, '
            f'n_cpe4={n_cpe4}, n_cpe4r={n_cpe4r}, n_coh2d4={n_coh2d4},\n'
            f'     n_cohax4={n_cohax4}, n_cax3={n_cax3}, n_cax4r={n_cax4r},'
            f' n_cps4r={n_cps4r},\n'
            f'     n_c3d10h={n_c3d10h}, n_c3d8r=n_c3d8r)\n'
        )
        nsets = list(self.node_sets.keys())
        esets = list(self.element_sets.keys())
        msg += '  Node Sets: %s\n' % nsets
        msg += '  Element Sets: %s\n' % esets
        for section in self.solid_sections:
            msg += str(section) + '\n'
        return msg

    @property
    def element_types(self):
        """simplified way to access all the elements as a dictionary"""
        element_types = {}
        element_types['r2d2'] = (self.r2d2_eids, self.r2d2)

        # plane strain        element_types['cpe3'] = (self.cpe3_eids, self.cpe3)
        element_types['cpe4'] = (self.cpe4_eids, self.cpe4)
        element_types['cpe4r'] = (self.cpe4r_eids, self.cpe4r)

        # plane stress
        element_types['cps3'] = (self.cps3_eids, self.cps3)
        element_types['cps4'] = (self.cps4_eids, self.cps4)
        element_types['cps4r'] = (self.cps4r_eids, self.cps4r)

        element_types['cohax4'] = (self.cohax4_eids, self.cohax4)
        element_types['coh2d4'] = (self.coh2d4_eids, self.coh2d4)
        element_types['cax3'] = (self.cax3_eids, self.cax3)
        element_types['cax4r'] = (self.cax4r_eids, self.cax4r)
        #element_types['cps4r'] = (self.cps4r_eids, self.cps4r)
        element_types['c3d10h'] = (self.c3d10h_eids, self.c3d10h)
        return element_types

    def write(self, abq_file, is_2d=False):
        """writes a Part"""
        #name, nids, nodes, element_types, node_sets, element_sets,
         #                solid_sections, log
        abq_file.write('*Part,name=%s\n' % write_name(self.name))

        abq_file.write('*Node\n')
        if is_2d:
            for nid, node in zip(self.nids, self.nodes):
                abq_file.write('%i,\t%s,\t%s,\t%s\n' % (nid, node[0], node[1], node[2]))
        else:
            for nid, node in zip(self.nids, self.nodes):
                abq_file.write('%i,\t%s,\t%s\n' % (nid, node[0], node[1]))

        for set_name, values in sorted(self.node_sets.items()):
            write_node_set_to_file(abq_file, set_name, values)

        for elem_type, (eids, elems) in self.element_types.items():
            if eids is None:
                continue
            abq_file.write('*Element,type=%s\n' % elem_type)
            nnodes = elems.shape[1]
            fmt = '%s,\t' * (nnodes - 1) + '%s\n'
            for eid, elem in zip(eids, elems):
                abq_file.write(fmt % tuple(elem))

        for set_name, values in sorted(self.element_sets.items()):
            write_element_set_to_file(abq_file, set_name, values)
        abq_file.write('*endpart\n')

def write_name(name):
    """Abaqus has odd rules for writing words without spaces vs. with spaces"""
    return '%r' % name if ' ' in name else '%s' % name

def write_element_set_to_file(abq_file, set_name, values_array):
    """writes an element set"""
    abq_file.write('*Elset, elset=%s\n' % write_name(set_name))
    write_set_to_file(abq_file, values_array)

def write_node_set_to_file(abq_file, set_name, values_array):
    """writes a node set"""
    abq_file.write('*Nset, nset=%s\n' % write_name(set_name))
    write_set_to_file(abq_file, values_array)

def write_set_to_file(abq_file, values_array):
    """writes 16 integer values per line to a set card"""
    assert isinstance(values_array, np.ndarray), type(values_array)
    nvalues = len(values_array)
    nrows = nvalues // 16
    nleftover = nvalues % 16
    if nrows:
        values_array_square = values_array[:nrows*16].reshape(nrows, 16)
        fmt = '%i,\t' * 16 + '\n'
        fmt2 = '%i,\t' * 15 + '%i\n'
        for row in values_array_square[:-1, :]:
            abq_file.write(fmt % tuple(row))
        abq_file.write(fmt2 % tuple(values_array_square[-1, :]))

    if nleftover:
        fmt = '%i,\t' * (nleftover - 1) +  '%i\n'
        leftover = values_array[nrows*16:]
        abq_file.write(fmt % tuple(leftover))
