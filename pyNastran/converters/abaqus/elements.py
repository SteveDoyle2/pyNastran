"""
General 3D solids
C3D4 (4-node linear tetrahedral element)
C3D6 (6-node linear triangular prism element)
C3D8 (3D 8-node linear isoparametric element)
C3D8R (the C3D8 element with reduced integration)
C3D10 (10-node quadratic tetrahedral element)
C3D15 (15-node quadratic triangular prism element)
C3D20 (3D 20-node quadratic isoparametric element)
C3D20R (the C3D20 element with reduced integration)
C3D20RI (incompressible C3D20 element with reduced integration)
``ABAQUS'' 3D solids for heat transfer (names are provided for compatibility)
DC3D4: identical to C3D4
DC3D6: identical to C3D6
DC3D8: identical to C3D8
DC3D10: identical to C3D10
DC3D15: identical to C3D15
DC3D20: identical to C3D20
Shell elements
S6 (6-node triangular shell element)
S8 (8-node quadratic shell element)
S8R (the S8 element with reduced integration)
Plane stress elements
  CPS6 (6-node triangular plane stress element)
  CPS8 (8-node quadratic plane stress element)
  CPS8R (the CPS8 element with reduced integration)
Plane strain elements
  CPE6 (6-node triangular plane strain element)
  CPE8 (8-node quadratic plane strain element)
  CPE8R (the CPS8 element with reduced integration)
Axisymmetric elements
  CAX6 (6-node triangular axisymmetric element)
  CAX8 (8-node quadratic axisymmetric element)
  CAX8R (the CAX8 element with reduced integration)
Beam elements
  B32 (3-node beam element)
  B32R (the B32 element with reduced integration)
Special elements
  D (3-node network element)
  GAPUNI (2-node unidirectional gap element)

http://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node194.html
"""
from typing import Dict
import numpy as np
from cpylog import SimpleLogger

allowed_element_types = [
    'r2d2', 'conn2d2', 'springa',
    'cpe3', 'cpe4', 'cpe4r', 'cpe8r',
    'cps3', 'cps4', 'cps4r', 'cps8r',

    'coh2d4', 'c3d10h', 'cohax4',
    'cax3', 'cax4r', 'mass', 'rotaryi', 't2d2', 'c3d8r',

    # 6/8 plates
    's8r',

    # solid
    'c3d4', 'c3d10',

    # lines
    'b31h', #2-node linear beam
]

class Elements:
    """a Part object is a series of nodes & elements (of various types)"""
    def __init__(self, element_types: Dict[str, np.ndarray], log: SimpleLogger):
        """
        creates a Part object

        Parameters
        ----------
        element_types : Dict[element_type] : node_ids
            element_type : str
                the element type
            bars:
                r2d2 : (nelements, 2) int ndarray
            beams:
                b31h : (nelements, 3) int ndarray - 2 nodes and g0
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
                c3d8r : (nelements, 8) int ndarray
                c3d10 : (nelements, 10) int ndarray
                c3d10h : (nelements, 10) int ndarray

        """
        self.log = log
        # bars
        self.r2d2 = None
        self.b31h = None

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

        # 6/8 shells
        self.s8r = None

        # solids
        self.c3d4 = None

        # solids
        self.c3d10 = None
        self.c3d10h = None
        self.c3d8r = None
        #-----------------------------------
        # eids
        self.r2d2_eids = None
        self.b31h_eids = None

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

        # 6/8 shells
        self.s8r_eids = None

        # solids
        self.c3d4_eids = None
        self.c3d10_eids = None

        # rigid elements
        self.c3d10h_eids = None
        self.c3d8r_eids = None
        self._store_elements(element_types)

    def _etypes_nnodes(self):
        """internal helper method"""
        etypes_nnodes = [
            ('r2d2', 2),  #  similar to a CBAR
            ('b31h', 3),  #  similar to a CBEAM?  TODO: why is this 3 nodes?

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

            ('s8r', 8), # CQUAD8

            #  solids
            ('c3d4', 4),  #  tet4
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
                if len(elements) == 0:
                    continue
                eids_elements = np.array(elements, dtype='int32')
                setattr(self, etype, eids_elements[:, 1:])  # r2d2
                setattr(self, etype_eids, eids_elements[:, 0]) #  r2d2_eids
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
                    #ieidi = ieid[0]
                    elems = getattr(self, etype)  # r2d2
                    elem = elems[ieid, :]
                    return etype, ieid, elem
        return None, None, None

    @property
    def nelements(self):
        """Gets the total number of elements"""
        n_r2d2 = self.r2d2.shape[0] if self.r2d2 is not None else 0  #  bar
        n_b31h = self.b31h.shape[0] if self.b31h is not None else 0  #  beam

        # plane strain
        n_cpe3 = self.cpe3.shape[0] if self.cpe3 is not None else 0
        n_cpe4 = self.cpe4.shape[0] if self.cpe4 is not None else 0
        n_cpe4r = self.cpe4r.shape[0] if self.cpe4r is not None else 0

        # plane stress
        n_cps3 = self.cps3.shape[0] if self.cps3 is not None else 0
        n_cps4 = self.cps4.shape[0] if self.cps4 is not None else 0
        n_cps4r = self.cps4r.shape[0] if self.cps4r is not None else 0

        n_coh2d4 = self.coh2d4.shape[0] if self.coh2d4 is not None else 0

        n_cohax4 = self.cohax4.shape[0] if self.cohax4 is not None else 0
        n_cax3 = self.cax3.shape[0] if self.cax3 is not None else 0
        n_cax4r = self.cax4r.shape[0] if self.cax4r is not None else 0

        # 6/8
        n_s8r = self.s8r.shape[0] if self.s8r is not None else 0
        # solids
        n_c3d8r = self.c3d8r.shape[0] if self.c3d8r is not None else 0
        n_c3d4 = self.c3d4.shape[0] if self.c3d4 is not None else 0
        n_c3d10h = self.c3d10h.shape[0] if self.c3d10h is not None else 0

        neids = (n_r2d2 + n_b31h +
                 n_cpe3 + n_cpe4 + n_cpe4r +  # plane strain
                 n_cps3 + n_cps4 + n_cps4r +  # plane stress
                 n_coh2d4 +
                 n_cohax4 + n_cax3 + n_cax4r +
                 #6 / 8 shells
                 n_s8r +
                 # solids
                 n_c3d8r + n_c3d4 + n_c3d10h)
        self.log.info(f'neids = {neids}')
        assert neids > 0, neids
        return neids

    def __repr__(self):
        """prints a summary for the part"""
        n_r2d2 = self.r2d2.shape[0] if self.r2d2 is not None else 0
        n_b31h = self.b31h.shape[0] if self.b31h is not None else 0

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

        # 6/8 shells
        n_s8r = self.s8r.shape[0] if self.s8r is not None else 0

        # solids
        n_c3d8r = self.c3d8r.shape[0] if self.c3d8r is not None else 0
        n_c3d4 = self.c3d4.shape[0] if self.c3d4 is not None else 0

        neids = (n_r2d2 + n_b31h +
                 n_cpe3 + n_cpe4 + n_cpe4r +  # plane strain
                 n_cps3 + n_cps4 + n_cps4r +  # plane stress
                 n_coh2d4 +
                 n_c3d10h + n_cohax4 + n_cax3 + n_cax4r +
                 n_s8r +
                 n_c3d8r + n_c3d4)
        assert neids == self.nelements, 'something is out of date...'
        msg = (
            f'Element(neids={neids:d},\n'
            f'        n_r2d2={n_r2d2}, n_b31h={n_b31h},\n'
            f'        n_cps3={n_cps3}, n_cpe3={n_cpe3},\n'
            f'        n_cpe4={n_cpe4}, n_cpe4r={n_cpe4r}, n_coh2d4={n_coh2d4},\n'
            f'        n_cohax4={n_cohax4}, n_cax3={n_cax3}, n_cax4r={n_cax4r},\n'
            f'        n_cps4r={n_cps4r}, n_c3d10h={n_c3d10h}, n_c3d8r=n_c3d8r)\n'
        )
        return msg

    @property
    def element_types(self):
        """simplified way to access all the elements as a dictionary"""
        element_types = {}
        element_types['r2d2'] = (self.r2d2_eids, self.r2d2)
        element_types['b31h'] = (self.b31h_eids, self.b31h)

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

        # 6/8 shells
        element_types['s8r'] = (self.s8r_eids, self.s8r)

        # solids
        element_types['c3d10'] = (self.c3d10_eids, self.c3d10)
        element_types['c3d10h'] = (self.c3d10h_eids, self.c3d10h)
        element_types['c3d4'] = (self.c3d4_eids, self.c3d4)
        return element_types

    def write(self, abq_file):
        """writes a Part"""
        for elem_type, (eids, elems) in self.element_types.items():
            if eids is None:
                continue
            neids, nnodes = elems.shape
            eids2 = eids.reshape(neids, 1)
            eids_elems = np.hstack([eids2, elems])

            abq_file.write('*Element,type=%s\n' % elem_type)
            fmt = '%d,\t' * (nnodes) + '%d\n'
            for eid_elem in eids_elems:
                abq_file.write(fmt % tuple(eid_elem))
