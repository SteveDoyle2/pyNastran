"""
defines:
 - SolidSection
 - Material
 - Part
"""
from __future__ import print_function
from six import iteritems
import numpy as np


class SolidSection(object):
    """a SolidSection defines depth and a material"""
    def __init__(self, param_map, data_lines, log):
        self.param_map = param_map
        self.data_lines = data_lines
        self.material = param_map['material']
        assert len(data_lines) == 1., data_lines
        line0 = data_lines[0]
        assert len(line0) == 1., data_lines

        self.thickness = float(line0[0])

        for line in data_lines:
            log.info('solid - %r' % line)

    def __repr__(self):
        """prints a summary for the solid section"""
        msg = 'SolidSection(\n'
        msg += '    param_map = %r,\n' % self.param_map
        msg += '    thickness = %s,\n' % self.thickness
        msg += ')\n'
        return msg

class Material(object):
    """a Material object is a series of nodes & elements (of various types)"""
    def __init__(self, name, sections):
        self.name = name
        if 'density' in sections:
            self.density = sections['density']
        if 'depvar' in sections:
            self.depvar = sections['depvar']
        if 'user_material' in sections:
            self.user_material = sections['user_material']
        self.sections = sections

    def __repr__(self):
        """prints a summary for the material"""
        msg = 'Material(\n'
        msg += '  name=%r,\n' % self.name
        for key, value in iteritems(self.sections):
            msg += '  %r : %r,\n' % (key, value)
        msg += ')\n'
        return msg

class Part(object):
    """a Part object is a series of nodes & elements (of various types)"""
    def __init__(self, name, nids, nodes, element_types, node_sets, element_sets,
                 solid_sections, log):
        """
        creates a Part object

        Parameters
        ----------
        name : str
            the name
        element_types : Dict[element_type] : nodes
            element_type : str
                the element type

            bars:
                r2d2 : (nelements, 2) int ndarray

            shells:
                cpe3 : (nelements, 3) int ndarray
                cpe4 : (nelements, 4) int ndarray
                cpe4r : (nelements, 4) int ndarray
                coh2d4 : (nelements, 4) int ndarray
                cohax4 : (nelements, 4) int ndarray
                cax3 : (nelements, 3) int ndarray
                cax4r : (nelements, 4) int ndarray
                cps4r : (nelements, 4) int ndarray

            solids:
                c3d10h : (nelements, 10) int ndarray
        """
        self.name = name
        self.log = log
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

        # shells
        self.cpe3 = None
        self.cpe4 = None
        self.cpe4r = None
        self.coh2d4 = None
        self.cohax4 = None
        self.cax3 = None
        self.cax4r = None
        self.cps4r = None

        # solids
        self.c3d10h = None
        #-----------------------------------
        # eids
        self.r2d2_eids = None
        self.cpe3_eids = None
        self.cpe4_eids = None
        self.cpe4r_eids = None
        self.coh2d4_eids = None
        self.cohax4_eids = None
        self.cax3_eids = None
        self.cax4r_eids = None

        # rigid elements
        self.c3d10h_eids = None
        self._store_elements(element_types)

    def _store_elements(self, element_types):
        """helper method for the init"""
        if 'r2d2' in element_types: # similar to CBAR
            elements = element_types['r2d2']
            self.r2d2 = np.array(elements, dtype='int32')
            self.r2d2_eids = self.r2d2[:, 0]
            assert self.r2d2.shape[1] == 3, self.r2d2.shape

        #-------------------------------------------------------
        # shells
        if 'cpe3' in element_types: # similar to CTRIA3
            elements = element_types['cpe3']
            self.cpe3 = np.array(elements, dtype='int32')
            self.cpe3_eids = self.cpe3[:, 0]
            assert self.cpe3.shape[1] == 4, self.cpe3.shape

        if 'cpe4' in element_types: # similar to CQUAD4
            elements = element_types['cpe4']
            self.cpe4 = np.array(elements, dtype='int32')
            self.cpe4_eids = self.cpe4[:, 0]
            assert self.cpe4.shape[1] == 5, self.cpe4.shape

        if 'cpe4r' in element_types: # similar to CQUAD4
            elements = element_types['cpe4r']
            self.cpe4r = np.array(elements, dtype='int32')
            self.cpe4r_eids = self.cpe4r[:, 0]
            assert self.cpe4r.shape[1] == 5, self.cpe4r.shape

        if 'coh2d4' in element_types:
            elements = element_types['coh2d4']
            self.coh2d4 = np.array(elements, dtype='int32')
            self.coh2d4_eids = self.coh2d4[:, 0]
            assert self.coh2d4.shape[1] == 5, self.coh2d4.shape

        if 'cohax4' in element_types:
            elements = element_types['cohax4']
            self.cohax4 = np.array(elements, dtype='int32')
            self.cohax4_eids = self.cohax4[:, 0]
            assert self.cohax4.shape[1] == 5, self.cohax4.shape

        if 'cax3' in element_types:
            elements = element_types['cax3']
            self.cax3 = np.array(elements, dtype='int32')
            self.cax3_eids = self.cax3[:, 0]
            assert self.cax3.shape[1] == 4, self.cax3.shape

        if 'cax4r' in element_types:
            elements = element_types['cax4r']
            self.cax4r = np.array(elements, dtype='int32')
            self.cax4r_eids = self.cax4r[:, 0]
            assert self.cax4r.shape[1] == 5, self.cax4r.shape

        if 'cps4r' in element_types:
            elements = element_types['cps4r']
            self.cps4r = np.array(elements, dtype='int32')
            self.cps4r_eids = self.cps4r[:, 0]
            assert self.cps4r.shape[1] == 5, self.cps4r.shape

        #-------------------------------------------------------
        # solids
        if 'c3d10h' in element_types: # similar to CTRIA3
            elements = element_types['c3d10h']
            self.c3d10h = np.array(elements, dtype='int32')
            self.c3d10h_eids = self.c3d10h[:, 0]
            assert self.c3d10h.shape[1] == 11, self.c3d10h.shape

    def element(self, eid):
        """gets a specific element of the part"""
        elem = None
        # bars
        if self.r2d2_eids is not None:
            ieid = np.where(eid == self.r2d2_eids)[0]
            self.log.info('ieid_r2d2 = %s, %s' % (ieid, len(ieid)))
            if len(ieid):
                ieid = ieid[0]
                etype = 'r2d2'
                elem = self.r2d2[ieid, :]
                return etype, ieid, elem

        #-------------------------------------------------------

         # shells
        if self.cpe3_eids is not None:
            ieid = np.where(eid == self.cpe3_eids)[0]
            self.log.debug('ieid_cpe3 = %s, %s' % (ieid, len(ieid)))
            if len(ieid):
                ieid = ieid[0]
                etype = 'cpe3'
                elem = self.cpe3[ieid, :]
                return etype, ieid, elem

        if self.cpe4_eids is not None:
            ieid = np.where(eid == self.cpe4_eids)[0]
            if len(ieid):
                ieid = ieid[0]
                etype = 'cpe4'
                elem = self.cpe4[ieid, :]
                return etype, ieid, elem

        if self.cpe4r_eids is not None:
            ieid = np.where(eid == self.cpe4r_eids)[0]
            if len(ieid):
                ieid = ieid[0]
                etype = 'cpe4r'
                elem = self.cpe4r[ieid, :]
                return etype, ieid, elem

        if self.coh2d4_eids is not None:
            ieid = np.where(eid == self.coh2d4_eids)[0]
            self.log.debug('ieid_coh2d4 = %s, %s' % (ieid, len(ieid)))
            if len(ieid):
                ieid = ieid[0]
                etype = 'coh2d4'
                elem = self.coh2d4[ieid, :]
                return etype, ieid, elem
            else:
                self.log.debug('ieid = %s' % ieid)

        if self.coh2d4_eids is not None:
            ieid = np.where(eid == self.coh2d4_eids)[0]
            print('ieid_coh2d4 = %s, %s' % (ieid, len(ieid)))
            if len(ieid):
                ieid = ieid[0]
                etype = 'coh2d4'
                elem = self.coh2d4[ieid, :]
                return etype, ieid, elem
            else:
                self.log.debug('ieid = %s' % ieid)

        if self.cohax4_eids is not None:
            ieid = np.where(eid == self.cohax4_eids)[0]
            print('ieid_cohax4 = %s, %s' % (ieid, len(ieid)))
            if len(ieid):
                ieid = ieid[0]
                etype = 'cohax4'
                elem = self.cohax4[ieid, :]
                return etype, ieid, elem
            else:
                self.log.debug('ieid = %s' % ieid)
        return None, None, None

    def check_materials(self, materials):
        """validates the materials"""
        for section in self.solid_sections:
            key = section.material
            if key in materials:
                self.log.debug('material=%r for part=%r exists' % (key, self.name))
            else:
                self.log.warning('key=%r is an invalid material' % key)


    def __repr__(self):
        """prints a summary for the part"""
        nnodes = self.nodes.shape[0]
        n_r2d2 = 0
        n_cpe3 = 0
        n_cpe4 = 0
        n_cpe4r = 0
        n_coh2d4 = 0
        n_c3d10h = 0

        n_cohax4 = 0
        n_cax3 = 0
        n_cax4r = 0
        n_cps4r = 0
        if self.r2d2 is not None:
            n_r2d2 = self.r2d2.shape[0]
        if self.cpe3 is not None:
            n_cpe3 = self.cpe3.shape[0]
        if self.cpe4 is not None:
            n_cpe4 = self.cpe4.shape[0]
        if self.cpe4r is not None:
            n_cpe4r = self.cpe4r.shape[0]
        if self.coh2d4 is not None:
            n_coh2d4 = self.coh2d4.shape[0]
        if self.c3d10h is not None:
            n_c3d10h = self.c3d10h.shape[0]

        if self.cohax4 is not None:
            n_cohax4 = self.cohax4.shape[0]
        if self.cax3 is not None:
            n_cax3 = self.cax3.shape[0]
        if self.cax4r is not None:
            n_cax4r = self.cax4r.shape[0]
        if self.cps4r is not None:
            n_cps4r = self.cps4r.shape[0]

        neids = (n_r2d2 + n_cpe3 + n_cpe4 + n_cpe4r + n_coh2d4 +
                 n_c3d10h + n_cohax4 + n_cax3 + n_cax4r + n_cps4r)
        msg = (
            'Part(name=%r, nnodes=%i, neids=%i,\n'
            '     n_r2d2=%i, n_cpe3=%i, n_cpe4=%i, n_cpe4r=%i, n_coh2d4=%i,\n'
            '     n_cohax4=%i, n_cax3=%i, n_cax4r=%i, n_cps4r=%i,\n'
            '     n_c3d10h=%i)\n' % (
                self.name, nnodes, neids,
                n_r2d2, n_cpe3, n_cpe4, n_cpe4r, n_coh2d4,
                n_cohax4, n_cax3, n_cax4r, n_cps4r,
                n_c3d10h,)
        )
        for section in self.solid_sections:
            msg += str(section) + '\n'
        return msg

    @property
    def element_types(self):
        """simplified way to access all the elements as a dictionary"""
        element_types = {}
        element_types['r2d2'] = (self.r2d2_eids, self.r2d2)
        element_types['cpe3'] = (self.cpe3_eids, self.cpe3)
        element_types['cpe4'] = (self.cpe4_eids, self.cpe4)
        element_types['cpe4r'] = (self.cpe4r_eids, self.cpe4r)
        element_types['cohax4'] = (self.cohax4_eids, self.cohax4)
        element_types['coh2d4'] = (self.coh2d4_eids, self.coh2d4)
        element_types['cax3'] = (self.cax3_eids, self.cax3)
        element_types['cax4r'] = (self.cax4r_eids, self.cax4r)
        #element_types['cps4r'] = (self.cps4r_eids, self.cps4r)
        element_types['c3d10h'] = (self.c3d10h_eids, self.c3d10h)
        return element_types

    def write(self, abq_file):
        """writes a Part"""
        #name, nids, nodes, element_types, node_sets, element_sets,
         #                solid_sections, log
        abq_file.write('*Part,name=%r\n' % self.name)
        abq_file.write('*Node\n')
        for inid, node in enumerate(self.nodes):
            abq_file.write('%i,%s,%s,%s\n' % (inid, node[0], node[1], node[2]))
        for elem_type, (eids, elems) in iteritems(self.element_types):
            if eids is None:
                continue
            abq_file.write('*Element,type=%s\n' % elem_type)
            for eid, elem in zip(eids, elems):
                #print(elem)
                nids_strs = (str(nid) for nid in elem.tolist())
                #abq_file.write(str(eid) + ','.join(nids_strs) + '\n')
                abq_file.write(','.join(nids_strs) + '\n')
        abq_file.write('*endpart')
