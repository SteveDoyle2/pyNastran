# pylint: disable=E1101,C0103
import warnings
from numpy import array


class BaseCardDeprecated(object):
    """
    Deprecated in:
      - version 0.7
    Removed in:
      - version 0.8
    """
    def rawFields(self):
        warnings.warn('deprecated; use raw_fields() instead of rawFields()',
                      DeprecationWarning, stacklevel=2)
        return self.raw_fields()

    def reprFields(self):
        warnings.warn('deprecated; use repr_fields() instead of reprFields()',
                      DeprecationWarning, stacklevel=2)
        return self.repr_fields()

    def reprCard(self):
        warnings.warn('deprecated; use print_repr_card(size, is_double) instead of reprCard()',
                      DeprecationWarning, stacklevel=2)
        return self.print_repr_card()

    def repr_card(self):
        warnings.warn('deprecated; use print_repr_card(size, is_double) instead of repr_card()',
                      DeprecationWarning, stacklevel=2)
        return self.print_repr_card()

    def printRawFields(self, size=8):
        warnings.warn('deprecated; use print_raw_card(size, is_double) instead of printRawFields(size)',
                      DeprecationWarning, stacklevel=2)
        return self.print_raw_card(size=size)

class CoordDeprecated(object):
    """defines deprecated methods for Coord"""
    def T(self):
        r"""
        Gets the 6 x 6 transformation

        .. math:: [\lambda] = [B_{ij}]

        .. math::
          [T] =
          \left[
            \begin{array}{cc}
            \lambda  & 0 \\
            0  & \lambda \\
            \end{array}
          \right]
        """
        warnings.warn('use beta_n(2) instead of T()', DeprecationWarning, stacklevel=2)
        return self.beta_n(2)

    def transformToGlobal(self, p, debug=False):
        """
        Transforms a node from the local frame to the global frame

        :param p: the xyz point in the local frame
        :param debug: debug flag (default=False; unused)

        :retval p2: the xyz point in the global frame
        :retval matrix: the transformation matrix
        """
        warnings.warn('deprecated; use p2=cid.transform_node_to_global(p);'
                      'M=cid.beta() instead of (p2,M)=cid.transformToGlobal(p)',
                      DeprecationWarning, stacklevel=2)
        if self.cid == 0:
            return p, array([[1., 0., 0.],
                             [0., 1., 0.],
                             [0., 0., 1.]], dtype='float64')

        p2 = self.transform_node_to_global(p)
        matrix = self.beta()
        return p2, matrix

    def transformVectorToGlobal(self, p):
        warnings.warn('deprecated; use transform_vector_to_global(p) instead of transformVectorToGlobal(p)',
                      DeprecationWarning, stacklevel=2)
        return self.transform_vector_to_global(p)

    def transformNodeToGlobal(self, p):
        warnings.warn('deprecated; use transform_node_to_global(p) instead of transformVectorToGlobal(p)',
                      DeprecationWarning, stacklevel=2)
        return self.transform_node_to_global(p)

    def transformToLocal(self, p, beta, debug=False):
        warnings.warn('deprecated; use transform_node_to_local(p) instead of transformToLocal(p)',
                      DeprecationWarning, stacklevel=2)
        return self.transform_node_to_local(p, debug)

    def transform_to_local(self, p, beta, debug=False):
        warnings.warn('deprecated; use transform_node_to_local(p) instead of transform_to_local(p)',
                      DeprecationWarning, stacklevel=2)
        return self.transform_node_to_local(p, debug)


class AeroDeprecated(object):
    """defines deprecated methods for Aero"""
    def IsSymmetricalXY(self):
        warnings.warn('deprecated; use is_symmetric_xy() instead of IsSymmetricalXY()',
                      DeprecationWarning, stacklevel=2)
        return self.is_symmetric_xy()

    def IsSymmetricalXZ(self):
        warnings.warn('deprecated; use is_symmetric_xz() instead of IsSymmetricalXZ()',
                      DeprecationWarning, stacklevel=2)
        return self.is_symmetric_xz()

    def IsAntiSymmetricalXY(self):
        warnings.warn('deprecated; use is_anti_symmetric_xy() instead of IsAntiSymmetricalXY()',
                      DeprecationWarning, stacklevel=2)
        return self.is_anti_symmetric_xy()

    def IsAntiSymmetricalXZ(self):
        warnings.warn('deprecated; use is_anti_symmetric_xz() instead of IsAntiSymmetricalXZ()',
                      DeprecationWarning, stacklevel=2)
        return self.is_anti_symmetric_xz()

    def EnableGroundEffect(self):
        warnings.warn('deprecated; use set_ground_effect(True) instead of EnableGroundEffect()',
                      DeprecationWarning, stacklevel=2)
        self.set_ground_effect(True)

    def DisableGroundEffect(self):
        warnings.warn('deprecated; use set_ground_effect(False) instead of DisableGroundEffect()',
                      DeprecationWarning, stacklevel=2)
        self.set_ground_effect(False)


class CAERO1Deprecated(object):
    """defines deprecated methods for CAERO1"""
    def Points(self):
        warnings.warn('deprecated; use get_points() instead of Points()',
                      DeprecationWarning, stacklevel=2)
        return self.get_points()

    def SetPoints(self, points):
        warnings.warn('deprecated; use set_points(points) instead of SetPoints(points)',
                      DeprecationWarning, stacklevel=2)
        return self.set_points(points)


class CAERO2Deprecated(object):
    """defines deprecated methods for CAERO2"""
    def Points(self):
        warnings.warn('deprecated; use get_points() instead of Points()',
                      DeprecationWarning, stacklevel=2)
        return self.get_points()

    def SetPoints(self, points):
        warnings.warn('deprecated; use set_points(points) instead of SetPoints(points)',
                      DeprecationWarning, stacklevel=2)
        return self.set_points(points)


class ElementDeprecated(object):
    """defines deprecated methods for Element"""
    def nodePositions(self, nodes=None):
        #warnings.warn('deprecated; use get_node_positions(nodes) instead of nodePositions(nodes)',
                      #DeprecationWarning, stacklevel=2)
        return self.get_node_positions(nodes=nodes)

    #def prepareNodeIDs(self, nids, allowEmptyNodes=False):
        #warnings.warn('deprecated; use prepare_node_ids(nids, allow_empty_nodes) instead of prepareNodeIDs(nids, allowEmptyNodes)',
                      #DeprecationWarning, stacklevel=2)
        #self.prepare_node_ids(nids, allow_empty_nodes=allowEmptyNodes)


class BDFMethodsDeprecated(object):
    """defines deprecated methods for BDFMethods"""
    def MassProperties(self):
        """
        .. seealso:: mass_properties
        """
        warnings.warn('deprecated; use mass_properties() instead of MassProperties()',
                      DeprecationWarning, stacklevel=2)
        return self.mass_properties()

    def Mass(self):
        """
        .. seealso:: mass
        """
        warnings.warn('deprecated; use mass_properties() instead of Mass()',
                      DeprecationWarning, stacklevel=2)
        mass, cg, I = self.mass_properties()
        return mass

    def mass(self, element_ids=None, sym_axis=None):
        """Calculates mass in the global coordinate system"""
        warnings.warn('deprecated; use mass_properties() instead of mass()',
                      DeprecationWarning, stacklevel=2)
        mass, cg, I = self.mass_properties(element_ids=element_ids,
                                           reference_point=None,
                                           sym_axis=sym_axis,
                                           num_cpus=1)
        return mass

    def resolveGrids(self, cid=0):
        """
        .. seealso:: resolve_grids
        .. deprecated: will be replaced in version 0.7 with resolve_grids
        """
        warnings.warn('deprecated; use resolve_grids(cid) instead of resolveGrids(cid)',
                      DeprecationWarning, stacklevel=2)
        return self.resolve_grids(cid)

    def unresolveGrids(self, femOld):
        """
        .. seealso:: unresolve_grids
        .. deprecated: will be replaced in version 0.7 with unresolve_grids
        """
        warnings.warn('deprecated; use unresolve_grids(femOld) instead of unresolveGrids(femOld)',
                      DeprecationWarning, stacklevel=2)
        return self.unresolve_grids(femOld)


class GetMethodsDeprecated(object):
    """defines deprecated methods for GetMethods"""
    def getElementIDsWithPID(self, pid):
        """
        Gets all the element IDs with a specific property ID

        :param pid: property ID
        :returns elementIDs: as a list

        .. deprecated:: will be removed in version 0.8

        The same functionality may be used by calling
          >>> self.getElementIDsWithPIDs([pid], mode='list')
        """
        warnings.warn('deprecated; use getElementIDsWithPIDs instead of '
                      'getElementIDsWithPID', DeprecationWarning, stacklevel=2)
        return self.getElementIDsWithPIDs([pid], mode='list')

    def Flfact(self, sid, msg):
        """
        .. deprecated:: will be removed in version 0.8
        """
        warnings.warn('deprecated; use FLFACT(sid) instead of Flfact(sid)',
                      DeprecationWarning, stacklevel=2)
        return self.FLFACT(sid, msg)

    def getNodes(self):
        """deprecated"""
        return self.get_nodes()

    def getNodeIDsWithElement(self, eid):
        """deprecated"""
        return self.get_node_ids_with_element(eid)

    def getNodeIDsWithElements(self, eids, msg=''):
        """deprecated"""
        return self.get_node_ids_with_elements(eids, msg=msg)

    #------------------------------
    def nNodes(self):
        """deprecated"""
        warnings.warn('deprecated; use nnodes instead of nNodes()',
                      DeprecationWarning, stacklevel=2)
        return self.nnodes

    def get_nnodes(self):
        """deprecated"""
        warnings.warn('deprecated; use nnodes instead of get_nnodes()',
                      DeprecationWarning, stacklevel=2)
        return len(self.nodes)

    def nodeIDs(self):
        """deprecated"""
        warnings.warn('deprecated; use node_ids instead of nodeIDs()',
                      DeprecationWarning, stacklevel=2)
        return self.get_node_ids()

    def get_node_ids(self):
        """deprecated"""
        warnings.warn('deprecated; use node_ids instead of get_node_ids()',
                      DeprecationWarning, stacklevel=2)
        return self.nodes.keys()

    #------------------------------
    def nElements(self):
        """deprecated"""
        warnings.warn('deprecated; use nelements instead of nElements()',
                      DeprecationWarning, stacklevel=2)
        return self.nelements

    def elementIDs(self):
        """deprecated"""
        warnings.warn('deprecated; use element_ids instead of elementIDs()',
                      DeprecationWarning, stacklevel=2)
        return self.element_ids

    def get_nelements(self):
        """deprecated"""
        warnings.warn('deprecated; use nelements instead of get_nelements()',
                      DeprecationWarning, stacklevel=2)
        return len(self.elements)

    def get_element_ids(self):
        """deprecated"""
        warnings.warn('deprecated; use element_ids instead of get_element_ids()',
                      DeprecationWarning, stacklevel=2)
        return self.elements.keys()


    #------------------------------
    #def nProperties(self):
        #"""deprecated"""
        #return self.nproperties

    #def elementIDs(self):
        #"""deprecated"""
        #adf
        #return self.property_ids

    #def get_nproperties(self):
        #"""deprecated"""
        #return len(self.properties)

    def get_property_ids(self):
        """deprecated"""
        warnings.warn('deprecated; use property_ids instead of get_property_ids()',
                      DeprecationWarning, stacklevel=2)
        return self.property_ids

    def propertyIDs(self):
        """deprecated"""
        warnings.warn('deprecated; use property_ids instead of propertyIDs()',
                      DeprecationWarning, stacklevel=2)
        return self.property_ids


    #------------------------------

    def get_coord_ids(self):
        warnings.warn('deprecated; use coord_ids instead of get_coord_ids()',
                      DeprecationWarning, stacklevel=2)
        return self.coord_ids


    def coordIDs(self):
        """deprecated"""
        warnings.warn('deprecated; use coord_ids instead of coordIDs()',
                      DeprecationWarning, stacklevel=2)
        return self.coord_ids

    #------------------------------

    def nCAeros(self):
        """deprecated"""
        warnings.warn('deprecated; use ncaeros instead of nCAeros()',
                      DeprecationWarning, stacklevel=2)
        return self.get_ncaeros()

    def get_ncaeros(self):
        warnings.warn('deprecated; use ncaeros instead of get_ncaeros()',
                      DeprecationWarning, stacklevel=2)
        return self.ncaeros

    #==================================================

    def getElementIDsWithPIDs(self, pids, mode='list'):
        """deprecated"""
        warnings.warn('deprecated; use get_element_ids_with_pids instead of '
                      'getElementIDsWithPIDs()',
                      DeprecationWarning, stacklevel=2)
        return self.get_element_ids_with_pids(pids, mode=mode)


    def _get_element_ids_with_pids(self, pids, mode='list'):
        if mode not in ['list', 'dict']:
            msg = "mode=%r is not supported.  Use 'list' or 'dict'\n" % mode
            raise ValueError(msg)
        if mode == 'list':
            return self.get_element_ids_list_with_pids(pids)
        else:
            return self.get_element_ids_dict_with_pids(pids)

    def getNodeIDToElementIDsMap(self):
        """deprecated"""
        warnings.warn('deprecated; use get_element_ids_with_pids instead of '
                      'getNodeIDToElementIDsMap()',
                      DeprecationWarning, stacklevel=2)
        return self.get_node_id_to_element_ids_map()

    def getPropertyIDToElementIDsMap(self):
        """deprecated"""
        warnings.warn('deprecated; use get_property_id_to_element_ids_map '
                      'instead of getPropertyIDToElementIDsMap()',
                      DeprecationWarning, stacklevel=2)
        return self.get_property_id_to_element_ids_map()

    def getMaterialIDToPropertyIDsMap(self):
        """deprecated"""
        warnings.warn('deprecated; use get_material_id_to_property_ids_map '
                      'instead of getMaterialIDToPropertyIDsMap()',
                      DeprecationWarning, stacklevel=2)
        return self.get_material_id_to_property_ids_map()

    def structuralMaterialIDs(self):
        """deprecated"""
        warnings.warn('deprecated; use get_structural_material_ids instead of '
                      'structuralMaterialIDs()',
                      DeprecationWarning, stacklevel=2)
        return self.get_structural_material_ids()

    def thermalMaterialIDs(self):
        """deprecated"""
        warnings.warn('deprecated; use get_thermal_material_ids instead of '
                      'thermalMaterialIDs()',
                      DeprecationWarning, stacklevel=2)
        return self.get_thermal_material_ids()

    def materialIDs(self):
        """deprecated"""
        warnings.warn('deprecated; use get_material_ids instead of '
                      'materialIDs()',
                      DeprecationWarning, stacklevel=2)
        return self.get_material_ids()


class ShellElementDeprecated(object):
    def Rho(self):
        """
        Returns the density
        """
        warnings.warn('deprecated; use pid.mid().rho instead of Rho()',
                      DeprecationWarning, stacklevel=2)
        #.. deprecated: will be removed in version 0.8
        #warnings.warn('replace with self.pid.mid().rho',
        #              DeprecationWarning, stacklevel=2)
        return self.pid.mid().rho


class DeprecatedCompositeShellProperty(object):
    """
    To be deprecated in:
      - Version 0.7
    To be removed in:
      - Version 0.8
    """
    def MassPerArea(self, iply='all', method='nplies'):
        return self.get_mass_per_area(iply, method)

    def Thickness(self, iply='all'):
        return self.get_thickness(iply)

    def nPlies(self):
        return self.get_nplies()

    def Nsm(self):
        return self.get_nonstructural_mass()

    def isSymmetrical(self):
        return self.is_symmetrical()

    def Rho(self, iply):
        return self.get_density(iply)

    def Theta(self, iply):
        return self.get_theta(iply)

    def sout(self, iply):
        return self.get_sout(iply)

class GridDeprecated(object):
    def Position(self, debug=False):
        return self.get_position(debug)

    def PositionWRT(self, model, cid, debug=False):
        return self.get_position_wrt(model, cid, debug=debug)

    def UpdatePosition(self, model, xyz, cid=0):
        return self.set_position(model, xyz, cid=cid)


class PointDeprecated(object):
    def Position(self, debug=False):
        return self.get_position(debug)

    def PositionWRT(self, model, cid, debug=False):
        return self.get_position_wrt(model, cid, debug)

    def UpdatePosition(self, model, xyz, cid=0):
        return self.set_position(model, xyz, cid=cid)


class SPOINTsDeprecated(object):
    def nDOF(self):
        return self.get_ndof()
