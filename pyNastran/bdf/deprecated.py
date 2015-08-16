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
        self.deprecated('rawFields()', 'raw_fields()', '0.7')
        return self.raw_fields()

    def reprFields(self):
        self.deprecated('reprFields()', 'repr_fields()', '0.7')
        return self.repr_fields()

    def reprCard(self):
        self.deprecated('reprCard()', 'print_repr_card(size, is_double)', '0.7')
        return self.print_repr_card()

    def repr_card(self):
        self.deprecated('repr_card()', 'print_repr_card(size, is_double)', '0.7')
        return self.print_repr_card()

    def printRawFields(self, size=8):
        self.deprecated('printRawFields(size)', 'print_raw_card(size, is_double)', '0.7')
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
        self.deprecated('T()', 'beta_n(2)', '0.7')
        return self.beta_n(2)

    def transformToGlobal(self, p, debug=False):
        """
        Transforms a node from the local frame to the global frame

        :param p: the xyz point in the local frame
        :param debug: debug flag (default=False; unused)

        :retval p2: the xyz point in the global frame
        :retval matrix: the transformation matrix
        """
        self.deprecated('(p2,M)=cid.transformToGlobal(p)', 'p2=cid.transform_node_to_global(p); M=cid.beta()', '0.7')
        if self.cid == 0:
            return p, array([[1., 0., 0.],
                             [0., 1., 0.],
                             [0., 0., 1.]], dtype='float64')

        p2 = self.transform_node_to_global(p)
        matrix = self.beta()
        return p2, matrix

    def transformVectorToGlobal(self, p):
        self.deprecated('transformVectorToGlobal(p)', 'transform_vector_to_global(p)', '0.7')
        return self.transform_vector_to_global(p)

    def transformNodeToGlobal(self, p):
        self.deprecated('transformNodeToGlobal(p)', 'transform_node_to_global(p)', '0.7')
        return self.transform_node_to_global(p)

    def transformToLocal(self, p, beta, debug=False):
        self.deprecated('transformToLocal(p)', 'transform_node_to_local(p)', '0.7')
        return self.transform_node_to_local(p, debug)

    def transform_to_local(self, p, beta, debug=False):
        self.deprecated('transform_to_local(p)', 'transform_node_to_local(p)', '0.7')
        return self.transform_node_to_local(p, debug)


class AeroDeprecated(object):
    """defines deprecated methods for Aero"""
    def IsSymmetricalXY(self):
        self.deprecated('IsSymmetricalXY()', 'is_symmetric_xy()', '0.7')
        return self.is_symmetric_xy()

    def IsSymmetricalXZ(self):
        self.deprecated('IsSymmetricalXZ()', 'is_symmetric_xz()', '0.7')
        return self.is_symmetric_xz()

    def IsAntiSymmetricalXY(self):
        self.deprecated('IsAntiSymmetricalXY()', 'is_anti_symmetric_xy()', '0.7')
        return self.is_anti_symmetric_xy()

    def IsAntiSymmetricalXZ(self):
        self.deprecated('IsAntiSymmetricalXZ()', 'is_anti_symmetric_xz()', '0.7')
        return self.is_anti_symmetric_xz()

    def EnableGroundEffect(self):
        self.deprecated('EnableGroundEffect()', 'set_ground_effect(True)', '0.7')
        self.set_ground_effect(True)

    def DisableGroundEffect(self):
        self.deprecated('DisableGroundEffect()', 'set_ground_effect(False)', '0.7')
        self.set_ground_effect(False)


class CAERO1Deprecated(object):
    """defines deprecated methods for CAERO1"""
    def Points(self):
        self.deprecated('Points()', 'get_points()', '0.7')
        return self.get_points()

    def SetPoints(self, points):
        self.deprecated('SetPoints(points)', 'set_points(points)', '0.7')
        return self.set_points(points)


class CAERO2Deprecated(object):
    """defines deprecated methods for CAERO2"""
    def Points(self):
        self.deprecated('Points()', 'get_points()', '0.7')
        return self.get_points()

    def SetPoints(self, points):
        self.deprecated('SetPoints(points)', 'set_points(points)', '0.7')
        return self.set_points(points)


class ElementDeprecated(object):
    """defines deprecated methods for Element"""
    def nodePositions(self, nodes=None):
        self.deprecated('nodePositions(nodes)', 'get_node_positions(nodes)', '0.8')
        return self.get_node_positions(nodes=nodes)


class BDFMethodsDeprecated(object):
    """defines deprecated methods for BDFMethods"""
    def MassProperties(self):
        """
        .. seealso:: mass_properties
        """
        self.deprecated('MassProperties()', 'mass_properties()', '0.7')
        return self.mass_properties()

    def Mass(self):
        """
        .. seealso:: mass
        """
        self.deprecated('Mass()', 'mass_properties()[0]', '0.7')
        mass, cg, I = self.mass_properties()
        return mass

    def mass(self, element_ids=None, sym_axis=None):
        """Calculates mass in the global coordinate system"""
        self.deprecated('mass()', 'mass_properties()[0]', '0.7')
        mass, cg, I = self.mass_properties(element_ids=element_ids,
                                           reference_point=None,
                                           sym_axis=sym_axis,
                                           num_cpus=1)
        return mass

    def resolveGrids(self, cid=0):
        """
        .. seealso:: resolve_grids
        """
        self.deprecated('resolveGrids(cid)', 'resolve_grids(cid)', '0.7')
        return self.resolve_grids(cid)

    def unresolveGrids(self, fem_old):
        """
        .. seealso:: unresolve_grids
        """
        self.deprecated('unresolveGrids(fem_old)', 'unresolve_grids(fem_old)', '0.7')
        return self.unresolve_grids(fem_old)


class GetMethodsDeprecated(object):
    """defines deprecated methods for GetMethods"""
    def getElementIDsWithPID(self, pid):
        """
        Gets all the element IDs with a specific property ID

        Parameters
        ----------
        pid : int
             property ID

        Returns
        -------
        elementIDs : List[int]
            as a list

        .. deprecated:: will be removed in version 0.8

        The same functionality may be used by calling
          >>> self.getElementIDsWithPIDs([pid], mode='list')
        """
        self.deprecated('getElementIDsWithPID(pid)', 'get_element_ids_list_with_pids([pid])', '0.8')
        return self.getElementIDsWithPIDs([pid], mode='list')

    def Flfact(self, sid, msg):
        """
        .. deprecated:: will be removed in version 0.8
        """
        self.deprecated('Flfact(sid)', 'FLFACT(sid)', '0.8')
        return self.FLFACT(sid, msg)

    def getNodes(self):
        """deprecated"""
        self.deprecated('getNodes', 'get_nodes()', '0.8')
        return self.get_nodes()

    def getNodeIDsWithElement(self, eid):
        """deprecated"""
        self.deprecated('getNodeIDsWithElement(eid)', 'get_node_ids_with_element(eid)', '0.8')
        return self.get_node_ids_with_element(eid)

    def getNodeIDsWithElements(self, eids, msg=''):
        """deprecated"""
        self.deprecated('getNodeIDsWithElements(eids)', 'get_node_ids_with_elements(eids)', '0.7')
        return self.get_node_ids_with_elements(eids, msg=msg)

    #------------------------------
    def nNodes(self):
        """deprecated"""
        self.deprecated('nNodes()', 'nnodes', '0.7')
        return self.nnodes

    def get_nnodes(self):
        """deprecated"""
        self.deprecated('get_nnodes()', 'nnodes', '0.7')
        return self.nnodes

    def nodeIDs(self):
        """deprecated"""
        self.deprecated('nodeIDs()', 'node_ids', '0.7')
        return self.node_ids

    def get_node_ids(self):
        """deprecated"""
        self.deprecated('get_node_ids()', 'node_ids', '0.7')
        return self.node_ids

    #------------------------------
    def nElements(self):
        """deprecated"""
        self.deprecated('nElements()', 'nelements', '0.7')
        return self.nelements

    def elementIDs(self):
        """deprecated"""
        self.deprecated('elementIDs()', 'element_ids', '0.7')
        return self.element_ids

    def get_nelements(self):
        """deprecated"""
        self.deprecated('get_nelements()', 'nelements', '0.7')
        return self.element_ids

    def get_element_ids(self):
        """deprecated"""
        self.deprecated('get_element_ids()', 'element_ids', '0.7')
        return self.element_ids

    #------------------------------
    #def nProperties(self):
        #"""deprecated"""
        #return self.nproperties

    #def elementIDs(self):
        #"""deprecated"""
        #return self.property_ids

    #def get_nproperties(self):
        #"""deprecated"""
        #return len(self.properties)

    def get_property_ids(self):
        """deprecated"""
        self.deprecated('get_property_ids()', 'property_ids', '0.7')
        return self.property_ids

    def propertyIDs(self):
        """deprecated"""
        self.deprecated('propertyIDs()', 'property_ids', '0.7')
        return self.property_ids

    #------------------------------

    def get_coord_ids(self):
        self.deprecated('get_coord_ids()', 'coord_ids', '0.7')
        return self.coord_ids


    def coordIDs(self):
        """deprecated"""
        self.deprecated('coordIDs()', 'coord_ids', '0.7')
        return self.coord_ids

    #------------------------------

    def nCAeros(self):
        """deprecated"""
        self.deprecated('nCAeros()', 'ncaeros', '0.7')
        return self.get_ncaeros()

    def get_ncaeros(self):
        self.deprecated('get_ncaeros()', 'ncaeros', '0.7')
        return self.ncaeros

    #==================================================

    def getElementIDsWithPIDs(self, pids, mode='list'):
        """deprecated"""
        self.deprecated('getElementIDsWithPIDs()', 'get_element_ids_with_pids(pids, mode)', '0.7')
        return self.get_element_ids_with_pids(pids, mode=mode)


    def _get_element_ids_with_pids(self, pids, mode='list'):
        if mode not in ['list', 'dict']:
            msg = "mode=%r is not supported.  Use 'list' or 'dict'\n" % mode
            raise ValueError(msg)
        if mode == 'list':
            self.deprecated('get_element_ids_with_pids(pids, "list")', 'get_element_ids_list_with_pids(pids)', '0.8')
            return self.get_element_ids_list_with_pids(pids)
        else:
            self.deprecated('get_element_ids_with_pids(pids, "dict")', 'get_element_ids_dict_with_pids(pids)', '0.8')
            return self.get_element_ids_dict_with_pids(pids)

    def getNodeIDToElementIDsMap(self):
        """deprecated"""
        self.deprecated('getNodeIDToElementIDsMap()', 'get_element_ids_with_pids()', '0.7')
        return self.get_node_id_to_element_ids_map()

    def getPropertyIDToElementIDsMap(self):
        """deprecated"""
        self.deprecated('getPropertyIDToElementIDsMap()', 'get_property_id_to_element_ids_map()', '0.7')
        return self.get_property_id_to_element_ids_map()

    def getMaterialIDToPropertyIDsMap(self):
        """deprecated"""
        self.deprecated('getMaterialIDToPropertyIDsMap()', 'get_material_id_to_property_ids_map()', '0.7')
        return self.get_material_id_to_property_ids_map()

    def structuralMaterialIDs(self):
        """deprecated"""
        self.deprecated('structuralMaterialIDs()', 'get_structural_material_ids()', '0.7')
        return self.get_structural_material_ids()

    def thermalMaterialIDs(self):
        """deprecated"""
        self.deprecated('thermalMaterialIDs()', 'get_thermal_material_ids()', '0.7')
        return self.get_thermal_material_ids()

    def materialIDs(self):
        """deprecated"""
        self.deprecated('materialIDs()', 'get_material_ids()', '0.7')
        return self.get_material_ids()


class ShellElementDeprecated(object):
    def Rho(self):
        """
        Returns the density
        """
        self.deprecated('Rho()', 'pid.mid().rho', '0.8')
        return self.pid.mid().rho


class DeprecatedCompositeShellProperty(object):
    """
    To be deprecated in:
      - Version 0.7
    To be removed in:
      - Version 0.8
    """
    def MassPerArea(self, iply='all', method='nplies'):
        self.deprecated('MassPerArea(iply, method)', 'get_mass_per_area(iply, method)', '0.8')
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
