import warnings

class BaseCardDeprecated(object):
    """
    Deprecated in:
      - version 0.7
    Removed in:
      - version 0.8
    """
    def rawFields(self):
        return self.raw_fields()

    def reprFields(self):
        return self.repr_fields()

    def printRawFields(self, size=8):
        return self.print_raw_fields(size=size)


class AeroDeprecated(object):
    def IsSymmetricalXY(self):
        return self.is_symmetric_xy()

    def IsSymmetricalXZ(self):
        return self.is_symmetric_xz()

    def IsAntiSymmetricalXY(self):
        return self.is_anti_symmetric_xy()

    def IsAntiSymmetricalXZ(self):
        return self.is_anti_symmetric_xz()

    def EnableGroundEffect(self):
        self.set_ground_effect(True)

    def DisableGroundEffect(self):
        self.set_ground_effect(False)


class CAERO1Deprecated(object):
    def Points(self):
        """deprecated"""
        return self.get_points()

    def SetPoints(self, points):
        """deprecated"""
        return self.set_points(points)


class CAERO2Deprecated(object):
    def Points(self):
        """deprecated"""
        return self.get_points()

    def SetPoints(self, points):
        """deprecated"""
        return self.set_points(points)


class ElementDeprecated(object):

    def nodePositions(self, nodes=None):
        """deprecated"""
        return self.get_node_positions(nodes=nodes)

    def prepareNodeIDs(self, nids, allowEmptyNodes=False):
        """deprecated"""
        self.prepare_node_ids(nids, allow_empty_nodes=allowEmptyNodes)


class BDFMethodsDeprecated(object):

    def MassProperties(self):
        """
        .. seealso:: mass_properties
        .. deprecated: will be replaced in version 0.7 with mass_properties
        """
        return self.mass_properties()

    def Mass(self):
        """
        .. seealso:: mass
        .. deprecated: will be replaced in version 0.7 with mass
        """
        return self.mass()


    def mass(self, element_ids=None, sym_axis=None):
        """Calculates mass in the global coordinate system"""
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
        return self.resolve_grids(cid)

    def unresolveGrids(self, femOld):
        """
        .. seealso:: unresolve_grids
        .. deprecated: will be replaced in version 0.7 with unresolve_grids
        """
        return self.unresolve_grids(femOld)

    #def sumForces(self):
        #"""
        #.. seealso:: sum_forces
        #.. deprecated: will be replaced in version 0.7 with sum_forces
        #"""
        #return self.sum_forces()

    #def sumMoments(self, p0):
        #"""
        #.. seealso:: sum_moments
        #.. deprecated: will be replaced in version 0.7 with sum_moments
        #"""
        #return self.sum_moments(p0)


class DeprecatedNastranMatrix(object):
    """
    Will be deprecated in:
      - v0.7
    Will be removed in:
      - v0.8
    """
    def isComplex(self):
        return self.is_complex()

    def isReal(self):
        return self.is_real()

    def isPolar(self):
        return self.is_polar()

    def getMatrix(self, isSparse=False, applySymmetry=True):
        return self.get_matrix(is_sparse=isSparse, apply_symmetry=applySymmetry)


class GetMethodsDeprecated(object):

    def getElementIDsWithPID(self, pid):
        """
        Gets all the element IDs with a specific property ID

        :param pid: property ID
        :returns elementIDs: as a list

        .. deprecated:: will be removed in version 0.8

        The same functionality may be used by calling
          >>> self.getElementIDsWithPIDs([pid], mode='list')
        """
        warnings.warn('getElementIDsWithPID has been deprecated; use '
                      'getElementIDsWithPIDs', DeprecationWarning, stacklevel=2)
        return self.getElementIDsWithPIDs([pid], mode='list')

    def Flfact(self, sid, msg):
        """
        .. deprecated:: will be removed in version 0.8
        """
        raise DeprecationWarning('use FLFACT(sid) instead of Flfact(sid)')
        return self.FLFACT(sid, msg)

    def nNodes(self):
        """deprecated"""
        return self.get_nnodes()

    def nodeIDs(self):
        """deprecated"""
        return self.get_node_ids()

    def getNodes(self):
        """deprecated"""
        return self.get_nodes()

    def getNodeIDsWithElement(self, eid):
        """deprecated"""
        return self.get_node_ids_with_element(eid)

    def getNodeIDsWithElements(self, eids, msg=''):
        """deprecated"""
        return self.get_node_ids_with_elements(self, eids, msg=msg)

    def nElements(self):
        """deprecated"""
        return self.get_nelements()

    def elementIDs(self):
        """deprecated"""
        return self.get_element_ids()

    def getElementIDsWithPIDs(self, pids, mode='list'):
        """deprecated"""
        return self.get_element_ids_with_pids(pids, mode=mode)

    def getNodeIDToElementIDsMap(self):
        """deprecated"""
        return self.get_node_id_to_element_ids_map()

    def getPropertyIDToElementIDsMap(self):
        """deprecated"""
        return self.get_property_id_to_element_ids_map()

    def getMaterialIDToPropertyIDsMap(self):
        """deprecated"""
        return self.get_material_id_to_property_ids_map()

    def propertyIDs(self):
        """deprecated"""
        return self.get_property_ids()

    def structuralMaterialIDs(self):
        """deprecated"""
        return self.get_structural_material_ids()

    def thermalMaterialIDs(self):
        """deprecated"""
        return self.get_thermal_material_ids()

    def materialIDs(self):
        """deprecated"""
        return self.get_material_ids()

    def coordIDs(self):
        """deprecated"""
        return self.get_coord_ids()

    def nCAeros(self):
        """deprecated"""
        return self.get_ncaeros()


class ShellElementDeprecated(object):
    def Rho(self):
        """
        Returns the density
        """
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

    def nDOF(self):
        return self.get_ndof()


class PointDeprecated(object):
    def Position(self, debug=False):
        return self.get_position(debug)

    def PositionWRT(self, model, cid, debug=False):
        return self.get_position_wrt(model, cid, debug)

    def UpdatePosition(self, model, xyz, cid=0):
        return self.set_position(model, xyz, cid=cid)

    def nDOF(self):
        return self.get_ndof()


class SPOINTsDeprecated(object):
    def nDOF(self):
        return self.get_ndof()
