"""
Links up the various cards in the BDF.
"""
# pylint: disable=E1101,C0103,R0902,R0904,R0914,W0611

from __future__ import print_function
from six import iteritems, itervalues
import warnings
import traceback


class XrefMesh(object):
    """
    Links up the various cards in the BDF.
    """
    def __init__(self):
        """
        The main BDF class defines all the parameters that are used.
        """
        pass

    def cross_reference(self, xref=True, xref_loads=True, xref_constraints=True):
        """
        Links up all the cards to the cards they reference
        """
        if xref:
            self.log.debug("Cross Referencing...")
            #for key,e in iteritems(self.elements):
                #print(e)

            self._cross_reference_nodes()
            self._cross_reference_coordinates()

            self._cross_reference_elements()
            self._cross_reference_properties()
            self._cross_reference_masses()
            self._cross_reference_materials()

            self._cross_reference_aero()
            if xref_constraints:
                self._cross_reference_constraints()
            if xref_loads:
                self._cross_reference_loads()
            #self.caseControlDeck.cross_reference(self)

    def _cross_reference_constraints(self):
        """
        Links the SPCADD, SPC, SPCAX, SPCD, MPCADD, MPC cards.
        """
        #self.spcObject.cross_reference(self)  # enable to output SPCs
        #self.mpcObject.cross_reference(self)  # enable to output MPCs

        #self.spcObject = constraintObject()
        for spcadd in itervalues(self.spcadds):
            self.spcObject.Add(spcadd)

        for spcs in itervalues(self.spcs):
            for spc in spcs:
                self.spcObject.append(spc)

        for mpcadd in itervalues(self.mpcadds):
            self.mpcObject.Add(mpcadd)

        for mpcs in itervalues(self.mpcs):
            for mpc in mpcs:
                self.mpcObject.append(mpc)
        #self.mpcObject = constraintObject()
        #self.spcObject.cross_reference(self)

    def _cross_reference_coordinates(self):
        """
        Links up all the coordinate cards to other coordinate cards and nodes
        """
        # CORD2x: links the rid to coordinate systems
        # CORD1x: links g1,g2,g3 to grid points
        for coord in itervalues(self.coords):
            coord.cross_reference(self)

        for coord in itervalues(self.coords):
            coord.setup()

    def _cross_reference_aero(self):
        """
        Links up all the aero cards
        """
        for caero in itervalues(self.caeros):
            caero.cross_reference(self)

        for spline in itervalues(self.splines):
            spline.cross_reference(self)

    def _cross_reference_nodes(self):
        """
        Links the nodes to coordinate systems
        """
        gridSet = self.gridSet
        for n in itervalues(self.nodes):
            try:
                n.cross_reference(self, gridSet)
            except:
                self.log.error("Couldn't cross reference GRID.\n%s" % (str(n)))
                raise

        if self.spoints:
            self.spointi = self.spoints.createSPOINTi()

    def _cross_reference_elements(self):
        """
        Links the elements to nodes, properties (and materials depending on
        the card).
        """
        for elem in itervalues(self.elements):
            try:
                elem.cross_reference(self)
            except Exception as e:
                if self._store_errors:
                    var = traceback.format_exception_only(type(e), e)
                    self._stored_errors.append(var)
                else:
                    msg = "Couldn't cross reference Element.\n%s" % str(elem)
                    self.log.error(msg)
                    raise

    def _cross_reference_masses(self):
        """
        Links the mass to nodes, properties (and materials depending on
        the card).
        """
        for mass in itervalues(self.masses):
            try:
                mass.cross_reference(self)
            except Exception as e:
                if self._store_errors:
                    var = traceback.format_exception_only(type(e), e)
                    self._stored_errors.append(var)
                else:
                    msg = "Couldn't cross reference Mass.\n%s" % str(mass)
                    self.log.error(msg)
                    raise

        for prop in itervalues(self.properties_mass):
            try:
                prop.cross_reference(self)
            except Exception as e:
                if self._store_errors:
                    var = traceback.format_exception_only(type(e), e)
                    self._stored_errors.append(var)
                else:
                    msg = "Couldn't cross reference PropertyMass.\n%s" % (str(prop))
                    self.log.error(msg)
                    raise

    def _cross_reference_properties(self):
        """
        Links the properties to materials
        """
        for prop in itervalues(self.properties):
            try:
                prop.cross_reference(self)
            except Exception as e:
                if self._store_errors:
                    var = traceback.format_exception_only(type(e), e)
                    self._stored_errors.append(var)
                else:
                    msg = "Couldn't cross reference Property.\n%s" % (str(prop))
                    self.log.error(msg)
                    raise

    def _cross_reference_materials(self):
        """
        Links the materials to materials (e.g. MAT1, CREEP)
        often this is a pass statement
        """
        for mat in itervalues(self.materials):  # MAT1
            try:
                mat.cross_reference(self)
            except Exception as e:
                if self._store_errors:
                    var = traceback.format_exception_only(type(e), e)
                    self._stored_errors.append(var)
                else:
                    msg = "Couldn't cross reference Material\n%s" % (str(mat))
                    self.log.error(msg)
                    raise

        # CREEP - depends on MAT1
        data = [self.MATS1, self.MATS3, self.MATS8,
                self.MATT1, self.MATT2, self.MATT3, self.MATT4, self.MATT5,
                self.MATT8, self.MATT9]
        for material_deps in data:
            for mat in itervalues(material_deps):
                try:
                    mat.cross_reference(self)
                except Exception as e:
                    if self._store_errors:
                        var = traceback.format_exception_only(type(e), e)
                        self._stored_errors.append(var)
                    else:
                        msg = "Couldn't cross reference Material\n%s" % (str(mat))
                        self.log.error(msg)
                        raise

    def _cross_reference_loads(self):
        """
        Links the loads to nodes, coordinate systems, and other loads.
        """
        for (lid, sid) in iteritems(self.loads):
            #self.log.debug("lid=%s sid=%s" %(lid,sid))
            for load in sid:
                try:
                    load.cross_reference(self)
                except Exception as e:
                    if self._store_errors:
                        var = traceback.format_exception_only(type(e), e)
                        self._stored_errors.append(var)
                    else:
                        self.log.error("lid=%s sid=%s" % (lid, sid))
                        msg = "Couldn't cross reference Load\n%s" % (str(load))
                        self.log.error(msg)
                        raise