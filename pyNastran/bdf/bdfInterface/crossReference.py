# pylint: disable=E1101,C0103,R0902,R0904,R0914,W0611

from __future__ import print_function
import sys
import warnings
#from pyNastran.bdf.cards.constraints import constraintObject2


class XrefMesh(object):
    """
    Links up the various cards in the BDF.
    """
    def __init__(self):
        """
        The main BDF class defines all the parameters that are used.
        """
        pass

    def crossReference(self, xref=True):
        """
        @see cross_reference
        """
        warnings.warn('crossReference is deprecated; use cross_reference')
        self.cross_reference(xref)

    def cross_reference(self, xref=True):
        """
        Links up all the cards to the cards they reference
        """
        if xref:
            self.log.debug("Cross Referencing...")
            #for key,e in self.elements.iteritems():
                #print(e)

            self._cross_reference_nodes()
            self._cross_reference_coordinates()

            self._cross_reference_elements()
            self._cross_reference_properties()
            self._cross_reference_materials()

            self._cross_reference_aero()
            self._cross_reference_constraints()
            self._cross_reference_loads()
            #self.caseControlDeck.cross_reference(self)

    def _cross_reference_constraints(self):
        """
        Links the SPCADD, SPC, SPCAX, SPCD, MPCADD, MPC cards.
        """
        #self.spcObject.cross_reference(self)  # enable to output SPCs
        #self.mpcObject.cross_reference(self)  # enable to output MPCs

        #self.spcObject2 = constraintObject2()
        for spcadd in self.spcadds.itervalues():
            self.spcObject2.Add(spcadd)

        for spcs in self.spcs.itervalues():
            for spc in spcs:
                self.spcObject2.append(spc)

        for mpcadd in self.mpcadds.itervalues():
            self.mpcObject2.Add(mpcadd)

        for mpcs in self.mpcs.itervalues():
            for mpc in mpcs:
                self.mpcObject2.append(mpc)
        #self.mpcObject2 = constraintObject2()
        #self.spcObject.cross_reference(self)

    def _cross_reference_coordinates(self):
        """
        Links up all the coordinate cards to other coordinate cards and nodes
        """
        # CORD2x: links the rid to coordinate systems
        # CORD1x: links g1,g2,g3 to grid points
        for coord in self.coords.itervalues():
            coord.cross_reference(self)

        # CORD1x: Since the grid points were already referenced,
        # we can now resolve the coordinate systems.
        # We couldnt do it in the previous step b/c
        # the grid's coordinate system might have been
        # unresolved
        for coord in self.coords.itervalues():
            coord.resolveCid()

    def _cross_reference_aero(self):
        """
        Links up all the aero cards
        """
        for caero in self.caeros.itervalues():
            caero.cross_reference(self)

        for spline in self.splines.itervalues():
            spline.cross_reference(self)

    def _cross_reference_nodes(self):
        """
        Links the nodes to coordinate systems
        """
        gridSet = self.gridSet
        for n in self.nodes.itervalues():
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
        for elem in self.elements.itervalues():
            try:
                elem.cross_reference(self)
            except:
                msg = "Couldn't cross reference Element.\n%s" % (str(elem))
                self.log.error(msg)
                raise

    def _cross_reference_properties(self):
        """
        Links the properties to materials
        """
        for prop in self.properties.itervalues():
            try:
                prop.cross_reference(self)
            except:
                msg = "Couldn't cross reference Property.\n%s" % (str(prop))
                self.log.error(msg)
                raise

    def _cross_reference_materials(self):
        """
        Links the materials to materials (e.g. MAT1, CREEP)
        often this is a pass statement
        """
        for mat in self.materials.itervalues():  # MAT1
            try:
                mat.cross_reference(self)
            except:
                msg = "Couldn't cross reference Material\n%s" % (str(mat))
                self.log.error(msg)
                raise

        for mat in self.materialDeps.itervalues():  # CREEP - depends on MAT1
            try:
                mat.cross_reference(self)
            except:
                msg = "Couldn't cross reference Material\n%s" % (str(mat))
                self.log.error(msg)
                raise

    def _cross_reference_loads(self):
        """
        Links the loads to nodes, coordinate systems, and other loads.
        """
        for (lid, sid) in self.loads.iteritems():
            #self.log.debug("lid=%s sid=%s" %(lid,sid))
            for load in sid:
                try:
                    load.cross_reference(self)
                except:
                    self.log.error("lid=%s sid=%s" % (lid, sid))
                    msg = "Couldn't cross reference Load\n%s" % (str(load))
                    self.log.error(msg)
                    raise
