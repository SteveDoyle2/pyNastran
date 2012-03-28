from __future__ import print_function
import sys
from pyNastran.bdf.cards.constraints import constraintObject2

class XrefMesh(object):
    def __init__(self):
        pass

    def crossReference(self,xref=True):
        """
        links up all the cards to the cards they reference
        """
        if xref:
            self.log.info("cross Referencing...")
            #for key,e in self.elements.items():
            #    print(e)

            self.crossReference_Nodes()
            self.crossReference_Coordinates()

            self.crossReference_Elements()
            self.crossReference_Properties()
            self.crossReference_Materials()

            self.crossReference_Aero()
            self.crossReference_Constraints()
            #self.crossReference_Loads()
            #self.caseControlDeck.crossReference(self)
        ###
    
    def crossReference_Constraints(self):
        #self.spcObject.crossReference(self)  # enable to output SPCs
        #self.mpcObject.crossReference(self)  # enable to output MPCs
        pass
        
        #self.spcObject2 = constraintObject2()
        for key,spcadd in sorted(self.spcadds.items()):
            self.spcObject2.Add(spcadd)

        for key,spcs in sorted(self.spcs.items()):
            #print spcs
            for spc in spcs:
                self.spcObject2.append(spc)
        ###

        #self.mpcObject2 = constraintObject2()
        for key,mpcadd in sorted(self.mpcadds.items()):
            self.mpcObject2.Add(mpcadd)

        for key,mpcs in sorted(self.mpcs.items()):
            #print spcs
            for mpc in mpcs:
                self.mpcObject2.append(mpc)
        ###
        #self.spcObject.crossReference(self)


    def crossReference_Coordinates(self):
        """
        links up all the coordinate cards to other coordinate cards and nodes
        """
        for cid,c in self.coords.items(): # CORD2x: links the rid to coordinate systems
            c.crossReference(self)        # CORD1x: links g1,g2,g3 to grid points
        ###
        for cid,c in self.coords.items(): # CORD1x: Since the grid points were already referenced,
            c.resolveCid()                # we can now resolve the coordinate systems.
        ###                               # We couldnt do it in the previous step b/c
                                          # the grid's coordinate system might have been
                                          # unresolved

    def crossReference_Aero(self):
        """
        links up all the aero cards
        """
        for ID,caero in self.caeros.items():
            caero.crossReference(self)
        ###
        for ID,spline in self.splines.items():
            spline.crossReference(self)
        ###

    def crossReference_Nodes(self):
        """
        links the nodes to coordinate systems
        """
        gridSet = self.gridSet
        for nid,n in self.nodes.items():
            try:
                n.crossReference(self,gridSet)
            except:
                sys.stderr.write('couldnt cross reference\n%s' %(str(n)))
                raise
        ###
        if self.spoints:
            self.spointi = self.spoints.createSPOINTi()
        ###

    def crossReference_Elements(self):
        """
        links the elements to nodes, properties (and materials depending on the card)
        """
        for eid,e in self.elements.items():
            try:
                e.crossReference(self)
            except:
                sys.stderr.write('couldnt cross reference\n%s' %(str(e)))
                raise
        ###

    def crossReference_Properties(self):
        """
        links the properties to materials
        """
        for pid,p in self.properties.items():
            #print p
            try:
                p.crossReference(self)
            except:
                sys.stderr.write('couldnt cross reference\n%s' %(str(p)))
                raise
        ###

    def crossReference_Materials(self):
        """
        links the materials to materials (e.g. CREEP)
        often this is a pass statement
        """
        for mid,m in self.materials.items():
            try:
                m.crossReference(self)
            except:
                sys.stderr.write('couldnt cross reference\n%s' %(str(m)))
                raise
        ###
        for mid,m in self.materialDeps.items():
            try:
                m.crossReference(self)
            except:
                sys.stderr.write('couldnt cross reference\n%s' %(str(m)))
                raise
        ###

    def crossReference_Loads(self):
        """
        links the loads to nodes, coordinate systems, and other loads
        """
        for lid,sid in self.loads.items():
            print("lid=%s sid=%s" %(lid,sid))
            for load in sid:
                try:
                    load.crossReference(self)
                except:
                    sys.stderr.write('couldnt cross reference\n%s' %(str(load)))
                    raise
        ###

