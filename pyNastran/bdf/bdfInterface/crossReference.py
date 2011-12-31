import sys

class XrefMesh(object):
    def __init__(self):
        pass

    def crossReference(self,xref=True):
        if xref:
            #print "cross Reference is a temp function"
            #for key,e in self.elements.items():
            #    print(e)

            self.crossReference_Nodes()
            self.crossReference_Coordinates()

            self.crossReference_Elements()
            self.crossReference_Properties()
            self.crossReference_Materials()

            self.crossReference_Aero()
            #self.crossReference_Loads()
            self.spcObject.crossReference(self)
            #self.mpcObject.crossReference(self)
            #self.caseControlDeck.crossReference(self)
        ###

    def crossReference_Coordinates(self):
        for cid,c in self.coords.items(): # CORD2x: links the rid to coordinate systems
            c.crossReference(self)        # CORD1x: links g1,g2,g3 to grid points
        ###
        for cid,c in self.coords.items(): # CORD1x: Since the grid points were already referenced,
            c.resolveCid()                # we can now resolve the coordinate systems.
        ###                               # We couldnt do it in the previous step b/c
                                          # the grid's coordinate system might have been
                                          # unresolved

    def crossReference_Aero(self):
        for ID,caero in self.caeros.items():
            caero.crossReference(self)
        ###
        for ID,spline in self.splines.items():
            spline.crossReference(self)
        ###

    def crossReference_Nodes(self):
        gridSet = self.gridSet
        for nid,n in self.nodes.items():
            #print "n.cid = ",n.cid
            #coord = self.Coord(n.cid)
            #print "*",str(coord)
            try:
                n.crossReference(self,gridSet)
            except:
                sys.stderr.write('couldnt cross reference\n%s' %(str(n)))
                raise
        ###

    def crossReference_Elements(self):
        for eid,e in self.elements.items():
            #print "n.cid = ",n.cid
            #coord = self.Coord(n.cid)
            #print "*",str(coord)
            try:
                e.crossReference(self)
            except:
                sys.stderr.write('couldnt cross reference\n%s' %(str(e)))
                raise
        ###

    def crossReference_Properties(self):
        for pid,p in self.properties.items():
            #print "n.cid = ",n.cid
            #coord = self.Coord(n.cid)
            #print "*",str(coord)
            #print p
            try:
                p.crossReference(self)
            except:
                sys.stderr.write('couldnt cross reference\n%s' %(str(p)))
                raise
        ###

    def crossReference_Materials(self):
        for mid,m in self.materials.items():
            #print "n.cid = ",n.cid
            #coord = self.Coord(n.cid)
            #print "*",str(coord)
            try:
                m.crossReference(self)
            except:
                sys.stderr.write('couldnt cross reference\n%s' %(str(m)))
                raise
        ###

    def crossReference_Loads(self):
        for lid,sid in self.loads.items():
        #    #print "n.cid = ",n.cid
        #    #coord = self.Coord(n.cid)
        #    #print "*",str(coord)
            l.crossReference(self)
        ###

