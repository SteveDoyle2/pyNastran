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
            #self.caseControlDeck.crossReference(self)
        ###

    def crossReference_Coordinates(self):
        for cid,c in self.coords.items():
            c.crossReference(self)
        ###
        for cid,c in self.coords.items():
            c.resolveCid()
        ###

    def crossReference_Aero(self):
        pass

    def crossReference_Nodes(self):
        gridSet = self.gridSet
        for nid,n in self.nodes.items():
            #print "n.cid = ",n.cid
            #coord = self.Coord(n.cid)
            #print "*",str(coord)
            n.crossReference(self,gridSet)
        ###

    def crossReference_Elements(self):
        for eid,e in self.elements.items():
            #print "n.cid = ",n.cid
            #coord = self.Coord(n.cid)
            #print "*",str(coord)
            e.crossReference(self)
        ###

    def crossReference_Properties(self):
        for pid,p in self.properties.items():
            #print "n.cid = ",n.cid
            #coord = self.Coord(n.cid)
            #print "*",str(coord)
            #print p
            p.crossReference(self)
        ###

    def crossReference_Materials(self):
        for mid,m in self.materials.items():
            #print "n.cid = ",n.cid
            #coord = self.Coord(n.cid)
            #print "*",str(coord)
            m.crossReference(self)
        ###

    def crossReference_Loads(self):
        for lid,sid in self.loads.items():
        #    #print "n.cid = ",n.cid
        #    #coord = self.Coord(n.cid)
        #    #print "*",str(coord)
            l.crossReference(self)
        ###

