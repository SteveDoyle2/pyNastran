class XrefMesh(object):
    def __init__(self):
        pass

    def crossReference(self):
        #print "cross Reference is a temp function"
        #for key,e in self.elements.items():
        #    print(e)
        
        #self.spcObject.crossReference(self)
        #self.caseControlDeck.crossReference(self)
        self.crossReference_Nodes()
        #self.crossReference_Elements()
        #self.crossReference_Properties()

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
            p.crossReference(self)
        ###
