import sys

from math import acos,pi
from numpy import degrees,dot
from numpy.linalg import norm

from pyNastran.bdf.bdf import BDF,ShellElement


class PanelBuckling(object):
    def __init__(self,bdfFileName):
        """
        Preliminary work on a panel buckling code
        Assumptions:
            1.  node connectivitity (no gaps)
            2.  large change in normal vector signifies change in panel
            4.  only CTRIA3, CQUAD4 elements
        """
        self.bdf = BDF(debug=True,log=None)
        self.bdf.readBDF(bdfFileName,xref=True)
        self.maxAngle = 5.

    def makePanels(self):
        self.shells = {}
        
        ## maps elements to edges
        Edges = {}
        for eid,element in sorted(self.bdf.elements.items()):
            if isinstance(element,ShellElement):
                Edges[eid] = []
                self.shells[eid] = element
            ###
        ###
        
        on = True
        ## maps edges to elements
        Edges2 = {}
        Normals = {}
        for eid,element in sorted(self.shells.items()):
            edges = []
            Normals[eid] = element.Normal()
            nodeIDs = element.nodeIDs()
            for i in range(len(nodeIDs)-1):
                key = tuple(sorted([nodeIDs[i],nodeIDs[i+1]]))
                Edges[eid].append(key)
                Edges2[key] = []
            key = tuple(sorted([nodeIDs[i+1],nodeIDs[0]]))
            Edges[eid].append(key)
            Edges2[key] = []
            print "Edges[eid=%s] = %s" %(eid,Edges[eid])
            if on:
                #print Edges[eid]
                on = False
            #print 'nodeIDs[%s] = %s' %(eid,element.nodeIDs())
            #Edges[eid] = edges
        del key,edges,element

        for eid,edges in sorted(Edges.items()):
            for edge in edges:
                Edges2[edge].append(eid)
                print "Edges2[edge=%s] = %s" %(edge,Edges2[edge])
            ###
        ###
        print "Edges[edge=%s] = %s" %(edge,Edges2[edge]) # edges to elements
        
        
        panelEids = {}
        self.panels = {}
        panel = []
        
    #def buildPanel(self,panel):
        maxPanelID = 0
        panelID = 0
        self.eidsDone = set([])
        for eid,edges in sorted(Edges.items()):
            print "eid = ",eid
            self.eidsDone = self.eidsDone.union([eid])
            # what panel are we working on
            if eid not in panelEids:
                panelEids[eid] = panelID
                maxPanelID +=1
            else:
                panelID = panelEids[eid]
            print "panelID = ",panelID
            
            # what elements are nearby
            nEid = Normals[eid]
            touchingElements = []
            for edge in edges:
                print "eid=%s edge=%s" %(eid,edge)
                # get only the unique entries
                touchingElements = list(set(Edges2[edge]))
                touchingElements = self.removeExtraElements(eid,touchingElements)

                if len(touchingElements)==0:
                    continue
                print "touchingElements = ",touchingElements
                for touchingElement in touchingElements:
                    nTouch = Normals[touchingElement]

                    angle = self.Angle(nEid,nTouch) # degrees
                    if angle<self.maxAngle:
                        panelEids[touchingElement] = panelID
                    print "    angle[%s,%s] = %g" %(eid,touchingElement,angle)
                ###
            ###
            
            print ""
            #sys.exit()
        del Edges
        
        reversedPanelEids = {}
        # i dont like initializing lists within dictionaries, messy
        for eid,panelID in panelEids.items():
            reversedPanelEids[panelID] = []

        for eid,panelID in panelEids.items():
            reversedPanelEids[panelID].append(eid)

        for panelID,eids in sorted(reversedPanelEids.items()):
            print "panelID=%s eids=%s" %(panelID,eids)
        print "done building panels!!!"

    def Angle(self,nEid,nTouch):
        """
        \f[ a \dot b = cos(\theta) * |a| * |b| \f]
        if |a|=|b|=1. (unit vectors)
        \f[ a \dot b = cos(\theta) \f]
        """
        d = dot(nEid,nTouch)
        if d>=1.0: # python is stupid...
            return 0.
        elif d<=-1.:
            return pi

        #try:
        angle = degrees(acos(d))
        #except ValueError:
        #    print "nEid=%s nTouch=%s d=%s" %(nEid,nTouch,d)
        #    raise
        return angle
    
    def angle2(self,other,arcos=False):
        """
        max(min(..)) because acos/asin sometimes are given values slightly outside [-1, 1]
        no idea how to use this...
        """
        pass
        if arcos:
            return acos(max(min((self * other)/(self.length * other.length),1),-1))           
        else:
            return math.asin(max(min((self * other.ortho())/(self.length * other.length),1),-1))

    def removeExtraElements(self,eid,touchingElements):
        #print "touchingElements1 = ",touchingElements
        for i,touchingElement in enumerate(touchingElements):
            #print "eidsDone = ",self.eidsDone
            if touchingElement in self.eidsDone:
                touchingElements.pop(i)
        # remove eid
        #touchingElements.pop(touchingElements.index(eid))
        #print "touchingElements2 = ",touchingElements
        return touchingElements


def main():
    p = PanelBuckling('aeroModel_2.bdf')
    p.makePanels()

if __name__=='__main__':
    main()


