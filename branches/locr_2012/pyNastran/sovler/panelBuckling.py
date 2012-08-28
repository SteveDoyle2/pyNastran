import sys
import copy

from math import acos, pi, asin
from numpy import degrees, dot
#from numpy.linalg import norm

from pyNastran.bdf.bdf import BDF, ShellElement


class PanelBuckling(object):
    def __init__(self, bdfFileName):
        """
        Preliminary work on a panel buckling code
        Assumptions:
        
        1. node connectivitity (no gaps)
        2. large change in normal vector signifies change in panel
        3. only CTRIA3, CQUAD4 elements
        """
        self.bdf = BDF(debug=True, log=None)
        self.bdf.readBDF(bdfFileName, xref=True)
        self.maxAngle = 10.

    def makePanels(self):
        self.shells = {}

        ## maps elements to edges
        Edges = {}
        for eid, element in sorted(self.bdf.elements.iteritems()):
            if isinstance(element, ShellElement):
                Edges[eid] = []
                self.shells[eid] = element

        on = True
        ## maps edges to elements
        Edges2 = {}
        Normals = {}
        for eid, element in sorted(self.shells.iteritems()):
            edges = []
            Normals[eid] = element.Normal()
            nodeIDs = element.nodeIDs()
            for i in xrange(len(nodeIDs) - 1):
                key = tuple(sorted([nodeIDs[i], nodeIDs[i + 1]]))
                Edges[eid].append(key)
                Edges2[key] = []
            key = tuple(sorted([nodeIDs[i + 1], nodeIDs[0]]))
            Edges[eid].append(key)
            Edges2[key] = []
            #print "Edges[eid=%s] = %s" %(eid,Edges[eid])
            if on:
                #print Edges[eid]
                on = False
            #print 'nodeIDs[%s] = %s' %(eid,element.nodeIDs())
            #Edges[eid] = edges
        del key, edges, element

        for eid, edges in sorted(Edges.iteritems()):
            for edge in edges:
                Edges2[edge].append(eid)
                #print "Edges2[edge=%s] = %s" %(edge,Edges2[edge])

        #print "Edges[edge=%s] = %s" %(edge,Edges2[edge]) # edges to elements

        self.buildEdgeCount(Edges, Edges2, Normals)
        #self.buildPanel2(Edges,Edges2,Normals)
        self.buildPanelMove(Edges, Edges2, Normals)

    def buildEdgeCount(self, Edges, Edges2, Normals):

        #Edges2 = {}  ## maps edges to elements
        #print "Edges[edge=%s] = %s" %(edge,Edges2[edge]) # edges to elements
        self.edgeCount = {}
        for edge, elements in sorted(Edges2.iteritems()):
            #elements = self.removeExtraElements2(elements)
            self.edgeCount[edge] = len(elements)
            print "edge=%s elements=%s edgeCount=%s" % (
                edge, elements, self.edgeCount[edge])


    def buildPanelMove(self, Edges, Edges2, Normals, eStart=None):
        self.unclaimed = Normals.keys()
        eStart = 22040
        if eStart is None:
            eStart = self.unclaimed[0]
        self.unclaimed = set(self.unclaimed)
        print "eStart = ", eStart

        #Edges2 = {}  ## maps edges to elements
        #print "Edges[edge=%s] = %s" %(edge,Edges2[edge]) # edges to elements

        patches = {}
        patch = [eStart]
        self.claimed = []
        iLevel = 0
        patch = self.buildPatch(patch, iLevel, Edges, Edges2, eStart)
        print "self.claimed = ", self.claimed

        # if a node is shared by *3* elements at the edge of a patch, it's a corner element
        # from 2, 1 is found, then 3 is found
        # from 1, 2 is found, then 3 is found
        # from 1,2,3 we can find  ?
        #
        #================
        #              ||
        #   1  |   ?   ||
        #      |       ||
        # -----*-------||
        #   2  |   3   ||
        #      |       ||

        # doesnt hand 2 element wide patches, but if we rego over all missed elements
        # a few times, we should be ok...

        print "patch = ", patch

    def buildPatch(self, patch, iLevel, Edges, Edges2, eStart):
        edges = Edges[eStart]
        print "------------"
        print "*eStart = ", eStart

        allEids = []
        for edge in edges:
            print "patch1 = ", patch
            eids = Edges2[edge]

            print "edge=%s elements=%s" % (edge, eids)
            #indexCurrentElement = eids.index(eStart)
            #eids.pop(indexCurrentElement)
            #Edges2[edge] = eids
            #allEids +=
            if self.edgeCount[edge] == 2:  # 2 unclaimed edge, 1 is symmetry
                #print "  unclaimed"
                allEids += copy.deepcopy(eids)
                self.claimed.append(allEids[0])
                self.unclaimed = self.unclaimed.difference(set([eStart]))
            elif self.edgeCount[edge] == 3:  # questionable panel
                pass

            print "patch2 = ", patch

        allEids = list(set(allEids))

        iLevel += 1
        print "allEids = ", allEids
        for eid in allEids:
            if eid in self.unclaimed:
                print "*eid=%s" % (eid)
                patch = self.buildPatch(patch, iLevel, Edges, Edges2, eid)
                print "patch = ", patch

        iLevel -= 1
        print "------------"
        return patch

            #sys.exit()
    def buildPanel2(self, Edges, Edges2, Normals):
        panelEids = {}
        self.panels = {}
        self.skippedEdges = []

        panel = []
        maxPanelID = -1

        self.eidsDone = set([])
        for eid, edges in sorted(Edges.iteritems()):
            print "eid = ", eid
            self.eidsDone = self.eidsDone.union([eid])
            # what panel are we working on
            if eid not in panelEids:
                maxPanelID += 1
                panelEids[eid] = maxPanelID
                panelID = copy.deepcopy(maxPanelID)
            else:
                panelID = panelEids[eid]
            print "panelID = ", panelID

            # what elements are nearby
            nEid = Normals[eid]
            touchingElements = []
            for edge in edges:
                print "edge=%s" % (str(edge))
                #print "eid=%s edge=%s" %(eid,edge)
                # get only the unique entries
                touchingElements = list(set(Edges2[edge]))
                touchingElements = self.removeExtraElements(
                    eid, touchingElements)

                if len(touchingElements) == 0:
                    continue
                print "touchingElements = ", touchingElements

                edgeCount = self.edgeCount[edge]
                if edgeCount == 1:
                    print "  symmetry..."
                elif edgeCount == 2:
                    for touchingElement in touchingElements:
                        if eid == touchingElement:
                            continue
                        nTouch = Normals[touchingElement]

                        angle = self.Angle(nEid, nTouch)  # degrees
                        if angle < self.maxAngle:
                            print "   *angle[%s,%s] = %g" % (eid, touchingElement, angle)
                            panelEids[touchingElement] = panelID
                        else:
                            print "    angle[%s,%s] = %g" % (eid, touchingElement, angle)

                else:
                    print "  edgeCount = ", edgeCount
                    self.skippedEdges.append(edge)

            print ""

        for edge in self.skippedEdges:
            elements = Edges2[edge]
            for eid in elements:
                if eid in panelEids:
                    del panelEids[eid]
        del Edges

        print "skippedEdges = ", self.skippedEdges

        reversedPanelEids = {}
        # i dont like initializing lists within dictionaries, messy
        for eid, panelID in panelEids.iteritems():
            reversedPanelEids[panelID] = []

        IDsToKeep = [74]
        for eid, panelID in panelEids.iteritems():
            reversedPanelEids[panelID].append(eid)

        eidsToRemove = []
        for panelID, eids in sorted(reversedPanelEids.iteritems()):
            print "panelID=%s nEids=%s eids=%s" % (panelID, len(eids), eids)
            if panelID not in IDsToKeep:
                print "****removing ", panelID
                eidsToRemove += eids

        print "eidsToRemove = ", eidsToRemove
        for eid in eidsToRemove:
            del self.bdf.elements[eid]

        #self.bdf.elements = shells2
        self.bdf.writeBDFAsPatran('fem.out.bdf')
        print "done building panels!!!"

    def Angle(self, nEid, nTouch):
        r"""
        \f[ a \cdot b = cos(\theta) * |a| * |b| \f]
        if \|a\|=\|b\|=1. (unit vectors)
        \f[ a \cdot b = cos(\theta) \f]
        """
        d = dot(nEid, nTouch)
        if d >= 1.0:  # python is stupid...
            return 0.
        elif d <= -1.:
            return pi

        #try:
        angle = degrees(acos(d))
        #except ValueError:
        #    print "nEid=%s nTouch=%s d=%s" %(nEid,nTouch,d)
        #    raise
        return angle

    def angle2(self, other, arcos=False):
        """
        max(min(..)) because acos/asin sometimes are given values slightly outside [-1, 1]
        no idea how to use this...
        """
        pass
        if arcos:
            return acos(max(min((self * other) / (self.length * other.length), 1), -1))
        else:
            return asin(max(min((self * other.ortho()) / (self.length * other.length), 1), -1))

    def removeExtraElements(self, eid, touchingElements):
        return list(set(touchingElements))

        #print "touchingElements1 = ",touchingElements
        for i, touchingElement in enumerate(touchingElements):
            #print "eidsDone = ",self.eidsDone
            if touchingElement in self.eidsDone:
                touchingElements.pop(i)
        # remove eid
        try:
            touchingElements.pop(touchingElements.index(eid))
        except:
            pass
        #print "touchingElements2 = ",touchingElements
        return touchingElements

    def removeExtraElements2(self, elements):
        return list(set(elements))


def main():
    p = PanelBuckling('aeroModel_2.bdf')
    p.makePanels()

if __name__ == '__main__':
    main()
