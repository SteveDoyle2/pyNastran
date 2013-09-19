from numpy import array,cross

# my code
from mathFunctions import Centroid,Triangle_AreaCentroidNormal,AreaNormal,ListPrint
from cart3d_reader import Cart3DReader
from FEM_Mesh import FEM_Mesh
from bdf import printCard

from logger import dummyLogger
loggerObj = dummyLogger()
log = loggerObj.startLog('debug') # or info

#------------------------------------------------------------------

class Model(object):
    def __init__(self):
        pass
    
    def getNodeIDLocations(self,nIDs):
        nodes = []
        for nid in nIDs:
            node = self.Node(nid)
            nodes.append(node)
        return nodes

    def Node(self,nid):
        raise Exception('overwrite this method...')

class StructuralModel(Model):
    """
    this class gets standard mesh parameters
    """
    def __init__(self,fem,pids,debug=False):
        Model.__init__(self)

        self.debug = debug
        if self.debug:
            print "*StructuralModel.init"
        self.fem = fem
        #nodes    = fem.getNodes()
        #elements = fem.getElements()
        nodeIDs    = fem.nodeIDs()
        elementIDs = fem.elementIDs()

        self.nNodes = len(nodeIDs)
        self.nElements = len(elementIDs)
        self.pids = pids
        #self.points   = points
        #self.elements = elements
        if self.debug:
            log().debug("***StructuralModel.init")

    def NodeIDs(self):
        return self.fem.nodeIDs()

    def getElementProperties(self,eid):
        """Returns area, centroid, normal"""
        e = self.fem.Element(eid)
        nodes = self.getElementNodes(eid)
        (area,centroid,normal) = e.AreaCentroidNormal(nodes)
        return (area,centroid,normal)

    def ElementIDs(self):
        #(elements,eids) = getElementsWithPIDs(self,properties)
        eids = self.fem.elementIDs()
        return eids
    
    def getElementNodeIDs(self,eid):
        e = self.fem.Element(eid)
        nIDs = e.nodes
        return nIDs

    def Node(self,nid):
        node = self.fem.Node(nid)
        return node.xyz

    def Element(self,eid):
        return self.fem.Element(eid)
    
    def etype(self,eid):
        return self.fem.Element(eid).type

    def getElementIDsWithPIDs(self):
        return self.fem.getElementIDsWithPIDs(self.pids)
        
    def getElementNodes(self,eid):
        nIDs = self.getElementNodeIDs(eid)
        #print "nIDs = ",nIDs
        nodes = []
        if nIDs:
            for nid in nIDs:
                node = self.Node(nid)
                nodes.append(node)
            return nodes
        else:
            return None

    def Centroid(self,eid):
        e = self.fem.Element(eid)
        nodes = self.getElementNodes(eid)
        centroid = e.Centroid(nodes)
        
        if isinstance(centroid,float):
            print "nodes = ",nodes
            print "*centroid[%s] = %s" %(eid,centroid)
            centroid = e.Centroid(nodes,debug=True)
            raise Exception('bad length')
        return centroid

    def Centroid_Area(self,eid,nodes):
        e = self.fem.Element(eid)
        raise Exception('not implemented')
        return (centroid,area)

    def Area(self,eid):
        nodes = self.fem.nodes
        e = self.fem.Element(eid)
        return e.Area(nodes)

    def Normal(self,eid):
        nodes = self.fem.nodes
        return self.fem.Normal(eid,nodes)

    def getElements(self): # dict
        return self.elements

    def Properties(self,eid):
        e = self.fem.Element(eid)
        (area,centroid,normal) = e.getAreaCentroidNormal()
        return (normal,centroid)

    def writeLoad(self,bdf,loadCase,nid,Fx,Fy,Fz,comment=''):
        """
        This function takes a:
           load case
           node ID
           Force
        and writes out a properly formatted card
        """
        cid = 0
        scaleFactor = 1.
        card = ['FORCE',loadCase,nid,cid,scaleFactor,Fx,Fy,Fz]
        #comment += " card=%s" %(card)
        #out = printCard(card)[:-1]+  '   $ %s\n' %(comment)
        out = printCard(card)
        bdf.write(out)

#------------------------------------------------------------------

class AeroModel(Model):
    def __init__(self,nodes,elements,Cp,pInf,qInf):
        """
        pInf - psi
        qInf - psi
        """
        Model.__init__(self)


        #Bref = 623.179569341
        #Cref = 623.179569341
        #Lref = 623.179569341
        #Sref = 1582876.10613
        #Xref = 268.217444245


        self.pInf = pInf # 499.3/144.
        self.qInf = qInf # 0.825=mach

        self.Sref = 1582876. # inches, bad value...
        self.Lref = 623.  # inch, bad value...
        self.xref = 268.

        self._nNodes    = len(nodes)
        self._nElements = len(elements)
        self.nodes      = nodes
        self.elements   = elements
        Cp = self.prepareCps(Cp)  # convert nodal Cp to centroidal Cp

        self.prepareCentroidAreaNormals()
        self.getMoments(Cp) # centroidal Cp

    def getMoments(self,Cp):
        pInf = self.pInf; qInf=self.qInf
        momentCenter = array([self.xref,0.,0.])
        sumMoments = array([0.,0.,0.])
        sumForces  = array([0.,0.,0.])

        for (key,cp) in self.Cps.items(): # centroidal based Cp
            area = self.areas[key]
            centroid = self.centroids[key]
            normal = self.normals[key]
            p = cp*qInf+pInf
            F = area*normal*p  # negative sign is b/c the normals are flipped...
            r = momentCenter-centroid

            sumForces  += F
            sumMoments += cross(r,F)
            #break
        ###
        log().info("pInf=%s [psi]; qInf= %s [psi]" %(pInf,qInf))
        
        log().info("sumForcesCFD  [lb]    = %s" %(ListPrint(sumForces)))
        log().info("sumMomentsCFD [ft-lb] = %s" %(ListPrint(sumMoments/12.)))
        Cf  = sumForces/(self.Sref*qInf)
        Cm = sumMoments/(self.Sref*qInf*self.Lref)*12.
        log().info("Cf = %s" %(ListPrint(Cf)))
        log().info("Cm = %s" %(ListPrint(Cm)))
        return (sumForces,sumMoments/12.)

    def prepareCentroidAreaNormals(self):
        self.centroids = {}
        self.areas     = {}
        self.normals   = {}
        for key,element in sorted(self.elements.items()):
            n1,n2,n3 = element
            n1 = self.nodes[n1]
            n2 = self.nodes[n2]
            n3 = self.nodes[n3]
            out = Triangle_AreaCentroidNormal([n1,n2,n3])
            #print "out = ",out
            (area,centroid,normal) = out
            if centroid is not None:
                assert len(centroid)==3,"eid=%s centroid=%s n1=%s n2=%s n3=%s" %(key,centroid,n1,n2,n3)
            
            self.areas[key]     = area
            self.centroids[key] = centroid
            self.normals[key]   = normal
        ###

    def prepareCps(self,Cp):
        """
        converts Cp applied to the node -> Cp applied on the element centroid
        """
        #self.Cps = Cp

        CpDict = {}
        for eid in self.ElementIDs():
            #print "eid = ",eid
            (n1,n2,n3) = self.getElementNodeIDs(eid)
            cp1 = Cp[n1]
            cp2 = Cp[n2]
            cp3 = Cp[n3]
            cp = (cp1+cp2+cp3)/3.
            CpDict[eid]=cp
        self.Cps = CpDict
        return CpDict

    def nNodes(self):
        return self._nNodes

    def Centroid(self,eid):
        return self.centroids[eid]

        nodes = self.getElementNodes(eid)
        return Centroid(*nodes)

    def nElements(self):
        return self._nElements
    
    def nCps(self):
        return self._nElements

    def Node(self,nid):
        return self.nodes[nid]

    def Element(self,eid):
        try:
            element = self.elements[eid]
        except IndexError:
            print "eid=%s len(elements)=%s" %(eid,len(self.elements))
            raise
        except TypeError:
            print "eid=%s len(elements)=%s" %(eid,len(self.elements))
            raise
        except KeyError:
            print "eid=%s len(elements)=%s" %(eid,len(self.elements))
            raise
        return element

    def Area(self,eid):
        return self.areas[eid]

        nodes = self.getElementNodes(eid)
        #print "nodes[%s]=%s" %(eid,nodes)
        (area,normal) = AreaNormal(nodes)
        return area

    def Normal(self,eid):
        return self.normals[eid]

        nodes = self.getElementNodes(eid)
        (area,normal) = AreaNormal(nodes)
        return normal

    def AreaNormalCentroid(self,eid):
        area     = self.areas[eid]
        normal   = self.normals[eid]
        centroid = self.centroids[eid]

        #nodes = self.getElementNodes(eid)
        #(area,normal,centroid) = Triangle_AreaNormalCentroid(nodes)
        return (area,normal,centroid)

    def NodeIDs(self):
        nNodes = self.nNodes()
        #print "nNodes = ",nNodes
        return range(nNodes)

    #def getElementIDsWithPIDs(self):
    #    return self.ElementIDs()

    def ElementIDs(self):
        nElements = self.nElements()
        #print "nElements = ",nElements
        return self.elements.keys()
        
    def Cp(self,eid):
        #print "Cp eid=%s" %eid
        #cp1 = self.Cps[n1]
        #cp2 = self.Cps[n2]
        #cp3 = self.Cps[n3]
        #cp = (cp1+cp2+cp3)/3.
        cp = self.Cps[eid]
        #print "eid=%s Cp=%s" %(eid,cp)
        return (cp)
        
        
        #try:
        #    cp = self.Cps[eid-1]
        #except:
        #    print "***error...Cp eid=%s" %eid
        #    raise
        #return cp

    #def getElementCentroid(self,eid):
    #    return None
    #def getElementArea(self,eid):
    #    return None
    #def getElementNormal(self,eid):
    #    return None

    def getElementNodeIDs(self,eid):
        #print "eid2 = ",eid
        e = self.Element(eid)
        (nid1,nid2,nid3) = e
        return (nid1,nid2,nid3)

    def getElementNodes(self,eid):
        (nid1,nid2,nid3) = self.getElementNodeIDs(eid)
        n1 = self.Node(nid1)
        n2 = self.Node(nid2)
        n3 = self.Node(nid3)
        #print "nids[%s] = %s %s %s" %(eid,nid1,nid2,nid3)
        #print "n[%s]=%s n[%s]=%s n[%s]=%s\n" %(nid1,n1,nid2,n2,nid3,n3)
        return (n1,n2,n3)
     
    def getElementProperties(self,eid):
        """
        Returns area, centroid, normal
        """
        nodes = self.getElementNodes(eid)
        (area,centroid,normal) = Triangle_AreaCentroidNormal(nodes)
        return (area,centroid,normal)

#------------------------------------------------------------------

if __name__=='__main__':
    pass
