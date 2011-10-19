import sys

from nastranParser.bdf.fieldWriter import printCard,setBlankIfDefault,setDefaultIfBlank

class BaseCard(object):

    def Is(self,typeCheck):
        if self.type==typeCheck:
            return True
        return False

    def printCard(self,fields):
        return printCard(fields)

    def setDefaultIfBlank(self,value,default):
        raise Exception('time to upgrade...')
        return setDefaultIfBlank(value,default)

    def setBlankIfDefault(self,value,default):
        return setBlankIfDefault(value,default)

    def expandThru(self,fields):
        """
        1,THRU,10
        1,3,THRU,19,15
        """
        out = []
        nFields = len(fields)
        i=0
        while(i<nFields):
            if fields[i]=='THRU':
                for j in range(fields[i-1],fields[i+1]):
                    out.append(j)
                ###
                i+=2
            else:
                out.append(i)
                i+=1
            ###
         ###
         #return list(set(out))
    
    def collapseThru(self,fields):
        """
        1,THRU,10
        1,3,THRU,19,15
        """
        fields = list(set(fields))
        
        #assumes sorted...

        pre,i = self._preCollapse(fields,dnMax=dnMax)
        mid = self._midCollapse(pre,dnMax=dnMax)
        #out = self._postCollapse(mid)
        #sys.exit()
        dnMax = 1

        out = []
        print "running post..."
        for data in mid:
            print "data = ",data
            nData = len(data)
            if nData == 1:
                out.append(data[0]) # 1 field only
            else:
                assert data[2]==1 # dn
                out += [data[0],'THRU',data[1]]
            ###
        ###
        print "dataOut = ",out
        return out
        ###

    def _midCollapse(self,preCollapse,dnMax=10000000):
        """
        input is lists of [[1,3,5,7],2]  dn=2
        dNmax = 2
        output is [1,7,2]
        """
        out = []
        print preCollapse
        for collapse in preCollapse:
            print "collapse = ",collapse
            (data,dn) = collapse
            print "data = ",data
            print "dn = ",dn
            if len(data)>1:
                if dn<=dnMax: # 1:11:2 - patran syntax
                    fields = [data[0],data[-1],dn]
                    out.append(fields)
                ###
                else: # bigger than dn
                    for field in data:
                        out.append(field)
                    ###
                ###
            else: # 1 item
                out.append([data[0]])
            ###
        return out
    
    def _preCollapse(self,fields,dnMax=10000000): # assumes sorted
        out = []
        nFields = len(fields)-1
        i=0
        while(i<nFields):
            dn = fields[i+1]-fields[i]
            print "preFields = ",fields[i:]
            (outFields,j) = self._subCollapse(fields[i:],dn,dnMax)
            print "outFields = ",outFields
            out.append([outFields,dn])
            i+=j
            ###
            #if i==nFields+1:
            #    out.append([[fields[i-1]],1])
            #    print "lastOut = ",out[-1]
            ###
        ###
        
        print "out = ",out,i
        print "--end of preCollapse"
        return (out,i)

    def _subCollapse(self,fields,dn,dnMax=10000000):
        """
        in  = [1,2,3,  7]
        out = [1,2,3]
        """
        # dn=1
        print "subIn = ",fields
        out = [fields[0]]
        nFields = len(fields)

        for i in range(1,nFields):
            dn = fields[i]-fields[i-1]
            print "i=%s field[%s]=%s fields[%s]=%s dn=%s dnMax=%s" %(i,i,fields[i],i-1,fields[i-1],dn,dnMax)
            if dn!=dnMax:
                #i+=1
                #out.append(fields[i])
                break
            out.append(fields[i])
            #if i==3:
            #    sys.exit()
        #i-=1
        print "subOut = ",out
        #i+=1
        print "iSubEnd = ",i,'\n'
        #sys.exit()
        return (out,i)
    ###
###

class Property(BaseCard):
    def __init__(self,card):
        #self.type = card[0]
        pass
        
    def __repr__(self):
        fields = [self.type,self.pid]
        return self.printCard(fields)

class Element(BaseCard):
    def __init__(self,card):
        #displayCard(card)
        self.nodes = None
        self.eid = int(card.field(1))
        self.id  = self.eid
        #self.nids = []
        pass

    def prepareNodeIDs(self,nids):
        self.nodes = []
        for nid in nids:
            self.nodes.append(int(nid))

    def Centroid(self,nodes,debug=False):
        return None

    def getNodeIDs(self):
        nids = self.nodes
        #nids = []
        #for node in nodes:
        #    nids.append(node.nid)
        #print "nids[%s] = %s" %(self.eid,nids)
        return nids

    def getNodes(self,nodes):
        """
        returns nodes...???
        """
        #print "self.type = ",self.type
        nids = self.nodes
        #print "nids = ",self.nodes
        nNodes = len(self.nodes)
        nodesList = []
        for nid in range(nNodes):
            nodesList.append(nodes[nid])
        return nodesList

    #def Normal(self,a,b):
    #    """finds the unit normal vector of 2 vectors"""
    #    return Normal(a,b)

    def CentroidTriangle(self,n1,n2,n3,debug=False):
        if debug:
            print "n1=%s \nn2=%s \nn3=%s" %(n1,n2,n3)
        centroid = (n1+n2+n3)/3.
        return centroid

    #def Area(self,a,b):
    #    return 0.5*numpy.linalg.norm(numpy.cross(a,b))

    def __repr__(self):
        fields = [self.type,self.eid,self.pid]+self.nodes
        return self.printCard(fields)

#dnMax = 2
if __name__=='__main__':
    card = BaseCard()

    """
    1,THRU,10
    1,3,THRU,19,15
    """
    card.collapseThru([1,2,3,4,5,10])
    card.collapseThru([1,3,4,5,6,17])

