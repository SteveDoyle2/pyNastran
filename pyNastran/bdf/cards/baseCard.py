import sys

import pyNastran
import pyNastran.bdf

from pyNastran.bdf.fieldWriter import printCard,setBlankIfDefault,setDefaultIfBlank
from pyNastran.bdf.BDF_Card import BDF_Card

class BaseCard(BDF_Card):
    #def __init__(self,card):
    #    pass

    def Is(self,typeCheck):
        if self.type==typeCheck:
            return True
        return False

    def removeTrailingNones(self,fields):
        self.wipeEmptyFields(fields)

    def printCard(self,fields):
        return printCard(fields)

    def setDefaultIfBlank(self,value,default):
        raise Exception('time to upgrade...')
        return setDefaultIfBlank(value,default)

    def setBlankIfDefault(self,value,default):
        return setBlankIfDefault(value,default)

    def crossReference(self,mesh):
        #self.mid = mesh.Material(self.mid)
        raise Exception('%s needs to implement this method' %(self.type))

   # def off_expandThru(self,fields):
   #     """
   #     not used...
   #     expands a list of values of the form [1,5,THRU,10,13]
   #     to be [1,5,6,7,8,9,10,13]
   #     """
   #     nFields = len(fields)
   #     fieldsOut = []
   #     for i in range(nFields):
   #         if fields[i]=='THRU':
   #             for j in range(fields[i-1],fields[i+1]):                    
   #                 fieldsOut.append(fields[j])
   #             ###
   #         else:
   #             fieldsOut.append(fields[i])
   #         ###
   #     ###
   #     return fieldsOut

    def expandThru(self,fields):
        """
        expands a list of values of the form [1,5,THRU,9,13]
        to be [1,5,6,7,8,9,13]
        """
        if len(fields)==1: return fields
        #print "expandThru"
        #print "fields = ",fields
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
                out.append(fields[i])
                i+=1
            ###
        ###
        #print "out = ",out,'\n'
        return list(set(out))
    
    def expandThruBy(self,fields):
        """
        expands a list of values of the form [1,5,THRU,9,13]
        to be [1,5,7,9,13]
        """
        if len(fields)==1: return fields
        #print "expandThruBy"
        #print "fields = ",fields
        out = []
        nFields = len(fields)
        i=0
        by = 1
        while(i<nFields):
            if fields[i]=='THRU':
                by = 1
                if i+2<nFields and fields[i+2]=='BY':
                    by = fields[i+3]
                    print "BY was found...untested..."
                for j in range(fields[i-1],fields[i+1],by):
                    out.append(j)
                ###
                if by>1:
                    i+=3
                else:
                    i+=2
                ###
            else:
                out.append(fields[i])
                i+=1
            ###
        ###
        #print "out = ",out,'\n'
        return list(set(out))

    def expandThruExclude(self,fields):
        """
        EXCLUDE isnt supported
        """
        return self.expandThru(fields)
        
        nFields = len(fields)
        for i in range(nFields):
            if fields[i]=='THRU':
                for j in range(fields[i-1],fields[i+1]):
                    fieldsOut.append(fields[j])
                ###
            else:
                fieldsOut.append(fields[i])
            ###
        ###


    def collapseThru(self,fields):
        return fields

    def collapseThruBy(self,fields):
        return fields

    def _collapseThru(self,fields):
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
    mid = 0 # ???
    def __init__(self,card):
        #self.type = card[0]
        pass

    def Mid(self):
        #print str(self)
        if isinstance(self.mid,int):
            return self.mid
        else:
            return self.mid.mid
        ###

    def crossReference(self,model):
        self.mid = model.Material(self.mid)
        
    def __repr__(self):
        fields = [self.type,self.Pid()]
        return self.printCard(fields)

class Element(BaseCard):
    pid = 0 # CONM2, rigid
    def __init__(self,card):
        #displayCard(card)
        self.nodes = None
        self.eid = int(card.field(1))
        #self.nids = []
        pass

    def Pid(self):
        if isinstance(self.pid,int):
            return self.pid
        else:
            return self.pid.pid
        ###

    def nodePositions(self):
        return [node.Position() for node in self.nodes]

    def nodeIDs(self):
        if isinstance(self.nodes[0],int):
            #print 'if'
            return [node     for node in self.nodes]
        else:
            #print 'else'
            return [node.nid for node in self.nodes]
        ###

    def prepareNodeIDs(self,nids,allowEmptyNodes=False):
        self.nodes = []
        for nid in nids:
            if isinstance(nid,int):
                self.nodes.append(int(nid))
            elif nid==None and allowEmptyNodes:
                self.nodes.append(nid)
            else: # string???
                self.nodes.append(int(nid))
                #raise Exception('this element may not have missing nodes...nids=%s allowEmptyNodes=False' %(nids))
            ###

    def Centroid(self,nodes,debug=False):
        return None

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
        fields = [self.type,self.eid,self.Pid()]+self.nodeIDs()
        return self.printCard(fields)

    def length(self):
        raise Exception('length not implemented in the %s class' %(self.type))
    def area(self):
        raise Exception('area not implemented in the %s class' %(self.type))
    def volume(self):
        raise Exception('volume not implemented in the %s class' %(self.type))
    def mass(self):
        raise Exception('mass not implemented in the %s class' %(self.type))

    def Jacobian(self):
        raise Exception('Jacobian not implemented for %s' %(self.type))
    def stiffnessMatrix(self):
        raise Exception('stiffnessMatrix not implemented in the %s class' %(self.type))
    def massMatrix(self):
        raise Exception('massMatrix not implemented in the %s class' %(self.type))


#dnMax = 2
if __name__=='__main__':
    card = BaseCard()

    """
    1,THRU,10
    1,3,THRU,19,15
    """
    card.collapseThru([1,2,3,4,5,10])
    card.collapseThru([1,3,4,5,6,17])

