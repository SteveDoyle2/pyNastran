language = 'english'
from pyNastran.bdf.bdf import BDF

class CodeAsterConverter(BDF):
    def __init__(self):
        BDF.__init__(self)
    
    def getElementsByPid(self):
        """builds a dictionary where the key is the property ID and the value is a list of element IDs"""
        props = {}
        for pid in self.properties:
            props[pid] = []
        for eid,element in self.elements.items():
            pid = element.Pid()
            props[pid].append(eid)
        ###
        return mats

    def getElementsByMid(self):
        """builds a dictionary where the key is the material ID and the value is a list of element IDs"""
        mats = {}
        for mid in self.materials:
            mats[mid] = []
        for eid,element in self.elements.items():
            mid = element.Mid()
            mats[mid].append(eid)
        ###
        return mats

    def getElementsByType(self):
        """builds a dictionary where the key is the element type and the value is a list of element IDs"""
        elems = {}
        #for eid,elements in self.elements:
            #elems[eid] = []
        for eid,element in self.elements.items():
            Type = element.type
            if Type not in elems:
                elems[Type] = []

            elems[Type].append(eid)
            #mid = element.Mid()
            #mats[mid].append(eid)
        ###
        return elems

    def CA_Executive(self):
        msg = ''
        if self.sol==101:
            #msg += 'MECA_STATIQUE % SOL 101 - linear statics\n'
            msg += 'stat(MECA_STATIQUE(MODELE=model,CHAM_MATER=material,CARA_ELEM=elemcar,\n'

            msg += 'ECIT=(_F(Charge=AllBoundaryConditions,),\n',
            msg += '      _F(Charge=AllLoads,),\n',
            msg += '      ),\n',

            msg += "TITRE='My Title'\n"
        return msg

    def CA_Nodes(self):
        if language=='english':
            msg = '% Grid Points\n'
        else:
            msg = ''

        msg += 'COORD_3D\n'
        form = '    grid%-'+str(self.maxNIDlen)+'s %8g %8g %8g\n'

        for nid,node in sorted(self.nodes.items()):
            p = node.Position()
            msg += form %(nid,p[0],p[1],p[2])
        ###
        msg += 'FINSF\n\n'
        return msg

    def CA_Elements(self):
        if language=='english':
            msg = '%Elements\n'
        else:
            msg = ''

        elems = self.getElementsByType()

        formE = '    elem%-'+str(self.maxEIDlen)+'s '
        formG =     'grid%-'+str(self.maxNIDlen)+'s '
        for Type,eids in sorted(elems.items()):
            msg += '%s\n' %(Type)
            for eid in eids:
                msg += formE %(eid)
                element = self.elements[eid]
                for nid in element.nodeIDs():
                    msg += formG %(nid)
                msg += '\n'
            msg += 'FINSF\n\n'
        return msg

    def CA_Properties(self):
        if language=='english':
            msg = '%Elements\n'
        else:
            msg = ''
        
        p = []
        for pid,prop in sorted(self.properties.items()):
            p.append('%s_%s' %(prop.type,pid))
        p = str(p)[1:-1] # chops the [] signs
        msg += "MODEL=AFFE_MODELE(MAILLAGE=MESH,\n"
        msg += "          AFFE=_F(GROUP_MA=(%s),\n" %(p)
        msg += "                  PHENOMENE='MECANIQUE',\n"
        msg += "                  MODELISATION=('POU_D_T'),),);\n\n"

        msg += "Prop = AFFE_CARA_ELEM(MODELE=FEMODL,),);\n"
        for pid,prop in sorted(self.properties.items()):
            msg += prop.writeCodeAster()
        msg = msg[:-1]
        msg += ');\nFINSF\n\n'
        return msg

    def CA_Materials(self):
        """
        might need to make this byPid instead...
        steel=DEFI_MATERIAU(ELAS=_F(E=210000.,NU=0.3,RHO=8e-9),);
        """
        if language=='english':
            msg = '%Materials\n'
        else:
            msg = ''
        mats = self.getElementsByMid()
        for mid,material in sorted(self.materials.items()):
            #msg += 'GROUP_MA name = %s_%s\n' %(material.type,mid)
            msg += material.writeCodeAster()

            eids = mats[mid]
            msg += '    '
            for eid in eids:
                msg += 'elem%s ' %(eid)
            msg = msg[:-1]
            msg += '\n'
        msg += 'FINSF\n\n'
        return msg

    def CA_SPCs(self):
        #for spcID,spcs in self.spcObject2.items():
        pass

    def buildMaxs(self):
        self.maxNIDlen = len(str(max(self.nodes)))
        self.maxEIDlen = len(str(max(self.elements)))
        self.maxPIDlen = len(str(max(self.properties)))
        self.maxMIDlen = len(str(max(self.materials)))

    def writeAsCodeAster(self,fname='fem.ca'):
        self.buildMaxs()

        msg = ''
        msg += '% BEGIN BULK\n'
        msg += 'DEBUT();\n\n'

        msg += self.CA_Nodes()
        msg += self.CA_Elements()
        msg += self.CA_Properties()
        msg += self.CA_Materials()

        msg += 'FIN()\n'
        msg += '%ENDDATA\n'
        f = open(fname,'wb')
        f.write(msg)
        f.close()
        

if __name__=='__main__':
    import sys
    ca = CodeAsterConverter()
    #model = 'solidBending'
    model = sys.argv[1]
    ca.readBDF(model+'.bdf')
    ca.writeAsCodeAster(model+'.comm')
    