
from pyNastran.bdf.bdf import BDF, PBARL, PBEAML  # PBAR,PBEAM,


class CodeAsterConverter(BDF):
    """
    Converts a BDF to Code Aster (comm/mail/py files).
    How:
    
     * Nodes/Coordinate Systems/Elements/Properties/Materials are
       directly extracted from the BDF.  All objects must reference
       each other properly.
     * Just like Nastran, extra materials/properties are allowed.
       No idea how Code_Aster handles SPOINTs or unassociated GRIDs.
     * Loads must be referenced by a single LOAD card in the Case Control deck.
       This is consistent with standard Nastran.

    Limitations:
    
     * All Case Control inputs must come from SUBCASE 1.
     * LOAD cards must bound FORCEx/MOMENTx/PLOAD4 cards in order for loads to be written
     * Only SOL 101 (Static)

    Supported Cards:
    
     * GRID, COORDx
     * LOAD, FORCEx, MOMENTx, PLOAD4
     * CBAR, CBEAM, CROD, CTUBE, CTETRA, CPENTA, CHEXA,CTRIA3/6, CQUAD4/8
     * PBAR, PBEAM, PROD, PTUBE, PSOLID, PSHELL
     * MAT1
     * GRAV (incorrect writing, but really easy to make it correct given proper format)

    @todo
      PCOMP,
      SPC, SPC1, MPC,
      RBE2, RBE3
    """
    def __init__(self, language='english'):
        self.language = 'english'
        BDF.__init__(self)

    def getElementsByPid(self):
        """builds a dictionary where the key is the property ID and the value is a list of element IDs"""
        props = {}
        for pid in self.properties:
            props[pid] = []
        for eid, element in self.elements.iteritems():
            pid = element.Pid()
            props[pid].append(eid)
        ###
        return mats

    def getElementsByMid(self):
        """builds a dictionary where the key is the material ID and the value is a list of element IDs"""
        mats = {0: []}

        for mid in self.materials:
            mats[mid] = []
        for eid, element in self.elements.iteritems():
            try:
                mid = element.Mid()
                mats[mid].append(eid)
            except:
                mats[0].append(eid)
        ###
        return mats

    def getElementsByType(self):
        """builds a dictionary where the key is the element type and the value is a list of element IDs"""
        elems = {}
        #for eid,elements in self.elements:
            #elems[eid] = []
        for eid, element in self.elements.iteritems():
            Type = element.asterType
            if Type not in elems:
                elems[Type] = []

            elems[Type].append(eid)
            #mid = element.Mid()
            #mats[mid].append(eid)
        ###
        return elems

    def getPropertiesByMid(self):
        """builds a dictionary where the key is the material ID and the value is a list of property IDs"""
        mats = {0: []}

        for mid in self.materials:
            mats[mid] = []
        for pid, property in self.properties.iteritems():
            try:
                mid = property.Mid()
                mats[mid].append(pid)
            except:
                mats[0].append(pid)
        ###
        return mats

    def CA_Executive(self):
        comm = ''
        if self.sol == 101:
            comm += 'MECA_STATIQUE % SOL 101 - linear statics\n'
            comm += 'stat(MECA_STATIQUE(MODELE=model,CHAM_MATER=material,CARA_ELEM=elemcar,\n'

            comm += 'ECIT=(_F(Charge=AllBoundaryConditions,),\n',
            comm += '      _F(Charge=AllLoads,),\n',
            comm += '      ),\n',

            comm += "TITRE='My Title'\n"

        if self.sol == 101:  # [K][U] = [F] #Kx=F
            pass
        elif self.sol == 103:  # [F] = [M][\ddot U] + [K][U] => U(s^2*M+K)=0  phi=det(M-lambda*K)
            pass
        elif self.sol == 129:  # [M][\ddot U] + [C][\dot U] + [K] [U] = [F]
            pass

        k = "#Calculate data for the stiffness Matrix\n"
        k += "StiffMtx=CALC_MATR_ELEM( OPTION='RIGI_MECA',MODELE=ModelDef,CHAM_MATER=MtrlFld,);\n\n"
        m = "#Calculate data for the Mass Matrix\n"
        m += "MassMtx=CALC_MATR_ELEM( OPTION='MASS_MECA',MODELE=ModelDef,CHAM_MATER=MtrlFld,);\n\n"

        K = "#Assign the Stiffness Matrix to the DOFs to be solved\n"
        K += "K=ASSE_MATRICE(MATR_ELEM=StiffMtx,NUME_DDL=NDOFs,);\n\n"
        M = "#Assign the Mass Matrix to the DOFs to be solved\n"
        M += "M=ASSE_MATRICE(MATR_ELEM=MassMtx,NUME_DDL=NDOFs,);\n"
        return comm

    def CA_Nodes(self, gridWord='grid'):
        mail = ''
        mail += '# CA_Nodes\n'
        if self.language == 'english':
            mail += '# Grid Points\n'
        else:
            mail += ''

        mail += 'COOR_3D\n'
        form = '    %s%-' + str(self.maxNIDlen) + 's %8s %8s %8s\n'

        for nid, node in sorted(self.nodes.iteritems()):
            p = node.Position()
            mail += form % (gridWord, nid, p[0], p[1], p[2])
        ###
        mail += 'FINSF\n\n'
        mail += self.breaker()
        return mail

    def CA_Elements(self, elemWord='elem', gridWord='grid'):
        mail = ''
        mail += '# CA_Elements\n'
        if self.language == 'english':
            mail += '# Elements\n'
        else:
            mail += ''

        elems = self.getElementsByType()

        formE = '    %s%-' + str(self.maxEIDlen) + 's '
        formG = '%s%-' + str(self.maxNIDlen) + 's '
        for Type, eids in sorted(elems.iteritems()):
            mail += '%s\n' % (Type)
            for eid in eids:
                mail += formE % (elemWord, eid)
                element = self.elements[eid]
                for nid in element.nodeIDs():
                    mail += formG % (gridWord, nid)
                mail += '\n'
            mail += 'FINSF\n\n'
        mail += self.breaker()
        return mail

    def CA_Properties(self):
        comm = ''
        comm += '# CA_Properties\n'
        if self.language == 'english':
            comm += '# Properties\n'
        else:
            comm += ''

        #p = []
        #for pid,prop in sorted(self.properties.iteritems()):
        #    p.append('%s_%s' %(prop.type,pid))
        #p = str(p)[1:-1] # chops the [] signs
        #comm += "MODEL=AFFE_MODELE(MAILLAGE=MESH,\n"
        #comm += "          AFFE=_F(GROUP_MA=(%s),\n" %(p)
        #comm += "                  PHENOMENE='MECANIQUE',\n"
        #comm += "                  MODELISATION=('POU_D_T'),),);\n\n"

        comm += "Prop = AFFE_CARA_ELEM(MODELE=FEMODL,\n"
        pyCA = ''
        iCut = 0
        iFace = 0
        iStart = 0
        for pid, prop in sorted(self.properties.iteritems()):
            if isinstance(prop, PBARL) or isinstance(prop, PBEAML):
                (pyCAi, iCut, iFace, iStart) = prop.writeCodeAster(
                    iCut, iFace, iStart)
                pyCA += pyCAi
                isSkipped = False
            else:
                prop = prop.writeCodeAster()
                isSkipped = False
                if 'skipped' in prop:
                    isSkipped = True
            ###
        if not isSkipped:
            comm = comm[:-2]
        comm += ');\n'
        #comm += ');\nFINSF\n\n'
        comm += self.breaker()
        return comm, pyCA

    def CA_Materials(self):
        """
        might need to make this by pid instead...
        steel=DEFI_MATERIAU(ELAS=_F(E=210000.,NU=0.3,RHO=8e-9),);
        """
        comm = ''
        comm += '# CA_Materials\n'
        if self.language == 'english':
            comm += '# Materials\n'
        else:
            comm += ''
        mats = self.getElementsByMid()
        for mid, material in sorted(self.materials.iteritems()):
            #comm += 'GROUP_MA name = %s_%s\n' %(material.type,mid)
            comm += material.writeCodeAster()

            eids = mats[mid]
            #comm += '    '
            #for eid in eids:
            #    comm += 'elem%s ' %(eid)
            #comm = comm[:-1]
            #comm += '\n'
        #comm = comm[:-2]
        #comm += '\n'
        #comm += ');\n'
        #comm += 'FINSF\n\n'
        comm += self.breaker()
        return comm

    def CA_MaterialField(self):
        """
        @code
        MtrlFld=AFFE_MATERIAU(MAILLAGE=MESH,
                              AFFE=(_F(GROUP_MA=('P32','P33','P42','P43','P46','P47','P48','P49','P61','P62','P63','P64','P65','P74',
                                                 'P75',),
                                       MATER=M3,),
                                    _F(GROUP_MA=('P11','P13','P14','P15','P55','P56','P59',),
                                       MATER=M6,),
        @endcode
        """
        comm = ''
        comm += '# CA_MaterialField\n'
        comm += 'MtrlFld=AFFE_MATERIAU(MAILLAGE=MESH,\n'
        comm += '                      AFFE=(\n'

        mat2Props = self.getPropertiesByMid()
        for mid, material in sorted(self.materials.iteritems()):
            comm += '                      _F(GROUP_MA=('
            pids = mat2Props[mid]
            #comm += "                      "
            for pid in pids:
                comm += "'P%s'," % (pid)
            comm = comm[:-1] + '),\n'
            comm += "                      MATER=M%s),\n" % (mid)
        ###
        comm = comm[:-1] + '));\n'

        comm += self.breaker()
        return comm

    def CA_Loads(self):
        """writes the load cards sorted by ID"""
        comm = '# CA_Loads\n'
        #if self.language=='english':
            #comm += '# Loads\n'
        #else:
            #comm += ''

        iSubcase = 1
        paramName = 'LOAD'

        #skippedLids = {}
        if self.loads:
            comm += '# LOADS\n'
            #loadKeys = self.loads.keys()

            key = self.caseControlDeck.getSubcaseParameter(
                iSubcase, paramName)[0]
            loadcase = self.loads[key]
            #print loadcase
            for i, load in enumerate(loadcase):
                comm += '# main LOAD lid=%s type=%s\n' % (loadcase[
                    i].lid, loadcase[i].__class__.__name__)

                #try:
                if 1:  # LOAD card
                    out = load.writeCodeAsterLoad(self, gridWord='N')
                    if len(out) == 3:  # LOAD card
                        (commi, loadIDs, loadTypes) = out
                        comm += commi
                    else:  # FORCEx, MOMENTx, GRAV
                        #skippedLids[(load.lid,load.type)] = out
                        comm += out
                #except:
                    #print 'failed printing load...type=%s key=%s' %(load.type,key)
                    #raise
                ###
            #loadcase.
            #for ID,grav in sorted(self.gravs.iteritems()):
            #    comm += grav.writeCodeAster(mag)
            ###
        ###

        #for lid_loadType,commi in sorted(skippedLids.iteritems()):
            #comm += commi

        comm += self.breaker()
        return comm

    def CA_SPCs(self):
        #for spcID,spcs in self.spcObject2.iteritems():
        comm = ''
        comm += '# CA_SPCs\n'
        comm += self.breaker()
        return comm

    def breaker(self):
        return '#-------------------------------------------------------------------------\n'

    def buildMaxs(self):
        self.maxNIDlen = len(str(max(self.nodes)))
        self.maxEIDlen = len(str(max(self.elements)))
        self.maxPIDlen = len(str(max(self.properties)))
        self.maxMIDlen = len(str(max(self.materials)))

    def writeAsCodeAster(self, fname='fem'):
        comm = ''
        mail = ''
        self.buildMaxs()  # gets number of nodes/elements/properties/materials

        comm += '# BEGIN BULK\n'
        comm += 'DEBUT();\n\n'

        comm += "#'Read the mesh' - we use the 'aster' file format here.\n"
        comm += 'mesh=LIRE_MAILLAGE(UNITE=20,\n'
        comm += "                   FORMAT='ASTER');\n\n"

        #comm += "#'MECA_STATIQUE' % SOL 101 - linear statics\n"
        comm += "# Assigning the model for which CA will calculate the results:\n"
        comm += "# 'Mecanique' - since we are dealing with a linear elastic model and '3D' since it's a 3D model.\n"
        comm += 'Meca=AFFE_MODELE(MAILLAGE=mesh,\n'
        comm += "                 AFFE=_F(TOUT='OUI',\n"
        comm += "                         PHENOMENE='MECANIQUE',\n"
        comm += "                         MODELISATION='3D',),);\n\n"
        comm += self.breaker()

        mail += self.CA_Nodes(gridWord='N')
        mail += self.CA_Elements(elemWord='E', gridWord='N')
        comm += self.CA_Materials()
        comm += self.CA_MaterialField()
        (commi, pyCA) = self.CA_Properties()
        comm += commi
        comm += self.CA_Loads()
        comm += self.CA_SPCs()

        comm += 'FIN();\n'
        comm += '# ENDDATA\n'
        print "writing fname=%s" % (fname + '.comm')
        f = open(fname + '.comm', 'wb')
        f.write(comm)
        f.close()

        if mail:
            f = open(fname + '.mail', 'wb')
            print "writing fname=%s" % (fname + '.mail')
            f.write(mail)
            f.close()

        if pyCA:
            f = open(fname + '.py', 'wb')
            print "writing fname=%s" % (fname + '.py')
            f.write(pyCA)
            f.close()
        ###


def main():
    import sys
    ca = CodeAsterConverter()
    #model = 'solidBending'
    model = sys.argv[1]
    ca.readBDF(model)
    ca.writeAsCodeAster(model)  # comm, py

if __name__ == '__main__':
    main()
