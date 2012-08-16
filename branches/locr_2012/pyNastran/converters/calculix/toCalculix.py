
from pyNastran.bdf.bdf import BDF, PBAR, PBARL, PBEAM, PBEAML


class CalculixConverter(BDF):
    """
    Converts a BDF to Calculix (inp/dat/py files).

    @warning Totally inaccurate....

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

    TODO:
      * PCOMP
      * SPC, SPC1, MPC
      * RBE2, RBE3
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
            Type = element.calculixType
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

    def Calculix_Executive(self):
        inp = ''
        if self.sol == 101:
            inp += 'MECA_STATIQUE % SOL 101 - linear statics\n'
            inp += 'stat(MECA_STATIQUE(MODELE=model,CHAM_MATER=material,CARA_ELEM=elemcar,\n'

            inp += 'ECIT=(_F(Charge=AllBoundaryConditions,),\n',
            inp += '      _F(Charge=AllLoads,),\n',
            inp += '      ),\n',

            inp += "TITRE='My Title'\n"
        return inp

    def Calculix_Nodes(self):
        dat = ''
        dat += '** Calculix_Nodes\n'
        dat += '*NODE, NSET=AllNodes\n'

        form = '%-' + str(self.maxNIDlen) + 's %8s,%8s,%8s\n'

        for nid, node in sorted(self.nodes.iteritems()):
            p = node.Position()
            dat += form % (nid, p[0], p[1], p[2])
        ###
        dat += '\n\n'
        dat += self.breaker()
        return dat

    def Calculix_Elements(self):
        """
        @todo sort elements by Type and Material ID
        """
        dat = ''
        dat += '** Calculix_Elements\n'

        elems = self.getElementsByType()
        #print elems
        formE = '%-' + str(self.maxEIDlen) + 's, '
        for Type, eids in sorted(elems.iteritems()):
            dat += '** eid,n1,n2,n3,etc... for a %s\n' % (
                self.elements[eids[0]].type)
            dat += '*ELEMENT, TYPE=%s, ELSET=Elements%i\n' % (Type, 999)
            for eid in eids:
                dat += formE % (eid)
                element = self.elements[eid]
                #print element
                for nid in element.nodeIDs():
                    #print nid
                    dat += '%s,' % (nid)
                dat = dat[:-1] + '\n'
        dat += self.breaker()
        #print dat
        return dat

    def Calculix_Properties(self):
        inp = ''
        inp += '** Calculix_Properties\n'

        #p = []
        #for pid,prop in sorted(self.properties.iteritems()):
        #    p.append('%s_%s' %(prop.type,pid))
        #p = str(p)[1:-1] # chops the [] signs
        #inp += "MODEL=AFFE_MODELE(MAILLAGE=MESH,\n"
        #inp += "          AFFE=_F(GROUP_MA=(%s),\n" %(p)
        #inp += "                  PHENOMENE='MECANIQUE',\n"
        #inp += "                  MODELISATION=('POU_D_T'),),);\n\n"

        inp += "Prop = AFFE_CARA_ELEM(MODELE=FEMODL,\n"
        pyCA = ''
        iCut = 0
        iFace = 0
        iStart = 0
        for pid, prop in sorted(self.properties.iteritems()):
            if isinstance(prop, PBARL) or isinstance(prop, PBEAML):
                (pyCAi, iCut, iFace, iStart) = prop.writeCalculix(
                    iCut, iFace, iStart)
                pyCA += pyCAi
                isSkipped = False
            else:
                try:
                    prop = prop.writeCalculix()
                except:
                    print prop
                    raise
                isSkipped = False
                if 'skipped' in prop:
                    isSkipped = True
            ###
        if not isSkipped:
            inp = inp[:-2]
        inp += ');\n'
        #inp += ');\nFINSF\n\n'
        inp += self.breaker()
        return inp, pyCA

    def Calculix_Materials(self):
        """
        might need to make this by pid instead...
        steel=DEFI_MATERIAU(ELAS=_F(E=210000.,NU=0.3,RHO=8e-9),);
        """
        inp = ''
        inp += '** Calculix_Materials\n'
        mats = self.getElementsByMid()
        for mid, material in sorted(self.materials.iteritems()):
            #inp += 'GROUP_MA name = %s_%s\n' %(material.type,mid)
            inp += material.writeCalculix()

            eids = mats[mid]
            #inp += '    '
            #for eid in eids:
            #    inp += 'elem%s ' %(eid)
            #inp = inp[:-1]
            #inp += '\n'
        #inp = inp[:-2]
        #inp += '\n'
        #inp += ');\n'
        #inp += 'FINSF\n\n'
        inp += self.breaker()
        return inp

    def Calculix_MaterialField(self):
        """
        MtrlFld=AFFE_MATERIAU(MAILLAGE=MESH,
                              AFFE=(_F(GROUP_MA=('P32','P33','P42','P43','P46','P47','P48','P49','P61','P62','P63','P64','P65','P74',
                                                 'P75',),
                                       MATER=M3,),
                                    _F(GROUP_MA=('P11','P13','P14','P15','P55','P56','P59',),
                                       MATER=M6,),
        """
        inp = ''
        inp += '** Calculix_MaterialField\n'
        inp += 'MtrlFld=AFFE_MATERIAU(MAILLAGE=MESH,\n'
        inp += '                      AFFE=(\n'

        mat2Props = self.getPropertiesByMid()
        for mid, material in sorted(self.materials.iteritems()):
            inp += '                      _F(GROUP_MA=('
            pids = mat2Props[mid]
            #inp += "                      "
            for pid in pids:
                inp += "'P%s'," % (pid)
            inp = inp[:-1] + '),\n'
            inp += "                      MATER=M%s),\n" % (mid)
        ###
        inp = inp[:-1] + '));\n'

        inp += self.breaker()
        return inp

    def Calculix_Loads(self):
        """writes the load cards sorted by ID"""
        inp = '** Calculix_Loads\n'
        #if self.language=='english':
            #inp += '** Loads\n'
        #else:
            #inp += ''

        iSubcase = 1
        paramName = 'LOAD'

        #skippedLids = {}
        if self.loads:
            inp += '** LOADS\n'
            loadKeys = self.loads.keys()
            if 1:
                key = self.caseControlDeck.getSubcaseParameter(
                    iSubcase, paramName)[0]
                loadcase = self.loads[key]
                #print loadcase
                for i, load in enumerate(loadcase):
                    inp += '** main LOAD lid=%s type=%s\n' % (loadcase[i].lid,
                                                              loadcase[i].__class__.__name__)

                    #try:
                    if 1:  # LOAD card
                        out = load.writeCalculixLoad(self)
                        if len(out) == 3:  # LOAD card
                            (inpi, loadIDs, loadTypes) = out
                            inp += inpi
                        else:  # FORCEx, MOMENTx, GRAV
                            #skippedLids[(load.lid,load.type)] = out
                            inp += out
                    #except:
                        #print 'failed printing load...type=%s key=%s' %(load.type,key)
                        #raise
                    ###
            #loadcase.
            #for ID,grav in sorted(self.gravs.iteritems()):
            #    inp += grav.writeCalculix(mag)
            ###
        ###

        #for lid_loadType,inpi in sorted(skippedLids.iteritems()):
            #inp += inpi

        inp += self.breaker()
        return inp

    def Calculix_SPCs(self):
        #for spcID,spcs in self.spcObject2.iteritems():
        inp = ''
        inp += '** Calculix_SPCs\n'
        inp += self.breaker()
        return inp

    def breaker(self):
        return '**-------------------------------------------------------------------------\n'

    def buildMaxs(self):
        self.maxNIDlen = len(str(max(self.nodes)))
        self.maxEIDlen = len(str(max(self.elements)))
        self.maxPIDlen = len(str(max(self.properties)))
        self.maxMIDlen = len(str(max(self.materials)))

    def writeAsCalculix(self, fname='fem'):
        inp = ''
        dat = ''
        self.buildMaxs()  # gets number of nodes/elements/properties/materials

        inp += '** BEGIN BULK\n'
        inp += 'DEBUT();\n\n'

        inp += "**'Read the mesh' - we use the 'aster' file format here.\n"
        inp += 'mesh=LIRE_MAILLAGE(UNITE=20,\n'
        inp += "                   FORMAT='ASTER');\n\n"

        #inp += "**'MECA_STATIQUE' % SOL 101 - linear statics\n"
        inp += "** Assigning the model for which CA will calculate the results:\n"
        inp += "** 'Mecanique' - since we are dealing with a linear elastic model and '3D' since it's a 3D model.\n"
        inp += 'Meca=AFFE_MODELE(MAILLAGE=mesh,\n'
        inp += "                 AFFE=_F(TOUT='OUI',\n"
        inp += "                         PHENOMENE='MECANIQUE',\n"
        inp += "                         MODELISATION='3D',),);\n\n"
        inp += self.breaker()

        dat += self.Calculix_Nodes()
        dat += self.Calculix_Elements()
        dat += self.Calculix_Materials()
        inp += self.Calculix_MaterialField()
        (inpi, pyCA) = self.Calculix_Properties()
        inp += inpi
        #inp += self.Calculix_Loads()
        inp += self.Calculix_SPCs()

        # Case Control Deck
        inp += '*STEP\n'
        inp += '*STATIC\n'
        inp += '*CLOAD\n'
        inp += 'LAST,2,1.\n'
        inp += '*NODE PRINT,NSET=NALL\n'
        inp += 'U\n'
        inp += '*EL PRINT,ELSET=EALL\n'
        inp += 'S\n'
        inp += '*END STEP\n'
        inp += '** END OF DATA\n'
        print "writing fname=%s" % (fname + '.inp')
        f = open(fname + '.inp', 'wb')
        f.write(inp)
        f.close()

        if dat:
            f = open(fname + '.dat', 'wb')
            print "writing fname=%s" % (fname + '.dat')
            f.write(dat)
            f.close()

        #if pyCA:
            #f = open(fname+'.py','wb')
            #print "writing fname=%s" %(fname+'.py')
            #f.write(pyCA)
            #f.close()
        ###


def main():
    import sys
    ca = CalculixConverter()
    #model = 'solidBending'
    model = sys.argv[1]
    ca.readBDF(model)
    ca.writeAsCalculix(model)  # inp, py

if __name__ == '__main__':
    main()
