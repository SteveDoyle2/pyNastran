from __future__ import print_function
#from codec import open as codec_open
import os
from six import iteritems
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

    def get_elements_by_pid(self):
        """
        builds a dictionary where the key is the property ID and the value is
        a list of element IDs
        """
        props = {}
        for pid in self.properties:
            props[pid] = []
        for eid, element in iteritems(self.elements):
            pid = element.Pid()
            props[pid].append(eid)
        mats = []
        return mats

    def get_elements_by_mid(self):
        """
        builds a dictionary where the key is the material ID and the value is
        a list of element IDs
        """
        mats = {0: []}

        for mid in self.materials:
            mats[mid] = []
        for eid, element in iteritems(self.elements):
            try:
                mid = element.Mid()
                mats[mid].append(eid)
            except:
                mats[0].append(eid)
        return mats

    def get_elements_by_type(self):
        """
        builds a dictionary where the key is the element type and the
        value is a list of element IDs
        """
        elems = {}
        #for eid,elements in self.elements:
            #elems[eid] = []
        for eid, element in iteritems(self.elements):
            if not hasattr(element, 'aster_type'):
                print('rejecting: %s' % element.type)
                continue
            Type = element.aster_type
            if Type not in elems:
                elems[Type] = []

            elems[Type].append(eid)
            #mid = element.Mid()
            #mats[mid].append(eid)
        return elems

    def get_properties_by_mid(self):
        """
        builds a dictionary where the key is the material ID and the
        value is a list of property IDs
        """
        mats = {0: []}

        for mid in self.materials:
            mats[mid] = []
        for pid, property in iteritems(self.properties):
            try:
                mid = property.Mid()
                mats[mid].append(pid)
            except:
                mats[0].append(pid)
        return mats

    def ca_executive(self):
        comm = ''
        if self.sol == 101:
            comm += 'MECA_STATIQUE % SOL 101 - linear statics\n'
            comm += 'stat(MECA_STATIQUE(MODELE=model, CHAM_MATER=material, CARA_ELEM=elemcar,\n'

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
        k += "StiffMtx = CALC_MATR_ELEM(OPTION='RIGI_MECA', MODELE=ModelDef, CHAM_MATER=MtrlFld);\n\n"
        m = "#Calculate data for the Mass Matrix\n"
        m += "MassMtx = CALC_MATR_ELEM(OPTION='MASS_MECA', MODELE=ModelDef, CHAM_MATER=MtrlFld);\n\n"

        K = "#Assign the Stiffness Matrix to the DOFs to be solved\n"
        K += "K = ASSE_MATRICE(MATR_ELEM=StiffMtx, NUME_DDL=NDOFs);\n\n"
        M = "#Assign the Mass Matrix to the DOFs to be solved\n"
        M += "M = ASSE_MATRICE(MATR_ELEM=MassMtx, NUME_DDL=NDOFs);\n"
        return comm

    def ca_nodes(self, grid_word='grid'):
        mail = ''
        mail += '# ca_nodes\n'
        if self.language == 'english':
            mail += '# Grid Points\n'
        else:
            mail += ''

        mail += 'COOR_3D\n'
        form = '    %s%-' + str(self.max_nid_len) + 's %8s %8s %8s\n'

        for nid, node in sorted(iteritems(self.nodes)):
            p = node.get_position()
            mail += form % (grid_word, nid, p[0], p[1], p[2])
        mail += 'FINSF\n\n'
        mail += self.breaker()
        return mail

    def ca_elements(self, elem_word='elem', grid_word='grid'):
        mail = ''
        mail += '# ca_elements\n'
        if self.language == 'english':
            mail += '# elements\n'
        else:
            mail += ''

        elems = self.get_elements_by_type()

        form_e = '    %s%-' + str(self.max_eid_len) + 's '
        form_g = '%s%-' + str(self.max_nid_len) + 's '
        for Type, eids in sorted(iteritems(elems)):
            mail += '%s\n' % (Type)
            for eid in eids:
                mail += form_e % (elem_word, eid)
                element = self.elements[eid]
                for nid in element.node_ids:
                    mail += form_g % (grid_word, nid)
                mail += '\n'
            mail += 'FINSF\n\n'
        mail += self.breaker()
        return mail

    def ca_properties(self):
        comm = ''
        comm += '# ca_properties\n'
        if self.language == 'english':
            comm += '# Properties\n'
        else:
            comm += ''

        p = []
        for pid, prop in sorted(iteritems(self.properties)):
            p.append('%s_%s' % (prop.type, pid))
        p = str(p)[1:-1] # chops the [] signs
        comm += "MODEL=AFFE_MODELE(MAILLAGE=MESH,\n"
        comm += "          AFFE=_F(GROUP_MA=(%s),\n" %(p)
        comm += "                  PHENOMENE='MECANIQUE',\n"
        comm += "                  MODELISATION=('POU_D_T'),),);\n\n"

        comm += "Prop = AFFE_CARA_ELEM(MODELE=FEMODL,\n"
        py_ca = ''
        icut = 0
        iface = 0
        istart = 0
        for pid, prop in sorted(iteritems(self.properties)):
            if isinstance(prop, (PBARL, PBEAML)):
                (py_cai, icut, iface, istart) = prop.write_code_aster(
                    iCut, iFace, istart)
                py_ca += py_cai
                isSkipped = False
            else:
                prop = prop.write_code_aster()
                is_skipped = False
                if 'skipped' in prop:
                    is_skipped = True

        if not is_skipped:
            comm = comm[:-2]
        comm += ');\n'
        #comm += ');\nFINSF\n\n'
        comm += self.breaker()
        return comm, py_ca

    def ca_materials(self):
        """
        might need to make this by pid instead...
        steel=DEFI_MATERIAU(ELAS=_F(E=210000.,NU=0.3,RHO=8e-9),);
        """
        comm = ''
        comm += '# ca_materials\n'
        if self.language == 'english':
            comm += '# materials\n'
        else:
            comm += ''
        mats = self.get_elements_by_mid()
        for mid, material in sorted(iteritems(self.materials)):
            #comm += 'GROUP_MA name = %s_%s\n' % (material.type, mid)
            comm += material.write_code_aster()

            eids = mats[mid]
            #comm += '    '
            #for eid in eids:
            #    comm += 'elem%s ' % (eid)
            #comm = comm[:-1]
            #comm += '\n'
        #comm = comm[:-2]
        #comm += '\n'
        #comm += ');\n'
        #comm += 'FINSF\n\n'
        comm += self.breaker()
        return comm

    def ca_material_field(self):
        """
        .. code-block:: ptyhon

          MtrlFld=AFFE_MATERIAU(MAILLAGE=MESH,
                              AFFE=(_F(GROUP_MA=('P32','P33','P42','P43','P46','P47','P48','P49','P61','P62','P63','P64','P65','P74',
                                                 'P75',),
                                       MATER=M3,),
                                    _F(GROUP_MA=('P11','P13','P14','P15','P55','P56','P59',),
                                       MATER=M6,),
        @endcode
        """
        comm = ''
        comm += '# ca_material_field\n'
        comm += 'MtrlFld=AFFE_MATERIAU(MAILLAGE=MESH,\n'
        comm += '                      AFFE=(\n'

        mat_to_props = self.get_properties_by_mid()
        for mid, material in sorted(iteritems(self.materials)):
            comm += '                      _F(GROUP_MA=('
            pids = mat_to_props[mid]
            #comm += "                      "
            for pid in pids:
                comm += "'P%s'," % (pid)
            comm = comm[:-1] + '),\n'
            comm += "                      MATER=M%s),\n" % (mid)

        comm = comm[:-1] + '));\n'
        comm += self.breaker()
        return comm

    def ca_loads(self):
        """writes the load cards sorted by ID"""
        comm = '# ca_loads\n'
        #if self.language=='english':
            #comm += '# loads\n'
        #else:
            #comm += ''

        isubcase = 1
        param_name = 'LOAD'

        #skippedLids = {}
        if self.loads:
            comm += '# LOADS\n'
            #loadKeys = self.loads.keys()

            key = self.case_control_deck.get_subcase_parameter(isubcase, param_name)[0]
            loadcase = self.loads[key]
            #print(loadcase)
            for i, load in enumerate(loadcase):
                loadi = loadcase[i]
                comm += '# main LOAD sid=%s type=%s\n' % (loadi.sid, loadi.__class__.__name__)
                if load.type != 'LOAD':
                    msg = 'LOAD card must be referenced in case control deck, not %s' % load.type
                    raise RuntimeError(msg)
                #try:
                if 1:  # LOAD card
                    out = load.write_code_aster_load(self, grid_word='N')
                    if len(out) == 3:  # LOAD card
                        (commi, load_ids, load_types) = out
                        comm += commi
                    else:  # FORCEx, MOMENTx, GRAV
                        #skipped_lids[(load.lid, load.type)] = out
                        comm += out
                #except:
                    #print('failed printing load...type=%s key=%s' % (load.type, key))
                    #raise
            #loadcase.
            #for ID, grav in sorted(iteritems(self.gravs)):
            #    comm += grav.write_code_aster(mag)

        #for lid_load_type, commi in sorted(iteritems(skipped_lids)):
            #comm += commi

        comm += self.breaker()
        return comm

    def ca_spcs(self):
        comm = ''
        comm += '# ca_spcs\n'
        for subcase_id, subcase in iteritems(self.subcases):
            if subcase_id == 0:
                continue
            if not subcase.has_parameter('SPC'):
                continue
            spc_id = subcase.get_parameter('SPC')[0]

            comm += '# SUBCASE %i; SPC=%s\n' % (subcase_id, spc_id)
            spcs = self.get_reduced_spcs(spc_id)
            for spc in spcs:
                for line in str(spc).split('\n'):
                    comm += '#%s\n' % line
        comm += self.breaker()
        return comm

    def ca_mpcs(self):
        comm = ''
        comm += '# ca_mpcs\n'
        for subcase_id, subcase in iteritems(self.subcases):
            if subcase_id == 0:
                continue
            if not 'MPC' in subcase:
                continue
            mpc_id = subcase.get_parameter('MPC')[0]

            comm += '# SUBCASE %i; MPC=%s\n' % (subcase_id, mpc_id)
            mpcs = self.get_reduced_spcs(mpc_id)
            for mpc in mpcs:
                for line in str(mpc).split('\n'):
                    comm += '#%s\n' % line
        comm += self.breaker()
        return comm

    def breaker(self):
        return '#-------------------------------------------------------------------------\n'

    def build_maxs(self):
        self.max_nid_len = len(str(max(self.nodes)))
        self.max_eid_len = len(str(max(self.elements)))
        self.max_pid_len = len(str(max(self.properties)))
        self.max_mid_len = len(str(max(self.materials)))

    def write_as_code_aster(self, model):
        comm = ''
        mail = ''
        self.build_maxs()  # gets number of nodes/elements/properties/materials

        comm += '# BEGIN BULK\n'
        comm += 'DEBUT();\n\n'

        comm += "#'Read the mesh' - we use the 'aster' file format here.\n"
        comm += 'mesh=LIRE_MAILLAGE(UNITE=20,\n'
        comm += "                   FORMAT='ASTER');\n\n"

        #comm += "#'MECA_STATIQUE' % SOL 101 - linear statics\n"
        comm += "# Assigning the model for which CA will calculate the results:\n"
        comm += "#   'Mecanique' - since we are dealing with a linear elastic model\n"
        comm += "#   '3D' since it's a 3D model.\n"
        comm += 'Meca = AFFE_MODELE(MAILLAGE=mesh,\n'
        comm += "                   AFFE=_F(TOUT='OUI',\n"
        comm += "                           PHENOMENE='MECANIQUE',\n"
        comm += "                           MODELISATION='3D'));\n\n"
        comm += self.breaker()

        mail += self.ca_nodes(grid_word='N')
        mail += self.ca_elements(elem_word='E', grid_word='N')
        comm += self.ca_materials()
        comm += self.ca_material_field()
        commi, py_ca = self.ca_properties()
        comm += commi
        comm += self.ca_loads()
        comm += self.ca_spcs()
        comm += self.ca_mpcs()

        comm += 'FIN();\n'
        comm += '# ENDDATA\n'
        assert comm != ''

        print('pwd=', os.getcwd())
        if comm:
            with open(model + '.comm', 'wb') as f:
                print("writing fname=%s" % (model + '.comm'))
                f.write(comm)

        #print(comm)
        #print(mail)
        #print(py_ca)
        if mail:
            print('mail')
            with open(model + '.mail', 'wb') as f:
                print("writing fname=%s" % (model + '.mail'))
                f.write(mail)

        if py_ca:
            print('py_ca')
            assert py_ca != ''
            with open(model + '.py', 'wb') as f:
                print("writing fname=%s" % (model + '.py'))
                f.write(py_ca)


def main():
    import sys
    import pyNastran
    from docopt import docopt
    msg = "Usage:\n"
    msg += "  nastranToCodeAster [-o] BDF_FILENAME\n" #
    msg += '  nastranToCodeAster -h | --help\n'
    msg += '  nastranToCodeAster -v | --version\n'
    msg += '\n'

    msg += "Positional Arguments:\n"
    msg += "  BDF_FILENAME   path to BDF/DAT/NAS file\n"
    msg += '\n'

    msg += 'Options:\n'
    msg += '  -o, --output    prints debug messages (default=False)\n'
    msg += '  -h, --help     show this help message and exit\n'
    msg += "  -v, --version  show program's version number and exit\n"

    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)

    for key, value in sorted(iteritems(data)):
        print("%-12s = %r" % (key.strip('--'), value))

    bdf_filename = data['BDF_FILENAME']
    fname_base = os.path.splitext(bdf_filename)[0]

    ca = CodeAsterConverter()
    ca.read_bdf(bdf_filename, encoding='ascii')
    ca.write_as_code_aster(fname_base)  # comm, py

if __name__ == '__main__':  # pragma: no cover
    main()
