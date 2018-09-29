"""
very old code aster converter that was never quite working
"""
from __future__ import print_function
import os
from codecs import open
from pyNastran.bdf.bdf import BDF


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
     * LOAD cards must bound FORCEx/MOMENTx/PLOAD4 cards in order for
       loads to be written
     * Only SOL 101 (Static)

    Supported Cards:

     * GRID, COORDx
     * LOAD, FORCEx, MOMENTx, PLOAD4
     * CBAR, CBEAM, CROD, CTUBE, CTETRA, CPENTA, CHEXA, CTRIA3/6, CQUAD4/8
     * PBAR, PBEAM, PROD, PTUBE, PSOLID, PSHELL
     * MAT1
     * GRAV (incorrect writing, but really easy to make it
             correct given proper format)

    @todo
      PCOMP,
      SPC, SPC1, MPC,
      RBE2, RBE3
    """
    def __init__(self, language='english'):
        self.language = language
        assert self.language in ['english']
        BDF.__init__(self)
        self.max_nid_len = 0
        self.max_eid_len = 0
        self.max_pid_len = 0
        self.max_mid_len = 0

    def get_elements_by_pid(self):
        """
        builds a dictionary where the key is the property ID and the value is
        a list of element IDs
        """
        props = {}
        for pid in self.properties:
            props[pid] = []
        for eid, element in self.elements.items():
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
        for eid, element in self.elements.items():
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
        for eid, element in self.elements.items():
            if element.type == 'CTRIA3':
                element_type = 'TRIA3'
            elif element.type == 'CTRIA6':
                element_type = 'TRIA6'
            elif element.type == 'CQUAD4':
                element_type = 'QUAD4'
            elif element.type == 'CQUAD8':
                element_type = 'QUAD8'

            elif element.type == 'CTETRA':
                element_type = 'TETRA4'
            elif element.type == 'CPENTA':
                element_type = 'PENTA6'
            elif element.type == 'CHEXA':
                element_type = 'HEXA8'
            elif element.type == 'CPYRAM':
                element_type = 'CPYRAM5'
            else:
                print('rejecting: %s' % element.type)
                continue
            if element_type not in elems:
                elems[element_type] = []

            elems[element_type].append(eid)
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
        for pid, prop in self.properties.items():
            try:
                mid = prop.Mid()
                mats[mid].append(pid)
            except:
                mats[0].append(pid)
        return mats

    def ca_executive(self):
        """
        Writes the executive control deck

        Only supports SOL 101
        """
        comm = ''
        if self.sol == 101:
            comm += ('ecit = (\n'
                     '    _F(Charge=AllBoundaryConditions,),\n',
                     '    _F(Charge=AllLoads,),)\n'
                     'MECA_STATIQUE % SOL 101 - linear statics\n'
                     'stat(MECA_STATIQUE(MODELE=model, CHAM_MATER=material, '
                     '                   CARA_ELEM=elemcar,\n'
                     '                   ECIT=ecit,\n',
                     "                   TITRE='My Title'\n")

        if self.sol == 101:  # [K][U] = [F] #Kx=F
            pass
        elif self.sol == 103:  # [F] = [M][\ddot U] + [K][U] => U(s^2*M+K)=0  phi=det(M-lambda*K)
            pass
        elif self.sol == 129:  # [M][\ddot U] + [C][\dot U] + [K] [U] = [F]
            pass

        k_msg = ("#Calculate data for the stiffness Matrix\n"
                 "StiffMtx = CALC_MATR_ELEM(OPTION='RIGI_MECA',\n"
                 "                          MODELE=ModelDef, \n"
                 "                          CHAM_MATER=MtrlFld);\n\n")
        mat_msg = ("#Calculate data for the Mass Matrix\n"
                   "MassMtx = CALC_MATR_ELEM(OPTION='MASS_MECA',\n"
                   "                         MODELE=ModelDef,"
                   "                         CHAM_MATER=MtrlFld);\n\n")

        k_msg = "#Assign the Stiffness Matrix to the DOFs to be solved\n"
        k_msg += "K = ASSE_MATRICE(MATR_ELEM=StiffMtx, NUME_DDL=NDOFs);\n\n"
        mat_msg = "#Assign the Mass Matrix to the DOFs to be solved\n"
        mat_msg += "M = ASSE_MATRICE(MATR_ELEM=MassMtx, NUME_DDL=NDOFs);\n"
        return comm

    def ca_nodes(self, grid_word='grid'):
        """writes the GRID cards"""
        mail = ''
        mail += '# ca_nodes\n'
        if self.language == 'english':
            mail += '# Grid Points\n'
        else:
            mail += ''

        mail += 'COOR_3D\n'
        form = '    %s%-' + str(self.max_nid_len) + 's %8s %8s %8s\n'

        for nid, node in sorted(self.nodes.items()):
            xyz = node.get_position()
            mail += form % (grid_word, nid, xyz[0], xyz[1], xyz[2])
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
        for etype, eids in sorted(elems.items()):
            mail += '%s\n' % etype
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
        """writes the properties"""
        comm = ''
        comm += '# ca_properties\n'
        if self.language == 'english':
            comm += '# Properties\n'
        else:
            comm += ''

        p = []
        for pid, prop in sorted(self.properties.items()):
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
        for pid, prop in sorted(self.properties.items()):
            ptype = prop.type
            if ptype == 'PBARL':
                msg = '# BAR Type=%s pid=%s\n' % (ptype, pid)
                msg2 = ''
                msg += code_aster_beam_section(prop, iface, istart, prop.dim)
                iface += 1
                istart += len(prop.dim)

                msg += code_aster_beam_section(prop, iface, istart, prop.dim)
                iface += 1
                msg2 += 'Cut_%s = geompy.MakeCut(Face_%i, Face_%i)\n' % (
                    icut + 1, iface + 1, iface + 2)
                msg2 += "geompy.addToStudy(Cut_%i,  'Cut_%i')\n" % (
                    icut + 1, icut + 1)
                istart += len(prop.dim)
                py_ca += (msg + msg2)
                continue

            elif ptype == 'PBEAML':
                msg = ''
                msg2 = 'Cut_%s = geompy.MakeCut(' % (icut + 1)
                for xxb, dim, nsm in zip(prop.xxb, prop.dim, prop.nsm):
                    msg += code_aster_beam_section(prop, iface, istart, prop.dim)
                    msg2 += 'Face_%i, ' % (iface + 1)
                    iface += 1
                    istart += len(prop.dim)
                msg2 = msg2[-2:]
                msg2 += ')\n'

                msg2 += "geompy.addToStudy(Cut_%i,  'Cut_%i')\n" % (icut + 1, icut + 1)
                icut += 1
                py_ca += (msg + msg2)
                continue


            elif ptype == 'PBAR':
                a = prop.Area()
                iy = prop.I11()
                iz = prop.I22()
                j = prop.J()
                msg = ''
                msg += "    POUTRE=_F(GROUP_MA='P%s', # PBAR\n" % (pid)
                msg += "              SECTION='GENERALE',\n"
                msg += "              CARA=('A','IY','IZ','JX')\n"
                msg += "              VALE=(%g,  %g,  %g,  %g,)\n" % (a, iy, iz, j)
                msg += "              ORIENTATION=(\n"
                msg += "                    CARA=('VECT_Y'),\n"
                msg += "                    VALE=(1.0,0.0,0.0,),),\n"

            elif ptype == 'PBEAM':
                a = prop.Area()
                iy = prop.I11()
                iz = prop.I22()
                j = prop.J()
                msg = "    POUTRE=_F(GROUP_MA='P%s', # PBEAM\n" % pid
                msg += "              SECTION='GENERALE',\n"
                msg += "              CARA=('A','IY','IZ','JX'), # area, moments of inertia\n"
                msg += "              VALE=(%g,  %g,  %g,  %g),\n" % (a, iy, iz, j)

                msg += "              ORIENTATION=_F( \n"
                ## .. todo:: is this correct
                msg += "                  CARA=('VECT_Y'), # direction of beam ???\n"
                msg += "                  VALE=(1.0,0.0,0.0,)"

                if [prop.n1a, prop.n1b] != [0., 0.]:
                    msg += "              \n),\n"
                    msg += "              CARA=('AX','AY'), # shear centers\n"
                    msg += "              VALE=(%g, %g),\n" % (prop.n1a, prop.n1b)
                    msg += "             ),\n"
                else:
                    msg += " )),\n"

            elif ptype == 'PROD':
                msg = "    POUTRE=_F(GROUP_MA='P%s', # PROD\n" % pid
                msg += "              SECTION='CERCLE',  # circular section\n"
                msg += "              CARA=('R')   # radius\n"
                #msg += "              VALE=(%g),),\n" % (prop.Radius())

                msg += "              SECTION='GENERALE',\n"
                msg += "              CARA=('A', 'JX')\n"
                msg += "              VALE=(%g, %g),\n"  %(prop.Area(), prop.J())
                msg += "                    CARA='VECT_Y'),\n"
                msg += "                    VALE=(1.0,0.0,0.0,),),\n"
            elif ptype == 'CELAS2': # TODO: not a property
                nodes = self.node_ids
                msg = 'DISCRET=_F( # CELAS2\n'
                if nodes[0]:
                    msg += "     CARA='K_T_D_N'\n"
                    msg += "     NOEUD=N%i,\n" % nodes[0]

                if nodes[1]:
                    msg += "     CARA='K_T_D_L'\n"
                    msg += "     NOEUD=N%i,\n" % nodes[1]
                    msg += "     AMOR_HYST=%g # ge - damping\n" % self.ge
                msg += "     )\n"
                msg += "\n"

                if self.c1 == 1:
                    msg += "VALE=(%g,0.,0.)\n" % self.k
                elif self.c1 == 2:
                    msg += "VALE=(0.,%g,0.)\n" % self.k
                elif self.c1 == 2:
                    msg += "VALE=(0.,0.,%g)\n" % self.k
                else:
                    raise ValueError('unsupported value of c1=%s' % self.c1)

            elif ptype == 'PSHELL':
                #
                # * http://www.caelinux.org/wiki/index.php/Contrib:KeesWouters/shell/static
                # * http://www.caelinux.org/wiki/index.php/Contrib:KeesWouters/platedynamics
                #
                # The angle_rep is a direction angle, use either angle(a,b) or
                # vecteur(x,y,z)
                # coque_ncou is the number of gauss nodes along the thickness, for
                # linear analysis one node is sufficient.
                #
                msg = "    COQUE=_F(GROUP_MA='P%s', # COQUE=PSHELL\n" % prop.pid
                msg += "              EPAIS=%g, # EPAIS=thickness\n" % prop.t
                msg += "              ANGL_REP=(0.,90.),  # ???\n"  #: .. todo:: what is this?
                #msg += "              VECTEUR=(1.0,0.0,0.0,)  #  local coordinate system\n"
                msg += "              EXCENTREMENT=%g,  # offset-Z1\n" % prop.z1
                msg += "              COQUE_NCOU=1,  # Number of Integration Layers\n"
                msg += "              CARA=('NSM'), # ???\n"  #: .. todo:: check
                msg += "              VALE=(%g),),\n" % prop.nsm
                return msg
            else:
                self.log.warning(msg)
                msg = '# skipping %s because write_code_aster is not implemented\n' % ptype
                continue
            comm += msg
        comm += ');\n'
        #comm += ');\nFINSF\n\n'
        comm += self.breaker()
        return comm, py_ca

    def ca_materials(self):
        """
        writes the material cards

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
        for mid, material in sorted(self.materials.items()):
            #comm += 'GROUP_MA name = %s_%s\n' % (material.type, mid)
            if material.type == 'MAT1':
                # MAT1
                pre = 'M%s = DEFI_MATRIAU(ELAS=_F(' % self.mid
                spaces = ' ' * len(pre)
                msg = '%sE=%g, # MAT1 mid=%s\n' % (pre, self.e, self.mid)
                #msg  = 'M%s = DEFI_MATRIAU(ELAS=_F( # MAT1\n' %(self.mid)
                #msg += spaces + 'E  =%g,\n'  %(self.e)
                msg += spaces + 'NU=%g,\n' % (self.nu)
                if '.' in '%g' % self.rho:
                    msg += spaces + 'RHO=%g));\n' % (self.rho)
                else:
                    msg += spaces + 'RHO=%.1f));\n' % (self.rho)
            else:
                raise NotImplementedError(material)
            comm += msg

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
                              AFFE=(_F(GROUP_MA=('P32','P33','P42','P43','P46','P47','P48','P49',
                                                 'P61','P62','P63','P64','P65','P74',
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
        for mid, material in sorted(self.materials.items()):
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

                # LOAD card
                out = write_code_aster_load(load, self, grid_word='N')
                if len(out) == 3:  # LOAD card
                    (commi, load_ids, load_types) = out
                    comm += commi
                else:  # FORCEx, MOMENTx, GRAV
                    #skipped_load_ids[(load.lid, load.type)] = out
                    comm += out

            #loadcase.
            #for ID, grav in sorted(self.gravs.items()):
            #    comm += grav.write_code_aster(mag)

        #for lid_load_type, commi in sorted(skipped_load_ids.items()):
            #comm += commi

        comm += self.breaker()
        return comm

    def ca_spcs(self):
        """creates SPC constraints"""
        comm = ''
        comm += '# ca_spcs\n'
        for subcase_id, subcase in self.subcases.items():
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
        """creates MPC constraints"""
        comm = ''
        comm += '# ca_mpcs\n'
        for subcase_id, subcase in self.subcases.items():
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
        """just a simple line break"""
        return '#-------------------------------------------------------------------------\n'

    def build_maxs(self):
        """
        gets the max node/element/property/material lengths in order to
        nicely format the text

        1000  -> 4
        10000 -> 5

        Now we can line up the node ids
        """
        self.max_nid_len = len(str(max(self.nodes)))
        self.max_eid_len = len(str(max(self.elements)))
        self.max_pid_len = len(str(max(self.properties)))
        self.max_mid_len = len(str(max(self.materials)))

    def write_as_code_aster(self, model):
        """writes a BDF as a code aster formatted file"""
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
            with open(model + '.comm', 'wb') as comm_file:
                print("writing fname=%s" % (model + '.comm'))
                comm_file.write(comm)

        #print(comm)
        #print(mail)
        #print(py_ca)
        if mail:
            print('mail')
            with open(model + '.mail', 'wb') as mail_file:
                print("writing fname=%s" % (model + '.mail'))
                mail_file.write(mail)

        if py_ca:
            print('py_ca')
            assert py_ca != ''
            with open(model + '.py', 'wb') as py_file:
                print("writing fname=%s" % (model + '.py'))
                py_file.write(py_ca)


def write_code_aster_load(load, model, grid_word='node'):
    """writes a BDF load card in CA format"""
    load_ids = load.get_load_ids()
    load_types = load.get_load_types()

    #msg = '# Loads\n'
    msg = ''
    (types_found, force_loads, moment_loads,
     force_constraints, moment_constraints,
     gravity_loads) = load.organize_loads(model)

    nids = []
    for nid in force_loads:
        nids.append(nid)
    for nid in moment_loads:
        nids.append(nid)

    if nids:
        msg += '# types_found = %s\n' % (list(types_found))
        msg += '# load_ids    = %s\n' % (load_ids)
        msg += "load_bc = AFFE_CHAR_MECA(MODELE=modmod,\n"
        #msg += "                        DDL_IMPO=(_F(GROUP_MA='Lleft',\n"
        msg += "                         FORCE_NODALE=(\n"

    #CHAR=AFFE_CHAR_MECA(MODELE=MODE,
    #             FORCE_NODALE=(
    #                     _F(NOEUD='N1',
    #                        FZ=-500.0),)

    spaces = "                           "
    for nid in sorted(nids):  # ,load in sorted(force_loads.items())
        msg += spaces + "_F(NOEUD='%s%s'," % (grid_word, nid)

        if nid in force_loads:
            force = force_loads[nid]
            if abs(force[0]) > 0.:
                msg += " FX=%s," % force[0]
            if abs(force[1]) > 0.:
                msg += " FY=%s," % force[1]
            if abs(force[2]) > 0.:
                msg += " FZ=%s," % force[2]

        if nid in moment_loads:
            moment = moment_loads[nid]
            if abs(moment[0]) > 0.:
                msg += " MX=%s," % moment[0]
            if abs(moment[1]) > 0.:
                msg += " MY=%s," % moment[1]
            if abs(moment[2]) > 0.:
                msg += " MZ=%s," % moment[2]
        #msg = msg[:-2]
        msg += '),\n'
        # finish the load

        #if moment in
        #msg += "                                   DX=0.0,\n"
        #msg += "                                   DY=0.0,\n"
        #msg += "                                   DZ=0.0,),\n"
        #msg += "                                _F(GROUP_MA='Lright',\n"
        #msg += "                                   DZ=0.0,),),\n"
    msg = msg[:-2]
    msg += ');\n'

    for gravity_load in gravity_loads:
        msg += 'CA_GRAVITY(%s);\n' % str(gravity_load)
    return msg, load_ids, load_types


def write_matrix(matrix):
    """
    was in NastranMatrix as write_code_aster
    assume set 1 = MAAX1,MAAX2, etc. and 100/n % on each
    """
    # for real combination
    comm = 'K_Mtx_AB=COMB_MATR_ASSE(COMB_R=(\n'
    comm += '    _F(MATR_ASSE = K_Mtx_A,COEF_R = 1.),\n'
    comm += '    _F(MATR_ASSE = K_Mtx_B,COEF_R = 1.)));\n'

    # for complex combination

    comm += "K_Mtx_AB=COMB_MATR_ASSE(COMB_C=(\n"
    comm += "_F(MATR_ASSE=K_Mtx_A,COEF_C=('RI',0.7,0.3,),)\n"
    comm += "_F(MATR_ASSE=K_Mtx_B,COEF_C=('RI',0.7,0.3,),),),);\n"
    comm = 'K_Mtx=ASSE_MATRICE(MATR_ELEM=ElMtx_K,NUME_DDL=%s,);'
    return comm


def write_conm2(elem):
    """writes a CONM2 card"""
    msg = "    DISCRET=_F(\n"
    msg += "             'CARA='M_T_D_N'\n"
    msg += "              NOEUD=N%s\n" % elem.Nid()
    msg += "              VALE=%g),\n" % elem.mass
    return msg

def code_aster_beam_section(prop, iface, istart, dims):
    """
    writes a PBARL/PBEAML cross-section

    ::

      ---msg1---
      H1 = 0.1
      W1 = 0.05

      ---msg2---
      Face_1 = geompy.MakeFaceHW(H1, W1, 1)
      geompy.addToStudy(Face_1, 'Face_1')

      ---msg---
      H1 = 0.1
      W1 = 0.05
      Face_1 = geompy.MakeFaceHW(H1, W1, 1)
      geompy.addToStudy(Face_1, 'Face_1')
    """
    msg1 = ''
    msg2 = 'Face_%s = geompy.MakeFaceHW(' % (iface + 1)
    for (i, dim) in enumerate(dims):
        msg1 += 'D%s = %s\n' % (istart + i, dim)
        msg2 += 'D%s,' % (istart + i)
    msg2 += '1)\n'
    msg2 += "geompy.addToStudy(Face_%i, 'Face_%i')\n" % (iface, iface)
    return msg1 + msg2



def main():
    """runs nastranToCodeAster"""
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

    for key, value in sorted(data.items()):
        print("%-12s = %r" % (key.strip('--'), value))

    bdf_filename = data['BDF_FILENAME']
    fname_base = os.path.splitext(bdf_filename)[0]

    model = CodeAsterConverter()
    model.read_bdf(bdf_filename, encoding='ascii')
    model.write_as_code_aster(fname_base)  # comm, py

if __name__ == '__main__':  # pragma: no cover
    main()
