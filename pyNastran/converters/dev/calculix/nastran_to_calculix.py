"""
defines:
 - CalculixConverter

"""
from collections import defaultdict

from numpy import array, zeros, cross
from numpy.linalg import norm  # type: ignore
from pyNastran.bdf.bdf import BDF, LOAD   # PBAR, PBARL, PBEAM, PBEAML,
from pyNastran.bdf.cards.loads.static_loads import Force, Moment


class CalculixConverter(BDF):
    """
    Converts a BDF to Calculix (inp/dat/py files).

    .. warning:: Totally inaccurate....

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

    .. todo::
      * PCOMP
      * SPC, SPC1, MPC
      * RBE2, RBE3
    """
    def __init__(self, language='english'):
        self.language = 'english'
        BDF.__init__(self)
        self.max_nid_len = None
        self.max_eid_len = None
        self.max_pid_len = None
        self.max_mid_len = None

    def get_elements_by_pid(self, element_ids=None):
        """
        builds a dictionary where the key is the property ID and the value
        is a list of element IDs
        """
        if element_ids is None:
            element_ids = self.elements.keys()

        props = defaultdict(list)
        for eid in element_ids:
            element = self.elements[eid]
            pid = element.Pid()
            props[pid].append(eid)
        return props

    def get_elements_by_mid(self):
        """
        builds a dictionary where the key is the material ID and the value
        is a list of element IDs
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

    def get_elements_by_type(self, element_ids=None):
        """
        builds a dictionary where the key is the element type and the value
        is a list of element IDs
        """
        if element_ids is None:
            element_ids = self.elements.keys()

        elems = defaultdict(list)
        for eid in element_ids:
            element = self.elements[eid]
            element_type = element.type
            elems[element_type].append(eid)
        return elems

    def get_properties_by_mid(self):
        """
        builds a dictionary where the key is the material ID and the value
        is a list of property IDs
        """
        mats = {0: []}

        for mid in self.materials:
            mats[mid] = []
        for pid, property in self.properties.items():
            try:
                mid = property.Mid()
                mats[mid].append(pid)
            except:
                mats[0].append(pid)
        return mats

    def calculix_executive(self):
        inp = ''
        if self.sol == 101:
            inp += 'MECA_STATIQUE % SOL 101 - linear statics\n'
            inp += 'stat(MECA_STATIQUE(MODELE=model,CHAM_MATER=material,CARA_ELEM=elemcar,\n'

            inp += 'ECIT=(_F(Charge=AllBoundaryConditions,),\n',
            inp += '      _F(Charge=AllLoads,),\n',
            inp += '      ),\n',

            inp += "TITRE='My Title'\n"
        return inp

    def calculix_nodes(self, fdat):
        """
        *NODE
         1, 0.000000, 0.000000, 0.000000
         2, 1.000000, 0.000000, 0.000000
         3, 1.000000, 1.000000, 0.000000
        """

        dat = ''
        dat += '** Calculix_Nodes\n'
        dat += '*NODE\n'
        fdat.write(dat)

        form = '%-' + str(self.max_nid_len) + 's %8s,%8s,%8s\n'

        for nid, node in sorted(self.nodes.items()):
            xyz = node.get_position()
            dat = form % (nid, xyz[0], xyz[1], xyz[2])
            fdat.write(dat)
        dat = '\n\n'
        dat += self.breaker()
        fdat.write(dat)

    def calculix_elements(self, fdat):
        """
        .. todo:: sort elements by Type and Material ID
        """
        dat = ''
        dat += '** Calculix_Elements\n'

        etype_map = {
            'CBAR'    : 'BR32R',
            'CBEAM'   : 'BR32R',

            'CTRIA3'  : 'C2D3',
            'CTRIA6'  : 'C2D6', # 'S6' ???
            'CQUAD4'  : 'C2D4', # 'S4' ???
            'CQUAD8'  : 'C2D8',
            'CSHEAR' : 'S4',

            'CTRIAX' : 'CAX6',
            'CQUADX' : 'CAX8',
            'CTRIAX6' : 'CAX6',
            'CQUADR' : 'CAX8',

            'CTETRA'  : 'C3D4',
            'CTETRA4'  : 'C3D4',
            'CPYRAM'  : 'C3D5',
            'CPYRAM5'  : 'C3D5',
            'CPYRAM13'  : 'C3D13',
            'CPENTA10' : 'C3D10',
            'CPENTA6'  : 'C3D6',
            'CPENTA15' : 'C3D15',
            'CHEXA'   : 'C3D8',
            'CHEXA8'   : 'C3D8',
            'CHEXA20'  : 'C3D20',
        }
        pid_eids = self.get_elements_by_pid(element_ids=None)
        form_elements = '%-' + str(self.nelements) + 's, '

        elsets = []
        for pid, eids in sorted(pid_eids.items()):
            elems = self.get_elements_by_type(eids)
            for etype, eids in sorted(elems.items()):
                calculix_type = etype_map[etype]
                elset = 'pid%i_Elements%s' % (pid, etype)
                elsets.append(elset)

                dat += '** eid,n1,n2,n3,etc... for a %s\n' % etype
                dat += '*ELEMENT, TYPE=%s, ELSET=%s\n' % (calculix_type, elset)
                for eid in eids:
                    dat += form_elements % eid
                    element = self.elements[eid]
                    for nid in element.node_ids:
                        dat += '%s,' % nid
                    dat = dat[:-1] + '\n'
        dat += self.breaker()
        #print(dat)
        fdat.write(dat)
        return elsets

    def calculix_properties(self, elsets):
        inp = ''
        inp += '** calculix_properties\n'

        for elset in elsets:
            #elset = 'pid%i_Elements%i' % (pid, etype)
            pid, etype = elset.lstrip('pid').split('_')
            pid = int(pid)
            etype = etype[8:] # element type
            prop = self.properties[pid]

            if prop.type == 'PSHELL':
                mid = prop.mid
                inp += '*SHELL SECTION,ELSET=%s,MATERIAL=MAT%i\n' % (elset, mid.mid)
                inp += '%s\n' % prop.t
                #def _write_calculix(self, marker='markerDummyProp',
                                    #element_set='ELsetDummyProp'):
                    #msg = '*SHELL SECTION,MATERIAL=M%s_%s,ELSET=%s,OFFSET=%s\n' % (
                        #marker, self.mid, element_set, self.z1)
                    #msg += '** THICKNESS\n'
                    #msg += '%s\n\n' % (self.t)
                    #return msg
            elif prop.type == 'PSOLID':
                mid = prop.mid
                inp += '*SOLID SECTION,ELSET=%s,MATERIAL=MAT%i\n' % (elset, mid.mid)
            elif prop.type == 'PBAR':
                mid = prop.mid
                inp += '*BEAM SECTION,ELSET=%s,MATERIAL=MAT%i\n' % (elset, mid.mid)
            elif prop.type == 'PBARL':
                mid = prop.mid
                #section_name = 'SQUARE'
                print("what is the section_name?")
                print(" ", sorted(prop.__dict__.keys()))
                inp += '*BEAM SECTION,ELSET=eids_pid%i,MATERIAL=MAT%i,SECTION=%s\n' % (
                    prop.pid, mid.mid, section_name)
                if section_name == 'SQUARE':
                    inp += '%s\n' % prop.dims[0]
                if section_name == 'RECT':
                    inp += '%s, %s\n' % tuple(prop.dims[0])
                else:
                    raise NotImplementedError(section_name)
            else:
                raise NotImplementedError(section_name)

        inp += self.breaker()
        return inp

    def calculix_materials(self):
        """
        might need to make this by pid instead...
        steel=DEFI_MATERIAU(ELAS=_F(E=210000.,NU=0.3,RHO=8e-9),);

        -----MAT1-----
        *MATERIAL,NAME=EL
         210000.0, .3
        *DENSITY
        7.8E-9
        *SOLID SECTION,MATERIAL=EL,ELSET=EALL
        """
        inp = '** calculix_materials\n'
        for mid, material in sorted(self.materials.items()):
            msg = '*MATERIAL,NAME=mid%i\n' % material.mid
            if material.type == 'MAT1':
                msg += '*ELASTIC\n%s, %s\n' % (material.E(), material.Nu())
                msg += '*DENSITY\n%s\n' % material.get_density()
                msg += '*SOLID SECTION,MATERIAL=EL,ELSET=EALL\n'
            elif material.type == 'MAT4':
                msg += '*ELASTIC\n%s, %s\n' % (material.E(), material.Nu())
                msg += '*DENSITY\n%s\n' % material.rho()
                msg += '*CONDUCTIVITY\n%s\n' % material.k
                msg += '*CONVECTION\n%s\n' % material.h
                msg += '*DENSITY\n%s\n' % material.get_density()
                msg += '*SOLID SECTION,MATERIAL=EL,ELSET=EALL\n'
            elif material.type == 'MAT2':
                raise NotImplementedError(material)
                #msg = '*ELASTIC,TYPE=ORTHO\n'
                #temperature = 0.  # default value - same properties for all values
                #msg += '%s,%s,%s\n' % (self.e, self.nu, temperature)
                #D = Dplate
                #D1111 = D[0, 0]
                #D1122 = 0.
                #D2222 = D[1, 1]
                #D1133 = D[0, 2]
                #D2233 = D[1, 2]
                #D3333 = D[2, 2]
                #D1212 = D[0, 1]
                #D1313 = D[0, 2]
                #msg += '%s,%s,%s,%s,%s,%s,%s,%s\n\n' % (
                    #D1111, D1122, D2222, D1133, D2233, D3333, D1212, D1313)

                #G23
                #temperature = self.tref
                #msg = '*ELASTIC,TYPE=ENGINEERING CONSTANTS  ** MAT2,mid=%s\n' % (
                    #self.mid)
                #msg += '** E1,E2,E3,NU12,NU13,NU23,G12,G13\n'
                #msg += '** G23,TEMPERATURE\n'
                #msg += '%s,%s,%s,%s,%s,%s,%s,%s\n' % (
                    #e1, e2, e3, nu12, nu13, nu23, g12, g13)
                #msg += '%s,%s\n' % (G23, temperature)
                #if self.rho > 0.:
                    #msg += '*DENSITY\n'
                    #msg += '%s\n' % (self.rho)
                #if self.a > 0:
                    #msg += '*EXPANSION,TYPE=ISO,ZERO=%s\n' % (self.tref)
                    #msg += '** ALPHA,ALPHA*TREF\n'
                    #msg += '%s,%s\n\n' % (self.a, self.a * self.tref)
                #return msg
            else:
                raise NotImplementedError(mid.type)
        inp += self.breaker()
        return inp

    def calculix_loads(self):
        """writes the load cards sorted by ID"""
        inp = '** calculix_loads\n'
        #if self.language=='english':
            #inp += '** Loads\n'
        #else:
            #inp += ''

        isubcase = 1
        param_name = 'LOAD'

        #skippedLids = {}
        if self.loads:
            inp += '** LOADS\n'
            #load_keys = self.loads.keys()
            #if isubcase in self.case_control_deck:
            if self.case_control_deck.has_subcase(isubcase):
                loadcase_id = self.case_control_deck.get_subcase_parameter(
                    isubcase, param_name)[0]
                #loadcase = self.loads[loadcase_id]
                self._write_loads_p0(loadcase_id) # bdf_file, size=8, is_double=False

        inp += self.breaker()
        return inp

    def _write_loads_p0(self, loadcase_id, p0=None):
        if not isinstance(loadcase_id, int):
            raise RuntimeError('loadcase_id must be an integer; loadcase_id=%r' % loadcase_id)
        if p0 is None:
            p = array([0., 0., 0.], dtype='float32')
        if isinstance(p0, int):
            p = self.nodes[p0].get_position()
        else:
            p = array(p0)

        load_case = self.loads[loadcase_id]
        #for (key, load_case) in self.loads.items():
            #if key != loadcase_id:
                #continue

        scale_factors2 = []
        loads2 = []
        for load in load_case:
            if isinstance(load, LOAD):
                scale_factors, loads = load.get_reduced_loads(
                    resolve_load_card=False, filter_zero_scale_factors=False)
                scale_factors2 += scale_factors
                loads2 += loads
            else:
                scale_factors2.append(1.)
                loads2.append(load)

        nnodes = self.nnodes
        print('nnodes = %s' % nnodes)
        force_moment = zeros((nnodes, 6), 'float64')
        #print(force_moment.shape)
        force = force_moment[:, :3]
        moment = force_moment[:, 3:]

        i = 0
        xyz = {}
        nid_to_i_map = {}
        for nid, node in self.nodes.items():
            nid_to_i_map[nid] = i
            xyz[nid] = node.get_position()

        unsupported_types = set()
        for load, scale in zip(loads2, scale_factors2):
            if isinstance(load, Force):  # FORCE, FORCE1, FORCE2
                forcei = load.mag * load.xyz
                self.log.info('%s %s' % (load, load.node))
                i = nid_to_i_map[load.node]
                force[i, :] += forcei * scale
            elif isinstance(load, Moment):  # MOMENT, MOMENT1, MOMENT2
                momenti = load.mag * load.xyz
                i = nid_to_i_map[load.node]
                moment[i, :] += momenti * scale
            elif load.type == 'PLOAD':
                nodes = load.node_ids
                nnodes = len(nodes)
                if nnodes == 3:
                    n1, n2, n3 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]]
                    axb = cross(n1 - n2, n1 - n3)
                    #centroid = (n1 + n2 + n3) / 3.
                elif nnodes == 4:
                    n1, n2, n3, n4 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]], xyz[nodes[3]]
                    axb = cross(n1 - n3, n2 - n4)
                    #centroid = (n1 + n2 + n3 + n4) / 4.
                else:
                    msg = 'invalid number of nodes on PLOAD card; nodes=%s' % str(nodes)
                    raise RuntimeError(msg)

                nunit = norm(axb)
                area = 0.5 * nunit
                try:
                    n = axb / nunit
                except FloatingPointError:
                    msg = ''
                    for i, nid in enumerate(nodes):
                        msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
                    msg += 'a x b = %s\n' % axb
                    msg += 'nunit = %s\n' % nunit
                    raise FloatingPointError(msg)
                forcei = load.pressure * area * n * scale / nnodes

                for nid in nodes:
                    i = nid_to_i_map[nid]
                    force[i, :] = forcei

            elif load.type == 'PLOAD1':
                elem = load.eid

            elif load.type == 'PLOAD2':
                # there are 4 pressures, but we assume p0
                pressure = load.pressures[0] * scale
                for eid in load.eids:
                    elem = self.elements[eid]
                    if elem.type in ['CTRIA3',
                                     'CQUAD4', 'CSHEAR']:
                        nodes = elem.node_ids
                        nnodes = len(nodes)
                        n = elem.Normal()
                        area = elem.Area()
                        forcei = pressure * n * area / nnodes
                        for nid in nodes:
                            i = nid_to_i_map[nid]
                            force[i, :] = forcei
                    else:
                        self.log.debug('case=%s etype=%r loadtype=%r not supported' % (
                            loadcase_id, elem.type, load.type))
            elif load.type == 'PLOAD4':
                # there are 4 possible pressures, but we assume p0
                pressure = load.pressures[0] * scale
                assert load.Cid() == 0, 'Cid() = %s' % (load.Cid())
                assert load.surf_or_line == 'SURF', 'surf_or_line = %s' % (load.surf_or_line)
                assert load.line_load_dir == 'NORM', 'line_load_dir = %s' % (load.line_load_dir)
                for elem in load.eids:
                    eid = elem.eid
                    if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIAR',]:
                        # triangles
                        nnodes = 3
                        nodes = elem.node_ids
                        n1, n2, n3 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]]
                        axb = cross(n1 - n2, n1 - n3)
                        nunit = norm(axb)
                        area = 0.5 * nunit
                        try:
                            n = axb / nunit
                        except FloatingPointError:
                            msg = ''
                            for i, nid in enumerate(nodes):
                                msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
                            msg += 'a x b = %s\n' % axb
                            msg += 'nunit = %s\n' % nunit
                            raise FloatingPointError(msg)
                        #centroid = (n1 + n2 + n3) / 3.
                    elif elem.type in ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        # quads
                        nnodes = 4
                        nodes = elem.node_ids
                        n1, n2, n3, n4 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]], xyz[nodes[3]]
                        axb = cross(n1 - n3, n2 - n4)
                        nunit = norm(axb)
                        area = 0.5 * nunit
                        try:
                            n = axb / nunit
                        except FloatingPointError:
                            msg = ''
                            for i, nid in enumerate(nodes):
                                msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
                            msg += 'a x b = %s\n' % axb
                            msg += 'nunit = %s\n' % nunit
                            raise FloatingPointError(msg)

                        centroid = (n1 + n2 + n3 + n4) / 4.
                    elif elem.type in ['CTETRA', 'CHEXA', 'CPENTA']:
                        area, centroid, normal = elem.get_face_area_centroid_normal(
                            load.g34_ref.nid, load.g1_ref.nid)
                        nnodes = None
                    else:
                        self.log.debug('case=%s eid=%s etype=%r loadtype=%r not supported' % (
                            loadcase_id, eid, elem.type, load.type))
                        continue
                    #r = centroid - p
                    forcei = pressure * area * n / nnodes
                    #m = cross(r, f)
                    for nid in nodes:
                        i = nid_to_i_map[nid]
                        force[i, :] = forcei
            elif load.type == 'GRAV':
                pass
                #def write_calculix_grav(self, gx, gy, gz):
                    #msg = '*DLOAD\n'
                    #msg += 'AllElements,GRAV,%s,%s,%s\n' % (gx, gy, gz)
                    #return msg

            else:
                # we collect them so we only get one print
                unsupported_types.add(load.type)

        for load_type in unsupported_types:
            self.log.debug('case=%s load_type=%r not supported' % (loadcase_id, load_type))
        force_moment.reshape((nnodes*6, 1))
        return force_moment

    def calculix_constraints(self):
        return self.calculix_spcs()

    def calculix_spcs(self):
        #for spc_id, spcs in self.spcObject2.items():
        inp = ''
        inp += '** Calculix_SPCs\n'
        inp += self.breaker()
        return inp

    def breaker(self):
        return '**-------------------------------------------------------------------------\n'

    def build_maxs(self):
        self.max_nid_len = len(str(max(self.nodes)))
        self.max_eid_len = len(str(max(self.elements)))
        self.max_pid_len = len(str(max(self.properties)))
        self.max_mid_len = len(str(max(self.materials)))

    def write_as_calculix(self, fname='fem'):
        inp = ''
        dat = ''
        self.build_maxs()  # gets number of nodes/elements/properties/materials

        inp += '** BEGIN BULK\n'
        #inp += 'DEBUT();\n\n'

        inp += (
            "**'Read the mesh' - we use the 'aster' file format here.\n"
            'mesh=LIRE_MAILLAGE(UNITE=20,\n'
            "                   FORMAT='ASTER');\n\n"

            "**'MECA_STATIQUE' % SOL 101 - linear statics\n"
            "** Assigning the model for which CA will calculate the results:\n"
            "** 'Mecanique' - since we are dealing with a linear elastic model "
            "and '3D' since it's a 3D model.\n"
            'Meca=AFFE_MODELE(MAILLAGE=mesh,\n'
            "                 AFFE=_F(TOUT='OUI',\n"
            "                         PHENOMENE='MECANIQUE',\n"
            "                         MODELISATION='3D',),);\n\n")
        inp += self.breaker()

        print("writing fname=%s" % (fname + '.dat'))
        with open(fname + '.dat', 'wb') as fdat:
            self.calculix_nodes(fdat)
            elsets = self.calculix_elements(fdat)
            dat = self.calculix_materials()
            fdat.write(dat)

        print("writing fname=%s" % (fname + '.inp'))

        with open(fname + '.inp', 'wb') as finp:
            inpi = self.calculix_properties(elsets)
            inp += inpi
            inp += self.calculix_loads()
            inp += self.calculix_constraints()
            finp.write(inp)

            # Case Control Deck
            inp = '*STEP\n'
            inp += '*STATIC\n'
            inp += '*CLOAD\n'
            inp += 'LAST,2,1.\n'
            inp += '*NODE PRINT,NSET=NALL\n'
            inp += 'U\n'
            inp += '*EL PRINT,ELSET=EALL\n'
            inp += 'S\n'
            inp += '*END STEP\n'
            inp += '** END OF DATA\n'
            finp.write(inp)

def _write_mat1(material, element_set='ELSetDummyMat'):
    # default value - same properties for all values
    temperature = material.tref
    msg = '*ELASTIC,TYPE=ISO,ELSET=%s\n' % (element_set)
    msg += '** E,NU,TEMPERATURE\n'
    msg += '%s,%s,%s\n' % (material.e, material.nu, temperature)

    if material.rho > 0.:
        msg += '*DENSITY\n'
        msg += '%s\n' % (material.rho)
    if material.a > 0:
        msg += '*EXPANSION,TYPE=ISO,ZERO=%s\n' % (material.tref)
        msg += '** ALPHA,ALPHA*TREF\n'
        msg += '%s,%s\n\n' % (material.a, material.a * material.tref)
    return msg

def main():  # pragma: no cover
    import sys
    code_aster = CalculixConverter()
    #bdf_filename = 'solidBending.bdf'
    bdf_filename = sys.argv[1]
    code_aster.read_bdf(bdf_filename, xref=False)
    code_aster.cross_reference()
    code_aster.write_as_calculix(bdf_filename + '.ca')  # inp, py

if __name__ == '__main__':  # pragma: no cover
    main()
