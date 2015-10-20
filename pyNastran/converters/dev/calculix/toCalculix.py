from __future__ import print_function
from six import iteritems
from six.moves import zip
from collections import defaultdict

from numpy import array, zeros
from pyNastran.bdf.bdf import BDF, PBAR, PBARL, PBEAM, PBEAML, LOAD
from pyNastran.bdf.cards.loads.staticLoads import Force


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

    def get_elements_by_pid(self, element_ids=None):
        """
        builds a dictionary where the key is the property ID and the value
        is a list of element IDs
        """
        if element_ids is None:
            element_ids = self.elements.iterkeys()

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
        for eid, element in iteritems(self.elements):
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
            element_ids = self.elements.iterkeys()

        elems = defaultdict(list)
        for eid in element_ids:
            element = self.elements[eid]
            Type = element.type
            elems[Type].append(eid)
        return elems

    def get_properties_by_mid(self):
        """
        builds a dictionary where the key is the material ID and the value
        is a list of property IDs
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

    def calculix_nodes(self, f):
        """
        *NODE
         1, 0.000000, 0.000000, 0.000000
         2, 1.000000, 0.000000, 0.000000
         3, 1.000000, 1.000000, 0.000000
        """

        dat = ''
        dat += '** Calculix_Nodes\n'
        dat += '*NODE\n'
        f.write(dat)

        form = '%-' + str(self.maxNIDlen) + 's %8s,%8s,%8s\n'

        for nid, node in sorted(iteritems(self.nodes)):
            p = node.get_position()
            dat = form % (nid, p[0], p[1], p[2])
            f.write(dat)
        dat = '\n\n'
        dat += self.breaker()
        f.write(dat)

    def calculix_elements(self, f):
        """
        @todo sort elements by Type and Material ID
        """
        dat = ''
        dat += '** Calculix_Elements\n'

        TypeMap = {
            'CBAR'    : 'BR32R',
            'CBEAM'   : 'BR32R',

            'CTRIA3'  : 'C2D3',
            'CTRIA6'  : 'C2D6',
            'CQUAD4'  : 'C2D4',
            'CQUAD8'  : 'C2D8',

            'CTETRA'  : 'C3D8',
            'CTETRA4'  : 'C3D8',
            'CPENTA10' : 'C3D10',
            'CPENTA5'  : 'C3D5',
            'CPENTA15' : 'C3D15',
            'CHEXA8'   : 'C3D8',
            'CHEXA20'  : 'C3D20',
        }
        pid_eids = self.get_elements_by_pid(element_ids=None)
        formE = '%-' + str(self.nelements) + 's, '

        elsets = []
        for pid, eids in sorted(iteritems(pid_eids)):
            elems = self.get_elements_by_type(eids)
            for Type, eids in sorted(iteritems(elems)):
                calculix_type = TypeMap[Type]
                elset = 'pid%i_Elements%s' % (pid, Type)
                elsets.append(elset)

                dat += '** eid,n1,n2,n3,etc... for a %s\n' % Type
                dat += '*ELEMENT, TYPE=%s, ELSET=%s\n' % (calculix_type, elset)
                for eid in eids:
                    dat += formE % eid
                    element = self.elements[eid]
                    for nid in element.nodeIDs():
                        dat += '%s,' % nid
                    dat = dat[:-1] + '\n'
        dat += self.breaker()
        #print(dat)
        f.write(dat)
        return elsets

    def calculix_properties(self, elsets):
        inp = ''
        inp += '** calculix_properties\n'

        for elset in elsets:
            #elset = 'pid%i_Elements%i' % (pid, elType)
            pid, elType = elset.lstrip('pid').split('_')
            pid = int(pid)
            elType = elType[8:] # element type
            prop = self.properties[pid]

            if prop.type == 'PSHELL':
                mid = prop.mid
                inp += '*SHELL SECTION,ELSET=%s,MATERIA	L=MAT%i\n' % (elset, mid.Mid())
                inp +='%s\n' % prop.t
            elif prop.type == 'PSOLID':
                mid = prop.mid
                inp += '*SOLID SECTION,ELSET=%s,MATERIAL=MAT%i\n' % (elset, mid.Mid())
            elif prop.type == 'PBAR':
                mid = prop.mid
                inp += '*BEAM SECTION,ELSET=%s,MATERIAL=MAT%i\n' % (elset, mid.Mid())
            elif prop.type == 'PBARL':
                mid = prop.mid
                #section_name = 'SQUARE'
                print("what is the section_name?")
                print(" ", sorted(prop.__dict__.keys()))
                msg += '*BEAM SECTION,ELSET=eids_pid%i,MATERIAL=MAT%i,SECTION=%s\n' % (prop.Pid(), mid.Mid(), section_name)
                if section_name == 'SQUARE':
                    msg += '%s\n' % prop.dims[0]
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
        for mid, material in sorted(iteritems(self.materials)):
            msg = '*MATERIAL,NAME=MAT%i\n'
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
            loadKeys = self.loads.keys()
            #if isubcase in self.case_control_deck:
            if self.case_control_deck.has_subcase(isubcase):
                loadcase_id = self.case_control_deck.get_subcase_parameter(isubcase, param_name)[0]
                #loadcase = self.loads[loadcase_id]
                self._write_loads(loadcase_id)

        inp += self.breaker()
        return inp

    def _write_loads(self, loadcase_id, p0=None):
        if not isinstance(loadcase_id, int):
            raise RuntimeError('loadcase_id must be an integer; loadcase_id=%r' % loadcase_id)
        if p0 is None:
            p = array([0., 0., 0.], dtype='float32')
        if isinstance(p0, int):
            p = self.model.nodes[p0].get_position()
        else:
            p = array(p0)

        loadCase = self.loads[loadcase_id]
        #for (key, loadCase) in iteritems(self.loads):
            #if key != loadcase_id:
                #continue

        scale_factors2 = []
        loads2 = []
        for load in loadCase:
            if isinstance(load, LOAD):
                scale_factors, loads = load.get_reduced_loads()
                scale_factors2 += scale_factors
                loads2 += loads
            else:
                scale_factors2.append(1.)
                loads2.append(load)

        nnodes = self.nnodes
        print('nnodes = %s' % nnodes)
        FM = zeros((nnodes, 6), 'float64')
        print(FM.shape)
        F = FM[:, :3]
        M = FM[:, 3:]

        i = 0
        xyz = {}
        nid_to_i_map = {}
        for nid, node in iteritems(self.nodes):
            nid_to_i_map[nid] = i
            xyz[nid] = node.get_position()

        unsupported_types = set([])
        for load, scale in zip(loads2, scale_factors2):
            if isinstance(load, Force):  # FORCE, FORCE1, FORCE2
                f = load.mag * load.xyz
                print(load, load.node)
                i = nid_to_i_map[load.node]
                F[i, :] += f * scale
            elif isinstance(load, Moment):  # MOMENT, MOMENT1, MOMENT2
                m = load.mag * load.xyz
                i = nid_to_i_map[load.node]
                M[i, :] += m * scale
            elif load.type == 'PLOAD':
                nodes = load.nodeIDs()
                nnodes = len(nodes)
                if nnodes == 3:
                    n1, n2, n3 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]]
                    axb = cross(n1 - n2, n1 - n3)
                    centroid = (n1 + n2 + n3) / 3.
                elif nnodes == 4:
                    n1, n2, n3, n4 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]], xyz[nodes[3]]
                    axb = cross(n1 - n3, n2 - n4)
                    centroid = (n1 + n2 + n3 + n4) / 4.
                else:
                    raise RuntimeError('invalid number of nodes on PLOAD card; nodes=%s' % str(nodes))

                nunit = norm(axb)
                A = 0.5 * nunit
                try:
                    n = axb / nunit
                except FloatingPointError:
                    msg = ''
                    for i, nid in enumerate(nodes):
                        msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
                    msg += 'a x b = %s\n' % axb
                    msg += 'nunit = %s\n' % nunit
                    raise FloatingPointError(msg)
                f = load.p * A * n * scale / nnodes

                for nid in nodes:
                    i = nid_to_i_map[nid]
                    F[i, :] = f

            elif load.type == 'PLOAD1':
                elem = load.eid

            elif load.type == 'PLOAD2':
                pressure = load.pressures[0] * scale  # there are 4 pressures, but we assume p0
                for eid in load.eids:
                    elem = self.elements[eid]
                    if elem.type in ['CTRIA3',
                                     'CQUAD4', 'CSHEAR']:
                        nodes = elem.nodeIDs()
                        nnodes = len(nodes)
                        n = elem.Normal()
                        A = elem.Area()
                        f = pressure * n * A / nnodes
                        for nid in nodes:
                            i = nid_to_i_map[nid]
                            F[i, :] = f
                    else:
                        self.log.debug('case=%s etype=%r loadtype=%r not supported' % (loadcase_id, elem.type, load.type))
            elif load.type == 'PLOAD4':
                pressure = load.pressures[0] * scale  # there are 4 possible pressures, but we assume p0
                assert load.Cid() == 0, 'Cid() = %s' % (load.Cid())
                assert load.sorl == 'SURF', 'sorl = %s' % (load.sorl)
                assert load.ldir == 'NORM', 'ldir = %s' % (load.ldir)
                for elem in load.eids:
                    eid = elem.eid
                    if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',]:
                        # triangles
                        nnodes = 3
                        nodes = elem.nodeIDs()
                        n1, n2, n3 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]]
                        axb = cross(n1 - n2, n1 - n3)
                        nunit = norm(axb)
                        A = 0.5 * nunit
                        try:
                            n = axb / nunit
                        except FloatingPointError:
                            msg = ''
                            for i, nid in enumerate(nodes):
                                msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
                            msg += 'a x b = %s\n' % axb
                            msg += 'nunit = %s\n' % nunit
                            raise FloatingPointError(msg)
                        centroid = (n1 + n2 + n3) / 3.
                    elif elem.type in ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        # quads
                        nnodes = 4
                        nodes = elem.nodeIDs()
                        n1, n2, n3, n4 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]], xyz[nodes[3]]
                        axb = cross(n1 - n3, n2 - n4)
                        nunit = norm(axb)
                        A = 0.5 * nunit
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
                        A, centroid, normal = elem.getFaceAreaCentroidNormal(load.g34.nid, load.g1.nid)
                        nnodes = None
                    else:
                        self.log.debug('case=%s eid=%s etype=%r loadtype=%r not supported' % (loadcase_id, eid, elem.type, load.type))
                        continue
                    #r = centroid - p
                    f = pressure * A * n / nnodes
                    #m = cross(r, f)
                    for nid in nodes:
                        i = nid_to_i_map[nid]
                        F[i, :] = f
            elif load.type == 'GRAV':
                pass
            else:
                # we collect them so we only get one print
                unsupported_types.add(load.type)

        for Type in unsupported_types:
            self.log.debug('case=%s loadtype=%r not supported' % (loadcase_id, Type))
        FM.reshape((nnodes*6,1))
        return FM

    def calculix_constraints(self):
        return self.calculix_SPCs()

    def calculix_SPCs(self):
        #for spcID,spcs in iteritems(self.spcObject2):
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

    def write_as_calculix(self, fname='fem'):
        inp = ''
        dat = ''
        self.buildMaxs()  # gets number of nodes/elements/properties/materials

        inp += '** BEGIN BULK\n'
        #inp += 'DEBUT();\n\n'

        #inp += "**'Read the mesh' - we use the 'aster' file format here.\n"
        #inp += 'mesh=LIRE_MAILLAGE(UNITE=20,\n'
        #inp += "                   FORMAT='ASTER');\n\n"

        #inp += "**'MECA_STATIQUE' % SOL 101 - linear statics\n"
        #inp += "** Assigning the model for which CA will calculate the results:\n"
        #inp += "** 'Mecanique' - since we are dealing with a linear elastic model and '3D' since it's a 3D model.\n"
        #inp += 'Meca=AFFE_MODELE(MAILLAGE=mesh,\n'
        #inp += "                 AFFE=_F(TOUT='OUI',\n"
        #inp += "                         PHENOMENE='MECANIQUE',\n"
        #inp += "                         MODELISATION='3D',),);\n\n"
        inp += self.breaker()

        fdat = open(fname + '.dat', 'wb')
        finp = open(fname + '.inp', 'wb')

        print("writing fname=%s" % (fname + '.dat'))
        print("writing fname=%s" % (fname + '.inp'))

        self.calculix_nodes(fdat)

        elsets = self.calculix_elements(fdat)

        dat = self.calculix_materials()
        fdat.write(dat)
        fdat.close()

        inpi = self.calculix_properties(elsets)
        inp += inpi
        inp += self.calculix_loads()
        inp += self.calculix_constraints()
        finp.write(inp)

        # Case Control Deck
        inp =  '*STEP\n'
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
        finp.close()


def main():
    import sys
    ca = CalculixConverter()
    #bdf_filename = 'solidBending.bdf'
    bdf_filename = sys.argv[1]
    ca.read_bdf(bdf_filename, xref=False)
    ca.cross_reference(xref=True, xref_elements=True, xref_properties=True,
                      xref_materials=True,
                      xref_loads=True,
                      xref_constraints=True,
                      xref_aero=True)
    ca.write_as_calculix(bdf_filename + '.ca')  # inp, py

if __name__ == '__main__':  # pragma: no cover
    main()
