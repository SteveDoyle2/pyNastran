from __future__ import annotations
import getpass
from typing import TYPE_CHECKING
from numpy import unique

from pyNastran.utils import object_attributes
from pyNastran.bdf.cards.base_card import deprecated
from pyNastran.bdf.case_control_deck import CaseControlDeck

from pyNastran.op2.op2_interface.internal_stats import get_op2_stats
from pyNastran.op2.result_objects.op2_results import Results

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2
    from cpylog import SimpleLogger

USER = getpass.getuser()

class OP2_F06_Common:
    def __init__(self):
        #: a dictionary that maps an integer of the subcaseName to the
        #: subcase_id
        self.isubcase_name_map = {}
        self.generalized_tables = {}
        self.case_control_deck = CaseControlDeck([], log=self.log)
        self.expected_times = {}
        self.params = {}
        self.table_names = []

        self.make_geom = False

        #: BDF Title
        self.title = None

        self.page_num = 1

        self.isubcases = []
        self.__objects_vector_init__()
        self.__objects_init__()
        self.__objects_common_init__()

    def has_result(self, result_name: str) -> bool:
        """checks to see if a result exists"""
        if '.' in result_name:
            sline = result_name.split('.')
            if len(sline) != 2:
                msg = 'result_name=%r has too many dots; only 2 are allowed' % result_name
                raise RuntimeError(msg)

            #obj_names = sline[:-1]
            #storage_obj = self
            #for obj_name in obj_names:
                #print(obj_name)
                #storage_obj = getattr(storage_obj, obj_name)

            obj_name, result_name = result_name.split('.')
            try:
                storage_obj = getattr(self.op2_results, obj_name)
            except AttributeError:
                return False
            return hasattr(storage_obj, result_name)
        else:
            return hasattr(self, result_name)

    def get_result(self, result_name: str) -> dict:
        """
        Getattr, but considers sub-objects

        Examples
        --------
        **Example 1**
        >>> self.eigenvectors = self.get_result('eigenvectors')

        **Example 2**
        >>> self.ato.displacements = self.get_result('ato.displacements')

        """
        if '.' in result_name:
            sline = result_name.split('.')
            if len(sline) != 2:
                raise RuntimeError('result_name=%r has too many dots; '
                                   'only 2 are allowed' % result_name)

            obj_name, result_name = result_name.split('.')

            storage_obj = getattr(self.op2_results, obj_name)
            storage_dict = getattr(storage_obj, result_name)
            return storage_dict
        else:
            try:
                storage_obj = getattr(self, result_name)
            except AttributeError:
                storage_obj = getattr(self.op2_results, result_name)
            return storage_obj

    def del_result(self, result_name: str) -> None:
        """
        delattr, but considers sub-objects
        """
        if '.' in result_name:
            sline = result_name.split('.')
            if len(sline) != 2:
                raise RuntimeError('result_name=%r has too many dots; '
                                   'only 2 are allowed' % result_name)
            obj_name, result_name = result_name.split('.')

            storage_obj = getattr(self.op2_results, obj_name)
            delattr(storage_obj, result_name)
        else:
            try:
                delattr(self, result_name)
            except AttributeError:
                delattr(self.op2_results, result_name)

    def deprecated(self, old_name: str, new_name: str, deprecated_version: str):
        """allows for simple OP2 vectorization"""
        return deprecated(old_name, new_name, deprecated_version, levels=[0, 1, 2])

    # ------------------------------------------------------------------
    # stress

    # ------------------------------------------------------------------
    # Strain Energy - Getter

    # ------------------------------------------------------
    # Strain Energy - Getter 2
    # ------------------------------------------------------------------
    # Strain Energy - Setter



    # ------------------------------------------------------------------
    # Stress - Getter

    #-------------------------------------------------------------------
    # Strain - Getter

    # ------------------------------------------------------------------
    # Force - Getter

    # ------------------------------------------------------------------
    # Force - Setter

        #(model.cquad4_strain_energy, 'CQUAD4', True),
        #(model.cquad8_strain_energy, 'CQUAD8', True),
        #(model.cquadr_strain_energy, 'CQUADR', True),
        #(model.cquadx_strain_energy, 'CQUADX', True),

        #(model.ctria3_strain_energy, 'CTRIA3', True),
        #(model.ctria6_strain_energy, 'CTRIA6', True),
        #(model.ctriar_strain_energy, 'CTRIAR', True),
        #(model.ctriax_strain_energy, 'CTRIAX', True),
        #(model.ctriax6_strain_energy, 'CTRIAX6', True),

        #(model.ctetra_strain_energy, 'CTETRA', True),
        #(model.cpenta_strain_energy, 'CPENTA', True),
        #(model.chexa_strain_energy, 'CHEXA', True),
        #(model.cpyram_strain_energy, 'CPYRAM', True),

        #(model.crod_strain_energy, 'CROD', True),
        #(model.ctube_strain_energy, 'CTUBE', True),
        #(model.conrod_strain_energy, 'CONROD', True),

        #(model.cbar_strain_energy, 'CBAR', True),
        #(model.cbeam_strain_energy, 'CBEAM', True),

        #(model.cgap_strain_energy, 'CGAP', True),
        #(model.celas1_strain_energy, 'CELAS1', True),
        #(model.celas2_strain_energy, 'CELAS2', True),
        #(model.celas3_strain_energy, 'CELAS3', True),
        #(model.celas4_strain_energy, 'CELAS4', True),
        #(model.cdum8_strain_energy, 'CDUM8', False),
        #(model.cbush_strain_energy, 'CBUSH', True),
        ##(model.chexa8fd_strain_energy, '', False),
        #(model.cbend_strain_energy, 'CBEND', False),
        #(model.dmig_strain_energy, 'DMIG', False),
        #(model.genel_strain_energy, 'GENEL', False),
        #(model.cshear_strain_energy, 'CSHEAR', True),
        #(model.conm2_strain_energy, 'CONM2', False),


    # --------------------------------------------------
    # Force - Setter

    # ------------------------------------------------------------------
    def __objects_vector_init__(self):
        """
        All OUG table is simple to vectorize, so we declere it in __objects_init__
        On the other hand, the rodForces object contains CROD/CTUBE/CONROD elements.
        It is difficult to handle initializing the CRODs/CONRODs given a
        mixed type case, so we split out the elements.

        """
        self.matrices = {}
        self.matdicts = {}
        #======================================================================
        # rods
        self.op2_results = Results()

        self.nonlinear_ctetra_stress_strain = {}
        self.nonlinear_cpenta_stress_strain = {}
        self.nonlinear_chexa_stress_strain = {}

        #======================================================================
        #: OES - CBEAM 94
        self.nonlinear_cbeam_stress = {}

        self.nonlinear_cbush_force_stress_strain = {}  # CBUSH 226
        self.nonlinear_cbush1d_stress_strain = {}

        #======================================================================

    def __objects_common_init__(self):
        #: the date the job was run on
        self.date = None

        #: Grid Point Weight Table
        #: create with:
        #:   PARAM   GRDPNT    0  (required for F06/OP2)
        #:   PARAM   POSTEXT YES  (required for OP2)
        self.grid_point_weight = {}
        self.oload_resultant = None

        #: LAMA
        self.eigenvalues = {}  # LAMA, LAMAS
        self.eigenvalues_fluid = {}  # LAMAF

        #self.convergence_history = {}
        #self.response1_table = {}

    def __objects_init__(self):
        """More variable declarations"""
        #: the date the job was run on
        self.date = None

        #: Grid Point Weight Table
        #: create with:
        #:   PARAM   GRDPNT    0  (required for F06/OP2)
        #:   PARAM   POSTEXT YES  (required for OP2)
        self.grid_point_weight = {}

        #: self.frequencies already exists as a BDF object
        #: but we need this for the FOL frequencies for the MONPNT1 and MONPNT3
        self._frequencies = None

        #: LAMA
        self.eigenvalues = {}  # LAMA, CLAMA, BLAMA, LAMAS
        self.eigenvalues_fluid = {}  # LAMAF

        #: OUG - displacement
        self.displacements = {}           # tCode=1 thermal=0

        #: OUG - velocity
        self.velocities = {}              # tCode=10 thermal=0

        #: OUG - acceleration
        self.accelerations = {}            # tCode=11 thermal=0

        #: OUG - temperatures
        self.temperatures = {}           # tCode=1 thermal=1

        #: OUG - eigenvectors
        self.eigenvectors = {}            # tCode=7 thermal=0
        self.eigenvectors_structure = {}  #  OUG1S
        self.eigenvectors_fluid = {}  #  OUG1F???

        # OES - tCode=5 thermal=0 s_code=0,1 (stress/strain)

        #: OES - nonlinear CROD/CONROD/CTUBE stress/strain
        self.nonlinear_crod_stress = {}
        self.nonlinear_crod_strain = {}

        self.nonlinear_ctube_stress = {}
        self.nonlinear_ctube_strain = {}

        self.nonlinear_conrod_stress = {}
        self.nonlinear_conrod_strain = {}

        self.nonlinear_cquad4_stress = {}
        self.nonlinear_ctria3_stress = {}

        self.nonlinear_cquad4_strain = {}
        self.nonlinear_ctria3_strain = {}

        #: OES - CELAS1 224, CELAS3 225,
        self.nonlinear_celas1_stress = {}
        self.nonlinear_celas3_stress = {}

        #: OES - GAPNL 86
        self.nonlinear_cgap_stress = {}

        self.nonlinear_cpyram_stress_strain = {}

        # OQG - spc/mpc forces
        self.spc_forces = {}  # OQG1, tCode=3?
        self.spc_forces_v = {} # OQGV1

        self.mpc_forces = {}  # tCode=39

        self.contact_forces = {} # OQGCF1
        self.contact_tractions_and_pressure = {}  # OBC1

        self.glue_forces = {} # OQGGF1

        # OQG - thermal forces
        self.thermal_gradient_and_flux = {}

        #: OGF - grid point forces
        self.grid_point_forces = {}  # tCode=19

        #: OGS1 - grid point stresses
        self.grid_point_surface_stresses = {}       # tCode=26
        self.grid_point_stresses_volume_direct = {}  # tCode=27
        self.grid_point_stresses_volume_principal = {} # tCode=28
        self.grid_point_stress_discontinuities = {}  # tCode=35

        # : OGSTR1 - grid point strains
        self.grid_point_surface_strains = {}       # tCode=26
        self.grid_point_strains_volume_direct = {}  # tCode=27
        self.grid_point_strains_volume_principal = {} # tCode=28
        self.grid_point_strain_discontinuities = {}  # tCode=35

        #: OPG - summation of loads for each element
        self.load_vectors = {}       # OPG1; tCode=2  thermal=0
        self.load_vectors_v = {}     # OPGV1
        self.thermal_load_vectors = {}  # tCode=2  thermal=1
        #self.applied_loads = {}       # tCode=19 thermal=0
        self.force_vectors = {}       # tCode=12 thermal=0

    def remove_empty_results(self) -> None:
        table_names = self.get_table_types() + ['grid_point_weight', 'oload_resultant']
        for table_name in table_names:
            obj = self.get_result(table_name)
            if obj is None:
                self.del_result(table_name)
            elif isinstance(obj, dict):
                if len(obj) == 0:
                    self.del_result(table_name)
                    #obj = None
                    #delattr(self, table_name)
                else:
                    print('saving %s dict' % table_name)
            else:
                print('saving %s' % table_name)

    def get_table_types(self):
        """Gets the names of the results."""
        base = self.op2_results.get_table_types()

        table_types = base + [
            'grid_point_weight',

            # OUG - displacement, temperatures, eigenvectors, velocity, acceleration
            'displacements',  # 'displacements_scaled',
            'temperatures',
            'eigenvectors', 'eigenvectors_structure', 'eigenvectors_fluid',
            'velocities',
            'accelerations',

            # OQG - spc/mpc forces
            'spc_forces', 'spc_forces_v',
            'mpc_forces',
            'contact_forces', 'contact_tractions_and_pressure',
            'glue_forces',
            'thermal_gradient_and_flux',

            # OGF - grid point forces
            'grid_point_forces',

            # OPG - summation of loads for each element '
            'load_vectors', 'load_vectors_v',
            'thermal_load_vectors',
            #'applied_loads',
            'force_vectors',
            # OES - isotropic CBEAM stress/strain
            'nonlinear_cbeam_stress',
            #'nonlinear_cbeam_strain',
            'nonlinear_cpyram_stress_strain',
        ]
        utables = unique(table_types)
        if len(table_types) != len(utables):
            msg = 'Non-unique tables: '
            for i, table_type in enumerate(table_types):
                if table_type in table_types[i+1:]:
                    msg += table_type + ', '
            raise AssertionError(msg)

        table_types += [
            # PVT/PVT0
            'params',

            # LAMA
            'eigenvalues',  # LAMA, CLAMA, BLAMA, LAMAS
            'eigenvalues_fluid',  # LAMAF

            # OUG1S
            #'eigenvectors_structure',

            # HISADD
            #'convergence_history',

            # R1TABRG
            #'response1_table',
        ]
        table_types += [
            #'cbush_stress',
            #'cbush_strain',
            #'cbush1d_stress_strain',

            # OES - nonlinear CROD/CONROD/CTUBE stress
            'nonlinear_crod_stress', 'nonlinear_crod_strain',
            'nonlinear_ctube_stress', 'nonlinear_ctube_strain',
            'nonlinear_conrod_stress', 'nonlinear_conrod_strain',

            # OESNLXR - CTRIA3/CQUAD4 stress
            'nonlinear_cbush_force_stress_strain',
            'nonlinear_cbush1d_stress_strain',
            'nonlinear_cgap_stress',
            'nonlinear_cquad4_stress', 'nonlinear_ctria3_stress',
            'nonlinear_cquad4_strain', 'nonlinear_ctria3_strain',

            # OESNLXR - solid
            'nonlinear_ctetra_stress_strain', 'nonlinear_cpenta_stress_strain', 'nonlinear_chexa_stress_strain',

            # OES - CEALS1 224, CELAS3 225
            'nonlinear_celas1_stress',
            'nonlinear_celas3_stress',

            # OGS1 - grid point stresses
            'grid_point_surface_stresses', # tCode=26
            'grid_point_stresses_volume_direct',  # tCode=27 # volume direct
            'grid_point_stresses_volume_principal', # tCode =28
            'grid_point_stress_discontinuities',  # tCode=35,

            'grid_point_surface_strains',
            'grid_point_strains_volume_direct',
            'grid_point_strains_volume_principal',
            'grid_point_strain_discontinuities',
        ]
        utables = unique(table_types)
        if len(table_types) != len(utables):
            msg = 'Non-unique tables: '
            for i, table_type in enumerate(table_types):
                if table_type in table_types[i+1:]:
                    msg += table_type + ', '
            raise AssertionError(msg)
        return table_types

    def _get_table_types_testing(self) -> list[str]:
        """testing method...don't use"""
        table_types = self.get_table_types()

        #if USER != 'sdoyle':
        return table_types
        stress = self.op2_results.stress
        strain = self.op2_results.strain
        force = self.op2_results.force
        strain_energy = self.op2_results.strain_energy
        skipped_attributes = [
            'card_count', 'data_code', 'element_mapper', 'isubcase_name_map',
            'labels', 'subtitles', 'additional_matrices', 'matrices', 'matdicts',
            'subcase_key', 'end_options', 'expected_times', 'generalized_tables',
            'op2_reader', 'table_count', 'table_mapper', ] + \
            stress.get_table_types(include_class=False) + strain.get_table_types(include_class=False) + \
            force.get_table_types(include_class=False) + strain_energy.get_table_types(include_class=False)

        tables = object_attributes(self, 'public', filter_properties=True)
        tables_with_properties = object_attributes(
            self, 'public', filter_properties=False,
            #keys_to_skip=skipped_attributes,
        )
        properties = set(tables_with_properties) - set(tables)
        assert len(properties) == 0, properties

        strain_energy = self.op2_results.strain_energy
        stress = self.op2_results.stress
        strain = self.op2_results.strain
        force = self.op2_results.force
        #thermal_load = self.op2_results.thermal_load
        #sum_objs = self.get_sum_objects()

        missing_attrs = []
        sum_objs = [strain_energy, stress, strain, force, ] # thermal_load
        for obj in sum_objs:
            for prop_name in obj.get_table_types():
                sprop_name = prop_name.split('.')[1]
                if sprop_name not in properties:
                    missing_attrs.append(prop_name)
        if missing_attrs:
            msg = 'missing the following getattr/setattr methods:  \n' + '  \n'.join(missing_attrs)
            raise RuntimeError(msg)
            #self.log.warning(msg)
        #for prop_name in properties:
            #if 'strain_energy' in prop_name:
                #getattr(strain_energy, prop_name)
            #print(prop_name)

        tables = [table for table in tables
                  if isinstance(getattr(self, table), dict)
                  and table not in skipped_attributes]

        if self.make_geom:
            return table_types

        for table in tables:
            #value = getattr(self, table)
            assert table in table_types, table
        return table_types

    #def get_f06_stats(self):
        #return self.get_op2_stats()

    def get_op2_stats(self, short: bool=False) -> str:
        """
        Gets info about the contents of the different attributes of the
        OP2 class.

        Examples
        --------
        ***Detailed OP2 Stats***
        >>> self.get_op2_stats()

        displacements[1]
          isubcase = 1
          type=RealDisplacementArray nnodes=72
          data: [t1, t2, t3, r1, r2, r3] shape=[1, 72, 6] dtype=float32
          gridTypes
          sort1
          lsdvmns = [1]

        spc_forces[1]
          isubcase = 1
          type=RealSPCForcesArray nnodes=72
          data: [t1, t2, t3, r1, r2, r3] shape=[1, 72, 6] dtype=float32
          gridTypes
          sort1
          lsdvmns = [1]

        ctetra_stress[1]
          type=RealSolidStressArray nelements=186 nnodes=930
          nodes_per_element=5 (including centroid)
          eType, cid
          data: [1, nnodes, 10] where 10=[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, von_mises]
          data.shape = (1, 930, 10)
          element name: CTETRA
          sort1
          lsdvmns = [1]

        ***Appreviated OP2 Stats***
        >>> self.get_op2_stats(short=True)
        displacements[1]; RealDisplacementArray; [1, 72, 6]; [t1, t2, t3, r1, r2, r3]

        """
        return get_op2_stats(self, short=short)

class Op2F06Attributes(OP2_F06_Common):
    def __init__(self):
        OP2_F06_Common.__init__(self)
