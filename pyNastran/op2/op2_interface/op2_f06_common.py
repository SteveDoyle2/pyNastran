from __future__ import annotations
import getpass
from typing import TYPE_CHECKING
from numpy import unique, int32, int64

from pyNastran import is_release
from pyNastran.f06.f06_formatting import get_key0
from pyNastran.utils import object_attributes
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import deprecated
from pyNastran.bdf.case_control_deck import CaseControlDeck

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
        >>> self.eigenvectors = get_result('eigenvectors')

        **Example 2**
        >>> self.ato.displacements = get_result('ato.displacements')

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

    def del_result(self, result_name):
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
                storage_obj = delattr(self.op2_results, result_name)

    def deprecated(self, old_name: str, new_name: str, deprecated_version: str):
        """allows for simple OP2 vectorization"""
        return deprecated(old_name, new_name, deprecated_version, levels=[0, 1, 2])

    # ------------------------------------------------------------------
    # stress
    @property
    def celas1_stress(self):
        self.deprecated('model.celas1_stress', 'model.op2_results.stress.celas1_stress', '1.4')
        return self.op2_results.stress.celas1_stress
    @property
    def celas2_stress(self):
        self.deprecated('model.celas2_stress', 'model.op2_results.stress.celas2_stress', '1.4')
        return self.op2_results.stress.celas2_stress
    @property
    def celas3_stress(self):
        self.deprecated('model.celas3_stress', 'model.op2_results.stress.celas3_stress', '1.4')
        return self.op2_results.stress.celas3_stress
    @property
    def celas4_stress(self):
        self.deprecated('model.celas4_stress', 'model.op2_results.stress.celas4_stress', '1.4')
        return self.op2_results.stress.celas4_stress

    @property
    def crod_stress(self):
        self.deprecated('model.crod_stress', 'model.op2_results.stress.crod_stress', '1.4')
        return self.op2_results.stress.crod_stress
    @property
    def ctube_stress(self):
        self.deprecated('model.ctube_stress', 'model.op2_results.stress.ctube_stress', '1.4')
        return self.op2_results.stress.crod_stress
    @property
    def conrod_stress(self):
        self.deprecated('model.conrod_stress', 'model.op2_results.stress.conrod_stress', '1.4')
        return self.op2_results.stress.conrod_stress

    @property
    def cbar_stress(self):
        self.deprecated('model.cbar_stress', 'model.op2_results.stress.cbar_stress', '1.4')
        return self.op2_results.stress.cbar_stress
    @property
    def cbeam_stress(self):
        self.deprecated('model.cbeam_stress', 'model.op2_results.stress.cbeam_stress', '1.4')
        return self.op2_results.stress.cbeam_stress
    @property
    def cbend_stress(self):
        self.deprecated('model.cbend_stress', 'model.op2_results.stress.cbend_stress', '1.4')
        return self.op2_results.stress.cbend_stress

    @property
    def cquad4_stress(self):
        self.deprecated('model.cquad4_stress', 'model.op2_results.stress.cquad4_stress', '1.4')
        return self.op2_results.stress.cquad4_stress
    @property
    def cquad8_stress(self):
        self.deprecated('model.cquad8_stress', 'model.op2_results.stress.cquad8_stress', '1.4')
        return self.op2_results.stress.cquad8_stress
    @property
    def cquadr_stress(self):
        self.deprecated('model.cquadr_stress', 'model.op2_results.stress.cquadr_stress', '1.4')
        return self.op2_results.stress.cquadr_stress

    @property
    def ctria3_stress(self):
        self.deprecated('model.ctria3_stress', 'model.op2_results.stress.ctria3_stress', '1.4')
        return self.op2_results.stress.ctria3_stress
    @property
    def ctria6_stress(self):
        self.deprecated('model.ctria6_stress', 'model.op2_results.stress.ctria6_stress', '1.4')
        return self.op2_results.stress.ctria6_stress
    @property
    def ctriar_stress(self):
        self.deprecated('model.ctriar_stress', 'model.op2_results.stress.ctriar_stress', '1.4')
        return self.op2_results.stress.ctriar_stress
    @property
    def ctriax_stress(self):
        return self.op2_results.stress.ctriax_stress

    @property
    def ctria3_composite_stress(self):
        self.deprecated('model.ctria3_composite_stress', 'model.op2_results.stress.ctria3_composite_stress', '1.4')
        return self.op2_results.stress.ctria3_composite_stress
    @property
    def ctria6_composite_stress(self):
        self.deprecated('model.ctria6_composite_stress', 'model.op2_results.stress.ctria6_composite_stress', '1.4')
        return self.op2_results.stress.ctria6_composite_stress
    @property
    def ctriar_composite_stress(self):
        self.deprecated('model.ctriar_composite_stress', 'model.op2_results.stress.ctriar_composite_stress', '1.4')
        return self.op2_results.stress.ctriar_composite_stress

    @property
    def cquad4_composite_stress(self):
        self.deprecated('model.cquad4_composite_stress', 'model.op2_results.stress.cquad4_composite_stress', '1.4')
        return self.op2_results.stress.cquad4_composite_stress
    @property
    def cquad8_composite_stress(self):
        self.deprecated('model.cquad8_composite_stress', 'model.op2_results.stress.cquad8_composite_stress', '1.4')
        return self.op2_results.stress.cquad8_composite_stress
    @property
    def cquadr_composite_stress(self):
        self.deprecated('model.cquadr_composite_stress', 'model.op2_results.stress.cquadr_composite_stress', '1.4')
        return self.op2_results.stress.cquadr_composite_stress

    @property
    def ctetra_stress(self):
        self.deprecated('model.ctetra_stress', 'model.op2_results.stress.ctetra_stress', '1.4')
        return self.op2_results.stress.ctetra_stress
    @property
    def chexa_stress(self):
        self.deprecated('model.chexa_stress', 'model.op2_results.stress.chexa_stress', '1.4')
        return self.op2_results.stress.chexa_stress
    @property
    def cpenta_stress(self):
        self.deprecated('model.cpenta_stress', 'model.op2_results.stress.cpenta_stress', '1.4')
        return self.op2_results.stress.cpenta_stress
    @property
    def cpyram_stress(self):
        self.deprecated('model.cpyram_stress', 'model.op2_results.stress.cpyram_stress', '1.4')
        return self.op2_results.stress.cpyram_stress

    @property
    def chexa_composite_stress(self):
        self.deprecated('model.chexa_composite_stress', 'model.op2_results.stress.chexa_composite_stress', '1.4')
        return self.op2_results.stress.chexa_composite_stress
    @property
    def cpenta_composite_stress(self):
        self.deprecated('model.cpenta_composite_stress', 'model.op2_results.stress.cpenta_composite_stress', '1.4')
        return self.op2_results.stress.cpenta_composite_stress

    @property
    def cshear_stress(self):
        self.deprecated('model.cshear_stress', 'model.op2_results.stress.cshear_stress', '1.4')
        return self.op2_results.stress.cshear_stress
    @property
    def cweld_stress(self):
        self.deprecated('model.cweld_stress', 'model.op2_results.stress.cweld_stress', '1.4')
        return self.op2_results.stress.cweld_stress
    @property
    def cbush_stress(self):
        self.deprecated('model.cbush_stress', 'model.op2_results.stress.cbush_stress', '1.4')
        return self.op2_results.stress.cbush_stress
    @property
    def cbar_stress_10nodes(self):
        self.deprecated('model.cbar_stress_10nodes', 'model.op2_results.stress.cbar_stress_10nodes', '1.4')
        return self.op2_results.stress.cbar_stress_10nodes
    @property
    def cfast_stress(self):
        self.deprecated('model.cfast_stress', 'model.op2_results.stress.cfast_stress', '1.4')
        return self.op2_results.stress.cfast_stress

    @property
    def cplstn3_stress(self):
        self.deprecated('model.cplstn3_stress', 'model.op2_results.stress.cplstn3_stress', '1.4')
        return self.op2_results.stress.cplstn3_stress
    @property
    def cplstn4_stress(self):
        self.deprecated('model.cplstn4_stress', 'model.op2_results.stress.cplstn4_stress', '1.4')
        return self.op2_results.stress.cplstn4_stress
    @property
    def cplstn6_stress(self):
        self.deprecated('model.cplstn6_stress', 'model.op2_results.stress.cplstn6_stress', '1.4')
        return self.op2_results.stress.cplstn6_stress
    @property
    def cplstn8_stress(self):
        self.deprecated('model.cplstn8_stress', 'model.op2_results.stress.cplstn8_stress', '1.4')
        return self.op2_results.stress.cplstn8_stress

    @property
    def cplsts3_stress(self):
        self.deprecated('model.cplsts3_stress', 'model.op2_results.stress.cplsts3_stress', '1.4')
        return self.op2_results.stress.cplsts3_stress
    @property
    def cplsts4_stress(self):
        self.deprecated('model.cplsts4_stress', 'model.op2_results.stress.cplsts4_stress', '1.4')
        return self.op2_results.stress.cplsts4_stress
    @property
    def cplsts6_stress(self):
        self.deprecated('model.cplsts6_stress', 'model.op2_results.stress.cplsts6_stress', '1.4')
        return self.op2_results.stress.cplsts6_stress
    @property
    def cplsts8_stress(self):
        self.deprecated('model.cplsts8_stress', 'model.op2_results.stress.cplsts8_stress', '1.4')
        return self.op2_results.stress.cplsts8_stress

    @property
    def hyperelastic_cquad4_stress(self):
        self.deprecated('model.hyperelastic_cquad4_stress', 'model.op2_results.stress.hyperelastic_cquad4_stress', '1.4')
        return self.op2_results.stress.hyperelastic_cquad4_stress
    @property
    def cbush1d_stress_strain(self):
        self.deprecated('model.cbush1d_stress_strain', 'model.op2_results.stress.cbush1d_stress_strain', '1.4')
        self.op2_results.stress.cbush1d_stress_strain

    ##@property
    ##def ctriax6_stress(self):
        ##return self.op2_results.stress.ctriax6_stress
    ##@property
    ##def cquadx_stress(self):
        ##return self.op2_results.stress.cquadx_stress
    ##@property
    ##def cbeam3_stress(self):
        ##return self.op2_results.stress.cbeam3_stress

    # ------------------------------------------------------------------
    # strain
    @property
    def celas1_strain(self):
        self.deprecated('model.celas1_strain', 'model.op2_results.strain.celas1_strain', '1.4')
        return self.op2_results.strain.celas1_strain
    @property
    def celas2_strain(self):
        self.deprecated('model.celas2_strain', 'model.op2_results.strain.celas2_strain', '1.4')
        return self.op2_results.strain.celas2_strain
    @property
    def celas3_strain(self):
        self.deprecated('model.celas3_strain', 'model.op2_results.strain.celas3_strain', '1.4')
        return self.op2_results.strain.celas3_strain
    @property
    def celas4_strain(self):
        self.deprecated('model.celas4_strain', 'model.op2_results.strain.celas4_strain', '1.4')
        return self.op2_results.strain.celas4_strain

    @property
    def crod_strain(self):
        self.deprecated('model.crod_strain', 'model.op2_results.strain.crod_strain', '1.4')
        return self.op2_results.strain.crod_strain
    @property
    def ctube_strain(self):
        self.deprecated('model.ctube_strain', 'model.op2_results.strain.ctube_strain', '1.4')
        return self.op2_results.strain.crod_strain
    @property
    def conrod_strain(self):
        self.deprecated('model.conrod_strain', 'model.op2_results.strain.conrod_strain', '1.4')
        return self.op2_results.strain.conrod_strain

    @property
    def cbar_strain(self):
        self.deprecated('model.cbar_strain', 'model.op2_results.strain.cbar_strain', '1.4')
        return self.op2_results.strain.cbar_strain
    @property
    def cbeam_strain(self):
        self.deprecated('model.cbeam_strain', 'model.op2_results.strain.cbeam_strain', '1.4')
        return self.op2_results.strain.cbeam_strain
    @property
    def cbend_strain(self):
        self.deprecated('model.cbend_strain', 'model.op2_results.strain.cbend_strain', '1.4')
        return self.op2_results.strain.cbend_strain

    @property
    def cquad4_strain(self):
        self.deprecated('model.cquad4_strain', 'model.op2_results.strain.cquad4_strain', '1.4')
        return self.op2_results.strain.cquad4_strain
    @property
    def cquad8_strain(self):
        self.deprecated('model.cquad8_strain', 'model.op2_results.strain.cquad8_strain', '1.4')
        return self.op2_results.strain.cquad8_strain
    @property
    def cquadr_strain(self):
        self.deprecated('model.cquadr_strain', 'model.op2_results.strain.cquadr_strain', '1.4')
        return self.op2_results.strain.cquadr_strain

    @property
    def ctria3_strain(self):
        self.deprecated('model.ctria3_strain', 'model.op2_results.strain.ctria3_strain', '1.4')
        return self.op2_results.strain.ctria3_strain
    @property
    def ctria6_strain(self):
        self.deprecated('model.ctria6_strain', 'model.op2_results.strain.ctria6_strain', '1.4')
        return self.op2_results.strain.ctria6_strain
    @property
    def ctriar_strain(self):
        self.deprecated('model.ctriar_strain', 'model.op2_results.strain.ctriar_strain', '1.4')
        return self.op2_results.strain.ctriar_strain
    @property
    def ctriax_strain(self):
        self.deprecated('model.ctriax_strain', 'model.op2_results.strain.ctriax_strain', '1.4')
        return self.op2_results.strain.ctriax_strain
    @property
    def ctriax6_strain(self):
        self.deprecated('model.ctriax6_strain', 'model.op2_results.strain.ctriax6_strain', '1.4')
        return self.op2_results.strain.ctriax6_strain

    @property
    def ctetra_strain(self):
        self.deprecated('model.ctetra_strain', 'model.op2_results.strain.ctetra_strain', '1.4')
        return self.op2_results.strain.ctetra_strain
    @property
    def cpenta_strain(self):
        self.deprecated('model.cpenta_strain', 'model.op2_results.strain.cpenta_strain', '1.4')
        return self.op2_results.strain.cpenta_strain
    @property
    def chexa_strain(self):
        self.deprecated('model.chexa_strain', 'model.op2_results.strain.chexa_strain', '1.4')
        return self.op2_results.strain.chexa_strain
    @property
    def cpyram_strain(self):
        self.deprecated('model.cpyram_strain', 'model.op2_results.strain.cpyram_strain', '1.4')
        return self.op2_results.strain.cpyram_strain

    @property
    def cplstn3_strain(self):
        self.deprecated('model.cplstn3_strain', 'model.op2_results.strain.cplstn3_strain', '1.4')
        return self.op2_results.strain.cplstn3_strain
    @property
    def cplstn4_strain(self):
        self.deprecated('model.cplstn4_strain', 'model.op2_results.strain.cplstn4_strain', '1.4')
        return self.op2_results.strain.cplstn4_strain
    @property
    def cplstn6_strain(self):
        self.deprecated('model.cplstn6_strain', 'model.op2_results.strain.cplstn6_strain', '1.4')
        return self.op2_results.strain.cplstn6_strain
    @property
    def cplstn8_strain(self):
        self.deprecated('model.cplstn8_strain', 'model.op2_results.strain.cplstn8_strain', '1.4')
        return self.op2_results.strain.cplstn8_strain

    @property
    def cplsts3_strain(self):
        self.deprecated('model.cplsts3_strain', 'model.op2_results.strain.cplsts3_strain', '1.4')
        return self.op2_results.strain.cplsts3_strain
    @property
    def cplsts4_strain(self):
        self.deprecated('model.cplsts4_strain', 'model.op2_results.strain.cplsts4_strain', '1.4')
        return self.op2_results.strain.cplsts4_strain
    @property
    def cplsts6_strain(self):
        self.deprecated('model.cplsts6_strain', 'model.op2_results.strain.cplsts6_strain', '1.4')
        return self.op2_results.strain.cplsts6_strain
    @property
    def cplsts8_strain(self):
        self.deprecated('model.cplsts8_strain', 'model.op2_results.strain.cplsts8_strain', '1.4')
        return self.op2_results.strain.cplsts8_strain

    @property
    def cshear_strain(self):
        self.deprecated('model.cshear_strain', 'model.op2_results.stress.cshear_strain', '1.4')
        return self.op2_results.strain.cshear_strain
    @property
    def cweld_strain(self):
        self.deprecated('model.cweld_strain', 'model.op2_results.stress.cweld_strain', '1.4')
        return self.op2_results.strain.cweld_strain
    @property
    def cbush_strain(self):
        self.deprecated('model.cbush_strain', 'model.op2_results.stress.cbush_strain', '1.4')
        return self.op2_results.strain.cbush_strain
    @property
    def cfast_strain(self):
        self.deprecated('model.cfast_strain', 'model.op2_results.stress.cfast_strain', '1.4')
        return self.op2_results.strain.cfast_strain
    @property
    def cbar_strain_10nodes(self):
        self.deprecated('model.cbar_strain_10nodes', 'model.op2_results.strain.cbar_strain_10nodes', '1.4')
        return self.op2_results.strain.cbar_strain_10nodes

    @property
    def chexa_composite_strain(self):
        self.deprecated('model.chexa_composite_strain', 'model.op2_results.strain.chexa_composite_strain', '1.4')
        return self.op2_results.strain.chexa_composite_strain
    @property
    def cpenta_composite_strain(self):
        self.deprecated('model.cpenta_composite_strain', 'model.op2_results.strain.cpenta_composite_strain', '1.4')
        return self.op2_results.strain.cpenta_composite_strain
    @property
    def hyperelastic_cquad4_strain(self):
        self.deprecated('model.hyperelastic_cquad4_strain', 'model.op2_results.strain.hyperelastic_cquad4_strain', '1.4')
        return self.op2_results.strain.hyperelastic_cquad4_strain

    @chexa_composite_strain.setter
    def chexa_composite_strain(self, chexa_composite_strain):
        self.deprecated('model.chexa_composite_strain', 'model.op2_results.strain.chexa_composite_strain', '1.4')
        self.op2_results.strain.chexa_composite_strain = chexa_composite_strain
    @cpenta_composite_strain.setter
    def cpenta_composite_strain(self, cpenta_composite_strain):
        self.deprecated('model.cpenta_composite_strain', 'model.op2_results.strain.cpenta_composite_strain', '1.4')
        self.op2_results.strain.cpenta_composite_strain = cpenta_composite_strain
    @hyperelastic_cquad4_strain.setter
    def hyperelastic_cquad4_strain(self, hyperelastic_cquad4_strain):
        self.deprecated('model.hyperelastic_cquad4_strain', 'model.op2_results.strain.hyperelastic_cquad4_strain', '1.4')
        self.op2_results.strain.hyperelastic_cquad4_strain = hyperelastic_cquad4_strain

    # ------------------------------------------------------------------
    # Strain Energy - Getter
    @property
    def celas1_strain_energy(self):
        self.deprecated('model.celas1_strain_energy', 'model.op2_results.strain_energy.celas1_strain_energy', '1.4')
        return self.op2_results.strain_energy.celas1_strain_energy
    @property
    def celas2_strain_energy(self):
        self.deprecated('model.celas2_strain_energy', 'model.op2_results.strain_energy.celas2_strain_energy', '1.4')
        return self.op2_results.strain_energy.celas2_strain_energy
    @property
    def celas3_strain_energy(self):
        self.deprecated('model.celas3_strain_energy', 'model.op2_results.strain_energy.celas3_strain_energy', '1.4')
        return self.op2_results.strain_energy.celas3_strain_energy
    @property
    def celas4_strain_energy(self):
        self.deprecated('model.celas4_strain_energy', 'model.op2_results.strain_energy.celas4_strain_energy', '1.4')
        return self.op2_results.strain_energy.celas4_strain_energy

    @property
    def crod_strain_energy(self):
        self.deprecated('model.crod_strain_energy', 'model.op2_results.strain_energy.crod_strain_energy', '1.4')
        return self.op2_results.strain_energy.crod_strain_energy
    @property
    def ctube_strain_energy(self):
        self.deprecated('model.ctube_strain_energy', 'model.op2_results.strain_energy.ctube_strain_energy', '1.4')
        return self.op2_results.strain_energy.crod_strain_energy
    @property
    def conrod_strain_energy(self):
        self.deprecated('model.conrod_strain_energy', 'model.op2_results.strain_energy.conrod_strain_energy', '1.4')
        return self.op2_results.strain_energy.conrod_strain_energy

    @property
    def cquad4_strain_energy(self):
        self.deprecated('model.cquad4_strain_energy', 'model.op2_results.strain_energy.cquad4_strain_energy', '1.4')
        return self.op2_results.strain_energy.cquad4_strain_energy
    @property
    def cquad8_strain_energy(self):
        self.deprecated('model.cquad8_strain_energy', 'model.op2_results.strain_energy.ctetra_strain_energy', '1.4')
        return self.op2_results.strain_energy.cquad8_strain_energy
    @property
    def cquadr_strain_energy(self):
        self.deprecated('model.cquadr_strain_energy', 'model.op2_results.strain_energy.cquadr_strain_energy', '1.4')
        return self.op2_results.strain_energy.cquadr_strain_energy
    @property
    def cquadx_strain_energy(self):
        self.deprecated('model.cquadx_strain_energy', 'model.op2_results.strain_energy.cquadx_strain_energy', '1.4')
        return self.op2_results.strain_energy.cquadx_strain_energy

    @property
    def ctria3_strain_energy(self):
        self.deprecated('model.ctria3_strain_energy', 'model.op2_results.strain_energy.ctria3_strain_energy', '1.4')
        return self.op2_results.strain_energy.ctria3_strain_energy
    @property
    def ctria6_strain_energy(self):
        self.deprecated('model.ctria6_strain_energy', 'model.op2_results.strain_energy.ctria6_strain_energy', '1.4')
        return self.op2_results.strain_energy.ctria6_strain_energy
    @property
    def ctriar_strain_energy(self):
        self.deprecated('model.ctriar_strain_energy', 'model.op2_results.strain_energy.ctriar_strain_energy', '1.4')
        return self.op2_results.strain_energy.ctriar_strain_energy
    @property
    def ctriax_strain_energy(self):
        self.deprecated('model.ctriax_strain_energy', 'model.op2_results.strain_energy.ctriax_strain_energy', '1.4')
        return self.op2_results.strain_energy.ctriax_strain_energy
    @property
    def ctriax6_strain_energy(self):
        self.deprecated('model.ctriax6_strain_energy', 'model.op2_results.strain_energy.ctriax6_strain_energy', '1.4')
        return self.op2_results.strain_energy.ctriax6_strain_energy

    @property
    def ctetra_strain_energy(self):
        self.deprecated('model.ctetra_strain_energy', 'model.op2_results.strain_energy.ctetra_strain_energy', '1.4')
        return self.op2_results.strain_energy.ctetra_strain_energy
    @property
    def cpenta_strain_energy(self):
        self.deprecated('model.cpenta_strain_energy', 'model.op2_results.strain_energy.cpenta_strain_energy', '1.4')
        return self.op2_results.strain_energy.cpenta_strain_energy
    @property
    def chexa_strain_energy(self):
        self.deprecated('model.chexa_strain_energy', 'model.op2_results.strain_energy.chexa_strain_energy', '1.4')
        return self.op2_results.strain_energy.chexa_strain_energy
    @property
    def cpyram_strain_energy(self):
        self.deprecated('model.cpyram_strain_energy', 'model.op2_results.strain_energy.cpyram_strain_energy', '1.4')
        return self.op2_results.strain_energy.cpyram_strain_energy

    @property
    def cbar_strain_energy(self):
        self.deprecated('model.cbar_strain_energy', 'model.op2_results.strain_energy.cbar_strain_energy', '1.4')
        return self.op2_results.strain_energy.cbar_strain_energy
    @property
    def cbeam_strain_energy(self):
        self.deprecated('model.cbeam_strain_energy', 'model.op2_results.strain_energy.cbeam_strain_energy', '1.4')
        return self.op2_results.strain_energy.cbeam_strain_energy
    @property
    def cbeam3_strain_energy(self):
        self.deprecated('model.cbeam3_strain_energy', 'model.op2_results.strain_energy.cbeam3_strain_energy', '1.4')
        return self.op2_results.strain_energy.cbeam3_strain_energy
    @property
    def cbend_strain_energy(self):
        self.deprecated('model.cbend_strain_energy', 'model.op2_results.strain_energy.cbend_strain_energy', '1.4')
        return self.op2_results.strain_energy.cbend_strain_energy

    # ------------------------------------------------------
    # Strain Energy - Getter 2
    @property
    def cgap_strain_energy(self):
        self.deprecated('model.cgap_strain_energy', 'model.op2_results.strain_energy.cgap_strain_energy', '1.4')
        return self.op2_results.strain_energy.cgap_strain_energy
    @property
    def cdamp1_strain_energy(self):
        self.deprecated('model.cdamp1_strain_energy', 'model.op2_results.strain_energy.cdamp1_strain_energy', '1.4')
        return self.op2_results.strain_energy.cdamp1_strain_energy
    @property
    def cdamp2_strain_energy(self):
        self.deprecated('model.cdamp2_strain_energy', 'model.op2_results.strain_energy.cdamp2_strain_energy', '1.4')
        return self.op2_results.strain_energy.cdamp2_strain_energy
    @property
    def cdamp3_strain_energy(self):
        self.deprecated('model.cdamp3_strain_energy', 'model.op2_results.strain_energy.cdamp3_strain_energy', '1.4')
        return self.op2_results.strain_energy.cdamp3_strain_energy
    @property
    def cdamp4_strain_energy(self):
        self.deprecated('model.cdamp4_strain_energy', 'model.op2_results.strain_energy.cdamp4_strain_energy', '1.4')
        return self.op2_results.strain_energy.cdamp4_strain_energy
    @property
    def cbush_strain_energy(self):
        self.deprecated('model.cbush_strain_energy', 'model.op2_results.strain_energy.cbush_strain_energy', '1.4')
        return self.op2_results.strain_energy.cbush_strain_energy
    @property
    def dmig_strain_energy(self):
        self.deprecated('model.dmig_strain_energy', 'model.op2_results.strain_energy.dmig_strain_energy', '1.4')
        return self.op2_results.strain_energy.dmig_strain_energy
    @property
    def genel_strain_energy(self):
        self.deprecated('model.genel_strain_energy', 'model.op2_results.strain_energy.genel_strain_energy', '1.4')
        return self.op2_results.strain_energy.genel_strain_energy
    @property
    def cshear_strain_energy(self):
        self.deprecated('model.cshear_strain_energy', 'model.op2_results.strain_energy.cshear_strain_energy', '1.4')
        return self.op2_results.strain_energy.cshear_strain_energy
    @property
    def conm2_strain_energy(self):
        self.deprecated('model.conm2_strain_energy', 'model.op2_results.strain_energy.conm2_strain_energy', '1.4')
        return self.op2_results.strain_energy.conm2_strain_energy
    @property
    def cdum8_strain_energy(self):
        self.deprecated('model.cdum8_strain_energy', 'model.op2_results.strain_energy.cdum8_strain_energy', '1.4')
        return self.op2_results.strain_energy.cdum8_strain_energy
    @property
    def rbe1_strain_energy(self):
        self.deprecated('model.rbe1_strain_energy', 'model.op2_results.strain_energy.rbe1_strain_energy', '1.4')
        return self.op2_results.strain_energy.rbe1_strain_energy
    @property
    def rbe3_strain_energy(self):
        self.deprecated('model.rbe3_strain_energy', 'model.op2_results.strain_energy.rbe3_strain_energy', '1.4')
        return self.op2_results.strain_energy.rbe3_strain_energy
    @property
    def cweld_strain_energy(self):
        self.deprecated('model.cweld_strain_energy', 'model.op2_results.strain_energy.cweld_strain_energy', '1.4')
        return self.op2_results.strain_energy.cweld_strain_energy
    @property
    def cfast_strain_energy(self):
        self.deprecated('model.cfast_strain_energy', 'model.op2_results.strain_energy.cfast_strain_energy', '1.4')
        return self.op2_results.strain_energy.cfast_strain_energy
    @property
    def cseam_strain_energy(self):
        self.deprecated('model.cseam_strain_energy', 'model.op2_results.strain_energy.cseam_strain_energy', '1.4')
        return self.op2_results.strain_energy.cseam_strain_energy
    # ------------------------------------------------------------------
    # Strain Energy - Setter

    @celas1_strain_energy.setter
    def celas1_strain_energy(self, celas1_strain_energy):
        self.deprecated('model.celas1_strain_energy', 'model.op2_results.strain_energy.celas1_strain_energy', '1.4')
        self.op2_results.strain_energy.celas1_strain_energy = celas1_strain_energy
    @celas2_strain_energy.setter
    def celas2_strain_energy(self, celas2_strain_energy):
        self.deprecated('model.celas2_strain_energy', 'model.op2_results.strain_energy.celas2_strain_energy', '1.4')
        self.op2_results.strain_energy.celas2_strain_energy = celas2_strain_energy
    @celas3_strain_energy.setter
    def celas3_strain_energy(self, celas3_strain_energy):
        self.deprecated('model.celas3_strain_energy', 'model.op2_results.strain_energy.celas3_strain_energy', '1.4')
        self.op2_results.strain_energy.celas3_strain_energy = celas3_strain_energy
    @celas4_strain_energy.setter
    def celas4_strain_energy(self, celas4_strain_energy):
        self.deprecated('model.celas4_strain_energy', 'model.op2_results.strain_energy.celas4_strain_energy', '1.4')
        self.op2_results.strain_energy.celas4_strain_energy = celas4_strain_energy


    @cdamp1_strain_energy.setter
    def cdamp1_strain_energy(self, cdamp1_strain_energy):
        self.deprecated('model.cdamp1_strain_energy', 'model.op2_results.strain_energy.cdamp1_strain_energy', '1.4')
        self.op2_results.strain_energy.cdamp1_strain_energy = cdamp1_strain_energy
    @cdamp2_strain_energy.setter
    def cdamp2_strain_energy(self, cdamp2_strain_energy):
        self.deprecated('model.cdamp2_strain_energy', 'model.op2_results.strain_energy.cdamp2_strain_energy', '1.4')
        self.op2_results.strain_energy.cdamp2_strain_energy = cdamp2_strain_energy
    @cdamp3_strain_energy.setter
    def cdamp3_strain_energy(self, cdamp3_strain_energy):
        self.deprecated('model.cdamp3_strain_energy', 'model.op2_results.strain_energy.cdamp3_strain_energy', '1.4')
        self.op2_results.strain_energy.cdamp3_strain_energy = cdamp3_strain_energy
    @cdamp4_strain_energy.setter
    def cdamp4_strain_energy(self, cdamp4_strain_energy):
        self.deprecated('model.cdamp4_strain_energy', 'model.op2_results.strain_energy.cdamp4_strain_energy', '1.4')
        self.op2_results.strain_energy.cdamp4_strain_energy = cdamp4_strain_energy

    @cgap_strain_energy.setter
    def cgap_strain_energy(self, cgap_strain_energy):
        self.deprecated('model.cgap_strain_energy', 'model.op2_results.strain_energy.cgap_strain_energy', '1.4')
        self.op2_results.strain_energy.cgap_strain_energy = cgap_strain_energy
    @cbush_strain_energy.setter
    def cbush_strain_energy(self, cbush_strain_energy):
        self.deprecated('model.cbush_strain_energy', 'model.op2_results.strain_energy.cbush_strain_energy', '1.4')
        self.op2_results.strain_energy.cbush_strain_energy = cbush_strain_energy

    @crod_strain_energy.setter
    def crod_strain_energy(self, crod_strain_energy):
        self.deprecated('model.crod_strain_energy', 'model.op2_results.strain_energy.crod_strain_energy', '1.4')
        self.op2_results.strain_energy.crod_strain_energy = crod_strain_energy
    @ctube_strain_energy.setter
    def ctube_strain_energy(self, ctube_strain_energy):
        self.deprecated('model.ctube_strain_energy', 'model.op2_results.strain_energy.ctube_strain_energy', '1.4')
        self.op2_results.strain_energy.crod_strain_energy = ctube_strain_energy
    @conrod_strain_energy.setter
    def conrod_strain_energy(self, conrod_strain_energy):
        self.deprecated('model.conrod_strain_energy', 'model.op2_results.strain_energy.conrod_strain_energy', '1.4')
        self.op2_results.strain_energy.conrod_strain_energy = conrod_strain_energy

    @cquad4_strain_energy.setter
    def cquad4_strain_energy(self, cquad4_strain_energy):
        self.deprecated('model.cquad4_strain_energy', 'model.op2_results.strain_energy.cquad4_strain_energy', '1.4')
        self.op2_results.strain_energy.cquad4_strain_energy = cquad4_strain_energy
    @cquad8_strain_energy.setter
    def cquad8_strain_energy(self, cquad8_strain_energy):
        self.deprecated('model.cquad8_strain_energy', 'model.op2_results.strain_energy.cquad8_strain_energy', '1.4')
        self.op2_results.strain_energy.cquad8_strain_energy = cquad8_strain_energy
    @cquadr_strain_energy.setter
    def cquadr_strain_energy(self, cquadr_strain_energy):
        self.deprecated('model.cquadr_strain_energy', 'model.op2_results.strain_energy.cquadr_strain_energy', '1.4')
        self.op2_results.strain_energy.cquadr_strain_energy = cquadr_strain_energy
    @cquadx_strain_energy.setter
    def cquadx_strain_energy(self, cquadx_strain_energy):
        self.deprecated('model.cquadx_strain_energy', 'model.op2_results.strain_energy.cquadx_strain_energy', '1.4')
        self.op2_results.strain_energy.cquadx_strain_energy = cquadx_strain_energy

    @ctetra_strain_energy.setter
    def ctria3_strain_energy(self, ctria3_strain_energy):
        self.deprecated('model.ctria3_strain_energy', 'model.op2_results.strain_energy.ctria3_strain_energy', '1.4')
        self.op2_results.strain_energy.ctria3_strain_energy = ctria3_strain_energy
    @ctria6_strain_energy.setter
    def ctria6_strain_energy(self, ctria6_strain_energy):
        self.deprecated('model.ctria6_strain_energy', 'model.op2_results.strain_energy.ctria6_strain_energy', '1.4')
        self.op2_results.strain_energy.ctria6_strain_energy = ctria6_strain_energy
    @ctriar_strain_energy.setter
    def ctriar_strain_energy(self, ctriar_strain_energy):
        self.deprecated('model.ctriar_strain_energy', 'model.op2_results.strain_energy.ctriar_strain_energy', '1.4')
        self.op2_results.strain_energy.ctriar_strain_energy = ctriar_strain_energy
    @ctetra_strain_energy.setter
    def ctriax_strain_energy(self, ctriax_strain_energy):
        self.deprecated('model.celas2_strain_energy', 'model.op2_results.strain_energy.celas2_strain_energy', '1.4')
        self.op2_results.strain_energy.ctriax_strain_energy = ctriax_strain_energy
    @ctetra_strain_energy.setter
    def ctriax6_strain_energy(self, ctriax6_strain_energy):
        self.deprecated('model.celas2_strain_energy', 'model.op2_results.strain_energy.celas2_strain_energy', '1.4')
        self.op2_results.strain_energy.ctriax6_strain_energy = ctriax6_strain_energy

    @ctetra_strain_energy.setter
    def ctetra_strain_energy(self, ctetra_strain_energy):
        self.deprecated('model.ctetra_strain_energy', 'model.op2_results.strain_energy.ctetra_strain_energy', '1.4')
        self.op2_results.strain_energy.ctetra_strain_energy = ctetra_strain_energy
    @chexa_strain_energy.setter
    def chexa_strain_energy(self, chexa_strain_energy):
        self.deprecated('model.chexa_strain_energy', 'model.op2_results.strain_energy.chexa_strain_energy', '1.4')
        self.op2_results.strain_energy.chexa_strain_energy = chexa_strain_energy
    @cpenta_strain_energy.setter
    def cpenta_strain_energy(self, cpenta_strain_energy):
        self.deprecated('model.cpenta_strain_energy', 'model.op2_results.strain_energy.cpenta_strain_energy', '1.4')
        self.op2_results.strain_energy.cpenta_strain_energy = cpenta_strain_energy
    @cpyram_strain_energy.setter
    def cpyram_strain_energy(self, cpyram_strain_energy):
        self.deprecated('model.cpyram_strain_energy', 'model.op2_results.strain_energy.cpyram_strain_energy', '1.4')
        self.op2_results.strain_energy.cpyram_strain_energy = cpyram_strain_energy
    @cfast_strain_energy.setter
    def cfast_strain_energy(self, cfast_strain_energy):
        self.deprecated('model.cfast_strain_energy', 'model.op2_results.strain_energy.cfast_strain_energy', '1.4')
        self.op2_results.strain_energy.cfast_strain_energy = cfast_strain_energy
    @cseam_strain_energy.setter
    def cseam_strain_energy(self, cseam_strain_energy):
        self.deprecated('model.cseam_strain_energy', 'model.op2_results.strain_energy.cseam_strain_energy', '1.4')
        self.op2_results.strain_energy.cseam_strain_energy = cseam_strain_energy
    @rbe3_strain_energy.setter
    def rbe3_strain_energy(self, rbe3_strain_energy):
        self.deprecated('model.rbe3_strain_energy', 'model.op2_results.strain_energy.rbe3_strain_energy', '1.4')
        self.op2_results.strain_energy.rbe3_strain_energy = rbe3_strain_energy
    # ------------------------------------------------------------------
    # Stress - Getter
    @celas1_stress.setter
    def celas1_stress(self, celas1_stress):
        self.deprecated('model.celas1_stress', 'model.op2_results.stress.celas1_stress', '1.4')
        self.op2_results.stress.celas1_stress = celas1_stress
    @celas2_stress.setter
    def celas2_stress(self, celas2_stress):
        self.deprecated('model.celas2_stress', 'model.op2_results.stress.celas2_stress', '1.4')
        self.op2_results.stress.celas2_stress = celas2_stress
    @celas3_stress.setter
    def celas3_stress(self, celas3_stress):
        self.deprecated('model.celas3_stress', 'model.op2_results.stress.celas3_stress', '1.4')
        self.op2_results.stress.celas3_stress = celas3_stress
    @celas4_stress.setter
    def celas4_stress(self, celas4_stress):
        self.deprecated('model.celas4_stress', 'model.op2_results.stress.celas4_stress', '1.4')
        self.op2_results.stress.celas4_stress = celas4_stress

    @cbush_stress.setter
    def cbush_stress(self, cbush_stress):
        self.deprecated('model.cbush_stress', 'model.op2_results.stress.cbush_stress', '1.4')
        self.op2_results.stress.cbush_stress = cbush_stress

    @crod_stress.setter
    def crod_stress(self, crod_stress):
        self.deprecated('model.crod_stress', 'model.op2_results.stress.crod_stress', '1.4')
        self.op2_results.stress.crod_stress = crod_stress
    @ctube_stress.setter
    def ctube_stress(self, ctube_stress):
        self.deprecated('model.ctube_stress', 'model.op2_results.stress.ctube_stress', '1.4')
        self.op2_results.stress.crod_stress = ctube_stress
    @conrod_stress.setter
    def conrod_stress(self, conrod_stress):
        self.deprecated('model.conrod_stress', 'model.op2_results.stress.conrod_stress', '1.4')
        self.op2_results.stress.conrod_stress = conrod_stress

    @cquad4_stress.setter
    def cquad4_stress(self, cquad4_stress):
        self.deprecated('model.cquad4_stress', 'model.op2_results.stress.cquad4_stress', '1.4')
        self.op2_results.stress.cquad4_stress = cquad4_stress
    @cquad8_stress.setter
    def cquad8_stress(self, cquad8_stress):
        self.deprecated('model.cquad8_stress', 'model.op2_results.stress.cquad8_stress', '1.4')
        self.op2_results.stress.cquad8_stress = cquad8_stress
    @cquadr_stress.setter
    def cquadr_stress(self, cquadr_stress):
        self.deprecated('model.cquadr_stress', 'model.op2_results.stress.cquadr_stress', '1.4')
        self.op2_results.stress.cquadr_stress = cquadr_stress

    @ctria3_stress.setter
    def ctria3_stress(self, ctria3_stress):
        self.deprecated('model.ctria3_stress', 'model.op2_results.stress.ctria3_stress', '1.4')
        self.op2_results.stress.ctria3_stress = ctria3_stress
    @ctria6_stress.setter
    def ctria6_stress(self, ctria6_stress):
        self.deprecated('model.ctria6_stress', 'model.op2_results.stress.ctria6_stress', '1.4')
        self.op2_results.stress.ctria6_stress = ctria6_stress
    @ctriar_stress.setter
    def ctriar_stress(self, ctriar_stress):
        self.deprecated('model.ctriar_stress', 'model.op2_results.stress.ctriar_stress', '1.4')
        self.op2_results.stress.ctriar_stress = ctriar_stress
    @ctriax_stress.setter
    def ctriax_stress(self, ctriax_stress):
        self.deprecated('model.ctriax_stress', 'model.op2_results.stress.ctriax_stress', '1.4')
        self.op2_results.stress.ctriax_stress = ctriax_stress

    @ctetra_stress.setter
    def ctetra_stress(self, ctetra_stress):
        self.deprecated('model.ctetra_stress', 'model.op2_results.stress.ctetra_stress', '1.4')
        self.op2_results.stress.ctetra_stress = ctetra_stress
    @chexa_stress.setter
    def chexa_stress(self, chexa_stress):
        self.deprecated('model.chexa_stress', 'model.op2_results.stress.chexa_stress', '1.4')
        self.op2_results.stress.chexa_stress = chexa_stress
    @cpenta_stress.setter
    def cpenta_stress(self, cpenta_stress):
        self.deprecated('model.cpenta_stress', 'model.op2_results.stress.cpenta_stress', '1.4')
        self.op2_results.stress.cpenta_stress = cpenta_stress
    @cpyram_stress.setter
    def cpyram_stress(self, cpyram_stress):
        self.deprecated('model.cpyram_stress', 'model.op2_results.stress.cpyram_stress', '1.4')
        self.op2_results.stress.cpyram_stress = cpyram_stress

    @hyperelastic_cquad4_stress.setter
    def hyperelastic_cquad4_stress(self, hyperelastic_cquad4_stress):
        self.deprecated('model.hyperelastic_cquad4_stress', 'model.op2_results.stress.hyperelastic_cquad4_stress', '1.4')
        self.op2_results.stress.hyperelastic_cquad4_stress = hyperelastic_cquad4_stress
    @cbush1d_stress_strain.setter
    def cbush1d_stress_strain(self, cbush1d_stress_strain):
        self.deprecated('model.cbush1d_stress_strain', 'model.op2_results.stress.cbush1d_stress_strain', '1.4')
        self.op2_results.stress.cbush1d_stress_strain = cbush1d_stress_strain

    #-------------------------------------------------------------------
    # Strain - Getter
    @celas1_strain.setter
    def celas1_strain(self, celas1_strain):
        self.deprecated('model.celas1_strain', 'model.op2_results.strain.celas1_strain', '1.4')
        self.op2_results.strain.celas1_strain = celas1_strain
    @celas2_strain.setter
    def celas2_strain(self, celas2_strain):
        self.deprecated('model.celas2_strain', 'model.op2_results.strain.celas2_strain', '1.4')
        self.op2_results.strain.celas2_strain = celas2_strain
    @celas3_strain.setter
    def celas3_strain(self, celas3_strain):
        self.deprecated('model.celas3_strain', 'model.op2_results.strain.celas3_strain', '1.4')
        self.op2_results.strain.celas3_strain = celas3_strain
    @celas4_strain.setter
    def celas4_strain(self, celas4_strain):
        self.deprecated('model.celas4_strain', 'model.op2_results.strain.celas4_strain', '1.4')
        self.op2_results.strain.celas4_strain = celas4_strain

    @cbush_strain.setter
    def cbush_strain(self, cbush_strain):
        self.deprecated('model.cbush_strain', 'model.op2_results.strain.cbush_strain', '1.4')
        self.op2_results.strain.cbush_strain = cbush_strain

    @crod_strain.setter
    def crod_strain(self, crod_strain):
        self.deprecated('model.crod_strain', 'model.op2_results.strain.crod_strain', '1.4')
        self.op2_results.strain.crod_strain = crod_strain
    @ctube_strain.setter
    def ctube_strain(self, ctube_strain):
        self.deprecated('model.ctube_strain', 'model.op2_results.strain.ctube_strain', '1.4')
        self.op2_results.strain.crod_strain = ctube_strain
    @conrod_strain.setter
    def conrod_strain(self, conrod_strain):
        self.deprecated('model.conrod_strain', 'model.op2_results.strain.conrod_strain', '1.4')
        self.op2_results.strain.conrod_strain = conrod_strain

    @cquad4_strain.setter
    def cquad4_strain(self, cquad4_strain):
        self.deprecated('model.cquad4_strain', 'model.op2_results.strain.cquad4_strain', '1.4')
        self.op2_results.strain.cquad4_strain = cquad4_strain
    @cquad8_strain.setter
    def cquad8_strain(self, cquad8_strain):
        self.deprecated('model.cquad8_strain', 'model.op2_results.strain.cquad8_strain', '1.4')
        self.op2_results.strain.cquad8_strain = cquad8_strain
    @cquadr_strain.setter
    def cquadr_strain(self, cquadr_strain):
        self.deprecated('model.cquadr_strain', 'model.op2_results.strain.cquadr_strain', '1.4')
        self.op2_results.strain.cquadr_strain = cquadr_strain

    @ctetra_strain.setter
    def ctria3_strain(self, ctria3_strain):
        self.deprecated('model.ctria3_strain', 'model.op2_results.strain.ctria3_strain', '1.4')
        self.op2_results.strain.ctria3_strain = ctria3_strain
    @ctria6_strain.setter
    def ctria6_strain(self, ctria6_strain):
        self.deprecated('model.ctria6_strain', 'model.op2_results.strain.ctria6_strain', '1.4')
        self.op2_results.strain.ctria6_strain = ctria6_strain
    @ctriar_strain.setter
    def ctriar_strain(self, ctriar_strain):
        self.deprecated('model.ctriar_strain', 'model.op2_results.strain.ctriar_strain', '1.4')
        self.op2_results.strain.ctriar_strain = ctriar_strain
    @ctriax_strain.setter
    def ctriax_strain(self, ctriax_strain):
        self.deprecated('model.ctriax_strain', 'model.op2_results.strain.ctriax_strain', '1.4')
        self.op2_results.strain.ctriax_strain = ctriax_strain
    @ctriax6_strain.setter
    def ctriax6_strain(self, ctriax6_strain):
        self.deprecated('model.ctriax6_strain', 'model.op2_results.strain.ctriax6_strain', '1.4')
        self.op2_results.strain.ctriax6_strain = ctriax6_strain

    @ctetra_strain.setter
    def ctetra_strain(self, ctetra_strain):
        self.deprecated('model.ctetra_strain', 'model.op2_results.strain.ctetra_strain', '1.4')
        self.op2_results.strain.ctetra_strain = ctetra_strain
    @chexa_strain.setter
    def chexa_strain(self, chexa_strain):
        self.deprecated('model.chexa_strain', 'model.op2_results.strain.chexa_strain', '1.4')
        self.op2_results.strain.chexa_strain = chexa_strain
    @cpenta_strain.setter
    def cpenta_strain(self, cpenta_strain):
        self.deprecated('model.cpenta_strain', 'model.op2_results.strain.cpenta_strain', '1.4')
        self.op2_results.strain.cpenta_strain = cpenta_strain
    @cpyram_strain.setter
    def cpyram_strain(self, cpyram_strain):
        self.deprecated('model.cpyram_strain', 'model.op2_results.strain.cpyram_strain', '1.4')
        self.op2_results.strain.cpyram_strain = cpyram_strain
    @cfast_strain.setter
    def cfast_strain(self, cfast_strain):
        self.deprecated('model.cfast_strain', 'model.op2_results.strain.cfast_strain', '1.4')
        self.op2_results.strain.cfast_strain = cfast_strain

    @property
    def ctria3_composite_strain(self):
        self.deprecated('model.ctria3_composite_strain', 'model.op2_results.strain.ctria3_composite_strain', '1.4')
        return self.op2_results.strain.ctria3_composite_strain
    @property
    def ctria6_composite_strain(self):
        self.deprecated('model.ctria6_composite_strain', 'model.op2_results.strain.ctria6_composite_strain', '1.4')
        return self.op2_results.strain.ctria6_composite_strain
    @property
    def ctriar_composite_strain(self):
        self.deprecated('model.ctriar_composite_strain', 'model.op2_results.strain.ctriar_composite_strain', '1.4')
        return self.op2_results.strain.ctriar_composite_strain

    @property
    def cquad4_composite_strain(self):
        self.deprecated('model.cquad4_composite_strain', 'model.op2_results.strain.cquad4_composite_strain', '1.4')
        return self.op2_results.strain.ctria3_composite_strain
    @property
    def cquad8_composite_strain(self):
        self.deprecated('model.cquad8_composite_strain', 'model.op2_results.strain.cquad8_composite_strain', '1.4')
        return self.op2_results.strain.cquad8_composite_strain
    @property
    def cquadr_composite_strain(self):
        self.deprecated('model.cquadr_composite_strain', 'model.op2_results.strain.cquadr_composite_strain', '1.4')
        return self.op2_results.strain.cquadr_composite_strain

    @ctria3_composite_strain.setter
    def ctria3_composite_strain(self, ctria3_composite_strain):
        self.deprecated('model.ctria3_composite_strain', 'model.op2_results.strain.ctria3_composite_strain', '1.4')
        self.op2_results.strain.ctria3_composite_strain = ctria3_composite_strain
    @ctria6_composite_strain.setter
    def ctria6_composite_strain(self, ctria6_composite_strain):
        self.deprecated('model.ctria6_composite_strain', 'model.op2_results.strain.ctria6_composite_strain', '1.4')
        self.op2_results.strain.ctria3_composite_strain = ctria6_composite_strain
    @ctriar_composite_strain.setter
    def ctriar_composite_strain(self, ctriar_composite_strain):
        self.deprecated('model.ctriar_composite_strain', 'model.op2_results.strain.ctriar_composite_strain', '1.4')
        self.op2_results.strain.ctriar_composite_strain = ctriar_composite_strain

    @cquad4_composite_strain.setter
    def cquad4_composite_strain(self, cquad4_composite_strain):
        self.deprecated('model.cquad4_composite_strain', 'model.op2_results.strain.cquad4_composite_strain', '1.4')
        self.op2_results.strain.cquad4_composite_strain = cquad4_composite_strain
    @cquad8_composite_strain.setter
    def cquad8_composite_strain(self, cquad8_composite_strain):
        self.deprecated('model.cquad8_composite_strain', 'model.op2_results.strain.cquad8_composite_strain', '1.4')
        self.op2_results.strain.cquad8_composite_strain = cquad8_composite_strain
    @ctriar_composite_strain.setter
    def cquadr_composite_strain(self, cquadr_composite_strain):
        self.deprecated('model.cquadr_composite_strain', 'model.op2_results.strain.cquadr_composite_strain', '1.4')
        self.op2_results.strain.cquadr_composite_strain = cquadr_composite_strain

    @cshear_strain.setter
    def cshear_strain(self, cshear_strain):
        self.deprecated('model.cshear_strain', 'model.op2_results.stress.cshear_strain', '1.4')
        self.op2_results.strain.cshear_strain = cshear_strain
    @cweld_strain.setter
    def cweld_strain(self, cweld_strain):
        self.deprecated('model.cweld_strain', 'model.op2_results.stress.cweld_strain', '1.4')
        self.op2_results.strain.cweld_strain = cweld_strain
    @cbar_strain_10nodes.setter
    def cbar_strain_10nodes(self, cbar_strain_10nodes):
        self.deprecated('model.cbar_strain_10nodes', 'model.op2_results.strain.cbar_strain_10nodes', '1.4')
        self.op2_results.strain.cbar_strain_10nodes = cbar_strain_10nodes
    @cbeam_strain.setter
    def cbeam_strain(self, cbeam_strain):
        self.deprecated('model.cbeam_strain', 'model.op2_results.strain.cbeam_strain', '1.4')
        self.op2_results.strain.cbeam_strain = cbeam_strain

    # ------------------------------------------------------------------
    # Force - Getter
    @property
    def celas1_force(self):
        self.deprecated('model.celas1_force', 'model.op2_results.force.celas1_force', '1.4')
        return self.op2_results.force.celas1_force
    @property
    def celas2_force(self):
        self.deprecated('model.celas2_force', 'model.op2_results.force.celas2_force', '1.4')
        return self.op2_results.force.celas2_force
    @property
    def celas3_force(self):
        self.deprecated('model.celas3_force', 'model.op2_results.force.celas3_force', '1.4')
        return self.op2_results.force.celas3_force
    @property
    def celas4_force(self):
        self.deprecated('model.celas4_force', 'model.op2_results.force.celas4_force', '1.4')
        return self.op2_results.force.celas4_force

    @property
    def cdamp1_force(self):
        self.deprecated('model.cdamp1_force', 'model.op2_results.force.cdamp1_force', '1.4')
        return self.op2_results.force.cdamp1_force
    @property
    def cdamp2_force(self):
        self.deprecated('model.cdamp2_force', 'model.op2_results.force.cdamp2_force', '1.4')
        return self.op2_results.force.cdamp2_force
    @property
    def cdamp3_force(self):
        self.deprecated('model.cdamp3_force', 'model.op2_results.force.cdamp3_force', '1.4')
        return self.op2_results.force.cdamp3_force
    @property
    def cdamp4_force(self):
        self.deprecated('model.cdamp4_force', 'model.op2_results.force.cdamp4_force', '1.4')
        return self.op2_results.force.cdamp4_force

    # ------------------------------------------------------------------
    # Force - Setter
    @celas1_force.setter
    def celas1_force(self, celas1_force):
        self.deprecated('model.cdamp1_force', 'model.op2_results.force.cdamp1_force', '1.4')
        self.op2_results.force.celas1_force = celas1_force
    @celas2_force.setter
    def celas2_force(self, celas2_force):
        self.deprecated('model.celas2_force', 'model.op2_results.force.celas2_force', '1.4')
        self.op2_results.force.celas2_force = celas2_force
    @celas3_force.setter
    def celas3_force(self, celas3_force):
        self.deprecated('model.celas3_force', 'model.op2_results.force.celas3_force', '1.4')
        self.op2_results.force.celas3_force = celas3_force
    @celas4_force.setter
    def celas4_force(self, celas4_force):
        self.deprecated('model.celas4_force', 'model.op2_results.force.celas4_force', '1.4')
        self.op2_results.force.celas4_force = celas4_force

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

    @cdamp1_force.setter
    def cdamp1_force(self, cdamp1_force):
        self.deprecated('model.cdamp1_force', 'model.op2_results.force.cdamp1_force', '1.4')
        self.op2_results.force.cdamp1_force = cdamp1_force
    @cdamp2_force.setter
    def cdamp2_force(self, cdamp2_force):
        self.deprecated('model.cdamp2_force', 'model.op2_results.force.cdamp2_force', '1.4')
        self.op2_results.force.cdamp2_force = cdamp2_force
    @cdamp3_force.setter
    def cdamp3_force(self, cdamp3_force):
        self.deprecated('model.cdamp3_force', 'model.op2_results.force.cdamp3_force', '1.4')
        self.op2_results.force.cdamp3_force = cdamp3_force
    @cdamp4_force.setter
    def cdamp4_force(self, cdamp4_force):
        self.deprecated('model.cdamp4_force', 'model.op2_results.force.cdamp4_force', '1.4')
        self.op2_results.force.cdamp4_force = cdamp4_force

    @property
    def cgap_force(self):
        self.deprecated('model.cgap_force', 'model.op2_results.force.cgap_force', '1.4')
        return self.op2_results.force.cgap_force
    @property
    def cbush_force(self):
        self.deprecated('model.cbush_force', 'model.op2_results.force.cbush_force', '1.4')
        return self.op2_results.force.cbush_force
    @property
    def cweld_force(self):
        self.deprecated('model.cweld_force', 'model.op2_results.force.cweld_force', '1.4')
        return self.op2_results.force.cweld_force
    @property
    def cfast_force(self):
        self.deprecated('model.cfast_force', 'model.op2_results.force.cfast_force', '1.4')
        return self.op2_results.force.cfast_force
    @property
    def cbear_force(self):
        self.deprecated('model.cbear_force', 'model.op2_results.force.cbear_force', '1.4')
        return self.op2_results.force.cbear_force

    @property
    def cshear_force(self):
        self.deprecated('model.cshear_force', 'model.op2_results.force.cshear_force', '1.4')
        return self.op2_results.force.cshear_force
    @property
    def cbend_force(self):
        self.deprecated('model.cbend_force', 'model.op2_results.force.cbend_force', '1.4')
        return self.op2_results.force.cbend_force
    @property
    def cconeax_force(self):
        self.deprecated('model.cconeax_force', 'model.op2_results.force.cconeax_force', '1.4')
        return self.op2_results.force.cconeax_force

    @property
    def ctetra_pressure_force(self):
        self.deprecated('model.ctetra_pressure_force', 'model.op2_results.force.ctetra_pressure_force', '1.4')
        return self.op2_results.force.ctetra_pressure_force
    @property
    def cpenta_pressure_force(self):
        self.deprecated('model.cpenta_pressure_force', 'model.op2_results.force.cpenta_pressure_force', '1.4')
        return self.op2_results.force.cpenta_pressure_force
    @property
    def chexa_pressure_force(self):
        self.deprecated('model.chexa_pressure_force', 'model.op2_results.force.chexa_pressure_force', '1.4')
        return self.op2_results.force.chexa_pressure_force
    @property
    def cpyram_pressure_force(self):
        self.deprecated('model.cpyram_pressure_force', 'model.op2_results.force.cpyram_pressure_force', '1.4')
        return self.op2_results.force.cpyram_pressure_force

    # --------------------------------------------------
    # Force - Setter
    @cgap_force.setter
    def cgap_force(self, cgap_force):
        self.deprecated('model.cgap_force', 'model.op2_results.force.cgap_force', '1.4')
        self.op2_results.force.cgap_force = cgap_force
    @cbush_force.setter
    def cbush_force(self, cbush_force):
        self.deprecated('model.cbush_force', 'model.op2_results.force.cbush_force', '1.4')
        self.op2_results.force.cbush_force = cbush_force
    @cweld_force.setter
    def cweld_force(self, cweld_force):
        self.deprecated('model.cweld_force', 'model.op2_results.force.cweld_force', '1.4')
        self.op2_results.force.cweld_force = cweld_force
    @cfast_force.setter
    def cfast_force(self, cfast_force):
        self.deprecated('model.cfast_force', 'model.op2_results.force.cfast_force', '1.4')
        self.op2_results.force.cfast_force = cfast_force
    @cbear_force.setter
    def cbear_force(self, cbear_force):
        self.deprecated('model.cbear_force', 'model.op2_results.force.cbear_force', '1.4')
        self.op2_results.force.cbear_force = cbear_force

    @cshear_force.setter
    def cshear_force(self, cshear_force):
        self.deprecated('model.cshear_force', 'model.op2_results.force.cshear_force', '1.4')
        self.op2_results.force.cshear_force = cshear_force
    @cbend_force.setter
    def cbend_force(self, cbend_force):
        self.deprecated('model.cbend_force', 'model.op2_results.force.cbend_force', '1.4')
        self.op2_results.force.cbend_force = cbend_force
    @cconeax_force.setter
    def cconeax_force(self, cconeax_force):
        self.deprecated('model.cconeax_force', 'model.op2_results.force.cconeax_force', '1.4')
        self.op2_results.force.cconeax_force = cconeax_force

    @property
    def cvisc_force(self):
        self.deprecated('model.cvisc_force', 'model.op2_results.force.cvisc_force', '1.4')
        return self.op2_results.force.cvisc_force
    @cvisc_force.setter
    def cvisc_force(self, cvisc_force):
        self.deprecated('model.cvisc_force', 'model.op2_results.force.cvisc_force', '1.4')
        self.op2_results.force.cvisc_force = cvisc_force

    @property
    def crod_force(self):
        self.deprecated('model.crod_force', 'model.op2_results.force.crod_force', '1.4')
        return self.op2_results.force.crod_force
    @property
    def conrod_force(self):
        self.deprecated('model.conrod_force', 'model.op2_results.force.conrod_force', '1.4')
        return self.op2_results.force.conrod_force
    @property
    def ctube_force(self):
        self.deprecated('model.ctube_force', 'model.op2_results.force.ctube_force', '1.4')
        return self.op2_results.force.ctube_force

    @property
    def cbeam_force(self):
        self.deprecated('model.cbeam_force', 'model.op2_results.force.cbeam_force', '1.4')
        return self.op2_results.force.cbeam_force
    @property
    def cbar_force(self):
        self.deprecated('model.cbar_force', 'model.op2_results.force.cbar_force', '1.4')
        return self.op2_results.force.cbar_force
    @property
    def cbar_force_10nodes(self):
        self.deprecated('model.cbar_force_10nodes', 'model.op2_results.force.cbar_force_10nodes', '1.4')
        return self.op2_results.force.cbar_force

    @property
    def ctria3_force(self):
        self.deprecated('model.cbar_force_10nodes', 'model.op2_results.force.cbar_force_10nodes', '1.4')
        return self.op2_results.force.ctria3_force
    @property
    def ctria6_force(self):
        self.deprecated('model.cbar_force_10nodes', 'model.op2_results.force.cbar_force_10nodes', '1.4')
        return self.op2_results.force.ctria6_force
    @property
    def ctriar_force(self):
        self.deprecated('model.cbar_force_10nodes', 'model.op2_results.force.cbar_force_10nodes', '1.4')
        return self.op2_results.force.ctriar_force
    @property
    def cquad4_force(self):
        self.deprecated('model.cbar_force_10nodes', 'model.op2_results.force.cbar_force_10nodes', '1.4')
        return self.op2_results.force.cquad4_force
    @property
    def cquad8_force(self):
        self.deprecated('model.cbar_force_10nodes', 'model.op2_results.force.cbar_force_10nodes', '1.4')
        return self.op2_results.force.cquad8_force
    @property
    def cquadr_force(self):
        self.deprecated('model.cbar_force_10nodes', 'model.op2_results.force.cbar_force_10nodes', '1.4')
        return self.op2_results.force.cquadr_force

    @crod_force.setter
    def crod_force(self, crod_force):
        self.deprecated('model.cbar_force_10nodes', 'model.op2_results.force.cbar_force_10nodes', '1.4')
        self.op2_results.force.crod_force = crod_force
    @conrod_force.setter
    def conrod_force(self, conrod_force):
        self.deprecated('model.cbar_force_10nodes', 'model.op2_results.force.cbar_force_10nodes', '1.4')
        self.op2_results.force.conrod_force = conrod_force
    @ctube_force.setter
    def ctube_force(self, ctube_force):
        self.deprecated('model.cbar_force_10nodes', 'model.op2_results.force.cbar_force_10nodes', '1.4')
        self.op2_results.force.ctube_force = ctube_force

    @cbeam_force.setter
    def cbeam_force(self, cbeam_force):
        self.deprecated('model.cbar_force_10nodes', 'model.op2_results.force.cbar_force_10nodes', '1.4')
        self.op2_results.force.cbeam_force = cbeam_force
    @cbar_force.setter
    def cbar_force(self, cbar_force):
        self.deprecated('model.cbar_force_10nodes', 'model.op2_results.force.cbar_force_10nodes', '1.4')
        self.op2_results.force.cbar_force = cbar_force
    @cbar_force_10nodes.setter
    def cbar_force_10nodes(self, cbar_force_10nodes):
        self.deprecated('model.cbar_force_10nodes', 'model.op2_results.force.cbar_force_10nodes', '1.4')
        self.op2_results.force.cbar_force = cbar_force_10nodes

    @ctria3_force.setter
    def ctria3_force(self, ctria3_force):
        self.deprecated('model.ctria3_force', 'model.op2_results.force.ctria3_force', '1.4')
        self.op2_results.force.ctria3_force = ctria3_force
    @ctria6_force.setter
    def ctria6_force(self, ctria6_force):
        self.deprecated('model.ctria6_force', 'model.op2_results.force.ctria6_force', '1.4')
        self.op2_results.force.ctria6_force = ctria6_force
    @ctriar_force.setter
    def ctriar_force(self, ctriar_force):
        self.deprecated('model.ctriar_force', 'model.op2_results.force.ctriar_force', '1.4')
        self.op2_results.force.ctriar_force = ctriar_force
    @cquad4_force.setter
    def cquad4_force(self, cquad4_force):
        self.deprecated('model.cquad4_force', 'model.op2_results.force.cquad4_force', '1.4')
        self.op2_results.force.cquad4_force = cquad4_force
    @cquad8_force.setter
    def cquad8_force(self, cquad8_force):
        self.deprecated('model.cquad8_force', 'model.op2_results.force.cquad8_force', '1.4')
        self.op2_results.force.cquad8_force = cquad8_force
    @cquadr_force.setter
    def cquadr_force(self, cquadr_force):
        self.deprecated('model.cquadr_force', 'model.op2_results.force.cquadr_force', '1.4')
        self.op2_results.force.cquadr_force = cquadr_force

    @ctetra_pressure_force.setter
    def ctetra_pressure_force(self, ctetra_pressure_force):
        self.deprecated('model.ctetra_pressure_force', 'model.op2_results.force.ctetra_pressure_force', '1.4')
        self.op2_results.force.ctetra_pressure_force = ctetra_pressure_force
    @cpenta_pressure_force.setter
    def cpenta_pressure_force(self, cpenta_pressure_force):
        self.deprecated('model.cpenta_pressure_force', 'model.op2_results.force.cpenta_pressure_force', '1.4')
        self.op2_results.force.cpenta_pressure_force = cpenta_pressure_force
    @chexa_pressure_force.setter
    def chexa_pressure_force(self, chexa_pressure_force):
        self.deprecated('model.chexa_pressure_force', 'model.op2_results.force.chexa_pressure_force', '1.4')
        self.op2_results.force.chexa_pressure_force = chexa_pressure_force
    @cpyram_pressure_force.setter
    def cpyram_pressure_force(self, cpyram_pressure_force):
        self.deprecated('model.cpyram_pressure_force', 'model.op2_results.force.cpyram_pressure_force', '1.4')
        self.op2_results.force.cpyram_pressure_force = cpyram_pressure_force

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

    def _get_result_length(self, res_types, res_key):
        """
        gets the length of the output data so we can line up:

          RealCRodStrain  - CROD
          RealCTubeStrain - CTUBE

        """
        res_length = 0
        for res_type in res_types:
            if not res_type:
                continue
            key0 = next(iter(res_type))
            if not isinstance(key0, integer_types) and not isinstance(res_key, integer_types):
                if not type(key0) == type(res_key):
                    msg = (
                        'bad compression check...\n'
                        'keys0=%s type(key0)=%s\n'
                        'res_key=%s type(res_key)=%s' % (
                            key0, type(key0), res_key, type(res_key))
                    )
                    raise RuntimeError(msg)

            #print('res_type.keys()=%s' % res_type.keys())
            # res_key_list = res_key[:-1] + [res_key[-1]]
            # res_key = tuple(res_key_list)

            if res_key in res_type:
                # the res_key is
                result = res_type[res_key]
                class_name = result.__class__.__name__
                res_length = max(len(class_name), res_length)
                #print('continue')
                #break
                continue
            elif len(res_type) != 0:
                #print('  not valid')
                # get the 0th key in the dictionary, where key0 is arbitrary
                key0 = get_key0(res_type)
                #print('  key0 = ', key0)

                # extract displacement[0]
                result = res_type[key0]

                # get the class name
                class_name = result.__class__.__name__
                res_length = max(len(class_name), res_length)

                if not is_release:
                    print(' %s - results not found...key=%s' % (class_name, res_key))
            else:  # empty result
                #print('else')
                pass
        #print('res_length =', res_length)
        return res_length

    def remove_empty_results(self):
        table_names = self.get_table_types() + ['grid_point_weight', 'oload_resultant']
        for table_name in table_names:
            obj = self.get_result(table_name)
            if obj is None:
                obj = self.del_result(table_name)
            elif isinstance(obj, dict):
                if len(obj) == 0:
                    obj = self.del_result(table_name)
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
            'displacements', # 'displacements_scaled',
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
            'eigenvalues', # LAMA, CLAMA, BLAMA, LAMAS
            'eigenvalues_fluid', # LAMAF

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

    def get_op2_stats(self, short=False):
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
        return _get_op2_stats(self, short=short)

class Op2F06Attributes(OP2_F06_Common):
    def __init__(self):
        OP2_F06_Common.__init__(self)


def _get_op2_stats(model: OP2, short=False):
    """see OP2.get_op2_stats(...)"""
    msg = []
    #msg += model.op2_results.responses.get_stats(short=short)

    msg.extend(_write_params(model.params))
    for key, weight in model.grid_point_weight.items():
        msg += weight.get_stats(key, short=short)

    msg += model.op2_results.psds.get_stats(short=short)

    table_types = model._get_table_types_testing()

    if short:
        msg += _get_op2_stats_short(model, table_types, model.log)
    else:
        msg += _get_op2_stats_full(model, table_types, model.log)

    if model.matrices:
        msg.append('matrices:\n')
        for unused_name, matrix in sorted(model.matrices.items()):
            #msg.append('matrices[%s].shape = %s\n' % (name, matrix.data.shape))
            msg.append('  ' + str(matrix) + '\n')

    if model.matdicts:
        msg.append('matdicts:\n')
        for unused_name, matrix_dict in sorted(model.matdicts.items()):
            #msg.append('matrices[%s].shape = %s\n' % (name, matrix.data.shape))
            msg.append('  ' + str(matrix_dict) + '\n')
    try:
        return ''.join(msg)
    except TypeError:
        for msgi in msg:
            print('TypeError...%r' % msgi.rstrip())
            assert isinstance(msgi, str), msgi
        raise
    except UnicodeDecodeError:
        for msgi in msg:
            print('UnicodeDecodeError...%r' % msgi.rstrip())
            assert isinstance(msgi, str), msgi
        raise

def _get_op2_stats_short(model: OP2, table_types: list[str],
                         log: SimpleLogger) -> list[str]:
    """helper for get_op2_stats(...)"""
    msg = []
    handled_previously = ['params', 'grid_point_weight', 'psds', 'cstm']
    no_data_classes = ['RealEigenvalues', 'ComplexEigenvalues', 'BucklingEigenvalues']
    for table_type in table_types:
        #table_type_print = ''
        if table_type in handled_previously:
            continue
        if table_type in ['gpdt', 'bgpdt', 'eqexin', 'monitor1', 'monitor3'] or table_type.startswith('responses.'):
            obj = model.get_result(table_type)
            if obj is None:
                continue
            elif isinstance(obj, dict):
                msg.extend(_get_op2_results_stats_dict(obj, table_type, short=True))
                continue
            stats = obj.get_stats(short=True)
            msg.extend(f'op2_results.{table_type}: ' + stats)  # TODO: a hack...not quite right...
            continue

        # # and not table_type.startswith('responses.')
        table_type_print = 'op2_results.' + table_type if '.' in table_type else table_type
        table = model.get_result(table_type)
        if table_type == 'superelement_tables':
            for key in table:
                msg.append(f'{table_type_print}[{key}]\n')
            continue

        try:
            sorted_tables = sorted(table.items(), key=_compare)
        except AttributeError:
            log.warning(f'table_type={table_type}; type(table)={type(table)}')
            raise

        for isubcase, subcase in sorted_tables:
            class_name = subcase.__class__.__name__
            if class_name in no_data_classes:
                msg.append('%s[%r]\n' % (table_type_print, isubcase))
            elif hasattr(subcase, 'data'):
                #data = subcase.data
                #shape = [int(i) for i in subcase.data.shape]
                #headers = subcase.get_headers()
                #headers_str = str(', '.join(headers))
                #msg.append('%s[%s]; %s; %s; [%s]\n' % (
                #table_type, isubcase, class_name, shape, headers_str))
                msg.append('%s[%s]\n' % (table_type_print, isubcase))
            elif table_type == 'params':  #  TODO: remove
                msgi = str(subcase)
            elif hasattr(subcase, 'get_stats'):
                msgi = '%s[%s] # unvectorized\n' % (table_type_print, isubcase)
                msg.append(msgi)
            else:
                msgi = 'skipping %r %s[%s]\n' % (class_name, table_type_print, isubcase)
                msg.append(msgi)
                #raise RuntimeError(msgi)
    return msg

def _get_op2_results_stats_dict(obj: dict[Any, Any],
                                table_type: str, short: bool) -> list[str]:
    msg = []
    for key, obji in obj.items():
        if isinstance(obji, list):
            for iobj, objii in enumerate(obji):
                stats = objii.get_stats(short=short)
                msg.extend(f'op2_results.{table_type}[{key}][{iobj}]: ' + stats)
        else:
            stats = obji.get_stats(short=short)
            msg.extend(f'op2_results.{table_type}[{key}]: ' + stats)
    return msg

def _get_op2_stats_full(model: OP2, table_types: list[str],
                        log: SimpleLogger) -> list[str]:
    """helper for get_op2_stats(...)"""
    msg = []
    handled_previously = ['params', 'grid_point_weight', 'psds']
    for table_type in table_types:
        table = model.get_result(table_type)
        if table_type in handled_previously:
            continue

        skip_results = ('gpdt', 'bgpdt', 'eqexin', 'monitor1', 'monitor3', 'cstm')
        if table_type in skip_results or table_type.startswith('responses.'):
            obj = model.get_result(table_type)
            if obj is None:
                continue
            elif isinstance(obj, dict):
                msg.extend(_get_op2_results_stats_dict(obj, table_type, short=False))
                continue

            stats = obj.get_stats(short=False)
            msg.extend(f'op2_results.{table_type}: ' + stats)  # TODO: a hack...not quite right...
            continue

        table_type_print = 'op2_results.' + table_type if '.' in table_type else table_type
        if table_type == 'superelement_tables':
            for key in table:
                msg.append(f'{table_type_print}[{key}]\n')
            continue

        try:
            for isubcase, subcase in sorted(table.items(), key=_compare):
                class_name = subcase.__class__.__name__
                if hasattr(subcase, 'get_stats'):
                    try:
                        stats = subcase.get_stats() # short=short
                    except Exception:
                        msgi = 'errored reading %s %s[%s]\n\n' % (
                            class_name, table_type_print, isubcase)
                        msg.append(msgi)
                        raise
                    else:
                        msg.append(f'{table_type_print}[{isubcase}]\n')
                        msg.extend(stats)
                        msg.append('\n')
                else:
                    msgi = 'skipping %s %s[%s]\n\n' % (class_name, table_type_print, isubcase)
                    msg.append(msgi)
                    raise RuntimeError(msgi)
        except Exception:
            log.warning(f'table_type={table_type}; type(table)={type(table)}')
            log.warning(str(table))
            raise
    return msg

def _write_params(params: dict[str, Any]):
    """helper for get_op2_stats(...)"""
    if not params:
        return []
    msg = ['params:\n']
    iparam = 0
    for key, param in sorted(params.items()):
        if len(param.values) == 1:
            msg.append(f'  {key} = {param.values[0]!r}\n')
        else:
            msg.append(f'  {key} = {param.values}\n')
        iparam += 1
        if iparam > 10:
            msg.append(f'  ...\n')
            break
    return msg

COMPARE_KEYS = (int, int32, int64, str, bytes)
def _compare(key_value):
    """helper for get_op2_stats(...)"""
    key = key_value[0]
    if isinstance(key, COMPARE_KEYS):
        return key
    #print('key=%s type=%s' % (key, type(key)))
    return key[0]
