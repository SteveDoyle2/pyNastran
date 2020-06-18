from __future__ import annotations
from typing import List, TYPE_CHECKING
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

    def has_result(self, result_name):
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

    def get_result(self, result_name):
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
    # Strain Energy - Getter
    @property
    def ctetra_strain_energy(self):
        self.deprecated('model.ctetra_strain_energy', 'model.op2_results.strain_energy.ctetra_strain_energy', '1.3')
        return self.op2_results.strain_energy.ctetra_strain_energy
    @property
    def chexa_strain_energy(self):
        return self.op2_results.strain_energy.chexa_strain_energy
    @property
    def cpenta_strain_energy(self):
        return self.op2_results.strain_energy.cpenta_strain_energy
    @property
    def cpyram_strain_energy(self):
        return self.op2_results.strain_energy.cpyram_strain_energy

    @property
    def celas1_strain_energy(self):
        return self.op2_results.strain_energy.celas1_strain_energy
    @property
    def celas2_strain_energy(self):
        return self.op2_results.strain_energy.celas2_strain_energy
    @property
    def celas3_strain_energy(self):
        return self.op2_results.strain_energy.celas3_strain_energy
    @property
    def celas4_strain_energy(self):
        return self.op2_results.strain_energy.celas4_strain_energy

    @property
    def crod_strain_energy(self):
        return self.op2_results.strain_energy.crod_strain_energy
    @property
    def ctube_strain_energy(self):
        return self.op2_results.strain_energy.crod_strain_energy
    @property
    def conrod_strain_energy(self):
        return self.op2_results.strain_energy.conrod_strain_energy

    @property
    def cquad4_strain_energy(self):
        return self.op2_results.strain_energy.cquad4_strain_energy
    @property
    def cquad8_strain_energy(self):
        return self.op2_results.strain_energy.cquad8_strain_energy
    @property
    def cquadr_strain_energy(self):
        return self.op2_results.strain_energy.cquadr_strain_energy
    @property
    def cquadx_strain_energy(self):
        return self.op2_results.strain_energy.cquadx_strain_energy

    @property
    def ctria3_strain_energy(self):
        return self.op2_results.strain_energy.ctria3_strain_energy
    @property
    def ctria6_strain_energy(self):
        return self.op2_results.strain_energy.ctria6_strain_energy
    @property
    def ctriar_strain_energy(self):
        return self.op2_results.strain_energy.ctriar_strain_energy
    @property
    def ctriax_strain_energy(self):
        return self.op2_results.strain_energy.ctriax_strain_energy
    @property
    def ctriax6_strain_energy(self):
        return self.op2_results.strain_energy.ctriax6_strain_energy

    @property
    def ctetra_strain_energy(self):
        return self.op2_results.strain_energy.ctetra_strain_energy
    @property
    def cpenta_strain_energy(self):
        return self.op2_results.strain_energy.cpenta_strain_energy
    @property
    def chexa_strain_energy(self):
        return self.op2_results.strain_energy.chexa_strain_energy
    @property
    def cpyram_strain_energy(self):
        return self.op2_results.strain_energy.cpyram_strain_energy

    @property
    def crod_strain_energy(self):
        return self.op2_results.strain_energy.crod_strain_energy
    @property
    def ctube_strain_energy(self):
        return self.op2_results.strain_energy.ctube_strain_energy
    @property
    def conrod_strain_energy(self):
        return self.op2_results.strain_energy.conrod_strain_energy

    @property
    def cbar_strain_energy(self):
        return self.op2_results.strain_energy.cbar_strain_energy
    @property
    def cbeam_strain_energy(self):
        return self.op2_results.strain_energy.cbeam_strain_energy
    @property
    def cbend_strain_energy(self):
        return self.op2_results.strain_energy.cbend_strain_energy

    # ------------------------------------------------------
    # Strain Energy - Getter 2
    @property
    def cgap_strain_energy(self):
        return self.op2_results.strain_energy.cgap_strain_energy
    @property
    def cdamp1_strain_energy(self):
        return self.op2_results.strain_energy.cdamp1_strain_energy
    @property
    def cdamp2_strain_energy(self):
        return self.op2_results.strain_energy.cdamp2_strain_energy
    @property
    def cdamp3_strain_energy(self):
        return self.op2_results.strain_energy.cdamp3_strain_energy
    @property
    def cdamp4_strain_energy(self):
        return self.op2_results.strain_energy.cdamp4_strain_energy
    @property
    def cbush_strain_energy(self):
        return self.op2_results.strain_energy.cbush_strain_energy
    @property
    def dmig_strain_energy(self):
        return self.op2_results.strain_energy.dmig_strain_energy
    @property
    def genel_strain_energy(self):
        return self.op2_results.strain_energy.genel_strain_energy
    @property
    def cshear_strain_energy(self):
        return self.op2_results.strain_energy.cshear_strain_energy
    @property
    def conm2_strain_energy(self):
        return self.op2_results.strain_energy.conm2_strain_energy
    @property
    def cdum8_strain_energy(self):
        return self.op2_results.strain_energy.cdum8_strain_energy
    @property
    def rbe1_strain_energy(self):
        return self.op2_results.strain_energy.rbe1_strain_energy
    @property
    def rbe3_strain_energy(self):
        return self.op2_results.strain_energy.rbe3_strain_energy
    @property
    def weldc_strain_energy(self):
        return self.op2_results.strain_energy.weldc_strain_energy

    # ------------------------------------------------------------------
    # Strain Energy - Setter

    @celas1_strain_energy.setter
    def celas1_strain_energy(self, celas1_strain_energy):
        self.op2_results.strain_energy.celas1_strain_energy = celas1_strain_energy
    @celas2_strain_energy.setter
    def celas2_strain_energy(self, celas2_strain_energy):
        self.op2_results.strain_energy.celas2_strain_energy = celas2_strain_energy
    @celas3_strain_energy.setter
    def celas3_strain_energy(self, celas3_strain_energy):
        self.op2_results.strain_energy.celas3_strain_energy = celas3_strain_energy
    @celas4_strain_energy.setter
    def celas4_strain_energy(self, celas4_strain_energy):
        self.op2_results.strain_energy.celas4_strain_energy = celas4_strain_energy


    @cdamp1_strain_energy.setter
    def cdamp1_strain_energy(self, cdamp1_strain_energy):
        self.op2_results.strain_energy.cdamp1_strain_energy = cdamp1_strain_energy
    @cdamp2_strain_energy.setter
    def cdamp2_strain_energy(self, cdamp2_strain_energy):
        self.op2_results.strain_energy.cdamp2_strain_energy = cdamp2_strain_energy
    @cdamp3_strain_energy.setter
    def cdamp3_strain_energy(self, cdamp3_strain_energy):
        self.op2_results.strain_energy.cdamp3_strain_energy = cdamp3_strain_energy
    @cdamp4_strain_energy.setter
    def cdamp4_strain_energy(self, cdamp4_strain_energy):
        self.op2_results.strain_energy.cdamp4_strain_energy = cdamp4_strain_energy

    @cgap_strain_energy.setter
    def cgap_strain_energy(self, cgap_strain_energy):
        self.op2_results.strain_energy.cgap_strain_energy = cgap_strain_energy
    @cbush_strain_energy.setter
    def cbush_strain_energy(self, cbush_strain_energy):
        self.op2_results.strain_energy.cbush_strain_energy = cbush_strain_energy

    @crod_strain_energy.setter
    def crod_strain_energy(self, crod_strain_energy):
        self.op2_results.strain_energy.crod_strain_energy = crod_strain_energy
    @ctube_strain_energy.setter
    def ctube_strain_energy(self, ctube_strain_energy):
        self.op2_results.strain_energy.crod_strain_energy = ctube_strain_energy
    @conrod_strain_energy.setter
    def conrod_strain_energy(self, conrod_strain_energy):
        self.op2_results.strain_energy.conrod_strain_energy = conrod_strain_energy

    @cquad4_strain_energy.setter
    def cquad4_strain_energy(self, cquad4_strain_energy):
        self.op2_results.strain_energy.cquad4_strain_energy = cquad4_strain_energy
    @cquad8_strain_energy.setter
    def cquad8_strain_energy(self, cquad8_strain_energy):
        self.op2_results.strain_energy.cquad8_strain_energy = cquad8_strain_energy
    @cquadr_strain_energy.setter
    def cquadr_strain_energy(self, cquadr_strain_energy):
        self.op2_results.strain_energy.cquadr_strain_energy = cquadr_strain_energy
    @cquadx_strain_energy.setter
    def cquadx_strain_energy(self, cquadx_strain_energy):
        self.op2_results.strain_energy.cquadx_strain_energy = cquadx_strain_energy

    @ctetra_strain_energy.setter
    def ctria3_strain_energy(self, ctria3_strain_energy):
        self.op2_results.strain_energy.ctria3_strain_energy = ctria3_strain_energy
    @ctria6_strain_energy.setter
    def ctria6_strain_energy(self, ctria6_strain_energy):
        self.op2_results.strain_energy.ctria6_strain_energy = ctria6_strain_energy
    @ctriar_strain_energy.setter
    def ctriar_strain_energy(self, ctriar_strain_energy):
        self.op2_results.strain_energy.ctriar_strain_energy = ctriar_strain_energy
    @ctetra_strain_energy.setter
    def ctriax_strain_energy(self, ctriax_strain_energy):
        self.op2_results.strain_energy.ctriax_strain_energy = ctriax_strain_energy
    @ctetra_strain_energy.setter
    def ctriax6_strain_energy(self, ctriax6_strain_energy):
        self.op2_results.strain_energy.ctriax6_strain_energy = ctriax6_strain_energy

    @ctetra_strain_energy.setter
    def ctetra_strain_energy(self, ctetra_strain_energy):
        self.op2_results.strain_energy.ctetra_strain_energy = ctetra_strain_energy
    @chexa_strain_energy.setter
    def chexa_strain_energy(self, chexa_strain_energy):
        self.op2_results.strain_energy.chexa_strain_energy = chexa_strain_energy
    @cpenta_strain_energy.setter
    def cpenta_strain_energy(self, cpenta_strain_energy):
        self.op2_results.strain_energy.cpenta_strain_energy = cpenta_strain_energy
    @cpyram_strain_energy.setter
    def cpyram_strain_energy(self, cpyram_strain_energy):
        self.op2_results.strain_energy.cpyram_strain_energy = cpyram_strain_energy

    # ------------------------------------------------------------------
    # Force - Getter
    @property
    def celas1_force(self):
        return self.op2_results.force.celas1_force
    @property
    def celas2_force(self):
        return self.op2_results.force.celas2_force
    @property
    def celas3_force(self):
        return self.op2_results.force.celas3_force
    @property
    def celas4_force(self):
        return self.op2_results.force.celas4_force

    @property
    def cdamp1_force(self):
        return self.op2_results.force.cdamp1_force
    @property
    def cdamp2_force(self):
        return self.op2_results.force.cdamp2_force
    @property
    def cdamp3_force(self):
        return self.op2_results.force.cdamp3_force
    @property
    def cdamp4_force(self):
        return self.op2_results.force.cdamp4_force

    # ------------------------------------------------------------------
    # Force - Setter
    @celas1_force.setter
    def celas1_force(self, celas1_force):
        self.op2_results.force.celas1_force = celas1_force
    @celas2_force.setter
    def celas2_force(self, celas2_force):
        self.op2_results.force.celas2_force = celas2_force
    @celas3_force.setter
    def celas3_force(self, celas3_force):
        self.op2_results.force.celas3_force = celas3_force
    @celas4_force.setter
    def celas4_force(self, celas4_force):
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
        self.op2_results.force.cdamp1_force = cdamp1_force
    @cdamp2_force.setter
    def cdamp2_force(self, cdamp2_force):
        self.op2_results.force.cdamp2_force = cdamp2_force
    @cdamp3_force.setter
    def cdamp3_force(self, cdamp3_force):
        self.op2_results.force.cdamp3_force = cdamp3_force
    @cdamp4_force.setter
    def cdamp4_force(self, cdamp4_force):
        self.op2_results.force.cdamp4_force = cdamp4_force

    @property
    def cgap_force(self):
        return self.op2_results.force.cgap_force
    @property
    def cbush_force(self):
        return self.op2_results.force.cbush_force
    @property
    def cweld_force(self):
        return self.op2_results.force.cweld_force
    @property
    def cfast_force(self):
        return self.op2_results.force.cfast_force
    @property
    def cbear_force(self):
        return self.op2_results.force.cbear_force

    @property
    def cshear_force(self):
        return self.op2_results.force.cshear_force
    @property
    def cbend_force(self):
        return self.op2_results.force.cbend_force
    @property
    def cconeax_force(self):
        return self.op2_results.force.cconeax_force

    @property
    def ctetra_pressure_force(self):
        return self.op2_results.force.ctetra_pressure_force
    @property
    def cpenta_pressure_force(self):
        return self.op2_results.force.cpenta_pressure_force
    @property
    def chexa_pressure_force(self):
        return self.op2_results.force.chexa_pressure_force
    @property
    def cpyram_pressure_force(self):
        return self.op2_results.force.cpyram_pressure_force

    # --------------------------------------------------
    # Force - Setter
    @cgap_force.setter
    def cgap_force(self, cgap_force):
        self.op2_results.force.cgap_force = cgap_force
    @cbush_force.setter
    def cbush_force(self, cbush_force):
        self.op2_results.force.cbush_force = cbush_force
    @cweld_force.setter
    def cweld_force(self, cweld_force):
        self.op2_results.force.cweld_force = cweld_force
    @cfast_force.setter
    def cfast_force(self, cfast_force):
        self.op2_results.force.cfast_force = cfast_force
    @cbear_force.setter
    def cbear_force(self, cbear_force):
        self.op2_results.force.cbear_force = cbear_force

    @cshear_force.setter
    def cshear_force(self, cshear_force):
        self.op2_results.force.cshear_force = cshear_force
    @cbend_force.setter
    def cbend_force(self, cbend_force):
        self.op2_results.force.cbend_force = cbend_force
    @cconeax_force.setter
    def cconeax_force(self, cconeax_force):
        self.op2_results.force.cconeax_force = cconeax_force

    @property
    def cvisc_force(self):
        return self.op2_results.force.vu_tria_force
    @property
    def vu_tria_force(self):
        return self.op2_results.force.vu_tria_force
    @property
    def vu_quad_force(self):
        return self.op2_results.force.vu_quad_force
    @property
    def cbeam_force_vu(self):
        return self.op2_results.force.cbeam_force_vu

    @cvisc_force.setter
    def cvisc_force(self, cvisc_force):
        self.op2_results.force.vu_tria_force = cvisc_force
    @vu_tria_force.setter
    def vu_tria_force(self, vu_tria_force):
        self.op2_results.force.vu_tria_force = vu_tria_force
    @vu_quad_force.setter
    def vu_quad_force(self, vu_quad_force):
        self.op2_results.force.vu_quad_force = vu_quad_force
    @cbeam_force_vu.setter
    def cbeam_force_vu(self, cbeam_force_vu):
        self.op2_results.force.cbeam_force_vu = cbeam_force_vu

    @property
    def crod_force(self):
        return self.op2_results.force.crod_force
    @property
    def conrod_force(self):
        return self.op2_results.force.conrod_force
    @property
    def ctube_force(self):
        return self.op2_results.force.ctube_force

    @property
    def cbeam_force(self):
        return self.op2_results.force.cbeam_force
    @property
    def cbar_force(self):
        return self.op2_results.force.cbar_force
    @property
    def cbar_force_10nodes(self):
        return self.op2_results.force.cbar_force

    @property
    def ctria3_force(self):
        return self.op2_results.force.ctria3_force
    @property
    def ctria6_force(self):
        return self.op2_results.force.ctria6_force
    @property
    def ctriar_force(self):
        return self.op2_results.force.ctriar_force
    @property
    def cquad4_force(self):
        return self.op2_results.force.cquad4_force
    @property
    def cquad8_force(self):
        return self.op2_results.force.cquad8_force
    @property
    def cquadr_force(self):
        return self.op2_results.force.cquadr_force

    @crod_force.setter
    def crod_force(self, crod_force):
        self.op2_results.force.crod_force = crod_force
    @conrod_force.setter
    def conrod_force(self, conrod_force):
        self.op2_results.force.conrod_force = conrod_force
    @ctube_force.setter
    def ctube_force(self, ctube_force):
        self.op2_results.force.ctube_force = ctube_force

    @cbeam_force.setter
    def cbeam_force(self, cbeam_force):
        self.op2_results.force.cbeam_force = cbeam_force
    @cbar_force.setter
    def cbar_force(self, cbar_force):
        self.op2_results.force.cbar_force = cbar_force
    @cbar_force_10nodes.setter
    def cbar_force_10nodes(self, cbar_force_10nodes):
        self.op2_results.force.cbar_force = cbar_force_10nodes

    @ctria3_force.setter
    def ctria3_force(self, ctria3_force):
        self.op2_results.force.ctria3_force = ctria3_force
    @ctria6_force.setter
    def ctria6_force(self, ctria6_force):
        self.op2_results.force.ctria6_force = ctria6_force
    @ctriar_force.setter
    def ctriar_force(self, ctriar_force):
        self.op2_results.force.ctriar_force = ctriar_force
    @cquad4_force.setter
    def cquad4_force(self, cquad4_force):
        self.op2_results.force.cquad4_force = cquad4_force
    @cquad8_force.setter
    def cquad8_force(self, cquad8_force):
        self.op2_results.force.cquad8_force = cquad8_force
    @cquadr_force.setter
    def cquadr_force(self, cquadr_force):
        self.op2_results.force.cquadr_force = cquadr_force

    @ctetra_pressure_force.setter
    def ctetra_pressure_force(self, ctetra_pressure_force):
        self.op2_results.force.ctetra_pressure_force = ctetra_pressure_force
    @cpenta_pressure_force.setter
    def cpenta_pressure_force(self, cpenta_pressure_force):
        self.op2_results.force.cpenta_pressure_force = cpenta_pressure_force
    @chexa_pressure_force.setter
    def chexa_pressure_force(self, chexa_pressure_force):
        self.op2_results.force.chexa_pressure_force = chexa_pressure_force
    @cpyram_pressure_force.setter
    def cpyram_pressure_force(self, cpyram_pressure_force):
        self.op2_results.force.cpyram_pressure_force = cpyram_pressure_force

    # ------------------------------------------------------------------
    # Stress - Getter
    @property
    def celas1_stress(self):
        return self.op2_results.stress.celas1_stress
    @property
    def celas2_stress(self):
        return self.op2_results.stress.celas2_stress
    @property
    def celas3_stress(self):
        return self.op2_results.stress.celas3_stress
    @property
    def celas4_stress(self):
        return self.op2_results.stress.celas4_stress

    @property
    def ctetra_stress(self):
        return self.op2_results.stress.ctetra_stress
    @property
    def cpenta_stress(self):
        return self.op2_results.stress.cpenta_stress
    @property
    def chexa_stress(self):
        return self.op2_results.stress.chexa_stress
    @property
    def cpyram_stress(self):
        return self.op2_results.stress.cpyram_stress
    # ------------------------------------------------------------------
    # Stress - Setter
    @celas1_stress.setter
    def celas1_strain(self, celas1_stress):
        self.op2_results.stress.celas1_stress = celas1_stress
    @celas2_stress.setter
    def celas2_stress(self, celas2_stress):
        self.op2_results.stress.celas2_stress = celas2_stress
    @celas3_stress.setter
    def celas3_stress(self, celas3_stress):
        self.op2_results.stress.celas3_stress = celas3_stress
    @celas4_stress.setter
    def celas4_stress(self, celas4_stress):
        self.op2_results.stress.celas4_stress = celas4_stress

    @ctetra_stress.setter
    def ctetra_stress(self, ctetra_stress):
        self.op2_results.stress.ctetra_stress = ctetra_stress
    @cpenta_stress.setter
    def cpenta_stress(self, cpenta_stress):
        self.op2_results.stress.cpenta_stress = cpenta_stress
    @chexa_stress.setter
    def chexa_stress(self, chexa_stress):
        self.op2_results.stress.chexa_stress = chexa_stress
    @cpyram_stress.setter
    def cpyram_stress(self, cpyram_stress):
        self.op2_results.stress.cpyram_stress = cpyram_stress
    # ------------------------------------------------------------------
    # Strain - Getter
    @property
    def celas1_strain(self):
        return self.op2_results.strain.celas1_strain
    @property
    def celas2_strain(self):
        return self.op2_results.strain.celas2_strain
    @property
    def celas3_strain(self):
        return self.op2_results.strain.celas3_strain
    @property
    def celas4_strain(self):
        return self.op2_results.strain.celas4_strain

    @property
    def ctetra_strain(self):
        return self.op2_results.strain.ctetra_strain
    @property
    def cpenta_strain(self):
        return self.op2_results.strain.cpenta_strain
    @property
    def chexa_strain(self):
        return self.op2_results.strain.chexa_strain
    @property
    def cpyram_strain(self):
        return self.op2_results.strain.cpyram_strain

    # ------------------------------------------------------------------
    # Strain - Setter
    @celas1_strain.setter
    def celas1_strain(self, celas1_strain):
        self.op2_results.strain.celas1_strain = celas1_strain
    @celas2_strain.setter
    def celas2_strain(self, celas2_strain):
        self.op2_results.strain.celas2_strain = celas2_strain
    @celas3_strain.setter
    def celas3_strain(self, celas3_strain):
        self.op2_results.strain.celas3_strain = celas3_strain
    @celas4_strain.setter
    def celas4_strain(self, celas4_strain):
        self.op2_results.strain.celas4_strain = celas4_strain

    @ctetra_strain.setter
    def ctetra_strain(self, ctetra_strain):
        self.op2_results.strain.ctetra_strain = ctetra_strain
    @cpenta_strain.setter
    def cpenta_strain(self, cpenta_strain):
        self.op2_results.strain.cpenta_strain = cpenta_strain
    @chexa_strain.setter
    def chexa_strain(self, chexa_strain):
        self.op2_results.strain.chexa_strain = chexa_strain
    @cpyram_strain.setter
    def cpyram_strain(self, cpyram_strain):
        self.op2_results.strain.cpyram_strain = cpyram_strain

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

        self.crod_stress = {}
        self.conrod_stress = {}
        self.ctube_stress = {}

        self.crod_strain = {}
        self.conrod_strain = {}
        self.ctube_strain = {}

        #======================================================================

        self.nonlinear_ctetra_stress_strain = {}
        self.nonlinear_cpenta_stress_strain = {}
        self.nonlinear_chexa_stress_strain = {}

        #======================================================================

        # bars/beams
        self.cbar_force_abs = {} # thermal=2
        self.cbar_force_srss = {} # thermal=4
        self.cbar_force_nrl = {} # thermal=8

        self.cbar_stress = {}
        self.cbar_strain = {}

        self.cbar_stress_10nodes = {}
        self.cbar_strain_10nodes = {}

        self.cbeam_stress = {}
        self.cbeam_strain = {}

        #======================================================================
        self.cbend_stress = {}
        self.cbend_strain = {}

        #======================================================================
        # shells

        self.ctria3_stress = {}
        self.ctria6_stress = {}
        self.cquad4_stress = {}
        self.cquad8_stress = {}
        self.cquadr_stress = {}
        self.ctriar_stress = {}

        self.ctria3_strain = {}
        self.ctria6_strain = {}
        self.cquad4_strain = {}
        self.cquad8_strain = {}
        self.cquadr_strain = {}
        self.ctriar_strain = {}

        self.cquad4_composite_stress = {}
        self.cquad8_composite_stress = {}
        self.cquadr_composite_stress = {}
        self.ctria3_composite_stress = {}
        self.ctria6_composite_stress = {}
        self.ctriar_composite_stress = {}

        self.cquad4_composite_strain = {}
        self.cquad8_composite_strain = {}
        self.cquadr_composite_strain = {}
        self.ctria3_composite_strain = {}
        self.ctria6_composite_strain = {}
        self.ctriar_composite_strain = {}

        self.cplstn3_stress = {}
        self.cplstn4_stress = {}
        self.cplstn6_stress = {}
        self.cplstn8_stress = {}
        self.cplsts3_stress = {}
        self.cplsts4_stress = {}
        self.cplsts6_stress = {}
        self.cplsts8_stress = {}

        self.cplstn3_strain = {}
        self.cplstn4_strain = {}
        self.cplstn6_strain = {}
        self.cplstn8_strain = {}
        self.cplsts3_strain = {}
        self.cplsts4_strain = {}
        self.cplsts6_strain = {}
        self.cplsts8_strain = {}

        self.cshear_stress = {}
        self.cshear_strain = {}

        #: OES - CBEAM 94
        self.nonlinear_cbeam_stress = {}

        # bushing
        self.cbush_stress = {}
        self.cbush_strain = {}
        self.cbush1d_stress_strain = {}

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
        self.displacements_scaled = {}    # tCode=1 thermal=8

        #: OUP

        self.displacement_scaled_response_spectra_nrl = {}  # thermal=8
        self.displacement_scaled_response_spectra_abs = {}  # thermal=2
        self.displacement_scaled_response_spectra_srss = {} # thermal=4
        #self.displacement_scaled_response_spectra_psd = {}
        #self.displacement_scaled_response_spectra_ato = {}
        #self.displacement_scaled_response_spectra_rms = {}
        #self.displacement_scaled_response_spectra_crm = {}
        #self.displacement_scaled_response_spectra_no = {}


        #: OUG - velocity
        self.velocities = {}              # tCode=10 thermal=0

        #self.velocity_scaled_response_spectra_nrl = {}
        self.velocity_scaled_response_spectra_abs = {}
        #self.velocity_scaled_response_spectra_psd = {}
        #self.velocity_scaled_response_spectra_ato = {}
        #self.velocity_scaled_response_spectra_rms = {}
        #self.velocity_scaled_response_spectra_crm = {}
        #self.velocity_scaled_response_spectra_no = {}

        #: OUG - acceleration
        self.accelerations = {}            # tCode=11 thermal=0

        self.acceleration_scaled_response_spectra_nrl = {}
        self.acceleration_scaled_response_spectra_abs = {}
        #self.acceleration_scaled_response_spectra_psd = {}
        #self.acceleration_scaled_response_spectra_ato = {}
        #self.acceleration_scaled_response_spectra_rms = {}
        #self.acceleration_scaled_response_spectra_crm = {}
        #self.acceleration_scaled_response_spectra_no = {}


        #: OUG - temperatures
        self.temperatures = {}           # tCode=1 thermal=1

        #: OUG - eigenvectors
        self.eigenvectors = {}            # tCode=7 thermal=0

        # OES - tCode=5 thermal=0 s_code=0,1 (stress/strain)

        #: OES - CTRIAX6
        self.ctriax_stress = {}
        self.ctriax_strain = {}

        #: OES - nonlinear CROD/CONROD/CTUBE stress/strain
        self.nonlinear_crod_stress = {}
        self.nonlinear_crod_strain = {}

        self.nonlinear_ctube_stress = {}
        self.nonlinear_ctube_strain = {}

        self.nonlinear_conrod_stress = {}
        self.nonlinear_conrod_strain = {}

        #: OESNLXR - CTRIA3/CQUAD4 strain
        self.hyperelastic_cquad4_strain = {}

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
        self.spc_forces_scaled_response_spectra_nrl = {}

        self.mpc_forces = {}  # tCode=39
        self.mpc_forces_RAQCONS = {}
        self.mpc_forces_RAQEATC = {}

        # OQG - thermal forces
        self.thermal_gradient_and_flux = {}

        #: OGF - grid point forces
        self.grid_point_forces = {}  # tCode=19

        #: OGS1 - grid point stresses
        self.grid_point_surface_stresses = {}       # tCode=26
        self.grid_point_stresses_volume_direct = {}  # tCode=27
        self.grid_point_stresses_volume_principal = {} # tCode=28
        self.grid_point_stress_discontinuities = {}  # tCode=35

        #: OPG - summation of loads for each element
        self.load_vectors = {}       # OPG1; tCode=2  thermal=0
        self.load_vectors_v = {}     # OPGV1
        self.thermal_load_vectors = {}  # tCode=2  thermal=1
        self.applied_loads = {}       # tCode=19 thermal=0
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
            'displacements', 'displacements_scaled',
            'temperatures',
            'eigenvectors',
            'velocities',
            'accelerations',

            # OQG - spc/mpc forces
            'spc_forces', 'spc_forces_v', 'spc_forces_scaled_response_spectra_nrl',
            'mpc_forces', 'mpc_forces_RAQCONS', 'mpc_forces_RAQEATC',
            'thermal_gradient_and_flux',

            # OGF - grid point forces
            'grid_point_forces',

            # OPG - summation of loads for each element '
            'load_vectors', 'load_vectors_v',
            'thermal_load_vectors',
            'applied_loads',
            'force_vectors',

            # OES - isotropic CROD/CONROD/CTUBE stress/strain
            'crod_stress', 'conrod_stress', 'ctube_stress',
            'crod_strain', 'conrod_strain', 'ctube_strain',

            # OES - isotropic CBAR stress/strain
            'cbar_stress', 'cbar_strain',
            #'cbar_force',
            'cbar_force_abs', 'cbar_force_nrl', 'cbar_force_srss',

            'cbar_stress_10nodes', 'cbar_strain_10nodes', # 'cbar_force_10nodes',

            # OES - isotropic CBEAM stress/strain
            'cbeam_stress', 'cbeam_strain',
            'nonlinear_cbeam_stress',
            #'nonlinear_cbeam_strain',
            'nonlinear_cpyram_stress_strain',

            # CBEND - isotropic CBEAM stress/strain
            'cbend_stress', 'cbend_strain', # 'cbend_force',

            # OES - isotropic CTRIA3/CQUAD4 stress/strain
            'ctria3_stress', 'ctriar_stress', 'ctria6_stress',
            'cquadr_stress', 'cquad4_stress', 'cquad8_stress',
            'ctria3_strain', 'ctriar_strain', 'ctria6_strain',
            'cquadr_strain', 'cquad4_strain', 'cquad8_strain',

            # OES - composite CTRIA3/CQUAD4 stress/strain
            'cquad4_composite_stress', 'cquad8_composite_stress', 'cquadr_composite_stress',
            'ctria3_composite_stress', 'ctria6_composite_stress', 'ctriar_composite_stress',
            'cquad4_composite_strain', 'cquad8_composite_strain', 'cquadr_composite_strain',
            'ctria3_composite_strain', 'ctria6_composite_strain', 'ctriar_composite_strain',

            # OES - CSHEAR stress/strain
            'cshear_stress', 'cshear_strain',

            'cplstn3_stress', 'cplstn4_stress', 'cplstn6_stress', 'cplstn8_stress',
            'cplsts3_stress', 'cplsts4_stress', 'cplsts6_stress', 'cplsts8_stress',

            'cplstn3_strain', 'cplstn4_strain', 'cplstn6_strain', 'cplstn8_strain',
            'cplsts3_strain', 'cplsts4_strain', 'cplsts6_strain', 'cplsts8_strain',
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

            # HISADD
            #'convergence_history',

            # R1TABRG
            #'response1_table',
        ]
        table_types += [
            # OES - CTRIAX6
            'ctriax_stress', 'ctriax_strain',
            'cbush_stress', 'cbush_strain',

            #'cbush_stress',
            #'cbush_strain',
            'cbush1d_stress_strain',

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

            'hyperelastic_cquad4_strain',

            # OES - CEALS1 224, CELAS3 225
            'nonlinear_celas1_stress',
            'nonlinear_celas3_stress',

            # OGS1 - grid point stresses
            'grid_point_surface_stresses', # tCode=26
            'grid_point_stresses_volume_direct',  # tCode=27 # volume direct
            'grid_point_stresses_volume_principal', # tCode =28
            'grid_point_stress_discontinuities',  # tCode=35,


            # unused?
            'displacement_scaled_response_spectra_nrl',
            'displacement_scaled_response_spectra_abs',
            'displacement_scaled_response_spectra_srss',
            'velocity_scaled_response_spectra_abs',
            'acceleration_scaled_response_spectra_nrl',
            'acceleration_scaled_response_spectra_abs',
        ]
        utables = unique(table_types)
        if len(table_types) != len(utables):
            msg = 'Non-unique tables: '
            for i, table_type in enumerate(table_types):
                if table_type in table_types[i+1:]:
                    msg += table_type + ', '
            raise AssertionError(msg)

        return table_types

    def _get_table_types_testing(self):
        """testing method...don't use"""
        skipped_attributes = [
            'card_count', 'data_code', 'element_mapper', 'isubcase_name_map',
            'labels', 'subtitles', 'additional_matrices', 'matrices', 'matdicts',
            'subcase_key', 'end_options', 'expected_times', 'generalized_tables',
            'op2_reader', 'table_count']

        table_types = self.get_table_types()
        tables = object_attributes(self, 'public', filter_properties=True)
        tables_with_properties = object_attributes(self, 'public', filter_properties=False)
        properties = set(tables_with_properties) - set(tables)

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
    except UnicodeDecodeError:
        for msgi in msg:
            print('UnicodeDecodeError...%r' % msgi.rstrip())
            assert isinstance(msgi, str), msgi
        raise

def _get_op2_stats_short(model: OP2, table_types: List[str], log) -> List[str]:
    """helper for get_op2_stats(...)"""
    msg = []
    handled_previously = ['params', 'grid_point_weight', 'psds']
    no_data_classes = ['RealEigenvalues', 'ComplexEigenvalues', 'BucklingEigenvalues']
    for table_type in table_types:
        #table_type_print = ''
        if table_type in handled_previously:
            continue
        if table_type in ['gpdt', 'bgpdt', 'eqexin'] or table_type.startswith('responses.'):
            obj = model.get_result(table_type)
            if obj is None:
                continue
            stats = obj.get_stats(short=True)
            msg.extend(f'op2_results.{table_type}: ' + stats)  # TODO: a hack...not quite right...
            continue

        # # and not table_type.startswith('responses.')
        table_type_print = 'op2_results.' + table_type if '.' in table_type else table_type
        table = model.get_result(table_type)
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

def _get_op2_stats_full(model: OP2, table_types: List[str], log):
    """helper for get_op2_stats(...)"""
    msg = []
    handled_previously = ['params', 'grid_point_weight', 'psds']
    for table_type in table_types:
        table = model.get_result(table_type)
        if table_type in handled_previously:
            continue
        if table_type in ['gpdt', 'bgpdt', 'eqexin'] or table_type.startswith('responses.'):
            obj = model.get_result(table_type)
            if obj is None:
                continue
            #elif isinstance(obj, dict):
                #print(obj)
            stats = obj.get_stats(short=False)
            msg.extend(f'op2_results.{table_type}: ' + stats)  # TODO: a hack...not quite right...
            continue

        table_type_print = 'op2_results.' + table_type if '.' in table_type else table_type
        try:
            for isubcase, subcase in sorted(table.items(), key=_compare):
                class_name = subcase.__class__.__name__
                if hasattr(subcase, 'get_stats'):
                    try:
                        stats = subcase.get_stats() # short=short
                    except:
                        msgi = 'errored reading %s %s[%s]\n\n' % (
                            class_name, table_type_print, isubcase)
                        msg.append(msgi)
                        raise
                    else:
                        msg.append('%s[%s]\n' % (table_type_print, isubcase))
                        msg.extend(stats)
                        msg.append('\n')
                else:
                    msgi = 'skipping %s %s[%s]\n\n' % (class_name, table_type_print, isubcase)
                    msg.append(msgi)
                    raise RuntimeError(msgi)
        except:
            log.warning(f'table_type={table_type}; type(table)={type(table)}')
            log.warning(str(table))
            raise
    return msg

def _write_params(params):
    """helper for get_op2_stats(...)"""
    if not params:
        return []
    msg = ['params:\n']
    for key, param in sorted(params.items()):
        if len(param.values) == 1:
            msg.append(f'  {key} = {param.values[0]!r}\n')
        else:
            msg.append(f'  {key} = {param.values}\n')
    return msg

COMPARE_KEYS = (int, int32, int64, str, bytes)
def _compare(key_value):
    """helper for get_op2_stats(...)"""
    key = key_value[0]
    if isinstance(key, COMPARE_KEYS):
        return key
    #print('key=%s type=%s' % (key, type(key)))
    return key[0]
