from typing import Any
import numpy as np

from pyNastran.op2.op2_interface.random_results import (
    RADCONS, RAECONS, RASCONS, RAPCONS, RAFCONS, RAGCONS, RANCONS, RARCONS, RAQCONS,
    RADEATC, RAEEATC, RASEATC, RAPEATC, RAFEATC, RAGEATC, RANEATC, RAREATC, RAQEATC,
    ROUGV1, ROQGM1,
    RADEFFM, SRSS, ABS, NRL,
    AutoCorrelationObjects, PowerSpectralDensityObjects, RootMeansSquareObjects,
    CumulativeRootMeansSquareObjects, NumberOfCrossingsObjects,
    PSDObjects,
)
from pyNastran.op2.result_objects.design_response import Responses

class Results:
    """storage object for even more op2_results (see op2.op2_results)"""
    def __init__(self):
        self.eqexin = None
        self.gpdt = None
        self.bgpdt = None
        self.cddata = []
        self.monitor1 = None
        self.monitor3 = None
        self.responses = Responses()

        self.separation_initial = {}
        self.separation_final = {}
        self.contact_slide_distance = {}
        self.glue_contact_slide_distance = {}
        self.contact_stress = {}
        self.contact_displacements = {}

        # bolts
        self.bolt_results = {}

        self.psds = PSDObjects()
        self.ato = AutoCorrelationObjects()
        self.psd = PowerSpectralDensityObjects()
        self.rms = RootMeansSquareObjects()
        self.no = NumberOfCrossingsObjects()
        self.crm = CumulativeRootMeansSquareObjects()

        self.acoustic = Acoustic()
        self.modal_contribution = ModalContribution()
        self.solution_set = SolutionSet()
        self.strength_ratio = StrengthRatio()
        self.failure_indices = FailureIndices()
        self.force = Force()
        self.thermal_load = ThermalLoad()
        self.stress = Stress('stress')
        self.strain = Strain()
        self.stressa = Stress('stressa')

        self.strain_energy = StrainEnergy()
        self.ROUGV1 = ROUGV1()   # relative disp/vel/acc/eigenvectors
        self.ROQGM1 = ROQGM1()   # relative mpc forces???

        self.RADEFFM = RADEFFM() # eigenvectors

        self.RADCONS = RADCONS() # eigenvectors
        self.RAFCONS = RAFCONS() # force
        self.RASCONS = RASCONS() # stress
        self.RAECONS = RAECONS() # strain
        self.RAGCONS = RAGCONS() # grid point forces
        self.RAPCONS = RAPCONS() # composite stress
        self.RANCONS = RANCONS() # strain energy
        self.RARCONS = RARCONS() # spc force
        self.RAQCONS = RAQCONS() # mpc force

        self.RADEATC = RADEATC() # eigenvectors
        self.RAFEATC = RAFEATC() # force
        self.RASEATC = RASEATC() # stress
        self.RAEEATC = RAEEATC() # strain
        self.RAGEATC = RAGEATC() # grid point forces
        self.RAPEATC = RAPEATC() # composite stress
        self.RANEATC = RANEATC() # strain energy
        self.RAREATC = RAREATC() # spc force
        self.RAQEATC = RAQEATC() # mpcforce
        self.srss = SRSS()
        self.abs = ABS()
        self.nrl = NRL()

        self.cstm = None
        self.trmbd = {}
        self.trmbu = {}
        self.superelement_tables = {}

    def _get_sum_objects_map(self):
        sum_objs = {
            'acoustic' : self.acoustic,
            'responses' : self.responses,
            'force' : self.force,
            'thermal_load' : self.thermal_load,
            'strain_energy' : self.strain_energy,
            'stress': self.stress,
            'strain': self.strain,
            'stressa': self.stressa,
            #self.ato,
            #self.psd,
            #self.rms,
            #self.no,
            #self.crm,
            #self.modal_contribution,
            #self.strength_ratio,
            #self.failure_indices,
            #self.solution_set,
            #self.ROUGV1,
            #self.RADEFFM,
            #self.RADCONS, self.RAFCONS, self.RASCONS, self.RAECONS, self.RAGCONS, self.RAPCONS, self.RANCONS,
            #self.RADEATC, self.RAFEATC, self.RASEATC, self.RAEEATC, self.RAGEATC, self.RAPEATC, self.RANEATC,
        }
        return sum_objs

    def _get_sum_objects(self) -> list[Any]:
        sum_objs = [
            self.acoustic,
            self.responses,
            self.force, self.thermal_load,
            self.stress, self.strain,
            self.stressa,
            self.strain_energy,
            self.ato, self.psd, self.rms, self.no, self.crm,
            self.modal_contribution, self.strength_ratio, self.failure_indices,
            self.solution_set,
            self.ROUGV1, self.ROQGM1,
            self.RADEFFM,
            self.RADCONS, self.RAFCONS, self.RASCONS, self.RAECONS, self.RAGCONS, self.RAPCONS, self.RANCONS, self.RARCONS, self.RAQCONS,
            self.RADEATC, self.RAFEATC, self.RASEATC, self.RAEEATC, self.RAGEATC, self.RAPEATC, self.RANEATC, self.RAREATC, self.RAQEATC,

            self.srss, self.abs, self.nrl,
        ]
        return sum_objs

    def _get_base_objects_map(self) -> dict[str, Any]:
        """gets only the objects that are do not contain sub-objects"""
        base_names = [
            'eqexin', 'gpdt', 'bgpdt', 'psds', 'monitor1', 'monitor3',
            'separation_initial', 'separation_final',
            'contact_slide_distance', 'glue_contact_slide_distance', 'contact_displacements',
            'superelement_tables',
            'cstm', 'trmbu', 'trmbd',
        ]
        base_objs_map = {}
        for base_name in base_names:
            obj = getattr(self, base_name)
            if obj:
                base_objs_map[base_name] = obj
        return base_objs_map

    def get_table_types(self) -> list[str]:
        """combines all the table_types from all objects and sub-objects"""
        base = [
            'eqexin', 'gpdt', 'bgpdt', 'psds', 'monitor1', 'monitor3',
            'separation_initial', 'separation_final',
            'contact_slide_distance', 'glue_contact_slide_distance', 'contact_displacements',
            'bolt_results',
            'superelement_tables',
            'cstm', 'trmbu', 'trmbd',
        ]
        sum_objs = self._get_sum_objects()
        for objs in sum_objs:
            base.extend(objs.get_table_types())
        return base

    def __repr__(self):
        msg = 'Results:\n'

        # all these objects have data
        base_obj_map = self._get_base_objects_map()
        sum_obj_map = self._get_sum_objects_map()

        for key, obj in base_obj_map.items():
            msg += f'  {key}\n'

        for key, obj in sum_obj_map.items():
            sub_results = obj.get_table_types()

            msgi = ''
            for sub_result in sub_results:
                unused_base, sub_result2 = sub_result.split('.')
                res = getattr(obj, sub_result2)
                if res is None or res == {}:
                    continue
                msgi += f'    {sub_result2}\n'
                #msg += f'  {key}\n'
            if msgi:
                msg += f'  {key}:\n'
                msg += msgi
        return msg


class Load:
    def __len__(self) -> int:
        names = self.get_table_types(include_class=False)
        length = sum((len(getattr(self, name)) for name in names))
        return length

    def get_table_types(self):
        raise NotImplementedError('get_table_types')


class SolutionSet(Load):
    def __init__(self):
        super().__init__()
        self.displacements = {}
        self.velocities = {}
        self.accelerations = {}
        self.eigenvectors = {}

    def get_table_types(self, include_class: bool=True) -> list[str]:
        tables = [
            'displacements', 'velocities', 'accelerations', 'eigenvectors',
        ]
        if include_class:
            return ['solution_set.' + table for table in tables]
        return tables

class Acoustic(Load):
    def __init__(self):
        super().__init__()
        self.displacements = {}

    def get_table_types(self, include_class: bool=True) -> list[str]:
        tables = [
            'displacements',
        ]
        if include_class:
            return ['acoustic.' + table for table in tables]
        return tables

class ModalContribution:
    def __init__(self):
        self.displacements = {}

        self.celas1_stress = {}
        self.celas2_stress = {}
        self.celas3_stress = {}
        self.celas4_stress = {}

        self.celas1_strain = {}
        self.celas2_strain = {}
        self.celas3_strain = {}
        self.celas4_strain = {}

        self.crod_stress = {}
        self.conrod_stress = {}
        self.ctube_stress = {}
        self.crod_strain = {}
        self.conrod_strain = {}
        self.ctube_strain = {}

        self.cbend_stress = {}

        self.ctetra_stress = {}
        self.cpenta_stress = {}
        self.chexa_stress = {}

        self.ctetra_strain = {}
        self.cpenta_strain = {}
        self.chexa_strain = {}

        self.cbar_stress = {}
        self.cbar_strain = {}
        self.cbeam_stress = {}
        self.cbeam_strain = {}

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

        self.cshear_stress = {}
        self.cshear_strain = {}
        self.cshear_force = {}
        self.cbush_stress = {}
        self.cbush_strain = {}

    def get_table_types(self, include_class: bool=True) -> list[str]:
        tables = [
            'displacements', # 'velocities', 'accelerations',
            #'load_vectors', 'spc_forces', 'mpc_forces',

            #'celas1_force', 'celas2_force', 'celas3_force', 'celas4_force',
            #'crod_force', 'conrod_force', 'ctube_force',
            #'cbar_force', 'cbeam_force',
            #'cquad4_force', 'cquad8_force', 'cquadr_force',
            #'ctria3_force', 'ctria6_force', 'ctriar_force',

            'celas1_stress', 'celas2_stress', 'celas3_stress', 'celas4_stress',
            'crod_stress', 'conrod_stress', 'ctube_stress',
            'cbar_stress', 'cbeam_stress',
            'ctria3_stress', 'ctriar_stress', 'ctria6_stress',
            'cquadr_stress', 'cquad4_stress', 'cquad8_stress',
            'ctetra_stress', 'cpenta_stress', 'chexa_stress',

            'celas1_strain', 'celas2_strain', 'celas3_strain', 'celas4_strain',
            'crod_strain', 'conrod_strain', 'ctube_strain',
            'cbar_strain', 'cbeam_strain',
            'ctria3_strain', 'ctriar_strain', 'ctria6_strain',
            'cquadr_strain', 'cquad4_strain', 'cquad8_strain',
            'ctetra_strain', 'cpenta_strain', 'chexa_strain',

            'cbend_stress', # 'cbend_strain', 'cbend_force',
            'cbush_stress', 'cbush_strain',
            'cshear_stress', 'cshear_strain', 'cshear_force',

            'cquad4_composite_stress', 'cquad8_composite_stress', 'cquadr_composite_stress',
            'ctria3_composite_stress', 'ctria6_composite_stress', 'ctriar_composite_stress',

            'cquad4_composite_strain', 'cquad8_composite_strain', 'cquadr_composite_strain',
            'ctria3_composite_strain', 'ctria6_composite_strain', 'ctriar_composite_strain',

            #'cbush_force',
            #'cdamp1_force', 'cdamp2_force', 'cdamp3_force', 'cdamp4_force',
            #'cvisc_force',
        ]
        if include_class:
            return ['modal_contribution.' + table for table in tables]
        return tables

class StrengthRatio:
    def __init__(self):
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

    def get_table_types(self, include_class: bool=True) -> list[str]:
        tables = [
            'cquad4_composite_stress', 'cquad8_composite_stress', 'cquadr_composite_stress',
            'ctria3_composite_stress', 'ctria6_composite_stress', 'ctriar_composite_stress',

            'cquad4_composite_strain', 'cquad8_composite_strain', 'cquadr_composite_strain',
            'ctria3_composite_strain', 'ctria6_composite_strain', 'ctriar_composite_strain',
        ]
        if include_class:
            return ['strength_ratio.' + table for table in tables]
        return tables

class FailureIndices:
    def __init__(self):
        self.cquad4_composite_force = {}
        self.cquad8_composite_force = {}
        self.cquadr_composite_force = {}
        self.ctria3_composite_force = {}
        self.ctria6_composite_force = {}
        self.ctriar_composite_force = {}

    def get_table_types(self, include_class: bool=True) -> list[str]:
        tables = [
            'cquad4_composite_force',
            'cquad8_composite_force',
            'cquadr_composite_force',
            'ctria3_composite_force',
            'ctria6_composite_force',
            'ctriar_composite_force',
        ]
        if include_class:
            return ['failure_indices.' + table for table in tables]
        return tables

class Force(Load):
    def __init__(self):
        super().__init__()
        self.celas1_force = {}
        self.celas2_force = {}
        self.celas3_force = {}
        self.celas4_force = {}

        self.cdamp1_force = {}
        self.cdamp2_force = {}
        self.cdamp3_force = {}
        self.cdamp4_force = {}

        self.crod_force = {}
        self.conrod_force = {}
        self.ctube_force = {}

        self.cbeam_force = {}
        self.cbar_force = {}
        self.cbar_force_10nodes = {}
        self.cbend_force = {}

        self.ctria3_force = {}
        self.ctria6_force = {}
        self.ctriar_force = {}
        self.cquad4_force = {}
        self.cquad8_force = {}
        self.cquadr_force = {}

        self.cvisc_force = {}
        self.cgap_force = {}
        self.cbear_force = {}
        self.cbush_force = {}
        self.cfast_force = {}
        self.cweld_force = {}
        self.cvisc_force = {}

        self.cshear_force = {}
        self.cconeax_force = {}

        # solidPressureForces
        self.chexa_pressure_force = {}
        self.cpenta_pressure_force = {}
        self.ctetra_pressure_force = {}
        self.cpyram_pressure_force = {}

    def get_table_types(self, include_class: bool=True) -> list[str]:
        tables = [
            # 0d
            'celas1_force', 'celas2_force', 'celas3_force', 'celas4_force',
            'cdamp1_force', 'cdamp2_force', 'cdamp3_force', 'cdamp4_force',
            'cvisc_force', 'cgap_force', 'cbush_force', 'cconeax_force',

            # 1d
            'crod_force', 'conrod_force', 'ctube_force',
            'cbar_force', 'cbeam_force', 'cbend_force', 'cbar_force_10nodes',
            'cfast_force', 'cweld_force', 'cbear_force',

            # 2d
            'ctria3_force', 'ctria6_force', 'ctriar_force',
            'cquad4_force', 'cquad8_force', 'cquadr_force',
            'cshear_force',

            # solid pressure forces
            'chexa_pressure_force', 'cpenta_pressure_force',
            'ctetra_pressure_force', 'cpyram_pressure_force',
        ]
        if include_class:
            return ['force.' + table for table in tables]
        return tables


class ThermalLoad(Load):
    def __init__(self):
        super().__init__()
        #OEF - Fluxes - tCode=4 thermal=1
        self.conv_thermal_load = {}

        #self.thermalLoad_CHBDY = {}
        self.chbdye_thermal_load = {}
        self.chbdyg_thermal_load = {}
        self.chbdyp_thermal_load = {}
        self.chbdye_thermal_load_flux = {}
        self.chbdyg_thermal_load_flux = {}
        self.chbdyp_thermal_load_flux = {}

        #self.thermalLoad_1D
        self.crod_thermal_load = {}
        self.cbeam_thermal_load = {}
        self.ctube_thermal_load = {}
        self.conrod_thermal_load = {}
        self.cbar_thermal_load = {}
        self.cbend_thermal_load = {}
        self.crod_thermal_load_flux = {}
        self.cbeam_thermal_load_flux = {}
        self.ctube_thermal_load_flux = {}
        self.conrod_thermal_load_flux = {}
        self.cbar_thermal_load_flux = {}
        self.cbend_thermal_load_flux = {}

        #self.thermalLoad_2D_3D
        self.cquad4_thermal_load = {}
        self.ctriax6_thermal_load = {}
        self.cquad8_thermal_load = {}
        self.ctria3_thermal_load = {}
        self.ctria6_thermal_load = {}
        self.ctetra_thermal_load = {}
        self.chexa_thermal_load = {}
        self.cpenta_thermal_load = {}

        self.cquad4_thermal_load_flux = {}
        self.ctriax6_thermal_load_flux = {}
        self.cquad8_thermal_load_flux = {}
        self.ctria3_thermal_load_flux = {}
        self.ctria6_thermal_load_flux = {}
        self.ctetra_thermal_load_flux = {}
        self.chexa_thermal_load_flux = {}
        self.cpenta_thermal_load_flux = {}
        #self.temperatureForces = {}

    def get_table_types(self, include_class: bool=True) -> list[str]:
        tables = [
            'conv_thermal_load',
            # flux
            'chbdye_thermal_load',
            'chbdyg_thermal_load',
            'chbdyp_thermal_load',
            'chbdye_thermal_load_flux',
            'chbdyg_thermal_load_flux',
            'chbdyp_thermal_load_flux',

            # 1D
            'crod_thermal_load',
            'cbeam_thermal_load',
            'ctube_thermal_load',
            'conrod_thermal_load',
            'cbar_thermal_load',
            'cbend_thermal_load',
            'crod_thermal_load_flux',
            'cbeam_thermal_load_flux',
            'ctube_thermal_load_flux',
            'conrod_thermal_load_flux',
            'cbar_thermal_load_flux',
            'cbend_thermal_load_flux',

            #self.thermalLoad_2D_3D
            'cquad4_thermal_load',
            'ctriax6_thermal_load',
            'cquad8_thermal_load',
            'ctria3_thermal_load',
            'ctria6_thermal_load',
            'ctetra_thermal_load',
            'chexa_thermal_load',
            'cpenta_thermal_load',

            # 2d/3d
            'cquad4_thermal_load_flux',
            'ctriax6_thermal_load_flux',
            'cquad8_thermal_load_flux',
            'ctria3_thermal_load_flux',
            'ctria6_thermal_load_flux',
            'ctetra_thermal_load_flux',
            'chexa_thermal_load_flux',
            'cpenta_thermal_load_flux',
        ]
        if include_class:
            return ['thermal_load.' + table for table in tables]
        return tables

class Stress:
    def __init__(self, word: str):
        """word : str
            'stress', 'stressa'
        """
        self.word = word
        self.celas1_stress = {}
        self.celas2_stress = {}
        self.celas3_stress = {}
        self.celas4_stress = {}

        # rods
        self.crod_stress = {}
        self.conrod_stress = {}
        self.ctube_stress = {}

        # bars/beams
        self.cbar_stress = {}
        self.cbar_stress_10nodes = {}
        self.cbeam_stress = {}
        self.cbend_stress = {}


        # other 1d/2d
        self.cshear_stress = {}
        self.cweld_stress = {}
        self.cfast_stress = {}
        self.cbush_stress = {}

        # shells
        self.ctria3_stress = {}
        self.ctria6_stress = {}
        self.cquad4_stress = {}
        self.cquad8_stress = {}
        self.cquadr_stress = {}
        self.ctriar_stress = {}

        self.cquad4_composite_stress = {}
        self.cquad8_composite_stress = {}
        self.cquadr_composite_stress = {}
        self.ctria3_composite_stress = {}
        self.ctria6_composite_stress = {}
        self.ctriar_composite_stress = {}

        self.cplstn3_stress = {}
        self.cplstn4_stress = {}
        self.cplstn6_stress = {}
        self.cplstn8_stress = {}
        self.cplsts3_stress = {}
        self.cplsts4_stress = {}
        self.cplsts6_stress = {}
        self.cplsts8_stress = {}

        # solids
        self.ctetra_stress = {}
        self.cpenta_stress = {}
        self.chexa_stress = {}
        self.cpyram_stress = {}

        # 269, 270
        self.chexa_composite_stress = {}
        self.cpenta_composite_stress = {}

        #: OES - CTRIAX6
        self.ctriax_stress = {}

        # bushing
        self.cbush1d_stress_strain = {}
        self.hyperelastic_cquad4_stress = {}

    def get_table_types(self, include_class: bool=True) -> list[str]:
        tables = [
            # OES - CELAS1/CELAS2/CELAS3/CELAS4 stress
            'celas1_stress', 'celas2_stress', 'celas3_stress', 'celas4_stress',

            # OES - isotropic CTETRA/CHEXA/CPENTA stress
            'ctetra_stress', 'cpenta_stress', 'chexa_stress', 'cpyram_stress',

            'chexa_composite_stress', 'cpenta_composite_stress',

            # OES - isotropic CROD/CONROD/CTUBE stress/strain
            'crod_stress', 'conrod_stress', 'ctube_stress',

            # OES - isotropic CBAR stress/strain
            'cbar_stress',
            'cbar_stress_10nodes',

            # OES - isotropic CBEAM stress/strain
            'cbeam_stress',

            # CBEND - isotropic CBEAM stress/strain
            'cbend_stress', # 'cbend_force',

            # OES - isotropic CTRIA3/CQUAD4 stress/strain
            'ctria3_stress', 'ctriar_stress', 'ctria6_stress',
            'cquadr_stress', 'cquad4_stress', 'cquad8_stress',
            'ctriax_stress',
            # OES - CTRIAX6
            'cbush_stress',

            # OES - composite CTRIA3/CQUAD4 stress/strain
            'cquad4_composite_stress', 'cquad8_composite_stress', 'cquadr_composite_stress',
            'ctria3_composite_stress', 'ctria6_composite_stress', 'ctriar_composite_stress',

            # OES - CSHEAR stress/strain
            'cshear_stress',

            'cplstn3_stress', 'cplstn4_stress', 'cplstn6_stress', 'cplstn8_stress',
            'cplsts3_stress', 'cplsts4_stress', 'cplsts6_stress', 'cplsts8_stress',

            'cweld_stress',
            'cfast_stress',
            'cbush1d_stress_strain',
            'hyperelastic_cquad4_stress',
        ]
        if include_class:
            return [f'{self.word}.' + table for table in tables]
        return tables

class Strain:
    def __init__(self):
        # springs
        self.celas1_strain = {}
        self.celas2_strain = {}
        self.celas3_strain = {}
        self.celas4_strain = {}

        # rods
        self.crod_strain = {}
        self.conrod_strain = {}
        self.ctube_strain = {}

        # bars/beams
        self.cbar_strain = {}
        self.cbar_strain_10nodes = {}
        self.cbeam_strain = {}
        self.cbend_strain = {}

        # other 1d
        self.cbush_strain = {}
        self.cfast_strain = {}
        self.cweld_strain = {}
        self.cshear_strain = {}

        # shells
        self.ctria3_strain = {}
        self.ctria6_strain = {}
        self.cquad4_strain = {}
        self.cquad8_strain = {}
        self.cquadr_strain = {}
        self.ctriar_strain = {}

        self.cquad4_composite_strain = {}
        self.cquad8_composite_strain = {}
        self.cquadr_composite_strain = {}
        self.ctria3_composite_strain = {}
        self.ctria6_composite_strain = {}
        self.ctriar_composite_strain = {}

        self.cplstn3_strain = {}
        self.cplstn4_strain = {}
        self.cplstn6_strain = {}
        self.cplstn8_strain = {}
        self.cplsts3_strain = {}
        self.cplsts4_strain = {}
        self.cplsts6_strain = {}
        self.cplsts8_strain = {}

        # solids
        self.ctetra_strain = {}
        self.cpenta_strain = {}
        self.chexa_strain = {}
        self.cpyram_strain = {}

        # 269, 270
        self.chexa_composite_strain = {}
        self.cpenta_composite_strain = {}

        #: OES - CTRIAX6
        self.ctriax_strain = {}
        self.ctriax6_strain = {}

        self.hyperelastic_cquad4_strain = {}


    def get_table_types(self, include_class: bool=True) -> list[str]:
        tables = [
            # OES - CELAS1/CELAS2/CELAS3/CELAS4 strain
            'celas1_strain', 'celas2_strain', 'celas3_strain', 'celas4_strain',

            # OES - isotropic CTETRA/CHEXA/CPENTA strain
            'ctetra_strain', 'cpenta_strain', 'chexa_strain', 'cpyram_strain',

            'chexa_composite_strain', 'cpenta_composite_strain',

            # OES - isotropic CROD/CONROD/CTUBE
            'crod_strain', 'conrod_strain', 'ctube_strain',

            # OES - isotropic CBAR
            'cbar_strain',          # OES - isotropic CBAR
            'cbar_strain_10nodes',  # OES - isotropic CBAR
            'cbeam_strain',         # OES - isotropic CBEAM
            'cbend_strain',         # OES - isotropic CBEND

            # OES - isotropic CTRIA3/CQUAD4
            'ctria3_strain', 'ctriar_strain', 'ctria6_strain',
            'cquadr_strain', 'cquad4_strain', 'cquad8_strain',
            'ctriax_strain', 'ctriax6_strain',

            # OES - composite CTRIA3/CQUAD4
            'cquad4_composite_strain', 'cquad8_composite_strain', 'cquadr_composite_strain',
            'ctria3_composite_strain', 'ctria6_composite_strain', 'ctriar_composite_strain',

            # OES - CSHEAR stress/strain
            'cplstn3_strain', 'cplstn4_strain', 'cplstn6_strain', 'cplstn8_strain',
            'cplsts3_strain', 'cplsts4_strain', 'cplsts6_strain', 'cplsts8_strain',

            # other 1d/2d
            'cshear_strain',
            'cweld_strain',
            'cfast_strain',
            'cbush_strain',
            'hyperelastic_cquad4_strain',


        ]
        if include_class:
            return ['strain.' + table for table in tables]
        return tables

class StrainEnergy:
    def __init__(self):
        """
        OEE - strain energy density; tCode=18
        """
        self.celas1_strain_energy = {}
        self.celas2_strain_energy = {}
        self.celas3_strain_energy = {}
        self.celas4_strain_energy = {}

        self.cdamp1_strain_energy = {}
        self.cdamp2_strain_energy = {}
        self.cdamp3_strain_energy = {}
        self.cdamp4_strain_energy = {}

        self.cquad4_strain_energy = {}
        self.cquad8_strain_energy = {}
        self.cquadr_strain_energy = {}
        self.cquadx_strain_energy = {}

        self.ctria3_strain_energy = {}
        self.ctria6_strain_energy = {}
        self.ctriar_strain_energy = {}
        self.ctriax_strain_energy = {}
        self.ctriax6_strain_energy = {}

        self.ctetra_strain_energy = {}
        self.cpenta_strain_energy = {}
        self.chexa_strain_energy = {}
        self.cpyram_strain_energy = {}

        self.crod_strain_energy = {}
        self.ctube_strain_energy = {}
        self.conrod_strain_energy = {}

        self.cbar_strain_energy = {}
        self.cbeam_strain_energy = {}
        self.cbend_strain_energy = {}
        self.cbeam3_strain_energy = {}

        self.cgap_strain_energy = {}
        self.cdum8_strain_energy = {}
        self.cbush_strain_energy = {}
        #self.chexa8fd_strain_energy = {}
        self.dmig_strain_energy = {}
        self.genel_strain_energy = {}
        self.cshear_strain_energy = {}
        self.conm2_strain_energy = {}
        self.rbe1_strain_energy = {}
        self.rbe3_strain_energy = {}
        self.cweld_strain_energy = {}
        self.cfast_strain_energy = {}
        self.cseam_strain_energy = {}

    def get_table_types(self, include_class: bool=True) -> list[str]:
        tables = [
            # OEE - strain energy density # tCode=18
            'cquad4_strain_energy', 'cquad8_strain_energy', 'cquadr_strain_energy',
            'cquadx_strain_energy',

            'ctria3_strain_energy', 'ctria6_strain_energy', 'ctriar_strain_energy',
            'ctriax_strain_energy', 'ctriax6_strain_energy',

            'cshear_strain_energy',

            'ctetra_strain_energy', 'cpenta_strain_energy',
            'chexa_strain_energy', 'cpyram_strain_energy',

            'crod_strain_energy', 'ctube_strain_energy', 'conrod_strain_energy',

            'cbar_strain_energy', 'cbeam_strain_energy', 'cbeam3_strain_energy',

            'cgap_strain_energy',
            'cbush_strain_energy',
            'celas1_strain_energy', 'celas2_strain_energy',
            'celas3_strain_energy', 'celas4_strain_energy',

            'cdamp1_strain_energy', 'cdamp2_strain_energy',
            'cdamp3_strain_energy', 'cdamp4_strain_energy',

            'cdum8_strain_energy',
            #'chexa8fd_strain_energy'
            'cbend_strain_energy',
            'dmig_strain_energy',
            'genel_strain_energy',
            'conm2_strain_energy',
            'rbe1_strain_energy', 'rbe3_strain_energy',
            'cweld_strain_energy', 'cfast_strain_energy', 'cseam_strain_energy',
        ]
        if include_class:
            return ['strain_energy.' + table for table in tables]
        return tables


class CSTM:
    def __init__(self):
        self.headers = {"cid": 0,
                        "cid_type": 1,
                        "unused_int_index": 2,
                        "unused_double_index": 3,
                        "ox": 4,
                        "oy": 5,
                        "oz": 6,
                        "T11": 7,
                        "T12": 8,
                        "T13": 9,
                        "T21": 10,
                        "T22": 11,
                        "T23": 12,
                        "T31": 13,
                        "T32": 14,
                        "T33": 15
                    }
        self.data = None  # type: Optional[np.ndarray]  # Coordinate Transformation Matrices from Native to Global

    def get_stats(self, short=False):
        return str(self) + '\n'

    def __repr__(self) -> str:
        msg = 'CSTM:\n'
        msg += f'  headers_str = {self.headers.keys()}\n'
        msg += f'  headers_ints = {self.headers.values()}\n'
        if self.data is not None:
            msg += f'  data.shape = {self.data.shape}'
        else:
            msg += '  data = None'
        return msg

class TRMBD:
    def __init__(self, **data: dict[str, Any]):
        self.isubcase = data['isubcase']
        self.analysis_code = data['analysis_code']
        self.device_code = data['device_code']
        self.sort_code = data['sort_code']
        self.table_name = data['table_name']
        self.table_code = data['table_code']
        self.title = data['title']
        self.subtitle = data['subtitle']
        self.label = data['label']

        self.times = np.array([], dtype='float64') # default type
        self.nodes: dict[str, np.ndarray] = {}
        self.eulersx: dict[str, np.ndarray] = {}
        self.eulersy: dict[str, np.ndarray] = {}
        self.eulersz: dict[str, np.ndarray] = {}

    def etypes(self) -> list[str]:
        etypes = list(self.nodes.keys())
        return etypes

    def get_stats(self, short=False):
        etypes = self.etypes()
        if short:
            return [f'op2_results.trmbd[{self.isubcase}]: TRMBD(time, nodes, eulersx, eulersy, eulersz)\n']
        else:
            return [
                f'op2_results.trmbd[{self.isubcase}]: TRMBD\n'
                f'  time = {self.time}\n'
                f'  etypes = {etypes}\n'
                f'  nodes = {self.nodes}\n'
                f'  eulersx, eulersy, eulersz\n']

    def __eq__(self, trmbd) -> bool:
        return True

    def __repr__(self) -> str:
        msg = 'op2_results.trmbd:\n'
        msg += f'  isubcase = {self.isubcase}\n'
        msg += f'  nodes, eulersx, eulersy, eulersz'
        return msg

class TRMBU:
    def __init__(self, ntimes: int, **data: dict[str, Any]):
        self.isubcase = data['isubcase']
        self.analysis_code = data['analysis_code']
        self.device_code = data['device_code']
        self.sort_code = data['sort_code']
        self.table_name = data['table_name']
        self.table_code = data['table_code']
        self.title = data['title']
        self.subtitle = data['subtitle']
        self.label = data['label']

        self.ntimes = ntimes
        self.times = np.array([], dtype='float64') # default dtype
        self.eulers: dict[str, np.ndarray] = {}

    def etypes(self) -> list[int]:
        etypes = list(self.eulers.keys())
        return etypes

    def get_stats(self, short=False):
        etypes = self.etypes()
        if short:
            return [f'op2_results.trmbu[{self.isubcase}]: TRMBU(time, eulers); etypes={etypes}\n']
        else:
            return [
                f'op2_results.trmbu[{self.isubcase}]: TRMBU\n'
                f'  time = {self.time}\n'
                f'  etypes = {etypes}\n'
                f'  eulers\n']

    def __eq__(self, trmbu) -> bool:
        return True

    def __repr__(self) -> str:
        msg = 'op2_results.trmbu:\n'
        msg += f'  isubcase = {self.isubcase}\n'
        msg += f'  eulers'
        return msg
