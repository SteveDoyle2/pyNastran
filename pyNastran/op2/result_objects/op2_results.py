from __future__ import print_function, unicode_literals, absolute_import

from pyNastran.op2.op2_interface.random_results import (
    RADCONS, RAECONS, RASCONS, RAPCONS, RAFCONS, RAGCONS, RANCONS,
    RADEATC, RAEEATC, RASEATC, RAPEATC, RAFEATC, RAGEATC, RANEATC,
    ROUGV1, RADEFFM,
    AutoCorrelationObjects, PowerSpectralDensityObjects, RootMeansSquareObjects,
    CumulativeRootMeansSquareObjects, NumberOfCrossingsObjects,
)

class Results(object):
    """storage object for even more op2_results (see op2.op2_results)"""
    def __init__(self):
        self.eqexin = None
        self.gpdt = None
        self.ato = AutoCorrelationObjects()
        self.psd = PowerSpectralDensityObjects()
        self.rms = RootMeansSquareObjects()
        self.no = NumberOfCrossingsObjects()
        self.crm = CumulativeRootMeansSquareObjects()

        self.modal_contribution = ModalContribution()
        self.strength_ratio = StrengthRatio()
        self.ROUGV1 = ROUGV1()   # relative disp/vel/acc/eigenvectors

        self.RADEFFM = RADEFFM() # eigenvectors

        self.RADCONS = RADCONS() # eigenvectors
        self.RAFCONS = RAFCONS() # force
        self.RASCONS = RASCONS() # stress
        self.RAECONS = RAECONS() # strain
        self.RAGCONS = RAGCONS() # grid point forces
        self.RAPCONS = RAPCONS() # composite stress
        self.RANCONS = RANCONS() # strain energy

        self.RADEATC = RADEATC() # eigenvectors
        self.RAFEATC = RAFEATC() # force
        self.RASEATC = RASEATC() # stress
        self.RAEEATC = RAEEATC() # strain
        self.RAGEATC = RAGEATC() # grid point forces
        self.RAPEATC = RAPEATC() # composite stress
        self.RANEATC = RANEATC() # strain energy

    def get_table_types(self):
        """combines all the table_types from all objects and sub-objects"""
        sum_objs = [
            self.ato, self.psd, self.rms, self.no, self.crm,
            self.modal_contribution, self.strength_ratio,
            self.ROUGV1,
            self.RADEFFM,
            self.RADCONS, self.RAFCONS, self.RASCONS, self.RAECONS, self.RAGCONS, self.RAPCONS, self.RANCONS,
            self.RADEATC, self.RAFEATC, self.RASEATC, self.RAEEATC, self.RAGEATC, self.RAPEATC, self.RANEATC,
        ]
        base = ['eqexin', 'gpdt']
        for objs in sum_objs:
            base.extend(objs.get_table_types())
        return base

class ModalContribution(object):
    def __init__(self):
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

    def get_table_types(self):
        tables = [
            #'displacements', 'velocities', 'accelerations',
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

            #'cbend_stress', 'cbend_strain', 'cbend_force',
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
        return ['modal_contribution.' + table for table in tables]

class StrengthRatio(object):
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

    def get_table_types(self):
        tables = [
            'cquad4_composite_stress', 'cquad8_composite_stress', 'cquadr_composite_stress',
            'ctria3_composite_stress', 'ctria6_composite_stress', 'ctriar_composite_stress',

            'cquad4_composite_strain', 'cquad8_composite_strain', 'cquadr_composite_strain',
            'ctria3_composite_strain', 'ctria6_composite_strain', 'ctriar_composite_strain',
        ]
        return ['strength_ratio.' + table for table in tables]
