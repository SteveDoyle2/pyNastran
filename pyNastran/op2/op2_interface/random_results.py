class RandomObjects(object):
    prefix = ''
    postfix = ''
    def __init__(self):
        self.displacements = {}
        self.velocities = {}
        self.accelerations = {}
        self.load_vectors = {}
        self.spc_forces = {}
        self.mpc_forces = {}

        self.crod_force = {}
        self.conrod_force = {}
        self.ctube_force = {}

        self.cbar_force = {}
        self.cbeam_force = {}

        self.cbush_stress = {}
        self.cbush_strain = {}

        self.crod_stress = {}
        self.conrod_stress = {}
        self.ctube_stress = {}
        self.cbar_stress = {}
        self.cbeam_stress = {}

        self.crod_strain = {}
        self.conrod_strain = {}
        self.ctube_strain = {}
        self.cbar_strain = {}
        self.cbeam_strain = {}

        self.ctetra_strain = {}
        self.cpenta_strain = {}
        self.chexa_strain = {}

        self.ctetra_stress = {}
        self.cpenta_stress = {}
        self.chexa_stress = {}

        self.celas1_stress = {}
        self.celas2_stress = {}
        self.celas3_stress = {}
        self.celas4_stress = {}

        self.celas1_strain = {}
        self.celas2_strain = {}
        self.celas3_strain = {}
        self.celas4_strain = {}

        self.celas1_force = {}
        self.celas2_force = {}
        self.celas3_force = {}
        self.celas4_force = {}

        self.ctria3_force = {}
        self.ctria6_force = {}
        self.ctriar_force = {}
        self.cquad4_force = {}
        self.cquad8_force = {}
        self.cquadr_force = {}

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

        self.cbend_stress = {}
        self.cbend_strain = {}
        self.cbend_force = {}

        self.cshear_stress = {}
        self.cshear_strain = {}
        self.cshear_force = {}

        self.cbush_force = {}
        self.cdamp1_force = {}
        self.cdamp2_force = {}
        self.cdamp3_force = {}
        self.cdamp4_force = {}
        self.cvisc_force = {}

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
            'displacements', 'velocities', 'accelerations',
            'load_vectors', 'spc_forces', 'mpc_forces',

            'celas1_force', 'celas2_force', 'celas3_force', 'celas4_force',
            'crod_force', 'conrod_force', 'ctube_force',
            'cbar_force', 'cbeam_force',
            'cquad4_force', 'cquad8_force', 'cquadr_force',
            'ctria3_force', 'ctria6_force', 'ctriar_force',

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

            'cquad4_composite_stress', 'cquad8_composite_stress', 'cquadr_composite_stress',
            'ctria3_composite_stress', 'ctria6_composite_stress', 'ctriar_composite_stress',

            'cquad4_composite_strain', 'cquad8_composite_strain', 'cquadr_composite_strain',
            'ctria3_composite_strain', 'ctria6_composite_strain', 'ctriar_composite_strain',

            'cbend_stress', 'cbend_strain', 'cbend_force',
            'cbush_stress', 'cbush_strain',
            'cshear_stress', 'cshear_strain', 'cshear_force',

            'cbush_force',
            'cdamp1_force', 'cdamp2_force', 'cdamp3_force', 'cdamp4_force',
            'cvisc_force',

        ]
        return [self.prefix + table + self.postfix for table in tables]

class AutoCorrelationObjects(RandomObjects):
    """storage class for the ATO objects"""
    prefix = 'ato.'
    #postfix = ''

class PowerSpectralDensityObjects(RandomObjects):
    """storage class for the PSD objects"""
    prefix = 'psd.'
    #postfix = ''

class RootMeansSquareObjects(RandomObjects):
    """storage class for the RMS objects"""
    prefix = 'rms.'
    #postfix = ''

class CumulativeRootMeansSquareObjects(RandomObjects):
    """storage class for the CRMS objects"""
    prefix = 'crm.'
    #postfix = ''

class NumberOfCrossingsObjects(RandomObjects):
    """storage class for the NO objects"""
    prefix = 'no.'
    #postfix = ''

class RAECONS(object):
    """storage class for the RAECONS objects"""
    def __init__(self):
        self.ctria3_strain = {}
        self.cquad4_strain = {}
        self.chexa_strain = {}

    def get_table_types(self):
        tables = [
            'chexa_strain',
            'ctria3_strain', 'cquad4_strain',
        ]
        return ['RAECONS.' + table for table in tables]

class RASCONS(object):
    """storage class for the RASCONS objects"""
    def __init__(self):
        self.ctetra_stress = {}
        self.cpenta_stress = {}
        self.chexa_stress = {}

        self.ctetra_strain = {}
        self.cpenta_strain = {}
        self.chexa_strain = {}

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

    def get_table_types(self):
        tables = [
            # OES - isotropic CTRIA3/CQUAD4 stress
            'ctria3_stress', 'ctriar_stress', 'ctria6_stress',
            'cquadr_stress', 'cquad4_stress', 'cquad8_stress',

            # OES - isotropic CTRIA3/CQUAD4 strain
            'ctria3_strain', 'ctriar_strain', 'ctria6_strain',
            'cquadr_strain', 'cquad4_strain', 'cquad8_strain',

            'ctetra_stress', 'chexa_stress', 'cpenta_stress',
            'ctetra_strain', 'chexa_strain', 'cpenta_strain',
        ]
        return ['RASCONS.' + table for table in tables]

class RAPCONS(object):
    """storage class for the RAPCONS objects"""
    def __init__(self):
        self.cquad4_composite_stress = {}
        self.cquad8_composite_stress = {}
        self.cquadr_composite_stress = {}
        self.ctria3_composite_stress = {}
        self.ctria6_composite_stress = {}
        self.ctriar_composite_stress = {}

    def get_table_types(self):
        tables = [
            'cquad4_composite_stress',
            'cquad8_composite_stress',
            'cquadr_composite_stress',
            'ctria3_composite_stress',
            'ctria6_composite_stress',
            'ctriar_composite_stress',
            #'cquad4_composite_strain',
            #'cquad8_composite_strain',
            #'cquadr_composite_strain',
            #'ctria3_composite_strain',
            #'ctria6_composite_strain',
            #'ctriar_composite_strain',
        ]
        return ['RAPCONS.' + table for table in tables]

class RAPEATC(object):
    """storage class for the RAPEATC objects"""
    def __init__(self):
        self.cquad4_composite_stress = {}
        self.cquad8_composite_stress = {}
        self.cquadr_composite_stress = {}
        self.ctria3_composite_stress = {}
        self.ctria6_composite_stress = {}
        self.ctriar_composite_stress = {}

    def get_table_types(self):
        tables = [
            'cquad4_composite_stress',
            'cquad8_composite_stress',
            'cquadr_composite_stress',
            'ctria3_composite_stress',
            'ctria6_composite_stress',
            'ctriar_composite_stress',

            #'cquad4_composite_strain',
            #'cquad8_composite_strain',
            #'cquadr_composite_strain',
            #'ctria3_composite_strain',
            #'ctria6_composite_strain',
            #'ctriar_composite_strain',
        ]
        return ['RAPEATC.' + table for table in tables]

class RAFCONS(object):
    """storage class for the RAFCONS objects"""
    def __init__(self):
        self.cbar_force = {}
        self.cquad4_force = {}
        self.cbush_force = {}

    def get_table_types(self):
        tables = [
            'cbar_force',
            'cquad4_force',
            'cbush_force',
        ]
        return ['RAFCONS.' + table for table in tables]

class RAGCONS(object):
    """storage class for the RAGCONS objects"""
    def __init__(self):
        self.grid_point_forces = {}
    def get_table_types(self):
        tables = [
            'grid_point_forces',
        ]
        return ['RAGCONS.' + table for table in tables]

class RAGEATC(object):
    """storage class for the RAGEATC objects"""
    def __init__(self):
        self.grid_point_forces = {}

    def get_table_types(self):
        tables = [
            'grid_point_forces',
        ]
        return ['RAGEATC.' + table for table in tables]


class RANCONS(object):
    """storage class for the RANCONS objects"""
    def __init__(self):
        self.cbar_strain_energy = {}
        self.cbush_strain_energy = {}
        self.chexa_strain_energy = {}
        self.ctria3_strain_energy = {}
        self.cquad4_strain_energy = {}

    def get_table_types(self):
        tables = [
            'cbar_strain_energy', 'cbush_strain_energy',
            'chexa_strain_energy',
            'ctria3_strain_energy', 'cquad4_strain_energy',
        ]
        return ['RANCONS.' + table for table in tables]

class RADEFFM(object):
    """storage class for the RADEFFM objects"""
    def __init__(self):
        self.eigenvectors = {}
    def get_table_types(self):
        tables = [
            'eigenvectors',
        ]
        return ['RADEFFM.' + table for table in tables]


class RADCONS(object):
    def __init__(self):
        self.eigenvectors = {}

    def get_table_types(self):
        tables = [
            'eigenvectors',
        ]
        return ['RADCONS.' + table for table in tables]


class RADEATC(object):
    """storage class for the RADEATC objects"""
    def __init__(self):
        self.eigenvectors = {}

    def get_table_types(self):
        tables = [
            'eigenvectors',
        ]
        return ['RADEATC.' + table for table in tables]


class RANEATC(object):
    """storage class for the RANEATC objects"""
    def __init__(self):
        self.cbar_strain_energy = {}
        self.cbush_strain_energy = {}
        self.chexa_strain_energy = {}
        self.ctria3_strain_energy = {}
        self.cquad4_strain_energy = {}

    def get_table_types(self):
        tables = [
            'cbar_strain_energy', 'cbush_strain_energy',
            'chexa_strain_energy',
            'ctria3_strain_energy', 'cquad4_strain_energy',
        ]
        return ['RANEATC.' + table for table in tables]


class ROUGV1(object):
    """storage class for the ROUGV1 objects"""
    def __init__(self):
        self.displacements = {}
        self.velocities = {}
        self.accelerations = {}
        self.eigenvectors = {}

    def get_table_types(self):
        tables = [
            'displacements', 'velocities', 'accelerations', 'eigenvectors',
        ]
        return ['ROUGV1.' + table for table in tables]

class RAFEATC(object):
    """storage class for the RAFEATC objects"""
    def __init__(self):
        self.cbar_force = {}
        self.cquad4_force = {}
        self.cbush_force = {}

    def get_table_types(self):
        tables = [
            'cbar_force',
            'cquad4_force',
            'cbush_force',
        ]
        return ['RAFEATC.' + table for table in tables]


class RASEATC(object):
    """storage class for the RASEATC objects"""
    def __init__(self):
        self.chexa_stress = {}
        self.cquad4_stress = {}

    def get_table_types(self):
        tables = [
            'chexa_stress',
            'cquad4_stress',
        ]
        return ['RASEATC.' + table for table in tables]

class RAEEATC(object):
    """storage class for the RAEEATC objects"""
    def __init__(self):
        self.chexa_strain = {}
        self.ctria3_strain = {}
        self.cquad4_strain = {}

    def get_table_types(self):
        tables = [
            'chexa_strain',
            'ctria3_strain', 'cquad4_strain',
        ]
        return ['RAEEATC.' + table for table in tables]
