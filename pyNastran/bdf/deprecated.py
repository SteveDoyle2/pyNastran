# pylint: disable=E1101,C0103

class DeprecatedCompositeShellProperty(object):
    """
    To be deprecated in:
      - Version 0.7
    To be removed in:
      - Version 0.8
    """
    def __init__(self):
        pass

    def MassPerArea(self, iply='all', method='nplies'):
        #self.deprecated('MassPerArea(iply, method)', 'get_mass_per_area(iply, method)', '0.8')
        return self.get_mass_per_area(iply, method)

    def Thickness(self, iply='all'):
        return self.get_thickness(iply)

    def nPlies(self):
        return self.get_nplies()

    def get_nplies(self):
        return self.nplies

    def get_material_ids(self):
        return self.material_ids

    def Nsm(self):
        return self.get_nonstructural_mass()

    def isSymmetrical(self):
        return self.is_symmetrical()

    def Rho(self, iply):
        return self.get_density(iply)

    def Theta(self, iply):
        return self.get_theta(iply)

    def sout(self, iply):
        return self.get_sout(iply)
