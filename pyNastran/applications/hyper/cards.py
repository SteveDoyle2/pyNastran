from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import double, double_or_blank, string, integer
from numpy import array, cross, sin, cos, dot, radians
from numpy.linalg import norm  # type: ignore


class FLOW(BaseCard):
    type = 'FLOW'
    def __init__(self, card, comment=''):
        self.flow_id = integer(card, 1, 'flow_id')
        self.mode = string(card, 2, 'mode')
        if 'VELO' == self.mode:
            self.Vx = double(card, 3, 'Vx')
            self.Vy = double(card, 4, 'Vy')
            self.Vz = double(card, 5, 'Vz')
        elif 'MACH' == self.mode:
            #self.Nx = double(card, 3, 'Nx')
            #self.Ny = double(card, 4, 'Ny')
            #self.Nz = double(card, 5, 'Nz')

            self.Ma = double(card, 3, 'Ma')
            self.alpha = double(card, 4, 'alpha')
            self.beta = double_or_blank(card, 5, 'beta', 0.)
            #self.phi = double_or_blank(card, 9, 'phi', 0.)
            self.phi = 0.0
        self.p = double_or_blank(card, 10, 'p', 0.)
        self.q = double_or_blank(card, 11, 'q', 0.)
        self.r = double_or_blank(card, 12, 'r', 0.)

    def get_omega(self):
        # p. 110
        return array([self.p, self.q, self.r], dtype='float32')

    def get_V(self):
        a = 1.0
        if 'VELO' == self.mode:
            V = array([self.Vx, self.Vy, self.Vz], dtype='float32')
            M = V / a
        elif 'MACH' == self.mode:
            # p. 111
            ca = cos(radians(self.alpha))
            cb = cos(radians(self.beta))
            cp = cos(radians(self.phi))

            sa = sin(radians(self.alpha))
            sb = sin(radians(self.beta))
            sp = sin(radians(self.phi))
            V = self.Ma * a * array([ca*cb,
                                     sp*sa*cb+cp*sb,
                                     cp*sa*cb-sp*sb], dtype='float32')
            #n = array([self.Nx, self.Ny, self.Nz])
            #n /= norm(n)
            #V = dot(V, n)
        return V, a

    def get_Ma(self):
        pass

    def raw_fields(self):
        if self.mode == 'VELO':
            cards = ['$FLOW', self.flow_id, self.mode, self.Vx, self.Vy, self.Vz,
                     None, None, None, self.p, self.q, self.r]
        else:
            cards = ['$FLOW', self.flow_id, self.mode, self.Nx, self.Ny, self.Nz,
                     self.Ma, self.alpha, self.beta, self.phi, self.p, self.q, self.r]
        return cards


class HYPER(BaseCard):
    type = 'HYPER'

    def raw_fields(self):
        return ['$HYPER', self.pid, self.set, self.Type, self.gamma, self._cp]

    def __init__(self, card, comment=''):
        self.pid = integer(card, 1, 'property_id')
        self.sets = []
        self.Types = []
        self.gammas = []
        self._cps = []
        #self.set = integer(card, 2, 'set_id')
        #self.Type = string(card, 3, 'Type')
        #if self.Type not in ['NEWTON','PRANDTL-MEYER', 'CP']:
        #    raise RuntimeError('Type=%r' % Type)
        #self.gamma = double_or_blank(card, 4, 'gamma', 1.4)

        i = 2
        while i < len(card):
            self.sets.append(integer(card, i, 'set_id'))
            Type = string(card, i+1, 'Type')
            self.Types.append(Type)
            #if self.Type not in ['NEWTON','PRANDTL-MEYER', 'CP']:
                #raise RuntimeError('Type=%r' % Type)
            self.gammas.append(double_or_blank(card, i+2, 'gamma', 1.4))

            _cp = None
            if Type == 'CP':
                _cp = double(card, i+3, 'Cp')
            elif Type == 'NEWTON':
                _cp = double_or_blank(card, i+3, 'Cp_nominal', 2.0)
            self._cps.append(_cp)
            i += 7

    def _Cp(self, Type, delta_radians):
        if self.Type == 'NEWTON':
            Cp = self._cp * sin(delta_radians)**2.0
        elif self.Type == 'TANWEDGE':
            gp1 = gamma + 1
            #Cp = 4.0/gp1 * (sin2_beta - 1/M1**2)
            Cp = 4.0/gp1 * sin2_beta
        return Cp

    def Cp(self, delta_radians):
        if self.Type == 'CP':
            return self._cp
        elif self.Type == 'NEWTON':
            Cp = self._cp * sin(delta_radians)**2.0
        elif self.Type == 'TANWEDGE':
            Cp = SELF._Cp('TANWEDGE', delta_radians)
        #elif self.Type == 'TAN_CONE':
            #pass
        #elif self.Type == 'INC_CONE':
            #pass
        #elif self.Type == 'VAN_DYKE':
            #pass
        #elif self.Type == 'FREEMOLC':
            #pass
        elif self.Type == 'HANKEY':
            Cp = 1.95 * sin(delta)**2 + 0.21*cos(delta)*sin(delta)
        elif self.Type == 'DALMBUCK':
            if delta <= 22.5:
                Cp = (1./sin(4*delta)**0.75 + 1.0) * sin(delta)**2
            else:
                Cp = self._Cp('TANWEDGE', delta_radians)
        return Cp

