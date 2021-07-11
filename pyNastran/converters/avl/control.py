import numpy as np

class Control:
    def __init__(self, name: str, gain: float, xhinge: float,
                 hinge_vector: np.ndarray, sign_deflection_duplicated: float):
        """
        Parameters
        ----------
        name : str
            name of control variable
        gain : float
            control deflection gain, units:  degrees deflection / control variable
        xhinge : float
            x/c location of hinge.
            If positive, control surface extent is Xhinge..1  (TE surface)
            If negative, control surface extent is 0..-Xhinge (LE surface)
        XYZhvec : (3,) float ndarray
            vector giving hinge axis about which surface rotates
            + deflection is + rotation about hinge vector by righthand rule
            Specifying XYZhvec = 0. 0. 0. puts the hinge vector along the hinge
        SgnDup : float
            sign of deflection for duplicated surface
            An elevator would have SgnDup = +1.0
            An aileron  would have SgnDup = -1.0

        """
        self.name = name
        self.gain = gain
        self.xhinge = xhinge
        self.hinge_vector = hinge_vector
        self.sign_deflection_duplicated = sign_deflection_duplicated
        assert isinstance(xhinge, float), xhinge
        assert isinstance(gain, float), gain
        assert isinstance(sign_deflection_duplicated, float), sign_deflection_duplicated

    def write(self) -> str:
        x, y, z = self.hinge_vector
        msg = (
            'CONTROL\n'
            '! name, gain, xhinge, [hinge_vector], sign_deflection_duplicated\n'
            f'{self.name} {self.gain} {self.xhinge} {x} {y} {z} {self.sign_deflection_duplicated}\n'
        )
        return msg

    def __repr__(self) -> str:
        msg = (f'Control(name={self.name!r}, gain={self.gain}, xhinge={self.xhinge}, '
               f'hinge_vector={self.hinge_vector}, sign_deflection_duplicated={self.sign_deflection_duplicated})')
        return msg

