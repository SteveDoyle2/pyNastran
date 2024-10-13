from typing import Union, Any
from pyNastran.utils.numpy_utils import integer_types

class RealMatrixCard:
    """
    interface for cards of the form:
       M2GG=MDMIG
       M2GG=MDMIG1, MDMIG2, MDMIG3
       M2GG=1.25*MDMIG1, 1.0*MDMIG2, 0.75*MDMIG3
       SET 100=M1, M2
       M2GG=100

    """
    type = 'RealMatrixCard'
    def __init__(self, value: int | list[tuple[float, str]]):
        """
        Creates an IntCard

        Parameters
        ----------
        value : int
            the value for the card

        """
        super().__init__()
        self.value = value

    def __iter__(self):
        """temporary method to emulate the old list access style"""
        value = self
        options = self.value
        param_type = 'OBJ-type'
        return iter([value, options, param_type])

    @classmethod
    def add_from_case_control(cls, line: str):
        """
        Creates a card from the Case Control Deck

        Parameters
        ----------
        line : str
            the line of the card

        """
        value = line.split('=')[1].strip()
        if value.isdigit():
            values2 = int(value) # SET id
        else:
            sline = value.split(',')
            values2 = []
            for val in sline:
                if '*' in val:
                    assert '*' in val, f'{cls.type} val={val!r} requires a *; line={line!r}'
                    scale_str, name = val.split('*')
                    scale = float(scale_str)
                else:
                    scale = 1.
                    name = val
                    #raise RuntimeError('val={val!r}')
                name = name.strip()
                values2.append((scale, name))
            #print(values2)

        out = cls(values2)
        return out

    def export_to_hdf5(self, h5_file, encoding: str) -> None:
        reals = []
        values = []
        for real, value in self.value:
            reals.append(real)
            values.append(value)
        h5_file.create_dataset('reals', data=reals)
        h5_file.create_dataset('values', data=values)

    @classmethod
    def load_hdf5(cls, h5_file, encoding: str) -> Any:
        from pyNastran.utils.dict_to_h5py import _cast
        reals = _cast(h5_file['reals'])
        values = _cast(h5_file['values'])
        values2 = [(key, value) for key, value in zip(reals, values)]
        return cls(values2)

    def _repr_rows(self) -> list[str]:
        if isinstance(self.value, integer_types):
            rows = [f'{self.type} = {self.value:d}']
        else:
            max_chars = 72 - 4 #  -4 for indentation
            rows = []
            msg = f'{self.type} = '
            for (scale, name) in self.value:
                msgi = f'{scale}*{name},'
                if len(msg) + len(msgi) < max_chars:
                    msg += msgi
                else:
                    rows.append(msg)
                    msg = ' ' + msgi
            if msg:
                rows.append(msg)
            rows[-1] = rows[-1].rstrip(',')
        return rows

    def __repr__(self) -> str:
        """writes a card"""
        rows = self._repr_rows()
        return '\n'.join(rows)

    def write(self, spaces: str):
        """writes a card with spaces"""
        rows = self._repr_rows()
        out = spaces + ('\n' + spaces).join(rows) + '\n'
        return out


class ImagMatrixCard:
    """
    interface for cards of the form:
       M2GG=MDMIG
       M2GG=MDMIG1, MDMIG2, MDMIG3
       M2GG=1.25*MDMIG1, 1.0*MDMIG2, 0.75*MDMIG3
       SET 100=M1, M2
       M2GG=100

    """
    type = 'ImagMatrixCard'
    def __init__(self, value: int | list[tuple[float, str]]):
        """
        Creates an IntCard
        K2PP = 1

        Parameters
        ----------
        value : int
            the value for the card

        """
        super().__init__()
        #assert isinstance(value, int), value
        self.value = value

    def __iter__(self):
        """temporary method to emulate the old list access style"""
        value = self
        options = self.value
        param_type = 'OBJ-type'
        return iter([value, options, param_type])

    @classmethod
    def add_from_case_control(cls, line: str):
        """
        Creates a card from the Case Control Deck

        Parameters
        ----------
        line : str
            the line of the card

        """
        value = line.split('=')[1].strip()
        if value.isdigit():
            values2 = int(value) # SET id
        else:
            #sline = value.split(',')
            sline = split_every_other_comma(value)

            values2 = []
            for val in sline:
                if '*' in val:
                    assert '*' in val, f'{cls.type} val={val!r} requires a *; line={line!r}'
                    scale, name = val.split('*')
                    assert '(' in scale and ')' in scale and ',' in scale, scale
                    scale2 = scale.strip()[1:-1]  #  remove ()
                    real_scale_str, imag_scale_str = scale2.split(',')
                    real_scale = float(real_scale_str)
                    imag_scale = float(imag_scale_str)
                else:
                    real_scale = 1.
                    imag_scale = 1.
                    name = val
                name = name.strip()
                values2.append((real_scale, imag_scale, name))
            #print(values2)

        out = cls(values2)
        return out

    def export_to_hdf5(self, h5_file, encoding: str) -> None:
        reals = []
        imags = []
        values = []
        #print(h5_file.keys(), h5_file['param_type'])
        for real, imag, value in self.value:
            reals.append(real)
            imags.append(imag)
            values.append(value)
        h5_file.create_dataset('reals', data=reals)
        h5_file.create_dataset('imags', data=imags)
        h5_file.create_dataset('values', data=values)

    @classmethod
    def load_hdf5(cls, h5_file, encoding: str) -> Any:
        from pyNastran.utils.dict_to_h5py import _cast
        reals = _cast(h5_file['reals'])
        imags = _cast(h5_file['imags'])
        values = _cast(h5_file['values'])
        values2 = [(real, imag, value) for real, imag, value in zip(reals, imags, values)]
        return cls(values2)

    def _repr_rows(self) -> list[str]:
        if isinstance(self.value, integer_types):
            rows = [f'{self.type} = {self.value:d}']
        else:
            max_chars = 72 - 4 #  -4 for indentation
            rows = []
            msg = f'{self.type} = '
            for (real_scale, imag_scale, name) in self.value:
                msgi = f'({real_scale},{imag_scale})*{name},'
                if len(msg) + len(msgi) < max_chars:
                    msg += msgi
                else:
                    rows.append(msg)
                    msg = ' ' + msgi
            if msg:
                rows.append(msg)
            rows[-1] = rows[-1].rstrip(',')
        return rows

    def __repr__(self) -> str:
        """writes a card"""
        rows = self._repr_rows()
        return '\n'.join(rows)

    def write(self, spaces: str):
        """writes a card with spaces"""
        rows = self._repr_rows()
        out = spaces + ('\n' + spaces).join(rows) + '\n'
        return out


def split_every_other_comma(chars: str) -> list[str]:
    i0 = 0
    ncomma = 0
    sline = []
    for i, char in enumerate(chars):
        #print(f'i={i:d} -> {char!r} ncomma={ncomma:d}')
        if char == ',' and ncomma == 0:
            ncomma += 1
        elif char == ',' and ncomma == 1:
            sline.append(chars[i0:i])
            i0 = i + 1
            ncomma = 0
            #print('adding')
        #else:
            #ncomma += 1
    if i0 < len(chars):
        #print(f'final add: {chars[i0:]!r}')
        sline.append(chars[i0:])
    return sline



class M2GG(RealMatrixCard):
    """
    M2GG=MDMIG
    M2GG=MDMIG1, MDMIG2, MDMIG3
    M2GG=1.25*MDMIG1, 1.0*MDMIG2, 0.75*MDMIG3
    SET 100=M1, M2
    M2GG=100

    """
    type = 'M2GG'

class A2GG(RealMatrixCard):
    type = 'A2GG'

class B2GG(RealMatrixCard):
    type = 'B2GG'

class K2GG(RealMatrixCard):
    type = 'K2GG'

class P2G(RealMatrixCard):
    type = 'P2G'

class K42GG(RealMatrixCard):
    type = 'K42GG'


class B2PP(ImagMatrixCard):
    type = 'B2PP'

class M2PP(ImagMatrixCard):
    type = 'M2PP'

class K2PP(ImagMatrixCard):
    type = 'K2PP'
