
class CaseControlCard:
    """basic card similar to the BaseCard class for the BDF"""
    def __iter__(self):
        """temporary method to emulate the old list access style"""
        value = self
        options = None
        param_type = 'OBJ-type'
        return iter([value, options, param_type])


class IntCard(CaseControlCard):
    """
    interface for cards of the form:
       NAME = 10

    """
    type = 'IntCard'
    def __init__(self, value):
        """
        Creates an IntCard

        Parameters
        ----------
        value : int
            the value for the card

        """
        super(IntCard, self).__init__()
        self.value = int(value)

    def __iter__(self):
        """temporary method to emulate the old list access style"""
        value = self
        options = []
        #param_type = 'STRESS-type'
        param_type = 'OBJ-type'
        return iter([value, options, param_type])

    @classmethod
    def add_from_case_control(cls, line: str, line_upper: str, lines: list[str], i: int):
        """
        Creates a card from the Case Control Deck

        Parameters
        ----------
        line : str
            the line of the card
        line_upper : str
            unused
        lines : list[str]
            unused
        i : int
            unused

        """
        value = line_upper.split('=')[1]
        try:
            out = cls(value)
        except ValueError:
            print(line)
            raise
        return out

    def export_to_hdf5(self, h5_file, encoding):
        h5_file.create_dataset('value', data=self.value)

    @classmethod
    def load_hdf5(cls, h5_file, encoding):
        from pyNastran.utils.dict_to_h5py import _cast
        value = h5_file['value']
        value2 = _cast(value)
        return cls(value2), []

    def __repr__(self):
        """writes a card"""
        return '%s = %i\n' % (self.type, self.value)

    def write(self, spaces):
        """writes a card with spaces"""
        return spaces + str(self)


class IntStrCard(IntCard):
    """
    interface for cards of the form:
       NAME = 10
       NAME = ALL

    """
    type = 'IntStrCard'
    allowed_strings: set[str] = set([])
    def __init__(self, value):
        """
        Creates an IntStrCard

        Parameters
        ----------
        value : int/str
            the value for the card

        """
        #super(IntStrCard, self).__init__()
        try:
            self.value = int(value)
        except ValueError:
            value = value.strip()
            if value not in self.allowed_strings:
                msg = 'value=%r not in [%s]' % (
                    value, ', '.join(self.allowed_strings))
                raise ValueError(msg)
            self.value = value

    def export_to_hdf5(self, h5_file, encoding):
        value_bytes = self.value.encode(encoding) if isinstance(self.value, str) else self.value
        #sub_group = h5_file.create_group(self.type)
        h5_file.create_dataset('value', data=value_bytes)

    @classmethod
    def load_hdf5(cls, h5_file, encoding):
        from pyNastran.utils.dict_to_h5py import _cast
        value = h5_file['value']

        casted_value = _cast(value)
        if isinstance(casted_value, int):
            value2 = casted_value
        else:
            value2 = casted_value.decode(encoding)# if isinstance(value, bytes) else value
        return cls(value2), []

    def __repr__(self):
        """writes a card"""
        return '%s = %s\n' % (self.type, self.value)


class StringCard(CaseControlCard):
    type = 'StringCard'
    allowed_values: list[str] = []
    def __init__(self, value, validate=True):
        super(StringCard, self).__init__()
        self.value = value.strip()
        if validate:
            self.validate()

    @classmethod
    def add_from_case_control(cls, line, line_upper, lines, i):
        """add method used by the CaseControl class"""
        value = line_upper.split('=')[1]
        return cls(value)

    def validate(self):
        if self.value not in self.allowed_values:
            msg = '%s: value=%r not in [%s]' % (
                self.type, self.value, ', '.join(self.allowed_values))
            raise ValueError(msg)

    def __repr__(self):
        """writes a card"""
        return '%s = %s\n' % (self.type, self.value)

    def write(self, spaces):
        return spaces + str(self)

    def export_to_hdf5(self, h5_file, encoding):
        value_bytes = self.value.encode(encoding)
        #sub_group = h5_file.create_group(self.type)
        h5_file.create_dataset('value', data=value_bytes)

    @classmethod
    def load_hdf5(cls, h5_file, encoding):
        from pyNastran.utils.dict_to_h5py import _cast
        value = h5_file['value']
        try:
            value2 = _cast(value).decode(encoding)
        except AttributeError:
            print(cls.type, _cast(value))
            raise
        return cls(value2), []
