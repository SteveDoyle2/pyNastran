"""Defines the BDFCard class that is passed into the various Nastran cards."""
from typing import List, Union, Optional, Any
from pyNastran.bdf.field_writer import print_card
from pyNastran.bdf.field_writer_16 import print_field_16
from pyNastran.bdf.cards.utils import wipe_empty_fields


class BDFCard:
    """A BDFCard is a list that has a default value of None for fields out of range."""
    def __init__(self, card: List[str], has_none: bool=True) -> None:
        """
        Parameters
        ----------
        card : List[str]
           the split values for the card
        has_none : bool; default=True
           helps with a special case to speed up runtime

        # definitely bad
        card = ['GRID', '1', '', '1.0', '2.0', '3.0']
        BDFCard(card, has_none=False)

        # definitely correct
        card = ['GRID', '1', None, '1.0', '2.0', '3.0']
        BDFCard(card, has_none=True)

        # ???
        card = ['GRID', '1', '', 1.0, 2.0, 3.0]
        BDFCard(card, has_none=True)
        """
        if has_none:
            long_fields = [print_field_16(field).strip() for field in card]
            card = wipe_empty_fields(long_fields)
        self.card = card  # type: List[Optional[str]]
        self.nfields = len(self.card)  # type: int

    def pop(self) -> Optional[str]:
        """card.pop()"""
        self.nfields -= 1
        return self.card.pop()

    def __setitem__(self, key: int, value: str) -> None:
        """card[4] = value"""
        self.card.__setitem__(key, value)

    def __getitem__(self, key: int) -> str:
        """print card[5]"""
        return self.card.__getitem__(key)

    def __getslice__(self, i: int, j: int) -> List[str]:
        """card[1:10]"""
        return self.card.__getslice__(i, j)

    def __setslice__(self, i: int, j: int, sequence: Any) -> Any:
        """card[1:10] = 2"""
        self.card.__setslice__(i, j, sequence)

    def index(self, value: str) -> int:
        """card.index(value)"""
        return self.card.index(value)

    def __repr__(self) -> str:
        """
        Prints the card as a list

        Returns
        -------
        msg : str
            the string representation of the card

        """
        #return str(self.card)
        return '%r' % self.card

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """prints the card in 8/16/16-double format"""
        return print_card(self.card, size=size, is_double=is_double)

    def __len__(self) -> int:
        """len(card)"""
        return self.nfields

    def fields(self, i: int=0, j: Optional[int]=None, defaults: Any=None) -> List[Any]:
        """
        Gets multiple fields on the card

        Parameters
        ----------
        i : int > 0
            the ith field on the card (following list notation)
        j : int / None
            int : the jth field on the card
            None : last field on the card
        defaults : List[int/float/str]
            the default value for the field (as a list)
            len(defaults)=i-j-1

        Returns
        -------
        values : List[varies]
            int/float/str
            the values on the ith-jth fields

        """
        if defaults is None:
            defaults = []
        if j is None:
            #if self.nfields is None:
                #return [None]
            j = self.nfields

        if defaults == []:
            defaults = [None] * (j - i + 1)
        out = []

        d = 0
        for n in range(i, j):
            value = self.field(n, defaults[d])
            out.append(value)
            d += 1
        return out

    def field(self, i: int,
              default: Optional[Union[int, float, str]]=None) -> Optional[Union[int, float, str]]:
        """
        Gets the ith field on the card

        Parameters
        ----------
        i : int
            the ith field on the card (following list notation)
        default : int/float/str/None
            the default value for the field

        Returns
        -------
        value : int/float/str/None
            the value on the ith field

        """
        if i < self.nfields and self.card[i] is not None and self.card[i] != '':
            return self.card[i]
        return default
