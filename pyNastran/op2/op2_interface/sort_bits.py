#from collections import namedtuple
#SortBits = namedtuple('SortBits', ['is_sort1', 'is_real', 'is_random'])
INV_MAP = {
    1: 0,
    0: 1,
}
class SortBits(list):
    """
    Bit Description
    === ===========
    1 Complex (on) flag
    2 SORT2 (on) flag
    3 Random (on) flag

    Value Sort type Data format Random
    ===== ========= =========== ======
    0     SORT1     Real        No
    1     SORT1     Complex     No
    2     SORT2     Real        No
    3     SORT2     Complex     No
    4     SORT1     Real        Yes
    6     SORT2     Real        Yes

    """
    @classmethod
    def add_from_sort_code(self, sort_code: int, is_table_1: bool) -> None:
        """
        Parameters
        ----------
        sort_code : int
            0-6 integer
        is_table_1 : bool
            is this a SORT1 table (e.g., OES1X)

        """
        i = 2
        bits = [0, 0, 0]
        while sort_code > 0:
            value = sort_code % 2
            sort_code = (sort_code - value) // 2
            bits[i] = value
            i -= 1

        # fixing bit[1]
        bits[1] = 0 if is_table_1 else 1
        return SortBits(bits)

    @property
    def is_real(self) -> int:
        return INV_MAP[self.is_complex]
    @property
    def is_complex(self) -> int:
        return self[0]
    @is_real.setter
    def is_real(self, value: int) -> int:
        self[0] = INV_MAP[value]
    @is_complex.setter
    def is_complex(self, value) -> int:
        self[0] = value

    @property
    def is_sort1(self) -> int:
        return INV_MAP[self[1]]
    @is_sort1.setter
    def is_sort1(self, value: int) -> int:
        self[1] = INV_MAP[value]
    @property
    def is_sort2(self) -> int:
        return self[1]
    @is_sort2.setter
    def is_sort2(self, value: int) -> int:
        self[1] = value

    @property
    def is_random(self) -> int:
        return self[2]
    @is_random.setter
    def is_random(self, value: int) -> int:
        self[2] = value

    def __repr__(self):
        return f'SortBits(is_complex={self[0]}, is_sort2={self[1]}, is_random={self[2]})'

sort_bits = SortBits([0, 0, 0])
assert sort_bits.is_real == 1, sort_bits
assert sort_bits.is_complex == 0, sort_bits

assert sort_bits.is_sort1 == 1, sort_bits
assert sort_bits.is_sort2 == 0, sort_bits

sort_bits.is_sort1 = 1
#sort_bits.is_sort2 = 0

assert sort_bits.is_sort1 == 1, sort_bits
assert sort_bits.is_sort2 == 0, sort_bits
