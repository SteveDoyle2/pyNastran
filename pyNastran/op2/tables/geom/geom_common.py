#pylint: disable=R0201,C0111
from struct import Struct

class SuppressLogging:
    def __init__(self):
        pass
    def debug(self, msg):
        pass
    def info(self, msg):
        pass
    def flush(self):
        pass

class SuppressFileIO:
    def __init__(self):
        pass
    def open(self, fname):
        pass
    def close(self):
        pass
    def write(self, msg):
        pass
    def flush(self):
        pass

class GeomCommon:
    def __init__(self):
        self.card_count = {}
        self.is_debug_file = False
        self._endian = b''
        self.struct_i = Struct('i')
        self.struct_2i = Struct('2i')
        self.binary_debug = SuppressFileIO()
        self.isuperelement = 0
        #self.log = SuppressLogging()

    def _read_fake(self, data: bytes, n: int) -> int:
        self.log.info(f'geom skipping {self.card_name} in {self.table_name}; ndata={len(data)-12}')
        #if (self.card_name == '' or '?' in self.card_name) and data:
            #self.show_data(data)
        #if self.table_name_str in ['GEOM3', 'DYNAMIC', 'DYNAMICS', 'GEOM4', 'EPT', 'MPT']: # 'GEOM2',
            #self.show_data(data)
            #aaa
        return len(data)

    def _read_fake_nx(self, data: bytes, n: int) -> int:
        """same as _read_fake, but casts to NX"""
        self.to_nx(f' because {self.card_name} was found')
        self.log.info(f'geom skipping {self.card_name} in {self.table_name}; ndata={len(data)-12}')
        return len(data)

    def increase_card_count(self, name: str, count_num: int=1):  # pragma: no cover
        msg = 'this should be overwritten; name=%s count_num=%s' % (name, count_num)
        raise NotImplementedError(msg)
