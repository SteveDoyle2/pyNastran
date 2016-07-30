from struct import Struct

class SuppressLogging(object):
    def __init__(self):
        pass
    def debug(self, msg):
        pass
    def info(self, msg):
        pass
    def flush(self):
        pass

class SuppressFileIO(object):
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

class GeomCommon(object):
    def __init__(self):
        self.card_count = {}
        self.is_debug_file = False
        self._endian = ''
        self.struct_i = Struct('i')
        self.struct_2i = Struct('2i')
        self.binary_debug = SuppressFileIO()
        self.log = SuppressLogging()

    def _read_fake(self, data, n):
        return len(data)

    def _increase_card_count(self, name, count_num=1):
        msg = 'this should be overwritten; name=%s count_num=%s' % (name, count_num)
        raise NotImplementedError(msg)

    def add_coord(self, coord, allow_overwrites=True):
        raise RuntimeError('this should be overwritten by the BDF class')

    def add_property(self, card, allow_overwrites=True):
        raise RuntimeError('this should be overwritten')

    def add_constraint_SPC(self, constraint):
        raise RuntimeError('this should be overwritten by the BDF class')

    def add_rigid_element(self, constraint):
        raise RuntimeError('this should be overwritten by the BDF class')

    def add_suport(self, constraint):
        raise RuntimeError('this should be overwritten by the BDF class')

    def add_thermal_load(self, load):
        raise RuntimeError('this should be overwritten by the BDF class')

    def add_load(self, load):
        raise RuntimeError('this should be overwritten by the BDF class')


    def add_TSTEPNL(self, card, allow_overwrites=True):
        raise RuntimeError('this should be overwritten by the BDF class')

    def add_NLPARM(self, card, allow_overwrites=True):
        raise RuntimeError('this should be overwritten by the BDF class')

    def add_material_dependence(self, material, allow_overwrites=True):
        raise RuntimeError('this should be overwritten by the BDF class')

    def add_creep_material(self, material, allow_overwrites=True):
        raise RuntimeError('this should be overwritten by the BDF class')

    def add_structural_material(self, material, allow_overwrites=True):
        raise RuntimeError('this should be overwritten by the BDF class')

    def add_thermal_material(self, material, allow_overwrites=True):
        raise RuntimeError('this should be overwritten by the BDF class')

    def add_gust(self, gust):
        raise RuntimeError('this should be overwritten by the BDF class')

    def add_table(self, table):
        raise RuntimeError('this should be overwritten by the BDF class')

