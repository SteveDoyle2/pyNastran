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

    def _read_fake(self, data, n):
        self.log.info(f'skipping {self.card_name} in {self.table_name}; ndata={len(data)-12}')
        #if (self.card_name == '' or '?' in self.card_name) and data:
            #self.show_data(data)
        #if self.table_name_str in ['GEOM3', 'DYNAMIC', 'DYNAMICS', 'GEOM4', 'EPT', 'MPT']: # 'GEOM2',
            #self.show_data(data)
            #aaa
        return len(data)

    def _read_fake_nx(self, data, n):
        """same as _read_fake, but casts to NX"""
        self.to_nx()
        self.log.info(f'skipping {self.card_name} in {self.table_name}; ndata={len(data)-12}')
        return len(data)

    def increase_card_count(self, name, count_num=1):  # pragma: no cover
        msg = 'this should be overwritten; name=%s count_num=%s' % (name, count_num)
        raise NotImplementedError(msg)

    def _add_coord_object(self, coord, allow_overwrites=True):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_property_object(self, card, allow_overwrites=True):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_constraint_spc_object(self, constraint):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_constraint_spcoff_object(self, constraint):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_rigid_element_object(self, constraint):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_suport_object(self, constraint):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_thermal_load_object(self, load):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_load_object(self, load):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_tstepnl_object(self, card, allow_overwrites=True):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_nlparm_object(self, card, allow_overwrites=True):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_material_dependence_object(self, material, allow_overwrites=True):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_creep_material_object(self, material, allow_overwrites=True):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_structural_material(self, material, allow_overwrites=True):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_thermal_material_object(self, material, allow_overwrites=True):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_gust_object(self, gust):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_table_object(self, table):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_aset_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_bset_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_cset_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_uset_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_seqset_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_constraint_spcadd_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_constraint_mpcadd_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_thermal_element_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_mass_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_spoint_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_plotel_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_load_combination_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_lseq_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_constraint_mpc_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_thermal_bc_object(self, obj, key):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_qset_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_element_object(self, obj):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_structural_material_object(self, allow_overwrites=True):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')

    def _add_hyperelastic_material_object(self, allow_overwrites=True):  # pragma: no cover
        raise RuntimeError('this should be overwritten by the BDF class')
