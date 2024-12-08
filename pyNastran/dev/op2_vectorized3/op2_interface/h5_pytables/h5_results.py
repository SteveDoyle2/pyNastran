from __future__ import annotations
import os
from typing import TYPE_CHECKING
#from itertools import count
#from collections import defaultdict

import numpy as np
#from numpy.lib import recfunctions as rfn  # random numpy package to merge arrays...
try:
    from tables import open_file, Group, Node, File
except ImportError:
    print('pytables was not found; no h5 support.  Run ">>> pip install tables"\n'
          'Do you have h5py installed?  That can cause conflicts.')
    #raise

#from .nodes import load_h5_node
#from .coords import load_h5_coord
#from .properties import load_h5_property
#from .materials import load_h5_material
#from .elements import load_h5_element
#from pyNastran.op2.op2_interface.op2_classes import (
    #ComplexDisplacementArray, ComplexLoadVectorArray, ComplexSPCForcesArray, ComplexMPCForcesArray
#)
from pyNastran.dev.bdf_vectorized3.bdf_interface.h5_pytables.utils import get_group_name, get_attributes
from .nodal import read_nodal_result
from .stress_strain import read_elemental_stress, read_elemental_strain, read_elemental_force
from pyNastran.utils import print_bad_path
#from .utils import get_name
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.op2_vectorized3.op2_geom import OP2, OP2Geom


USE_PANDAS = True
if USE_PANDAS:
    import pandas as pd

# -------------------------------------------------------------------------------------------
class MyDataFrame:
    def __init__(self, data, names):
        self.data = data
        self.names = list(names)
        data2 = {name: self.data[name] for name in self.names}
        x = 1
    def __getattr__(self, name):
        if name == 'loc':
            return self
        else:
            raise NotImplementedError(name)
        x = 1
    def __getitem__(self, value):
        if isinstance(value, str):
            i = self.names.index(value)
            return self.data[value]
        elif isinstance(value, np.ndarray):
            dtype_name = value.dtype.name
            if dtype_name == 'bool':
                data2 = {name: self.data[name][value] for name in self.names}
                return data2
            else:
                raise RuntimeError(value)
        else:
            raise TypeError((value, type(value)))
        raise RuntimeError((value, type(value)))
        #return self.data[value]

# -------------------------------------------------------------------------------------------

def read_h5_geometry_result(model: OP2Geom, h5_filename: str, root_path: str='/'):
    assert os.path.exists(h5_filename), print_bad_path(h5_filename)
    from pyNastran.dev.bdf_vectorized3.bdf_interface.h5_pytables.h5_geometry import read_geometry_from_h5
    with open_file(h5_filename, mode="r", title="", root_uep="/", filters=None) as h5_file:
        read_geometry_from_h5(model, h5_file, '/NASTRAN/INPUT/')
        read_result_from_h5(model, h5_file, '/NASTRAN/', '/INDEX/NASTRAN/')

def read_h5_result(model: OP2, h5_filename: str, root_path: str='/'):
    assert os.path.exists(h5_filename), print_bad_path(h5_filename)
    with open_file(h5_filename, mode="r", title="", root_uep="/", filters=None) as h5_file:
        read_result_from_h5(model, h5_file,
                            '/NASTRAN/', '/INDEX/NASTRAN/')

def h5py_to_dataframe(group: h5py._hl.dataset.Dataset) -> pd.DataFrame:
    data = {}
    for key in group.dtype.names:
        data[key] = group[key]
    return pd.DataFrame(data)

def get_result_domains(h5_file: File, node: Node):
    data = node.read()
    #dtype([('ID', '<i8'), ('SUBCASE', '<i8'), ('STEP', '<i8'), ('ANALYSIS', '<i8'),
    # ('TIME_FREQ_EIGR', '<f8'), ('EIGI', '<f8'), ('MODE', '<i8'), ('DESIGN_CYCLE', '<i8'),
    # ('RANDOM', '<i8'), ('SE', '<i8'), ('AFPM', '<i8'), ('TRMC', '<i8'), ('INSTANCE', '<i8'),
    # ('MODULE', '<i8'), ('SUBSTEP', '<i8'), ('IMPFID', '<i8')])
    #z = np.zeros((2,2), dtype='U2')
    #o = np.ones((2,1), dtype='O')
    #np.hstack([o, z])

    if USE_PANDAS:
        domains = h5py_to_dataframe(data)
    else:
        domains = MyDataFrame(data, data.dtype.names)
    # domains =  rfn.merge_arrays([data[key] for key in data.dtype.names])
    #domains = np.hstack([idi, subcase, step, analysis, time_freq_eigr, eigi, mode, design_cycle,
                         #random, se, afpm, trmc, instance, module, substep, impfid])
    attributes = get_attributes(node)
    assert len(attributes) == 1, attributes
    version = attributes['version'][0]  # 20200
    return domains


def read_elemental_result(model: OP2, domains: np.ndarray,
                          elemental: Node, elemental_index: Node):
    for h5_node_ in elemental._f_iter_nodes():
        if isinstance(h5_node_, Group):
            name = get_group_name(h5_node_)
            if name == 'ELEMENT_FORCE':
                read_elemental_force(model, domains,
                                     elemental[name], elemental_index[name])
            elif name == 'ENERGY':
                print(h5_node_)
            elif name == 'STRAIN':
                read_elemental_strain(model, domains,
                                      elemental[name], elemental_index[name])
            elif name == 'STRESS':
                read_elemental_stress(model, domains,
                                      elemental[name], elemental_index[name])
            else:
                print(h5_node_)
                raise NotImplementedError(name)

        elif isinstance(h5_node_, Node):
            data = h5_node_.read()
            assert isinstance(data, np.ndarray), data
            name = h5_node_.name
            group = elemental[name].read()
            index = elemental_index[name].read()
            if name == 'APPLIED_LOAD_CPLX':
                fake
            else:
                print(h5_node_)
                raise NotImplementedError(name)
        else:
            print(h5_node_)
            raise NotImplementedError(h5_node_)
    x = 1

def _set_date(model: OP2) -> None:
    #'Tue Dec 29 13:07:02 2020 (UTC-8)'
    time = model.time

    base, utc = time[:-1].split('(')
    #base
    #'Tue Dec 29 13:07:02 2020 '
    #utc
    #'UTC-8'
    day_of_week, month_str, day_str, time, year_str = base.strip().split()
    year = int(year_str)
    day = int(day_str)
    month_str_to_int = {
        'Jan': 1,
        'Feb': 2,
        'Mar': 3,
        'Apr': 4,
        'May': 5,
        'Jun': 6,
        'Jul': 7,
        'Aug': 8,
        'Sep': 9,
        'Oct': 10,
        'Nov': 11,
        'Dec': 12,
    }
    month_int = month_str_to_int[month_str]
    model.date = (month_int, day, year)

def read_result_from_h5(model: OP2, h5_file: File,
                        nastran_path: str,
                        index_nastran_path: str, ):

    nastran = h5_file.get_node(nastran_path)
    attributes = get_attributes(nastran)
    assert len(attributes) == 7, attributes
    model.bdf_filename = attributes['INPUT'].split('\n')[0].strip()# .decode('latin1')
    model.time = attributes['TIME']# .decode('latin1')
    model.nastran_format = attributes['VERSION'] # .decode('latin1')

    _set_date(model)

    result_path = nastran_path + 'RESULT/'
    index_result_path = index_nastran_path + 'RESULT/'
    index_result = h5_file.get_node(index_result_path)
    result = h5_file.get_node(result_path)

    domains_node = h5_file.get_node(result_path + '/DOMAINS')
    domains = get_result_domains(h5_file, domains_node)

    elemental_index = h5_file.get_node(index_result_path + '/ELEMENTAL')
    elemental = h5_file.get_node(result_path + '/ELEMENTAL')
    read_elemental_result(model, domains, elemental, elemental_index)

    nodal_index = h5_file.get_node(index_result_path + '/NODAL')
    nodal = h5_file.get_node(result_path + '/NODAL')
    read_nodal_result(model, domains, nodal, nodal_index)

    for h5_node_ in result._f_iter_nodes():
        if isinstance(h5_node_, Group):
            name = get_group_name(h5_node_)
            if name == 'ELEMENTAL':
                pass
                #read_elemental_result(model, domains, h5_node_)
            elif name == 'NODAL':
                pass
            elif name in {'AERODYNAMIC', 'MATRIX', 'MONITOR', 'SUMMARY'}:
                pass
                #read_nodal_result(model, domains, h5_node_)
            else:
                print(h5_node_)
                raise NotImplementedError(name)
        elif isinstance(h5_node_, Node):
            data = h5_node_.read()
            assert isinstance(data, np.ndarray), data
            name = h5_node_.name
            if name == 'DOMAINS':
                pass
                #get_result_domains(h5_file, h5_node_)
            else:
                print(h5_node_)
                raise NotImplementedError(name)
        else:
            print(h5_node_)
            raise NotImplementedError(h5_node_)
        print(f' - name {name!r}')
    x = 1
