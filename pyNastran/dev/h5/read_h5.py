"""
https://public.kitware.com/pipermail/paraview/2016-January/035894.html
https://discourse.paraview.org/t/paraviewweb-visualizer/1268/5
https://discourse.vtk.org/t/appending-data-field-to-a-vtk-file-python-vtk/3220/3
https://www.paraview.org/Wiki/Python_Programmable_Filter
https://stackoverflow.com/questions/54603267/how-to-show-vtkunstructuredgrid-in-python-script-based-on-paraview/54633793#54633793
https://vtkpythonpackage.readthedocs.io/en/latest/index.html
https://cvw.cac.cornell.edu/ParaViewAdv/pythonshell
"""
from itertools import count
#from pprint import pprint
from typing import Optional, Any, TYPE_CHECKING

import os
import sys
import vtk # vtk > 9
#import vtkmodules
#vtkmodules.vtkCommonCore.
import h5py
import numpy as np
import pandas as pd
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from pyNastran.utils import print_bad_path # , object_attributes, object_methods, object_stats
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.dev.h5.h5_utils import get_tree, h5py_to_dataframe
from pyNastran.dev.h5.h5_result_objects import RealVectorTable, RealVectorTableOptistruct, RealStrainEnergyOptistruct

from pyNastran.dev.h5.geometry.h5_case_control import load_case_control, load_parameters
from pyNastran.dev.h5.geometry.h5_elements import element_map
from pyNastran.dev.h5.geometry.h5_constraints import constraint_map
from pyNastran.dev.h5.geometry.h5_loads import load_map
from pyNastran.dev.h5.geometry.h5_tables import table_map
from pyNastran.dev.h5.geometry.h5_geometry import (
    coord_map, node_map,
    material_map,
    matrix_map, design_map, dynamic_map, partition_map,
    load_geometry_block)
from pyNastran.dev.h5.geometry.h5_properties import property_map
#from pyNastran.gui.utils.vtk.base_utils import numpy_to_vtk

Function = Any
from cpylog import SimpleLogger

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

class BDF2:
    def __init__(self):
        self.log = SimpleLogger(level='debug', encoding='utf-8')
        self.card_count = {}
        self.CTRIA3 = None
        self.CQUAD4 = None
        self.CTRIA6 = None
        self.CQUAD8 = None
        self.CTETRA = None
        self.CPENTA = None
        self.CPYRAM = None
        self.CHEXA = None

def break_domain_by_case(domains_df: pd.DataFrame, INDEX_DOMAIN) -> pd.DataFrame:
    #subcase = domains_df['SUBCASE']
    #analysis_code = domains_df['ANALYSIS']
    #eigr = mycase['TIME_FREQ_EIGR']
    #eigi = mycase['EIGI']
    #mode = mycase['MODE']

    subcase = domains_df['SUBCASE']
    isubcase = subcase != 0

    keys = ['SUBCASE','ANALYSIS', 'STEP', 'DESIGN_CYCLE', 'RANDOM', 'SE', 'AFPM', 'TRMC', 'INSTANCE', 'MODULE']
    domain_no_zero = domains_df.loc[isubcase]
    grouped = domain_no_zero.groupby(keys)
    #for g in grouped:
        #print(g)
    grouped_df = grouped.size().reset_index().rename(columns={0:'count'})[keys]
    #step = mycase['STEP']
    #design_cycle = mycase['DESIGN_CYCLE']
    #random = mycase['RANDOM']
    #se = mycase['SE']
    #afpm = mycase['AFPM']
    #trmc = mycase['TRMC']
    #inst = mycase['INSTANCE']
    #module = mycase['MODULE']
    return grouped # _df


class pyNastranH5:
    def __init__(self,
                 add_aero: bool=True,
                 add_constraints: bool=True,
                 add_results: bool=True,
                 subcases=None):
        self.filename = None
        self.geom_model = None
        self.add_aero = add_aero
        self.add_constraints = add_constraints
        self.add_results = add_results
        self.flags = {
            'aero': add_aero,
            'constraint': add_constraints,
            'results': add_results,
            'subcases': subcases,
        }
        self.subcases = subcases
        self.results = {}
        self.log = SimpleLogger(level='debug', encoding='utf-8')

    def read_h5_nastran(self, h5_filename: str,
                        subcases: Optional[list[int]]=None) -> None:
        self.filename = h5_filename
        print(f'opening {h5_filename}')
        assert os.path.exists(h5_filename), print_bad_path(h5_filename)
        hdf5_file = h5py.File(h5_filename, 'r')
        # ['INDEX', 'NASTRAN']

        """
        {'INDEX': {
            'NASTRAN': {
                'RESULT': {
                    'ELEMENTAL': {
                        'STRESS': {'QUAD4': None}},
                    'NODAL': {
                        'EIGENVECTOR': None}}}},
        'NASTRAN': {
            'INPUT': {
                'DOMAINS': None,
                'DYNAMIC': {'EIGRL': {'IDENTITY': None}},
                'ELEMENT': {'CQUAD4': None},
                'MATERIAL': {'MAT1': None},
                'NODE': {'GRID': None},
                'PARAMETER': {'CASECC': {'SUBCASE': None},
                'MDLPRM': None,
                'PVT': {'INT': None}},
                'PROPERTY': {'PSHELL': None}
            },
            'RESULT': {
                'DOMAINS': None,
                'ELEMENTAL': {'STRESS': {'QUAD4': None}},
                'NODAL': {'EIGENVECTOR': None},
                'SUMMARY': {'EIGENVALUE': None}}}}
        """
        tree = get_tree(hdf5_file)
        #pprint(tree)

        self._load_geometry(hdf5_file)
        #return
        self._load_results(hdf5_file, subcases=subcases)

    def _load_geometry(self, hdf5_file: h5py.File) -> None:
        base_str = 'NASTRAN/INPUT'

        inputs = hdf5_file.get(base_str)
        #['CONSTRAINT', 'COORDINATE_SYSTEM', 'DESIGN', 'DOMAINS', 'ELEMENT', 'LOAD', 'MATERIAL', 'NODE', 'PARAMETER', 'PROPERTY', 'TABLE']
        quiet_skip_inputs = ['DOMAINS', 'PARAMETER']
        geom_model = BDF()
        geom_model.flags = self.flags

        log = geom_model.log
        for geom_name in list(inputs):
            if geom_name in quiet_skip_inputs:
                continue

            geom_group = inputs.get(geom_name)
            if geom_name == 'CONSTRAINT':
                #if not self.add_constraints:
                    #continue
                load_geometry_block(geom_group, constraint_map, geom_model)

            elif geom_name == 'COORDINATE_SYSTEM':
                load_geometry_block(geom_group, coord_map, geom_model)
            elif geom_name == 'DESIGN':
                load_geometry_block(geom_group, design_map, geom_model)
            #elif geom_name == 'DOMAINS':
            elif geom_name == 'DYNAMIC':
                load_geometry_block(geom_group, dynamic_map, geom_model)
            elif geom_name == 'ELEMENT':
                load_geometry_block(geom_group, element_map, geom_model)
            elif geom_name == 'LOAD':
                load_geometry_block(geom_group, load_map, geom_model)
            elif geom_name == 'MATERIAL':
                #if not self.add_material:
                    #continue
                load_geometry_block(geom_group, material_map, geom_model)
            elif geom_name == 'MATRIX':
                load_geometry_block(geom_group, matrix_map, geom_model)
            elif geom_name == 'NODE':
                load_geometry_block(geom_group, node_map, geom_model)
            elif geom_name == 'PARTITION':
                load_geometry_block(geom_group, partition_map, geom_model)
            elif geom_name == 'PROPERTY':
                load_geometry_block(geom_group, property_map, geom_model)
            elif geom_name == 'TABLE':
                load_geometry_block(geom_group, table_map, geom_model)
            else:
                log.warning(f'skipping {geom_name}')

        load_parameters(hdf5_file, geom_model)

        cc_str = 'NASTRAN/INPUT/PARAMETER/CASECC/SUBCASE'
        cc = hdf5_file.get(cc_str)
        assert cc is not None, cc
        load_case_control(geom_model, cc)

        #geom_model.write_bdf(r'C:\NASA\m4\formats\git\pyNastran\models\msc\spike.bdf')
        finish_geometry(geom_model)
        self.geom_model = geom_model

    def _load_results(self, hdf5_file: h5py.File,
                      subcases: Optional[list[int]]=None) -> None:
        geom_model = self.geom_model
        node_ids = geom_model._node_ids
        element_ids = geom_model._element_ids
        domains_str = '/NASTRAN/RESULT/DOMAINS'
        domains = hdf5_file.get(domains_str)

        #    ID  SUBCASE  STEP  ANALYSIS  TIME_FREQ_EIGR  EIGI  MODE  DESIGN_CYCLE  RANDOM  SE  AFPM  TRMC  INSTANCE  MODULE
        # 0   1        0     0         0    0.000000e+00   0.0     0             0       0   0     0     0         0       0
        # 1   2        1     0         2   -3.087735e-10   0.0     1             0       0   0     0     0         0       0
        # 2   3        1     0         2   -2.082743e-10   0.0     2             0       0   0     0     0         0       0
        # 3   4        1     0         2   -1.514309e-10   0.0     3             0       0   0     0     0         0       0


        domains_df = h5py_to_dataframe(domains)
        assert domains is not None, domains
        assert domains_df is not None, domains_df

        result_str = '/NASTRAN/RESULT'
        result_index_str = '/INDEX/NASTRAN/RESULT'
        #stress = hdf5_file.get('/INDEX/NASTRAN/RESULT/ELEMENTAL/STRESS')
        #pprint(stress)
        #pprint(get_tree(stress))
        results_dict = {}
        iresult = 0
        model = OP2()
        log = model.log
        results = hdf5_file.get(result_str)
        for res_name in results:
            #results_group = hdf5_file.get(name)
            if res_name in ['NODAL', 'ELEMENTAL',
                            'AERODYNAMIC', 'DOMAINS', 'MATRIX', 'MONITOR', 'OPTIMIZATION', 'SUMMARY']:
                pass
            else:
                log.warning(f'skipping result {res_name}')

        #iresult = self._load_nodal_results(
            #iresult, results_dict,
            #geom_model, model, node_ids,
            #domains_df, hdf5_file, subcases=subcases)
        #return
        # --------------------------------------------------------------
        iresult = self._load_elemental_results(
            iresult,
            result_str, result_index_str,
            results_dict,
            geom_model, model, element_ids,
            domains_df, hdf5_file, subcases=subcases)
        self.results = results_dict

        #geom_model.write_bdf(r'C:\NASA\m4\formats\git\pyNastran\models\msc\6+element-nastran-sol103.bdf')
        #self.results = results  # type: dict[int, RealVectorTableOptistruct]
        #self.results_model = model

    def _load_elemental_results(self, iresult: int,
                                result_str: str, result_index_str: str,
                                results: dict[int, Any],
                                geom_model: BDF,
                                model: OP2,
                                element_ids: np.ndarray,
                                domains_df: pd.DataFrame,
                                hdf5_file: h5py.File,
                                subcases: Optional[list[int]]=None) -> int:

        elemental_str = result_str + '/ELEMENTAL'
        elemental_index_str = result_index_str + '/ELEMENTAL'
        assert domains_df is not None
        element = hdf5_file.get(elemental_str)
        element_index = hdf5_file.get(elemental_index_str)
        #pprint(get_tree(element))
        #pprint(get_tree(element_index))

        basename = ''

        for ires, name in enumerate(element):
            group = element.get(name)
            index = element_index.get(name)
            if name == 'ENERGY':
                iresult = load_strain_energy(basename, iresult, results,
                                             domains_df, group, index,
                                             element_ids,
                                             geom_model, model,
                                             subcases=subcases)
            elif name in ['STRESS', 'STRAIN']:
                is_stress = name == 'STRESS'
                iresult = load_stress_strain(basename, iresult, results,
                                             domains_df, group, index,
                                             element_ids,
                                             geom_model, model,
                                             is_stress=is_stress,
                                             subcases=subcases)
            elif name in ['ELEMENT_FORCE', 'FAILURE_INDEX']:
                pass
                #load_eigenvector(results, i, domains_df, group, index, model)
            else:
                raise NotImplementedError(name)

        #element = hdf5_file.get(element_str)
        #pprint(get_tree(element))
        #load_geometry_block(element, element_map, geom_model)

    def _load_nodal_results(self, iresult: int,
                            results: dict[int, Any],
                            geom_model: BDF,
                            model: OP2,
                            node_ids: np.ndarray,
                            domains_df: pd.DataFrame,
                            hdf5_file: h5py.File,
                            subcases: Optional[list[int]]=None) -> int:
        assert domains_df is not None
        nodal_str = '/NASTRAN/RESULT/NODAL'
        nodal_index_str = '/INDEX/NASTRAN/RESULT/NODAL'

        name_map = {
            'EIGENVECTOR': 'Eigenvector',
            'EIGENVECTOR_CPLX': 'Eigenvector',

            'DISPLACEMENT': 'Displacement',
            'DISPLACEMENT_CPLX': 'Displacement',

            'VELOCITY': 'Velocity',

            'APPLIED_LOAD': 'Applied Load',
            'APPLIED_LOAD_CPLX': 'Applied Load',

            'SPC_FORCE': 'SPC Force',
            'SPC_FORCE_CPLX': 'SPC Force',

            'MPC_FORCE': 'MPC Force',
            'MPC_FORCE_CPLX': 'MPC Force',

            'GRID_FORCE': 'Grid Point Force',
            'KINETIC_ENERGY': 'Kinetic Energy',
            'GRID_WEIGHT': 'Grid Point Weight',
            'TEMPERATURE': 'Temperature',
            #'KINETIC_ENERGY': 'Eigenvector',
        }

        node = hdf5_file.get(nodal_str)
        node_index = hdf5_file.get(nodal_index_str)
        #pprint(get_tree(node))
        #pprint(get_tree(node_index))
        for ires, name in enumerate(node):
            group = node.get(name)
            index = node_index.get(name)
            basename = name_map[name]
            if name == 'EIGENVECTOR':
                iresult = load_eigenvector(basename, iresult, results, domains_df, group, index,
                                           node_ids, geom_model, model, subcases=subcases)
            elif name in ['EIGENVECTOR_CPLX', 'APPLIED_LOAD_CPLX', 'DISPLACEMENT_CPLX', 'SPC_FORCE_CPLX', 'MPC_FORCE_CPLX']:
                iresult = load_eigenvector_complex(basename, iresult, results, domains_df, group, index,
                                                   node_ids, geom_model, model, subcases=subcases)
            elif name in ['DISPLACEMENT', 'VELOCITY']:
                iresult = load_eigenvector(basename, iresult, results, domains_df, group, index,
                                           node_ids, geom_model, model, subcases=subcases)
            elif name in ['SPC_FORCE', 'MPC_FORCE', 'APPLIED_LOAD']:
                iresult = load_eigenvector(basename, iresult, results, domains_df, group, index,
                                           node_ids, geom_model, model, subcases=subcases)
            elif name in ['GRID_FORCE', 'KINETIC_ENERGY', 'GRID_WEIGHT', 'TEMPERATURE']:
                # GRID_POINT_FORCE
                pass
            else:
                raise NotImplementedError(name)
        #load_geometry_block(node, node_map, model)
        return iresult

def finish_geometry(geom_model: BDF):
    from pyNastran.dev.h5.fill_unstructured_grid import _load_nodes
    nodes, node_ids, nid_map, idtype = _load_nodes(geom_model)

    geom_model._nodes = nodes
    geom_model._node_ids = node_ids
    geom_model._element_ids = None
    geom_model._nid_map = nid_map
    geom_model._idtype = idtype

def load_strain_energy(basename_orig: str,
                       iresult: int,
                       results: dict[int, Function],

                       domains_df: pd.DataFrame,
                       element_group: h5py._hl.group.Group,
                       element_index: h5py._hl.group.Group,

                       ids: np.ndarray,
                       geom_model: BDF,
                       model: OP2,
                       subcases: Optional[list[int]]=None) -> int:
    return iresult
    basename = 'Strain Energy'
    #assert ids is not None
    for ires, element_group_name in enumerate(element_group):
        group = element_group.get(element_group_name)
        index = element_index.get(element_group_name)
        if element_group_name == 'IDENT':
            #group = <class 'h5py._hl.dataset.Dataset'>
            #('IDENT', 'ELNAME', 'ETOTAL', 'CVALRES', 'ESUBT', 'ETOTPOS', 'ETOTNEG')
            #IDENT = array([ 1,  2,  3,  4,  5])
            #ELNAME = array([b'BAR     ', b'HEXA    ', b'PENTA   '])
            pass
        elif element_group_name == 'STRAIN_ELEM':
            #group = <class 'h5py._hl.dataset.Dataset'>
            #('ID', 'ENERGY', 'PCT', 'DEN', 'IDENT', 'DOMAIN_ID')

            #index = <class 'h5py._hl.dataset.Dataset'>
            #('DOMAIN_ID', 'POSITION', 'LENGTH')
            INDEX_DOMAIN = index['DOMAIN_ID']
            INDEX_POSITION = index['POSITION']
            INDEX_LENGTH = index['LENGTH']

            EID_BASE = group['ID']
            ieid = ~(EID_BASE == 100000000)
            EID = EID_BASE[ieid]
            ENERGY = group['ENERGY', ieid]
            PCT = group['PCT', ieid]
            DEN = group['DEN', ieid]
            IDENT = group['IDENT', ieid]
            DOMAIN = group['DOMAIN_ID', ieid]

            grouped_df = break_domain_by_case(domains_df, INDEX_DOMAIN)
            subcases = geom_model.subcases
            for tuplei in grouped_df:
                indexi, dfi = tuplei
                DOMAINs = dfi['ID']

                # dfi
                #       ID  SUBCASE  STEP  ANALYSIS  TIME_FREQ_EIGR  EIGI  MODE  DESIGN_CYCLE  RANDOM  SE  AFPM  TRMC  INSTANCE  MODULE
                # 0      1        1     0         2   -2.572570e-07   0.0     1             0       0   0     0     0         0       0
                # 1      2        1     0         2    8.711537e-07   0.0     2             0       0   0     0     0         0       0
                # 2      3        1     0         2    2.950070e-06   0.0     3             0       0   0     0     0         0       0
                # 3      4        1     0         2    3.668822e-06   0.0     4             0       0   0     0     0         0       0
                # 4      5        1     0         2    4.037721e-06   0.0     5             0       0   0     0     0         0       0
                # 5      6        1     0         2    5.389803e-06   0.0     6             0       0   0     0     0         0       0
                # 6      7        1     0         2    9.181562e+03   0.0     7             0       0   0     0     0         0       0
                # 7      8        1     0         2    2.184327e+04   0.0     8             0       0   0     0     0         0       0
                # 8      9        1     0         2    2.399637e+04   0.0     9             0       0   0     0     0         0       0
                # 9     10        1     0         2    2.903699e+04   0.0    10             0       0   0     0     0         0       0
                # 10    11        1     0         2    3.472390e+04   0.0    11             0       0   0     0     0         0       0
                # 11    12        1     0         2    4.077224e+04   0.0    12             0       0   0     0     0         0       0
                # 12    13        1     0         2    6.833794e+04   0.0    13             0       0   0     0     0         0       0
                # 13    14        1     0         2    8.648623e+04   0.0    14             0       0   0     0     0         0       0
                # 14    15        1     0         2    9.260647e+04   0.0    15             0       0   0     0     0         0       0
                # 15    16        1     0         2   6.013470e-154   0.0     1             0       0   0     0     0         0       0
                # 16    17        1     0         2   6.013471e-154   0.0     1             0       0   0     0     0         0       0
                # 17    18        1     0         2   6.013471e-154   0.0     1             0       0   0     0     0         0       0
                # 18    19        1     0         2   6.013471e-154   0.0     1             0       0   0     0     0         0       0
                # 19    20        1     0         2   6.013471e-154   0.0     1             0       0   0     0     0         0       0
                # 20    21        1     0         2   6.013471e-154   0.0     1             0       0   0     0     0         0       0
                # 21    22        1     0         2   6.013470e-154   0.0     2             0       0   0     0     0         0       0
                # 22    23        1     0         2   6.013471e-154   0.0     2             0       0   0     0     0         0       0
                # 23    24        1     0         2   6.013471e-154   0.0     2             0       0   0     0     0         0       0
                # 24    25        1     0         2   6.013471e-154   0.0     2             0       0   0     0     0         0       0
                # 25    26        1     0         2   6.013471e-154   0.0     2             0       0   0     0     0         0       0
                # 26    27        1     0         2   6.013471e-154   0.0     2             0       0   0     0     0         0       0
                # 27    28        1     0         2   6.013470e-154   0.0     3             0       0   0     0     0         0       0
                # 28    29        1     0         2   6.013471e-154   0.0     3             0       0   0     0     0         0       0
                # 29    30        1     0         2   6.013471e-154   0.0     3             0       0   0     0     0         0       0
                # 30    31        1     0         2   6.013471e-154   0.0     3             0       0   0     0     0         0       0
                # 31    32        1     0         2   6.013471e-154   0.0     3             0       0   0     0     0         0       0
                # 32    33        1     0         2   6.013471e-154   0.0     3             0       0   0     0     0         0       0
                # 33    34        1     0         2   6.013470e-154   0.0     4             0       0   0     0     0         0       0
                # 34    35        1     0         2   6.013471e-154   0.0     4             0       0   0     0     0         0       0
                # 35    36        1     0         2   6.013471e-154   0.0     4             0       0   0     0     0         0       0
                # 36    37        1     0         2   6.013471e-154   0.0     4             0       0   0     0     0         0       0
                # 37    38        1     0         2   6.013471e-154   0.0     4             0       0   0     0     0         0       0
                # 38    39        1     0         2   6.013471e-154   0.0     4             0       0   0     0     0         0       0
                # 39    40        1     0         2   6.013470e-154   0.0     5             0       0   0     0     0         0       0
                # 40    41        1     0         2   6.013471e-154   0.0     5             0       0   0     0     0         0       0
                # 41    42        1     0         2   6.013471e-154   0.0     5             0       0   0     0     0         0       0
                # 42    43        1     0         2   6.013471e-154   0.0     5             0       0   0     0     0         0       0
                # 43    44        1     0         2   6.013471e-154   0.0     5             0       0   0     0     0         0       0
                # 44    45        1     0         2   6.013471e-154   0.0     5             0       0   0     0     0         0       0
                # 45    46        1     0         2   6.013470e-154   0.0     6             0       0   0     0     0         0       0
                # 46    47        1     0         2   6.013471e-154   0.0     6             0       0   0     0     0         0       0
                # 47    48        1     0         2   6.013471e-154   0.0     6             0       0   0     0     0         0       0
                # 48    49        1     0         2   6.013471e-154   0.0     6             0       0   0     0     0         0       0
                # 49    50        1     0         2   6.013471e-154   0.0     6             0       0   0     0     0         0       0
                # 50    51        1     0         2   6.013471e-154   0.0     6             0       0   0     0     0         0       0
                # 51    52        1     0         2   6.013470e-154   0.0     7             0       0   0     0     0         0       0
                # 52    53        1     0         2   6.013471e-154   0.0     7             0       0   0     0     0         0       0
                # 53    54        1     0         2   6.013471e-154   0.0     7             0       0   0     0     0         0       0
                # 54    55        1     0         2   6.013471e-154   0.0     7             0       0   0     0     0         0       0
                # 55    56        1     0         2   6.013471e-154   0.0     7             0       0   0     0     0         0       0
                # 56    57        1     0         2   6.013471e-154   0.0     7             0       0   0     0     0         0       0
                # 57    58        1     0         2   6.013470e-154   0.0     8             0       0   0     0     0         0       0
                # 58    59        1     0         2   6.013471e-154   0.0     8             0       0   0     0     0         0       0
                # 59    60        1     0         2   6.013471e-154   0.0     8             0       0   0     0     0         0       0
                # 60    61        1     0         2   6.013471e-154   0.0     8             0       0   0     0     0         0       0
                # 61    62        1     0         2   6.013471e-154   0.0     8             0       0   0     0     0         0       0
                # 62    63        1     0         2   6.013471e-154   0.0     8             0       0   0     0     0         0       0
                # 63    64        1     0         2   6.013470e-154   0.0     9             0       0   0     0     0         0       0
                # 64    65        1     0         2   6.013471e-154   0.0     9             0       0   0     0     0         0       0
                # 65    66        1     0         2   6.013471e-154   0.0     9             0       0   0     0     0         0       0
                # 66    67        1     0         2   6.013471e-154   0.0     9             0       0   0     0     0         0       0
                # 67    68        1     0         2   6.013471e-154   0.0     9             0       0   0     0     0         0       0
                # 68    69        1     0         2   6.013471e-154   0.0     9             0       0   0     0     0         0       0
                # 69    70        1     0         2   6.013470e-154   0.0    10             0       0   0     0     0         0       0
                # 70    71        1     0         2   6.013471e-154   0.0    10             0       0   0     0     0         0       0
                # 71    72        1     0         2   6.013471e-154   0.0    10             0       0   0     0     0         0       0
                # 72    73        1     0         2   6.013471e-154   0.0    10             0       0   0     0     0         0       0
                # 73    74        1     0         2   6.013471e-154   0.0    10             0       0   0     0     0         0       0
                # 74    75        1     0         2   6.013471e-154   0.0    10             0       0   0     0     0         0       0
                # 75    76        1     0         2   6.013470e-154   0.0    11             0       0   0     0     0         0       0
                # 76    77        1     0         2   6.013471e-154   0.0    11             0       0   0     0     0         0       0
                # 77    78        1     0         2   6.013471e-154   0.0    11             0       0   0     0     0         0       0
                # 78    79        1     0         2   6.013471e-154   0.0    11             0       0   0     0     0         0       0
                # 79    80        1     0         2   6.013471e-154   0.0    11             0       0   0     0     0         0       0
                # 80    81        1     0         2   6.013471e-154   0.0    11             0       0   0     0     0         0       0
                # 81    82        1     0         2   6.013470e-154   0.0    12             0       0   0     0     0         0       0
                # 82    83        1     0         2   6.013471e-154   0.0    12             0       0   0     0     0         0       0
                # 83    84        1     0         2   6.013471e-154   0.0    12             0       0   0     0     0         0       0
                # 84    85        1     0         2   6.013471e-154   0.0    12             0       0   0     0     0         0       0
                # 85    86        1     0         2   6.013471e-154   0.0    12             0       0   0     0     0         0       0
                # 86    87        1     0         2   6.013471e-154   0.0    12             0       0   0     0     0         0       0
                # 87    88        1     0         2   6.013470e-154   0.0    13             0       0   0     0     0         0       0
                # 88    89        1     0         2   6.013471e-154   0.0    13             0       0   0     0     0         0       0
                # 89    90        1     0         2   6.013471e-154   0.0    13             0       0   0     0     0         0       0
                # 90    91        1     0         2   6.013471e-154   0.0    13             0       0   0     0     0         0       0
                # 91    92        1     0         2   6.013471e-154   0.0    13             0       0   0     0     0         0       0
                # 92    93        1     0         2   6.013471e-154   0.0    13             0       0   0     0     0         0       0
                # 93    94        1     0         2   6.013470e-154   0.0    14             0       0   0     0     0         0       0
                # 94    95        1     0         2   6.013471e-154   0.0    14             0       0   0     0     0         0       0
                # 95    96        1     0         2   6.013471e-154   0.0    14             0       0   0     0     0         0       0
                # 96    97        1     0         2   6.013471e-154   0.0    14             0       0   0     0     0         0       0
                # 97    98        1     0         2   6.013471e-154   0.0    14             0       0   0     0     0         0       0
                # 98    99        1     0         2   6.013471e-154   0.0    14             0       0   0     0     0         0       0
                # 99   100        1     0         2   6.013470e-154   0.0    15             0       0   0     0     0         0       0
                # 100  101        1     0         2   6.013471e-154   0.0    15             0       0   0     0     0         0       0
                # 101  102        1     0         2   6.013471e-154   0.0    15             0       0   0     0     0         0       0
                # 102  103        1     0         2   6.013471e-154   0.0    15             0       0   0     0     0         0       0
                # 103  104        1     0         2   6.013471e-154   0.0    15             0       0   0     0     0         0       0
                # 104  105        1     0         2   6.013471e-154   0.0    15             0       0   0     0     0         0       0

                #print(indexi)
                #print(dfi)
                #print('---------------------')
                idomain = np.searchsorted(DOMAINs, INDEX_DOMAIN)
                exists = idomain < len(INDEX_DOMAIN)
                if not np.all(exists):
                    if not np.any(exists):
                        continue
                    idomain = idomain[exists]
                    print('partial domain...')

                #idomain = []
            #for domain, position, length in zip(INDEX_DOMAIN, INDEX_POSITION, INDEX_LENGTH):
                #idomain = (DOMAIN_ID == domain)
                #mycase = domains_df.loc[idomain]
                position = INDEX_POSITION[idomain]
                length = INDEX_LENGTH[idomain]
                # ID=1; Analysis=0 -> eigenvalues
                #    ID  SUBCASE  STEP  ANALYSIS  TIME_FREQ_EIGR  EIGI  MODE  DESIGN_CYCLE  RANDOM  SE  AFPM  TRMC  INSTANCE  MODULE
                # 0   1        0     0         0    0.000000e+00   0.0     0             0       0   0     0     0         0       0
                # 1   2        1     0         2   -3.087735e-10   0.0     1             0       0   0     0     0         0       0
                # 2   3        1     0         2   -2.082743e-10   0.0     2             0       0   0     0     0         0       0
                # 3   4        1     0         2   -1.514309e-10   0.0     3             0       0   0     0     0         0       0
                subcase = dfi['SUBCASE'][idomain] # .values[0]

                analysis_code = dfi['ANALYSIS'][idomain]
                eigr = dfi['TIME_FREQ_EIGR'][idomain]# .values[0]
                eigi = dfi['EIGI'][idomain]
                mode = dfi['MODE'][idomain]# .values[0]

                #SEID Integer, default=0 Super element id of the data block
                #AFPMID Integer, default=0 Acoustic field point mesh id
                #TRIMID Trim id, default=0 Trim component id

                step = dfi['STEP'][idomain].values[0]
                design_cycle = dfi['DESIGN_CYCLE'][idomain].values[0]
                random = dfi['RANDOM'][idomain].values[0]
                se = dfi['SE'][idomain].values[0]
                afpm = dfi['AFPM'][idomain].values[0]
                trmc = dfi['TRMC'][idomain].values[0]
                inst = dfi['INSTANCE'][idomain].values[0]
                module = dfi['MODULE'][idomain].values[0]
                assert step in [0, 1], step
                #assert design_cycle == 0, design_cycle
                assert random == 0, random
                assert se == 0, se
                assert afpm == 0, afpm
                assert trmc == 0, trmc
                assert inst == 0, inst
                assert module == 0, module

                iresults = np.full(len(idomain), np.nan, dtype='int32')
                is_freq = np.abs(eigr).max() != 0 or np.abs(eigi).max() != 0
                is_modes = mode.max() != 0
                #step_type = ''
                for itime, subcasei, analysis_codei, modei, eigri, eigii, idomaini, positioni, lengthi in zip(count(), subcase, analysis_code, mode, eigr, eigi, idomain, position, length):
                    name = _get_eigenvector_name(basename, subcasei, analysis_codei, modei, eigri, eigii,
                                                 is_modes, is_freq)
                    if subcasei in subcases:
                        #print('found subcase')
                        pass
                    i0 = positioni
                    i1 = positioni + lengthi
                    iresults[itime] = iresult
                    _domain = DOMAIN[i0]

                    eid = EID[i0:i1]
                    energy = ENERGY[i0:i1]
                    percent = PCT[i0:i1]
                    density = DEN[i0:i1]
                    ident = IDENT[i0:i1] # TODO: what is this????
                    print(name, 'ident', list(np.unique(ident)))
                    results[iresult] = RealStrainEnergyOptistruct(
                        name, itime, iresult, iresults,
                        _domain, position, length,
                        eid,
                        energy, percent, density,
                        _domain, location='element')
                    iresult += 1
        #if name in ['BAR', 'QUAD4_COMP', 'TRIA3_COMP']:
            #pass
        else:
            raise NotImplementedError(name)
    return iresult

def load_stress_strain(basename_orig: str,
                       iresult: int,
                       results: dict[int, Function],

                       domains_df: pd.DataFrame,
                       element_group: h5py._hl.group.Group,
                       element_index: h5py._hl.group.Group,

                       ids: np.ndarray,
                       geom_model: BDF,
                       model: OP2,
                       is_stress: bool=True,
                       subcases: Optional[list[int]]=None) -> int:
    for ires, name in enumerate(element_group):
        group = element_group.get(name)
        index = element_index.get(name)
        if name in ['ELAS1', 'ELAS2', 'ELAS3',
                    'CONROD', 'ROD', 'TUBE', # 1D-basic
                    'BAR', 'BARS', 'BEAM', 'SHEAR',
                    'HEXA', 'PENTA', 'TETRA',
                    'QUAD4', 'QUAD4_COMP',
                    'QUAD8',
                    'QUADR_COMP',
                    'QUAD_CN',
                    'TRIA3', 'TRIA3_COMP',
                    'TRIA6',
                    'TRIAR_COMP',
                    ]:
            # real
            pass
        elif name in ['ELAS1_CPLX', 'ELAS2_CPLX', 'ELAS3_CPLX',
                      'BAR_CPLX', 'BEAM_CPLX',
                      'HEXA_CPLX', 'PENTA_CPLX', 'TETRA_CPLX',
                      'ROD_CPLX', 'CONROD_CPLX', 'TUBE_CPLX',
                      'SHEAR_CPLX',
                      'QUAD4_COMP_CPLX', 'QUAD_CN_CPLX',
                      'QUAD8_CPLX',
                      'QUADR_CPLX', 'QUADR_COMP_CPLX',
                      'TRIA3_CPLX', 'TRIA3_COMP_CPLX',
                      'TRIA6_CPLX',
                      'TRIAR_CPLX', 'TRIAR_COMP_CPLX',
                      ]:
            # complex
            pass
        else:
            raise NotImplementedError(name)
    return


def load_eigenvector(basename_orig: str,
                     iresult: int,
                     results: dict[int, Function],
                     domains_df: pd.DataFrame,
                     group: h5py._hl.dataset.Dataset,
                     index: h5py._hl.dataset.Dataset,
                     ids: np.ndarray,
                     geom_model: BDF,
                     model: OP2,
                     subcases: Optional[list[int]]=None) -> int:
    """
    <HDF5 dataset "EIGENVECTOR": shape (147,), type "|V64">
    Dataset:
    attrs  : <Attributes of HDF5 object at 2156696909096>
    chunks : (510,)
    compression : 'gzip'
    compression_opts : 1
    dims   : <Dimensions of HDF5 object at 2156696909096>
    dtype  : dtype([('ID', '<i8'), ('X', '<f8'), ('Y', '<f8'), ('Z', '<f8'), ('RX', '<f8'), ('RY', '<f8'), ('RZ', '<f8'), ('DOMAIN_ID', '<i8')])
    external : None
    file   : <HDF5 file "6+element-nastran-sol103.h5" (mode r)>
    fillvalue : (0, 0., 0., 0., 0., 0., 0., 0)
    fletcher32 : False
    id     : <h5py.h5d.DatasetID object at 0x000001F625273528>
    is_virtual : False
    maxshape : (None,)
    name   : '/NASTRAN/RESULT/NODAL/EIGENVECTOR'
    nbytes : 9408
    ndim   : 1
    parent : <HDF5 group "/NASTRAN/RESULT/NODAL" (1 members)>
    ref    : <HDF5 object reference>
    regionref : <h5py._hl.base._RegionProxy object at 0x000001F625247DC8>
    scaleoffset : None
    shape  : (147,)
    shuffle : True
    size   : 147
    -------------------------------------------------
    Dataset:
    attrs  : <Attributes of HDF5 object at 1418106158632>
    chunks : None
    compression : None
    compression_opts : None
    dims   : <Dimensions of HDF5 object at 1418106158632>
    dtype  : dtype([('DOMAIN_ID', '<i8'), ('POSITION', '<i8'), ('LENGTH', '<i8')])
    external : None
    file   : <HDF5 file "6+element-nastran-sol103.h5" (mode r)>
    fillvalue : (0, 0, 0)
    fletcher32 : False
    id     : <h5py.h5d.DatasetID object at 0x0000014A2DB6BE28>
    is_virtual : False
    maxshape : (3,)
    name   : '/INDEX/NASTRAN/RESULT/NODAL/EIGENVECTOR'
    nbytes : 72
    ndim   : 1
    parent : <HDF5 group "/INDEX/NASTRAN/RESULT/NODAL" (1 members)>
    ref    : <HDF5 object reference>
    regionref : <h5py._hl.base._RegionProxy object at 0x0000014A2DBB1108>
    scaleoffset : None
    shape  : (3,)
    shuffle : False
    size   : 3
    """
    basename_orig += ' (real)'
    # TODO: check real/imaginary or magnitude/phase
    #'ID', 'X', 'Y', 'Z', 'RX', 'RY', 'RZ', 'DOMAIN_ID',
    #'DOMAIN_ID', 'POSITION', 'LENGTH'
    INDEX_DOMAIN = index['DOMAIN_ID']
    INDEX_POSITION = index['POSITION']
    INDEX_LENGTH = index['LENGTH']

    NID = group['ID']
    names = group.dtype.names
    if 'X' in names:
        TX = group['X']
        TY = group['Y']
        TZ = group['Z']
        RX = group['RX']
        RY = group['RY']
        RZ = group['RZ']
        solver = 'msc'
    elif 'VALUE' in names:
        VALUE = group['VALUE']
        #(20214885, 6)
        solver = 'optistruct'
    else:
        raise RuntimeError(str(names))
    DOMAIN = group['DOMAIN_ID']
    # what are the subcases...
    #udomains = np.unique(DOMAIN)
    #all_data = np.stack([X, Y, Z, RX, RY, RZ], axis=1, out=None)
    DOMAIN_ID = domains_df['ID']
    grouped_df = break_domain_by_case(domains_df, INDEX_DOMAIN)
    subcases = geom_model.subcases
    for tuplei in grouped_df:
        indexi, dfi = tuplei
        DOMAINs = dfi['ID']
        #print(indexi)
        #print(dfi)
        #print('---------------------')
        idomain = np.searchsorted(DOMAINs, INDEX_DOMAIN)
        exists = idomain < len(INDEX_DOMAIN)
        if not np.all(exists):
            if np.any(exists):
                raise RuntimeError(idomain)
            continue

        #idomain = []
    #for domain, position, length in zip(INDEX_DOMAIN, INDEX_POSITION, INDEX_LENGTH):
        #idomain = (DOMAIN_ID == domain)
        #mycase = domains_df.loc[idomain]
        position = INDEX_POSITION[idomain]
        length = INDEX_LENGTH[idomain]
        # ID=1; Analysis=0 -> eigenvalues
        #    ID  SUBCASE  STEP  ANALYSIS  TIME_FREQ_EIGR  EIGI  MODE  DESIGN_CYCLE  RANDOM  SE  AFPM  TRMC  INSTANCE  MODULE
        # 0   1        0     0         0    0.000000e+00   0.0     0             0       0   0     0     0         0       0
        # 1   2        1     0         2   -3.087735e-10   0.0     1             0       0   0     0     0         0       0
        # 2   3        1     0         2   -2.082743e-10   0.0     2             0       0   0     0     0         0       0
        # 3   4        1     0         2   -1.514309e-10   0.0     3             0       0   0     0     0         0       0
        subcase = dfi['SUBCASE']# .values[0]
        analysis_code = dfi['ANALYSIS']
        eigr = dfi['TIME_FREQ_EIGR']# .values[0]
        eigi = dfi['EIGI']
        mode = dfi['MODE']# .values[0]

        #SEID Integer, default=0 Super element id of the data block
        #AFPMID Integer, default=0 Acoustic field point mesh id
        #TRIMID Trim id, default=0 Trim component id

        step = dfi['STEP'].values[0]
        design_cycle = dfi['DESIGN_CYCLE'].values[0]
        random = dfi['RANDOM'].values[0]
        se = dfi['SE'].values[0]
        afpm = dfi['AFPM'].values[0]
        trmc = dfi['TRMC'].values[0]
        inst = dfi['INSTANCE'].values[0]
        module = dfi['MODULE'].values[0]
        assert step in [0, 1], step
        #assert design_cycle == 0, design_cycle
        assert random == 0, random
        assert se == 0, se
        assert afpm == 0, afpm
        assert trmc == 0, trmc
        assert inst == 0, inst
        assert module == 0, module

        iresults = np.full(len(idomain), np.nan, dtype='int32')
        is_freq = np.abs(eigr).max() != 0 or np.abs(eigi).max() != 0
        is_modes = mode.max() != 0
        step_type = ''
        for itime, subcasei, analysis_codei, modei, eigri, eigii, idomaini, positioni, lengthi in zip(count(), subcase, analysis_code, mode, eigr, eigi, idomain, position, length):
            if subcasei in subcases:
                subcase = subcases[subcasei]
                if 'ANALYSIS' in subcase:
                    analysisi = subcase['ANALYSIS'][0]
                    step_type = f' ({analysisi})'
                else:
                    print('preload step?')
            else:
                print(f'what is subcase {subcasei}?')
            i0 = positioni
            i1 = positioni + lengthi
            #name = f'{basename}_subcase={subcasei:d}; mode={modei:d}; freq={eigri}'

            basename = f'{basename_orig}'
            name = _get_eigenvector_name(basename, subcasei, analysis_codei, modei, eigri, eigii,
                                         is_modes, is_freq)
                #name = f'{basename}: subcase={subcasei:d}; mode={modei:d}; freq={eigri}'
            print(name)
            #print(itime, domain, positioni, lengthi)

            # make it so we can determine the other "times"
            iresults[itime] = iresult
            _domain = DOMAIN[i0]
            __domain = _domain # # DOMAIN[i0:i1]
            _nids = NID[i0:i1]
            assert len(_nids) == len(ids)
            if solver == 'msc':
                results[iresult] = RealVectorTable(
                    name, itime, iresult, iresults,
                    _domain, position, length,
                    _nids,
                    TX[i0:i1], TY[i0:i1], TZ[i0:i1],
                    RX[i0:i1], RY[i0:i1], RZ[i0:i1],
                    _domain, 'node')
            elif solver == 'optistruct':
                results[iresult] = RealVectorTableOptistruct(
                    name, itime, iresult, iresults,
                    _domain, position, length,
                    _nids,
                    VALUE[i0:i1, :],
                    #VALUE[i0:i1, 3:],
                    #VALUE[i0:i1, :3],
                    _domain, location='node')
            else:
                raise NotImplementedError(solver)

            iresult += 1
    return iresult

def load_eigenvector_complex(basename: str,
                             iresult: int,
                             results: dict[int, Function],
                             domains_df: pd.DataFrame,
                             group: h5py._hl.dataset.Dataset,
                             index: h5py._hl.dataset.Dataset,
                             ids: np.ndarray,
                             geom_model: BDF,
                             model: OP2,
                             subcases: Optional[list[int]]=None) -> int:
    basename += ' (complex)'
    # TODO: check real/imaginary or magnitude/phase
    #'ID', 'X', 'Y', 'Z', 'RX', 'RY', 'RZ', 'DOMAIN_ID',
    #'DOMAIN_ID', 'POSITION', 'LENGTH'
    INDEX_DOMAIN = index['DOMAIN_ID']
    INDEX_POSITION = index['POSITION']
    INDEX_LENGTH = index['LENGTH']

    NID = group['ID']
    names = group.dtype.names
    if 'X' in names:
        TX = group['XR'] + group['XI'] * 1j
        TY = group['YR'] + group['YI'] * 1j
        TZ = group['ZR'] + group['ZI'] * 1j
        RX = group['RXR'] + group['RXI'] * 1j
        RY = group['RYR'] + group['RYI'] * 1j
        RZ = group['RZR'] + group['RZI'] * 1j
        solver = 'msc'
    #elif 'VALUE' in names:
        #asdf
    else:
        raise NotImplementedError(str(names))
    DOMAIN = group['DOMAIN_ID']
    # what are the subcases...
    #udomains = np.unique(DOMAIN)
    #all_data = np.stack([X, Y, Z, RX, RY, RZ], axis=1, out=None)
    DOMAIN_ID = domains_df['ID']
    grouped_df = break_domain_by_case(domains_df, INDEX_DOMAIN)
    for tuplei in grouped_df:
        indexi, dfi = tuplei
        DOMAINs = dfi['ID']
        #print(indexi)
        #print(dfi)
        #print('---------------------')
        idomain = np.searchsorted(DOMAINs, INDEX_DOMAIN)
        exists = idomain < len(INDEX_DOMAIN)
        if not np.all(exists):
            if np.any(exists):
                raise RuntimeError(idomain)
            continue

        #idomain = []
    #for domain, position, length in zip(INDEX_DOMAIN, INDEX_POSITION, INDEX_LENGTH):
        #idomain = (DOMAIN_ID == domain)
        #mycase = domains_df.loc[idomain]
        position = INDEX_POSITION[idomain]
        length = INDEX_LENGTH[idomain]
        # ID=1; Analysis=0 -> eigenvalues
        #    ID  SUBCASE  STEP  ANALYSIS  TIME_FREQ_EIGR  EIGI  MODE  DESIGN_CYCLE  RANDOM  SE  AFPM  TRMC  INSTANCE  MODULE
        # 0   1        0     0         0    0.000000e+00   0.0     0             0       0   0     0     0         0       0
        # 1   2        1     0         2   -3.087735e-10   0.0     1             0       0   0     0     0         0       0
        # 2   3        1     0         2   -2.082743e-10   0.0     2             0       0   0     0     0         0       0
        # 3   4        1     0         2   -1.514309e-10   0.0     3             0       0   0     0     0         0       0
        subcase = dfi['SUBCASE']
        analysis_code = dfi['ANALYSIS']
        eigr = dfi['TIME_FREQ_EIGR']
        eigi = dfi['EIGI']
        mode = dfi['MODE']

        step = dfi['STEP'].values[0]
        design_cycle = dfi['DESIGN_CYCLE'].values[0]
        random = dfi['RANDOM'].values[0]
        se = dfi['SE'].values[0]
        afpm = dfi['AFPM'].values[0]
        trmc = dfi['TRMC'].values[0]
        inst = dfi['INSTANCE'].values[0]
        module = dfi['MODULE'].values[0]
        assert step in [0, 1], step
        assert design_cycle == 0, design_cycle
        assert random == 0, random
        assert se == 0, se
        assert afpm == 0, afpm
        assert trmc == 0, trmc
        assert inst == 0, inst
        assert module == 0, module

        iresults = np.full(len(idomain), np.nan, dtype='int32')
        is_freq = np.abs(eigr).max() != 0 or np.abs(eigi).max() != 0
        for itime, subcasei, analysis_codei, modei, eigri, eigii, idomaini, positioni, lengthi in zip(
                count(), subcase, analysis_code, mode, eigr, eigi, idomain, position, length):
            i0 = positioni
            i1 = positioni + lengthi
            name = _get_name(basename, is_freq, subcasei, analysis_codei, modei, eigri, eigii)
            print(name)
            #print(itime, domain, positioni, lengthi)

            # make it so we can determine the other "times"
            iresults[itime] = iresult
            _domain = DOMAIN[i0]
            if solver == 'msc':
                results[iresult] = RealVectorTable(
                    name, itime, iresult, iresults,
                    _domain, position, length,
                    NID[i0:i1],
                    TX[i0:i1], TY[i0:i1], TZ[i0:i1],
                    RX[i0:i1], RY[i0:i1], RZ[i0:i1],
                    DOMAIN[i0:i1], location='node')
            elif solver == 'optistruct':
                results[iresult] = RealVectorTableOptistruct(
                    name, itime, iresult, iresults,
                    _domain, position, length,
                    NID[i0:i1],
                    VALUE[i0:i1, :],
                    #VALUE[i0:i1, 3:],
                    #VALUE[i0:i1, :3],
                    DOMAIN[i0:i1], location='node')
            else:
                raise NotImplementedError(solver)
            iresult += 1
    return iresult

def _get_eigenvector_name(basename, subcasei, analysis_codei, modei, eigri, eigii,
                          is_modes, is_freq):
    if modei == 0: # static?
        if is_freq:
            name = f'{basename}: subcase={subcasei:d}; freq={eigri:g}'
        else:
            name = f'{basename}: subcase={subcasei:d}'
    elif analysis_codei == 0:  # ???
        name = f'{basename}: subcase={subcasei:d}; static? mode={modei:d}; eigr={eigri:g}; eigi={eigii:g}'
    elif analysis_codei == 2:  # modes?
        if is_modes:
            name = f'{basename}: subcase={subcasei:d}; mode={modei:d}; freq={eigri:g}'
        else:
            raise NotImplementedError(analysis_codei)
    elif analysis_codei == 8:  # buckling (pre?/post?)
        # we left off eigr/eigi...(eigr=0; eigi!=0)
        name = f'{basename}: subcase={subcasei:d}; buckling mode={modei:d}; load_factor={eigri:g}; eigi={eigii:g}'
    elif analysis_codei == 9:  # loadstep?
        # we left off eigr/eigi...(eigr=0; eigi!=0)
        name = f'{basename}: subcase={subcasei:d}; loadstep? mode={modei:d}; eigi={eigii:g}'
    else:
        raise NotImplementedError(analysis_codei)
    return name

def _get_name(basename: str, is_freq: bool, subcasei, analysis_codei, modei, eigri, eigii):
    if modei == 0: # static
        if is_freq:
            name = f'{basename}: subcase={subcasei:d}; freq={eigri:g}; eigi={eigii:g}'
        else:
            name = f'{basename}: subcase={subcasei:d}'
            raise NotImplementedError(analysis_codei)
    elif analysis_codei == 9:
        # we left off eigr/eigi...(eigr=0; eigi!=0)
        name = f'{basename}: subcase={subcasei:d}; loadstep? mode={modei:d} eigi={eigii:g}'
        #raise NotImplementedError(analysis_codei)
    else:
        raise NotImplementedError(analysis_codei)
        #name = f'{basename}: subcase={subcasei:d}; mode={modei:d}; freq={eigri}'
    #name = f'Eigenvector_subcase={subcasei:d}; mode={modei:d}; freq={eigri:g} eigi={eigii:}'
    return name

class vtkNastranReader(vtk.vtkPolyDataAlgorithm):
    """
    References
    ----------
    https://kitware.github.io/paraview-docs/latest/python/paraview.simple.VisItNASTRANReader.html
    https://github.com/Kitware/VTK/blob/master/IO/Geometry/vtkSTLReader.cxx
    https://vtk.org/Wiki/VTK/Tutorials/New_Pipeline
    https://vtk.org/Wiki/ParaView/Examples/Plugins/Reader
    https://vtk.org/Wiki/VTK/Examples/Python/STLReader
    https://gitlab.kitware.com/paraview/paraview/-/blob/master/Examples/Plugins/PythonAlgorithm/PythonAlgorithmExamples.py
    https://gitlab.kitware.com/paraview/visitbridge/-/blob/master/databases/readers/NASTRAN/NASTRANPluginInfo.C
    https://vtk.org/Wiki/VisIt_Database_Bridge
    https://blog.kitware.com/developing-hdf5-readers-using-vtkpythonalgorithm/
    """
    def __init__(self):
        vtk.vtkPolyDataAlgorithm.__init__(self)
        self.filename = None
        self.model = pyNastranH5()
    def SetFileName(self, char_filename: str):
        self.filename = char_filename
        self.model.read_h5_nastran(self.filename)
    def GetFileName(self) -> str:
        return self.filename
    def IsFilePolyData(self) -> int:
        return 0
    def IsFileStructuredGrid(self) -> int:
        return 0
    def IsFileUnstructuredGrid(self) -> int:
        return 1
    def IsFileRectilinearGrid(self) -> int:
        return 0

    def RequestData(unused_request: vtk.vtkInformation,
                    unused_inputVector: vtk.vtkInformationVector,
                    unused_outputVector: vtk.vtkInformationVector) -> int:
        #vtkInformation* outInfo = outputVector->GetInformationObject(0)
        #vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()))

        ## All of the data in the first piece.
        #if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0):
            #return 0
        raise NotImplementedError()

    def GetSTLFileType(filename: str) -> int:
        ft = vtk.vtksys.SystemTools.DetectFileType(filename)
        #ft = vtksys::SystemTools::DetectFileType(filename)
        #vtksys::SystemTools::FileTypeEnum
        if ft == 0: # vtksys::SystemTools::FileTypeBinary:
            return VTK_BINARY
        elif ft == 1: # vtksys::SystemTools::FileTypeText:
            return VTK_ASCII
        elif ft2 == 2: # vtksys::SystemTools::FileTypeUnknown:
            vtkWarningMacro("File type not recognized; attempting binary")
            return VTK_BINARY
        else:
            raise RuntimeError("Case not handled, file type is %s" % ft)
        return VTK_BINARY  # should not happen

    def PrintSelf(file_obj, indent: vtk.vtkIndent):
        self.Superclass.PrintSelf(file_obj, indent)

        msg = indent + "Merging: " + ("On\n" if self.Merging else "Off\n")
        msg += indent + "ScalarTags: " + ("On\n" if self.ScalarTags else "Off\n")
        msg += indent + "Locator: "
        file_obj.write(msg)
        if self.Locator:
            asfd
            #self.Locator.PrintSelf(os << endl, indent.GetNextIndent());
        else:
            print('(none)\n')

    def ProcessRequest(request: vtk.vtkInformation,
                       inputVector: vtk.vtkInformationVector,
                       outputVector: vtk.vtkInformationVector) -> int:
        raise NotImplementedError()

    def ReadHeader(char_fname: str='') -> int:
        """Read the header of a vtk data file. More..."""
        raise NotImplementedError()

    def ReadCellData(ds: vtk.vtkDataSet, numCells: vtk.vtkIdTypeArray):
        """Read the cell data of a vtk data file. More..."""
        raise NotImplementedError()

    def ReadPointData(ds: vtk.vtkDataSet, numPts: vtk.vtkIdTypeArray) -> int:
        """Read the point data of a vtk data file. More..."""
        raise NotImplementedError()

    def ReadPointCoordinates(vtkPointSet_ps, numCells: vtk.vtkIdTypeArray) -> int:
        """Read point coordinates. More..."""
        raise NotImplementedError()

    def ReadPointCoordinates(g: vtk.vtkGraph, numPts: vtk.vtkIdTypeArray) -> int:
        """Read point coordinates. More..."""
        raise NotImplementedError()

    def ReadVertexData(g: vtk.vtkGraph, numVertices: vtk.vtkIdTypeArray) -> int:
        """Read the vertex data of a vtk data file. More..."""
        raise NotImplementedError()

    def ReadEdgeData(g: vtk.vtkGraph, numEdges: vtk.vtkIdTypeArray) -> int:
        """Read the edge data of a vtk data file. More..."""
        raise NotImplementedError()

    def ReadRowData(t: vtk.vtkTable, numEdges: vtk.vtkIdTypeArray) -> int:
        """Read the row data of a vtk data file. More..."""
        raise NotImplementedError()

    def ReadCells(vtkCellArray, cellArray) -> int:
        raise NotImplementedError()

    # ------------------------------------------------------------
    # Data Descriptors
    def CellArrayStatus(self):
        """This property lists which cell-centered arrays to read."""
        raise NotImplementedError()
    #def FileName(self):
    #"""The list of files to be read by the reader."""
    def MaterialStatus(self):
        """Select the materials to be loaded from the dataset, if any."""
        raise NotImplementedError()
    def MeshStatus(self):
        """Select the meshes to be loaded from the dataset, if any."""
        raise NotImplementedError()
    def PointArrayStatus(self):
        """This property lists which point-centered arrays to read."""
        raise NotImplementedError()
    def TimestepValues(self):
        """Available timestep values."""
        raise NotImplementedError()
    # ------------------------------------------------------------
    # Data Descriptors
    def CellData(self):
        """Returns cell data information"""
        raise NotImplementedError()
    def FieldData(self):
        """Returns field data information"""
        raise NotImplementedError()
    def PointData(self):
        """Returns point data information"""
        raise NotImplementedError()
    # ------------------------------------------------------------
    # Methods
    def FileNameChanged(self):
        """Called when the filename of a source proxy is changed."""
        raise NotImplementedError()
    def GetCellDataInformation(self):
        """Returns the associated cell data information."""
        raise NotImplementedError()
    def GetDataInformation(self, idx=None):
        """This method returns a DataInformation wrapper around a vtkPVDataInformation"""
        raise NotImplementedError()
    def GetFieldDataInformation(self):
        """Returns the associated cell data information."""
        raise NotImplementedError()
    def GetPointDataInformation(self):
        """Returns the associated point data information."""
        raise NotImplementedError()
    def UpdatePipeline(self, time=None):
        """This method updates the server-side VTK pipeline and the associated data information. Make sure to update a source to validate the output meta-data."""
        raise NotImplementedError()
    def UpdatePipelineInformation(self):
        """This method updates the meta-data of the server-side VTK pipeline and the associated information properties"""
        raise NotImplementedError()

def run():
    reader = vtk.vtkDataSetReader()
    reader.SetFileName("bunny-0.1.vtk")
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()

    data = reader.GetOutput()

    calc = vtk.vtkArrayCalculator()
    calc.SetInputData(data)

    calc.SetFunction("5")
    calc.SetResultArrayName("MyResults")
    calc.Update()

    # Gives: AttributeError: 'NoneType' object has no attribute 'GetPointData'
    #print(calc.GetPolyDataOutput().GetPointData().GetArray("MyResults").getValue(10))

    #writer = vtk.vtkUnstructuredGridWriter()
    #writer.SetInputData(data)
    #writer.SetFileName("Output.vtk")
    #writer.Write()

def run2(hdf5_filename: str):
    """
    https://github.com/Kitware/VTK/blob/master/IO/Geometry/vtkSTLReader.cxx
    https://vtk.org/Wiki/VTK/Examples/Python/STLReader
    """
    #vtkSTLReader
    reader = vtkNastranReader()
    print(reader)
    reader.SetFileName(hdf5_filename)

    mapper = vtk.vtkPolyDataMapper()
    #if vtk.VTK_MAJOR_VERSION <= 5:
        #mapper.SetInput(reader.GetOutput())
    #else:
    mapper.SetInputConnection(reader.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Create a rendering window and renderer
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)

    # Create a renderwindowinteractor
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    # Assign actor to the renderer
    ren.AddActor(actor)

    # Enable user interface interactor
    iren.Initialize()
    renWin.Render()
    iren.Start()

if __name__ == '__main__':  # pragma: no cover
    # tasks
    # - figure out nastran hdf5
    #   - done
    # - figure out how to make a vtk reader
    #   -
    # day 1: 3.5
    # day 2: 7
    # day 3: 3
    # day 4: 7
    # day 5: 3
    # 1/1: day 6: 6
    # ---------------------> 29.5 (for basic hdf5 and )
    # 1/3: 5  # gui works
    #
    #
    # tasks:
    #  - superelements?
    #  - solutions?
    #  - spoints
    #  - elements
    #  - is there enough RAM to store objects?
    #  - representative model
    #  - what are your goals once in Paraview/Paravis
    model = pyNastranH5()
    h5_filenames = [
        r'C:\NASA\m4\formats\git\pyNastran\models\sol_101_elements\buckling_solid_shell_bar.h5',
        r'C:\NASA\m4\formats\git\pyNastran\models\sol_101_elements\buckling2_solid_shell_bar.h5',
        r'C:\NASA\m4\formats\git\pyNastran\models\sol_101_elements\buckling_solid_shell_bar.h5',
        r'C:\NASA\m4\formats\git\pyNastran\models\aero\2_mode_flutter\0012_flutter.h5',
        r'C:\NASA\m4\formats\git\pyNastran\models\aero\aerobeam.h5',
        r'C:\NASA\m4\formats\git\pyNastran\models\aero\cpmopt.h5',
        r'C:\NASA\m4\formats\git\pyNastran\models\msc\6+element-nastran-sol103.h5',
        r'C:\NASA\m4\formats\git\pyNastran\models\msc\mode_echo.h5',
        r'C:\NASA\m4\formats\git\pyNastran\models\msc\ex1.h5',
        r'C:\NASA\m4\formats\git\pyNastran\models\sol_101_elements\static_solid_shell_bar.h5',

        r'C:\NASA\m4\formats\git\pyNastran\models\elements\static_elements.h5',
        r'C:\NASA\m4\formats\git\pyNastran\models\elements\modes_elements.h5',
        r'C:\NASA\m4\formats\git\pyNastran\models\elements\time_elements.h5',
        r'C:\NASA\m4\formats\git\pyNastran\models\elements\modes_complex_elements.h5',
        r'C:\NASA\m4\formats\git\pyNastran\models\elements\freq_elements.h5',
        r'C:\NASA\m4\formats\git\pyNastran\models\elements\freq_elements2.h5',
        ##r'C:\NASA\m4\formats\git\pyNastran\models\elements\loadstep_elements.h5',  # no nonlinear examples
        r'C:\NASA\m4\formats\git\pyNastran\models\elements\time_thermal_elements.h5',

        r'C:\NASA\m4\formats\git\pyNastran\models\bwb\bwb_saero_saved.h5',
        r'C:\NASA\m4\formats\git\pyNastran\models\sol200\d200obus.h5',
    ]
    for h5_filename in h5_filenames:
        model.read_h5_nastran(h5_filename)
    #model = pyNastranH5()
    #model.read_h5_nastran(h5_filename)
    #run2(h5_filename)
