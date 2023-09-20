from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
from numpy.lib import recfunctions as rfn  # random numpy package to merge arrays...
try:
    from tables import open_file, Group, Node, File
except ImportError:
    print('pytables was not found; no h5 support.  Run ">>> pip install tables"\n'
          'Do you have h5py installed?  That can cause conflicts.')
    raise

from pyNastran.dev.bdf_vectorized3.bdf_interface.h5_pytables.h5_nodes import load_h5_node
from pyNastran.dev.bdf_vectorized3.bdf_interface.h5_pytables.h5_coords import load_h5_coord
#from pyNastran.dev.bdf_vectorized3.bdf_interface.h5_pytables.h5_parameter import load_h5_parameter
#from pyNastran.dev.bdf_vectorized3.bdf_interface.h5_pytables.h5_partition import load_h5_partition
#from pyNastran.dev.bdf_vectorized3.bdf_interface.h5_pytables.h5_properties import load_h5_property
from pyNastran.dev.bdf_vectorized3.bdf_interface.h5_pytables.h5_materials import load_h5_material
from pyNastran.dev.bdf_vectorized3.bdf_interface.h5_pytables.h5_elements import load_h5_element
#from pyNastran.dev.bdf_vectorized3.bdf_interface.h5_pytables.h5_dynamic import load_h5_dynamic
from pyNastran.dev.bdf_vectorized3.bdf_interface.h5_pytables.utils import get_group_name
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


def read_h5_geometry(model: BDF, h5_filename: str,
                     root_path: str='/'):
    with open_file(h5_filename, mode="r", title="", root_uep="/", filters=None) as h5_file:
        read_geometry_from_h5(model, h5_file, '/NASTRAN/INPUT/')


def read_geometry_from_h5(model: BDF, h5_file: File, input_path: str):
    inputi = h5_file.get_node(input_path)
    log = model.log
    for h5_node_ in inputi._f_iter_nodes():
        #print(f'node {h5_node_}')
        #h5_node_._c_classid = 'GROUP'

        #if isinstance(h5_node_, np.ndarray):
            #h5_node = h5_node_.read()
            #name = h5_node_.name
        if isinstance(h5_node_, Group):
            #name = h5_node.name
            ## TODO: janky way to get the name
            name = get_group_name(h5_node_)

            if name in {'DYNAMIC', }:
                log.warning(f' - name {name!r}')
                #load_h5_dynamic(model, h5_node_)
                continue
            elif name == 'ELEMENT':
                load_h5_element(model, h5_node_)
                continue
            elif name == 'NODE':
                load_h5_node(model, h5_node_)
                continue
            elif name == 'PROPERTY':
                continue
                load_h5_property(model, h5_node_)
                continue
            elif name == 'MATERIAL':
                load_h5_material(model, h5_node_)
                continue
            elif name == 'COORDINATE_SYSTEM':
                load_h5_coord(model, h5_node_)
                continue
            elif name == 'CONSTRAINT':
                log.warning(f'skipping h5group name={name}')
                continue
            elif name == 'DESIGN':
                log.warning(f'skipping h5group name={name}')
                continue
            elif name == 'LOAD':
                log.warning(f'skipping h5group name={name}')
                continue
            #elif name == 'PARAMETER':
                #load_h5_parameter(model, h5_node_)
                #continue
            #elif name == 'PARTITION':
                #load_h5_partition(model, h5_node_)
                #continue
            log.warning(f'skipping h5group name={name!r}')
            continue
        elif isinstance(h5_node_, Node):
            h5_node = h5_node_.read()
            assert isinstance(h5_node, np.ndarray), h5_node
            name = h5_node_.name
        else:
            raise NotImplementedError(h5_node_)
        log.warning(f' - name {name!r}')
        # DOMAINS

    #  this only gets the folders, but not the files/nodes/tables/arrays
    #  (these terms all mean the same thing)
    #
    # the problem is NODE, which is a folder, doesn't show up in this...
    if 0:
        for (input_group_key, input_group) in inputi._v_children.items():  #  good
        #for node in inputi.walk_groups():
            print(f'child {input_group_key!r}')
            #if input_group_key == 'NODE':
                #x = 1
            if input_group_key == 'ELEMENT':
                #_load_h5_element(model, input_group) # good
                pass

            #elif input_group_key == 'MATERIAL':

            #elif input_group_key == 'PROPERTY':
            #elif input_group_key == 'DYNAMIC':
                #pass
            #elif input_group_key == 'PARAMETER':
                #pass
            #elif input_group_key == 'DOMAINS':
                #pass
            else:
                print(f'INPUT input_group_key={input_group_key!r}')
            #'PARAMETER' (Group), 'PROPERTY' (Group), 'DOMAINS' (Table)]
    model.coord.setup()
    print('finished loading h5 file')
