"""
defines:
    bdf merge        (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME]\n'
    bdf equivalence  IN_BDF_FILENAME EQ_TOL\n'
    bdf renumber     IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'
    bdf mirror       IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--plane PLANE] [--tol TOL]\n'
    bdf export_mcids IN_BDF_FILENAME [-o OUT_GEOM_FILENAME]\n'
    bdf solid_dof    IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--spc SPC]\n'
    bdf split_cbars_by_pin_flags IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'
    bdf flutter UNITS [-o OUT_BDF_FILENAME]

"""
from __future__ import annotations
import sys
from typing import TYPE_CHECKING

import pyNastran
from pyNastran.bdf.mesh_utils.bdf_renumber_exclude import bdf_renumber_exclude
from pyNastran.bdf.mesh_utils.pierce_shells import pierce_shell_model

# test imports
# if something is imported and tested, it should be removed from here
from pyNastran.bdf.mesh_utils.shift import update_nodes
from pyNastran.bdf.mesh_utils.mirror_mesh import write_bdf_symmetric
#from pyNastran.bdf.mesh_utils.collapse_bad_quads import convert_bad_quads_to_tris
from pyNastran.bdf.mesh_utils.delete_bad_elements import delete_bad_shells, get_bad_shells
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused
from pyNastran.bdf.mesh_utils.free_faces import write_skin_solid_faces
from pyNastran.bdf.mesh_utils.get_oml import get_oml_eids

from pyNastran.bdf.mesh_utils.run_jobs import cmd_line_run_jobs
from pyNastran.bdf.mesh_utils.host_jobs import cmd_line_host_jobs
from .solid_dof import solid_dof
from .flip_shell_normals import flip_shell_normals

from .cmd_line.bdf_cmd_line import cmd_line
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


if __name__ == '__main__':  # pragma: no cover
    # for the exe, we pass all the args, but we hack them to have the bdf prefix
    from copy import deepcopy
    argv_root = deepcopy(sys.argv)
    argv_root[0] = 'bdf'
    cmd_line(argv=argv_root)
