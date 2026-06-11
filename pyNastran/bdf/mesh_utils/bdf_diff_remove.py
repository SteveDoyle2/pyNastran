import os
from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.mesh_utils.cmd_line.utils_bdf import read_lax_obj
from .bdf_diff import scalar_obj_keys, dict_cards


def apply_diff(bdf_filename: PathLike,
               remove_filename: PathLike, encoding=None) -> None:
    obj_filename = os.path.splitext(bdf_filename)[0] + ".obj"
    is_obj = True
    model = read_lax_obj(bdf_filename, obj_filename, is_obj, xref=False, save_file_structure=True)
    log = model.log
    # else:
    #     model = BDF()
    #     model.read_bdf(
    #         bdf_filename, xref=False, punch=False, encoding=encoding, save_file_structure=True)

    model_to_remove = BDF()
    model_to_remove.read_bdf(
        remove_filename, xref=False, punch=False, encoding=encoding, save_file_structure=True)

    # dict_attrs = ['elements', 'properties', 'masses', 'properties_mass', 'rigid_elements']
    ifiles_used = set([])
    for attr in dict_cards:
        model_attr = getattr(model, attr)  # model.elements
        remove_attr = getattr(model_to_remove, attr)
        for eid, elem_remove in remove_attr.items():
            elem = model_attr[eid]
            ifiles_used.add(elem.ifile)
            del model_attr[eid]

    # out_files_map : dict[source_bdf, out_bdf]
    #     source_bdf : str
    #         the name of the original bdf
    #     out_bdf : str
    #         the name of the output bdf
    out_files_map = {}
    print(f'active_filenames = {model.active_filenames}')
    out_files_map[model.active_filenames[0]] = bdf_filename[:-4] + "_NEW" + bdf_filename[-4:]
    for ifile, include_filenames in model.include_filenames.items():
        log.debug(f'ifile={ifile} -> include_filenames={include_filenames}')
        if ifile not in ifiles_used:
            continue
        for include_filename in include_filenames:
            base, ext = os.path.splitext(include_filename)
            new_filename = base + "_NEW" + ext
            out_files_map[include_filename] = new_filename

    relative_dirname = ''
    model.write_bdfs(
        out_files_map,
        relative_dirname=relative_dirname,
        encoding=None,
        size=8, is_double=False,
        enddata=None, close=True,
        is_windows=None)
