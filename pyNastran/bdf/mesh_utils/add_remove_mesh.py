from typing import Optional
import os
from pathlib import Path
from pyNastran.utils import PathLike, print_bad_path
from pyNastran.bdf.bdf import BDF

def add_remove_mesh(bdf_filename: PathLike,
                    bdf_add_remove_filenames: list[tuple[str, PathLike]]=None,) -> BDF:
    # check that files exist
    assert os.path.exists(bdf_filename), print_bad_path(bdf_filename)
    for add_remove, bdf_filenamei in bdf_add_remove_filenames:
        assert add_remove in {'add', 'remove'}, (add_remove, bdf_filenamei)
        assert os.path.exists(bdf_filenamei), print_bad_path(bdf_filenamei)

    model_root = BDF()
    model_root._parse = False
    model_root.read_bdf(bdf_filename, xref=False)
    slot_to_type_map = model_root._slot_to_type_map
    slot_names = list(slot_to_type_map.keys())

    string_groups = [
        # dict[str]
        'params', 'dmig', 'dmi', 'dmiax', 'dti',
    ]

    dict_list_groups = [
        # dict[spc_id][list_index]
        'loads',
        'spcs', 'spcoffs',
        'mpcs',
    ]

    group_names = string_groups + dict_list_groups
    for group_name in group_names:
        # temp to check special cards
        assert group_name in slot_names, (group_name, slot_names)

    for add_remove, bdf_filenamei in bdf_add_remove_filenames:
        modeli = BDF()
        modeli._parse = False
        modeli.read_bdf(bdf_filenamei, xref=False)
        log = modeli.log

        #print(modeli.elements)
        for slot_name in slot_names:
            if not hasattr(modeli, slot_name):  # panlsts
                continue

            root_group = getattr(model_root, slot_name)
            modeli_group = getattr(modeli, slot_name)
            if root_group is None and modeli_group is None:
                continue
            if isinstance(root_group, (list, dict)):
                if len(root_group) == 0 and len(modeli_group) == 0:
                    continue
            log.debug(slot_name)

            if add_remove == 'add':
                if slot_name in dict_list_groups:
                    # spcs, mpcs
                    log.warning(f'adding {slot_name}')
                    # for idi, card_lines in modeli_group.items():
                    #     del root_group[idi]

                elif slot_name in string_groups:
                    for idi, card_lines in modeli_group.items():
                        log.debug(slot_name, card_lines[0])
                        root_group[idi] = card_lines
                else:
                    # CQUAD4, CONM2
                    for idi, card_lines in modeli_group.items():
                        log.debug(slot_name, card_lines[0])
                        root_group[idi] = card_lines
            else:
                # remove
                if slot_name in dict_list_groups:
                    # spcs, mpcs
                    log.warning(f'removing {slot_name}')
                    for idi, card_lines in modeli_group.items():
                        del root_group[idi]

                elif slot_name in string_groups:
                    for idi, card_lines in modeli_group.items():
                        del root_group[idi]
                else:
                    # CQUAD4, CONM2
                    for idi, card_lines in modeli_group.items():
                        if not isinstance(card_lines, list):  # CORD2R
                            continue
                        print(card_lines)
                        msgi = str((slot_name, card_lines[0]))
                        log.debug(msgi)
                        del root_group[idi]
    #for

def main():
    import pyNastran
    pkg_path = Path(pyNastran.__path__[0])
    print(pkg_path)
    model_path = pkg_path / '..' / 'models'
    assert model_path.exists(), print_bad_path(model_path)
    bdf_filename = model_path / 'solid_bending' / 'solid_bending.bdf'
    add_bdf_filenames = []
    add_remove_bdf_filenames = [
        ('remove', model_path / 'solid_bending' / 'solid_bending.bdf'),
    ]
    add_remove_mesh(
        bdf_filename,
        add_remove_bdf_filenames,
    )


if __name__ == '__main__':  # pragma: no cover
    main()
