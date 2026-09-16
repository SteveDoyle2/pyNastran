import os
from typing import TextIO
import numpy as np
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.field_writer_8 import print_card_8


def cmd_line_caero_to_wkk():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('bdf_filename', type=str)
    # parser.add_argument('output_filename', type=str, default='')
    dmi_dmik_group = parser.add_mutually_exclusive_group()
    dmi_dmik_group.add_argument('--dmi', action='store_true')
    dmi_dmik_group.add_argument('--dmik', action='store_true')
    args = parser.parse_args()
    print('args', args)

    dmi_type = 'DMIK'
    if args.dmi:
        dmi_type = 'DMI'
    dmi_type_lower = dmi_type.lower()

    bdf_filename = args.bdf_filename
    if 1: # args.output_filename != '':
        base, ext = os.path.splitext(bdf_filename)
        # bdf_filename_out = base + '_' + 'out' + ext
        bdf_filename_out = f'{dmi_type_lower}_wkk.blk'
    else:  # pragma: no cover
        bdf_filename_out = args.output_filename
    caero_to_wkk(bdf_filename, bdf_filename_out=bdf_filename_out,
                 dmi_type=dmi_type)

def caero_to_wkk(bdf_filename: str, bdf_filename_out: str='',
                 dmi_type: str='DMIK'):
    assert dmi_type in {'DMI', 'DMIK'}, dmi_type
    skip_cards = [
        'CQUAD4', 'CTRIA3', 'CQUAD8', 'CTRIA6',
        'CTETRA', 'CHEXA', 'CPYRAM', 'CPENTA',
        'CBUSH', 'CONM1', 'CONM2',
        'RBE2', 'RBE3',
        'PSHELL', 'PCOMP', 'PCOMPG', 'PSOLID', 'PBUSH',
        'SPC', 'SPC1', 'MPC',
    ]
    model = read_bdf(bdf_filename, skip_cards=skip_cards)
    tin = 1 # 'float32'
    tout = 0 # same as tin
    # form = 'column'
    form = 1  # square
    name = 'WKK'
    GCj_list = []
    GCi_list = []
    reals_list = []
    # DMI          WKK       0       3       1       1            6704       1
    # DMI          WKK       1       1    1.00
    panels_initial = list(model.caeros)
    for eid, caero in model.caeros.items():
        ids = caero.box_ids.ravel()
        list_ids = ids.tolist()
        assert isinstance(list_ids, list), list_ids
        # model.add_dmij(name, form, tin, tout, )
        # print(list_ids)
        npaneli = len(list_ids)
        GCj = list(list_ids)
        GCi = list(list_ids)
        reals = np.ones(npaneli).tolist()
        reals_list.append(reals)
        GCi_list.append(GCi)
        GCj_list.append(GCj)

    # print('GCj_list', GCj_list)
    # print('GCi_list', GCi_list)
    GCj = np.hstack(GCj_list)
    panels = GCj
    panels.sort()
    npanel = len(panels)
    if 0:  # pragma: no cover
        panels2 = np.column_stack([panels, panels])
        three_five = np.column_stack([
            np.full(npanel, 3, dtype='int32'),
            np.full(npanel, 5, dtype='int32'),
        ])
        GCj = np.column_stack([panels2, three_five])
        GCi = np.column_stack([panels2, three_five])
    # GCi = np.hstack(GCi_list)
    Real = np.hstack(reals_list)

    # def add_dmik(self, name: str, ifo: int,
    #              tin: int, ncols: int,
    #              GCj: np.ndarray, GCi: np.ndarray,
    #              Real: np.ndarray, Complex=None,
    #              tout: int=0, polar: int=0,
    #              comment: str='') -> DMIK:

    ncols = len(Real)
    # model.add_dmik(
    #     name, form, tin, ncols,
    #     GCj, GCi, Real, Complex=None,
    #     tout=tout, polar=0,
    #     comment='')

    # def add_dmi(self, name: str, form: int | str,
    #             tin: int | str,
    #             nrows: int, ncols: int,
    #             GCj: np.ndarray, GCi: np.ndarray,
    #             Real: np.ndarray, Complex=None,
    #             tout: int | str=0,
    #             comment: str='') -> DMI:
    nrows = ncols
    model_out = BDF(log=model.log)

    if 0:  # pragma: no cover
        model_out.add_dmi(
            name, form, tin, nrows, ncols,
            panels, panels, Real, Complex=None,
            tout=0,
            comment='')
        if bdf_filename_out:
            model_out.write_bdf(bdf_filename_out)
    else:
        model.log.info(f'writing {bdf_filename_out}')
        with open(bdf_filename_out, 'w') as bdf_file:
            bdf_file.write(f'$ form = {form} (square)\n')
            if dmi_type == 'DMI':
                write_dmi_wkk(
                    bdf_file,
                    model, form, tin, tout,
                    npanel, panels_initial, panels)
            elif dmi_type == 'DMIK':
                write_dmik_wkk(
                    bdf_file,
                    model, form, tin, tout,
                    npanel, panels_initial, panels)
            else:  # pragma: no cover
                raise NotImplementedError(dmi_type)


def write_dmi_wkk(bdf_file: TextIO,
                  model: BDF,
                  form: int, tin: int, tout: int, npanel: int,
                  panels_initial: list[int], panels: np.ndarray):
    name = '%-8s' % 'WKK'
    npanel2 = npanel * 2
    real = 1.0

    card0 = ['$ DMI', 'name', '0', 'form', 'tin', 'tout', '', 'npanel2', 'npanel2',]
    bdf_file.write(print_card_8(card0))
    card = ['DMI', name, '0', form, tin, tout, '', npanel2, npanel2,]
    bdf_file.write(print_card_8(card))
    bdf_file.write('$ ---------- force ----------\n')
    for ipanel, panel in enumerate(panels):
        if panel in panels_initial:
            bdf_file.write(f'$ force {panel}\n')
            caero_card = str(model.caeros[panel]).split('\n')
            bdf_file.write('$' + '\n$ '.join(caero_card) + '\n')

        # panels start at 1, but ipanel starts at 0
        j = ipanel * 2 + 1
        i = j
        # DMI          WKK       1       1    1.00
        # DMI          WKK       2       2    1.00
        card = ['DMI', name, j, i, real]
        bdf_file.write(print_card_8(card))

    bdf_file.write('$ ---------- moment ----------\n')
    for ipanel, panel in enumerate(panels):
        if panel in panels_initial:
            bdf_file.write(f'$ force {panel}\n')
        # panels start at 1, but ipanel starts at 0
        # offset by 1 for moment
        j = ipanel * 2 + 2
        i = j
        # bdf_file.write(f'$ moment {ipanel}\n')
        card = ['DMI', name, j, i, real]
        bdf_file.write(print_card_8(card))


def write_dmik_wkk(bdf_file: TextIO,
                   model: BDF,
                   form: int, tin: int, tout: int, npanel: int,
                   panels_initial: list[int], panels: np.ndarray):
    name = '%-8s' % 'WKK'
    npanel2 = npanel * 2
    real = 1.0
    assert form == 1, (form, 'square')
    polar = ''

    card0 = ['$ DMIK', 'name', '0', 'form', 'tin', 'tout', 'polar', '', 'ncol']
    bdf_file.write(print_card_8(card0))
    card = ['DMIK', name, '0', form, tin, tout, polar, '', npanel2,]
    bdf_file.write(print_card_8(card))
    bdf_file.write('$ ---------- force ----------\n')
    dof = 3
    for ipanel, panel in enumerate(panels):
        if panel in panels_initial:
            bdf_file.write(f'$ force {panel}\n')
            caero_card = str(model.caeros[panel]).split('\n')
            bdf_file.write('$' + '\n$ '.join(caero_card) + '\n')

        # DMI          WKK       1       1    1.00
        # DMI          WKK       2       2    1.00
        card = ['DMIK', name, panel, dof, '', panel, dof, real]
        bdf_file.write(print_card_8(card))

    bdf_file.write('$ ---------- moment ----------\n')
    dof = 5
    for ipanel, panel in enumerate(panels):
        if panel in panels_initial:
            bdf_file.write(f'$ moment {panel}\n')
        card = ['DMIK', name, panel, dof, '', panel, dof, real]
        bdf_file.write(print_card_8(card))

if __name__ == '__main__':  # pragma: no cover
    cmd_line_caero_to_wkk()
