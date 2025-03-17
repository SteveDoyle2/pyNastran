"""
myname	htail act. Stiffness	vtail act stiff	mass	wing skin	htail skin
id	99	100	101	102,103,104	201
type	PELAS, PROP	PBUSH	CONM2	PCOMP	PCOMP
layer	N/A	N/A	-1	ALL	1-3,5-6
variable name	K	K1	M	T	T
value	1.00E+06	1.00E+07	100	0.2	0.1

"""
import copy
from pathlib import Path
import numpy as np
from pyNastran.bdf.bdf import BDF, PBUSH

properties_list = [
    'PROP', 'PBUSH', 'PSHELL', 'PCOMP', 'PCOMPG',
]
masses_list = ['CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4']

elements_list = [
    'CQUAD4', 'CTRIA3',
    'CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM'
    'CBUSH', 'CBAR', 'CBEAM', 'CBEND',
]
pbush_field_map = {idi: field for field, idi in PBUSH.pname_map.items()}


def cmd_line_setup_trade(argv=None, quiet: bool=False) -> None:
    base_model = BDF()
    base_model.add_pbush(100, k=[1., 2., 3., 4., 5., 6.])
    base_model.add_pbush(200, k=[10., 20., 30., 40., 50., 60.])
    #base_model.validate()

    meta_data = [
        ('stiff1', 'PBUSH', 100, 'K1,K2', 1000.),  # nan is blank
        ('stiff2', 'PBUSH', 200, 'K2', 2000.),
    ]
    data = np.array([
        [101.,   205.],
        [np.nan, 508.],  # nan is default, not 0.0
        [np.nan, np.nan],  # nan is default, not 0.0
    ])
    row0 = data[0]
    base_cards = get_base_cards(base_model, meta_data, row0)
    for irow, row in enumerate(data):
        build_model(base_cards, meta_data, row, irow)


def get_base_cards(base_model: BDF,
                   meta_data,
                   row: np.ndarray) -> list:
    base_cards = []
    assert len(meta_data) == len(row)
    for meta_datai in meta_data:
        param_type = meta_datai[1]
        pid_eid = meta_datai[2]
        if param_type in properties_list:
            param_group = 'properties'
            card = base_model.properties[pid_eid]
            del base_model.properties[pid_eid]
        elif param_type in masses_list:
            param_group = 'masses'
            card = base_model.masses[pid_eid]
            del base_model.masses[pid_eid]
        elif param_type in elements_list:
            param_group = 'elements'
            card = base_model.elements[pid_eid]
            del base_model.elements[pid_eid]
        else:  # pragma: no cover
            raise RuntimeError(f'{param_type} is not supported')
        #print(f'card:\n{str(card)}')
        response_types = meta_datai[3].split(',')
        value = meta_datai[4]
        base_cards.append((param_group, card, response_types, value))
    return base_cards

def _validate_model(meta_data, row: np.ndarray) -> None:
    # validation
    for meta_datai, value in zip(meta_data, row):
        name = meta_datai[0]
        param_type = meta_datai[1]
        pid_eid = meta_datai[2]
        response_types = meta_datai[3].split(',')
        for response_type in response_types:
            if param_type == 'PBUSH':
                assert response_type in pbush_field_map, response_type
            else:  # pragma: no cover
                raise RuntimeError(f'{param_type}-{response_type} is not a PBUSH')

def build_model(base_cards: list,
                meta_data: list,
                row: np.ndarray,
                irow: int) -> None:
    dirname = Path('.')
    size = 8
    assert len(meta_data) == len(row)
    #model = BDF()

    _validate_model(meta_data, row)
    #-------------------------------------------------------------------
    # remove cards to overwrite from model and save them in a backup file
    # drop the ENDDATA

    #-------------------------------------------------------------------
    # create a bunch of files to replace
    group_cards = []
    base_filename = dirname / f'run_{irow}.blk'

    print(f'{base_filename}')
    with open(base_filename, 'w') as blk_file:
        for base_card_data, datai in zip(base_cards, row):
            (param_group, card_bkp, response_types, value) = base_card_data
            #print('response_types = ', response_types)
            #print(card_bkp)
            card = copy.deepcopy(card_bkp)
            for response_type in response_types:
                card.update_by_pname_fid(response_type, value)
            print(card.rstrip())
            blk_file.write(card.write_card(size=size))
    print('-----------------')
def main():
    cmd_line_setup_trade()


if __name__ == '__main__':
    main()
