from pyNastran.dev.bdf_vectorized3.bdf import BDF


def convert(model: BDF, units_to: list[str], units: list[str]) -> None:
    xyz_scale = 2.0
    mass_scale = 3.0
    time_scale = 1.0
    gravity_scale = 1.0
    temperature_scale = 1.0

    area_scale = xyz_scale ** 2
    volume_scale = xyz_scale ** 3
    force_scale = xyz_scale ** 2
    velocity_scale = xyz_scale / time_scale
    accel_scale = xyz_scale / time_scale ** 2
    area_inertia_scale = xyz_scale ** 4
    mass_inertia_scale = mass_scale * xyz_scale ** 2

    density_scale = force_scale / volume_scale
    pressure_scale = force_scale / area_scale
    scales_dict = {
        'xyz_scale': xyz_scale,
        'area_scale': area_scale,
        'volume_scale': volume_scale,

        'time_scale': time_scale,
        'gravity_scale': gravity_scale,

        'velocity_scale': velocity_scale,
        'accel_scale': accel_scale,
        'force_scale': force_scale,
        'pressure_scale': pressure_scale,

        'density_scale': density_scale,
        'mass_inertia_scale': mass_inertia_scale,
        'area_inertia_scale': area_inertia_scale,
        'nsm_per_length_scale': mass_scale / xyz_scale,
        'nsm_per_area_scale': mass_scale / area_scale,
    }
    CARDS_TO_SKIP = {
        'SPOINT', 'EPOINT', 'CTUBE', 'CROD',
        'SET1', 'ASET', 'BSET', 'CSET', 'OMIT', 'USET', 'SUPORT',
        'PLOTEL', 'SPCADD', 'MPCADD', 'DCONADD',
        'CQUAD', 'CTETRA', 'CPYRAM', 'CPENTA', 'CHEXA',
    }
    SUPPORTED_CARDS = {
        'GRID',  'POINT', 'CONM1', 'CONM2',
        'CBAR', 'PBAR', 'PBARL',
        'CBEAM', 'PBEAM', 'PBEAML',
        'PROD', 'PTUBE', 'CONROD',
        'RBE1', 'RBE2', 'RBE3', 'RBAR', 'RBAR1', 'RROD',
        'CTRIA3', 'CQUAD4', 'CTRIAR', 'CQUADR', 'CTRIA6', 'CQUAD8',
        'PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE',
    }
    cards = [card for card in model._cards_to_setup
             if card.n and card.type not in CARDS_TO_SKIP]
    skipped_cards = []
    for card in cards:
        if card.type in SUPPORTED_CARDS:
            card.convert(**scales_dict)
        elif hasattr(card, 'convert'):  # pragma: no cover
            card.convert(**scales_dict)
            raise RuntimeError(card.type)
        else:
            skipped_cards.append(card.type)
    if skipped_cards:
        model.log.debug('cant convert {skipped_cards}')

