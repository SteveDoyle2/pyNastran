from pyNastran.bdf.bdf import BDF


def flip_shell_normals(model: BDF, zero_zoffset: float) -> None:
    log = model.log
    skip_elements = {
        # nothing to convert (verified)
        'CCONEAX',
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
        'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',
        'CVISC', 'CBUSH', 'CBUSH1D', 'CBUSH2D',
        'CROD', 'CTUBE',
        'CBAR', 'CBEAM',
        'CSHEAR', 'CQUADX', 'CTRIAX', 'CTRIAX6',
        'CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM',
        'CRAC2D', 'CRAC3D',
        'CHBDYG', 'CHBDYE', 'CHBDYP',

        # TODO: NX-verify
        'CTRAX3', 'CTRAX6',
        'CPLSTN3', 'CPLSTN6', 'CPLSTN4', 'CPLSTN8',
        'CQUADX4', 'CQUADX8',

        # acoustic
        'CHACAB',
    }

    shells = {
        'CTRIA3', 'CTRIA6', 'CTRIAR',
        'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUAD',
    }
    for eid, elem in model.elements.items():
        elem_type = elem.type
        if elem_type in shells:
            elem.flip_normal()
            if zero_zoffset:
                elem.zoffset = 0.
        elif elem_type in skip_elements:
            pass
        else:
            log.warning(f'cannot flip {elem_type}')
