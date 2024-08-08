def assign(op2_filename: str,
           unit: int=12,
           is_binary: bool=True) -> str:
    _format = 'UNFORMATTED' if is_binary else 'FORMATTED'
    out = f'ASSIGN OUTPUT2=static.op2,{_format},UNIT={unit:d},DELETE\n'

def export_kelm(sol: int) -> str:
    if sol == 101:
        msg = [
            'DIAG 8',
            'COMPILE SEMG NOLIST',
            'ALTER 153',
            #'MATPRT KELM',
            "OUTPUT2 KELM,,,,///12///'KELM'",
        ]
    else:
        raise RuntimeError(sol)
    out = '\n'.join(msg)
    return out
