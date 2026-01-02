"""
https://iopscience.iop.org/article/10.1088/1742-6596/1925/1/012058/pdf
http://www.p-over-a.co.uk/toolbox/Huth.aspx
https://help.altair.com/inspire/en_us/topics/shared/glossary/gloss_connection_stiffness_c.htm
"""
import numpy as np

def cis_ass_to_nab(shear_type: str='single',
                   connection_type: str='rivet',
                   mat1_type: str='metallic',
                   mat2_type: str='metallic',
                   ) -> tuple[int, float, float]:
    assert shear_type in {'single', 'double'}, f'shear_type={shear_type!r}'
    assert connection_type in {'rivet', 'bolt'}, f'connection_type={connection_type!r}'
    assert mat1_type in {'metallic', 'composite'}, f'mat1_type={mat1_type!r}'
    assert mat2_type in {'metallic', 'composite'}, f'mat2_type={mat2_type!r}'

    words_dict = {
        ('rivet', 'metallic', 'metallic'): 'metallic_metallic_riveted',
        ('bolt', 'metallic', 'metallic'): 'metallic_metallic_bolted',
        ('bolt', 'metallic', 'composite'): 'metallic_composite_bolted',
        ('bolt', 'composite', 'metallic'): 'metallic_composite_bolted',
        ('bolt', 'composite', 'composite'): 'composite_composite_bolted',
    }
    word = words_dict[(connection_type, mat1_type, mat2_type)]

    if shear_type == 'single':
        n = 1  # single shear
    else: # shear_type has already been checked
        n = 2

    if word == 'metallic_metallic_riveted':
        a = 2 / 5.
        b = 2.2
    elif word == 'metallic_metallic_bolted':
        a = 2 / 3.
        b = 3.0
    elif word == 'composite_composite_bolted':  # graphite epoxy
        a = 2 / 3.
        b = 4.2
    elif word == 'metallic_composite_bolted':  # graphite epoxy
        a = 2 / 3.  # same for both
        b = (4.2 + 3.0) / 2 # assumed and it's exactly what they use
    else:  # pragma: no cover
        raise NotImplementedError(word)
    return n, a, b

def get_huth(diameter: float,
             t1: float, t2: float,
             E1: float, E2: float, Ef: float,
             shear_type: str='single',
             connection_type: str='rivet',
             mat1_type: str='metallic',
             mat2_type: str='metallic',
             ) -> float:
    """
    1 = Plate 1 (central one in double shear)
    2 = Plate 2 (outer ones in double shear)
    f = fastener

    #t1 = tp (the thicker member)
    #t2 = ts (the thinner member)
    # This gives the lower compliance and thus a more
    # conservative pin load than the opposite.
    #
    # Therefore I believe t2 refers to the outer sheets,
    # and it only makes sense to distinguish between
    # t1 and t2 in the double shear case.

    """
    n, a, b = cis_ass_to_nab(
        shear_type=shear_type,
        connection_type=connection_type,
        mat1_type=mat1_type, mat2_type=mat2_type)
    part2 = 1/(t1*E1) + 1/(n*t2*E2) + 1/(2*t1*Ef) + 1/(2*n*t2*Ef)
    c = ((t1 + t2) / (2 * diameter)) ** a * (b / n) * part2
    k_huth = 1 / c
    return k_huth
