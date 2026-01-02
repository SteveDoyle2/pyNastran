"""
https://iopscience.iop.org/article/10.1088/1742-6596/1925/1/012058/pdf
http://www.p-over-a.co.uk/toolbox/Huth.aspx
https://help.altair.com/inspire/en_us/topics/shared/glossary/gloss_connection_stiffness_c.htm
"""
import numpy as np

def cis_ass_to_nab(cis: str, ass: str) -> tuple[int, float, float]:
    if cis == 'single_shear':
        n = 1  # single shear
    else:
        assert cis == 'double_shear', cis
        n = 2

    if ass == 'metallic_metallic_riveted':
        a = 2 / 5.
        b = 2.2
    elif ass == 'metallic_metallic_bolted':
        a = 2 / 3.
        b = 3.0
    elif ass == 'composite_composite_bolted':  # graphite epoxy
        a = 2 / 3.
        b = 4.2
    elif ass == 'metallic_composite_bolted':  # graphite epoxy
        a = 2 / 3.  # same for both
        b = (4.2 + 3.0) / 2 # assumed and it's exactly what they use
    else:  # pragma: no cover
        raise NotImplementedError(ass)
    return n, a, b

def get_huth(diameter: float,
             t1: float, t2: float,
             E1: float, E2: float, Ef: float,
             cis: str, ass: str) -> float:
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
    n, a, b = cis_ass_to_nab(cis, ass)
    part2 = 1/(t1*E1) + 1/(n*t2*E2) + 1/(2*t1*Ef) + 1/(2*n*t2*Ef)
    c = ((t1 + t2) / (2 * diameter)) ** a * (b / n) * part2
    k_huth = 1 / c
    return k_huth

def main():
    diameter = 0.1875 # inch
    t1 = 0.040 # inch
    t2 = 0.063 # inch
    E1 = 10500000 # psi
    E2 = 10500000 # psi
    E3 = 16000000 # psi

    # ----------------------------------------------
    x = 1

    # ----------------------------------------------------------
    # metallic_metallic_riveted
    cis = 'single_shear'
    ass = 'metallic_metallic_riveted'
    k_huth = get_huth(diameter, t1, t2, E1, E2, E3, cis, ass)
    print(f'k_huth = {k_huth:g}')  # 147,424  # psi
    assert np.allclose(k_huth, 147424), k_huth
    x = 1

    cis = 'double_shear'
    ass = 'metallic_metallic_riveted'
    k_huth = get_huth(diameter, t1, t2, E1, E2, E3, cis, ass)
    print(f'k_huth = {k_huth:g}')  # 365,986  # psi
    assert np.allclose(k_huth, 365896.0556), k_huth
    # --------------------------
    #  metallic_metallic_bolted
    cis = 'single_shear'
    ass = 'metallic_metallic_bolted'
    k_huth = get_huth(diameter, t1, t2, E1, E2, E3, cis, ass)
    print(f'k_huth = {k_huth:g}')  # 152,589  # psi
    assert np.allclose(k_huth, 152589), k_huth

    cis = 'double_shear'
    ass = 'metallic_metallic_bolted'
    k_huth = get_huth(diameter, t1, t2, E1, E2, E3, cis, ass)
    print(f'k_huth = {k_huth:g}')  # 378,714  # psi
    assert np.allclose(k_huth, 378714), k_huth

    # -------------------------------------------
    # metalic_composite_bolted
    cis = 'single_shear'
    ass = 'composite_composite_bolted'
    k_huth = get_huth(diameter, t1, t2, E1, E2, E3, cis, ass)
    print(f'k_huth = {k_huth:g}')
    assert np.allclose(k_huth, 108992), k_huth

    cis = 'double_shear'  # n=2
    ass = 'composite_composite_bolted'
    k_huth = get_huth(diameter, t1, t2, E1, E2, E3, cis, ass)
    print(f'k_huth = {k_huth:g}')
    assert np.allclose(k_huth,  270510), k_huth
    # -------------------------------------------
    # metalic_composite_bolted
    cis = 'single_shear'
    ass = 'metallic_composite_bolted'
    k_huth = get_huth(diameter, t1, t2, E1, E2, E3, cis, ass)
    print(f'k_huth = {k_huth:g}')
    assert np.allclose(k_huth, 127157), k_huth

    cis = 'double_shear'  # n=2
    ass = 'metallic_composite_bolted'
    k_huth = get_huth(diameter, t1, t2, E1, E2, E3, cis, ass)
    print(f'k_huth = {k_huth:g}')
    assert np.allclose(k_huth, 315595), k_huth
    # -------------------------------------------
    x = 1
