from typing import Any
import numpy as np
from pyNastran.dev.bdf_vectorized3.cards.base_card import searchsorted_filter

def get_density_from_property(property_id: np.ndarray,
                              allowed_properties: list[Any]) -> np.ndarray:
    assert len(allowed_properties) > 0, allowed_properties
    nproperties = len(property_id)
    rho = np.full(nproperties, np.nan, dtype='float64')
    for prop in allowed_properties:
        i_lookup, i_all = searchsorted_filter(prop.property_id, property_id, msg='')
        if len(i_lookup) == 0:
            continue

        # we're at least using some properties
        rhoi = prop.rho()
        rho[i_lookup] = rhoi[i_all]
    return rho

def get_density_from_material(material_id: np.ndarray,
                              allowed_materials: list[Any],
                              debug: bool=False) -> np.ndarray:
    assert len(allowed_materials) > 0, allowed_materials
    nmaterials = len(material_id)
    if nmaterials == 0:
        raise RuntimeError(f'material_id={material_id}')

    material_id_check = np.zeros(nmaterials, dtype='int32')
    rho = np.full(nmaterials, np.nan, dtype='float64')
    if debug:  # pragma: no cover
        print(f'material_id = {material_id}')
    for mat in allowed_materials:
        mat_material_ids = mat.material_id
        if debug:  # pragma: no cover
            print('mat.material_id = ', mat.material_id)
            print('mat.rho = ', mat.rho)
            #print('i_lookup = ', i_lookup)

        i_lookup, i_all = searchsorted_filter(mat_material_ids, material_id)
        if len(i_all) == 0:
            continue
        if debug:  # pragma: no cover
            print('i_lookup = ', i_lookup)
            print('i_all = ', i_all)

        material_id_check[i_lookup] = mat_material_ids[i_all]
        if debug:  # pragma: no cover
            print('material_id_check = ', material_id_check)
            print('mat_material_ids = ', mat_material_ids)

        rho[i_lookup] = mat.rho[i_all]
        #if debug:  # pragma: no cover
            #print('material_id_check = ', material_id_check)
            #print('mat_material_ids = ', mat_material_ids)
            #x = 1
            #  Exception: AssertionError: mass/area=[0.26646691]
            #  PCOMP          2
            #      1.4330029      0.     YES       11.205997     45.     YES
            #      11.999867    -45.     YES       11.463921     90.     YES
            #      2      1.      0.     YES       2      1.      0.     YES
            #      1      1.     90.     YES       1      1.    -45.     YES
            #      1      1.     45.     YES       1      1.      0.     YES
            #
            #  mids=[1 1 1 1 2 2 1 1 1 1]
    return rho

def get_mass_from_property(property_id: np.ndarray,
                           allowed_properties: list[Any]) -> np.ndarray:
    assert len(allowed_properties) > 0, allowed_properties
    nproperties = len(property_id)

    mass = np.full(nproperties, np.nan, dtype='float64')
    for prop in allowed_properties:
        i_lookup, i_all = searchsorted_filter(prop.property_id, property_id)
        if len(i_all) == 0:
            continue

        # PBUSH can't have nan mass if it's defined
        prop_mass_og = prop.mass()
        prop_mass = prop_mass_og.copy()
        inan = np.where(np.isnan(prop_mass))
        prop_mass[inan] = 0.0

        # we're at least using some properties
        mass[i_lookup] = prop_mass[i_all]
    return mass

def basic_mass_material_id(property_id: np.ndarray,
                           allowed_properties: list[Any], msg: str='') -> np.ndarray:
    """Intended for elements that can only have a single material id"""
    assert len(allowed_properties) > 0, allowed_properties
    material_id = np.full(len(property_id), np.nan, dtype='float64')
    for prop in allowed_properties:
        i_lookup, i_all = searchsorted_filter(prop.property_id, property_id, msg=msg)
        if len(i_lookup) == 0:
            continue
        material_id[i_lookup] = prop.material_id[i_all]
    return material_id

def expanded_mass_material_id(element_id: np.ndarray,
                              property_id: np.ndarray,
                              allowed_properties: list[Any],
                              msg: str='') -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Intended for elements that can have multiple material ids (e.g., PCOMP)"""
    assert len(allowed_properties) > 0, allowed_properties
    element_ids_list = []
    property_ids_list = []
    material_ids_list = []

    for prop in allowed_properties:
        i_lookup, i_all = searchsorted_filter(prop.property_id, property_id, msg=msg)
        if len(i_lookup) == 0:
            continue
        nproperties = len(prop.property_id)

        element_idi = element_id[i_lookup]
        property_idi1 = property_id[i_lookup]
        property_idi, mass_midi = prop.expanded_mass_material_id()
        if len(mass_midi) == nproperties:
            # [eid, pid, mid]
            #material_id[i_lookup]
            property_idi2 = property_idi[i_all]
            assert np.array_equal(property_idi1, property_idi2)

            material_idi = mass_midi[i_all]
        else:
            element_idi, property_idi, mass_midi = prop.expanded_mass_material_id_lookup(
                element_idi, property_idi1,
            )
            material_idi = mass_midi

        element_ids_list.append(element_idi)
        property_ids_list.append(property_idi)
        material_ids_list.append(material_idi)
        #x = 1

    element_id = np.hstack(element_ids_list)
    property_id = np.hstack(property_ids_list)
    material_id = np.hstack(material_ids_list)
    return element_id, property_id, material_id
