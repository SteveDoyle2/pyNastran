from itertools import count
from typing import List, Any
import numpy as np


def build_obj(obj):
    """
    there are some cases in build objects that set things that aren't consistent,
    so this exists to combine those
    """
    if not obj.is_built:
        obj.build()

def apply_mag_phase(floats: Any, is_magnitude_phase: bool,
                    isave1: List[int], isave2: List[int]) -> Any:
    """converts mag/phase data to real/imag"""
    if is_magnitude_phase:
        mag = floats[:, isave1]
        phase = floats[:, isave2]
        rtheta = np.radians(phase)
        real_imag = mag * (np.cos(rtheta) + 1.j * np.sin(rtheta))
    else:
        real = floats[:, isave1]
        imag = floats[:, isave2]
        real_imag = real + 1.j * imag
    return real_imag

def to_mag_phase(real_imag: Any, is_magnitude_phase: bool) -> Any:
    """converts real/imag data to mag/phase"""
    if is_magnitude_phase:
        mag = np.abs(real_imag)
        phase = np.angle(real_imag, deg=True)
        return mag, phase
    else:
        return real_imag.real, real_imag.imag

def get_superelement_adaptivity_index(subtitle: str, superelement: str) -> str:
    """determines the SUPERELEMENT/ADAPTIVITY_INDEX from the subtitle"""
    superelement_adaptivity_index = ''
    if 'SUPERELEMENT' in superelement:
        # 'SUPERELEMENT 0'

        # F:\work\pyNastran\examples\Dropbox\move_tpl\opt7.op2
        # 'SUPERELEMENT 0       ,   1'
        split_superelement = superelement.split()
        if len(split_superelement) == 2:
            word, value1 = split_superelement
            assert word == 'SUPERELEMENT', 'split_superelement=%s' % split_superelement
            subtitle = '%s; SUPERELEMENT %s' % (subtitle, value1)
            value1 = int(value1)

            if superelement_adaptivity_index:
                superelement_adaptivity_index = '%s; SUPERELEMENT %s' % (
                    superelement_adaptivity_index, value1)
            else:
                superelement_adaptivity_index = 'SUPERELEMENT %s' % value1
        elif len(split_superelement) == 4:
            word, value1, comma, value2 = split_superelement
            assert word == 'SUPERELEMENT', 'split_superelement=%s' % split_superelement
            value1 = int(value1)
            value2 = int(value2)

            if superelement_adaptivity_index:
                superelement_adaptivity_index = '%s; SUPERELEMENT %s,%s' % (
                    superelement_adaptivity_index, value1, value2)
            else:
                superelement_adaptivity_index = 'SUPERELEMENT %s,%s' % (value1, value2)
        else:
            raise RuntimeError(split_superelement)
    return superelement_adaptivity_index

def update_subtitle_with_adaptivity_index(subtitle, superelement_adaptivity_index,
                                          adpativity_index):
    """
    Parameters
    ----------
    subtitle : str
        the subtitle
    superelement_adaptivity_index : str
        the superelement marker
    adpativity_index : str
        the mesh adaptivity index
        'ADAPTIVITY INDEX=      1'

    Returns
    -------
    subtitle : str
        the new subtitle
        title + 'SUPERELEMENT 0'
        title + 'SUPERELEMENT 0, 1'
        title + 'ADAPTIVITY_INDEX=1'
        title + 'SUPERELEMENT 0 (id=2000)'
    """
    if adpativity_index.strip():
        def hexi(src):
            return ''.join(str(ord(c)) for c in src)
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\cqrsee101q3.op2
        adpativity_index2 = (
            adpativity_index \
            .replace('\x02\x00\x00\x00', hexi('\x02\x00\x00\x00')) \
            .replace('\x03\x00\x00\x00', hexi('\x03\x00\x00\x00')) \
            .replace('\x04\x00\x00\x00', hexi('\x04\x00\x00\x00')) \
            .replace('\x05\x00\x00\x00', hexi('\x05\x00\x00\x00')) \

            # form feed; \f
            .replace('\x0c\x00\x00\x00', hexi('\x0c\x00\x00\x00')) \

            #  carriage return
            .replace('\r\x00\x00\x00', hexi('\r\x00\x00\x00')) \

            .replace('\x0f\x00\x00\x00', hexi('\x0f\x00\x00\x00')) \
            # null
            .replace('\x00\x00\x00\x00', hexi('\x00\x00\x00\x00')) \
            .strip()
        )
        #print('adpativity_index = %r' % adpativity_index.strip())
        if adpativity_index2.isdigit():
            #print(adpativity_index2)
            subtitle2 = '%s; %s' % (subtitle, adpativity_index2)
            superelement_adaptivity_index = adpativity_index2
        else:
            if 'ADAPTIVITY INDEX=' not in adpativity_index2:
                msg = f'subtitle={subtitle!r}\nsuperelement_adaptivity_index={superelement_adaptivity_index!r} adpativity_index={adpativity_index2!r}'
                raise AssertionError(msg)

            # F:\work\pyNastran\examples\Dropbox\move_tpl\pet1018.op2
            #'ADAPTIVITY INDEX=      1'
            split_adpativity_index = adpativity_index2.split()
            assert len(split_adpativity_index) == 3, split_adpativity_index
            word1, word2, adpativity_index_value = split_adpativity_index
            assert word1 == 'ADAPTIVITY', 'split_adpativity_index=%s' % split_adpativity_index
            assert word2 == 'INDEX=', 'split_adpativity_index=%s' % split_adpativity_index

            adpativity_index_value = int(adpativity_index_value)
            subtitle2 = '%s; ADAPTIVITY_INDEX=%s' % (subtitle, adpativity_index_value)
            if superelement_adaptivity_index:
                superelement_adaptivity_index = '%s; ADAPTIVITY_INDEX=%s' % (
                    superelement_adaptivity_index, adpativity_index_value)
            else:
                superelement_adaptivity_index = 'ADAPTIVITY_INDEX=%s' % adpativity_index_value
    else:
        subtitle2 = subtitle
    return subtitle2, superelement_adaptivity_index

def update_label2(label2, isubcase):
    """strips off SUBCASE from the label2 to simplfify the output keys (e.g., displacements)"""
    # strip off any comments
    # 'SUBCASE  1 $ STAT'
    # 'SUBCASE  1 $ 0.900 P'
    label2 = label2.split('$')[0].strip()

    if label2:
        subcase_expected = 'SUBCASE %i' % isubcase
        subcase_equal_expected = 'SUBCASE = %i' % isubcase
        if subcase_expected == label2:
            label2 = ''
        elif label2 == 'NONLINEAR':
            pass
        elif subcase_expected in label2:
            # 'SUBCASE 10' in 'NONLINEAR    SUBCASE 10'
            nchars = len(subcase_expected)
            ilabel_1 = label2.index(subcase_expected)
            ilabel_2 = ilabel_1 + nchars
            label2_prime = label2[:ilabel_1] + label2[ilabel_2:]
            label2 = label2_prime.strip()
        elif subcase_equal_expected in label2:
            # 'SUBCASE = 10'
            slabel = label2.split('=')
            assert len(slabel) == 2, slabel
            label2 = ''
        elif label2.startswith('NONLINEAR '):
            # 'NONLINEAR    SUBCASE   1'
            # sline =['', '   SUBCASE   1']
            sline = label2.split('NONLINEAR ', 1)
            label2 = 'NONLINEAR ' + sline[1].strip()
        elif 'PVAL ID=' in label2 and 'SUBCASE=' in label2:
            # 'PVAL ID=       1                       SUBCASE=       1'
            # '    PVAL ID=       1                       SUBCASE=       1'
            ilabel2 = label2.index('SUBCASE')
            slabel = label2[:ilabel2].strip().split('=')
            assert slabel[0] == 'PVAL ID', slabel
            label2 = slabel[0].strip() + '=' + slabel[1].strip()
        elif 'SUBCASE' in label2:
            # 'SUBCASE    10'
            # 'SUBCASE = 10'
            # 'SUBCASE = 1    SEGMENT = 1'
            # 'SUBCASE = 1    HARMONIC = 0 ,C'
            slabel = label2.split('$')[0].strip().split()

            # 'SUBCASE    10'
            # 'SUBCASE = 10'
            # 'SUBCASE = 1    SEGMENT = 1'
            # 'SUBCASE = 1    HARMONIC = 0 ,C'
            if len(slabel) == 2:
                label2 = ''
            elif len(slabel) == 3 and slabel[1] == '=':
                label2 = ''
            else:
                assert slabel[0] == 'SUBCASE', slabel

                # 'SEGMENT = 1'
                label2 = slabel[3] + '=' + slabel[5]

        elif 'SUBCOM' in label2:
            subcom, isubcase = label2.split()
            label2 = ''
        elif 'SYM' in label2 or 'REPCASE' in label2:
            # 'SYM 401'
            # 'REPCASE 108'
            pass
        #else:
            #print('label2   = %r' % label2)
            #print('subcasee = %r' % subcase_expected)
    return label2

def grids_comp_array_to_index(grids1, comps1, grids2, comps2,
                              make_matrix_symmetric):
    """maps the dofs"""
    #from pyNastran.femutils.utils import unique2d
    ai = np.vstack([grids1, comps1]).T
    bi = np.vstack([grids2, comps2]).T
    #print('grids2 =', grids2)
    #print('comps2 =', comps2)
    #c = np.vstack([a, b])
    #assert c.shape[1] == 2, c.shape
    #urows = unique2d(c)
    #urows = c

    nid_comp_to_dof_index = {}
    j = 0
    a_keys = set()
    for nid_dof in ai:
        #nid_dof = (int(nid), int(dof))
        nid_dof = tuple(nid_dof)
        if nid_dof not in a_keys:
            a_keys.add(nid_dof)
            if nid_dof not in nid_comp_to_dof_index:
                nid_comp_to_dof_index[nid_dof] = j
                j += 1
    nja = len(a_keys)
    del a_keys

    b_keys = set()
    for nid_dof in bi:
        nid_dof = tuple(nid_dof)
        if nid_dof not in b_keys:
            b_keys.add(nid_dof)
        if nid_dof not in nid_comp_to_dof_index:
            nid_comp_to_dof_index[nid_dof] = j
            j += 1
    njb = len(b_keys)
    del b_keys


    nj = len(nid_comp_to_dof_index)
    if make_matrix_symmetric:
        ja = np.zeros(nj, dtype='int32')
        for i, nid_dof in zip(count(), ai):
            j[i] = nid_comp_to_dof_index[tuple(nid_dof)]
        for i, nid_dof in zip(count(), bi):
            j[i] = nid_comp_to_dof_index[tuple(nid_dof)]
        return j, j, nj, nj, nj
    else:
        ja = np.zeros(grids1.shape, dtype='int32')
        for i, nid_dof in zip(count(), ai.tolist()):
            ja[i] = nid_comp_to_dof_index[tuple(nid_dof)]

        jb = np.zeros(grids2.shape, dtype='int32')
        for i, nid_dof in zip(count(), bi.tolist()):
            jb[i] = nid_comp_to_dof_index[tuple(nid_dof)]

        return ja, jb, nja, njb, nj
