from typing import Any
import numpy as np
from pyNastran.op2.op2_helper import polar_to_real_imag


def reshape_bytes_block(block: bytes) -> bytes:
    """
    Converts the nonsense 64-bit string to 32-bit format.

    Note
    ----
    Requires a multiple of 8 characters.
    """
    nwords = len(block) // 2
    block2 = b''.join([block[8*i:8*i+4] for i in range(nwords)])
    return block2

def reshape_bytes_block_size(name_bytes: bytes, size: int=4) -> bytes:
    """
    Converts the nonsense 64-bit string to 32-bit format.  Right strips
    the output

    Example
    -------
    >>> reshape_bytes_block_size('1234567 ', size=4)
    '1234567'
    >>> reshape_bytes_block_size('1234    567     ', size=8)
    '1234567'
    """
    if size == 4:
        name_str = name_bytes.decode('latin1').rstrip()
    else:
        name_str = reshape_bytes_block(name_bytes).decode('latin1').rstrip()
    return name_str

def reshape_bytes_block_strip(name_bytes: bytes, size: int=4) -> str:
    """
    Converts the nonsense 64-bit string to 32-bit format.  Strips
    the output.

    Example
    -------
    >>> reshape_bytes_block_strip('1234567 ', size=4)
    '1234567'
    >>> reshape_bytes_block_strip('1234    567     ', size=8)
    '1234567'
    """
    if size == 4:
        name_str = name_bytes.decode('latin1').strip()
    else:
        name_str = reshape_bytes_block(name_bytes).decode('latin1').strip()
    return name_str

def mapfmt(fmt: bytes, size: int) -> bytes:
    """
    Changes the type of of a format string from 32-bit to 64-bit.

    WARNING: doesn't handle strings.

    Example
    -------
    >>> mapfmt(b'2i 6f')
    b'2q 6d'
    """
    if size == 4:
        return fmt
    return fmt.replace(b'i', b'q').replace(b'f', b'd')

def mapfmt_str(fmt: str, size: int) -> str:
    """Same as mapfmt, but works on strings instead of bytes."""
    if size == 4:
        return fmt
    return fmt.replace('i', 'q').replace('f', 'd')

def build_obj(obj):
    """
    there are some cases in build objects that set things that aren't consistent,
    so this exists to combine those
    """
    if not obj.is_built:
        try:
            obj.build()
        except (KeyboardInterrupt, SystemError, NameError, SyntaxError, AttributeError, MemoryError, TypeError):  # pragma: no cover
            raise
        except Exception as e:
            raise RuntimeError(str(obj)) from e
        obj.is_built = True

def apply_mag_phase(floats: Any, is_magnitude_phase: bool,
                    isave_real: list[int], isave_imag: list[int]) -> Any:
    """converts mag/phase data to real/imag"""
    if is_magnitude_phase:
        mag = floats[:, isave_real]
        phase = floats[:, isave_imag]
        real_imag = polar_to_real_imag(mag, phase)
    else:
        real = floats[:, isave_real]
        imag = floats[:, isave_imag]
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
        if ',' not in superelement:
            split_superelement = superelement.split()
            word, value1 = split_superelement
            assert word == 'SUPERELEMENT', f"split_superelement={split_superelement}; expected something of the form 'SUPERELEMENT 0'"
            subtitle = f'{subtitle}; SUPERELEMENT {value1}'
            value1 = int(value1)

            if superelement_adaptivity_index:
                superelement_adaptivity_index = '%s; SUPERELEMENT %s' % (
                    superelement_adaptivity_index, value1)
            else:
                superelement_adaptivity_index = f'SUPERELEMENT {value1}'
        else:
            split_superelement = superelement.split(',')
            assert len(split_superelement) == 2, f"split_superelement={split_superelement}; expected something of the form 'SUPERELEMENT 0 , 1'"
            word_value1, value2 = split_superelement
            word, value1 = word_value1.split()
            assert word == 'SUPERELEMENT', f"word_value1={word_value1}; expected something of the form 'SUPERELEMENT 0 , 1'"
            value1 = int(value1)
            value2 = int(value2)

            if superelement_adaptivity_index:
                superelement_adaptivity_index = '%s; SUPERELEMENT %s,%s' % (
                    superelement_adaptivity_index, value1, value2)
            else:
                superelement_adaptivity_index = f'SUPERELEMENT {value1},{value2}'
    return superelement_adaptivity_index

def update_subtitle_with_adaptivity_index(subtitle: str, superelement_adaptivity_index: str,
                                          adpativity_index: str) -> str:
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
