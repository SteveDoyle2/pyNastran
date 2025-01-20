from cpylog import SimpleLogger
from pyNastran.op2.op2_interface.nx_tables import NX_VERSIONS
from pyNastran.op2.op2_interface.utils import reshape_bytes_block


MSC_VERSIONS = [b'V2005R3B', b'XXXXXXXX']
MSC_LONG_VERSION = [
    b'XXXXXXXX20140', b'XXXXXXXX20141', b'XXXXXXXX20142',
    b'XXXXXXXX20150', b'XXXXXXXX20151', b'XXXXXXXX20152',
    b'XXXXXXXX20160', b'XXXXXXXX20161', b'XXXXXXXX20162',
    b'XXXXXXXX20170', b'XXXXXXXX20171', b'XXXXXXXX20172',
    b'XXXXXXXX20180', b'XXXXXXXX20181', b'XXXXXXXX20182',
    b'XXXXXXXX20190', b'XXXXXXXX20191', b'XXXXXXXX20192',
    b'XXXXXXXX20200', b'XXXXXXXX20201', b'XXXXXXXX20202',
    b'XXXXXXXX20210', b'XXXXXXXX20211', b'XXXXXXXX20212', b'XXXXXXXX20214',
    b'XXXXXXXX20220', b'XXXXXXXX20221', b'XXXXXXXX20222',
    b'XXXXXXXX20230', b'XXXXXXXX20231', b'XXXXXXXX20232', b'XXXXXXXX20234',
    b'XXXXXXXX20240', b'XXXXXXXX20241', b'XXXXXXXX20242',  # not checked
]

OPTISTRUCT_VERSIONS = [
    b'OS11XXXX', b'OS12.210', b'OS14.210',
    b'OS2017', b'OS2017.1', b'OS2017.2',          # OS2017 not checked
    b'OS2018', b'OS2018.1', b'OS2018.2',          # OS2018, OS2018.2 not checked

    # Version between 2019 and 2024.1 checked against Altair documentation
    b'OS2019', b'OS2019.1', b'OS2019.2',               # OS2019 not tested
    b'OS2019.1.1', b'OS2019.2.1', b'OS2019.2.2',       # Subversions in documentation, but not tested

    b'OS2020', b'OS2020.1',                            # OS2020.2 removed
    b'OS2020.1.1',                                     # Subversion in documentation, but not tested

    b'OS2021', b'OS2021.1', b'OS2021.2',               # OS2021 not tested
    b'OS2021.2.1',

    b'OS2022', b'OS2022.1', b'OS2022.2', b'OS2022.3',  # OS2022.3 added, OS2022 not tested

    b'OS2023', b'OS2023.1',                            # OS2023.2 removed, OS2023.1 not tested

    b'OS2024', b'OS2024.1', b'OS2024.2',  # OS2024.2 not checked
    b'OS2025', b'OS2025.1', b'OS2025.2',  # not checked
    b'OS2026', b'OS2026.1', b'OS2026.2',  # not checked
]

AUTODESK_VERSIONS = [
    b'NE  0824',  # this means NEi Nastran...
]


def parse_nastran_version(data: bytes, version: bytes, encoding: bytes,
                           log: SimpleLogger) -> tuple[str, str]:
    """parses a Nastran version string"""
    version_str = ''
    if len(data) == 32:
        #self.show_data(data[:16], types='ifsdqlILQ', endian=None)
        #self.show_data(data[16:], types='ifsdqlILQ', endian=None)
        if data[:16].strip() in MSC_LONG_VERSION:
            version2 = reshape_bytes_block(data)
            version_str = version2.decode(encoding).strip()
            # 'XXXXXXXX20140   0   \x00\x00\x00\x00        '
            # 'XXXXXXXX20141   0   \x00\x00\x00\x00        '
            mode = 'msc'
        else:
            raise NotImplementedError(f'check={data[:16].strip()} data={data!r}; '
                                      f'len(data)={len(data)}')
    elif len(data) == 8:
        mode, version_str = _parse_nastran_version_8(data, version, encoding, log)
    elif len(data) == 16:
        mode, version_str = _parse_nastran_version_16(data, version, encoding, log)
    else:  # pragma: no cover
        raise NotImplementedError(f'version={version!r}; n={len(data)}')
    return mode, version_str

def _parse_nastran_version_16(data: bytes, version: bytes, encoding: str,
                              log) -> tuple[str, str]:
    """parses an 8 character version string"""
    version2 = reshape_bytes_block(version, is_interlaced_block=True)
    if version2[:2] == b'NX':
        version_str = version2[2:].decode('latin1')
        if version_str in NX_VERSIONS:
            mode = 'nx'
        else:  # pragma: no cover
            raise RuntimeError(f'unknown version={version_str}')
    elif version2 in MSC_VERSIONS:
        version_str = version2.decode('latin1')
        mode = 'msc'
    else:  # pragma: no cover
        raise RuntimeError(f'unknown version={version!r}; version2={version2!r}')
    return mode, version_str

def _parse_nastran_version_8(data: bytes, version: bytes, encoding: str,
                             log) -> tuple[str, str]:
    """parses an 8 character version string"""
    version_str = version.strip().decode(encoding)
    if version.startswith(b'NX'):
        mode = 'nx'
        version_str = version[2:].strip().decode(encoding)
        if version_str not in NX_VERSIONS:
            log.warning(f'nx version={version_str!r} is not supported')
    elif version.startswith(b'MODEP'):
        # TODO: why is this separate?
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_ac11103.op2
        #print('found NX table?...')
        #log.warning('Assuming NX Nastran')
        mode = 'nx'
    elif version.startswith(b'AEROFREQ'):
        # TODO: why is this separate?
        # C:\Users\Steve\Dropbox\pyNastran_examples\move_tpl\loadf.op2
        #print('found MSC table?...')
        #log.warning('Assuming MSC Nastran')
        mode = 'msc'
    elif version.startswith(b'AEROTRAN'):
        # TODO: why is this separate?
        # C:\Users\Steve\Dropbox\pyNastran_examples\move_tpl\loadf.op2
        #log.warning('Assuming MSC Nastran')
        mode = 'msc'
    elif version in MSC_VERSIONS:
        mode = 'msc'
    #elif version in [b'XXXXXXXX']:
        ##log.warning('Assuming MSC Nastran')
        #mode = 'msc'
    elif version in OPTISTRUCT_VERSIONS:
        # should this be called optistruct or radioss?
        mode = 'optistruct'
    elif version in AUTODESK_VERSIONS:
        mode = 'autodesk'
    #elif data[:20] == b'XXXXXXXX20141   0   ':
        #self.set_as_msc()
        #self.set_table_type()
    elif version == b'ADINAOUT':
        mode = 'adina'
        #print('version_str = ',version_str)
    elif version == b'NASA95':
        mode = 'nasa95'
    else:  # pragma: no cover
        raise RuntimeError(f'unknown version={version!r}')
    return mode, version_str
