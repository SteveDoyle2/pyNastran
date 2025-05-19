"""
defines readers for BDF objects in the OP2 EPT/EPTS table
"""
#pylint: disable=C0103,R0914
from __future__ import annotations
from struct import unpack, Struct
from functools import partial
from typing import TYPE_CHECKING

import numpy as np

#from pyNastran import is_release
#from pyNastran.bdf.cards.properties.mass import PMASS, NSM, NSML
#from pyNastran.bdf.cards.properties.bars import PBAR, PBARL, PBEND, PBEAM3
#from pyNastran.bdf.cards.properties.beam import PBEAM, PBEAML, PBCOMP
#from pyNastran.bdf.cards.properties.bush import PBUSH, PBUSHT
#from pyNastran.bdf.cards.properties.damper import PDAMP, PVISC
#from pyNastran.bdf.cards.properties.properties import PFAST, PGAP
#from pyNastran.bdf.cards.properties.rods import PROD, PTUBE
from pyNastran.bdf.cards.properties.shell import map_failure_theory_int # PSHEAR, PSHELL, PCOMP
#from pyNastran.bdf.cards.properties.springs import PELAS, PELAST
from pyNastran.bdf.cards.properties.beam import pbeam_op2_data_to_init, pbeaml_op2_data_to_init

#from pyNastran.bdf.cards.thermal.thermal import PCONV, PHBDY, PCONVM
# PCOMPG, PBUSH1D, PBEAML, PBEAM3
from pyNastran.op2.op2_interface.op2_reader import (
    mapfmt, reshape_bytes_block_size) # reshape_bytes_block,
from .utils import get_minus1_start_end, get_ints, get_ints_floats, get_ints_strings, get_ints_floats_strings
from .geom2 import PID_CAP, MID_CAP, DoubleCardError
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.op2_vectorized3.op2_geom import OP2Geom
    from pyNastran.dev.bdf_vectorized3.cards.elements.solid import PSOLID



class EPT:
    """defines methods for reading op2 properties"""

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    def _read_fake(self, data: bytes, n: int) -> int:
        return self.op2._read_fake(data, n)

    def read_ept_4(self, data: bytes, ndata: int):
        return self.op2._read_geom_4(self.ept_map, data, ndata)

    def __init__(self, op2: OP2Geom):
        self.op2 = op2
        self.ept_map = {
            (3201, 32, 55): ['NSM', self._read_nsm],          # record 2
            (52, 20, 181): ['PBAR', self._read_pbar],         # record 11 - buggy
            (9102, 91, 52): ['PBARL', self._read_pbarl],      # record 12 - almost there...
            (2706, 27, 287): ['PCOMP', self._read_pcomp],     # record 22 - buggy
            (302, 3, 46): ['PELAS', self._read_pelas],        # record 39
            (2102, 21, 121): ['PGAP', self._read_pgap],       # record 42
            (902, 9, 29): ['PROD', self._read_prod],          # record 49
            (1002, 10, 42): ['PSHEAR', self._read_pshear],    # record 50
            (2402, 24, 281): ['PSOLID', self._read_psolid],   # record 51
            (2302, 23, 283): ['PSHELL', self._read_pshell],   # record 52
            (1602, 16, 30): ['PTUBE', self._read_ptube],      # record 56

            (5402, 54, 262): ['PBEAM', self._read_pbeam],      # record 14 - not done
            (9202, 92, 53): ['PBEAML', self._read_pbeaml],     # record 15
            (2502, 25, 248): ['PBEND', self._read_pbend],      # record 16 - not done
            (1402, 14, 37): ['PBUSH', self._read_pbush],       # record 19 - not done
            (3101, 31, 219): ['PBUSH1D', self._read_pbush1d],  # record 20 - not done
            (152, 19, 147): ['PCONEAX', self._read_pconeax],   # record 24 - not done
            (11001, 110, 411): ['PCONV', self._read_pconv],    # record 25 - not done
            # record 26
            (202, 2, 45): ['PDAMP', self._read_pdamp],      # record 27 - not done
            (2802, 28, 236): ['PHBDY', self._read_phbdy],   # record 43 - not done
            (402, 4, 44): ['PMASS', self._read_pmass],      # record 48
            (1802, 18, 31): ['PVISC', self._read_pvisc],    # record 59
            (10201, 102, 400): ['PVAL', self._read_pval],   # record 58 - not done
            (2606, 26, 289): ['VIEW', self._read_view],     # record 62 - not done
            (3201, 32, 991) : ['NSM', self._read_nsm_2],  # record
            (3301, 33, 992) : ['NSM1', self._read_nsm1],  # record
            (3701, 37, 995) : ['NSML1', self._read_nsml1_nx],    # record
            (3601, 36, 62): ['NSML1', self._read_nsml1_msc],  # record 7
            (15006, 150, 604): ['PCOMPG', self._read_pcompg],  # record

            (702, 7, 38): ['PBUSHT', self._read_pbusht],  # record 1
            (3301, 33, 56): ['NSM1', self._read_fake],  # record 3
            (3401, 34, 57) : ['NSMADD', self._read_fake],    # record 5
            (3501, 35, 58): ['NSML', self._read_fake],  # record 6
            (3501, 35, 994) : ['NSML', self._read_nsml],
            (1502, 15, 36): ['PAABSF', self._read_fake],  # record 8
            (8300, 83, 382): ['PACABS', self._read_fake],  # record 9
            (8500, 85, 384): ['PACBAR', self._read_fake],  # record 10
            (5403, 55, 349): ['PBCOMP', self._read_pbcomp],  # record 13
            (13301, 133, 509): ['PBMSECT', self._read_fake],  # record 17
            (2902, 29, 420): ['PCONVM', self._read_pconvm],  # record 26
            (1202, 12, 33): ['PDAMPT', self._read_pdampt],  # record 28
            (8702, 87, 412): ['PDAMP5', self._read_pdamp5],  # record 29
            (6802, 68, 164): ['PDUM8', self._read_fake],  # record 37
            (6902, 69, 165): ['PDUM9', self._read_fake],  # record 38
            (1302, 13, 34): ['PELAST', self._read_pelast],  # record 41
            (12001, 120, 480): ['PINTC', self._read_fake],  # record 44
            (12101, 121, 484): ['PINTS', self._read_fake],  # record 45
            (4606, 46, 375): ['PLPLANE', self._read_plplane],  # record 46
            (4706, 47, 376): ['PLSOLID', self._read_plsolid],  # record 47
            (10301, 103, 399): ['PSET', self._read_pset],  # record 57
            (3002, 30, 415): ['VIEW3D', self._read_fake],  # record 63

            (13501, 135, 510) : ['PFAST', self._read_pfast_msc],  # MSC-specific
            (3601, 36, 55) : ['PFAST', self._read_pfast_nx],  # NX-specific
            (3801, 38, 979) : ['PPLANE', self._read_pplane],
            (11801, 118, 560) : ['PWELD', self._read_fake],
            (3401, 34, 993) : ['NSMADD', self._read_nsmadd],
            (9300, 93, 684) : ['ELAR', self._read_fake],
            (9400, 94, 685) : ['ELAR2', self._read_fake],
            (16006, 160, 903) : ['PCOMPS', self._read_fake],

            # MSC-specific
            (14602, 146, 692): ['PSLDN1', self._read_fake],
            (16502, 165, 916): ['PAXSYMH', self._read_fake],
            (13201, 132, 513): ['PBRSECT', self._read_fake],

            (13701, 137, 638): ['PWSEAM', self._read_fake],
            (7001, 70, 632): ['???', self._read_fake],
            (15106, 151, 953): ['PCOMPG1', self._read_fake],
            (3901, 39, 969): ['PSHL3D', self._read_fake],
            (17006, 170, 901): ['MATCID', self._read_fake],

            (9601, 96, 691): ['PJOINT', self._read_fake],
            (16502, 165, 916): ['???', self._read_fake],

            (9701, 97, 692): ['PJOINT2', self._read_fake],
            (13401, 134, 611): ['PBEAM3', self._read_pbeam3],
            (8901, 89, 905): ['PSOLCZ', self._read_fake],
            (9801, 98, 698): ['DESC', self._read_desc],
            #(9701, 97, 692): ['???', self._read_fake],
            #(9701, 97, 692): ['???', self._read_fake],
            #(9701, 97, 692): ['???', self._read_fake],

        }

    def _add_op2_property(self, prop):
        """helper method for op2"""
        return
        op2 = self.op2
        #if prop.pid > 100000000:
            #raise RuntimeError('bad parsing; pid > 100000000...%s' % str(prop))
        #print(str(prop)[:-1])
        ntables = op2.table_names.count(b'EPT') + op2.table_names.count(b'EPTS')
        if prop is None:
            return
        pid = prop.pid
        allow_overwrites = (
            ntables > 1 and
            pid in op2.properties and
            op2.properties[pid].type == prop.type)
        op2._add_methods.add_property_object(prop, allow_overwrites=allow_overwrites)

    def _add_op2_property_mass(self, prop):
        """helper method for op2"""
        op2 = self.op2
        #if prop.pid > 100000000:
            #raise RuntimeError('bad parsing; pid > 100000000...%s' % str(prop))
        #print(str(prop)[:-1])
        ntables = op2.table_names.count(b'EPT') + op2.table_names.count(b'EPTS')
        pid = prop.pid
        allow_overwrites = (
            ntables > 1 and
            pid in op2.properties_mass and
            op2.properties_mass[pid].type == prop.type)
        op2._add_methods.add_property_mass_object(prop, allow_overwrites=allow_overwrites)

    def _add_pconv(self, prop: PCONV) -> None:
        if prop.pconid > 100000000:
            raise RuntimeError('bad parsing pconid > 100000000...%s' % str(prop))
        self.op2._add_methods.add_convection_property_object(prop)

# HGSUPPR

    def _read_desc(self, data: bytes, n: int) -> int:
        """
        RECORD – DESC(9801,98,698)

        Word Name Type Description
        1 DID        I Description identification number
        2 NWORDS     I Number of words for the description string
        3 DESC   CHAR4 Description
        Words 3 repeats NWORDS times

        data = (1, 14, 'FACE CONTACT(1)                                         ')
        """
        op2 = self.op2
        assert self.size == 4, 'DESC size={self.size} is not supported'
        #op2.show_data(data[n:], types='ifs')
        struct_2i = Struct(op2._endian + b'2i')
        while n < len(data):

            datai = data[n:n+8]
            desc_id, nwords = struct_2i.unpack(datai)
            #print(desc_id, nwords)
            ndatai = 8 + nwords * 4
            word_bytes = data[n+8:n+ndatai]
            word = word_bytes.decode('ascii').rstrip()
            assert len(word_bytes) == nwords * 4
            #print('word_bytes =', word_bytes)
            op2.log.warning(f'geom skipping DESC={desc_id}: {word!r}')
            n += ndatai
        assert n == len(data), n
        return n

    def _read_nsml(self, data: bytes, n: int) -> int:
        """
        NX 2019.2
        RECORD – NSML(3501, 35, 994)

        Defines a set of lumped nonstructural mass by ID.
        Word Name Type Description
        1 SID         I Set identification number
        2 PROP(2) CHAR4 Set of properties or elements
        4 ID          I Property of element identification number
        5 VALUE      RS Lumped nonstructural mass value
        Words 4 and 5 repeat until -1 occurs

          ints    = (3, ELEMENT, 0,   200, 0.7, -1, 4, PSHELL, 0, 6401, 4.2, -1)
          floats  = (3, ELEMENT, 0.0, 200, 0.7, -1, 4, PSHELL, 0.0, 6401, 4.2, -1)

        """
        op2 = self.op2
        n0 = n
        #op2.show_data(data[n:])
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        istart, iend = get_minus1_start_end(ints)

        ncards = 0
        size = self.size
        for (i0, i1) in zip(istart, iend):
            #data = (4, ELEMENT, 2.1, 1, 3301, -1, -2)
            assert ints[i1] == -1, ints[i1]
            sid = ints[i0]
            prop_bytes = data[n0+(i0+1)*size:n0+(i0+3)*size]
            #print(sid, prop_bytes)
            ids = ints[i0+4:i1:2].tolist()
            values = floats[i0+5:i1:2].tolist()
            #print(ids, values)
            assert len(ids) == len(values)
            nsm_type = prop_bytes.decode('latin1').rstrip()
            nsml = op2.add_nsml(sid, nsm_type, ids, values)
            #print(nsml)
            str(nsml)
            n += (i1 - i0 + 1) * size
            ncards += 1
        op2.card_count['NSML'] = ncards
        return n

    def _read_nsmadd(self, data: bytes, n: int) -> int:
        """
        NX 2019.2
        (3401, 34, 993)

        RECORD – NSMADD(3401,34,993)
        Combines the nonstructural mass inputs.

        Word Name Type Description
        1 SID I Set identification number
        2 ID  I Set of properties or elements
        Word 2 repeats until End of Record

        (1, 2, 3, 4, -1)
        """
        op2 = self.op2
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        istart, iend = get_minus1_start_end(ints)

        ncards = 0
        istart = [0] + list(iend + 1)
        size = self.size
        for (i0, i1) in zip(istart, iend):
            assert ints[i1] == -1, ints[i1]
            sid, *nsms = ints[i0:i1]
            nsmadd = op2.add_nsmadd(sid, nsms)
            #print(nsmadd)
            str(nsmadd)
            n += (i1 - i0 + 1) * size
            ncards += 1
        op2.card_count['NSMADD'] = ncards
        return n

    def _read_nsml1_nx(self, data: bytes, n: int) -> int:
        """
        NSML1(3701, 37, 995)
        Alternate form of NSML entry. Defines lumped nonstructural mass entries by VALUE, ID list.

        Word Name Type Description
        1 SID      I Set identification number
        2 PROP CHAR4 Set of properties
        3 TYPE CHAR4 Set of elements
        4 VALUE   RS Lumped nonstructural mass value
        5 SPECOPT  I Specification option
        SPECOPT=1 By IDs
          6 ID I Property of element identification number
          Word 6 repeats until -1 occurs
        SPECOPT=2 All
          6 ALL(2) CHAR4 Keyword ALL
          Words 6 and 7 repeat until -1 occurs
        SPECOPT=3 Thru range
          6 ID1         I Starting identification number
          7 THRU(2) CHAR4 Keyword THRU
          9 ID2         I Ending identification number
          Words 6 through 9 repeat until -1 occurs
        SPECOPT=4 Thru range with by
          6 ID1         I Starting identification number
          7 THRU(2) CHAR4 Keyword THRU
          9 ID2         I Ending identification number
          10 BY(2)  CHAR4 Keyword BY
          12 N I Increment
          Words 6 through 12 repeat until -1 occurs

        data = (
            3701, 37, 995,
            1, ELEMENT, 466.2,
            3, 249311, THRU, 250189, -1,
            3, 250656, THRU, 251905, -1,
            3, 270705, THRU, 275998, -1,
            3, 332687, THRU, 334734, -1,
            -2,

            2, ELEMENT, 77.7,
            3, 225740, THRU 227065, -1,
            3, 227602, THRU, 228898, -1,
            3, 229435, THRU, 230743, -1,
            3, 231280, THRU, 233789, -1,
            3, 233922, THRU, 235132, -1,
            3, 235265, THRU, 236463, -1,
            3, 338071, THRU, 341134, -1, -2)
        """
        #ints    = (1, ELEMENT, 466.2,
        #           3, 249311, THRU, 250189, -1,
        #           3, 250656, THRU, 251905, -1,
        #           3, 270705, THRU, 275998, -1,
        #           3, 332687, THRU, 334734, -1,
        #           -2,
        #
        #           2, ELEMENT, 77.7,
        #           3, 225740, THRU 227065, -1,
        #           3, 227602, THRU, 228898, -1,
        #           3, 229435, THRU, 230743, -1,
        #           3, 231280, THRU, 233789, -1,
        #           3, 233922, THRU, 235132, -1,
        #           3, 235265, THRU, 236463, -1,
        #           3, 338071, THRU, 341134, -1, -2)
        op2 = self.op2
        n0 = n
        #op2.show_data(data[n:])
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        iminus2 = np.where(ints == -2)[0]
        istart = [0] + list(iminus2[:-1] + 1)
        iend = iminus2
        #print(istart, iend)
        assert len(data[n:]) > 12, data[n:]
        #op2.show_data(data[n:], types='ifs')

        ncards = 0
        istart = [0] + list(iend + 1)
        size = self.size
        for (i0, i1) in zip(istart, iend):
            #data = (4, ELEMENT, 2.1, 1, 3301, -1, -2)
            assert ints[i1] == -2, ints[i1]
            sid = ints[i0]
            nsm_type = data[n0+(i0+1)*size:n0+(i0+2)*size].decode('latin1').rstrip()
            value = float(floats[i0+3])
            #print(f'sid={sid} nsm_type={nsm_type} value={value}')

            iminus1 = i0 + np.where(ints[i0:i1] == -1)[0]
            #print('-1', iminus1)
            #print('-2', iminus2)
            istart2 = [i0 + 4] + list(iminus1[:-1] + 1)
            iend2 = iminus1
            #print(istart2, iend2)

            for istarti, iendi in zip(istart2, iend2):
                #print(istarti, iendi)
                spec_opt = ints[istarti] # 4
                #print(f'ints[{istarti}] = spec_opt = {spec_opt}')
                if spec_opt == 1:
                    # 6 ID I Property of element identification number

                    ivalues = list(range(istarti, iendi))
                    #print('ivalues =', ivalues)
                    pid_eids = ints[ivalues].tolist()
                    #print('pid_eids =', pid_eids)
                elif spec_opt == 3:
                    # datai = (3, 249311, 'THRU    ', 250189)
                    #print(f'i0={i0}')
                    #datai = data[n0+(i0+6)*size:n0+i1*size]
                    #op2.show_data(datai)
                    ids = ints[istarti:iendi]
                    istart = ids[1]
                    iend = ids[-1]
                    pid_eids = list(range(istart, iend+1))
                else:
                    raise NotImplementedError(spec_opt)

            if nsm_type == 'ELEM':
                nsm_type = 'ELEMENT'
            #for pid_eid in pid_eids:
            #nsml = op2.add_nsml1(sid, nsm_type, pid_eids, [value])
            assert len(pid_eids) > 0, pid_eids
            nsml1 = op2.add_nsml1(sid, nsm_type, value, pid_eids)
            #print(nsml1)
            str(nsml1)
            n += (i1 - i0 + 1) * size
            ncards += 1
        op2.card_count['NSML'] = ncards
        return n

    def _read_nsml1_msc(self, data: bytes, n: int) -> int:
        r"""
        NSML1(3601, 36, 62)

        Word Name Type Description
        1 SID      I Set identification number
        2 PROP CHAR4 Set of property or elements
        3 VALUE   RS Lumped nonstructural mass value
        4 SPECOPT  I Specification option
        SPECOPT=1 By IDs
          5 IDs , =FLG1LIST in ixidlst.prm
          6 ID I Property or element ID
          Word 6 repeats until End of Record
        SPECOPT=2 means ALL, =FLG1ALL in ixidlst.prm
          5 ALL(2) CHAR4 Keyword ALL
          Words 5 through 6 repeat until End of Record
        SPECOPT=3 means THRU range, =FLG1THRU in ixidlst.prm
          5 ID1 I Starting ID
          6 THRU(2) CHAR4 Keyword THRU
          8 ID2 I Ending ID
          Words 5 through 8 repeat until End of Record
        SPECOPT=4 means THRU range with BY, =FLG1THBY in ixidlst.prm
          5 ID1 I Starting ID
          6 THRU(2) CHAR4 Keyword THRU
          8 ID2 I Ending ID
          9 BY(2) CHAR4 Keyword BY
          11 N I Increment
          Words 5 through 11 repeat until End of Record
        End SPECOPT
        Words 4 through max repeat until End of Record

        C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\elsum15.op2

        data = (4, ELEMENT, 2.1, 1, 3301, -1, -2)

        """
        op2 = self.op2
        op2.log.info(f'geom skipping NSML1 in {op2.table_name}; ndata={len(data)-12}')
        #op2.show_data(data[n:], types='ifs')
        #bbb
        return len(data)

    def _read_nsm1(self, data: bytes, n: int) -> int:
        """
        NSM1(3301, 33, 992)

        Defines the properties of a nonstructural mass.
        Word Name Type Description
        1 SID      I Set identification number
        2 PROP CHAR4 Set of properties
        3 TYPE CHAR4 Set of elements
        4 ORIGIN   I Entry origin
        5 VALUE   RS Nonstructural mass value
        6 SPECOPT  I Specification option
        SPECOPT=1 By IDs
          7 ID I
          Word 7 repeats until -1 occurs
        SPECOPT=2 All
          7 ALL(2) CHAR4
          Words 7 and 8 repeat until -1 occurs
        SPECOPT=3 Thru range
          7 ID          I
          8 THRU(2) CHAR4
          10 ID         I
          Words 7 through 10 repeat until -1 occurs
        SPECOPT=4 Thru range with by
          7 ID I
          8 THRU(2) CHAR4
          10 ID        I
          11 BY(2) CHAR4
          13 N         I
          Words 7 through 13 repeat until -1 occurs

        data = (3, PCOMP,   0, 0.37, 2, ALL, -1,
                4, ELEMENT, 2, 2.1, 1, 3301, -1)

        """
        op2 = self.op2
        #op2.show_data(data[n:], types='ifs')
        n0 = n
        #op2.show_data(data[n:])
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        istart, iend = get_minus1_start_end(ints)

        ncards = 0
        size = self.size
        for (i0, i1) in zip(istart, iend):
            assert ints[i1] == -1, ints[i1]
            # 1 SID      I Set identification number
            sid = ints[i0]

            # 2 PROP CHAR4 Set of properties
            # 3 TYPE CHAR4 Set of elements
            # 4 ORIGIN   I Entry origin
            # 5 VALUE   RS Nonstructural mass value
            # 6 SPECOPT  I Specification option
            nsm_type = data[n0+(i0+1)*size:n0+(i0+3)*size].decode('latin1').rstrip()
            zero_two = ints[i0+3]
            value = float(floats[i0+4])
            spec_opt = ints[i0+5]
            assert zero_two in [0, 2], zero_two
            #nii = 6
            #print(ints[i0+nii:i1])
            #print(floats[i0+nii:i1])
            #print(sid, nsm_type, value, spec_opt)

            iminus1 = i0 + np.where(ints[i0:i1] == -1)[0]
            #print('-1', iminus1)
            #print('-2', iminus2)
            istart2 = [i0 + 5] + list(iminus1[:-1] + 1)
            iend2 = iminus1
            #print(istart2, iend2)

            if spec_opt == 1:
                # 7 ID I
                ids = ints[i0+6:i1]
            elif spec_opt == 2:
                word = data[n0+(i0+6)*size:n0+i1*size]
                ids = word
            elif spec_opt == 3:  # thru
                # datai = (249311, 'THRU    ', 250189)
                #datai = data[n0+(i0+6)*size:n0+i1*size]
                ids = ints[i0+6:i1]
                istart = ids[0]
                iend = ids[-1]
                ids = list(range(istart, iend+1))
            else:
                raise NotImplementedError(spec_opt)
            #print(sid, nsm_type, zero_two, value, ids)
            #if nsm_type == 'ELEM':
                #nsm_type = 'ELEMENT'
            #for pid_eid in pid_eids:
            #nsml = self.add_nsml1(sid, nsm_type, pid_eids, [value])
            nsm1 = op2.add_nsm1(sid, nsm_type, value, ids)
            #print(nsm1)
            str(nsm1)
            n += (i1 - i0 + 1) * size
            ncards += 1
        op2.card_count['NSM1'] = ncards
        return n

    def _read_nsm(self, data: bytes, n: int) -> int:
        """NSM"""
        op2 = self.op2
        n = op2.reader_geom2._read_dual_card(
            data, n,
            self._read_nsm_nx, self._read_nsm_msc,
            'NSM', None)
        return n

    def _read_nsm_2(self, data: bytes, n: int) -> int:
        """
        NX 2019.2
        NSM(3201, 32, 991)

        RECORD – NSM(3201,32,991)
        Defines the properties of a nonstructural mass.

        Word Name Type Description
        1 SID      I Set identification number
        2 PROP CHAR4 Set of properties
        3 TYPE CHAR4 Set of elements   <---- not right...it's an integer and not used...
        4 ID       I Property or element identification number
        5 VALUE   RS Nonstructural mass value
        Words 5 through 6 repeat until End of Record

        NSM,2,conrod,1007,0.3

        data    = (2, CONROD,  0, 1007, 0.3, -1,
                   2, ELEMENT, 0,  200, 0.20, -1,
                   3, PSHELL,  0, 3301, 0.20, -1,
                   3, ELEMENT, 2,  200, 1.0, -1,
                   4, PSHELL,  2, 6401, 4.2, -1)
        """
        op2 = self.op2
        n0 = n
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        istart, iend = get_minus1_start_end(ints)

        ncards = 0
        size = self.size
        for (i0, i1) in zip(istart, iend):
            #data = (4, ELEMENT, 2.1, 1, 3301, -1, -2)
            assert ints[i1] == -1, ints[i1]
            sid = ints[i0]
            prop_type = data[n0+(i0+1)*size:n0+(i0+3)*size]
            elem_type = data[n0+(i0+3)*size:n0+(i0+4)*size]
            nsm_type = prop_type.decode('latin1').rstrip()
            dunno_int = ints[i0+3]
            #print(ints[i0+4:i1])
            #print(floats[i0+4:i1])
            ids = ints[i0+4:i1:2].tolist()
            values = floats[i0+5:i1:2].tolist()
            assert len(ids) == len(values)
            assert dunno_int in [0, 2], (sid, prop_type, (ints[i0+3], floats[i0+4]), ids, values)
            #print(sid, prop_type, (ints[i0+3], floats[i0+4]), ids, values)
            nsm = op2.add_nsm(sid, nsm_type, ids, values)
            #print(nsm[0])
            str(nsm)
            n += (i1 - i0 + 1) * size
            ncards += 1
        op2.card_count['NSM'] = ncards
        return n

    def _read_nsm_msc(self, data: bytes, n: int) -> int:
        """
        NSM(3201,32,55) - the marker for Record 2

        MSC
        1 SID       I Set identification number
        2 PROP  CHAR4 Set of property or elements
        3 ID        I Property or element identification number
        4 VALUE    RS Nonstructural mass value
        ORIGIN=0 NSM Bulk Data entry
          5 ID I Property or element ID
          6 VALUE RS Nonstructural mass value
          Words 5 through 6 repeat until End of Record
        ORIGIN=2 NSML Bulk Data entry
          5 ID I Property or element ID
          6 VALUE RS Nonstructural mass value
          Words 5 through 6 repeat until End of Record
        Words 3 through 4 repeat until End of Record
        """
        op2 = self.op2
        model = op2
        properties = []
        struct1 = Struct(op2._endian + b'i 4s if')
        ndelta = 16

        i = 0
        op2.show_data(data)
        ints = np.frombuffer(data[n:], op2.idtype).copy()
        floats = np.frombuffer(data[n:], op2.fdtype).copy()

        nsm = model.nsm
        while n < len(data):
            edata = data[n:n+ndelta]
            out = struct1.unpack(edata)
            (sid, prop_set, pid, value) = out
            #            538976312
            assert pid < 100000000, f'sid={sid} prop_set={prop_set} pid={pid} value={value}'
            i += 4
            n += ndelta

            prop_set = prop_set.decode('utf8').rstrip(' ') # \x00
            values = [value]
            #print('ints[i:]=', ints[i:])
            while ints[i] != -1:
                value2 = floats[i]
                values.append(value2)
                n += 4
                i += 1
            op2.log.info("MSC: NSM-sid=%s prop_set=%s pid=%s values=%s" % (
                sid, prop_set, pid, values))

            cards.append()
            card = [sid, prop_set, pid, values]
            cards.append(card)
            #op2._add_methods.add_nsm_object(prop)
            #properties.append(prop)

            # handle the trailing -1
            i += 1
            n += 4

        for card in cards:
            nsm.add(*card)
        nsm.parse_cards()
        properties = []
        return n, properties

    def _read_nsm_nx(self, data: bytes, n: int) -> int:
        """
        NSM(3201,32,55) - the marker for Record 2

        1 SID         I Set identification number
        2 PROP(2) CHAR4 Set of properties or elements
        4 ORIGIN      I  Entry origin
        5 ID          I  Property or element identification number
        6 VALUE      RS Nonstructural mass value
        Words 5 through 6 repeat until End of Record
        """
        op2 = self.op2
        model = op2
        nsm = model.nsm
        nsml = model.nsml
        properties = []

        #NX: C:\Users\sdoyle\Dropbox\move_tpl\nsmlcr2s.op2
        struct1 = Struct(op2._endian + b'i 8s ii f')
        ndelta = 24
        #op2.show_data(data[12:], 'ifs')

        i = 0
        ints = np.frombuffer(data[n:], op2.idtype).copy()
        floats = np.frombuffer(data[n:], op2.fdtype).copy()

        unused_packs = break_by_minus1(ints)
        #for pack in packs:
            #print(pack)

        #ipack = 0
        cards = []
        while n < len(data):
            #print('ints[i:]=', ints[i:].tolist())
            #i1, i2 = packs[ipack]
            #print('idata=%s' % idata[i1:i2])
            #print('fdata=%s' % fdata[i1:i2])
            #print(idata[i1:i2])
            edata = data[n:n+ndelta]
            out = struct1.unpack(edata)
            (sid, prop_set, origin, pid, value) = out
            #            538976312
            assert pid < 100000000
            i += 6
            n += ndelta

            prop_set = prop_set.decode('utf8').rstrip(' ') # \x00
            pids = [pid]
            values = [value]
            #print('ints[i:]=', ints[i:].tolist())
            while ints[i] != -1:
                pid = ints[i]
                value2 = floats[i+1]
                assert pid != -1
                pids.append(pid)
                values.append(value2)
                n += 8
                i += 2
            card = (origin, sid, prop_set, pids, values)
            assert origin in [0, 1], card
            #if origin == 0:
                ##print(f'NX NSM: sid={sid} prop_set={prop_set} pids={pids} values={values}')
                #nsm.add(sid, prop_set, pids, values)
            #elif origin == 1:
                ##print(f'NX NSML: sid={sid} prop_set={prop_set} pids={pids} values={values}')
                #nsml.add(sid, prop_set, pids, values)
            #else:
                #raise NotImplementedError((sid, prop_set, pids, values))
            cards.append(card)
            #print('----')

            # handle the trailing -1
            i += 1
            n += 4
            #ipack += 1
        for (origin, sid, prop_set, pids, values) in cards:
            if origin == 1:
                obj = nsm
            else:
                obj = nsml
            obj.add(sid, prop_set, pids, values)

        nsm.parse_cards()
        nsml.parse_cards()
        return n, properties

# NSM1
# NSML1
# NSMADD
# NSML
# NSML1
# PAABSF
# PACABS
# PACBAR

    def _read_pbar(self, data: bytes, n: int) -> int:
        """
        PBAR(52,20,181) - the marker for Record 11
        .. warning:: this makes a funny property...

        MSC 2016/NX10

        Word Name Type Description
        1  PID  I Property identification number
        2  MID  I Material identification number
        3  A   RS Area
        4  I1  RS Area moment of inertia in plane 1
        5  I2  RS Area moment of inertia in plane 2
        6  J   RS Torsional constant
        7  NSM RS Nonstructural mass per unit length
        8  FE  RS
        9  C1  RS Stress recovery location at point C in element y-axis
        10 C2  RS Stress recovery location at point C in element z-axis
        11 D1  RS Stress recovery location at point D in element y-axis
        12 D2  RS Stress recovery location at point D in element z-axis
        13 E1  RS Stress recovery location at point E in element y-axis
        14 E2  RS Stress recovery location at point E in element z-axis
        15 F1  RS Stress recovery location at point F in element y-axis
        16 F2  RS Stress recovery location at point F in element z-axis
        17 K1  RS Area factor for shear in plane 1
        18 K2  RS Area factor for shear in plane 2
        19 I12 RS Area product of inertia for plane 1 and 2
        """
        op2 = self.op2
        ntotal = 76 * self.factor  # 19*4
        nproperties = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(data, n, nproperties, 19,
                                          size=op2.size, endian=op2._endian)

        prop = op2.pbar
        assert ints.shape[1] == 19
        property_id = ints[:, 0]
        material_id = ints[:, 1]
        assert len(property_id) == len(np.unique(property_id))
        A = floats[:, 2]
        I = floats[:, [3, 4, 18]]
        J = floats[:, 5]
        nsm = floats[:, 6]
        c = floats[:, [8, 9]]
        d = floats[:, [10, 11]]
        e = floats[:, [12, 13]]
        f = floats[:, [14, 15]]
        k = floats[:, [16, 17]]
        #prop.n = nproperties
        prop._save(property_id, material_id, A, J, c, d, e, f, I, k, nsm, force=True)
        prop.write()

        #struct1 = Struct(mapfmt(op2._endian + b'2i17f', self.size))
        #for unused_i in range(nentries):
            #edata = data[n:n+ntotal]
            #out = struct1.unpack(edata)
            ##(pid, mid, a, I1, I2, J, nsm, fe, c1, c2, d1, d2,
             ##e1, e2, f1, f2, k1, k2, I12) = out
            #prop = PBAR.add_op2_data(out)
            #self._add_op2_property(prop)
            #n += ntotal
        op2.card_count['PBAR'] = nproperties
        return n

    def _read_pbarl(self, data: bytes, n: int) -> int:
        """
        PBARL(9102,91,52) - the marker for Record 12
        TODO: buggy
        It's possible to have a PBARL and a PBAR at the same time.
        NSM is at the end of the element.
        """
        op2 = self.op2
        pbar = op2.pbar
        pbar.parse_cards()
        pbar_pids_to_remove = set()

        valid_types = {
            'ROD': 1,
            'TUBE': 2,
            'TUBE2': 2,
            'I': 6,
            'CHAN': 4,
            'T': 4,
            'BOX': 4,
            'BAR': 2,
            'CROSS': 4,
            'H': 4,
            'T1': 4,
            'I1': 4,
            'CHAN1': 4,
            'Z': 4,
            'CHAN2': 4,
            "T2": 4,
            'BOX1': 6,
            'HEXA': 3,
            'HAT': 4,
            'HAT1': 5,
            'DBOX': 10,  # was 12
            #'MLO TUBE' : 2,
        }  # for GROUP="MSCBML0"

        size = self.size
        ntotal = 28 * self.factor  # 7*4 - ROD - shortest entry...could be buggy... # TODO fix this
        if size == 4:
            struct1 = Struct(op2._endian + b'2i 8s 8s f')
        else:
            struct1 = Struct(op2._endian + b'2q 16s 16s d')

        #nentries = (len(data) - n) // ntotal
        #print(self.show_ndata(80))
        ndata = len(data)

        nproperties = 0
        while ndata - n > ntotal:
            edata = data[n:n+ntotal]
            n += ntotal

            out = struct1.unpack(edata)
            (pid, mid, group, bar_type_bytes, value) = out
            if pid > 100000000 or pid < 1:
                op2.log.debug("  pid=%s mid=%s group=%r bar_type=%r value=%s" % (
                    pid, mid, group, bar_type_bytes, value))
                raise RuntimeError('bad parsing...')

            bar_type = reshape_bytes_block_size(bar_type_bytes, size=size)
            group = reshape_bytes_block_size(group, size=size)
            data_in = [pid, mid, group, bar_type, value]

            expected_length = valid_types[bar_type]
            iformat = op2._endian + b'%if' % expected_length
            #print(iformat)

            ndelta = expected_length * 4
            dims_nsm = list(unpack(iformat, data[n:n+ndelta]))
            data_in += dims_nsm
            #print("  pid=%s mid=%s group=%r bar_type=%r value=%s dims_nsm=%s" % (
                #pid, mid, group, bar_type, value, dims_nsm))

            # TODO why do i need the +4???
            #      is that for the nsm?
            #min_len =  expected_length * 4 + 4
            #if len(data)
            #data = data[n + expected_length * 4 + 4:]
            n += ndelta

            #prin( "len(out) = ",len(out)))
            #print("PBARL = %s" % data_in)
            if pid in pbar.property_id:
                pbar_pids_to_remove.add(pid)

            if 1:
                dim = data_in[4:-1]
                #dim = dims_nsm[:-1]
                nsm = dims_nsm[-1]
                #print("  pid=%s mid=%s group=%r bar_type=%r dim=%s nsm=%s" % (
                    #pid, mid, group, bar_type, dim, nsm))
                op2.add_pbarl(pid, mid, bar_type, dim, group=group, nsm=nsm, comment='')
            else:
                asdf
                prop = PBARL.add_op2_data(data_in)  # last value is nsm
                pid = prop.pid
                if pid in op2.properties:
                    #op2.log.debug("removing:\n%s" % op2.properties[pid])
                    op2._type_to_id_map['PBAR'].remove(pid)
                    del op2.properties[pid]
                self._add_op2_property(prop)

            #op2.properties[pid] = prop
            #print(prop.get_stats())
            #print(op2.show_data(data[n-8:-100]))

            # the PBARL ends with a -1 flag
            #value, = unpack(op2._endian + b'i', data[n:n+4])
            nproperties = 1
            n += 4 * self.factor

        # TODO: filter PBARs
        #if len(op2._type_to_id_map['PBAR']) == 0 and 'PBAR' in op2.card_count:
            #del op2._type_to_id_map['PBAR']
            #del op2.card_count['PBAR']

        op2.pbarl.parse_cards()
        if pbar_pids_to_remove:
            pbar_pids_to_remove = np.array(list(pbar_pids_to_remove), dtype='int32')
            op2.pbar = pbar.remove_card_by_id(pbar_pids_to_remove)

        if nproperties:
            op2.card_count['PBARL'] = nproperties
        #op2.increase_card_count('PBARL', count_num=1)
        #assert len(data) == n
        if self.size == 8:
            n += 16
            #n += 8  # same for 32/64 bit - not 100% that it's always active
        return n

    def _read_pbcomp(self, data: bytes, n: int) -> int:
        """
        PBCOMP(5403, 55, 349)

                    pid      mid  A      I1      I2         I12 J           NSM
        PBCOMP         3       2 2.00E-4 6.67E-9 1.67E-9    0.0 4.58E-9     0.0 +
                 pid mid
        floats = (3, 2, 0.0002, 6.67e-09, 1.67e-09, 0.0, 4.58e-09, 0.0, 1.0, 1.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        ints   = (3, 2, 0.0002, 6.67E-9, 1.67E-9, 0, 4.58E-9, 0, 1.0, 1.0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

        """
        op2 = self.op2
        model = op2
        struct1 = Struct(mapfmt(op2._endian + b'2i 12f i', self.size))
        struct2 = Struct(mapfmt(op2._endian + b'3f 2i', self.size))
        nproperties = 0
        ntotal1 = 60 * self.factor  # 4*15
        ntotal2 = 20 * self.factor

        ndata = len(data)
        #print(ntotal1, ntotal2)
        if self.factor == 2:
            op2.show_data(data[12*self.factor:], types='qd')
        #print(len(data[12*self.factor:]))
        while n < ndata:
            #op2.log.debug(f"n={n} ndata={ndata}")
            edata = data[n:n+ntotal1]
            #if len(edata) == ntotal1:
            data1 = struct1.unpack(edata)
            #else:
                #op2.show_data(edata, types='qdi')
                #n += ntotal2
                #continue
            nsections = data1[-1]
            if op2.is_debug_file:
                (pid, mid, a, i1, i2, i12, j, nsm, k1, k2,
                 m1, m2, n1, n2, unused_nsections) = data1
                op2.log.info(f'PBCOMP pid={pid} mid={mid} nsections={nsections} '
                              f'k1={k1} k2={k2} m=({m1},{m2}) n=({n1},{n2})\n')
            #if pid > 0 and nsections == 0:
                #print('n1')
                #n += ntotal1
                #continue
            #if pid == 0 and nsections == 0:
                #print('n2')
                #n += ntotal2
                #continue

            data2 = []
            n += ntotal1
            if nsections in [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]:
                # 16 Y   RS    Lumped area location along element's y-axis
                # 17 Z   RS    Lumped area location along element's z-axis
                # 18 C   RS    Fraction of the total area for the lumped area
                # 19 MID I     Material identification number
                # 20     UNDEF None
                # Words 16 through 20 repeat NSECT times
                for unused_i in range(nsections):
                    datai = data[n:n+ntotal2]
                    xi, yi, ci, mid, unused_null = struct2.unpack(datai)
                    data2.append((xi, yi, ci, mid))
                    n += ntotal2
            else:
                op2.log.error(f'PBCOMP={data1[0]} has no sections; check your bdf')
                return n
                #raise NotImplementedError('PBCOMP nsections=%r' % nsections)

            if op2.is_debug_file:
                op2.binary_debug.write('     PBCOMP: %s\n' % str([data1, data2]))
                msg = (
                    '    i=%-2s so=%s xxb=%.1f a=%g i1=%g i2=%g i12=%g j=%g nsm=%g '
                    'c=[%s,%s] d=[%s,%s] e=[%s,%s] f=[%s,%s]' % (
                        nsections, None, -9999., a, i1, i2, i12, j, nsm,
                        None, None, None, None, None, None, None, None,)
                )
                op2.log.debug(msg)
            #op2.log.debug(data1)
            #op2.log.debug(data2)

            data_in = [data1, data2]
            prop = model.add_pbcomp(data_in)
            pid = data1[0]
            if pid in op2.properties:
                op2._type_to_id_map['PBEAM'].remove(pid)
                del op2.properties[pid]

            self._add_op2_property(prop)
            nproperties += 1
        #print(f"n={n} ndata={ndata}")
        assert nproperties > 0, f'PBCOMP nproperties={nproperties:d}'
        if len(op2._type_to_id_map['PBEAM']) == 0 and 'PBEAM' in op2.card_count:
            del op2._type_to_id_map['PBEAM']
            del op2.card_count['PBEAM']
        op2.card_count['PBCOMP'] = nproperties
        return n

    def _read_pbeam(self, data: bytes, n: int) -> int:
        """
        PBEAM(5402,54,262) - the marker for Record 14
        .. todo:: add object
        """
        op2 = self.op2
        model = op2
        #op2.log.warning('skipping PBEAM')
        #return len(data)
        cross_section_type_map = {
            0 : 'variable',
            1 : 'constant',
            2 : '???',
        }

        struct1 = Struct(mapfmt(op2._endian + b'4if', self.size))
        struct2 = Struct(mapfmt(op2._endian + b'16f', self.size))
        struct3 = Struct(mapfmt(op2._endian + b'16f', self.size))
        unused_ntotal = 768 # 4*(5+16*12)
        #nproperties = (len(data) - n) // ntotal
        #assert nproperties > 0, 'ndata-n=%s n=%s datai\n%s' % (len(data)-n, n, op2.show_data(data[n:100+n]))
        ndata = len(data)
        #op2.show_data(data[12:], 'if')
        #assert ndata % ntotal == 0, 'ndata-n=%s n=%s ndata%%ntotal=%s' % (len(data)-n, n, ndata % ntotal)
        nproperties = 0

        ntotal1 = 20 * self.factor
        ntotal2 = 64 * self.factor
        while n < ndata:
        #while 1: #for i in range(nproperties):
            edata = data[n:n+ntotal1]
            n += ntotal1
            data_in = list(struct1.unpack(edata))
            #if op2.is_debug_file:
                #op2.log.info('PBEAM pid=%s mid=%s nsegments=%s ccf=%s cweld=%s\n' % tuple(data_in))
            (pid, unused_mid, unused_nsegments, ccf, cweld) = data_in
            #op2.log.info('PBEAM pid=%s mid=%s nsegments=%s ccf=%s cweld=%s' % tuple(data_in))

            # Constant cross-section flag: 1=yes and 0=no
            # what is 2?
            if ccf not in [0, 1, 2]:
                msg = ('  PBEAM pid=%s mid=%s nsegments=%s ccf=%s cweld=%s; '
                       'ccf must be in [0, 1, 2]\n' % tuple(data_in))
                raise ValueError(msg)

            cross_section_type = cross_section_type_map[ccf]
            #print('cross_section_type = %s' % cross_section_type)

            is_pbcomp = False
            is_bad_so = False

            so = []
            xxb = []
            for i in range(11):
                edata = data[n:n+ntotal2]
                if len(edata) != ntotal2:
                    endpack = []
                    raise RuntimeError(f'PBEAM unexpected length i={i:d}...')
                n += ntotal2
                pack = struct2.unpack(edata)
                (soi, xxbi, a, i1, i2, i12, j, nsm, c1, c2,
                 d1, d2, e1, e2, f1, f2) = pack
                xxb.append(xxbi)
                so.append(soi)

                if soi == 0.0:
                    so_str = 'NO'
                elif soi == 1.0:
                    so_str = 'YES'
                else:
                    so_str = str(soi)
                    is_bad_so = True
                    #msg = 'PBEAM pid=%s i=%s x/xb=%s soi=%s; soi not in 0.0 or 1.0' % (
                        #pid, i, xxb, soi)
                    #raise NotImplementedError(msg)

                #if xxb != 0.0:
                    #msg = 'PBEAM pid=%s i=%s x/xb=%s soi=%s; xxb not in 0.0 or 1.0' % (
                        #pid, i, xxb, soi)
                    #raise NotImplementedError(msg)

                pack2 = (so_str, xxbi, a, i1, i2, i12, j, nsm, c1, c2,
                         d1, d2, e1, e2, f1, f2)
                data_in.append(pack2)
                if op2.is_debug_file:
                    op2.binary_debug.write(f'     {pack}\n')
                    msg = (
                        '    i=%-2s' % i + ' so=%s xxb=%.1f a=%g i1=%g i2=%g i12=%g j=%g nsm=%g '
                        'c=[%s,%s] d=[%s,%s] e=[%s,%s] f=[%s,%s]' % (tuple(pack2))
                    )
                    op2.binary_debug.write(msg)
                #msg = (
                    #'    i=%-2s' % i + ' so=%s xxb=%.1f a=%g i1=%g i2=%g i12=%g j=%g nsm=%g '
                    #'c=[%s,%s] d=[%s,%s] e=[%s,%s] f=[%s,%s]' % (tuple(pack2))
                #)
                #print(msg)

            edata = data[n:n+ntotal2]
            if len(edata) != ntotal2:
                endpack = []
                raise RuntimeError('PBEAM unexpected length 2...')
            endpack = struct3.unpack(edata)
            n += ntotal2

            assert len(endpack) == 16, endpack
            #(k1, k2, s1, s2, nsia, nsib, cwa, cwb, # 8
             #m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b) = endpack # 8 -> 16
            if op2.is_debug_file:
                op2.binary_debug.write('    k=[%s,%s] s=[%s,%s] nsi=[%s,%s] cw=[%s,%s] '
                                        'ma=[%s,%s] mb=[%s,%s] na=[%s,%s] nb=[%s,%s]' % (
                                            tuple(endpack)))
            data_in.append(endpack)

            if is_bad_so:
            #if soi < 0.:
                xxb_str = ', '.join(['%g' % xxbi for xxbi in xxb])
                so_str = ', '.join(['%g' % soi for soi in so])
                msg = (f'PBEAM pid={pid} i={i} soi=[{so_str}]; '
                       'soi not 0.0 or 1.0; assuming PBCOMP & dropping')
                op2.log.error(msg)
                is_pbcomp = True

            if min(xxb) < 0.0 or max(xxb) > 1.0:
                xxb_str = ', '.join(['%g' % xxbi for xxbi in xxb])
                msg = (f'PBEAM pid={pid} i={i} x/xb=[{xxb_str}]; '
                       'x/xb must be between 0.0 and 1.0; assuming PBCOMP & dropping')
                op2.log.error(msg)
                is_pbcomp = True

            if is_pbcomp:
                continue
            if pid in op2.pbcomp.property_id:
                continue

            (pid, mid, xxb, so, area, i1, i2, i12, j, nsm,
             c1, c2, d1, d2, e1, e2, f1, f2,
             k1, k2, s1, s2,
             nsia, nsib, cwa, cwb,
             m1a, m2a, m1b, m2b,
             n1a, n2a, n1b, n2b) = pbeam_op2_data_to_init(data_in)
            pbeam = model.add_pbeam(
                pid, mid, xxb, so, area, i1, i2, i12, j, nsm,
                c1, c2, d1, d2, e1, e2, f1, f2,
                k1=k1, k2=k2, s1=s1, s2=s2,
                nsia=nsia, nsib=nsib, cwa=cwa, cwb=cwb,
                m1a=m1a, m2a=m2a, m1b=m1b, m2b=m2b,
                n1a=n1a, n2a=n2a, n1b=n1b, n2b=n2b)

            #prop = PBEAM.add_op2_data(data_in)
            nproperties += 1
            #self._add_op2_property(prop)
        if nproperties:
            op2.card_count['PBEAM'] = nproperties
        pbeam = model.pbeam
        pbeam.parse_cards()
        pbeam.remove_duplicates()
        print(pbeam.property_id)
        return n

    def _read_pbeaml(self, data: bytes, n: int) -> int:
        """
        PBEAML(9202,92,53)

        Word Name Type Description
        1 PID        I   Property identification number
        2 MID        I   Material identification number
        3 GROUP(2) CHAR4 Cross-section group name
        5 TYPE(2)  CHAR4 Cross section type
        7 VALUE      RS  Cross section values for XXB, SO, NSM, and dimensions
        Word 7 repeats until (-1) occurs
        """
        op2 = self.op2
        #strs = numpy.core.defchararray.reshapesplit(data, sep=",")
        #ints = np.frombuffer(data[n:], self._uendian + 'i').copy()
        #floats = np.frombuffer(data[n:], self._uendian + 'f').copy()
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        istart, iend = get_minus1_start_end(ints)

        size = self.size
        nproperties = len(istart)
        if size == 4:
            struct1 = Struct(op2._endian + b'2i 8s 8s')
        else:
            struct1 = Struct(op2._endian + b'2q 16s 16s')

        for unused_i, (istarti, iendi) in enumerate(zip(istart, iend)):
            idata = data[n+istarti*size : n+(istarti+6)*size]
            pid, mid, group, beam_type = struct1.unpack(idata)

            group = reshape_bytes_block_size(group, size=size)
            beam_type = reshape_bytes_block_size(beam_type, size=size)
            fvalues = floats[istarti+6: iendi]
            if op2.is_debug_file:
                op2.binary_debug.write('     %s\n' % str(fvalues))
                op2.log.debug(f'pid={pid:d} mid={mid:d} group={group} beam_type={beam_type}')
                op2.log.debug(fvalues)
            #op2.log.debug(f'pid={pid:d} mid={mid:d} group={group} beam_type={beam_type}')
            data_in = [pid, mid, group, beam_type, fvalues]

            pid, mid, beam_type, xxb, dims, group, so, nsm = pbeaml_op2_data_to_init(
                data_in, op2.pbeaml.valid_types)
            op2.add_pbeaml(pid, mid, beam_type, xxb, dims,
                           so=so, nsm=nsm, group=group)
            if 0:
                prop = PBEAML.add_op2_data(data_in)
                if pid in op2.properties:
                    # this is a fake PSHELL
                    propi = op2.properties[pid]
                    assert propi.type in ['PBEAM'], propi.get_stats()
                    nproperties -= 1
                    continue
                self._add_op2_property(prop)
        if nproperties:
            op2.card_count['PBEAML'] = nproperties
        return len(data)

    def _read_pbend(self, data: bytes, n: int) -> int:
        """PBEND"""
        self.op2.log.warning('skipping PBEND')
        return len(data)
        op2 = self.op2
        n = op2.reader_geom2._read_dual_card(
            data, n,
            self._read_pbend_nx, self._read_pbend_msc,
            'PBEND', op2._add_methods.add_property_object)
        return n

    def _read_pbend_msc(self, data: bytes, n: int) -> int:
        """
        PBEND

        1 PID     I  Property identification number
        2 MID     I  Material identification number
        3 A       RS Area
        4 I1      RS Area moment of inertia in plane 1
        5 I2      RS Area moment of inertia in plane 2
        6 J       RS Torsional constant
        7 FSI     I  flexibility and stress intensification factors
        8 RM      RS Mean cross-sectional radius of the curved pipe
        9 T       RS Wall thickness of the curved pipe
        10 P      RS Internal pressure
        11 RB     RS Bend radius of the line of centroids
        12 THETAB RS Arc angle of element
        13 C1     RS Stress recovery location at point C in element y-axis
        14 C2     RS Stress recovery location at point C in element z-axis
        15 D1     RS Stress recovery location at point D in element y-axis
        16 D2     RS Stress recovery location at point D in element z-axis
        17 E1     RS Stress recovery location at point E in element y-axis
        18 E2     RS Stress recovery location at point E in element z-axis
        19 F1     RS Stress recovery location at point F in element y-axis
        20 F2     RS Stress recovery location at point F in element z-axis
        21 K1     RS Area factor for shear in plane 1
        22 K2     RS Area factor for shear in plane 2
        23 NSM    RS Nonstructural mass per unit length
        24 RC     RS Radial offset of the geometric centroid
        25 ZC     RS Offset of the geometric centroid
        26 DELTAN  I Radial offset of the neutral axis from the geometric
                     centroid
        """
        op2 = self.op2
        ntotal = 104  # 26*4
        struct1 = Struct(op2._endian + b'2i 4f i 18f f')  # delta_n is a float, not an integer
        nproperties = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nproperties > 0, 'table=%r len=%s' % (op2.table_name, len(data) - n)
        properties = []
        for unused_i in range(nproperties):
            edata = data[n:n+104]
            out = struct1.unpack(edata)
            (pid, mid, area, i1, i2, j, fsi, rm, t, p, rb, theta_b,
             c1, c2, d1, d2, e1, e2, f1, f2, k1, k2, nsm, rc, zc,
             delta_n) = out
            beam_type = fsi

            if (area, rm, t, p) == (0., 0., 0., 0.):
                area = None
                rm = None
                t = None
                p = None
                delta_n = None
                beam_type = 2
            if delta_n == 0:
                #: Radial offset of the neutral axis from the geometric
                #: centroid, positive is toward the center of curvature
                delta_n = None
            pbend = PBEND(pid, mid, beam_type, area, i1, i2, j,
                          c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                          nsm, rc, zc, delta_n, fsi, rm, t, p, rb, theta_b)
            #print(pbend)
            pbend.validate()

            properties.append(pbend)
            n += ntotal
        return n, properties

    def _read_pbend_nx(self, data: bytes, n: int) -> int:
        """
        PBEND

        1 PID     I  Property identification number
        2 MID     I  Material identification number
        3 A       RS Area
        4 I1      RS Area moment of inertia in plane 1
        5 I2      RS Area moment of inertia in plane 2
        6 J       RS Torsional constant
        7 FSI     I  Flexibility and stress intensification factors
        8 RM      RS Mean cross-sectional radius of the curved pipe
        9 T       RS Wall thickness of the curved pipe
        10 P      RS Internal pressure
        11 RB     RS Bend radius of the line of centroids
        12 THETAB RS Arc angle of element
        13 C1     RS Stress recovery location at point C in element y-axis
        14 C2     RS Stress recovery location at point C in element z-axis
        15 D1     RS Stress recovery location at point D in element y-axis
        16 D2     RS Stress recovery location at point D in element z-axis
        17 E1     RS Stress recovery location at point E in element y-axis
        18 E2     RS Stress recovery location at point E in element z-axis
        19 F1     RS Stress recovery location at point F in element y-axis
        20 F2     RS Stress recovery location at point F in element z-axis
        21 K1     RS Area factor for shear in plane 1
        22 K2     RS Area factor for shear in plane 2
        23 NSM    RS Nonstructural mass per unit length
        24 RC     RS Radial offset of the geometric centroid
        25 ZC     RS Offset of the geometric centroid
        26 DELTAN RS Radial offset of the neutral axis from the geometric
                     centroid
        27 SACL   RS Miter spacing at center line.
        28 ALPHA  RS One-half angle between the adjacent miter axis
                     (Degrees).
        29 FLANGE I  For FSI=5, defines the number of flanges attached.
        30 KX     RS For FSI=6, the user defined flexibility factor for the
                  torsional moment.
        31 KY     RS For FSI=6, the user defined flexibility factor for the
                  out-of-plane bending moment.
        32 KZ     RS For FSI=6, the user defined flexbility factor for the
                  in-plane bending moment.
        33 Not used
        """
        op2 = self.op2
        #op2.log.info('geom skipping PBEND in EPT')
        #return len(data)
        ntotal = 132  # 33*4
        struct1 = Struct(op2._endian + b'2i 4f i 21f i 4f')
        nproperties = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nproperties > 0, 'table=%r len=%s' % (op2.table_name, len(data) - n)
        properties = []
        for unused_i in range(nproperties):
            edata = data[n:n+132]
            out = struct1.unpack(edata)
            (pid, mid, area, i1, i2, j, fsi, rm, t, p, rb, theta_b,
             c1, c2, d1, d2, e1, e2, f1, f2, k1, k2, nsm, rc, zc,
             delta_n, unused_sacl, unused_alpha, unused_flange,
             unused_kx, unused_ky, unused_kz, unused_junk,) = out
            beam_type = fsi

            pbend = PBEND(pid, mid, beam_type, area, i1, i2, j,
                          c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                          nsm, rc, zc, delta_n, fsi, rm, t, p, rb, theta_b)
            pbend.validate()
            properties.append(pbend)
            n += ntotal
        return n, properties

# PBMSECT
# PBRSECT

    def _read_pbush(self, data: bytes, n: int) -> int:
        """
        The PBUSH card is different between MSC and NX Nastran.

        DMAP NX 11
        ----------
        NX has 23 fields in NX 11-NX 2019.2 (same as MSC 2005)
        NX has 18 fields in the pre-2001 format

        DMAP MSC 2005
        -------------
        MSC has 23 fields in 2005
        MSC has 18 fields in the pre-2001 format

        DMAP MSC 2016
        -------------
        MSC has 24 fields in 2016.1
        MSC has 18 fields in the pre-2001 format

        DMAP MSC 2021
        -------------
        MSC has 27 fields in 2021

        """
        op2 = self.op2
        card_name = 'PBUSH'
        card_obj = self.op2.pbush
        methods = {
            72 : self._read_pbush_nx_72,  # 72=4*18
            92 : self._read_pbush_msc_92, # 92=4*23
            96 : self._read_pbush_msc_96, # 96=4*24
            108 : self._read_pbush_msc_108, # 108=4*27
        }
        try:
            n = op2.reader_geom2._read_double_card(
                card_name, card_obj, self._add_op2_property,
                methods, data, n)
        except DoubleCardError:
            nx_method = partial(self._read_pbush_nx_72, card_obj)
            msc_method = partial(self._read_pbush_msc_92, card_obj)
            n = op2.reader_geom2._read_dual_card(
                data, n,
                nx_method, msc_method,
                card_name, self._add_op2_property)

        # we're listing nx twice because NX/MSC used to be consistent
        # the new form for MSC is not supported
        #n = self._read_dual_card(data, n, self._read_pbush_nx, self._read_pbush_msc,
                                 #'PBUSH', self._add_op2_property)
        return n

    def _read_pbush_nx_72(self, card_obj: PBUSH, data: bytes, n: int) -> tuple[int, list[PBUSH]]:
        """
        PBUSH(1402,14,37) - 18 fields
        legacy MSC/NX format
        """
        op2 = self.op2
        pbush = op2.pbush
        ntotal = 72 * self.factor  # 18*4
        struct1 = Struct(mapfmt(op2._endian + b'i17f', self.size))
        ndata = len(data) - n
        nproperties = ndata // ntotal
        assert nproperties > 0, f'table={op2.table_name} len={ndata}'
        assert ndata % ntotal == 0, f'table={op2.table_name} leftover = {ndata} % {ntotal} = {ndata % ntotal}'

        n, ints, floats = get_ints_floats(data, n, nproperties, 18,
                                          size=op2.size, endian=op2._endian)
        property_id = ints[:, 0]
        assert property_id.min() > 0, property_id
        k_fields = floats[:, [1, 2, 3, 4, 5, 6]]
        b_fields = floats[:, [7, 8, 9, 10, 11, 12]]
        g1 = floats[:, 13]
        sa = floats[:, 14]
        st = floats[:, 15]
        ea = floats[:, 16]
        et = floats[:, 17]

        ge_fields = g1
        #rcv = [sa, st, ea, et]
        rcv_fields = np.column_stack([sa, st, ea, et])

        _mass = None
        alpha = None
        tref = None
        coincident_length = None
        pbush._save(property_id, k_fields, b_fields, ge_fields, rcv_fields,
                    _mass, alpha, tref, coincident_length)
        #prop.set_from_op2(property_id, k, b, sa, st)
        #op2.log.warning('Skipping PBUSH-72')

        props = []
        #for unused_i in range(nproperties):
            #edata = data[n:n+ntotal]
            #out = struct1.unpack(edata)
            #(pid,
             #k1, k2, k3, k4, k5, k6,
             #b1, b2, b3, b4, b5, b6,
             #g1, sa, st, ea, et) = out
            ##op2.log.debug(out)
            #assert pid > 0, pid
            #g2 = g3 = g4 = g5 = g6 = g1
            #data_in = (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
                       #g1, g2, g3, g4, g5, g6, sa, st, ea, et)
            #prop = PBUSH.add_op2_data(data_in)
            #props.append(prop)
            #n += ntotal
        return n, props

    def _read_pbush_msc_92(self, card_obj: PBUSH, data: bytes, n: int) -> tuple[int, list[PBUSH]]:
        """PBUSH(1402,14,37) - 23 fields

        MSC 2005r2 to <MSC 2016
        """
        op2 = self.op2
        ntotal = 92 * self.factor # 23*4

        ndata = len(data) - n
        nproperties = ndata // ntotal
        assert nproperties > 0, f'table={op2.table_name} len={ndata}'
        assert ndata % ntotal == 0, f'table={op2.table_name} leftover = {ndata} % {ntotal} = {ndata % ntotal}'

        n, ints, floats = get_ints_floats(data, n, nproperties, 23,
                                          size=op2.size, endian=op2._endian)
        property_id = ints[:, 0]
        assert property_id.min() > 0, property_id

        #(pid, k1, k2, k3, k4, k5, k6,
        #      b1, b2, b3, b4, b5, b6,
         #g1, g2, g3, g4, g5, g6, sa, st, ea, et) = out
        k = floats[:, [1,  2,  3,  4,  5,  6]]
        b = floats[:, [7,  8,  9,  10, 11, 12]]
        g = floats[:, [13, 14, 15, 16, 17, 18]]
        sa = floats[:, 14]
        st = floats[:, 15]
        ea = floats[:, 16]
        et = floats[:, 17]
        #prop.set_from_op2(property_id, k, b, g, sa, st)
        op2.log.warning('Skipping PBUSH-92')

        props = []
        #struct1 = Struct(mapfmt(op2._endian + b'i22f', self.size))
        #for unused_i in range(nproperties):
            #edata = data[n:n+ntotal]
            #out = struct1.unpack(edata)
            ##(pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
             ##g1, g2, g3, g4, g5, g6, sa, st, ea, et) = out
            #pid = out[0]
            #assert pid > 0, pid
            #prop = PBUSH.add_op2_data(out)
            #props.append(prop)
            #n += ntotal
        return n, props

    def _read_pbush_msc_96(self, card_obj: PBUSH, data: bytes, n: int) -> tuple[int, list[PBUSH]]:
        """PBUSH(1402,14,37) - 24 fields

        MSC 2016.1? to 2020
        """
        op2 = self.op2
        ntotal = 96 * self.factor # 24*4
        struct1 = Struct(mapfmt(op2._endian + b'i22f f', self.size))

        ndata = len(data) - n
        nentries = ndata // ntotal
        assert nentries > 0, f'table={op2.table_name} len={ndata}'
        assert ndata % ntotal == 0, f'table={op2.table_name} leftover = {ndata} % {ntotal} = {ndata % ntotal}'

        props = []
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
             g1, g2, g3, g4, g5, g6, sa, st, ea, et, mass) = out
            pid = out[0]
            assert pid > 0, pid

            k = [k1, k2, k3, k4, k5, k6]
            b = [b1, b2, b3, b4, b5, b6]
            ge = [g1, g2, g3, g4, g5, g6]
            rcv = [sa, st, ea, et]
            prop = model.add_pbush(pid, k, b, ge, rcv=rcv, mass=mass)
            props.append(prop)
            n += ntotal
        return n, props

    def _read_pbush_msc_108(self, card_obj: PBUSH, data: bytes, n: int) -> tuple[int, list[PBUSH]]:
        """
        PBUSH(1402,14,37) - 27 fields
        MSC 2021 to current

        ints    = (1402, 14, 37, 2, 100000.0, 200000.0, 300000.0, 0.15, 0.25, 0.35, 1000.0, 2000.0, 3000.0, 0.0015, 0.0025, 0.0035, 0,
                   -1577048263, -1577048263, -1577048263, -1577048263, -1577048263, 1065353216, 1065353216, 1065353216, 1065353216, 0, 0, 0, 0)
        floats  = (1402, 14, 37,
                   2, 100000.0, 200000.0, 300000.0, 0.15, 0.25, 0.35, 1000.0, 2000.0, 3000.0, 0.0015, 0.0025, 0.0035, 0.0,
                   -1.7367999061094683e-18, -1.7367999061094683e-18, -1.7367999061094683e-18, -1.7367999061094683e-18, -1.7367999061094683e-18, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0)
        """
        op2 = self.op2
        ntotal = 108 * self.factor # 27*4
        struct1 = Struct(mapfmt(op2._endian + b'i22f 4f', self.size))
        #op2.show_data(data, types='ifs')

        ndata = len(data) - n
        nentries = ndata // ntotal
        assert nentries > 0, f'table={op2.table_name} len={ndata}'
        assert ndata % ntotal == 0, f'table={op2.table_name} leftover = {ndata} % {ntotal} = {ndata % ntotal}'

        props = []
        model = op2 # .model
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)

            (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
             g1, g2, g3, g4, g5, g6, sa, st, ea, et,
             mass, alpha, tref, coinl) = out
            #model.log.warning(f'skipping PBUSH pid={pid} alpha={alpha} tref={tref} coinl={coinl}')
            pid = out[0]
            assert pid > 0, pid
            k = [k1, k2, k3, k4, k5, k6]
            b = [b1, b2, b3, b4, b5, b6]
            ge = [g1, g2, g3, g4, g5, g6]
            rcv = [sa, st, ea, et]
            prop = model.add_pbush(pid, k, b, ge, rcv=rcv, mass=mass,
                                   alpha=alpha, tref=tref, coincident_length=coinl)
            str(prop)
            props.append(prop)
            n += ntotal
        return n, props

    def _read_pbush1d(self, data: bytes, n: int) -> int:
        """
        Record 18 -- PBUSH1D(3101,31,219)

        1  PID    I  Property identification number
        2  K      RS Stiffness
        3  C      RS Viscous Damping
        4  M      RS Mass
        5  ALPHA  RS Temperature coefficient
        6  SA     RS Stress recovery coefficient
        7  EA/SE  RS Strain recovery coefficient

        8  TYPEA  I  Shock data type:0=Null, 1=Table, 2=Equation
        9  CVT    RS Coefficient of translation velocity tension
        10 CVC    RS Coefficient of translation velocity compression
        11 EXPVT  RS Exponent of velocity tension
        12 EXPVC  RS Exponent of velocity compression
        13 IDTSU  I  TABLEDi or DEQATN entry identification number for scale factor vs displacement
        14 IDTCU  I  DEQATN entry identification number for scale factor vs displacement
        15 IDTSUD I  DEQATN entry identification number for derivative tension
        16 IDCSUD I  DEQATN entry identification number for derivative compression

        17 TYPES  I  Spring data type: 0=Null, 1=Table, 2=Equation
        18 IDTS   I  TABLEDi or DEQATN entry identification number for tension compression
        19 IDCS   I  DEQATN entry identification number for compression
        20 IDTDU  I  DEQATN entry identification number for scale factor vs displacement
        21 IDCDU  I  DEQATN entry identification number for force vs displacement

        22 TYPED  I  Damper data type: 0=Null, 1=Table, 2=Equation
        23 IDTD   I  TABLEDi or DEQATN entry identification number for tension compression
        24 IDCD   I  DEQATN entry identification number for compression
        25 IDTDV  I  DEQATN entry identification number for scale factor versus velocity
        26 IDCDV  I  DEQATN entry identification number for force versus velocity

        27 TYPEG  I  General data type: 0=Null, 1=Table, 2=Equation
        28 IDTG   I  TABLEDi or DEQATN entry identification number for tension compression
        29 IDCG   I  DEQATN entry identification number for compression
        30 IDTDU  I  DEQATN entry identification number for scale factor versus displacement
        31 IDCDU  I  DEQATN entry identification number for force versus displacement
        32 IDTDV  I  DEQATN entry identification number for scale factor versus velocity
        33 IDCDV  I  DEQATN entry identification number for force vs velocity

        34 TYPEF  I  Fuse data type: 0=Null, 1=Table
        35 IDTF   I  TABLEDi entry identification number for tension
        36 IDCF   I  TABLEDi entry identification number for compression

        37 UT     RS Ultimate tension
        38 UC     RS Ultimate compression
        """
        op2 = self.op2
        type_map = {
            0 : None,  # NULL
            1 : 'TABLE',
            2 : 'EQUAT',
        }
        ntotal = 152 * self.factor  # 38*4
        struct1 = Struct(mapfmt(op2._endian + b'i 6f i 4f 24i 2f', self.size))
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            (pid, k, c, m, unused_alpha, sa, se,
             typea, cvt, cvc, expvt, expvc, idtsu, idtcu, idtsud, idcsud,
             types, idts, idcs, idtdus, idcdus,
             typed, idtd, idcd, idtdvd, idcdvd,
             typeg, idtg, idcg, idtdug, idcdug, idtdvg, idcdvg,
             typef, idtf, idcf,
             unused_ut, unused_uc) = out
            #  test_op2_other_05
            #pbush1d, 204, 1.e+5, 1000., , , , , , +pb1
            #+pb1, spring, table, 205, , , , , , +pb2
            #+pb2, damper, table, 206
            #pid=204 k=100000.0 c=1000.0 m=0.0 sa=nan se=nan


            msg = f'PBUSH1D pid={pid} k={k} c={c} m={m} sa={sa} se={se}'
            optional_vars = {}
            typea_str = type_map[typea]
            types_str = type_map[types]
            typed_str = type_map[typed]
            unused_typeg_str = type_map[typeg]
            unused_typef_str = type_map[typef]

            if min([typea, types, typed, typeg, typef]) < 0:
                raise RuntimeError(f'typea={typea} types={types} typed={typed} typeg={typeg} typef={typef}')
            if typea in [1, 2]:  # SHOCKA?
                #pbush1d, 204, 1.e+5, 1000., , , , , , +pb4
                #+pb4, shocka, table, 1000., , 1., , 214, , +pb41
                #+pb41, spring, table, 205

                idts = idtsu # if typea_str == 'TABLE' else 0
                idets = idtsu # if typea_str == 'EQUAT' else 0
                optional_vars['SHOCKA'] = [typea_str, cvt, cvc, expvt, expvc,
                                           idts, idets, idtcu, idtsud, idcsud]
                #(shock_type, shock_cvt, shock_cvc, shock_exp_vt, shock_exp_vc,
                 #shock_idts, shock_idets, shock_idecs, shock_idetsd, shock_idecsd
                #)
                #print('shock_idts, shock_idets', typea_str, idtsu, idtsu)
                msg += (
                    f'  SHOCKA type={typea} cvt={cvt} cvc={cvc} expvt={expvt} expvc={expvc}\n'
                    f'    idtsu={idtsu} (idts={idts} idets={idets}) idtcu={idtcu} idtsud={idtsud} idcsud={idcsud}')
            if types in [1, 2]: # SPRING: Spring data type: 0=Null, 1=Table, 2=Equation
                #(spring_type, spring_idt, spring_idc, spring_idtdu, spring_idcdu) = values
                # SPRING, TYPE IDT IDC IDTDU IDCDU
                optional_vars['SPRING'] = [types_str, idts, idcs, idtdus, idcdus]
                msg += f'  SPRING type={types} idt={idts} idc={idcs} idtdu={idtdus} idcdu={idcdus}'
            if typed in [1, 2]: # Damper data type: 0=Null, 1=Table, 2=Equation
                optional_vars['DAMPER'] = [typed_str, idtd, idcd, idtdvd, idcdvd]
                msg += f'  DAMPER type={typed} idt={idtd} idc={idtd} idtdv={idtdvd} idcdv={idcdvd}'
            if typeg in [1, 2]: # general, GENER?: 0=Null, 1=Table 2=Equation
                # C:\NASA\m4\formats\git\examples\move_tpl\ar29scbt.bdf
                #pbush1d, 206, 1.e+3, 10., , , , , , +pb6
                #+pb6, gener, equat, 315, , 3015, , 3016
                msg += f'  GENER  type={typeg} idt={idtg} idc={idcg} idtdu={idtdug} idcdu={idcdug} idtdv={idtdvg} idcdv={idcdvg}'
                optional_vars['GENER'] = [idtg, idcg, idtdug, idcdug, idtdvg, idcdvg]
            if typef in [1, 2]: # Fuse data type: 0=Null, 1=Table
                raise NotImplementedError(f'typef={typef} idtf={idtf} idcf={idcf}')

            if op2.is_debug_file:
                op2.binary_debug.write(msg)

            pbush1d = op2.add_pbush1d(pid, k=k, c=c, m=m, sa=sa, se=se,
                                      optional_vars=optional_vars,)
            str(pbush1d)
            n += ntotal
        op2.card_count['PBUSH1D'] = nentries
        return n

    #def _read_pbusht(self, data: bytes, n: int) -> int:
        #"""reads the PBUSHT(702, 7, 38)"""
        #n, props = self._read_pbusht_nx(data, n)
        #for prop in props:
            ##print(prop)
            #op2._add_pbusht_object(prop)
        #return n

    def _read_pbusht(self, data: bytes, n: int) -> int:
        """
        NX 12 / MSC 2005
        Word Name Type Description
        1 PID       I Property identification number
        2 TKID(6)   I TABLEDi entry identification numbers for stiffness
        8 TBID(6)   I TABLEDi entry identification numbers for viscous damping
        14 TGEID(6) I TABLEDi entry identification number for structural damping
        20 TKNID(6) I TABLEDi entry identification numbers for force versus deflection

        old style
        Word Name Type Description
        1 PID       I Property identification number
        2 TKID(6)   I TABLEDi entry identification numbers for stiffness
        8 TBID(6)   I TABLEDi entry identification numbers for viscous damping
        14 TGEID    I TABLEDi entry identification number for structural damping
        15 TKNID(6) I TABLEDi entry IDs for force versus deflection
        """
        op2 = self.op2
        card_name = 'PBUSHT'
        card_obj = op2.pbusht
        methods = {
            80 : self._read_pbusht_80,
            100 : self._read_pbusht_100,
            136 : self._read_pbusht_136,
        }
        try:
            n = op2.reader_geom2._read_double_card(
                card_name, card_obj, None,
                methods, data, n)
        except DoubleCardError:
            raise
            op2.log.warning(f'try-except {card_name}')
            #n = self._read_split_card(data, n,
                                      #self._read_cquad8_current, self._read_cquad8_v2001,
                                      #card_name, self.add_op2_element)
        #nelements = op2.card_count['CQUAD8']
        #op2.log.debug(f'nCQUAD8 = {nelements}')

        #n = self._read_dual_card(data, n, self._read_ctriax_8, self._read_ctriax_9,
                                 #'CTRIAX', self.add_op2_element)
        return n

    def _read_pbusht_nx_old(self, data: bytes, n: int) -> int:
        op2 = self.op2
        #op2.show_data(data[12:])
        ndata = (len(data) - n) // self.factor

        if ndata % 100 == 0 and ndata % 80 == 0:
            op2.log.warning(f"skipping PBUSHT in EPT because nfields={ndata//4}, which is "
                             'nproperties*25 or nproperties*20')
            return len(data), []
        if ndata % 100 == 0:
            n, props = self._read_pbusht_100(data, n)
        elif ndata % 80 == 0:
            n, props = self._read_pbusht_80(data, n)
        else:
            # C:\MSC.Software\msc_nastran_runs\mbsh14.op2
            # ints = (1,
            #         51, 51, 0, 0, 0, 0,
            #         61, 61, 0, 0, 0, 0,
            #         0,  0,  0, 0, 0, 0,
            #         0, '', '', 0, 0, '', '', 0, 0, 925353388, 0, 0, 0, 0, 0,
            #         7,
            #         51, 51, 0, 0, 0, 0,
            #         61, 61, 0, 0, 0, 0,
            #         0,  0,  0, 0, 0, 0,
            #         0, '', '', 0, 0, '', '', 0, 0, 925353388, 0, 0, 0, 0, 0)
            # strings = (b"1 51 51 \x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00=\x00\x00\x00=\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00        \x00\x00\x00\x00\x00\x00\x00\x00        \x00\x00\x00\x00\x00\x00\x00\x00\xac\xc5'7\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x07\x00\x00\x003\x00\x00\x003\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00=\x00\x00\x00=\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00        \x00\x00\x00\x00\x00\x00\x00\x00        \x00\x00\x00\x00\x00\x00\x00\x00\xac\xc5'7\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00",)
            # ints    = (1, 51, 51, 0,   0,   0,   0,   61, 61, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   '    ', '    ', 0,   0,   '    ', '    ', 0,   0,   1e-5, 0,   0,   0,   0  , 0,
            #
            # 7, 51, 51, 0,   0,   0,   0,   61, 61, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, '    ', '    ', 0, 0, '    ', '    ', 0, 0, 1e-5, 0, 0, 0, 0, 0)
            #op2.show_data(data[n:], types='is')
            raise NotImplementedError('You have blank lines in your PBUSHT')
        return n, props

    def _read_pbusht_80(self, card_obj, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 PID     I Property identification number
        2 TKID(6) I TABLEDi entry identification numbers for stiffness
        8 TBID(6) I TABLEDi entry identification numbers for viscous damping
        14 TGEID  I TABLEDi entry identification number for structural damping
        15 TKNID(6) I TABLEDi entry identification numbers for force versus deflection
        16,17,18,19,20
        ???
        """
        op2 = self.op2
        ntotal = 80 * self.factor
        struct1 = Struct(op2._endian + b'20i')
        nproperties = (len(data) - n) // ntotal
        assert nproperties > 0, 'table=%r len=%s' % (op2.table_name, len(data) - n)

        n, ints = get_ints(
            data, n, nproperties, 20,
            size=op2.size, endian=op2._endian)

        #(pid,
         #k1, k2, k3, k4, k5, k6,
         #b1, b2, b3, b4, b5, b6,
         #g1, g2, g3, g4, g5, g6,
         #n1, n2, n3, n4, n5, n6) = out
        property_id = ints[:, 0]
        k_tables = ints[:, [1, 2, 3, 4, 5, 6]]
        b_tables = ints[:, [7, 8, 9, 10, 11, 12]]

        ge_tablesi = ints[:, [13]]
        ge_tables = np.tile(ge_tablesi, 6)
        assert ge_tables.shape == (nproperties, 6), ge_tables.shape
        kn_tables = ints[:, [14, 15, 16, 17, 18, 19]]
        #print(k_tables)

        #print(ints[:, 25:])
        prop = op2.pbusht
        prop._save(property_id, k_tables, b_tables, ge_tables, kn_tables)
        prop.write()


        props = []
        #for unused_i in range(nproperties):
            #edata = data[n:n+ntotal]
            #out = struct1.unpack(edata)
            ##(pid,
             ##k1, k2, k3, k4, k5, k6,
             ##b1, b2, b3, b4, b5, b6,
             ##g1, sa, st, ea, et) = out
            #(pid,
             #k1, k2, k3, k4, k5, k6,
             #b1, b2, b3, b4, b5, b6,
             #g1,
             #n1, n2, n3, n4, n5, n6) = out
            #g2 = g3 = g4 = g5 = g6 = g1
            #k_tables = [k1, k2, k3, k4, k5, k6]
            #b_tables = [b1, b2, b3, b4, b5, b6]
            #ge_tables = [g1, g2, g3, g4, g5, g6]
            #kn_tables = [n1, n2, n3, n4, n5, n6]
            #prop = PBUSHT(pid, k_tables, b_tables, ge_tables, kn_tables)
            #props.append(prop)
            #n += ntotal
        return n, props

    def _read_pbusht_100(self, card_obj, data: bytes, n: int) -> int:
        op2 = self.op2
        props = []
        ntotal = 100 * self.factor
        nproperties = (len(data) - n) // ntotal
        assert nproperties > 0, 'table=%r len=%s' % (op2.table_name, len(data) - n)

        n, ints = get_ints(
            data, n, nproperties, 25,
            size=op2.size, endian=op2._endian)

        #(pid,
         #k1, k2, k3, k4, k5, k6,
         #b1, b2, b3, b4, b5, b6,
         #g1, g2, g3, g4, g5, g6,
         #n1, n2, n3, n4, n5, n6) = out
        property_id = ints[:, 0]
        k_tables = ints[:, [1, 2, 3, 4, 5, 6]]
        b_tables = ints[:, [7, 8, 9, 10, 11, 12]]
        ge_tables = ints[:, [13, 14, 15, 16, 17, 18]]
        kn_tables = ints[:, [19, 20, 21, 22, 23, 24]]
        #print(k_tables)

        #print(ints[:, 25:])
        prop = op2.pbusht
        prop._save(property_id, k_tables, b_tables, ge_tables, kn_tables)
        prop.write()

        #struct1 = Struct(mapfmt(op2._endian + b'25i', self.size))
        #for unused_i in range(nproperties):
            #edata = data[n:n+ntotal]
            #out = struct1.unpack(edata)
            #(pid,
             #k1, k2, k3, k4, k5, k6,
             #b1, b2, b3, b4, b5, b6,
             #g1, g2, g3, g4, g5, g6,
             #n1, n2, n3, n4, n5, n6) = out
            #k_tables = [k1, k2, k3, k4, k5, k6]
            #b_tables = [b1, b2, b3, b4, b5, b6]
            #ge_tables = [g1, g2, g3, g4, g5, g6]
            #kn_tables = [n1, n2, n3, n4, n5, n6]
            #prop = PBUSHT(pid, k_tables, b_tables, ge_tables, kn_tables)
            #props.append(prop)
            #n += ntotal
        return n, props

    def _read_pbusht_136(self, card_obj, data: bytes, n: int) -> int:
        r"""not 100%

        1  PID           I Property identification number
        2  TKID(6)       I TABLEDi entry identification numbers for stiffness
        8  TBID(6)       I TABLEDi entry identification numbers for viscous damping
        14 TGEID(6)      I TABLEDi entry identification number for structural damping
        20 TKNID(6)      I TABLEDi entry IDs for force vs. deflection
        26 FDC(2)    CHAR4 Force deflection curve rule
        28 FUSE          I Failure level
        29 DIR           I Fuse direction
        30 OPTION(2) CHAR4 Failure mode
        32 LOWER        RS Lower failure bound
        33 UPPER        RS Upper failure bound
        34 FRATE        RS FACTOR of scales the stiffness
        35 LRGR          I Controls large rotation
        36 UNDEF(4)        none

        # C:\MSC.Software\msc_nastran_runs\mbsh14.op2
        PBUSHT	1	 K	51	51
                 B	61	61
        PBUSHT	7	 K	51	51
                 B	61	61

        538976288 = '    '
        ints    = (
            702, 7, 38,
            1, (51, 51, 0, 0, 0, 0), (61, 61, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0), 0, 538976288, 538976288, 0, 0, 538976288, 538976288, 0, 0, 925353388, 0, 0, 0, 0, 0,
            7, (51, 51, 0, 0, 0, 0), (61, 61, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0), 0, 538976288, 538976288, 0, 0, 538976288, 538976288, 0, 0, 925353388, 0, 0, 0, 0, 0)
        floats  = (
            702, 7, 38,
            1, 51, 51, 0.0, 0.0, 0.0, 0.0, 61, 61, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 538976288, 538976288, 0.0, 0.0, 538976288, 538976288, 0.0, 0.0, 1.e-7, 0.0, 0.0, 0.0, 0.0, 0.0,
            7, 51, 51, 0.0, 0.0, 0.0, 0.0, 61, 61, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 538976288, 538976288, 0.0, 0.0, 538976288, 538976288, 0.0, 0.0, 1.e-7, 0.0, 0.0, 0.0, 0.0, 0.0)
        """
        op2 = self.op2
        props = []
        ntotal = 136 * self.factor              #  k  b  g  n  fdc
        #struct1 = Struct(mapfmt(op2._endian + b'i 6i 6i 6i 6i 4s  2i i 5i', self.size))
        nproperties = (len(data) - n) // ntotal
        assert nproperties > 0, 'table=%r len=%s' % (op2.table_name, len(data) - n)

        n, ints, strings = get_ints_strings(
            data, n, nproperties, 34,
            size=op2.size, endian='<')
        property_id = ints[:, 0]
        k_tables = ints[:, [1, 2, 3, 4, 5, 6]]
        b_tables = ints[:, [7, 8, 9, 10, 11, 12]]
        ge_tables = ints[:, [13, 14, 15, 16, 17, 18]]
        kn_tables = ints[:, [19, 20, 21, 22, 23, 24]]
        word1 = strings[:, 25]
        a_word2_c = ints[:, [26, 27, 28]]
        other = ints[:, 29:]
        assert other.shape == (nproperties, 5), other.shape
        assert other.sum() == 0, other
        #print(word1)
        #print(a_word2_c)

        k_tables_blank = strings[:, [1, 2, 3, 4, 5, 6]]
        b_tables_blank = strings[:, [7, 8, 9, 10, 11, 12]]
        ge_tables_blank = strings[:, [13, 14, 15, 16, 17, 18]]
        kn_tables_blank = strings[:, [19, 20, 21, 22, 23, 24]]

        k_tables[np.where(k_tables_blank == b'    ')] = 0
        b_tables[np.where(b_tables_blank == b'    ')] = 0
        ge_tables[np.where(ge_tables_blank == b'    ')] = 0
        kn_tables[np.where(kn_tables_blank == b'    ')] = 0
        #print(i)
        #print(k_tables)

        #print(ints[:, 25:])
        #print(strings[:, 25:])
        prop = op2.pbusht
        prop._save(property_id, k_tables, b_tables, ge_tables, kn_tables)
        prop.write()

        #ntotal = 136 * self.factor              #  k  b  g  n  fdc
        #struct1 = Struct(mapfmt(op2._endian + b'i 6i 6i 6i 6i 4s  2i i 5i', self.size))
        #for unused_i in range(nproperties):
            #edata = data[n:n+ntotal]
            #out = struct1.unpack(edata)
            #(pid,
             #k1, k2, k3, k4, k5, k6,
             #b1, b2, b3, b4, b5, b6,
             #g1, g2, g3, g4, g5, g6,
             #n1, n2, n3, n4, n5, n6,
             #word1, a, word2, c, *other) = out
            #assert len(other) == 5, other


            #k_tables = [ki if ki != 538976288 else 0
                        #for ki in [k1, k2, k3, k4, k5, k6]]

            #b_tables = [bi if bi != 538976288 else 0
                        #for bi in [b1, b2, b3, b4, b5, b6]]
            #ge_tables = [gei if gei != 538976288 else 0
                        #for gei in [g1, g2, g3, g4, g5, g6]]
            #kn_tables = [kni if kni != 538976288 else 0
                        #for kni in [n1, n2, n3, n4, n5, n6]]
            #op2.log.warning(
                #f'PBUSHT: pid={pid} '
                #f'k={k_tables} '
                #f'b={b_tables} '
                #f'ge={ge_tables} '
                #f'n={kn_tables} ' +
                #'words=' + str([word1, a, word2, c]) +
                #f' other={other}')
            #assert sum(other) == 0, other
            #prop = PBUSHT(pid, k_tables, b_tables, ge_tables, kn_tables)
            #props.append(prop)
            #n += ntotal
        return n, props

    def _read_pcomp(self, data: bytes, n: int) -> int:
        r"""
        PCOMP(2706,27,287) - the marker for Record 22

        standard:
          EPTS; 64-bit: C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\cqrdbxdra3lg.op2

        optistruct:
          ints    = (2706, 27, 287,
                     5,
                     3, -2.75, 0, 0, 1, 0,   0,
                     2,  0.25, 0, 2, # why is sout=2?
                     3,  5.0,  0, 3, # why is sout=3?
                     2,  0.25, 0, 2, # why is sout=2?

                     6, 5, -3.0, 0, 0, 1, 0, 0,
                     2, 0.25, 0, 2,
                     2, 0.25, 0, 2,
                     3, 5.0,  0, 3,
                     2, 0.25, 0, 2,
                     2, 0.25, 0, 2, 7, 7, -1068498944, 0, 0, 1, 0, 0, 2, 0.25, 0, 2, 2, 0.25, 0, 2, 2, 0.25, 0, 2, 3, 5.0, 0, 3, 2, 0.25, 0, 2, 2, 0.25, 0, 2, 2, 0.25, 0, 2)
          floats  = (2706, 27, 287,
                     5, 3, -2.75, 0.0, 0.0, 1, 0.0, 0.0, 2, 0.25, 0.0, 2, 3, 5.0, 0.0, 3, 2, 0.25, 0.0, 2, 6, 5, -3.0, 0.0, 0.0, 1, 0.0, 0.0, 2, 0.25, 0.0, 2, 2, 0.25, 0.0, 2, 3, 5.0, 0.0, 3, 2, 0.25, 0.0, 2, 2, 0.25, 0.0, 2, 9.80908925027372e-45, 9.80908925027372e-45, -3.25, 0.0, 0.0, 1, 0.0, 0.0, 2, 0.25, 0.0, 2, 2, 0.25, 0.0, 2, 2, 0.25, 0.0, 2, 3, 5.0, 0.0, 3, 2, 0.25, 0.0, 2, 2, 0.25, 0.0, 2, 2, 0.25, 0.0, 2)
        """
        op2 = self.op2
        if self.size == 4:
            n2, props = self._read_pcomp_32_bit(data, n)
            nproperties = len(props)
            for prop in props:
                self._add_op2_property(prop)
            op2.card_count['PCOMP'] = nproperties
        else:
            n2 = op2.reader_geom2._read_dual_card(
                data, n, self._read_pcomp_32_bit,
                self._read_pcomp_64_bit,
                'PCOMP', self._add_op2_property)
        op2.pcomp.parse_cards()
        return n2

    def _read_pcomp_64_bit(self, data: bytes, n: int) -> tuple[int, list[PCOMP]]:
        r"""
        PCOMP(2706,27,287) - the marker for Record 22

        1  PID   I  Property identification number
        2  N(C)  I  Number of plies
        3  Z0    RS Distance from the reference plane to the bottom surface
        4  NSM   RS Nonstructural mass per unit area
        5  SB    RS Allowable shear stress of the bonding material
        6  FT    I  Failure theory
        7  TREF  RS Reference temperature
        8  GE    RS Damping coefficient

        9  MID   I  Material identification number
        10 T     RS Thicknesses of the ply
        11 THETA RS Orientation angle of the longitudinal direction of the ply
        12 SOUT  I Stress or strain output request of the ply
        Words 9 through 12 repeat N times

        TODO:
           64-bit bug: why is the number of plies 0???

          doubles (float64) = (
          1, 0.0, 1.7368e-18, 0.0, 1.0, 1.5e-323, 0.0, 0.0,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
          -1, -1, -1, -1,
          21, 0.0, 1.7368e-18, 0.0, 1.0, 1.5e-323, 0.0, 0.0,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
          -1, -1, -1, -1)
          long long (int64) = (
          1, 0,   1.7368e-18, 0,   1.0, 3, 0, 0, 1, 4592590756007337001, 0, 1,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
          -1, -1, -1, -1,
          21, 0, 4341475431749739292, 0, 4607182418800017408, 3, 0, 0,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
          -1, -1, -1, -1)

          doubles (float64) = (5e-324, 0.0, -0.005, 0.0, 0.0, 0.0, 0.0, 0.0,
                               4e-323, 0.005, 0.0, 5e-324,
                               4e-323, 0.005, 0.0, 5e-324,
                               nan, nan, nan, nan)
          long long (int64) = (1, 0, -4650957407178058629, 0, 0, 0, 0, 0,
                                  8, 4572414629676717179, 0, 1,
                                  8, 4572414629676717179, 0, 1,
                               -1, -1, -1, -1)

        C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\dbxdr12lg.op2
        data = (3321, 2, -0.5, 0.0, 1.0, 4, 0.0, 0.0,
                3, 0.5, 0, 1,
                3, 0.5, 0, 1)
        """
        op2 = self.op2
        op2.to_nx(' because PCOMP-64 was found')
        nproperties = 0
        s1 = Struct(mapfmt(op2._endian + b'2i3fi2f', self.size))
        ntotal1 = 32 * self.factor
        s2 = Struct(mapfmt(op2._endian + b'i2fi', self.size))

        four_minus1 = Struct(mapfmt(op2._endian + b'4i', self.size))
        ndata = len(data)
        ntotal2 = 16 * self.factor
        props = []
        while n < (ndata - ntotal1):
            out = s1.unpack(data[n:n+ntotal1])
            (pid, nlayers, z0, nsm, sb, ft, tref, ge) = out
            assert pid > 0
            if op2.binary_debug:
                op2.binary_debug.write(f'PCOMP pid={pid} nlayers={nlayers} z0={z0} nsm={nsm} '
                                        f'sb={sb} ft={ft} Tref={tref} ge={ge}')
            assert isinstance(nlayers, int), out
            #print(f'PCOMP pid={pid} nlayers={nlayers} z0={z0} nsm={nsm} '
                  #f'sb={sb} ft={ft} Tref={tref} ge={ge}')
            n += ntotal1

            # None, 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE', 'NO'
            is_symmetrical = 'NO'
            #if nlayers < 0:
                #is_symmetrical = 'SYM'
                #nlayers = abs(nlayers)

            mids = []
            T = []
            thetas = []
            souts = []
            edata2 = data[n:n+ntotal2]
            idata = four_minus1.unpack(edata2)
            while idata != (-1, -1, -1, -1):
                (mid, t, theta, sout) = s2.unpack(edata2)
                mids.append(mid)
                T.append(t)
                thetas.append(theta)
                souts.append(sout)
                if op2.is_debug_file:
                    op2.binary_debug.write(f'      mid={mid} t={t} theta={theta} sout={sout}\n')
                n += ntotal2
                #print(f'      mid={mid} t={t} theta={theta} sout={sout}')
                edata2 = data[n:n+ntotal2]
                if n == ndata:
                    op2.log.warning('  no (-1, -1, -1, -1) flag was found to close the PCOMPs')
                    break
                idata = four_minus1.unpack(edata2)

            if self.size == 4:
                assert 0 < nlayers < 400, 'pid=%s nlayers=%s z0=%s nms=%s sb=%s ft=%s Tref=%s ge=%s' % (
                    pid, nlayers, z0, nsm, sb, ft, tref, ge)
            else:
                assert nlayers == 0, nlayers
                nlayers = len(mids)

            data_in = [
                pid, z0, nsm, sb, ft, tref, ge,
                is_symmetrical, mids, T, thetas, souts]

            lam = ''
            if is_symmetrical:
                lam = 'SYM'

            op2.add_pcomp(pid, mids, T, thetas=thetas, souts=souts,
                          nsm=nsm, sb=sb, ft=ft, tref=tref, ge=ge, lam=lam,
                          z0=z0, comment='')
            #prop = PCOMP.add_op2_data(data_in)
            nproperties += 1
            n += ntotal2
            #props.append(prop)
        return n, props

    def _read_pcomp_32_bit(self, data: bytes, n: int) -> tuple[int, list[PCOMP]]:  # pragma: no cover
        """PCOMP(2706,27,287) - the marker for Record 22"""
        op2 = self.op2
        model = op2
        nproperties = 0
        s1 = Struct(mapfmt(op2._endian + b'2i3fi2f', self.size))
        ntotal1 = 32 * self.factor
        s2 = Struct(mapfmt(op2._endian + b'i2fi', self.size))

        sout_mapper = {
            0: 'NO',
            1: 'YES', }
        ndata = len(data)
        ntotal2 = 16 * self.factor
        props = []
        while n < (ndata - ntotal1):
            out = s1.unpack(data[n:n+ntotal1])
            (pid, nlayers, z0, nsm, sb, ft_int, tref, ge) = out
            assert pid > 0

            if op2.binary_debug:
                op2.binary_debug.write(f'PCOMP pid={pid} nlayers={nlayers} z0={z0} nsm={nsm} '
                                        f'sb={sb} ft={ft_int} Tref={tref} ge={ge}')
            assert isinstance(nlayers, int), out
            #print(f'PCOMP pid={pid} nlayers={nlayers} z0={z0} nsm={nsm} '
                  #f'sb={sb} ft={ft} Tref={tref} ge={ge}')
            n += ntotal1

            mids = []
            T = []
            thetas = []
            souts = []

            # None, 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE', 'NO'
            is_symmetrical = 'NO'
            if nlayers < 0:
                is_symmetrical = 'SYM'
                nlayers = abs(nlayers)
            assert nlayers > 0, out

            assert 0 < nlayers < 400, 'pid=%s nlayers=%s z0=%s nsm=%s sb=%s ft=%s Tref=%s ge=%s' % (
                pid, nlayers, z0, nsm, sb, ft_int, tref, ge)

            if op2.is_debug_file:
                op2.binary_debug.write('    pid=%s nlayers=%s z0=%s nsm=%s sb=%s ft=%s Tref=%s ge=%s\n' % (
                    pid, nlayers, z0, nsm, sb, ft_int, tref, ge))
            #if op2._nastran_format == 'optistruct':
                #print('    pid=%s nlayers=%s z0=%s nsm=%s sb=%s ft=%s Tref=%s ge=%s' % (
                    #pid, nlayers, z0, nsm, sb, ft, tref, ge))
            for unused_ilayer in range(nlayers):
                (mid, t, theta, sout) = s2.unpack(data[n:n+ntotal2])
                if op2._nastran_format == 'optistruct':
                    #print(f'      mid={mid} t={t} theta={theta} sout={sout}')
                    if sout in [2, 3]: # TODO: Why is this 2/3?
                        sout = 1 # YES

                mids.append(mid)
                assert mid > 0

                T.append(t)
                thetas.append(theta)
                souts.append(sout)
                if op2.is_debug_file:
                    op2.binary_debug.write(f'      mid={mid} t={t} theta={theta} sout={sout}\n')
                n += ntotal2

            data_in = [
                pid, z0, nsm, sb, ft_int, tref, ge,
                is_symmetrical, mids, T, thetas, souts]
            lam = 'SYM' if is_symmetrical else ''
            try:
                ft = map_failure_theory_int(ft_int)
            except NotImplementedError:  # pragma: no cover
                raise RuntimeError(f'unsupported ft.  pid={pid} ft={ft_int!r}.'
                               f'\nPCOMP = {data_in}')

            souts_str = [sout_mapper[sout] for sout in souts]
            prop = model.add_pcomp(
                pid, mids, T, thetas, souts_str,
                nsm=nsm, sb=sb, ft=ft, tref=tref, ge=ge, lam=lam,
                z0=z0)
            #print(prop)
            props.append(prop)
            nproperties += 1
        return n, props

    def _read_pcompg(self, data: bytes, n: int) -> int:
        """
        PCOMP(2706,27,287)

        1 PID      I  Property identification number
        2 LAMOPT   I  Laminate option
        3 Z0       RS Distance from the reference plane to the bottom surface
        4 NSM      RS Nonstructural mass per unit area
        5 SB       RS Allowable shear stress of the bonding material
        6 FT       I  Failure theory
        7 TREF     RS Reference temperature
        8 GE       RS Damping coefficient

        9  GPLYIDi I  Global ply IDs.
        10 MID     I  Material identification number
        11 T       RS Thicknesses of the ply
        12 THETA   RS Orientation angle of the longitudinal direction of the ply
        13 SOUT    I  Stress or strain output request of the ply
        Words 9 through 13 repeat N times (until -1, -1, -1, -1, -1 as Nplies doesn't exist...)

        float = (15006, 150, 604,
                 5, 0.0, 1.7368e-18, 0.0, 0.0, 0.0, 20.0, 0.0,
                     5e-324, 5e-324, 2.0, 0.0, 0.0,
                     1e-323, 1e-323, 3.0, 0.0, 0.0,
                     1.5e-323, 1e-323, 3.0, 0.0, 0.0,
                     2e-323, 5e-324, 2.0, 0.0, 0.0,
                     nan, nan, nan, nan, nan)
        int   = (15006, 150, 604,
                 5, 0,   1.7368e-18, 0,   0,   0,   20.0, 0,
                     1, 1, 4611686018427387904, 0, 0,
                     2, 2, 4613937818241073152, 0, 0,
                     3, 2, 4613937818241073152, 0, 0,
                     4, 1, 4611686018427387904, 0, 0,
                     -1, -1, -1, -1, -1)

        """
        op2 = self.op2  # type: OP2Geom
        #op2.log.warning('skipping PCOMPG')
        #return len(data)
        nproperties = 0
        s1 = Struct(mapfmt(op2._endian + b'2i 3f i 2f', self.size))
        s2 = Struct(mapfmt(op2._endian + b'2i 2f i', self.size))
        struct_i5 = Struct(mapfmt(op2._endian + b'5i', self.size))

        # lam - SYM, MEM, BEND, SMEAR, SMCORE, None
        lam_map = {
            0 : None,
            # MEM
            # BEND
            # SMEAR
            # SMCORE
        }

        # ft - HILL, HOFF, TSAI, STRN, None
        ft_map = {
            0 : None,
            # HILL
            # HOFF
            3 : 'TSAI',
            # STRN
        }
        # sout - YES, NO
        sout_map = {
            0 : 'NO',
            1 : 'YES',
        }
        pcomp_pids_to_remove = set()
        ndata = len(data)
        #op2.show_data(data, types='qd')
        ntotal1 = 32 * self.factor
        ntotal2 = 20 * self.factor
        while n < (ndata - ntotal1):
            out = s1.unpack(data[n:n+ntotal1])
            (pid, lam_int, z0, nsm, sb, ft_int, tref, ge) = out
            if op2.binary_debug:
                op2.binary_debug.write(f'PCOMPG pid={pid} lam_int={lam_int} z0={z0} nsm={nsm} '
                                        f'sb={sb} ft_int={ft_int} tref={tref} ge={ge}')
            #print(f'PCOMPG pid={pid} lam_int={lam_int} z0={z0} nsm={nsm} sb={sb} '
                  #f'ft_int={ft_int} tref={tref} ge={ge}')
            assert isinstance(lam_int, int), out
            assert pid > -1, out
            n += ntotal1

            mids = []
            thicknesses = []
            thetas = []
            souts = []
            global_ply_ids = []

            # None, 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE', 'NO'
            #is_symmetrical = 'NO'
            #if nlayers < 0:
                #is_symmetrical = 'SYM'
                #nlayers = abs(nlayers)
            #assert nlayers > 0, out

            #assert 0 < nlayers < 400, 'pid=%s nlayers=%s z0=%s nms=%s sb=%s ft=%s tref=%s ge=%s' % (
                #pid, nlayers, z0, nsm, sb, ft, tref, ge)

            #if op2.is_debug_file:
                #op2.binary_debug.write('    pid=%s nlayers=%s z0=%s nms=%s sb=%s ft=%s tref=%s ge=%s\n' % (
                    #pid, nlayers, z0, nsm, sb, ft, tref, ge))
            ilayer = 0
            while ilayer < 1000:
                ints5 = struct_i5.unpack(data[n:n+ntotal2])
                if ints5 == (-1, -1, -1, -1, -1):
                    if op2.is_debug_file:
                        op2.binary_debug.write('      global_ply=%-1 mid=%-1 t=%-1 theta=%-1 sout=-1\n')
                    break
                (global_ply, mid, t, theta, sout_int) = s2.unpack(data[n:n+ntotal2])
                #print('  ', (global_ply, mid, t, theta, sout_int))
                try:
                    sout = sout_map[sout_int]
                except KeyError:
                    op2.log.error('cant parse global_ply=%s sout=%s; assuming 0=NO' % (
                        global_ply, sout_int))
                    sout = 'NO'

                global_ply_ids.append(global_ply)
                mids.append(mid)
                thicknesses.append(t)
                thetas.append(theta)
                souts.append(sout)
                if op2.is_debug_file:
                    op2.binary_debug.write('      global_ply=%s mid=%s t=%s theta=%s sout_int=%s sout=%r\n' % (
                        global_ply, mid, t, theta, sout_int, sout))
                n += ntotal2
                ilayer += 1
            n += ntotal2

            try:
                ft = ft_map[ft_int]
            except KeyError:
                op2.log.error('pid=%s cant parse ft=%s; should be HILL, HOFF, TSAI, STRN'
                               '...skipping' % (pid, ft_int))
                continue

            try:
                lam = lam_map[lam_int]
            except KeyError:
                op2.log.error('pid=%s cant parse lam=%s; should be HILL, HOFF, TSAI, STRN'
                               '...skipping' % (pid, lam_int))
                continue

            # apparently Nastran makes duplicate property ids...
            if pid in op2.pcomp.property_id:
                pcomp_pids_to_remove.add(pid)

            op2.add_pcompg(pid, global_ply_ids, mids, thicknesses, thetas=thetas, souts=souts,
                           nsm=nsm, sb=sb, ft=ft, tref=tref, ge=ge, lam=lam, z0=z0, comment='')
            nproperties += 1

        op2.pcompg.parse_cards()
        if pcomp_pids_to_remove:
            pcomp_pids_to_remove = np.array(list(pcomp_pids_to_remove), dtype='int32')
            op2.pcomp = op2.pcomp.remove_card_by_id(pcomp_pids_to_remove)

        if nproperties:
            op2.card_count['PCOMPG'] = nproperties
        return n

# PCOMPA

    def _read_pconeax(self, data: bytes, n: int) -> int:
        """
        (152,19,147) - Record 24
        """
        self.op2.log.info('geom skipping PCONEAX in EPT')
        return len(data)

    def _read_pconv(self, data: bytes, n: int) -> int:
        """common method for reading PCONVs"""
        op2 = self.op2
        #n = self._read_dual_card(data, n, self._read_pconv_nx, self._read_pconv_msc,
                                 #'PCONV', self._add_pconv)

        card_name = 'PCONV'
        card_obj = op2.pconv
        methods = {
            16 : self._read_pconv_nx_16,  # 16=4*4
            56 : self._read_pconv_msc_56, # 56=4*14
        }
        try:
            n, elements = op2.reader_geom2._read_double_card_load(
                card_name, card_obj,
                methods, data, n)
        except DoubleCardError:
            nx_method = partial(self._read_pconv_nx_16, card_obj)
            msc_method = partial(self._read_pconv_msc_56, card_obj)
            n, elements = op2._read_dual_card_load(
                data, n,
                nx_method, msc_method,
                card_name, self._add_op2_property)

        nelements = len(elements)
        for prop in elements:
            key = prop.pconid
            if key in op2.convection_properties:
                prop_old = op2.convection_properties[key]
                if prop != prop_old:
                    op2.log.warning(prop.raw_fields())
                    op2.log.warning(prop_old.raw_fields())
                    op2.log.warning(f'PCONV pconid={key}; old, new\n{prop_old}{prop}')
                    # this will fail due to a duplicate id
                    self._add_pconv(prop)
                #else:
                    # already exists
            else:
                self._add_pconv(prop)
        op2.card_count['PCONV'] = nelements

        return n

    def _read_pconv_nx_16(self, card_obj: PCONV, data: bytes, n: int) -> int:
        """
        (11001,110,411)- NX version
        """
        op2 = self.op2
        ntotal = 16  # 4*4
        nproperties = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0

        n, ints, floats = get_ints_floats(data, n, nproperties, 4,
                                          size=op2.size, endian=op2._endian)
        prop = op2.pconv
        pconv_id = ints[:, 0]
        material_id = ints[:, 1]
        form = ints[:, 2]
        expf = floats[:, 3]
        prop._save_nx(pconv_id, material_id, form, expf)
        props = []
        #struct_3if = Struct(op2._endian + b'3if')
        #for unused_i in range(nentries):
            #out = struct_3if.unpack(data[n:n+ntotal])
            #(pconid, mid, form, expf) = out
            #ftype = tid = chlen = gidin = ce = e1 = e2 = e3 = None
            #data_in = (pconid, mid, form, expf, ftype, tid, chlen,
                       #gidin, ce, e1, e2, e3)

            #prop = PCONV.add_op2_data(data_in)
            #props.append(prop)
            #n += ntotal
        return n, props

    def _read_pconv_msc_56(self, card_obj: PCONV, data: bytes, n: int) -> int:
        """
        (11001,110,411)- MSC version - Record 25
        """
        op2 = self.op2
        ntotal = 56  # 14*4
        s = Struct(op2._endian + b'3if 4i fii 3f')
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        props = []
        for unused_i in range(nentries):
            out = s.unpack(data[n:n+ntotal])
            (pconid, mid, form, expf, ftype, tid, unused_undef1, unused_undef2, chlen,
             gidin, ce, e1, e2, e3) = out
            data_in = (pconid, mid, form, expf, ftype, tid, chlen,
                       gidin, ce, e1, e2, e3)

            prop = PCONV.add_op2_data(data_in)
            props.append(prop)
            n += ntotal
        return n, props

    def _read_pconvm(self, data: bytes, n: int) -> int:
        """Record 24 -- PCONVM(2902,29,420)

        1 PID    I Property identification number
        2 MID    I Material identification number
        3 FORM   I Type of formula used for free convection
        4 FLAG   I Flag for mass flow convection
        5 COEF  RS Constant coefficient used for forced convection
        6 EXPR  RS Reynolds number convection exponent
        7 EXPPI RS Prandtl number convection exponent into the working fluid
        8 EXPPO RS Prandtl number convection exponent out of the working fluid
        """
        op2 = self.op2
        ntotal = 32  # 8*4
        nproperties = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(data, n, nproperties, 8,
                                          size=op2.size, endian=op2._endian)
        property_id = ints[:, 0]
        material_id = ints[:, 1]
        form = ints[:, 2]
        flag = ints[:, 3]
        coeff = floats[:, 4]
        expr = floats[:, 5]
        expri = floats[:, 6]
        exppo = floats[:, 7]
        i0 = (property_id == 0)
        pconvm = op2.pconvm
        pconvm._save(property_id, material_id, form, flag, coeff, expr, expri, exppo)
        if i0.sum():
            zsum = np.abs(ints[i0, :]).sum(axis=1)
            assert zsum.max() == 0, ints[i0, :]
            op2.log.warning('filtering PCONVM with property_ids=0')
            i = np.arange(nproperties, dtype='int32')[property_id != 0]
            pconvm.__apply_slice__(pconvm, i)
        assert pconvm.pconvm_id.min() > 0, pconvm.pconvm_id
        #_save(self, pconvm_id, material_id, flag, coefficient, expr, exppi, exppo)
        #structi = Struct(op2._endian + b'4i 4f')
        #for unused_i in range(nproperties):
            #out = structi.unpack(data[n:n+ntotal])
            #if out != (0, 0, 0, 0, 0., 0., 0., 0.):
                #(pconid, mid, form, flag, coeff, expr, expri, exppo) = out
                ##print(out)
                #prop = PCONVM(pconid, mid, coeff, form=form, flag=flag,
                              #expr=expr, exppi=expri, exppo=exppo, comment='')
                #op2._add_methods.add_convection_property_object(prop)
            #n += ntotal
        op2.card_count['PCONVM'] = nproperties
        return n

    def _read_pdamp(self, data: bytes, n: int) -> int:
        """
        PDAMP(202,2,45) - the marker for Record ???
        """
        op2 = self.op2
        ntotal = 8 * self.factor # 2*4
        nentries = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(data, n, nentries, 2,
                                          size=op2.size, endian=op2._endian)
        prop = op2.pdamp
        property_id = ints[:, 0]
        b = floats[:, 1]
        #prop.n = nentries
        prop._save(property_id, b)
        prop.write()

        #struct_if = Struct(mapfmt(op2._endian + b'if', self.size))
        #for unused_i in range(nentries):
            #out = struct_if.unpack(data[n:n+ntotal])
            ##(pid, b) = out
            #prop = PDAMP.add_op2_data(out)
            #self._add_op2_property(prop)
            #n += ntotal
        op2.card_count['PDAMP'] = nentries
        return n

    def _read_pdampt(self, data: bytes, n: int) -> int:  # 26
        self.op2.log.info('geom skipping PDAMPT in EPT')
        return len(data)

    def _read_pdamp5(self, data: bytes, n: int) -> int:  # 26
        self.op2.log.info('geom skipping PDAMP5 in EPT')
        return len(data)

# PDUM1
# PDUM2
# PDUM3
# PDUM4
# PDUM5
# PDUM6
# PDUM7
# PDUM8
# PDUM9

    def _read_pelas(self, data: bytes, n: int) -> int:
        """PELAS(302,3,46) - the marker for Record 39"""
        op2 = self.op2
        ntotal = 16 * self.factor # 4*4
        nproperties = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(data, n, nproperties, 4,
                                          size=op2.size, endian=op2._endian)
        prop = op2.pelas
        property_id = ints[:, 0]
        k = floats[:, 1]
        ge = floats[:, 2]
        s = floats[:, 3]
        prop._save(property_id, k, ge, s)
        prop.write()

        #struct_i3f = Struct(mapfmt(op2._endian + b'i3f', self.size))
        #for unused_i in range(nproperties):
            #edata = data[n:n+ntotal]
            #out = struct_i3f.unpack(edata)
            ##(pid, k, ge, s) = out
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PELAS=%s\n' % str(out))
            #prop = PELAS.add_op2_data(out)
            #self._add_op2_property(prop)
            #n += ntotal
        op2.card_count['PELAS'] = nproperties
        return n

    def _read_pfast_msc(self, data: bytes, n: int) -> int:
        r"""
        Word Name Type Description
        1 PID       I Property identification number
        2 MID       I Material property identification number
        3 D        RS Diameter of the fastener
        4 CONNBEH   I Connection behavior (0=FF/F, 1=FR, 10=RF/R, 11=RR)
        5 CONNTYPE  I Connection type (0=clamp, 1=hinge, 2=bolt)
        6 EXTCON    I External constraint flag (0=off, 1=on)
        7 CONDTYPE  I Condition type (0=rigid, 1=equivalent)
        8 WELDTYPE  I Weld type (0=spot weld, 1=but seam, 2=T-seam)

        9 MINLEN   RS Minimum length of spot weld
        10 MAXLEN  RS Maximum length of spot weld
        11 GMCHK    I Perform geometry check
        12 SPCGS    I SPC the master grid GS
        13 CMASS   RS Concentrated mass
        14 GE      RS Structureal Damping

        15 UNDEF(3) none Not used
        18 MCID    I Element stiffness coordinate system
        19 MFLAG   I Defined the coordinate system type
        20 KT(3)  RS Stiffness values in direction 1
        23 KR(3)  RS Rotation stiffness values in direction 1

        C:\MSC.Software\msc_nastran_runs\cfmass.op2
                   pid mid  D    con  con  ext  cond weld min max  chk  spc  cmass ge  und  und  und  mcid mfag kt1      kt2       kt3       kr1    kr2      kr3
        ints    = (99, 0,   0.1, 0,   0,   0,   0,   -1, 0.2, 5.0, 0,   0,   7.9, 0,   0,   0,   0,   -1, 0,   471200.0, 181200.0, 181200.0, 226.6, 45610.0, 45610.0)
        floats  = (99, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, -1, 0.2, 5.0, 0.0, 0.0, 7.9, 0.0, 0.0, 0.0, 0.0, -1, 0.0, 471200.0, 181200.0, 181200.0, 226.6, 45610.0, 45610.0)
        """
        op2 = self.op2
        #op2.show_data(data[n:], types='ifs')
        #ntotal = 92 * self.factor # 26*4
        #struct1 = Struct(op2._endian + b'ifii 3f')

        ntotal = 100 * self.factor # 25*4
        ndatai = len(data) - n
        nproperties = ndatai // ntotal
        delta = ndatai % ntotal
        assert delta == 0, 'len(data)-n=%s n=%s' % (ndatai, ndatai / 100.)

        n, ints, floats = get_ints_floats(data, n, nproperties, 25,
                                          size=op2.size, endian=op2._endian)
        property_id = ints[:, 0]
        assert property_id.min() > 0, property_id
        d = ints[:, 1]

        mcid = floats[:, 2] # ????
        #unused_connbeh  3 i
        #unused_conntype 4 i
        #unused_extcon   5 i
        #unused_condtype 6 i
        #unused_weldtype 7 i
        #unused_minlen 8 # f
        #unused_maxlen 9 # f
        #unused_gmcheck 10 # i
        #unused_spcgs 11   # i
        mass = floats[:, 12]
        ge = floats[:, 13]
        #aa 14
        #bb 15
        #cc 16
        mcid = ints[:, 17]
        # mflag 18
        kt = floats[:, [19, 20, 21]]
        kr = floats[:, [22, 23, 24]]
        op2.log.warning('skipping PFAST-MSC')

        #struct1 = Struct(op2._endian + b'2if 5i 2f2i2f 3i 2i 6f')
        #for unused_i in range(nproperties):
            #edata = data[n:n+ntotal]
            #out = struct1.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PFAST=%s\n' % str(out))
            #(pid, d, mcid, unused_connbeh, unused_conntype, unused_extcon,
             #unused_condtype, unused_weldtype, unused_minlen, unused_maxlen,
             #unused_gmcheck, unused_spcgs, mass, ge,
             #unused_aa, unused_bb, unused_cc, mcid, mflag,
             #kt1, kt2, kt3, kr1, kr2, kr3) = out

            #data_in = (pid, d, mcid, mflag, kt1, kt2, kt3,
                       #kr1, kr2, kr3, mass, ge)
            #prop = PFAST.add_op2_data(data_in)
            #str(prop)
            ##print(prop)
            #self._add_op2_property(prop)
            #n += ntotal
        op2.card_count['PFAST'] = nproperties
        return n

    def _read_pfast_nx(self, data: bytes, n: int) -> int:
        """
        PFAST(3601,36,55)
        NX only
        """
        op2 = self.op2
        ntotal = 48
        nproperties = (len(data) - n) // ntotal
        delta = (len(data) - n) % ntotal
        assert delta == 0, 'len(data)-n=%s n=%s' % (len(data) - n, (len(data) - n) / 48.)

        n, ints, floats = get_ints_floats(data, n, nproperties, 12,
                                          size=op2.size, endian=op2._endian)
        property_id = ints[:, 0]
        assert property_id.min() > 0, property_id

        #(pid, d, mcid, mflag, kt1, kt2, kt3, kr1, kr2, kr3, mass, ge) = out
        property_id = ints[:, 0]
        d = floats[:, 1]
        mcid = ints[:, 2]
        mflag = ints[:, 3]
        kt = floats[:, [4, 5, 6]]
        kr = floats[:, [7, 8, 9]]
        mass = floats[:, 10]
        ge = floats[:, 11]
        op2.log.warning('skipping PFAST-NX')
        #struct1 = Struct(op2._endian + b'ifii 8f')
        #for unused_i in range(nproperties):
            #edata = data[n:n+ntotal]
            #out = struct1.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PFAST=%s\n' % str(out))
            #(pid, d, mcid, mflag, kt1, kt2, kt3, kr1, kr2, kr3, mass, ge) = out

            #data_in = (pid, d, mcid, mflag, kt1, kt2, kt3,
                       #kr1, kr2, kr3, mass, ge)
            #prop = PFAST.add_op2_data(data_in)
            #self._add_op2_property(prop)
            #n += ntotal
        op2.card_count['PFAST'] = nproperties
        op2.to_nx(' because PFAST-NX was found')
        return n

    def _read_pelast(self, data: bytes, n: int) -> int:
        """
        Record 41 -- PELAST(1302,13,34)

        1 PID   I Property identification number
        2 TKID  I TABLEDi entry identification number for stiffness
        3 TGEID I TABLEDi entry identification number for structural
                  damping
        4 TKNID I TABLEDi entry
        """
        op2 = self.op2
        ntotal = 16 * self.factor
        nproperties = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(data, n, nproperties, 4,
                                          size=op2.size, endian=op2._endian)
        prop = op2.pelast
        property_id = ints[:, 0]
        table_k = ints[:, 1]
        table_ge = ints[:, 2]
        table_k_nonlinear = ints[:, 3]
        prop._save(property_id, table_k, table_ge, table_k_nonlinear)
        prop.write()

        #struct_4i = Struct(mapfmt(op2._endian + b'4i', self.size))
        #for unused_i in range(nproperties):
            #edata = data[n:n+ntotal]
            #out = struct_4i.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PELAST=%s\n' % str(out))
            ##(pid, tkid, tgeid, tknid) = out
            #prop = PELAST.add_op2_data(out)
            #op2._add_methods.add_pelast_object(prop)
            #n += ntotal
        op2.card_count['PELAST'] = nproperties
        return n

    def _read_pgap(self, data: bytes, n: int) -> int:
        """
        PGAP(2102,21,121) - the marker for Record 42
        """
        op2 = self.op2
        ntotal = 44 * self.factor
        struct_i10f = Struct(mapfmt(op2._endian + b'i10f', self.size))
        nproperties = (len(data) - n) // ntotal
        for unused_i in range(nproperties):
            edata = data[n:n+ntotal]
            out = struct_i10f.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  PGAP=%s\n' % str(out))
            (pid, u0, f0, ka, kb, kt, mu1, mu2, tmax, mar, trmin) = out
            op2.add_pgap(pid, u0=u0, f0=f0, ka=ka, kb=kb, mu1=mu1, kt=kt, mu2=mu2,
                         tmax=tmax, mar=mar, trmin=trmin, comment='')
            #prop = PGAP.add_op2_data(out)
            #self._add_op2_property(prop)
            n += ntotal
        op2.card_count['PGAP'] = nproperties
        return n

    def _read_phbdy(self, data: bytes, n: int) -> int:
        """
        PHBDY(2802,28,236) - the marker for Record 43
        """
        op2 = self.op2
        nproperties = (len(data) - n) // 16

        n, ints, floats = get_ints_floats(data, n, nproperties, 4,
                                          size=op2.size, endian=op2._endian)
        property_id = ints[:, 0]
        area_factor = floats[:, 1]
        diameter = floats[:, [2, 3]]
        phbdy = op2.phbdy
        phbdy._save(property_id, area_factor, diameter)

        #struct_i3f = Struct(op2._endian + b'ifff')
        #for unused_i in range(nproperties):
            #edata = data[n:n+16]
            #out = struct_i3f.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PHBDY=%s\n' % str(out))
            ##(pid, af, d1, d2) = out
            #prop = PHBDY.add_op2_data(out)
            #op2._add_methods.add_phbdy_object(prop)
            #n += 16
        op2.card_count['PHBDY'] = nproperties

        return n

    def _read_pintc(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping PINTC in EPT')
        return len(data)

    def _read_pints(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping PINTS in EPT')
        return len(data)

    def _read_pbeam3(self, data: bytes, n: int) -> int:
        op2 = self.op2
        if not hasattr(self, 'pbeam3'):
            op2.log.warning('skipping PBEAM3')
            return len(data)
        card_name = 'PBEAM3'
        card_obj = PBEAM3
        methods = {
            264 : self._read_pbeam3_264,
            456 : self._read_pbeam3_456,
        }
        try:
            n = op2.reader_geom2._read_double_card(
                card_name, card_obj, self._add_op2_property,
                methods, data, n)
        except DoubleCardError:
            raise
            op2.log.warning(f'try-except {card_name}')
        return n

    def _read_pbeam3_456(self, card_obj, data: bytes, n: int) -> int:
        r"""

        # per C:\MSC.Software\msc_nastran_runs\b3plod3.op2
        ints    = (2201, 1, 1.0, 0.1833, 0.0833, 0, -1.0, 0, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, -0.5,
                         2, 1.0, 0.1833, 0.0833, 0, -1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         2, 1.0, 0.1833, 0.0833, 0, -1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         1.0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   2901, 2, 0.1, 0.1, 0.1, 0, 0.2, 0, 0.5, 0, 0, 0.5, -0.5, 0, 0, -0.5,
                         2, 0.1, 0.1, 0.1, 0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         2, 0.1, 0.1, 0.1, 0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         1.0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        floats  = (2201, 1, 1.0, 0.1833, 0.0833, 0.0, -1.0, 0.0, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, -0.5,
                         2, 1.0, 0.1833, 0.0833, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         2, 1.0, 0.1833, 0.0833, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   2901, 2, 0.1, 0.1, 0.1, 0.0, 0.2, 0.0, 0.5, 0.0, 0.0, 0.5, -0.5, 0.0, 0.0, -0.5,
                         2, 0.1, 0.1, 0.1, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         2, 0.1, 0.1, 0.1, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        """
        op2 = self.op2
        #op2.show_data(data[n:])
        ntotal = 456 * self.factor # 114*4
        #
        struct1 = Struct(mapfmt(op2._endian +
                                b'2i' # pid, mid
                                b'3f' # A, Iy, Iz
                                b'5f'  # # a, b, c, d, e
                                b'5f fi  14f i' #fj ki  14f i
                                b'2i3f' #aa-ee - good
                                b'5f'   #ff-jj
                                b'5f'   #kk-oo
                                b'5f'   #pp-tt
                                b'6f'   #uu-zz
                                b'5f'   #aaa-eee
                                b'4i'   #fff-iii
                                # jjj-ooo
                                b'2f iii f'
                                # ppp-ttt
                                b'5f'
                                # uuu-zzz
                                b'6f'
                                b'30f', self.size))

        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert ndatai % ntotal == 0

        props = []
        for unused_i in range(nentries):
            #print(n, ntotal)
            datai = data[n:n+ntotal]
            #op2.show_data(datai, types='ifqd')
            n += ntotal

            (pid, mid, A, iz, iy,
             a, b, c, d, e,
             f, g, h, i, j,
             k, inta, l, m, ni, o, p, q, r, s, t, u, v, w, x, y, z,
             aa, bb, cc, dd, ee,
             ff, gg, hh, ii, jj,
             kk, ll, mm, nn, oo,
             pp, qq, rr, ss, tt,
             uu, vv, ww, xx, yy, zz,
             aaa, bbb, ccc, ddd, eee,
             fff, ggg, hhh, iii,
             jjj, kkk, lll, mmm, nnn, ooo,
             ppp, qqq, rrr, sss, ttt,
             uuu, vvv, www, xxx, yyy, zzz,
             *other) = struct1.unpack(datai)
            #print(pid, mid, A, iz, iy)
            #print('a-e', (a, b, c, d, e))
            #print('f-j', (f, g, h, i, j))
            #print(k, inta, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z)
            #print('aa-ee', (aa, bb, cc, dd, ee))
            #print('ff-jj', (ff, gg, hh, ii, jj))
            #print('kk-oo', (kk, ll, mm, nn, oo))
            #print('pp-tt', (pp, qq, rr, ss, tt))
            #print('uu-zz', (uu, vv, ww, xx, yy, zz))
            #print('aaa-eee', (aaa, bbb, ccc, ddd, eee))
            #print('fff-jjj', (fff, ggg, hhh, iii))
            #print('jjj-ooo', (jjj, kkk, lll, mmm, nnn, ooo))
            #print('ppp-ttt', (ppp, qqq, rrr, sss, ttt))
            #print('uuu-zzz', (uuu, vvv, www, xxx, yyy, zzz))

            if mid == 0:
                continue
            #assert sum(other) < 100, other
            prop = PBEAM3(
                pid, mid, A, iz, iy, iyz=None, j=None, nsm=0.,
                so=None,
                cy=None, cz=None,
                dy=None, dz=None,
                ey=None, ez=None,
                fy=None, fz=None,
                ky=1., kz=1.,
                ny=None, nz=None, my=None, mz=None,
                nsiy=None, nsiz=None, nsiyz=None,
                cw=None, stress='GRID',
                w=None, wy=None, wz=None, comment='')
            assert pid > 0, prop.get_stats()
            assert mid > 0, prop.get_stats()
            str(prop)
            props.append(prop)
            #self._add_op2_property(prop)
        #op2.card_count['PBEAM3'] = nentries
        return n, props

    def _read_pbeam3_264(self, card_obj, data: bytes, n: int) -> int:
        """
        TODO: partial
        # per test_cbeam_cbeam3???
        ints    = (2901, 2, 0.1, 0.1, 0.1, 0,   0.02, 0,   0.5, 0,   0,   0.5, -0.5, 0,   0,   -0.5, 2, 0.1, 0.1, 0.1, 0,   0.02, 0,   0,   0,   0,   0,   0,   0,   0,   0,   2, 0.1, 0.1, 0.1, 0,   0.02, 0,   0,   0,   0,   0,   0,   0,   0,   0,   1.0, 1.0,   0,   0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   -2,   0,   0)
        floats  = (2901, 2, 0.1, 0.1, 0.1, 0.0, 0.02, 0.0, 0.5, 0.0, 0.0, 0.5, -0.5, 0.0, 0.0, -0.5, 2, 0.1, 0.1, 0.1, 0.0, 0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2, 0.1, 0.1, 0.1, 0.0, 0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nan, 0.0, 0.0)
        """
        op2 = self.op2
        ntotal = 264 * self.factor # 66*4
        #                                       p/m ayz ae fj ki  14f i
        struct1 = Struct(mapfmt(op2._endian + b'2i 3f  5f 5f fi  14f i 30f 4i', self.size))

        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert ndatai % ntotal == 0

        props = []
        for unused_i in range(nentries):
            pid, mid, A, iz, iy, a, b, c, d, e, f, g, h, i, j, k, inta, *other = struct1.unpack(data[n:n+ntotal])
            #print(pid, mid, A, iz, iy)
            #print((a, b, c, d, e))
            #print((f, g, h, i, j))
            #print(k, inta)
            assert sum(other) < 100, other
            prop = PBEAM3(
                pid, mid, A, iz, iy, iyz=None, j=None, nsm=0.,
                so=None,
                cy=None, cz=None,
                dy=None, dz=None,
                ey=None, ez=None,
                fy=None, fz=None,
                ky=1., kz=1.,
                ny=None, nz=None, my=None, mz=None,
                nsiy=None, nsiz=None, nsiyz=None,
                cw=None, stress='GRID',
                w=None, wy=None, wz=None, comment='')
            assert pid > 0, prop.get_stats()
            assert mid > 0, prop.get_stats()
            str(prop)
            props.append(prop)
            n += ntotal
        return n, props

    def _read_pplane(self, data: bytes, n: int) -> int:
        """
        RECORD – PPLANE(3801,38,979)
        Word Name Type Description
        1 PID    I Property identification number
        2 MID    I Material identification number
        3 T     RS Default membrane thickness for Ti on the connection entry
        4 NSM   RS Nonstructural mass per unit area
        5 FOROPT I Formulation option number
        6 CSOPT  I Reserved for coordinate system definition of plane
        7 UNDEF(2) None

        ints    = (1, 1, 1.0, 0,   0,   0,   0,   0,   2, 2, 1.0, 0, 0, 0, 0, 0)
        floats  = (1, 1, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2, 2, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        """
        op2 = self.op2
        ntotal = 32 * self.factor # 8*4
        struct1 = Struct(mapfmt(op2._endian + b'2i 2f 4i', self.size))

        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert ndatai % ntotal == 0
        for unused_i in range(nentries):
            out = struct1.unpack(data[n:n+ntotal])
            pid, mid, t, nsm, foropt, csopt = out[:6]
            #print(out)
            assert csopt == 0, csopt
            pplane = op2.add_pplane(pid, mid, t=t, nsm=nsm,
                                     formulation_option=foropt)
            #pplane.validate()
            #print(pplane)
            #str(pplane)
            n += ntotal
        op2.card_count['PLPLANE'] = nentries
        return n

    def _read_plplane(self, data: bytes, n: int) -> int:
        """
        PLPLANE(4606,46,375)

        NX 10
        1 PID     I Property identification number
        2 MID     I Material identification number
        3 CID     I Coordinate system identification number
        4 STR CHAR4 Location of stress and strain output
        5 T      RS Default membrane thickness for Ti on the connection entry
        6 CSOPT  I  Reserved for coordinate system definition of plane
        7 UNDEF(5) None

        MSC 2016
        PID       I Property identification number
        2 MID     I Material identification number
        3 CID     I Coordinate system identification number
        4 STR CHAR4 Location of stress and strain output
        5 UNDEF(7 ) none Not used

        .. warning:: CSOPT ad T are not supported
        """
        op2 = self.op2
        ntotal = 44 * self.factor  # 4*11
        if self.size == 4:
            s = Struct(op2._endian + b'3i 4s f 6i')
        else:
            s = Struct(op2._endian + b'3q 8s d 6q')
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            out = s.unpack(data[n:n+ntotal])
            pid, mid, cid, location, unused_t, unused_csopt = out[:6]
            location = location.decode('latin1')
            #op2.show_data(data[n:n+ntotal], 'ifs')
            op2.add_plplane(pid, mid, cid=cid, stress_strain_output_location=location)
            n += ntotal
        op2.card_count['PLPLANE'] = nentries
        return n

    def _read_plsolid(self, data: bytes, n: int) -> int:
        """
        MSC 2016
        1 PID I Property identification number
        2 MID I Material identification number
        3 STR CHAR4 Location of stress and strain output
        4 UNDEF(4 ) none Not used

        NX 10
        1 PID I Property identification number
        2 MID I Material identification number
        3 STR CHAR4 Location of stress and strain output
        4 CSOPT I Reserved for coordinate system definition of plane
        5 UNDEF(3) None

        .. warning:: CSOPT is not supported
        """
        op2 = self.op2
        ntotal = 28 * self.factor  # 4*7
        if self.size == 4:
            struct1 = Struct(op2._endian + b'2i 4s 4i')
        else:
            struct1 = Struct(op2._endian + b'2q 8s 4q')
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            out = struct1.unpack(data[n:n+ntotal])
            pid, mid, location, unused_csopt, unused_null_a, unused_null_b, unused_null_c = out
            location = location.decode('latin1')
            #op2.show_data(data[n:n+ntotal], 'ifs')
            op2.add_plsolid(pid, mid, stress_strain=location, ge=0.)
            n += ntotal
        op2.card_count['PLSOLID'] = nentries
        return n

    def _read_pmass(self, data: bytes, n: int) -> int:
        """
        PMASS(402,4,44) - the marker for Record 48
        """
        op2 = self.op2
        ntotal = 8 * self.factor # 2*4
        nproperties = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(data, n, nproperties, 2,
                                          size=op2.size, endian=op2._endian)
        prop = op2.pmass
        property_id = ints[:, 0]
        mass = floats[:, 1]
        prop._save(property_id, mass)
        prop.write()

        #struct_if = Struct(mapfmt(op2._endian + b'if', self.size))
        #for unused_i in range(nproperties):
            #edata = data[n:n + ntotal]
            #out = struct_if.unpack(edata)
            ##out = (pid, mass)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PMASS=%s\n' % str(out))
            #prop = PMASS.add_op2_data(out)
            #self._add_op2_property_mass(prop)
            #n += ntotal
        return n

    def _read_prod(self, data: bytes, n: int) -> int:
        """
        PROD(902,9,29) - the marker for Record 49
        """
        op2 = self.op2
        ntotal = 24 * self.factor  # 6*4
        nproperties = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(data, n, nproperties, 6,
                                          size=op2.size, endian=op2._endian)
        prop = op2.prod
        property_id = ints[:, 0]
        material_id = ints[:, 1]
        A = floats[:, 2]
        J = floats[:, 3]
        c = floats[:, 4]
        nsm = floats[:, 5]
        prop._save(property_id, material_id, A, J, c, nsm)
        filter_large_property_ids(prop)
        prop.write()

        #struct_2i4f = Struct(mapfmt(op2._endian + b'2i4f', self.size))
        #for unused_i in range(nproperties):
            #edata = data[n:n+ntotal]
            #out = struct_2i4f.unpack(edata)
            ##(pid, mid, a, j, c, nsm) = out
            #prop = PROD.add_op2_data(out)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PROD=%s\n' % str(out))
            #self._add_op2_property(prop)
            #n += ntotal
        op2.card_count['PROD'] = nproperties
        return n

    def _read_pshear(self, data: bytes, n: int) -> int:
        """
        PSHEAR(1002,10,42) - the marker for Record 50
        """
        op2 = self.op2
        ntotal = 24 * self.factor
        nproperties = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(data, n, nproperties, 6,
                                          size=op2.size, endian=op2._endian)
        prop = op2.pshear
        property_id = ints[:, 0]
        material_id = ints[:, 1]
        t = floats[:, 2]
        nsm = floats[:, 3]
        f1 = floats[:, 4]
        f2 = floats[:, 5]
        prop._save(property_id, material_id, t, nsm, f1, f2)
        prop.write()

        #struct_2i4f = Struct(mapfmt(op2._endian + b'2i4f', self.size))
        #for unused_i in range(nproperties):
            #edata = data[n:n+ntotal]
            #out = struct_2i4f.unpack(edata)
            ##(pid, mid, t, nsm, f1, f2) = out
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PSHEAR=%s\n' % str(out))
            #prop = PSHEAR.add_op2_data(out)
            #self._add_op2_property(prop)
            #n += ntotal
        op2.card_count['PSHEAR'] = nproperties
        return n

    def _read_pshell(self, data: bytes, n: int) -> int:
        """
        PSHELL(2302,23,283) - the marker for Record 51
        """
        op2 = self.op2  # type: OP2Geom
        ntotal = 44 * self.factor  # 11*4
        nproperties = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(data, n, nproperties, 11,
                                          size=op2.size, endian=op2._endian)
        prop = op2.pshell

        property_id = ints[:, 0]

        #prop.property_id = property_id
        material_id = ints[:, [1, 3, 5, 10]]
        t = floats[:, 2]
        twelveIt3 = floats[:, 4]
        tst = floats[:, 6]
        nsm = floats[:, 7]
        z = floats[:, [8, 9]]
        #prop.n = nproperties
        prop._save(property_id, material_id, t, twelveIt3, tst, nsm, z)

        op2.pshell.parse_cards()

        pcomp_pids = op2.pcomp.property_id
        iall = np.arange(nproperties, dtype='int32')
        isave = [i for i, pid in enumerate(property_id)
                 if pid not in pcomp_pids]
        if isave:
            ifiltered = np.setdiff1d(iall, isave)
            if len(ifiltered):
                ibig = np.where(prop.property_id >= PID_CAP)
                removed_properties = prop.property_id[ibig]
                op2.log.warning(f'removing PCOMP property ids\n'
                                f'property_ids={removed_properties}\n')
                prop.slice_card_by_index(isave)

        filter_large_property_ids(prop)

        if prop.material_id.max() >= MID_CAP:
            ibig = np.where(prop.material_id >= PID_CAP)
            #i = np.where(prop.material_id < PID_CAP)[0]

            ibig_pid = np.unique(ibig[0])
            iall = np.arange(len(prop.property_id), dtype='int32')
            i = np.setdiff1d(iall, ibig_pid)  # subtract

            removed_materials = prop.material_id[ibig_pid, :]
            removed_properties = prop.property_id[ibig_pid]
            op2.log.warning(f'removing {prop.type} property ids with material_ids > {MID_CAP}\n'
                            f'removed_property_ids={removed_properties}\n'
                            f'removed_material_ids={removed_materials}')
            prop.__apply_slice__(prop, i)
        #print('prop.property_id', prop.property_id)
        #print('prop.material_id', prop.material_id)
        prop._filter_pcomp_pshells()
        prop.write()

        #s = Struct(mapfmt(op2._endian + b'iififi4fi', self.size))
        #for unused_i in range(nproperties):
            #edata = data[n:n+ntotal]
            #out = s.unpack(edata)
            #(pid, mid1, unused_t, mid2, unused_bk, mid3, unused_ts,
             #unused_nsm, unused_z1, unused_z2, mid4) = out
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PSHELL=%s\n' % str(out))
            #prop = PSHELL.add_op2_data(out)
            #n += ntotal

            #if pid in op2.properties:
                ## this is a fake PSHELL
                #propi = op2.properties[pid]
                #if prop == propi:
                    #op2.log.warning(f'Fake PSHELL {pid:d} (skipping):\n{propi}')
                    #nproperties -= 1
                    #continue
                ##assert propi.type in ['PCOMP', 'PCOMPG'], propi.get_stats()
                #op2.log.error(f'PSHELL {pid:d} is also {propi.type} (skipping PSHELL):\n{propi}{prop}')
                #nproperties -= 1
                #continue
            ##continue
            ##if max(pid, mid1, mid2, mid3, mid4) > 1e8:
                ##self.big_properties[pid] = prop
            ##else:
            #self._add_op2_property(prop)
        nproperties = len(prop.property_id)
        if nproperties:
            op2.card_count['PSHELL'] = nproperties
        return n

    def _read_psolid(self, data: bytes, n: int) -> int:
        """
        PSOLID(2402,24,281) - the marker for Record 52
        """
        op2 = self.op2
        #print("reading PSOLID")
        ntotal = 28 * self.factor # 7*4

        nproperties = (len(data) - n) // ntotal
        prop: PSOLID = op2.psolid
        n, ints, floats, strings = get_ints_floats_strings(
            data, n, nproperties, 7, size=op2.size, endian=op2._uendian)

        #  adding an H...make it SMECH instead of SMEC
        fctn = strings[:, 6].astype('|U8')
        i_smech = np.where(fctn == 'SMEC')[0]
        fctn[i_smech] = 'SMECH'

        ##data_in = [pid, mid, cid, inp, stress, isop, fctn]
        property_id = ints[:, 0]
        material_id = ints[:, 1]
        coord_id = ints[:, 2]
        integ = ints[:, 3]
        stress = ints[:, 4]
        isop = ints[:, 5]
        prop._save(property_id, material_id, coord_id,
                   integ, stress, isop, fctn)

        ifake = (fctn == 'FAKE')
        if ifake.sum():
            ibool = ~ifake
            i = np.arange(nproperties, dtype='int32')[ibool]
            removed_property_ids = prop.property_id[ifake]
            op2.log.warning(f'filtering {prop.type} property_id={removed_property_ids} because fctn="FAKE"')
            prop.__apply_slice__(prop, i)

        filter_large_property_ids(prop)
        #filter_large_material_ids(prop)
        prop.write()


        #if self.size == 4:
            #struct_6i4s = Struct(op2._endian + b'6i4s')
        #else:
            #struct_6i4s = Struct(op2._endian + b'6q8s')

        #nproperties_found = 0
        #for unused_i in range(nproperties):
            #edata = data[n:n+ntotal]
            #out = struct_6i4s.unpack(edata)
            ##(pid, mid, cid, inp, stress, isop, fctn) = out
            ##data_in = [pid, mid, cid, inp, stress, isop, fctn]
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PSOLID=%s\n' % str(out))

            #n += ntotal
            #fctn = out[-1]
            #if fctn == b'FAKE':
                #op2.log.warning('    PSOLID=%s; is this a PCOMPLS?' % str(out))
                #continue
            #prop = PSOLID.add_op2_data(out)
            #self._add_op2_property(prop)
            #nproperties_found += 1
        op2.card_count['PSOLID'] = len(prop.property_id)
        return n

# PSOLIDL
# PTRIA6
# PTSHELL

    def _read_ptube(self, data: bytes, n: int) -> int:
        """
        PTUBE(1602,16,30) - the marker for Record 56

        .. todo:: OD2 only exists for heat transfer...
                  how do i know if there's heat transfer at this point?
                  I could store all the tubes and add them later,
                  but what about themal/non-thermal subcases?

        .. warning:: assuming OD2 is not written (only done for thermal)
        """
        op2 = self.op2
        ntotal = 20 * self.factor # 5*4
        nproperties = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(data, n, nproperties, 5,
                                          size=op2.size, endian=op2._endian)
        prop = op2.ptube
        property_id = ints[:, 0]
        material_id = ints[:, 1]
        diameter = np.full((nproperties, 2), np.nan, dtype=floats.dtype)
        diameter[:, 0] = floats[:, 2]
        t = floats[:, 3]
        nsm = floats[:, 4]
        prop._save(property_id, material_id, diameter, t, nsm)
        prop.write()

        #struct_2i3f = Struct(op2._endian + b'2i3f')
        #for unused_i in range(nproperties):
            #edata = data[n:n+20]  # or 24???
            #out = struct_2i3f.unpack(edata)
            #(pid, mid, OD, t, nsm) = out
            #data_in = [pid, mid, OD, t, nsm]
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PTUBE=%s\n' % str(out))
            #prop = PTUBE.add_op2_data(data_in)
            #self._add_op2_property(prop)
            #n += 20
        op2.card_count['PTUBE'] = nproperties
        return n

    def _read_pset(self, data: bytes, n: int) -> int:
        op2 = self.op2
        if not hasattr(self.op2, 'pval'):
            op2.log.warning('skipping PVAL')
            return len(data)

        struct_5i4si = Struct(op2._endian + b'5i4si')
        nentries = 0
        while  n < len(data):
            edata = data[n:n+28]
            out = struct_5i4si.unpack(edata)
            #print(out)
            idi, poly1, poly2, poly3, cid, typei, typeid = out
            typei = typei.rstrip().decode('latin1')
            assert typei in ['SET', 'ELID'], (idi, poly1, poly2, poly3, cid, typei, typeid)
            if op2.is_debug_file:
                op2.binary_debug.write('  PVAL=%s\n' % str(out))
            #print(idi, poly1, poly2, poly3, cid, typei, typeid)
            typeids = []
            n += 28
            while typeid != -1:
                typeids.append(typeid)
                typeid, = op2.struct_i.unpack(data[n:n+4])
                n += 4
                #print(val)
            #print(typeids)
            # PSET ID POLY1 POLY2 POLY3 CID SETTYP ID
            if len(typeids) == 1:
                typeids = typeids[0]
            op2.add_pset(idi, poly1, poly2, poly3, cid, typei, typeids)
        op2.card_count['PSET'] = nentries
        return n

    def _read_pval(self, data: bytes, n: int) -> int:
        """
        PVAL(10201,102,400)

        Word Name Type Description
        1 ID       I p-value set identification number
        2 POLY1    I Polynomial order in 1 direction of the CID system
        3 POLY2    I Polynomial order in 2 direction of the CID system
        4 POLY3    I Polynomial order in 2 direction of the CID system
        5 CID      I Coordinate system identification number
        6 TYPE CHAR4 Type of set provided: "SET" or "ELID"
        7 TYPEID   I SET identification number or element identification
                     number with this p-value specification.
        Words 1 through 7 repeat until End of Record
        """
        op2 = self.op2
        if not hasattr(self.op2, 'pval'):
            op2.log.warning('skipping PVAL')
            return len(data)
        #op2.show_data(data[n:])
        if self.size == 4:
            struct_5i4si = Struct(op2._endian + b'5i 4s i')
            struct_i = op2.struct_i
        else:
            struct_5i4si = Struct(op2._endian + b'5q 8s q')
            struct_i = op2.struct_q

        nentries = 0
        ntotal = 28 * self.factor
        size = self.size
        while  n < len(data):
            edata = data[n:n+ntotal]
            out = struct_5i4si.unpack(edata)
            #print(out)
            idi, poly1, poly2, poly3, cid, typei, typeid = out
            typei = typei.rstrip().decode('latin1')
            assert typei in ['SET', 'ELID'], f'idi={idi} poly1={poly1} poly2={poly2} poly3={poly3} cid={cid} typei={typei} typeid={typeid}'
            if op2.is_debug_file:
                op2.binary_debug.write('  PVAL=%s\n' % str(out))
            #print(idi, poly1, poly2, poly3, cid, typei, typeid)
            typeids = []
            n += ntotal
            while typeid != -1:
                typeids.append(typeid)
                typeid, = struct_i.unpack(data[n:n+size])
                n += size
                #print(val)
            #print(typeids)
            # PVAL ID POLY1 POLY2 POLY3 CID SETTYP ID
            op2.add_pval(idi, poly1, poly2, poly3, cid, typei, typeids)
        op2.card_count['PVAL'] = nentries
        return n

    def _read_pvisc(self, data: bytes, n: int) -> int:
        """PVISC(1802,18,31) - the marker for Record 39"""
        op2 = self.op2
        nproperties = (len(data) - n) // 12

        n, ints, floats = get_ints_floats(data, n, nproperties, 3,
                                          size=op2.size, endian=op2._endian)
        prop = op2.pvisc
        property_id = ints[:, 0]
        ce = floats[:, 1]
        cr = floats[:, 2]
        prop._save(property_id, ce, cr)
        prop.write()

        #struct_i2f = Struct(op2._endian + b'i2f')
        #for unused_i in range(nproperties):
            #edata = data[n:n+12]
            #out = struct_i2f.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PVISC=%s\n' % str(out))
            ##(pid, ce, cr) = out
            #prop = PVISC.add_op2_data(out)
            #self._add_op2_property(prop)
            #n += 12
        op2.card_count['PVISC'] = nproperties
        return n

# PWELD
# PWSEAM
    def _read_view(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping VIEW in EPT')
        return len(data)

    def _read_view3d(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping VIEW3D in EPT')
        return len(data)

def filter_large_property_ids(prop):
    if len(prop.property_id) == 0 or prop.property_id.max() < PID_CAP:
        return
    ibig = np.where(prop.property_id >= PID_CAP)[0]
    i = np.where(prop.property_id < PID_CAP)[0]
    removed_properties = prop.property_id[ibig]
    prop.model.log.warning(f'{prop.type}: removing property ids > {PID_CAP}\n'
                           f'property_ids={removed_properties}\n')
    #prop.slice_card_by_index(i)
    prop.__apply_slice__(prop, i)

def filter_large_material_ids(prop):
    if len(prop.material_id) == 0 or prop.material_id.max() < MID_CAP:
        return
    ibig = np.where(prop.material_id >= MID_CAP)[0]
    i = np.where(prop.material_id < MID_CAP)[0]
    removed_properties = prop.property_id[ibig]
    removed_materials = prop.material_id[ibig]
    prop.model.log.warning(f'{prop.type}: removing property ids with material_id > {MID_CAP}\n'
                           f'property_ids={removed_properties}\n'
                           f'material_ids={removed_materials}\n'
                           )
    #prop.slice_card_by_index(i)
    prop.__apply_slice__(prop, i)

def break_by_minus1(idata):
    """helper for ``read_nsm_nx``"""
    i1 = 0
    i = 0
    i2 = None
    packs = []
    for idatai in idata:
        #print('data[i:] = ', data[i:])
        if idatai == -1:
            i2 = i
            packs.append((i1, i2))
            i1 = i2 + 1
            i += 1
            continue
        i += 1
    #print(packs)
    return packs
