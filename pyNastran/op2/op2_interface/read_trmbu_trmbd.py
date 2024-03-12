from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.op2.result_objects.op2_results import TRMBD, TRMBU
if TYPE_CHECKING:
    from pyNastran.op2.op2_interface.op2_reader import OP2Reader
    from pyNastran.op2.op2 import OP2


#def trmbx_tsteps_to_tstep_index_mapper(tsteps_set: set[float]) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    ## Create time step to tid mapper per subcase
    #tsteps = np.array(list(tsteps_set))

    ## Sort because set does not retain ordering
    #tstep_indices = np.argsort(tsteps)
    #time_steps = tsteps[tstep_indices]

    #ntimes = time_steps.shape[0]
    #tstep_to_index_mapper = {
        #time_steps[itime]: itime for itime in range(ntimes)}


    ## Create time step to tid mapper per subcase
    ##tsteps = np.array(list(time_steps))
    ##tstep_indices = np.argsort(tsteps)  # Sort because set does not retain ordering
    ##time_steps = tsteps[tstep_indices]

    ##ntimes = time_steps.shape[0]
    ##tstep_to_index_mapper = {
        ##time_steps[itime]: itime for itime in range(ntimes)
    ##}
    #return time_steps, tstep_to_index_mapper

def read_trmbu(self: OP2Reader) -> None:
    """
    Reads the TRMBU table

    Table of Euler Angles for transformation from material to basic coordinate system
    in the undeformed configuration

    Record - HEADER
    +------+---------+-------+---------------------------------+
    | Word | Name    | Type  | Description                     |
    +======+=========+=======+=================================+
    |  1   | NAME(2) | CHAR4 | Data block name                 |
    +------+---------+-------+---------------------------------+
    |  3   | WORD    | I     | No Def or Month, Year, One, One |
    +------+---------+-------+---------------------------------+

    Record IDENT
    +---------+------------+-------+---------------------------------------------------+
    | Word    | Name       | Type  | Description                                       |
    +=========+============+=======+===================================================+
    |  1      | ACODE      | I     | Device code + 10*Approach code = 60 + iand        |
    |         |            |       | (print, plot)                                     |
    +---------+------------+-------+---------------------------------------------------+
    |  2      | TCODE      | I     | Table code                                        |
    +---------+------------+-------+---------------------------------------------------+
    |  3      | ELTYPE(C)  | I     | Element type number                               |
    +---------+------------+-------+---------------------------------------------------+
    |  4      | SUBCASE    | I     | Subcase number                                    |
    +---------+------------+-------+---------------------------------------------------+
    | 5       | TIME       | RS    | Time Step                                         |
    +---------+------------+-------+---------------------------------------------------+
    | 6-9     | UNDEF(4)   | None  |                                                   |
    +---------+------------+-------+---------------------------------------------------+
    | 10      | NUMWDE     | I     | Number of words per entry in DATA record          |
    +---------+------------+-------+---------------------------------------------------+
    | 11-50   | UNDEF(40)  | None  |                                                   |
    +---------+------------+-------+---------------------------------------------------+
    | 51-82   | TITLE(32)  | CHAR4 | Title character string (TITLE)                    |
    +---------+------------+-------+---------------------------------------------------+
    | 83-114  | SUBTITL(32)| CHAR4 | Subtitle character string (SUBTITLE)              |
    +---------+------------+-------+---------------------------------------------------+
    | 115-146 | LABEL(32)  | CHAR4 | LABEL character string (LABEL)                    |
    +---------+------------+-------+---------------------------------------------------+

    Record DATA
    +------+------------+-------+---------------------------------------------------+
    | Word | Name       | Type  | Description                                       |
    +======+============+=======+===================================================+
    |  1   | EID        | I     | Element ID * 10 + device code                     |
    +------+------------+-------+---------------------------------------------------+
    |  3   | AX         | RS    | Euler angle X                                     |
    +------+------------+-------+---------------------------------------------------+
    |  3   | AY         | RS    | Euler angle Y                                     |
    +------+------------+-------+---------------------------------------------------+
    |  3   | AZ         | RS    | Euler angle Z                                     |
    +------+------------+-------+---------------------------------------------------+

    Record TRAILER
    +------+------------+-------+---------------------------------------------------+
    | Word | Name       | Type  | Description                                       |
    +======+============+=======+===================================================+
    |  1   | UNDEF(6)   | None  |                                                   |
    +------+------------+-------+---------------------------------------------------+
    """
    #trmbu
    op2: OP2 = self.op2
    #is_geometry = op2.is_geometry
    unused_table_name = self._read_table_name(rewind=False)
    self.read_markers([-1])
    data = self._read_record()
    # print(self.show_data(data, types='ifsqd'))
    # 101, 466286, 15,      1, 1, 180, 0
    # 101, 466286, ncoords, 1, 1, 180, 0
    factor = self.factor
    #if self.size == 4:
        #idtype = 'int32'
        #fdtype = 'float32'
    #else:
        #idtype = 'int64'
        #fdtype = 'float64'
    assert len(data) == 28 * factor, len(data)

    self.read_3_markers([-2, 1, 0])
    data = self._read_record() # CSTM
    # print(self.show_data(data, types='s'))
    # assert len(data) == 8 * factor, len(data)

    self.read_3_markers([-3, 1, 0])
    itable = -4

    read_record = self._read_record if self.read_mode == 2 else self._skip_record
    blocks = []
    while 1:
        markers = self.get_nmarkers(1, rewind=True)
        if markers == [0]:
            break
        data = read_record()
        blocks.append(data)
        self.read_markers([itable, 1, 0])
        itable -= 1
    markers = self.get_nmarkers(1, rewind=False)

    if self.read_mode == 1:
        return
    nblocks = len(blocks)

    assert(nblocks % 2 == 0)  # Should be even, first block is IDENT, second block is DATA

    # Get time steps per subcase

    # Set because time step can be repeated
    #time_steps_set = set([])  # {t1, t2, ...}
    #subcases_set = set([])

    # Read data
    element_type_to_str_map = {
        300: 'CHEXA',
        301: 'CPENTA',
        302: 'CTETRA',
        303: 'CPYRAM',

        312: 'CTRAX3',
        314: 'CTRAX6',

        313: 'CQUADX4',
        315: 'CQUADX8',

        320: 'CPLSTS3',
        322: 'CPLSTS6',

        343: 'CTRIA6',
        344: 'CQUAD8',

        345: 'CTRIAR',
        346: 'CQUADR',

        347: 'CBAR',
        348: 'CBEAM',

        363: 'CROD',
    }
    trmbus = op2.op2_results.trmbu

    # Set because time step can be repeated
    #time_steps_set = set([])  # {t1, t2, ...}
    subcases_set = set([])

    size = op2.size
    ntimes = nblocks // 2
    for i in range(0, nblocks, 2):
        itime = i // 2
        block0 = blocks[i]
        identifiers_int = np.frombuffer(block0, dtype=op2.idtype8)
        identifiers_float = np.frombuffer(block0, dtype=op2.fdtype8)

        #acode = identifiers_int[0]
        #tcode = identifiers_int[1]
        #eltype = identifiers_int[2]
        #isubcase = identifiers_int[3]
        #time_step = identifiers_float[4]

        #subcases_set.add(isubcase)
        #time_steps_set.add(time_step)

        identifiers_int = np.frombuffer(block0, dtype=op2.idtype8)
        identifiers_float = np.frombuffer(block0, dtype=op2.fdtype8)

        approach_code = identifiers_int[0]
        tcode = identifiers_int[1]
        element_type = identifiers_int[2]
        isubcase = identifiers_int[3]
        time_step = identifiers_float[4]
        subcases_set.add(isubcase)

        op2.data_code = {}
        op2.subtable_name = ''
        tCode = tcode
        int3 = element_type
        op2.element_type = element_type
        op2._set_approach_code(approach_code, tCode, int3, isubcase)
        #------------------------------------------------------------------------

        if isubcase not in trmbus:
            op2._read_title_helper(block0)
            trmbus[isubcase] = TRMBU(ntimes, **op2.data_code)
            time = np.zeros(ntimes, dtype=op2.fdtype8)
            trmbus[isubcase].time = time
        trmbu = trmbus[isubcase]
        time = trmbu.time
        time[itime] = time_step

        eid_euler_data = blocks[i+1]
        int_data = np.frombuffer(eid_euler_data, dtype=op2.idtype8)
        float_data = np.frombuffer(eid_euler_data, dtype=op2.fdtype8)
        numwide = 4
        if element_type in element_type_to_str_map:
            element = element_type_to_str_map[element_type]

        elif element_type == 316:
            element = 'CPLSTN3'
        #elif element_type == 316:
            #numwide = 4
            #element = 'CPLSTS3'
        #elif element_type == 318:
            #nnodes = 6
            #element = 'CPLSTS6'
        elif element_type == 317:
            element = 'CPLSTN4'
        elif element_type == 318:
            element = 'CPLSTN6'
        elif element_type == 319:
            element = 'CPLSTN8'

        elif element_type == 321:
            element = 'CPLSTS4'
        elif element_type == 323:
            element = 'CPLSTS8'
        else:
            print(int_data[:10])
            print(op2.code_information())
            raise NotImplementedError(f"Element type {element_type} is not implemented.\n")

        assert numwide == 4, numwide
        nelements = int_data.shape[0] // numwide  # element_id + 3 euler angles
        assert int_data.shape[0] % numwide == 0, 'failed eid+3 euler'
        int_data = int_data.reshape(nelements, numwide)
        float_data = float_data.reshape(nelements, numwide)

        eids = (int_data[:, 0] - op2.device_code) // 10
        eulers = float_data[:, 1:]

        if element not in trmbu.eulers:
            trmbu.eulers[element] = np.empty([ntimes, nelements, 4])

        #itime = tstep_to_index_mapper[time_step]
        trmbu.eulers[element][itime, :, 0] = eids
        trmbu.eulers[element][itime, :, 1:] = eulers
    assert len(subcases_set) == 1, subcases_set
    asdf
    return

def read_trmbd(self: OP2Reader) -> None:
    """
    Reads the TRMBD table

    Table of Euler Angles for transformation from material to basic coordinate system in the deformed configuration

    Record - HEADER
    +------+---------+-------+---------------------------------+
    | Word | Name    | Type  | Description                     |
    +======+=========+=======+=================================+
    |  1   | NAME(2) | CHAR4 | Data block name                 |
    +------+---------+-------+---------------------------------+
    |  3   | WORD    | I     | No Def or Month, Year, One, One |
    +------+---------+-------+---------------------------------+

    Record IDENT
    +------+------------+-------+---------------------------------------------------+
    | Word | Name       | Type  | Description                                       |
    +======+============+=======+===================================================+
    |  1   | ACODE(C)   | I     | Device code + 10*Approach code = 60 + iand        |
    |      |            |       | (print, plot)                                     |
    +------+------------+-------+---------------------------------------------------+
    |  2   | TCODE      | I     | Table code                                        |
    +------+------------+-------+---------------------------------------------------+
    |  3   | ELTYPE(C)  | I     | Element type number                               |
    +------+------------+-------+---------------------------------------------------+
    |  4   | SUBCASE    | I     | Subcase number                                    |
    +------+------------+-------+---------------------------------------------------+
    | TCODE,1 = 1       | Sort 1                                                    |
    +------+------------+-------+---------------------------------------------------+
    | ACODE,4 = 02      | Real eigenvalues                                          |
    +------+------------+-------+---------------------------------------------------+
    | 5    | UNDEF(5)   | None  |                                                   |
    +------+------------+-------+---------------------------------------------------+
    | ACODE,4 = 06      | Transient                                                 |
    +------+------------+-------+---------------------------------------------------+
    | 5    | TIME       | RS    | Current time step                                 |
    +------+------------+-------+---------------------------------------------------+
    | 6    | UNDEF(4)   | None  |                                                   |
    +------+------------+-------+---------------------------------------------------+
    | ACODE,4 = 10      | Nonlinear statics                                         |
    +------+------------+-------+---------------------------------------------------+
    | 5    | LFTSFQ     | RS    | Load step for SOL401 arc-length only              |
    +------+------------+-------+---------------------------------------------------+
    | 6    | TIME       | RS    | Time for SOL401 arc-length only                   |
    +------+------------+-------+---------------------------------------------------+
    | 7    | AL_TOTAL   | RS    | Accumulated arc-length for SOL401 arc-length only |
    +------+------------+-------+---------------------------------------------------+
    | 8    | UNDEF(2)   | None  |                                                   |
    +------+------------+-------+---------------------------------------------------+
    | End ACODE,4       |       |                                                   |
    +------+------------+-------+---------------------------------------------------+
    | TCODE,1 = 02      | Sort 2                                                    |
    +------+------------+-------+---------------------------------------------------+
    | End TCODE,1       |       |                                                   |
    +------+------------+-------+---------------------------------------------------+
    | 10   | NUMWDE     | I     | Number of words per entry in DATA record          |
    +------+------------+-------+---------------------------------------------------+
    | 11   | UNDEF(40)  | None  |                                                   |
    +------+------------+-------+---------------------------------------------------+
    | 51   | TITLE(32)  | CHAR4 | Title character string (TITLE)                    |
    +------+------------+-------+---------------------------------------------------+
    | 83   | SUBTITL(32)| CHAR4 | Subtitle character string (SUBTITLE)              |
    +------+------------+-------+---------------------------------------------------+
    | 115  | LABEL(32)  | CHAR4 | LABEL character string (LABEL)                    |
    +------+------------+-------+---------------------------------------------------+

    Record DATA
    +------+------------+-------+---------------------------------------------------+
    | Word | Name       | Type  | Description                                       |
    +======+============+=======+===================================================+
    |  1   | ELID       | I     | Element ID * 10 + device code                     |
    +------+------------+-------+---------------------------------------------------+
    |  2   | GRID       | I     | External grid ID                                  |
    +------+------------+-------+---------------------------------------------------+
    |  3   | AX         | RS    | Euler angle X                                     |
    +------+------------+-------+---------------------------------------------------+
    |  3   | AY         | RS    | Euler angle Y                                     |
    +------+------------+-------+---------------------------------------------------+
    |  3   | AZ         | RS    | Euler angle Z                                     |
    +------+------------+-------+---------------------------------------------------+
    | Words 2 thru 5 repeat for each end of corner grid point                       |
    | (for example, 8 grids for CHEXA, 6 grids for XPENTA, 2 grids for CBAR,        |
    | 4 grids for CQUADX4                                                           |
    +------+------------+-------+---------------------------------------------------+

    Record TRAILER
    +------+------------+-------+---------------------------------------------------+
    | Word | Name       | Type  | Description                                       |
    +======+============+=======+===================================================+
    |  1   | UNDEF(6)   | None  |                                                   |
    +------+------------+-------+---------------------------------------------------+
    """
    op2: OP2 = self.op2
    #is_geometry = op2.is_geometry
    unused_table_name = self._read_table_name(rewind=False)
    self.read_markers([-1])
    data = self._read_record()
    #print(self.show_data(data, types='ifsqd'))
    # 101, 466286, 15,      1, 1, 180, 0
    # 101, 466286, ncoords, 1, 1, 180, 0
    header_int = np.frombuffer(data, dtype=op2.idtype8)
    #header_float = np.frombuffer(data, dtype=op2.fdtype8)
    ncoords = header_int[2]
    #print('data =', header_int)
    #print('ncoords =', ncoords)
    factor = self.factor
    #if self.size == 4:
        #idtype = 'int32'
        #fdtype = 'float32'
    #else:
        #idtype = 'int64'
        #fdtype = 'float64'
    assert len(data) == 28 * factor, len(data)

    self.read_3_markers([-2, 1, 0])
    data = self._read_record() # CSTM
    #assert len(data) == 8 * factor, len(data)

    self.read_3_markers([-3, 1, 0])

    itable = -4

    blocks = []
    read_record = self._read_record if self.read_mode == 2 else self._skip_record
    while 1:
        markers = self.get_nmarkers(1, rewind=True)
        if markers == [0]:
            break
        data = read_record()
        blocks.append(data)
        self.read_markers([itable, 1, 0])
        itable -= 1
    markers = self.get_nmarkers(1, rewind=False)

    if self.read_mode == 1:
        return
    nblocks = len(blocks)

    # Should be even, first block is IDENT, second block is DATA
    assert(nblocks % 2 == 0)


    element_type_to_str_map = {
        300: ('CHEXA', 8),
        301: ('CPENTA', 6),
        302: ('CTETRA', 4),
        303: ('CPYRAM', 5),

        312: ('CTRAX3', 3),
        314: ('CTRAX6', 3),

        313: ('CQUADX4', 4),
        315: ('CQUADX8', 4),

        320: ('CPLSTS3', 3),
        322: ('PLSTS6', 3),

        343: ('CTRIA6', 3),
        344: ('CQUAD8', 4),

        345: ('CTRIAR', 3),
        346: ('CQUADR', 4),

        347: ('CBAR', 2),
        348: ('CBEAM', 2), # ???

        363: ('CROD', 2),
    }
    # Set because time step can be repeated
    time_steps_set = set([])  # {t1, t2, ...}
    subcases_set = set([])

    # Get time steps per subcase
    size = op2.size
    ntimes = nblocks // 2
    trmbds = op2.op2_results.trmbd

    for i in range(0, nblocks, 2):
        itime = i // 2
        block0 = blocks[i]
        identifiers_int = np.frombuffer(block0, dtype=op2.idtype8)
        identifiers_float = np.frombuffer(block0, dtype=op2.fdtype8)

        approach_code = identifiers_int[0]
        tCode = identifiers_int[1]
        int3 = identifiers_int[2]

        element_type = identifiers_int[2]
        isubcase = identifiers_int[3]
        time_step = identifiers_float[4]

        subcases_set.add(isubcase)
        time_steps_set.add(time_step)


        #------------------------------------------------------------------------
        # Read data
    #for i in range(0, nblocks, 2):
        #block0 = blocks[i]
        #identifiers_int = np.frombuffer(block0, dtype=op2.idtype8)
        #identifiers_float = np.frombuffer(block0, dtype=op2.fdtype8)

        #acode = identifiers_int[0]
        #tcode = identifiers_int[1]
        #element_type = identifiers_int[2]
        #isubcase = identifiers_int[3]
        #time_step = identifiers_float[4]

        op2.data_code = {}
        op2.subtable_name = ''
        int3 = element_type
        op2.element_type = element_type
        op2._set_approach_code(approach_code, tCode, int3, isubcase)

        if isubcase not in trmbds:
            op2._read_title_helper(block0)
            trmbds[isubcase] = TRMBD(**op2.data_code)
            trmbds[isubcase].time = np.zeros(ntimes, dtype=op2.fdtype8)
        trmbd = trmbds[isubcase]
        time = trmbd.time
        time[itime] = time_step
        #------------------------------------------------------------------------
        block1 = blocks[i+1]
        int_data = np.frombuffer(block1, dtype=op2.idtype8)
        float_data = np.frombuffer(block1, dtype=op2.fdtype8)

        #grids = int_data[:, np.array([1, 5, 9, 13, 17, 21])]
        if element_type in element_type_to_str_map:
            element, nnodes = element_type_to_str_map[element_type]

        elif element_type == 316:
            nnodes = 3
            element = 'CPLSTN3'
        elif element_type == 317:
            nnodes = 4
            element = 'CPLSTN4'
        elif element_type == 318:
            nnodes = 3
            element = 'CPLSTN6'
        elif element_type == 319:
            nnodes = 4
            element = 'CPLSTN8'

        #elif element_type == 323:
            #nnodes = 6
            #element = 'CPLSTS6'
        elif element_type == 321:
            nnodes = 4
            element = 'CPLSTS4'
        elif element_type == 323:
            nnodes = 4
            element = 'CPLSTS8'
        else:
            #print(int_data[:20])
            print(op2.code_information())
            raise NotImplementedError(f"Element type {element_type:d} is not implemented.\n")
        nelements, int_data, float_data = reshape_trmbd(
            element, nnodes, int_data, float_data)

        eids = (int_data[:, 0] - op2.device_code) // 10
        #nelements = len(eids)
        #print('ndatai =', int_data.size)
        #print('nelements =', nelements)
        index = np.arange(1, 1+nnodes*4, 4)
        grids = int_data[:, index]
        eulers_x = float_data[:, index+1]
        eulers_y = float_data[:, index+2]
        eulers_z = float_data[:, index+3]

        if element not in trmbd.eulersx:
            #isubcase = time_steps[subcase]
            ntimes = time_steps.shape[0]
            trmbd.eulersx[element] = np.empty([ntimes, nelements, eulers_x.shape[1]])
            trmbd.eulersy[element] = np.empty([ntimes, nelements, eulers_y.shape[1]])
            trmbd.eulersz[element] = np.empty([ntimes, nelements, eulers_z.shape[1]])
        if element not in trmbd.nodes:
            nodes = np.empty([nelements, grids.shape[1] + 1], dtype='int32')
            nodes[:, 0] = eids
            nodes[:, 1:] = grids
            trmbd.nodes[element] = nodes

        itime = tstep_to_index_mapper[time_step]
        trmbd.eulersx[element][itime, :, :] = eulers_x
        trmbd.eulersy[element][itime, :, :] = eulers_y
        trmbd.eulersz[element][itime, :, :] = eulers_z

    assert len(subcases_set) == 1, subcases_set
    asdf
    return

def reshape_trmbd(element_name: str, nnodes: int, int_data, float_data):
    ndata_per_element = 1 + nnodes + 3 * nnodes  # 1+4*(nnnodes) = 1+4*2 = 9
    n_elements = int_data.shape[0] // ndata_per_element  # elid + 2 grid + 3*2 euler angles
    assert int_data.shape[0] % ndata_per_element == 0, f'{element_name}: nnodes={nnodes} int_data.shape={int_data.shape}'

    int_data = int_data.reshape(n_elements, ndata_per_element)
    float_data = float_data.reshape(n_elements, ndata_per_element)
    return n_elements, int_data, float_data

