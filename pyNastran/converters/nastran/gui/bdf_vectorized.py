from io import StringIO
import numpy as np
import vtk

#VTK_TRIANGLE = 5
#VTK_QUADRATIC_TRIANGLE = 22

#VTK_QUAD = 9
#VTK_QUADRATIC_QUAD = 23

#VTK_TETRA = 10
#VTK_QUADRATIC_TETRA = 24

#VTK_WEDGE = 13
#VTK_QUADRATIC_WEDGE = 26

#VTK_HEXAHEDRON = 12
#VTK_QUADRATIC_HEXAHEDRON = 25

def add_vectorized_elements(model, nelements: int, idtype: str, log):
    dim_array = np.full(nelements, -1, dtype='int32')
    pids_array = np.zeros(nelements, 'int32')
    nnodes_array = np.full(nelements, -1, dtype='int32')
    #mids = np.zeros(nelements, 'int32')
    #material_coord = np.zeros(nelements, 'int32')
    #min_interior_angle = np.zeros(nelements, 'float32')
    #max_interior_angle = np.zeros(nelements, 'float32')
    #dideal_theta = np.zeros(nelements, 'float32')
    #max_skew_angle = np.zeros(nelements, 'float32')
    #max_warp_angle = np.zeros(nelements, 'float32')
    #max_aspect_ratio = np.zeros(nelements, 'float32')
    #area = np.zeros(nelements, 'float32')
    #area_ratio = np.zeros(nelements, 'float32')
    #taper_ratio = np.zeros(nelements, 'float32')
    #min_edge_length = np.zeros(nelements, 'float32')

    if len(model.ctria3):
        model.ctria3.quality()
    #if len(model.tria6):
        #model.tria6.quality()
    #if len(model.quad4):
        #model.quad4.quality()
    #if len(model.cquad8):
        #model.cquad8.quality()
    #if len(model.cquad):
        #model.cquad.quality()

    nids_list = []  # type: list[int]
    unused_ieid = 0
    unused_cell_offset = 0

    eids_array = np.zeros(nelements, dtype=idtype)
    cell_types_array = np.zeros(nelements, dtype=idtype)
    cell_offsets_array = np.zeros(nelements, dtype=idtype)

    results = {
        'pid' : pids_array,
        'eid' : eids_array,
        'nnodes' : nnodes_array,
        'dim' : dim_array,
    }

    cell_type_point = vtk.vtkVertex().GetCellType()
    cell_type_line = vtk.vtkLine().GetCellType()
    cell_type_tri3 = 5
    cell_type_tri6 = 22
    cell_type_quad4 = 9
    cell_type_quad8 = 23
    cell_type_tetra4 = 10
    cell_type_tetra10 = 24
    cell_type_pyram5 = vtk.vtkPyramid().GetCellType()
    cell_type_pyram13 = vtk.vtkQuadraticPyramid().GetCellType()
    cell_type_hexa8 = 12
    cell_type_hexa20 = 25
    cell_type_penta6 = 13
    cell_type_penta15 = 26

    unused_all_eids = model.elements2.eids
    #print('type(eids) =', type(all_eids)) # list
    #print('all_eids =', all_eids)

    #ncelas1 = len(model.celas1)
    #ncelas2 = len(model.celas2)
    #ncelas3 = len(model.celas3)
    #ncelas4 = len(model.celas4)

    #ncdamp1 = len(model.cdamp1)
    #ncdamp2 = len(model.cdamp2)
    #ncdamp3 = len(model.cdamp3)
    #ncdamp4 = len(model.cdamp4)
    #ncdamp5 = len(model.cdamp5)

    #nconrod = len(model.conrod)
    #ncrod = len(model.crod)
    #nctube = len(model.ctube)

    #ncbar = len(model.cbar)
    #ncbeam = len(model.cbeam)

    #ncshear = len(model.cshear)

    #nctria3 = len(model.ctria3)
    #ncquad4 = len(model.cquad4)
    #nctria6 = len(model.ctria6)
    #ncquad8 = len(model.cquad8)
    #ncquad = len(model.cquad)

    #nctetra4 = len(model.ctetra4)
    #ncpenta6 = len(model.cpenta6)
    #nchexa8 = len(model.chexa8)
    #ncpyram5 = len(model.cpyram5)
    #nsolids = nctetra4 + ncpenta6 + nchexa8

    ieid0 = 0
    cell_offset0 = 0
    nids_list = []

    nodes = model.nodes

    ieid0, cell_offset0 = map_elements_vectorized_fill_spring(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.celas1, cell_type_line, cell_type_point)
    ieid0, cell_offset0 = map_elements_vectorized_fill_spring(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.celas2, cell_type_line, cell_type_point)
    ieid0, cell_offset0 = map_elements_vectorized_fill_spring(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.celas3, cell_type_line, cell_type_point)
    ieid0, cell_offset0 = map_elements_vectorized_fill_spring(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.celas4, cell_type_line, cell_type_point)

    ieid0, cell_offset0 = map_elements_vectorized_fill_spring(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cdamp1, cell_type_line, cell_type_point)
    ieid0, cell_offset0 = map_elements_vectorized_fill_spring(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cdamp2, cell_type_line, cell_type_point)
    ieid0, cell_offset0 = map_elements_vectorized_fill_spring(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cdamp3, cell_type_line, cell_type_point)
    ieid0, cell_offset0 = map_elements_vectorized_fill_spring(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cdamp4, cell_type_line, cell_type_point)

    ieid0, cell_offset0 = map_elements_vectorized_fill_spring(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cvisc, cell_type_line, cell_type_point)

    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.plotel, cell_type_line, nnodesi=2, dimi=2)

    ieid0, cell_offset0 = map_elements_vectorized_fill_spring(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cbush, cell_type_line, cell_type_point)

    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.conrod, cell_type_line, nnodesi=2, dimi=2)
    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.crod, cell_type_line, nnodesi=2, dimi=2)
    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.ctube, cell_type_line, nnodesi=2, dimi=2)

    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cbar, cell_type_line, nnodesi=2, dimi=1)

    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cbeam, cell_type_line, nnodesi=2, dimi=1)

    #model.cbend
    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cshear, cell_type_quad4, nnodesi=4, dimi=2)

    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.ctria3, cell_type_tri3, nnodesi=3, dimi=2)
    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cquad4, cell_type_quad4, nnodesi=4, dimi=2)

    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.ctria6, cell_type_tri6, nnodesi=6, dimi=2,
        allow0=True, cell_type_allow=cell_type_tri3)
    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cquad8, cell_type_quad8, nnodesi=8, dimi=2,
        allow0=True, cell_type_allow=cell_type_quad8)
    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cquad, cell_type_quad8, nnodesi=8, dimi=2,
        allow0=True, cell_type_allow=cell_type_quad8)

    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.ctriar, cell_type_tri6, nnodesi=6, dimi=2,
        allow0=True, cell_type_allow=cell_type_tri3)
    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cquadr, cell_type_tri6, nnodesi=6, dimi=2,
        allow0=True, cell_type_allow=cell_type_quad4)

    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.ctetra4, cell_type_tetra4, nnodesi=4, dimi=3)
    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cpenta6, cell_type_penta6, nnodesi=8, dimi=3)
    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.chexa8, cell_type_hexa8, nnodesi=8, dimi=3)
    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cpyram5, cell_type_pyram5, nnodesi=8, dimi=3)

    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.ctetra10, cell_type_tetra10, nnodesi=4, dimi=3,
        allow0=True, cell_type_allow=cell_type_tetra4)
    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cpenta15, cell_type_penta15, nnodesi=8, dimi=3,
        allow0=True, cell_type_allow=cell_type_penta6)
    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.chexa20, cell_type_hexa20, nnodesi=8, dimi=3,
        allow0=True, cell_type_allow=cell_type_hexa8)
    ieid0, cell_offset0 = map_elements_vectorized_fill(
        log,
        ieid0, cell_offset0,
        nodes, nids_list,
        eids_array, pids_array, nnodes_array, dim_array,
        cell_types_array, cell_offsets_array,
        model.cpyram13, cell_type_pyram13, nnodesi=8, dimi=3,
        allow0=True, cell_type_allow=cell_type_pyram5)

    # model.chbdyg
    # model.chbdye
    # model.chbdyp
    return cell_types_array, cell_offsets_array, nids_list, eids_array, results

def map_elements_vectorized_fill(log,
                                 ieid0: int, cell_offset0: int,
                                 nodes, nids_list: list[int],
                                 eids_array, pids_array, nnodes_array, dim_array,
                                 cell_types_array, cell_offsets_array,
                                 model_obj, cell_type: int, nnodesi=None, dimi=None,
                                 allow0=False, cell_type_allow=None) -> tuple[int, int]:
    """helper method for ``map_elements_vectorized``"""
    assert nnodesi is not None
    assert dimi is not None

    nelems = len(model_obj)
    #print('nelements =', nelems)
    if nelems:
        ieid0_in = ieid0
        elem = model_obj
        elem.make_current()
        try:
            nnodes = model_obj.nids.shape[1]
        except IndexError:
            log.info('%s; nelem=%s nnodes=???' % (
                model_obj.card_name, nelems))
            raise

        log.debug('%s; nelem=%s nnodes=%s' % (
            model_obj.card_name, nelems, nnodes))
        eid = elem.eid
        ieid = np.arange(ieid0, ieid0 + nelems)
        dim = np.full(nelems, dimi, dtype='int32')
        dim_array[ieid] = dim

        nids = elem.nids
        #print('nids =', nids)
        if hasattr(elem, 'pid'):
            pid = elem.pid
            pids_array[ieid] = pid

        try:
            inids = nodes.get_node_index(nids, allow0=allow0)
        except Exception:
            bdf_file = StringIO()
            model_obj.write_card(bdf_file=bdf_file)
            print(bdf_file.getvalue())
            raise
        if inids.min() == -1:
            nrequired = elem.nrequired
            #print('nrequired =', nrequired)
            #print('nids.shape =', nids.shape)
            #print('inids.shape =', inids.shape)
            inids = inids[:, :nrequired]
            nnodesi = nrequired
            assert cell_type_allow is not None
            unused_eid_type = np.full(nelems, cell_type_allow, dtype='int32')
        else:
            unused_eid_type = np.full(nelems, cell_type, dtype='int32')

        nnodes = np.full((nelems, 1), nnodesi, dtype='int32')
        if len(ieid) != len(nnodes):
            msg = 'len(ieid)=%s len(nnodes)=%s' % (len(ieid), len(nnodes))
            raise RuntimeError(msg)

        # nnodes is a 2D column matrix and nnodes_array is a 1D array,
        # so we ravel it
        nnodes_array[ieid] = nnodes.ravel()
        assert -1 not in inids, inids
        nnodes_inids = np.hstack([nnodes, inids])
        nids_list.append(nnodes_inids)
        log.debug('  ieid = %s' % ieid)
        assert len(ieid) == nelems

        nnodes[0] = -1
        cumsum = cell_offset0 + np.cumsum(nnodes + 1)
        assert len(ieid) == len(cumsum)
        try:
            eids_array[ieid] = eid
        except IndexError:
            # this hits when an element doesn't correctly implement make_current()
            log.error('nelems = %s' % nelems)
            log.error('card_count = %s' % model_obj.model.card_count)
            log.error('eids_array = %s' % eids_array)
            log.error('eids_array.shape = %s' % str(eids_array.shape))
            log.error('ieid.max() = %s' % ieid.max())
            raise
        cell_types_array[ieid] = cell_type
        cell_offsets_array[ieid] = cumsum
        ieid0 += nelems
        cell_offset0 += cumsum[-1]
        #print('nids_list[%s] = %s' % (model_obj.card_name, nids_list[-1]))
        #print('cell_offsets_array =', cell_offsets_array)
        log.debug("  ieid0_in=%s ieid0_out=%s" % (ieid0_in, ieid0))
    return ieid0, cell_offset0

def map_elements_vectorized_fill_spring(log,
                                        ieid0: int, cell_offset0: int,
                                        nodes, nids_list,
                                        eids_array, pids_array, nnodes_array, dim_array,
                                        cell_types_array, cell_offsets_array,
                                        model_obj, cell_type_line: int, cell_type_point: int) -> tuple[int, int]:
    """helper method for ``map_elements_vectorized``"""
    nelems = len(model_obj)
    if nelems:
        ieid0_in = ieid0
        dimi = 1
        nnodesi = 2
        elem = model_obj
        elem.make_current()
        log.debug('%s; nelem=%s nnodes=%s' % (
            model_obj.card_name, nelems, model_obj.nids.shape[1]))

        ieid = np.arange(ieid0, ieid0 + nelems)
        eid_type = np.full(nelems, cell_type_line, dtype='int32')
        dim = np.full(nelems, dimi, dtype='int32')
        nnodes = np.full((nelems, 1), nnodesi, dtype='int32')
        dim_array[ieid] = dim

        eid = elem.eid
        eids_array[ieid] = eid

        nids = elem.nids
        inonzero = np.where(nids[:, 1] > 0)[0]
        izero = np.where(nids[:, 1] == 0)[0]

        if len(inonzero):
            cell_types_array[ieid[inonzero]] = cell_type_line

        nzero = len(izero)
        if nzero:
            #raise NotImplementedError('izero')
            eid_type[izero] = cell_type_point
            nnodes[izero] = 1
            cell_types_array[ieid[izero]] = cell_type_point
        assert len(ieid) == len(nnodes), 'len(ieid)=%s len(nnodes)=%s' % (len(ieid), nnodes)
        nnodes_array[ieid] = nnodes.ravel()

        if hasattr(elem, 'pid'):
            pid = elem.pid
            pids_array[ieid] = pid
        inids = nodes.get_node_index(nids, allow0=True)

        nnodes_inids = np.hstack([nnodes, inids])#.ravel()
        #print("  nnodes_inids=", nnodes_inids)
        if nzero:
            nnodes_inids_ravel = nnodes_inids.flatten()
            index_minus1 = np.where(nnodes_inids_ravel != -1)[0]
            nnodes_inids_cleaned = nnodes_inids_ravel[index_minus1]
            assert -1 not in nnodes_inids_cleaned, nnodes_inids_cleaned
            nids_list.append(nnodes_inids_cleaned)
        else:
            nids_list.append(nnodes_inids.ravel())

        nnodes[0] = -1
        cumsum = cell_offset0 + np.cumsum(nnodes + 1)
        #print("nnodes =", nnodes)
        #print('offset0=%s cumsum=%s' % (cell_offset0, cumsum))
        log.debug('  ieid = %s' % ieid)
        assert len(ieid) == nelems
        assert len(ieid) == len(cumsum)
        cell_offsets_array[ieid] = cumsum
        ieid0 += nelems
        log.debug("  ieid0_in=%s ieid0_out=%s" % (ieid0_in, ieid0))
        cell_offset0 += cumsum[-1]
        #print('nids_list[%s] = %s' % (model_obj.card_name, nids_list[-1]))
        #print('cell_offsets_array =', cell_offsets_array)
    return ieid0, cell_offset0
