from io import StringIO
import numpy as np
#from typing import List

def map_elements_vectorized_fill(log,
                                 ieid0: int, cell_offset0: int,
                                 nodes, nids_list,
                                 eids_array, pids_array, nnodes_array, dim_array,
                                 cell_types_array, cell_offsets_array,
                                 model_obj, cell_type, nnodesi=None, dimi=None,
                                 allow0=False, cell_type_allow=None):
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
        except:
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
                                        ieid0, cell_offset0,
                                        nodes, nids_list,
                                        eids_array, pids_array, nnodes_array, dim_array,
                                        cell_types_array, cell_offsets_array,
                                        model_obj, cell_type_line, cell_type_point):
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
