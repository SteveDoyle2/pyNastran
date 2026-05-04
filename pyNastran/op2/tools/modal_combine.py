import os
from typing import Optional

import numpy as np
from cpylog import SimpleLogger
from pyNastran.f06.parse_flutter import make_flutter_response, FlutterResponse
from pyNastran.utils import PathLike  #, print_bad_path
from pyNastran.op2.op2 import OP2
from pyNastran.op2.tables.oug.oug_eigenvectors import ComplexEigenvectorArray
from .solution_combination import _read_op2


def run_flutter_combination(op2_filename: PathLike | OP2,
                            f06_filename: PathLike,
                            f06_units: str,
                            out_units: str,
                            op2_filename_out: PathLike,
                            save_individual: bool=False,
                            nmodes_to_keep: Optional[int]=None) -> None:

    responses, mass_lama = make_flutter_response(
        f06_filename, f06_units=f06_units, out_units=out_units,
        use_rhoref=False, read_flutter=True, make_alt=True, log=None)
    # print(responses)

    subcases = list(responses)
    assert len(subcases) == 1, subcases
    subcase0 = subcases[0]
    response: FlutterResponse = responses[subcase0]
    nvel = response.eigenvector.shape[0] # (2, 81, 81)
    log = SimpleLogger(level='info')

    # log.info(f'response.eigenvector={response.eigenvector}')
    if len(response.eigenvector) == 0:
        raise RuntimeError(f'f06={str(f06_filename)} has no eigenvectors')

    # print('response.eigr_eigi_velocity.shape = ', response.eigr_eigi_velocity.shape)
    vels_out = response.eigr_eigi_velocity[:, 0, 2]
    vel_units_out = response.out_units['eas']
    vels_out = vels_out.round()
    log.info(f'nvel={nvel}; vel={vels_out} {vel_units_out}')
    log.info(f'nmodes_to_keep = {nmodes_to_keep}')
    assert nvel > 0, nvel

    mode = None
    model = _read_op2(op2_filename,
                      mode, log,
                      include_results=['eigenvectors'],
                      include_complex_modes=False,
                      subcases=subcases)
    mode = model._nastran_format
    # for obj in model.eigenvectors.values():
    #     print(''.join(obj.get_stats()))

    model_out = None
    # base, ext = os.path.splitext(str(f06_filename))
    dirname = os.path.dirname(str(f06_filename))

    base, ext = os.path.splitext(op2_filename_out)

    for ivel, vel in enumerate(vels_out):
        log.info(f'creating ivel={ivel}/{nvel}; vel={vels_out} {vel_units_out}')
        # (nmodes, nmodes) where each column is a complex mode shape
        eigenvector = response.eigenvector[ivel, :nmodes_to_keep, :].T

        eigrs_eigis = response.eigr_eigi_velocity[ivel, :nmodes_to_keep, :]
        csv_filename1 = os.path.join(dirname, f'ivel{ivel}_eigenvector.csv')
        csv_filename2 = os.path.join(dirname, f'ivel{ivel}_eigr_eigi.csv')
        np.savetxt(csv_filename1, eigenvector)
        np.savetxt(csv_filename2, eigrs_eigis)

        # looks like the eigenvalue order is pseudo-random :(
        eigrs = eigrs_eigis[:, 0]
        eigis = eigrs_eigis[:, 1]
        # print(eigrs_eigis.shape)
        # print(eigenvector)
        if save_individual:
            # print('creating model_out')
            model_out = OP2(mode=mode)
            model_out._nastran_format = mode
            model_out.set_revision_from_model(model)
        model_out = run_modal_combination(
            model, eigenvector, eigrs, eigis,
            subcase0=subcase0,
            log=log, model_out=model_out)
        log.debug(f'ivel={ivel} vel={vel} {vel_units_out}')
        if save_individual:
            op2_filename_outi = f'{base}_ivel{ivel}_{vel:.0f}_{vel_units_out}{ext}'
            model_out.log.info(f'writing {str(op2_filename_outi)}')
            model_out.write_op2(op2_filename_outi)

    assert model_out is not None, model_out
    if not save_individual:
        model_out.log.info(f'writing {str(op2_filename_out)}')
        model_out.set_revision_from_model(model)
        model_out.write_op2(op2_filename_out)


def run_modal_combination(op2_filename: PathLike | OP2,
                          combination_factors: PathLike | np.ndarray,
                          eigrs: np.ndarray, eigis: np.ndarray,
                          model_out: Optional[OP2]=None,
                          mode: Optional[str]=None,
                          subcase0: int=0,
                          subcase_out: int=0,
                          log: Optional[SimpleLogger]=None) -> OP2:
    """
      _count : 0
      _encoding : 'utf8'
      _nnodes : 100
      _ntotals : [...]
      _times : array([...], dtype=int32)
      _times_dtype : 'float32'
      acoustic_flag : 0
      analysis_code : 9
      analysis_fmt : 'int32'
      approach_code : 92
      class_name : 'ComplexEigenvectorArray'
      data   : array(..., dtype=float64)
      data_code : {
          '_encoding': 'utf8', 'load_as_h5': False, 'size': 4, 'is_msc': True,
          'table_name': 'BOUGV1', 'approach_code': 92, 'isubcase': 3, 'table_code': 7,
          'tCode': 7, 'sort_code': 0, 'sort_method': 1, 'device_code': 2,
          'analysis_code': 9, 'sort_bits': SortBits(is_complex=1, is_sort2=0,
          is_random=0), 'random_code': 0, 'format_code': 2, 'num_wide': 14,
          'acoustic_flag': 0, 'thermal': 0, 'nonlinear_factor': np.int32(80),
          'name': 'mode', 'mode': 80, 'eigr': 0.0, 'eigi': 436.6650085449219,
          'data_names': ['mode', 'eigr', 'eigi'], '_times_dtype': 'float32',
          'thermal_bits': [0, 0, 0, 0, 0], 'subtitle': '', 'subtitle_original': '',
          'pval_step': '', 'superelement_adaptivity_index': '', 'label': '',
          'title': '', 'result_name': 'eigenvectors', '_count': 0}
      data_frame : None
      data_names : ['mode', 'eigr', 'eigi']
      dataframe : None
      device_code : 2
      dt     : None
      eigi   : 436.6650085449219
      eigis  : [0.0,  ...]
      eigr   : 0.0
      eigrs  : [0.0, ...]
      format_code : 2
      gridtype_str : array(..., dtype='<U1')
      h5_file : None
      headers : ['t1', 't2', 't3', 'r1', 'r2', 'r3']
      is_built : True
      is_cid : False
      is_complex : True
      is_msc : True
      is_real : False
      is_sort1 : True
      is_sort2 : False
      isubcase : 3
      itime  : 81
      itotal : 0
      label  : ''
      load_as_h5 : False
      mode   : 80
      modes  : array(..., dtype=int32)
      name   : 'mode'
      node_gridtype : array([], shape=(0, 2), dtype=float64)
      nonlinear_factor : np.int32(80)
      ntimes : 81
      ntotal : ...
      num_wide : 14
      ogs    : None
      pval_step : ''
      random_code : 0
      result_name : 'eigenvectors'
      size   : 4
      sort_bits : SortBits(is_complex=1, is_sort2=0, is_random=0)
      sort_code : 0 -> 1
      sort_method : 1
      subtitle : ''
      subtitle_original : ''
      superelement_adaptivity_index : ''
      tCode  : 7 -> 1007
      table_code : 7
      table_name : 'BOUGV1'
      table_name_str : 'BOUGV1'
      thermal : 0
      thermal_bits : [0, 0, 0, 0, 0]
      title  : ''

    Parameters
    ----------
    op2_filename : PathLike | OP2
        the op2 filename
    combination_factors : np.ndarray or file
        the (nmodes, nmodes) combinations
    eigrs : np.ndarray
        the real eigenvalues for the new modes
    eigis : np.ndarray
        the imaginary eigenvalues for the new modes
    model_out : OP2; default=None -> OP2
        lets you create a combined op2
    mode : str | None; default=None
        msc, nx, optistruct
    subcase0 : int; default=0 -> first subcase
        the modal eigenvectors subcase
    subcase_out : int; default=0 -> max+1
        the output subcase id
    """
    subcases = None if subcase0 == 0 else [subcase0]
    include_results = ['eigenvectors']

    if isinstance(combination_factors, PathLike):
        mpfs = np.loadtxt(combination_factors, delimiter=',')
    else:
        # nmodes, 2
        assert isinstance(combination_factors, np.ndarray), combination_factors
        mpfs = combination_factors
    assert len(mpfs.shape) == 2, mpfs.shape

    # assert A.shape[1] == 1, mpfs.shape
    nmodes_source, nmodes_out = mpfs.shape
    assert mpfs.dtype == np.complex128, mpfs.dtype
    #-----------------------------------------------

    model = _read_op2(op2_filename,
                      mode, log,
                      # exclude_results=None,
                      include_results=include_results,
                      include_complex_modes=False,
                      subcases=subcases)
    log = model.log
    if model_out is None:
        mode = model._nastran_format
        assert mode is not None, mode
        model_out = OP2(mode=mode, log=log)
        assert model_out._nastran_format is not None, mode

    log.level = 'debug'
    if subcase0 == 0:
        subcases = list(model.eigenvectors)
        assert len(subcases) == 1, subcases
        subcase0 = subcases[0]

    if subcase_out == 0:
        if len(model_out.eigenvectors) > 0:
            subcase_out = max(model_out.eigenvectors) + 1
        else:
            subcase_out = max(model.eigenvectors) + 1

    obj0 = model.eigenvectors[subcase0]

    if obj0.table_name == 'BOPHIG':
        table_name = 'BOUGV1'
    elif obj0.table_name == 'BOPHIG1':
        table_name = 'BOUGV1'
    else:
        raise RuntimeError(f'Unknown table name: {obj0.table_name}')
    # table_name = obj0.table_name
    assert obj0.analysis_code == 2, obj0.get_stats()

    phi = obj0.get_phi()  # (ndof, nmodes)
    assert phi.shape[1] == nmodes_source

    # nmodes
    modes_out = np.arange(nmodes_out, dtype='int32')
    modal_disp = phi @ mpfs
    data = obj0.phi_to_data(modal_disp)
    assert isinstance(data[0, 0, 0], np.complex128), data.dtype

    data = data.astype(np.complex64)
    eigrs = np.asarray(eigrs, dtype='float32')
    eigis = np.asarray(eigis, dtype='float32')

    obj = ComplexEigenvectorArray.add_complex_modes_case(
        table_name, obj0.node_gridtype, data, subcase_out,
        modes_out, eigrs, eigis,
        table_code=7,  # modes,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=obj0.title, subtitle='', label='')
    # print(obj.data.shape)
    obj.subcase = subcase_out
    # print(''.join(obj.get_stats()))
    model_out.eigenvectors[subcase_out] = obj
    return model_out

def cmd_line_flutter_combine(argv=None, quiet: bool=False, log=None):
    if argv is None:
        argv = sys.argv

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('modal_combine', help='')
    parser.add_argument('op2_filename', help='input op2')
    parser.add_argument('f06_filename', help='input f06', default=None)
    parser.add_argument('-o', '--out', help='output op2', default=None)
    parser.add_argument('--in_units', help='input units (english_in, english_ft, english_kt, si, si_mm)', default='si')
    parser.add_argument('--out_units', help='output units (english_in, english_ft, english_kt, si, si_mm)', default='si')
    parser.add_argument('--nmodes', help='number of modes to keep', default=None)
    parser.add_argument('--split', help='split the op2 into separate files to reduce RAM usage', action='store_true')
    # parser.add_argument('--out_units', help='output units (english_in, english_ft, english_kt, , si, si_mm)', default='si')

    args = parser.parse_args(args=argv[1:])
    if not quiet:  # pragma: no cover
        print(args)
    op2_filename = args.op2_filename
    f06_filename = args.f06_filename
    op2_filename_out = args.out
    f06_in_units = args.in_units
    out_units = args.out_units
    save_individual = args.split
    if args.nmodes is not None:
        nmodes_to_keep = int(args.nmodes)

    base, ext = os.path.splitext(op2_filename)
    if f06_filename is None:
        f06_filename = f'{base}.f06'
    if op2_filename_out is None:
        op2_filename_out = f'{base}.flutter_combined{ext}'

    print(op2_filename, op2_filename_out)

    run_flutter_combination(
        op2_filename, f06_filename,
        f06_in_units, out_units, op2_filename_out,
        nmodes_to_keep=nmodes_to_keep,
        save_individual=save_individual,
    )
