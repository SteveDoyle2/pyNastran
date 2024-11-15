test_bdf demo
=============

In this demo, we’ll show off test_bdf

.. code:: ipython3

    from IPython.display import HTML as html_print
    from pyNastran.bdf.bdf import BDF, read_bdf, CaseControlDeck


.. code:: ipython3

    model = BDF()

    # add_grid(nid, xyz, cp=0, cd=0, ps='', seid=0)
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [1., 0., 0.])
    model.add_grid(3, [1., 1., 0.])
    model.add_grid(4, [0., 1., 0.])

    eid = 10
    pid = 100
    mid = 1000
    # add_cbar(eid, pid, nids, x, g0, offt='GGG', pa=0, pb=0, wa=None, wb=None)
    model.add_cbar(eid,   pid, [1, 2], [0., 0., 1.], None, offt='GGG')
    model.add_cbar(eid+1, pid, [2, 3], [0., 0., 1.], None, offt='GGG')
    model.add_cbar(eid+2, pid, [3, 4], [0., 0., 1.], None, offt='GGG')
    model.add_cbar(eid+3, pid, [4, 1], [0., 0., 1.], None, offt='GGG')

    eid_cquad4 = 15
    pid_pshell = 101
    # add_pshell(pid, mid1=None, t=None, mid2=None, twelveIt3=1.0,
    #            mid3=None, tst=0.833333,
    #            nsm=0.0, z1=None, z2=None, mid4=None)
    model.add_pshell(pid_pshell, mid1=mid, t=0.1, mid2=mid, mid3=mid)
    model.add_cquad4(eid_cquad4, pid_pshell, [1, 2, 3, 4])

    dim = [3., 3., 1., 1.] # TODO: should be [1., 2., 3., 4.]
    # add_pbarl(pid, mid, Type, dim, group='MSCBML0', nsm=0.0)
    pbarl = model.add_pbarl(pid, mid, 'BOX', dim, nsm=0.0)
    pbarl.validate()

    E = 3.e7
    G = None
    nu = 0.3
    mat = model.add_mat1(mid, E, G, nu)

    spc_id = 1
    nids = 1
    # add_spc1(conid, components, nodes
    model.add_spc1(spc_id, 123456, nids)


    dresp_id = 100
    label = 'resp1'
    response_type = 'STRESS'
    property_type = 'PSHELL'
    pid = 3
    atta = 9 # von mises upper surface stress
    region = None
    attb = None
    atti = [pid_pshell]
    # add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)
    model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

    dresp_id += 1
    atta = 17 # von mises lower surface stress
    model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

    # add_dconstr(oid, dresp_id, lid=-1e+20, uid=1e+20, lowfq=0.0, highfq=1e+20)
    dconstr_id = 10000
    model.add_dconstr(dconstr_id, dresp_id, lid=-35000., uid=35000.)

    dresp_id += 1
    dresp = model.add_dresp1(dresp_id, 'WT', 'WEIGHT', None, None, None, None, None)
    dresp.validate()

    oid = 1000
    dvids = 1
    coeffs = 1.
    # add_dvprel1(oid, prop_type, pid, pname_fid, dvids, coeffs,
    #             p_min=None, p_max=1e+20, c0=0.0)
    model.add_dvprel1(oid, 'PSHELL', pid_pshell, 'T', dvids, coeffs)

    # add_desvar(desvar_id, label, xinit, xlb=-1e+20, xub=1e+20,
    #            delx=None, ddval=None)
    model.add_desvar(1, 'DIM1', 0.1, xlb=1e-5)
    model.add_desvar(2, 'DIM2', 0.2, xlb=1e-5)
    model.add_desvar(3, 'DIM3', 0.3, xlb=1e-5)
    model.add_desvar(4, 'DIM4', 0.4, xlb=1e-5)
    model.add_desvar(5, 'DV5', 0.1, xlb=1e-5)

    #model.add_dlink(6)

    eid = 10 # TODO: remove
    load_id = 1
    # add_pload4(sid, eids, pressures, g1=None, g34=None,
    #            cid=0, nvector=None, surf_or_line='SURF', line_load_dir='NORM')
    pload4 = model.add_pload4(load_id, [eid_cquad4], [1., None, None, None],
                              comment=' load')
    #print(pload4.get_stats())

    eid = 10 # TODO: should be 100
    scale = 'LE' # TODO: should be 100.
    # add_pload1(sid, eid, load_type, scale, x1, p1, x2=None, p2=None)
    model.add_pload1(load_id, eid, 'FZ', scale, 0., 1.)  # TODO: change atti to None

    # add_eigrl(sid, v1=None, v2=None, nd=None, msglvl=0, maxset=None, shfscl=None,
    #           norm=None, options=None, values=None)
    eigrl = model.add_eigrl(42, nd=42)

    model.sol = 103  # start=103
    cc = CaseControlDeck([
        'DESOBJ = 102',  # DRESP1
        'DESSUB = %s' % dconstr_id,  # DCONSTR
        'SUBCASE 1',
        '  METHOD = 42',  # TODO: remove
        '  LOAD = %s' % load_id,  # TODO: remove
        '  SPC = %s' % spc_id,
        #'  TRIM = 42',  # TODO: add
        'ANALYSIS = SAERO',
    ])
    #print(cc)
    model.case_control_deck = cc
    model.validate()

    # rerun between each change
    # 1. change SOL=103 -> SOL=144
    model.sol = 144
    # 2. add the trim in the case control deck
    #help(model.add_trim)
    mach = 0.8
    q = 100.
    labels = ['Z']
    uxs = [2.5]
    trim = model.add_trim(42, mach, q, labels, uxs, aeqr=1.0, trim_type=1)
    print(trim)
    # 3. add a trim card
    # x. change to SOL=200
    print(model.trims)
    model.write_bdf('junk.bdf')
    !cat junk.bdf
    print('----------------------------------------------------------------------------------------------------')


.. parsed-literal::

    TRIM          42      .8    100.       Z     2.5

    {42: TRIM          42      .8    100.       Z     2.5
    }



.. raw:: html

    <text style=color:blue>DEBUG:   write_mesh.py:145            ---starting BDF.write_bdf of junk.bdf---
    </text>


.. parsed-literal::

    $pyNastran: version=msc
    $pyNastran: punch=False
    $pyNastran: encoding=utf-8
    $pyNastran: nnodes=4
    $pyNastran: nelements=5
    $EXECUTIVE CONTROL DECK
    SOL 144
    CEND
    $CASE CONTROL DECK
    DESOBJ = 102
    DESSUB = 10000
    SUBCASE 1
        ANALYSIS = SAERO
        LOAD = 1
        METHOD = 42
        SPC = 1
    BEGIN BULK
    $NODES
    GRID           1              0.      0.      0.
    GRID           2              1.      0.      0.
    GRID           3              1.      1.      0.
    GRID           4              0.      1.      0.
    $ELEMENTS
    CBAR          10     100       1       2      0.      0.      1.
    CBAR          11     100       2       3      0.      0.      1.
    CBAR          12     100       3       4      0.      0.      1.
    CBAR          13     100       4       1      0.      0.      1.
    CQUAD4        15     101       1       2       3       4
    $PROPERTIES
    PBARL        100    1000             BOX
                  3.      3.      1.      1.      0.
    PSHELL       101    1000      .1    1000            1000
    $MATERIALS
    MAT1        1000    3.+7              .3
    $LOADS
    $ load
    PLOAD4         1      15      1.
                   0      0.      0.      0.    SURF
    PLOAD1         1      10      FZ      LE      0.      1.      0.      1.
    $DYNAMIC
    EIGRL         42                      42
    $STATIC AERO
    TRIM          42      .8    100.       Z     2.5
    $SPCs
    SPC1           1  123456       1
    $OPTIMIZATION
    DCONSTR    10000     101 -35000.  35000.
    DESVAR         1   DIM1       .1  .00001
    DESVAR         2   DIM2       .2  .00001
    DESVAR         3   DIM3       .3  .00001
    DESVAR         4   DIM4       .4  .00001
    DESVAR         5    DV5       .1  .00001
    DRESP1       100   resp1  STRESS  PSHELL               9             101
    DRESP1       101   resp1  STRESS  PSHELL              17             101
    DRESP1       102      WT  WEIGHT                                     ALL
    DVPREL1     1000  PSHELL     101       T
                   1      1.
    ----------------------------------------------------------------------------------------------------


.. parsed-literal::

    c:\python37\lib\site-packages\IPython\utils\_process_win32.py:131: ResourceWarning: unclosed file <_io.BufferedWriter name=5>
      return process_handler(cmd, _system_body)
    ResourceWarning: Enable tracemalloc to get the object allocation traceback
    c:\python37\lib\site-packages\IPython\utils\_process_win32.py:131: ResourceWarning: unclosed file <_io.BufferedReader name=6>
      return process_handler(cmd, _system_body)
    ResourceWarning: Enable tracemalloc to get the object allocation traceback
    c:\python37\lib\site-packages\IPython\utils\_process_win32.py:131: ResourceWarning: unclosed file <_io.BufferedReader name=7>
      return process_handler(cmd, _system_body)
    ResourceWarning: Enable tracemalloc to get the object allocation traceback


.. code:: ipython3

    from pyNastran.bdf.test.test_bdf import run_bdf as test_bdf
    model.write_bdf('junk.bdf')
    test_bdf('.', 'junk.bdf')



.. raw:: html

    <text style=color:blue>DEBUG:   write_mesh.py:145            ---starting BDF.write_bdf of junk.bdf---
    </text>


.. parsed-literal::

    debug = False
    bdf_model = junk.bdf



.. raw:: html

    <text style=color:green>INFO:    test_bdf.py:374              starting fem1
    </text>



.. raw:: html

    <text style=color:green>INFO:    test_bdf.py:841              starting fem2
    </text>



.. raw:: html

    <text style=color:orange>WARNING: test_bdf.py:863              PARAM,POST,0 is not supported by the OP2 reader
    </text>



.. raw:: html

    <text style=color:red>ERROR:   test_bdf.py:1162             A TRIM or DIVERG card is required for STATIC AERO - SOL 144
    SUBCASE 1
        ANALYSIS = SAERO
        DESOBJ = 102
        DESSUB = 10000
        LOAD = 1
        METHOD = 42
        SPC = 1

    </text>


::


    ---------------------------------------------------------------------------

    RuntimeError                              Traceback (most recent call last)

    <ipython-input-22-746e54da5a3e> in <module>
          1 from pyNastran.bdf.test.test_bdf import run_bdf as test_bdf
          2 model.write_bdf('junk.bdf')
    ----> 3 test_bdf('.', 'junk.bdf')


    c:\nasa\m4\formats\git\pynastran\pyNastran\bdf\test\test_bdf.py in run_bdf(folder, bdf_filename, debug, xref, check, punch, mesh_form, is_folder, print_stats, encoding, sum_load, size, is_double, hdf5, stop, nastran, post, dynamic_vars, quiet, dumplines, dictsort, run_extract_bodies, run_skin_solids, save_file_structure, nerrors, dev, crash_cards, safe_xref, run_pickle, stop_on_failure, log)
        342         run_pickle=run_pickle,
        343         stop_on_failure=stop_on_failure,
    --> 344         log=log,
        345     )
        346     return fem1, fem2, diff_cards


    c:\nasa\m4\formats\git\pynastran\pyNastran\bdf\test\test_bdf.py in run_and_compare_fems(bdf_model, out_model, debug, xref, check, punch, mesh_form, print_stats, encoding, sum_load, size, is_double, save_file_structure, stop, nastran, post, hdf5, dynamic_vars, quiet, dumplines, dictsort, nerrors, dev, crash_cards, safe_xref, run_extract_bodies, run_skin_solids, run_pickle, stop_on_failure, log)
        412                         encoding=encoding, debug=debug, quiet=quiet,
        413                         ierror=ierror, nerrors=nerrors,
    --> 414                         stop_on_failure=stop_on_failure, log=log)
        415
        416         diff_cards = compare(fem1, fem2, xref=xref, check=check,


    c:\nasa\m4\formats\git\pynastran\pyNastran\bdf\test\test_bdf.py in run_fem2(bdf_model, out_model, xref, punch, sum_load, size, is_double, mesh_form, safe_xref, encoding, debug, quiet, stop_on_failure, ierror, nerrors, log)
        878                 fem2, p0, sol_base, subcase_keys, subcases, sol_200_map,
        879                 ierror=ierror, nerrors=nerrors,
    --> 880                 stop_on_failure=stop_on_failure)
        881
        882     if mesh_form is not None:


    c:\nasa\m4\formats\git\pynastran\pyNastran\bdf\test\test_bdf.py in validate_case_control(fem2, p0, sol_base, subcase_keys, subcases, unused_sol_200_map, stop_on_failure, ierror, nerrors)
        922         ierror = check_case(
        923             sol_base, subcase, fem2, p0, isubcase, subcases,
    --> 924             ierror=ierror, nerrors=nerrors, stop_on_failure=stop_on_failure)
        925     return ierror
        926


    c:\nasa\m4\formats\git\pynastran\pyNastran\bdf\test\test_bdf.py in check_case(sol, subcase, fem2, p0, isubcase, subcases, ierror, nerrors, stop_on_failure)
       1106
       1107     elif sol == 144:
    -> 1108         ierror = _check_static_aero_case(fem2, log, sol, subcase, ierror, nerrors)
       1109     elif sol == 145:
       1110         ierror = _check_flutter_case(fem2, log, sol, subcase, ierror, nerrors)


    c:\nasa\m4\formats\git\pynastran\pyNastran\bdf\test\test_bdf.py in _check_static_aero_case(fem2, log, sol, subcase, ierror, nerrors)
       1161             sol, subcase)
       1162         log.error(msg)
    -> 1163         ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)
       1164     if fem2.aeros is None:
       1165         msg = 'An AEROS card is required for STATIC AERO - SOL %i; AEROS=%s' % (sol, fem2.aeros)


    c:\nasa\m4\formats\git\pynastran\pyNastran\bdf\test\test_bdf.py in stop_if_max_error(msg, error, ierror, nerrors)
        954     """if the error count is greater than nerrors, stop"""
        955     if ierror == nerrors:
    --> 956         raise error(msg)
        957     ierror += 1
        958     return ierror


    RuntimeError: A TRIM or DIVERG card is required for STATIC AERO - SOL 144
    SUBCASE 1
        ANALYSIS = SAERO
        DESOBJ = 102
        DESSUB = 10000
        LOAD = 1
        METHOD = 42
        SPC = 1
