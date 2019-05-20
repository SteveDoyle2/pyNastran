
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
    
    model.sol = 200  # start with 103
    cc = CaseControlDeck([
        'DESOBJ = 102',  # DRESP1
        'DESSUB = %s' % dconstr_id,  # DCONSTR
        'SUBCASE 1',
        '  METHOD = 42',  # TODO: remove
        '  LOAD = %s' % load_id,  # TODO: remove
        '  SPC = %s' % spc_id,
        '  TRIM = 42',  # TODO: add
        'ANALYSIS = SAERO',
    ])
    #print(cc)
    model.case_control_deck = cc
    model.validate()
    
    
    model.write_bdf('junk.bdf')
    #!cat junk.bdf
    print('----------------------------------------------------------------------------------------------------')
    



.. raw:: html

    <text style=color:blue>DEBUG:   write_mesh.py:245            ---starting BDF.write_bdf of junk.bdf---
    </text>


.. parsed-literal::

    ----------------------------------------------------------------------------------------------------
    

.. code:: ipython3

    from pyNastran.bdf.test.test_bdf import run_bdf as test_bdf
    model.write_bdf('junk.bdf')
    test_bdf('.', 'junk.bdf')



.. raw:: html

    <text style=color:blue>DEBUG:   write_mesh.py:245            ---starting BDF.write_bdf of junk.bdf---
    </text>


.. parsed-literal::

    debug = False
    bdf_model = junk.bdf
    


.. raw:: html

    <text style=color:green>INFO:    test_bdf.py:347              starting fem1
    </text>



.. raw:: html

    <text style=color:green>INFO:    test_bdf.py:797              starting fem2
    </text>



.. raw:: html

    <text style=color:orange>WARNING: test_bdf.py:819              PARAM,POST,0 is not supported by the OP2 reader
    </text>



.. raw:: html

    <text style=color:red>ERROR:   test_bdf.py:1106             An AEROS card is required for STATIC AERO - SOL 144; AEROS=None
    </text>



.. raw:: html

    <text style=color:red>ERROR:   test_bdf.py:1110             An CAEROx card is required for STATIC AERO - SOL 144
    </text>



.. raw:: html

    <text style=color:red>ERROR:   test_bdf.py:1114             An SPLINEx card is required for STATIC AERO - SOL 144
    </text>


::


    ---------------------------------------------------------------------------

    RuntimeError                              Traceback (most recent call last)

    <ipython-input-3-746e54da5a3e> in <module>
          1 from pyNastran.bdf.test.test_bdf import run_bdf as test_bdf
          2 model.write_bdf('junk.bdf')
    ----> 3 test_bdf('.', 'junk.bdf')
    

    c:\nasa\m4\formats\git\pynastran\pyNastran\bdf\test\test_bdf.py in run_bdf(folder, bdf_filename, debug, xref, check, punch, mesh_form, is_folder, print_stats, encoding, sum_load, size, is_double, hdf5, stop, nastran, post, dynamic_vars, quiet, dumplines, dictsort, run_extract_bodies, save_file_structure, nerrors, dev, crash_cards, safe_xref, pickle_obj, stop_on_failure, log)
        317         pickle_obj=pickle_obj,
        318         stop_on_failure=stop_on_failure,
    --> 319         log=log,
        320     )
        321     return fem1, fem2, diff_cards
    

    c:\nasa\m4\formats\git\pynastran\pyNastran\bdf\test\test_bdf.py in run_and_compare_fems(bdf_model, out_model, debug, xref, check, punch, mesh_form, print_stats, encoding, sum_load, size, is_double, save_file_structure, stop, nastran, post, hdf5, dynamic_vars, quiet, dumplines, dictsort, nerrors, dev, crash_cards, safe_xref, run_extract_bodies, pickle_obj, stop_on_failure, log)
        384                         encoding=encoding, debug=debug, quiet=quiet,
        385                         ierror=ierror, nerrors=nerrors,
    --> 386                         stop_on_failure=stop_on_failure, log=log)
        387 
        388         diff_cards = compare(fem1, fem2, xref=xref, check=check,
    

    c:\nasa\m4\formats\git\pynastran\pyNastran\bdf\test\test_bdf.py in run_fem2(bdf_model, out_model, xref, punch, sum_load, size, is_double, mesh_form, safe_xref, encoding, debug, quiet, stop_on_failure, ierror, nerrors, log)
        834                 fem2, p0, sol_base, subcase_keys, subcases, sol_200_map,
        835                 ierror=ierror, nerrors=nerrors,
    --> 836                 stop_on_failure=stop_on_failure)
        837 
        838     if mesh_form is not None:
    

    c:\nasa\m4\formats\git\pynastran\pyNastran\bdf\test\test_bdf.py in validate_case_control(fem2, p0, sol_base, subcase_keys, subcases, unused_sol_200_map, stop_on_failure, ierror, nerrors)
        872         ierror = check_case(
        873             sol_base, subcase, fem2, p0, isubcase, subcases,
    --> 874             ierror=ierror, nerrors=nerrors, stop_on_failure=stop_on_failure)
        875     return ierror
        876 
    

    c:\nasa\m4\formats\git\pynastran\pyNastran\bdf\test\test_bdf.py in check_case(sol, subcase, fem2, p0, isubcase, subcases, ierror, nerrors, stop_on_failure)
       1079 
       1080     elif sol == 200:
    -> 1081         _check_case_sol_200(sol, subcase, fem2, p0, isubcase, subcases, log)
       1082     elif sol in [114, 115, 116, 118]:
       1083         # cyclic statics, modes, buckling, frequency
    

    c:\nasa\m4\formats\git\pynastran\pyNastran\bdf\test\test_bdf.py in _check_case_sol_200(sol, subcase, fem2, p0, isubcase, subcases, log)
       1270     elif analysis in ['SAERO', 'DIVERG', 'DIVERGE']:
       1271         solution = 144
    -> 1272         check_case(solution, subcase, fem2, p0, isubcase, subcases)
       1273     elif analysis == 'FLUTTER':
       1274         solution = 145
    

    c:\nasa\m4\formats\git\pynastran\pyNastran\bdf\test\test_bdf.py in check_case(sol, subcase, fem2, p0, isubcase, subcases, ierror, nerrors, stop_on_failure)
       1092         subcase, fem2, p0, isubcase, sol,
       1093         ierror=ierror, nerrors=nerrors,
    -> 1094         stop_on_failure=stop_on_failure)
       1095     return ierror
       1096 
    

    c:\nasa\m4\formats\git\pynastran\pyNastran\bdf\test\test_bdf.py in _check_case_parameters(subcase, fem2, p0, isubcase, sol, ierror, nerrors, stop_on_failure)
       1339                 'trims=%s\n'
       1340                 'subcase:\n%s' % (trim_id, str(fem2.trims), str(subcase)))
    -> 1341             raise RuntimeError(msg)
       1342         trim = fem2.trims[trim_id]
       1343 
    

    RuntimeError: TRIM = 42
    trims={}
    subcase:
    SUBCASE 1
        ANALYSIS = SAERO
        DESOBJ = 102
        DESSUB = 10000
        LOAD = 1
        METHOD = 42
        SPC = 1
        TRIM = 42
    

