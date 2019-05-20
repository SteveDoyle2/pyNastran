
BDF Introduction
================

The Jupyter notebook for this demo can be found in: -
docs/quick_start/demo/bdf_demo.ipynb -
https://github.com/SteveDoyle2/pyNastran/tree/master/docs/quick_start/demo/bdf_demo.ipynb

Import pyNastran

.. code:: ipython3

    import os
    import pyNastran
    print (pyNastran.__file__)
    print (pyNastran.__version__)
    pkg_path = pyNastran.__path__[0]
    
    from pyNastran.bdf.bdf import BDF, read_bdf
    from pyNastran.utils import object_attributes, object_methods
    
    print("pkg_path = %s" % pkg_path)

Let’s load the iSat model into the pyNastranGUI
-----------------------------------------------

it’s a .dat file, so instead of:

::

   >>> pyNastranGUI -i bdf_filename

we need to include the format:

::

   >>> pyNastranGUI -f nastran -i bdf_filename

Alternatively, we could load the model and the results, but in this demo
we’re just showing off the geometry. To do that instead:

::

   >>> pyNastranGUI -f nastran -i bdf_filename -o op2_filename

.. code:: ipython3

    bdf_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_Rgd.dat'))
    #bdf_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Dploy_Sm.dat'))
    print(bdf_filename)
    
    # look at the model
    !pyNastranGUI -f nastran -i {bdf_filename} > junk.out

Loading a BDF
-------------

There are two ways to load a BDF; the long way or the short way.

The short way instantiates the **``BDF``** class and the short way uses
the **``read_bdf``** function. As this demo was written for the Jupyter
Notebook, we’ll use **``read_bdf``** and then mention the other method.
The class-based method allows finer control over things like: - what
cards should be loaded - OpenMDAO dynamic syntax support

The class-based method
~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    print(bdf_filename)
    
    # create the BDF object
    bdf = BDF()
    
    # read the file from the GUI
    # don't cross-reference
    bdf.read_bdf(bdf_filename, xref=False)


.. parsed-literal::

    c:\nasa\m4\formats\git\pynastran_1.2\models\iSat\ISat_Launch_Sm_Rgd.dat
    


.. raw:: html

    <text style=color:blue>DEBUG:   bdf.py:1107                  ---starting BDF.read_bdf of c:\nasa\m4\formats\git\pynastran_1.2\models\iSat\ISat_Launch_Sm_Rgd.dat---
    </text>



.. raw:: html

    <text style=color:blue>DEBUG:   pybdf.py:542                 opening 'c:\\nasa\\m4\\formats\\git\\pynastran_1.2\\models\\iSat\\ISat_Launch_Sm_Rgd.dat'
    </text>



.. raw:: html

    <text style=color:blue>DEBUG:   bdf.py:1155                  ---finished BDF.read_bdf of c:\nasa\m4\formats\git\pynastran_1.2\models\iSat\ISat_Launch_Sm_Rgd.dat---
    </text>


The function-based method
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    bdf = read_bdf(bdf_filename, xref=False)



.. raw:: html

    <text style=color:blue>DEBUG:   bdf.py:1107                  ---starting BDF.read_bdf of c:\nasa\m4\formats\git\pynastran_1.2\models\iSat\ISat_Launch_Sm_Rgd.dat---
    </text>



.. raw:: html

    <text style=color:blue>DEBUG:   pybdf.py:542                 opening 'c:\\nasa\\m4\\formats\\git\\pynastran_1.2\\models\\iSat\\ISat_Launch_Sm_Rgd.dat'
    </text>



.. raw:: html

    <text style=color:blue>DEBUG:   bdf.py:1155                  ---finished BDF.read_bdf of c:\nasa\m4\formats\git\pynastran_1.2\models\iSat\ISat_Launch_Sm_Rgd.dat---
    </text>


For simplicity of using the demo, we’ll again use the ``read_bdf``
method

.. code:: ipython3

    #bdf_filename = r'D:\work\pynastran_0.8.0_py27\models\iSat\ISat_Launch_Sm_Rgd.dat'
    bdf_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_Rgd.dat'))
    
    # read the file as a path
    bdf_xref = read_bdf(bdf_filename, xref=True)



.. raw:: html

    <text style=color:blue>DEBUG:   bdf.py:1107                  ---starting BDF.read_bdf of c:\nasa\m4\formats\git\pynastran_1.2\models\iSat\ISat_Launch_Sm_Rgd.dat---
    </text>



.. raw:: html

    <text style=color:blue>DEBUG:   pybdf.py:542                 opening 'c:\\nasa\\m4\\formats\\git\\pynastran_1.2\\models\\iSat\\ISat_Launch_Sm_Rgd.dat'
    </text>



.. raw:: html

    <text style=color:blue>DEBUG:   cross_reference.py:159       Cross Referencing...
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2360                PSHELL pid=1 midsurface: z1=0.400000006 z2=-0.400000006 t=0.035999998 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2360                PSHELL pid=2 midsurface: z1=0.400000006 z2=-0.400000006 t=0.054000005 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2360                PSHELL pid=3 midsurface: z1=0.400000006 z2=-0.400000006 t=0.017999999 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2360                PSHELL pid=7 midsurface: z1=0.418000013 z2=-0.418000013 t=0.035999998 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2360                PSHELL pid=34 midsurface: z1=0.194000006 z2=-0.194000006 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2360                PSHELL pid=38 midsurface: z1=0.284000009 z2=-0.284000009 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2360                PSHELL pid=46 midsurface: z1=0.199000001 z2=-0.199000001 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2360                PSHELL pid=37 midsurface: z1=0.308999985 z2=-0.308999985 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:blue>DEBUG:   bdf.py:1155                  ---finished BDF.read_bdf of c:\nasa\m4\formats\git\pynastran_1.2\models\iSat\ISat_Launch_Sm_Rgd.dat---
    </text>


Interrogating the BDF object
----------------------------

IDE’s like WingIDE, PyCharm, Spyder and “Python Tools for Visual Studio”
make it very easy to program with their object introspection ability.
Unfortunately, because pyNastran has so many functions, it can be
difficult to learn the code.

**Some handy object introspection methods were created that will work on
all pyNastran objects and even non-pyNastran objects**. By convention,
private data members/functions start with an underscore \_, and public
ones do not.

We can use the generic object attributes/methods functions

.. code:: ipython3

    print(object_attributes(bdf))
    print(object_methods(bdf))


.. parsed-literal::

    ['MATS1', 'MATS3', 'MATS8', 'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9', 'active_filename', 'active_filenames', 'aecomps', 'aefacts', 'aelinks', 'aelists', 'aeparams', 'aero', 'aeros', 'aestats', 'aesurf', 'aesurfs', 'ao_element_flags', 'asets', 'axic', 'axif', 'baror', 'bconp', 'bcrparas', 'bcs', 'bctadds', 'bctparas', 'bctsets', 'bdf_filename', 'beamor', 'blseg', 'bsets', 'bsurf', 'bsurfs', 'cMethods', 'caero_ids', 'caeros', 'card_count', 'cards_to_read', 'case_control_deck', 'case_control_lines', 'convection_properties', 'coord_ids', 'coords', 'creep_materials', 'csets', 'csschds', 'csuper', 'csupext', 'dareas', 'dconadds', 'dconstrs', 'ddvals', 'debug', 'delays', 'dequations', 'desvars', 'divergs', 'dlinks', 'dload_entries', 'dloads', 'dmigs', 'dmijis', 'dmijs', 'dmiks', 'dmis', 'doptprm', 'dphases', 'dresps', 'dscreen', 'dtable', 'dti', 'dumplines', 'dvcrels', 'dvgrids', 'dvmrels', 'dvprels', 'echo', 'element_ids', 'elements', 'epoints', 'executive_control_lines', 'flfacts', 'flutters', 'force_echo_off', 'frequencies', 'grdset', 'gridb', 'gusts', 'hyperelastic_materials', 'include_dir', 'include_filenames', 'initial_superelement_models', 'is_bdf_vectorized', 'is_long_ids', 'is_msc', 'is_nx', 'is_superelements', 'is_zona', 'load_combinations', 'loads', 'log', 'masses', 'material_ids', 'materials', 'methods', 'mkaeros', 'monitor_points', 'mpcadds', 'mpcs', 'nastran_format', 'ncaeros', 'ncoords', 'nelements', 'nid_map', 'nlparms', 'nlpcis', 'nmaterials', 'nnodes', 'node_ids', 'nodes', 'normals', 'npoints', 'nproperties', 'nsmadds', 'nsms', 'nxstrats', 'omits', 'paeros', 'params', 'pbusht', 'pdampt', 'pelast', 'phbdys', 'plotels', 'point_ids', 'points', 'properties', 'properties_mass', 'property_ids', 'punch', 'qsets', 'radcavs', 'radmtx', 'radset', 'random_tables', 'read_includes', 'reject_cards', 'reject_count', 'reject_lines', 'rigid_elements', 'ringaxs', 'ringfl', 'rotors', 'rsolmap_to_str', 'save_file_structure', 'se_bsets', 'se_csets', 'se_qsets', 'se_sets', 'se_suport', 'se_usets', 'sebndry', 'sebulk', 'seconct', 'seelt', 'seexcld', 'selabel', 'seload', 'seloc', 'sempln', 'senqset', 'seqgp', 'setree', 'sets', 'sol', 'sol_iline', 'sol_method', 'spcadds', 'spcoffs', 'spcs', 'special_cards', 'splines', 'spoints', 'subcases', 'superelement_models', 'suport', 'suport1', 'system_command_lines', 'tables', 'tables_d', 'tables_m', 'tables_sdamping', 'tempds', 'thermal_materials', 'tics', 'transfer_functions', 'trims', 'tstepnls', 'tsteps', 'type_slot_str', 'units', 'usets', 'values_to_skip', 'view3ds', 'views', 'wtmass', 'zona']
    ['AEFact', 'AELIST', 'AELink', 'AEList', 'AEParam', 'AEStat', 'AESurf', 'Acsid', 'Aero', 'Aeros', 'CAero', 'CMethod', 'Coord', 'DAREA', 'DConstr', 'DDVal', 'DELAY', 'DEQATN', 'DLoad', 'DMIG', 'DPHASE', 'DResp', 'DVcrel', 'DVmrel', 'DVprel', 'Desvar', 'Element', 'Elements', 'EmptyNode', 'EmptyNodes', 'FLFACT', 'Flutter', 'Gust', 'HyperelasticMaterial', 'Load', 'MPC', 'Mass', 'Material', 'Materials', 'Method', 'NLParm', 'NSM', 'Node', 'Nodes', 'PAero', 'Phbdy', 'Point', 'Points', 'Properties', 'Property', 'PropertyMass', 'RandomTable', 'RigidElement', 'SET1', 'SPC', 'Set', 'Spline', 'StructuralMaterial', 'Table', 'TableD', 'TableM', 'ThermalMaterial', 'add_accel', 'add_accel1', 'add_acsrce', 'add_aecomp', 'add_aecompl', 'add_aefact', 'add_aelink', 'add_aelist', 'add_aeparm', 'add_aero', 'add_aeros', 'add_aestat', 'add_aesurf', 'add_aesurfs', 'add_aset', 'add_aset1', 'add_axic', 'add_bcrpara', 'add_bctadd', 'add_bctpara', 'add_bctset', 'add_bset', 'add_bset1', 'add_bsurf', 'add_bsurfs', 'add_caero1', 'add_caero2', 'add_caero3', 'add_caero4', 'add_caero5', 'add_card', 'add_card_fields', 'add_card_ifile', 'add_card_lines', 'add_cbar', 'add_cbarao', 'add_cbeam', 'add_cbeam3', 'add_cbend', 'add_cbush', 'add_cbush1d', 'add_cbush2d', 'add_cconeax', 'add_cdamp1', 'add_cdamp2', 'add_cdamp3', 'add_cdamp4', 'add_cdamp5', 'add_celas1', 'add_celas2', 'add_celas3', 'add_celas4', 'add_cfast', 'add_cfluid2', 'add_cfluid3', 'add_cfluid4', 'add_cgap', 'add_cgen', 'add_chbdye', 'add_chbdyg', 'add_chbdyp', 'add_chexa', 'add_cihex1', 'add_cihex2', 'add_cmass1', 'add_cmass2', 'add_cmass3', 'add_cmass4', 'add_cmfree', 'add_conm1', 'add_conm2', 'add_conrod', 'add_conv', 'add_convm', 'add_cord1c', 'add_cord1r', 'add_cord1s', 'add_cord2c', 'add_cord2r', 'add_cord2s', 'add_cpenta', 'add_cplstn3', 'add_cplstn4', 'add_cplstn6', 'add_cplstn8', 'add_cplsts3', 'add_cpyram', 'add_cquad', 'add_cquad4', 'add_cquad8', 'add_cquadr', 'add_cquadx', 'add_cquadx4', 'add_cquadx8', 'add_crac2d', 'add_crac3d', 'add_creep', 'add_crod', 'add_cset', 'add_cset1', 'add_cshear', 'add_csschd', 'add_csuper', 'add_csupext', 'add_ctetra', 'add_ctrax3', 'add_ctrax6', 'add_ctria3', 'add_ctria6', 'add_ctriar', 'add_ctriax', 'add_ctriax6', 'add_ctube', 'add_cvisc', 'add_darea', 'add_dconadd', 'add_dconstr', 'add_ddval', 'add_deform', 'add_delay', 'add_deqatn', 'add_desvar', 'add_diverg', 'add_dlink', 'add_dload', 'add_dmi', 'add_dmig', 'add_dmig_uaccel', 'add_dmij', 'add_dmiji', 'add_dmik', 'add_doptprm', 'add_dphase', 'add_dresp1', 'add_dresp2', 'add_dresp3', 'add_dscreen', 'add_dtable', 'add_dti', 'add_dvcrel1', 'add_dvcrel2', 'add_dvgrid', 'add_dvmrel1', 'add_dvmrel2', 'add_dvprel1', 'add_dvprel2', 'add_eigb', 'add_eigc', 'add_eigp', 'add_eigr', 'add_eigrl', 'add_epoint', 'add_flfact', 'add_flutter', 'add_force', 'add_force1', 'add_force2', 'add_freq', 'add_freq1', 'add_freq2', 'add_freq3', 'add_freq4', 'add_freq5', 'add_genel_flexibility', 'add_genel_stiffness', 'add_gmload', 'add_gmspc', 'add_grav', 'add_grdset', 'add_grid', 'add_gust', 'add_load', 'add_loadcyn', 'add_lseq', 'add_mat1', 'add_mat10', 'add_mat11', 'add_mat2', 'add_mat3', 'add_mat3d', 'add_mat4', 'add_mat5', 'add_mat8', 'add_mat9', 'add_matg', 'add_mathe', 'add_mathp', 'add_mats1', 'add_matt1', 'add_matt2', 'add_matt3', 'add_matt4', 'add_matt5', 'add_matt8', 'add_matt9', 'add_mkaero1', 'add_mkaero2', 'add_moment', 'add_moment1', 'add_moment2', 'add_monpnt1', 'add_monpnt2', 'add_monpnt3', 'add_mpc', 'add_mpcadd', 'add_nlparm', 'add_nlpci', 'add_nsm', 'add_nsm1', 'add_nsmadd', 'add_nsml', 'add_nsml1', 'add_nxstrat', 'add_omit1', 'add_paero1', 'add_paero2', 'add_paero3', 'add_paero4', 'add_paero5', 'add_param', 'add_pbar', 'add_pbarl', 'add_pbcomp', 'add_pbeam', 'add_pbeam3', 'add_pbeaml', 'add_pbend', 'add_pbmsect', 'add_pbrsect', 'add_pbush', 'add_pbush1d', 'add_pbusht', 'add_pcomp', 'add_pcompg', 'add_pcomps', 'add_pconeax', 'add_pconv', 'add_pconvm', 'add_pdamp', 'add_pdamp5', 'add_pdampt', 'add_pelas', 'add_pelast', 'add_pfast', 'add_pgap', 'add_phbdy', 'add_pihex', 'add_pload', 'add_pload1', 'add_pload2', 'add_pload4', 'add_ploadx1', 'add_plotel', 'add_plplane', 'add_plsolid', 'add_pmass', 'add_point', 'add_pointax', 'add_pplane', 'add_prac2d', 'add_prac3d', 'add_presax', 'add_prod', 'add_pshear', 'add_pshell', 'add_psolid', 'add_ptube', 'add_pvisc', 'add_qbdy1', 'add_qbdy2', 'add_qbdy3', 'add_qhbdy', 'add_qset', 'add_qset1', 'add_qvect', 'add_qvol', 'add_radbc', 'add_radm', 'add_randps', 'add_randt1', 'add_rbar', 'add_rbar1', 'add_rbe1', 'add_rbe2', 'add_rbe3', 'add_rforce', 'add_rforce1', 'add_rgyro', 'add_ringax', 'add_rload1', 'add_rload2', 'add_rotord', 'add_rotorg', 'add_rrod', 'add_rspint', 'add_rspline', 'add_rsscon', 'add_sebndry', 'add_sebset', 'add_sebset1', 'add_sebulk', 'add_seconct', 'add_secset', 'add_secset1', 'add_seelt', 'add_seexcld', 'add_selabel', 'add_seload', 'add_seloc', 'add_sempln', 'add_senqset', 'add_seqgp', 'add_seqset', 'add_seqset1', 'add_seset', 'add_sesup', 'add_set1', 'add_set3', 'add_setree', 'add_sload', 'add_spc', 'add_spc1', 'add_spcadd', 'add_spcax', 'add_spcd', 'add_spline1', 'add_spline2', 'add_spline3', 'add_spline4', 'add_spline5', 'add_spoint', 'add_suport', 'add_suport1', 'add_tabdmp1', 'add_tabled1', 'add_tabled2', 'add_tabled3', 'add_tabled4', 'add_tablem1', 'add_tablem2', 'add_tablem3', 'add_tablem4', 'add_tables1', 'add_tablest', 'add_tabrnd1', 'add_tabrndg', 'add_temp', 'add_tempax', 'add_tempd', 'add_tf', 'add_tic', 'add_tload1', 'add_tload2', 'add_trim', 'add_tstep', 'add_tstep1', 'add_tstepnl', 'add_uset', 'add_uset1', 'clear_attributes', 'create_card_object', 'cross_reference', 'deprecated', 'disable_cards', 'export_hdf5_file', 'export_hdf5_filename', 'geom_check', 'get_MPCx_node_ids', 'get_MPCx_node_ids_c1', 'get_SPCx_node_ids', 'get_SPCx_node_ids_c1', 'get_area_breakdown', 'get_bdf_cards', 'get_bdf_cards_dict', 'get_bdf_stats', 'get_card_ids_by_card_types', 'get_cards_by_card_types', 'get_custom_types', 'get_dependent_nid_to_components', 'get_displacement_index', 'get_displacement_index_xyz_cp_cd', 'get_dload_entries', 'get_element_faces', 'get_element_ids_dict_with_pids', 'get_element_ids_list_with_pids', 'get_element_nodes_by_element_type', 'get_elements_nodes_by_property_type', 'get_elements_properties_nodes_by_element_type', 'get_encoding', 'get_h5attrs', 'get_length_breakdown', 'get_load_arrays', 'get_mass_breakdown', 'get_material_id_to_property_ids_map', 'get_material_ids', 'get_mklist', 'get_mpcs', 'get_nid_map', 'get_node_id_to_element_ids_map', 'get_node_id_to_elements_map', 'get_node_ids_with_elements', 'get_pid_to_node_ids_and_elements_array', 'get_point_grids', 'get_pressure_array', 'get_property_id_to_element_ids_map', 'get_reduced_dloads', 'get_reduced_loads', 'get_reduced_mpcs', 'get_reduced_nsms', 'get_reduced_spcs', 'get_rigid_elements_with_node_ids', 'get_rslot_map', 'get_spcs', 'get_structural_material_ids', 'get_thermal_material_ids', 'get_volume_breakdown', 'get_xyz_in_coord', 'get_xyz_in_coord_array', 'get_xyz_in_coord_no_xref', 'include_zip', 'increase_card_count', 'is_reject', 'load', 'load_hdf5_file', 'load_hdf5_filename', 'mass_properties', 'mass_properties_no_xref', 'mass_properties_nsm', 'object_attributes', 'object_methods', 'pop_parse_errors', 'pop_xref_errors', 'read_bdf', 'reject_card_lines', 'replace_cards', 'reset_errors', 'reset_rslot_map', 'safe_acsid', 'safe_aefact', 'safe_aelist', 'safe_caero', 'safe_coord', 'safe_cross_reference', 'safe_element', 'safe_elements', 'safe_empty_nodes', 'safe_get_elements', 'safe_get_nodes', 'safe_get_points', 'safe_material', 'safe_paero', 'safe_property', 'safe_property_mass', 'safe_tabled', 'safe_tableh', 'save', 'saves', 'set_as_msc', 'set_as_nx', 'set_as_zona', 'set_cards', 'set_dynamic_syntax', 'set_error_storage', 'set_param', 'sum_forces_moments', 'sum_forces_moments_elements', 'superelement_nodes', 'transform_xyzcp_to_xyz_cid', 'uncross_reference', 'update_card', 'update_model_by_desvars', 'update_solution', 'validate', 'write_bdf', 'write_bdfs', 'write_skin_solid_faces']
    

Let’s clean that up a bit
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    print("attributes = [%s]\n" % ', '.join(bdf.object_attributes()))
    print("methods = [%s]\n" % ', '.join(bdf.object_methods()))


.. parsed-literal::

    attributes = [MATS1, MATS3, MATS8, MATT1, MATT2, MATT3, MATT4, MATT5, MATT8, MATT9, active_filename, active_filenames, aecomps, aefacts, aelinks, aelists, aeparams, aero, aeros, aestats, aesurf, aesurfs, ao_element_flags, asets, axic, axif, baror, bconp, bcrparas, bcs, bctadds, bctparas, bctsets, bdf_filename, beamor, blseg, bsets, bsurf, bsurfs, cMethods, caeros, card_count, cards_to_read, case_control_deck, case_control_lines, convection_properties, coords, creep_materials, csets, csschds, csuper, csupext, dareas, dconadds, dconstrs, ddvals, debug, delays, dequations, desvars, divergs, dlinks, dload_entries, dloads, dmigs, dmijis, dmijs, dmiks, dmis, doptprm, dphases, dresps, dscreen, dtable, dti, dumplines, dvcrels, dvgrids, dvmrels, dvprels, echo, elements, epoints, executive_control_lines, flfacts, flutters, force_echo_off, frequencies, grdset, gridb, gusts, hyperelastic_materials, include_dir, include_filenames, initial_superelement_models, is_bdf_vectorized, is_msc, is_nx, is_superelements, is_zona, load_combinations, loads, masses, materials, methods, mkaeros, monitor_points, mpcadds, mpcs, nastran_format, nid_map, nlparms, nlpcis, nodes, normals, npoints, nsmadds, nsms, nxstrats, omits, paeros, params, pbusht, pdampt, pelast, phbdys, plotels, points, properties, properties_mass, punch, qsets, radcavs, radmtx, radset, random_tables, read_includes, reject_cards, reject_count, reject_lines, rigid_elements, ringaxs, ringfl, rotors, rsolmap_to_str, save_file_structure, se_bsets, se_csets, se_qsets, se_sets, se_suport, se_usets, sebndry, sebulk, seconct, seelt, seexcld, selabel, seload, seloc, sempln, senqset, seqgp, setree, sets, sol, sol_iline, sol_method, spcadds, spcoffs, spcs, special_cards, splines, spoints, superelement_models, suport, suport1, system_command_lines, tables, tables_d, tables_m, tables_sdamping, tempds, thermal_materials, tics, transfer_functions, trims, tstepnls, tsteps, type_slot_str, units, usets, values_to_skip, view3ds, views, wtmass, zona]
    
    methods = [AEFact, AELIST, AELink, AEList, AEParam, AEStat, AESurf, Acsid, Aero, Aeros, CAero, CMethod, Coord, DAREA, DConstr, DDVal, DELAY, DEQATN, DLoad, DMIG, DPHASE, DResp, DVcrel, DVmrel, DVprel, Desvar, Element, Elements, EmptyNode, EmptyNodes, FLFACT, Flutter, Gust, HyperelasticMaterial, Load, MPC, Mass, Material, Materials, Method, NLParm, NSM, Node, Nodes, PAero, Phbdy, Point, Points, Properties, Property, PropertyMass, RandomTable, RigidElement, SET1, SPC, Set, Spline, StructuralMaterial, Table, TableD, TableM, ThermalMaterial, add_accel, add_accel1, add_acsrce, add_aecomp, add_aecompl, add_aefact, add_aelink, add_aelist, add_aeparm, add_aero, add_aeros, add_aestat, add_aesurf, add_aesurfs, add_aset, add_aset1, add_axic, add_bcrpara, add_bctadd, add_bctpara, add_bctset, add_bset, add_bset1, add_bsurf, add_bsurfs, add_caero1, add_caero2, add_caero3, add_caero4, add_caero5, add_card, add_card_fields, add_card_ifile, add_card_lines, add_cbar, add_cbarao, add_cbeam, add_cbeam3, add_cbend, add_cbush, add_cbush1d, add_cbush2d, add_cconeax, add_cdamp1, add_cdamp2, add_cdamp3, add_cdamp4, add_cdamp5, add_celas1, add_celas2, add_celas3, add_celas4, add_cfast, add_cfluid2, add_cfluid3, add_cfluid4, add_cgap, add_cgen, add_chbdye, add_chbdyg, add_chbdyp, add_chexa, add_cihex1, add_cihex2, add_cmass1, add_cmass2, add_cmass3, add_cmass4, add_cmfree, add_conm1, add_conm2, add_conrod, add_conv, add_convm, add_cord1c, add_cord1r, add_cord1s, add_cord2c, add_cord2r, add_cord2s, add_cpenta, add_cplstn3, add_cplstn4, add_cplstn6, add_cplstn8, add_cplsts3, add_cpyram, add_cquad, add_cquad4, add_cquad8, add_cquadr, add_cquadx, add_cquadx4, add_cquadx8, add_crac2d, add_crac3d, add_creep, add_crod, add_cset, add_cset1, add_cshear, add_csschd, add_csuper, add_csupext, add_ctetra, add_ctrax3, add_ctrax6, add_ctria3, add_ctria6, add_ctriar, add_ctriax, add_ctriax6, add_ctube, add_cvisc, add_darea, add_dconadd, add_dconstr, add_ddval, add_deform, add_delay, add_deqatn, add_desvar, add_diverg, add_dlink, add_dload, add_dmi, add_dmig, add_dmig_uaccel, add_dmij, add_dmiji, add_dmik, add_doptprm, add_dphase, add_dresp1, add_dresp2, add_dresp3, add_dscreen, add_dtable, add_dti, add_dvcrel1, add_dvcrel2, add_dvgrid, add_dvmrel1, add_dvmrel2, add_dvprel1, add_dvprel2, add_eigb, add_eigc, add_eigp, add_eigr, add_eigrl, add_epoint, add_flfact, add_flutter, add_force, add_force1, add_force2, add_freq, add_freq1, add_freq2, add_freq3, add_freq4, add_freq5, add_genel_flexibility, add_genel_stiffness, add_gmload, add_gmspc, add_grav, add_grdset, add_grid, add_gust, add_load, add_loadcyn, add_lseq, add_mat1, add_mat10, add_mat11, add_mat2, add_mat3, add_mat3d, add_mat4, add_mat5, add_mat8, add_mat9, add_matg, add_mathe, add_mathp, add_mats1, add_matt1, add_matt2, add_matt3, add_matt4, add_matt5, add_matt8, add_matt9, add_mkaero1, add_mkaero2, add_moment, add_moment1, add_moment2, add_monpnt1, add_monpnt2, add_monpnt3, add_mpc, add_mpcadd, add_nlparm, add_nlpci, add_nsm, add_nsm1, add_nsmadd, add_nsml, add_nsml1, add_nxstrat, add_omit1, add_paero1, add_paero2, add_paero3, add_paero4, add_paero5, add_param, add_pbar, add_pbarl, add_pbcomp, add_pbeam, add_pbeam3, add_pbeaml, add_pbend, add_pbmsect, add_pbrsect, add_pbush, add_pbush1d, add_pbusht, add_pcomp, add_pcompg, add_pcomps, add_pconeax, add_pconv, add_pconvm, add_pdamp, add_pdamp5, add_pdampt, add_pelas, add_pelast, add_pfast, add_pgap, add_phbdy, add_pihex, add_pload, add_pload1, add_pload2, add_pload4, add_ploadx1, add_plotel, add_plplane, add_plsolid, add_pmass, add_point, add_pointax, add_pplane, add_prac2d, add_prac3d, add_presax, add_prod, add_pshear, add_pshell, add_psolid, add_ptube, add_pvisc, add_qbdy1, add_qbdy2, add_qbdy3, add_qhbdy, add_qset, add_qset1, add_qvect, add_qvol, add_radbc, add_radm, add_randps, add_randt1, add_rbar, add_rbar1, add_rbe1, add_rbe2, add_rbe3, add_rforce, add_rforce1, add_rgyro, add_ringax, add_rload1, add_rload2, add_rotord, add_rotorg, add_rrod, add_rspint, add_rspline, add_rsscon, add_sebndry, add_sebset, add_sebset1, add_sebulk, add_seconct, add_secset, add_secset1, add_seelt, add_seexcld, add_selabel, add_seload, add_seloc, add_sempln, add_senqset, add_seqgp, add_seqset, add_seqset1, add_seset, add_sesup, add_set1, add_set3, add_setree, add_sload, add_spc, add_spc1, add_spcadd, add_spcax, add_spcd, add_spline1, add_spline2, add_spline3, add_spline4, add_spline5, add_spoint, add_suport, add_suport1, add_tabdmp1, add_tabled1, add_tabled2, add_tabled3, add_tabled4, add_tablem1, add_tablem2, add_tablem3, add_tablem4, add_tables1, add_tablest, add_tabrnd1, add_tabrndg, add_temp, add_tempax, add_tempd, add_tf, add_tic, add_tload1, add_tload2, add_trim, add_tstep, add_tstep1, add_tstepnl, add_uset, add_uset1, clear_attributes, create_card_object, cross_reference, deprecated, disable_cards, export_hdf5_file, export_hdf5_filename, geom_check, get_MPCx_node_ids, get_MPCx_node_ids_c1, get_SPCx_node_ids, get_SPCx_node_ids_c1, get_area_breakdown, get_bdf_cards, get_bdf_cards_dict, get_bdf_stats, get_card_ids_by_card_types, get_cards_by_card_types, get_custom_types, get_dependent_nid_to_components, get_displacement_index, get_displacement_index_xyz_cp_cd, get_dload_entries, get_element_faces, get_element_ids_dict_with_pids, get_element_ids_list_with_pids, get_element_nodes_by_element_type, get_elements_nodes_by_property_type, get_elements_properties_nodes_by_element_type, get_encoding, get_h5attrs, get_length_breakdown, get_load_arrays, get_mass_breakdown, get_material_id_to_property_ids_map, get_material_ids, get_mklist, get_mpcs, get_nid_map, get_node_id_to_element_ids_map, get_node_id_to_elements_map, get_node_ids_with_elements, get_pid_to_node_ids_and_elements_array, get_point_grids, get_pressure_array, get_property_id_to_element_ids_map, get_reduced_dloads, get_reduced_loads, get_reduced_mpcs, get_reduced_nsms, get_reduced_spcs, get_rigid_elements_with_node_ids, get_rslot_map, get_spcs, get_structural_material_ids, get_thermal_material_ids, get_volume_breakdown, get_xyz_in_coord, get_xyz_in_coord_array, get_xyz_in_coord_no_xref, include_zip, increase_card_count, is_reject, load, load_hdf5_file, load_hdf5_filename, mass_properties, mass_properties_no_xref, mass_properties_nsm, pop_parse_errors, pop_xref_errors, read_bdf, reject_card_lines, replace_cards, reset_errors, reset_rslot_map, safe_acsid, safe_aefact, safe_aelist, safe_caero, safe_coord, safe_cross_reference, safe_element, safe_elements, safe_empty_nodes, safe_get_elements, safe_get_nodes, safe_get_points, safe_material, safe_paero, safe_property, safe_property_mass, safe_tabled, safe_tableh, save, saves, set_as_msc, set_as_nx, set_as_zona, set_cards, set_dynamic_syntax, set_error_storage, set_param, sum_forces_moments, sum_forces_moments_elements, superelement_nodes, transform_xyzcp_to_xyz_cid, uncross_reference, update_card, update_model_by_desvars, update_solution, validate, write_bdf, write_bdfs, write_skin_solid_faces]
    
    

Some other very handy methods that will be used later by ``test_bdf``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    print(bdf.get_bdf_stats())
    print("card_count = %s\n" % bdf.card_count)
    print("reject_count = %s" % bdf.reject_count)


.. parsed-literal::

    ---BDF Statistics---
    SOL 103
    
    bdf.spcs[1]
      SPC:     1
    
    bdf.params
      PARAM    : 8
    
    bdf.nodes
      GRID     : 5380
    
    bdf.elements
      CBAR     : 827
      CBUSH    : 104
      CHEXA    : 25
      CQUAD4   : 4580
      CTRIA3   : 32
    
    bdf.rigid_elements
      RBE2     : 44
    
    bdf.properties
      PBAR     : 1
      PBARL    : 18
      PBUSH    : 2
      PSHELL   : 8
      PSOLID   : 4
    
    bdf.masses
      CONM2    : 15
    
    bdf.materials
      MAT1     : 14
      MAT8     : 8
    
    bdf.coords
      CORD2R   : 75
    
    bdf.methods
      EIGRL    : 1
    
    bdf.usets
      USET     : 1
    
    
    card_count = {'ENDDATA': 1, 'PARAM': 8, 'SPC': 1, 'USET': 1, 'EIGRL': 1, 'CORD2R': 75, 'GRID': 5380, 'CQUAD4': 4580, 'CBAR': 827, 'CHEXA': 25, 'RBE2': 44, 'CTRIA3': 32, 'CBUSH': 104, 'CONM2': 15, 'MAT1': 14, 'MAT8': 8, 'PSHELL': 8, 'PBARL': 18, 'PSOLID': 4, 'PBAR': 1, 'PBUSH': 2}
    
    reject_count = {}
    

Cross-referencing
-----------------

Cross-referencing a BDF allows improved usability of the **``BDF``**
class. It comes with some negative side effects, but in general is a
very useful thing. It dramatically minimizes the amount of code you need
to write, greatly simplifies future operations, and is highly
recommended.

The major downside is it slows down the code.

Without Cross-Referencing (xref=False)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here the raw values of the the data objects are returned to us

.. code:: ipython3

    cquad = bdf.elements[1]
    print(cquad)
    nid1 = cquad.nodes[0]
    print("nid1 = %s" % nid1)
    n1 = bdf.nodes[nid1]
    cd4 = n1.cd
    c4 = bdf.coords[cd4]
    print("i (xref=False) = %s" % str(c4.i))
    #print object_attributes(c4)


.. parsed-literal::

    $*
    $*  ELEMENT CARDS
    $*
    CQUAD4         1       1       1       2       4       3
    
    nid1 = 1
    i (xref=False) = [1. 0. 0.]
    

Cross-Referenced (xref=True)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we can trace the referenced objects very easily.

A cross-referenced attribute is indicated with the **``*_ref``** suffix:
\* ``cquad4_element.nodes`` : not cross referenced \*
``cquad4_element.nodes_ref`` : cross referenced

.. code:: ipython3

    print("i (xref=True) = %s" % bdf_xref.elements[1].nodes_ref[0].cd_ref.i)


.. parsed-literal::

    i (xref=True) = [1. 0. 0.]
    

So how is this done?

.. code:: ipython3

    cquad.nodes_ref = []
    cquad.nodes_ref.append(n1)
    print(cquad.nodes_ref[0])


.. parsed-literal::

    $*
    $*  GRID CARDS
    $*
    GRID           1       4    -4.5    -7.5    -14.       4
    
    

Let’s show off the GRID card
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    # some Grid methods
    n1 = bdf_xref.nodes[1]
    print(n1)
    
    # the comment
    c1 = bdf_xref.nodes[1].comment
    c2 = bdf_xref.nodes[2].comment
    print("c1=%r" % c1)
    print("c2=%r" % c2)
    
    
    # get the position of a node
    # in the local cooordinate system
    print("xyz = %s" % n1.xyz)
    
    # in the global frame
    print("position = %s" % n1.get_position())
    
    # in an arbitrary frame
    print("wrt5 = %s" % n1.get_position_wrt(bdf, 5))
    print("wrt4 = %s" % n1.get_position_wrt(bdf, 4))
    


.. parsed-literal::

    $*
    $*  GRID CARDS
    $*
    GRID           1       4    -4.5    -7.5    -14.       4
    
    c1='$*\n$*  GRID CARDS\n$*\n'
    c2=''
    xyz = [ -4.5  -7.5 -14. ]
    position = [ -4.5  -7.5 -14. ]
    wrt5 = [  2.12132034  14.         -26.59188309]
    wrt4 = [ -4.5  -7.5 -14. ]
    

Now let’s modify the **``GRID``** card and write it out

.. code:: ipython3

    n1 = bdf_xref.nodes[1]
    n1.xyz[1] = -7.5
    print("repr  = %s" % n1.repr_fields())
    print("raw   = %s" % n1.repr_fields())
    
    # n1.xyz[1] = 100000000000.
    print("repr2 = %s" % n1.repr_fields())
    print(n1)
    print(n1.write_card(size=8))
    print(n1.write_card(size=16, is_double=False))
    print(n1.write_card(size=16, is_double=True))


.. parsed-literal::

    repr  = ['GRID', 1, 4, -4.5, -7.5, -14.0, 4, '', None]
    raw   = ['GRID', 1, 4, -4.5, -7.5, -14.0, 4, '', None]
    repr2 = ['GRID', 1, 4, -4.5, -7.5, -14.0, 4, '', None]
    $*
    $*  GRID CARDS
    $*
    GRID           1       4    -4.5    -7.5    -14.       4
    
    $*
    $*  GRID CARDS
    $*
    GRID           1       4    -4.5    -7.5    -14.       4                
    
    $*
    $*  GRID CARDS
    $*
    GRID*                  1               4            -4.5            -7.5
    *                   -14.               4                                
    
    $*
    $*  GRID CARDS
    $*
    GRID*                  1               4-4.500000000D+00-7.500000000D+00
    *       -1.400000000D+01               4                                
    
    

Calculating the mass of the structure
-------------------------------------

You can also calculate the mass of individual groups

.. code:: ipython3

    mass, cg, I = bdf_xref.mass_properties()
    print("mass = %s\n" % mass)
    
    bdf_xref.mass_properties()
    area_breakdown = bdf_xref.get_area_breakdown(property_ids=None, sum_bar_area=True)
    table_lines = ['%-3s: %s\n' % (k, v) for k, v in sorted(area_breakdown.items())]
    print('area_breakdown:\n%s\n' % ''.join(table_lines))
    
    pids_to_mass, mass_type_to_mass = bdf_xref.get_mass_breakdown(property_ids=None, stop_if_no_mass=True)
    table_lines = ['%-3s: %s\n' % (k, v) for k, v in sorted(pids_to_mass.items())]
    print('mass_breakdown properties:\n%s\n' % ''.join(table_lines))
    print('mass_breakdown masses:\n%s\n' % mass_type_to_mass)
    
    volume_breakdown = bdf_xref.get_volume_breakdown(property_ids=None)
    table_lines = ['%-3s: %s\n' % (k, v) for k, v in sorted(volume_breakdown.items())]
    print('volume_breakdown:\n%s' % ''.join(table_lines))
    


.. parsed-literal::

    mass = 1.7746011578443164
    
    area_breakdown:
    1  : 2807.9999999999873
    2  : 3775.821643408031
    3  : 3126.7012947840108
    4  : 30.79009842432805
    7  : 2815.099327279977
    9  : 14.828303331033284
    10 : 0.47123889803846714
    12 : 0.5364976066492582
    13 : 0.8885208097843912
    14 : 0.7040464383245798
    15 : 0.5826161160904479
    16 : 0.6148182523366601
    19 : 27.937341701307243
    21 : 1.8849560445432187
    22 : 0.13412440166231454
    23 : 0.17601160958114495
    33 : 4.523893217594095
    34 : 726.2061840339735
    35 : 7.488919259999999
    36 : 8.228384167359806
    37 : 3527.999999999994
    38 : 1270.3425529271904
    39 : 2.6317113721498187
    41 : 1.1313811432273508
    42 : 10.990017829889638
    43 : 18.115180916725283
    46 : 741.00586611462
    
    
    mass_breakdown properties:
    1  : 0.027277887741788746
    2  : 0.047992762103254816
    3  : 0.02099754503346675
    4  : 0.012215670441971119
    5  : 0.33015780000000006
    7  : 0.027813499178780537
    8  : 0.08158356674999993
    9  : 0.07764247249451299
    10 : 0.0002359723127060847
    11 : 0.041699540901982614
    12 : 0.0004572896419732434
    13 : 0.003885126000791797
    14 : 0.00035284798842724434
    15 : 0.0036261119992419943
    16 : 0.0
    19 : 0.017748849515790373
    20 : 0.163081991669256
    21 : 0.0036250296551127333
    22 : 0.0
    23 : 0.0
    33 : 0.0013462746328861252
    34 : 0.0035610911559462153
    35 : 0.0
    36 : 0.007196651333033216
    37 : 0.0945659013652082
    38 : 0.007602230793126197
    39 : 0.0024328318514722706
    41 : 0.0007353854868221845
    42 : 0.0088541268554014
    43 : 0.01224145559837195
    46 : 0.003671235343011555
    
    
    mass_breakdown masses:
    {'CONM2': 0.7720000099999998}
    
    volume_breakdown:
    1  : 101.08799438399943
    2  : 203.8943876231405
    3  : 56.280620179411706
    4  : 68.35241635875428
    5  : 2100.0
    7  : 101.34357015187966
    8  : 110.24999999999991
    9  : 48.56269340913401
    10 : 0.9110618695410359
    11 : 466.6577492124528
    12 : 2.5587504307014672
    13 : 1.3327812149482106
    14 : 1.9743503011887258
    15 : 0.5826161159686574
    16 : 1.9983305902008948
    19 : 99.3131533594663
    20 : 2376.0
    21 : 4.7123901113580455
    22 : 0.15084661593460977
    23 : 0.16128631010351485
    33 : 5.197811005397997
    34 : 13.507435023031892
    35 : 120.88885224756376
    36 : 40.26864596920936
    37 : 65.62080000000002
    38 : 23.62837148444574
    39 : 13.612837415073471
    41 : 4.11482736197198
    42 : 49.5430003771425
    43 : 68.496696425457
    46 : 13.78270910973193
    
    

Examples of xref on elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    eid100 = bdf_xref.elements[100]
    print(eid100)
    print("nodes = %s" % eid100.nodes)
    print("--node0--\n%s" % eid100.nodes_ref[0])
    print("--cd--\n%s" % eid100.nodes_ref[0].cd)
    print("cd.cid = %s" % eid100.nodes_ref[0].cd_ref.cid)
    
    print("area = %s" % eid100.Area())
    print("mass = %s" % eid100.Mass())
    print("--pid--\n%s" % eid100.pid)
    print("pid.pid = %s" % eid100.pid_ref.pid)
    print("pid.Pid() = %s" % eid100.Pid())
    
    print(eid100.pid_ref.mid1_ref)
    print("type = %s" % eid100.pid_ref.mid1_ref.type)
    print("nu12 = %s" % eid100.pid_ref.mid1_ref.nu12)
    print("mass = %s" % eid100.Mass())


.. parsed-literal::

    CQUAD4       100       1     149     152     161     160
    
    nodes = [149, 152, 161, 160]
    --node0--
    GRID         149       4      3.     7.5   -16.5       4
    
    --cd--
    4
    cd.cid = 4
    area = 3.75
    mass = 3.642880307396999e-05
    --pid--
    1
    pid.pid = 1
    pid.Pid() = 1
    $*
    $*  I-DEAS Material: 6  name: BUS_CFRP_PW_ORTHO
    $* M46J PW ETW
    MAT8           6   1.7+7   1.7+7     .98 340000. 180000. 180000..0001712
                               71.33
    
    type = MAT8
    nu12 = 0.98
    mass = 3.642880307396999e-05
    

Write the modified deck
-----------------------

Let’s first switch to the desktop to make the file easy to find

.. code:: ipython3

    import getpass
    name = getpass.getuser()
    os.chdir(os.path.join(r'C:\Users', name, 'Desktop'))

.. code:: ipython3

    pwd




.. parsed-literal::

    'C:\\Users\\sdoyle\\Desktop'



There are two ways to write a deck - **``interspersed``** : alternate
properties and elements (similar to how Patran writes decks) -
**``not-interspersed (default)``** : much faster

We can also use 8 or 16 character field width as well as double
precision.

Note that double precision only works for certain cards (e.g. ``GRID``,
``COORD``, ``DMIG``) and not much else.

.. code:: ipython3

    bdf_xref.write_bdf('fem.bdf', interspersed=False, size=8, is_double=False)
    !tail -n 5 "fem.bdf"
    
    bdf_xref.write_bdf('fem.bdf', interspersed=True, size=16, is_double=False)
    !tail "fem.bdf"
    
    bdf_xref.write_bdf('fem.bdf', interspersed=True, size=16, is_double=True)
    !tail "fem.bdf"



.. raw:: html

    <text style=color:blue>DEBUG:   write_mesh.py:157            ---starting BDF.write_bdf of fem.bdf---
    </text>


.. parsed-literal::

    CORD2R        75        1.355-13-2.19-15    -40.1.355-13-2.19-15      0.
                 40.-2.19-15    -40.
    CORD2R        76        1.355-13-2.19-15    -40.1.355-13-2.19-15      0.
                 40.-2.19-15    -40.
    ENDDATA
    


.. raw:: html

    <text style=color:blue>DEBUG:   write_mesh.py:157            ---starting BDF.write_bdf of fem.bdf---
    </text>


.. parsed-literal::

    *
    CORD2R*               75                 1.3549966049-13-2.1854783949-15
    *                   -40. 1.3549966049-13-2.1854783949-15              0.
    *                    40.-2.1854783949-15            -40.
    *
    CORD2R*               76                 1.3549966049-13-2.1854783949-15
    *                   -40. 1.3549966049-13-2.1854783949-15              0.
    *                    40.-2.1854783949-15            -40.
    *
    ENDDATA
    


.. raw:: html

    <text style=color:blue>DEBUG:   write_mesh.py:157            ---starting BDF.write_bdf of fem.bdf---
    </text>


.. parsed-literal::

    *
    CORD2R*               75                1.3549966049D-13-2.185478395D-15
    *       -4.000000000D+011.3549966049D-13-2.185478395D-150.0000000000D+00
    *       4.0000000000D+01-2.185478395D-15-4.000000000D+01
    *
    CORD2R*               76                1.3549966049D-13-2.185478395D-15
    *       -4.000000000D+011.3549966049D-13-2.185478395D-150.0000000000D+00
    *       4.0000000000D+01-2.185478395D-15-4.000000000D+01
    *
    ENDDATA
    

.. code:: ipython3

    bdf_filename




.. parsed-literal::

    'c:\\nasa\\m4\\formats\\git\\pynastran_1.2\\models\\iSat\\ISat_Launch_Sm_Rgd.dat'



pyNastranGUI
------------

.. code:: ipython3

    print(bdf_filename)
    %echo {bdf_filename}
    #!pyNastranGUI -f nastran -i {bdf_filename}
    
    solid_bending_bdf = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending.bdf'))
    solid_bending_op2 = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending.op2'))
    
    !pyNastranGUI -f nastran -i {solid_bending_bdf} -o {solid_bending_op2}  > junk.out
    print("done")


.. parsed-literal::

    c:\nasa\m4\formats\git\pynastran_1.2\models\iSat\ISat_Launch_Sm_Rgd.dat
    c:\nasa\m4\formats\git\pynastran_1.2\models\iSat\ISat_Launch_Sm_Rgd.dat
    done
    

We can also script the GUI!
===========================

.. code:: ipython3

    solid_bending_bdf = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending.bdf'))
    solid_bending_op2 = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending.op2'))
    
    if os.path.exists('wireframe_solid_bending.png'):
        os.remove('wireframe_solid_bending.png')
    
    with open('script.py', 'w') as f:
        f.write('self.on_wireframe()\n')
        picture_filename = os.path.join(os.getcwd(), 'wireframe_solid_bending.png')
        f.write("self.on_take_screenshot(%r)\n" % picture_filename)
        f.write('sys.exit()')
    
    !pwd
    !pyNastranGUI -f nastran -i {solid_bending_bdf} -o {solid_bending_op2} --postscript script.py > junk.out
    
    # display in a popup
    !wireframe_solid_bending.png
    
    from IPython.display import Image
    from IPython.display import display
    assert os.path.exists('wireframe_solid_bending.png')
    
    # display in iPython
    i = Image(filename='wireframe_solid_bending.png')
    display(i)
    print("the picture is visible")


.. parsed-literal::

    /cygdrive/c/Users/sdoyle/Desktop
    


.. image:: bdf_demo_files%5Cbdf_demo_41_1.png


.. parsed-literal::

    the picture is visible
    


