
BDF Demo
========

The iPython notebook for this demo can be found in: -
docs:raw-latex:`\quick`\_start:raw-latex:`\demo`:raw-latex:`\bdf`\_demo.ipynb
-
https://github.com/SteveDoyle2/pyNastran/tree/master/docs/quick\_start/demo/bdf\_demo.ipynb

Import pyNastran

.. code:: python

    import os
    import pyNastran
    print (pyNastran.__file__)
    print (pyNastran.__version__)
    pkg_path = pyNastran.__path__[0]
    
    from pyNastran.bdf.bdf import BDF, read_bdf
    from pyNastran.utils import object_attributes, object_methods
    
    print("pkg_path = %s" % pkg_path)


.. parsed-literal::

    f:\work\pynastran\pynastran\master3\pyNastran\__init__.pyc
    0.8.0+dev.e79582b
    pkg_path = f:\work\pynastran\pynastran\master3\pyNastran
    

Let's load the iSat model into the pyNastranGUI
-----------------------------------------------

it's a .dat file, so instead of:

::

    >>> pyNastranGUI -i bdf_filename

we need to include the format:

::

    >>> pyNastranGUI -f nastran -i bdf_filename

Alternatively, we could load the model and the results, but in this demo
we're just showing off the geometry. To do that instead:

::

    >>> pyNastranGUI -f nastran -i bdf_filename -o op2_filename

.. code:: python

    bdf_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_Rgd.dat'))
    print(bdf_filename)
    
    # look at the model
    !pyNastranGUI -f nastran -i {bdf_filename} > junk.out


.. parsed-literal::

    f:\work\pynastran\pynastran\master3\models\iSat\ISat_Launch_Sm_Rgd.dat
    

.. parsed-literal::

    QObject::startTimer: QTimer can only be used with threads started with QThread
    QObject::startTimer: QTimer can only be used with threads started with QThread
    QObject::startTimer: QTimer can only be used with threads started with QThread
    QObject::startTimer: QTimer can only be used with threads started with QThread
    

Loading a BDF
-------------

There are two ways to load a BDF; the long way or the short way.

The short way instantiates the ``BDF`` class and the short way uses the
``read_bdf`` function. As this demo was written for the Jupyter
Notebook, we'll use read\_bdf and then mention the other method. The
class-based method allows finer control over things like: - what cards
should be loaded - OpenMDAO dynamic syntax support

The Class-based method
~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    print(bdf_filename)
    
    # create the BDF object
    bdf = BDF()
    
    # read the file from the GUI
    # don't cross-reference
    bdf.read_bdf(bdf_filename, xref=False)


.. parsed-literal::

    f:\work\pynastran\pynastran\master3\models\iSat\ISat_Launch_Sm_Rgd.dat
    DEBUG:     fname=bdf.py                    lineNo=1008   ---starting BDF.read_bdf of f:\work\pynastran\pynastran\master3\models\iSat\ISat_Launch_Sm_Rgd.dat---
    DEBUG:     fname=bdf.py                    lineNo=2999   opening 'f:\\work\\pynastran\\pynastran\\master3\\models\\iSat\\ISat_Launch_Sm_Rgd.dat'
    DEBUG:     fname=bdf.py                    lineNo=1047   ---finished BDF.read_bdf of f:\work\pynastran\pynastran\master3\models\iSat\ISat_Launch_Sm_Rgd.dat---
    

The function-based method
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    bdf = read_bdf(bdf_filename, xref=False)


.. parsed-literal::

    DEBUG:     fname=bdf.py                    lineNo=1008   ---starting BDF.read_bdf of f:\work\pynastran\pynastran\master3\models\iSat\ISat_Launch_Sm_Rgd.dat---
    DEBUG:     fname=bdf.py                    lineNo=2999   opening 'f:\\work\\pynastran\\pynastran\\master3\\models\\iSat\\ISat_Launch_Sm_Rgd.dat'
    DEBUG:     fname=bdf.py                    lineNo=1047   ---finished BDF.read_bdf of f:\work\pynastran\pynastran\master3\models\iSat\ISat_Launch_Sm_Rgd.dat---
    

For simplicity of using the demo, we'll again use the ``read_bdf``
method

.. code:: python

    #bdf_filename = r'D:\work\pynastran_0.8.0_py27\models\iSat\ISat_Launch_Sm_Rgd.dat'
    bdf_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_Rgd.dat'))
    
    # read the file as a path
    bdf_xref = read_bdf(bdf_filename, xref=True)


.. parsed-literal::

    DEBUG:     fname=bdf.py                    lineNo=1008   ---starting BDF.read_bdf of f:\work\pynastran\pynastran\master3\models\iSat\ISat_Launch_Sm_Rgd.dat---
    DEBUG:     fname=bdf.py                    lineNo=2999   opening 'f:\\work\\pynastran\\pynastran\\master3\\models\\iSat\\ISat_Launch_Sm_Rgd.dat'
    DEBUG:     fname=cross_reference.py        lineNo=379    Cross Referencing...
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=1 midsurface: z1=0.400000006 z2=-0.400000006 t=0.035999998 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=2 midsurface: z1=0.400000006 z2=-0.400000006 t=0.054000005 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=3 midsurface: z1=0.400000006 z2=-0.400000006 t=0.017999999 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=7 midsurface: z1=0.418000013 z2=-0.418000013 t=0.035999998 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=34 midsurface: z1=0.194000006 z2=-0.194000006 t=0.0186 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=37 midsurface: z1=0.308999985 z2=-0.308999985 t=0.0186 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=38 midsurface: z1=0.284000009 z2=-0.284000009 t=0.0186 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=46 midsurface: z1=0.199000001 z2=-0.199000001 t=0.0186 not in range of -1.5t < zi < 1.5t
    DEBUG:     fname=bdf.py                    lineNo=1047   ---finished BDF.read_bdf of f:\work\pynastran\pynastran\master3\models\iSat\ISat_Launch_Sm_Rgd.dat---
    

Interrogating the BDF object
----------------------------

IDE's like WingIDE, PyCharm, Spyder and "Python Tools for Visual Studio"
make it very easy to program with their object introspection ability.
Unfortunately, because pyNastran has so many functions, it can be
difficult to learn the code.

Some handy object introspection methods were created that will work on
all pyNastran objects and even non-pyNastran objects. By convention,
private data members/functions start with an underscore \_, and public
ones do not.

We can use the generic object attributes/methods functions

.. code:: python

    print(object_attributes(bdf))
    print(object_methods(bdf))


.. parsed-literal::

    ['MATS1', 'MATS3', 'MATS8', 'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9', 'active_filename', 'active_filenames', 'aecomps', 'aefacts', 'aelinks', 'aelists', 'aeparams', 'aero', 'aeros', 'aestats', 'aesurfs', 'asets', 'bcrparas', 'bcs', 'bctadds', 'bctparas', 'bctsets', 'bdf_filename', 'bsets', 'bsurf', 'bsurfs', 'cMethods', 'caero_ids', 'caeros', 'card_count', 'cards_to_read', 'case_control_deck', 'case_control_lines', 'convection_properties', 'coord_ids', 'coords', 'creep_materials', 'csets', 'csschds', 'dareas', 'dconadds', 'dconstrs', 'ddvals', 'debug', 'delays', 'dequations', 'desvars', 'divergs', 'dlinks', 'dload_entries', 'dloads', 'dmigs', 'dmijis', 'dmijs', 'dmiks', 'dmis', 'doptprm', 'dphases', 'dresps', 'dscreen', 'dtable', 'dumplines', 'dvcrels', 'dvmrels', 'dvprels', 'echo', 'element_ids', 'elements', 'epoints', 'executive_control_lines', 'flfacts', 'flutters', 'frequencies', 'gridSet', 'gusts', 'hyperelastic_materials', 'iSolLine', 'include_dir', 'is_long_ids', 'is_msc', 'is_nx', 'loads', 'log', 'masses', 'material_ids', 'materials', 'methods', 'mkaeros', 'monitor_points', 'mpcadds', 'mpcs', 'nastran_format', 'ncaeros', 'ncoords', 'nelements', 'nlparms', 'nlpcis', 'nmaterials', 'nnodes', 'node_ids', 'nodes', 'nproperties', 'paeros', 'params', 'pbusht', 'pdampt', 'pelast', 'phbdys', 'plotels', 'point_ids', 'points', 'properties', 'properties_mass', 'property_ids', 'punch', 'qsets', 'random_tables', 'read_includes', 'reject_cards', 'reject_count', 'reject_lines', 'rejects', 'rigid_elements', 'rsolmap_toStr', 'se_bsets', 'se_csets', 'se_qsets', 'se_sets', 'se_suport', 'se_usets', 'sets', 'sol', 'solMethod', 'spcadds', 'spcs', 'special_cards', 'splines', 'spoints', 'subcases', 'suport', 'suport1', 'tables', 'tables_sdamping', 'tempds', 'thermal_materials', 'tics', 'transfer_functions', 'trims', 'tstepnls', 'tsteps', 'units', 'usets', 'values_to_skip']
    ['AEFact', 'AELIST', 'AELink', 'AEList', 'AEParam', 'AEStat', 'Acsid', 'Aero', 'Aeros', 'CAero', 'CMethod', 'Coord', 'DConstr', 'DDVal', 'DELAY', 'DEQATN', 'DLoad', 'DMIG', 'DPHASE', 'DResp', 'DVcrel', 'DVmrel', 'DVprel', 'Desvar', 'Element', 'Elements', 'FLFACT', 'Flfact', 'Flutter', 'Gust', 'HyperelasticMaterial', 'Load', 'MPC', 'Mass', 'Material', 'Materials', 'Method', 'NLParm', 'Node', 'Nodes', 'PAero', 'Phbdy', 'Properties', 'Property', 'PropertyMass', 'RandomTable', 'RigidElement', 'SET1', 'SPC', 'Set', 'SetSuper', 'Spline', 'StructuralMaterial', 'Table', 'ThermalMaterial', 'add_AECOMP', 'add_AEFACT', 'add_AELINK', 'add_AELIST', 'add_AEPARM', 'add_AERO', 'add_AEROS', 'add_AESTAT', 'add_AESURF', 'add_ASET', 'add_BCRPARA', 'add_BCTADD', 'add_BCTPARA', 'add_BCTSET', 'add_BSET', 'add_BSURF', 'add_BSURFS', 'add_CAERO', 'add_CSET', 'add_CSSCHD', 'add_DAREA', 'add_DCONSTR', 'add_DDVAL', 'add_DELAY', 'add_DEQATN', 'add_DESVAR', 'add_DIVERG', 'add_DLINK', 'add_DMI', 'add_DMIG', 'add_DMIJ', 'add_DMIJI', 'add_DMIK', 'add_DPHASE', 'add_DRESP', 'add_DTABLE', 'add_DVCREL', 'add_DVMREL', 'add_DVPREL', 'add_FLFACT', 'add_FLUTTER', 'add_FREQ', 'add_LSEQ', 'add_MKAERO', 'add_MONPNT', 'add_NLPARM', 'add_NLPCI', 'add_PAERO', 'add_PARAM', 'add_PBUSHT', 'add_PDAMPT', 'add_PELAST', 'add_PHBDY', 'add_QSET', 'add_SEBSET', 'add_SECSET', 'add_SEQSET', 'add_SESET', 'add_SET', 'add_SEUSET', 'add_SPLINE', 'add_TEMPD', 'add_TF', 'add_TRIM', 'add_TSTEP', 'add_TSTEPNL', 'add_USET', 'add_card', 'add_card_class', 'add_card_fields', 'add_card_lines', 'add_cmethod', 'add_constraint', 'add_constraint_MPC', 'add_constraint_SPC', 'add_convection_property', 'add_coord', 'add_creep_material', 'add_damper', 'add_dload', 'add_dload_entry', 'add_element', 'add_epoint', 'add_gust', 'add_hyperelastic_material', 'add_load', 'add_mass', 'add_material_dependence', 'add_method', 'add_node', 'add_plotel', 'add_property', 'add_property_mass', 'add_random_table', 'add_rigid_element', 'add_sesuport', 'add_spoint', 'add_structural_material', 'add_suport', 'add_suport1', 'add_table', 'add_table_sdamping', 'add_thermal_BC', 'add_thermal_element', 'add_thermal_load', 'add_thermal_material', 'auto_reject_bdf', 'create_card_object', 'create_card_object_fields', 'create_card_object_list', 'cross_reference', 'deprecated', 'disable_cards', 'echo_bdf', 'fill_dmigs', 'geom_check', 'getElementIDsWithPID', 'getNodes', 'get_bdf_cards', 'get_bdf_cards_dict', 'get_bdf_stats', 'get_card_ids_by_card_types', 'get_cards_by_card_types', 'get_dependent_nid_to_components', 'get_displacement_index_transforms', 'get_dload_entries', 'get_element_ids_dict_with_pids', 'get_element_ids_list_with_pids', 'get_encoding', 'get_material_id_to_property_ids_map', 'get_material_ids', 'get_mpcs', 'get_node_id_to_element_ids_map', 'get_node_id_to_elements_map', 'get_node_ids_with_element', 'get_node_ids_with_elements', 'get_property_id_to_element_ids_map', 'get_reduced_mpcs', 'get_reduced_spcs', 'get_rigid_elements_with_node_ids', 'get_solid_skin_faces', 'get_spcs', 'get_structural_material_ids', 'get_thermal_material_ids', 'get_xyz_in_coord', 'is_reject', 'load', 'mass_properties', 'object_attributes', 'object_methods', 'pop_parse_errors', 'pop_xref_errors', 'process_card', 'read_bdf', 'remove_unassociated_properties', 'replace_cards', 'resolve_grids', 'safe_cross_reference', 'save', 'set_as_msc', 'set_as_nx', 'set_dynamic_syntax', 'set_error_storage', 'skin_solid_elements', 'sum_forces_moments', 'sum_forces_moments_elements', 'uncross_reference', 'unresolve_grids', 'update_solution', 'validate', 'write_bdf', 'write_bdf_symmetric', 'write_caero_model', 'write_skin_solid_faces']
    

Let's clean that up a bit
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    print("attributes = [%s]\n" % ', '.join(bdf.object_attributes()))
    print("methods = [%s]\n" % ', '.join(bdf.object_methods()))


.. parsed-literal::

    attributes = [MATS1, MATS3, MATS8, MATT1, MATT2, MATT3, MATT4, MATT5, MATT8, MATT9, active_filename, active_filenames, aecomps, aefacts, aelinks, aelists, aeparams, aero, aeros, aestats, aesurfs, asets, bcrparas, bcs, bctadds, bctparas, bctsets, bdf_filename, bsets, bsurf, bsurfs, cMethods, caeros, card_count, cards_to_read, case_control_deck, case_control_lines, convection_properties, coords, creep_materials, csets, csschds, dareas, dconadds, dconstrs, ddvals, debug, delays, dequations, desvars, divergs, dlinks, dload_entries, dloads, dmigs, dmijis, dmijs, dmiks, dmis, doptprm, dphases, dresps, dscreen, dtable, dumplines, dvcrels, dvmrels, dvprels, echo, elements, epoints, executive_control_lines, flfacts, flutters, frequencies, gridSet, gusts, hyperelastic_materials, iSolLine, include_dir, is_msc, is_nx, loads, masses, materials, methods, mkaeros, monitor_points, mpcadds, mpcs, nastran_format, nlparms, nlpcis, nodes, paeros, params, pbusht, pdampt, pelast, phbdys, plotels, points, properties, properties_mass, punch, qsets, random_tables, read_includes, reject_cards, reject_count, reject_lines, rejects, rigid_elements, rsolmap_toStr, se_bsets, se_csets, se_qsets, se_sets, se_suport, se_usets, sets, sol, solMethod, spcadds, spcs, special_cards, splines, spoints, suport, suport1, tables, tables_sdamping, tempds, thermal_materials, tics, transfer_functions, trims, tstepnls, tsteps, units, usets, values_to_skip]
    
    methods = [AEFact, AELIST, AELink, AEList, AEParam, AEStat, Acsid, Aero, Aeros, CAero, CMethod, Coord, DConstr, DDVal, DELAY, DEQATN, DLoad, DMIG, DPHASE, DResp, DVcrel, DVmrel, DVprel, Desvar, Element, Elements, FLFACT, Flfact, Flutter, Gust, HyperelasticMaterial, Load, MPC, Mass, Material, Materials, Method, NLParm, Node, Nodes, PAero, Phbdy, Properties, Property, PropertyMass, RandomTable, RigidElement, SET1, SPC, Set, SetSuper, Spline, StructuralMaterial, Table, ThermalMaterial, add_AECOMP, add_AEFACT, add_AELINK, add_AELIST, add_AEPARM, add_AERO, add_AEROS, add_AESTAT, add_AESURF, add_ASET, add_BCRPARA, add_BCTADD, add_BCTPARA, add_BCTSET, add_BSET, add_BSURF, add_BSURFS, add_CAERO, add_CSET, add_CSSCHD, add_DAREA, add_DCONSTR, add_DDVAL, add_DELAY, add_DEQATN, add_DESVAR, add_DIVERG, add_DLINK, add_DMI, add_DMIG, add_DMIJ, add_DMIJI, add_DMIK, add_DPHASE, add_DRESP, add_DTABLE, add_DVCREL, add_DVMREL, add_DVPREL, add_FLFACT, add_FLUTTER, add_FREQ, add_LSEQ, add_MKAERO, add_MONPNT, add_NLPARM, add_NLPCI, add_PAERO, add_PARAM, add_PBUSHT, add_PDAMPT, add_PELAST, add_PHBDY, add_QSET, add_SEBSET, add_SECSET, add_SEQSET, add_SESET, add_SET, add_SEUSET, add_SPLINE, add_TEMPD, add_TF, add_TRIM, add_TSTEP, add_TSTEPNL, add_USET, add_card, add_card_class, add_card_fields, add_card_lines, add_cmethod, add_constraint, add_constraint_MPC, add_constraint_SPC, add_convection_property, add_coord, add_creep_material, add_damper, add_dload, add_dload_entry, add_element, add_epoint, add_gust, add_hyperelastic_material, add_load, add_mass, add_material_dependence, add_method, add_node, add_plotel, add_property, add_property_mass, add_random_table, add_rigid_element, add_sesuport, add_spoint, add_structural_material, add_suport, add_suport1, add_table, add_table_sdamping, add_thermal_BC, add_thermal_element, add_thermal_load, add_thermal_material, auto_reject_bdf, create_card_object, create_card_object_fields, create_card_object_list, cross_reference, deprecated, disable_cards, echo_bdf, fill_dmigs, geom_check, getElementIDsWithPID, getNodes, get_bdf_cards, get_bdf_cards_dict, get_bdf_stats, get_card_ids_by_card_types, get_cards_by_card_types, get_dependent_nid_to_components, get_displacement_index_transforms, get_dload_entries, get_element_ids_dict_with_pids, get_element_ids_list_with_pids, get_encoding, get_material_id_to_property_ids_map, get_material_ids, get_mpcs, get_node_id_to_element_ids_map, get_node_id_to_elements_map, get_node_ids_with_element, get_node_ids_with_elements, get_property_id_to_element_ids_map, get_reduced_mpcs, get_reduced_spcs, get_rigid_elements_with_node_ids, get_solid_skin_faces, get_spcs, get_structural_material_ids, get_thermal_material_ids, get_xyz_in_coord, is_reject, load, mass_properties, pop_parse_errors, pop_xref_errors, process_card, read_bdf, remove_unassociated_properties, replace_cards, resolve_grids, safe_cross_reference, save, set_as_msc, set_as_nx, set_dynamic_syntax, set_error_storage, skin_solid_elements, sum_forces_moments, sum_forces_moments_elements, uncross_reference, unresolve_grids, update_solution, validate, write_bdf, write_bdf_symmetric, write_caero_model, write_skin_solid_faces]
    
    

Some other very handy methods that will be used later by ``test_bdf``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    print(bdf.get_bdf_stats())
    print("card_count = %s\n" % bdf.card_count)
    print("reject_count = %s" % bdf.reject_count)


.. parsed-literal::

    ---BDF Statistics---
    SOL 103
    
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
    
    bdf.materials
      MAT1     : 14
      MAT8     : 8
    
    bdf.coords
      CORD2R   : 75
    
    bdf.methods
      EIGRL    : 1
    
    
    card_count = {u'ENDDATA': 1, u'CQUAD4': 4580, u'PSOLID': 4, u'PARAM': 8, u'PBUSH': 2, u'CBUSH': 104, u'PBARL': 18, u'MAT1': 14, u'CTRIA3': 32, u'CORD2R': 75, u'PSHELL': 8, u'USET': 1, u'GRID': 5380, u'EIGRL': 1, u'CBAR': 827, u'PBAR': 1, u'CHEXA': 25, u'SPC': 1, u'CONM2': 15, u'RBE2': 44, u'MAT8': 8}
    
    reject_count = {}
    

Cross-referencing
-----------------

Cross-referencing a BDF allows improved usability of the BDF class. It
comes with some negative side effects, but in general is a very useful
thing. It dramatically minimizes the amount of code you need to write,
greatly simplifies future operations, and is highly recommended.

The major downside is it prevents decks from being saved to object files
for faster loading.

Without Cross-Referencing (xref=False)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here the raw values of the the data objects are returned to us

.. code:: python

    cquad = bdf.elements[1]
    print(cquad)
    nid1 = cquad.nodes[0]
    print("nid1 = %s" % nid1)
    n1 = bdf.nodes[nid1]
    cd4 = n1.cd
    c4 = bdf.coords[cd4]
    print("i xref=False %s" % str(c4.i))
    #print object_attributes(c4)


.. parsed-literal::

    $*
    $*  ELEMENT CARDS
    $*
    CQUAD4         1       1       1       2       4       3
    
    nid1 = 1
    i xref=False [ 1.  0.  0.]
    

Cross-Referenced (xref=True)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we can trace the referenced objects very easily

.. code:: python

    print("i xref=True %s" % bdf_xref.elements[1].nodes[0].cd.i)


.. parsed-literal::

    i xref=True [ 1.  0.  0.]
    

So how is this done?

.. code:: python

    cquad.nodes[0] = n1
    print(cquad.nodes[0])


.. parsed-literal::

    $*
    $*  GRID CARDS
    $*
    GRID           1       4    -4.5    -7.5    -14.       4
    
    

Let's show off the GRID card
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

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
    
    c1=u'$*\n$*  GRID CARDS\n$*\n'
    c2=u''
    xyz = [ -4.5  -7.5 -14. ]
    position = [ -4.5  -7.5 -14. ]
    wrt5 = [  2.12132034  14.         -26.59188309]
    wrt4 = [ -4.5  -7.5 -14. ]
    

Now let's modify the GRID card and write it out

.. code:: python

    n1 = bdf_xref.nodes[1]
    n1.xyz[1] = -7.5
    print("repr  = %s" % n1.repr_fields())
    print("raw   = %s" % n1.repr_fields())
    
    #n1.xyz[1] = 100000000000.
    print("repr2 = %s" % n1.repr_fields())
    print(n1)
    print(n1.write_card(size=8))
    print(n1.write_card(size=16, is_double=False))
    print(n1.write_card(size=16, is_double=True))


.. parsed-literal::

    repr  = [u'GRID', 1, 4, -4.5, -7.5, -14.0, 4, u'', None]
    raw   = [u'GRID', 1, 4, -4.5, -7.5, -14.0, 4, u'', None]
    repr2 = [u'GRID', 1, 4, -4.5, -7.5, -14.0, 4, u'', None]
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

.. code:: python

    mass, cg, I = bdf_xref.mass_properties()
    print("mass = %s" % mass)


.. parsed-literal::

    DEBUG:     fname=bdf_methods.py            lineNo=524    Mass/MOI sym_axis = []
    mass = 1.70202557344
    

Examples of xref on elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    eid100 = bdf_xref.elements[100]
    print(eid100)
    print("nodes = %s" % eid100.nodes)
    print("--node0--\n%s" % eid100.nodes[0])
    print("--cd--\n%s" % eid100.nodes[0].cd)
    print("cd.cid = %s" % eid100.nodes[0].cd.cid)
    
    print("area = %s" % eid100.Area())
    print("mass = %s" % eid100.Mass())
    print("--pid--\n%s" % eid100.pid)
    print("pid.pid = %s" % eid100.pid.pid)
    print("pid.Pid() = %s" % eid100.Pid())
    
    print(eid100.pid.mid1)
    print("type = %s" % eid100.pid.mid1.type)
    print("nu12 = %s" % eid100.pid.mid1.nu12)
    print("mass = %s" % eid100.Mass())


.. parsed-literal::

    CQUAD4       100       1     149     152     161     160
    
    nodes = [GRID         149       4      3.     7.5   -16.5       4
    , GRID         152       4     1.5     7.5   -16.5       4
    , GRID         161       4     1.5     7.5    -14.       4
    , GRID         160       4      3.     7.5    -14.       4
    ]
    --node0--
    GRID         149       4      3.     7.5   -16.5       4
    
    --cd--
    CORD2R         4              0.      0.      0.      0.      0.      1.
                  1.      0.      0.
    
    cd.cid = 4
    area = 3.75
    mass = 3.6428803074e-05
    --pid--
    $*
    $*  PROPERTY CARDS
    $*
    $*
    $*  I-DEAS property: 1  name: BUS PNL HCMB 2PLY
    PSHELL         1       6    .036       61415.815       7         3.551-6
                  .4     -.4
    
    pid.pid = 1
    pid.Pid() = 1
    $*
    $*  I-DEAS Material: 6  name: BUS_CFRP_PW_ORTHO
    $* M46J PW ETW
    MAT8           6   1.7+7   1.7+7     .98 340000. 180000. 180000..0001712
                               71.33
    
    type = MAT8
    nu12 = 0.98
    mass = 3.6428803074e-05
    

Write the modified deck
-----------------------

Let's first switch to the desktop to make the file easy to find

.. code:: python

    import getpass
    name = getpass.getuser()
    os.chdir(os.path.join(r'C:\Users', name, 'Desktop'))

.. code:: python

    pwd




.. parsed-literal::

    u'C:\\Users\\Steve\\Desktop'



There are two ways to write a deck - interspersed : alternate properties
and elements (similar to how Patran writes decks) - not-interspersed
(default) : much faster

We can also use 8 or 16 character field width as well as double
precision.

Note that double precision only works for certain cards (e.g. GRID,
COORD, DMIG) and not much else.

.. code:: python

    bdf_xref.write_bdf('fem.bdf', interspersed=False, size=8, is_double=False)
    !tail -n 5 "fem.bdf"
    
    bdf_xref.write_bdf('fem.bdf', interspersed=True, size=16, is_double=False)
    !tail "fem.bdf"
    
    bdf_xref.write_bdf('fem.bdf', interspersed=True, size=16, is_double=True)
    !tail "fem.bdf"


.. parsed-literal::

    CORD2R        75        1.355-13-2.19-15    -40.1.355-13-2.19-15      0.
                 40.-2.19-15    -40.
    CORD2R        76        1.355-13-2.19-15    -40.1.355-13-2.19-15      0.
                 40.-2.19-15    -40.
    ENDDATA
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
    

.. code:: python

    bdf_filename




.. parsed-literal::

    'f:\\work\\pynastran\\pynastran\\master3\\models\\iSat\\ISat_Launch_Sm_Rgd.dat'



pyNastranGUI
------------

.. code:: python

    print(bdf_filename)
    %echo {bdf_filename}
    #!pyNastranGUI -f nastran -i {bdf_filename}
    
    solid_bending_bdf = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending.bdf'))
    solid_bending_op2 = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending.op2'))
    
    !pyNastranGUI -f nastran -i {solid_bending_bdf} -o {solid_bending_op2}  > junk.out
    print("done")


.. parsed-literal::

    f:\work\pynastran\pynastran\master3\models\iSat\ISat_Launch_Sm_Rgd.dat
    f:\work\pynastran\pynastran\master3\models\iSat\ISat_Launch_Sm_Rgd.dat
    done
    

.. parsed-literal::

    QObject::startTimer: QTimer can only be used with threads started with QThread
    QObject::startTimer: QTimer can only be used with threads started with QThread
    QObject::startTimer: QTimer can only be used with threads started with QThread
    QObject::startTimer: QTimer can only be used with threads started with QThread
    

We can also script the GUI!
===========================

.. code:: python

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

    /cygdrive/c/Users/Steve/Desktop
    

.. parsed-literal::

    QObject::startTimer: QTimer can only be used with threads started with QThread
    QObject::startTimer: QTimer can only be used with threads started with QThread
    QObject::startTimer: QTimer can only be used with threads started with QThread
    QObject::startTimer: QTimer can only be used with threads started with QThread
    


.. image:: bdf_demo_files%5Cbdf_demo_41_2.png


.. parsed-literal::

    the picture is visible
    


