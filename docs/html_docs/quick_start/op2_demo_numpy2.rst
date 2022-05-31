
OP2: Numpy Demo #2 (Composite Plate Stress)
===========================================

The Jupyter notebook for this demo can be found in: -
docs/quick_start/demo/op2_demo_numpy1.ipynb -
https://github.com/SteveDoyle2/pyNastran/tree/main/docs/quick_start/demo/op2_demo_numpy1.ipynb

It’s recommended that you first go through: -
https://github.com/SteveDoyle2/pyNastran/tree/main/docs/quick_start/demo/op2_intro.ipynb
-
https://github.com/SteveDoyle2/pyNastran/tree/main/docs/quick_start/demo/op2_demo.ipynb
-
https://github.com/SteveDoyle2/pyNastran/tree/main/docs/quick_start/demo/op2_demo_numpy1.ipynb

In this tutorial, composite plate stresses will be covered.

Load the model
--------------

If the BWB example OP2 doesn’t exist, we’ll run Nastran to create it.

.. code:: ipython3

    import os
    import copy
    import numpy as np
    np.set_printoptions(precision=2, threshold=20, linewidth=100, suppress=True)

    import pyNastran
    from pyNastran.op2.op2 import read_op2
    from pyNastran.utils.nastran_utils import run_nastran
    pkg_path = pyNastran.__path__[0]
    model_path = os.path.join(pkg_path, '..', 'models')

    bdf_filename = os.path.join(model_path, 'bwb', 'bwb_saero.bdf')
    op2_filename = os.path.join(model_path, 'bwb', 'bwb_saero.op2')
    if not os.path.exists(op2_filename):
        keywords = ['scr=yes', 'bat=no', 'old=no']
        run_nastran(bdf_filename, nastran_cmd='nastran', keywords=keywords, run=True)
        import shutil
        op2_filename2 = os.path.join('bwb_saero.op2')
        shutil.move(op2_filename2, op2_filename)

    assert os.path.exists(op2_filename), print_bad_path(op2_filename)
    model = read_op2(op2_filename, build_dataframe=False, debug=False)

    print(model.get_op2_stats(short=True))



.. raw:: html

    <text style=color:green>INFO:    op2_scalar.py:1469           op2_filename = 'c:\\nasa\\m4\\formats\\git\\pynastran_1.2\\pyNastran\\..\\models\\bwb\\bwb_saero.op2'
    </text>


.. parsed-literal::

    displacements[1]
    spc_forces[1]
    grid_point_forces[1]
    cbar_stress[1]
    cbar_strain[1]
    cquad4_composite_stress[1]
    ctria3_composite_stress[1]
    cquad4_composite_strain[1]
    ctria3_composite_strain[1]



Accessing the Composite Stress
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    isubcase = 1
    stress = model.cquad4_composite_stress[isubcase]
    print(stress)
    headers = stress.get_headers()
    imax = headers.index('major')


.. parsed-literal::

      type=RealCompositePlateStressArray nelements=9236 ntotal=92360
      data: [1, ntotal, 9] where 9=[o11, o22, t12, t1z, t2z, angle, major, minor, max_shear]
      element_layer.shape = (92360, 2)
      data.shape = (1, 92360, 9)
      element type: QUAD4LC-composite
      sort1
      lsdvmns = [1]



Composite Stress/Strain data is tricky to access as there is not a good way to index the data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let’s cheat a bit using the element ids and layers to make a pivot
table. - **table** is (ntimes, nelements, nlayers, ndata) -
**max_principal_stress_table** is (nelements, nlayers)

.. code:: ipython3

    from pyNastran.femutils.utils import pivot_table

    eids = stress.element_layer[:, 0]
    layers = stress.element_layer[:, 1]

    ## now pivot the stress
    table, rows_new = pivot_table(stress.data, eids, layers)

    # now access the max principal stress for the static result
    # table is (itime, nelements, nlayers, data)
    itime = 0
    max_principal_stress_table = table[itime,:,:,imax]
    ueids = np.unique(eids)
    print('max_principal_stress_table:\n%s' % max_principal_stress_table)


.. parsed-literal::

    max_principal_stress_table:
    [[ 239.3   163.91   98.41 ...  -35.77  -34.6   -19.86]
     [  18.61   78.52   25.52 ...  -63.92  -62.48  -12.99]
     [   2.99  105.48   49.37 ... -137.74 -127.07  -41.14]
     ...
     [ 157.    170.3   112.79 ...   44.56   47.13   38.9 ]
     [ 123.96  143.01   97.41 ...   40.99   44.06   42.47]
     [  90.04  109.97   79.86 ...   33.18   36.12   24.04]]


More realistic pivot table
--------------------------

All the elements have 10 layers. Let’s remove the last 5 layers.

By having empty layers, the pivot table now has nan data in it.

.. code:: ipython3

    # drop out 5 layers
    eids2 = stress.element_layer[:-5, 0]
    layers2 = stress.element_layer[:-5, 1]
    data2 = stress.data[:, :-5, :]

    # now pivot the stress
    table, rows_new = pivot_table(data2, eids2, layers2)

    # access the table data
    # table is (itime, nelements, nlayers, data)
    itime = 0
    max_principal_stress_table2 = table[itime,:,:,imax]
    print('max_principal_stress_table2:\n%s' % max_principal_stress_table2)


.. parsed-literal::

    max_principal_stress_table2:
    [[ 239.3   163.91   98.41 ...  -35.77  -34.6   -19.86]
     [  18.61   78.52   25.52 ...  -63.92  -62.48  -12.99]
     [   2.99  105.48   49.37 ... -137.74 -127.07  -41.14]
     ...
     [ 157.    170.3   112.79 ...   44.56   47.13   38.9 ]
     [ 123.96  143.01   97.41 ...   40.99   44.06   42.47]
     [  90.04  109.97   79.86 ...     nan     nan     nan]]

