OP2 Introduction
================

The Jupyter notebook for this demo can be found in: -
docs/quick_start/demo/op2_demo.ipynb -
https://github.com/SteveDoyle2/pyNastran/tree/master/docs/quick_start/demo/op2_demo.ipynb

Why use the OP2? Why not use the F06/PCH file?
----------------------------------------------

Most people are comfortable with the F06. However, it’s: - Ironically, a
lot harder to parse. The OP2 is very structured. - Much, much, much
slower. We can read entire blocks of arrays with a single call. The data
is already typed. - Much, much more memory inefficient because we aren’t
appending strings onto lists and turning that into a numpy array.

F06 parsers get ridiculously hard when you start do complicated results,
like: - single subcase buckling - superelements - SOL 200 optimization
with sub-optimization - SPOINTs

The pyNastran OP2 Reader is fast, highly validated, and it supports most
result types. The data in the OP2 is also more accurate because there is
no rounding.

Validating an OP2
-----------------

The ``test_op2`` script is created when you run
``python setup.py develop`` or ``python setup.py install`` on pyNastran.
Assuming it’s on your path (it’ll be in Python27:raw-latex:`\Scripts `or
something similar), you can run:

::

   >>> test_op2 -f solid_bending.op2

The ``-f`` tells us to print out ``solid_bending.test_op2.f06``, which
can be compared to your F06 for a small file to build confidence in the
reader. It’s also useful when you want an F06 of your model without
rerunning Nastran just to see what’s in it.

If you have a large model, you can make ``test_op2`` run much, much
faster. The ``-c`` flag disables double-reading of the OP2. By default,
``test_op2`` uses two different read methods (the old method and new
method) to ensure that results are read in properly. When running the
code, this is turned off, but is turned on for ``test_op2``.

::

   >>> test_op2 -fc solid_bending.op2

Import the packages
-------------------

.. code:: ipython3

    import os
    import copy
    import numpy as np
    np.set_printoptions(precision=2, threshold=20, suppress=True)
    
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    
    from pyNastran.utils import print_bad_path
    from pyNastran.op2.op2 import read_op2
    from pyNastran.utils import object_methods, object_attributes
    from pyNastran.utils.nastran_utils import run_nastran
    
    import pandas as pd

Sets default precision of real numbers for pandas output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    pd.set_option('precision', 3)
    np.set_printoptions(precision=3, threshold=20)

As with the BDF, we can use the long form and the short form. However,
the long form for the ``OP2`` doesn’t really add anything. So, let’s
just use the short form.

In addition to the default numpy support, there is also **``pandas``**
dataframe support.

.. code:: ipython3

    #op2_filename = r'D:\work\pynastran_0.8.0\models\iSat\ISat_Launch_Sm_Rgd.op2'
    #op2_filename = r'D:\work\pynastran_0.8.0\models\iSat\ISat_Launch_Sm_4pt.op2'
    bdf_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_4pt.dat'))
    op2_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_4pt.op2'))
    
    if 1:
        from pyNastran.bdf.bdf import read_bdf
        op2_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_4pt2.op2'))
        model = read_bdf(bdf_filename, debug=None)
        model.set_param('POSTEXT', 'YES')
        model.set_param('POST', -2)
        bdf_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_4pt2.bdf'))
        model.write_bdf(bdf_filename)
    
    if not os.path.exists(op2_filename):
        run_nastran(bdf_filename)
    assert os.path.exists(op2_filename), print_bad_path(op2_filename)
    
    # define the input file with a file path
    op2 = read_op2(op2_filename, build_dataframe=True, debug=False)



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=1 midsurface: z1=0.4 z2=-0.4 t=0.036 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=2 midsurface: z1=0.4 z2=-0.4 t=0.054 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=3 midsurface: z1=0.4 z2=-0.4 t=0.018 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=7 midsurface: z1=0.418 z2=-0.418 t=0.036 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=34 midsurface: z1=0.194 z2=-0.194 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=38 midsurface: z1=0.284 z2=-0.284 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=46 midsurface: z1=0.199 z2=-0.199 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=37 midsurface: z1=0.309 z2=-0.309 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:green>INFO:    op2_scalar.py:1588           op2_filename = 'c:\\nasa\\m4\\formats\\git\\pynastran\\models\\iSat\\ISat_Launch_Sm_4pt2.op2'
    </text>


OP2 Introspection
-----------------

The ``get_op2_stats()`` function lets you quickly understand what in an
op2.

.. code:: ipython3

    print(op2.get_op2_stats())


.. parsed-literal::

    params:
      AUTOSPC = 'YES'
      GRDPNT = 0
      K6ROT = 100.0
      OMODES = 11
      POST = -2
      POSTEXT = 'YES'
    GridPointWeight['']:  reference_point=0
      mass=[    1.7746     1.7746     1.7746]
      cg  =[-6.02244e-18    -2.5306   -18.4677]
           [-0.0338514 -1.01609e-17   -18.4677]
           [-0.0338514    -2.5306 -6.56299e-20]
    
      IS  =[    705.69   -1.56673   0.141188]
           [  -1.56673    621.837    135.836]
           [  0.141188    135.836    415.862]
    
      IQ  =[   689.184                      ]
           [              348.385           ]
           [                         705.821]
    
      Q  = [ 0.0884613 0.00159687   0.996078]
           [ -0.892013  -0.444887  0.0799325]
           [   0.44327  -0.895585 -0.0379308]
    op2_results.eqexin: EQEXIN(nid, ndof, doftype); nnodes=5379
    op2_results.bgpdt: BGPDT(cd, xyz); nnodes=5379
    op2_results.force.cbar_force[1]
      type=RealCBarForceArray ntimes=167 nelements=827; table_name='OEF1X'
      data: [ntimes, nnodes, 8] where 8=[bending_moment_a1, bending_moment_a2, bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque]
      data.shape = (167, 827, 8)
      element.shape = (827,)
      element name: CBAR-34
      sort1
      modes = [  1   2   3 ... 165 166 167]
      eigns = [   2757.98     3568.25     9686.269 ... 6162759.    6169884.5
     6229575.5  ]
      cycles = [  8.358   9.507  15.664 ... 395.101 395.329 397.237]
    
    eigenvectors[1]
      isubcase = 1
      type=RealEigenvectorArray ntimes=167 nnodes=5379, table_name=OUGV1
      data: [t1, t2, t3, r1, r2, r3] shape=[167, 5379, 6] dtype=float32
      node_gridtype.shape = (5379, 2)
      sort1
      modes = [  1   2   3 ... 165 166 167]
      eigns = [   2757.98     3568.25     9686.269 ... 6162759.    6169884.5
     6229575.5  ]
      mode_cycles = [  8.358   9.507  15.664 ... 395.101 395.329 397.237]
    
    ctria3_stress[1]
      type=RealPlateStressArray ntimes=167 nelements=32 nnodes_per_element=1 nlayers=2 ntotal=64
      data: [ntimes, ntotal, 8] where 8=[fiber_distance, oxx, oyy, txy, angle, omax, omin, von_mises]
      element_node.shape = (64, 2)
      data.shape=(167, 64, 8)
      element type: CTRIA3-74
      s_code: 1
      sort1
      modes = [  1   2   3 ... 165 166 167]
      eigns = [   2757.98     3568.25     9686.269 ... 6162759.    6169884.5
     6229575.5  ]
      mode2s = [0 0 0 ... 0 0 0]
      cycles = [  8.358   9.507  15.664 ... 395.101 395.329 397.237]
    
    cquad4_stress[1]
      type=RealPlateStressArray ntimes=167 nelements=4580 nnodes_per_element=1 nlayers=2 ntotal=9160
      data: [ntimes, ntotal, 8] where 8=[fiber_distance, oxx, oyy, txy, angle, omax, omin, von_mises]
      element_node.shape = (9160, 2)
      data.shape=(167, 9160, 8)
      element type: CQUAD4-33
      s_code: 1
      sort1
      modes = [  1   2   3 ... 165 166 167]
      eigns = [   2757.98     3568.25     9686.269 ... 6162759.    6169884.5
     6229575.5  ]
      mode2s = [0 0 0 ... 0 0 0]
      cycles = [  8.358   9.507  15.664 ... 395.101 395.329 397.237]
    
    eigenvalues[ISAT_SM_LAUNCH_4PT MODES TO 400 HZ]
      type=RealEigenvalues neigenvalues=167
      title, extraction_order, eigenvalues, radians, cycles, generalized_mass, generalized_stiffness
    
    
    

If that’s too long…
~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    print(op2.get_op2_stats(short=True))


.. parsed-literal::

    params:
      AUTOSPC = 'YES'
      GRDPNT = 0
      K6ROT = 100.0
      OMODES = 11
      POST = -2
      POSTEXT = 'YES'
    GridPointWeight['']: ref_point=0 mass=1.7746; [reference_point, M0, S, mass, cg, IS, IQ, Q]
    op2_results.eqexin: EQEXIN(nid, ndof, doftype); nnodes=5379
    op2_results.bgpdt: BGPDT(cd, xyz); nnodes=5379
    op2_results.force.cbar_force[1]
    eigenvectors[1]
    ctria3_stress[1]
    cquad4_stress[1]
    eigenvalues['ISAT_SM_LAUNCH_4PT MODES TO 400 HZ']
    
    

Accessing the Eigenvectors object
----------------------------------

Eigenvectors are the simplest object. They use the same class as for
displacements, velocity, acceleration, SPC Forces, MPC Forces, Applied
Loads, etc. These are all node-based tables with TX, TY, TZ, RX, RY, RZ.
Results are in the analysis coordinate frame (CD), which is defined by
the GRID card.

Numpy-based Approach
~~~~~~~~~~~~~~~~~~~~

We’ll first show off the standard ``numpy`` based results on a transient
case. Static results are the same, except that you’ll always use the 0th
index for the “time” index.

The tutorial is intetionally just accessing the objects in a very clear,
though inefficient way. The OP2 objects can take full advantage of the
numpy operations.

.. code:: ipython3

    # what modes did we analyze:  1 to 167
    eigenvector_keys = list(op2.eigenvectors.keys())
    print("loadcases = %s" % eigenvector_keys)
    
    # get subcase 1
    eig1 = op2.eigenvectors[1]
    
    modes = eig1.modes
    times = eig1._times #  the generic version of modes
    print("modes = %s\n" % modes)
    print("times = %s\n" % times)
    
    imode2 = 1 # corresponds to mode 2
    mode2 = eig1.data[imode2, :, :]
    
    print('first 10 nodes and grid types\nNid Gridtype\n%s' % eig1.node_gridtype[:10, :])
    node_ids = eig1.node_gridtype[:, 0]
    
    index_node10 = np.where(node_ids == 10)[0]  # we add the [0] because it's 1d
    mode2_node10 = mode2[index_node10]
    print("translation mode2_node10 = %s" % eig1.data[imode2, index_node10, :3].ravel())
    print("rotations mode2_node10 = %s" % eig1.data[imode2, index_node10, 3:].ravel())


.. parsed-literal::

    loadcases = [1]
    modes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167]
    
    times = [  1.   2.   3. ... 165. 166. 167.]
    
    first 10 nodes and grid types
    Nid Gridtype
    [[ 1  1]
     [ 2  1]
     [ 3  1]
     [ 4  1]
     [ 5  1]
     [ 6  1]
     [ 7  1]
     [ 8  1]
     [ 9  1]
     [10  1]]
    translation mode2_node10 = [0.    0.008 0.002]
    rotations mode2_node10 = [-0.  0. -0.]
    

Pandas-based Approach
~~~~~~~~~~~~~~~~~~~~~

If you like pandas, you can access all the OP2 objects, which is very
useful within the Jupyter Notebook. Different objects will look
differently, but you can change the layout.

If you’re trying to learn pandas, there are many tutorials online, such
as: http://pandas.pydata.org/pandas-docs/stable/10min.html

or a very long, but good video:

.. code:: ipython3

    from IPython.display import YouTubeVideo
    YouTubeVideo('5JnMutdy6Fw')
    #https://www.youtube.com/watch?v=5JnMutdy6Fw




.. raw:: html

    
    <iframe
        width="400"
        height="300"
        src="https://www.youtube.com/embed/5JnMutdy6Fw"
        frameborder="0"
        allowfullscreen
    ></iframe>
    



.. code:: ipython3

    # get subcase 1
    eig1 = op2.eigenvectors[1]
    
    eig1.data_frame




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead tr th {
            text-align: left;
        }
    
        .dataframe thead tr:last-of-type th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr>
          <th></th>
          <th>Mode</th>
          <th>1</th>
          <th>2</th>
          <th>3</th>
          <th>4</th>
          <th>5</th>
          <th>6</th>
          <th>7</th>
          <th>8</th>
          <th>9</th>
          <th>10</th>
          <th>...</th>
          <th>158</th>
          <th>159</th>
          <th>160</th>
          <th>161</th>
          <th>162</th>
          <th>163</th>
          <th>164</th>
          <th>165</th>
          <th>166</th>
          <th>167</th>
        </tr>
        <tr>
          <th></th>
          <th>Freq</th>
          <th>8.358</th>
          <th>9.507</th>
          <th>15.664</th>
          <th>20.229</th>
          <th>20.306</th>
          <th>20.548</th>
          <th>21.500</th>
          <th>21.701</th>
          <th>21.716</th>
          <th>28.444</th>
          <th>...</th>
          <th>382.715</th>
          <th>385.300</th>
          <th>387.258</th>
          <th>390.518</th>
          <th>390.989</th>
          <th>391.049</th>
          <th>393.165</th>
          <th>395.101</th>
          <th>395.329</th>
          <th>397.237</th>
        </tr>
        <tr>
          <th></th>
          <th>Eigenvalue</th>
          <th>2.758e+03</th>
          <th>3.568e+03</th>
          <th>9.686e+03</th>
          <th>1.615e+04</th>
          <th>1.628e+04</th>
          <th>1.667e+04</th>
          <th>1.825e+04</th>
          <th>1.859e+04</th>
          <th>1.862e+04</th>
          <th>3.194e+04</th>
          <th>...</th>
          <th>5.782e+06</th>
          <th>5.861e+06</th>
          <th>5.921e+06</th>
          <th>6.021e+06</th>
          <th>6.035e+06</th>
          <th>6.037e+06</th>
          <th>6.103e+06</th>
          <th>6.163e+06</th>
          <th>6.170e+06</th>
          <th>6.230e+06</th>
        </tr>
        <tr>
          <th></th>
          <th>Radians</th>
          <th>52.516</th>
          <th>59.735</th>
          <th>98.419</th>
          <th>127.102</th>
          <th>127.585</th>
          <th>129.107</th>
          <th>135.087</th>
          <th>136.351</th>
          <th>136.445</th>
          <th>178.722</th>
          <th>...</th>
          <th>2404.668</th>
          <th>2420.914</th>
          <th>2433.217</th>
          <th>2453.694</th>
          <th>2456.658</th>
          <th>2457.035</th>
          <th>2470.327</th>
          <th>2482.490</th>
          <th>2483.925</th>
          <th>2495.912</th>
        </tr>
        <tr>
          <th>NodeID</th>
          <th>Item</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th rowspan="5" valign="top">1</th>
          <th>t1</th>
          <td>5.548e-03</td>
          <td>4.669e-06</td>
          <td>1.816e-04</td>
          <td>-5.670e-02</td>
          <td>1.721e-04</td>
          <td>4.175e-02</td>
          <td>-8.661e-05</td>
          <td>1.341e-03</td>
          <td>1.582e-03</td>
          <td>-2.439e-01</td>
          <td>...</td>
          <td>-5.721e-02</td>
          <td>5.368e-02</td>
          <td>3.839e-02</td>
          <td>-0.133</td>
          <td>1.974e-02</td>
          <td>-0.028</td>
          <td>-0.033</td>
          <td>0.104</td>
          <td>0.069</td>
          <td>1.901e-02</td>
        </tr>
        <tr>
          <th>t2</th>
          <td>-2.133e-04</td>
          <td>5.699e-03</td>
          <td>-2.393e-02</td>
          <td>5.802e-04</td>
          <td>-1.812e-04</td>
          <td>-1.971e-04</td>
          <td>6.490e-05</td>
          <td>3.562e-02</td>
          <td>-3.164e-02</td>
          <td>-1.291e-02</td>
          <td>...</td>
          <td>3.090e-01</td>
          <td>3.746e-01</td>
          <td>5.836e-02</td>
          <td>-0.024</td>
          <td>5.890e-02</td>
          <td>-0.015</td>
          <td>-0.177</td>
          <td>0.010</td>
          <td>-0.053</td>
          <td>-1.187e-01</td>
        </tr>
        <tr>
          <th>t3</th>
          <td>8.469e-04</td>
          <td>1.512e-03</td>
          <td>-7.038e-03</td>
          <td>-8.160e-03</td>
          <td>-1.385e-03</td>
          <td>6.209e-03</td>
          <td>-1.005e-04</td>
          <td>9.286e-03</td>
          <td>-7.856e-03</td>
          <td>-3.757e-02</td>
          <td>...</td>
          <td>4.531e-02</td>
          <td>-1.270e-01</td>
          <td>-2.550e-01</td>
          <td>-0.179</td>
          <td>1.140e-03</td>
          <td>-0.042</td>
          <td>0.037</td>
          <td>-0.263</td>
          <td>-0.213</td>
          <td>1.474e-01</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>8.399e-06</td>
          <td>-2.241e-04</td>
          <td>1.035e-03</td>
          <td>-4.509e-05</td>
          <td>6.317e-05</td>
          <td>9.634e-06</td>
          <td>-2.505e-06</td>
          <td>-1.322e-03</td>
          <td>1.172e-03</td>
          <td>5.433e-04</td>
          <td>...</td>
          <td>3.061e-02</td>
          <td>9.862e-04</td>
          <td>-2.992e-02</td>
          <td>-0.035</td>
          <td>-1.128e-04</td>
          <td>-0.007</td>
          <td>0.053</td>
          <td>0.004</td>
          <td>-0.024</td>
          <td>3.404e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>2.507e-04</td>
          <td>1.228e-06</td>
          <td>8.731e-06</td>
          <td>-2.571e-03</td>
          <td>6.177e-06</td>
          <td>1.767e-03</td>
          <td>-3.825e-06</td>
          <td>5.682e-05</td>
          <td>5.614e-05</td>
          <td>-1.009e-02</td>
          <td>...</td>
          <td>1.174e-02</td>
          <td>-1.239e-03</td>
          <td>-1.026e-02</td>
          <td>-0.031</td>
          <td>4.137e-03</td>
          <td>-0.011</td>
          <td>-0.026</td>
          <td>-0.010</td>
          <td>-0.007</td>
          <td>9.107e-04</td>
        </tr>
        <tr>
          <th>...</th>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th rowspan="5" valign="top">5633</th>
          <th>t2</th>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>...</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000e+00</td>
        </tr>
        <tr>
          <th>t3</th>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>...</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000e+00</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>-3.006e-07</td>
          <td>5.476e-05</td>
          <td>-6.343e-04</td>
          <td>6.336e-06</td>
          <td>2.494e-06</td>
          <td>-2.716e-06</td>
          <td>-5.488e-07</td>
          <td>2.376e-04</td>
          <td>-2.019e-04</td>
          <td>-6.017e-05</td>
          <td>...</td>
          <td>4.288e-04</td>
          <td>3.523e-03</td>
          <td>-9.686e-04</td>
          <td>-0.007</td>
          <td>-9.864e-04</td>
          <td>-0.001</td>
          <td>0.014</td>
          <td>0.008</td>
          <td>-0.028</td>
          <td>2.645e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>-1.723e-06</td>
          <td>1.278e-06</td>
          <td>1.805e-06</td>
          <td>1.940e-04</td>
          <td>3.376e-07</td>
          <td>8.449e-06</td>
          <td>-3.548e-08</td>
          <td>4.728e-05</td>
          <td>4.650e-05</td>
          <td>-2.129e-04</td>
          <td>...</td>
          <td>-2.084e-02</td>
          <td>-1.173e-03</td>
          <td>6.237e-03</td>
          <td>-0.013</td>
          <td>1.097e-03</td>
          <td>-0.005</td>
          <td>-0.006</td>
          <td>-0.005</td>
          <td>-0.005</td>
          <td>-6.870e-03</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>7.271e-06</td>
          <td>3.394e-06</td>
          <td>2.716e-06</td>
          <td>-1.478e-04</td>
          <td>4.099e-07</td>
          <td>5.572e-05</td>
          <td>-2.954e-07</td>
          <td>1.383e-05</td>
          <td>-1.663e-05</td>
          <td>-6.515e-04</td>
          <td>...</td>
          <td>-2.914e-02</td>
          <td>-1.307e-03</td>
          <td>6.513e-03</td>
          <td>-0.014</td>
          <td>1.144e-03</td>
          <td>-0.005</td>
          <td>-0.008</td>
          <td>-0.008</td>
          <td>-0.007</td>
          <td>-8.940e-03</td>
        </tr>
      </tbody>
    </table>
    <p>32274 rows × 167 columns</p>
    </div>



Accessing the plate stress/strain
---------------------------------

Results are stored on a per element type basis.

The OP2 is the same as an F06, so CQUAD4 elements have centroidal-based
results or centroidal-based as well as the results at the 4 corner
nodes.

Be careful about what you’re accessing.

.. code:: ipython3

    # element forces/stresses/strains are by element type consistent with the F06, so...
    plate_stress = op2.cquad4_stress[1]
    print("plate_stress_obj = %s" % type(plate_stress))
    
    # the set of variables in the RealPlateStressArray
    print("plate_stress = %s\n" % plate_stress.__dict__.keys())
    
    # list of parameters that define the object (e.g. what is the nonlinear variable name
    print("data_code_keys = %s\n" % plate_stress.data_code.keys())
    
    # nonlinear variable name
    name = plate_stress.data_code['name']
    print("name = %r" % plate_stress.data_code['name'])
    
    print("list-type variables = %s" % plate_stress.data_code['data_names'])
    
    # the special loop parameter
    # for modal analysis, it's "modes"
    # for transient, it's "times"
    # or be lazy and use "_times"
    print("modes = %s" % plate_stress.modes) # name + 's'
    
    
    # extra list-type parameter for modal analysis; see data_names
    #print("mode_cycles =", plate_stress.mode_cycles)


.. parsed-literal::

    plate_stress_obj = <class 'pyNastran.op2.tables.oes_stressStrain.real.oes_plates.RealPlateStressArray'>
    plate_stress = dict_keys(['element_type', 'element_name', 'nonlinear_factor', '_times', 'result_name', 'approach_code', 'analysis_code', 'data', 'isubcase', 'ogs', 'pval_step', 'name', 'superelement_adaptivity_index', '_count', 'is_built', 'format_code', 'sort_code', 'table_code', 'title', 'subtitle', 'label', 'num_wide', 'device_code', 'table_name', 'data_frame', 'size', 'dt', 'ntimes', 'ntotal', '_ntotals', 'load_as_h5', 'h5_file', 'data_code', 'ielement', 'nelements', 'nnodes', '_encoding', '_times_dtype', 'cycle', 'data_names', 'eign', 'is_msc', 'is_nasa95', 'is_strain_flag', 'is_stress_flag', 'load_set', 'mode', 'mode2', 's_code', 'sort_bits', 'sort_method', 'stress_bits', 'subtitle_original', 'tCode', 'thermal', 'thermal_bits', 'modes', 'eigns', 'mode2s', 'cycles', 'itotal', 'itime', 'element_node', 'words'])
    
    data_code_keys = dict_keys(['_encoding', 'load_as_h5', 'size', 'is_msc', 'is_nasa95', 'table_name', 'approach_code', 'isubcase', 'table_code', 'tCode', 'sort_code', 'sort_method', 'device_code', 'analysis_code', 'sort_bits', 'element_type', 'load_set', 'format_code', 'num_wide', 's_code', 'thermal', 'nonlinear_factor', 'name', 'mode', 'eign', 'mode2', 'cycle', 'data_names', '_times_dtype', 'thermal_bits', 'element_name', 'subtitle', 'subtitle_original', 'pval_step', 'superelement_adaptivity_index', 'label', 'title', 'stress_bits', 'is_stress_flag', 'is_strain_flag', 'result_name', '_count'])
    
    name = 'mode'
    list-type variables = ['mode', 'eign', 'mode2', 'cycle']
    modes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167]
    

Similar to the BDF, we can use object_attributes/methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    #print "attributes =", object_attributes(plate_stress)
    print("methods = %s\n" % object_methods(plate_stress))
    print('methods2= %s\n' % plate_stress.object_methods())
    print("headers = %s\n" % plate_stress.get_headers())


.. parsed-literal::

    methods = ['add_new_eid_sort1', 'add_new_eid_sort2', 'add_sort1', 'add_sort2', 'apply_data_code', 'approach_code_str', 'build', 'build_dataframe', 'cast_grid_type', 'code_information', 'eid_to_element_node_index', 'export_to_hdf5', 'finalize', 'get_data_code', 'get_disp_temp', 'get_element_index', 'get_element_type', 'get_force_flux', 'get_headers', 'get_nnodes_bilinear', 'get_stats', 'get_table_code_name', 'get_unsteady_value', 'is_bilinear', 'is_magnitude_phase', 'is_sort1_new', 'is_thermal', 'object_attributes', 'object_methods', 'print_data_members', 'print_table_code', 'recast_gridtype_as_string', 'set_as_sort1', 'set_table_type', 'update_data_code', 'update_dt', 'update_t_code', 'write_f06', 'write_op2']
    
    methods2= ['add_new_eid_sort1', 'add_new_eid_sort2', 'add_sort1', 'add_sort2', 'apply_data_code', 'approach_code_str', 'build', 'build_dataframe', 'cast_grid_type', 'code_information', 'eid_to_element_node_index', 'export_to_hdf5', 'finalize', 'get_data_code', 'get_disp_temp', 'get_element_index', 'get_element_type', 'get_force_flux', 'get_headers', 'get_nnodes_bilinear', 'get_stats', 'get_table_code_name', 'get_unsteady_value', 'is_bilinear', 'is_magnitude_phase', 'is_sort1_new', 'is_thermal', 'print_data_members', 'print_table_code', 'recast_gridtype_as_string', 'set_as_sort1', 'set_table_type', 'update_data_code', 'update_dt', 'update_t_code', 'write_f06', 'write_op2']
    
    headers = ['fiber_distance', 'oxx', 'oyy', 'txy', 'angle', 'omax', 'omin', 'von_mises']
    
    

Number of Nodes on a CQUAD4
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  For CENT, there is 1 centroidal stress at two locations
-  For BILIN, there are 5 stresses at two locations (4 nodes +
   centroidal)
-  node_id=0 indicates a centroidal quantity
-  CTRIA3s are always centroidal

What sets this?
^^^^^^^^^^^^^^^

::

   STRESS(real, sort1, BILIN) = ALL   # centroid + 4 corner nodes
   STRESS(real, sort1, CENT) = ALL    # centroid

   STRAIN(real, sort1, BILIN) = ALL   # centroid + 4 corner nodes
   STRAIN(real, sort1, CENT) = ALL    # centroid

How do we know if we’re bilinear?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   print("is_bilinear = %s\n" % plate_stress.is_bilinear())

What locations are chosen?
^^^^^^^^^^^^^^^^^^^^^^^^^^

That depends on fiber distance/fiber curvature… - fiber_curvature - mean
stress (:math:`\sigma_{alt}`) & slope (:math:`\sigma_{mean}`)

::

   $$ \sigma_{top} = \sigma_{alt} + \frac{t}{2} \sigma_{mean}$$

   $$ \sigma_{btm} = \sigma_{alt} + \frac{t}{2} \sigma_{mean}$$

-  fiber_distance - upper and lower surface stress (o_top; o_btm)
-  If you have stress, fiber_distance is always returned regardless of
   your option.

.. _what-sets-this-1:

What sets this?
^^^^^^^^^^^^^^^

::

   STRAIN(real, sort1, FIBER) = ALL   # fiber distance/default
   STRAIN(real, sort1, STRCUR) = ALL  # strain curvature

How do we know if we’re using fiber_distance?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   print("is_fiber_distance = %s" % plate_stress.is_fiber_distance)

Accessing results
-----------------

Note that this is intentionally done iinefficiently to access specific entries in order to explain the data structure.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    # element forces/stresses/strains are by element type consistent
    # with the F06, so...
    
    def abs_max_min(vals):
        absvals = list(abs(vals))
        maxval = max(absvals)
        i = absvals.index(maxval)
        return vals[i]
    
    #-----------------------------
    # again, we have linear quads, so two locations per element
    print("element_node[:10, :] =\n%s..." % plate_stress.element_node[:10, :])
    
    # lets get the stress for the first 3 CQUAD4 elements
    eids = plate_stress.element_node[:, 0]
    ueids = np.unique(eids)
    print('ueids = %s' % ueids[:3])
    
    # get the first index of the first 5 elements
    ieids = np.searchsorted(eids, ueids[:3])
    print('ieids = %s' % ieids)
    
    # the easy way to slice data for linear plates
    ieids5 = np.vstack([ieids, ieids + 1]).ravel()
    ieids5.sort()
    
    print('verify5:\n%s' % ieids5)
    
    #-----------------------------
    itime = 0 # static analysis / mode 1
    if plate_stress.is_von_mises:  # True
        ovm = plate_stress.data[itime, :, 7]
        print('we have von mises data; ovm=%s\n' % ovm)
    else:
        omax_shear = plate_stress.data[itime, :, 7]
        print('we have max shear data; omax_shear=%s\n' % omax_shear)
    
    
    print("[layer1, layer2, ...] = %s" % ovm[ieids5])
    
    ieid1000 = np.where(eids == 1000)[0]
    print('ieid1000 = %s' % ieid1000)
    ovm_mode6_eid1000 = ovm[ieid1000]
    print("ovm_mode6_eid1000 = %s -> %s" % (ovm_mode6_eid1000, abs_max_min(ovm_mode6_eid1000)))


.. parsed-literal::

    element_node[:10, :] =
    [[1 0]
     [1 0]
     [2 0]
     [2 0]
     [3 0]
     [3 0]
     [4 0]
     [4 0]
     [5 0]
     [5 0]]...
    ueids = [1 2 3]
    ieids = [0 2 4]
    verify5:
    [0 1 2 3 4 5]
    we have von mises data; ovm=[54.223  5.041 13.143 ...  2.34   6.146  7.368]
    
    [layer1, layer2, ...] = [54.223  5.041 13.143 21.223 78.546 17.91 ]
    ieid1000 = [1998 1999]
    ovm_mode6_eid1000 = [90.618 94.093] -> 94.09257
    

.. code:: ipython3

    # see the difference between "transient"/"modal"/"frequency"-style results
    # and "nodal"/"elemental"-style results
    # just change imode
    
    imode = 5  # mode 6; could just as easily be dt
    iele = 10  # element 10
    ilayer = 1
    
    ieid10 = np.where(eids == iele)[0][ilayer]
    print('ieid10 = %s' % ieid10)
    print(plate_stress.element_node[ieid10, :])
    
    
    # headers = [u'fiber_distance', u'oxx', u'oyy', u'txy', u'angle', u'omax', u'omin', u'von_mises']
    print("ps.modes = %s" % plate_stress.modes[imode])
    print("ps.cycles = %s" % plate_stress.cycles[imode])
    print("oxx = %s" % plate_stress.data[imode, ieid10, 1])
    print("oyy = %s" % plate_stress.data[imode, ieid10, 2])
    print("txy = %s" % plate_stress.data[imode, ieid10, 3])
    print("omax = %s" % plate_stress.data[imode, ieid10, 5])
    print("omin = %s" % plate_stress.data[imode, ieid10, 6])
    print("ovm/max_shear = %s" % plate_stress.data[imode, ieid10, 7])
    
    if plate_stress.is_fiber_distance:
        print("fiber_distance = %s" % plate_stress.data[imode, ieid10, 0])
    else:
        print("curvature = %s" % plate_stress.data[imode, ieid10, 0])


.. parsed-literal::

    ieid10 = 19
    [10  0]
    ps.modes = 6
    ps.cycles = 20.548039949617547
    oxx = 18.868874
    oyy = 20.159315
    txy = 8.309595
    omax = 27.848701
    omin = 11.179487
    ovm/max_shear = 24.273378
    fiber_distance = -0.4
    

.. code:: ipython3

    from pyNastran.bdf.bdf import read_bdf
    from pyNastran.bdf.mesh_utils.mass_properties import mass_properties
    bdf_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_4pt.dat'))
    model = read_bdf(bdf_filename, debug=False)
    mass, cg, I = mass_properties(model)



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=1 midsurface: z1=0.4 z2=-0.4 t=0.036 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=2 midsurface: z1=0.4 z2=-0.4 t=0.054 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=3 midsurface: z1=0.4 z2=-0.4 t=0.018 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=7 midsurface: z1=0.418 z2=-0.418 t=0.036 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=34 midsurface: z1=0.194 z2=-0.194 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=38 midsurface: z1=0.284 z2=-0.284 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=46 midsurface: z1=0.199 z2=-0.199 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=37 midsurface: z1=0.309 z2=-0.309 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>


Let’s print out the actual mass properties from the OP2 and get the same result as the F06
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We need ``PARAM,POSTEXT,YES`` in out BDF to get the Grid Point Weight
Table

.. code:: ipython3

    gpw = op2.grid_point_weight
    #print(gpw)
    if gpw:
        gpwi = gpw['']
        #print(gpw.object_attributes())
        #print(gpwi)
        gpwi.object_methods()
        #print(gpwi.object_attributes())
        #gpw.write_f06?
        print(gpwi.get_stats())
        print('M0 = ', gpwi.MO)
        print('S  = ', gpwi.S)
        print('mass = ', gpwi.mass)
        print('cg = ', gpwi.cg)
        print('IS = ', gpwi.IS)
        print('IQ = ', gpwi.IQ)
        print('Q  = ', gpwi.Q)


.. parsed-literal::

    GridPointWeight['']: ref_point=0 mass=1.7746; [reference_point, M0, S, mass, cg, IS, IQ, Q]
    
    M0 =  [[   1.775   -0.      -0.      -0.     -32.773    4.491]
     [  -0.       1.775    0.      32.773   -0.      -0.06 ]
     [  -0.       0.       1.775   -4.491    0.06    -0.   ]
     [  -0.      32.773   -4.491 1322.291    1.415   -1.251]
     [ -32.773   -0.       0.06     1.415 1227.076 -218.771]
     [   4.491   -0.06    -0.      -1.251 -218.771  427.228]]
    S  =  [[1. 0. 0.]
     [0. 1. 0.]
     [0. 0. 1.]]
    mass =  [1.775 1.775 1.775]
    cg =  [[ -0.     -2.531 -18.468]
     [ -0.034  -0.    -18.468]
     [ -0.034  -2.531  -0.   ]]
    IS =  [[705.69   -1.567   0.141]
     [ -1.567 621.837 135.836]
     [  0.141 135.836 415.862]]
    IQ =  [689.184 348.385 705.821]
    Q  =  [[ 0.088  0.002  0.996]
     [-0.892 -0.445  0.08 ]
     [ 0.443 -0.896 -0.038]]
    

We can also write the full ``F06``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    import getpass
    name = getpass.getuser()
    os.chdir(os.path.join(r'C:\Users', name, 'Desktop'))
    
    # write the F06 with Real/Imaginary or Magnitude/Phase
    # only matters for complex results
    op2.write_f06('isat.f06', is_mag_phase=False)
    
    !head -n 40 isat.f06


.. parsed-literal::

                               O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R
    0                                                     REFERENCE POINT =        0
                                                                    M O
                          *  1.774604E+00 -2.303930E-19 -5.452775E-20 -1.068744E-17 -3.277277E+01  4.490806E+00 *
                          * -2.303930E-19  1.774604E+00  1.829591E-19  3.277277E+01 -1.803164E-17 -6.007288E-02 *
                          * -5.452775E-20  1.829591E-19  1.774604E+00 -4.490806E+00  6.007288E-02 -1.164670E-19 *
                          * -1.068744E-17  3.277277E+01 -4.490806E+00  1.322291E+03  1.414705E+00 -1.250593E+00 *
                          * -3.277277E+01 -1.803164E-17  6.007288E-02  1.414705E+00  1.227076E+03 -2.187709E+02 *
                          *  4.490806E+00 -6.007288E-02 -1.164670E-19 -1.250593E+00 -2.187709E+02  4.272284E+02 *
                                                                     S
                                               *  1.000000E+00  0.000000E+00  0.000000E+00 *
                                               *  0.000000E+00  1.000000E+00  0.000000E+00 *
                                               *  0.000000E+00  0.000000E+00  1.000000E+00 *
                                   DIRECTION
                              MASS AXIS SYSTEM (S)     MASS              X-C.G.        Y-C.G.        Z-C.G.
                                      X            1.774604E+00     -6.022436E-18 -2.530596E+00 -1.846765E+01
                                      Y            1.774604E+00     -3.385143E-02 -1.016094E-17 -1.846765E+01
                                      Z            1.774604E+00     -3.385143E-02 -2.530596E+00 -6.562987E-20
                                                                    I(S)
                                               *  7.056902E+02 -1.566725E+00  1.411882E-01 *
                                               * -1.566725E+00  6.218375E+02  1.358363E+02 *
                                               *  1.411882E-01  1.358363E+02  4.158619E+02 *
                                                                    I(Q)
                                               *  6.891835E+02                             *
                                               *                3.483848E+02               *
                                               *                              7.058213E+02 *
                                                                     Q
                                               *  8.846130E-02  1.596867E-03  9.960783E-01 *
                                               * -8.920127E-01 -4.448868E-01  7.993249E-02 *
                                               *  4.432697E-01 -8.955854E-01 -3.793084E-02 *
    
    1    ISAT_SM_LAUNCH_4PT MODES TO 400 HZ                                    JULY      16, 2020  pyNastran v1.4.0+dev.cc0dbf554  PAGE     1
    
    1    ISAT_SM_LAUNCH_4PT MODES TO 400 HZ                                    JULY      16, 2020  pyNastran v1.4.0+dev.cc0dbf554  PAGE     2
         DEFAULT                                                                                                                        
    
                                                  R E A L   E I G E N V A L U E S
                                                 ISAT_SM_LAUNCH_4PT MODES TO 400 HZ
       MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED
        NO.       ORDER                                                                       MASS              STIFFNESS
    

.. code:: ipython3

    #from IPython.display import display, Math, Latex

The mass results are different as pyNastran’s mass assumes point masses

.. math:: m_{plates} = A (\rho t + nsm)

.. math:: m_{solid} = V \rho

.. math:: m_{bars} = L (\rho A + nsm)

.. math:: I = m r^2

The larger your model is and the further from the origin, the more
accurate the result. For some applications (e.g. a weight breakdown),
this is probably be fine.

.. code:: ipython3

    print('cg =\n%s' % gpw[''].cg)
    print('cg = %s' % cg)


.. parsed-literal::

    cg =
    [[ -0.     -2.531 -18.468]
     [ -0.034  -0.    -18.468]
     [ -0.034  -2.531  -0.   ]]
    cg = [ -0.034  -2.531 -18.468]
    

It’s not like Nastran is perfect either.
----------------------------------------

Limitations
~~~~~~~~~~~

1. You cannot do weight statements in Nastran by
   component/property/material.

2. Everything is always summmed up (e.g. you can have different geometry
   in Subcase 2 and MPCs connecting physical geometry, with other parts
   flying off into space).

These are things that pyNastran ``can`` do.

.. code:: ipython3

    from pyNastran.bdf.bdf import read_bdf
    bdf_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_4pt.dat'))
    model = read_bdf(bdf_filename, debug=False)



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=1 midsurface: z1=0.4 z2=-0.4 t=0.036 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=2 midsurface: z1=0.4 z2=-0.4 t=0.054 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=3 midsurface: z1=0.4 z2=-0.4 t=0.018 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=7 midsurface: z1=0.418 z2=-0.418 t=0.036 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=34 midsurface: z1=0.194 z2=-0.194 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=38 midsurface: z1=0.284 z2=-0.284 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=46 midsurface: z1=0.199 z2=-0.199 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>



.. raw:: html

    <text style=color:orange>WARNING: shell.py:2544                PSHELL pid=37 midsurface: z1=0.309 z2=-0.309 t=0.0186 not in range of -1.5t < zi < 1.5t
    </text>


Weight Statement
~~~~~~~~~~~~~~~~

Let’s get the breakdown by property ID

.. code:: ipython3

    #help(model.mass_properties)
    
    pid_to_eids_map = model.get_element_ids_dict_with_pids()
    #print(pid_to_eids_map.keys())
    print('pid, mass, cg, [ixx, iyy, izz, ixy, ixz, iyz]')
    for pid, eids in sorted(pid_to_eids_map.items()):
        mass, cg, inertia = mass_properties(model, element_ids=eids, mass_ids=[], reference_point=[0., 0., 0.])
        print('%-6s %-.6f %-38s %s' % (pid, mass, cg, inertia))
    
    mass_ids = list(model.masses.keys())
    mass, cg, inertia = mass_properties(model, element_ids=[], mass_ids=mass_ids, reference_point=[0., 0., 0.])
    print('%-6s %-.6f %-38s %s' % ('mass', mass, cg, inertia))


.. parsed-literal::

    pid, mass, cg, [ixx, iyy, izz, ixy, ixz, iyz]
    1      0.027278 [  0.   0. -20.]                       [3.699 6.553 4.384 0.    0.    0.   ]
    2      0.047993 [ -0.   0. -20.]                       [18.033 18.033 12.454 -0.    -0.     0.   ]
    3      0.020998 [  0.  -0. -20.]                       [5.881 3.907 5.27  0.    0.    0.   ]
    4      0.012216 [  0.043   0.438 -19.702]              [2.346 3.23  2.019 0.01  0.005 0.052]
    5      0.330158 [  0.    2.2 -20. ]                    [63.317 28.366 41.752  0.     0.     0.   ]
    7      0.027813 [  0.  -0. -20.]                       [ 8.141  8.141  9.438 -0.    -0.    -0.   ]
    8      0.081584 [  0.   0. -20.]                       [15.087 15.087 30.174 -0.     0.     0.   ]
    9      0.077642 [  0.   0. -20.]                       [17.017 17.017 18.911 -0.    -0.     0.   ]
    10     0.000236 [  0.  -0. -20.]                       [ 0.035  0.035  0.057 -0.    -0.    -0.   ]
    11     0.041700 [ -1.025  23.773 -12.016]              [ 0.666  0.988  0.348 -0.037  0.263 -0.056]
    12     0.000457 [ 0.    -5.92  20.506]                 [0.013 0.013 0.    0.    0.    0.001]
    13     0.003885 [ 0.    -6.949  9.892]                 [ 0.002  0.     0.002  0.     0.    -0.   ]
    14     0.000353 [-0.     0.    14.391]                 [ 0.003  0.012  0.009  0.    -0.     0.   ]
    15     0.003626 [0.    0.    7.867]                    [ 0.     0.092  0.091  0.    -0.     0.   ]
    16     0.000000 [0. 0. 0.]                             [0. 0. 0. 0. 0. 0.]
    19     0.017749 [ -0.23    6.021 -35.642]              [ 1.77   4.395  5.817 -0.053  0.005 -0.145]
    20     0.163082 [  0.      0.    -18.545]              [ 9.01 34.77 25.76  0.    0.    0.  ]
    21     0.003625 [ -0.  -0. -20.]                       [ 0.728  0.728  1.41  -0.    -0.    -0.   ]
    22     0.000000 [0. 0. 0.]                             [0. 0. 0. 0. 0. 0.]
    23     0.000000 [0. 0. 0.]                             [0. 0. 0. 0. 0. 0.]
    33     0.001346 [-0.    -2.175  0.369]                 [ 0.077  0.085  0.162  0.    -0.    -0.001]
    34     0.003561 [-0.    -0.    14.833]                 [ 0.271  0.271  0.067 -0.    -0.    -0.   ]
    35     0.000000 [0. 0. 0.]                             [0. 0. 0. 0. 0. 0.]
    36     0.007197 [  0.      0.    -14.783]              [ 0.835  3.605  3.02   0.     0.    -0.   ]
    37     0.094566 [ -0.      0.    -19.499]              [ 8.975 52.72  46.49   0.    -0.     0.   ]
    38     0.007602 [ 0.    -9.329 27.311]                 [0.681 1.214 0.562 0.    0.    0.07 ]
    39     0.002433 [ 0.    -8.954  4.04 ]                 [ 0.022  0.045  0.024  0.    -0.    -0.003]
    41     0.000735 [-0.    -0.     2.193]                 [ 0.009  0.057  0.062  0.     0.    -0.   ]
    42     0.008854 [ -1.554  20.121 -19.007]              [ 1.166  1.459  0.293 -0.001  0.056  0.   ]
    43     0.012241 [  0.      0.    -19.499]              [2.228 7.588 6.32  0.    0.    0.   ]
    46     0.003671 [ 0.    0.   15.28]                    [ 0.178  0.348  0.335  0.    -0.     0.   ]
    60     0.000000 [0. 0. 0.]                             [0. 0. 0. 0. 0. 0.]
    61     0.000000 [0. 0. 0.]                             [0. 0. 0. 0. 0. 0.]
    mass   0.772000 [  0.     -8.256 -18.238]              [392.813 338.699 118.704  -0.     -0.    138.698]
    
