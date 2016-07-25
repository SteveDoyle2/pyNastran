
OP2 Demo
========

The iPython notebook for this demo can be found in: -
docs:raw-latex:`\quick`\_start:raw-latex:`\demo`:raw-latex:`\op`2\_demo.ipynb
-
https://github.com/SteveDoyle2/pyNastran/tree/master/docs/quick\_start/demo/op2\_demo.ipynb

Why use the OP2? Why not use the F06/PCH file?
----------------------------------------------

Most people are comfortable with the F06. However, it's: - Ironically, a
lot harder to parse. The OP2 is very structured. - Much, much, much
slower. We can read entire blocks of arrays with a single call. The data
is already typed. - Much, much more memory inefficient because we aren't
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
Assuming it's on your path (it'll be in Python27:raw-latex:`\Scripts `or
something similar), you can run:

::

    >>> test_op2 -f solid_bending.op2

The ``-f`` tells us to print out ``solid_bending.test_op2.f06``, which
can be compared to your F06 for a small file to build confidence in the
reader. It's also useful when you want an F06 of your model without
rerunning Nastran just to see what's in it.

If you have a large model, you can make ``test_op2`` run much, much
faster. The ``-c`` flag disables double-reading of the OP2. By default,
``test_op2`` uses two different read methods (the old method and new
method) to ensure that results are read in properly. When running the
code, this is turned off, but is turned on for ``test_op2``.

::

    >>> test_op2 -fc solid_bending.op2

Import the packages
-------------------

.. code:: python

    import os
    import copy
    import numpy as np
    
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    
    from pyNastran.utils import print_bad_path
    from pyNastran.op2.op2 import read_op2
    from pyNastran.utils import object_methods, object_attributes
    
    import pandas as pd

Sets default precision of real numbers for pandas output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    pd.set_option('precision', 3)
    np.set_printoptions(precision=3, threshold=20)

As with the BDF, we can use the long form and the short form. However,
the long form for the ``OP2`` doesn't really add anything. So, let's
just use the short form.

Besides massive speed improvements in the OP2 relative to v0.7, this
version adds ``pandas`` dataframe support.

.. code:: python

    #op2_filename = r'D:\work\pynastran_0.8.0\models\iSat\ISat_Launch_Sm_Rgd.op2'
    #op2_filename = r'D:\work\pynastran_0.8.0\models\iSat\ISat_Launch_Sm_4pt.op2'
    op2_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_4pt.op2'))
    assert os.path.exists(op2_filename), print_bad_path(op2_filename)
    
    # define the input file with a file path
    op2 = read_op2(op2_filename, build_dataframe=True, debug=False)


.. parsed-literal::

    INFO:      fname=op2_scalar.py             lineNo=1176   op2_filename = 'f:\\work\\pynastran\\pynastran\\master3\\models\\iSat\\ISat_Launch_Sm_4pt.op2'
    

OP2 Introspection
-----------------

The ``get_op2_stats()`` function lets you quickly understand what in an
op2.

.. code:: python

    print(op2.get_op2_stats())


.. parsed-literal::

    GridPointWeight:  reference_point=0
      mass=[    1.7746     1.7746     1.7746]
      cg  =[-1.02627e-17   -2.53061   -18.4676]
           [-0.0338509 -1.63349e-17   -18.4676]
           [-0.0338509   -2.53061 -2.18142e-19]
    
      IS  =[    705.69   -1.56671   0.141187]
           [  -1.56671    621.838    135.836]
           [  0.141187    135.836    415.861]
    
      IQ  =[   689.183                      ]
           [              348.384           ]
           [                         705.821]
    
      Q  = [ 0.0884636 0.00159685   0.996078]
           [ -0.892013  -0.444886  0.0799345]
           [  0.443269  -0.895586 -0.0379318]
    eigenvectors[1]
      isubcase = 1
      type=RealEigenvectorArray ntimes=167 nnodes=5379
      data: [t1, t2, t3, r1, r2, r3] shape=[167, 5379, 6] dtype=float32
      gridTypes
      sort1
      modes = [  1   2   3 ..., 165 166 167]
      eigrs = [  2.758e+03   3.568e+03   9.686e+03 ...,   6.163e+06   6.170e+06
       6.230e+06]
      mode_cycles = [1112674317 1114566525 1120195719 ..., 1159407589 1159413465 1159462558]
    
    cbar_force[1]
      type=RealCBarForceArray ntimes=167 nelements=827
      data: [ntimes, nnodes, 8] where 8=[bending_moment_a1, bending_moment_a2, bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque]
      data.shape = (167, 827, 8)
      element name: CBAR-34
      sort1
      modes = [  1   2   3 ..., 165 166 167]
      eigrs = [  2.758e+03   3.568e+03   9.686e+03 ...,   6.163e+06   6.170e+06
       6.230e+06]
      cycles = [   8.358    9.507   15.663 ...,  395.101  395.329  397.237]
    
    ctria3_stress[1]
      type=RealPlateStressArray ntimes=167 nelements=32 nnodes_per_element=1 nlayers=2 ntotal=64
      data: [ntimes, ntotal, 8] where 8=[fiber_distance, oxx, oyy, txy, angle, omax, omin, von_mises]
      data.shape=(167L, 64L, 8L)
      element type: CTRIA3
      s_code: 1
      sort1
      modes = [  1   2   3 ..., 165 166 167]
      eigrs = [  2.758e+03   3.568e+03   9.686e+03 ...,   6.163e+06   6.170e+06
       6.230e+06]
      mode2s = [         2          3          4 ...,        166        167 1159462558]
      cycles = [  2.803e-45   4.204e-45   5.605e-45 ...,   2.326e-43   2.340e-43
       2.496e+03]
    
    cquad4_stress[1]
      type=RealPlateStressArray ntimes=167 nelements=4580 nnodes_per_element=1 nlayers=2 ntotal=9160
      data: [ntimes, ntotal, 8] where 8=[fiber_distance, oxx, oyy, txy, angle, omax, omin, von_mises]
      data.shape=(167L, 9160L, 8L)
      element type: CQUAD4
      s_code: 1
      sort1
      modes = [  1   2   3 ..., 165 166 167]
      eigrs = [  2.758e+03   3.568e+03   9.686e+03 ...,   6.163e+06   6.170e+06
       6.230e+06]
      mode2s = [         2          3          4 ...,        166        167 1159462558]
      cycles = [  2.803e-45   4.204e-45   5.605e-45 ...,   2.326e-43   2.340e-43
       2.496e+03]
    
    eigenvalues[ISAT_SM_LAUNCH_4PT MODES TO 400 HZ]
      type=RealEigenvalues neigenvalues=167
      title, extraction_order, eigenvalues, radians, cycles, generalized_mass, generalized_stiffness
    
    
    

If that's too long...
~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    print(op2.get_op2_stats(short=True))


.. parsed-literal::

    GridPointWeight: ref_point=0 mass=1.7746; [reference_point, M0, S, mass, cg, IS, IQ, Q]
    eigenvectors[1]
    cbar_force[1]
    ctria3_stress[1]
    cquad4_stress[1]
    eigenvalues[u'ISAT_SM_LAUNCH_4PT MODES TO 400 HZ']
    
    

Acccessing the Eigenvectors object
----------------------------------

Eigenvectors are the simplest object. They use the same class as for
displacements, velocity, acceleration, SPC Forces, MPC Forces, Applied
Loads, etc. These are all node-based tables with TX, TY, TZ, RX, RY, RZ.
Results are in the analysis coordinate frame (CD), which is defined by
the GRID card.

Numpy-based Approach
~~~~~~~~~~~~~~~~~~~~

We'll first show off the standard ``numpy`` based results on a transient
case. Static results are the same, except that you'll always use the 0th
index for the "time" index.

The tutorial is intetionally just accessing the objects in a very clear,
though inefficient way. The OP2 objects can take full advantage of the
numpy operations.

.. code:: python

    # what modes did we analyze:  1 to 167
    print("loadcases = %s" % op2.eigenvectors.keys())
    
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
    
    times = [   1.    2.    3. ...,  165.  166.  167.]
    
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
    translation mode2_node10 = [  1.696e-05   7.937e-03   1.510e-03]
    rotations mode2_node10 = [ -2.241e-04   1.228e-06  -1.187e-06]
    

Pandas-based Approach
~~~~~~~~~~~~~~~~~~~~~

If you like pandas, you can access all the OP2 objects, which is very
useful within the Jupyter Notebook. Different objects will look
differently, but you can change the layout.

If you're trying to learn pandas, there are many tutorials online, such
as: http://pandas.pydata.org/pandas-docs/stable/10min.html

or a very long, but good video:

.. code:: python

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
            



.. code:: python

    # get subcase 1
    eig1 = op2.eigenvectors[1]
    
    eig1.data_frame




.. raw:: html

    <div>
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
          <th>EigenvalueReal</th>
          <th>2757.896</th>
          <th>3568.136</th>
          <th>9685.530</th>
          <th>16154.883</th>
          <th>16278.047</th>
          <th>16668.713</th>
          <th>18248.492</th>
          <th>18591.637</th>
          <th>18617.254</th>
          <th>31930.465</th>
          <th>...</th>
          <th>5782436.500</th>
          <th>5860846.500</th>
          <th>5920603.000</th>
          <th>6020617.500</th>
          <th>6035178.000</th>
          <th>6037030.000</th>
          <th>6102521.500</th>
          <th>6162773.500</th>
          <th>6169898.500</th>
          <th>6229584.500</th>
        </tr>
        <tr>
          <th></th>
          <th>Freq</th>
          <th>8.358</th>
          <th>9.507</th>
          <th>15.663</th>
          <th>20.229</th>
          <th>20.306</th>
          <th>20.548</th>
          <th>21.500</th>
          <th>21.701</th>
          <th>21.716</th>
          <th>28.440</th>
          <th>...</th>
          <th>382.715</th>
          <th>385.301</th>
          <th>387.260</th>
          <th>390.518</th>
          <th>390.990</th>
          <th>391.050</th>
          <th>393.165</th>
          <th>395.101</th>
          <th>395.329</th>
          <th>397.237</th>
        </tr>
        <tr>
          <th></th>
          <th>Radians</th>
          <th>52.516</th>
          <th>59.734</th>
          <th>98.415</th>
          <th>127.102</th>
          <th>127.585</th>
          <th>129.107</th>
          <th>135.087</th>
          <th>136.351</th>
          <th>136.445</th>
          <th>178.691</th>
          <th>...</th>
          <th>2404.670</th>
          <th>2420.919</th>
          <th>2433.229</th>
          <th>2453.695</th>
          <th>2456.660</th>
          <th>2457.037</th>
          <th>2470.328</th>
          <th>2482.493</th>
          <th>2483.928</th>
          <th>2495.914</th>
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
          <th rowspan="6" valign="top">1</th>
          <th>t1</th>
          <td>5.548e-03</td>
          <td>4.671e-06</td>
          <td>-1.818e-04</td>
          <td>-5.670e-02</td>
          <td>1.722e-04</td>
          <td>-4.175e-02</td>
          <td>8.632e-05</td>
          <td>-1.341e-03</td>
          <td>1.582e-03</td>
          <td>-2.428e-01</td>
          <td>...</td>
          <td>5.723e-02</td>
          <td>-5.369e-02</td>
          <td>-3.838e-02</td>
          <td>-1.326e-01</td>
          <td>-1.973e-02</td>
          <td>0.028</td>
          <td>0.033</td>
          <td>-0.104</td>
          <td>-6.919e-02</td>
          <td>-1.904e-02</td>
        </tr>
        <tr>
          <th>t2</th>
          <td>-2.133e-04</td>
          <td>5.699e-03</td>
          <td>2.392e-02</td>
          <td>5.801e-04</td>
          <td>-1.812e-04</td>
          <td>1.971e-04</td>
          <td>-6.526e-05</td>
          <td>-3.563e-02</td>
          <td>-3.164e-02</td>
          <td>-1.292e-02</td>
          <td>...</td>
          <td>-3.090e-01</td>
          <td>-3.746e-01</td>
          <td>-5.840e-02</td>
          <td>-2.385e-02</td>
          <td>-5.889e-02</td>
          <td>0.015</td>
          <td>0.177</td>
          <td>-0.010</td>
          <td>5.252e-02</td>
          <td>1.187e-01</td>
        </tr>
        <tr>
          <th>t3</th>
          <td>8.469e-04</td>
          <td>1.512e-03</td>
          <td>7.038e-03</td>
          <td>-8.160e-03</td>
          <td>-1.385e-03</td>
          <td>-6.209e-03</td>
          <td>1.004e-04</td>
          <td>-9.286e-03</td>
          <td>-7.856e-03</td>
          <td>-3.743e-02</td>
          <td>...</td>
          <td>-4.535e-02</td>
          <td>1.271e-01</td>
          <td>2.550e-01</td>
          <td>-1.792e-01</td>
          <td>-1.136e-03</td>
          <td>0.042</td>
          <td>-0.037</td>
          <td>0.263</td>
          <td>2.141e-01</td>
          <td>-1.473e-01</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>8.399e-06</td>
          <td>-2.241e-04</td>
          <td>-1.035e-03</td>
          <td>-4.509e-05</td>
          <td>6.317e-05</td>
          <td>-9.634e-06</td>
          <td>2.518e-06</td>
          <td>1.322e-03</td>
          <td>1.172e-03</td>
          <td>5.440e-04</td>
          <td>...</td>
          <td>-3.061e-02</td>
          <td>-9.829e-04</td>
          <td>2.993e-02</td>
          <td>-3.527e-02</td>
          <td>1.148e-04</td>
          <td>0.007</td>
          <td>-0.053</td>
          <td>-0.004</td>
          <td>2.357e-02</td>
          <td>-3.403e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>2.507e-04</td>
          <td>1.228e-06</td>
          <td>-8.742e-06</td>
          <td>-2.571e-03</td>
          <td>6.180e-06</td>
          <td>-1.767e-03</td>
          <td>3.812e-06</td>
          <td>-5.683e-05</td>
          <td>5.614e-05</td>
          <td>-1.004e-02</td>
          <td>...</td>
          <td>-1.174e-02</td>
          <td>1.241e-03</td>
          <td>1.025e-02</td>
          <td>-3.112e-02</td>
          <td>-4.135e-03</td>
          <td>0.011</td>
          <td>0.026</td>
          <td>0.009</td>
          <td>7.311e-03</td>
          <td>-9.083e-04</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>5.261e-05</td>
          <td>-1.187e-06</td>
          <td>-1.986e-04</td>
          <td>-1.310e-04</td>
          <td>-2.860e-05</td>
          <td>-4.676e-05</td>
          <td>-1.092e-07</td>
          <td>-1.774e-04</td>
          <td>1.806e-04</td>
          <td>1.011e-03</td>
          <td>...</td>
          <td>4.109e-05</td>
          <td>2.184e-02</td>
          <td>2.495e-03</td>
          <td>8.832e-02</td>
          <td>1.660e-02</td>
          <td>-0.030</td>
          <td>-0.100</td>
          <td>0.022</td>
          <td>2.547e-02</td>
          <td>-5.581e-03</td>
        </tr>
        <tr>
          <th rowspan="6" valign="top">2</th>
          <th>t1</th>
          <td>5.548e-03</td>
          <td>4.671e-06</td>
          <td>-1.818e-04</td>
          <td>-5.670e-02</td>
          <td>1.722e-04</td>
          <td>-4.175e-02</td>
          <td>8.632e-05</td>
          <td>-1.341e-03</td>
          <td>1.582e-03</td>
          <td>-2.428e-01</td>
          <td>...</td>
          <td>5.723e-02</td>
          <td>-5.369e-02</td>
          <td>-3.838e-02</td>
          <td>-1.326e-01</td>
          <td>-1.973e-02</td>
          <td>0.028</td>
          <td>0.033</td>
          <td>-0.104</td>
          <td>-6.919e-02</td>
          <td>-1.904e-02</td>
        </tr>
        <tr>
          <th>t2</th>
          <td>-1.081e-04</td>
          <td>5.696e-03</td>
          <td>2.353e-02</td>
          <td>3.180e-04</td>
          <td>-2.384e-04</td>
          <td>1.036e-04</td>
          <td>-6.548e-05</td>
          <td>-3.598e-02</td>
          <td>-3.128e-02</td>
          <td>-1.090e-02</td>
          <td>...</td>
          <td>-3.090e-01</td>
          <td>-3.309e-01</td>
          <td>-5.341e-02</td>
          <td>1.528e-01</td>
          <td>-2.568e-02</td>
          <td>-0.045</td>
          <td>-0.022</td>
          <td>0.034</td>
          <td>1.035e-01</td>
          <td>1.075e-01</td>
        </tr>
        <tr>
          <th>t3</th>
          <td>3.455e-04</td>
          <td>1.510e-03</td>
          <td>7.055e-03</td>
          <td>-3.018e-03</td>
          <td>-1.398e-03</td>
          <td>-2.676e-03</td>
          <td>9.274e-05</td>
          <td>-9.172e-03</td>
          <td>-7.968e-03</td>
          <td>-1.735e-02</td>
          <td>...</td>
          <td>-2.187e-02</td>
          <td>1.246e-01</td>
          <td>2.345e-01</td>
          <td>-1.170e-01</td>
          <td>7.135e-03</td>
          <td>0.020</td>
          <td>-0.090</td>
          <td>0.244</td>
          <td>1.995e-01</td>
          <td>-1.454e-01</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>8.399e-06</td>
          <td>-2.241e-04</td>
          <td>-1.035e-03</td>
          <td>-4.509e-05</td>
          <td>6.317e-05</td>
          <td>-9.634e-06</td>
          <td>2.518e-06</td>
          <td>1.322e-03</td>
          <td>1.172e-03</td>
          <td>5.440e-04</td>
          <td>...</td>
          <td>-3.061e-02</td>
          <td>-9.829e-04</td>
          <td>2.993e-02</td>
          <td>-3.527e-02</td>
          <td>1.148e-04</td>
          <td>0.007</td>
          <td>-0.053</td>
          <td>-0.004</td>
          <td>2.357e-02</td>
          <td>-3.403e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>2.507e-04</td>
          <td>1.228e-06</td>
          <td>-8.742e-06</td>
          <td>-2.571e-03</td>
          <td>6.180e-06</td>
          <td>-1.767e-03</td>
          <td>3.812e-06</td>
          <td>-5.683e-05</td>
          <td>5.614e-05</td>
          <td>-1.004e-02</td>
          <td>...</td>
          <td>-1.174e-02</td>
          <td>1.241e-03</td>
          <td>1.025e-02</td>
          <td>-3.112e-02</td>
          <td>-4.135e-03</td>
          <td>0.011</td>
          <td>0.026</td>
          <td>0.009</td>
          <td>7.311e-03</td>
          <td>-9.083e-04</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>5.261e-05</td>
          <td>-1.187e-06</td>
          <td>-1.986e-04</td>
          <td>-1.310e-04</td>
          <td>-2.860e-05</td>
          <td>-4.676e-05</td>
          <td>-1.092e-07</td>
          <td>-1.774e-04</td>
          <td>1.806e-04</td>
          <td>1.011e-03</td>
          <td>...</td>
          <td>4.109e-05</td>
          <td>2.184e-02</td>
          <td>2.495e-03</td>
          <td>8.832e-02</td>
          <td>1.660e-02</td>
          <td>-0.030</td>
          <td>-0.100</td>
          <td>0.022</td>
          <td>2.547e-02</td>
          <td>-5.581e-03</td>
        </tr>
        <tr>
          <th rowspan="6" valign="top">3</th>
          <th>t1</th>
          <td>6.169e-03</td>
          <td>7.911e-06</td>
          <td>-2.160e-04</td>
          <td>-6.310e-02</td>
          <td>1.897e-04</td>
          <td>-4.617e-02</td>
          <td>9.580e-05</td>
          <td>-1.466e-03</td>
          <td>1.704e-03</td>
          <td>-2.679e-01</td>
          <td>...</td>
          <td>2.695e-02</td>
          <td>-6.243e-02</td>
          <td>-6.576e-03</td>
          <td>-2.369e-01</td>
          <td>-3.571e-02</td>
          <td>0.065</td>
          <td>0.135</td>
          <td>-0.079</td>
          <td>-5.365e-02</td>
          <td>-2.056e-02</td>
        </tr>
        <tr>
          <th>t2</th>
          <td>-2.295e-04</td>
          <td>6.255e-03</td>
          <td>2.639e-02</td>
          <td>6.019e-04</td>
          <td>-2.806e-04</td>
          <td>1.856e-04</td>
          <td>-7.132e-05</td>
          <td>-3.892e-02</td>
          <td>-3.453e-02</td>
          <td>-1.453e-02</td>
          <td>...</td>
          <td>-1.994e-01</td>
          <td>-3.102e-01</td>
          <td>-1.168e-01</td>
          <td>1.054e-01</td>
          <td>-4.058e-02</td>
          <td>-0.023</td>
          <td>0.200</td>
          <td>0.023</td>
          <td>1.385e-02</td>
          <td>1.724e-01</td>
        </tr>
        <tr>
          <th>t3</th>
          <td>8.457e-04</td>
          <td>1.512e-03</td>
          <td>7.034e-03</td>
          <td>-8.138e-03</td>
          <td>-1.386e-03</td>
          <td>-6.198e-03</td>
          <td>1.003e-04</td>
          <td>-9.286e-03</td>
          <td>-7.856e-03</td>
          <td>-3.736e-02</td>
          <td>...</td>
          <td>-4.493e-02</td>
          <td>1.259e-01</td>
          <td>2.542e-01</td>
          <td>-1.777e-01</td>
          <td>-1.024e-03</td>
          <td>0.042</td>
          <td>-0.038</td>
          <td>0.262</td>
          <td>2.140e-01</td>
          <td>-1.467e-01</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>8.883e-06</td>
          <td>-2.240e-04</td>
          <td>-1.036e-03</td>
          <td>-5.241e-05</td>
          <td>6.649e-05</td>
          <td>-6.664e-06</td>
          <td>2.507e-06</td>
          <td>1.321e-03</td>
          <td>1.174e-03</td>
          <td>5.725e-04</td>
          <td>...</td>
          <td>-3.039e-02</td>
          <td>2.253e-04</td>
          <td>2.894e-02</td>
          <td>-3.716e-02</td>
          <td>-4.793e-04</td>
          <td>0.008</td>
          <td>-0.047</td>
          <td>-0.006</td>
          <td>2.042e-02</td>
          <td>-3.308e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>2.507e-04</td>
          <td>1.229e-06</td>
          <td>-8.748e-06</td>
          <td>-2.571e-03</td>
          <td>6.181e-06</td>
          <td>-1.767e-03</td>
          <td>3.812e-06</td>
          <td>-5.683e-05</td>
          <td>5.614e-05</td>
          <td>-1.004e-02</td>
          <td>...</td>
          <td>-1.174e-02</td>
          <td>1.238e-03</td>
          <td>1.025e-02</td>
          <td>-3.111e-02</td>
          <td>-4.135e-03</td>
          <td>0.011</td>
          <td>0.026</td>
          <td>0.009</td>
          <td>7.308e-03</td>
          <td>-9.065e-04</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>4.657e-05</td>
          <td>2.289e-06</td>
          <td>-8.563e-06</td>
          <td>-2.151e-05</td>
          <td>8.189e-06</td>
          <td>1.310e-05</td>
          <td>-1.439e-07</td>
          <td>-1.571e-04</td>
          <td>1.600e-04</td>
          <td>1.240e-03</td>
          <td>...</td>
          <td>6.409e-03</td>
          <td>-2.870e-02</td>
          <td>6.276e-03</td>
          <td>-3.529e-02</td>
          <td>-1.277e-02</td>
          <td>0.016</td>
          <td>0.091</td>
          <td>-0.024</td>
          <td>-3.187e-02</td>
          <td>1.810e-03</td>
        </tr>
        <tr>
          <th rowspan="6" valign="top">4</th>
          <th>t1</th>
          <td>6.169e-03</td>
          <td>7.956e-06</td>
          <td>-2.157e-04</td>
          <td>-6.310e-02</td>
          <td>1.907e-04</td>
          <td>-4.617e-02</td>
          <td>9.580e-05</td>
          <td>-1.466e-03</td>
          <td>1.704e-03</td>
          <td>-2.679e-01</td>
          <td>...</td>
          <td>2.726e-02</td>
          <td>-6.325e-02</td>
          <td>-6.637e-03</td>
          <td>-2.366e-01</td>
          <td>-3.568e-02</td>
          <td>0.065</td>
          <td>0.135</td>
          <td>-0.079</td>
          <td>-5.359e-02</td>
          <td>-2.070e-02</td>
        </tr>
        <tr>
          <th>t2</th>
          <td>-1.295e-04</td>
          <td>6.253e-03</td>
          <td>2.619e-02</td>
          <td>4.724e-04</td>
          <td>-3.533e-04</td>
          <td>1.577e-04</td>
          <td>-7.179e-05</td>
          <td>-3.925e-02</td>
          <td>-3.419e-02</td>
          <td>-1.184e-02</td>
          <td>...</td>
          <td>-1.877e-01</td>
          <td>-3.177e-01</td>
          <td>-1.326e-01</td>
          <td>1.642e-01</td>
          <td>-3.435e-02</td>
          <td>-0.035</td>
          <td>0.187</td>
          <td>0.006</td>
          <td>6.929e-04</td>
          <td>1.884e-01</td>
        </tr>
        <tr>
          <th>t3</th>
          <td>3.469e-04</td>
          <td>1.510e-03</td>
          <td>7.059e-03</td>
          <td>-3.040e-03</td>
          <td>-1.396e-03</td>
          <td>-2.688e-03</td>
          <td>9.276e-05</td>
          <td>-9.173e-03</td>
          <td>-7.968e-03</td>
          <td>-1.742e-02</td>
          <td>...</td>
          <td>-2.215e-02</td>
          <td>1.256e-01</td>
          <td>2.349e-01</td>
          <td>-1.179e-01</td>
          <td>7.068e-03</td>
          <td>0.020</td>
          <td>-0.089</td>
          <td>0.243</td>
          <td>1.992e-01</td>
          <td>-1.457e-01</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>7.731e-06</td>
          <td>-2.241e-04</td>
          <td>-1.037e-03</td>
          <td>-3.840e-05</td>
          <td>6.177e-05</td>
          <td>-1.181e-05</td>
          <td>2.531e-06</td>
          <td>1.322e-03</td>
          <td>1.169e-03</td>
          <td>5.129e-04</td>
          <td>...</td>
          <td>-3.086e-02</td>
          <td>-3.467e-03</td>
          <td>3.029e-02</td>
          <td>-3.335e-02</td>
          <td>6.161e-04</td>
          <td>0.007</td>
          <td>-0.058</td>
          <td>-0.003</td>
          <td>2.650e-02</td>
          <td>-3.497e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>2.507e-04</td>
          <td>1.229e-06</td>
          <td>-8.746e-06</td>
          <td>-2.571e-03</td>
          <td>6.177e-06</td>
          <td>-1.767e-03</td>
          <td>3.812e-06</td>
          <td>-5.682e-05</td>
          <td>5.614e-05</td>
          <td>-1.004e-02</td>
          <td>...</td>
          <td>-1.174e-02</td>
          <td>1.238e-03</td>
          <td>1.026e-02</td>
          <td>-3.112e-02</td>
          <td>-4.135e-03</td>
          <td>0.011</td>
          <td>0.026</td>
          <td>0.010</td>
          <td>7.313e-03</td>
          <td>-9.084e-04</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>4.712e-05</td>
          <td>2.923e-07</td>
          <td>5.697e-05</td>
          <td>-2.570e-05</td>
          <td>3.632e-06</td>
          <td>1.221e-05</td>
          <td>-1.684e-07</td>
          <td>-1.412e-04</td>
          <td>1.769e-04</td>
          <td>1.298e-03</td>
          <td>...</td>
          <td>1.457e-02</td>
          <td>-1.664e-02</td>
          <td>4.624e-03</td>
          <td>-3.608e-02</td>
          <td>-1.077e-02</td>
          <td>0.016</td>
          <td>0.083</td>
          <td>-0.026</td>
          <td>-3.403e-02</td>
          <td>3.453e-04</td>
        </tr>
        <tr>
          <th rowspan="6" valign="top">5</th>
          <th>t1</th>
          <td>6.801e-03</td>
          <td>1.081e-05</td>
          <td>-2.255e-04</td>
          <td>-6.955e-02</td>
          <td>2.031e-04</td>
          <td>-5.058e-02</td>
          <td>1.054e-04</td>
          <td>-1.626e-03</td>
          <td>1.863e-03</td>
          <td>-2.930e-01</td>
          <td>...</td>
          <td>-1.463e-03</td>
          <td>-4.748e-02</td>
          <td>1.290e-02</td>
          <td>-2.882e-01</td>
          <td>-4.040e-02</td>
          <td>0.084</td>
          <td>0.165</td>
          <td>-0.056</td>
          <td>-3.263e-02</td>
          <td>-2.358e-02</td>
        </tr>
        <tr>
          <th>t2</th>
          <td>-2.553e-04</td>
          <td>6.819e-03</td>
          <td>2.910e-02</td>
          <td>8.055e-04</td>
          <td>-4.971e-04</td>
          <td>2.453e-04</td>
          <td>-7.785e-05</td>
          <td>-4.223e-02</td>
          <td>-3.750e-02</td>
          <td>-1.564e-02</td>
          <td>...</td>
          <td>-1.560e-01</td>
          <td>-3.697e-01</td>
          <td>-2.080e-01</td>
          <td>1.525e-01</td>
          <td>-5.946e-02</td>
          <td>-0.022</td>
          <td>0.442</td>
          <td>0.012</td>
          <td>-6.531e-02</td>
          <td>2.889e-01</td>
        </tr>
        <tr>
          <th>t3</th>
          <td>8.469e-04</td>
          <td>1.512e-03</td>
          <td>7.038e-03</td>
          <td>-8.160e-03</td>
          <td>-1.385e-03</td>
          <td>-6.209e-03</td>
          <td>1.004e-04</td>
          <td>-9.286e-03</td>
          <td>-7.856e-03</td>
          <td>-3.743e-02</td>
          <td>...</td>
          <td>-4.535e-02</td>
          <td>1.271e-01</td>
          <td>2.550e-01</td>
          <td>-1.792e-01</td>
          <td>-1.136e-03</td>
          <td>0.042</td>
          <td>-0.037</td>
          <td>0.263</td>
          <td>2.141e-01</td>
          <td>-1.473e-01</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>8.399e-06</td>
          <td>-2.241e-04</td>
          <td>-1.035e-03</td>
          <td>-4.509e-05</td>
          <td>6.317e-05</td>
          <td>-9.634e-06</td>
          <td>2.518e-06</td>
          <td>1.322e-03</td>
          <td>1.172e-03</td>
          <td>5.440e-04</td>
          <td>...</td>
          <td>-3.061e-02</td>
          <td>-9.829e-04</td>
          <td>2.993e-02</td>
          <td>-3.527e-02</td>
          <td>1.148e-04</td>
          <td>0.007</td>
          <td>-0.053</td>
          <td>-0.004</td>
          <td>2.357e-02</td>
          <td>-3.403e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>2.507e-04</td>
          <td>1.228e-06</td>
          <td>-8.742e-06</td>
          <td>-2.571e-03</td>
          <td>6.180e-06</td>
          <td>-1.767e-03</td>
          <td>3.812e-06</td>
          <td>-5.683e-05</td>
          <td>5.614e-05</td>
          <td>-1.004e-02</td>
          <td>...</td>
          <td>-1.174e-02</td>
          <td>1.241e-03</td>
          <td>1.025e-02</td>
          <td>-3.112e-02</td>
          <td>-4.135e-03</td>
          <td>0.011</td>
          <td>0.026</td>
          <td>0.009</td>
          <td>7.311e-03</td>
          <td>-9.083e-04</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>5.261e-05</td>
          <td>-1.187e-06</td>
          <td>-1.986e-04</td>
          <td>-1.310e-04</td>
          <td>-2.860e-05</td>
          <td>-4.676e-05</td>
          <td>-1.092e-07</td>
          <td>-1.774e-04</td>
          <td>1.806e-04</td>
          <td>1.011e-03</td>
          <td>...</td>
          <td>4.109e-05</td>
          <td>2.184e-02</td>
          <td>2.495e-03</td>
          <td>8.832e-02</td>
          <td>1.660e-02</td>
          <td>-0.030</td>
          <td>-0.100</td>
          <td>0.022</td>
          <td>2.547e-02</td>
          <td>-5.581e-03</td>
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
          <th rowspan="6" valign="top">5629</th>
          <th>t1</th>
          <td>-7.413e-05</td>
          <td>-8.245e-05</td>
          <td>-3.907e-04</td>
          <td>3.482e-03</td>
          <td>3.748e-05</td>
          <td>2.988e-04</td>
          <td>-2.694e-06</td>
          <td>2.441e-05</td>
          <td>1.075e-03</td>
          <td>2.742e-03</td>
          <td>...</td>
          <td>4.748e-04</td>
          <td>5.261e-03</td>
          <td>-4.662e-02</td>
          <td>5.886e-02</td>
          <td>-6.227e-03</td>
          <td>-0.017</td>
          <td>0.162</td>
          <td>-0.557</td>
          <td>6.614e-01</td>
          <td>-1.042e-01</td>
        </tr>
        <tr>
          <th>t2</th>
          <td>-4.452e-05</td>
          <td>-2.089e-04</td>
          <td>-5.165e-03</td>
          <td>2.748e-04</td>
          <td>-1.754e-04</td>
          <td>4.173e-04</td>
          <td>-3.617e-06</td>
          <td>-1.361e-04</td>
          <td>1.100e-04</td>
          <td>8.914e-04</td>
          <td>...</td>
          <td>2.047e-01</td>
          <td>6.117e-02</td>
          <td>-5.444e-02</td>
          <td>-1.529e-02</td>
          <td>-1.469e-02</td>
          <td>0.023</td>
          <td>0.183</td>
          <td>0.290</td>
          <td>-3.938e-01</td>
          <td>3.587e-01</td>
        </tr>
        <tr>
          <th>t3</th>
          <td>-1.283e-04</td>
          <td>1.048e-03</td>
          <td>8.982e-03</td>
          <td>5.709e-04</td>
          <td>-1.808e-04</td>
          <td>1.258e-03</td>
          <td>5.577e-06</td>
          <td>-5.649e-03</td>
          <td>-5.097e-03</td>
          <td>7.869e-03</td>
          <td>...</td>
          <td>-1.857e-01</td>
          <td>-2.785e-02</td>
          <td>6.353e-02</td>
          <td>5.410e-02</td>
          <td>2.473e-02</td>
          <td>-0.030</td>
          <td>-0.234</td>
          <td>0.069</td>
          <td>6.091e-02</td>
          <td>-3.214e-01</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>-3.005e-07</td>
          <td>5.476e-05</td>
          <td>6.343e-04</td>
          <td>6.332e-06</td>
          <td>2.491e-06</td>
          <td>2.715e-06</td>
          <td>5.464e-07</td>
          <td>-2.376e-04</td>
          <td>-2.019e-04</td>
          <td>-6.031e-05</td>
          <td>...</td>
          <td>-4.279e-04</td>
          <td>-3.524e-03</td>
          <td>9.710e-04</td>
          <td>-6.896e-03</td>
          <td>9.867e-04</td>
          <td>0.001</td>
          <td>-0.014</td>
          <td>-0.008</td>
          <td>2.789e-02</td>
          <td>-2.645e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>1.195e-05</td>
          <td>-1.468e-05</td>
          <td>-9.874e-05</td>
          <td>2.887e-07</td>
          <td>7.293e-06</td>
          <td>-1.234e-04</td>
          <td>1.826e-07</td>
          <td>7.492e-05</td>
          <td>1.152e-04</td>
          <td>-9.475e-04</td>
          <td>...</td>
          <td>2.369e-02</td>
          <td>1.095e-03</td>
          <td>-8.118e-03</td>
          <td>-1.289e-02</td>
          <td>-1.785e-03</td>
          <td>0.005</td>
          <td>0.015</td>
          <td>-0.021</td>
          <td>3.344e-02</td>
          <td>3.849e-03</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>-2.865e-06</td>
          <td>1.522e-05</td>
          <td>6.912e-05</td>
          <td>-4.279e-06</td>
          <td>-4.743e-06</td>
          <td>2.949e-05</td>
          <td>1.857e-07</td>
          <td>-1.044e-04</td>
          <td>-6.757e-05</td>
          <td>-1.060e-04</td>
          <td>...</td>
          <td>2.703e-02</td>
          <td>1.362e-03</td>
          <td>-5.113e-03</td>
          <td>-1.492e-02</td>
          <td>-6.335e-04</td>
          <td>0.005</td>
          <td>0.001</td>
          <td>0.027</td>
          <td>-1.418e-02</td>
          <td>1.118e-02</td>
        </tr>
        <tr>
          <th rowspan="6" valign="top">5630</th>
          <th>t1</th>
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
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
        </tr>
        <tr>
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
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000e+00</td>
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
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>-1.815e-05</td>
          <td>-9.454e-05</td>
          <td>-3.223e-04</td>
          <td>-3.568e-05</td>
          <td>1.340e-05</td>
          <td>-3.384e-05</td>
          <td>1.329e-06</td>
          <td>7.127e-04</td>
          <td>4.621e-04</td>
          <td>-6.382e-04</td>
          <td>...</td>
          <td>-3.555e-02</td>
          <td>-8.501e-03</td>
          <td>1.420e-02</td>
          <td>-1.374e-02</td>
          <td>1.390e-04</td>
          <td>0.004</td>
          <td>-0.017</td>
          <td>0.004</td>
          <td>-2.480e-04</td>
          <td>-2.458e-03</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>-1.174e-04</td>
          <td>8.335e-07</td>
          <td>-1.801e-05</td>
          <td>1.328e-03</td>
          <td>2.448e-05</td>
          <td>7.252e-04</td>
          <td>-3.178e-07</td>
          <td>-1.708e-05</td>
          <td>-1.350e-05</td>
          <td>3.852e-03</td>
          <td>...</td>
          <td>6.103e-03</td>
          <td>-6.432e-03</td>
          <td>1.082e-02</td>
          <td>-3.451e-02</td>
          <td>-4.203e-03</td>
          <td>0.011</td>
          <td>0.018</td>
          <td>-0.003</td>
          <td>7.885e-03</td>
          <td>2.321e-02</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>1.512e-05</td>
          <td>3.817e-05</td>
          <td>2.898e-04</td>
          <td>-7.734e-06</td>
          <td>-1.064e-06</td>
          <td>-1.914e-06</td>
          <td>-6.212e-07</td>
          <td>-2.275e-04</td>
          <td>-1.247e-04</td>
          <td>5.015e-04</td>
          <td>...</td>
          <td>1.508e-02</td>
          <td>-1.308e-03</td>
          <td>3.008e-03</td>
          <td>-9.727e-03</td>
          <td>3.836e-04</td>
          <td>0.001</td>
          <td>-0.003</td>
          <td>-0.030</td>
          <td>-3.103e-02</td>
          <td>-2.712e-02</td>
        </tr>
        <tr>
          <th rowspan="6" valign="top">5631</th>
          <th>t1</th>
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
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
        </tr>
        <tr>
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
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000e+00</td>
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
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>-9.862e-07</td>
          <td>5.862e-05</td>
          <td>5.579e-04</td>
          <td>1.046e-05</td>
          <td>6.905e-05</td>
          <td>5.601e-06</td>
          <td>-1.679e-06</td>
          <td>-2.394e-04</td>
          <td>-2.043e-04</td>
          <td>-3.901e-05</td>
          <td>...</td>
          <td>1.138e-03</td>
          <td>-1.261e-02</td>
          <td>1.119e-02</td>
          <td>-1.439e-02</td>
          <td>1.245e-03</td>
          <td>0.004</td>
          <td>-0.024</td>
          <td>-0.012</td>
          <td>4.048e-03</td>
          <td>1.465e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>8.388e-06</td>
          <td>-1.919e-06</td>
          <td>-7.635e-06</td>
          <td>-2.048e-04</td>
          <td>-1.957e-07</td>
          <td>-2.855e-04</td>
          <td>5.311e-07</td>
          <td>6.254e-05</td>
          <td>-5.671e-05</td>
          <td>-2.159e-03</td>
          <td>...</td>
          <td>-2.994e-02</td>
          <td>-4.564e-03</td>
          <td>1.167e-02</td>
          <td>-1.208e-02</td>
          <td>-2.319e-03</td>
          <td>0.004</td>
          <td>0.012</td>
          <td>-0.005</td>
          <td>-5.452e-03</td>
          <td>-5.399e-03</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>-4.235e-05</td>
          <td>3.105e-06</td>
          <td>1.132e-06</td>
          <td>3.700e-04</td>
          <td>3.678e-07</td>
          <td>2.318e-04</td>
          <td>-3.299e-07</td>
          <td>-1.454e-05</td>
          <td>-9.195e-06</td>
          <td>8.113e-04</td>
          <td>...</td>
          <td>7.605e-03</td>
          <td>-3.327e-03</td>
          <td>1.359e-02</td>
          <td>8.849e-04</td>
          <td>-7.085e-04</td>
          <td>0.002</td>
          <td>0.011</td>
          <td>-0.018</td>
          <td>-1.781e-02</td>
          <td>-1.326e-02</td>
        </tr>
        <tr>
          <th rowspan="6" valign="top">5632</th>
          <th>t1</th>
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
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
        </tr>
        <tr>
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
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000e+00</td>
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
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>1.756e-05</td>
          <td>-9.628e-05</td>
          <td>-3.117e-04</td>
          <td>4.014e-05</td>
          <td>1.268e-05</td>
          <td>3.502e-05</td>
          <td>1.054e-06</td>
          <td>5.821e-04</td>
          <td>6.400e-04</td>
          <td>9.411e-04</td>
          <td>...</td>
          <td>3.064e-02</td>
          <td>-2.242e-03</td>
          <td>5.439e-04</td>
          <td>-1.809e-02</td>
          <td>1.961e-04</td>
          <td>0.005</td>
          <td>-0.016</td>
          <td>0.004</td>
          <td>1.558e-04</td>
          <td>-4.253e-03</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>-1.170e-04</td>
          <td>-2.698e-07</td>
          <td>2.598e-05</td>
          <td>1.325e-03</td>
          <td>-3.278e-05</td>
          <td>7.227e-04</td>
          <td>-2.756e-06</td>
          <td>1.174e-05</td>
          <td>1.113e-05</td>
          <td>3.844e-03</td>
          <td>...</td>
          <td>1.025e-02</td>
          <td>4.245e-03</td>
          <td>2.363e-03</td>
          <td>-2.781e-02</td>
          <td>-4.249e-03</td>
          <td>0.009</td>
          <td>0.026</td>
          <td>0.024</td>
          <td>1.089e-02</td>
          <td>-1.333e-02</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>1.548e-05</td>
          <td>-4.294e-05</td>
          <td>-2.770e-04</td>
          <td>-1.257e-05</td>
          <td>-3.928e-06</td>
          <td>-5.062e-06</td>
          <td>3.049e-07</td>
          <td>1.836e-04</td>
          <td>2.254e-04</td>
          <td>5.881e-04</td>
          <td>...</td>
          <td>2.334e-02</td>
          <td>-1.135e-03</td>
          <td>5.283e-03</td>
          <td>-1.865e-03</td>
          <td>-3.915e-03</td>
          <td>0.002</td>
          <td>0.024</td>
          <td>-0.032</td>
          <td>-2.891e-02</td>
          <td>-1.447e-02</td>
        </tr>
        <tr>
          <th rowspan="6" valign="top">5633</th>
          <th>t1</th>
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
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
        </tr>
        <tr>
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
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000e+00</td>
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
          <td>0.000e+00</td>
          <td>0.000e+00</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000</td>
          <td>0.000e+00</td>
          <td>0.000e+00</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>-3.006e-07</td>
          <td>5.476e-05</td>
          <td>6.343e-04</td>
          <td>6.334e-06</td>
          <td>2.491e-06</td>
          <td>2.716e-06</td>
          <td>5.464e-07</td>
          <td>-2.376e-04</td>
          <td>-2.019e-04</td>
          <td>-6.030e-05</td>
          <td>...</td>
          <td>-4.279e-04</td>
          <td>-3.524e-03</td>
          <td>9.710e-04</td>
          <td>-6.896e-03</td>
          <td>9.867e-04</td>
          <td>0.001</td>
          <td>-0.014</td>
          <td>-0.008</td>
          <td>2.789e-02</td>
          <td>-2.645e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>-1.723e-06</td>
          <td>1.278e-06</td>
          <td>-1.805e-06</td>
          <td>1.940e-04</td>
          <td>3.380e-07</td>
          <td>-8.450e-06</td>
          <td>3.548e-08</td>
          <td>-4.728e-05</td>
          <td>4.650e-05</td>
          <td>-2.113e-04</td>
          <td>...</td>
          <td>2.084e-02</td>
          <td>1.171e-03</td>
          <td>-6.235e-03</td>
          <td>-1.349e-02</td>
          <td>-1.096e-03</td>
          <td>0.005</td>
          <td>0.006</td>
          <td>0.005</td>
          <td>4.639e-03</td>
          <td>6.872e-03</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>7.271e-06</td>
          <td>3.394e-06</td>
          <td>-2.722e-06</td>
          <td>-1.478e-04</td>
          <td>4.108e-07</td>
          <td>-5.572e-05</td>
          <td>2.948e-07</td>
          <td>-1.384e-05</td>
          <td>-1.663e-05</td>
          <td>-6.516e-04</td>
          <td>...</td>
          <td>2.914e-02</td>
          <td>1.305e-03</td>
          <td>-6.509e-03</td>
          <td>-1.448e-02</td>
          <td>-1.144e-03</td>
          <td>0.005</td>
          <td>0.008</td>
          <td>0.008</td>
          <td>7.160e-03</td>
          <td>8.942e-03</td>
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

Be careful about what you're accessing.

.. code:: python

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
    plate_stress = ['_add_new_eid', '_add_new_node', 'subtitle', 'words', 's_code', 'is_built', 'stress_bits', 'load_set', 'itotal', '_add', '_ntotals', 'nnodes', 'sort_bits', 'isubcase', 'element_name', 'itime', 'nonlinear_factor', 'title', '_times', 'ntotal', 'approach_code', 'is_stress_flag', 'label', 'element_node', 'is_msc', 'num_wide', 'mode', 'format_code', 'device_code', 'modes', '_times_dtype', 'thermal_bits', 'mode2s', 'data_frame', 'mode2', 'dt', 'is_strain_flag', 'data', 'cycle', 'name', 'nelements', 'eigr', 'ielement', 'thermal', 'analysis_code', 'eigrs', 'table_code', 'element_type', 'table_name', 'data_code', 'isTransient', 'sort_code', 'cycles', 'ntimes', 'data_names']
    
    data_code_keys = [u'subtitle', u'stress_bits', u'load_set', u'thermal', u's_code', u'sort_bits', u'isubcase', u'element_name', u'mode2', u'title', u'approach_code', u'is_stress_flag', u'label', u'is_msc', u'num_wide', u'format_code', u'device_code', u'_times_dtype', u'thermal_bits', u'nonlinear_factor', u'is_strain_flag', u'cycle', u'name', u'eigr', u'analysis_code', u'table_code', u'element_type', u'table_name', u'mode', u'sort_code', u'data_names']
    
    name = u'mode'
    list-type variables = [u'mode', u'eigr', u'mode2', u'cycle']
    modes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167]
    

Similar to the BDF, we can use object\_attributes/methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    #print "attributes =", object_attributes(plate_stress)
    print("methods = %s\n" % object_methods(plate_stress))
    print('methods2= %s\n' % plate_stress.object_methods())
    print("headers = %s\n" % plate_stress.get_headers())


.. parsed-literal::

    methods = ['apply_data_code', 'approach_code_str', 'build', 'build_dataframe', 'cast_grid_type', 'code_information', 'eid_to_element_node_index', 'get_data_code', 'get_element_index', 'get_element_type', 'get_headers', 'get_nnodes_bilinear', 'get_stats', 'get_unsteady_value', 'is_bilinear', 'is_complex', 'is_curvature', 'is_fiber_distance', 'is_magnitude_phase', 'is_max_shear', 'is_real', 'is_sort1', 'is_sort2', 'is_strain', 'is_stress', 'is_thermal', 'is_von_mises', 'object_attributes', 'object_methods', 'print_data_members', 'print_table_code', 'recast_gridtype_as_string', 'set_table_type', 'update_data_code', 'update_dt', 'write_f06']
    
    methods2= ['apply_data_code', 'approach_code_str', 'build', 'build_dataframe', 'cast_grid_type', 'code_information', 'eid_to_element_node_index', 'get_data_code', 'get_element_index', 'get_element_type', 'get_headers', 'get_nnodes_bilinear', 'get_stats', 'get_unsteady_value', 'is_bilinear', 'is_complex', 'is_curvature', 'is_fiber_distance', 'is_magnitude_phase', 'is_max_shear', 'is_real', 'is_sort1', 'is_sort2', 'is_strain', 'is_stress', 'is_thermal', 'is_von_mises', 'print_data_members', 'print_table_code', 'recast_gridtype_as_string', 'set_table_type', 'update_data_code', 'update_dt', 'write_f06']
    
    headers = [u'fiber_distance', u'oxx', u'oyy', u'txy', u'angle', u'omax', u'omin', u'von_mises']
    
    

Number of Nodes on a CQUAD4
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  For linear CQUAD4s, there is 1 centroidal stress at two locations
-  For bilinear quads, there are 5 stresses at two locations (4 nodes +
   centroidal)
-  node\_id=0 indicates a centroidal quantity
-  CTRIA3s are always centroidal

What sets this?
^^^^^^^^^^^^^^^

::

    STRESS(real, sort1, BILIN) = ALL   # bilinear cquad
    STRESS(real, sort1, CENT) = ALL    # linear quad

    STRAIN(real, sort1, BILIN) = ALL   # bilinear cquad
    STRAIN(real, sort1, CENT) = ALL    # linear quad

How do we know if we're bilinear?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    print("is_bilinear = %s\n" % plate_stress.is_bilinear())

What locations are chosen?
^^^^^^^^^^^^^^^^^^^^^^^^^^

That depends on fiber distance/fiber curvature... - fiber\_curvature -
mean stress (oa) & slope (om)

::

    $$ \sigma_{top} = \sigma_{alt} + \frac{t}{2} \sigma_{mean}$$

    $$ \sigma_{btm} = \sigma_{alt} + \frac{t}{2} \sigma_{mean}$$

-  fiber\_distance - upper and lower surface stress (o\_top; o\_btm)
-  If you have stress, fiber\_distance is always returned regardless of
   your option.

What sets this?
^^^^^^^^^^^^^^^

::

    STRAIN(real, sort1, FIBER) = ALL   # fiber distance/default
    STRAIN(real, sort1, STRCUR) = ALL  # strain curvature

How do we know if we're using fiber\_distance?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    print("is_fiber_distance = %s" % plate_stress.is_fiber_distance())

Accessing results
-----------------

Note that this is intentionally done iinefficiently to access specific entries in order to explain the data structure.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

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
    if plate_stress.is_von_mises():  # True
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
    we have von mises data; ovm=[ 54.222   5.041  13.143 ...,   2.34    6.146   7.368]
    
    [layer1, layer2, ...] = [ 54.222   5.041  13.143  21.222  78.544  17.91 ]
    ieid1000 = [1998 1999]
    ovm_mode6_eid1000 = [ 90.618  94.09 ] -> 94.0905
    

.. code:: python

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
    
    if plate_stress.is_fiber_distance():
        print("fiber_distance = %s" % plate_stress.data[imode, ieid10, 0])
    else:
        print("curvature = %s" % plate_stress.data[imode, ieid10, 0])


.. parsed-literal::

    ieid10 = 19
    [10  0]
    ps.modes = 6
    ps.cycles = 9.80908925027e-45
    oxx = -18.8701
    oyy = -20.1605
    txy = -8.30956
    omax = -11.1807
    omin = -27.8499
    ovm/max_shear = 24.2743
    fiber_distance = -0.4
    

.. code:: python

    from pyNastran.bdf.bdf import read_bdf
    bdf_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_4pt.dat'))
    model = read_bdf(bdf_filename, debug=False)
    mass, cg, I = model.mass_properties()


.. parsed-literal::

    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=1 midsurface: z1=0.400000006 z2=-0.400000006 t=0.035999998 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=2 midsurface: z1=0.400000006 z2=-0.400000006 t=0.054000005 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=3 midsurface: z1=0.400000006 z2=-0.400000006 t=0.017999999 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=7 midsurface: z1=0.418000013 z2=-0.418000013 t=0.035999998 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=34 midsurface: z1=0.194000006 z2=-0.194000006 t=0.0186 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=37 midsurface: z1=0.308999985 z2=-0.308999985 t=0.0186 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=38 midsurface: z1=0.284000009 z2=-0.284000009 t=0.0186 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=46 midsurface: z1=0.199000001 z2=-0.199000001 t=0.0186 not in range of -1.5t < zi < 1.5t
    

Let's print out the actual mass properties from the OP2 and get the same result as the F06
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We need ``PARAM,POSTEXT,YES`` in out BDF to get the Grid Point Weight
Table

.. code:: python

    gpw = op2.grid_point_weight
    #print(gpw.object_attributes())
    
    print(gpw)


.. parsed-literal::

                               O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R
    0                                                     REFERENCE POINT =        0
                                                                    M O
                          *  1.774601E+00  1.402827E-19  2.212874E-19 -1.821217E-17 -3.277270E+01  4.490826E+00 *
                          *  1.402827E-19  1.774601E+00 -4.675622E-19  3.277270E+01 -2.898787E-17 -6.007191E-02 *
                          *  2.212874E-19 -4.675622E-19  1.774601E+00 -4.490826E+00  6.007191E-02 -3.871152E-19 *
                          * -1.821217E-17  3.277270E+01 -4.490826E+00  1.322289E+03  1.414696E+00 -1.250574E+00 *
                          * -3.277270E+01 -2.898787E-17  6.007191E-02  1.414696E+00  1.227074E+03 -2.187713E+02 *
                          *  4.490826E+00 -6.007191E-02 -3.871152E-19 -1.250574E+00 -2.187713E+02  4.272278E+02 *
                                                                     S
                                               *  1.000000E+00  0.000000E+00  0.000000E+00 *
                                               *  0.000000E+00  1.000000E+00  0.000000E+00 *
                                               *  0.000000E+00  0.000000E+00  1.000000E+00 *
                                   DIRECTION
                              MASS AXIS SYSTEM (S)     MASS              X-C.G.        Y-C.G.        Z-C.G.
                                      X            1.774601E+00     -1.026268E-17 -2.530611E+00 -1.846764E+01
                                      Y            1.774601E+00     -3.385094E-02 -1.633486E-17 -1.846764E+01
                                      Z            1.774601E+00     -3.385094E-02 -2.530611E+00 -2.181421E-19
                                                                    I(S)
                                               *  7.056896E+02 -1.566714E+00  1.411869E-01 *
                                               * -1.566714E+00  6.218375E+02  1.358363E+02 *
                                               *  1.411869E-01  1.358363E+02  4.158613E+02 *
                                                                    I(Q)
                                               *  6.891835E+02                             *
                                               *                3.483842E+02               *
                                               *                              7.058207E+02 *
                                                                     Q
                                               *  8.846355E-02  1.596853E-03  9.960781E-01 *
                                               * -8.920128E-01 -4.448861E-01  7.993453E-02 *
                                               *  4.432690E-01 -8.955857E-01 -3.793179E-02 *
    
    PAGE 1
    
    

We can also write the full F06
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    import getpass
    name = getpass.getuser()
    os.chdir(os.path.join(r'C:\Users', name, 'Desktop'))
    
    # write the F06 with Real/Imaginary or Magnitude/Phase
    # only matters for complex results
    op2.write_f06('isat.f06', is_mag_phase=False)
    
    !head -n 40 isat.f06


.. parsed-literal::

    F06:
     grid_point_weight
    RealEigenvalues    case=u'ISAT_SM_LAUNCH_4PT MODES TO 400 HZ'
     RealEigenvectorArray SUBCASE=1 SUBTITLE=
     RealCBarForceArray   SUBCASE=1 SUBTITLE=  - CBAR-34
     RealPlateStressArray SUBCASE=1 SUBTITLE=  - CQUAD4
     RealPlateStressArray SUBCASE=1 SUBTITLE=  - CTRIA3
                               O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R
    0                                                     REFERENCE POINT =        0
                                                                    M O
                          *  1.774601E+00  1.402827E-19  2.212874E-19 -1.821217E-17 -3.277270E+01  4.490826E+00 *
                          *  1.402827E-19  1.774601E+00 -4.675622E-19  3.277270E+01 -2.898787E-17 -6.007191E-02 *
                          *  2.212874E-19 -4.675622E-19  1.774601E+00 -4.490826E+00  6.007191E-02 -3.871152E-19 *
                          * -1.821217E-17  3.277270E+01 -4.490826E+00  1.322289E+03  1.414696E+00 -1.250574E+00 *
                          * -3.277270E+01 -2.898787E-17  6.007191E-02  1.414696E+00  1.227074E+03 -2.187713E+02 *
                          *  4.490826E+00 -6.007191E-02 -3.871152E-19 -1.250574E+00 -2.187713E+02  4.272278E+02 *
                                                                     S
                                               *  1.000000E+00  0.000000E+00  0.000000E+00 *
                                               *  0.000000E+00  1.000000E+00  0.000000E+00 *
                                               *  0.000000E+00  0.000000E+00  1.000000E+00 *
                                   DIRECTION
                              MASS AXIS SYSTEM (S)     MASS              X-C.G.        Y-C.G.        Z-C.G.
                                      X            1.774601E+00     -1.026268E-17 -2.530611E+00 -1.846764E+01
                                      Y            1.774601E+00     -3.385094E-02 -1.633486E-17 -1.846764E+01
                                      Z            1.774601E+00     -3.385094E-02 -2.530611E+00 -2.181421E-19
                                                                    I(S)
                                               *  7.056896E+02 -1.566714E+00  1.411869E-01 *
                                               * -1.566714E+00  6.218375E+02  1.358363E+02 *
                                               *  1.411869E-01  1.358363E+02  4.158613E+02 *
                                                                    I(Q)
                                               *  6.891835E+02                             *
                                               *                3.483842E+02               *
                                               *                              7.058207E+02 *
                                                                     Q
                                               *  8.846355E-02  1.596853E-03  9.960781E-01 *
                                               * -8.920128E-01 -4.448861E-01  7.993453E-02 *
                                               *  4.432690E-01 -8.955857E-01 -3.793179E-02 *
    
    1    ISAT_SM_LAUNCH_4PT MODES TO 400 HZ                                     JANUARY   4, 2016  pyNastran v0.8.0+dev.82cefee  PAGE     1
    
    1    ISAT_SM_LAUNCH_4PT MODES TO 400 HZ                                     JANUARY   4, 2016  pyNastran v0.8.0+dev.82cefee  PAGE     2
         DEFAULT                                                                                                                        
    
                                                  R E A L   E I G E N V A L U E S
                                                 ISAT_SM_LAUNCH_4PT MODES TO 400 HZ
       MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED
        NO.       ORDER                                                                       MASS              STIFFNESS
    

.. code:: python

    #from IPython.display import display, Math, Latex

The mass results are different as pyNastran's mass assumes point masses

.. math:: m_{plates} = A * (rho * t + nsm)

.. math:: m_{solid} = V * rho

.. math:: m_{bars} = L * (rho * A + nsm)

.. math:: I = m*r^2

The larger your model is and the further from the origin, the more
accurate the result. For some applications (e.g. a weight breakdown),
this is probably be fine.

.. code:: python

    print('cg =\n%s' % gpw.cg)
    print('cg = %s' % cg)


.. parsed-literal::

    cg =
    [[ -1.026e-17  -2.531e+00  -1.847e+01]
     [ -3.385e-02  -1.633e-17  -1.847e+01]
     [ -3.385e-02  -2.531e+00  -2.181e-19]]
    cg = [ -0.035  -2.623 -18.53 ]
    

It's not like Nastran is perfect either.
----------------------------------------

Limitations
~~~~~~~~~~~

1. You cannot do weight statements in Nastran by
   component/property/material.

2. Everything is always summmed up (e.g. you can have different geometry
   in Subcase 2 and MPCs connecting physical geomtry, with other parts
   flying off into space).

These are things that pyNastran ``can`` do.

.. code:: python

    from pyNastran.bdf.bdf import read_bdf
    bdf_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_4pt.dat'))
    model = read_bdf(bdf_filename, debug=False)


.. parsed-literal::

    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=1 midsurface: z1=0.400000006 z2=-0.400000006 t=0.035999998 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=2 midsurface: z1=0.400000006 z2=-0.400000006 t=0.054000005 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=3 midsurface: z1=0.400000006 z2=-0.400000006 t=0.017999999 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=7 midsurface: z1=0.418000013 z2=-0.418000013 t=0.035999998 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=34 midsurface: z1=0.194000006 z2=-0.194000006 t=0.0186 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=37 midsurface: z1=0.308999985 z2=-0.308999985 t=0.0186 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=38 midsurface: z1=0.284000009 z2=-0.284000009 t=0.0186 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.py                  lineNo=1434   PSHELL pid=46 midsurface: z1=0.199000001 z2=-0.199000001 t=0.0186 not in range of -1.5t < zi < 1.5t
    

Let's get the breakdown by property ID
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    from six import iteritems
    #help(model.mass_properties)
    
    pid_to_eids_map = model.get_element_ids_dict_with_pids()
    #print(pid_to_eids_map.keys())
    print('pid, mass, cg, [ixx, iyy, izz, ixy, ixz]')
    for pid, eids in sorted(iteritems(pid_to_eids_map)):
        mass, cg, inertia = model.mass_properties(element_ids=eids, reference_point=[0., 0., 0.])
        print('%-3s %-.6f %-38s %s' % (pid, mass, cg, inertia))


.. parsed-literal::

    pid, mass, cg, [ixx, iyy, izz, ixy, ixz]
    1   0.027278 [  4.297e-15   2.980e-15  -2.000e+01]  [  1.461e+01   1.746e+01   4.384e+00   4.098e-17   1.450e-15   4.619e-16]
    2   0.047993 [ -6.506e-16   2.530e-16  -2.000e+01]  [  3.723e+01   3.723e+01   1.245e+01  -3.469e-17  -5.690e-16   2.047e-16]
    3   0.020998 [  6.842e-17  -4.699e-16  -2.000e+01]  [  1.428e+01   1.231e+01   5.270e+00   2.992e-17   5.281e-16   1.011e-15]
    4   0.012216 [  0.043   0.438 -19.702]              [  7.090e+00   7.972e+00   2.021e+00   1.057e-02  -5.272e-03  -5.334e-02]
    5   0.330158 [  0.    2.2 -20. ]                    [  1.970e+02   1.604e+02   4.335e+01   0.000e+00   4.441e-16  -1.453e+01]
    7   0.027813 [  6.578e-17  -1.143e-14  -2.000e+01]  [  1.927e+01   1.927e+01   9.438e+00  -3.408e-15  -2.069e-16  -1.751e-15]
    8   0.081584 [  0.000e+00   6.804e-16  -2.000e+01]  [  4.772e+01   4.772e+01   3.017e+01  -8.882e-16   0.000e+00  -8.882e-16]
    9   0.012578 [  0.000e+00  -7.585e-16  -2.000e+01]  [  7.788e+00   7.788e+00   3.064e+00   5.551e-17  -4.857e-17  -1.180e-16]
    10  0.000236 [  0.000e+00  -4.595e-16  -2.000e+01]  [  1.290e-01   1.290e-01   5.747e-02  -8.674e-19  -1.735e-18  -8.674e-19]
    11  0.041700 [ -1.025  23.773 -12.016]              [ 30.253   7.053  23.957  -1.053   0.776 -11.967]
    12  0.000457 [  0.     -5.92   20.506]              [ 0.221  0.205  0.016  0.     0.    -0.054]
    13  0.000000 [ 0.  0.  0.]                          [ 0.  0.  0.  0.  0.  0.]
    14  0.000353 [ -2.305e-16   0.000e+00   1.439e+01]  [  7.655e-02   8.537e-02   8.821e-03   0.000e+00  -2.168e-18   0.000e+00]
    15  0.000000 [ 0.  0.  0.]                          [ 0.  0.  0.  0.  0.  0.]
    16  0.000000 [ 0.  0.  0.]                          [ 0.  0.  0.  0.  0.  0.]
    19  0.017749 [ -0.23    6.021 -35.642]              [ 24.961  26.944   6.461  -0.078   0.151  -3.954]
    20  0.163082 [  6.808e-16   0.000e+00  -1.855e+01]  [  6.510e+01   9.086e+01   2.576e+01   0.000e+00  -1.776e-15   0.000e+00]
    21  0.003625 [ -1.196e-16  -1.077e-15  -2.000e+01]  [  2.178e+00   2.178e+00   1.410e+00   0.000e+00   0.000e+00  -2.776e-17]
    22  0.000000 [ 0.  0.  0.]                          [ 0.  0.  0.  0.  0.  0.]
    23  0.000000 [ 0.  0.  0.]                          [ 0.  0.  0.  0.  0.  0.]
    33  0.001346 [ -3.612e-14  -2.175e+00   3.691e-01]  [  8.358e-02   8.533e-02   1.683e-01   1.802e-16  -3.670e-17  -1.832e-03]
    34  0.003561 [ -5.899e-17  -1.903e-18   1.483e+01]  [  1.054e+00   1.054e+00   6.701e-02  -5.828e-19  -4.120e-18  -8.674e-19]
    35  0.000000 [ 0.  0.  0.]                          [ 0.  0.  0.  0.  0.  0.]
    36  0.007197 [  3.676e-15   0.000e+00  -1.478e+01]  [  2.408e+00   5.178e+00   3.020e+00   6.939e-18   8.153e-17  -1.041e-17]
    37  0.094566 [ -4.439e-15   0.000e+00  -1.950e+01]  [  4.493e+01   8.867e+01   4.649e+01   5.454e-14  -2.776e-16   0.000e+00]
    38  0.007602 [  6.830e-06  -9.329e+00   2.731e+01]  [  7.013e+00   6.884e+00   1.223e+00  -2.922e-07   1.641e-06  -1.867e+00]
    39  0.002433 [  1.348e-13  -8.954e+00   4.040e+00]  [  2.568e-01   8.454e-02   2.187e-01  -2.920e-15   1.207e-15  -9.123e-02]
    41  0.000735 [ -9.583e-16  -1.843e-17   2.193e+00]  [  1.226e-02   6.032e-02   6.220e-02   5.127e-16   8.674e-19  -5.421e-20]
    42  0.008854 [ -1.554  20.121 -19.007]              [ 7.95   4.679  3.9   -0.277  0.318 -3.386]
    43  0.012241 [  6.191e-15   2.214e-18  -1.950e+01]  [  6.882e+00   1.224e+01   6.320e+00   7.003e-15  -4.749e-17   1.762e-18]
    46  0.003671 [  3.544e-15   7.383e-18   1.528e+01]  [  1.035e+00   1.205e+00   3.350e-01   1.239e-07   8.544e-17   4.120e-18]
    60  0.000000 [ 0.  0.  0.]                          [ 0.  0.  0.  0.  0.  0.]
    61  0.000000 [ 0.  0.  0.]                          [ 0.  0.  0.  0.  0.  0.]
    
