
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
    op2_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_4pt.op2'))
    assert os.path.exists(op2_filename), print_bad_path(op2_filename)
    
    # define the input file with a file path
    op2 = read_op2(op2_filename, build_dataframe=True, debug=False)



.. raw:: html

    <text style=color:green>INFO:    op2_scalar.py:1469           op2_filename = 'c:\\nasa\\m4\\formats\\git\\pynastran_1.2\\models\\iSat\\ISat_Launch_Sm_4pt.op2'
    </text>


.. parsed-literal::

    c:\nasa\m4\formats\git\pynastran_1.2\pyNastran\op2\op2.py:752: FutureWarning: 
    Panel is deprecated and will be removed in a future version.
    The recommended way to represent these types of 3-dimensional data are with a MultiIndex on a DataFrame, via the Panel.to_frame() method
    Alternatively, you can use the xarray package http://xarray.pydata.org/en/stable/.
    Pandas provides a `.to_xarray()` method to help automate this conversion.
    
      obj.build_dataframe()
    

OP2 Introspection
-----------------

The ``get_op2_stats()`` function lets you quickly understand what in an
op2.

.. code:: ipython3

    print(op2.get_op2_stats())


.. parsed-literal::

    eigenvectors[1]
      isubcase = 1
      type=RealEigenvectorArray ntimes=167 nnodes=5379, table_name=OUGV1
      data: [t1, t2, t3, r1, r2, r3] shape=[167, 5379, 6] dtype=float32
      node_gridtype.shape = (5379, 2)
      sort1
      modes = [  1   2   3 ... 165 166 167]
      eigns = [   2757.896    3568.136    9686.188 ... 6162773.5   6169898.5
     6229583.   ]
      mode_cycles = [  8.358   9.507  15.664 ... 395.101 395.329 397.237]
    
    cbar_force[1]
      type=RealCBarForceArray ntimes=167 nelements=827; table_name='OEF1X'
      data: [ntimes, nnodes, 8] where 8=[bending_moment_a1, bending_moment_a2, bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque]
      data.shape = (167, 827, 8)
      element.shape = (827,)
      element name: CBAR-34
      sort1
      modes = [  1   2   3 ... 165 166 167]
      eigns = [   2757.896    3568.136    9686.188 ... 6162773.5   6169898.5
     6229583.   ]
      cycles = [  8.358   9.507  15.664 ... 395.101 395.329 397.237]
    
    ctria3_stress[1]
      type=RealPlateStressArray ntimes=167 nelements=32 nnodes_per_element=1 nlayers=2 ntotal=64
      data: [ntimes, ntotal, 8] where 8=[fiber_distance, oxx, oyy, txy, angle, omax, omin, von_mises]
      element_node.shape = (64, 2)
      data.shape=(167, 64, 8)
      element type: CTRIA3
      s_code: 1
      sort1
      modes = [  1   2   3 ... 165 166 167]
      eigns = [   2757.896    3568.136    9686.188 ... 6162773.5   6169898.5
     6229583.   ]
      mode2s = [0 0 0 ... 0 0 0]
      cycles = [  8.358   9.507  15.664 ... 395.101 395.329 397.237]
    
    cquad4_stress[1]
      type=RealPlateStressArray ntimes=167 nelements=4580 nnodes_per_element=1 nlayers=2 ntotal=9160
      data: [ntimes, ntotal, 8] where 8=[fiber_distance, oxx, oyy, txy, angle, omax, omin, von_mises]
      element_node.shape = (9160, 2)
      data.shape=(167, 9160, 8)
      element type: CQUAD4
      s_code: 1
      sort1
      modes = [  1   2   3 ... 165 166 167]
      eigns = [   2757.896    3568.136    9686.188 ... 6162773.5   6169898.5
     6229583.   ]
      mode2s = [0 0 0 ... 0 0 0]
      cycles = [  8.358   9.507  15.664 ... 395.101 395.329 397.237]
    
    
    

If that’s too long…
~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    print(op2.get_op2_stats(short=True))


.. parsed-literal::

    eigenvectors[1]
    cbar_force[1]
    ctria3_stress[1]
    cquad4_stress[1]
    
    

Acccessing the Eigenvectors object
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

    loadcases = dict_keys([1])
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
          <th>59.734</th>
          <th>98.418</th>
          <th>127.102</th>
          <th>127.585</th>
          <th>129.107</th>
          <th>135.087</th>
          <th>136.351</th>
          <th>136.445</th>
          <th>178.721</th>
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
          <th>2495.913</th>
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
          <td>-5.548e-03</td>
          <td>4.671e-06</td>
          <td>-1.816e-04</td>
          <td>-5.670e-02</td>
          <td>-1.721e-04</td>
          <td>-4.175e-02</td>
          <td>-8.632e-05</td>
          <td>-1.341e-03</td>
          <td>1.582e-03</td>
          <td>2.438e-01</td>
          <td>...</td>
          <td>5.723e-02</td>
          <td>-5.369e-02</td>
          <td>-3.837e-02</td>
          <td>1.326e-01</td>
          <td>-1.973e-02</td>
          <td>-0.028</td>
          <td>0.033</td>
          <td>-0.104</td>
          <td>6.919e-02</td>
          <td>1.904e-02</td>
        </tr>
        <tr>
          <th>t2</th>
          <td>2.133e-04</td>
          <td>5.699e-03</td>
          <td>2.393e-02</td>
          <td>5.803e-04</td>
          <td>1.812e-04</td>
          <td>1.971e-04</td>
          <td>6.526e-05</td>
          <td>-3.562e-02</td>
          <td>-3.164e-02</td>
          <td>1.291e-02</td>
          <td>...</td>
          <td>-3.091e-01</td>
          <td>-3.746e-01</td>
          <td>-5.840e-02</td>
          <td>2.385e-02</td>
          <td>-5.889e-02</td>
          <td>-0.015</td>
          <td>0.177</td>
          <td>-0.011</td>
          <td>-5.252e-02</td>
          <td>-1.187e-01</td>
        </tr>
        <tr>
          <th>t3</th>
          <td>-8.469e-04</td>
          <td>1.512e-03</td>
          <td>7.038e-03</td>
          <td>-8.160e-03</td>
          <td>1.385e-03</td>
          <td>-6.210e-03</td>
          <td>-1.004e-04</td>
          <td>-9.286e-03</td>
          <td>-7.856e-03</td>
          <td>3.756e-02</td>
          <td>...</td>
          <td>-4.534e-02</td>
          <td>1.271e-01</td>
          <td>2.550e-01</td>
          <td>1.792e-01</td>
          <td>-1.136e-03</td>
          <td>-0.042</td>
          <td>-0.037</td>
          <td>0.263</td>
          <td>-2.141e-01</td>
          <td>1.472e-01</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>-8.399e-06</td>
          <td>-2.241e-04</td>
          <td>-1.035e-03</td>
          <td>-4.509e-05</td>
          <td>-6.317e-05</td>
          <td>-9.635e-06</td>
          <td>-2.518e-06</td>
          <td>1.322e-03</td>
          <td>1.172e-03</td>
          <td>-5.433e-04</td>
          <td>...</td>
          <td>-3.061e-02</td>
          <td>-9.825e-04</td>
          <td>2.993e-02</td>
          <td>3.527e-02</td>
          <td>1.148e-04</td>
          <td>-0.007</td>
          <td>-0.053</td>
          <td>-0.004</td>
          <td>-2.357e-02</td>
          <td>3.403e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>-2.507e-04</td>
          <td>1.228e-06</td>
          <td>-8.730e-06</td>
          <td>-2.571e-03</td>
          <td>-6.174e-06</td>
          <td>-1.767e-03</td>
          <td>-3.812e-06</td>
          <td>-5.683e-05</td>
          <td>5.614e-05</td>
          <td>1.008e-02</td>
          <td>...</td>
          <td>-1.174e-02</td>
          <td>1.241e-03</td>
          <td>1.025e-02</td>
          <td>3.112e-02</td>
          <td>-4.135e-03</td>
          <td>-0.011</td>
          <td>0.026</td>
          <td>0.009</td>
          <td>-7.311e-03</td>
          <td>9.082e-04</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>-5.261e-05</td>
          <td>-1.187e-06</td>
          <td>-1.986e-04</td>
          <td>-1.310e-04</td>
          <td>2.861e-05</td>
          <td>-4.677e-05</td>
          <td>1.092e-07</td>
          <td>-1.774e-04</td>
          <td>1.806e-04</td>
          <td>-1.008e-03</td>
          <td>...</td>
          <td>4.063e-05</td>
          <td>2.184e-02</td>
          <td>2.495e-03</td>
          <td>-8.831e-02</td>
          <td>1.660e-02</td>
          <td>0.030</td>
          <td>-0.100</td>
          <td>0.022</td>
          <td>-2.547e-02</td>
          <td>5.581e-03</td>
        </tr>
        <tr>
          <th rowspan="6" valign="top">2</th>
          <th>t1</th>
          <td>-5.548e-03</td>
          <td>4.671e-06</td>
          <td>-1.816e-04</td>
          <td>-5.670e-02</td>
          <td>-1.721e-04</td>
          <td>-4.175e-02</td>
          <td>-8.632e-05</td>
          <td>-1.341e-03</td>
          <td>1.582e-03</td>
          <td>2.438e-01</td>
          <td>...</td>
          <td>5.723e-02</td>
          <td>-5.369e-02</td>
          <td>-3.837e-02</td>
          <td>1.326e-01</td>
          <td>-1.973e-02</td>
          <td>-0.028</td>
          <td>0.033</td>
          <td>-0.104</td>
          <td>6.919e-02</td>
          <td>1.904e-02</td>
        </tr>
        <tr>
          <th>t2</th>
          <td>1.081e-04</td>
          <td>5.696e-03</td>
          <td>2.353e-02</td>
          <td>3.182e-04</td>
          <td>2.384e-04</td>
          <td>1.036e-04</td>
          <td>6.548e-05</td>
          <td>-3.598e-02</td>
          <td>-3.128e-02</td>
          <td>1.090e-02</td>
          <td>...</td>
          <td>-3.090e-01</td>
          <td>-3.309e-01</td>
          <td>-5.341e-02</td>
          <td>-1.528e-01</td>
          <td>-2.568e-02</td>
          <td>0.045</td>
          <td>-0.022</td>
          <td>0.034</td>
          <td>-1.035e-01</td>
          <td>-1.075e-01</td>
        </tr>
        <tr>
          <th>t3</th>
          <td>-3.455e-04</td>
          <td>1.510e-03</td>
          <td>7.056e-03</td>
          <td>-3.018e-03</td>
          <td>1.398e-03</td>
          <td>-2.676e-03</td>
          <td>-9.274e-05</td>
          <td>-9.172e-03</td>
          <td>-7.968e-03</td>
          <td>1.739e-02</td>
          <td>...</td>
          <td>-2.187e-02</td>
          <td>1.246e-01</td>
          <td>2.345e-01</td>
          <td>1.170e-01</td>
          <td>7.135e-03</td>
          <td>-0.020</td>
          <td>-0.090</td>
          <td>0.244</td>
          <td>-1.995e-01</td>
          <td>1.454e-01</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>-8.399e-06</td>
          <td>-2.241e-04</td>
          <td>-1.035e-03</td>
          <td>-4.509e-05</td>
          <td>-6.317e-05</td>
          <td>-9.635e-06</td>
          <td>-2.518e-06</td>
          <td>1.322e-03</td>
          <td>1.172e-03</td>
          <td>-5.433e-04</td>
          <td>...</td>
          <td>-3.061e-02</td>
          <td>-9.825e-04</td>
          <td>2.993e-02</td>
          <td>3.527e-02</td>
          <td>1.148e-04</td>
          <td>-0.007</td>
          <td>-0.053</td>
          <td>-0.004</td>
          <td>-2.357e-02</td>
          <td>3.403e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>-2.507e-04</td>
          <td>1.228e-06</td>
          <td>-8.730e-06</td>
          <td>-2.571e-03</td>
          <td>-6.174e-06</td>
          <td>-1.767e-03</td>
          <td>-3.812e-06</td>
          <td>-5.683e-05</td>
          <td>5.614e-05</td>
          <td>1.008e-02</td>
          <td>...</td>
          <td>-1.174e-02</td>
          <td>1.241e-03</td>
          <td>1.025e-02</td>
          <td>3.112e-02</td>
          <td>-4.135e-03</td>
          <td>-0.011</td>
          <td>0.026</td>
          <td>0.009</td>
          <td>-7.311e-03</td>
          <td>9.082e-04</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>-5.261e-05</td>
          <td>-1.187e-06</td>
          <td>-1.986e-04</td>
          <td>-1.310e-04</td>
          <td>2.861e-05</td>
          <td>-4.677e-05</td>
          <td>1.092e-07</td>
          <td>-1.774e-04</td>
          <td>1.806e-04</td>
          <td>-1.008e-03</td>
          <td>...</td>
          <td>4.063e-05</td>
          <td>2.184e-02</td>
          <td>2.495e-03</td>
          <td>-8.831e-02</td>
          <td>1.660e-02</td>
          <td>0.030</td>
          <td>-0.100</td>
          <td>0.022</td>
          <td>-2.547e-02</td>
          <td>5.581e-03</td>
        </tr>
        <tr>
          <th rowspan="6" valign="top">3</th>
          <th>t1</th>
          <td>-6.169e-03</td>
          <td>7.911e-06</td>
          <td>-2.157e-04</td>
          <td>-6.310e-02</td>
          <td>-1.896e-04</td>
          <td>-4.617e-02</td>
          <td>-9.580e-05</td>
          <td>-1.466e-03</td>
          <td>1.704e-03</td>
          <td>2.690e-01</td>
          <td>...</td>
          <td>2.695e-02</td>
          <td>-6.243e-02</td>
          <td>-6.576e-03</td>
          <td>2.369e-01</td>
          <td>-3.571e-02</td>
          <td>-0.065</td>
          <td>0.135</td>
          <td>-0.079</td>
          <td>5.365e-02</td>
          <td>2.056e-02</td>
        </tr>
        <tr>
          <th>t2</th>
          <td>2.295e-04</td>
          <td>6.255e-03</td>
          <td>2.639e-02</td>
          <td>6.021e-04</td>
          <td>2.805e-04</td>
          <td>1.856e-04</td>
          <td>7.132e-05</td>
          <td>-3.892e-02</td>
          <td>-3.453e-02</td>
          <td>1.452e-02</td>
          <td>...</td>
          <td>-1.994e-01</td>
          <td>-3.102e-01</td>
          <td>-1.168e-01</td>
          <td>-1.054e-01</td>
          <td>-4.058e-02</td>
          <td>0.023</td>
          <td>0.200</td>
          <td>0.023</td>
          <td>-1.385e-02</td>
          <td>-1.724e-01</td>
        </tr>
        <tr>
          <th>t3</th>
          <td>-8.457e-04</td>
          <td>1.512e-03</td>
          <td>7.034e-03</td>
          <td>-8.137e-03</td>
          <td>1.386e-03</td>
          <td>-6.198e-03</td>
          <td>-1.003e-04</td>
          <td>-9.286e-03</td>
          <td>-7.856e-03</td>
          <td>3.749e-02</td>
          <td>...</td>
          <td>-4.493e-02</td>
          <td>1.259e-01</td>
          <td>2.542e-01</td>
          <td>1.777e-01</td>
          <td>-1.024e-03</td>
          <td>-0.042</td>
          <td>-0.038</td>
          <td>0.262</td>
          <td>-2.140e-01</td>
          <td>1.467e-01</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>-8.883e-06</td>
          <td>-2.240e-04</td>
          <td>-1.036e-03</td>
          <td>-5.241e-05</td>
          <td>-6.649e-05</td>
          <td>-6.665e-06</td>
          <td>-2.507e-06</td>
          <td>1.321e-03</td>
          <td>1.174e-03</td>
          <td>-5.719e-04</td>
          <td>...</td>
          <td>-3.039e-02</td>
          <td>2.256e-04</td>
          <td>2.894e-02</td>
          <td>3.716e-02</td>
          <td>-4.793e-04</td>
          <td>-0.008</td>
          <td>-0.047</td>
          <td>-0.006</td>
          <td>-2.042e-02</td>
          <td>3.308e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>-2.507e-04</td>
          <td>1.229e-06</td>
          <td>-8.736e-06</td>
          <td>-2.571e-03</td>
          <td>-6.175e-06</td>
          <td>-1.767e-03</td>
          <td>-3.812e-06</td>
          <td>-5.683e-05</td>
          <td>5.614e-05</td>
          <td>1.008e-02</td>
          <td>...</td>
          <td>-1.174e-02</td>
          <td>1.238e-03</td>
          <td>1.025e-02</td>
          <td>3.111e-02</td>
          <td>-4.135e-03</td>
          <td>-0.011</td>
          <td>0.026</td>
          <td>0.009</td>
          <td>-7.308e-03</td>
          <td>9.064e-04</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>-4.657e-05</td>
          <td>2.289e-06</td>
          <td>-8.570e-06</td>
          <td>-2.151e-05</td>
          <td>-8.187e-06</td>
          <td>1.310e-05</td>
          <td>1.439e-07</td>
          <td>-1.571e-04</td>
          <td>1.600e-04</td>
          <td>-1.241e-03</td>
          <td>...</td>
          <td>6.409e-03</td>
          <td>-2.870e-02</td>
          <td>6.276e-03</td>
          <td>3.529e-02</td>
          <td>-1.277e-02</td>
          <td>-0.016</td>
          <td>0.091</td>
          <td>-0.024</td>
          <td>3.187e-02</td>
          <td>-1.809e-03</td>
        </tr>
        <tr>
          <th rowspan="6" valign="top">4</th>
          <th>t1</th>
          <td>-6.169e-03</td>
          <td>7.956e-06</td>
          <td>-2.155e-04</td>
          <td>-6.310e-02</td>
          <td>-1.906e-04</td>
          <td>-4.617e-02</td>
          <td>-9.580e-05</td>
          <td>-1.466e-03</td>
          <td>1.704e-03</td>
          <td>2.690e-01</td>
          <td>...</td>
          <td>2.726e-02</td>
          <td>-6.325e-02</td>
          <td>-6.636e-03</td>
          <td>2.366e-01</td>
          <td>-3.568e-02</td>
          <td>-0.065</td>
          <td>0.135</td>
          <td>-0.079</td>
          <td>5.359e-02</td>
          <td>2.070e-02</td>
        </tr>
        <tr>
          <th>t2</th>
          <td>1.295e-04</td>
          <td>6.253e-03</td>
          <td>2.619e-02</td>
          <td>4.726e-04</td>
          <td>3.533e-04</td>
          <td>1.577e-04</td>
          <td>7.179e-05</td>
          <td>-3.925e-02</td>
          <td>-3.419e-02</td>
          <td>1.183e-02</td>
          <td>...</td>
          <td>-1.877e-01</td>
          <td>-3.177e-01</td>
          <td>-1.326e-01</td>
          <td>-1.642e-01</td>
          <td>-3.435e-02</td>
          <td>0.035</td>
          <td>0.187</td>
          <td>0.006</td>
          <td>-6.914e-04</td>
          <td>-1.884e-01</td>
        </tr>
        <tr>
          <th>t3</th>
          <td>-3.469e-04</td>
          <td>1.510e-03</td>
          <td>7.059e-03</td>
          <td>-3.040e-03</td>
          <td>1.396e-03</td>
          <td>-2.688e-03</td>
          <td>-9.276e-05</td>
          <td>-9.173e-03</td>
          <td>-7.968e-03</td>
          <td>1.746e-02</td>
          <td>...</td>
          <td>-2.215e-02</td>
          <td>1.256e-01</td>
          <td>2.349e-01</td>
          <td>1.179e-01</td>
          <td>7.068e-03</td>
          <td>-0.020</td>
          <td>-0.089</td>
          <td>0.243</td>
          <td>-1.992e-01</td>
          <td>1.457e-01</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>-7.731e-06</td>
          <td>-2.241e-04</td>
          <td>-1.037e-03</td>
          <td>-3.841e-05</td>
          <td>-6.177e-05</td>
          <td>-1.181e-05</td>
          <td>-2.531e-06</td>
          <td>1.322e-03</td>
          <td>1.169e-03</td>
          <td>-5.122e-04</td>
          <td>...</td>
          <td>-3.086e-02</td>
          <td>-3.467e-03</td>
          <td>3.029e-02</td>
          <td>3.335e-02</td>
          <td>6.161e-04</td>
          <td>-0.007</td>
          <td>-0.058</td>
          <td>-0.003</td>
          <td>-2.650e-02</td>
          <td>3.497e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>-2.507e-04</td>
          <td>1.229e-06</td>
          <td>-8.734e-06</td>
          <td>-2.571e-03</td>
          <td>-6.171e-06</td>
          <td>-1.767e-03</td>
          <td>-3.812e-06</td>
          <td>-5.682e-05</td>
          <td>5.614e-05</td>
          <td>1.008e-02</td>
          <td>...</td>
          <td>-1.174e-02</td>
          <td>1.238e-03</td>
          <td>1.026e-02</td>
          <td>3.112e-02</td>
          <td>-4.135e-03</td>
          <td>-0.011</td>
          <td>0.026</td>
          <td>0.010</td>
          <td>-7.314e-03</td>
          <td>9.083e-04</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>-4.712e-05</td>
          <td>2.923e-07</td>
          <td>5.696e-05</td>
          <td>-2.570e-05</td>
          <td>-3.632e-06</td>
          <td>1.221e-05</td>
          <td>1.684e-07</td>
          <td>-1.412e-04</td>
          <td>1.769e-04</td>
          <td>-1.299e-03</td>
          <td>...</td>
          <td>1.457e-02</td>
          <td>-1.664e-02</td>
          <td>4.624e-03</td>
          <td>3.608e-02</td>
          <td>-1.077e-02</td>
          <td>-0.016</td>
          <td>0.083</td>
          <td>-0.026</td>
          <td>3.404e-02</td>
          <td>-3.449e-04</td>
        </tr>
        <tr>
          <th rowspan="6" valign="top">5</th>
          <th>t1</th>
          <td>-6.801e-03</td>
          <td>1.081e-05</td>
          <td>-2.253e-04</td>
          <td>-6.955e-02</td>
          <td>-2.029e-04</td>
          <td>-5.058e-02</td>
          <td>-1.054e-04</td>
          <td>-1.626e-03</td>
          <td>1.863e-03</td>
          <td>2.943e-01</td>
          <td>...</td>
          <td>-1.462e-03</td>
          <td>-4.749e-02</td>
          <td>1.290e-02</td>
          <td>2.882e-01</td>
          <td>-4.040e-02</td>
          <td>-0.084</td>
          <td>0.165</td>
          <td>-0.056</td>
          <td>3.263e-02</td>
          <td>2.358e-02</td>
        </tr>
        <tr>
          <th>t2</th>
          <td>2.553e-04</td>
          <td>6.819e-03</td>
          <td>2.910e-02</td>
          <td>8.057e-04</td>
          <td>4.970e-04</td>
          <td>2.453e-04</td>
          <td>7.785e-05</td>
          <td>-4.223e-02</td>
          <td>-3.750e-02</td>
          <td>1.563e-02</td>
          <td>...</td>
          <td>-1.560e-01</td>
          <td>-3.697e-01</td>
          <td>-2.080e-01</td>
          <td>-1.525e-01</td>
          <td>-5.946e-02</td>
          <td>0.022</td>
          <td>0.442</td>
          <td>0.012</td>
          <td>6.531e-02</td>
          <td>-2.889e-01</td>
        </tr>
        <tr>
          <th>t3</th>
          <td>-8.469e-04</td>
          <td>1.512e-03</td>
          <td>7.038e-03</td>
          <td>-8.160e-03</td>
          <td>1.385e-03</td>
          <td>-6.210e-03</td>
          <td>-1.004e-04</td>
          <td>-9.286e-03</td>
          <td>-7.856e-03</td>
          <td>3.756e-02</td>
          <td>...</td>
          <td>-4.534e-02</td>
          <td>1.271e-01</td>
          <td>2.550e-01</td>
          <td>1.792e-01</td>
          <td>-1.136e-03</td>
          <td>-0.042</td>
          <td>-0.037</td>
          <td>0.263</td>
          <td>-2.141e-01</td>
          <td>1.472e-01</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>-8.399e-06</td>
          <td>-2.241e-04</td>
          <td>-1.035e-03</td>
          <td>-4.509e-05</td>
          <td>-6.317e-05</td>
          <td>-9.635e-06</td>
          <td>-2.518e-06</td>
          <td>1.322e-03</td>
          <td>1.172e-03</td>
          <td>-5.433e-04</td>
          <td>...</td>
          <td>-3.061e-02</td>
          <td>-9.825e-04</td>
          <td>2.993e-02</td>
          <td>3.527e-02</td>
          <td>1.148e-04</td>
          <td>-0.007</td>
          <td>-0.053</td>
          <td>-0.004</td>
          <td>-2.357e-02</td>
          <td>3.403e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>-2.507e-04</td>
          <td>1.228e-06</td>
          <td>-8.730e-06</td>
          <td>-2.571e-03</td>
          <td>-6.174e-06</td>
          <td>-1.767e-03</td>
          <td>-3.812e-06</td>
          <td>-5.683e-05</td>
          <td>5.614e-05</td>
          <td>1.008e-02</td>
          <td>...</td>
          <td>-1.174e-02</td>
          <td>1.241e-03</td>
          <td>1.025e-02</td>
          <td>3.112e-02</td>
          <td>-4.135e-03</td>
          <td>-0.011</td>
          <td>0.026</td>
          <td>0.009</td>
          <td>-7.311e-03</td>
          <td>9.082e-04</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>-5.261e-05</td>
          <td>-1.187e-06</td>
          <td>-1.986e-04</td>
          <td>-1.310e-04</td>
          <td>2.861e-05</td>
          <td>-4.677e-05</td>
          <td>1.092e-07</td>
          <td>-1.774e-04</td>
          <td>1.806e-04</td>
          <td>-1.008e-03</td>
          <td>...</td>
          <td>4.063e-05</td>
          <td>2.184e-02</td>
          <td>2.495e-03</td>
          <td>-8.831e-02</td>
          <td>1.660e-02</td>
          <td>0.030</td>
          <td>-0.100</td>
          <td>0.022</td>
          <td>-2.547e-02</td>
          <td>5.581e-03</td>
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
          <td>7.413e-05</td>
          <td>-8.245e-05</td>
          <td>-3.908e-04</td>
          <td>3.482e-03</td>
          <td>-3.748e-05</td>
          <td>2.988e-04</td>
          <td>2.694e-06</td>
          <td>2.440e-05</td>
          <td>1.075e-03</td>
          <td>-2.721e-03</td>
          <td>...</td>
          <td>4.755e-04</td>
          <td>5.261e-03</td>
          <td>-4.662e-02</td>
          <td>-5.886e-02</td>
          <td>-6.227e-03</td>
          <td>0.017</td>
          <td>0.162</td>
          <td>-0.557</td>
          <td>-6.614e-01</td>
          <td>1.042e-01</td>
        </tr>
        <tr>
          <th>t2</th>
          <td>4.452e-05</td>
          <td>-2.089e-04</td>
          <td>-5.166e-03</td>
          <td>2.748e-04</td>
          <td>1.754e-04</td>
          <td>4.173e-04</td>
          <td>3.617e-06</td>
          <td>-1.361e-04</td>
          <td>1.100e-04</td>
          <td>-9.064e-04</td>
          <td>...</td>
          <td>2.047e-01</td>
          <td>6.117e-02</td>
          <td>-5.444e-02</td>
          <td>1.529e-02</td>
          <td>-1.469e-02</td>
          <td>-0.023</td>
          <td>0.183</td>
          <td>0.290</td>
          <td>3.938e-01</td>
          <td>-3.587e-01</td>
        </tr>
        <tr>
          <th>t3</th>
          <td>1.283e-04</td>
          <td>1.048e-03</td>
          <td>8.983e-03</td>
          <td>5.709e-04</td>
          <td>1.808e-04</td>
          <td>1.258e-03</td>
          <td>-5.577e-06</td>
          <td>-5.649e-03</td>
          <td>-5.097e-03</td>
          <td>-7.903e-03</td>
          <td>...</td>
          <td>-1.857e-01</td>
          <td>-2.785e-02</td>
          <td>6.353e-02</td>
          <td>-5.410e-02</td>
          <td>2.473e-02</td>
          <td>0.030</td>
          <td>-0.234</td>
          <td>0.069</td>
          <td>-6.092e-02</td>
          <td>3.214e-01</td>
        </tr>
        <tr>
          <th>r1</th>
          <td>3.005e-07</td>
          <td>5.476e-05</td>
          <td>6.343e-04</td>
          <td>6.334e-06</td>
          <td>-2.493e-06</td>
          <td>2.715e-06</td>
          <td>-5.464e-07</td>
          <td>-2.376e-04</td>
          <td>-2.019e-04</td>
          <td>6.017e-05</td>
          <td>...</td>
          <td>-4.279e-04</td>
          <td>-3.524e-03</td>
          <td>9.711e-04</td>
          <td>6.896e-03</td>
          <td>9.867e-04</td>
          <td>-0.001</td>
          <td>-0.014</td>
          <td>-0.008</td>
          <td>-2.789e-02</td>
          <td>2.645e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>-1.195e-05</td>
          <td>-1.468e-05</td>
          <td>-9.874e-05</td>
          <td>2.889e-07</td>
          <td>-7.292e-06</td>
          <td>-1.234e-04</td>
          <td>-1.826e-07</td>
          <td>7.492e-05</td>
          <td>1.152e-04</td>
          <td>9.511e-04</td>
          <td>...</td>
          <td>2.369e-02</td>
          <td>1.095e-03</td>
          <td>-8.118e-03</td>
          <td>1.289e-02</td>
          <td>-1.785e-03</td>
          <td>-0.005</td>
          <td>0.015</td>
          <td>-0.021</td>
          <td>-3.344e-02</td>
          <td>-3.849e-03</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>2.865e-06</td>
          <td>1.522e-05</td>
          <td>6.913e-05</td>
          <td>-4.280e-06</td>
          <td>4.744e-06</td>
          <td>2.949e-05</td>
          <td>-1.857e-07</td>
          <td>-1.044e-04</td>
          <td>-6.757e-05</td>
          <td>1.044e-04</td>
          <td>...</td>
          <td>2.703e-02</td>
          <td>1.362e-03</td>
          <td>-5.113e-03</td>
          <td>1.492e-02</td>
          <td>-6.335e-04</td>
          <td>-0.005</td>
          <td>0.001</td>
          <td>0.027</td>
          <td>1.418e-02</td>
          <td>-1.118e-02</td>
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
          <td>1.815e-05</td>
          <td>-9.454e-05</td>
          <td>-3.224e-04</td>
          <td>-3.568e-05</td>
          <td>-1.340e-05</td>
          <td>-3.384e-05</td>
          <td>-1.329e-06</td>
          <td>7.127e-04</td>
          <td>4.621e-04</td>
          <td>6.393e-04</td>
          <td>...</td>
          <td>-3.555e-02</td>
          <td>-8.501e-03</td>
          <td>1.420e-02</td>
          <td>1.374e-02</td>
          <td>1.390e-04</td>
          <td>-0.004</td>
          <td>-0.017</td>
          <td>0.004</td>
          <td>2.479e-04</td>
          <td>2.458e-03</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>1.174e-04</td>
          <td>8.335e-07</td>
          <td>-1.801e-05</td>
          <td>1.328e-03</td>
          <td>-2.449e-05</td>
          <td>7.252e-04</td>
          <td>3.178e-07</td>
          <td>-1.708e-05</td>
          <td>-1.350e-05</td>
          <td>-3.866e-03</td>
          <td>...</td>
          <td>6.104e-03</td>
          <td>-6.432e-03</td>
          <td>1.082e-02</td>
          <td>3.451e-02</td>
          <td>-4.203e-03</td>
          <td>-0.011</td>
          <td>0.018</td>
          <td>-0.003</td>
          <td>-7.885e-03</td>
          <td>-2.321e-02</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>-1.512e-05</td>
          <td>3.817e-05</td>
          <td>2.898e-04</td>
          <td>-7.733e-06</td>
          <td>1.063e-06</td>
          <td>-1.915e-06</td>
          <td>6.212e-07</td>
          <td>-2.275e-04</td>
          <td>-1.247e-04</td>
          <td>-5.015e-04</td>
          <td>...</td>
          <td>1.508e-02</td>
          <td>-1.308e-03</td>
          <td>3.008e-03</td>
          <td>9.727e-03</td>
          <td>3.836e-04</td>
          <td>-0.001</td>
          <td>-0.003</td>
          <td>-0.030</td>
          <td>3.103e-02</td>
          <td>2.712e-02</td>
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
          <td>9.862e-07</td>
          <td>5.862e-05</td>
          <td>5.580e-04</td>
          <td>1.046e-05</td>
          <td>-6.905e-05</td>
          <td>5.601e-06</td>
          <td>1.679e-06</td>
          <td>-2.394e-04</td>
          <td>-2.043e-04</td>
          <td>3.883e-05</td>
          <td>...</td>
          <td>1.138e-03</td>
          <td>-1.261e-02</td>
          <td>1.119e-02</td>
          <td>1.439e-02</td>
          <td>1.245e-03</td>
          <td>-0.004</td>
          <td>-0.024</td>
          <td>-0.012</td>
          <td>-4.048e-03</td>
          <td>-1.465e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>-8.388e-06</td>
          <td>-1.919e-06</td>
          <td>-7.634e-06</td>
          <td>-2.048e-04</td>
          <td>1.955e-07</td>
          <td>-2.855e-04</td>
          <td>-5.311e-07</td>
          <td>6.254e-05</td>
          <td>-5.671e-05</td>
          <td>2.168e-03</td>
          <td>...</td>
          <td>-2.994e-02</td>
          <td>-4.564e-03</td>
          <td>1.167e-02</td>
          <td>1.208e-02</td>
          <td>-2.319e-03</td>
          <td>-0.004</td>
          <td>0.012</td>
          <td>-0.005</td>
          <td>5.452e-03</td>
          <td>5.399e-03</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>4.235e-05</td>
          <td>3.105e-06</td>
          <td>1.133e-06</td>
          <td>3.700e-04</td>
          <td>-3.676e-07</td>
          <td>2.318e-04</td>
          <td>3.299e-07</td>
          <td>-1.454e-05</td>
          <td>-9.195e-06</td>
          <td>-8.160e-04</td>
          <td>...</td>
          <td>7.604e-03</td>
          <td>-3.327e-03</td>
          <td>1.359e-02</td>
          <td>-8.851e-04</td>
          <td>-7.085e-04</td>
          <td>-0.002</td>
          <td>0.011</td>
          <td>-0.018</td>
          <td>1.781e-02</td>
          <td>1.326e-02</td>
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
          <td>-1.756e-05</td>
          <td>-9.628e-05</td>
          <td>-3.117e-04</td>
          <td>4.014e-05</td>
          <td>-1.268e-05</td>
          <td>3.502e-05</td>
          <td>-1.054e-06</td>
          <td>5.821e-04</td>
          <td>6.400e-04</td>
          <td>-9.421e-04</td>
          <td>...</td>
          <td>3.064e-02</td>
          <td>-2.242e-03</td>
          <td>5.440e-04</td>
          <td>1.809e-02</td>
          <td>1.961e-04</td>
          <td>-0.005</td>
          <td>-0.016</td>
          <td>0.004</td>
          <td>-1.559e-04</td>
          <td>4.253e-03</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>1.170e-04</td>
          <td>-2.698e-07</td>
          <td>2.597e-05</td>
          <td>1.325e-03</td>
          <td>3.278e-05</td>
          <td>7.228e-04</td>
          <td>2.756e-06</td>
          <td>1.174e-05</td>
          <td>1.113e-05</td>
          <td>-3.858e-03</td>
          <td>...</td>
          <td>1.025e-02</td>
          <td>4.245e-03</td>
          <td>2.363e-03</td>
          <td>2.781e-02</td>
          <td>-4.249e-03</td>
          <td>-0.009</td>
          <td>0.026</td>
          <td>0.024</td>
          <td>-1.089e-02</td>
          <td>1.333e-02</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>-1.548e-05</td>
          <td>-4.294e-05</td>
          <td>-2.770e-04</td>
          <td>-1.258e-05</td>
          <td>3.928e-06</td>
          <td>-5.063e-06</td>
          <td>-3.048e-07</td>
          <td>1.836e-04</td>
          <td>2.254e-04</td>
          <td>-5.880e-04</td>
          <td>...</td>
          <td>2.334e-02</td>
          <td>-1.135e-03</td>
          <td>5.283e-03</td>
          <td>1.865e-03</td>
          <td>-3.915e-03</td>
          <td>-0.002</td>
          <td>0.024</td>
          <td>-0.032</td>
          <td>2.892e-02</td>
          <td>1.447e-02</td>
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
          <td>3.006e-07</td>
          <td>5.476e-05</td>
          <td>6.343e-04</td>
          <td>6.336e-06</td>
          <td>-2.493e-06</td>
          <td>2.716e-06</td>
          <td>-5.464e-07</td>
          <td>-2.376e-04</td>
          <td>-2.019e-04</td>
          <td>6.017e-05</td>
          <td>...</td>
          <td>-4.280e-04</td>
          <td>-3.524e-03</td>
          <td>9.711e-04</td>
          <td>6.896e-03</td>
          <td>9.867e-04</td>
          <td>-0.001</td>
          <td>-0.014</td>
          <td>-0.008</td>
          <td>-2.789e-02</td>
          <td>2.645e-02</td>
        </tr>
        <tr>
          <th>r2</th>
          <td>1.723e-06</td>
          <td>1.278e-06</td>
          <td>-1.805e-06</td>
          <td>1.940e-04</td>
          <td>-3.377e-07</td>
          <td>-8.446e-06</td>
          <td>-3.548e-08</td>
          <td>-4.728e-05</td>
          <td>4.650e-05</td>
          <td>2.129e-04</td>
          <td>...</td>
          <td>2.084e-02</td>
          <td>1.171e-03</td>
          <td>-6.235e-03</td>
          <td>1.349e-02</td>
          <td>-1.096e-03</td>
          <td>-0.005</td>
          <td>0.006</td>
          <td>0.005</td>
          <td>-4.639e-03</td>
          <td>-6.872e-03</td>
        </tr>
        <tr>
          <th>r3</th>
          <td>-7.271e-06</td>
          <td>3.394e-06</td>
          <td>-2.717e-06</td>
          <td>-1.478e-04</td>
          <td>-4.097e-07</td>
          <td>-5.573e-05</td>
          <td>-2.948e-07</td>
          <td>-1.383e-05</td>
          <td>-1.663e-05</td>
          <td>6.515e-04</td>
          <td>...</td>
          <td>2.914e-02</td>
          <td>1.305e-03</td>
          <td>-6.509e-03</td>
          <td>1.448e-02</td>
          <td>-1.144e-03</td>
          <td>-0.005</td>
          <td>0.008</td>
          <td>0.008</td>
          <td>-7.160e-03</td>
          <td>-8.942e-03</td>
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
    plate_stress = dict_keys(['element_type', 'element_name', 'nonlinear_factor', '_times', 'result_name', 'approach_code', 'analysis_code', 'data', 'isubcase', 'ogs', 'pval_step', 'name', 'superelement_adaptivity_index', '_count', 'is_built', 'format_code', 'sort_code', 'table_code', 'title', 'subtitle', 'label', 'num_wide', 'device_code', 'table_name', 'data_frame', 'dt', 'ntimes', 'ntotal', '_ntotals', 'load_as_h5', 'h5_file', 'data_code', 'ielement', 'nelements', 'nnodes', '_encoding', '_times_dtype', 'cycle', 'data_names', 'eign', 'is_msc', 'is_strain_flag', 'is_stress_flag', 'load_set', 'mode', 'mode2', 's_code', 'sort_bits', 'sort_method', 'stress_bits', 'subtitle_original', 'tCode', 'thermal', 'thermal_bits', 'modes', 'eigns', 'mode2s', 'cycles', 'itotal', 'itime', 'element_node', 'words'])
    
    data_code_keys = dict_keys(['_encoding', 'load_as_h5', 'is_msc', 'table_name', 'approach_code', 'isubcase', 'table_code', 'tCode', 'sort_code', 'sort_method', 'device_code', 'analysis_code', 'sort_bits', 'element_type', 'load_set', 'format_code', 'num_wide', 's_code', 'thermal', 'nonlinear_factor', 'name', 'mode', 'eign', 'mode2', 'cycle', 'data_names', '_times_dtype', 'thermal_bits', 'element_name', 'subtitle', 'subtitle_original', 'pval_step', 'superelement_adaptivity_index', 'label', 'title', 'stress_bits', 'is_stress_flag', 'is_strain_flag', 'result_name', '_count'])
    
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

    methods = ['add_new_eid_sort1', 'add_new_node_sort1', 'add_sort1', 'apply_data_code', 'approach_code_str', 'build', 'build_dataframe', 'cast_grid_type', 'code_information', 'eid_to_element_node_index', 'export_to_hdf5', 'finalize', 'get_data_code', 'get_element_index', 'get_element_type', 'get_headers', 'get_nnodes_bilinear', 'get_stats', 'get_unsteady_value', 'is_bilinear', 'is_magnitude_phase', 'is_sort1_new', 'is_thermal', 'object_attributes', 'object_methods', 'print_data_members', 'print_table_code', 'recast_gridtype_as_string', 'set_as_sort1', 'set_table_type', 'update_data_code', 'update_dt', 'update_t_code', 'write_f06', 'write_op2']
    
    methods2= ['add_new_eid_sort1', 'add_new_node_sort1', 'add_sort1', 'apply_data_code', 'approach_code_str', 'build', 'build_dataframe', 'cast_grid_type', 'code_information', 'eid_to_element_node_index', 'export_to_hdf5', 'finalize', 'get_data_code', 'get_element_index', 'get_element_type', 'get_headers', 'get_nnodes_bilinear', 'get_stats', 'get_unsteady_value', 'is_bilinear', 'is_magnitude_phase', 'is_sort1_new', 'is_thermal', 'print_data_members', 'print_table_code', 'recast_gridtype_as_string', 'set_as_sort1', 'set_table_type', 'update_data_code', 'update_dt', 'update_t_code', 'write_f06', 'write_op2']
    
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

   print("is_fiber_distance = %s" % plate_stress.is_fiber_distance())

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
    we have von mises data; ovm=[54.222  5.041 13.143 ...  2.34   6.146  7.368]
    
    [layer1, layer2, ...] = [54.222  5.041 13.143 21.222 78.545 17.91 ]
    ieid1000 = [1998 1999]
    ovm_mode6_eid1000 = [90.618 94.091] -> 94.09056
    

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
    ps.cycles = 20.548073657198046
    oxx = -18.872536
    oyy = -20.16303
    txy = -8.309847
    omax = -11.182922
    omin = -27.852644
    ovm/max_shear = 24.276606
    fiber_distance = -0.4
    

.. code:: ipython3

    from pyNastran.bdf.bdf import read_bdf
    bdf_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_4pt.dat'))
    model = read_bdf(bdf_filename, debug=False)
    mass, cg, I = model.mass_properties()



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


Let’s print out the actual mass properties from the OP2 and get the same result as the F06
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We need ``PARAM,POSTEXT,YES`` in out BDF to get the Grid Point Weight
Table

.. code:: ipython3

    gpw = op2.grid_point_weight
    #print(gpw.object_attributes())
    
    print(gpw)
    gpw.object_methods()
    #gpw.write_f06?
    print(gpw.get_stats())


.. parsed-literal::

    
    
    

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

    1    ISAT_SM_LAUNCH_4PT MODES TO 400 HZ                                    FEBRUARY  14, 2018  pyNastran v1.2.0       PAGE     1
         
    0                                                                                                            SUBCASE 1
                                             R E A L   E I G E N V E C T O R   N O .          1
     
          POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
                 1      G     -5.547863E-03   2.133077E-04  -8.469186E-04  -8.399206E-06  -2.506956E-04  -5.261146E-05
                 2      G     -5.547863E-03   1.080848E-04  -3.455275E-04  -8.399206E-06  -2.506956E-04  -5.261146E-05
                 3      G     -6.169366E-03   2.295251E-04  -8.457433E-04  -8.882780E-06  -2.506924E-04  -4.657186E-05
                 4      G     -6.169462E-03   1.295465E-04  -3.468718E-04  -7.731095E-06  -2.506902E-04  -4.712092E-05
                 5      G     -6.801341E-03   2.553037E-04  -8.469186E-04  -8.399206E-06  -2.506956E-04  -5.261146E-05
                 6      G     -6.801341E-03   1.500808E-04  -3.455275E-04  -8.399206E-06  -2.506956E-04  -5.261146E-05
                 7      G     -7.420456E-03   2.779656E-04  -8.458295E-04  -9.019739E-06  -2.506924E-04  -6.147318E-05
                 8      G     -7.420568E-03   1.595261E-04  -3.468334E-04  -7.642398E-06  -2.506902E-04  -6.249677E-05
                 9      G     -8.054819E-03   2.972998E-04  -8.469186E-04  -8.399206E-06  -2.506956E-04  -5.261146E-05
                10      G     -8.054819E-03   1.920769E-04  -3.455275E-04  -8.399206E-06  -2.506956E-04  -5.261146E-05
                11      G     -5.547945E-03  -2.034865E-04   7.615836E-04   9.685779E-06  -2.505747E-04  -5.059179E-05
                12      G     -6.170889E-03  -2.232806E-04   7.602651E-04   1.005450E-05  -2.505695E-04  -4.706052E-05
                13      G     -6.800818E-03  -2.519154E-04   7.615836E-04   9.685779E-06  -2.505747E-04  -5.059179E-05
                14      G     -7.421867E-03  -2.748493E-04   7.598438E-04   1.041645E-05  -2.505710E-04  -6.034294E-05
                15      G     -8.053692E-03  -3.003443E-04   7.615836E-04   9.685779E-06  -2.505747E-04  -5.059179E-05
                16      G     -5.547945E-03  -1.023029E-04   2.604342E-04   9.685779E-06  -2.505747E-04  -5.059179E-05
                17      G     -6.170872E-03  -1.254381E-04   2.618640E-04   9.119666E-06  -2.505694E-04  -4.744373E-05
                18      G     -6.800818E-03  -1.507318E-04   2.604342E-04   9.685779E-06  -2.505747E-04  -5.059179E-05
                19      G     -7.421871E-03  -1.645201E-04   2.621645E-04   8.841090E-06  -2.505694E-04  -6.108151E-05
                20      G     -8.053692E-03  -1.991607E-04   2.604342E-04   9.685779E-06  -2.505747E-04  -5.059179E-05
                21      G     -2.821161E-03  -1.092790E-04   7.579715E-04   7.133808E-06  -2.468418E-04  -3.410970E-05
                22      G     -3.443996E-03  -1.233896E-04   7.561403E-04   7.799596E-06  -2.468382E-04  -3.001459E-05
                23      G     -4.055370E-03  -1.449481E-04   7.579715E-04   7.133808E-06  -2.468418E-04  -3.410970E-05
                24      G     -4.670648E-03  -1.632725E-04   7.558515E-04   7.958040E-06  -2.468395E-04  -4.140726E-05
                25      G     -5.289579E-03  -1.806171E-04   7.579715E-04   7.133808E-06  -2.468418E-04  -3.410970E-05
                26      G     -2.821161E-03  -4.105964E-05   2.642879E-04   7.133808E-06  -2.468418E-04  -3.410970E-05
                27      G     -3.443890E-03  -6.412427E-05   2.660819E-04   6.400805E-06  -2.468367E-04  -2.950787E-05
                28      G     -4.055370E-03  -7.672868E-05   2.642879E-04   7.133808E-06  -2.468418E-04  -3.410970E-05
                29      G     -4.670699E-03  -9.216915E-05   2.663240E-04   6.263966E-06  -2.468367E-04  -4.123877E-05
                30      G     -5.289579E-03  -1.123977E-04   2.642879E-04   7.133808E-06  -2.468418E-04  -3.410970E-05
                31      G     -1.459303E-04   2.173024E-05   2.685585E-04   7.402167E-06  -2.425398E-04  -1.785639E-05
                32      G     -1.459303E-04  -1.398253E-05   7.536381E-04   7.402167E-06  -2.425398E-04  -1.785639E-05
                33      G     -7.555047E-04  -8.368386E-06   2.703055E-04   6.154360E-06  -2.425350E-04  -8.907110E-06
                34      G     -7.555627E-04  -2.869231E-05   7.518338E-04   8.680357E-06  -2.425363E-04  -9.971128E-06
    

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

    print('cg =\n%s' % gpw.cg)
    print('cg = %s' % cg)


.. parsed-literal::

    cg =
    None
    cg = [ -0.034  -2.531 -18.468]
    

It’s not like Nastran is perfect either.
----------------------------------------

Limitations
~~~~~~~~~~~

1. You cannot do weight statements in Nastran by
   component/property/material.

2. Everything is always summmed up (e.g. you can have different geometry
   in Subcase 2 and MPCs connecting physical geomtry, with other parts
   flying off into space).

These are things that pyNastran ``can`` do.

.. code:: ipython3

    from pyNastran.bdf.bdf import read_bdf
    bdf_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_4pt.dat'))
    model = read_bdf(bdf_filename, debug=False)



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


Weight Statement
~~~~~~~~~~~~~~~~

Let’s get the breakdown by property ID

.. code:: ipython3

    #help(model.mass_properties)
    
    pid_to_eids_map = model.get_element_ids_dict_with_pids()
    #print(pid_to_eids_map.keys())
    print('pid, mass, cg, [ixx, iyy, izz, ixy, ixz, iyz]')
    for pid, eids in sorted(pid_to_eids_map.items()):
        mass, cg, inertia = model.mass_properties(element_ids=eids, mass_ids=[], reference_point=[0., 0., 0.])
        print('%-6s %-.6f %-38s %s' % (pid, mass, cg, inertia))
    
    mass_ids = list(model.masses.keys())
    mass, cg, inertia = model.mass_properties(element_ids=[], mass_ids=mass_ids, reference_point=[0., 0., 0.])
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
    
