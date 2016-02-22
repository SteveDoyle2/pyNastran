
OP2 Demo
========

Most people are comfortable with the F06. However, it's: - a lot harder
to parse - much, much, much slower - much, much more memory inefficient

The pyNastran OP2 Reader is fast, highly validated, and it supports most
result types.

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

    INFO:      fname=op2_scalar.pyc            lineNo=1020   op2_filename = 'f:\\work\\pynastran\\pynastran\\master3\\models\\iSat\\ISat_Launch_Sm_4pt.op2'
    

OP2 Introspection
-----------------

The ``get_op2_stats()`` function lets you quickly understand what in an
op2.

.. code:: python

    print(op2.get_op2_stats())


.. parsed-literal::

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
    
    ctria3_stress[1]
      type=RealPlateStressArray ntimes=167 nelements=32 nnodes_per_element=1 nlayers=2 ntotal=64
      data: [ntimes, ntotal, 8] where 8=[fiber_distance, oxx, oyy, txy, angle, omax, omin, von_mises]
      data.shape=(167L, 64L, 8L)
      element type: CTRIA3
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
    
    print(eig1.data_frame)


.. parsed-literal::

    Mode               1          2          3          4          5          6    \
    Freq           8.358      9.507      15.663     20.229     20.306     20.548    
    NodeID Item                                                                     
    1      t1    5.548e-03  4.671e-06 -1.818e-04 -5.670e-02  1.722e-04 -4.175e-02   
           t2   -2.133e-04  5.699e-03  2.392e-02  5.801e-04 -1.812e-04  1.971e-04   
           t3    8.469e-04  1.512e-03  7.038e-03 -8.160e-03 -1.385e-03 -6.209e-03   
           r1    8.399e-06 -2.241e-04 -1.035e-03 -4.509e-05  6.317e-05 -9.634e-06   
           r2    2.507e-04  1.228e-06 -8.742e-06 -2.571e-03  6.180e-06 -1.767e-03   
           r3    5.261e-05 -1.187e-06 -1.986e-04 -1.310e-04 -2.860e-05 -4.676e-05   
    2      t1    5.548e-03  4.671e-06 -1.818e-04 -5.670e-02  1.722e-04 -4.175e-02   
           t2   -1.081e-04  5.696e-03  2.353e-02  3.180e-04 -2.384e-04  1.036e-04   
           t3    3.455e-04  1.510e-03  7.055e-03 -3.018e-03 -1.398e-03 -2.676e-03   
           r1    8.399e-06 -2.241e-04 -1.035e-03 -4.509e-05  6.317e-05 -9.634e-06   
           r2    2.507e-04  1.228e-06 -8.742e-06 -2.571e-03  6.180e-06 -1.767e-03   
           r3    5.261e-05 -1.187e-06 -1.986e-04 -1.310e-04 -2.860e-05 -4.676e-05   
    3      t1    6.169e-03  7.911e-06 -2.160e-04 -6.310e-02  1.897e-04 -4.617e-02   
           t2   -2.295e-04  6.255e-03  2.639e-02  6.019e-04 -2.806e-04  1.856e-04   
           t3    8.457e-04  1.512e-03  7.034e-03 -8.138e-03 -1.386e-03 -6.198e-03   
           r1    8.883e-06 -2.240e-04 -1.036e-03 -5.241e-05  6.649e-05 -6.664e-06   
           r2    2.507e-04  1.229e-06 -8.748e-06 -2.571e-03  6.181e-06 -1.767e-03   
           r3    4.657e-05  2.289e-06 -8.563e-06 -2.151e-05  8.189e-06  1.310e-05   
    4      t1    6.169e-03  7.956e-06 -2.157e-04 -6.310e-02  1.907e-04 -4.617e-02   
           t2   -1.295e-04  6.253e-03  2.619e-02  4.724e-04 -3.533e-04  1.577e-04   
           t3    3.469e-04  1.510e-03  7.059e-03 -3.040e-03 -1.396e-03 -2.688e-03   
           r1    7.731e-06 -2.241e-04 -1.037e-03 -3.840e-05  6.177e-05 -1.181e-05   
           r2    2.507e-04  1.229e-06 -8.746e-06 -2.571e-03  6.177e-06 -1.767e-03   
           r3    4.712e-05  2.923e-07  5.697e-05 -2.570e-05  3.632e-06  1.221e-05   
    5      t1    6.801e-03  1.081e-05 -2.255e-04 -6.955e-02  2.031e-04 -5.058e-02   
           t2   -2.553e-04  6.819e-03  2.910e-02  8.055e-04 -4.971e-04  2.453e-04   
           t3    8.469e-04  1.512e-03  7.038e-03 -8.160e-03 -1.385e-03 -6.209e-03   
           r1    8.399e-06 -2.241e-04 -1.035e-03 -4.509e-05  6.317e-05 -9.634e-06   
           r2    2.507e-04  1.228e-06 -8.742e-06 -2.571e-03  6.180e-06 -1.767e-03   
           r3    5.261e-05 -1.187e-06 -1.986e-04 -1.310e-04 -2.860e-05 -4.676e-05   
    ...                ...        ...        ...        ...        ...        ...   
    5629   t1   -7.413e-05 -8.245e-05 -3.907e-04  3.482e-03  3.748e-05  2.988e-04   
           t2   -4.452e-05 -2.089e-04 -5.165e-03  2.748e-04 -1.754e-04  4.173e-04   
           t3   -1.283e-04  1.048e-03  8.982e-03  5.709e-04 -1.808e-04  1.258e-03   
           r1   -3.005e-07  5.476e-05  6.343e-04  6.332e-06  2.491e-06  2.715e-06   
           r2    1.195e-05 -1.468e-05 -9.874e-05  2.887e-07  7.293e-06 -1.234e-04   
           r3   -2.865e-06  1.522e-05  6.912e-05 -4.279e-06 -4.743e-06  2.949e-05   
    5630   t1    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00   
           t2    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00   
           t3    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00   
           r1   -1.815e-05 -9.454e-05 -3.223e-04 -3.568e-05  1.340e-05 -3.384e-05   
           r2   -1.174e-04  8.335e-07 -1.801e-05  1.328e-03  2.448e-05  7.252e-04   
           r3    1.512e-05  3.817e-05  2.898e-04 -7.734e-06 -1.064e-06 -1.914e-06   
    5631   t1    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00   
           t2    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00   
           t3    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00   
           r1   -9.862e-07  5.862e-05  5.579e-04  1.046e-05  6.905e-05  5.601e-06   
           r2    8.388e-06 -1.919e-06 -7.635e-06 -2.048e-04 -1.957e-07 -2.855e-04   
           r3   -4.235e-05  3.105e-06  1.132e-06  3.700e-04  3.678e-07  2.318e-04   
    5632   t1    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00   
           t2    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00   
           t3    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00   
           r1    1.756e-05 -9.628e-05 -3.117e-04  4.014e-05  1.268e-05  3.502e-05   
           r2   -1.170e-04 -2.698e-07  2.598e-05  1.325e-03 -3.278e-05  7.227e-04   
           r3    1.548e-05 -4.294e-05 -2.770e-04 -1.257e-05 -3.928e-06 -5.062e-06   
    5633   t1    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00   
           t2    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00   
           t3    0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00  0.000e+00   
           r1   -3.006e-07  5.476e-05  6.343e-04  6.334e-06  2.491e-06  2.716e-06   
           r2   -1.723e-06  1.278e-06 -1.805e-06  1.940e-04  3.380e-07 -8.450e-06   
           r3    7.271e-06  3.394e-06 -2.722e-06 -1.478e-04  4.108e-07 -5.572e-05   
    
    Mode               7          8          9          10     ...            158  \
    Freq           21.500     21.701     21.716     28.440     ...        382.715   
    NodeID Item                                                ...                  
    1      t1    8.632e-05 -1.341e-03  1.582e-03 -2.428e-01    ...      5.723e-02   
           t2   -6.526e-05 -3.563e-02 -3.164e-02 -1.292e-02    ...     -3.090e-01   
           t3    1.004e-04 -9.286e-03 -7.856e-03 -3.743e-02    ...     -4.535e-02   
           r1    2.518e-06  1.322e-03  1.172e-03  5.440e-04    ...     -3.061e-02   
           r2    3.812e-06 -5.683e-05  5.614e-05 -1.004e-02    ...     -1.174e-02   
           r3   -1.092e-07 -1.774e-04  1.806e-04  1.011e-03    ...      4.109e-05   
    2      t1    8.632e-05 -1.341e-03  1.582e-03 -2.428e-01    ...      5.723e-02   
           t2   -6.548e-05 -3.598e-02 -3.128e-02 -1.090e-02    ...     -3.090e-01   
           t3    9.274e-05 -9.172e-03 -7.968e-03 -1.735e-02    ...     -2.187e-02   
           r1    2.518e-06  1.322e-03  1.172e-03  5.440e-04    ...     -3.061e-02   
           r2    3.812e-06 -5.683e-05  5.614e-05 -1.004e-02    ...     -1.174e-02   
           r3   -1.092e-07 -1.774e-04  1.806e-04  1.011e-03    ...      4.109e-05   
    3      t1    9.580e-05 -1.466e-03  1.704e-03 -2.679e-01    ...      2.695e-02   
           t2   -7.132e-05 -3.892e-02 -3.453e-02 -1.453e-02    ...     -1.994e-01   
           t3    1.003e-04 -9.286e-03 -7.856e-03 -3.736e-02    ...     -4.493e-02   
           r1    2.507e-06  1.321e-03  1.174e-03  5.725e-04    ...     -3.039e-02   
           r2    3.812e-06 -5.683e-05  5.614e-05 -1.004e-02    ...     -1.174e-02   
           r3   -1.439e-07 -1.571e-04  1.600e-04  1.240e-03    ...      6.409e-03   
    4      t1    9.580e-05 -1.466e-03  1.704e-03 -2.679e-01    ...      2.726e-02   
           t2   -7.179e-05 -3.925e-02 -3.419e-02 -1.184e-02    ...     -1.877e-01   
           t3    9.276e-05 -9.173e-03 -7.968e-03 -1.742e-02    ...     -2.215e-02   
           r1    2.531e-06  1.322e-03  1.169e-03  5.129e-04    ...     -3.086e-02   
           r2    3.812e-06 -5.682e-05  5.614e-05 -1.004e-02    ...     -1.174e-02   
           r3   -1.684e-07 -1.412e-04  1.769e-04  1.298e-03    ...      1.457e-02   
    5      t1    1.054e-04 -1.626e-03  1.863e-03 -2.930e-01    ...     -1.463e-03   
           t2   -7.785e-05 -4.223e-02 -3.750e-02 -1.564e-02    ...     -1.560e-01   
           t3    1.004e-04 -9.286e-03 -7.856e-03 -3.743e-02    ...     -4.535e-02   
           r1    2.518e-06  1.322e-03  1.172e-03  5.440e-04    ...     -3.061e-02   
           r2    3.812e-06 -5.683e-05  5.614e-05 -1.004e-02    ...     -1.174e-02   
           r3   -1.092e-07 -1.774e-04  1.806e-04  1.011e-03    ...      4.109e-05   
    ...                ...        ...        ...        ...    ...            ...   
    5629   t1   -2.694e-06  2.441e-05  1.075e-03  2.742e-03    ...      4.748e-04   
           t2   -3.617e-06 -1.361e-04  1.100e-04  8.914e-04    ...      2.047e-01   
           t3    5.577e-06 -5.649e-03 -5.097e-03  7.869e-03    ...     -1.857e-01   
           r1    5.464e-07 -2.376e-04 -2.019e-04 -6.031e-05    ...     -4.279e-04   
           r2    1.826e-07  7.492e-05  1.152e-04 -9.475e-04    ...      2.369e-02   
           r3    1.857e-07 -1.044e-04 -6.757e-05 -1.060e-04    ...      2.703e-02   
    5630   t1    0.000e+00  0.000e+00  0.000e+00  0.000e+00    ...      0.000e+00   
           t2    0.000e+00  0.000e+00  0.000e+00  0.000e+00    ...      0.000e+00   
           t3    0.000e+00  0.000e+00  0.000e+00  0.000e+00    ...      0.000e+00   
           r1    1.329e-06  7.127e-04  4.621e-04 -6.382e-04    ...     -3.555e-02   
           r2   -3.178e-07 -1.708e-05 -1.350e-05  3.852e-03    ...      6.103e-03   
           r3   -6.212e-07 -2.275e-04 -1.247e-04  5.015e-04    ...      1.508e-02   
    5631   t1    0.000e+00  0.000e+00  0.000e+00  0.000e+00    ...      0.000e+00   
           t2    0.000e+00  0.000e+00  0.000e+00  0.000e+00    ...      0.000e+00   
           t3    0.000e+00  0.000e+00  0.000e+00  0.000e+00    ...      0.000e+00   
           r1   -1.679e-06 -2.394e-04 -2.043e-04 -3.901e-05    ...      1.138e-03   
           r2    5.311e-07  6.254e-05 -5.671e-05 -2.159e-03    ...     -2.994e-02   
           r3   -3.299e-07 -1.454e-05 -9.195e-06  8.113e-04    ...      7.605e-03   
    5632   t1    0.000e+00  0.000e+00  0.000e+00  0.000e+00    ...      0.000e+00   
           t2    0.000e+00  0.000e+00  0.000e+00  0.000e+00    ...      0.000e+00   
           t3    0.000e+00  0.000e+00  0.000e+00  0.000e+00    ...      0.000e+00   
           r1    1.054e-06  5.821e-04  6.400e-04  9.411e-04    ...      3.064e-02   
           r2   -2.756e-06  1.174e-05  1.113e-05  3.844e-03    ...      1.025e-02   
           r3    3.049e-07  1.836e-04  2.254e-04  5.881e-04    ...      2.334e-02   
    5633   t1    0.000e+00  0.000e+00  0.000e+00  0.000e+00    ...      0.000e+00   
           t2    0.000e+00  0.000e+00  0.000e+00  0.000e+00    ...      0.000e+00   
           t3    0.000e+00  0.000e+00  0.000e+00  0.000e+00    ...      0.000e+00   
           r1    5.464e-07 -2.376e-04 -2.019e-04 -6.030e-05    ...     -4.279e-04   
           r2    3.548e-08 -4.728e-05  4.650e-05 -2.113e-04    ...      2.084e-02   
           r3    2.948e-07 -1.384e-05 -1.663e-05 -6.516e-04    ...      2.914e-02   
    
    Mode               159        160        161        162     163     164  \
    Freq           385.301    387.260    390.518    390.990 391.050 393.165   
    NodeID Item                                                               
    1      t1   -5.369e-02 -3.838e-02 -1.326e-01 -1.973e-02   0.028   0.033   
           t2   -3.746e-01 -5.840e-02 -2.385e-02 -5.889e-02   0.015   0.177   
           t3    1.271e-01  2.550e-01 -1.792e-01 -1.136e-03   0.042  -0.037   
           r1   -9.829e-04  2.993e-02 -3.527e-02  1.148e-04   0.007  -0.053   
           r2    1.241e-03  1.025e-02 -3.112e-02 -4.135e-03   0.011   0.026   
           r3    2.184e-02  2.495e-03  8.832e-02  1.660e-02  -0.030  -0.100   
    2      t1   -5.369e-02 -3.838e-02 -1.326e-01 -1.973e-02   0.028   0.033   
           t2   -3.309e-01 -5.341e-02  1.528e-01 -2.568e-02  -0.045  -0.022   
           t3    1.246e-01  2.345e-01 -1.170e-01  7.135e-03   0.020  -0.090   
           r1   -9.829e-04  2.993e-02 -3.527e-02  1.148e-04   0.007  -0.053   
           r2    1.241e-03  1.025e-02 -3.112e-02 -4.135e-03   0.011   0.026   
           r3    2.184e-02  2.495e-03  8.832e-02  1.660e-02  -0.030  -0.100   
    3      t1   -6.243e-02 -6.576e-03 -2.369e-01 -3.571e-02   0.065   0.135   
           t2   -3.102e-01 -1.168e-01  1.054e-01 -4.058e-02  -0.023   0.200   
           t3    1.259e-01  2.542e-01 -1.777e-01 -1.024e-03   0.042  -0.038   
           r1    2.253e-04  2.894e-02 -3.716e-02 -4.793e-04   0.008  -0.047   
           r2    1.238e-03  1.025e-02 -3.111e-02 -4.135e-03   0.011   0.026   
           r3   -2.870e-02  6.276e-03 -3.529e-02 -1.277e-02   0.016   0.091   
    4      t1   -6.325e-02 -6.637e-03 -2.366e-01 -3.568e-02   0.065   0.135   
           t2   -3.177e-01 -1.326e-01  1.642e-01 -3.435e-02  -0.035   0.187   
           t3    1.256e-01  2.349e-01 -1.179e-01  7.068e-03   0.020  -0.089   
           r1   -3.467e-03  3.029e-02 -3.335e-02  6.161e-04   0.007  -0.058   
           r2    1.238e-03  1.026e-02 -3.112e-02 -4.135e-03   0.011   0.026   
           r3   -1.664e-02  4.624e-03 -3.608e-02 -1.077e-02   0.016   0.083   
    5      t1   -4.748e-02  1.290e-02 -2.882e-01 -4.040e-02   0.084   0.165   
           t2   -3.697e-01 -2.080e-01  1.525e-01 -5.946e-02  -0.022   0.442   
           t3    1.271e-01  2.550e-01 -1.792e-01 -1.136e-03   0.042  -0.037   
           r1   -9.829e-04  2.993e-02 -3.527e-02  1.148e-04   0.007  -0.053   
           r2    1.241e-03  1.025e-02 -3.112e-02 -4.135e-03   0.011   0.026   
           r3    2.184e-02  2.495e-03  8.832e-02  1.660e-02  -0.030  -0.100   
    ...                ...        ...        ...        ...     ...     ...   
    5629   t1    5.261e-03 -4.662e-02  5.886e-02 -6.227e-03  -0.017   0.162   
           t2    6.117e-02 -5.444e-02 -1.529e-02 -1.469e-02   0.023   0.183   
           t3   -2.785e-02  6.353e-02  5.410e-02  2.473e-02  -0.030  -0.234   
           r1   -3.524e-03  9.710e-04 -6.896e-03  9.867e-04   0.001  -0.014   
           r2    1.095e-03 -8.118e-03 -1.289e-02 -1.785e-03   0.005   0.015   
           r3    1.362e-03 -5.113e-03 -1.492e-02 -6.335e-04   0.005   0.001   
    5630   t1    0.000e+00  0.000e+00  0.000e+00  0.000e+00   0.000   0.000   
           t2    0.000e+00  0.000e+00  0.000e+00  0.000e+00   0.000   0.000   
           t3    0.000e+00  0.000e+00  0.000e+00  0.000e+00   0.000   0.000   
           r1   -8.501e-03  1.420e-02 -1.374e-02  1.390e-04   0.004  -0.017   
           r2   -6.432e-03  1.082e-02 -3.451e-02 -4.203e-03   0.011   0.018   
           r3   -1.308e-03  3.008e-03 -9.727e-03  3.836e-04   0.001  -0.003   
    5631   t1    0.000e+00  0.000e+00  0.000e+00  0.000e+00   0.000   0.000   
           t2    0.000e+00  0.000e+00  0.000e+00  0.000e+00   0.000   0.000   
           t3    0.000e+00  0.000e+00  0.000e+00  0.000e+00   0.000   0.000   
           r1   -1.261e-02  1.119e-02 -1.439e-02  1.245e-03   0.004  -0.024   
           r2   -4.564e-03  1.167e-02 -1.208e-02 -2.319e-03   0.004   0.012   
           r3   -3.327e-03  1.359e-02  8.849e-04 -7.085e-04   0.002   0.011   
    5632   t1    0.000e+00  0.000e+00  0.000e+00  0.000e+00   0.000   0.000   
           t2    0.000e+00  0.000e+00  0.000e+00  0.000e+00   0.000   0.000   
           t3    0.000e+00  0.000e+00  0.000e+00  0.000e+00   0.000   0.000   
           r1   -2.242e-03  5.439e-04 -1.809e-02  1.961e-04   0.005  -0.016   
           r2    4.245e-03  2.363e-03 -2.781e-02 -4.249e-03   0.009   0.026   
           r3   -1.135e-03  5.283e-03 -1.865e-03 -3.915e-03   0.002   0.024   
    5633   t1    0.000e+00  0.000e+00  0.000e+00  0.000e+00   0.000   0.000   
           t2    0.000e+00  0.000e+00  0.000e+00  0.000e+00   0.000   0.000   
           t3    0.000e+00  0.000e+00  0.000e+00  0.000e+00   0.000   0.000   
           r1   -3.524e-03  9.710e-04 -6.896e-03  9.867e-04   0.001  -0.014   
           r2    1.171e-03 -6.235e-03 -1.349e-02 -1.096e-03   0.005   0.006   
           r3    1.305e-03 -6.509e-03 -1.448e-02 -1.144e-03   0.005   0.008   
    
    Mode            165        166        167  
    Freq        395.101    395.329    397.237  
    NodeID Item                                
    1      t1    -0.104 -6.919e-02 -1.904e-02  
           t2    -0.010  5.252e-02  1.187e-01  
           t3     0.263  2.141e-01 -1.473e-01  
           r1    -0.004  2.357e-02 -3.403e-02  
           r2     0.009  7.311e-03 -9.083e-04  
           r3     0.022  2.547e-02 -5.581e-03  
    2      t1    -0.104 -6.919e-02 -1.904e-02  
           t2     0.034  1.035e-01  1.075e-01  
           t3     0.244  1.995e-01 -1.454e-01  
           r1    -0.004  2.357e-02 -3.403e-02  
           r2     0.009  7.311e-03 -9.083e-04  
           r3     0.022  2.547e-02 -5.581e-03  
    3      t1    -0.079 -5.365e-02 -2.056e-02  
           t2     0.023  1.385e-02  1.724e-01  
           t3     0.262  2.140e-01 -1.467e-01  
           r1    -0.006  2.042e-02 -3.308e-02  
           r2     0.009  7.308e-03 -9.065e-04  
           r3    -0.024 -3.187e-02  1.810e-03  
    4      t1    -0.079 -5.359e-02 -2.070e-02  
           t2     0.006  6.929e-04  1.884e-01  
           t3     0.243  1.992e-01 -1.457e-01  
           r1    -0.003  2.650e-02 -3.497e-02  
           r2     0.010  7.313e-03 -9.084e-04  
           r3    -0.026 -3.403e-02  3.453e-04  
    5      t1    -0.056 -3.263e-02 -2.358e-02  
           t2     0.012 -6.531e-02  2.889e-01  
           t3     0.263  2.141e-01 -1.473e-01  
           r1    -0.004  2.357e-02 -3.403e-02  
           r2     0.009  7.311e-03 -9.083e-04  
           r3     0.022  2.547e-02 -5.581e-03  
    ...             ...        ...        ...  
    5629   t1    -0.557  6.614e-01 -1.042e-01  
           t2     0.290 -3.938e-01  3.587e-01  
           t3     0.069  6.091e-02 -3.214e-01  
           r1    -0.008  2.789e-02 -2.645e-02  
           r2    -0.021  3.344e-02  3.849e-03  
           r3     0.027 -1.418e-02  1.118e-02  
    5630   t1     0.000  0.000e+00  0.000e+00  
           t2     0.000  0.000e+00  0.000e+00  
           t3     0.000  0.000e+00  0.000e+00  
           r1     0.004 -2.480e-04 -2.458e-03  
           r2    -0.003  7.885e-03  2.321e-02  
           r3    -0.030 -3.103e-02 -2.712e-02  
    5631   t1     0.000  0.000e+00  0.000e+00  
           t2     0.000  0.000e+00  0.000e+00  
           t3     0.000  0.000e+00  0.000e+00  
           r1    -0.012  4.048e-03  1.465e-02  
           r2    -0.005 -5.452e-03 -5.399e-03  
           r3    -0.018 -1.781e-02 -1.326e-02  
    5632   t1     0.000  0.000e+00  0.000e+00  
           t2     0.000  0.000e+00  0.000e+00  
           t3     0.000  0.000e+00  0.000e+00  
           r1     0.004  1.558e-04 -4.253e-03  
           r2     0.024  1.089e-02 -1.333e-02  
           r3    -0.032 -2.891e-02 -1.447e-02  
    5633   t1     0.000  0.000e+00  0.000e+00  
           t2     0.000  0.000e+00  0.000e+00  
           t3     0.000  0.000e+00  0.000e+00  
           r1    -0.008  2.789e-02 -2.645e-02  
           r2     0.005  4.639e-03  6.872e-03  
           r3     0.008  7.160e-03  8.942e-03  
    
    [32274 rows x 167 columns]
    

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
    
    
    # extra list-type parameter for modal analysis; see dataNames
    #print("mode_cycles =", plate_stress.mode_cycles)


.. parsed-literal::

    plate_stress_obj = <class 'pyNastran.op2.tables.oes_stressStrain.real.oes_plates.RealPlateStressArray'>
    plate_stress = ['_add_new_eid', '_add_new_node', 'subtitle', 'words', 's_code', 'is_built', 'stress_bits', 'load_set', 'itotal', '_add', '_ntotals', 'nnodes', 'isubcase', 'element_name', 'itime', 'nonlinear_factor', 'title', '_times', 'ntotal', 'approach_code', 'sort_bits', 'label', 'element_node', 'is_msc', 'num_wide', 'mode', 'format_code', 'device_code', 'modes', '_times_dtype', 'thermal_bits', 'mode2s', 'data_frame', 'mode2', 'dt', 'data', 'cycle', 'name', 'nelements', 'eigr', 'ielement', 'thermal', 'analysis_code', 'eigrs', 'table_code', 'element_type', 'table_name', 'data_code', 'isTransient', 'sort_code', 'cycles', 'ntimes', 'data_names']
    
    data_code_keys = [u'subtitle', u'stress_bits', u'load_set', u'thermal', u's_code', u'isubcase', u'element_name', u'mode2', u'title', u'approach_code', u'sort_bits', u'label', u'is_msc', u'num_wide', u'format_code', u'device_code', u'_times_dtype', u'thermal_bits', u'nonlinear_factor', u'cycle', u'name', u'eigr', u'analysis_code', u'table_code', u'element_type', u'table_name', u'mode', u'sort_code', u'data_names']
    
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
    model = read_bdf(bdf_filename)
    mass, cg, I = model.mass_properties()


.. parsed-literal::

    DEBUG:     fname=bdf.pyc                   lineNo=723    ---starting BDF.read_bdf of f:\work\pynastran\pynastran\master3\models\iSat\ISat_Launch_Sm_4pt.dat---
    DEBUG:     fname=bdf.pyc                   lineNo=2658   opening 'f:\\work\\pynastran\\pynastran\\master3\\models\\iSat\\ISat_Launch_Sm_4pt.dat'
    DEBUG:     fname=bdf.pyc                   lineNo=2658   opening u'f:\\work\\pynastran\\pynastran\\master3\\models\\iSat\\write_modes.v2005'
    DEBUG:     fname=cross_reference.pyc       lineNo=527    Cross Referencing...
    WARNING:   fname=shell.pyc                 lineNo=1387   PSHELL pid=1 midsurface: z1=0.400000006 z2=-0.400000006 t=0.035999998 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.pyc                 lineNo=1387   PSHELL pid=2 midsurface: z1=0.400000006 z2=-0.400000006 t=0.054000005 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.pyc                 lineNo=1387   PSHELL pid=3 midsurface: z1=0.400000006 z2=-0.400000006 t=0.017999999 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.pyc                 lineNo=1387   PSHELL pid=7 midsurface: z1=0.418000013 z2=-0.418000013 t=0.035999998 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.pyc                 lineNo=1387   PSHELL pid=34 midsurface: z1=0.194000006 z2=-0.194000006 t=0.0186 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.pyc                 lineNo=1387   PSHELL pid=37 midsurface: z1=0.308999985 z2=-0.308999985 t=0.0186 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.pyc                 lineNo=1387   PSHELL pid=38 midsurface: z1=0.284000009 z2=-0.284000009 t=0.0186 not in range of -1.5t < zi < 1.5t
    WARNING:   fname=shell.pyc                 lineNo=1387   PSHELL pid=46 midsurface: z1=0.199000001 z2=-0.199000001 t=0.0186 not in range of -1.5t < zi < 1.5t
    DEBUG:     fname=bdf.pyc                   lineNo=749    ---finished BDF.read_bdf of f:\work\pynastran\pynastran\master3\models\iSat\ISat_Launch_Sm_4pt.dat---
    

.. code:: python

    import getpass
    name = getpass.getuser()
    os.chdir(os.path.join(r'C:\Users', name, 'Desktop'))
    
    # write the F06 with Real/Imaginary or Magnitude/Phase
    # only matters for complex results
    #op2.write_f06('isat.f06', is_mag_phase=False)
    
    !head -n 40 isat.f06
    gpw = op2.grid_point_weight
    print(gpw.object_attributes())
    
    # The mass results are different as pyNastran's mass assumes point masses.
    # The larger your model is and the further from the origin, the more accurate
    # the result.
    # For some applications (e.g. a weight breakdown), this may be fine.
    print('cg =\n%s' % gpw.cg)
    print('cg = %s' % cg)


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
    
    1    ISAT_SM_LAUNCH_4PT MODES TO 400 HZ                                     JANUARY   4, 2016  pyNastran v0.8.0+dev.f897603  PAGE   506
    
    1    ISAT_SM_LAUNCH_4PT MODES TO 400 HZ                                     JANUARY   4, 2016  pyNastran v0.8.0+dev.f897603  PAGE   507
         DEFAULT                                                                                                                        
    
                                                  R E A L   E I G E N V A L U E S
                                                 ISAT_SM_LAUNCH_4PT MODES TO 400 HZ
       MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED
        NO.       ORDER                                                                       MASS              STIFFNESS
    ['IQ', 'IS', 'MO', 'Q', 'S', 'cg', 'mass', 'reference_point']
    cg =
    [[ -1.026e-17  -2.531e+00  -1.847e+01]
     [ -3.385e-02  -1.633e-17  -1.847e+01]
     [ -3.385e-02  -2.531e+00  -2.181e-19]]
    cg = [ -0.035  -2.623 -18.53 ]
    
