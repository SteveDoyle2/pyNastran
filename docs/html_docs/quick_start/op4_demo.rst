
OP4 Demo
--------

The Jupyter notebook for this demo can be found in: -
docs/quick_start/demo/op4_demo.ipynb -
https://github.com/SteveDoyle2/pyNastran/tree/master/docs/quick_start/demo/op4_demo.ipynb

The OP4 is a Nastran input/output format that can store matrices.

The OP2 can as well, but is less validated in regards to matrices.

Import pyNastran
^^^^^^^^^^^^^^^^

.. code:: ipython3

    import os
    
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    
    from pyNastran.utils import print_bad_path
    from pyNastran.op4.op4 import read_op4
    import numpy as np
    from numpy import float32, float64, int32, int64, product
    
    # decrease output precision
    np.set_printoptions(precision=3, threshold=20)

Print the docstring
^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    help(read_op4)


.. parsed-literal::

    Help on function read_op4 in module pyNastran.op4.op4:
    
    read_op4(op4_filename=None, matrix_names=None, precision='default', debug=False, log=None)
        Reads a NASTRAN OUTPUT4 file, and stores the
        matrices as the output arguments.  The number of
        matrices read is defined by the list matrix_names.  By default, all
        matrices will be read.  The resulting output is a dictionary of
        matrices that are accessed by their name.
        
        .. code-block:: python
        
          >>> from pyNastran.op4.op4 import OP4
          >>> op4 = OP4()
        
          # get all the matrices
          >>> matrices = op4.read_op4(op4_filename)
          >>> (formA, A) = matrices['A']
          >>> (formB, B) = matrices['B']
          >>> (formC, C) = matrices['C']
        
          # or to reduce memory usage
          >>> matrices = op4.read_op4(op4_filename, matrix_names=['A', 'B'])
          >>> (formA, A) = matrices['A']
          >>> (formB, B) = matrices['B']
        
          # or because you only want A
          >>> matrices = op4.read_op4(op4_filename, matrix_names='A')
          >>> (formA, A) = matrices['A']
        
          # get all the matrices, but select the file using a file dialog
          >>> matrices = op4.read_op4()
          >>>
        
        Parameters
        ----------
        op4_filename : str / None
            an OP4 filename.  Type=STRING.
        matrix_names : List[str], str / None
            matrix name(s) (None -> all)
        precision : str; {'default', 'single', 'double'}
            specifies if the matrices are in single or double precsion
            which means the format will be whatever the file is in
        
        Returns
        -------
        matricies : dict[str] = (int, matrix)
            dictionary of matrices where the key is the name and the value is [form, matrix]
        
            +------+----------------+
            | Form |   Definition   |
            +======+================+
            |  1   | Square         |
            +------+----------------+
            |  2   | Rectangular    |
            +------+----------------+
            |  3   | Diagonal       |
            +------+----------------+
            |  6   | Symmetric      |
            +------+----------------+
            |  8   | Id entity      |
            +------+----------------+
            |  9   | Pseudoidentity |
            +------+----------------+
        
            +--------+-------------------------+
            |  Type  | Object                  |
            +========+=========================+
            | Dense  | NUMPY.NDARRAY           |
            +--------+-------------------------+
            | Sparse | SCIPY.SPARSE.COO_MATRIX |
            +--------+-------------------------+
        
        .. note:: based off the MATLAB code SAVEOP4 developed by ATA-E and
                  later UCSD.
        .. note:: it's strongly recommended that you convert sparse matrices to
                  another format before doing math on them.  This is standard
                  with sparse matrices.
    
    

So as you can see, Nastran has many matrix formats.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # read the op4, will pop open a dialog box
    #matrices = read_op4()

.. code:: ipython3

    op4_filename = os.path.join(pkg_path, '..', 'models', 'iSat', 'ISat_Launch_Sm_4pt.op4')
    assert os.path.exists(op4_filename), print_bad_path(op4_filename)
    
    #specify the file
    matrices = read_op4(op4_filename)

There are more ways to read an OP4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # only 1 matrix
    matrices = read_op4(op4_filename, matrix_names='FLAMA', debug=False)
    
    # 1 or more matrices
    matrices = read_op4(op4_filename, matrix_names=['FLAMA','UGEXT'])

.. code:: ipython3

    # extract a matrix
    form, flama = matrices['FLAMA']
    print("form = %s" % form)
    print("type = %s" % type(flama))


.. parsed-literal::

    form = 2
    type = <class 'numpy.ndarray'>
    

.. code:: ipython3

    print("keys = %s" % matrices.keys())


.. parsed-literal::

    keys = dict_keys(['FLAMA', 'UGEXT'])
    

.. code:: ipython3

    print(matrices.keys())
    form_flama, flama = matrices['FLAMA']
    print("shape = %s" % str(flama.shape))
    print("flamat nvals = %s" % flama.size)
    
    form_ugext, ugext = matrices['UGEXT']
    print("form_ugext=%s type=%s" % (form_ugext, type(ugext[0,0])))
    #print "ugext", ugext
    print("ugext.shape = %s" % str(ugext.shape))
    print("ugext nvals = %s" % ugext.size)


.. parsed-literal::

    dict_keys(['FLAMA', 'UGEXT'])
    shape = (3, 167)
    flamat nvals = 501
    form_ugext=2 type=<class 'numpy.float64'>
    ugext.shape = (32274, 167)
    ugext nvals = 5389758
    

.. code:: ipython3

    print(ugext[:,:])
    #print(flama)


.. parsed-literal::

    [[-5.548e-03  4.671e-06 -1.816e-04 ... -1.037e-01  6.919e-02  1.904e-02]
     [ 2.133e-04  5.699e-03  2.393e-02 ... -1.050e-02 -5.252e-02 -1.187e-01]
     [-8.469e-04  1.512e-03  7.038e-03 ...  2.626e-01 -2.141e-01  1.472e-01]
     ...
     [ 3.006e-07  5.476e-05  6.343e-04 ... -8.222e-03 -2.789e-02  2.645e-02]
     [ 1.723e-06  1.278e-06 -1.805e-06 ...  4.866e-03 -4.639e-03 -6.872e-03]
     [-7.271e-06  3.394e-06 -2.717e-06 ...  7.772e-03 -7.160e-03 -8.942e-03]]
    
