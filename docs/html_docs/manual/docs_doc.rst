====
docs
====

------------
Introduction
------------
This is meant as a tutorial on how to use the pyNastran's reference documentation

The easiest way to access the documentation is to click:

Note that a static model is a SOL 101 or SOL 144.  A dynamic/"transient" solution is any transient/modal/load step/frequency based solution (e.g. 103, 109, 145).


The **head**/**tail**/**file_slice** methods can be found at:

    https://github.com/SteveDoyle2/pyNastran/blob/v0.7/docs_sphinx/manual/py_docs/bdf_doc.py

These examples can be found at:

    https://github.com/SteveDoyle2/pyNastran/blob/v0.7/docs_sphinx/manual/py_docs/op2_doc.py

---------------------
Example 1: Read Write
---------------------
This example will demonstate:

 - reading the OP2

 - getting some basic information

 - writing the F06

our model

.. code-block:: python

    >>> import pyNastran
    >>> pkg_path = pyNastran.__path__[0]
    >>> test_path = os.path.join(pkg_path, '..', 'models', 'solid_bending')
    >>> op2_filename = os.path.join(test_path, 'solid_bending.op2')
    >>> f06_filename = os.path.join(test_path, 'solid_bending_out.f06')

instantiate the model

.. code-block:: python

    >>> from pyNastran.op2.op2 import OP2
    >>> model = OP2()
    >>> model.read_op2(op2_filename)
    >>> print(model.get_op2_stats())
    op2.displacements[1]
      type=RealDisplacementArray nnodes=72
      data: [t1, t2, t3, r1, r2, r3] shape=[1, 72, 6] dtype=float32
      gridTypes
      lsdvmns = [1]
    
    op2.spc_forces[1]
      type=RealSPCForcesArray nnodes=72
      data: [t1, t2, t3, r1, r2, r3] shape=[1, 72, 6] dtype=float32
      gridTypes
      lsdvmns = [1]
    
    op2.ctetra_stress[1]
      type=RealSolidStressArray nelements=186 nnodes=930
      nodes_per_element=5 (including centroid)
      eType, cid
      data: [1, nnodes, 10] where 10=[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, von_mises]
      data.shape = (1, 930, 10)
      element types: CTETRA
      lsdvmns = [1]
    >>> model.write_f06(f06_filename)
    F06:
     RealDisplacementArray SUBCASE=1
     RealSPCForcesArray    SUBCASE=1
     RealSolidStressArray  SUBCASE=1 - CTETRA
    >>> print(tail(f06_filename, 21))
    0       186           0GRID CS  4 GP
    0                CENTER  X   9.658666E+02  XY  -2.978357E+01   A   2.559537E+04  LX-0.02 0.20 0.98  -1.094517E+04    2.288671E+04
                             Y   7.329372E+03  YZ   5.895411E+02   B  -7.168877E+01  LY-1.00-0.03-0.01
                             Z   2.454026E+04  ZX  -5.050599E+03   C   7.311813E+03  LZ 0.03-0.98 0.20
    0                     8  X   9.658666E+02  XY  -2.978357E+01   A   2.559537E+04  LX-0.02 0.20 0.98  -1.094517E+04    2.288671E+04
                             Y   7.329372E+03  YZ   5.895411E+02   B  -7.168877E+01  LY-1.00-0.03-0.01
                             Z   2.454026E+04  ZX  -5.050599E+03   C   7.311813E+03  LZ 0.03-0.98 0.20
    0                    62  X   9.658666E+02  XY  -2.978357E+01   A   2.559537E+04  LX-0.02 0.20 0.98  -1.094517E+04    2.288671E+04
                             Y   7.329372E+03  YZ   5.895411E+02   B  -7.168877E+01  LY-1.00-0.03-0.01
                             Z   2.454026E+04  ZX  -5.050599E+03   C   7.311813E+03  LZ 0.03-0.98 0.20
    0                     4  X   9.658666E+02  XY  -2.978357E+01   A   2.559537E+04  LX-0.02 0.20 0.98  -1.094517E+04    2.288671E+04
                             Y   7.329372E+03  YZ   5.895411E+02   B  -7.168877E+01  LY-1.00-0.03-0.01
                             Z   2.454026E+04  ZX  -5.050599E+03   C   7.311813E+03  LZ 0.03-0.98 0.20
    0                    58  X   9.658666E+02  XY  -2.978357E+01   A   2.559537E+04  LX-0.02 0.20 0.98  -1.094517E+04    2.288671E+04
                             Y   7.329372E+03  YZ   5.895411E+02   B  -7.168877E+01  LY-1.00-0.03-0.01
                             Z   2.454026E+04  ZX  -5.050599E+03   C   7.311813E+03  LZ 0.03-0.98 0.20
    1    MSC.NASTRAN JOB CREATED ON 28-JAN-12 AT 12:52:32                       JANUARY  28, 2012  pyNastran v0.7.1       PAGE     3
    
    1                                        * * * END OF JOB * * *

--------------------------------
Example 2: Displacement (static)
--------------------------------
This example will demonstate:

 - calculating total deflection of the nodes for a static case for an OP2
 - calculate von mises stress and max shear


.. math:: \sqrt\left(T_x^2 + T_y^2 + T_z^2\right)

our model

.. code-block:: python

    >>> import pyNastran
    >>> pkg_path = pyNastran.__path__[0]
    >>> test_path = os.path.join(pkg_path, '..', 'models', 'solid_bending')
    >>> op2_filename = os.path.join(test_path, 'solid_bending.op2')
    >>> out_filename = os.path.join(test_path, 'solid_bending.out')

instantiate the model

.. code-block:: python

    >>> from pyNastran.op2.op2 import OP2
    >>> model = OP2()
    >>> model.read_op2(op2_filename)
    >>> print(model.get_op2_stats())

we're analyzing a static problem, so itime=0

we're also assuming subcase 1

.. code-block:: python

    >>> itime = 0
    >>> isubcase = 1

get the displacement object

.. code-block:: python

    >>> disp = model.displacements[isubcase]

displacement is an array

.. code-block:: python

    # data = [tx, ty, tz, rx, ry, rz]
    # for some itime
    # all the nodes -> :
    # get [tx, ty, tz] -> :3
    >>> txyz = disp.data[itime, :, :3]

calculate the total deflection of the vector

.. code-block:: python

    >>> from numpy.linalg import norm
    >>> total_xyz = norm(txyz, axis=1)

since norm's axis parameter can be tricky, we'll double check the length

.. code-block:: python

    >>> nnodes = disp.data.shape[1]
    >>> assert len(total_xyz) == nnodes

we could also have found nnodes by using the attribute.

It has an underscore because the object is also used for elements.

.. code-block:: python

    >>> nnodes2 = disp._nnodes
    >>> assert nnodes == nnodes2
    >>> assert nnodes == 72

additionally we know we have 72 nodes from the shape:

.. code-block:: python

    op2.displacements[1]
      type=RealDisplacementArray nnodes=72
      data: [t1, t2, t3, r1, r2, r3] shape=[1, 72, 6] dtype=float32
      gridTypes
      lsdvmns = [1]

now we'll loop over the nodes and print the total deflection

.. code-block:: python

    >>> msg = 'nid, gridtype, tx, ty, tz, txyz'
    >>> print(msg)
    >>> for (nid, grid_type), txyz, total_xyzi in zip(disp.node_gridtype, txyz, total_xyz):
    >>>     msg = '%s, %s, %s, %s, %s, %s' % (nid, grid_type, txyz[0], txyz[1], txyz[2], total_xyzi)
    >>>     print(msg)

    nid, gridtype, tx, ty, tz, txyz
    1, 1, 0.00764469, 4.01389e-05, 0.000111137, 0.00764561
    2, 1, 0.00762899, 5.29171e-05, 0.000142154, 0.0076305
    3, 1, 0.00944763, 6.38675e-05, 7.66179e-05, 0.00944816
    4, 1, 0.00427092, 2.62277e-05, 7.27848e-05, 0.00427162
    5, 1, 0.00152884, 1.71054e-05, -3.47525e-06, 0.00152894
    ...

----------------------------------
Example 3: Eigenvector (transient)
----------------------------------
This example will demonstate:

 - calculate von mises stress and max shear for solid elements for a static case for an OP2


.. math:: \sqrt\left(T_x^2 + T_y^2 + T_z^2\right)

our model

.. code-block:: python

    >>> import pyNastran
    >>> pkg_path = pyNastran.__path__[0]
    >>> test_path = os.path.join(pkg_path, '..', 'models', 'solid_bending')
    >>> op2_filename = os.path.join(test_path, 'solid_bending.op2')
    >>> out_filename = os.path.join(test_path, 'solid_bending.out')

instantiate the model

.. code-block:: python

    >>> from pyNastran.op2.op2 import OP2
    >>> model = OP2()
    >>> model.read_op2(op2_filename)
    >>> print(model.get_op2_stats())

    op2.ctetra_stress[1]
      type=RealSolidStressArray nelements=186 nnodes=930
      nodes_per_element=5 (including centroid)
      eType, cid
      data: [1, nnodes, 10] where 10=[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, von_mises]
      data.shape = (1, 930, 10)
      element types: CTETRA
      lsdvmns = [1]

we're analyzing a static problem, so itime=0

we're also assuming subcase 1

.. code-block:: python

    >>> itime = 0
    >>> isubcase = 1

get the stress object (there is also cpenta_stress and chexa_stress as well as ctetra_strain/cpenta_strain/chexa_strain)

.. code-block:: python

    >>> stress = model.ctetra_stress[isubcase]

The stress/strain data can often be von_mises/max_shear (same for fiber_distance/curvature), so check!

.. code-block:: python

     #data = [oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, von_mises]
    >>> o1 = stress.data[itime, :, 6]
    >>> o3 = stress.data[itime, :, 8]
    >>> if stress.is_von_mises():
    >>>     max_shear = (o1 - o3) / 2.
    >>>     von_mises = stress.data[itime, :, 9]
    >>> else:
    >>>     from numpy import sqrt
    >>>     o2 = data[itime, :, 8]
    >>>     von_mises = sqrt(0.5*((o1-o2)**2 + (o2-o3)**2, (o3-o1)**2))
    >>>     max_shear = stress.data[itime, :, 9]
    >>> for (eid, node), vm, ms in zip(stress.element_node, von_mises, max_shear):
    >>>     print(eid, 'CEN/4' if node == 0 else node, vm, ms)

    1 CEN/4 15900.2 2957.35
    1 8     15900.2 2957.35
    1 13    15900.2 2957.35
    1 67    15900.2 2957.35
    1 33    15900.2 2957.35
    2 CEN/4 16272.3 6326.18
    2 8     16272.3 6326.18
    2 7     16272.3 6326.18
    2 62    16272.3 6326.18
    2 59    16272.3 6326.18

Note that because element_node is an integer array, the centroid is 0.  We renamed it to CEN/4 when we wrote it

--------------------------------
Example 4: Solid Stress (static)
--------------------------------
This example will demonstate:

 - calculating total deflection of the nodes for a dynamic case for an OP2


.. math:: \sqrt\left(T_x^2 + T_y^2 + T_z^2\right)

our model

.. code-block:: python

    >>> import pyNastran
    >>> pkg_path = pyNastran.__path__[0]
    >>> test_path = os.path.join(pkg_path, '..', 'models', 'plate_py')
    >>> op2_filename = os.path.join(test_path, 'plate_py.op2')

ut_filename = os.path.join(test_path, 'solid_bending.out')

instantiate the model

.. code-block:: python

    >>> from pyNastran.op2.op2 import OP2
    >>> model = OP2()
    >>> model.read_op2(op2_filename)
    >>> print(model.get_op2_stats())

    op2.eigenvectors[1]
      type=RealEigenvectorArray ntimes=10 nnodes=231
      data: [t1, t2, t3, r1, r2, r3] shape=[10, 231, 6] dtype=float32
      gridTypes
      modes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    eigrs = [-0.00037413835525512695, -0.00022113323211669922, -0.0001882314682006836, -0.00010025501251220703, 0.0001621246337890625, 0.00
    07478296756744385, 1583362560.0, 2217974016.0, 10409966592.0, 11627085824.0]
    mode_cycles = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    >>> isubcase = 1
    >>> eigenvector = model.eigenvectors[isubcase]

"time/mode/frequency are stored by id, so to get mode 5:

.. code-block:: python

    >>> modes = eigenvector._times  # it may not be "time" so we don't use the name "time"
    >>> from numpy import where
    >>> imode5 = where(modes == 5)[0]
    >>> txyz = eigenvector.data[imode5, :, :3]

calculate the total deflection of the vector

.. code-block:: python

    >>> from numpy.linalg import norm
    >>> total_xyz = norm(txyz, axis=1)

get the eigenvalue

.. code-block:: python

    >>> print('eigr5 = %s' % eigenvector.eigrs[imode5])
    eigr5 = 0.000162124633789

------------------------------------------
Example 5: Isotropic Plate Stress (static)
------------------------------------------
This example will demonstate:

 - print the fiber distance and the max principal stress for a static case for an OP2

our model

.. code-block:: python

    >>> import pyNastran
    >>> pkg_path = pyNastran.__path__[0]
    >>> test_path = os.path.join(pkg_path, '..', 'models', 'sol_101_elements')
    >>> op2_filename = os.path.join(test_path, 'static_solid_shell_bar.op2')

instantiate the model

.. code-block:: python

    >>> from pyNastran.op2.op2 import OP2
    >>> model = OP2()
    >>> model.read_op2(op2_filename)
    >>> print(model.get_op2_stats())

    op2.cquad4_stress[1]
      type=RealPlateStressArray nelements=2 nnodes_per_element=5 nlayers=2 ntotal=20
      data: [1, ntotal, 8] where 8=[fiber_distance, oxx, oyy, txy, angle, omax, omin, von_mises]
      data.shape=(1L, 20L, 8L)
      element types: CQUAD4
      lsdvmns = [1]
    >>> isubcase = 1
    >>> itime = 0 # this is a static case
    >>> stress = model.cquad4_stress[isubcase]
    >>> assert stress.nnodes == 5, 'this is a bilinear quad'

write the data

.. code-block:: python

    #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
    >>> eids = stress.element_node[:, 0]
    >>> nids = stress.element_node[:, 1]
    >>> if stress.is_fiber_distance():
    >>>     fiber_dist = stress.data[itime, :, 0]
    >>> else:
    >>>     raise RuntimeError('found fiber curvature; expected fiber distance')
    >>> maxp = stress.data[itime, :, 5]
    >>> for (eid, nid, fdi, maxpi) in zip(eids, nids, fiber_dist, maxp):
    >>>     print(eid, 'CEN/4' if nid == 0 else nid, fdi, maxpi)

    6 CEN/4 -0.125 8022.26
    6 CEN/4  0.125 12015.9
    6 4     -0.125 7580.84
    6 4      0.125 11872.9
    6 1     -0.125 8463.42
    6 1      0.125 12158.9
    6 14    -0.125 8463.69
    6 14     0.125 12158.9
    6 15    -0.125 7581.17
    6 15     0.125 11872.9
    7 CEN/4 -0.125 10016.3
    7 CEN/4  0.125 10019.5
    7 3     -0.125 10307.1
    7 3      0.125 10311.0
    7 2     -0.125 9725.54
    7 2      0.125 9727.9
    7 17    -0.125 9725.54
    7 17     0.125 9728.06
    7 16    -0.125 10307.1
    7 16     0.125 10311.1

note we have 2 layers (upper and lower surface) for any PSHELL-based elements

------------------------------------------
Example 6: Composite Plate Stress (static)
------------------------------------------
This example will demonstate:

 - print the fiber distance and the max principal stress for a static case for an OP2

our model

.. code-block:: python

    >>> import pyNastran
    >>> pkg_path = pyNastran.__path__[0]
    >>> test_path = os.path.join(pkg_path, '..', 'models', 'sol_101_elements')
    >>> op2_filename = os.path.join(test_path, 'static_solid_shell_bar.op2')

instantiate the model

.. code-block:: python

    >>> from pyNastran.op2.op2 import OP2
    >>> model = OP2()
    >>> model.read_op2(op2_filename)
    >>> print(model.get_op2_stats())
    op2.ctria3_composite_stress[1]
      type=RealCompositePlateStressArray nelements=4 ntotal=18
      data: [1, ntotal, 9] where 9=[o11, o22, t12, t1z, t2z, angle, major, minor, max_shear]
      data.shape = (1, 18, 9)
      element types: CTRIA3
      lsdvmns = [1]
    >>> isubcase = 1
    >>> itime = 0 # this is a static case
    >>> stress = model.ctria3_composite_stress[isubcase]

In the previous example, we had an option for a variable number of nodes for the CQUAD4s (1/5), but only nnodes=1 for the CTRIA3s.

In this example, we have 4 layers on one element and 5 on another, but they're all at the centroid.

.. code-block:: python

 #[o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
    >>> eids = stress.element_layer[:, 0]
    >>> layers = stress.element_layer[:, 1]
    >>> maxp = stress.data[itime, :, 6]
    >>> if stress.is_fiber_distance():
    >>>     fiber_dist = stress.data[itime, :, 0]
    >>> else:
    >>>     raise RuntimeError('found fiber curvature; expected fiber distance')
    >>> maxp = stress.data[itime, :, 5]
    >>> for (eid, layer, maxpi) in zip(eids, layers, maxp):
    >>>     print(eid, 'CEN/4', layer, maxpi)

    7  CEN/4 1  89.3406
    7  CEN/4 2  89.3745
    7  CEN/4 3  89.4313
    7  CEN/4 4  89.5115
    8  CEN/4 1 -85.6691
    8  CEN/4 2 -85.6121
    8  CEN/4 3 -85.5193
    8  CEN/4 4 -85.3937
    8  CEN/4 5 -85.2394
    9  CEN/4 1  86.3663
    9  CEN/4 2  86.6389
    9  CEN/4 3  87.0977
    9  CEN/4 4  87.7489
    10 CEN/4 1 -87.6962
    10 CEN/4 2 -87.4949
    10 CEN/4 3 -87.1543
    10 CEN/4 4 -86.6662
    10 CEN/4 5 -86.0192

