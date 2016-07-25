
Static & Transient DataFrames in PyNastran
==========================================

The iPython notebook for this demo can be found in: -
docs:raw-latex:`\quick`\_start:raw-latex:`\demo`:raw-latex:`\op`2\_pandas\_multi\_case.ipynb
-
https://github.com/SteveDoyle2/pyNastran/tree/master/docs/quick\_start/demo/op2\_pandas\_multi\_case.ipynb

.. code:: python

    import os
    import pandas as pd
    import pyNastran
    from pyNastran.op2.op2 import read_op2
    
    pkg_path = pyNastran.__path__[0]
    model_path = os.path.join(pkg_path, '..', 'models')

Solid Bending
-------------

Let's show off ``combine=True/False``. We'll talk about the keys soon.

.. code:: python

    solid_bending_op2 = os.path.join(model_path, 'solid_bending', 'solid_bending.op2')
    solid_bending = read_op2(solid_bending_op2, combine=False, debug=False)
    print(solid_bending.displacements.keys())


.. parsed-literal::

    INFO:      fname=op2_scalar.py             lineNo=1176   op2_filename = 'f:\\work\\pynastran\\pynastran\\master3\\pyNastran\\..\\models\\solid_bending\\solid_bending.op2'
    [(1, 1, 1, 0, u'DEFAULT')]
    

.. code:: python

    solid_bending_op2 = os.path.join(model_path, 'solid_bending', 'solid_bending.op2')
    solid_bending2 = read_op2(solid_bending_op2, combine=True, debug=False)
    print(solid_bending2.displacements.keys())


.. parsed-literal::

    INFO:      fname=op2_scalar.py             lineNo=1176   op2_filename = 'f:\\work\\pynastran\\pynastran\\master3\\pyNastran\\..\\models\\solid_bending\\solid_bending.op2'
    [1]
    

Single Subcase Buckling Example
-------------------------------

The keys cannot be "combined" despite us telling the program that it was
OK. We'll get the following values that we need to handle. ####
isubcase, analysis\_code, sort\_method, count, subtitle \* isubcase ->
the same key that you're used to accessing \* sort\_method -> 1 (SORT1),
2 (SORT2) \* count -> the optimization count \* subtitle -> the analysis
subtitle (changes for superlements) \* analysis code -> the "type" of
solution

### Partial code for calculating analysis code:

::

       if trans_word == 'LOAD STEP':  # nonlinear statics
          analysis_code = 10
      elif trans_word in ['TIME', 'TIME STEP']:  # TODO check name
          analysis_code = 6
      elif trans_word == 'EIGENVALUE':  # normal modes
          analysis_code = 2
      elif trans_word == 'FREQ':  # TODO check name
          analysis_code = 5
      elif trans_word == 'FREQUENCY':
          analysis_code = 5
      elif trans_word == 'COMPLEX EIGENVALUE':
          analysis_code = 9
      else:
          raise NotImplementedError('transient_word=%r is not supported...' % trans_word)

Let's look at an odd case:
~~~~~~~~~~~~~~~~~~~~~~~~~~

You can do buckling as one subcase or two subcases (makes parsing it a
lot easier!).

However, you **have** to do this once you start messing around with
superelements or multi-step optimization.

For optimization, sometimes Nastran will downselect elements and do an
optimization on that and print out a subset of the elements. At the end,
it will rerun an analysis to double check the constraints are satisfied.
It does not always do multi-step optimization.

.. code:: python

    op2_filename = os.path.join(model_path, 'sol_101_elements', 'buckling_solid_shell_bar.op2')
    model = read_op2(op2_filename, combine=True, debug=False, build_dataframe=True)


.. parsed-literal::

    INFO:      fname=op2_scalar.py             lineNo=1176   op2_filename = 'f:\\work\\pynastran\\pynastran\\master3\\pyNastran\\..\\models\\sol_101_elements\\buckling_solid_shell_bar.op2'
    

.. code:: python

    stress_keys = model.cquad4_stress.keys()
    print (stress_keys)
    
    # isubcase, analysis_code, sort_method, count, subtitle
    key0 = (1, 1, 1, 0, 'DEFAULT1')
    key1 = (1, 8, 1, 0, 'DEFAULT1')


.. parsed-literal::

    [(1, 1, 1, 0, u'DEFAULT1'), (1, 8, 1, 0, u'DEFAULT1')]
    

Keys: \* key0 is the "static" key \* key1 is the "buckling" key

Similarly: \* Transient solutions can have preload \* Frequency
solutions can have loadsets (???)

Moving onto the data frames
---------------------------

-  The static case is the initial deflection state
-  The buckling case is "transient", where the modes (called load steps
   or lsdvmn here) represent the "times"

pyNastran reads these tables differently and handles them differently
internally. They look very similar though.

.. code:: python

    stress_static = model.cquad4_stress[key0].data_frame
    stress_transient = model.cquad4_stress[key1].data_frame
    
    # The final calculated factor:
    #   Is it a None or not?
    # This defines if it's static or transient
    print('stress_static.nonlinear_factor = %s' % model.cquad4_stress[key0].nonlinear_factor)
    print('stress_transient.nonlinear_factor = %s' % model.cquad4_stress[key1].nonlinear_factor)
    
    print('data_names  = %s' % model.cquad4_stress[key1].data_names)
    print('loadsteps   = %s' % model.cquad4_stress[key1].lsdvmns)
    print('eigenvalues = %s' % model.cquad4_stress[key1].eigrs)
    


.. parsed-literal::

    stress_static.nonlinear_factor = None
    stress_transient.nonlinear_factor = 4
    data_names  = [u'lsdvmn', u'eigr']
    loadsteps   = [1, 2, 3, 4]
    eigenvalues = [-49357660160.0, -58001940480.0, -379750744064.0, -428462538752.0]
    

Static Table
------------

.. code:: python

    # Sets default precision of real numbers for pandas output\n"
    pd.set_option('precision', 2)
    
    stress_static.head(20)




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th></th>
          <th></th>
          <th>index</th>
          <th>fiber_distance</th>
          <th>oxx</th>
          <th>oyy</th>
          <th>txy</th>
          <th>angle</th>
          <th>omax</th>
          <th>omin</th>
          <th>von_mises</th>
        </tr>
        <tr>
          <th>ElementID</th>
          <th>NodeID</th>
          <th>Location</th>
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
          <th rowspan="10" valign="top">6</th>
          <th rowspan="2" valign="top">CEN</th>
          <th>Top</th>
          <td>0</td>
          <td>-0.12</td>
          <td>5.85e-07</td>
          <td>9.73e-06</td>
          <td>-1.36e-07</td>
          <td>-89.15</td>
          <td>9.73e-06</td>
          <td>5.83e-07</td>
          <td>9.46e-06</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>1</td>
          <td>0.12</td>
          <td>4.71e-07</td>
          <td>9.44e-06</td>
          <td>-1.61e-07</td>
          <td>-88.97</td>
          <td>9.44e-06</td>
          <td>4.69e-07</td>
          <td>9.21e-06</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">4</th>
          <th>Top</th>
          <td>2</td>
          <td>-0.12</td>
          <td>-6.50e-07</td>
          <td>9.48e-06</td>
          <td>-1.36e-07</td>
          <td>-89.23</td>
          <td>9.48e-06</td>
          <td>-6.52e-07</td>
          <td>9.82e-06</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>3</td>
          <td>0.12</td>
          <td>-8.37e-07</td>
          <td>9.11e-06</td>
          <td>-1.61e-07</td>
          <td>-89.08</td>
          <td>9.12e-06</td>
          <td>-8.39e-07</td>
          <td>9.56e-06</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">1</th>
          <th>Top</th>
          <td>4</td>
          <td>-0.12</td>
          <td>-6.50e-07</td>
          <td>9.98e-06</td>
          <td>-1.36e-07</td>
          <td>-89.27</td>
          <td>9.99e-06</td>
          <td>-6.51e-07</td>
          <td>1.03e-05</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>5</td>
          <td>0.12</td>
          <td>-8.37e-07</td>
          <td>9.76e-06</td>
          <td>-1.61e-07</td>
          <td>-89.13</td>
          <td>9.76e-06</td>
          <td>-8.39e-07</td>
          <td>1.02e-05</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">14</th>
          <th>Top</th>
          <td>6</td>
          <td>-0.12</td>
          <td>1.82e-06</td>
          <td>9.98e-06</td>
          <td>-1.36e-07</td>
          <td>-89.05</td>
          <td>9.99e-06</td>
          <td>1.82e-06</td>
          <td>9.21e-06</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>7</td>
          <td>0.12</td>
          <td>1.78e-06</td>
          <td>9.76e-06</td>
          <td>-1.61e-07</td>
          <td>-88.85</td>
          <td>9.76e-06</td>
          <td>1.78e-06</td>
          <td>9.01e-06</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">15</th>
          <th>Top</th>
          <td>8</td>
          <td>-0.12</td>
          <td>1.82e-06</td>
          <td>9.48e-06</td>
          <td>-1.36e-07</td>
          <td>-88.98</td>
          <td>9.48e-06</td>
          <td>1.82e-06</td>
          <td>8.72e-06</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>9</td>
          <td>0.12</td>
          <td>1.78e-06</td>
          <td>9.11e-06</td>
          <td>-1.61e-07</td>
          <td>-88.75</td>
          <td>9.12e-06</td>
          <td>1.78e-06</td>
          <td>8.37e-06</td>
        </tr>
        <tr>
          <th rowspan="10" valign="top">7</th>
          <th rowspan="2" valign="top">CEN</th>
          <th>Top</th>
          <td>10</td>
          <td>-0.12</td>
          <td>7.16e-07</td>
          <td>1.02e-05</td>
          <td>1.22e-07</td>
          <td>89.26</td>
          <td>1.02e-05</td>
          <td>7.14e-07</td>
          <td>9.82e-06</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>11</td>
          <td>0.12</td>
          <td>7.31e-07</td>
          <td>1.04e-05</td>
          <td>1.53e-07</td>
          <td>89.10</td>
          <td>1.04e-05</td>
          <td>7.29e-07</td>
          <td>1.01e-05</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">3</th>
          <th>Top</th>
          <td>12</td>
          <td>-0.12</td>
          <td>-7.30e-07</td>
          <td>1.04e-05</td>
          <td>1.22e-07</td>
          <td>89.37</td>
          <td>1.04e-05</td>
          <td>-7.31e-07</td>
          <td>1.08e-05</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>13</td>
          <td>0.12</td>
          <td>-8.05e-07</td>
          <td>1.07e-05</td>
          <td>1.53e-07</td>
          <td>89.24</td>
          <td>1.07e-05</td>
          <td>-8.07e-07</td>
          <td>1.12e-05</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">2</th>
          <th>Top</th>
          <td>14</td>
          <td>-0.12</td>
          <td>-7.30e-07</td>
          <td>9.90e-06</td>
          <td>1.22e-07</td>
          <td>89.34</td>
          <td>9.90e-06</td>
          <td>-7.31e-07</td>
          <td>1.03e-05</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>15</td>
          <td>0.12</td>
          <td>-8.05e-07</td>
          <td>1.01e-05</td>
          <td>1.53e-07</td>
          <td>89.20</td>
          <td>1.01e-05</td>
          <td>-8.07e-07</td>
          <td>1.05e-05</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">17</th>
          <th>Top</th>
          <td>16</td>
          <td>-0.12</td>
          <td>2.16e-06</td>
          <td>9.90e-06</td>
          <td>1.22e-07</td>
          <td>89.10</td>
          <td>9.90e-06</td>
          <td>2.16e-06</td>
          <td>9.02e-06</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>17</td>
          <td>0.12</td>
          <td>2.27e-06</td>
          <td>1.01e-05</td>
          <td>1.53e-07</td>
          <td>88.88</td>
          <td>1.01e-05</td>
          <td>2.26e-06</td>
          <td>9.18e-06</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">16</th>
          <th>Top</th>
          <td>18</td>
          <td>-0.12</td>
          <td>2.16e-06</td>
          <td>1.04e-05</td>
          <td>1.22e-07</td>
          <td>89.15</td>
          <td>1.04e-05</td>
          <td>2.16e-06</td>
          <td>9.52e-06</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>19</td>
          <td>0.12</td>
          <td>2.27e-06</td>
          <td>1.07e-05</td>
          <td>1.53e-07</td>
          <td>88.96</td>
          <td>1.07e-05</td>
          <td>2.26e-06</td>
          <td>9.79e-06</td>
        </tr>
      </tbody>
    </table>
    </div>



Transient Table
---------------

.. code:: python

    # Sets default precision of real numbers for pandas output\n"
    pd.set_option('precision', 3)
    #import numpy as np
    #np.set_printoptions(formatter={'all':lambda x: '%g'})
    
    stress_transient.head(20)




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr>
          <th></th>
          <th></th>
          <th>LoadStep</th>
          <th>Item</th>
          <th>1</th>
          <th>2</th>
          <th>3</th>
          <th>4</th>
        </tr>
        <tr>
          <th></th>
          <th></th>
          <th>EigenvalueReal</th>
          <th></th>
          <th>-49357660160.0</th>
          <th>-58001940480.0</th>
          <th>-3.79750744064e+11</th>
          <th>-4.28462538752e+11</th>
        </tr>
        <tr>
          <th></th>
          <th></th>
          <th>Freq</th>
          <th></th>
          <th>35358.7915137</th>
          <th>38330.227181</th>
          <th>98077.5138317</th>
          <th>104178.13059</th>
        </tr>
        <tr>
          <th></th>
          <th></th>
          <th>Radians</th>
          <th></th>
          <th>222165.839318</th>
          <th>240835.920244</th>
          <th>616239.193872</th>
          <th>654570.499451</th>
        </tr>
        <tr>
          <th>ElementID</th>
          <th>NodeID</th>
          <th>Location</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th rowspan="20" valign="top">6</th>
          <th rowspan="16" valign="top">CEN</th>
          <th>Top</th>
          <td>fiber_distance</td>
          <td>-0.125</td>
          <td>-0.125</td>
          <td>-0.125</td>
          <td>-0.125</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>oxx</td>
          <td>-36570.457</td>
          <td>-158687.391</td>
          <td>-149706.203</td>
          <td>1068952.125</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>oyy</td>
          <td>206374.969</td>
          <td>1083602.750</td>
          <td>403245.969</td>
          <td>6158211.500</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>txy</td>
          <td>229.650</td>
          <td>-12673.086</td>
          <td>4394314.500</td>
          <td>-357167.656</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>angle</td>
          <td>89.946</td>
          <td>-89.416</td>
          <td>46.800</td>
          <td>-86.005</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>omax</td>
          <td>206375.188</td>
          <td>1083732.125</td>
          <td>4529773.000</td>
          <td>6183155.500</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>omin</td>
          <td>-36570.672</td>
          <td>-158816.656</td>
          <td>-4276233.500</td>
          <td>1044008.062</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>von_mises</td>
          <td>226881.938</td>
          <td>1171244.000</td>
          <td>7627279.000</td>
          <td>5732896.500</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>fiber_distance</td>
          <td>0.125</td>
          <td>0.125</td>
          <td>0.125</td>
          <td>0.125</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>oxx</td>
          <td>-28156.799</td>
          <td>-95551.906</td>
          <td>-194234.062</td>
          <td>-488197.969</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>oyy</td>
          <td>140208.719</td>
          <td>732509.188</td>
          <td>7016.848</td>
          <td>-278514.844</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>txy</td>
          <td>74085.039</td>
          <td>-35219.672</td>
          <td>4534850.000</td>
          <td>-353332.000</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>angle</td>
          <td>69.325</td>
          <td>-87.569</td>
          <td>45.636</td>
          <td>-53.263</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>omax</td>
          <td>168165.734</td>
          <td>734004.500</td>
          <td>4442357.500</td>
          <td>-14798.063</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>omin</td>
          <td>-56113.816</td>
          <td>-97047.195</td>
          <td>-4629575.000</td>
          <td>-751914.750</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>von_mises</td>
          <td>202150.672</td>
          <td>787028.500</td>
          <td>7857081.500</td>
          <td>744626.000</td>
        </tr>
        <tr>
          <th rowspan="4" valign="top">4</th>
          <th>Top</th>
          <td>fiber_distance</td>
          <td>-0.125</td>
          <td>-0.125</td>
          <td>-0.125</td>
          <td>-0.125</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>oxx</td>
          <td>-99755.844</td>
          <td>-580174.062</td>
          <td>-292532.719</td>
          <td>793623.688</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>oyy</td>
          <td>-1101563.000</td>
          <td>1460770.000</td>
          <td>-3137639.000</td>
          <td>6441436.000</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>txy</td>
          <td>229.650</td>
          <td>-12673.086</td>
          <td>4394314.500</td>
          <td>-357167.656</td>
        </tr>
      </tbody>
    </table>
    </div>


