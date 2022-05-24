
Static & Transient DataFrames in PyNastran
==========================================

The Jupyter notebook for this demo can be found in: -
docs/quick_start/demo/op2_pandas_multi_case.ipynb -
https://github.com/SteveDoyle2/pyNastran/tree/main/docs/quick_start/demo/op2_pandas_multi_case.ipynb

.. code:: ipython3

    import os
    import pandas as pd
    import pyNastran
    from pyNastran.op2.op2 import read_op2
    
    pkg_path = pyNastran.__path__[0]
    model_path = os.path.join(pkg_path, '..', 'models')

Solid Bending
-------------

Let’s show off ``combine=True/False``. We’ll talk about the keys soon.

.. code:: ipython3

    solid_bending_op2 = os.path.join(model_path, 'solid_bending', 'solid_bending.op2')
    solid_bending = read_op2(solid_bending_op2, combine=False, debug=False)
    print(solid_bending.displacements.keys())



.. raw:: html

    <text style=color:green>INFO:    op2_scalar.py:1459           op2_filename = 'c:\\nasa\\m4\\formats\\git\\pynastran\\pyNastran\\..\\models\\solid_bending\\solid_bending.op2'
    </text>


.. parsed-literal::

    dict_keys([(1, 1, 1, 0, 0, '', '')])
    

.. parsed-literal::

    c:\nasa\m4\formats\git\pynastran\pyNastran\op2\op2.py:740: FutureWarning: 
    Panel is deprecated and will be removed in a future version.
    The recommended way to represent these types of 3-dimensional data are with a MultiIndex on a DataFrame, via the Panel.to_frame() method
    Alternatively, you can use the xarray package http://xarray.pydata.org/en/stable/.
    Pandas provides a `.to_xarray()` method to help automate this conversion.
    
      obj.build_dataframe()
    

.. code:: ipython3

    solid_bending_op2 = os.path.join(model_path, 'solid_bending', 'solid_bending.op2')
    solid_bending2 = read_op2(solid_bending_op2, combine=True, debug=False)
    print(solid_bending2.displacements.keys())



.. raw:: html

    <text style=color:green>INFO:    op2_scalar.py:1459           op2_filename = 'c:\\nasa\\m4\\formats\\git\\pynastran\\pyNastran\\..\\models\\solid_bending\\solid_bending.op2'
    </text>


.. parsed-literal::

    dict_keys([1])
    

Single Subcase Buckling Example
-------------------------------

The keys cannot be “combined” despite us telling the program that it was
OK. We’ll get the following values that we need to handle. ####
isubcase, analysis_code, sort_method, count, subtitle \* isubcase -> the
same key that you’re used to accessing \* sort_method -> 1 (SORT1), 2
(SORT2) \* count -> the optimization count \* subtitle -> the analysis
subtitle (changes for superlements) \* analysis code -> the “type” of
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

Let’s look at an odd case:
~~~~~~~~~~~~~~~~~~~~~~~~~~

You can do buckling as one subcase or two subcases (makes parsing it a
lot easier!).

However, you **have** to do this once you start messing around with
superelements or multi-step optimization.

For optimization, sometimes Nastran will downselect elements and do an
optimization on that and print out a subset of the elements. At the end,
it will rerun an analysis to double check the constraints are satisfied.
It does not always do multi-step optimization.

.. code:: ipython3

    op2_filename = os.path.join(model_path, 'sol_101_elements', 'buckling_solid_shell_bar.op2')
    model = read_op2(op2_filename, combine=True, debug=False, build_dataframe=True)



.. raw:: html

    <text style=color:green>INFO:    op2_scalar.py:1459           op2_filename = 'c:\\nasa\\m4\\formats\\git\\pynastran\\pyNastran\\..\\models\\sol_101_elements\\buckling_solid_shell_bar.op2'
    </text>


.. code:: ipython3

    stress_keys = model.cquad4_stress.keys()
    print (stress_keys)
    
    # subcase, analysis_code, sort_method, count, isuperelmemnt_adaptivity_index, pval_step
    key0 = (1, 1, 1, 0, 0, '', '')
    key1 = (1, 8, 1, 0, 0, '', '')


.. parsed-literal::

    dict_keys([(1, 1, 1, 0, 0, '', ''), (1, 8, 1, 0, 0, '', '')])
    

Keys: \* key0 is the “static” key \* key1 is the “buckling” key

Similarly: \* Transient solutions can have preload \* Frequency
solutions can have loadsets (???)

Moving onto the data frames
---------------------------

-  The static case is the initial deflection state
-  The buckling case is “transient”, where the modes (called load steps
   or lsdvmn here) represent the “times”

pyNastran reads these tables differently and handles them differently
internally. They look very similar though.

.. code:: ipython3

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

    stress_static.nonlinear_factor = nan
    stress_transient.nonlinear_factor = 4
    data_names  = ['lsdvmn', 'eigr']
    loadsteps   = [1, 2, 3, 4]
    eigenvalues = [-49357660160.0, -58001940480.0, -379750744064.0, -428462538752.0]
    

Static Table
------------

.. code:: ipython3

    # Sets default precision of real numbers for pandas output\n"
    pd.set_option('precision', 2)
    
    stress_static.head(20)




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
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

.. code:: ipython3

    # Sets default precision of real numbers for pandas output\n"
    pd.set_option('precision', 3)
    #import numpy as np
    #np.set_printoptions(formatter={'all':lambda x: '%g'})
    
    stress_transient.head(20)




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
          <th>-379750744064.0</th>
          <th>-428462538752.0</th>
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
          <td>-1.250e-01</td>
          <td>-1.250e-01</td>
          <td>-1.250e-01</td>
          <td>-1.250e-01</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>oxx</td>
          <td>-3.657e+04</td>
          <td>-1.587e+05</td>
          <td>-1.497e+05</td>
          <td>1.069e+06</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>oyy</td>
          <td>2.064e+05</td>
          <td>1.084e+06</td>
          <td>4.032e+05</td>
          <td>6.158e+06</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>txy</td>
          <td>2.296e+02</td>
          <td>-1.267e+04</td>
          <td>4.394e+06</td>
          <td>-3.572e+05</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>angle</td>
          <td>8.995e+01</td>
          <td>-8.942e+01</td>
          <td>4.680e+01</td>
          <td>-8.601e+01</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>omax</td>
          <td>2.064e+05</td>
          <td>1.084e+06</td>
          <td>4.530e+06</td>
          <td>6.183e+06</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>omin</td>
          <td>-3.657e+04</td>
          <td>-1.588e+05</td>
          <td>-4.276e+06</td>
          <td>1.044e+06</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>von_mises</td>
          <td>2.269e+05</td>
          <td>1.171e+06</td>
          <td>7.627e+06</td>
          <td>5.733e+06</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>fiber_distance</td>
          <td>1.250e-01</td>
          <td>1.250e-01</td>
          <td>1.250e-01</td>
          <td>1.250e-01</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>oxx</td>
          <td>-2.816e+04</td>
          <td>-9.555e+04</td>
          <td>-1.942e+05</td>
          <td>-4.882e+05</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>oyy</td>
          <td>1.402e+05</td>
          <td>7.325e+05</td>
          <td>7.017e+03</td>
          <td>-2.785e+05</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>txy</td>
          <td>7.409e+04</td>
          <td>-3.522e+04</td>
          <td>4.535e+06</td>
          <td>-3.533e+05</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>angle</td>
          <td>6.933e+01</td>
          <td>-8.757e+01</td>
          <td>4.564e+01</td>
          <td>-5.326e+01</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>omax</td>
          <td>1.682e+05</td>
          <td>7.340e+05</td>
          <td>4.442e+06</td>
          <td>-1.480e+04</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>omin</td>
          <td>-5.611e+04</td>
          <td>-9.705e+04</td>
          <td>-4.630e+06</td>
          <td>-7.519e+05</td>
        </tr>
        <tr>
          <th>Bottom</th>
          <td>von_mises</td>
          <td>2.022e+05</td>
          <td>7.870e+05</td>
          <td>7.857e+06</td>
          <td>7.446e+05</td>
        </tr>
        <tr>
          <th rowspan="4" valign="top">4</th>
          <th>Top</th>
          <td>fiber_distance</td>
          <td>-1.250e-01</td>
          <td>-1.250e-01</td>
          <td>-1.250e-01</td>
          <td>-1.250e-01</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>oxx</td>
          <td>-9.976e+04</td>
          <td>-5.802e+05</td>
          <td>-2.925e+05</td>
          <td>7.936e+05</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>oyy</td>
          <td>-1.102e+06</td>
          <td>1.461e+06</td>
          <td>-3.138e+06</td>
          <td>6.441e+06</td>
        </tr>
        <tr>
          <th>Top</th>
          <td>txy</td>
          <td>2.296e+02</td>
          <td>-1.267e+04</td>
          <td>4.394e+06</td>
          <td>-3.572e+05</td>
        </tr>
      </tbody>
    </table>
    </div>


