Manipulating the Pandas DataFrame
=================================

The Jupyter notebook for this demo can be found in: -
docs/quick_start/demo/op2_pandas_unstack.ipynb -
https://github.com/SteveDoyle2/pyNastran/tree/main/docs/quick_start/demo/op2_pandas_unstack.ipynb

This example will use pandas unstack
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The unstack method on a DataFrame moves on index level from rows to
columns. First let’s read in some data:

.. code:: ipython3

    import os
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    from pyNastran.op2.op2 import read_op2
    import pandas as pd
    pd.set_option('precision', 2)

    op2_filename = os.path.join(pkg_path, '..', 'models', 'iSat', 'iSat_launch_100Hz.op2')
    from pyNastran.op2.op2 import read_op2
    isat = read_op2(op2_filename, build_dataframe=True, debug=False, skip_undefined_matrices=True)



.. raw:: html

    <text style=color:green>INFO:    op2_scalar.py:1588           op2_filename = 'c:\\nasa\\m4\\formats\\git\\pynastran\\pyNastran\\..\\models\\iSat\\iSat_launch_100Hz.op2'
    </text>


.. parsed-literal::

    self.cannot apply column_names=['Mode', 'Freq'] to RealStrainEnergyArray: 'QUAD4'
    self.cannot apply column_names=['Mode', 'Freq'] to RealStrainEnergyArray: 'TRIA3'
    self.cannot apply column_names=['Mode', 'Freq'] to RealStrainEnergyArray: 'HEXA'
    self.cannot apply column_names=['Mode', 'Freq'] to RealStrainEnergyArray: 'BAR'
    self.cannot apply column_names=['Mode', 'Freq'] to RealStrainEnergyArray: 'BUSH'


.. code:: ipython3

    cbar = isat.cbar_force[1].data_frame

.. code:: ipython3

    cbar.head()




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
          <th>24</th>
          <th>25</th>
          <th>26</th>
          <th>27</th>
          <th>28</th>
          <th>29</th>
          <th>30</th>
          <th>31</th>
          <th>32</th>
          <th>33</th>
        </tr>
        <tr>
          <th></th>
          <th>Freq</th>
          <th>8.36</th>
          <th>9.51</th>
          <th>15.67</th>
          <th>20.24</th>
          <th>20.31</th>
          <th>20.55</th>
          <th>21.50</th>
          <th>21.71</th>
          <th>21.72</th>
          <th>28.54</th>
          <th>...</th>
          <th>80.08</th>
          <th>86.49</th>
          <th>88.17</th>
          <th>88.48</th>
          <th>89.93</th>
          <th>94.29</th>
          <th>94.37</th>
          <th>96.04</th>
          <th>98.70</th>
          <th>98.89</th>
        </tr>
        <tr>
          <th></th>
          <th>Eigenvalue</th>
          <th>2758.15</th>
          <th>3568.63</th>
          <th>9689.98</th>
          <th>16168.04</th>
          <th>16278.16</th>
          <th>16679.71</th>
          <th>18248.43</th>
          <th>18600.70</th>
          <th>18632.55</th>
          <th>32159.89</th>
          <th>...</th>
          <th>253141.17</th>
          <th>295300.94</th>
          <th>306886.00</th>
          <th>309040.66</th>
          <th>319267.09</th>
          <th>350984.50</th>
          <th>351566.19</th>
          <th>364166.31</th>
          <th>384601.34</th>
          <th>386090.47</th>
        </tr>
        <tr>
          <th></th>
          <th>Radians</th>
          <th>52.52</th>
          <th>59.74</th>
          <th>98.44</th>
          <th>127.15</th>
          <th>127.59</th>
          <th>129.15</th>
          <th>135.09</th>
          <th>136.38</th>
          <th>136.50</th>
          <th>179.33</th>
          <th>...</th>
          <th>503.13</th>
          <th>543.42</th>
          <th>553.97</th>
          <th>555.91</th>
          <th>565.04</th>
          <th>592.44</th>
          <th>592.93</th>
          <th>603.46</th>
          <th>620.16</th>
          <th>621.36</th>
        </tr>
        <tr>
          <th>ElementID</th>
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
          <th rowspan="5" valign="top">3323</th>
          <th>bending_moment_a1</th>
          <td>-0.16</td>
          <td>-0.23</td>
          <td>-1.33</td>
          <td>-2.32e+00</td>
          <td>-1.88</td>
          <td>-0.80</td>
          <td>-1.34e-03</td>
          <td>1.42</td>
          <td>1.47</td>
          <td>4.65</td>
          <td>...</td>
          <td>-43.42</td>
          <td>63.36</td>
          <td>-43.07</td>
          <td>-3.35</td>
          <td>11.08</td>
          <td>-14.38</td>
          <td>0.75</td>
          <td>29.36</td>
          <td>0.49</td>
          <td>-4.56</td>
        </tr>
        <tr>
          <th>bending_moment_a2</th>
          <td>0.19</td>
          <td>0.05</td>
          <td>0.18</td>
          <td>5.58e-03</td>
          <td>-0.11</td>
          <td>-0.42</td>
          <td>-4.19e-03</td>
          <td>-1.11</td>
          <td>0.10</td>
          <td>-1.57</td>
          <td>...</td>
          <td>-4.50</td>
          <td>5.33</td>
          <td>1.63</td>
          <td>4.86</td>
          <td>2.15</td>
          <td>0.09</td>
          <td>-1.27</td>
          <td>-10.58</td>
          <td>-0.67</td>
          <td>3.48</td>
        </tr>
        <tr>
          <th>bending_moment_b1</th>
          <td>0.17</td>
          <td>0.21</td>
          <td>2.01</td>
          <td>2.66e+00</td>
          <td>1.88</td>
          <td>0.73</td>
          <td>2.29e-03</td>
          <td>-1.38</td>
          <td>-1.31</td>
          <td>-3.98</td>
          <td>...</td>
          <td>34.70</td>
          <td>-74.02</td>
          <td>35.13</td>
          <td>3.54</td>
          <td>-15.03</td>
          <td>10.97</td>
          <td>-0.67</td>
          <td>-17.69</td>
          <td>-0.63</td>
          <td>6.39</td>
        </tr>
        <tr>
          <th>bending_moment_b2</th>
          <td>-0.19</td>
          <td>-0.05</td>
          <td>-0.18</td>
          <td>-3.54e-03</td>
          <td>0.11</td>
          <td>0.43</td>
          <td>4.18e-03</td>
          <td>1.11</td>
          <td>-0.10</td>
          <td>1.57</td>
          <td>...</td>
          <td>4.50</td>
          <td>-5.34</td>
          <td>-1.62</td>
          <td>-4.86</td>
          <td>-2.15</td>
          <td>-0.08</td>
          <td>1.27</td>
          <td>10.56</td>
          <td>0.67</td>
          <td>-3.49</td>
        </tr>
        <tr>
          <th>shear1</th>
          <td>-0.13</td>
          <td>-0.18</td>
          <td>-1.33</td>
          <td>-1.99e+00</td>
          <td>-1.50</td>
          <td>-0.61</td>
          <td>-1.45e-03</td>
          <td>1.12</td>
          <td>1.11</td>
          <td>3.45</td>
          <td>...</td>
          <td>-31.25</td>
          <td>54.95</td>
          <td>-31.28</td>
          <td>-2.76</td>
          <td>10.44</td>
          <td>-10.14</td>
          <td>0.57</td>
          <td>18.82</td>
          <td>0.44</td>
          <td>-4.38</td>
        </tr>
      </tbody>
    </table>
    <p>5 rows × 33 columns</p>
    </div>



First I’m going to pull out a small subset to work with

.. code:: ipython3

    csub = cbar.loc[3323:3324,1:2]
    csub




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
        </tr>
        <tr>
          <th></th>
          <th>Freq</th>
          <th>8.36</th>
          <th>9.51</th>
        </tr>
        <tr>
          <th></th>
          <th>Eigenvalue</th>
          <th>2758.15</th>
          <th>3568.63</th>
        </tr>
        <tr>
          <th></th>
          <th>Radians</th>
          <th>52.52</th>
          <th>59.74</th>
        </tr>
        <tr>
          <th>ElementID</th>
          <th>Item</th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th rowspan="8" valign="top">3323</th>
          <th>bending_moment_a1</th>
          <td>-0.16</td>
          <td>-0.23</td>
        </tr>
        <tr>
          <th>bending_moment_a2</th>
          <td>0.19</td>
          <td>0.05</td>
        </tr>
        <tr>
          <th>bending_moment_b1</th>
          <td>0.17</td>
          <td>0.21</td>
        </tr>
        <tr>
          <th>bending_moment_b2</th>
          <td>-0.19</td>
          <td>-0.05</td>
        </tr>
        <tr>
          <th>shear1</th>
          <td>-0.13</td>
          <td>-0.18</td>
        </tr>
        <tr>
          <th>shear2</th>
          <td>0.15</td>
          <td>0.04</td>
        </tr>
        <tr>
          <th>axial</th>
          <td>0.80</td>
          <td>-0.21</td>
        </tr>
        <tr>
          <th>torque</th>
          <td>-0.04</td>
          <td>0.06</td>
        </tr>
        <tr>
          <th rowspan="8" valign="top">3324</th>
          <th>bending_moment_a1</th>
          <td>0.14</td>
          <td>0.29</td>
        </tr>
        <tr>
          <th>bending_moment_a2</th>
          <td>-0.19</td>
          <td>-0.05</td>
        </tr>
        <tr>
          <th>bending_moment_b1</th>
          <td>-0.15</td>
          <td>-0.26</td>
        </tr>
        <tr>
          <th>bending_moment_b2</th>
          <td>0.19</td>
          <td>0.05</td>
        </tr>
        <tr>
          <th>shear1</th>
          <td>0.12</td>
          <td>0.22</td>
        </tr>
        <tr>
          <th>shear2</th>
          <td>-0.15</td>
          <td>-0.04</td>
        </tr>
        <tr>
          <th>axial</th>
          <td>-0.80</td>
          <td>0.21</td>
        </tr>
        <tr>
          <th>torque</th>
          <td>0.04</td>
          <td>-0.06</td>
        </tr>
      </tbody>
    </table>
    </div>



I happen to like the way that’s organized, but let’s say that I want the
have the item descriptions in columns and the mode ID’s and element
numbers in rows. To do that, I’ll first move the element ID’s up to the
columns using a .unstack(level=0) and the transpose the result:

.. code:: ipython3

    csub.unstack(level=0).T




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
          <th></th>
          <th>Item</th>
          <th>axial</th>
          <th>bending_moment_a1</th>
          <th>bending_moment_a2</th>
          <th>bending_moment_b1</th>
          <th>bending_moment_b2</th>
          <th>shear1</th>
          <th>shear2</th>
          <th>torque</th>
        </tr>
        <tr>
          <th>Mode</th>
          <th>Freq</th>
          <th>Eigenvalue</th>
          <th>Radians</th>
          <th>ElementID</th>
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
          <th rowspan="2" valign="top">1</th>
          <th rowspan="2" valign="top">8.36</th>
          <th rowspan="2" valign="top">2758.15</th>
          <th rowspan="2" valign="top">52.52</th>
          <th>3323</th>
          <td>0.80</td>
          <td>-0.16</td>
          <td>0.19</td>
          <td>0.17</td>
          <td>-0.19</td>
          <td>-0.13</td>
          <td>0.15</td>
          <td>-0.04</td>
        </tr>
        <tr>
          <th>3324</th>
          <td>-0.80</td>
          <td>0.14</td>
          <td>-0.19</td>
          <td>-0.15</td>
          <td>0.19</td>
          <td>0.12</td>
          <td>-0.15</td>
          <td>0.04</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">2</th>
          <th rowspan="2" valign="top">9.51</th>
          <th rowspan="2" valign="top">3568.63</th>
          <th rowspan="2" valign="top">59.74</th>
          <th>3323</th>
          <td>-0.21</td>
          <td>-0.23</td>
          <td>0.05</td>
          <td>0.21</td>
          <td>-0.05</td>
          <td>-0.18</td>
          <td>0.04</td>
          <td>0.06</td>
        </tr>
        <tr>
          <th>3324</th>
          <td>0.21</td>
          <td>0.29</td>
          <td>-0.05</td>
          <td>-0.26</td>
          <td>0.05</td>
          <td>0.22</td>
          <td>-0.04</td>
          <td>-0.06</td>
        </tr>
      </tbody>
    </table>
    </div>



unstack requires unique row indices so I can’t work with CQUAD4 stresses
as they’re currently output, but I’ll work with CHEXA stresses. Let’s
pull out the first two elements and first two modes:

.. code:: ipython3

    chs = isat.chexa_stress[1].data_frame.loc[3684:3685,1:2]
    chs




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
          <th>Mode</th>
          <th>1</th>
          <th>2</th>
        </tr>
        <tr>
          <th></th>
          <th></th>
          <th>Freq</th>
          <th>8.36</th>
          <th>9.51</th>
        </tr>
        <tr>
          <th></th>
          <th></th>
          <th>Eigenvalue</th>
          <th>2758.15</th>
          <th>3568.63</th>
        </tr>
        <tr>
          <th></th>
          <th></th>
          <th>Radians</th>
          <th>52.52</th>
          <th>59.74</th>
        </tr>
        <tr>
          <th>ElementID</th>
          <th>NodeID</th>
          <th>Item</th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th rowspan="5" valign="top">3684</th>
          <th rowspan="5" valign="top">0</th>
          <th>oxx</th>
          <td>1.22e-12</td>
          <td>-3.41e-13</td>
        </tr>
        <tr>
          <th>oyy</th>
          <td>-3.35e-12</td>
          <td>-2.27e-13</td>
        </tr>
        <tr>
          <th>ozz</th>
          <td>1.25e-12</td>
          <td>4.55e-13</td>
        </tr>
        <tr>
          <th>txy</th>
          <td>-3.27e-13</td>
          <td>1.82e-12</td>
        </tr>
        <tr>
          <th>tyz</th>
          <td>2.84e-14</td>
          <td>3.98e-13</td>
        </tr>
        <tr>
          <th>...</th>
          <th>...</th>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th rowspan="5" valign="top">3685</th>
          <th rowspan="5" valign="top">1037</th>
          <th>txz</th>
          <td>-2.84e-13</td>
          <td>-1.82e-12</td>
        </tr>
        <tr>
          <th>omax</th>
          <td>-7.47e-15</td>
          <td>2.08e-12</td>
        </tr>
        <tr>
          <th>omid</th>
          <td>-1.15e-13</td>
          <td>-2.71e-13</td>
        </tr>
        <tr>
          <th>omin</th>
          <td>-1.00e-12</td>
          <td>-1.70e-12</td>
        </tr>
        <tr>
          <th>von_mises</th>
          <td>9.43e-13</td>
          <td>3.30e-12</td>
        </tr>
      </tbody>
    </table>
    <p>180 rows × 2 columns</p>
    </div>



Now I want to put ElementID and the Node ID in the rows along with the
Load ID, and have the items in the columns:

.. code:: ipython3

    cht = chs.unstack(level=[0,1]).T
    cht




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
          <th></th>
          <th></th>
          <th>Item</th>
          <th>omax</th>
          <th>omid</th>
          <th>omin</th>
          <th>oxx</th>
          <th>oyy</th>
          <th>ozz</th>
          <th>txy</th>
          <th>txz</th>
          <th>tyz</th>
          <th>von_mises</th>
        </tr>
        <tr>
          <th>Mode</th>
          <th>Freq</th>
          <th>Eigenvalue</th>
          <th>Radians</th>
          <th>ElementID</th>
          <th>NodeID</th>
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
          <th rowspan="5" valign="top">8.36</th>
          <th rowspan="5" valign="top">2758.15</th>
          <th rowspan="5" valign="top">52.52</th>
          <th rowspan="5" valign="top">3684</th>
          <th>0</th>
          <td>1.48e-12</td>
          <td>1.02e-12</td>
          <td>-3.38e-12</td>
          <td>1.22e-12</td>
          <td>-3.35e-12</td>
          <td>1.25e-12</td>
          <td>-3.27e-13</td>
          <td>-2.27e-13</td>
          <td>2.84e-14</td>
          <td>4.64e-12</td>
        </tr>
        <tr>
          <th>55</th>
          <td>4.81e-12</td>
          <td>1.92e-13</td>
          <td>-3.57e-13</td>
          <td>4.53e-12</td>
          <td>-2.42e-13</td>
          <td>3.55e-13</td>
          <td>-6.54e-13</td>
          <td>-9.30e-13</td>
          <td>2.20e-15</td>
          <td>4.92e-12</td>
        </tr>
        <tr>
          <th>51</th>
          <td>2.32e-12</td>
          <td>1.49e-13</td>
          <td>-1.41e-12</td>
          <td>-1.39e-12</td>
          <td>2.32e-12</td>
          <td>1.35e-13</td>
          <td>-1.30e-13</td>
          <td>-1.46e-13</td>
          <td>7.51e-15</td>
          <td>3.25e-12</td>
        </tr>
        <tr>
          <th>778</th>
          <td>-1.38e-12</td>
          <td>-3.27e-12</td>
          <td>-6.12e-12</td>
          <td>-6.08e-12</td>
          <td>-1.38e-12</td>
          <td>-3.31e-12</td>
          <td>-5.81e-14</td>
          <td>-3.41e-13</td>
          <td>-1.97e-14</td>
          <td>4.14e-12</td>
        </tr>
        <tr>
          <th>758</th>
          <td>5.79e-12</td>
          <td>4.11e-12</td>
          <td>7.57e-14</td>
          <td>5.68e-12</td>
          <td>1.14e-13</td>
          <td>4.18e-12</td>
          <td>-4.55e-13</td>
          <td>-3.41e-13</td>
          <td>-3.91e-14</td>
          <td>5.09e-12</td>
        </tr>
        <tr>
          <th>...</th>
          <th>...</th>
          <th>...</th>
          <th>...</th>
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
        </tr>
        <tr>
          <th rowspan="5" valign="top">2</th>
          <th rowspan="5" valign="top">9.51</th>
          <th rowspan="5" valign="top">3568.63</th>
          <th rowspan="5" valign="top">59.74</th>
          <th rowspan="5" valign="top">3685</th>
          <th>1015</th>
          <td>1.19e-12</td>
          <td>-9.78e-14</td>
          <td>-6.96e-13</td>
          <td>3.41e-13</td>
          <td>-5.68e-14</td>
          <td>1.14e-13</td>
          <td>-2.27e-13</td>
          <td>-9.09e-13</td>
          <td>1.13e-13</td>
          <td>1.67e-12</td>
        </tr>
        <tr>
          <th>50</th>
          <td>4.88e-13</td>
          <td>-1.28e-13</td>
          <td>-1.21e-12</td>
          <td>-3.98e-13</td>
          <td>1.14e-13</td>
          <td>-5.68e-13</td>
          <td>-5.68e-13</td>
          <td>-4.90e-13</td>
          <td>-2.27e-13</td>
          <td>1.49e-12</td>
        </tr>
        <tr>
          <th>46</th>
          <td>6.22e-13</td>
          <td>2.26e-14</td>
          <td>-5.87e-13</td>
          <td>-2.84e-13</td>
          <td>3.41e-13</td>
          <td>0.00e+00</td>
          <td>4.54e-13</td>
          <td>9.59e-14</td>
          <td>-2.27e-13</td>
          <td>1.05e-12</td>
        </tr>
        <tr>
          <th>1031</th>
          <td>2.32e-12</td>
          <td>-6.90e-13</td>
          <td>-1.63e-12</td>
          <td>-2.27e-13</td>
          <td>-6.82e-13</td>
          <td>9.09e-13</td>
          <td>4.55e-13</td>
          <td>-1.82e-12</td>
          <td>-2.27e-13</td>
          <td>3.57e-12</td>
        </tr>
        <tr>
          <th>1037</th>
          <td>2.08e-12</td>
          <td>-2.71e-13</td>
          <td>-1.70e-12</td>
          <td>4.55e-13</td>
          <td>-3.41e-13</td>
          <td>0.00e+00</td>
          <td>-3.98e-13</td>
          <td>-1.82e-12</td>
          <td>-1.14e-13</td>
          <td>3.30e-12</td>
        </tr>
      </tbody>
    </table>
    <p>68 rows × 10 columns</p>
    </div>



Maybe I’d like my rows organized with the modes on the inside. I can do
that by swapping levels:

We actually need to get rid of the extra rows using dropna():

.. code:: ipython3

    cht = cht.dropna()
    cht




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
          <th></th>
          <th></th>
          <th>Item</th>
          <th>omax</th>
          <th>omid</th>
          <th>omin</th>
          <th>oxx</th>
          <th>oyy</th>
          <th>ozz</th>
          <th>txy</th>
          <th>txz</th>
          <th>tyz</th>
          <th>von_mises</th>
        </tr>
        <tr>
          <th>Mode</th>
          <th>Freq</th>
          <th>Eigenvalue</th>
          <th>Radians</th>
          <th>ElementID</th>
          <th>NodeID</th>
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
          <th rowspan="18" valign="top">1</th>
          <th rowspan="18" valign="top">8.36</th>
          <th rowspan="18" valign="top">2758.15</th>
          <th rowspan="18" valign="top">52.52</th>
          <th rowspan="9" valign="top">3684</th>
          <th>0</th>
          <td>1.48e-12</td>
          <td>1.02e-12</td>
          <td>-3.38e-12</td>
          <td>1.22e-12</td>
          <td>-3.35e-12</td>
          <td>1.25e-12</td>
          <td>-3.27e-13</td>
          <td>-2.27e-13</td>
          <td>2.84e-14</td>
          <td>4.64e-12</td>
        </tr>
        <tr>
          <th>55</th>
          <td>4.81e-12</td>
          <td>1.92e-13</td>
          <td>-3.57e-13</td>
          <td>4.53e-12</td>
          <td>-2.42e-13</td>
          <td>3.55e-13</td>
          <td>-6.54e-13</td>
          <td>-9.30e-13</td>
          <td>2.20e-15</td>
          <td>4.92e-12</td>
        </tr>
        <tr>
          <th>51</th>
          <td>2.32e-12</td>
          <td>1.49e-13</td>
          <td>-1.41e-12</td>
          <td>-1.39e-12</td>
          <td>2.32e-12</td>
          <td>1.35e-13</td>
          <td>-1.30e-13</td>
          <td>-1.46e-13</td>
          <td>7.51e-15</td>
          <td>3.25e-12</td>
        </tr>
        <tr>
          <th>778</th>
          <td>-1.38e-12</td>
          <td>-3.27e-12</td>
          <td>-6.12e-12</td>
          <td>-6.08e-12</td>
          <td>-1.38e-12</td>
          <td>-3.31e-12</td>
          <td>-5.81e-14</td>
          <td>-3.41e-13</td>
          <td>-1.97e-14</td>
          <td>4.14e-12</td>
        </tr>
        <tr>
          <th>758</th>
          <td>5.79e-12</td>
          <td>4.11e-12</td>
          <td>7.57e-14</td>
          <td>5.68e-12</td>
          <td>1.14e-13</td>
          <td>4.18e-12</td>
          <td>-4.55e-13</td>
          <td>-3.41e-13</td>
          <td>-3.91e-14</td>
          <td>5.09e-12</td>
        </tr>
        <tr>
          <th>60</th>
          <td>2.88e-12</td>
          <td>1.71e-12</td>
          <td>-4.27e-12</td>
          <td>2.63e-12</td>
          <td>1.90e-12</td>
          <td>-4.21e-12</td>
          <td>-4.26e-13</td>
          <td>-6.34e-13</td>
          <td>8.53e-14</td>
          <td>6.64e-12</td>
        </tr>
        <tr>
          <th>56</th>
          <td>1.66e-12</td>
          <td>-1.65e-12</td>
          <td>-5.92e-12</td>
          <td>-5.87e-12</td>
          <td>1.65e-12</td>
          <td>-1.69e-12</td>
          <td>-3.07e-13</td>
          <td>-4.22e-13</td>
          <td>8.53e-14</td>
          <td>6.59e-12</td>
        </tr>
        <tr>
          <th>880</th>
          <td>2.63e-12</td>
          <td>2.11e-12</td>
          <td>-4.74e-12</td>
          <td>-4.72e-12</td>
          <td>2.10e-12</td>
          <td>2.61e-12</td>
          <td>-3.41e-13</td>
          <td>0.00e+00</td>
          <td>8.53e-14</td>
          <td>7.12e-12</td>
        </tr>
        <tr>
          <th>952</th>
          <td>-8.35e-15</td>
          <td>-1.79e-12</td>
          <td>-3.00e-12</td>
          <td>-1.73e-12</td>
          <td>-1.14e-13</td>
          <td>-2.96e-12</td>
          <td>-4.26e-13</td>
          <td>2.27e-13</td>
          <td>5.68e-14</td>
          <td>2.61e-12</td>
        </tr>
        <tr>
          <th rowspan="9" valign="top">3685</th>
          <th>0</th>
          <td>5.11e-13</td>
          <td>1.44e-13</td>
          <td>-5.41e-13</td>
          <td>4.97e-13</td>
          <td>1.56e-13</td>
          <td>-5.40e-13</td>
          <td>6.75e-14</td>
          <td>-2.84e-14</td>
          <td>1.42e-14</td>
          <td>9.25e-13</td>
        </tr>
        <tr>
          <th>45</th>
          <td>8.09e-13</td>
          <td>-5.47e-13</td>
          <td>-8.09e-13</td>
          <td>6.79e-13</td>
          <td>-7.89e-13</td>
          <td>-4.37e-13</td>
          <td>6.39e-14</td>
          <td>-4.02e-13</td>
          <td>5.15e-14</td>
          <td>1.50e-12</td>
        </tr>
        <tr>
          <th>41</th>
          <td>7.67e-13</td>
          <td>-1.49e-13</td>
          <td>-4.37e-13</td>
          <td>-1.49e-13</td>
          <td>7.21e-13</td>
          <td>-3.91e-13</td>
          <td>1.92e-13</td>
          <td>1.14e-13</td>
          <td>5.33e-14</td>
          <td>1.09e-12</td>
        </tr>
        <tr>
          <th>1021</th>
          <td>7.75e-13</td>
          <td>1.86e-13</td>
          <td>-5.38e-13</td>
          <td>1.56e-13</td>
          <td>1.60e-13</td>
          <td>1.07e-13</td>
          <td>1.94e-13</td>
          <td>-6.25e-13</td>
          <td>5.47e-14</td>
          <td>1.14e-12</td>
        </tr>
        <tr>
          <th>1015</th>
          <td>8.67e-13</td>
          <td>1.56e-13</td>
          <td>-6.65e-13</td>
          <td>2.98e-13</td>
          <td>1.53e-13</td>
          <td>-9.24e-14</td>
          <td>2.84e-14</td>
          <td>-7.39e-13</td>
          <td>4.51e-14</td>
          <td>1.33e-12</td>
        </tr>
        <tr>
          <th>50</th>
          <td>9.64e-13</td>
          <td>7.11e-13</td>
          <td>-6.57e-14</td>
          <td>-3.91e-14</td>
          <td>7.11e-13</td>
          <td>9.38e-13</td>
          <td>-4.26e-14</td>
          <td>-1.57e-13</td>
          <td>1.42e-14</td>
          <td>9.29e-13</td>
        </tr>
        <tr>
          <th>46</th>
          <td>2.26e-12</td>
          <td>7.91e-13</td>
          <td>7.55e-14</td>
          <td>2.25e-12</td>
          <td>7.96e-13</td>
          <td>7.82e-14</td>
          <td>8.37e-14</td>
          <td>-7.01e-14</td>
          <td>1.42e-14</td>
          <td>1.93e-12</td>
        </tr>
        <tr>
          <th>1031</th>
          <td>-8.72e-13</td>
          <td>-1.24e-12</td>
          <td>-1.87e-12</td>
          <td>-1.68e-12</td>
          <td>-8.81e-13</td>
          <td>-1.42e-12</td>
          <td>8.53e-14</td>
          <td>-2.84e-13</td>
          <td>4.26e-14</td>
          <td>8.75e-13</td>
        </tr>
        <tr>
          <th>1037</th>
          <td>-7.47e-15</td>
          <td>-1.15e-13</td>
          <td>-1.00e-12</td>
          <td>-9.95e-14</td>
          <td>-1.14e-13</td>
          <td>-9.09e-13</td>
          <td>-7.11e-15</td>
          <td>-2.84e-13</td>
          <td>2.84e-14</td>
          <td>9.43e-13</td>
        </tr>
        <tr>
          <th rowspan="18" valign="top">2</th>
          <th rowspan="18" valign="top">9.51</th>
          <th rowspan="18" valign="top">3568.63</th>
          <th rowspan="18" valign="top">59.74</th>
          <th rowspan="9" valign="top">3684</th>
          <th>0</th>
          <td>2.22e-12</td>
          <td>-1.06e-13</td>
          <td>-2.23e-12</td>
          <td>-3.41e-13</td>
          <td>-2.27e-13</td>
          <td>4.55e-13</td>
          <td>1.82e-12</td>
          <td>1.14e-12</td>
          <td>3.98e-13</td>
          <td>3.85e-12</td>
        </tr>
        <tr>
          <th>55</th>
          <td>3.64e-12</td>
          <td>9.55e-13</td>
          <td>-1.33e-12</td>
          <td>2.33e-12</td>
          <td>-5.68e-14</td>
          <td>9.95e-13</td>
          <td>1.82e-12</td>
          <td>4.26e-13</td>
          <td>1.12e-12</td>
          <td>4.31e-12</td>
        </tr>
        <tr>
          <th>51</th>
          <td>7.48e-13</td>
          <td>-3.94e-13</td>
          <td>-1.46e-12</td>
          <td>-3.41e-13</td>
          <td>5.68e-14</td>
          <td>-8.24e-13</td>
          <td>1.57e-13</td>
          <td>3.85e-13</td>
          <td>9.24e-13</td>
          <td>1.91e-12</td>
        </tr>
        <tr>
          <th>778</th>
          <td>1.24e-12</td>
          <td>-5.76e-13</td>
          <td>-2.37e-12</td>
          <td>-2.27e-12</td>
          <td>7.96e-13</td>
          <td>-2.27e-13</td>
          <td>1.25e-13</td>
          <td>4.55e-13</td>
          <td>7.73e-13</td>
          <td>3.13e-12</td>
        </tr>
        <tr>
          <th>758</th>
          <td>1.03e-12</td>
          <td>-1.04e-12</td>
          <td>-3.57e-12</td>
          <td>-5.68e-13</td>
          <td>-2.33e-12</td>
          <td>-6.82e-13</td>
          <td>1.82e-12</td>
          <td>4.55e-13</td>
          <td>9.08e-13</td>
          <td>3.99e-12</td>
        </tr>
        <tr>
          <th>60</th>
          <td>1.02e-12</td>
          <td>-2.20e-12</td>
          <td>-2.85e-12</td>
          <td>-1.48e-12</td>
          <td>-3.41e-13</td>
          <td>-2.22e-12</td>
          <td>1.82e-12</td>
          <td>2.64e-13</td>
          <td>1.14e-13</td>
          <td>3.59e-12</td>
        </tr>
        <tr>
          <th>56</th>
          <td>5.04e-13</td>
          <td>-3.52e-13</td>
          <td>-2.37e-12</td>
          <td>-1.82e-12</td>
          <td>-7.96e-13</td>
          <td>3.98e-13</td>
          <td>-9.10e-13</td>
          <td>-3.23e-13</td>
          <td>1.14e-13</td>
          <td>2.55e-12</td>
        </tr>
        <tr>
          <th>880</th>
          <td>1.79e-12</td>
          <td>6.09e-13</td>
          <td>-2.51e-12</td>
          <td>1.25e-12</td>
          <td>-2.27e-12</td>
          <td>9.09e-13</td>
          <td>-9.09e-13</td>
          <td>-4.55e-13</td>
          <td>3.41e-13</td>
          <td>3.85e-12</td>
        </tr>
        <tr>
          <th>952</th>
          <td>1.30e-12</td>
          <td>1.93e-13</td>
          <td>-2.63e-12</td>
          <td>-9.09e-13</td>
          <td>-6.82e-13</td>
          <td>4.55e-13</td>
          <td>1.82e-12</td>
          <td>4.55e-13</td>
          <td>2.27e-13</td>
          <td>3.51e-12</td>
        </tr>
        <tr>
          <th rowspan="9" valign="top">3685</th>
          <th>0</th>
          <td>1.18e-12</td>
          <td>-3.91e-15</td>
          <td>-9.47e-13</td>
          <td>-1.14e-13</td>
          <td>1.14e-13</td>
          <td>2.27e-13</td>
          <td>4.97e-13</td>
          <td>-9.09e-13</td>
          <td>-1.71e-13</td>
          <td>1.84e-12</td>
        </tr>
        <tr>
          <th>45</th>
          <td>8.96e-13</td>
          <td>-2.21e-13</td>
          <td>-1.07e-12</td>
          <td>2.27e-13</td>
          <td>-5.68e-14</td>
          <td>-5.68e-13</td>
          <td>-2.84e-13</td>
          <td>-6.52e-13</td>
          <td>5.54e-13</td>
          <td>1.71e-12</td>
        </tr>
        <tr>
          <th>41</th>
          <td>8.01e-13</td>
          <td>-3.97e-13</td>
          <td>-8.59e-13</td>
          <td>-2.84e-13</td>
          <td>0.00e+00</td>
          <td>-1.71e-13</td>
          <td>-4.62e-13</td>
          <td>-6.12e-13</td>
          <td>3.56e-13</td>
          <td>1.48e-12</td>
        </tr>
        <tr>
          <th>1021</th>
          <td>9.04e-13</td>
          <td>4.58e-13</td>
          <td>-3.01e-12</td>
          <td>-1.48e-12</td>
          <td>5.68e-13</td>
          <td>-7.39e-13</td>
          <td>4.50e-13</td>
          <td>-1.82e-12</td>
          <td>9.06e-14</td>
          <td>3.71e-12</td>
        </tr>
        <tr>
          <th>1015</th>
          <td>1.19e-12</td>
          <td>-9.78e-14</td>
          <td>-6.96e-13</td>
          <td>3.41e-13</td>
          <td>-5.68e-14</td>
          <td>1.14e-13</td>
          <td>-2.27e-13</td>
          <td>-9.09e-13</td>
          <td>1.13e-13</td>
          <td>1.67e-12</td>
        </tr>
        <tr>
          <th>50</th>
          <td>4.88e-13</td>
          <td>-1.28e-13</td>
          <td>-1.21e-12</td>
          <td>-3.98e-13</td>
          <td>1.14e-13</td>
          <td>-5.68e-13</td>
          <td>-5.68e-13</td>
          <td>-4.90e-13</td>
          <td>-2.27e-13</td>
          <td>1.49e-12</td>
        </tr>
        <tr>
          <th>46</th>
          <td>6.22e-13</td>
          <td>2.26e-14</td>
          <td>-5.87e-13</td>
          <td>-2.84e-13</td>
          <td>3.41e-13</td>
          <td>0.00e+00</td>
          <td>4.54e-13</td>
          <td>9.59e-14</td>
          <td>-2.27e-13</td>
          <td>1.05e-12</td>
        </tr>
        <tr>
          <th>1031</th>
          <td>2.32e-12</td>
          <td>-6.90e-13</td>
          <td>-1.63e-12</td>
          <td>-2.27e-13</td>
          <td>-6.82e-13</td>
          <td>9.09e-13</td>
          <td>4.55e-13</td>
          <td>-1.82e-12</td>
          <td>-2.27e-13</td>
          <td>3.57e-12</td>
        </tr>
        <tr>
          <th>1037</th>
          <td>2.08e-12</td>
          <td>-2.71e-13</td>
          <td>-1.70e-12</td>
          <td>4.55e-13</td>
          <td>-3.41e-13</td>
          <td>0.00e+00</td>
          <td>-3.98e-13</td>
          <td>-1.82e-12</td>
          <td>-1.14e-13</td>
          <td>3.30e-12</td>
        </tr>
      </tbody>
    </table>
    </div>



.. code:: ipython3

    # mode, eigr, freq, rad, eids, nids # initial
    # nids, eids, eigr, freq, rad, mode # final

    cht.swaplevel(0,4).swaplevel(1,5).swaplevel(2,5).swaplevel(4, 5)




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
          <th></th>
          <th></th>
          <th>Item</th>
          <th>omax</th>
          <th>omid</th>
          <th>omin</th>
          <th>oxx</th>
          <th>oyy</th>
          <th>ozz</th>
          <th>txy</th>
          <th>txz</th>
          <th>tyz</th>
          <th>von_mises</th>
        </tr>
        <tr>
          <th>ElementID</th>
          <th>NodeID</th>
          <th>Freq</th>
          <th>Radians</th>
          <th>Eigenvalue</th>
          <th>Mode</th>
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
          <th rowspan="9" valign="top">3684</th>
          <th>0</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>1.48e-12</td>
          <td>1.02e-12</td>
          <td>-3.38e-12</td>
          <td>1.22e-12</td>
          <td>-3.35e-12</td>
          <td>1.25e-12</td>
          <td>-3.27e-13</td>
          <td>-2.27e-13</td>
          <td>2.84e-14</td>
          <td>4.64e-12</td>
        </tr>
        <tr>
          <th>55</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>4.81e-12</td>
          <td>1.92e-13</td>
          <td>-3.57e-13</td>
          <td>4.53e-12</td>
          <td>-2.42e-13</td>
          <td>3.55e-13</td>
          <td>-6.54e-13</td>
          <td>-9.30e-13</td>
          <td>2.20e-15</td>
          <td>4.92e-12</td>
        </tr>
        <tr>
          <th>51</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>2.32e-12</td>
          <td>1.49e-13</td>
          <td>-1.41e-12</td>
          <td>-1.39e-12</td>
          <td>2.32e-12</td>
          <td>1.35e-13</td>
          <td>-1.30e-13</td>
          <td>-1.46e-13</td>
          <td>7.51e-15</td>
          <td>3.25e-12</td>
        </tr>
        <tr>
          <th>778</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>-1.38e-12</td>
          <td>-3.27e-12</td>
          <td>-6.12e-12</td>
          <td>-6.08e-12</td>
          <td>-1.38e-12</td>
          <td>-3.31e-12</td>
          <td>-5.81e-14</td>
          <td>-3.41e-13</td>
          <td>-1.97e-14</td>
          <td>4.14e-12</td>
        </tr>
        <tr>
          <th>758</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>5.79e-12</td>
          <td>4.11e-12</td>
          <td>7.57e-14</td>
          <td>5.68e-12</td>
          <td>1.14e-13</td>
          <td>4.18e-12</td>
          <td>-4.55e-13</td>
          <td>-3.41e-13</td>
          <td>-3.91e-14</td>
          <td>5.09e-12</td>
        </tr>
        <tr>
          <th>60</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>2.88e-12</td>
          <td>1.71e-12</td>
          <td>-4.27e-12</td>
          <td>2.63e-12</td>
          <td>1.90e-12</td>
          <td>-4.21e-12</td>
          <td>-4.26e-13</td>
          <td>-6.34e-13</td>
          <td>8.53e-14</td>
          <td>6.64e-12</td>
        </tr>
        <tr>
          <th>56</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>1.66e-12</td>
          <td>-1.65e-12</td>
          <td>-5.92e-12</td>
          <td>-5.87e-12</td>
          <td>1.65e-12</td>
          <td>-1.69e-12</td>
          <td>-3.07e-13</td>
          <td>-4.22e-13</td>
          <td>8.53e-14</td>
          <td>6.59e-12</td>
        </tr>
        <tr>
          <th>880</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>2.63e-12</td>
          <td>2.11e-12</td>
          <td>-4.74e-12</td>
          <td>-4.72e-12</td>
          <td>2.10e-12</td>
          <td>2.61e-12</td>
          <td>-3.41e-13</td>
          <td>0.00e+00</td>
          <td>8.53e-14</td>
          <td>7.12e-12</td>
        </tr>
        <tr>
          <th>952</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>-8.35e-15</td>
          <td>-1.79e-12</td>
          <td>-3.00e-12</td>
          <td>-1.73e-12</td>
          <td>-1.14e-13</td>
          <td>-2.96e-12</td>
          <td>-4.26e-13</td>
          <td>2.27e-13</td>
          <td>5.68e-14</td>
          <td>2.61e-12</td>
        </tr>
        <tr>
          <th rowspan="9" valign="top">3685</th>
          <th>0</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>5.11e-13</td>
          <td>1.44e-13</td>
          <td>-5.41e-13</td>
          <td>4.97e-13</td>
          <td>1.56e-13</td>
          <td>-5.40e-13</td>
          <td>6.75e-14</td>
          <td>-2.84e-14</td>
          <td>1.42e-14</td>
          <td>9.25e-13</td>
        </tr>
        <tr>
          <th>45</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>8.09e-13</td>
          <td>-5.47e-13</td>
          <td>-8.09e-13</td>
          <td>6.79e-13</td>
          <td>-7.89e-13</td>
          <td>-4.37e-13</td>
          <td>6.39e-14</td>
          <td>-4.02e-13</td>
          <td>5.15e-14</td>
          <td>1.50e-12</td>
        </tr>
        <tr>
          <th>41</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>7.67e-13</td>
          <td>-1.49e-13</td>
          <td>-4.37e-13</td>
          <td>-1.49e-13</td>
          <td>7.21e-13</td>
          <td>-3.91e-13</td>
          <td>1.92e-13</td>
          <td>1.14e-13</td>
          <td>5.33e-14</td>
          <td>1.09e-12</td>
        </tr>
        <tr>
          <th>1021</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>7.75e-13</td>
          <td>1.86e-13</td>
          <td>-5.38e-13</td>
          <td>1.56e-13</td>
          <td>1.60e-13</td>
          <td>1.07e-13</td>
          <td>1.94e-13</td>
          <td>-6.25e-13</td>
          <td>5.47e-14</td>
          <td>1.14e-12</td>
        </tr>
        <tr>
          <th>1015</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>8.67e-13</td>
          <td>1.56e-13</td>
          <td>-6.65e-13</td>
          <td>2.98e-13</td>
          <td>1.53e-13</td>
          <td>-9.24e-14</td>
          <td>2.84e-14</td>
          <td>-7.39e-13</td>
          <td>4.51e-14</td>
          <td>1.33e-12</td>
        </tr>
        <tr>
          <th>50</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>9.64e-13</td>
          <td>7.11e-13</td>
          <td>-6.57e-14</td>
          <td>-3.91e-14</td>
          <td>7.11e-13</td>
          <td>9.38e-13</td>
          <td>-4.26e-14</td>
          <td>-1.57e-13</td>
          <td>1.42e-14</td>
          <td>9.29e-13</td>
        </tr>
        <tr>
          <th>46</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>2.26e-12</td>
          <td>7.91e-13</td>
          <td>7.55e-14</td>
          <td>2.25e-12</td>
          <td>7.96e-13</td>
          <td>7.82e-14</td>
          <td>8.37e-14</td>
          <td>-7.01e-14</td>
          <td>1.42e-14</td>
          <td>1.93e-12</td>
        </tr>
        <tr>
          <th>1031</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>-8.72e-13</td>
          <td>-1.24e-12</td>
          <td>-1.87e-12</td>
          <td>-1.68e-12</td>
          <td>-8.81e-13</td>
          <td>-1.42e-12</td>
          <td>8.53e-14</td>
          <td>-2.84e-13</td>
          <td>4.26e-14</td>
          <td>8.75e-13</td>
        </tr>
        <tr>
          <th>1037</th>
          <th>8.36</th>
          <th>52.52</th>
          <th>2758.15</th>
          <th>1</th>
          <td>-7.47e-15</td>
          <td>-1.15e-13</td>
          <td>-1.00e-12</td>
          <td>-9.95e-14</td>
          <td>-1.14e-13</td>
          <td>-9.09e-13</td>
          <td>-7.11e-15</td>
          <td>-2.84e-13</td>
          <td>2.84e-14</td>
          <td>9.43e-13</td>
        </tr>
        <tr>
          <th rowspan="9" valign="top">3684</th>
          <th>0</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>2.22e-12</td>
          <td>-1.06e-13</td>
          <td>-2.23e-12</td>
          <td>-3.41e-13</td>
          <td>-2.27e-13</td>
          <td>4.55e-13</td>
          <td>1.82e-12</td>
          <td>1.14e-12</td>
          <td>3.98e-13</td>
          <td>3.85e-12</td>
        </tr>
        <tr>
          <th>55</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>3.64e-12</td>
          <td>9.55e-13</td>
          <td>-1.33e-12</td>
          <td>2.33e-12</td>
          <td>-5.68e-14</td>
          <td>9.95e-13</td>
          <td>1.82e-12</td>
          <td>4.26e-13</td>
          <td>1.12e-12</td>
          <td>4.31e-12</td>
        </tr>
        <tr>
          <th>51</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>7.48e-13</td>
          <td>-3.94e-13</td>
          <td>-1.46e-12</td>
          <td>-3.41e-13</td>
          <td>5.68e-14</td>
          <td>-8.24e-13</td>
          <td>1.57e-13</td>
          <td>3.85e-13</td>
          <td>9.24e-13</td>
          <td>1.91e-12</td>
        </tr>
        <tr>
          <th>778</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>1.24e-12</td>
          <td>-5.76e-13</td>
          <td>-2.37e-12</td>
          <td>-2.27e-12</td>
          <td>7.96e-13</td>
          <td>-2.27e-13</td>
          <td>1.25e-13</td>
          <td>4.55e-13</td>
          <td>7.73e-13</td>
          <td>3.13e-12</td>
        </tr>
        <tr>
          <th>758</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>1.03e-12</td>
          <td>-1.04e-12</td>
          <td>-3.57e-12</td>
          <td>-5.68e-13</td>
          <td>-2.33e-12</td>
          <td>-6.82e-13</td>
          <td>1.82e-12</td>
          <td>4.55e-13</td>
          <td>9.08e-13</td>
          <td>3.99e-12</td>
        </tr>
        <tr>
          <th>60</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>1.02e-12</td>
          <td>-2.20e-12</td>
          <td>-2.85e-12</td>
          <td>-1.48e-12</td>
          <td>-3.41e-13</td>
          <td>-2.22e-12</td>
          <td>1.82e-12</td>
          <td>2.64e-13</td>
          <td>1.14e-13</td>
          <td>3.59e-12</td>
        </tr>
        <tr>
          <th>56</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>5.04e-13</td>
          <td>-3.52e-13</td>
          <td>-2.37e-12</td>
          <td>-1.82e-12</td>
          <td>-7.96e-13</td>
          <td>3.98e-13</td>
          <td>-9.10e-13</td>
          <td>-3.23e-13</td>
          <td>1.14e-13</td>
          <td>2.55e-12</td>
        </tr>
        <tr>
          <th>880</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>1.79e-12</td>
          <td>6.09e-13</td>
          <td>-2.51e-12</td>
          <td>1.25e-12</td>
          <td>-2.27e-12</td>
          <td>9.09e-13</td>
          <td>-9.09e-13</td>
          <td>-4.55e-13</td>
          <td>3.41e-13</td>
          <td>3.85e-12</td>
        </tr>
        <tr>
          <th>952</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>1.30e-12</td>
          <td>1.93e-13</td>
          <td>-2.63e-12</td>
          <td>-9.09e-13</td>
          <td>-6.82e-13</td>
          <td>4.55e-13</td>
          <td>1.82e-12</td>
          <td>4.55e-13</td>
          <td>2.27e-13</td>
          <td>3.51e-12</td>
        </tr>
        <tr>
          <th rowspan="9" valign="top">3685</th>
          <th>0</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>1.18e-12</td>
          <td>-3.91e-15</td>
          <td>-9.47e-13</td>
          <td>-1.14e-13</td>
          <td>1.14e-13</td>
          <td>2.27e-13</td>
          <td>4.97e-13</td>
          <td>-9.09e-13</td>
          <td>-1.71e-13</td>
          <td>1.84e-12</td>
        </tr>
        <tr>
          <th>45</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>8.96e-13</td>
          <td>-2.21e-13</td>
          <td>-1.07e-12</td>
          <td>2.27e-13</td>
          <td>-5.68e-14</td>
          <td>-5.68e-13</td>
          <td>-2.84e-13</td>
          <td>-6.52e-13</td>
          <td>5.54e-13</td>
          <td>1.71e-12</td>
        </tr>
        <tr>
          <th>41</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>8.01e-13</td>
          <td>-3.97e-13</td>
          <td>-8.59e-13</td>
          <td>-2.84e-13</td>
          <td>0.00e+00</td>
          <td>-1.71e-13</td>
          <td>-4.62e-13</td>
          <td>-6.12e-13</td>
          <td>3.56e-13</td>
          <td>1.48e-12</td>
        </tr>
        <tr>
          <th>1021</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>9.04e-13</td>
          <td>4.58e-13</td>
          <td>-3.01e-12</td>
          <td>-1.48e-12</td>
          <td>5.68e-13</td>
          <td>-7.39e-13</td>
          <td>4.50e-13</td>
          <td>-1.82e-12</td>
          <td>9.06e-14</td>
          <td>3.71e-12</td>
        </tr>
        <tr>
          <th>1015</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>1.19e-12</td>
          <td>-9.78e-14</td>
          <td>-6.96e-13</td>
          <td>3.41e-13</td>
          <td>-5.68e-14</td>
          <td>1.14e-13</td>
          <td>-2.27e-13</td>
          <td>-9.09e-13</td>
          <td>1.13e-13</td>
          <td>1.67e-12</td>
        </tr>
        <tr>
          <th>50</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>4.88e-13</td>
          <td>-1.28e-13</td>
          <td>-1.21e-12</td>
          <td>-3.98e-13</td>
          <td>1.14e-13</td>
          <td>-5.68e-13</td>
          <td>-5.68e-13</td>
          <td>-4.90e-13</td>
          <td>-2.27e-13</td>
          <td>1.49e-12</td>
        </tr>
        <tr>
          <th>46</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>6.22e-13</td>
          <td>2.26e-14</td>
          <td>-5.87e-13</td>
          <td>-2.84e-13</td>
          <td>3.41e-13</td>
          <td>0.00e+00</td>
          <td>4.54e-13</td>
          <td>9.59e-14</td>
          <td>-2.27e-13</td>
          <td>1.05e-12</td>
        </tr>
        <tr>
          <th>1031</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>2.32e-12</td>
          <td>-6.90e-13</td>
          <td>-1.63e-12</td>
          <td>-2.27e-13</td>
          <td>-6.82e-13</td>
          <td>9.09e-13</td>
          <td>4.55e-13</td>
          <td>-1.82e-12</td>
          <td>-2.27e-13</td>
          <td>3.57e-12</td>
        </tr>
        <tr>
          <th>1037</th>
          <th>9.51</th>
          <th>59.74</th>
          <th>3568.63</th>
          <th>2</th>
          <td>2.08e-12</td>
          <td>-2.71e-13</td>
          <td>-1.70e-12</td>
          <td>4.55e-13</td>
          <td>-3.41e-13</td>
          <td>0.00e+00</td>
          <td>-3.98e-13</td>
          <td>-1.82e-12</td>
          <td>-1.14e-13</td>
          <td>3.30e-12</td>
        </tr>
      </tbody>
    </table>
    </div>



Alternatively I can do that by first using reset_index to move all the
index columns into data, and then using set_index to define the order of
columns I want as my index:

.. code:: ipython3

    cht.reset_index().set_index(['ElementID','NodeID','Mode','Freq']).sort_index()




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
          <th>Item</th>
          <th>Eigenvalue</th>
          <th>Radians</th>
          <th>omax</th>
          <th>omid</th>
          <th>omin</th>
          <th>oxx</th>
          <th>oyy</th>
          <th>ozz</th>
          <th>txy</th>
          <th>txz</th>
          <th>tyz</th>
          <th>von_mises</th>
        </tr>
        <tr>
          <th>ElementID</th>
          <th>NodeID</th>
          <th>Mode</th>
          <th>Freq</th>
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
          <th rowspan="18" valign="top">3684</th>
          <th rowspan="2" valign="top">0</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>1.48e-12</td>
          <td>1.02e-12</td>
          <td>-3.38e-12</td>
          <td>1.22e-12</td>
          <td>-3.35e-12</td>
          <td>1.25e-12</td>
          <td>-3.27e-13</td>
          <td>-2.27e-13</td>
          <td>2.84e-14</td>
          <td>4.64e-12</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>2.22e-12</td>
          <td>-1.06e-13</td>
          <td>-2.23e-12</td>
          <td>-3.41e-13</td>
          <td>-2.27e-13</td>
          <td>4.55e-13</td>
          <td>1.82e-12</td>
          <td>1.14e-12</td>
          <td>3.98e-13</td>
          <td>3.85e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">51</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>2.32e-12</td>
          <td>1.49e-13</td>
          <td>-1.41e-12</td>
          <td>-1.39e-12</td>
          <td>2.32e-12</td>
          <td>1.35e-13</td>
          <td>-1.30e-13</td>
          <td>-1.46e-13</td>
          <td>7.51e-15</td>
          <td>3.25e-12</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>7.48e-13</td>
          <td>-3.94e-13</td>
          <td>-1.46e-12</td>
          <td>-3.41e-13</td>
          <td>5.68e-14</td>
          <td>-8.24e-13</td>
          <td>1.57e-13</td>
          <td>3.85e-13</td>
          <td>9.24e-13</td>
          <td>1.91e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">55</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>4.81e-12</td>
          <td>1.92e-13</td>
          <td>-3.57e-13</td>
          <td>4.53e-12</td>
          <td>-2.42e-13</td>
          <td>3.55e-13</td>
          <td>-6.54e-13</td>
          <td>-9.30e-13</td>
          <td>2.20e-15</td>
          <td>4.92e-12</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>3.64e-12</td>
          <td>9.55e-13</td>
          <td>-1.33e-12</td>
          <td>2.33e-12</td>
          <td>-5.68e-14</td>
          <td>9.95e-13</td>
          <td>1.82e-12</td>
          <td>4.26e-13</td>
          <td>1.12e-12</td>
          <td>4.31e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">56</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>1.66e-12</td>
          <td>-1.65e-12</td>
          <td>-5.92e-12</td>
          <td>-5.87e-12</td>
          <td>1.65e-12</td>
          <td>-1.69e-12</td>
          <td>-3.07e-13</td>
          <td>-4.22e-13</td>
          <td>8.53e-14</td>
          <td>6.59e-12</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>5.04e-13</td>
          <td>-3.52e-13</td>
          <td>-2.37e-12</td>
          <td>-1.82e-12</td>
          <td>-7.96e-13</td>
          <td>3.98e-13</td>
          <td>-9.10e-13</td>
          <td>-3.23e-13</td>
          <td>1.14e-13</td>
          <td>2.55e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">60</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>2.88e-12</td>
          <td>1.71e-12</td>
          <td>-4.27e-12</td>
          <td>2.63e-12</td>
          <td>1.90e-12</td>
          <td>-4.21e-12</td>
          <td>-4.26e-13</td>
          <td>-6.34e-13</td>
          <td>8.53e-14</td>
          <td>6.64e-12</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>1.02e-12</td>
          <td>-2.20e-12</td>
          <td>-2.85e-12</td>
          <td>-1.48e-12</td>
          <td>-3.41e-13</td>
          <td>-2.22e-12</td>
          <td>1.82e-12</td>
          <td>2.64e-13</td>
          <td>1.14e-13</td>
          <td>3.59e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">758</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>5.79e-12</td>
          <td>4.11e-12</td>
          <td>7.57e-14</td>
          <td>5.68e-12</td>
          <td>1.14e-13</td>
          <td>4.18e-12</td>
          <td>-4.55e-13</td>
          <td>-3.41e-13</td>
          <td>-3.91e-14</td>
          <td>5.09e-12</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>1.03e-12</td>
          <td>-1.04e-12</td>
          <td>-3.57e-12</td>
          <td>-5.68e-13</td>
          <td>-2.33e-12</td>
          <td>-6.82e-13</td>
          <td>1.82e-12</td>
          <td>4.55e-13</td>
          <td>9.08e-13</td>
          <td>3.99e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">778</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>-1.38e-12</td>
          <td>-3.27e-12</td>
          <td>-6.12e-12</td>
          <td>-6.08e-12</td>
          <td>-1.38e-12</td>
          <td>-3.31e-12</td>
          <td>-5.81e-14</td>
          <td>-3.41e-13</td>
          <td>-1.97e-14</td>
          <td>4.14e-12</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>1.24e-12</td>
          <td>-5.76e-13</td>
          <td>-2.37e-12</td>
          <td>-2.27e-12</td>
          <td>7.96e-13</td>
          <td>-2.27e-13</td>
          <td>1.25e-13</td>
          <td>4.55e-13</td>
          <td>7.73e-13</td>
          <td>3.13e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">880</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>2.63e-12</td>
          <td>2.11e-12</td>
          <td>-4.74e-12</td>
          <td>-4.72e-12</td>
          <td>2.10e-12</td>
          <td>2.61e-12</td>
          <td>-3.41e-13</td>
          <td>0.00e+00</td>
          <td>8.53e-14</td>
          <td>7.12e-12</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>1.79e-12</td>
          <td>6.09e-13</td>
          <td>-2.51e-12</td>
          <td>1.25e-12</td>
          <td>-2.27e-12</td>
          <td>9.09e-13</td>
          <td>-9.09e-13</td>
          <td>-4.55e-13</td>
          <td>3.41e-13</td>
          <td>3.85e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">952</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>-8.35e-15</td>
          <td>-1.79e-12</td>
          <td>-3.00e-12</td>
          <td>-1.73e-12</td>
          <td>-1.14e-13</td>
          <td>-2.96e-12</td>
          <td>-4.26e-13</td>
          <td>2.27e-13</td>
          <td>5.68e-14</td>
          <td>2.61e-12</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>1.30e-12</td>
          <td>1.93e-13</td>
          <td>-2.63e-12</td>
          <td>-9.09e-13</td>
          <td>-6.82e-13</td>
          <td>4.55e-13</td>
          <td>1.82e-12</td>
          <td>4.55e-13</td>
          <td>2.27e-13</td>
          <td>3.51e-12</td>
        </tr>
        <tr>
          <th rowspan="18" valign="top">3685</th>
          <th rowspan="2" valign="top">0</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>5.11e-13</td>
          <td>1.44e-13</td>
          <td>-5.41e-13</td>
          <td>4.97e-13</td>
          <td>1.56e-13</td>
          <td>-5.40e-13</td>
          <td>6.75e-14</td>
          <td>-2.84e-14</td>
          <td>1.42e-14</td>
          <td>9.25e-13</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>1.18e-12</td>
          <td>-3.91e-15</td>
          <td>-9.47e-13</td>
          <td>-1.14e-13</td>
          <td>1.14e-13</td>
          <td>2.27e-13</td>
          <td>4.97e-13</td>
          <td>-9.09e-13</td>
          <td>-1.71e-13</td>
          <td>1.84e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">41</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>7.67e-13</td>
          <td>-1.49e-13</td>
          <td>-4.37e-13</td>
          <td>-1.49e-13</td>
          <td>7.21e-13</td>
          <td>-3.91e-13</td>
          <td>1.92e-13</td>
          <td>1.14e-13</td>
          <td>5.33e-14</td>
          <td>1.09e-12</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>8.01e-13</td>
          <td>-3.97e-13</td>
          <td>-8.59e-13</td>
          <td>-2.84e-13</td>
          <td>0.00e+00</td>
          <td>-1.71e-13</td>
          <td>-4.62e-13</td>
          <td>-6.12e-13</td>
          <td>3.56e-13</td>
          <td>1.48e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">45</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>8.09e-13</td>
          <td>-5.47e-13</td>
          <td>-8.09e-13</td>
          <td>6.79e-13</td>
          <td>-7.89e-13</td>
          <td>-4.37e-13</td>
          <td>6.39e-14</td>
          <td>-4.02e-13</td>
          <td>5.15e-14</td>
          <td>1.50e-12</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>8.96e-13</td>
          <td>-2.21e-13</td>
          <td>-1.07e-12</td>
          <td>2.27e-13</td>
          <td>-5.68e-14</td>
          <td>-5.68e-13</td>
          <td>-2.84e-13</td>
          <td>-6.52e-13</td>
          <td>5.54e-13</td>
          <td>1.71e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">46</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>2.26e-12</td>
          <td>7.91e-13</td>
          <td>7.55e-14</td>
          <td>2.25e-12</td>
          <td>7.96e-13</td>
          <td>7.82e-14</td>
          <td>8.37e-14</td>
          <td>-7.01e-14</td>
          <td>1.42e-14</td>
          <td>1.93e-12</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>6.22e-13</td>
          <td>2.26e-14</td>
          <td>-5.87e-13</td>
          <td>-2.84e-13</td>
          <td>3.41e-13</td>
          <td>0.00e+00</td>
          <td>4.54e-13</td>
          <td>9.59e-14</td>
          <td>-2.27e-13</td>
          <td>1.05e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">50</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>9.64e-13</td>
          <td>7.11e-13</td>
          <td>-6.57e-14</td>
          <td>-3.91e-14</td>
          <td>7.11e-13</td>
          <td>9.38e-13</td>
          <td>-4.26e-14</td>
          <td>-1.57e-13</td>
          <td>1.42e-14</td>
          <td>9.29e-13</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>4.88e-13</td>
          <td>-1.28e-13</td>
          <td>-1.21e-12</td>
          <td>-3.98e-13</td>
          <td>1.14e-13</td>
          <td>-5.68e-13</td>
          <td>-5.68e-13</td>
          <td>-4.90e-13</td>
          <td>-2.27e-13</td>
          <td>1.49e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">1015</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>8.67e-13</td>
          <td>1.56e-13</td>
          <td>-6.65e-13</td>
          <td>2.98e-13</td>
          <td>1.53e-13</td>
          <td>-9.24e-14</td>
          <td>2.84e-14</td>
          <td>-7.39e-13</td>
          <td>4.51e-14</td>
          <td>1.33e-12</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>1.19e-12</td>
          <td>-9.78e-14</td>
          <td>-6.96e-13</td>
          <td>3.41e-13</td>
          <td>-5.68e-14</td>
          <td>1.14e-13</td>
          <td>-2.27e-13</td>
          <td>-9.09e-13</td>
          <td>1.13e-13</td>
          <td>1.67e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">1021</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>7.75e-13</td>
          <td>1.86e-13</td>
          <td>-5.38e-13</td>
          <td>1.56e-13</td>
          <td>1.60e-13</td>
          <td>1.07e-13</td>
          <td>1.94e-13</td>
          <td>-6.25e-13</td>
          <td>5.47e-14</td>
          <td>1.14e-12</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>9.04e-13</td>
          <td>4.58e-13</td>
          <td>-3.01e-12</td>
          <td>-1.48e-12</td>
          <td>5.68e-13</td>
          <td>-7.39e-13</td>
          <td>4.50e-13</td>
          <td>-1.82e-12</td>
          <td>9.06e-14</td>
          <td>3.71e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">1031</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>-8.72e-13</td>
          <td>-1.24e-12</td>
          <td>-1.87e-12</td>
          <td>-1.68e-12</td>
          <td>-8.81e-13</td>
          <td>-1.42e-12</td>
          <td>8.53e-14</td>
          <td>-2.84e-13</td>
          <td>4.26e-14</td>
          <td>8.75e-13</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>2.32e-12</td>
          <td>-6.90e-13</td>
          <td>-1.63e-12</td>
          <td>-2.27e-13</td>
          <td>-6.82e-13</td>
          <td>9.09e-13</td>
          <td>4.55e-13</td>
          <td>-1.82e-12</td>
          <td>-2.27e-13</td>
          <td>3.57e-12</td>
        </tr>
        <tr>
          <th rowspan="2" valign="top">1037</th>
          <th>1</th>
          <th>8.36</th>
          <td>2758.15</td>
          <td>52.52</td>
          <td>-7.47e-15</td>
          <td>-1.15e-13</td>
          <td>-1.00e-12</td>
          <td>-9.95e-14</td>
          <td>-1.14e-13</td>
          <td>-9.09e-13</td>
          <td>-7.11e-15</td>
          <td>-2.84e-13</td>
          <td>2.84e-14</td>
          <td>9.43e-13</td>
        </tr>
        <tr>
          <th>2</th>
          <th>9.51</th>
          <td>3568.63</td>
          <td>59.74</td>
          <td>2.08e-12</td>
          <td>-2.71e-13</td>
          <td>-1.70e-12</td>
          <td>4.55e-13</td>
          <td>-3.41e-13</td>
          <td>0.00e+00</td>
          <td>-3.98e-13</td>
          <td>-1.82e-12</td>
          <td>-1.14e-13</td>
          <td>3.30e-12</td>
        </tr>
      </tbody>
    </table>
    </div>



