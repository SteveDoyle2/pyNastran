r"""
defines:
 - rainflow_from_arrays(A, features, casenames=None)
 - rainflow_from_csv(input_csv, casenames, features, write_csvs=True, delimiter=',',
                     xmax=None, legend_alpha=1.0)

Fatigue Overview
================
Fatigue cycling for a single feature (e.g., element or part) is:

    ^ stress for bar under axial load; 100k cycles
    |
  4 |         /\
  3 |  /\    /  \
    | /  \  /    \
  1 |/    \/      \
  0 +--------------\------->
    |               \  /
 -2 |                \/

Pretty basic so far, but we also may have:

   ^ stress for bar under torsional load; 200k cycles
   |
   |
 3 |  /\
 2 | /  \  /\
   |/    \/  \
 0 +-------------->

So for the same bar (a "feature"), we have two loading spectra that do
not happen at the same time.  The first "case" is "axial", and the
second case is "torison".

We need to reduce this into min/max stresses, which can be used to
compute damage for each case.  We sum this up using Miner's rule and
can compute a fatigue margin.


Fatigue Implementation
======================
We only handle the stress cycling aspect here, so not the damage or
the fatigue margin calculation.

We'll define the first case per the above graphs:

Bar Axial
=========
    0
    3
    1
    4
   -2
    0 (implicit if not explicitly defined)

Bar Torsion
===========
   0
   3
   0
   2
   0

So, we can vstack these together and define cases:

 iRow     Bar Stress
 =====    ==========
   0          0
   1          3
   2          1
   3          4
   4          2
   5          0
   6          0  <---- start of torsion
   7          3
   8          0
   9          2
  10          0

A is either (10, ) or (10, 1)

# TODO: check for off by 1 errors
cases = [
    ('axial', 0, 5),
    ('torison', 6, 10),
]

features = {
    # the key indicates the column of A
    0 : 'bar',
}

If we wanted a second feature, we'd just make A (10, 2) and add a new
entry to "features".


We choose to make the values implicit, so we can define impulsive cases as:

Bar Impulse
===========
    0
   10
"""
# encoding: utf-8
from __future__ import print_function, unicode_literals
from codecs import open

from six.moves import range
from six import iteritems, iterkeys
import numpy as np


def rainflow(icase, stress_in):
    """
    Does rainflow counting based on stress (not nominal stress).
    Works with a non-minimum first value.

    Parameters
    ----------
    icase : int
        unused; simply a flag for debugging
    stress_in : (N, ) ndarray
        the stresses to cycle count
    """
    stress = reorganize_load(icase, stress_in)
    stress, cycles = ASTM_E1049_rainflow(stress)
    max_stress = []
    min_stress = []

    for i in range(len(cycles)):
        max_stress.append(max(cycles[i][0], cycles[i][1]))
        min_stress.append(min(cycles[i][0], cycles[i][1]))
    return (max_stress, min_stress)


def reorganize_load(icase, stress_in):
    """
    Reorganize Order History - ASTM E1049 5.4.5.2 (1)

    We want to start at the minimum value and then loop around
    until we get back to where we started

    Parameters
    ----------
    icase : int
        unused; simply a flag for debugging
    stress_in : (N, ) ndarray
        the stresses to cycle count
    """
    data = np.asarray(stress_in)
    imin = np.where(data == data.min())[0]
    imin0 = imin[0]

    # .. todo:: the double type conversion (list & array) is less than ideal
    x = list(stress_in)
    if imin0 != 0:
        x = x[imin0:] + x[:imin0]

    # remove repeats
    i = 0
    y = [x[0]]
    while i < len(x) - 1:
        if x[i+1] != x[i]:
            y.append(x[i+1])
        i += 1

    # handles single-cycle impulse loading
    # and 0 to 359 degrees with an assumed 360 degrees loop
    if y[0] != y[-1]:
        #print('  y[0]=%s y[-1]=%s; adding y[0]' % (y[0], y[-1]))
        y.append(y[0])
    return y

def ASTM_E1049_rainflow(stress_in):
    """
    Does rainflow counting based on stress (not nominal stress).
    Works with a non-minimum first value.

    From ASTM E1049 5.4.5.2 (1)

    The ASTM spec is very Fortran heavy, so rather than Python-ifying it,
    we just copy the spec.

    .. note:: Assumes the minimum value is stress_in[0] and
                      the final value is a repeat of stress_in[0].
    """
    x = 0.0
    y = 0.0

    stress = {}
    cycles = {}
    i = 0
    icyc = 0
    while stress_in:
        temp_stress = stress_in.pop(0)
        stress[i] = temp_stress
        if i > 1:
            if((stress[i] - stress[i-1] >= 0 and stress[i-1] - stress[i-2] >= 0) or
               (stress[i] - stress[i-1] <= 0 and stress[i-1] - stress[i-2] < 0)):
                del stress[i]
                i -= 1
                stress[i] = temp_stress
        i += 1

    i = -1
    istop = 0
    goto = 2
    while istop == 0:
        if goto == 2:
            i += 1
            if i >= len(stress):
                istop = 1  # rainflow counting finished
            goto = 3
        elif goto == 3:
            if i < 2:
                goto = 2
            else:
                y = abs(stress[i-1] - stress[i-2])
                x = abs(stress[i]   - stress[i-1])
                goto = 4
        elif goto == 4:
            if x < y:
                goto = 2
            else:
                goto = 5
        elif goto == 5:
            cycles[icyc] = [min(stress[i-1], stress[i-2]),
                            max(stress[i-1], stress[i-2]), y]
            icyc += 1
            jend = len(stress) - 1
            for j in range(i, jend+1):
                stress[j-2] = stress[j]
            del stress[jend]
            del stress[jend - 1]
            i -= 2
            goto = 3
    return (stress, cycles)


def _fix_input(A, casenames=None, features=None, nfeatures_max=None, filter_stresses=True):
    """
    Prepares the data for rainflow counting

    Parameters
    ----------
    A : (n, ) or (n, m) float ndarray
       the matrix of stresses/forces; (n, 0) is the data of interest
    casenames : [(casename, min_index, max_index), ...]; default=None
        Let's you analyze a subsection of steps
        None : consider all load steps as one case
    features : dict[count] = name; default=None
        dict:
            count : int
                the column index of the A matrix
            name : str
                unused
        None : auto-generates features with count=icolumn
    nfeatures_max : int; default=None
        limit the number of cases to the max N values
        None : all

    Returns
    -------
    A : (n, m) float ndarray
        the matrix of stresses/forces
    casenames : [(casename, min_index, max_index), ...]
        Let's you analyze a subsection of steps
    features : dict[count] = name
        count : int
            the column index of the A matrix
        name : str
            unused
    features_to_analyze : List[int]
        the indicies from A to analyze
    filter_stresses : bool; default=True
        flag to filter stresses that are equal and don't create cycles
    """
    if len(A.shape) == 1:
        A = A.reshape(len(A), 1)

    nfeatures = A.shape[1]
    if casenames is None:
        casename = 'case1'
        min_index = 0
        max_index = A.shape[0]
        casenames = [(casename, min_index, max_index)]

    if features is None:
        features = {}
        for i in range(nfeatures):
            features[i] = 'feature%i' % i

    if nfeatures_max is not None:
        maxs = np.abs(A).max(axis=0)

        smaxs = np.sort(maxs)
        nmax = len(maxs)

        # the highest is at the end of the array
        icutoff = nmax - nfeatures_max - 1
        assert icutoff >= 0, icutoff
        try:
            max_cutoff = smaxs[icutoff]
            get_features = True
        except IndexError:
            features_to_analyze = iterkeys(features)
            get_features = False
        if get_features:
            features_to_analyze = np.where(maxs > max_cutoff)[0]
    else:
        features_to_analyze = np.array(list(iterkeys(features)))

    if filter_stresses:
        delta = A[:, features_to_analyze].max(axis=0) - A[:, features_to_analyze].min(axis=0)
        features_to_analyze = np.where(delta > 0)[0]

    return A, casenames, features, features_to_analyze

def rainflow_from_arrays(A, features=None, casenames=None, nfeatures_max=None, filter_stresses=True):
    """
    calculates rainflow from an array

    A : (nfeatures, nload_steps) float ndarray
        the matrix of forces/stresses
    features : dict[count] = name; default=None
        dict:
            count : int
                the column index of the A matrix
            name : str
                unused
        None : auto-generates features with count=icolumn
    casenames : [(casename, min_index, max_index), ...]; default=None
        Let's you analyze a subsection of steps
        None : consider all load steps as one case
    nfeatures_max : int; default=None
        limit the number of cases to the max N values
        None : all
    filter_stresses : bool; default=True
        flag to filter stresses that are equal and don't create cycles

    Returns
    -------
    min_max_stresses_dict : dict[count] : feature_stresses
        count : int
            the column index of the A matrix
        feature_stresses : List[min_max_stress1, min_max_stress2, ...]
            the stresses for each feature
            min_max_stressi : List[min_stress, max_stress]
                the stresses for each case
                n = len(min_stress) == len(max_stress), but varies for each min_stress

                min_stress : List[stress_mini, ...]; length=n
                    the min stresses
                max_stress : List[stress_maxi, ...]; length=n
                    the max stresses
    """
    A, casenames, features, features_to_analyze = _fix_input(
        A, casenames, features, nfeatures_max, filter_stresses)

    icase = 0
    min_max_stresses_dict = {}
    for ifeature in sorted(features_to_analyze):
        min_max_stresses = []
        for case_name, min_index, max_index in casenames:
            stress_case = A[min_index:max_index, ifeature]
            min_max_stress = rainflow(icase, stress_case)
            #min_stress, max_stress = min_max_stress
            min_max_stresses.append(min_max_stress)
            icase += 1
        min_max_stresses_dict[ifeature] = min_max_stresses
    return min_max_stresses_dict

def write_min_max_stress_dict(stress_filename, min_max_stresses_dict, eids):
    """
    Writes the data in the following form:

    # nfeatures = 8484
    # ncases = 1
    # icase = 0
    # Eid, stress_min_cycle1, stress_max_cycle1, stress_min_cycle2, stress_max_cycle2, ...
    0, 6.29373e+06, 6.28326e+06, 6.34438e+06, 6.29244e+06
    1, 570227.0, 569266.0,
    2, 577549.0, 576540.0, 582015.0, 577219.0,
    10, 577798.0, 576761.0, 582172.0, 577387.0, 582789.0, 582463.0,
    20,
    21, 570227.0, 569266.0,
    """
    nfeatures = len(min_max_stresses_dict)
    key0 = next(iterkeys(min_max_stresses_dict))
    feature0 = min_max_stresses_dict[key0]
    ncases = len(feature0)
    msg =  '# nfeatures = %i\n' % nfeatures
    msg += '# ncases = %i\n' % ncases

    for icase in range(ncases):
        msg += '# icase = %i\n' % icase
        msg += '# Eid, stress_min_cycle1, stress_max_cycle1, stress_min_cycle2, stress_max_cycle2, ...\n'
        for ieid, cases in sorted(iteritems(min_max_stresses_dict)):
            eid = eids[ieid]
            line = '%s, ' % eid
            case = cases[icase]
            min_stresses, max_stresses = case
            #print_count = 0
            for min_stress, max_stress in zip(min_stresses, max_stresses):
                if not np.allclose(min_stress, max_stress):
                    line += '%s, %s, ' % (min_stress, max_stress)
                    #print_count += 1
            msg += line + '\n'
    with open(stress_filename, 'w') as stress_file:
        stress_file.write(msg)


def rainflow_from_csv(input_csv, casenames, features,
                      write_csvs=True, delimiter=',',
                      xmax=None, legend_alpha=1.0, debug=True):
    """
    Rainflow counts from csv files.

    This supports multiple features as separate columns.

    Parameters
    ----------
    input_csv : str
        a csv file as described below
    casenames : str
        tuple of cases; see below
    features : dict
       key : int; column id to parse
       value : str; name

    xmax : float
        the max value for the x (cycle) axis; helps to change the legend
    delimiter : str; default=','
        the delimiter for the output file (doesn't apply to input)
    legend_alpha : float; default=1.0
        the transparency
        1=solid
        0=transparent
    debug : bool; default=True
        prints debug messages

    Returns
    -------
    files : list[str]
        filenames are of the form icase_icase_name.csv

    Input_csv
    ---------
      # name1_stress, name2_stress, ...
      0.00, 0.0 # case 0 - line 1
      20.0, 1.0 # case 1
      50.0, 2.0 # case 2
      etc.
      0.00, 0.0 # case 0
    casenames = (
       # (casename, irow_start, irow_stop)
       ('normal',  0,  62),
       ('impulse', 63, 65),
       etc.
    )
    features = {  # the indicies are column numbers
        0 : 'fillet',
        1 : 'groove',
    }
    features = ['feature1', 'feature2']
    so we get:
       feature0_normal_fillet.csv
       feature0_impulse_fillet.csv
       feature1_normal_groove.csv
       feature1_impulse_groove.csv

    We'll also get corresponding png files. of the form:
       fillet.png
       groove.png

    that show our cycling.
    """
    import matplotlib.pyplot as plt
    A = np.loadtxt(input_csv, delimiter=',', skiprows=1)
    A, casenames, features, features_to_analyze = _fix_input(
        A, casenames, features, nfeatures_max=None)

    icase = 0
    for ifeature, feature_name in sorted(iteritems(features)):
        plt.figure(ifeature)
        legend = []
        for case_name, min_index, max_index in casenames:
            csv_out = 'feature%i_%s_%s.csv'  % (ifeature, case_name, feature_name)
            if debug: # pragma: no cover
                print(csv_out)

            stress_case = A[min_index:max_index, ifeature]
            min_stress, max_stress = rainflow(icase, stress_case)
            if len(min_stress) == 0:
                min_stress = [A[min_index, ifeature]]
                max_stress = [A[max_index - 1, ifeature]]

            B = np.vstack([min_stress, max_stress]).T
            header = 'max stress%smin_stress\n' % delimiter
            np.savetxt(csv_out, B, delimiter=delimiter, header=header)
            plt.plot(range(min_index, max_index), stress_case)
            legend.append(case_name)
            icase += 1
        # add the legend in the middle of the plot
        leg = plt.legend(legend, fancybox=True)
        # set the alpha value of the legend: it will be translucent
        leg.get_frame().set_alpha(legend_alpha)

        plt.title(feature_name)
        if xmax:
            plt.xlim([0, xmax])
        plt.xlabel('Cycle Number')
        plt.ylabel('Stress (ksi)')
        plt.grid(True)
        plt.savefig('%s.png' % feature_name)
    #plt.show()

if __name__ == '__main__':
    main()

