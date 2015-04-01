# this needs some serious TLC

import sys
from six.moves import range
from six import iteritems
from numpy import where, array, vstack, savetxt, loadtxt


def rainflow(icase, stress_in):
    """
    Does rainflow counting based on stress (not nominal stress).
    Works with a non-minimum first value.
    """
    stress = reorganize_load(icase, stress_in)
    #print('stress[%i] = %s' % (icase, stress))
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
    """
    data = array(stress_in)
    imin = where(data == data.min())[0]
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
    if y[0] != y[-1]:
        print('  y[0]=%s y[-1]=%s; adding y[0]' % (y[0], y[-1]))
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
            if ((stress[i] - stress[i-1] >= 0 and stress[i-1] - stress[i-2] >= 0) or
                (stress[i] - stress[i-1] <= 0 and stress[i-1] - stress[i-2] <  0)):
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


def rainflow_from_csv(input_csv, casenames, features,
                      write_csvs=True, delimiter=',',
                      xmax=None, legend_alpha=1.0):
    """
    Rainflow counts from csv files.

    This supports multiple features as separate columns.

    :param fname:        a file as described below
    :param casenames:    allows for case splitting
    :param features:     columns to parse

    :param xmax:         the max value for the x (cycle) axis; helps to change the legend
    :param delimiter:    the delimiter for the output file (doesn't apply to input); default=','
    :param legend_alpha: the transparency (1=solid, 0=transparent; default=1.0)
    :returns files:      filenames are of the form icase_icase_name.csv

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
    A = loadtxt(input_csv, delimiter=',', skiprows=1)
    if len(A.shape) == 1:
        A = A.reshape(len(A), 1)
    icase = 0

    for ifeature, feature_name in sorted(iteritems(features)):
        plt.figure(ifeature)
        legend = []
        for case_name, min_index, max_index in casenames:
            csv_out = 'feature%i_%s_%s.csv'  %(ifeature, case_name, feature_name)
            print(csv_out)

            stress_case = A[min_index:max_index, ifeature]
            min_stress, max_stress = rainflow(icase, stress_case)
            if len(min_stress) == 0:
                min_stress = [A[min_index, ifeature]]
                max_stress = [A[max_index - 1, ifeature]]

            B = vstack([min_stress, max_stress]).T
            f = open(csv_out, 'wb')
            f.write('# max stress%smin_stress\n' % delimiter)
            savetxt(f, B, delimiter=delimiter)
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


def rainflow_from_csv2(input_csv, names, write_csvs=True):
    """
    Rainflow counts from csv files

    :param fname:   a file as described below
    :param names:   data to append to the filename
    :retval out_stress: dictionary of name:[N, 2] where
                        2=[min_stress, max_stress] and N=nMinMax
    :returns files: of the form icase_icase_name.csv

    .. note:: This is not the ideal method as it only handles one case,
              and the input has to be hand created with min/max stress
              for each row, so it will be removed at some point.

    Input_csv
    ---------
      # case, min_stress, max_stress
      1, 0.00, 20.0
      1, 20.0, 50.0
      1, 50.0, 0.0
      2, etc.
      3, etc.
      12, etc.

    names = {
        1 : 'A',
        2 : 'B',
        3 : 'C',
        12 : 'D',
    }

    so we get:
       icase_1_A.csv
       icase_2_B.csv
       icase_3_C.csv
       icase_12_D.csv
    """
    A = loadtxt(input_csv, delimiter=',', skiprows=1)
    nrows = A.shape[0]
    ncols = A.shape[1]
    assert ncols == 3, ncols  # icycle, stress_min, smax
    cases = {}
    case = []

    row_id_old = None
    for row in A:
        row_id = int(row[0])
        if row_id == row_id_old:
            case.append(row[1])
            case.append(row[2])
        else:
            if row_id_old is not None:
                if min(case) == max(case):
                    print('  case[%i] = %s'  %(row_id_old, case))
                else:
                    if row_id_old in cases:
                        raise KeyError('row_id=%s exists' % row_id_old)
                    cases[row_id_old] = case

            case = []
            case.append(row[1])
            case.append(row[2])
            row_id_old = int(row[0])
            #print(row_id_old)
    if min(case) == max(case):
        print('  case[%i] = %s'  %(row_id_old, case))
    else:
        if row_id_old in cases:
            raise KeyError('row_id=%s already exists' % row_id_old)
        cases[row_id_old] = case

    out_stress = {}
    for icase, case in sorted(iteritems(cases)):
        name = names[icase]
        max_stress, min_stress = rainflow(icase, case)
        A = vstack([min_stress, max_stress]).T

        if write_csvs:
            print('  max[%i] = %s' % (icase, max_stress))
            print('  min[%i] = %s\n' % (icase, min_stress))

            print(A)
            fname = 'icase_%s_%s.out'% (icase, name)
            f = open(fname, 'wb')
            f.write('# max stress,min_stress\n')
            savetxt(f, A, delimiter=',')
        out_stress[name] = A
    return out_stress

if __name__ == '__main__':
    main()

