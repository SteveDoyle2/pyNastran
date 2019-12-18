"""
CBUSH          6       3       7       2                               0
XYPUNCH,ELFORCE/ 6 (2)
XYPUNCH,ELFORCE/ 6 (3)
XYPUNCH,ELFORCE/ 6 (4)
XYPUNCH,ELFORCE/ 6 (5)
XYPUNCH,ELFORCE/ 6 (6)
XYPUNCH,ELFORCE/ 6 (7)

CBAR           3       1       5       6      0.      0.      1.
XYPUNCH,ELFORCE/ 3 (2)
XYPUNCH,ELFORCE/ 3 (3)
XYPUNCH,ELFORCE/ 3 (4)
XYPUNCH,ELFORCE/ 3 (5)
XYPUNCH,ELFORCE/ 3 (6)
XYPUNCH,ELFORCE/ 3 (7)

CBAR           5       1       2       3      0.      0.      1.
XYPUNCH,ELFORCE/ 5 (2)
XYPUNCH,ELFORCE/ 5 (3)
XYPUNCH,ELFORCE/ 5 (4)
XYPUNCH,ELFORCE/ 5 (5)
XYPUNCH,ELFORCE/ 5 (6)
XYPUNCH,ELFORCE/ 5 (7)

"""
from pyNastran.op2.op2 import OP2
from pyNastran.f06.f06_formatting import write_float_12e

op2 = OP2()
op2.read_op2('good_sine.op2')
acc = op2.accelerations[103]
nnids = [1, 9]  # 1, 5, 9, 1000

cases = [
    (acc, nnids, [1, 2], ['real', 'imag', 'mag', 'phase']),
]
#$SUBCASE     101                                                               1

iline = 1
pch_filename = 'test.pch'
with open(pch_filename, 'w') as pch_file:
    for case_line in cases:
        (obj, nids, comps, output) = case_line
        title = obj.subtitle
        times = obj._times
        for nid in nids:
            for comp in comps:
                for outputi in output:
                    msg = '%-59s %20i\n' % (title, iline)
                    iline += 1
                    nid = [nid]
                    data = obj.extract_xyplot(nid, comp, outputi)
                    data = data.reshape(len(data))
                    for idata, (time, datai) in enumerate(zip(times, data)):
                        msg += '%10i                  %12s        %12s%20i\n' % (
                            idata + 1, write_float_12e(time), write_float_12e(datai), iline)
                        iline += 1
                    pch_file.write(msg)
                    # print(msg)
                    # break
                # break
            # break

#acc.extract_xyplot(nids, 1, 'real')
#acc.extract_xyplot(nids, 1, 'imag')
#acc.extract_xyplot(nids, 1, 'mag')
#acc.extract_xyplot(nids, 1, 'phase')

