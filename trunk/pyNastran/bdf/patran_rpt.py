def rpt_read():
    lines = open('patran.rpt', 'r').readlines()[14:]

    nheader = 2
    headers = ''
    for n in xrange(nheader):
        header = lines[n].strip('\n\r\t')
        headers += header
    headers = headers.split('-')

    headers2 = []
    for header in headers:
        if header:
            headers2.append(header)
    headers = headers2
    #headers2 = [header if header.strip() for header in headers]
    print headers

    i = 2
    j = 0
    res = 0
    results = {0:[], 1:[]}
    while i < len(lines):
        data = []

        for n in xrange(nheader):
            data += lines[i+n].strip().split()
        if "MSC.Patran" in data:
            res += 1
            i += 11
            continue
        if len(data) == 0:
            break
        results[res].append(data)
        #print "i=%s j=%s data =%s" %(i, j, data)
        i += nheader
        j += 1
        if int(float(data[0])) == 0:
            asdf

    key_map = {}
    for ikey, key in enumerate(headers):
        key_map[key] = ikey
    return headers, results, key_map


def main():
    headers, results, key_map = rpt_read()

    for key, rows in results.iteritems():
        data2 = {}
        iz = key_map['X Location']
        iozz = key_map['Z Component']
        csv_filename = 'rpt_%i.csv' % key
        f = open(csv_filename, 'wb')
        f.write('#x,ozz\n')
        for row in rows:
            z = float(row[iz])
            ozz = float(row[iozz])
            print "z=%s ozz=%s" % (z, ozz)
            f.write('%g,%g\n' % (z, ozz))
        f.close()

        ix = 0
        dx = 0.1
        iname = 1
        csv_simplify(csv_filename, None, ix, iname, tol=dx)
        #break

from numpy import loadtxt, argsort, allclose, vstack, array, abs, where

def csv_simplify(csv_filename, x0, ix, iname, tol=0.05):
    A = loadtxt(csv_filename, delimiter=',')

    f = open(csv_filename + '2', 'wb')

    isort = argsort(A[:, ix])
    X = A[isort, ix]
    aResponse = A[isort, iname]

    f.write('#z,response\n')
    update_flag = 1
    tol = 0.001
    response = None

    for i, x in enumerate(X):
        # reset the max RSS per subcase
        if update_flag:
            x_orig = X[i]
            update_flag = False

        # if the points are the same in the direction of interest, find the max
        print "x=%s x_orig=%s" % (x, x_orig)
        if allclose(x, x_orig, atol=tol):
            if response is None:
                response = aResponse[i]
            else:
                response = vstack((response, aResponse[i]))
            update_flag = False
        else:
            print "response =", response
            if len(response.shape) == 2:
                # max
                response = response[0:].max(axis=0)

                # min
                #response = response[0].max(axis=0)

                # abs(max), abs(min) + sign
                if 0:
                    values2 = array([response.max(),
                                     response.min()])

                    # we figure out the absolute max/min
                    abs_vals = abs(values2)
                    abs_val = abs_vals.max()

                    # find the location of the absolute max value
                    # 1.  we take the first value (the where[0]) to chop the return value
                    #     since there is no else conditional
                    # 2.  we take the first value (the where[0][0]) to only get the max
                    #     value if 2+ values are returned
                    j = where(abs_val == abs_vals)[0][0]

                    # get the raw value from the absoluted value, so:
                    # value = abs(raw_value)
                    response = response[j]

            #print response1, response2
            #f.write( "%g,%g,%g,%g\n" % (x, response1, response2, response1/1e6, response2/1e6))
            f.write( "%g,%g,%g\n" % (x, response, response/1e6))
            #f.write("-------------------\n")
            x_orig = X[i]
            response = aResponse[i]

if __name__ == "__main__":
    main()
