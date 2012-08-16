import scipy.sparse as ss
from numpy import array, zeros

#row = array(

A = ss.lil_matrix((3, 3), dtype='float32')

rows = []
cols = []
data = []
for i in range(3):
    for j in range(3):
        k = float((i + 1) * (j + 1))
        rows.append(i)
        cols.append(j)
        data.append(k)
        A[i, j] = k

B = ss.coo_matrix((data, (rows, cols)), shape=(3, 3))
print A.todense()
print B.todense()
