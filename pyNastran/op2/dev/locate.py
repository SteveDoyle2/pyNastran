from numpy import where, array

def findvals(vec, ref_ids):
    print('vec =', vec)
    print('ref_ids =', ref_ids)
    # if ref_ids == 0:
        # return array([-1], dtype='int32')
    i = where(ref_ids == vec)[0]
    print(i)
    return i
