import scipy
if scipy.__version__ >= '1.14': # TODO: pick a better version
    # tested on 1.14.1, 1.16.3
    from scipy.sparse import csc_matrix, dok_matrix
else:
    from scipy.sparse._dok import dok_matrix
    from scipy.sparse.csc import csc_matrix
