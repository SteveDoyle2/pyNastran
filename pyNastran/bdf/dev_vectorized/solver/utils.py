from six.moves import zip
from numpy import dot, ndarray, zeros

def partition_sparse(Is, Js, Vs):
    I2 = []
    J2 = []
    V2 = []
    for (i, j, v) in (Is, Js, Vs):
        if abs(v) >= 1e-8:
            I2.append(i)
            J2.append(j)
            V2.append(v)
    return(I2, J2, V2)


def getDOF_Set(nAll, dofs):
    dofsAll = set([i for i in xrange(nAll)])
    dofs = list(dofsAll.difference(set(dofs)))
    return dofs

def remove_dofs(dofsAll, dofs_remove):
    dofs = list(dofsAll.difference(set(dofs_remove)))
    return dofs


def partition_dense_symmetric(A, dofs_in):
    nAll = A.shape[0]
    dofs = getDOF_Set(nAll, dofs_in)
    dofs.sort()
    n = len(dofs)
    A2 = zeros((n, n), 'float64')
    for (i, dofI) in enumerate(dofs):
        for (j, dofJ) in enumerate(dofs):
            v = A[dofI, dofJ]
            if abs(v) >= 1e-8:
                A2[i, j] = v
    return (A2, dofs)


def partition_dense_vector(F, dofs_in):
    nAll = F.shape[0]
    #print("dofs = %s" % dofs)
    dofs = getDOF_Set(nAll, dofs_in)
    dofs.sort()
    #print("dofs = %s" % dofs)
    n = len(dofs)
    F2 = zeros(n, 'float64')
    for (i, dofI) in enumerate(dofs):
        v = F[dofI]
        if abs(v) >= 1e-8:
            F2[i] = v
    return (F2, dofs)


def partition_sparse_vector(F, dofs):
    dofs.sort()
    #n = len(dofs)
    F2i = []
    F2v = []
    for (i, dofI) in enumerate(dofs):
        v = F2v[dofI]
        if abs(v) >= 1e-8:
            F2i.append(i)
            F2v.append(v)
    return(F2i, F2v)


def departition_dense_vector(n, IsVs):
    V = zeros(n)
    for IV in IsVs:
        (Is, Vs) = IV
        for (i, v) in zip(Is, Vs):
            V[i] = v
    return(V)


def reverse_dict(A):
    B = {}
    for (key, value) in A.iteritems():
        B[value] = key
    return B


def triple(A, B):
    """
    A.T @ B @ A

    A   [n x m]
    A.T [m x n]
    T [m x m] = A.T [m x n] @ B [n x n] @ A [n x m]
    """
    assert isinstance(A, ndarray), type(A)
    assert isinstance(A, ndarray), type(B)
    return dot(A.T, dot(B, A))

