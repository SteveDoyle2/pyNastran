from numpy import ndarray, zeros

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


def get_dof_set(nall, dofs):
    dofs_all = {i for i in range(nall)}
    dofs = list(dofs_all.difference(set(dofs)))
    return dofs

def remove_dofs(dofs_all, dofs_remove):
    dofs = list(dofs_all.difference(set(dofs_remove)))
    return dofs


def partition_dense_symmetric(A, dofs_in):
    nall = A.shape[0]
    dofs = get_dof_set(nall, dofs_in)
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
    nall = F.shape[0]
    #print("partition_dense_vector:  dofs_in = %s" % sorted(dofs_in))
    dofs = get_dof_set(nall, dofs_in)
    dofs.sort()
    #print("partition_dense_vector:  dofs = %s" % dofs)
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
    return V


def reverse_dict(A):
    return {v: k for k, v in A.items()}


def triple(A, B):
    """
    A.T @ B @ A

    A   [n x m]
    A.T [m x n]
    T [m x m] = A.T [m x n] @ B [n x n] @ A [n x m]
    """
    assert isinstance(A, ndarray), type(A)
    assert isinstance(A, ndarray), type(B)
    return A.T @ B @ A

