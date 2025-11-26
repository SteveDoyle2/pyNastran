from itertools import count
import numpy as np


def get_mesh(mesh, i, j, k):
    print('-------------------')
    print(i, j, k)
    print('mesh', mesh.shape)
    imesh = mesh[i, :, :]
    print('imesh', imesh.shape)
    jmesh = imesh[:, j, :]
    kmesh = jmesh[:, :, k]
    print('kmesh', kmesh.shape)
    return kmesh.flatten()


def main():
    # define the number of nodes in x, y, and z
    nx = 3
    ny = 2
    nz = 2
    # Define the ranges for x, y, and z coordinates
    x_coords = np.linspace(-5, 5, nx)  # 10 points from -5 to 5 for x
    y_coords = np.linspace(-5, 5, ny)  # 10 points from -5 to 5 for y
    z_coords = np.linspace(-5, 5, nz)  # 10 points from -5 to 5 for z
    print(x_coords)
    i0 = np.arange(nx-1)
    i1 = np.arange(ny-1)
    i2 = np.arange(nz-1)
    print('i0,i1,i2', i0, i1, i2)

    ngrid = nx * ny * nz
    mesh = np.arange(1, ngrid+1, dtype='int32').reshape(nx, ny, nz)


    n1 = get_mesh(mesh, i0, i1, i2)
    n2 = get_mesh(mesh, i0+1, i1, i2)
    n3 = get_mesh(mesh, i0+1, i1+1, i2)
    n4 = get_mesh(mesh, i0, i1+1, i2)

    n5 = get_mesh(mesh, i0, i1, i2+1)
    n6 = get_mesh(mesh, i0+1, i1, i2+1)
    n7 = get_mesh(mesh, i0+1, i1+1, i2+1)
    n8 = get_mesh(mesh, i0, i1+1, i2+1)

    # Create the 3D meshgrid
    X, Y, Z = np.meshgrid(x_coords, y_coords, z_coords)
    x = X.ravel()
    y = Y.ravel()
    z = Z.ravel()

    from pyNastran.bdf.mesh_utils.mesh import create_structured_3d
    X, Y, Z, elements = create_structured_3d(
        x_coords, y_coords, z_coords,
        nx, ny, nz)
    x = X.ravel()
    y = Y.ravel()
    z = Z.ravel()

    # X, Y, and Z are now 3D arrays, where:
    # X contains the x-coordinates for each point in the grid
    # Y contains the y-coordinates for each point in the grid
    # Z contains the z-coordinates for each point in the grid

    # You can access the coordinates of a specific point, for example:
    #print(X[0, 0, 0], Y[0, 0, 0], Z[0, 0, 0])
    #print(X.shape, Y.shape, Z.shape)
    pid = 100
    mid = 101
    with open('fem.pch', 'w') as bdf_file:
        bdf_file.write(f'$ pyNastran: punch=True\n')
        bdf_file.write(f'MAT1,{mid},3.0e7,,0.3\n')
        bdf_file.write(f'PSOLID,{pid},{mid}\n')
        for i, xi, yi, zi in zip(count(), x, y, z):
            bdf_file.write(f'GRID,{i+1},,{xi},{yi},{zi}\n')

        for eid, n1i, n2i, n3i, n4i, n5i, n6i, n7i, n8i in zip(count(), n1, n2, n3, n4, n5, n6, n7, n8):
            bdf_file.write(f'CHEXA   {eid+1:8d}{pid:8d}{n1i:8d}{n2i:8d}{n3i:8d}{n4i:8d}{n5i:8d}{n6i:8d}\n'
                           f'        {n7i:8d}{n8i:8d}\n')



if __name__ == '__main__':
    main()
