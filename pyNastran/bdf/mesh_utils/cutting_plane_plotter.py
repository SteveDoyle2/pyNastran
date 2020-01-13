"""Plots a model cutting plane"""
import numpy as np
from pyNastran.utils import int_version
try:
    import matplotlib.pyplot as plt
    import matplotlib
    MPL_VER = int_version('matplotlib', matplotlib.__version__)
    IS_MATPLOTLIB = True
except ImportError:  # pragma: no cover
    IS_MATPLOTLIB = False
    MPL_VER = None

from pyNastran.bdf.mesh_utils.cut_model_by_plane import (
    cut_edge_model_by_coord, cut_face_model_by_coord, export_face_cut)


def cut_and_plot_model(title, p1, p2, zaxis,
                       model, coord, nodal_result,
                       log, ytol,
                       plane_atol=1e-5,
                       csv_filename=None,
                       invert_yaxis=False,
                       cut_type='edge', plot=True, show=True):
    """
    Cuts a Nastran model with a cutting plane

    Parameters
    ----------
    title : str
        the title for the plot
    p1 : ???
        ???
    p2 : ???
        ???
    zaxis : ???
        ???
    model : str / BDF
        str : the bdf filename
        model : a properly configurated BDF object
    coord : Coord
        the coordinate system to cut the model with
    nodal_result : (nelements, ) float np.ndarray
        the result to cut the model with
    log : logger
        a logging object
    ytol : float
        the tolerance to filter edges (using some large value) to prevent
        excessive computations
    plane_atol : float; default=1e-5
        the tolerance for a line that's located on the y=0 local plane
    csv_filename : str; default=None
        None : don't write a csv
        str : write a csv
    """
    #assert cut_type.lower() in ['edge'], 'cut_type=%r and must be edge' % cut_type

    cut_type = cut_type.lower()
    if cut_type == 'edge':
        # bar cutting version
        csv_filename_edge = None
        #if csv_filename:
            #csv_filename_edge = csv_filename + '_edge.csv'
        local_points_array, global_points_array, result_array = cut_edge_model_by_coord(
            model, coord, ytol,
            nodal_result, plane_atol=plane_atol, csv_filename=csv_filename_edge)
        if len(local_points_array) == 0:
            log.error('No elements were piereced.  Check your cutting plane.')
            return
        if plot or csv_filename:
            plot_cutting_plane_edges(title, p1, p2, zaxis,
                                     local_points_array, global_points_array, result_array,
                                     csv_filename=csv_filename, invert_yaxis=invert_yaxis,
                                     plot=plot, show=show)
    elif cut_type == 'face':
        csv_filename_face = None
        if csv_filename:
            csv_filename_face = csv_filename + '_face.csv'
        geometry_arrays, results_arrays, unused_rods_array = cut_face_model_by_coord(
            model, coord, ytol,
            nodal_result, plane_atol=plane_atol,
            csv_filename=csv_filename_face)

        plot_cutting_plane_faces(title, p1, p2, zaxis,
                                 geometry_arrays, results_arrays,
                                 csv_filename=csv_filename, invert_yaxis=invert_yaxis,
                                 show=show)
    else:  # pragma: no cover
        raise NotImplementedError(f'cut_type={cut_type!r} and must be [edge, face]')


def plot_cutting_plane_faces(title, p1, p2, zaxis,
                             geometry_arrays, result_arrays,
                             csv_filename=None, invert_yaxis=False, show=True):
    """for faces"""
    try:
        local_x = result_arrays[:, 0]
    except TypeError:
        print('results_arrays', result_arrays)
        raise

    #x = global_points_array[:, 0]
    #y = global_points_array[:, 1]
    #z = global_points_array[:, 2]
    # bl
    x = result_arrays[:, -4]
    y = result_arrays[:, -3]
    z = result_arrays[:, -2]

    #x = local_points_array[:, 0]
    #y = local_points_array[:, 1]
    cp = result_arrays[:, -1]
    plt.close()
    fig = plt.figure(1)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(local_x, cp, '-o')
    ax.set_xlabel('Local X')
    ax.set_ylabel(title)
    ax.grid(True)
    if invert_yaxis:
        ax.invert_yaxis()

    fig = plt.figure(2)
    plt.plot(x, y, 'o')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid(True)

    fig = plt.figure(3)
    plt.plot(x, z, 'o')
    plt.xlabel('X')
    plt.ylabel('Z')
    plt.grid(True)
    if show:
        plt.show()

    csv_filename = '%s.csv' % title

    header = (
        #'# cut 1\n'
        'p1: %s\n'
        'p2: %s\n'
        'zaxis manual: %s\n'
        'x, y, z, %s\n' % (str(p1), str(p2), str(zaxis), title)
    )
    #for body in bodies:
    #with open(csv_filename, 'w') as csv_file:
    #np.savetxt(csv_filename, result_arrays[0], delimiter=',', header=header, comments='# ',
    #           encoding=None)
    export_face_cut(csv_filename, geometry_arrays, result_arrays, header=header)

def plot_cutting_plane_edges(title, p1, p2, zaxis,
                             local_points_array, global_points_array, result_array,
                             csv_filename=None, invert_yaxis=False, plot=True, show=True):
    """for edges"""
    local_x = local_points_array[:, 0]

    is_complex = np.iscomplexobj(result_array)

    global_x = global_points_array[:, 0]
    global_y = global_points_array[:, 1]
    global_z = global_points_array[:, 2]

    #x = local_points_array[:, 0]
    #y = local_points_array[:, 1]
    nrows, ncols = result_array.shape
    assert ncols == 8, result_array.shape


    ncp = ncols - 7
    nresult_cols = ncp
    if is_complex:
        nresult_cols = 2 * ncp

    cp = result_array[:, -1]
    if plot:
        assert IS_MATPLOTLIB, IS_MATPLOTLIB
        plt.close()
        fig = plt.figure(1)
        ax = fig.add_subplot(1, 1, 1)

        colors = ['C0<', 'C1>'] if MPL_VER >= [2, 1] else ['b>', 'r<']
        ax.plot(local_x, cp.real, colors[0], label='real')
        if is_complex:
            ax.plot(local_x, cp.imag, colors[0], label='imag')
        ax.set_xlabel('Local X')
        ax.set_ylabel(title)
        ax.legend()
        ax.grid(True)
        if invert_yaxis:
            ax.invert_yaxis()

        fig = plt.figure(2)
        plt.plot(global_x, global_y, 'o')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.grid(True)

        fig = plt.figure(3)
        plt.plot(global_x, global_z, 'o')
        plt.xlabel('X')
        plt.ylabel('Z')
        plt.grid(True)

    if csv_filename is None:
        if plot and show:
            plt.show()
        return

    #for body in bodies:
    #with open(csv_filename, 'w') as csv_file:
    #Cp = result_array[:, 3:]

    nd_msg = ''

    # x, y, z, Cp
    result_array2 = np.zeros((nrows, 6 + nresult_cols), dtype='float64')

    # xyz
    result_array2[:, :3] = global_points_array
    result_array2[:, 3:6] = local_points_array

    if is_complex:
        if ncp > 1:
            nd_msg = ' (real, imag for each component)'

        # convert to float64
        header = (
            'cut 1\n'
            'p1: %s\n'
            'p2: %s\n'
            'zaxis manual: %s\n'
            'x, y, z, x_local, y_local, z_local, %s_real, %s_imag%s\n' % (
                str(p1), str(p2), str(zaxis),
                title, title, nd_msg)
        )
        # Cp
        for i in range(ncp):
            result_array2[:, 6 + 2*i] = cp.real
            result_array2[:, 7 + 2*i] = cp.imag
    else:
        if ncp > 1:
            nd_msg = ' (each component)'

        result_array2[:, -ncp] = result_array[:, -ncp]
        header = (
            'cut 1\n'
            'p1: %s\n'
            'p2: %s\n'
            'zaxis manual: %s\n'
            'x, y, z, x_local, y_local, z_local, %s%s\n' % (
                str(p1), str(p2), str(zaxis), title, nd_msg)
        )

    # encoding=None - added in numpy 1.14
    np.savetxt(csv_filename, result_array2, delimiter=',', header=header, comments='# ',)
    if plot and show:
        plt.show()
