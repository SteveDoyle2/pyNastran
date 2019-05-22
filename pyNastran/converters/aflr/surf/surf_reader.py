import os
from numpy import zeros, cross, vstack, array
from numpy.linalg import norm  # type: ignore


def combine_surfs(surf_filenames, surf_out_filename=None):
    """
    Combines multiple SURFs into a single file

    Parameters
    ----------
    surf_filenames : List[str]
        list of surf filenames
    surf_out_filename : str; default=None -> no writing
        string of stl output filename

    Returns
    -------
    surf : SurfReader
        the surf object
    """
    nodes = []
    tris = []
    quads = []
    tri_props = []
    quad_props = []

    n0 = 0
    for fname in surf_filenames:
        surf = SurfReader()
        surf.read_surf(fname)
        nnodes = surf.nodes.shape[0]
        nodes.append(surf.nodes)
        tris.append(surf.tris + n0)
        quads.append(surf.quads + n0)
        tri_props.append(surf.tri_props)
        quad_props.append(surf.quad_props)
        n0 += nnodes

    surf.nodes = vstack(nodes)
    surf.tris = vstack(tris)
    surf.quads = vstack(quads)
    surf.tri_props = vstack(tri_props)
    surf.quad_props = vstack(quad_props)

    if surf_out_filename is not None:
        surf.write_surf(surf_out_filename)
    return surf


class SurfReader:
    def __init__(self, log=None, debug=False):
        """
        Initializes the SurfReader object

        Parameters
        ----------
        debug : bool/None; default=True
            used to set the logger if no logger is passed in
                True:  logs debug/info/error messages
                False: logs info/error messages
                None:  logs error messages
        log : logging module object / None
            if log is set, debug is ignored and uses the
            settings the logging object has
        """
        self.log = log
        self.debug = debug

        self.nodes = None
        self.node_props = None
        #self.nodes_failed = None
        self.tris = None
        self.tri_props = None
        self.quads = None
        self.quad_props = None

        self.nodes_failed = None

    def read_surf(self, surf_filename):
        """
        # AFLR3,14.18.115,28.206.154,WIN/doc/ug_io/3d_input_output_grids.html

        +---------+--------------------------------------------------------------+
        | Grid BC |                         Description                          |
        +=========+==============================================================+
        |    -6   | internal embedded/transparent surface that will be converted |
        |         | to internal/interior/volume faces with BL volume grid        |
        +---------+--------------------------------------------------------------+
        |    -5   | embedded/transparent surface with BL volume grid             |
        +---------+--------------------------------------------------------------+
        |    -1   | standard surface with BL volume grid; **wall**               |
        +---------+--------------------------------------------------------------+
        |     0   | standard surface                                             |
        +---------+--------------------------------------------------------------+
        |     1   | standard surface; **farfield**                               |
        +---------+--------------------------------------------------------------+
        |     2   | standard surface that intersects the BL region;              |
        |         | **boundary layer**                                           |
        +---------+--------------------------------------------------------------+
        |     3   | embedded/transparent surface or source surface that will be  |
        |         | converted to source nodes; **source**                        |
        +---------+--------------------------------------------------------------+
        |     4   | embedded/transparent surface that intersects the BL region   |
        +---------+--------------------------------------------------------------+
        |     5   | embedded/transparent surface                                 |
        +---------+--------------------------------------------------------------+
        |     6   | internal embedded/transparent surface that will be converted |
        |         | to internal/interior/volume faces                            |
        +---------+--------------------------------------------------------------+
        |     7   | fixed surface that intersects and directly connects to the   |
        |         | BL region; **transparent**                                   |
        +---------+--------------------------------------------------------------+
        """
        #basename = os.path.splitext(surf_filename)[0]
        #fail_surf_filename = basename +
        #print(fail_surf_filename)


        is_failed_surf = True if '.FAIL.surf' in surf_filename else False
        # if '.FAIL.' in surf_filename:
            # is_failed_surf = True

        with open(surf_filename, 'r') as surf_file:
            ntris, nquads, nnodes = surf_file.readline().strip().split()
            ntris = int(ntris)
            nquads = int(nquads)
            nnodes = int(nnodes)

            nodes = zeros((nnodes, 3), dtype='float64')
            node_props = zeros((nnodes, 2), dtype='float64')

            tris = zeros((ntris, 3), dtype='int32')
            quads = zeros((nquads, 4), dtype='int32')

            tri_props = zeros((ntris, 3), dtype='int32')
            quad_props = zeros((nquads, 3), dtype='int32')


            # read nodes
            if is_failed_surf:
                for inode in range(nnodes):
                    (x, y, z) = surf_file.readline().strip().split()
                    nodes[inode, :] = [x, y, z]
            else:
                for inode in range(nnodes):
                    x_y_z_initalnormalspacing_blthickness = surf_file.readline().strip().split()
                    if len(x_y_z_initalnormalspacing_blthickness) == 4:
                        (x, y, z, initial_normal_spacing) = x_y_z_initalnormalspacing_blthickness
                        bl_thickness = 0.0
                    elif len(x_y_z_initalnormalspacing_blthickness) == 5:
                        (x, y, z, initial_normal_spacing,
                         bl_thickness) = x_y_z_initalnormalspacing_blthickness
                    else:
                        raise NotImplementedError(x_y_z_initalnormalspacing_blthickness)
                    nodes[inode, :] = [x, y, z]
                    node_props[inode, :] = [initial_normal_spacing, bl_thickness]

            # read tris
            for itri in range(ntris):
                (n1, n2, n3, surface_id, reconnection_flag,
                 grid_bc_flag) = surf_file.readline().strip().split()
                tris[itri, :] = [n1, n2, n3]
                tri_props[itri, :] = [surface_id, reconnection_flag, grid_bc_flag]

            # read quads
            for iquad in range(nquads):
                (n1, n2, n3, n4, surface_id, reconnection_flag,
                 grid_bc_flag) = surf_file.readline().strip().split()
                quads[iquad, :] = [n1, n2, n3, n4]
                quad_props[iquad, :] = [surface_id, reconnection_flag, grid_bc_flag]

        self.nodes = nodes
        self.node_props = node_props
        self.tris = tris
        self.tri_props = tri_props
        self.quads = quads
        self.quad_props = quad_props

    def read_surf_failnode(self, surf_filename):
        #print(surf_filename)
        basename = os.path.splitext(surf_filename)[0]
        fail_node_filename = basename + '.FAIL.node'

        # nnodes = self.nodes.shape[0]
        nodes_failed = []
        if os.path.exists(fail_node_filename):
            with open(fail_node_filename, 'r') as fail_file:
                lines = fail_file.readlines()
                for line in lines:
                    nid = int(line.split()[0])
                    nodes_failed.append(nid)
            nodes_failed = array(nodes_failed, dtype='int32')
        self.nodes_failed = nodes_failed

    def get_normals(self):
        ntris = self.tris.shape[0]
        nquads = self.quads.shape[0]
        nelements = ntris + nquads
        normals = zeros((nelements, 3), dtype='float64')
        n1 = self.tris[:, 0] - 1
        n2 = self.tris[:, 1] - 1
        n3 = self.tris[:, 2] - 1
        a = self.nodes[n2, :] - self.nodes[n1, :]
        b = self.nodes[n3, :] - self.nodes[n1, :]
        normal = cross(a, b)
        mag = norm(normal, axis=1)
        assert len(mag) == ntris, 'axis is wrong'
        normal /= mag[:, None]
        normals[:ntris, :] = normal

        n1 = self.quads[:, 0] - 1
        n2 = self.quads[:, 1] - 1
        n3 = self.quads[:, 2] - 1
        n4 = self.quads[:, 3] - 1
        #print(n1)
        #print(n3)
        a = self.nodes[n3, :] - self.nodes[n1, :]
        b = self.nodes[n4, :] - self.nodes[n2, :]
        normal = cross(a, b)
        mag = norm(normal, axis=1)
        assert len(mag) == nquads, 'axis is wrong'
        normal /= mag[:, None]
        normals[ntris:, :] = normal
        return normals


class TagReader:
    def __init__(self, log=None, debug=False):
        self.log = log
        self.debug = debug

    def read_tag_filename(self, tag_filename):
        with open(tag_filename, 'r') as tag_file:
            lines = tag_file.readlines()

        lines = [
            line.strip().split('#')[0].strip() for line in lines
            if line.strip().split('#')[0].strip()
        ]
        #print('\n'.join(lines))

        #ID Group  Visc  Recon  Rebuild Fixed  Source  Trans  Delete  Spacing Thcknss Layers
        data = {}
        for line in lines:
            sline = line.split()
            if len(sline) == 12:
                (ID, name, is_visc, is_recon, is_rebuild, is_fixed, is_source,
                 is_trans, is_delete, bl_spacing, bl_thickness, nlayers) = sline
                ID = int(ID)
                is_visc = int(is_visc)
                is_recon = int(is_recon)
                is_rebuild = int(is_rebuild)
                is_source = int(is_source)
                is_trans = int(is_trans)
                is_delete = int(is_delete)
                bl_spacing = float(bl_spacing)
                bl_thickness = float(bl_thickness)
                nlayers = int(nlayers)
            else:
                msg = 'invalid data line in tags file=%r\n' % tag_filename
                msg += 'split line = %s\n' % sline
                msg += 'expected 12 values; got %s\n' % (len(sline))
                try:
                    msg += '  ID           = %s\n' % sline[0]
                    msg += '  name         = %s\n' % sline[1]
                    msg += '  is_visc      = %s\n' % sline[2]
                    msg += '  is_recon     = %s\n' % sline[3]
                    msg += '  is_rebuild   = %s\n' % sline[4]
                    msg += '  is_fixed     = %s\n' % sline[5]
                    msg += '  is_source    = %s\n' % sline[6]
                    msg += '  is_trans     = %s\n' % sline[7]
                    msg += '  is_delete    = %s\n' % sline[8]
                    msg += '  bl_spacing   = %s\n' % sline[9]
                    msg += '  bl_thickness = %s\n' % sline[10]
                    msg += '  nlayers = %s\n' % sline[11]
                except IndexError:
                    pass
                raise RuntimeError(msg)

            data[ID] = [name, is_visc, is_recon, is_rebuild, is_fixed, is_source,
                        is_trans, is_delete, bl_spacing, bl_thickness, nlayers]
        return data
