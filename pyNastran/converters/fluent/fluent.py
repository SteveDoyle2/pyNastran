import os
from typing import Optional

import numpy as np
from cpylog import SimpleLogger, get_logger2
from pyNastran.utils import PathLike, print_bad_path

class Fluent:
    def __init__(self, auto_read_write_h5: bool=True,
                 log=None, debug: bool=True):
        """
        Parameters
        ----------
        auto_read_write_h5: bool; default=True
            the vrt/cel/daten files are really slow, so save an h5 file
            and load it instead
        log: SimpleLogger or None
            the log
        debug: bool; default=False
            logging debug
        """
        self.auto_read_write_h5 = auto_read_write_h5
        self.log = get_logger2(log=log, debug=debug)

        # vrt
        self.node_id = np.array([], dtype='int32')
        self.xyz = np.zeros((0, 3), dtype='int32')

        # cel
        self.quads = np.zeros((0, 6), dtype='int32')
        self.tris = np.zeros((0, 5), dtype='int32')

        # daten
        self.titles = np.zeros((0, ), dtype='|U8')
        self.element_id = np.array([], dtype='int32')
        self.results = np.zeros((0, 0), dtype='float64')

    @classmethod
    def from_data(self,
                  node_id: np.ndarray, xyz: np.ndarray,
                  tris: np.ndarray, quads: np.ndarray,
                  element_id: np.ndarray, titles: np.ndarray, results: np.ndarray,
                  element_ids: np.ndarray, region: np.ndarray,
                  auto_read_write_h5: bool=True,
                 log=None, debug: bool=True):
        model = Fluent(
            auto_read_write_h5=auto_read_write_h5,
            log=log, debug=debug)
        model.node_id = node_id
        model.xyz = xyz
        model.element_id = element_id
        model.quads = quads
        model.tris = tris
        model.titles = titles
        model.results = results
        model.element_ids = element_ids  # TODO: consider removing
        model.region = region            # TODO: consider removing
        return model

    def _get_h5_filename(self, h5_filename: PathLike) -> str:
        if h5_filename == '':
            base, ext = os.path.splitext(self.fluent_filename)
            h5_filename = base + '.h5'
        return h5_filename

    def write_h5(self, h5_filename: PathLike='',
                 require_write: bool=True) -> bool:
        is_written = False
        try:
            import tables
        except ImportError:
            if require_write:
                raise
            return is_written
        h5_filename = self._get_h5_filename(h5_filename)

        h5file = tables.open_file(h5_filename, 'w', driver='H5FD_CORE')
        names = ['node_id', 'xyz', 'element_id', 'titles', 'results',
                 'quads', 'tris', 'element_ids', 'region']
        h5file.create_array(h5file.root, 'format', ['fluent'])
        for name in names:
            datai = getattr(self, name)
            h5file.create_array(h5file.root, name, datai)

        h5file.close()
        self.fluent_filename = h5_filename
        is_written = True
        return is_written

    def read_h5(self, h5_filename: PathLike='',
                require_load: bool=True) -> bool:
        is_loaded = False
        try:
            import tables
        except ImportError:
            if require_load:
                raise
            return is_loaded
        names = ['node_id', 'xyz', 'element_id', 'titles', 'results',
                 'quads', 'tris', 'element_ids', 'region']
        h5_filename = self._get_h5_filename(h5_filename)

        assert os.path.exists(h5_filename), print_bad_path(h5_filename)
        h5file = tables.open_file(h5_filename, 'r', driver='H5FD_CORE')

        try:
            table = h5file.get_node('/format')
            formati = table.read()
            assert formati == [b'fluent'], f'format={formati!r}'
        except Exception:
            print('---error no format?---')
            print(h5file)

        for name in names:
            table = h5file.get_node(f'/{name}')
            datai = table.read()
            if name == 'titles':
                # strings are stored as bytes so cast them back to strings
                data_list = [title.decode('latin1') for title in datai]
                data_array = np.array(data_list)
                setattr(self, name, data_array)
            else:
                setattr(self, name, datai)
        h5file.close()
        is_loaded = True
        return is_loaded

    def read_fluent(self, fluent_filename: PathLike) -> None:
        base, ext = os.path.splitext(fluent_filename)
        vrt_filename = base + '.vrt'
        daten_filename = base + '.daten'
        cell_filename = base + '.cel'
        h5_filename = base + '.h5'
        self.fluent_filename = fluent_filename
        if self.auto_read_write_h5 and os.path.exists(h5_filename):
            # fails if pytables isn't installed
            is_loaded = self.read_h5(h5_filename, require_load=False)
            if is_loaded:
                return
        assert os.path.exists(vrt_filename), print_bad_path(vrt_filename)
        assert os.path.exists(daten_filename), print_bad_path(daten_filename)
        assert os.path.exists(cell_filename), print_bad_path(cell_filename)

        (quads, tris), (element_ids, region) = read_cell(cell_filename)
        node, xyz = read_vrt(vrt_filename)
        element_id, titles, results = read_daten(daten_filename, scale=1.0)

        #tri_centroid, tri_pressure = tri_split(xyz, tris, element_id, results)
        #quad_centroid, quad_pressure = quad_split(xyz, quads, element_id, results)
        self.node_id = node
        self.xyz = xyz
        self.element_id = element_id  # result element ids
        self.titles = titles
        self.results = results

        self.quads = np.zeros((0, 6), dtype='int32')
        self.tris = np.zeros((0, 5), dtype='int32')
        if len(quads):
            self.quads = quads
        if len(tris):
            self.tris = tris

        self.element_ids = element_ids  # TODO: consider removing
        self.region = region            # TODO: consider removing
        if self.auto_read_write_h5:
            # fails if pytables isn't installed
            self.write_h5(h5_filename, require_write=False)

    def write_fluent(self, fluent_filename: str) -> None:
        assert len(self.node_id) > 0, self.node_id
        assert len(self.node_id) == len(self.xyz)
        assert self.tris.shape[1] == 5, self.tris.shape
        assert self.quads.shape[1] == 6, self.quads.shape
        assert len(self.element_id) > 0, self.element_id

        base, ext = os.path.splitext(fluent_filename)
        vrt_filename = base + '.vrt'
        daten_filename = base + '.daten'
        cell_filename = base + '.cel'

        #(quads, tris), (element_ids, region)
        # node, xyz = read_vrt(vrt_filename)
        # element_id, titles, results = read_daten(daten_filename, scale=1.0)
        write_cell(cell_filename, self.quads, self.tris)
        write_vrt(vrt_filename, self.node_id, self.xyz)

        elements_list = []
        results_list = []
        if len(self.tris):
            itri = np.searchsorted(self.element_id, self.tris[:, 0])
            elements_list.append(self.element_id[itri])
            results_list.append(self.results[itri, :])

        if len(self.quads):
            iquad = np.searchsorted(self.element_id, self.quads[:, 0])
            elements_list.append(self.element_id[iquad])
            results_list.append(self.results[iquad, :])
        element_id = np.hstack(elements_list)
        results = np.vstack(results_list)
        write_daten(daten_filename, element_id, self.titles, results)

    def get_area_centroid_normal(self, tris: np.ndarray,
                                 quads: np.ndarray) -> tuple[np.ndarray, np.ndarray,
                                                             np.ndarray, np.ndarray,
                                                             np.ndarray, np.ndarray]:
        tri_nodes = tris[:, 2:]
        quad_nodes = quads[:, 2:]
        assert tri_nodes.shape[1] == 3, tri_nodes.shape
        assert quad_nodes.shape[1] == 4, quad_nodes.shape

        utri_nodes = np.unique(tri_nodes.ravel())
        uquad_nodes = np.unique(quad_nodes.ravel())
        assert len(np.setdiff1d(uquad_nodes, self.node_id)) == 0
        assert len(np.setdiff1d(utri_nodes, self.node_id)) == 0

        itri = np.searchsorted(self.node_id, tri_nodes)
        tri_xyz1 = self.xyz[itri[:, 0]]
        tri_xyz2 = self.xyz[itri[:, 1]]
        tri_xyz3 = self.xyz[itri[:, 2]]
        tri_centroid = (tri_xyz1 + tri_xyz2 + tri_xyz3) / 3

        tri_normal = np.cross(tri_xyz2 - tri_xyz1, tri_xyz3 - tri_xyz1)
        tri_areai = np.linalg.norm(tri_normal, axis=1)
        tri_area = 0.5 * tri_areai
        tri_normal /= tri_areai[:, np.newaxis]

        iquad = np.searchsorted(self.node_id, quad_nodes)
        xyz1 = self.xyz[iquad[:, 0]]
        xyz2 = self.xyz[iquad[:, 1]]
        xyz3 = self.xyz[iquad[:, 2]]
        xyz4 = self.xyz[iquad[:, 3]]

        quad_centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4
        quad_normal = np.cross(xyz3 - xyz1, xyz4 - xyz2)
        quad_areai = np.linalg.norm(quad_normal, axis=1)
        quad_area = 0.5 * quad_areai
        quad_normal /= quad_areai[:, np.newaxis]
        assert len(tri_centroid) == len(tri_area)
        assert len(quad_centroid) == len(quad_area)
        return tri_area, quad_area, tri_centroid, quad_centroid, tri_normal, quad_normal

    def get_area_centroid(self, tris: np.ndarray,
                          quads: np.ndarray) -> tuple[np.ndarray, np.ndarray,
                                                      np.ndarray, np.ndarray]:
        out = self.get_area_centroid_normal(tris, quads)
        tri_area, quad_area, tri_centroid, quad_centroid, tri_normal, quad_normal = out
        del tri_normal, quad_normal
        return tri_area, quad_area, tri_centroid, quad_centroid


    def get_filtered_data(self,
                          regions_to_remove: Optional[list[int]]=None,
                          regions_to_include: Optional[list[int]]=None,
                          ) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                                     np.ndarray, np.ndarray]:
        if regions_to_remove is None:
            regions_to_remove = []
        if regions_to_include is None:
            regions_to_include = []
        region_split = bool(len(regions_to_remove) + len(regions_to_include))
        if region_split:
            element_id, tris, quads, region, results, quad_results, tri_results = filter_by_region(
                self, regions_to_remove, regions_to_include)
            if 0:
                element_ids = np.unique(np.hstack([
                    quads[:, 0], tris[:, 0],
                ]))
                model2 = self.from_data(
                    self.node_id, self.xyz,
                    tris, quads,
                    element_id, self.titles, results,
                    element_ids, region,
                    auto_read_write_h5=True,
                    log=self.log)
                assert isinstance(model2, Fluent)
        else:
            tris = self.tris
            quads = self.quads
            tri_regions = tris[:, 1]
            quad_regions = quads[:, 1]

            # we reordered the tris/quads to be continuous to make them easier to add
            iquad = np.searchsorted(self.element_id, quads[:, 0])
            itri = np.searchsorted(self.element_id, tris[:, 0])
            quad_results = self.results[iquad, :]
            tri_results = self.results[itri, :]
            element_id = self.element_id
            region = np.hstack([quad_regions, tri_regions])
            results = np.vstack([quad_results, tri_results])
        return element_id, tris, quads, region, results #, quad_results, tri_results

def filter_by_region(model: Fluent,
                     regions_to_remove: list[int],
                     regions_to_include: list[int]) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                                                             np.ndarray, np.ndarray]:
    tris = model.tris
    quads = model.quads
    results = model.results

    tri_regions = tris[:, 1]
    quad_regions = quads[:, 1]

    # we reordered the tris/quads to be continuous to make them easier to add
    iquad = np.searchsorted(model.element_id, quads[:, 0])
    itri = np.searchsorted(model.element_id, tris[:, 0])
    #-----------------------------
    is_remove = (len(regions_to_remove) == 0)
    is_include = (len(regions_to_include) == 0)
    assert (is_remove and not is_include) or (not is_remove and is_include)
    if regions_to_remove:
        itri_regions = np.logical_and.reduce([(tri_regions != regioni) for regioni in regions_to_remove])
        iquad_regions = np.logical_and.reduce([(quad_regions != regioni) for regioni in regions_to_remove])
    else:
        itri_regions = np.logical_or.reduce([(tri_regions == regioni) for regioni in regions_to_include])
        iquad_regions = np.logical_or.reduce([(quad_regions == regioni) for regioni in regions_to_include])

    quad_results = results[iquad, :][iquad_regions, :]
    tri_results = results[itri, :][itri_regions, :]

    tris = tris[itri_regions, :]
    quads = quads[iquad_regions, :]
    tri_eids = tris[:, 0]
    quad_eids = quads[:, 0]
    element_id = np.unique(np.hstack([quad_eids, tri_eids]))
    tri_regions = tris[:, 1]
    quad_regions = quads[:, 1]
    region = np.hstack([quad_regions, tri_regions])
    results = np.vstack([quad_results, tri_results])
    return element_id, tris, quads, region, results, quad_results, tri_results


def write_daten(daten_filename: PathLike,
                element_id: np.ndarray,
                titles: np.ndarray,
                results: np.ndarray) -> None:
    assert len(titles) > 1, titles
    titles = titles[1:]
    ntitles = len(titles)
    title_str = '# Shell Id, ' + ', '.join(list(titles)) + '\n'
    fmt = '%-8d' + ' %15s' * ntitles + '\n'
    with open(daten_filename, 'w') as daten_file:
        daten_file.write(title_str)
        for eid, res in zip(element_id, results.tolist()):
            daten_file.write(fmt % tuple([eid] + res))

def read_daten(daten_filename: PathLike,
               scale: float=1.0) -> tuple[np.ndarray, np.ndarray]:
    """
    # Shell Id, Cf Components X, Cf Components Y, Cf Components Z, Pressure Coefficient
         1       0.00158262927085      -0.00013470219276       0.00033890815696      -0.09045402854837
   4175456       0.00005536470364      -0.00050913009274       0.00226408640311      -0.44872838534457

    """
    with open(daten_filename, 'r') as vrt_file:
        lines = vrt_file.readlines()

    element_id_list = []
    data_list = []
    for i, line in enumerate(lines[1:]):
        sline = line.strip().split()
        eid, *datai = sline
        element_id_list.append(eid)
        data_list.append(datai)
    element_id = np.array(element_id_list, dtype='int32')

    titles_sline1 = lines[0].split('#')[1].strip().split(',')
    titles_sline2 = [val.strip() for val in titles_sline1]
    titles = np.array(titles_sline2)

    results = np.array(data_list, dtype='float64')
    if scale != 1.0:
        results *= scale

    # we drop the element_id from the check
    assert results.shape[1] == len(titles)-1, f'shape={str(results.shape)}; titles={titles}'
    return element_id, titles, results

def write_vrt(vrt_filename: PathLike, node_id: str, xyz: np.ndarray) -> None:
    title_str = 'PROSTAR_VERTEX\n'
    title_str += '4000         0         0         0         0         0         0         0\n'
    fmt = '%-8d' + ' %15s' * 3 + '\n'
    with open(vrt_filename, 'w') as vrt_file:
        vrt_file.write(title_str)
        for nid, (x, y, z) in zip(node_id, xyz):
            vrt_file.write(fmt % (nid, x, y, z))
    return

def read_vrt(vrt_filename: PathLike) -> tuple[np.ndarray, np.ndarray]:
    """
     PROSTAR_VERTEX
     4000         0         0         0         0         0         0         0
             1     10.86546244027553      0.15709716073078      1.16760397346681
       2118954      8.91578650948743      1.88831278193532      1.27715047467296
    """
    with open(vrt_filename, 'r') as vrt_file:
        lines = vrt_file.readlines()

    node_id_list = []
    xyz_list = []
    for i, line in enumerate(lines[2:]):
        sline = line.strip().split()
        nid, *xyzi = sline
        node_id_list.append(nid)
        xyz_list.append(xyzi)
    node_id = np.array(node_id_list, dtype='int32')
    xyz = np.array(xyz_list, dtype='float64')
    return node_id, xyz

def write_cell(vrt_filename: PathLike, quads: np.ndarray, tris: np.ndarray) -> None:
    title_str = 'PROSTAR_CELL\n'
    title_str += '4000         0         0         0         0         0         0         0\n'
    dim_fmt  = '%-8d %8s %8s %8d %8s\n'
    tri_fmt  = '%-8d %8s %8s %8s\n'
    quad_fmt = '%-8d %8s %8s %8s %8s\n'

    ntri, ntri_col = tris.shape
    nquad, nquad_col = quads.shape
    assert ntri_col == 5, tris.shape
    assert nquad_col == 6, quads.shape
    with open(vrt_filename, 'w') as cel_file:
        # row1 = [eid, dim, pid, 3, 4]
        # row2 = [eid, n1, n2, n3, n4]
        cel_file.write(title_str)
        for eid, pid, n1, n2, n3 in zip(tris[:, 0], tris[:, 1],
                                        tris[:, 2], tris[:, 3], tris[:, 4]):
            # {eid}      {dim}      {pid}          3          4
            # {eid}        407        406        472        473  #  nnodes = dim

            #cel_file.write(f'{eid}          3            {dim}         {pid}         4\n'
            #               f'{eid}          {n1}          {n2}          {n3}\n')
            cel_file.write(dim_fmt % (eid, '3', '3', pid, '4'))
            cel_file.write(tri_fmt % (eid, n1, n2, n3))

        for eid, pid, n1, n2, n3, n4 in zip(quads[:, 0], quads[:, 1],
                                            quads[:, 2], quads[:, 3], quads[:, 4], quads[:, 5]):
            cel_file.write(dim_fmt % (eid, '3', '4', pid, '4'))
            cel_file.write(quad_fmt % (eid, n1, n2, n3, n4))
    return

def read_cell(cell_filename: PathLike) -> tuple[tuple[np.ndarray, np.ndarray],
                                                tuple[np.ndarray, np.ndarray]]:
    """
PROSTAR_CELL
     4000         0         0         0         0         0         0         0
         1          3          3          3          4
         1          3          4          5
   ...
   4175456          3          3         27          4
   4175456    2118943    2118954    2118941

    """
    with open(cell_filename, 'r') as cell_file:
        lines = cell_file.readlines()

    element_ids_list = []
    regions_list = []

    quad_list = []
    tri_list = []
    i = 2
    while i < len(lines):
        sline1 = lines[i].strip().split()
        sline2 = lines[i+1].strip().split()

        idi1, *ids1 = sline1
        idi2, *ids2 = sline2
        assert idi1 == idi2, (sline1, sline2)
        assert len(ids1) == 4, ids1

        #790          3          3          3          4
        #790        406        403        472
        #791          3          4          3          4
        #791        407        406        472        473
        #
        #{eid}      {dim}      {pid}          3          4
        #{eid}        407        406        472        473  #  nnodes = dim
        dim = int(ids1[1])
        pid = ids1[2]
        assert ids1[0] == '3', (idi1, ids1)
        assert ids1[1] in '34', (idi1, ids1)
        #assert ids1[2] in '34', (idi1, ids1)
        assert ids1[3] == '4', (idi1, ids1)
        nids1 = len(ids1)
        nids2 = len(ids2)
        assert nids2 == dim, (nids2, dim)
        ids = [idi1, pid] + list(ids2)
        element_ids_list.append(idi1)
        regions_list.append(pid)
        if dim == 3:
            tri_list.append(ids)
        elif dim == 4:
            quad_list.append(ids)
        else:  # pragma: no cover
            raise RuntimeError(dim)
        assert nids1 in {4}, ids1
        assert nids2 in {3, 4}, ids2
        i += 2

    quads = np.array(quad_list, dtype='int32')
    tris = np.array(tri_list, dtype='int32')
    if len(tri_list):
        assert tris.shape[1] == 2+3, (tris.shape, quads.shape)
    if len(quad_list):
        assert quads.shape[1] == 2+4, (tris.shape, quads.shape)

    element_ids = np.array(element_ids_list, dtype='int32')
    regions = np.array(regions_list, dtype='int32')
    return (quads, tris), (element_ids, regions)

#def tri_split(xyz: np.ndarray,
              #tris: np.ndarray,
              #element_id: np.ndarray,
              #results: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    # [eid, pid, n1, n2, n3]
    #eids = tris[:, 0]
    #n1 = tris[:, 2] - 1
    #n2 = tris[:, 3] - 1
    #n3 = tris[:, 4] - 1
    #xyz1 = xyz[n1, :]
    #xyz2 = xyz[n2, :]
    #xyz3 = xyz[n3, :]
    #tri_centroid = (xyz1 + xyz2 + xyz3) / 3
    #ieid = np.searchsorted(element_id, eids)
    #tri_pressure = results[ieid]
    #return tri_centroid, tri_pressure

#def quad_split(xyz: np.ndarray,
               #quads: np.ndarray,
               #element_id: np.ndarray,
               #results: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    #eids = quads[:, 0]
    #n1 = quads[:, 2] - 1
    #n2 = quads[:, 3] - 1
    #n3 = quads[:, 4] - 1
    #n4 = quads[:, 5] - 1
    #xyz1 = xyz[n1, :]
    #xyz2 = xyz[n2, :]
    #xyz3 = xyz[n3, :]
    #xyz4 = xyz[n4, :]
    #quad_centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4
    #ieid = np.searchsorted(element_id, eids)
    #quad_pressure = results[ieid]
    #return quad_centroid, quad_pressure

def read_fluent(fluent_filename: PathLike,
                auto_read_write_h5: bool=True,
                log=None, debug: bool=False) -> Fluent:
    model = Fluent(auto_read_write_h5=auto_read_write_h5,
                   log=log, debug=debug)
    model.read_fluent(fluent_filename)
    return model
