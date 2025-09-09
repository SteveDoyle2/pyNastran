import os
from typing import Optional

import numpy as np
from pyNastran.utils import PathLike, print_bad_path
from pyNastran.converters.fluent.utils import (
    read_vrt, read_cell, read_daten,
    write_vrt, write_cell, write_daten,
    filter_by_region,
)
from cpylog import __version__ as CPYLOG_VERSION, SimpleLogger
if CPYLOG_VERSION > '1.6.0':
    from cpylog import get_logger
else:  # pragma: no cover
    from cpylog import get_logger2 as get_logger

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
        self.log = get_logger(log, debug)

        # vrt
        self.node_id = np.array([], dtype='int32')
        self.xyz = np.zeros((0, 3), dtype='int32')

        # cel
        self.quads = np.zeros((0, 6), dtype='int32')
        self.tris = np.zeros((0, 5), dtype='int32')

        # daten
        self.titles = np.zeros((0, ), dtype='|U8')
        self.result_element_id = np.array([], dtype='int32')
        self.results = np.zeros((0, 0), dtype='float64')

        # TODO: may be removed
        self.element_ids = np.array([], dtype='int32')

    @property
    def element_id(self) -> np.ndarray:
        return self.result_element_id

    @classmethod
    def from_data(self,
                  node_id: np.ndarray, xyz: np.ndarray,
                  tris: np.ndarray, quads: np.ndarray,
                  result_element_id: np.ndarray, titles: np.ndarray,
                  quad_results: np.ndarray, tri_results: np.ndarray,
                  auto_read_write_h5: bool=True,
                  log=None, debug: bool=True):
        model = Fluent(
            auto_read_write_h5=auto_read_write_h5,
            log=log, debug=debug)
        assert len(node_id) > 0
        assert len(result_element_id) > 0
        model.node_id = node_id
        model.xyz = xyz
        model.result_element_id = result_element_id
        model.quads = quads
        model.tris = tris
        model.titles = titles
        model.results = np.vstack([quad_results, tri_results])

        element_ids = np.unique(np.hstack([
            quads[:, 0], tris[:, 0],
        ]))
        model.element_ids = element_ids  # TODO: consider removing
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
        names = ['node_id', 'xyz',
                 'titles', 'results', 'result_element_id', #'element_id',
                 'quads', 'tris', 'element_ids']
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
        names = [
            'node_id', 'xyz',
            'titles', 'results', 'result_element_id', #'element_id',
            'quads', 'tris', 'element_ids',
        ]
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
        assert len(self.result_element_id) > 0, self.result_element_id
        assert len(self.element_ids) > 0, self.element_ids
        is_loaded = True
        return is_loaded

    def read_fluent(self, fluent_filename: PathLike) -> None:
        base, ext = os.path.splitext(fluent_filename)
        vrt_filename = base + '.vrt'
        daten_filename = base + '.daten'
        cell_filename = base + '.cel'
        h5_filename = base + '.h5'
        assert os.path.exists(vrt_filename), print_bad_path(vrt_filename)
        assert os.path.exists(daten_filename), print_bad_path(daten_filename)
        assert os.path.exists(cell_filename), print_bad_path(cell_filename)
        self.fluent_filename = fluent_filename
        if self.auto_read_write_h5 and os.path.exists(h5_filename):
            # fails if pytables isn't installed
            self.log.debug(f'reading fluent {h5_filename}')
            is_loaded = self.read_h5(h5_filename, require_load=False)
            if is_loaded:
                return
        assert os.path.exists(vrt_filename), print_bad_path(vrt_filename)
        assert os.path.exists(daten_filename), print_bad_path(daten_filename)
        assert os.path.exists(cell_filename), print_bad_path(cell_filename)

        (quads, tris), element_ids = read_cell(cell_filename, self.log)
        node, xyz = read_vrt(vrt_filename, self.log)
        result_element_id, titles, results = read_daten(daten_filename, scale=1.0)

        self.node_id = node
        self.xyz = xyz

        self.result_element_id = result_element_id
        self.titles = titles
        self.results = results

        self.quads = np.zeros((0, 6), dtype='int32')
        self.tris = np.zeros((0, 5), dtype='int32')
        if len(quads):
            self.quads = quads
        if len(tris):
            self.tris = tris

        # ids of the quad/tri elements in the order
        # respects interspersed quads/tris
        self.element_ids = element_ids  # TODO: consider removing

        assert len(result_element_id) > 0, result_element_id
        assert len(element_ids) > 0, element_ids
        if self.auto_read_write_h5:
            # fails if pytables isn't installed
            self.write_h5(h5_filename, require_write=False)

    def write_fluent(self, fluent_filename: str) -> None:
        assert len(self.node_id) > 0, self.node_id
        assert len(self.node_id) == len(self.xyz)
        assert self.tris.shape[1] == 5, self.tris.shape
        assert self.quads.shape[1] == 6, self.quads.shape
        assert len(self.result_element_id) > 0, self.result_element_id

        base, ext = os.path.splitext(fluent_filename)
        vrt_filename = base + '.vrt'
        daten_filename = base + '.daten'
        cell_filename = base + '.cel'

        #(quads, tris), (element_ids, region)
        # node, xyz = read_vrt(vrt_filename)
        # result_element_id, titles, results = read_daten(daten_filename, scale=1.0)
        write_cell(cell_filename, self.quads, self.tris)
        write_vrt(vrt_filename, self.node_id, self.xyz)

        results_elements_list = []
        results_list = []
        if len(self.tris):
            itri = np.searchsorted(self.result_element_id, self.tris[:, 0])
            results_elements_list.append(self.result_element_id[itri])
            results_list.append(self.results[itri, :])

        if len(self.quads):
            iquad = np.searchsorted(self.result_element_id, self.quads[:, 0])
            results_elements_list.append(self.result_element_id[iquad])
            results_list.append(self.results[iquad, :])
        result_element_id = np.hstack(results_elements_list)
        results = np.vstack(results_list)
        write_daten(daten_filename, result_element_id, self.titles, results)

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
        out = (
            tri_area, quad_area,
            tri_centroid, quad_centroid,
            tri_normal, quad_normal,
        )
        return out

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
                          return_model: bool=False,
                          #deepcopy: bool=True,
                          ) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                                     np.ndarray, np.ndarray]:
        # False
        #assert np.array_equal(self.result_element_id, np.unique(self.result_element_id))
        if regions_to_remove is None:
            regions_to_remove = []
        if regions_to_include is None:
            regions_to_include = []
        region_split = bool(len(regions_to_remove) + len(regions_to_include))
        assert len(self.result_element_id) > 0
        if region_split:
            self.log.info(f'regions_to_remove={regions_to_remove}; '
                          f'regions_to_include={regions_to_include}')
            result_element_id, tris, quads, quad_results, tri_results = filter_by_region(
                self, regions_to_remove, regions_to_include)
            assert np.array_equal(result_element_id, np.unique(result_element_id))
            if len(result_element_id) == 0 and 0:
                # error state to not crash gui
                region_id = self.region[0]
                self.log.warning(f'no elements remaining; including region_id={region_id}')
                return self.get_filtered_data(
                    regions_to_include=[region_id],
                    return_model=return_model)
            assert len(result_element_id) > 0, 'no elements remaining'

            model2 = self.from_data(
                self.node_id, self.xyz,
                tris, quads,
                result_element_id, self.titles,
                quad_results, tri_results,
                auto_read_write_h5=True,
                log=self.log)
            assert isinstance(model2, Fluent)
            str(model2)
            assert len(model2.result_element_id) > 0
            region = model2.region
            results = model2.results
        else:
            self.log.debug('no filtering of regions')
            tris = self.tris
            quads = self.quads

            # we reordered the tris/quads to be continuous to make them easier to add
            assert np.array_equal(quads[:, 0], np.unique(quads[:, 0]))
            assert np.array_equal(tris[:, 0], np.unique(tris[:, 0]))

            iquad = np.searchsorted(self.result_element_id, quads[:, 0])
            itri = np.searchsorted(self.result_element_id, tris[:, 0])
            quad_results = self.results[iquad, :]
            tri_results = self.results[itri, :]
            result_element_id = self.result_element_id
            assert np.array_equal(result_element_id, np.unique(result_element_id))

            region = self.region
            results = np.vstack([quad_results, tri_results])
            self.results = results
            assert len(result_element_id) > 0, 'no elements remaining2'
            assert len(self.result_element_id) > 0
            model2 = self
        if return_model:
            return model2
        return result_element_id, tris, quads, region, results

    @property
    def region(self) -> np.ndarray:
        return np.hstack([self.quads[:, 1], self.tris[:, 1]])

    #@property
    #def results(self) -> np.ndarray:
        #return np.vstack([self.quad_results, self.tri_results])

    def __repr__(self) -> str:
        msg = (
            'Fluent:\n'
            f' - node_id = {self.node_id}\n'
            f' - xyz =\n{self.xyz}\n'
            f' - result_element_id = {self.result_element_id}\n'
            f' - quads = {self.quads}\n'
            f' - tris  = {self.tris}\n'
            f' - titles = {self.titles}\n'
            f' - results =\n{self.results}\n'
            f' - element_ids = {self.element_ids}\n'
            f' - region = {self.region}\n'
        )
        return msg

def read_fluent(fluent_filename: PathLike,
                auto_read_write_h5: bool=True,
                log=None, debug: bool=False) -> Fluent:
    model = Fluent(auto_read_write_h5=auto_read_write_h5,
                   log=log, debug=debug)
    model.read_fluent(fluent_filename)
    return model
