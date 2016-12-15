from numpy import vstack, amax, amin, arange, ones, zeros, where, abs, degrees
from pyNastran.converters.aflr2.aflr2 import AFLR2

from six import iteritems
from six.moves import range
import os
#from numpy import zeros, arange, mean, amax, amin, array, where

import vtk
from vtk import vtkLine
from pyNastran.utils import print_bad_path


class BEDGE_IO(object):
    def __init__(self):
        pass

    def get_bedge_wildcard_geometry_results_functions(self):
        data = (
            'AFLR3 BEDGE',
            'AFLR3 BEDGE (*.bedge)', self.load_bedge_geometry,
             None, None)
        return data

    def load_bedge_geometry(self, bedge_filename, dirname, plot=True):
        #skipReading = self.remove_old_openfoam_geometry(openfoam_filename)
        #if skipReading:
        #    return

        model = AFLR2()

        self.modelType = 'bedge'
        print('bedge_filename = %s' % bedge_filename)

        model.read_bedge(bedge_filename)
        nnodes = model.nodes.shape[0]
        nbars = model.bars.shape[0]
        nelements = nbars

        nodes = model.nodes
        self.nElements = nelements
        self.nNodes = nnodes

        print("nNodes = %s" % self.nNodes)
        print("nElements = %s" % self.nElements)
        assert nelements > 0, nelements

        self.grid.Allocate(self.nElements, 1000)
        self.grid2.Allocate(1, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)

        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.update_axes_length(dim_max)
        self.log.info('max = %s' % mmax)
        self.log.info('min = %s' % mmin)

        for inode, node in enumerate(nodes):
            points.InsertPoint(inode, node)

        #tris = model.tris - 1
        #quads = model.quads - 1
        bars = model.bars

        #print('bars!!!!!!')
        for eid, element in enumerate(bars):
            #print(element)
            elem = vtk.vtkLine()
            n1, n2 = element
            try:
                elem.GetPointIds().SetId(0, n1)
                elem.GetPointIds().SetId(1, n2)
            except KeyError:
                print("nodeIDs =", element)
                print(str(element))
                continue

            self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

        self.nElements = nelements
        self.grid.SetPoints(points)
        self.grid.Modified()
        #self.grid2.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()
            #self.grid2.Update()
        #print("updated grid")

        # loadCart3dResults - regions/loads
        self.TurnTextOn()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        self.iSubcaseNameMap = {1: ['AFLR BEDGE', '']}
        cases = {}
        ID = 1

        form, cases = self._fill_bedge_case(bedge_filename, cases, ID, nnodes, nelements, model)
        if plot:
            self._finish_results_io2(form, cases)

    def clear_bedge(self):
        pass

    def _load_bedge_results(self, openfoam_filename, dirname):
        model = Cart3DReader(log=self.log, debug=False)
        #self.modelType = model.modelType
        #(nodes, elements, regions, loads) = model.read_cart3d(cart3dFileName)

        model.infilename = cart3d_filename
        if is_binary(infilename):
            model.infile = open(cart3d_filename, 'rb')
            (model.nPoints, model.nElements) = self.read_header_binary()
            points = model.read_points_binary(self.nPoints)
            elements = model.read_elements_binary(self.nElements)
            regions = model.read_regions_binary(self.nElements)
            #loads = {}
        else:
            model.infile = open(cart3d_filename, 'r')
            model.read_header_ascii()
            points = model.read_points_ascii(bypass=True)
            elements = model.read_elements_ascii(bypass=True)
            regions = model.read_regions_ascii(bypass=True)
            loads = model.read_results_ascii(0, model.infile, result_names=result_names)
        self.load_cart3d_geometry(cart3d_filename, dirname)


    def _fill_bedge_case(self, bedge_filename, cases, ID, nnodes, nelements, model):
        print("is_centroidal=%s isNodal=%s" % (self.is_centroidal, self.is_nodal))
        assert self.is_centroidal != self.is_nodal

        base, ext = os.path.splitext(bedge_filename)
        assert ext == '.bedge', bedge_filename

        tag_filename = base + '.tags'

        cases_new = []
        has_tag_data = False
        results_form = []
        if self.is_centroidal:
            geometry_form = [
                #('Region', 0, []),
                ('ElementID', 0, []),
                ('CurveID', 1, []),
                ('SubcurveID', 2, []),
                ('GridBC', 3, []),
                #('TurnAngle', 3, [])
            ]

            eids = arange(1, nelements + 1)

            if 0:
                surf_ids = element_props[:, 0]
                recon_flags = element_props[:, 1]
                grid_bcs = element_props[:, 2]

            cases[(ID, 0, 'ElementID',  1, 'centroid', '%i')] = eids
            cases[(ID, 1, 'CurveID',    1, 'centroid', '%i')] = model.curves
            cases[(ID, 2, 'SubcurveID', 1, 'centroid', '%i')] = model.subcurves
            cases[(ID, 3, 'GridBC',     1, 'centroid', '%i')] = model.grid_bcs

            if hasattr(model, 'turn_angle'):
                gf = ('TurnAngle', 4, [])
                geometry_form.append(gf)
                cases[(ID, 4, 'TurnAngle', 1, 'centroid', '%.2f')] = degrees(abs(model.turn_angle))


        elif self.is_nodal:
            geometry_form = [
                ('NodeID', 0, []),
                #('normSpacing', 1, []),
                #('BL_thick', 2, []),
            ]
            nids = arange(1, nnodes+1)
            #norm_spacing = model.node_props[:, 0]
            #bl_thickness = model.node_props[:, 1]
            cases[(ID, 0, 'NodeID',      1, 'node', '%i')] = nids
            #cases[(ID, 1, 'normSpacing', 1, 'node', '%.3e')] = norm_spacing
            #cases[(ID, 2, 'BL_thick',    1, 'node', '%.3e')] = bl_thickness

        form = [
            ('Geometry', None, geometry_form),
        ]
        if has_tag_data:
            form.append(('Tag Data', None, tag_form),)

        results_form = []
        if len(results_form):
            form.append(('Results', None, results_form))
        return form, cases
