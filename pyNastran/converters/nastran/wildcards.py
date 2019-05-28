try:
    import h5py
    IS_H5PY = True
except ImportError:
    IS_H5PY = False

bdf_h5 = ''
if IS_H5PY:
    bdf_h5 = '*.h5; '
GEOM_METHODS_BDF = ('Nastran Geometry - BDF (*.bdf; *.dat; *.nas; *.ecd; '
                    '*.op2; *.pch; %s*.obj)' % bdf_h5)
GEOM_BDF_SAVE = ('Nastran Geometry - BDF (*.bdf; *.dat; *.nas; *.ecd; '
                 '*.pch; %s*.obj)' % bdf_h5)
