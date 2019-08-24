import os

from cpylog import get_logger2
from pyNastran.converters.tecplot.tecplot import read_tecplot

def main():
    log = get_logger2(debug=True)

    dirnames = ['ascii', 'binary']
    filenames = [

        #'ascii/3dgeom.dat', #  good; multi-zone, geometry
        #'ascii/block_febrick_3d.dat', # 3d unstructured block; good
        #'ascii/block_fetet_3d.dat', # bad; no decimal values
        #'ascii/channel.dat', # 2d structured point; good
        #'ascii/cylinder_slice.dat', # 3d structured point; good
        #'ascii/cylindrical.dat',  # 3d structured empty lines; good
        #'ascii/ell.dat', # 2d; good
        #'ascii/humanoid_quad.dat', # good
        #'ascii/humanoid_tri.dat', # good

        #'ascii/movie.dat',  # csv -> bad
        'ascii/multzn2d.dat',  #  2d structured; good
        #'ascii/plane_slice.dat',  # 2d structured multi-line; good
        #'ascii/point_febrick_3d_02.dat',  # difficult header and funny write bug; bad
        #'ascii/point_fequad_2d.dat',  # 2d; good
        #'ascii/point_fetet_3d.dat',  # good
        #'ascii/point_fetri_2d_01.dat',  # good
        #'ascii/point_fetri_2d_02.dat',  # good
        #'ascii/point_fetri_2d_03.dat',  # good

        #'ascii/simp3dbk.dat',  # 3d structured block - bad
        #'ascii/simp3dpt.dat', #  good
        #'ascii/simpscat.dat', #  bad -> text
        #'ascii/simpxy.dat',  # no xyz; it's a plot -> bad
        #'ascii/simpxy2.dat',  # no xyz; it's a plot -> bad
        #'ascii/tiny.dat',  # good
    ]
    run_filenames(filenames, log)
    return
    dirnames = ['binary']
    for dirname in dirnames:
        fnames = [os.path.join(dirname, fname) for fname in os.listdir(dirname)
                  if not fname.endswith('.png')]
        run_filenames(fnames, log)

def run_filenames(fnames, log):
    for fname in fnames:
        #print(fname)
        #continue
        try:
            model = read_tecplot(fname, log=log)
            log.info('read %r' % fname)
        except Exception as error:
            log.warning('failed reading %r' % fname)
            log.error(error)
            print('')
            raise
            continue
        print(model)
        model.write_tecplot('junk.plt', res_types=None, adjust_nids=True)

if __name__ == '__main__':  # pragma: no cover
    main()

