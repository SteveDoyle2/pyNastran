import sys
import glob

from pyNastran.converters.dev.tecplot.tecplot_reader import merge_tecplot_files


def run():
    files = sys.argv[1]
    tecplot_filenames = glob.glob(files)

    tecplot_filename_out = sys.argv[2]
    merge_tecplot_files(tecplot_filenames, tecplot_filename_out=tecplot_filename_out)


if __name__ == '__main__':
    run()