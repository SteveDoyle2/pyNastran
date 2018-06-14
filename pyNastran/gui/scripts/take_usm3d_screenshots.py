import os
from pyNastran.converters.usm3d.time_accurate_results import get_flo_files_from_n, get_n_list


# disable logging
info = self.show_info
debug = self.show_debug
command = self.show_command

self.show_info = False
self.show_debug = False
self.show_command = False

#===========================
# this could be dynamic...
model_name = 'bay'

dirname = os.getcwd()
n_list = get_n_list(dirname, model_name)#[:20]

#n_list = [65495, 65475]
# take pictures every N steps
N = 440
nlist = [i for i in range(max(n_list)) if i%N==0 ]
flo_filenames = get_flo_files_from_n(dirname, model_name, n_list, include_dirname_in_path=True)


self.on_load_geometry(infile_name=model_name+'.cogsg', geometry_format='usm3d')
self.zoom(1.6)

for ni, flo_filename in zip(n_list, flo_filenames):
    self.on_load_results(out_filename=flo_filename)
    png_name = 'Cp_%s_%s.png' %(model_name, ni)
    self.on_take_screenshot(png_name)

#===========================
# reenable the logging
self.show_info = info
self.show_debug = debug
self.show_command = command
