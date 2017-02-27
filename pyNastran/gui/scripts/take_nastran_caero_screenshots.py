from pyNastran.bdf.bdf import read_bdf

# find the list of CAERO cards
bdf_filename = self.infilename
model = read_bdf(bdf_filename, debug=False)
caeros = model.caeros.keys()
#print caeros
del model

# make a bdf for each CAERO card
for caeroi in caeros:
    model2 = read_bdf(bdf_filename, debug=False)

    for eid in caeros:
        if eid != caeroi:
            del model2.caeros[eid]
    bdf_filename2 = 'caero_%s.bdf' % caeroi
    model2.write_bdf(bdf_filename2)

# load each model and take a picture
for caeroi in caeros:
    bdf_filename = 'caero_%i.bdf' % caeroi
    self.on_load_geometry(infile_name=bdf_filename, geometry_format='nastran')
    self.zoom(1.6)
    png_name = 'caero_%s.png' %(caeroi)
    self.on_take_screenshot(png_name)
