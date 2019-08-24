from pyNastran.converters.dev.lsdyna.lsdyna import read_lsdyna

if __name__ == '__main__':
    # https://www.nhtsa.gov/DOT/NHTSA/NVS/Biomechanics%20&%20Trauma/Thor-Lx%20Finite%20Element%20Model/thor-no-skin-right-non-sae.key
    model = read_lsdyna('thor_no_skin_right_non_sae.key')
