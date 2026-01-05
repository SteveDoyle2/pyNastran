import os
from pathlib import Path
import pyNastran
from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import read_bdf, BDF, DOPTPRM, DVPREL1, DVPREL2, DESVAR
from pyNastran.op2.op2 import read_op2
from pyNastran.utils.nastran_utils import run_nastran
PKG_PATH = Path(pyNastran.__path__[0])


def FullyStressedDesign:
    def __init__(self, bdf_filename: PathLike, mode: str='nx'):
        model = BDF(mode=mode)
        self.model = model
        self.model.read_bdf(bdf_filename)

    def get_updated_filenames(base: str, n0: int=0) -> tuple[str, str, str]:
        """bdf f06 op2"""
        n = n0
        while 1:
            bdf_filename = f'{base}.{n}.bdf'
            f06_filename = f'{base}.{n}.f06'
            op2_filename = f'{base}.{n}.op2'
            if (os.path.exists(f06_filename) or
                os.path.exists(f06_filename) or
                os.path.exists(f06_filename)):
                n += 1
                continue
            break
        return bdf_filename, f06_filename, op2_filename

    def setup_optimization(self):
        model: BDF = self.model
        desglb_id = 0
        for subcase_id, subcase in model.subcases.items():
            if subcase_id == 0:
                continue
            assert subcase_id == 1, subcase_id
            dessub_id = subcase.get_int_parameter('DESSUB')
            desobj_id = subcase.get_int_parameter('DESOBJ')
            #desglb_id = subcase.get_int_parameter('DESGLB')

        # case control
        doptprm: DOPTPRM = model.doptprm
        ncycles = doptprm.params['FSDMAX']
        dvprels = model.dvprels[dessub_id]
        for dvprel in dvprels:
            if dvprel.type == 'DVPREL1':
                dvprel = cast(DVPREL1, dvprel)
                asdf
            else:
                assert dvprel.type == 'DVPREL2', dvprel

            if dvprel.
        'CONV1' : 0.0001,  # NX 2019.2
        'CONV2' : 1e-20,
        'CONVDV' : 0.001,  # 0.0001 for topology optimization

        # DESSUB Selects the design constraints to be used
        # in a design optimization task for the current subcase.


    def run(self):
        self.supress_outputs()
        self.setup_optimization()
        model: BDF = self.model
        bdf_filename_base = model.bdf_filename
        base = os.path.splitext(bdf_filename_base)[0]
        bdf_filename = os.path.basename(bdf_filename)
        op2_filename = base + '.op2'

        out = self.get_updated_filenames(base, n0=0)
        bdf_filenamei, f06_filenamei, op2_filenamei = out
        run_nastran(bdf_filename, nastran_cmd=nastran_cmd, keywords=[], run=True)
        op2_model = read_op2(op2_filename)


    def suppress_outputs(self, word: str='PLOT', overwrite: bool=False):
        model: BDF = self.model
        options_to_remove = ['PRINT', 'PUNCH']
        keys = ['STRESS', 'STRAIN', 'FORCE', 'GPFORCE', 'DISPLACEMENT']
        for subcase in model.subcases.items():
            for key in subcase.keys():
                if key not in keys:
                    continue
                value, options = subcase.get_parameter(key)
                if overwrite:
                    for option in options_to_remove:
                        if option in options:
                            options_to_remove.remove(option)
                else:
                    if word not in options:
                        options.append(word)

def main():
    bdf_filename = PKG_PATH / '..' / 'models' / 'bwb' / 'bwb_saero.bdf'
    fsd = FullyStressedDesign(bdf_filename, mode=)
    self.
    fsd.run()


if __name__ == '__main__':  # pragma: no cover
    main()
