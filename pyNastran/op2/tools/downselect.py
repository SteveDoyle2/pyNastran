from cpylog.test_log import PKG_PATH
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2


def downselect(
        bdf_filename: BDF,
        op2_filename: OP2,
        eids=None, nids=None,
        percent_eids_target: float=0.90):
    pass

def main():
    import pyNastran
    from pathlib import Path
    PKG_PATH = Path(pyNastran.__file__[0]).parent
    model_path = PKG_PATH / 'models'
    bdf_filename = model_path / 'model.bdf'
    op2_filename = model_path / 'model.op2'
    downselect(bdf_filename, op2_filename)

if __name__ == __main():
    main()
