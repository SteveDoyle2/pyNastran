from cpylog import SimpleLogger
from pyNastran.op2.op2 import OP2
from pyNastran.utils import PathLike

class ResultsEnvelope:
    def __init__(self, log=None):
        log = SimpleLogger(level='debug')
        self.op2_filenames = []
        self.result_names = None
        self.model = OP2(log=log)
        self.nfiles_max = 1
        self.log = log

    def add_op2_filenames(self, op2_filenames: list[PathLike]):
        self.op2_filenames = op2_filenames

    def set_result_names(self, result_names: list[str]):
        self.result_names = result_names

    def run(self) -> None:
        log = self.log
        assert len(self.op2_filenames), self.op2_filenames
        assert self.result_names is not None, self.result_names
        from collections import defaultdict
        gres_dict = defaultdict(list)
        for op2_filename in self.op2_filenames:
            model = OP2(log=log)
            model.read_op2(op2_filename)
            for result_name in self.result_names:
                res_dict_total = self.model.get_result(result_name)
                res_dict = model.get_result(result_name)
                # quick skip
                # if len(res_dict_total) == 0:
                #     model.set_result(result_name, res_dict)
                #     continue

                assert self.nfiles_max == 1, self.nfiles_max
                if isinstance(res_dict, dict):
                    for subcase, resi in res_dict.items():
                        key = (result_name, subcase)
                        gres_dict[key].append(resi)
                # elif isinstance(resi, TableArray):
                #     new_data = []
                else:
                    raise RuntimeError(result_name)
        keys = list(gres_dict.keys())
        print(keys)

def main():
    import pyNastran
    from pathlib import Path
    PKG_PATH = Path(pyNastran.__path__[0])
    assert PKG_PATH.exists()
    op2_filename = PKG_PATH / '..' / 'models' / 'solid_bending' /  'solid_bending.op2'
    env = ResultsEnvelope(log=None)
    env.set_result_names(['stress.ctetra_stress'])
    env.add_op2_filenames([op2_filename, op2_filename])
    env.run()

if __name__ == '__main__':
    main()
