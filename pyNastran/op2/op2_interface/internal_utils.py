from __future__ import annotations
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.f06.f06_formatting import get_key0
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from pyNastran.op2.op2 import OP2


def get_result_length(op2: OP2, res_types: list[dict], res_key: str) -> int:
    """
    gets the length of the output data so we can line up:

      RealCRodStrain  - CROD
      RealCTubeStrain - CTUBE

    """
    res_length = 0
    for res_type in res_types:
        if not res_type:
            continue
        key0 = next(iter(res_type))
        if not isinstance(key0, integer_types) and not isinstance(res_key, integer_types):
            if not type(key0) == type(res_key):
                msg = (
                    'bad compression check...\n'
                    'keys0=%s type(key0)=%s\n'
                    'res_key=%s type(res_key)=%s' % (
                        key0, type(key0), res_key, type(res_key))
                )
                raise RuntimeError(msg)

        #print('res_type.keys()=%s' % res_type.keys())
        # res_key_list = res_key[:-1] + [res_key[-1]]
        # res_key = tuple(res_key_list)

        if res_key in res_type:
            # the res_key is
            result = res_type[res_key]
            class_name = result.__class__.__name__
            res_length = max(len(class_name), res_length)
            #print('continue')
            #break
            continue
        elif len(res_type) != 0:
            #print('  not valid')
            # get the 0th key in the dictionary, where key0 is arbitrary
            key0 = get_key0(res_type)
            #print('  key0 = ', key0)

            # extract displacement[0]
            result = res_type[key0]

            # get the class name
            class_name = result.__class__.__name__
            res_length = max(len(class_name), res_length)

            if not op2.warn_on_op2_missed_table:
                print(' %s - results not found...key=%s' % (class_name, res_key))
        else:  # empty result
            #print('else')
            pass
    #print('res_length =', res_length)
    return res_length
