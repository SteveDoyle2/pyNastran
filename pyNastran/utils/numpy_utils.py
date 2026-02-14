"""Interface to various numpy utilities"""
import numpy as np

integer_types = (int, np.int32, np.int64)
integer_string_types = (int, np.int32, np.int64, bytes, str)
integer_float_types = (int, np.int32, np.int64, float, np.float32)
float_types = (float, np.float32, np.float64)
bytes_type = (bytes, np.bytes_)

def zip_strict(*arrays):
   lengths = [len(array) for array in arrays]
   assert min(lengths) == max(lengths), f'lengths={lengths} should be the same'
   assert len(lengths) > 0, lengths
   return zip(*arrays)

def empty_array(shape, dtype: str, default_int: int=-1) -> np.ndarray:
   """creates a null int/float array"""
   if dtype in {'int32', 'int64'}:
      out = np.full(shape, default_int, dtype=dtype)
   else:
      out = np.full(shape, np.nan, dtype=dtype)
   return out

def cast_ints(ints: list[int], dtype: str='int32') -> np.ndarray:
   if dtype == 'int32':
      try:
         int_array = np.array(ints, dtype=dtype)
      except OverflowError:
         int_array = np.array(ints, dtype='int64')
   elif dtype == 'int64':
      int_array = np.array(ints, dtype=dtype)
   else:
      raise RuntimeError(dtype)
   return int_array
