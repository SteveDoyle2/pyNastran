"""Interface to various numpy utilities"""
import numpy as np

integer_types = (int, np.int32, np.int64)
integer_string_types = (int, np.int32, np.int64, bytes, str)
integer_float_types = (int, np.int32, np.int64, float, np.float32)
float_types = (float, np.float32, np.float64)
bytes_type = (bytes, np.bytes_)
