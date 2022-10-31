import numpy as np

def get_minus1_start_end(ints: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    iminus1 = np.where(ints == -1)[0]
    istart = [0] + list(iminus1[:-1] + 1)
    iend = iminus1
    return istart, iend
