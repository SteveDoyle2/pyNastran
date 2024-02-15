import numpy as np

def volume_ctetra(n1: np.ndarray, n2: np.ndarray,
                  n3: np.ndarray, n4: np.ndarray) -> np.ndarray:
    #volume = -dot(n1 - n4, cross(n2 - n4, n3 - n4)) / 6.
    n14 = n1 - n4
    n24 = n2 - n4
    n34 = n3 - n4
    n234 = np.cross(n24, n34, axis=1)

    # dot product
    volume = np.einsum("ij, ij->i", n14, n234) / -6.
    return volume

def volume_cpyram(n1: np.ndarray, n2: np.ndarray,
                  n3: np.ndarray, n4: np.ndarray,
                  n5: np.ndarray) -> np.ndarray:
    vpyramid1 = volume_ctetra(n1, n2, n3, n5)
    vpyramid2 = volume_ctetra(n1, n3, n4, n5)
    volume = vpyramid1 + vpyramid2
    return volume

def volume_cpenta(n1: np.ndarray, n2: np.ndarray,
                  n3: np.ndarray, n4: np.ndarray,
                  n5: np.ndarray, n6: np.ndarray) -> np.ndarray:
    nelement = n1.shape[0]
    length1 = np.linalg.norm(n1 - n4, axis=1)
    length2 = np.linalg.norm(n2 - n5, axis=1)
    length3 = np.linalg.norm(n3 - n6, axis=1)
    assert len(length1) == nelement
    avg_length = (length1 + length2 + length3) / 3

    #c1 = (n1 + n2 + n3) / 3.
    #c2 = (n4 + n5 + n6) / 3.
    #avg_length = np.linalg.norm(c1-c2, axis=1)
    #assert len(avg_length) == nelements

    a1 = np.cross(n2-n1, n3-n1, axis=1)
    a2 = np.cross(n5-n4, n6-n4, axis=1)
    assert a1.shape == (nelement, 3)

    # need to divide A1 and A2 by 0.5
    # additional 0.5 factor for the average area
    a_avg = 0.25 * (np.linalg.norm(a1, axis=1) + np.linalg.norm(a2, axis=1))
    assert len(a_avg) == nelement
    #vpyramid1 = volume4_array(n1, n2, n3, n4)
    v_triangular_prism = avg_length * a_avg
    assert len(v_triangular_prism) == nelement
    return v_triangular_prism

def volume_chexa(n1: np.ndarray, n2: np.ndarray,
                 n3: np.ndarray, n4: np.ndarray,
                 n5: np.ndarray, n6: np.ndarray,
                 n7: np.ndarray, n8: np.ndarray) -> np.ndarray:
    nelement = n1.shape[0]
    #vpyramid1 = volume4_array(n1, n2, n3, n5)
    #vpyramid2 = volume4_array(n1, n3, n4, n5)
    #vpyramid = vpyramid1 + vpyramid2

    # https://www.osti.gov/servlets/purl/632793/
    #volume = (
        #det3(x7 - x0, x1 - x0, x3 - x5) +
        #det3(x7 - x0, x4 - x0, x5 - x6) +
        #det3(x7 - x0, x2 - x0, x6 - x3)
    #) / 6.
    #  swap points
    # x2 / x3
    # x6 / x7
    def det3(a, b, c):
        stack = np.dstack([a, b, c])
        d3 = np.linalg.det(stack)
        #except FloatingPointError:
            # recasting it to a float64 array gets rid of the underflow
            #if stack.dtype.name == 'float64':
                #stack2 = stack.astype('float64')
                #d3 = np.linalg.det(stack2)
            #else:
                #raise NotImplementedError(stack.dtype.name)
            #abs_stack = np.abs(stack-stack2)
            #print(f'dstack = {abs_stack.max():.3e}')
        return d3

    #volume = (
        #det3(x6 - x0, x1 - x0, x2 - x5) +
        #det3(x6 - x0, x4 - x0, x5 - x7) +
        #det3(x6 - x0, x3 - x0, x7 - x2)
    #) / 6.
    # add 1
    n71 = n7 - n1
    v1 = det3(n71, n2 - n1, n3 - n6)
    v2 = det3(n71, n5 - n1, n6 - n8)
    v3 = det3(n71, n4 - n1, n8 - n3)
    assert len(v1) == nelement
    #volume = (
        #det3(n71, n2 - n1, n3 - n6) +
        #det3(n71, n5 - n1, n6 - n8) +
        #det3(n71, n4 - n1, n8 - n3)
    #) / 6.
    volume = (v1 + v2 + v3) / 6.
    return volume
