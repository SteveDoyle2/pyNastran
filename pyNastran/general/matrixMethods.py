from math import sqrt


def gauss(nGauss):
    """
    http://en.wikipedia.org/wiki/Gaussian_quadrature
    """
    if nGauss == 1:
        P = [0.]
        W = [2]
    elif nGauss == 2:

        p = 1. / sqrt(3)
        P = [w, -w]
    elif nGauss == 3:
        p = sqrt(3 / 5.)
        P = [p, 0., -p]
    elif nGauss == 4:
        p1 = (3 - 2. * sqrt(6 / 5)) / 7.
        p2 = (3 + 2. * sqrt(6 / 5)) / 7.
        P = [-p1, p1, -p2, p2]
        w1 = (18 + sqrt(30)) / 36.
        w2 = (18 - sqrt(30)) / 36.
        W = [w1, w1, w2, w2]
    elif nGauss == 5:
        p1 = 1 / 3. * sqrt(5 - 2 * sqrt(10. / 7.))
        p2 = 1 / 3. * sqrt(5 + 2 * sqrt(10. / 7.))
        w1 = (322 + 13 * sqrt(70)) / 900.
        w2 = (322 - 13 * sqrt(70)) / 900.
        P = [0., -p1, p1, -p2, p2]
        W = [128. / 225., w1, w1, w2, w2]
    else:
        raise NotImplementedError('the code only supports up to 5 quadrature '
                                  'points')
    return (P, W)
