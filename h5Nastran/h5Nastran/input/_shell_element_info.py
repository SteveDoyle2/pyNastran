from __future__ import print_function, absolute_import

import numpy as np
import tables
from numpy import cross

from pyNastran.op2.data_in_material_coord import is_mcid, angle2vec, check_theta, calc_imat


shell_element_info_dtype = np.dtype([
    ('EID', '<i8', ()),
    ('CENTER', '<f8', (3,)),
    ('THETA_RAD', '<f8', ())
])

shell_element_info_format = tables.descr_from_dtype(shell_element_info_dtype)[0]


class QuadShell(object):
    def get_shell_info(self, bdf, cards):
        elms = cards.get(self.__class__.__name__, [])
        elm_count = len(elms)
        result = np.empty(elm_count, dtype=shell_element_info_dtype)
    
        if len(result) == 0:
            return result
    
        _eids = result['EID']
        _centers = result['CENTER']
        _theta_rad = result['THETA_RAD']
    
        eids = np.array(sorted(elms.keys()))
        elems = np.array([elms[eid] for eid in eids])
        mcid = np.array([is_mcid(e) for e in elems])
        elemsmcid = elems[mcid]
        elemstheta = elems[~mcid]
    
        thetadeg = np.zeros(elems.shape)
        thetadeg[~mcid] = np.array([check_theta(e) for e in elemstheta])
        thetarad = np.deg2rad(thetadeg)
    
        thisquad = np.array([True] * len(elemstheta), dtype=bool)
    
        if len(thisquad) > 0:
            quadelems = elemstheta[thisquad]
            corner = np.array([e.get_node_positions() for e in quadelems])
            g1 = corner[:, 0, :]
            g2 = corner[:, 1, :]
            g3 = corner[:, 2, :]
            g4 = corner[:, 3, :]
            betarad = angle2vec(g3 - g1, g2 - g1)
            gammarad = angle2vec(g4 - g2, g1 - g2)
            alpharad = (betarad + gammarad) / 2.
            tmp = thetarad[~mcid]
            tmp[thisquad] -= betarad
            tmp[thisquad] += alpharad
            thetarad[~mcid] = tmp
    
        thisquad = np.array([True] * len(elemsmcid), dtype=bool)
    
        if len(thisquad) > 0:
            quadelems = elemsmcid[thisquad]
            corner = np.array([e.get_node_positions() for e in quadelems])
            g1 = corner[:, 0, :]
            g2 = corner[:, 1, :]
            g3 = corner[:, 2, :]
            g4 = corner[:, 3, :]
            normals = np.array([e.Normal() for e in quadelems])
            csysi = np.array([bdf.coords[e.theta_mcid].i for e in quadelems])
            imat = calc_imat(normals, csysi)
            tmp = thetarad[mcid]
            tmp[thisquad] = angle2vec(g2 - g1, imat)
            # getting sign of THETA
            check_normal = cross(g2 - g1, imat)
            tmp[thisquad] *= np.sign((check_normal * normals).sum(axis=1))
            betarad = angle2vec(g3 - g1, g2 - g1)
            gammarad = angle2vec(g4 - g2, g1 - g2)
            alpharad = (betarad + gammarad) / 2.
            tmp[thisquad] -= betarad
            tmp[thisquad] += alpharad
            thetarad[mcid] = tmp
    
        _eids[:] = eids
        _centers[:] = [elm.Centroid() for elm in elems]
        _theta_rad[:] = thetarad
    
        return result


class TriaShell(object):
    def get_shell_info(self, bdf, cards):
        elms = cards.get(self.__class__.__name__, [])
        elm_count = len(elms)
        result = np.empty(elm_count, dtype=shell_element_info_dtype)

        if len(result) == 0:
            return result

        _eids = result['EID']
        _centers = result['CENTER']
        _theta_rad = result['THETA_RAD']

        eids = np.array(sorted(elms.keys()))
        elems = np.array([elms[eid] for eid in eids])
        mcid = np.array([is_mcid(e) for e in elems])
        elemsmcid = elems[mcid]
        elemstheta = elems[~mcid]

        thetadeg = np.zeros(elems.shape)
        thetadeg[~mcid] = np.array([check_theta(e) for e in elemstheta])
        thetarad = np.deg2rad(thetadeg)

        thistria = np.array([True] * len(elemsmcid), dtype=bool)

        if len(thistria) > 0:
            triaelems = elemsmcid[thistria]
            corner = np.array([e.get_node_positions() for e in triaelems])
            g1 = corner[:, 0, :]
            g2 = corner[:, 1, :]
            # g3 = corner[:, 2, :]
            normals = np.array([e.Normal() for e in triaelems])
            csysi = np.array([bdf.coords[e.theta_mcid].i for e in triaelems])
            imat = calc_imat(normals, csysi)
            tmp = thetarad[mcid]
            tmp[thistria] = angle2vec(g2 - g1, imat)
            # getting sign of THETA
            check_normal = cross(g2 - g1, imat)
            tmp[thistria] *= np.sign((check_normal * normals).sum(axis=1))
            thetarad[mcid] = tmp

        _eids[:] = eids
        _centers[:] = [elm.Centroid() for elm in elems]
        _theta_rad[:] = thetarad

        return result
