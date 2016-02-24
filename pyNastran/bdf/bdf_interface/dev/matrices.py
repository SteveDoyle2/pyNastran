from six import iteritems


def AJJ(model, xyz=None):
    if xyz is None:
        xyz = {}
        for nid, node in iteritems(model.nodes):
            xyz[nid] = node.get_position()
    for caero_id, caero in iteritems(model.caeros):
        c = caero.get_centroids()

    for spline_id, spline in iteritems(model.splines):
        s = spline.spline_nodes


    return Ajj