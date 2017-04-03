# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals


try:
    from scipy.spatial.distance import cdist as _dist
except:
    from ...utils.matrices import distance as _dist


def distance(protein1, protein2):
    return _dist(
        protein1.atom_positions(),
        protein2.atom_positions(),
    )
