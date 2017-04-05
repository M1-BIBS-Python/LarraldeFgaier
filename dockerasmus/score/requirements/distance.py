# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals


try:
    from scipy.spatial.distance import cdist as _dist
except ImportError:
    from ...utils.matrices import distance as _dist


def distance(protein1, protein2):
    """The euclidean distances of ``protein1`` atoms to ``protein2`` atoms.

    .. hint::
        If available, the distance matrix will be computed using
        `scipy.spatial.distance.cdist`. *Speedup is substantial !*
    """
    return _dist(
        protein1.atom_positions(),
        protein2.atom_positions(),
    )
