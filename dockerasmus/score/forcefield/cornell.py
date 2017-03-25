# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import numpy

from ...utils.matrices import compose, distance


def score(protein1, protein2):
    """Matrix computation of the VDW term of Cornell's force field
    """

    atoms_1 = list(protein1.iteratoms())
    atoms_2 = list(protein2.iteratoms())

    mx_epsilon = compose(
        lambda e1, e2: numpy.sqrt(e1*e2),
        numpy.array([a.epsilon for a in atoms_1]),
        numpy.array([a.epsilon for a in atoms_2]),
    )

    mx_radius = compose(
        lambda r1, r2: r1+r2,
        numpy.array([a.radius for a in atoms_1]),
        numpy.array([a.radius for a in atoms_2]),
    )

    mx_A = mx_epsilon * mx_radius**12
    mx_B = 2 * mx_epsilon * mx_radius**6

    mx_distance = distance(
        numpy.array([a.pos for a in atoms_1]),
        numpy.array([a.pos for a in atoms_2]),
    )

    mx_q = compose(
        lambda q1, q2: q1*q2,
        numpy.array([a.charge for a in atoms_1]),
        numpy.array([a.charge for a in atoms_2]),
    )

    return numpy.sum(numpy.triu(k=1, m=\
        mx_A/(mx_distance**12)-mx_B/(mx_distance**6))+mx_q/(mx_epsilon*mx_distance)
    )
