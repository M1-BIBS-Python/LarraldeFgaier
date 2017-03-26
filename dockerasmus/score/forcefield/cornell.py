# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import numpy

from ...utils.matrices import compose, distance


def score(protein1, protein2):
    """Matrix computation of the VDW term of Cornell's force field
    """

    mx_epsilon = compose(
        lambda e1, e2: numpy.sqrt(e1*e2),
        numpy.array([a.epsilon for a in protein1.iteratoms()]),
        numpy.array([a.epsilon for a in protein2.iteratoms()]),
    )

    mx_radius_6 = compose(
        lambda r1, r2: r1+r2,
        numpy.array([a.radius for a in protein1.iteratoms()]),
        numpy.array([a.radius for a in protein2.iteratoms()]),
    )**6

    #mx_A = mx_epsilon * mx_radius_12
    #mx_B = 2 * mx_epsilon * mx_radius_6
    mx_B = 2 * mx_epsilon * mx_radius_6
    mx_A = 0.5 * mx_B * mx_radius_6

    mx_distance = distance(
        numpy.array([a.pos for a in protein1.iteratoms()]),
        numpy.array([a.pos for a in protein2.iteratoms()]),
    )
    mx_distance_6 = mx_distance**6

    mx_q = compose(
        lambda q1, q2: q1*q2,
        numpy.array([a.charge for a in protein1.iteratoms()]),
        numpy.array([a.charge for a in protein2.iteratoms()]),
    )

    return numpy.sum(numpy.triu(k=1, m=\
        mx_A/(mx_distance_6**2)-mx_B/(mx_distance_6)+mx_q/(80*mx_distance))
    )
