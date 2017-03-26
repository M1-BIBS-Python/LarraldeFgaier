# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import numpy

from ...utils.matrices import compose, distance


def score(protein1, protein2, diel=65.0):
    """Matrix computation of the VDW term of Cornell's force field

    Arguments:
        protein1 (Protein): the reference protein (receptor)
        protein2 (Protein): the moving protein (ligand)

    Keyword Arguments:
        diel (float): the value of the dielectric constant for the medium.
            Default value uses results obtained after simulations in
            `Cossins et al. <https://dx.doi.org/10.1371/journal.pcbi.1002066>`_.
            for the cytoplasm of E. Coli. [default: 65.0]

    .. seealso::
        Source concept of the forcefield described in
        `Cornell et al. <https://dx.doi.org/10.1021/ja00124a002>`_
    """

    # Matrix of Van der Waals
    mx_epsilon = compose(
        lambda e1, e2: numpy.sqrt(e1*e2),
        numpy.array([a.epsilon for a in protein1.iteratoms()]),
        numpy.array([a.epsilon for a in protein2.iteratoms()]),
    )

    # Matrix of Van der Waals optimal radius for each aminoacid,
    # to the power of 6
    mx_radius_6 = compose(
        lambda r1, r2: r1+r2,
        numpy.array([a.radius for a in protein1.iteratoms()]),
        numpy.array([a.radius for a in protein2.iteratoms()]),
    )**6

    # Matrices A and B from Cornell
    mx_B = 2 * mx_epsilon * mx_radius_6
    mx_A = 0.5 * mx_B * mx_radius_6

    # Matrix of aminoacid-wise distance
    mx_distance = distance(
        numpy.array([a.pos for a in protein1.iteratoms()]),
        numpy.array([a.pos for a in protein2.iteratoms()]),
    )
    mx_distance_6 = mx_distance**6

    # Matrix of aminoacid-wise charge
    mx_q = compose(
        lambda q1, q2: q1*q2,
        numpy.array([a.charge for a in protein1.iteratoms()]),
        numpy.array([a.charge for a in protein2.iteratoms()]),
    )

    # Final matrix (using numpy.triu(k=1, ...) means we're only
    # summing for i<j, where i is the receptor and j the ligand).
    return numpy.sum(numpy.triu(k=1, m=\
        mx_A/(mx_distance_6**2)-mx_B/(mx_distance_6)+mx_q/(80*mx_distance))
    )
