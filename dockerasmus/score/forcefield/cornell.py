# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import warnings
import functools

from ..base import ScoringFunction
from ...utils.matrices import compose, distance


__all__ = ["score"]


class _CornellScoringFunction(ScoringFunction):
    """Matrix computation of the VDW term of Cornell's force field

    Arguments:
        protein1 (Protein): the reference protein (receptor)
        protein2 (Protein): the moving protein (ligand)

    Keyword Arguments:
        diel (float): the value of the dielectric constant for the medium.
            Default value uses results obtained after simulations in
            `Cossins et al. <https://dx.doi.org/10.1371/journal.pcbi.1002066>`_.
            for the cytoplasm of E. Coli. [default: 65.0]

    Exemple:
        >>> cornell.score(barnase, barstar)
        -22.4...

    .. note::
        Available backends: `theano`, `numpy`

    .. seealso::
        Source concept of the forcefield described in
        `Cornell et al. <https://dx.doi.org/10.1021/ja00124a002>`_
    """
    _backends = ['theano', 'numpy']

    def _setup_theano(self, theano):

        ### Dielectric constant
        diel = theano.tensor.dscalar('diel')

        ### Epsilon matrix from protein vectors
        v_eps1, v_eps2 = theano.tensor.dvectors('v_eps1', 'v_eps2')
        mx_eps = theano.tensor.sqrt(theano.tensor.outer(v_eps1, v_eps2))

        ### Radius matrix from protein vectors
        v_rad1, v_rad2 = theano.tensor.dvectors('v_rad1', 'v_rad2')
        # z = theano.tensor.zeros((v_rad2.shape[0], v_rad1.shape[0]))
        # mx_rad = theano.tensor.transpose(v_rad1+z)+v_rad2
        mx_rad = theano.tensor.transpose(theano.tensor.repeat(
            theano.tensor.reshape(
                v_rad1, (1, v_rad1.shape[0])
            ), v_rad2.shape[0], axis=0,
        )) + v_rad2

        ### Charge matrix from protein vectors
        v_q1, v_q2 = theano.tensor.dvectors('v_q1', 'v_q2')
        mx_q = theano.tensor.outer(v_q1, v_q2)

        ### Euclidian distance matrix from protein matrices
        mx_pos1, mx_pos2 = theano.tensor.dmatrices('mx_pos1', 'mx_pos2')
        # Arrange each vector for pairwise euclidian norm computation
        mx_arr_pos1 = theano.tensor.repeat(mx_pos1, mx_pos2.shape[0], axis=0)
        mx_arr_pos2 = theano.tensor.tile(mx_pos2, (mx_pos1.shape[0], 1))
        # Calculate vector of euclidian distance
        v_d = theano.tensor.sqrt(theano.tensor.sum((mx_arr_pos1-mx_arr_pos2)**2, axis=1))
        # Rearrange vector as a matrix where mx_d[i,j] is the distance between
        # the i-th atom of protein 1 and the j-th atom of protein 2
        mx_d = theano.tensor.reshape(v_d, (mx_pos1.shape[0], mx_pos2.shape[0]))

        ### Van der Waals constants
        mx_A = mx_eps * mx_rad**12
        mx_B = 2 * mx_eps * mx_rad**6

        ### Final scoring function
        self._score = theano.function(
            [v_eps1, v_eps2, v_rad1, v_rad2, v_q1, v_q2, mx_pos1, mx_pos2, diel],
            theano.tensor.sum(theano.tensor.triu(
                mx_A/(mx_d**12)-mx_B/(mx_d**6)+mx_q/(diel*mx_d), k=1,
            ))
        )

    def _setup_numpy(self, numpy):

        def _score(v_eps1, v_eps2, v_rad1, v_rad2, v_q1, v_q2, mx_pos1, mx_pos2, diel):
            # Matrix of Van der Waals epsilon
            mx_epsilon = compose(lambda e1, e2: numpy.sqrt(e1*e2), v_eps1, v_eps2)

            # Matrix of Van der Waals optimal radius for each aminoacid,
            # to the power of 6
            mx_radius_6 = compose(lambda r1, r2: r1+r2, v_rad1, v_rad2)**6

            # Matrices A and B from Cornell
            mx_B = 2 * mx_epsilon * mx_radius_6
            mx_A = 0.5 * mx_B * mx_radius_6

            # Matrix of aminoacid-wise distance
            mx_distance = distance(mx_pos1, mx_pos2)
            mx_distance_6 = mx_distance**6

            # Matrix of aminoacid-wise charge
            mx_q = compose(lambda q1, q2: q1*q2, v_q1, v_q2)

            # Final matrix (using numpy.triu(k=1, ...) means we're only
            # summing for i<j, where i is the receptor and j the ligand).
            return numpy.sum(numpy.triu(k=1, m=\
                mx_A/(mx_distance_6**2)-mx_B/(mx_distance_6)+mx_q/(80*mx_distance))
            )
        self._score = _score

    def __call__(self, receptor, ligand, diel=65.0):
        return float(self._score(
            receptor.atom_epsilon(),
            ligand.atom_epsilon(),
            receptor.atom_radius(),
            ligand.atom_radius(),
            receptor.atom_charges(),
            ligand.atom_charges(),
            #mx_distance,
            receptor.atom_positions(),
            ligand.atom_positions(),
            diel
        ))


score = _CornellScoringFunction()
