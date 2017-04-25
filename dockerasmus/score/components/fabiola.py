# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import division

import math

from .base import BaseComponent


class Fabiola(BaseComponent):
    """A scoring component modeling hydrogen bonds.

    Reference:
        `Fabiola, F., Bertram, R., Korostelev, A. & Chapman, M. S.
        "An improved hydrogen bond potential: Impact on medium resolution
        protein structures". Protein Sci 11, 1415–1423 (2002).
        <https://dx.doi.org/10.1110%2Fps.4890102>`_
    """

    backends = ["theano", "numpy"]

    ## FIXME: mx_arccos returns angle+pi (that's why numpy.pi is
    ## substracted), something is wrong with algorithm (since same
    ## issue occurs in both backend implementations)

    def _setup_theano(self, theano):
        from ...utils.tensors import normalized, distance

        ### Parameters
        r_null = theano.tensor.dscalar('r_null')
        sigma = r_null * theano.tensor.sqrt(2/3)
        theta_low, theta_high = theano.tensor.dscalars('theta_low', 'theta_high')
        m = theano.tensor.dscalar('m')

        # Matrice des coordonnes des vecteurs O->C, avec les lignes
        # repetees le nombre de fois necessaires
        mx_pos_c, mx_pos_o, mx_pos_n = theano.tensor.dmatrices('mx_pos_c', 'mx_pos_o', 'mx_pos_n')

        mx_vec_o_to_c = theano.tensor.repeat(
            theano.tensor.reshape(
                normalized(mx_pos_c - mx_pos_o), (mx_pos_o.shape[0], 1, 3)
            ), mx_pos_n.shape[0], axis=1,
        )
        # Matrice des coordonnees des vecteurs O->N, avec i,j la distance
        # entre le i-eme O et le j-eme N
        mx_vec_o_to_n = normalized(
            theano.tensor.repeat(mx_pos_o, mx_pos_n.shape[0], axis=0)
            - theano.tensor.tile(mx_pos_n, (mx_pos_o.shape[0], 1))
        ).reshape((mx_pos_o.shape[0], mx_pos_n.shape[0], 3))
        # Produit scalaire interne des coordonnees de chaque vecteur
        # O->N et O->C
        mx_cos = theano.tensor.sum(mx_vec_o_to_c*mx_vec_o_to_n, axis=-1)
        mx_angles = math.pi - theano.tensor.arccos(mx_cos)

        # Matrice des deviations angulaires en prenant a chaque fois
        # le minimum de deviation avec theta_high et theta_low
        mx_theta_rel = theano.tensor.min([
                abs(mx_angles - theano.tensor.deg2rad(theta_high)),
                abs(mx_angles - theano.tensor.deg2rad(theta_low)),
            ], axis=0,
        )

        # Matrice des distances entre le ieme O et le jeme N
        mx_distance = distance(mx_pos_o, mx_pos_n)

        ### Final function
        self._call = theano.function(
            [mx_pos_o, mx_pos_c, mx_pos_n, m, r_null, theta_low, theta_high],
            theano.tensor.sum(
                ((sigma/mx_distance)**6-(sigma/mx_distance)**4)*theano.tensor.cos(mx_theta_rel)**m
            ),
        )

    def _setup_numpy(self, numpy):
        from ...utils.matrices import normalized, distance
        def call(mx_pos_o, mx_pos_c, mx_pos_n, m, r_null, theta_low, theta_high):
            # Parameter
            sigma = r_null * numpy.sqrt(2/3)
            # Matrice des coordonnes des vecteurs O->C, avec les lignes
            # repetees le nombre de fois necessaires
            mx_vec_o_to_c = normalized(mx_pos_c - mx_pos_o).reshape((mx_pos_o.shape[0], 1, 3)).repeat(mx_pos_n.shape[0], axis=1)
            # Matrice des coordonnees des vecteurs O->N, avec i,j la distance
            # entre le i-eme O et le j-eme N
            mx_vec_o_to_n = normalized(
                numpy.repeat(mx_pos_o, mx_pos_n.shape[0], axis=0)
                - numpy.tile(mx_pos_n, (mx_pos_o.shape[0], 1))
            ).reshape((mx_pos_o.shape[0], mx_pos_n.shape[0], 3))

            # Produit scalaire interne des coordonnees de chaque vecteur
            # O->N et O->C
            mx_cos = numpy.sum(mx_vec_o_to_c*mx_vec_o_to_n, axis=-1)
            mx_angles = numpy.pi - numpy.arccos(mx_cos)

            # Matrice des deviations angulaires en prenant a chaque fois
            # le minimum de deviation avec theta_high et theta_low
            mx_theta_rel = numpy.min([
                    numpy.abs(mx_angles-numpy.radians(theta_high)),
                    numpy.abs(mx_angles-numpy.radians(theta_low)),
                ], axis=0,
            )
            # Matrice des distances entre le ieme O et le jeme N
            mx_distance = distance(mx_pos_o, mx_pos_n)
            ######
            return numpy.sum(
                ((sigma/mx_distance)**6-(sigma/mx_distance)**4)*numpy.cos(mx_theta_rel)**m
            )
        self._call = call

    def __call__(self, ocn_atoms_positions, m=4, r_null=3, theta_low=115, theta_high=155):
        return self._call(
            ocn_atoms_positions[0],
            ocn_atoms_positions[2],
            ocn_atoms_positions[5],
            m, r_null, theta_low, theta_high,
        ) + self._call(
            ocn_atoms_positions[1],
            ocn_atoms_positions[3],
            ocn_atoms_positions[4],
            m, r_null, theta_low, theta_high,
        )
