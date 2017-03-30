# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import division

import theano
import numpy
from ...utils.tensors import distance, normalized


class TheanoComponent(object):
    pass


class _LennardJones(TheanoComponent):

    def __init__(self):

        print("Compiling Lennard Jones")

        ### Epsilon matrix from protein vectors
        v_eps1, v_eps2 = theano.tensor.dvectors('v_eps1', 'v_eps2')
        mx_eps = theano.tensor.sqrt(theano.tensor.outer(v_eps1, v_eps2))

        ### Radius matrix from protein vectors
        v_rad1, v_rad2 = theano.tensor.dvectors('v_rad1', 'v_rad2')
        mx_rad = theano.tensor.transpose(theano.tensor.repeat(
            theano.tensor.reshape(
                v_rad1, (1, v_rad1.shape[0])
            ), v_rad2.shape[0], axis=0,
        )) + v_rad2

        ### Atomwise distance matrix
        mx_distance = theano.tensor.dmatrix('mx_distance')

        ### Van der Waals constants
        mx_A = mx_eps * mx_rad**12
        mx_B = 2 * mx_eps * mx_rad**6

        self._call = theano.function(
            [v_eps1, v_eps2, v_rad1, v_rad2, mx_distance],
            theano.tensor.sum(theano.tensor.triu(
                mx_A/(mx_distance**12) - mx_B/(mx_distance**6), k=1,
            ))
        )

    def __call__(self, v_eps1, v_eps2, v_rad1, v_rad2, mx_distance):
        return self._call(v_eps1, v_eps2, v_rad1, v_rad2, mx_distance)


class _Coulomb(TheanoComponent):

    def __init__(self):

        print("Compiling Coulomb")

        ### Dielectric constant
        diel = theano.tensor.dscalar('diel')

        ### Charge matrix from protein vectors
        v_q1, v_q2 = theano.tensor.dvectors('v_q1', 'v_q2')
        mx_q = theano.tensor.outer(v_q1, v_q2)

        ### Atomwise distance matrix
        mx_distance = theano.tensor.dmatrix('mx_distance')

        self._call = theano.function(
            [v_q1, v_q2, mx_distance, diel],
            theano.tensor.sum(theano.tensor.triu(
                mx_q/(diel*mx_distance)
            ))
        )

    def __call__(self, v_q1, v_q2, mx_distance, diel):
        return self._call(v_q1, v_q2, mx_distance, diel)


class _HBonds(TheanoComponent):

    def __init__(self):

        print("Compiling HBonds")

        r_null = theano.tensor.dscalar('r_null')
        sigma = r_null * theano.tensor.sqrt(2/3)

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
        mx_angles = theano.tensor.arccos(mx_cos)

        theta_low, theta_high = theano.tensor.dscalars('theta_low', 'theta_high')
        # Matrice des deviations angulaires en prenant a chaque fois
        # le minimum de deviation avec theta_high et theta_low
        mx_theta_rel = theano.tensor.min([
                abs(mx_angles - theano.tensor.deg2rad(theta_high)),
                abs(mx_angles - theano.tensor.deg2rad(theta_low)),
            ], axis=0,
        )
        #mx_theta_rel = abs(mx_angles - theano.tensor.deg2rad(150))

        # Matrice des distances entre le ieme O et le jeme N
        mx_distance = distance(mx_pos_o, mx_pos_n)


        m = theano.tensor.dscalar('m')
        ######
        self._call = theano.function(
            [mx_pos_o, mx_pos_c, mx_pos_n, m, r_null, theta_low, theta_high],
            theano.tensor.sum(
                ((sigma/mx_distance)**6-(sigma/mx_distance)**4)*numpy.cos(mx_theta_rel)**m
            ),
        )

    def __call__(self, mx_pos_o, mx_pos_c, mx_pos_n, m, r_null, theta_low, theta_high):
        return self._call(
            mx_pos_o, mx_pos_c, mx_pos_n,
            m, r_null, theta_low, theta_high,
        )


lennard_jones = _LennardJones()
coulomb = _Coulomb()
h_bonds = _HBonds()
