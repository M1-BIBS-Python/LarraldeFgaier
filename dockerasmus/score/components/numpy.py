# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import numpy

from ...utils.matrices import normalized, distance, compose


def lennard_jones(v_eps1, v_eps2, v_rad1, v_rad2, mx_distance):
    # Matrix of van der waals epsilon for each aminoacid atom
    mx_epsilon = numpy.sqrt(numpy.outer(v_eps1, v_eps2))
    # Matrix of Van der Waals optimal radius for each aminoacid atom
    # to the power of 6
    mx_radius_6 = compose(lambda r1, r2: r1+r2, v_rad1, v_rad2)**6
    # Matrices A and B from Cornell
    mx_B = 2 * mx_epsilon * mx_radius_6
    mx_A = 0.5 * mx_B * mx_radius_6
    # Matrix of aminoacid-wise distance
    mx_distance_6 = mx_distance**6
    #######
    return numpy.sum(numpy.triu(mx_A/(mx_distance_6**2)-mx_B/(mx_distance_6), 1))


def coulomb(v_q1, v_q2, mx_distance, diel=65.0, upper=True):
    # Matrix of aminoacid-wise charge
    mx_q = compose(lambda q1, q2: q1*q2, v_q1, v_q2)
    #####
    return numpy.sum(numpy.triu(mx_q/(diel*mx_distance), 1))


def h_bonds(mx_pos_o, mx_pos_c, mx_pos_n, m=4, r_null=3, theta_low=115, theta_high=155):
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
    mx_cos = numpy.sum(mx_vec_o_to_c*mx_vec_o_to_n, axis=2)
    mx_angles = numpy.arccos(mx_cos)
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


def rmsd(mx_pos1, mx_pos2):
    N = mx_pos_1.shape[0]
    mx_d = distance(mx_pos_1, mx_pos_2)
    return numpy.sqrt(numpy.sum(mx_d**2)/N)
