# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import division

import math
import unittest
import numpy
import collections
import warnings

from dockerasmus.pdb import Protein, Chain, Residual, Atom
from dockerasmus.score import ScoringFunction
from dockerasmus.score.components import Fabiola
from dockerasmus.utils.matrices import distance

from ...utils import mock, suppress_tf_log


class TestFabiola(unittest.TestCase):

    #   TEST SITUATION: 2N on the ligand, 2C + 2O on the receptor
    #
    #
    #   N1 (0, 3)   N2 (1, 3)
    #    x         x
    #
    #
    #   O1 (0, 1)
    #    x                   x  O2 (3, 1)
    #
    #    x             x
    #   C1 (0, 0)     C2 (2, 0)
    #
    #
    #
    #   angle(C1, O1, N1) = 180°
    #   angle(C1, O1, N2) = 153.4°
    #   angle(C2, O2, N1) = 78.7°
    #   angle(C2, O2, N2) = 90°
    #
    #
    #   ||O1 N1|| = 2
    #   ||O1 N2|| = √5
    #   ||O2 N1|| = √13
    #   ||O2 N2|| = √8
    #
    #
    #   With default values: m=4, θ_low=115, θ_high=155 r_0=3 (=> σ=√6)
    #
    #   FAB = [(√6 /   2)⁶ - (√6 /   2)⁴] cos(180        - θ_high)⁴
    #       + [(√6 / √13)⁶ - (√6 / √13)⁴] cos( 78.690067 - θ_low )⁴
    #       + [(√6 /  √5)⁶ - (√6 /  √5)⁴] cos(153.434948 - θ_high)⁴
    #       + [(√6 /  √8)⁶ - (√6 /  √8)⁴] cos( 90        - θ_low )⁴

    expected_default = (
        # C1 - O1 ... N1
        ((6**.5/2)**6 - (6**.5/2)**4) * math.cos(math.radians(180-155))**4
        # C2 - O2 ... N1
      + ((6**.5/13**.5)**6 - (6**.5/13**.5)**4) * math.cos(math.radians(78.690067-115))**4
        # C1 - O1 ... N2
      + ((6**.5/5**.5)**6 - (6**.5/5**.5)**4) * math.cos(math.radians(153.434948-155))**4
        # C2 - O2 ... N2
      + ((6**.5/8**.5)**6 - (6**.5/8**.5)**4) * math.cos(math.radians(90-115))**4
    )


    def setUp(self):
        self.ocn_atoms_positions = [
            numpy.zeros((0, 3)),                     # O lig
            numpy.array([[0, 1, 0], [3, 1, 0]]), # O rec
            numpy.zeros((0, 3)),                     # C lig
            numpy.array([[0, 0, 0], [2, 0, 0]]), # C rec
            numpy.array([[0, 3, 0], [1, 3, 0]]), # N lig
            numpy.zeros((0, 3)),                     # N ref
        ]


    def test_numpy(self):
        fa = Fabiola(force_backend='numpy')
        self.assertAlmostEqual(
            float(fa(self.ocn_atoms_positions)),
            self.expected_default,
        )

    def test_theano(self):
        fa = Fabiola(force_backend='theano')
        self.assertAlmostEqual(
            float(fa(self.ocn_atoms_positions)),
            self.expected_default,
        )


    def test_wrapped(self):

        r_rec = Residual(1)
        r_rec['C1'] = Atom(0, 0, 0, 1, 'C', r_rec)
        r_rec['C2'] = Atom(2, 0, 0, 2, "C", r_rec)
        r_rec['O1'] = Atom(0, 1, 0, 3, "O", r_rec)
        r_rec['O2'] = Atom(3, 1, 0, 4, "O", r_rec)

        r_lig = Residual(2)
        r_lig['N1'] = Atom(0, 3, 0, 5, "N", r_lig)
        r_lig['N2'] = Atom(1, 3, 0, 6, "N", r_lig)

        
        c_rec = Chain('A', residuals=collections.OrderedDict([(1, r_rec)]))
        c_lig = Chain('A', residuals=collections.OrderedDict([(2, r_lig)]))

        rec = Protein(chains=collections.OrderedDict([('A', c_rec)]))
        lig = Protein(chains=collections.OrderedDict([('A', c_lig)]))

        score_lj = ScoringFunction(Fabiola)
        self.assertAlmostEqual(score_lj(rec, lig), self.expected_default)


def setUpModule():
    warnings.simplefilter('ignore', category=ImportWarning)
    warnings.simplefilter('ignore', category=DeprecationWarning)

def tearDownModule():
    warnings.simplefilter(warnings.defaultaction)
