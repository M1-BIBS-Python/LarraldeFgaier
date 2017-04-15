# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import division

import unittest
import numpy
import collections

from dockerasmus.pdb import Protein, Chain, Residual, Atom
from dockerasmus.score import ScoringFunction
from dockerasmus.score.components import Coulomb
from dockerasmus.utils.matrices import distance

from ...utils import mock


class TestCoulomb(unittest.TestCase):

    #   TEST SITUATION: 4 atoms
    #   A & B from prot 1, C & D from prot 2
    #
    #   A, C charged +1
    #      B charged +2
    #      D charged -1
    #
    #
    #  A(id:1, x=0, y=0) B(id:2, x=1, y=0)
    #   x================x
    #   .     PROT 1     .
    #   .                .
    #   .                .
    #   .     PROT 2     .
    #   x================x
    #  C(id:3, x=0, y=1) D(id:4, x=1, y=1)
    #
    #     distance =  1 √2     q = +1 +2
    #                √2  1         -1 -2
    #
    #   CL = 1/(diel*1)  - 1/(diel*√2) + 2/(diel*√2) - 2/(diel*1)
    #      = 1/(diel*√2) - 1/(diel*1)
    #      = √2/(diel*2) - 2/(diel*2)
    #      = (√2-2)/(2*diel)

    diel = 65.0
    expected = (2**.5 - 2) / (2*diel)

    def setUp(self):
        self.charges = numpy.array([[1, 2], [1, -1]])
        self.positions = [numpy.array([[0,0], [0,1]]),
                          numpy.array([[1,0], [1,1]])]
        self.distances = distance(*self.positions)

    def test_numpy(self):
        cl = Coulomb(force_backend='numpy')
        self.assertAlmostEqual(
            float(cl(self.charges, self.distances)),
            self.expected,
        )

    def test_mxnet(self):
        cl = Coulomb(force_backend='mxnet')
        self.assertAlmostEqual(
            float(cl(self.charges, self.distances)),
            self.expected,
        )

    def test_theano(self):
        cl = Coulomb(force_backend='theano')
        self.assertAlmostEqual(
            float(cl(self.charges, self.distances)),
            self.expected,
        )


    def test_wrapped(self):
        # mock ligand
        prot1 = mock.MagicMock()
        prot1.atom_charges.return_value = self.charges[0]
        prot1.atom_positions.return_value = self.positions[0]
        # mock receptor
        prot2 = mock.MagicMock()
        prot2.atom_charges.return_value = self.charges[1]
        prot2.atom_positions.return_value = self.positions[1]
        # test scoring
        score_lj = ScoringFunction(Coulomb)
        self.assertAlmostEqual(score_lj(prot1, prot2), self.expected)
