# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import division

import unittest
import numpy
import collections

from dockerasmus.pdb import Protein, Chain, Residual, Atom
from dockerasmus.score import ScoringFunction
from dockerasmus.score.components import LennardJones
from dockerasmus.utils.matrices import distance

from ...utils import mock


class TestLennardJones(unittest.TestCase):

    #   TEST SITUATION: 4 atoms of the same kind
    #   A & B from prot 1, C & D from prot 2
    #
    #
    #   A(id:1, x=0, y=0) B(id:2, x=1, y=0)
    #     x===============x
    #     .     PROT 1    .
    #     .               .
    #     .               .
    #     .     PROT 2    .
    #     x===============X
    #   C(id:1, x=0, y=1) D(id:2, x=1, y=1)
    #
    #         pwd =  1  1     radius = .5  .5  (so they sum to 1)
    #                1  1              .5  .5
    #
    #    distance =  1 √2   (since A->D = (1,1), B->C = (1,1)
    #               √2  1      and A->C = (0,1), B->D = (0,1))
    #
    # Since we're only summing for i>j, we only are
    # computing Van der Waals on B & C
    #
    # a =     pwd * radius¹² = 1
    # b = 2 * pwd * radius⁶  = 2
    # d = √2
    #
    # LJ = a/d^12  - b/d^6
    #    = 1/64 - 2/8



    def setUp(self):
        self.eps = numpy.array([[1, 1], [1, 1]])
        self.vdw_radius = numpy.array([[.5, .5], [.5, .5]])
        self.distance = distance(numpy.array([[0,0], [0,1]]),
                                 numpy.array([[1,0], [1,1]]))

    def test_numpy(self):
        lj = LennardJones(force_backend='numpy')
        self.assertAlmostEqual(
            float(lj(self.eps, self.vdw_radius, self.distance)),
            1/64 - 2/8
        )

    def test_theano(self):
        lj = LennardJones(force_backend='theano')
        self.assertAlmostEqual(
            float(lj(self.eps, self.vdw_radius, self.distance)),
            1/64 - 2/8
        )

    def test_wrapped(self):
        a = Atom(0, 0, 0, 1)
        b = Atom(1, 0, 0, 2)
        c = Atom(0, 1, 0, 3)
        d = Atom(1, 1, 0, 4)

        r1 = Residual(1, atoms=collections.OrderedDict([('A', a), ('B', b)]))
        r2 = Residual(2, atoms=collections.OrderedDict([('C', c), ('D', d)]))

        c1 = Chain('A', residuals=collections.OrderedDict([(1, r1)]))
        c2 = Chain('A', residuals=collections.OrderedDict([(2, r2)]))

        prot1 = Protein(chains=collections.OrderedDict([('A', c1)]))
        prot2 = Protein(chains=collections.OrderedDict([('A', c2)]))

        with mock.patch('dockerasmus.pdb.Atom.pwd',
                        new_callable=mock.PropertyMock,
                        return_value=1):
            with mock.patch('dockerasmus.pdb.Atom.radius',
                            new_callable=mock.PropertyMock,
                            return_value=.5):
                score_lj = ScoringFunction(LennardJones)
                self.assertAlmostEqual(score_lj(prot1, prot2), 1/64 - 2/8)
