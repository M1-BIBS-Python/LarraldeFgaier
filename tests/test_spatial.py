# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import unittest
import os
import math
import numpy

from dockerasmus.pdb import Protein
from dockerasmus import spatial

from .utils import DATADIR


class TestArrays(unittest.TestCase):

    def assertArrayAlmostEqual(self, actual, desired, decimal=7):
        numpy.testing.assert_almost_equal(actual, desired, decimal=decimal)


class TestMatrices(TestArrays):

    def test_rotation_matrix(self):
        # No rotation
        no_rotation = spatial.RotationMatrix()
        self.assertArrayAlmostEqual(no_rotation, numpy.identity(4))
        # Rotation around x axis
        point = numpy.array([0, 1, 1, 1])
        rotated_point = spatial.RotationMatrix(math.pi/2, 0, 0).dot(point)
        self.assertArrayAlmostEqual(rotated_point, [0., -1., 1., 1.])
        # Rotation around z axis
        point = numpy.array([1, 1, 0, 1])
        rotated_point = spatial.RotationMatrix(0, 0, math.pi).dot(point)
        self.assertArrayAlmostEqual(rotated_point, [-1., -1., 0., 1.])

    def test_translation_matrix(self):
        # No translation
        no_translation = spatial.TranslationMatrix()
        self.assertArrayAlmostEqual(no_translation, numpy.identity(4))
        # with translation
        point = numpy.array([0, 0, 0, 1])
        translated_point = spatial.TranslationMatrix(1, 1, 1).dot(point)
        self.assertArrayAlmostEqual(translated_point, [1, 1, 1, 1])


class TestTransform(TestArrays):

    @classmethod
    def setUpClass(cls):
        cls.arginine = Protein.from_pdb_file(
            os.path.join(DATADIR, 'arginine.pdb')
        )

    def test_transform_cartesian(self):
        new_arginine = spatial.transform_cartesian(
            self.arginine, x=2, y=2, z=1
        )

        for new_atom in new_arginine['A'][-3].values():
            self.assertArrayAlmostEqual(
                new_atom.pos,
                self.arginine['A'][-3][new_atom.name].pos + numpy.array([2, 2, 1])
            )

    def test_transform_spherical(self):
        new_arginine = spatial.transform_spherical(
            self.arginine, r=1, phi=0, theta=0
        )

        for new_atom in new_arginine['A'][-3].values():
            self.assertArrayAlmostEqual(
                new_atom.pos,
                self.arginine['A'][-3][new_atom.name].pos + numpy.array([0, 0, 1])
            )

        new_arginine_2 = spatial.transform_spherical(
            self.arginine, r=2, phi=math.pi/2, theta=math.pi/2
        )

        for new_atom_2 in new_arginine_2['A'][-3].values():
            self.assertArrayAlmostEqual(
                new_atom_2.pos,
                self.arginine['A'][-3][new_atom_2.name].pos + numpy.array([0, 2, 0])
            )
