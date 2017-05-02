# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import unittest

from dockerasmus.pdb import Atom, Residue, Protein

from ..utils import DATADIR


class AtomTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        arginine = Residue(id=1, name="ARG")
        arginine["CA"] = cls.carbon = Atom(0, 0, 0, 1, "CA", arginine)
        cls.oxygen = Atom(1, 1, 1, 2, "O")
        cls.unnamed = Atom(2, 2, 2, 3)


class TestProperties(AtomTestCase):

    def test_mass(self):
        self.assertEqual(self.carbon.mass, 12)
        self.assertEqual(self.oxygen.mass, 15.99491461957)

    def test_mass_unnamed(self):
        with self.assertRaises(ValueError):
            _ = self.unnamed.mass

    def test_pos(self):
        self.assertEqual(list(self.carbon.pos), [0, 0, 0])
        self.assertEqual(list(self.oxygen.pos), [1, 1, 1])

    def test_radius(self):
        self.assertEqual(self.carbon.radius, 1.908)

    def test_radius_unnamed(self):
        with self.assertRaises(ValueError):
            _ = self.unnamed.radius

    def test_radius_no_residue(self):
        with self.assertRaises(ValueError):
            _ = self.oxygen.radius

    def test_charge(self):
        self.assertEqual(self.carbon.charge, -0.2637)

    def test_charge_unnamed(self):
        with self.assertRaises(ValueError):
            _ = self.unnamed.charge

    def test_charge_no_residue(self):
        with self.assertRaises(ValueError):
            _ = self.oxygen.charge


class TestMagicMethods(AtomTestCase):

    def test_repr(self):
        self.assertEqual(str(self.carbon), "Atom 1(0, 0, 0)")
        self.assertEqual(str(self.oxygen), "Atom 2(1, 1, 1)")

    def test_eq(self):
        self.assertEqual(self.oxygen, Atom(1, 1, 1, 2, "O"))
        self.assertNotEqual(self.oxygen, 42)
        self.assertNotEqual(self.oxygen, Atom(1, 1, 1, 3, "O"))


class TestMethods(AtomTestCase):

    def setUp(self):
        self.arginine = Protein.from_pdb_file(
            os.path.join(DATADIR, 'arginine.pdb')
        )["A"][-3]

    def test_distance_to(self):
        self.assertEqual(self.carbon.distance_to([0, 0, 1]), 1)
        self.assertEqual(self.carbon.distance_to([0, 4, 3]), 5)
        self.assertEqual(self.oxygen.distance_to([1, 2, 3]), 5**.5)


    def test_distance_to_invalid_length(self):
        with self.assertRaises(ValueError):
            _ = self.oxygen.distance_to([1, 2])

    def test_distance_to_other_type(self):
        with self.assertRaises(TypeError):
            _ = self.oxygen.distance_to(42)

    def test_nearest(self):
        self.assertEqual(
            self.arginine['CA'].nearest('CB'),
            self.arginine['CB']
        )
        self.assertEqual(
            self.arginine['CA'].nearest('O'),
            self.arginine['O'],
        )

    def test_nearest_no_residue(self):
        with self.assertRaises(ValueError):
            _ = self.oxygen.nearest("C")

    def test_nearest_not_found(self):
        with self.assertRaises(ValueError):
            _ = self.arginine['CA'].nearest('CE3')
