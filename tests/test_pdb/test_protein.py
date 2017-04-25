# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import unittest
import collections

from dockerasmus.pdb import Protein, Chain, Atom, Residual
from dockerasmus.constants import ATOMIC_MASSES

from ..utils import DATADIR


class TestProtein(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        ## A protein made only of one arginine
        cls.arginine_prot = Protein.from_pdb_file(
            os.path.join(DATADIR, 'arginine.pdb')
        )
        ## An unrealistic protein, which properties
        ## are easy to compute manually
        cls.prot = Protein(chains={
            'A': Chain('A', residuals={
                1: Residual(1, atoms={
                    'C1': Atom(0, 0, 0, 1, 'C1'),
                    'C2': Atom(0, 0, 1, 2, 'C2')
                })
            })
        })

    def setUp(self):
        ## Custom made protein with elements accessible
        ## outside of prot.__getitem__
        self.atom_1 = Atom(0,0,0,1)
        self.atom_2 = Atom(1,0,0,2)
        self.atom_3 = Atom(0,1,0,3)
        self.atom_4 = Atom(0,0,1,4)
        self.res_1 = Residual(1, atoms={'C1': self.atom_1, 'C2': self.atom_2})
        self.res_2 = Residual(2, atoms={'C1': self.atom_3})
        self.res_3 = Residual(3, atoms={'C1': self.atom_4})
        self.chain_b = Chain('B', residuals={1: self.res_1})
        self.chain_c = Chain('C', residuals={2: self.res_2})
        self.chain_d = Chain('D', residuals={3: self.res_3})
        self.prot2 = Protein(chains={'B': self.chain_b, 'C': self.chain_c,
                                     'D': self.chain_d})


class TestProperties(TestProtein):

    def test_mass(self):
        self.assertEqual(self.prot.mass, 2*ATOMIC_MASSES['C'])

    def test_mass_center(self):
        self.assertEqual(list(self.prot.mass_center), [0, 0, .5])

    def test_radius(self):
        self.assertEqual(self.prot.radius, .5)


class TestMagicMethods(TestProtein):
    ## TODO: test_getitem

    def test_getitem_simple(self):
        self.assertEqual(self.prot2['B'], self.chain_b)
        self.assertEqual(self.prot2['D'], self.chain_d)
        with self.assertRaises(KeyError):
            _ = self.prot2['A']

    def test_getitem_wordslice(self):
        prot_1 = Protein(chains=collections.OrderedDict([
            ('B', self.chain_b), ('C', self.chain_c), ('D', self.chain_d)
        ]))
        prot_2 = Protein(chains=collections.OrderedDict([
            ('B', self.chain_b), ('C', self.chain_c),
        ]))
        prot_3 = Protein(chains=collections.OrderedDict([
            ('C', self.chain_c), ('D', self.chain_d)
        ]))
        prot_4 = Protein(chains=collections.OrderedDict([
            ('B', self.chain_b),
        ]))

        self.assertEqual(self.prot2[:], prot_1)
        self.assertEqual(self.prot2[:'D'], prot_2)
        self.assertEqual(self.prot2['C':], prot_3)
        self.assertEqual(self.prot2['B':'C'], prot_4)

    def test_contains_chain(self):
        self.assertIn(self.arginine_prot['A'], self.arginine_prot)
        self.assertNotIn(Chain('B'), self.arginine_prot)

    def test_contains_residual(self):
        self.assertIn(self.arginine_prot['A'][-3], self.arginine_prot)
        self.assertNotIn(Residual(-4), self.arginine_prot)

    def test_contains_atom(self):
        self.assertIn(self.arginine_prot['A'][-3]['O'], self.arginine_prot)
        self.assertNotIn(Atom(0,0,0,4), self.arginine_prot)

    def test_contains_text(self):
        self.assertIn('A', self.arginine_prot)
        self.assertNotIn('B', self.arginine_prot)

    def test_contains_othertype(self):
        with self.assertRaises(TypeError):
            _ = [1] in self.arginine_prot
        with self.assertRaises(TypeError):
            _ = set() in self.arginine_prot
        with self.assertRaises(TypeError):
            _ = object() in self.arginine_prot

    def test_iadd_protein(self):
        test_prot = Protein()
        test_prot += self.prot2
        self.assertEqual(set(test_prot), set(self.prot2))

    def test_iadd_chain(self):
        test_prot = Protein()
        test_prot += self.chain_b
        self.assertEqual(set(test_prot), {'B'})

    def test_iadd_same_chain(self):
        test_prot = Protein(chains={'B': self.chain_b})
        with self.assertRaises(ValueError):
            test_prot += self.chain_b

    def test_iadd_same_protein(self):
        test_prot = Protein(chains={'B': self.chain_b})
        with self.assertRaises(ValueError):
            test_prot += self.prot2

    def test_iadd_othertype(self):
        test_prot = Protein()
        with self.assertRaises(TypeError):
            test_prot += 1
        with self.assertRaises(TypeError):
            test_prot += "a"




class TestMethods(TestProtein):

    def test_nearest_atom(self):
        pos1 = [11.296, 86.721, 94.521] # Coordinates of Atom 32
        self.assertEqual(self.arginine_prot.nearest_atom(pos1).id, 32)

        pos2 = [12.924, 87.357, 96.42]
        self.assertEqual(self.arginine_prot.nearest_atom(pos2).id, 38)

    def test_copy(self):
        arginine1 = self.arginine_prot
        arginine2 = arginine1.copy()

        for chain1, chain2 in zip(arginine1.itervalues(), arginine2.itervalues()):
            for res1, res2 in zip(chain1.itervalues(), chain2.itervalues()):
                for atom_id in res1:
                    self.assertEqual(res1[atom_id], res2[atom_id])
                    self.assertIsNot(res1[atom_id], res2[atom_id])
                    res1[atom_id].x += 1
                    self.assertNotEqual(res1[atom_id], res2[atom_id])

    def test_iteritems(self):
        self.assertEqual(
            list(self.arginine_prot.iteritems()),
            [('A', self.arginine_prot['A'])],
        )

    def test_itervalues(self):
        self.assertEqual(
            list(self.arginine_prot.itervalues()),
            [self.arginine_prot['A']],
        )

    def test_iteratoms(self):
        self.assertEqual(
            sorted(self.prot.iteratoms(), key=lambda a: a.id),
            [self.prot['A'][1]['C1'],
             self.prot['A'][1]['C2']]
        )

    def test_rmsd_ref(self):
        ## rmsd = sqrt(0.5[(0-0)²+(0-0)²+(0-0)²+(0-0)²+(0-0)²+(1-0)²])
        ##      = sqrt(0.5 * 1²) = sqrt(0.5)
        self.assertEqual(self.prot.rmsd([0,0,0]), (.5)**.5)

    def test_rmsd_prot(self):
        self.assertEqual(self.prot.rmsd(self.prot), 0)

    def test_rmsd_othertype(self):
        with self.assertRaises(TypeError):
            _ = self.prot.rmsd(1)
        with self.assertRaises(TypeError):
            _ = self.prot.rmsd("test")

    def test_contact_map_nearest(self):
        self.assertEqual( ## Nearest atom are identical
            self.prot.contact_map(self.prot, mode='nearest')[1,1],
            0,
        )

    def test_contact_map_farthest(self):
        self.assertEqual(
            self.prot.contact_map(self.prot, mode='farthest')[1,1],
            1,
        )

    def test_contact_map_mass_center(self):
        self.assertEqual( ## Same mass center
            self.prot.contact_map(self.prot, mode='mass_center')[1,1],
            0,
        )

    def test_contact_map_othermode(self):
        with self.assertRaises(ValueError):
            _ = self.prot.contact_map(self.prot, mode='rubbish')

    def test_contact_map_othertype(self):
        with self.assertRaises(TypeError):
            _ = self.prot.contact_map(1)

    def test_atom_int(self):
        self.assertEqual(self.prot2.atom(1), self.atom_1)
        self.assertEqual(self.prot2.atom(4), self.atom_4)
        with self.assertRaises(KeyError):
            _ = self.prot2.atom(40)

    def test_atom_other(self):
        self.assertEqual(self.prot2.atom("1"), self.atom_1)
        self.assertEqual(self.prot2.atom(1.0), self.atom_1)
        with self.assertRaises(TypeError):
            _ = self.prot2.atom("hello")

    def test_residual_int(self):
        self.assertEqual(self.prot2.residual(1), self.res_1)
        self.assertEqual(self.prot2.residual(3), self.res_3)
        with self.assertRaises(KeyError):
            _ = self.prot2.residual(42)

    def test_residual_other(self):
        self.assertEqual(self.prot2.residual("1"), self.res_1)
        self.assertEqual(self.prot2.residual(2.0), self.res_2)
        with self.assertRaises(TypeError):
            _ = self.prot2.residual("hi!")
