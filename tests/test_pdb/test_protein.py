# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import unittest
import collections

from dockerasmus.pdb import Protein, Chain, Atom, Residue
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
        res = Residue(1, 'LEU')
        res['C1'] = Atom(0, 0, 0, 1, 'CA', res)
        res['C2'] = Atom(0, 0, 1, 2, 'CB', res)
        cls.prot = Protein(chains={'A': Chain('A', residues={1:res})})

        ## Custom made protein with elements accessible
        ## outside of prot.__getitem__
        cls.atom_1 = Atom(0,0,0,1)
        cls.atom_2 = Atom(1,0,0,2)
        cls.atom_3 = Atom(0,1,0,3)
        cls.atom_4 = Atom(0,0,1,4)
        cls.res_1 = Residue(1, atoms={'C1': cls.atom_1, 'C2': cls.atom_2})
        cls.res_2 = Residue(2, atoms={'C1': cls.atom_3})
        cls.res_3 = Residue(3, atoms={'C1': cls.atom_4})
        cls.chain_b = Chain('B', residues={1: cls.res_1})
        cls.chain_c = Chain('C', residues={2: cls.res_2})
        cls.chain_d = Chain('D', residues={3: cls.res_3})

    def setUp(self):
        ## Reset self.prot2 before each test since it can be
        ## mutated by __iadd__
        self.prot2 = Protein(chains={'B': self.chain_b, 'C': self.chain_c,
                                     'D': self.chain_d})


class TestProperties(TestProtein):

    def test_mass(self):
        self.assertEqual(self.prot.mass, 2*ATOMIC_MASSES['C'])

    def test_mass_center(self):
        self.assertEqual(list(self.prot.mass_center), [0, 0, .5])

    def test_radius(self):
        self.assertEqual(self.prot.radius, .5)

    def test_atom_charges(self):
        self.assertEqual(list(self.prot.atom_charges()), [-0.0518, -0.1102])


class TestMagicMethods(TestProtein):

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

    def test_contains_residue(self):
        self.assertIn(self.arginine_prot['A'][-3], self.arginine_prot)
        self.assertNotIn(Residue(-4), self.arginine_prot)

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

    def test_add_protein(self):
        prot = Protein() + self.prot2
        self.assertEqual(set(prot), set(self.prot2))

    def test_add_chain(self):
        prot = Protein() + self.chain_b
        self.assertEqual(set(prot), {'B'})

    def test_add_same_chain(self):
        with self.assertRaises(ValueError):
            _ = Protein() + self.chain_b + self.chain_b

    def test_add_same_protein(self):
        with self.assertRaises(ValueError):
            _ = self.prot2 + self.prot2

    def test_iadd_othertype(self):
        with self.assertRaises(TypeError):
            _ = self.prot2 + 1
        with self.assertRaises(TypeError):
            _ = self.prot2 + ["a", "b", "c"]

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
        self.atom_1 = Atom(0,0,0,1)
        self.atom_2 = Atom(1,0,0,2)
        self.atom_3 = Atom(0,1,0,3)
        self.atom_4 = Atom(0,0,1,4)
        self.res_1 = Residue(1, atoms={'C1': self.atom_1, 'C2': self.atom_2})
        self.res_2 = Residue(2, atoms={'C1': self.atom_3})
        self.res_3 = Residue(3, atoms={'C1': self.atom_4})
        self.chain_b = Chain('B', residues={1: self.res_1})
        self.chain_c = Chain('C', residues={2: self.res_2})
        self.chain_d = Chain('D', residues={3: self.res_3})
        self.prot2 = Protein(chains={'B': self.chain_b, 'C': self.chain_c,
                                     'D': self.chain_d})
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

    def test_atom_othertype(self):
        self.assertEqual(self.prot2.atom("1"), self.atom_1)
        self.assertEqual(self.prot2.atom(1.0), self.atom_1)
        with self.assertRaises(TypeError):
            _ = self.prot2.atom("hello")

    def test_residue_int(self):
        self.assertEqual(self.prot2.residue(1), self.res_1)
        self.assertEqual(self.prot2.residue(3), self.res_3)
        with self.assertRaises(KeyError):
            _ = self.prot2.residue(42)

    def test_residue_othertype(self):
        self.assertEqual(self.prot2.residue("1"), self.res_1)
        self.assertEqual(self.prot2.residue(2.0), self.res_2)
        with self.assertRaises(TypeError):
            _ = self.prot2.residue("hi!")

    def test_interface_identical_proteins(self):
        # Since the 2 proteins are the same,
        # every residue has a very close interfacing
        # residue: itself in the other protein !
        self.assertEqual(
            list(self.arginine_prot.interface(self.arginine_prot)),
            [(self.arginine_prot['A'][-3],)*2],
        )

    def test_interface_othertype(self):
        with self.assertRaises(TypeError):
            for r1, r2 in self.prot2.interface(1):
                pass
