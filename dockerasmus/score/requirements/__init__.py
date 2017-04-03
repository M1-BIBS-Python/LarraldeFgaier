# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

from .distance import distance

__all__ = [
    "epsilon", "distance", "vdw_radius", "charge", "ocn_atoms_positions",
]


def epsilon(protein1, protein2):
    return protein1.atom_epsilon(), protein2.atom_epsilon()


def vdw_radius(protein1, protein2):
    return protein1.atom_radius(), protein2.atom_radius()


def charge(protein1, protein2):
    return protein1.atom_charges(), protein2.atom_charges()


def ocn_atoms_positions(protein1, protein2):
    import numpy

    return (
        # Position of O atoms
        numpy.array([atom.pos for atom in protein1.iteratoms() if atom.name.startswith('O')]),
        numpy.array([atom.pos for atom in protein2.iteratoms() if atom.name.startswith('O')]),

        # Positions of C atoms linked to each O atom
        numpy.array([atom.nearest("C").pos for atom in protein1.iteratoms() if atom.name.startswith('O')]),
        numpy.array([atom.nearest("C").pos for atom in protein2.iteratoms() if atom.name.startswith('O')]),

        # Positions of N atoms
        numpy.array([atom.pos for atom in protein1.iteratoms() if atom.name.startswith('N')]),
        numpy.array([atom.pos for atom in protein2.iteratoms() if atom.name.startswith('N')]),
    )
