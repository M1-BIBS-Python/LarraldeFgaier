# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

from .distance import distance

__all__ = [
    "potential_well_depth", "distance", "vdw_radius", "charge",
    "ocn_atoms_positions",
]


def potential_well_depth(protein1, protein2):
    """The :math:`\epsilon` of the atoms of ``protein1`` and ``protein2``.
    """
    return protein1.atom_pwd(), protein2.atom_pwd()


def vdw_radius(protein1, protein2):
    """The Van der Waals radius of the atoms of ``protein1`` and ``protein2``.
    """
    return protein1.atom_radius(), protein2.atom_radius()


def charge(protein1, protein2):
    """The charge of the atoms of ``protein1`` and ``protein2``.
    """
    return protein1.atom_charges(), protein2.atom_charges()


def ocn_atoms_positions(protein1, protein2):
    """The positions of *O*, *C* and *N* atoms in ``protein1`` and ``protein2``.
    """
    import numpy

    return (
        # Position of O atoms
        numpy.array([atom.pos for atom in protein1.iteratoms() if atom.name.startswith('O')] or numpy.zeros((0,3))),
        numpy.array([atom.pos for atom in protein2.iteratoms() if atom.name.startswith('O')] or numpy.zeros((0,3))),

        # Positions of C atoms linked to each O atom
        numpy.array([atom.nearest("C").pos for atom in protein1.iteratoms() if atom.name.startswith('O')] or numpy.zeros((0,3))),
        numpy.array([atom.nearest("C").pos for atom in protein2.iteratoms() if atom.name.startswith('O')] or numpy.zeros((0,3))),

        # Positions of N atoms
        numpy.array([atom.pos for atom in protein1.iteratoms() if atom.name.startswith('N')] or numpy.zeros((0,3))),
        numpy.array([atom.pos for atom in protein2.iteratoms() if atom.name.startswith('N')] or numpy.zeros((0,3))),
    )
