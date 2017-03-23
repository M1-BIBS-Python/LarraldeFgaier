# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import six
import numpy


class Atom(object):

    _atomic_mass = {
        'C': 12.0107, 'H': 1.00794, 'N': 14.0067, 'O': 15.9994, 'S': 32.065,
    }

    def __init__(self, x, y, z, atom_id, atom_name=None):
        self.id = atom_id
        self.name = atom_name
        self.x = x
        self.y = y
        self.z = z

    @property
    def pos(self):
        """The position of the atom
        """
        return numpy.array([self.x, self.y, self.z])

    def __repr__(self):
        return "Atom {}({}, {}, {})".format(self.id, self.x, self.y, self.z)

    def __eq__(self, other):
        return all(
            getattr(self, attr, None)==getattr(other, attr, None)
                for attr in ('x', 'y', 'z', 'id')
        )

    if six.PY3:
        def itervalues(self):
            return six.itervalues(self)

        def iteritems(self):
            return six.iteritems(self)

    def distance_to(self, other):
        """Computes the distance to ``other``

        Arguments:
            other (numpy.array): the position to compute the
                distance to (must be array-like of dimension 3)
        """
        if len(other) != 3:
            raise ValueError("")
        return numpy.linalg.norm(self.pos - other)

    @property
    def mass(self):
        """The mass of the atom
        """
        if self.name is None:
            raise ValueError("Cannot find atom type !")
        return self._atomic_mass[self.name[0]]
