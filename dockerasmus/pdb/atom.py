# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import six
import numpy
import functools

from .. import constants
from .utils import method_requires


class Atom(object):

    def __init__(self, x, y, z, atom_id, atom_name=None, residual=None):
        self.id = atom_id
        self.name = atom_name
        self.x = x
        self.y = y
        self.z = z
        self.residual = residual

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
            raise ValueError
        return numpy.linalg.norm(self.pos - other)

    @property
    @method_requires(["name"], "Cannot find atom type !")
    def mass(self):
        """The mass of the atom
        """
        return constants.ATOMIC_MASSES[self.name[0]]

    @property
    @method_requires(["name", "residual"], "Cannot find atom type !")
    def epsilon(self):
        return constants.AMINOACID_EPSILON[self.residual][self.name]

    @property
    @method_requires(["name", "residual"], "Cannot find atom type !")
    def radius(self):
        return constants.AMINOACID_RADIUS[self.residual][self.name]

    @property
    @method_requires(["name", "residual"], "Cannot find atom type !")
    def charge(self):
        return constants.AMINOACID_CHARGES[self.residual][self.name]
