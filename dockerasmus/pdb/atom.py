# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import six
import numpy
import functools

from .. import constants
from ..utils.decorators import method_requires


class Atom(object):
    __slots__ = ("id", "name", "x", "y", "z", "residual")

    def __init__(self, x, y, z, id, name=None, residual=None):
        self.id = id
        self.name = name
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
    def epsilon(self):
        return self._read_from_constants(constants.AMINOACID_EPSILON)

    @property
    def radius(self):
        return self._read_from_constants(constants.AMINOACID_RADIUS)

    @property
    def charge(self):
        return self._read_from_constants(constants.AMINOACID_CHARGES)

    @method_requires(["name", "residual"], "Cannot find atom type !")
    def _read_from_constants(self, table):
        try:
            return table[self.residual.name][self.name]
        except KeyError:
            if self.residual.nter:
                return table["NTER"][self.name]
            elif self.residual.cter:
                return table["CTER"][self.name]
            else:
                raise KeyError("{}, {} ({})".format(self.residual.name, self.name, self.id))
