# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import collections
import six
import numpy

from .atom import Atom


class Residual(dict):

    def __init__(self, id, name=None, atoms=None):
        super(Residual, self).__init__(atoms or {})
        self.id = id
        self.name = name

    def __contains__(self, item):
        """Checks if `item` is contained in the residual.

        Arguments:
            item: either an atom_id (`int`) or an `Atom` object
                to check if present within the residual.
        """
        if isinstance(item, six.text_type):
            return super(Residual, self).__contains__(item)
        elif isinstance(item, Atom):
            return any(item==atom for atom in self.itervalues())
        elif isinstance(item, int):
            return any(item == atom.id for atom in self.itervalues())

        else:
            raise TypeError(
                "'in <Residual>' requires Atom or {}"
                " as left operand, not {}".format(
                    six.text_type.__name__,type(item).__name__)
            )

    if six.PY3:
        def itervalues(self):
            return six.itervalues(self)

        def iteritems(self):
            return six.iteritems(self)

    @property
    def mass(self):
        """The mass of the residual.

        Computed as sum of the masses of the non-hydrogen
        atoms of the residual.
        """
        return sum(atom.mass for atom in self.itervalues())

    @property
    def mass_center(self):
        """The position of the mass center of the residual.

        Computed as the barycenter of the positions of the
        atoms weighted by their atomic masses.
        """
        mass = self.mass
        return sum((atom.mass/mass)*atom.pos for atom in self.itervalues())

    def distance_to(self, other):
        """The distance of the mass_center of the residual to `other`
        """
        if len(other) != 3:
            raise ValueError("")
        return numpy.linalg.norm(self.mass_center - other)

    def rmsd(self, ref):
        """The RMSD of the atoms of the residual, with ref as reference.

        Arguments:
            ref (numpy.array): the x,y,z coordinates of the reference
                position.
        """
        return numpy.sqrt(
            (1/len(self))*sum(atom.distance_to(ref) for atom in self)
        )
