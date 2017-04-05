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
        """Instantiate a new `Atom` object.

        Arguments:
            x (`int`): the position of the atom on the x-axis.
            y (`int`): the position of the atom on the y-axis.
            z (`int`): the position of the atom on the z-axis.
            id (`int`): the id of the atom in the protein.
            name (`str`): the name of the atom element ('C',
                'CA', 'O', etc.). Giving a name to the Atom is
                required to access to the `mass` property.
            residual (`Residual`): a reference to the residual
                this atom is part of. Giving a reference to
                the residual of the Atom is required to access
                the `charge`, `epsilon` and `radius` properties.
        """
        self.id = id
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.residual = residual

    def __repr__(self):
        return "Atom {}({}, {}, {})".format(self.id, self.x, self.y, self.z)

    def __eq__(self, other):
        return all(
            getattr(self, attr, None)==getattr(other, attr, None)
                for attr in ('x', 'y', 'z', 'id', 'name')
        )

    @property
    def charge(self):
        """The electrostatic charge of the atom.
        """
        return self._read_from_constants(constants.AMINOACID_CHARGES)

    @property
    def epsilon(self):
        r"""The :math:`\epsilon` of the atom (*empirical*).

        See Also:
            `dockerasmus.score.components.LennardJones`
        """
        return self._read_from_constants(constants.AMINOACID_EPSILON)

    @property
    @method_requires(["name"], "Cannot find atom type !")
    def mass(self):
        """The mass of the atom.
        """
        return constants.ATOMIC_MASSES[self.name[0]]

    @property
    def pos(self):
        """The position of the atom.
        """
        return numpy.array([self.x, self.y, self.z])

    @property
    def radius(self):
        """The optimal Van der Waals radius of the atom (*empirical*).
        """
        return self._read_from_constants(constants.AMINOACID_RADIUS)

    if six.PY3:
        def itervalues(self):
            return six.itervalues(self)

        def iteritems(self):
            return six.iteritems(self)

    def distance_to(self, other):
        """Computes the distance to ``other``.

        Arguments:
            other (`numpy.ndarray`): the position to compute the
                distance to (must be array-like of dimension 3)
        """
        if len(other) != 3:
            raise ValueError
        return numpy.linalg.norm(self.pos - other)

    @method_requires(["name", "residual"], "Cannot find atom residual !")
    def nearest(self, other_atom):
        """The nearest ``other_atom`` in ``self.residual``.

        Arguments:
            other_atom (`str`): the name of the other atom
                to find (can be a generic atom such as 'C', 'O', 'N', etc.
                or a specific atom such as 'CA', 'OH1', etc.).

        Example: Carbon atom nearest to the oxygen atom in an arginine
            >>> oxygen = arginine["O"]
            >>> oxygen.nearest("C")
            Atom 36(13.559, 86.257, 95.222)
            >>> arginine["C"]
            Atom 36(13.559, 86.257, 95.222)
        """
        return min([
            atom for atom in self.residual.itervalues()
                if atom.name.startswith(other_atom)],
            key=lambda atom: self.distance_to(atom.pos)
        )

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
