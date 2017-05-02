# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import six
import numpy

from .. import constants
from ..utils.decorators import method_requires


class Atom(object):
    __slots__ = ("id", "name", "x", "y", "z", "residue")

    def __init__(self, x, y, z, id, name=None, residue=None):
        """Instantiate a new `Atom` object.

        Arguments:
            x (`int`): the position of the atom on the x-axis.
            y (`int`): the position of the atom on the y-axis.
            z (`int`): the position of the atom on the z-axis.
            id (`int`): the id of the atom in the protein.
            name (`str`): the name of the atom element ('C',
                'CA', 'O', etc.). Giving a name to the Atom is
                required to access to the `mass` property.
            residue (`Residue`): a reference to the residue
                this atom is part of. Giving a reference to
                the residue of the Atom is required to access
                the `charge`, `epsilon` and `radius` properties.
        """
        self.id = id
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.residue = residue

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
    def pwd(self):
        """The well depth of the Lennard-Jones potential of the atom.

        See Also:
            `dockerasmus.score.components.LennardJones`
        """
        return self._read_from_constants(constants.AMINOACID_POTENTIAL_WELL_DEPTH)

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

    def distance_to(self, other):
        """Computes the distance to ``other``.

        Arguments:
            other (`numpy.ndarray`): the position to compute the
                distance to (must be array-like of dimension 3)

        Raises:
            ValueError: when ``other`` is not of dimension 3.
            TypeError: when ``other`` is not a sequence.
        """
        if len(other) != 3:
            raise ValueError("Invalid position: {}".format(other))
        return numpy.linalg.norm(self.pos - other)

    @method_requires(["name", "residue"], "Cannot find atom residue !")
    def nearest(self, other_atom):
        """The nearest ``other_atom`` in ``self.residue``.

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

        Raises:
            ValueError: when no residue is set or when ``other_atom``
                cannot be found in the residue.
        """
        try:
            return min([
                atom for atom in self.residue.itervalues()
                    if atom.name.startswith(other_atom)
                    and atom != self],
                key=lambda atom: self.distance_to(atom.pos)
            )
        except ValueError:
            err = ValueError(
                "Could not find atom named '{}' in "
                "residue {}".format(other_atom, self.residue)
            )
            six.raise_from(err, None)

    @method_requires(["name", "residue"], "Cannot find atom type !")
    def _read_from_constants(self, table):
        try:
            return table[self.residue.name][self.name]
        except KeyError:
            if self.residue.nter:
                return table["NTER"][self.name]
            elif self.residue.cter:
                return table["CTER"][self.name]
            else:
                err = KeyError("{}, {} ({})".format(self.residue.name, self.name, self.id))
                six.raise_from(err, None)
