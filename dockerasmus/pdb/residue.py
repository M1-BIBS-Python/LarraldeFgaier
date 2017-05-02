# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import six
import numpy

from .atom import Atom


class Residue(dict):
    __slots__ = ("id", "_name")

    CTER_ATOMS = frozenset({"OXT"})
    NTER_ATOMS = frozenset({"H1", "H2", "H3"})

    def __init__(self, id, name=None, atoms=None):
        super(Residue, self).__init__(atoms or {})
        self.id = id
        self._name = name

    def __contains__(self, item):
        """Check if `item` is contained in the residue.

        Arguments:
            item: either an atom_id (`int`) or an `Atom` object
                to check if present within the residue.
        """
        if isinstance(item, six.text_type):
            return super(Residue, self).__contains__(item)
        elif isinstance(item, Atom):
            return any(item == atom for atom in self.itervalues())
        elif isinstance(item, int):
            return any(item == atom.id for atom in self.itervalues())

        else:
            raise TypeError(
                "'in <Residue>' requires Atom or {}"
                " as left operand, not {}".format(
                    six.text_type.__name__,type(item).__name__)
            )

    def __hash__(self):
        return hash(hash(frozenset(self)) + hash(self.id) + hash(self.name))

    @property
    def cter(self):
        return not self.CTER_ATOMS.isdisjoint(self)

    @property
    def mass(self):
        """The mass of the residue.

        Computed as sum of the masses of the non-hydrogen
        atoms of the residue.
        """
        return sum(atom.mass for atom in self.itervalues())

    @property
    def mass_center(self):
        """The position of the mass center of the residue.

        Computed as the barycenter of the positions of the
        atoms weighted by their atomic masses.
        """
        mass = self.mass
        return sum((atom.mass/mass)*atom.pos for atom in self.itervalues())

    @property
    def name(self):
        if self._name is not None:
            if self._name == "HIS":
                return "HID" if "HD1" in self else "HIE"
            return self._name

    @property
    def nter(self):
        return not self.NTER_ATOMS.isdisjoint(self)

    if six.PY3:
        def itervalues(self):
            return six.itervalues(self)

        def iteritems(self):
            return six.iteritems(self)

    def distance_to(self, other):
        """The distance of the mass center of the residue to ``other``
        """
        if len(other) != 3:
            raise ValueError("") #TODO
        return numpy.linalg.norm(self.mass_center - other)

    def rmsd(self, ref):
        """The RMSD of the atoms of the residue, with ref as reference.

        Arguments:
            ref (`numpy.ndarray` or `list`): the x,y,z
                coordinates of the reference position.
        """
        return numpy.sqrt(
            (1/len(self))*sum(atom.distance_to(ref) for atom in self)
        )
