# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import collections
import six

from .residue import Residue
from .atom import Atom



class Chain(collections.OrderedDict):
    __slots__ = ("id", "name")

    def __init__(self, id, name=None, residues=None):
        super(Chain, self).__init__(residues or [])
        self.id = id
        self.name = name

    def __contains__(self, item):
        if isinstance(item, int):
            return super(Chain, self).__contains__(item)
        elif isinstance(item, Atom):
            return any(item in res for res in self.itervalues())
        elif isinstance(item, Residue):
            return any(item == res for res in self.itervalues())
        else:
            raise TypeError(
                "'in <Chain>' requires Residue, Atom or int"
                " as left operand, not {}".format(type(item).__name__)
            )

    @property
    def mass(self):
        """The mass of the chain.

        Warning:
            Computed as the sum of the masses of the residues
            of the chain (it does not take the masses of the atoms
            in the peptidic bound into account).
        """
        return sum(res.mass for res in self.itervalues())

    if six.PY3:
        def itervalues(self):
            return six.itervalues(self)

        def iteritems(self):
            return six.iteritems(self)
