# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import collections
import six

from .residual import Residual
from .atom import Atom



class Chain(collections.OrderedDict):
    __slots__ = ("id", "name")

    def __init__(self, id, name=None, residuals=None):
        super(Chain, self).__init__(residuals or [])
        self.id = id
        self.name = name

    def __contains__(self, item):
        if isinstance(item, Residual):
            return super(Chain, self).__contains__(item)
        elif isinstance(item, Atom):
            return any(item in res for res in self)
        elif isinstance(item, int):
            return any(item == res.id for res in self.itervalues())
        else:
            raise TypeError(
                "'in <Chain>' requires Residual, Atom or int"
                " as left operand, not {}".format(type(item).__name__)
            )

    @property
    def mass(self):
        """The mass of the chain

        Computed as the sum of the masses of the residuals
        of the chain (it does not take the masses of the aminous
        acids in the peptidic bound into account).
        """
        return sum(res.mass for res in self.itervalues())

    if six.PY3:
        def itervalues(self):
            return six.itervalues(self)

        def iteritems(self):
            return six.iteritems(self)
