# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import itertools
import string


def infinitewords(start='A', stop=None):
    """An infinite word iterator in alphabetic order.

    Exemple:
        >>> for word in infinitewords('AA', 'AC'):
        ...     print(word)
        AA
        AB
    """
    uppercase = string.ascii_uppercase
    k = 0
    STARTED = False

    while True:
        k += 1
        for word in (''.join(x) for x in itertools.product(*[uppercase]*k)):
            if word == start:
                STARTED = True
            if word == stop:
                return
            if STARTED:
                yield word


def nth(iterable, n, default=None):
        """Returns the nth item of an interator or a default value
        """
        return next(itertools.islice(iterable, n, None), default)
