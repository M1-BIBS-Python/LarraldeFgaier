# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import itertools
import string


def wordrange(start='A', stop=None):
    """A word range.

    Exemple:
        >>> for word in wordrange('AA', 'AC'):
        ...     print(word)
        AA
        AB
    """
    STARTED = False
    uppercase = string.ascii_uppercase
    if six.PY2:
        uppercase = uppercase.decode('utf-8')

    if stop is not None and start > stop:
        return

    for k in itertools.count(start=1):
        for word in (''.join(x) for x in itertools.product(*[uppercase]*k)):
            STARTED |= word == start
            if word == stop:
                return
            if STARTED:
                yield word


def nth(iterable, n, default=None):
    """Returns the nth item of an interator or a default value

    Exemple:
        >>> nth(range(10), 2)
        2
        >>> nth(range(10), 11, "oops")
        'oops'
    """
    return next(itertools.islice(iterable, n, None), default)
