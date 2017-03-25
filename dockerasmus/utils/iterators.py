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
    STARTED = False

    if stop is not None and start > stop:
        raise ValueError(
            "'{}' does not come before '{}' in alphabetic order !".format(start, stop)
        )

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
