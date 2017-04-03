# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import itertools
import string
import six


def wordrange(start='A', stop=None):
    """A word range.

    Example:
        >>> for word in wordrange('AA', 'AC'):
        ...     print(word)
        AA
        AB
    """
    # TODO: check that stop can actually be reached
    STARTED = False
    uppercase = string.ascii_uppercase
    if six.PY2:
        uppercase = uppercase.decode('utf-8')

    for k in itertools.count(start=1):
        for word in (''.join(x) for x in itertools.product(*[uppercase]*k)):
            STARTED |= word == start
            if word == stop:
                return
            if STARTED:
                yield word


def nth(iterable, n, default=None):
    """Return the nth item of an interable or a default value.

    Example:
        >>> nth(range(10), 2)
        2
        >>> nth(range(10), 11, "oops")
        'oops'
    """
    return next(itertools.islice(iterable, n, None), default)
