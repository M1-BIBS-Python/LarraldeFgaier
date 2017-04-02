# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import functools


def method_requires(attributes, msg):
    """Raise a value error with `msg` if any of `attributes` is ``None``.

    Example:
        >>> class Car(object):
        ...     def __init__(self, wheel=None):
        ...         self.wheel = None
        ...         self.direction = "forward"
        ...     @method_requires(["wheel"], "Can't turn if there's no wheel !")
        ...     def turn(self, direction="right"):
        ...         self.direction = direction
        >>> car = Car()
        >>> car.turn("left")
        Traceback (most recent call last):
            raise ValueError(msg)
        ValueError: Can't turn if there's no wheel !
    """
    def decorator(func):
        @functools.wraps(func)
        def new_func(self, *args, **kwargs):
            for attr in attributes:
                if getattr(self, attr, None) is None:
                    raise ValueError(msg)
            return func(self, *args, **kwargs)
        return new_func
    return decorator
