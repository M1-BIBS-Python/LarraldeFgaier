# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import importlib

from . import decorators
from . import iterators
from . import matrices


__all__ = ["decorators", "iterators", "matrices", "maybe_import"]


def maybe_import(module_name):
    """Import a module if available

    Exemple:
        >>> maybe_import("math")            # Returns the actual module
        <module 'math' from ...>
        >>> maybe_import("nomodulehere")    # Returns None
    """
    try:
        return importlib.import_module(module_name)
    except ImportError:
        return None
