# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import abc
import six

from ... import utils


@six.add_metaclass(abc.ABCMeta)
class BaseComponent(object):
    """The individual component of a scoring function

    For instance, the scoring function used in AutoDock4
    can be decomposed as several independent components:

    .. math::

       V = W_{vdw} V_{vdw}   + W_{hb} V_{hb}
         + W_{elec} V_{elec} + W_{sol} V_{sol}

    where :math:`V_{i}` is an individual scoring component
    (Van der Waals interactions, Hydrogen bonds, electro
    -static forces or desolvation) and :math:`W_{i}`
    the weight associated to that component.

    Warning:
        Components can rely on several different backends, i.e.
        libraries that handle linear algebra computations. Those
        can and will likely handle floating point precision
        differently, resulting in seemingly different results
        from the same inputs. However, Components are tested
        against manually computed cases, and code is verified
        to be the same across all backends, so these errors
        should only occur in Components with rapidly growing
        functions (:math:`exp` and such).
    """

    backends = []


    @classmethod
    def _make_argspec(cls):
        spec = utils.getargspec(cls.__call__)
        if spec.defaults is None:
            cls._args = args = spec.args[1:]
            cls._kwargs = kwargs = []
        else:
            cls._args = args = spec.args[1:-len(spec.defaults)]
            cls._kwargs = kwargs = spec.args[-len(spec.defaults):]
        return args, kwargs

    @classmethod
    def kwargs(cls):
        """The list of __call__ keyword arguments.
        """
        try:
            kwargs = cls._kwargs
        except AttributeError:
            _, kwargs = cls._make_argspec()
        return kwargs

    @classmethod
    def args(cls):
        """The list of __call__ positional arguments.
        """
        try:
            args = cls._args
        except AttributeError:
            args, _ = cls._make_argspec()
        return args

    def __init__(self, force_backend=None):

        # try forcing a specific backend if required
        if force_backend is not None:
            if force_backend not in self.backends:
                raise ValueError("Unknown backend: {}".format(force_backend))
            backend_module = utils.maybe_import(force_backend)
            if backend_module is None:
                raise ValueError("Unavailable backend: {}".format(force_backend))
            else:
                backend = force_backend
                setup_function = getattr(self, "_setup_{}".format(force_backend))
                setup_function(backend_module)
                self.backend = force_backend

        # find the best available backend (first in the backends list)
        else:
            unavailable_backends = []
            for backend in self.backends:
                backend_module = utils.maybe_import(backend)
                # the backend was imported
                if backend_module is not None:
                    setup_function = getattr(self, "_setup_{}".format(backend))
                    setup_function(backend_module)
                    self.backend = backend
                    if unavailable_backends: # if we are not using the best backend
                        warnings.warn("Unavailable backends: {}, using {}".format(
                            ', '.join(unavailable_backends), backend,
                        ), UserWarning)
                    break
                # the backend is unavailable
                else:
                    unavailable_backends.append(backend)

        # raise an error if no available backend was found.
        if backend_module is None:
            raise RuntimeError(
                "Could not find any available backend for {} "
                "among: {}".format(
                    type(self).__name__, ', '.join(self.backends)
                )
            )

    @abc.abstractmethod
    def __call__(self):
        """Compute the score component
        """
        pass
